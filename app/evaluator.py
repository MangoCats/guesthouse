"""Parametric formula evaluator for database-driven element positioning.

Evaluates JSON formula specs in topological dependency order to compute
element polygons, bounding boxes, and positions.  Extends the existing
_resolve_point_spec / _resolve_dir_spec vocabulary from engine.py with
element-reference resolution, constant-valued distances, and wall/item
polygon construction.

Phase 12a: core evaluator, topological sort, dependency extraction.
"""
import json
import math
from collections import defaultdict, deque

from shared.geometry import seg_vecs, offset_pt, line_isect


# ---------------------------------------------------------------------------
# Exceptions
# ---------------------------------------------------------------------------

class CycleError(Exception):
    """Raised when the dependency graph contains a cycle."""

    def __init__(self, cycle_nodes):
        self.cycle_nodes = cycle_nodes
        super().__init__(f"Circular dependency: {cycle_nodes}")


class FormulaError(Exception):
    """Raised when a formula cannot be evaluated."""


# ---------------------------------------------------------------------------
# Face / polygon helpers (mirrors engine.py conventions)
# ---------------------------------------------------------------------------

_FACE_MAP = {"south": (0, 1), "east": (1, 2), "north": (2, 3), "west": (3, 0)}


def _face_midpoint(poly, face):
    """Midpoint of a named face of a 4-vertex polygon."""
    indices = _FACE_MAP.get(face)
    if not indices or len(poly) < 4:
        return None
    i, j = indices
    a, b = poly[i], poly[j]
    return [(a[0] + b[0]) / 2, (a[1] + b[1]) / 2]


def _face_vertices(poly, face):
    """Return (vertex_i, vertex_j) for a named face."""
    indices = _FACE_MAP.get(face)
    if not indices or len(poly) < 4:
        return None, None
    return poly[indices[0]], poly[indices[1]]


def _bbox_from_poly(poly):
    """Bounding box dict from a polygon (list of [E, N])."""
    es = [p[0] for p in poly]
    ns = [p[1] for p in poly]
    return {"w": min(es), "s": min(ns), "e": max(es), "n": max(ns)}


def _line_poly_intersections(origin, direction, poly):
    """Find intersections of a ray from origin along direction with polygon edges.

    Returns list of [E, N] intersection points sorted by distance from origin.
    """
    results = []
    n = len(poly)
    for i in range(n):
        a = poly[i]
        b = poly[(i + 1) % n]
        # Segment direction
        seg_dx = b[0] - a[0]
        seg_dy = b[1] - a[1]
        det = direction[0] * seg_dy - direction[1] * seg_dx
        if abs(det) < 1e-12:
            continue
        da_x = a[0] - origin[0]
        da_y = a[1] - origin[1]
        t = (da_x * seg_dy - da_y * seg_dx) / det
        u = (da_x * direction[1] - da_y * direction[0]) / det
        if 0 <= u <= 1 and t > 1e-9:
            pt = [origin[0] + t * direction[0], origin[1] + t * direction[1]]
            results.append((t, pt))
    results.sort(key=lambda x: x[0])
    return [pt for _, pt in results]


# ---------------------------------------------------------------------------
# Toilet plan-view polygon (SVG coordinates, from gen_floorplan.py)
# ---------------------------------------------------------------------------
_TOILET_SVG = [
    (-1.905, 0), (-1.905, 2.032), (-0.841, 2.032),
    (-1.078, 2.224), (-1.267, 2.455), (-1.408, 2.719),
    (-1.495, 3.005), (-1.524, 3.302),
    (-1.732, 5.461), (-1.699, 5.799), (-1.600, 6.124),
    (-1.440, 6.423), (-1.225, 6.686), (-0.962, 6.901),
    (-0.663, 7.061), (-0.338, 7.160), (0, 7.193),
    (0.338, 7.160), (0.663, 7.061), (0.962, 6.901),
    (1.225, 6.686), (1.440, 6.423), (1.600, 6.124),
    (1.699, 5.799), (1.732, 5.461),
    (1.524, 3.302), (1.495, 3.005), (1.408, 2.719),
    (1.267, 2.455), (1.078, 2.224), (0.847, 2.035),
    (0.841, 2.032), (1.905, 2.032), (1.905, 0),
]
_SVG_TO_FT = 10.0 / 30.48


# ---------------------------------------------------------------------------
# FormulaEvaluator
# ---------------------------------------------------------------------------

class FormulaEvaluator:
    """Topological DAG evaluator for parametric element positioning.

    Usage:
        ev = FormulaEvaluator(constants, base_points, inner_poly, radii)
        ev.add_formula("IW1", "position", formula_json)
        ev.topo_sort()
        ev.evaluate_all()
        result = ev.elements["IW1"]  # {"poly": [...], "bbox": {...}}
    """

    def __init__(self, constants, base_points, inner_poly, radii):
        """Initialise with base geometry context.

        Args:
            constants: dict name → float value
            base_points: dict name → (E, N) or [E, N]
            inner_poly: list of [E, N] or (E, N) — inner wall polygon
            radii: dict name → float radius values
        """
        self.constants = dict(constants)
        self.points = {k: list(v) for k, v in base_points.items()}
        self.inner_poly = [list(p) for p in inner_poly] if inner_poly else []
        self.radii = dict(radii) if radii else {}

        # Computed element results: name → {"poly": [...], "bbox": {...}}
        self.elements = {}

        # Formula storage: (element_name, param_name) → formula dict
        self.formulas = {}

        # Lock state: (element_name, param_name) → locked_value (parsed JSON)
        self.locked = {}

        # Dependencies: (element_name, param_name) → set of dep keys
        self.deps = {}

        # Topological order
        self.eval_order = []

    def add_formula(self, element_name, param_name, formula_json, locked=False,
                    locked_value=None):
        """Register a formula for evaluation."""
        key = (element_name, param_name)
        if isinstance(formula_json, str):
            formula_json = json.loads(formula_json)
        self.formulas[key] = formula_json
        if locked and locked_value is not None:
            if isinstance(locked_value, str):
                locked_value = json.loads(locked_value)
            self.locked[key] = locked_value
        # Extract dependencies
        self.deps[key] = extract_deps(formula_json)

    def load_formulas_from_db(self, db_path=None, variant=None):
        """Load all formulas from the database."""
        from app.database import get_all_formulas
        rows = get_all_formulas(variant=variant, db_path=db_path)
        for row in rows:
            fj = row["formula_json"]
            if isinstance(fj, str):
                fj = json.loads(fj)
            locked = bool(row.get("locked"))
            lv = row.get("locked_value")
            if lv and isinstance(lv, str):
                lv = json.loads(lv)
            self.add_formula(
                row["element_name"], row["param_name"], fj,
                locked=locked, locked_value=lv if locked else None,
            )

    def topo_sort(self):
        """Topological sort using Kahn's algorithm.  Raises CycleError on cycles."""
        # Build adjacency: dep_key → set of formula keys that depend on it
        in_degree = defaultdict(int)
        adjacency = defaultdict(set)

        # Initialise all formula nodes
        for key in self.formulas:
            if key not in in_degree:
                in_degree[key] = 0

        for key, dep_set in self.deps.items():
            for dep in dep_set:
                # dep is (dep_type, dep_name)
                # Find formula keys that this dep satisfies
                dep_type, dep_name = dep
                if dep_type == "element":
                    # Look for any formula key with this element_name
                    for fkey in self.formulas:
                        if fkey[0] == dep_name and fkey != key:
                            adjacency[fkey].add(key)
                            in_degree[key] += 1
                # Constants and points are leaves — no formula dependency

        queue = deque(k for k in self.formulas if in_degree[k] == 0)
        order = []
        while queue:
            node = queue.popleft()
            order.append(node)
            for dependent in adjacency[node]:
                in_degree[dependent] -= 1
                if in_degree[dependent] == 0:
                    queue.append(dependent)

        if len(order) != len(self.formulas):
            remaining = set(self.formulas) - set(order)
            raise CycleError(remaining)

        self.eval_order = order

    def evaluate_all(self):
        """Evaluate all formulas in topological order."""
        if not self.eval_order:
            self.topo_sort()
        for key in self.eval_order:
            elem_name, param_name = key
            # Check lock
            if key in self.locked:
                self.elements[elem_name] = self.locked[key]
                continue
            formula = self.formulas[key]
            result = self._evaluate_formula(elem_name, formula)
            if result is not None:
                self.elements[elem_name] = result

    def _evaluate_formula(self, elem_name, formula):
        """Evaluate a single formula spec and return result dict."""
        ftype = formula.get("type") if isinstance(formula, dict) else None
        if ftype == "wall_rect":
            return self._eval_wall_rect(formula)
        if ftype == "item_rect":
            return self._eval_item_rect(formula)
        if ftype == "item_circle":
            return self._eval_item_circle(formula)
        if ftype == "wall_opening":
            return self._eval_wall_opening(formula)
        if ftype == "four_corner":
            return self._eval_four_corner(formula)
        if ftype == "toilet_shape":
            return self._eval_toilet_shape(formula)
        if ftype == "bath_sink_shape":
            return self._eval_bath_sink_shape(formula)
        if ftype == "dining_triangle":
            return self._eval_dining_triangle(elem_name, formula)
        if ftype == "dining_chair":
            return self._eval_dining_chair(formula)
        if ftype == "ellipse_rect":
            return self._eval_ellipse_rect(formula)
        return None

    # -------------------------------------------------------------------
    # Formula type evaluators
    # -------------------------------------------------------------------

    def _eval_wall_rect(self, formula):
        """Evaluate a wall_rect formula → {"poly": [...], "bbox": {...}}."""
        anchor = self.resolve_point(formula.get("anchor"))
        along = self.resolve_dir(formula.get("along"))
        thick_dir = self.resolve_dir(formula.get("thickness_dir"))
        thickness = self.resolve_length(formula.get("thickness"))
        if anchor is None or along is None or thick_dir is None or thickness is None:
            return None

        end_mode = formula.get("end_mode", "fixed")

        if end_mode == "intersect":
            # Find wall extent by intersecting with target polygon
            target = formula.get("end_target", "inner_poly")
            poly = self.inner_poly if target == "inner_poly" else []
            if not poly:
                return None

            # Find intersections along the wall direction from anchor
            isects = _line_poly_intersections(anchor, along, poly)
            if not isects:
                return None
            # Take the nearest intersection as the far end
            select = formula.get("select", "nearest")
            far_pt = isects[0] if select == "nearest" else isects[-1]

            # Also check the opposite direction for the near end
            neg_along = [-along[0], -along[1]]
            neg_isects = _line_poly_intersections(anchor, neg_along, poly)

            if neg_isects:
                near_pt = neg_isects[0]
            else:
                near_pt = anchor

            # Build 4-corner rectangle: SW, SE, NE, NW
            # Convention: along direction goes from near to far
            # thick_dir goes from outer face to inner face
            sw = near_pt
            se = far_pt
            ne = [far_pt[0] + thickness * thick_dir[0],
                  far_pt[1] + thickness * thick_dir[1]]
            nw = [near_pt[0] + thickness * thick_dir[0],
                  near_pt[1] + thickness * thick_dir[1]]
            poly_out = [sw, se, ne, nw]
        else:
            # Fixed length
            length = self.resolve_length(formula.get("length"))
            if length is None:
                return None

            sw = anchor
            se = [anchor[0] + length * along[0],
                  anchor[1] + length * along[1]]
            ne = [se[0] + thickness * thick_dir[0],
                  se[1] + thickness * thick_dir[1]]
            nw = [anchor[0] + thickness * thick_dir[0],
                  anchor[1] + thickness * thick_dir[1]]
            poly_out = [sw, se, ne, nw]

        return {"poly": poly_out, "bbox": _bbox_from_poly(poly_out)}

    def _eval_item_rect(self, formula):
        """Evaluate an item_rect formula → {"poly": [...], "bbox": {...}}."""
        anchor = self.resolve_point(formula.get("anchor"))
        along = self.resolve_dir(formula.get("along"))
        across = self.resolve_dir(formula.get("across"))
        width = self.resolve_length(formula.get("width"))
        depth = self.resolve_length(formula.get("depth"))
        if any(v is None for v in (anchor, along, across, width, depth)):
            return None

        corner = formula.get("anchor_corner", "sw")
        # Compute SW corner from anchor + corner type
        if corner == "sw":
            sw = anchor
        elif corner == "se":
            sw = [anchor[0] - width * along[0], anchor[1] - width * along[1]]
        elif corner == "nw":
            sw = [anchor[0] - depth * across[0], anchor[1] - depth * across[1]]
        elif corner == "ne":
            sw = [anchor[0] - width * along[0] - depth * across[0],
                  anchor[1] - width * along[1] - depth * across[1]]
        else:
            sw = anchor

        se = [sw[0] + width * along[0], sw[1] + width * along[1]]
        ne = [se[0] + depth * across[0], se[1] + depth * across[1]]
        nw = [sw[0] + depth * across[0], sw[1] + depth * across[1]]
        poly = [sw, se, ne, nw]
        return {"poly": poly, "bbox": _bbox_from_poly(poly)}

    def _eval_item_circle(self, formula):
        """Evaluate an item_circle formula → {"center": [...], "radius": R, "poly": [...], "bbox": {...}}."""
        center = self.resolve_point(formula.get("center"))
        radius = self.resolve_length(formula.get("radius"))
        if center is None or radius is None:
            return None
        # Approximate circle as polygon (24 segments)
        n_pts = formula.get("n_pts", 24)
        poly = []
        for i in range(n_pts):
            angle = 2 * math.pi * i / n_pts
            poly.append([center[0] + radius * math.cos(angle),
                         center[1] + radius * math.sin(angle)])
        return {
            "center": center,
            "radius": radius,
            "poly": poly,
            "bbox": _bbox_from_poly(poly),
        }

    def _eval_four_corner(self, formula):
        """Evaluate a four_corner formula → {"poly": [...], "bbox": {...}}.

        Each corner is an independent point spec:
          {"type": "four_corner", "sw": spec, "se": spec, "ne": spec, "nw": spec}
        """
        sw = self.resolve_point(formula.get("sw"))
        se = self.resolve_point(formula.get("se"))
        ne = self.resolve_point(formula.get("ne"))
        nw = self.resolve_point(formula.get("nw"))
        if any(v is None for v in (sw, se, ne, nw)):
            return None
        poly = [sw, se, ne, nw]
        return {"poly": poly, "bbox": _bbox_from_poly(poly)}

    def _eval_wall_opening(self, formula):
        """Evaluate a wall_opening formula → {"poly": [...], "bbox": {...}}.

        Computes a 4-point polygon spanning from outer to inner wall face.
        Parametric t-values are computed on the outer segment and applied to
        both outer and inner segments.

        Positioning modes:

        1. Gap from endpoint (default):
           "from_end": true, "gap": dist_from_end, "width": opening_width
           "from_end": false, "gap": dist_from_start, "width": opening_width

        2. Reference point projection:
           "ref_point": point_spec, "ref_offset": length_spec, "width": w
           Projects ref_point onto outer segment, offsets by ref_offset,
           then places opening of given width starting from that position.

        3. Centered:
           "centered": true, "center_t": float (0-1), "width": w
           Centers opening at the given parametric position (0.5 = midpoint).

        4. Centered between two projections:
           "center_refs": [point_spec, point_spec], "width": w
           Projects both points onto the segment, centers opening between them.

        Result poly = [outer_start, outer_end, inner_end, inner_start].
        """
        os = self.resolve_point(formula.get("outer_start"))
        oe = self.resolve_point(formula.get("outer_end"))
        ins = self.resolve_point(formula.get("inner_start"))
        ine = self.resolve_point(formula.get("inner_end"))
        width = self.resolve_length(formula.get("width"))
        if any(v is None for v in (os, oe, ins, ine, width)):
            return None

        # Outer segment vector
        dx, dy = oe[0] - os[0], oe[1] - os[1]
        seg_len = math.sqrt(dx * dx + dy * dy)
        if seg_len < 1e-12:
            return None

        # Determine t_start and t_end (parametric on outer segment)
        if "ref_point" in formula:
            # Mode 2: reference point projection + offset
            ref_pt = self.resolve_point(formula["ref_point"])
            ref_offset = self.resolve_length(formula.get("ref_offset", 0))
            if ref_pt is None or ref_offset is None:
                return None
            # Project ref_point onto outer segment
            rx, ry = ref_pt[0] - os[0], ref_pt[1] - os[1]
            t_ref = (rx * dx + ry * dy) / (dx * dx + dy * dy)
            t_start = t_ref + ref_offset / seg_len
            t_end = t_start + width / seg_len
        elif "center_refs" in formula:
            # Mode 4: centered between two projected points
            refs = formula["center_refs"]
            p1 = self.resolve_point(refs[0])
            p2 = self.resolve_point(refs[1])
            if p1 is None or p2 is None:
                return None
            t1 = ((p1[0] - os[0]) * dx + (p1[1] - os[1]) * dy) / (dx * dx + dy * dy)
            t2 = ((p2[0] - os[0]) * dx + (p2[1] - os[1]) * dy) / (dx * dx + dy * dy)
            t_center = (t1 + t2) / 2
            half_w = width / (2 * seg_len)
            t_start = t_center - half_w
            t_end = t_center + half_w
        elif formula.get("centered"):
            # Mode 3: centered at fixed t
            t_center = formula.get("center_t", 0.5)
            half_w = width / (2 * seg_len)
            t_start = t_center - half_w
            t_end = t_center + half_w
        else:
            # Mode 1: gap from endpoint
            gap = self.resolve_length(formula.get("gap", 0))
            if gap is None:
                return None
            from_end = formula.get("from_end", True)
            if from_end:
                t_end = 1.0 - gap / seg_len
                t_start = t_end - width / seg_len
            else:
                t_start = gap / seg_len
                t_end = t_start + width / seg_len

        # Interpolate on outer segment
        outer_start = [os[0] + t_start * dx, os[1] + t_start * dy]
        outer_end = [os[0] + t_end * dx, os[1] + t_end * dy]

        # Interpolate on inner segment (same t values)
        idx, idy = ine[0] - ins[0], ine[1] - ins[1]
        inner_end = [ins[0] + t_end * idx, ins[1] + t_end * idy]
        inner_start = [ins[0] + t_start * idx, ins[1] + t_start * idy]

        # Default poly order: [outer_start, outer_end, inner_end, inner_start]
        # Other orderings supported via "poly_order" field:
        #   "inner_first": [inner_start, inner_end, outer_end, outer_start]
        #   "outer_reversed": [outer_end, outer_start, inner_start, inner_end]
        #   "face_pair": [outer_start, inner_start, inner_end, outer_end]
        order = formula.get("poly_order", "outer_first")
        if order == "inner_first":
            poly = [inner_start, inner_end, outer_end, outer_start]
        elif order == "outer_reversed":
            poly = [outer_end, outer_start, inner_start, inner_end]
        elif order == "face_pair":
            poly = [outer_start, inner_start, inner_end, outer_end]
        else:
            poly = [outer_start, outer_end, inner_end, inner_start]
        return {"poly": poly, "bbox": _bbox_from_poly(poly)}

    def _eval_toilet_shape(self, formula):
        """Evaluate a toilet_shape formula → {"poly": [...], "bbox": {...}}.

        Transforms the fixed toilet SVG polygon to building coordinates.
        """
        center = self.resolve_point(formula.get("center"))
        facing = self.resolve_dir(formula.get("facing_dir"))
        width_dir = self.resolve_dir(formula.get("width_dir"))
        if any(v is None for v in (center, facing, width_dir)):
            return None
        poly = [(center[0] + dx * _SVG_TO_FT * width_dir[0]
                            + dy * _SVG_TO_FT * facing[0],
                 center[1] + dx * _SVG_TO_FT * width_dir[1]
                            + dy * _SVG_TO_FT * facing[1])
                for dx, dy in _TOILET_SVG]
        poly = [[p[0], p[1]] for p in poly]
        return {"poly": poly, "bbox": _bbox_from_poly(poly)}

    def _eval_bath_sink_shape(self, formula):
        """Evaluate a bath_sink_shape formula → {"poly": [...], "bbox": {...}}.

        Semicircular bulge polygon: rectangular base with arc on outward face.
        """
        anchor = self.resolve_point(formula.get("anchor"))
        al = self.resolve_dir(formula.get("along"))
        out = self.resolve_dir(formula.get("outward"))
        length = self.resolve_length(formula.get("length"))
        depth = self.resolve_length(formula.get("depth"))
        if any(v is None for v in (anchor, al, out, length, depth)):
            return None

        half_len = length / 2
        quarter_len = length / 4
        rect_depth = depth - quarter_len
        n_arc = formula.get("n_arc", 32)

        def _pt(along, outward):
            return [anchor[0] + along * al[0] + outward * out[0],
                    anchor[1] + along * al[1] + outward * out[1]]

        pts = []
        pts.append(_pt(half_len, 0))
        pts.append(_pt(-half_len, 0))
        pts.append(_pt(-half_len, rect_depth))
        for i in range(n_arc + 1):
            t = math.pi - math.pi * i / n_arc
            pts.append(_pt(math.cos(t) * quarter_len,
                           rect_depth + math.sin(t) * quarter_len))
        pts.append(_pt(half_len, rect_depth))
        return {"poly": pts, "bbox": _bbox_from_poly(pts)}

    def _eval_dining_triangle(self, elem_name, formula):
        """Evaluate a dining_triangle formula → {"poly": [...], "bbox": {...}}.

        Triangle table with apex arc and fillet corners.  Stores side endpoint
        data for dining chair positioning.
        """
        bc = self.resolve_point(formula.get("base_center"))
        to_apex = self.resolve_dir(formula.get("toward_apex"))
        along = self.resolve_dir(formula.get("along_base"))
        base_w = self.resolve_length(formula.get("base_width"))
        height = self.resolve_length(formula.get("height"))
        apex_r = self.resolve_length(formula.get("apex_radius"))
        fillet_r = self.resolve_length(formula.get("fillet_radius"))
        if any(v is None for v in (bc, to_apex, along, base_w, height,
                                    apex_r, fillet_r)):
            return None

        to_base = [-to_apex[0], -to_apex[1]]
        ne = [bc[0] + base_w / 2 * along[0], bc[1] + base_w / 2 * along[1]]
        nw = [bc[0] - base_w / 2 * along[0], bc[1] - base_w / 2 * along[1]]
        apex_tip = [bc[0] + height * to_apex[0], bc[1] + height * to_apex[1]]
        arc_c = [apex_tip[0] + apex_r * to_base[0],
                 apex_tip[1] + apex_r * to_base[1]]

        def _sym(p):
            vx, vy = p[0] - bc[0], p[1] - bc[1]
            v_dot = vx * to_apex[0] + vy * to_apex[1]
            return [bc[0] + 2 * v_dot * to_apex[0] - vx,
                    bc[1] + 2 * v_dot * to_apex[1] - vy]

        d_base_ne = [-along[0], -along[1]]
        dx_r = ne[0] - arc_c[0]
        dn_r = ne[1] - arc_c[1]
        dist_r = math.sqrt(dx_r**2 + dn_r**2)
        angle_cp = math.atan2(dn_r, dx_r)
        delta = math.acos(max(-1, min(1, apex_r / dist_r)))
        alpha_r = angle_cp - delta
        t_right = [arc_c[0] + apex_r * math.cos(alpha_r),
                    arc_c[1] + apex_r * math.sin(alpha_r)]
        t_left = _sym(t_right)

        dtr = [t_right[0] - ne[0], t_right[1] - ne[1]]
        dtr_len = math.sqrt(dtr[0]**2 + dtr[1]**2)
        d_tang_ne = [dtr[0] / dtr_len, dtr[1] / dtr_len]
        cos_th = d_base_ne[0] * d_tang_ne[0] + d_base_ne[1] * d_tang_ne[1]
        half_angle = math.acos(max(-1, min(1, cos_th))) / 2
        fillet_dist = fillet_r / math.sin(half_angle)
        bis_ne = [d_base_ne[0] + d_tang_ne[0], d_base_ne[1] + d_tang_ne[1]]
        bis_ne_len = math.sqrt(bis_ne[0]**2 + bis_ne[1]**2)
        bis_ne = [bis_ne[0] / bis_ne_len, bis_ne[1] / bis_ne_len]
        fc_ne = [ne[0] + fillet_dist * bis_ne[0],
                 ne[1] + fillet_dist * bis_ne[1]]
        v_ne = [fc_ne[0] - ne[0], fc_ne[1] - ne[1]]
        t_proj = v_ne[0] * d_tang_ne[0] + v_ne[1] * d_tang_ne[1]
        f_ne_tang = [ne[0] + t_proj * d_tang_ne[0],
                     ne[1] + t_proj * d_tang_ne[1]]
        f_nw_tang = _sym(f_ne_tang)

        fc_ne_d = ((fc_ne[0] - bc[0]) * along[0] +
                    (fc_ne[1] - bc[1]) * along[1])
        f_ne_base = [bc[0] + fc_ne_d * along[0], bc[1] + fc_ne_d * along[1]]
        f_nw_base = _sym(f_ne_base)
        fc_nw = _sym(fc_ne)

        ne_fil_a0 = math.atan2(f_ne_base[1] - fc_ne[1],
                                f_ne_base[0] - fc_ne[0])
        ne_fil_a1 = math.atan2(f_ne_tang[1] - fc_ne[1],
                                f_ne_tang[0] - fc_ne[0])
        apex_a0 = math.atan2(t_right[1] - arc_c[1], t_right[0] - arc_c[0])
        apex_a1 = math.atan2(t_left[1] - arc_c[1], t_left[0] - arc_c[0])
        nw_fil_a0 = math.atan2(f_nw_tang[1] - fc_nw[1],
                                f_nw_tang[0] - fc_nw[0])
        nw_fil_a1 = math.atan2(f_nw_base[1] - fc_nw[1],
                                f_nw_base[0] - fc_nw[0])

        def _norm_sweep(a0, a1):
            d = (a1 - a0) % (2 * math.pi)
            if d > math.pi:
                d -= 2 * math.pi
            return d

        def _arc_pts(center, radius, a0, sweep, n=8):
            pts = []
            for i in range(1, n):
                t = a0 + sweep * i / n
                pts.append([center[0] + radius * math.cos(t),
                            center[1] + radius * math.sin(t)])
            return pts

        poly = [f_nw_base, f_ne_base]
        sweep = _norm_sweep(ne_fil_a0, ne_fil_a1)
        poly += _arc_pts(fc_ne, fillet_r, ne_fil_a0, sweep, 8)
        poly.append(f_ne_tang)
        poly.append(t_right)
        sweep = _norm_sweep(apex_a0, apex_a1)
        poly += _arc_pts(arc_c, apex_r, apex_a0, sweep, 16)
        poly.append(t_left)
        poly.append(f_nw_tang)
        sweep = _norm_sweep(nw_fil_a0, nw_fil_a1)
        poly += _arc_pts(fc_nw, fillet_r, nw_fil_a0, sweep, 8)

        result = {"poly": poly, "bbox": _bbox_from_poly(poly)}
        # Store side endpoints for dining chair positioning
        result["ne_side"] = [f_ne_tang, t_right]
        result["nw_side"] = [t_left, f_nw_tang]
        result["base_center"] = bc
        return result

    def _eval_dining_chair(self, formula):
        """Evaluate a dining_chair formula → {"poly": [...], "bbox": {...}}.

        Positioned at midpoint of a dining table side, facing outward.
        """
        table_name = formula.get("table")
        side_key = formula.get("side")  # "ne_side" or "nw_side"
        ch_short = self.resolve_length(formula.get("chair_width"))
        ch_long = self.resolve_length(formula.get("chair_depth"))
        gap = self.resolve_length(formula.get("gap", 2.0 / 12.0))
        if table_name is None or side_key is None or ch_short is None or ch_long is None:
            return None
        table = self.elements.get(table_name)
        if not table or side_key not in table:
            return None

        side_start, side_end = table[side_key]
        bc = table.get("base_center")
        if bc is None:
            return None

        mid_e = (side_start[0] + side_end[0]) / 2
        mid_n = (side_start[1] + side_end[1]) / 2
        se_d = (side_end[0] - side_start[0], side_end[1] - side_start[1])
        sl = math.sqrt(se_d[0]**2 + se_d[1]**2)
        su = (se_d[0] / sl, se_d[1] / sl)
        sn = (-su[1], su[0])
        to_ctr = (bc[0] - mid_e, bc[1] - mid_n)
        if sn[0] * to_ctr[0] + sn[1] * to_ctr[1] > 0:
            sn = (-sn[0], -sn[1])
        cc_e = mid_e + sn[0] * (ch_long / 2 + gap)
        cc_n = mid_n + sn[1] * (ch_long / 2 + gap)
        corners = []
        for ds, dn in [(-1, -1), (1, -1), (1, 1), (-1, 1)]:
            ce = cc_e + su[0] * ds * ch_short / 2 + sn[0] * dn * ch_long / 2
            cn = cc_n + su[1] * ds * ch_short / 2 + sn[1] * dn * ch_long / 2
            corners.append([ce, cn])
        return {"poly": corners, "bbox": _bbox_from_poly(corners)}

    def _eval_ellipse_rect(self, formula):
        """Evaluate an ellipse_rect formula → bounding rect of an ellipse.

        anchor: center on wall face
        along: direction along wall
        outward: direction away from wall
        rx: half-width along wall
        ry: half-depth away from wall
        """
        anchor = self.resolve_point(formula.get("anchor"))
        al = self.resolve_dir(formula.get("along"))
        out = self.resolve_dir(formula.get("outward"))
        rx = self.resolve_length(formula.get("rx"))
        ry = self.resolve_length(formula.get("ry"))
        if any(v is None for v in (anchor, al, out, rx, ry)):
            return None

        sw = [anchor[0] - rx * al[0] - ry * out[0],
              anchor[1] - rx * al[1] - ry * out[1]]
        se = [anchor[0] + rx * al[0] - ry * out[0],
              anchor[1] + rx * al[1] - ry * out[1]]
        ne = [anchor[0] + rx * al[0] + ry * out[0],
              anchor[1] + rx * al[1] + ry * out[1]]
        nw = [anchor[0] - rx * al[0] + ry * out[0],
              anchor[1] - rx * al[1] + ry * out[1]]
        poly = [sw, se, ne, nw]
        return {"poly": poly, "bbox": _bbox_from_poly(poly)}

    # -------------------------------------------------------------------
    # Resolution functions
    # -------------------------------------------------------------------

    def resolve_point(self, spec):
        """Resolve a point spec to [E, N] coordinates.

        Extends engine.py's _resolve_point_spec with:
          - {"element": name, "corner": idx_or_face} → element polygon corner
          - {"line_poly_isect": {...}} → line-polygon intersection
          - constant-valued dist in offset specs
        """
        if spec is None:
            return None

        # List/tuple passthrough: already [E, N]
        if isinstance(spec, (list, tuple)) and len(spec) == 2:
            return [float(spec[0]), float(spec[1])]

        # String: named point
        if isinstance(spec, str):
            pt = self.points.get(spec)
            return list(pt) if pt else None

        if not isinstance(spec, dict):
            return None

        # Element corner reference
        if "element" in spec:
            elem_name = spec["element"]
            elem = self.elements.get(elem_name)
            if not elem or "poly" not in elem:
                return None
            corner = spec.get("corner")
            poly = elem["poly"]
            if isinstance(corner, int):
                if 0 <= corner < len(poly):
                    return list(poly[corner])
                return None
            if isinstance(corner, str):
                # Named corner for 4-vertex polygons
                corner_map = {"SW": 0, "SE": 1, "NE": 2, "NW": 3,
                              "sw": 0, "se": 1, "ne": 2, "nw": 3}
                idx = corner_map.get(corner)
                if idx is not None and idx < len(poly):
                    return list(poly[idx])
                return None
            return None

        # Element centroid
        if "element_centroid" in spec:
            name = spec["element_centroid"]
            elem = self.elements.get(name)
            if elem and "poly" in elem:
                poly = elem["poly"]
                n = len(poly)
                if n == 0:
                    return None
                cx = sum(p[0] for p in poly) / n
                cy = sum(p[1] for p in poly) / n
                return [cx, cy]
            return None

        # Face midpoint of an element
        if "face_mid" in spec:
            name = spec["face_mid"]
            face = spec.get("face")
            # Try evaluated elements first, then fall back to base context
            elem = self.elements.get(name)
            if elem and "poly" in elem:
                return _face_midpoint(elem["poly"], face)
            return None

        # Midpoint of two specs
        if "midpoint" in spec:
            pair = spec["midpoint"]
            if len(pair) != 2:
                return None
            a = self.resolve_point(pair[0])
            b = self.resolve_point(pair[1])
            if a and b:
                return [(a[0] + b[0]) / 2, (a[1] + b[1]) / 2]
            return None

        # Offset: base + dist * direction
        if "offset" in spec:
            base = self.resolve_point(spec["offset"])
            direction = self.resolve_dir(spec.get("dir"))
            dist = self.resolve_length(spec.get("dist", 0))
            if base and direction and dist is not None:
                return [base[0] + dist * direction[0],
                        base[1] + dist * direction[1]]
            return None

        # Line intersection
        if spec.get("type") == "line_intersection":
            p1 = self.resolve_point(spec.get("line1_point"))
            d1 = self.resolve_dir(spec.get("line1_dir"))
            p2 = self.resolve_point(spec.get("line2_point"))
            d2 = self.resolve_dir(spec.get("line2_dir"))
            if p1 and d1 and p2 and d2:
                det = d1[0] * d2[1] - d1[1] * d2[0]
                if abs(det) < 1e-12:
                    return None
                t = ((p2[0] - p1[0]) * d2[1] - (p2[1] - p1[1]) * d2[0]) / det
                return [p1[0] + t * d1[0], p1[1] + t * d1[1]]
            return None

        # Line-polygon intersection
        if "line_poly_isect" in spec:
            lpi = spec["line_poly_isect"]
            origin = self.resolve_point(lpi.get("origin"))
            direction = self.resolve_dir(lpi.get("dir"))
            poly_name = lpi.get("poly", "inner_poly")
            poly = self.inner_poly if poly_name == "inner_poly" else []
            select = lpi.get("select", "nearest")
            if origin and direction and poly:
                isects = _line_poly_intersections(origin, direction, poly)
                if isects:
                    return isects[0] if select == "nearest" else isects[-1]
            return None

        # Arc point (mirror of engine.py)
        if "arc_point" in spec:
            ap = spec["arc_point"]
            center = self.resolve_point(ap.get("center"))
            radius = self.radii.get(ap.get("radius_key"))
            ref = self.resolve_point(ap.get("reference"))
            side = ap.get("side", "east")
            if center and radius and ref:
                dn = ref[1] - center[1]
                disc = radius**2 - dn**2
                if disc < 0:
                    return None
                de = math.sqrt(disc)
                if side == "west":
                    de = -de
                return [center[0] + de, center[1] + dn]
            return None

        # Ray-circle intersection: find point on ray at given distance from target
        # {"ray_circle_isect": {"origin": pt, "dir": dir, "target": pt, "distance": len}}
        if "ray_circle_isect" in spec:
            rci = spec["ray_circle_isect"]
            origin = self.resolve_point(rci.get("origin"))
            direction = self.resolve_dir(rci.get("dir"))
            target = self.resolve_point(rci.get("target"))
            distance = self.resolve_length(rci.get("distance"))
            if origin and direction and target and distance is not None:
                # Solve |origin + t*dir - target|² = distance²
                dx = origin[0] - target[0]
                dy = origin[1] - target[1]
                d_al = dx * direction[0] + dy * direction[1]
                d2 = dx**2 + dy**2
                disc = distance**2 - d2 + d_al**2
                if disc < 0:
                    return None
                select = rci.get("select", "farthest")
                if select == "nearest":
                    t = -d_al - math.sqrt(disc)
                else:
                    t = -d_al + math.sqrt(disc)
                return [origin[0] + t * direction[0],
                        origin[1] + t * direction[1]]
            return None

        return None

    def resolve_dir(self, spec):
        """Resolve a direction spec to [dE, dN] unit vector."""
        if spec is None:
            return None

        if isinstance(spec, (list, tuple)) and len(spec) == 2:
            return [float(spec[0]), float(spec[1])]

        if isinstance(spec, str):
            if spec == "east":
                return [1.0, 0.0]
            if spec == "west":
                return [-1.0, 0.0]
            if spec == "north":
                return [0.0, 1.0]
            if spec == "south":
                return [0.0, -1.0]
            return None

        if not isinstance(spec, dict):
            return None

        # Face along direction
        if "face_along" in spec:
            name = spec["face_along"]
            face = spec.get("face")
            elem = self.elements.get(name)
            if elem and "poly" in elem:
                a, b = _face_vertices(elem["poly"], face)
                if a is not None:
                    dx, dy = b[0] - a[0], b[1] - a[1]
                    length = math.sqrt(dx * dx + dy * dy)
                    if length < 1e-12:
                        return None
                    return [dx / length, dy / length]
            return None

        # Face perpendicular
        if "face_perp" in spec:
            name = spec["face_perp"]
            face = spec.get("face")
            elem = self.elements.get(name)
            if elem and "poly" in elem:
                a, b = _face_vertices(elem["poly"], face)
                if a is not None:
                    dx, dy = b[0] - a[0], b[1] - a[1]
                    length = math.sqrt(dx * dx + dy * dy)
                    if length < 1e-12:
                        return None
                    return [dy / length, -dx / length]
            return None

        # Segment direction (between two points)
        if "segment" in spec:
            pair = spec["segment"]
            if len(pair) != 2:
                return None
            a = self.resolve_point(pair[0])
            b = self.resolve_point(pair[1])
            if a and b:
                dx, dy = b[0] - a[0], b[1] - a[1]
                length = math.sqrt(dx * dx + dy * dy)
                if length < 1e-12:
                    return None
                return [dx / length, dy / length]
            return None

        # Segment perpendicular
        if "segment_perp" in spec:
            pair = spec["segment_perp"]
            if len(pair) != 2:
                return None
            a = self.resolve_point(pair[0])
            b = self.resolve_point(pair[1])
            if a and b:
                dx, dy = b[0] - a[0], b[1] - a[1]
                length = math.sqrt(dx * dx + dy * dy)
                if length < 1e-12:
                    return None
                return [dy / length, -dx / length]
            return None

        # Negate direction
        if "neg" in spec:
            d = self.resolve_dir(spec["neg"])
            if d:
                return [-d[0], -d[1]]
            return None

        # Rotate direction by angle (radians)
        if "rotated" in spec:
            d = self.resolve_dir(spec["rotated"])
            angle = spec.get("angle", 0)
            if d and isinstance(angle, (int, float)):
                cos_a = math.cos(angle)
                sin_a = math.sin(angle)
                return [d[0] * cos_a - d[1] * sin_a,
                        d[0] * sin_a + d[1] * cos_a]
            return None

        # Perpendicular of direction (rotate 90° CW: [dy, -dx])
        if "perp" in spec:
            d = self.resolve_dir(spec["perp"])
            if d:
                return [d[1], -d[0]]
            return None

        return None

    def resolve_length(self, spec):
        """Resolve a length expression to a float value.

        Length specs:
          float/int             → literal value
          {"const": name}       → constant lookup
          {"neg": spec}         → negation
          {"add": [a, b, ...]}  → sum
          {"sub": [a, b]}       → a - b
          {"mul": [a, b]}       → product
          {"dist": [pt1, pt2]}  → distance between two points
        """
        if spec is None:
            return None
        if isinstance(spec, (int, float)):
            return float(spec)
        if isinstance(spec, dict):
            if "const" in spec:
                val = self.constants.get(spec["const"])
                return float(val) if val is not None else None
            if "radius_key" in spec:
                val = self.radii.get(spec["radius_key"])
                return float(val) if val is not None else None
            if "neg" in spec:
                val = self.resolve_length(spec["neg"])
                return -val if val is not None else None
            if "add" in spec:
                vals = [self.resolve_length(s) for s in spec["add"]]
                if any(v is None for v in vals):
                    return None
                return sum(vals)
            if "sub" in spec:
                pair = spec["sub"]
                if len(pair) != 2:
                    return None
                a = self.resolve_length(pair[0])
                b = self.resolve_length(pair[1])
                if a is None or b is None:
                    return None
                return a - b
            if "mul" in spec:
                vals = [self.resolve_length(s) for s in spec["mul"]]
                if any(v is None for v in vals):
                    return None
                result = 1.0
                for v in vals:
                    result *= v
                return result
            if "proj" in spec:
                # Signed projection of (target - anchor) onto direction
                proj = spec["proj"]
                target = self.resolve_point(proj.get("target"))
                anchor = self.resolve_point(proj.get("anchor"))
                direction = self.resolve_dir(proj.get("dir"))
                if target and anchor and direction:
                    dx = target[0] - anchor[0]
                    dy = target[1] - anchor[1]
                    return dx * direction[0] + dy * direction[1]
                return None
            if "dist" in spec:
                pair = spec["dist"]
                if len(pair) != 2:
                    return None
                a = self.resolve_point(pair[0])
                b = self.resolve_point(pair[1])
                if a is None or b is None:
                    return None
                return math.sqrt((b[0] - a[0])**2 + (b[1] - a[1])**2)
        return None

    # -------------------------------------------------------------------
    # Dependency queries
    # -------------------------------------------------------------------

    def get_dependents(self, elem_name):
        """Return set of element names that depend on elem_name (transitive)."""
        direct = set()
        for key, dep_set in self.deps.items():
            for dep_type, dep_name in dep_set:
                if dep_type == "element" and dep_name == elem_name:
                    direct.add(key[0])
        # Transitive closure
        result = set()
        queue = deque(direct)
        while queue:
            name = queue.popleft()
            if name in result:
                continue
            result.add(name)
            # Find what depends on this
            for key, dep_set in self.deps.items():
                for dt, dn in dep_set:
                    if dt == "element" and dn == name and key[0] not in result:
                        queue.append(key[0])
        return result

    def get_dependencies(self, elem_name):
        """Return set of element names that elem_name depends on (transitive)."""
        direct = set()
        for key in self.formulas:
            if key[0] == elem_name:
                for dep_type, dep_name in self.deps.get(key, set()):
                    if dep_type == "element":
                        direct.add(dep_name)
        # Transitive closure
        result = set()
        queue = deque(direct)
        while queue:
            name = queue.popleft()
            if name in result:
                continue
            result.add(name)
            for key in self.formulas:
                if key[0] == name:
                    for dt, dn in self.deps.get(key, set()):
                        if dt == "element" and dn not in result:
                            queue.append(dn)
        return result


# ---------------------------------------------------------------------------
# Dependency extraction (static analysis of formula JSON)
# ---------------------------------------------------------------------------

def extract_deps(formula):
    """Extract dependency set from a formula JSON spec.

    Returns set of (dep_type, dep_name) tuples where dep_type is
    'element', 'point', or 'constant'.
    """
    deps = set()
    _extract_deps_recursive(formula, deps)
    return deps


def _extract_deps_recursive(spec, deps):
    """Recursively walk a formula spec and collect dependencies."""
    if spec is None:
        return
    if isinstance(spec, (int, float, bool)):
        return
    if isinstance(spec, str):
        # A bare string is a named point reference
        deps.add(("point", spec))
        return
    if isinstance(spec, (list, tuple)):
        for item in spec:
            _extract_deps_recursive(item, deps)
        return
    if not isinstance(spec, dict):
        return

    # Element reference
    if "element" in spec:
        deps.add(("element", spec["element"]))

    # Element centroid
    if "element_centroid" in spec:
        deps.add(("element", spec["element_centroid"]))

    # Dining chair table reference
    if "table" in spec:
        deps.add(("element", spec["table"]))

    # Face midpoint of element
    if "face_mid" in spec:
        deps.add(("element", spec["face_mid"]))

    # Face along/perp of element
    if "face_along" in spec:
        deps.add(("element", spec["face_along"]))
    if "face_perp" in spec:
        deps.add(("element", spec["face_perp"]))

    # Constant reference
    if "const" in spec:
        deps.add(("constant", spec["const"]))

    # Recurse into known sub-specs
    for key in ("anchor", "along", "across", "thickness_dir", "dir", "dist",
                "offset", "center", "reference", "line1_point", "line1_dir",
                "line2_point", "line2_dir", "thickness", "length", "width",
                "depth", "radius", "origin", "neg", "perp", "rotated",
                "sw", "se", "ne", "nw",
                "outer_start", "outer_end", "inner_start", "inner_end",
                "gap", "ref_point", "ref_offset",
                "facing_dir", "width_dir", "outward", "rx", "ry",
                "base_center", "toward_apex", "along_base", "base_width",
                "height", "apex_radius", "fillet_radius",
                "chair_width", "chair_depth", "distance", "target"):
        if key in spec:
            _extract_deps_recursive(spec[key], deps)

    # Length arithmetic
    for key in ("add", "sub", "mul"):
        if key in spec:
            _extract_deps_recursive(spec[key], deps)

    # Projection
    if "proj" in spec:
        _extract_deps_recursive(spec["proj"], deps)

    # Distance between points
    if "dist" in spec and isinstance(spec["dist"], (list, tuple)):
        _extract_deps_recursive(spec["dist"], deps)

    # Midpoint pair
    if "midpoint" in spec:
        _extract_deps_recursive(spec["midpoint"], deps)

    # Segment pair
    if "segment" in spec:
        _extract_deps_recursive(spec["segment"], deps)
    if "segment_perp" in spec:
        _extract_deps_recursive(spec["segment_perp"], deps)

    # Arc point sub-dict
    if "arc_point" in spec:
        _extract_deps_recursive(spec["arc_point"], deps)

    # Line-polygon intersection
    if "line_poly_isect" in spec:
        _extract_deps_recursive(spec["line_poly_isect"], deps)

    # Ray-circle intersection
    if "ray_circle_isect" in spec:
        _extract_deps_recursive(spec["ray_circle_isect"], deps)

    # Center refs (pair of point specs)
    if "center_refs" in spec:
        _extract_deps_recursive(spec["center_refs"], deps)


# ---------------------------------------------------------------------------
# IW formula definitions (Phase 12c)
# ---------------------------------------------------------------------------
# Direction shorthands for formula specs
_W18W1_AL = {"segment": ["W18", "W1"]}      # along W18→W1 (≈ west)
_W18W1_IN = {"segment_perp": ["W18", "W1"]} # inward (≈ north)
_W2W5_AL = {"segment": ["W2", "W5"]}        # along W2→W5 (≈ north)
_W2W5_IN = {"segment_perp": ["W2", "W5"]}   # inward (≈ east)
_W9W10_AL = {"segment": ["W9", "W10"]}      # along W9→W10 (≈ east)
_W9W10_IN = {"segment_perp": ["W9", "W10"]} # inward (≈ south)
_W6W7_AL = {"segment": ["W6", "W7"]}        # along W6→W7 (≈ east)
_W6W7_IN = {"segment_perp": ["W6", "W7"]}   # inward (≈ south)

# Helper: line intersection point spec
def _li(p1, d1, p2, d2):
    return {"type": "line_intersection",
            "line1_point": p1, "line1_dir": d1,
            "line2_point": p2, "line2_dir": d2}

# Helper: offset point spec
def _off(base, dist, direction):
    return {"offset": base, "dist": dist, "dir": direction}

# Helper: line-polygon intersection
def _lpi(origin, direction, select="farthest"):
    return {"line_poly_isect": {"origin": origin, "dir": direction,
                                "poly": "inner_poly", "select": select}}

# _w2w5_ref = line_isect(W2, w2w5_al, W1, w18w1_al)
_W2W5_REF = _li("W2", _W2W5_AL, "W1", _W18W1_AL)

# _iw2_w_anchor = offset(W2, IW2_DIST_W2W5, w2w5_in)
_IW2_W_ANCHOR = _off("W2", {"const": "IW2_DIST_W2W5"}, _W2W5_IN)

# Virtual W2 ref: offset(W7, IW_W2_REF_DIST, -w6w7_al)
_IW_W2_REF = _off("W7", {"const": "IW_W2_REF_DIST"}, {"neg": _W6W7_AL})

# _iw_al = perp(w6w7_al) rotated: (-w6w7_al[1], w6w7_al[0])
# seg_vecs gives along=(dx/L, dy/L), inward=(dy/L, -dx/L)
# _iw_al = (-inward[1], inward[0]) = (-(-dx/L), dy/L) = (dx/L, dy/L)??
# Actually: _iw_al = (-_w6w7_al[1], _w6w7_al[0]) = left perp of w6w7_al
# That's the same as neg(w6w7_in) since w6w7_in = (dy/L, -dx/L)
# _iw_al = (-dy/L, dx/L) = neg of inward? No...
# w6w7_al = (dx/L, dy/L), so (-al[1], al[0]) = (-dy/L, dx/L)
# w6w7_in = (dy/L, -dx/L)
# So _iw_al = (-dy/L, dx/L) which is neg(w6w7_in).
# Actually neg(w6w7_in) = (-dy/L, dx/L). Yes!
_IW_AL = {"neg": _W6W7_IN}   # ≈ north (left perp of W6→W7)

# _iw_in = w6w7_al
_IW_IN = _W6W7_AL  # ≈ east

# _iw2s_w_anchor = offset(_iw_w2_ref, IW2S_W2REF_OFFSET, _iw_in)
_IW2S_W_ANCHOR = _off(_IW_W2_REF, {"const": "IW2S_W2REF_OFFSET"}, _IW_IN)


def get_iw_formulas():
    """Return dict of {element_name: formula_json} for all 13 IW walls."""
    return {
        # --- IW1: west end at IW2 west face, extends east to inner poly ---
        # NW = line_isect(offset(W9, IW1_OFFSET, w9w10_in), w9w10_al,
        #                 offset(W2, IW2_DIST, w2w5_in), w2w5_al)
        # SW = offset(NW, WALL_6IN, w9w10_in)
        # NE = line_poly_isect from NW along w9w10_al, farthest
        # SE = line_poly_isect from SW along w9w10_al, farthest
        "IW1": {
            "type": "four_corner",
            "nw": _li(
                _off("W9", {"const": "IW1_OFFSET_FROM_W9"}, _W9W10_IN),
                _W9W10_AL,
                _IW2_W_ANCHOR,
                _W2W5_AL,
            ),
            "sw": _off(
                _li(
                    _off("W9", {"const": "IW1_OFFSET_FROM_W9"}, _W9W10_IN),
                    _W9W10_AL,
                    _IW2_W_ANCHOR,
                    _W2W5_AL,
                ),
                {"const": "WALL_6IN"},
                _W9W10_IN,
            ),
            "ne": _lpi(
                _li(
                    _off("W9", {"const": "IW1_OFFSET_FROM_W9"}, _W9W10_IN),
                    _W9W10_AL,
                    _IW2_W_ANCHOR,
                    _W2W5_AL,
                ),
                _W9W10_AL,
            ),
            "se": _lpi(
                _off(
                    _li(
                        _off("W9", {"const": "IW1_OFFSET_FROM_W9"}, _W9W10_IN),
                        _W9W10_AL,
                        _IW2_W_ANCHOR,
                        _W2W5_AL,
                    ),
                    {"const": "WALL_6IN"},
                    _W9W10_IN,
                ),
                _W9W10_AL,
            ),
        },

        # --- IW3: fixed-height wall near south wall ---
        # SW = offset(_w2w5_ref, -IW3_DIST_W2W5, w18w1_al)
        # SE = offset(SW, -WALL_4IN, w18w1_al)
        # height = IW7_OFFSET_FROM_W18W1 + WALL_4IN
        # NE = offset(SE, height, w18w1_in)
        # NW = offset(SW, height, w18w1_in)
        "IW3": {
            "type": "four_corner",
            "sw": _off(_W2W5_REF, {"neg": {"const": "IW3_DIST_W2W5"}}, _W18W1_AL),
            "se": _off(
                _off(_W2W5_REF, {"neg": {"const": "IW3_DIST_W2W5"}}, _W18W1_AL),
                {"neg": {"const": "WALL_4IN"}},
                _W18W1_AL,
            ),
            "ne": _off(
                _off(
                    _off(_W2W5_REF, {"neg": {"const": "IW3_DIST_W2W5"}}, _W18W1_AL),
                    {"neg": {"const": "WALL_4IN"}},
                    _W18W1_AL,
                ),
                {"add": [{"const": "IW7_OFFSET_FROM_W18W1"}, {"const": "WALL_4IN"}]},
                _W18W1_IN,
            ),
            "nw": _off(
                _off(_W2W5_REF, {"neg": {"const": "IW3_DIST_W2W5"}}, _W18W1_AL),
                {"add": [{"const": "IW7_OFFSET_FROM_W18W1"}, {"const": "WALL_4IN"}]},
                _W18W1_IN,
            ),
        },

        # --- IW9: from IW3 east face, extends to IW1 south face ---
        # SW = offset(IW3.SE, -IW3_OFFSET_IW9, w18w1_al)
        # SE = offset(SW, -WALL_4IN, w18w1_al)
        # NE = line_isect(SE, w18w1_in, IW1.SW, IW1_south_dir)
        # NW = line_isect(SW, w18w1_in, IW1.SW, IW1_south_dir)
        "IW9": {
            "type": "four_corner",
            "sw": _off(
                {"element": "IW3", "corner": "SE"},
                {"neg": {"const": "IW3_OFFSET_IW9"}},
                _W18W1_AL,
            ),
            "se": _off(
                _off(
                    {"element": "IW3", "corner": "SE"},
                    {"neg": {"const": "IW3_OFFSET_IW9"}},
                    _W18W1_AL,
                ),
                {"neg": {"const": "WALL_4IN"}},
                _W18W1_AL,
            ),
            "ne": _li(
                _off(
                    _off(
                        {"element": "IW3", "corner": "SE"},
                        {"neg": {"const": "IW3_OFFSET_IW9"}},
                        _W18W1_AL,
                    ),
                    {"neg": {"const": "WALL_4IN"}},
                    _W18W1_AL,
                ),
                _W18W1_IN,
                {"element": "IW1", "corner": "SW"},
                {"face_along": "IW1", "face": "south"},
            ),
            "nw": _li(
                _off(
                    {"element": "IW3", "corner": "SE"},
                    {"neg": {"const": "IW3_OFFSET_IW9"}},
                    _W18W1_AL,
                ),
                _W18W1_IN,
                {"element": "IW1", "corner": "SW"},
                {"face_along": "IW1", "face": "south"},
            ),
        },

        # --- IW7: spans IW3.SE → IW9.SW, 6' north of W18-W1 ---
        # SW = offset(IW3.SE, IW7_OFFSET, w18w1_in)
        # SE = offset(IW9.SW, IW7_OFFSET, w18w1_in)
        # NW = offset(SW, WALL_4IN, w18w1_in)
        # NE = offset(SE, WALL_4IN, w18w1_in)
        "IW7": {
            "type": "four_corner",
            "sw": _off({"element": "IW3", "corner": "SE"},
                       {"const": "IW7_OFFSET_FROM_W18W1"}, _W18W1_IN),
            "se": _off({"element": "IW9", "corner": "SW"},
                       {"const": "IW7_OFFSET_FROM_W18W1"}, _W18W1_IN),
            "ne": _off(
                _off({"element": "IW9", "corner": "SW"},
                     {"const": "IW7_OFFSET_FROM_W18W1"}, _W18W1_IN),
                {"const": "WALL_4IN"}, _W18W1_IN,
            ),
            "nw": _off(
                _off({"element": "IW3", "corner": "SE"},
                     {"const": "IW7_OFFSET_FROM_W18W1"}, _W18W1_IN),
                {"const": "WALL_4IN"}, _W18W1_IN,
            ),
        },

        # --- IW11: positioned by IW9-IW11 gap ---
        # _w2w5_ref_s = line_isect(W2, w2w5_al, W1, w18w1_al)  (same as _W2W5_REF)
        # _iw9_se_pos = offset(_w2w5_ref_s,
        #     -(IW3_DIST_W2W5 + 2*WALL_4IN + IW3_OFFSET_IW9), w18w1_al)
        # iw11_sw = offset(_iw9_se_pos, -IW9_IW11_GAP, w18w1_al)
        # iw11_se = offset(iw11_sw, -WALL_4IN, w18w1_al)
        # iw11_ne = line_isect(iw11_se, w18w1_in, iw1_sw, w9w10_al)
        # iw11_nw = line_isect(iw11_sw, w18w1_in, iw1_sw, w9w10_al)
        "IW11": {
            "type": "four_corner",
            "sw": _off(
                _off(_W2W5_REF,
                     {"neg": {"add": [{"const": "IW3_DIST_W2W5"},
                                      {"mul": [2, {"const": "WALL_4IN"}]},
                                      {"const": "IW3_OFFSET_IW9"}]}},
                     _W18W1_AL),
                {"neg": {"const": "IW9_IW11_GAP"}},
                _W18W1_AL,
            ),
            "se": _off(
                _off(
                    _off(_W2W5_REF,
                         {"neg": {"add": [{"const": "IW3_DIST_W2W5"},
                                          {"mul": [2, {"const": "WALL_4IN"}]},
                                          {"const": "IW3_OFFSET_IW9"}]}},
                         _W18W1_AL),
                    {"neg": {"const": "IW9_IW11_GAP"}},
                    _W18W1_AL,
                ),
                {"neg": {"const": "WALL_4IN"}},
                _W18W1_AL,
            ),
            "ne": _li(
                _off(
                    _off(
                        _off(_W2W5_REF,
                             {"neg": {"add": [{"const": "IW3_DIST_W2W5"},
                                              {"mul": [2, {"const": "WALL_4IN"}]},
                                              {"const": "IW3_OFFSET_IW9"}]}},
                             _W18W1_AL),
                        {"neg": {"const": "IW9_IW11_GAP"}},
                        _W18W1_AL,
                    ),
                    {"neg": {"const": "WALL_4IN"}},
                    _W18W1_AL,
                ),
                _W18W1_IN,
                {"element": "IW1", "corner": "SW"},
                _W9W10_AL,
            ),
            "nw": _li(
                _off(
                    _off(_W2W5_REF,
                         {"neg": {"add": [{"const": "IW3_DIST_W2W5"},
                                          {"mul": [2, {"const": "WALL_4IN"}]},
                                          {"const": "IW3_OFFSET_IW9"}]}},
                         _W18W1_AL),
                    {"neg": {"const": "IW9_IW11_GAP"}},
                    _W18W1_AL,
                ),
                _W18W1_IN,
                {"element": "IW1", "corner": "SW"},
                _W9W10_AL,
            ),
        },

        # --- IW4: parallel to IW11, west face IW4_GAP_IW11 east of IW11.SE ---
        # iw4_sw = offset(IW11.SE, -IW4_GAP_IW11, w18w1_al)
        # iw4_se = offset(iw4_sw, -WALL_4IN, w18w1_al)
        # iw4_nw = line_isect(iw4_sw, w18w1_in, IW12.NW, IW12_n_dir)
        # iw4_ne = line_isect(iw4_se, w18w1_in, IW12.NW, IW12_n_dir)
        "IW4": {
            "type": "four_corner",
            "sw": _off({"element": "IW11", "corner": "SE"},
                       {"neg": {"const": "IW4_GAP_IW11"}}, _W18W1_AL),
            "se": _off(
                _off({"element": "IW11", "corner": "SE"},
                     {"neg": {"const": "IW4_GAP_IW11"}}, _W18W1_AL),
                {"neg": {"const": "WALL_4IN"}}, _W18W1_AL,
            ),
            "ne": _li(
                _off(
                    _off({"element": "IW11", "corner": "SE"},
                         {"neg": {"const": "IW4_GAP_IW11"}}, _W18W1_AL),
                    {"neg": {"const": "WALL_4IN"}}, _W18W1_AL,
                ),
                _W18W1_IN,
                {"element": "IW12", "corner": "NW"},
                {"face_along": "IW12", "face": "north"},
            ),
            "nw": _li(
                _off({"element": "IW11", "corner": "SE"},
                     {"neg": {"const": "IW4_GAP_IW11"}}, _W18W1_AL),
                _W18W1_IN,
                {"element": "IW12", "corner": "NW"},
                {"face_along": "IW12", "face": "north"},
            ),
        },

        # --- IW12: S face IW12_S_OFFSET north of W18-W1 ---
        # spans IW11.SE → IW4.SW (IW4.SW inlined to break cycle)
        # IW4.SW = offset(IW11.SE, -IW4_GAP_IW11, w18w1_al)
        "IW12": {
            "type": "four_corner",
            "sw": _li(
                _off("W18", {"const": "IW12_S_OFFSET_W18W1"}, _W18W1_IN),
                {"neg": _W18W1_AL},
                {"element": "IW11", "corner": "SE"},
                _W18W1_IN,
            ),
            "se": _li(
                _off("W18", {"const": "IW12_S_OFFSET_W18W1"}, _W18W1_IN),
                {"neg": _W18W1_AL},
                _off({"element": "IW11", "corner": "SE"},
                     {"neg": {"const": "IW4_GAP_IW11"}}, _W18W1_AL),
                _W18W1_IN,
            ),
            "ne": _li(
                _off(
                    _off("W18", {"const": "IW12_S_OFFSET_W18W1"}, _W18W1_IN),
                    {"const": "WALL_4IN"}, _W18W1_IN,
                ),
                {"neg": _W18W1_AL},
                _off({"element": "IW11", "corner": "SE"},
                     {"neg": {"const": "IW4_GAP_IW11"}}, _W18W1_AL),
                _W18W1_IN,
            ),
            "nw": _li(
                _off(
                    _off("W18", {"const": "IW12_S_OFFSET_W18W1"}, _W18W1_IN),
                    {"const": "WALL_4IN"}, _W18W1_IN,
                ),
                {"neg": _W18W1_AL},
                {"element": "IW11", "corner": "SE"},
                _W18W1_IN,
            ),
        },

        # --- IW2: lower section, 6'6" from W2-W5 ---
        # _iw2_e_anchor = offset(_iw2_w_anchor, WALL_6IN, w2w5_in)
        # iw2_sw = line_isect(_iw2_w_anchor, w2w5_al, IW1.NW, w9w10_al)
        # iw2_se = line_isect(_iw2_e_anchor, w2w5_al, IW1.NW, w9w10_al)
        # iw2_nw = offset(iw2_sw, IW2_LENGTH, w2w5_al)
        # iw2_ne = offset(iw2_se, IW2_LENGTH, w2w5_al)
        "IW2": {
            "type": "four_corner",
            "sw": _li(_IW2_W_ANCHOR, _W2W5_AL,
                      {"element": "IW1", "corner": "NW"}, _W9W10_AL),
            "se": _li(
                _off("W2", {"add": [{"const": "IW2_DIST_W2W5"}, {"const": "WALL_6IN"}]},
                     _W2W5_IN),
                _W2W5_AL,
                {"element": "IW1", "corner": "NW"}, _W9W10_AL,
            ),
            "ne": _off(
                _li(
                    _off("W2", {"add": [{"const": "IW2_DIST_W2W5"}, {"const": "WALL_6IN"}]},
                         _W2W5_IN),
                    _W2W5_AL,
                    {"element": "IW1", "corner": "NW"}, _W9W10_AL,
                ),
                {"const": "IW2_LENGTH"}, _W2W5_AL,
            ),
            "nw": _off(
                _li(_IW2_W_ANCHOR, _W2W5_AL,
                    {"element": "IW1", "corner": "NW"}, _W9W10_AL),
                {"const": "IW2_LENGTH"}, _W2W5_AL,
            ),
        },

        # --- IW2S: upper/shower section ---
        # iw2s_nw = line_isect(_iw2s_w_anchor, _iw_al, W6, w6w7_al)
        # iw2s_ne = line_isect(_iw2s_e_anchor, _iw_al, W6, w6w7_al)
        # iw2s_sw = offset(iw2s_nw, -IW2S_LENGTH, _iw_al)
        # iw2s_se = offset(iw2s_ne, -IW2S_LENGTH, _iw_al)
        "IW2S": {
            "type": "four_corner",
            "nw": _li(_IW2S_W_ANCHOR, _IW_AL, "W6", _W6W7_AL),
            "ne": _li(
                _off(_IW_W2_REF,
                     {"add": [{"const": "IW2S_W2REF_OFFSET"}, {"const": "WALL_6IN"}]},
                     _IW_IN),
                _IW_AL, "W6", _W6W7_AL,
            ),
            "sw": _off(
                _li(_IW2S_W_ANCHOR, _IW_AL, "W6", _W6W7_AL),
                {"neg": {"const": "IW2S_LENGTH"}}, _IW_AL,
            ),
            "se": _off(
                _li(
                    _off(_IW_W2_REF,
                         {"add": [{"const": "IW2S_W2REF_OFFSET"}, {"const": "WALL_6IN"}]},
                         _IW_IN),
                    _IW_AL, "W6", _W6W7_AL,
                ),
                {"neg": {"const": "IW2S_LENGTH"}}, _IW_AL,
            ),
        },

        # --- IW2O: oblique connector from IW2 north to IW2S south ---
        # Midpoints of IW2 north face and IW2S south face
        # _iw2_n_mid = midpoint(IW2.NW, IW2.NE)
        # _iw2s_s_mid = midpoint(IW2S.SW, IW2S.SE)
        # Direction: _iw2s_s_mid - _iw2_n_mid, normalized
        # Perpendicular: left perp of that direction
        # half_thickness = IW2O_THICKNESS / 2
        # SW = mid_start - half_t * perp
        # SE = mid_start + half_t * perp
        # NW = mid_end - half_t * perp
        # NE = mid_end + half_t * perp
        # perp gives right perp [dy,-dx]; layout uses left perp [-dy,dx]
        # So: layout's -half_t*left_perp = +half_t*right_perp (SW)
        #     layout's +half_t*left_perp = -half_t*right_perp (SE)
        "IW2O": {
            "type": "four_corner",
            "sw": _off(
                {"midpoint": [{"element": "IW2", "corner": "NW"},
                              {"element": "IW2", "corner": "NE"}]},
                {"mul": [0.5, {"const": "IW2O_THICKNESS"}]},
                {"perp": {"segment": [
                    {"midpoint": [{"element": "IW2", "corner": "NW"},
                                  {"element": "IW2", "corner": "NE"}]},
                    {"midpoint": [{"element": "IW2S", "corner": "SW"},
                                  {"element": "IW2S", "corner": "SE"}]},
                ]}},
            ),
            "se": _off(
                {"midpoint": [{"element": "IW2", "corner": "NW"},
                              {"element": "IW2", "corner": "NE"}]},
                {"neg": {"mul": [0.5, {"const": "IW2O_THICKNESS"}]}},
                {"perp": {"segment": [
                    {"midpoint": [{"element": "IW2", "corner": "NW"},
                                  {"element": "IW2", "corner": "NE"}]},
                    {"midpoint": [{"element": "IW2S", "corner": "SW"},
                                  {"element": "IW2S", "corner": "SE"}]},
                ]}},
            ),
            "nw": _off(
                {"midpoint": [{"element": "IW2S", "corner": "SW"},
                              {"element": "IW2S", "corner": "SE"}]},
                {"mul": [0.5, {"const": "IW2O_THICKNESS"}]},
                {"perp": {"segment": [
                    {"midpoint": [{"element": "IW2", "corner": "NW"},
                                  {"element": "IW2", "corner": "NE"}]},
                    {"midpoint": [{"element": "IW2S", "corner": "SW"},
                                  {"element": "IW2S", "corner": "SE"}]},
                ]}},
            ),
            "ne": _off(
                {"midpoint": [{"element": "IW2S", "corner": "SW"},
                              {"element": "IW2S", "corner": "SE"}]},
                {"neg": {"mul": [0.5, {"const": "IW2O_THICKNESS"}]}},
                {"perp": {"segment": [
                    {"midpoint": [{"element": "IW2", "corner": "NW"},
                                  {"element": "IW2", "corner": "NE"}]},
                    {"midpoint": [{"element": "IW2S", "corner": "SW"},
                                  {"element": "IW2S", "corner": "SE"}]},
                ]}},
            ),
        },

        # --- IW8: perpendicular to W2-W5, centered between W18-W1 and W6-W7 ---
        # _d = proj(W6, W1, w18w1_in)  → distance W1 to W6 along w18w1_in
        # _mid = d / 2
        # _iw8_s_anchor = offset(W18, _mid - WALL_6IN/2, w18w1_in)
        # _iw8_n_anchor = offset(_iw8_s_anchor, WALL_6IN, w2w5_al)
        # iw8_sw = line_isect(_s_anchor, w2w5_in, W2, w2w5_al)
        # iw8_nw = line_isect(_n_anchor, w2w5_in, W2, w2w5_al)
        # iw8_se = line_isect(_s_anchor, w2w5_in, _iw2_w_anchor, w2w5_al)
        # iw8_ne = line_isect(_n_anchor, w2w5_in, _iw2_w_anchor, w2w5_al)
        "IW8": {
            "type": "four_corner",
            "sw": _li(
                _off("W18",
                     {"sub": [{"mul": [0.5, {"proj": {"target": "W6", "anchor": "W1", "dir": _W18W1_IN}}]},
                              {"mul": [0.5, {"const": "WALL_6IN"}]}]},
                     _W18W1_IN),
                _W2W5_IN, "W2", _W2W5_AL,
            ),
            "nw": _li(
                _off("W18",
                     {"add": [{"sub": [{"mul": [0.5, {"proj": {"target": "W6", "anchor": "W1", "dir": _W18W1_IN}}]},
                                       {"mul": [0.5, {"const": "WALL_6IN"}]}]},
                              {"const": "WALL_6IN"}]},
                     _W18W1_IN),
                _W2W5_IN, "W2", _W2W5_AL,
            ),
            "se": _li(
                _off("W18",
                     {"sub": [{"mul": [0.5, {"proj": {"target": "W6", "anchor": "W1", "dir": _W18W1_IN}}]},
                              {"mul": [0.5, {"const": "WALL_6IN"}]}]},
                     _W18W1_IN),
                _W2W5_IN, _IW2_W_ANCHOR, _W2W5_AL,
            ),
            "ne": _li(
                _off("W18",
                     {"add": [{"sub": [{"mul": [0.5, {"proj": {"target": "W6", "anchor": "W1", "dir": _W18W1_IN}}]},
                                       {"mul": [0.5, {"const": "WALL_6IN"}]}]},
                              {"const": "WALL_6IN"}]},
                     _W18W1_IN),
                _W2W5_IN, _IW2_W_ANCHOR, _W2W5_AL,
            ),
        },

        # --- IW5: S face 30" south of IW1 S face ---
        # _iw5_s_anchor = offset(IW1.SW, IW5_S_OFFSET, w9w10_in)
        # _iw5_n_anchor = offset(_s_anchor, -WALL_3IN, w9w10_in)
        # iw5_ne = offset(_n_anchor, proj(W15, _n_anchor, w9w10_al), w9w10_al)
        # iw5_se = offset(_s_anchor, proj(W15, _s_anchor, w9w10_al), w9w10_al)
        # iw5_nw = line_isect(_n_anchor, w9w10_al, IW11.SE, w18w1_in)
        # iw5_sw = line_isect(_s_anchor, w9w10_al, IW11.SE, w18w1_in)
        "IW5": {
            "type": "four_corner",
            "sw": _li(
                _off({"element": "IW1", "corner": "SW"},
                     {"const": "IW5_S_OFFSET_FROM_IW1"}, _W9W10_IN),
                _W9W10_AL,
                {"element": "IW11", "corner": "SE"},
                _W18W1_IN,
            ),
            "se": _lpi(
                _off({"element": "IW1", "corner": "SW"},
                     {"const": "IW5_S_OFFSET_FROM_IW1"}, _W9W10_IN),
                _W9W10_AL,
            ),
            "ne": _lpi(
                _off(
                    _off({"element": "IW1", "corner": "SW"},
                         {"const": "IW5_S_OFFSET_FROM_IW1"}, _W9W10_IN),
                    {"neg": {"const": "WALL_3IN"}}, _W9W10_IN,
                ),
                _W9W10_AL,
            ),
            "nw": _li(
                _off(
                    _off({"element": "IW1", "corner": "SW"},
                         {"const": "IW5_S_OFFSET_FROM_IW1"}, _W9W10_IN),
                    {"neg": {"const": "WALL_3IN"}}, _W9W10_IN,
                ),
                _W9W10_AL,
                {"element": "IW11", "corner": "SE"},
                _W18W1_IN,
            ),
        },

        # --- IW6: 1" partition, 5'6" from W6 ---
        # _iw6_n_anchor = offset(W6, IW6_OFFSET, w6w7_in)
        # _iw6_s_anchor = offset(_n_anchor, IW6_THICKNESS, w6w7_in)
        # iw6_ne = line_isect(_n_anchor, w6w7_al, _iw2s_w_anchor, _iw_al)
        # iw6_se = line_isect(_s_anchor, w6w7_al, _iw2s_w_anchor, _iw_al)
        # iw6_nw = line_poly_isect(_n_anchor, -w6w7_al, inner_poly)
        # iw6_sw = line_poly_isect(_s_anchor, -w6w7_al, inner_poly)
        "IW6": {
            "type": "four_corner",
            "ne": _li(
                _off("W6", {"const": "IW6_OFFSET_FROM_W6"}, _W6W7_IN),
                _W6W7_AL,
                _IW2S_W_ANCHOR,
                _IW_AL,
            ),
            "se": _li(
                _off(
                    _off("W6", {"const": "IW6_OFFSET_FROM_W6"}, _W6W7_IN),
                    {"const": "IW6_THICKNESS"}, _W6W7_IN,
                ),
                _W6W7_AL,
                _IW2S_W_ANCHOR,
                _IW_AL,
            ),
            "nw": _lpi(
                _off("W6", {"const": "IW6_OFFSET_FROM_W6"}, _W6W7_IN),
                {"neg": _W6W7_AL},
            ),
            "sw": _lpi(
                _off(
                    _off("W6", {"const": "IW6_OFFSET_FROM_W6"}, _W6W7_IN),
                    {"const": "IW6_THICKNESS"}, _W6W7_IN,
                ),
                {"neg": _W6W7_AL},
            ),
        },
    }


def get_layout_item_formulas():
    """Return dict of {element_name: formula_json} for core layout items.

    These are the 5 items computed by floorplan/layout.py: dryer, washer,
    counter, dresser, shelves.  The bed is also here but depends on O9
    opening parametrics (IW11-anchored).
    """
    # --- Dryer ---
    # SW = line_isect(offset(W2, APPLIANCE_OFFSET_FROM_W2, w2w5_in), w2w5_al,
    #                 offset(W1, APPLIANCE_OFFSET_FROM_W1, w2w5_al), w9w10_al)
    # SE = offset(SW, APPLIANCE_WIDTH, w2w5_in)
    # NW = offset(SW, APPLIANCE_DEPTH, w2w5_al)
    _dryer_sw = _li(
        _off("W2", {"const": "APPLIANCE_OFFSET_FROM_W2"}, _W2W5_IN),
        _W2W5_AL,
        _off("W1", {"const": "APPLIANCE_OFFSET_FROM_W1"}, _W2W5_AL),
        _W9W10_AL,
    )

    # --- Washer ---
    # SW = offset(dryer.NW, APPLIANCE_GAP, w2w5_al)
    _washer_sw = _off(
        {"element": "DRYER", "corner": "NW"},
        {"const": "APPLIANCE_GAP"},
        _W2W5_AL,
    )

    # --- Counter ---
    # SW anchor = offset(dryer.SE, COUNTER_GAP, w2w5_in)
    # counter has 4 corners, each a line intersection
    # ctr_sw = line_isect(ctr_sw_anchor, w2w5_al, W1, w9w10_al)
    # ctr_se = line_isect(ctr_se_anchor, w2w5_al, W1, w9w10_al)
    # ctr_nw = offset(ctr_sw, COUNTER_LENGTH, w2w5_al)
    # ctr_ne = offset(ctr_se, COUNTER_LENGTH, w2w5_al)
    # where ctr_se_anchor = offset(ctr_sw_anchor, COUNTER_DEPTH, w2w5_in)
    _ctr_sw_anchor = _off(
        {"element": "DRYER", "corner": "SE"},
        {"const": "COUNTER_GAP"},
        _W2W5_IN,
    )
    _ctr_se_anchor = _off(
        _ctr_sw_anchor,
        {"const": "COUNTER_DEPTH"},
        _W2W5_IN,
    )

    # --- Dresser ---
    # NE = line_isect(offset(IW11.NW, DRESSER_GAP_IW15, w18w1_al), w18w1_in,
    #                 offset(IW1.SW, DRESSER_GAP_IW1, w9w10_in), w9w10_al)
    # along = w9w10_al (east), across = w9w10_in (south)
    # width goes westward from NE: NW = offset(NE, -DRESSER_WIDTH, w9w10_al)
    # depth goes southward from NE: SE = offset(NE, DRESSER_DEPTH, w9w10_in)
    _dresser_ne = _li(
        _off({"element": "IW11", "corner": "NW"}, {"const": "DRESSER_GAP_IW15"},
             _W18W1_AL),
        _W18W1_IN,
        _off({"element": "IW1", "corner": "SW"}, {"const": "DRESSER_GAP_IW1"},
             _W9W10_IN),
        _W9W10_AL,
    )

    # --- Shelves ---
    # NE = line_isect(offset(IW1.SW, SHELVES_GAP_IW1, w9w10_in), w9w10_al,
    #                 offset(IW9.NW, SHELVES_GAP_IW9, w18w1_al), w18w1_in)
    # NW = offset(NE, SHELVES_LENGTH, w18w1_al)  [westward]
    # SE = offset(NE, SHELVES_DEPTH, w9w10_in)   [southward]
    _shelves_ne = _li(
        _off({"element": "IW1", "corner": "SW"}, {"const": "SHELVES_GAP_IW1"},
             _W9W10_IN),
        _W9W10_AL,
        _off({"element": "IW9", "corner": "NW"}, {"const": "SHELVES_GAP_IW9"},
             _W18W1_AL),
        _W18W1_IN,
    )

    # --- Bed ---
    # The bed is anchored relative to opening O9, which is anchored to IW11.
    # O9 east end on F18-F1: te9 = t_iw11_sw + O9_OFFSET/seg_len + O9_WIDTH/seg_len
    # Bed SE on wall = offset(W18, bed_t, (dE, dN)) where bed_t = te9 + BED_GAP_O9/seg_len
    # But the formula system works with W-series points (inner wall face), while
    # O9 uses F-series (outer face).  The bed SE is:
    #   bed_se_wall = W18 + (te9 + BED_GAP_O9/seg_len) * (W1 - W18)
    #   bed_se = offset(bed_se_wall, BED_WALL_GAP, w18w1_in)
    # te9 is computed from IW11.SW projected onto F18-F1.  Rather than replicate
    # this parametric chain, we use a four_corner formula with proj-based offsets
    # from IW11.SW along the south wall.
    #
    # Procedural code computes:
    #   _t_sw9 = proj(IW11.SW, F18, seg_F18_F1_unit)
    #   _ts9 = _t_sw9 + O9_OFFSET_IW11 / seg9_len
    #   _te9 = _ts9 + 2 * O9_HALF_WIDTH / seg9_len
    #   _bed_t = _te9 + BED_GAP_O9 / _seg_len   (on W18-W1, not F18-F1!)
    #   _bed_se_wall = W18 + _bed_t * (W1 - W18)  [along W18-W1 direction]
    #
    # The key subtlety: _bed_t uses _te9 (F18-F1 parametric) applied to
    # W18-W1 segment length.  Since F18-F1 and W18-W1 are parallel,
    # _te9 = dist_along / F18F1_len, and _bed_t = _te9 + BED_GAP_O9/W18W1_len.
    # These have different denominators (F vs W segment lengths), which is
    # a quirk of the procedural code.
    #
    # For formula equivalence, we compute the bed anchor directly:
    #   bed_along_dist = proj(IW11.SW, W18, w18w1_neg_al) + O9_OFFSET_IW11
    #                    + 2*O9_HALF_WIDTH + BED_GAP_O9
    # where all distances use W18-W1 direction but proj uses F18-F1 for IW11.
    # This won't be bit-identical because F18≠W18 offset matters.
    #
    # Instead, replicate the exact procedural formula using a four_corner type
    # with explicit offset computations matching the F18-F1 parametric approach.
    # We'll compute the bed position as:
    #   bed_se_wall = offset(W18, total_dist_along_w18w1, neg(w18w1_al))
    #   where total_dist_along_w18w1 = proj(IW11.SW, F18, F18F1_dir)/F18F1_len * W18W1_len
    #                                 + O9_OFFSET_IW11/F18F1_len * W18W1_len
    #                                 + 2*O9_HALF_WIDTH/F18F1_len * W18W1_len
    #                                 + BED_GAP_O9
    # This is getting complex. Let's use four_corner and match exactly.
    #
    # Actually: the procedural code uses seg_vec(F18, F1) for the parametric
    # but seg_vec(W18, W1) for _bed_t denominator on the GAP term. So:
    #   _bed_se_wall = W18 + (_te9 + BED_GAP_O9/_seg_w18w1_len) * seg_w18w1_vec
    # where _te9 = (proj_iw11_sw_on_f18f1 + O9_OFFSET_IW11 + O9_WIDTH) / seg_f18f1_len
    #
    # For bit-identical results, we need the exact same computation.
    # The simplest approach: compute bed_se_wall as
    #   W18 + proj(IW11.SW, F18, f18f1_dir) * w18w1_unit
    #        + (O9_OFFSET_IW11 + 2*O9_HALF_WIDTH) / f18f1_len * w18w1_vec
    #        + BED_GAP_O9 * w18w1_unit
    # But we don't have f18f1_len in the formula system (F-series points aren't
    # in the evaluator context).
    #
    # For now, skip the bed — it requires F-series points which the evaluator
    # doesn't have access to yet.  The bed will be migrated when F-series
    # support is added.

    return {
        "DRYER": {
            "type": "item_rect",
            "anchor": _dryer_sw,
            "along": _W2W5_IN,
            "across": _W2W5_AL,
            "width": {"const": "APPLIANCE_WIDTH"},
            "depth": {"const": "APPLIANCE_DEPTH"},
            "anchor_corner": "sw",
        },
        "WASHER": {
            "type": "item_rect",
            "anchor": _washer_sw,
            "along": _W2W5_IN,
            "across": _W2W5_AL,
            "width": {"const": "APPLIANCE_WIDTH"},
            "depth": {"const": "APPLIANCE_DEPTH"},
            "anchor_corner": "sw",
        },
        "COUNTER": {
            "type": "four_corner",
            "sw": _li(_ctr_sw_anchor, _W2W5_AL, "W1", _W9W10_AL),
            "se": _li(_ctr_se_anchor, _W2W5_AL, "W1", _W9W10_AL),
            "ne": _off(
                _li(_ctr_se_anchor, _W2W5_AL, "W1", _W9W10_AL),
                {"const": "COUNTER_LENGTH"},
                _W2W5_AL,
            ),
            "nw": _off(
                _li(_ctr_sw_anchor, _W2W5_AL, "W1", _W9W10_AL),
                {"const": "COUNTER_LENGTH"},
                _W2W5_AL,
            ),
        },
        "DRESSER": {
            "type": "item_rect",
            "anchor": _dresser_ne,
            "along": _W9W10_AL,
            "across": {"neg": _W9W10_IN},
            "width": {"const": "DRESSER_WIDTH"},
            "depth": {"const": "DRESSER_DEPTH"},
            "anchor_corner": "ne",
        },
        "SHELVES": {
            "type": "item_rect",
            "anchor": _shelves_ne,
            "along": {"neg": _W18W1_AL},
            "across": {"neg": _W9W10_IN},
            "width": {"const": "SHELVES_LENGTH"},
            "depth": {"const": "SHELVES_DEPTH"},
            "anchor_corner": "ne",
        },
    }


def _add(a, b):
    """Helper: length addition spec."""
    return {"add": [a, b]}


def get_outer_opening_formulas():
    """Return dict of {name: formula_json} for outer openings O1-O11, O8a.

    Each formula is a wall_opening type with parametric positioning on
    paired outer/inner wall segments.
    """
    # Cumulative gap from F5 for F2-F5 openings:
    # O3: gap = O3_GAP_F5
    # O2: gap = O3_GAP_F5 + O3_WIDTH + O2_GAP_O3
    # O1: gap = O3_GAP_F5 + O3_WIDTH + O2_GAP_O3 + O2_WIDTH + O1_GAP_O2
    _o3_gap = {"const": "O3_GAP_F5"}
    _o2_gap = _add(_add(_o3_gap, {"const": "O3_WIDTH"}), {"const": "O2_GAP_O3"})
    _o1_gap = _add(_add(_o2_gap, {"const": "O2_WIDTH"}), {"const": "O1_GAP_O2"})

    # F18-F1 south wall openings (O9 anchored to IW11.SW projection):
    # O9: ref_point = IW11.SW on F18-F1, ref_offset = O9_OFFSET_IW11
    # O8a: positioned relative to O9 start
    #   O9_start_t = t_ref + O9_OFFSET_IW11/seg_len
    #   O8a_end_t = O9_start_t - O8A_GAP_O9/seg_len
    #   O8a_start_t = O8a_end_t - 2*O8A_HALF_WIDTH/seg_len
    # O10: O9_end_t + O9_O10_WALL/seg_len
    # O11: O10_end_t + O10_O11_WALL/seg_len
    # For O8a/O10/O11, express as ref_point + cumulative ref_offset.
    _o9_start_offset = {"const": "O9_OFFSET_IW11"}
    _o9_width = {"mul": [2, {"const": "O9_HALF_WIDTH"}]}
    _o10_start_offset = _add(_add(_o9_start_offset, _o9_width),
                             {"const": "O9_O10_WALL"})
    _o10_width = {"mul": [2, {"const": "O10_HALF_WIDTH"}]}
    _o11_start_offset = _add(_add(_o10_start_offset, _o10_width),
                             {"const": "O10_O11_WALL"})
    # O8a is before O9: end_offset = O9_OFFSET_IW11 - O8A_GAP_O9
    # start_offset = end_offset - 2*O8A_HALF_WIDTH
    _o8a_start_offset = {"sub": [_o9_start_offset,
                                 _add({"const": "O8A_GAP_O9"},
                                      {"mul": [2, {"const": "O8A_HALF_WIDTH"}]})]}

    return {
        # --- F2-F5 segment (south-east wall) ---
        "O3": {
            "type": "wall_opening",
            "outer_start": "F2", "outer_end": "F5",
            "inner_start": "W2", "inner_end": "W5",
            "from_end": True,
            "gap": {"const": "O3_GAP_F5"},
            "width": {"const": "O3_WIDTH"},
        },
        "O2": {
            "type": "wall_opening",
            "outer_start": "F2", "outer_end": "F5",
            "inner_start": "W2", "inner_end": "W5",
            "from_end": True,
            "gap": _o2_gap,
            "width": {"const": "O2_WIDTH"},
        },
        "O1": {
            "type": "wall_opening",
            "outer_start": "F2", "outer_end": "F5",
            "inner_start": "W2", "inner_end": "W5",
            "from_end": True,
            "gap": _o1_gap,
            "width": {"const": "O1_WIDTH"},
        },

        # --- F6-F7 segment (east wall) ---
        "O4": {
            "type": "wall_opening",
            "outer_start": "F6", "outer_end": "F7",
            "inner_start": "W6", "inner_end": "W7",
            "centered": True,
            "center_t": 0.5,
            "width": {"mul": [2, {"const": "O4_HALF_WIDTH"}]},
            "poly_order": "inner_first",
        },

        # --- F9-F10 segment (north-east wall) ---
        "O5": {
            "type": "wall_opening",
            "outer_start": "F9", "outer_end": "F10",
            "inner_start": "W9", "inner_end": "W10",
            "ref_point": {"face_mid": "IW2S", "face": "east"},
            "ref_offset": {"sub": [{"const": "O5_OFFSET_FROM_IW2"},
                                   {"const": "O5_WIDTH"}]},
            "width": {"const": "O5_WIDTH"},
            "poly_order": "inner_first",
        },
        "O6": {
            "type": "wall_opening",
            "outer_start": "F9", "outer_end": "F10",
            "inner_start": "W9", "inner_end": "W10",
            "from_end": True,
            "gap": {"const": "O6_GAP_F10"},
            "width": {"const": "O6_WIDTH"},
            "poly_order": "inner_first",
        },

        # --- F12-F13 segment (north-west wall, diagonal) ---
        "O7": {
            "type": "wall_opening",
            "outer_start": "F12", "outer_end": "F13",
            "inner_start": "W12", "inner_end": "W13",
            "from_end": False,
            "gap": {"const": "O7_NW_GAP"},
            "width": {"mul": [2, {"const": "O7_HALF_WIDTH"}]},
        },

        # --- F14-F15 segment (south-west wall) ---
        "O8": {
            "type": "wall_opening",
            "outer_start": "F14", "outer_end": "F15",
            "inner_start": "W14", "inner_end": "W15",
            "center_refs": [
                {"face_mid": "IW5", "face": "south"},
                "F15",
            ],
            "width": {"mul": [2, {"const": "O8_HALF_WIDTH"}]},
            "poly_order": "outer_reversed",
        },

        # --- F18-F1 segment (south wall) ---
        "O9": {
            "type": "wall_opening",
            "outer_start": "F18", "outer_end": "F1",
            "inner_start": "W18", "inner_end": "W1",
            "ref_point": {"element": "IW11", "corner": "SW"},
            "ref_offset": _o9_start_offset,
            "width": _o9_width,
        },
        "O8a": {
            "type": "wall_opening",
            "outer_start": "F18", "outer_end": "F1",
            "inner_start": "W18", "inner_end": "W1",
            "ref_point": {"element": "IW11", "corner": "SW"},
            "ref_offset": _o8a_start_offset,
            "width": {"mul": [2, {"const": "O8A_HALF_WIDTH"}]},
        },
        "O10": {
            "type": "wall_opening",
            "outer_start": "F18", "outer_end": "F1",
            "inner_start": "W18", "inner_end": "W1",
            "ref_point": {"element": "IW11", "corner": "SW"},
            "ref_offset": _o10_start_offset,
            "width": _o10_width,
        },
        "O11": {
            "type": "wall_opening",
            "outer_start": "F18", "outer_end": "F1",
            "inner_start": "W18", "inner_end": "W1",
            "ref_point": {"element": "IW11", "corner": "SW"},
            "ref_offset": _o11_start_offset,
            "width": {"mul": [2, {"const": "O11_HALF_WIDTH"}]},
        },
    }


def get_rough_opening_formulas():
    """Return dict of {name: formula_json} for rough openings RO1-RO7.

    RO1-RO4, RO6-RO7 use wall_opening on rectangular IW wall faces.
    RO5 uses four_corner (IW6 is trapezoidal).
    """
    return {
        # RO1: in IW1, ref = IW2 east face midpoint
        # poly = [SW+start*unit, SW+end*unit, NW+end*unit, NW+start*unit]
        # = outer_first (south face outer, north face inner)
        "RO1": {
            "type": "wall_opening",
            "outer_start": {"element": "IW1", "corner": "SW"},
            "outer_end": {"element": "IW1", "corner": "SE"},
            "inner_start": {"element": "IW1", "corner": "NW"},
            "inner_end": {"element": "IW1", "corner": "NE"},
            "ref_point": {"face_mid": "IW2", "face": "east"},
            "ref_offset": {"const": "RO1_OFFSET_FROM_IW2"},
            "width": {"const": "IW1_RO_WIDTH"},
        },

        # RO2: in IW11, centered between IW12.NW and IW5.SW projections
        # _ro_poly_bbox(SE, SW, unit(SE→NE), ...) → face_pair order
        "RO2": {
            "type": "wall_opening",
            "outer_start": {"element": "IW11", "corner": "SE"},
            "outer_end": {"element": "IW11", "corner": "NE"},
            "inner_start": {"element": "IW11", "corner": "SW"},
            "inner_end": {"element": "IW11", "corner": "NW"},
            "center_refs": [
                {"element": "IW12", "corner": "NW"},
                {"element": "IW5", "corner": "SW"},
            ],
            "width": {"const": "IW4_RO_WIDTH"},
            "poly_order": "face_pair",
        },

        # RO3: in IW9, south edge = IW7.NE projection + RO3_IW7_GAP
        # _ro_poly_bbox(IW9.SE, IW9.SW, unit(SE→NE), start, start+width)
        "RO3": {
            "type": "wall_opening",
            "outer_start": {"element": "IW9", "corner": "SE"},
            "outer_end": {"element": "IW9", "corner": "NE"},
            "inner_start": {"element": "IW9", "corner": "SW"},
            "inner_end": {"element": "IW9", "corner": "NW"},
            "ref_point": {"element": "IW7", "corner": "NE"},
            "ref_offset": {"const": "RO3_IW7_GAP"},
            "width": {"const": "RO3_WIDTH"},
            "poly_order": "face_pair",
        },

        # RO4: in IW2O, centered along length
        # _ro_poly_bbox(IW2O.SW, IW2O.SE, unit(SW→NW), ...)
        "RO4": {
            "type": "wall_opening",
            "outer_start": {"element": "IW2O", "corner": "SW"},
            "outer_end": {"element": "IW2O", "corner": "NW"},
            "inner_start": {"element": "IW2O", "corner": "SE"},
            "inner_end": {"element": "IW2O", "corner": "NE"},
            "centered": True,
            "center_t": 0.5,
            "width": {"const": "IW2_RO_WIDTH"},
            "poly_order": "face_pair",
        },

        # RO5: in IW6 (trapezoidal), four_corner needed
        # Each face independently projects IW2S west face midpoint,
        # then offsets by -IW6_RO_OFFSET_W to get end, -width to get start
        # South face: IW6.SW → IW6.SE, ref = IW2S west mid
        # North face: IW6.NW → IW6.NE, ref = IW2S west mid
        "RO5": {
            "type": "four_corner",
            "sw": {
                "offset": {"element": "IW6", "corner": "SW"},
                "dir": {"face_along": "IW6", "face": "south"},
                "dist": {"sub": [
                    {"sub": [
                        {"proj": {"target": {"face_mid": "IW2S", "face": "west"},
                                  "anchor": {"element": "IW6", "corner": "SW"},
                                  "dir": {"face_along": "IW6", "face": "south"}}},
                        {"const": "IW6_RO_OFFSET_W"},
                    ]},
                    {"const": "IW6_RO_WIDTH"},
                ]},
            },
            "se": {
                "offset": {"element": "IW6", "corner": "SW"},
                "dir": {"face_along": "IW6", "face": "south"},
                "dist": {"sub": [
                    {"proj": {"target": {"face_mid": "IW2S", "face": "west"},
                              "anchor": {"element": "IW6", "corner": "SW"},
                              "dir": {"face_along": "IW6", "face": "south"}}},
                    {"const": "IW6_RO_OFFSET_W"},
                ]},
            },
            "ne": {
                "offset": {"element": "IW6", "corner": "NW"},
                "dir": {"neg": {"face_along": "IW6", "face": "north"}},
                "dist": {"sub": [
                    {"proj": {"target": {"face_mid": "IW2S", "face": "west"},
                              "anchor": {"element": "IW6", "corner": "NW"},
                              "dir": {"neg": {"face_along": "IW6", "face": "north"}}}},
                    {"const": "IW6_RO_OFFSET_W"},
                ]},
            },
            "nw": {
                "offset": {"element": "IW6", "corner": "NW"},
                "dir": {"neg": {"face_along": "IW6", "face": "north"}},
                "dist": {"sub": [
                    {"sub": [
                        {"proj": {"target": {"face_mid": "IW2S", "face": "west"},
                                  "anchor": {"element": "IW6", "corner": "NW"},
                                  "dir": {"neg": {"face_along": "IW6", "face": "north"}}}},
                        {"const": "IW6_RO_OFFSET_W"},
                    ]},
                    {"const": "IW6_RO_WIDTH"},
                ]},
            },
        },

        # RO6: in IW11, centered between IW12.SW and W18 projections
        "RO6": {
            "type": "wall_opening",
            "outer_start": {"element": "IW11", "corner": "SE"},
            "outer_end": {"element": "IW11", "corner": "NE"},
            "inner_start": {"element": "IW11", "corner": "SW"},
            "inner_end": {"element": "IW11", "corner": "NW"},
            "center_refs": [
                {"element": "IW12", "corner": "SW"},
                "W18",
            ],
            "width": {"const": "IW11_RO_WIDTH"},
            "poly_order": "face_pair",
        },

        # RO7: in IW9, centered between IW9.SE (origin) and IW7.SW
        "RO7": {
            "type": "wall_opening",
            "outer_start": {"element": "IW9", "corner": "SE"},
            "outer_end": {"element": "IW9", "corner": "NE"},
            "inner_start": {"element": "IW9", "corner": "SW"},
            "inner_end": {"element": "IW9", "corner": "NW"},
            "center_refs": [
                {"element": "IW9", "corner": "SE"},
                {"element": "IW7", "corner": "SW"},
            ],
            "width": {"const": "IW9_RO_WIDTH"},
            "poly_order": "face_pair",
        },
    }


# ---------------------------------------------------------------------------
# Variant item formula definitions (Phase 12e)
# ---------------------------------------------------------------------------
# Direction shorthands for variant item formulas
_W11W12_AL = {"segment": ["W11", "W12"]}
_W11W12_IN = {"segment_perp": ["W11", "W12"]}
_W12W13_AL = {"segment": ["W12", "W13"]}
_W16W17_AL = {"segment": ["W16", "W17"]}
_W16W17_IN = {"segment_perp": ["W16", "W17"]}

# IW8 face directions
_IW8_AL = {"face_along": "IW8", "face": "south"}
_IW8_IN = {"segment_perp": [{"element": "IW8", "corner": "SW"},
                              {"element": "IW8", "corner": "SE"}]}
_IW8_OUT = {"neg": _IW8_IN}
_IW8_N_AL = {"face_along": "IW8", "face": "north"}

# IW2s east face directions
_IW2S_E_AL = {"segment": [{"element": "IW2S", "corner": "SE"},
                            {"element": "IW2S", "corner": "NE"}]}
_IW2S_E_OUT = {"segment_perp": [{"element": "IW2S", "corner": "SE"},
                                  {"element": "IW2S", "corner": "NE"}]}

# IW1 north face directions
_IW1_N_AL = {"segment": [{"element": "IW1", "corner": "NW"},
                           {"element": "IW1", "corner": "NE"}]}
_IW1_N_CW = {"segment_perp": [{"element": "IW1", "corner": "NW"},
                                {"element": "IW1", "corner": "NE"}]}
_IW1_N_OUT = {"neg": _IW1_N_CW}
_IW1_W = {"neg": _IW1_N_AL}

# IW12 corner: intersection of IW2 east face extended with IW1 north face
_IW12_CORNER = _li(
    {"element": "IW2", "corner": "SE"},
    {"segment": [{"element": "IW2", "corner": "SE"},
                  {"element": "IW2", "corner": "NE"}]},
    {"element": "IW1", "corner": "NW"},
    _IW1_N_AL,
)

# IW4/IW1 corner: intersection of IW4 west face extended with IW1 north face
_IW41_CORNER = _li(
    {"element": "IW4", "corner": "NW"},
    {"segment": [{"element": "IW4", "corner": "NW"},
                  {"element": "IW4", "corner": "SW"}]},
    {"element": "IW1", "corner": "NW"},
    _IW1_N_AL,
)

# IW2s NE corner distance along north wall from W9
_IW2S_NE_D = {"proj": {"target": {"element": "IW2S", "corner": "NE"},
                         "anchor": "W9",
                         "dir": _W9W10_AL}}


def _nwp(d_along, d_inward=0):
    """North wall point spec: W9 + d_along * w9w10_al + d_inward * w9w10_in."""
    base = _off("W9", d_along, _W9W10_AL)
    if isinstance(d_inward, (int, float)) and d_inward == 0:
        return base
    return _off(base, d_inward, _W9W10_IN)


def _iwp(d_e, d_n=0):
    """IW12 corner + d_e * iw1_n_al + d_n * iw1_n_out."""
    base = _off(_IW12_CORNER, d_e, _IW1_N_AL)
    if isinstance(d_n, (int, float)) and d_n == 0:
        return base
    return _off(base, d_n, _IW1_N_OUT)


def _lwp(d_w=0, d_n=0):
    """IW41 corner + d_w * (-iw1_n_al) + d_n * iw1_n_out."""
    base = _off(_IW41_CORNER, d_w, _IW1_W) if not (isinstance(d_w, (int, float)) and d_w == 0) else _IW41_CORNER
    if isinstance(d_n, (int, float)) and d_n == 0:
        return base
    return _off(base, d_n, _IW1_N_OUT)


def _c(name):
    """Shorthand for constant reference."""
    return {"const": name}


def _in(val):
    """Convert inches to feet (literal)."""
    return val / 12.0


def get_variant_item_formulas():
    """Return list of (element_name, variant, formula_json) for variant items.

    variant is None for items common to all non-bare/sf variants,
    or a specific variant name for variant-specific items.
    """
    formulas = []

    def _f(name, formula, variant=None):
        formulas.append((name, variant, formula))

    # IW2s NE distance along north wall (used by many kitchen items)
    _iw2_d = _IW2S_NE_D

    # Kitchen chain distances along north wall
    # _st_d = _iw2_d + NORTH_CTR_LENGTH + KITCHEN_APPL_GAP
    _st_d = {"add": [_iw2_d, _c("NORTH_CTR_LENGTH"), _c("KITCHEN_APPL_GAP")]}
    # _ks_d = _st_d + STOVE_WIDTH + KITCHEN_APPL_GAP + 2/12
    _ks_d = {"add": [_st_d, _c("STOVE_WIDTH"), _c("KITCHEN_APPL_GAP"), _in(2)]}
    # _dw_d = _ks_d + KITCHEN_SINK_WIDTH + KITCHEN_APPL_GAP
    _dw_d = {"add": [_ks_d, _c("KITCHEN_SINK_WIDTH"), _c("KITCHEN_APPL_GAP")]}

    # ===================================================================
    # HAMPER — all variants, positioned relative to washer NW
    # ===================================================================
    # hamper sw = offset(offset(W2, _washer_nw_d + 2/12, w2w5_al), 2/12, w2w5_in)
    # _washer_nw_d = proj(WASHER.NW, W2, w2w5_al)
    _washer_nw_d = {"proj": {"target": {"element": "WASHER", "corner": "NW"},
                              "anchor": "W2", "dir": _W2W5_AL}}
    _hm_sw = _off(
        _off("W2", {"add": [_washer_nw_d, _in(2)]}, _W2W5_AL),
        _in(2), _W2W5_IN,
    )
    _f("hamper", {
        "type": "item_rect",
        "anchor": _hm_sw,
        "along": _W2W5_IN,
        "across": _W2W5_AL,
        "width": _c("HAMPER_W"),
        "depth": _c("HAMPER_D"),
        "anchor_corner": "sw",
    })

    # ===================================================================
    # WATER HEATER — all variants, circle tangent to arc at C7
    # ===================================================================
    # _wh_ref = offset(IW2S.NE, WH_RADIUS, _iw2_e_out)
    _wh_ref = _off({"element": "IW2S", "corner": "NE"}, _c("WH_RADIUS"), _IW2S_E_OUT)
    # wh_tangent_r = R_a7 - WALL_OUTER - WH_RADIUS
    _wh_tangent_r = {"sub": [{"sub": [{"radius_key": "R_a7"}, _c("WALL_OUTER")]},
                              _c("WH_RADIUS")]}
    _f("water_heater", {
        "type": "item_circle",
        "center": {"ray_circle_isect": {
            "origin": _wh_ref,
            "dir": _IW2S_E_AL,
            "target": "C7",
            "distance": _wh_tangent_r,
            "select": "farthest",
        }},
        "radius": _c("WH_RADIUS"),
    })

    # ===================================================================
    # TOILETS — all variants
    # ===================================================================
    # South toilet: center on IW8 south face, aligned with dryer centroid
    # _d_dryer_al = proj(DRYER_centroid, IW8.SW, iw8_al)
    _dryer_cx = {"element_centroid": "DRYER"}
    _d_dryer_al = {"proj": {"target": _dryer_cx,
                             "anchor": {"element": "IW8", "corner": "SW"},
                             "dir": _IW8_AL}}
    _toilet_s_center = _off({"element": "IW8", "corner": "SW"},
                             {"sub": [_d_dryer_al, _in(4)]}, _IW8_AL)
    _f("toilet_s", {
        "type": "toilet_shape",
        "center": _toilet_s_center,
        "facing_dir": _IW8_IN,
        "width_dir": _IW8_AL,
    })

    # North toilet: center on IW8 north face
    _d_dryer_al_n = {"proj": {"target": _dryer_cx,
                               "anchor": {"element": "IW8", "corner": "NW"},
                               "dir": _IW8_AL}}
    _toilet_n_center = _off({"element": "IW8", "corner": "NW"},
                             {"sub": [_d_dryer_al_n, _in(4)]}, _IW8_AL)
    _f("toilet_n", {
        "type": "toilet_shape",
        "center": _toilet_n_center,
        "facing_dir": _IW8_OUT,
        "width_dir": _IW8_AL,
    })

    # ===================================================================
    # UTIL SINK — all variants, on IW8 south face
    # ===================================================================
    # _sink_mid = midpoint(dryer_centroid, counter_centroid)
    # _d_sink_al = proj(_sink_mid, IW8.SW, iw8_al)
    _ctr_cx = {"element_centroid": "COUNTER"}
    _sink_mid = {"midpoint": [_dryer_cx, _ctr_cx]}
    _d_sink_al = {"proj": {"target": _sink_mid,
                            "anchor": {"element": "IW8", "corner": "SW"},
                            "dir": _IW8_AL}}
    _sink_s_pt = _off({"element": "IW8", "corner": "SW"}, _d_sink_al, _IW8_AL)
    # ellipse center = offset(_sink_s_pt, SINK_RY, iw8_in)
    _sink_center = _off(_sink_s_pt, _c("SINK_RY"), _IW8_IN)
    _f("util_sink", {
        "type": "ellipse_rect",
        "anchor": _sink_center,
        "along": _IW8_AL,
        "outward": _IW8_IN,
        "rx": _c("SINK_RX"),
        "ry": _c("SINK_RY"),
    })

    # ===================================================================
    # BATH SINK — all variants, on IW8 north face
    # ===================================================================
    # _d_iw2_al = proj(IW2.SW, IW8.NW, iw8_al)
    # _bath_sink_east_d = _d_iw2_al - 9/12
    # _bath_sink_ctr_d = _bath_sink_east_d - BATH_SINK_LENGTH/2
    _d_iw2_al = {"proj": {"target": {"element": "IW2", "corner": "SW"},
                            "anchor": {"element": "IW8", "corner": "NW"},
                            "dir": _IW8_AL}}
    _bs_ctr_d = {"sub": [{"sub": [_d_iw2_al, _in(9)]},
                          {"mul": [0.5, _c("BATH_SINK_LENGTH")]}]}
    _bs_anchor = _off({"element": "IW8", "corner": "NW"}, _bs_ctr_d, _IW8_AL)
    _f("bath_sink", {
        "type": "bath_sink_shape",
        "anchor": _bs_anchor,
        "along": _IW8_AL,
        "outward": _IW8_OUT,
        "length": _c("BATH_SINK_LENGTH"),
        "depth": _c("BATH_SINK_DEPTH"),
    })

    # ===================================================================
    # KITCHEN SINK — all variants
    # ===================================================================
    _f("kitchen_sink", {
        "type": "item_rect",
        "anchor": _nwp(_ks_d, _c("KITCHEN_SINK_DEPTH")),
        "along": _W9W10_AL,
        "across": {"neg": _W9W10_IN},
        "width": _c("KITCHEN_SINK_WIDTH"),
        "depth": _c("KITCHEN_SINK_DEPTH"),
        "anchor_corner": "sw",
    })

    # ===================================================================
    # STOVE — standard + daybed only
    # ===================================================================
    _stove_formula = {
        "type": "item_rect",
        "anchor": _nwp(_st_d, {"add": [_c("KITCHEN_APPL_GAP"), _c("STOVE_DEPTH")]}),
        "along": _W9W10_AL,
        "across": {"neg": _W9W10_IN},
        "width": _c("STOVE_WIDTH"),
        "depth": _c("STOVE_DEPTH"),
        "anchor_corner": "sw",
    }
    _f("stove", _stove_formula, "standard")
    _f("stove", _stove_formula, "daybed")

    # ===================================================================
    # DISHWASHER — standard + daybed only
    # ===================================================================
    _dw_formula = {
        "type": "item_rect",
        "anchor": _nwp(_dw_d, _c("DW_DEPTH")),
        "along": _W9W10_AL,
        "across": {"neg": _W9W10_IN},
        "width": _c("DW_WIDTH"),
        "depth": _c("DW_DEPTH"),
        "anchor_corner": "sw",
    }
    _f("dishwasher", _dw_formula, "standard")
    _f("dishwasher", _dw_formula, "daybed")

    # ===================================================================
    # FRIDGE — different positioning per variant
    # ===================================================================
    # Standard fridge: near IW12 corner
    _std_fridge = {
        "type": "four_corner",
        "sw": _iwp(_c("KITCHEN_APPL_GAP"), _c("KITCHEN_APPL_GAP")),
        "se": _iwp({"add": [_c("KITCHEN_APPL_GAP"), _c("STD_FRIDGE_W")]},
                    _c("KITCHEN_APPL_GAP")),
        "ne": _iwp({"add": [_c("KITCHEN_APPL_GAP"), _c("STD_FRIDGE_W")]},
                    {"add": [_c("KITCHEN_APPL_GAP"), _c("STD_FRIDGE_D")]}),
        "nw": _iwp(_c("KITCHEN_APPL_GAP"),
                    {"add": [_c("KITCHEN_APPL_GAP"), _c("STD_FRIDGE_D")]}),
    }
    _f("fridge", _std_fridge, "standard")
    _f("fridge", _std_fridge, "daybed")

    # Minik fridge: on north wall
    _fr_mk_d = {"add": [_ks_d, _c("KITCHEN_SINK_WIDTH"), _in(3)]}
    _fr_mk_i = _in(3)
    _f("fridge", {
        "type": "item_rect",
        "anchor": _nwp(_fr_mk_d, {"add": [_fr_mk_i, _c("MINIK_FRIDGE_D")]}),
        "along": _W9W10_AL,
        "across": {"neg": _W9W10_IN},
        "width": _c("MINIK_FRIDGE_W"),
        "depth": _c("MINIK_FRIDGE_D"),
        "anchor_corner": "sw",
    }, "minik")

    # ===================================================================
    # ICE MAKER — all variants, different positions
    # ===================================================================
    # Standard: _dw_d + DW_WIDTH + 6/12, gap = KITCHEN_APPL_GAP
    _ice_d_std = {"add": [_dw_d, _c("DW_WIDTH"), _in(6)]}
    _ice_i_std = _c("KITCHEN_APPL_GAP")
    _f("ice_maker", {
        "type": "item_rect",
        "anchor": _nwp(_ice_d_std, {"add": [_ice_i_std, _c("ICE_DEPTH")]}),
        "along": _W9W10_AL,
        "across": {"neg": _W9W10_IN},
        "width": _c("ICE_WIDTH"),
        "depth": _c("ICE_DEPTH"),
        "anchor_corner": "sw",
    }, "standard")

    # Daybed: _dw_d + DW_WIDTH + 2/12, gap = KITCHEN_APPL_GAP
    _ice_d_db = {"add": [_dw_d, _c("DW_WIDTH"), _in(2)]}
    _f("ice_maker", {
        "type": "item_rect",
        "anchor": _nwp(_ice_d_db, {"add": [_ice_i_std, _c("ICE_DEPTH")]}),
        "along": _W9W10_AL,
        "across": {"neg": _W9W10_IN},
        "width": _c("ICE_WIDTH"),
        "depth": _c("ICE_DEPTH"),
        "anchor_corner": "sw",
    }, "daybed")

    # Minik: after fridge
    _ice_d_mk = {"add": [_fr_mk_d, _c("MINIK_FRIDGE_W"), _in(3)]}
    _ice_i_mk = _in(3)
    _f("ice_maker", {
        "type": "item_rect",
        "anchor": _nwp(_ice_d_mk, {"add": [_ice_i_mk, _c("ICE_DEPTH")]}),
        "along": _W9W10_AL,
        "across": {"neg": _W9W10_IN},
        "width": _c("ICE_WIDTH"),
        "depth": _c("ICE_DEPTH"),
        "anchor_corner": "sw",
    }, "minik")

    # ===================================================================
    # WORK COUNTER — standard + daybed only
    # ===================================================================
    _wc_d2 = {"add": [_c("KITCHEN_APPL_GAP"), _c("STD_FRIDGE_W"),
                       _c("KITCHEN_APPL_GAP")]}
    _wc_formula = {
        "type": "four_corner",
        "sw": _iwp(_wc_d2, 0),
        "se": _iwp({"add": [_wc_d2, _c("WORK_CTR_W")]}, 0),
        "ne": _iwp({"add": [_wc_d2, _c("WORK_CTR_W")]}, _c("WORK_CTR_D")),
        "nw": _iwp(_wc_d2, _c("WORK_CTR_D")),
    }
    _f("work_counter", _wc_formula, "standard")
    _f("work_counter", _wc_formula, "daybed")

    # ===================================================================
    # MICROWAVE — different per variant
    # ===================================================================
    # Standard: on work counter
    _mw_d2_std = {"add": [_c("KITCHEN_APPL_GAP"), _c("STD_FRIDGE_W"),
                           _c("KITCHEN_APPL_GAP"), _in(2)]}
    _mw_d1_std = _in(2)
    _mw_std = {
        "type": "four_corner",
        "sw": _iwp(_mw_d2_std, _mw_d1_std),
        "se": _iwp({"add": [_mw_d2_std, _c("MICROWAVE_W")]}, _mw_d1_std),
        "ne": _iwp({"add": [_mw_d2_std, _c("MICROWAVE_W")]},
                    {"add": [_mw_d1_std, _c("MICROWAVE_D")]}),
        "nw": _iwp(_mw_d2_std, {"add": [_mw_d1_std, _c("MICROWAVE_D")]}),
    }
    _f("microwave", _mw_std, "standard")
    _f("microwave", _mw_std, "daybed")

    # Minik: on north wall
    _mw_mk_d = {"add": [_iw2_d, _in(2)]}
    _mw_mk_i = _in(3)
    _f("microwave", {
        "type": "item_rect",
        "anchor": _nwp(_mw_mk_d, {"add": [_mw_mk_i, _c("MICROWAVE_D")]}),
        "along": _W9W10_AL,
        "across": {"neg": _W9W10_IN},
        "width": _c("MICROWAVE_W"),
        "depth": _c("MICROWAVE_D"),
        "anchor_corner": "sw",
    }, "minik")

    # ===================================================================
    # KITCHEN COUNTER (minik) / NORTH COUNTER (standard/daybed)
    # ===================================================================
    _f("kitchen_counter", {
        "type": "item_rect",
        "anchor": _nwp(_iw2_d, _c("KITCHEN_CTR_DEPTH")),
        "along": _W9W10_AL,
        "across": {"neg": _W9W10_IN},
        "width": _c("KITCHEN_CTR_LENGTH"),
        "depth": _c("KITCHEN_CTR_DEPTH"),
        "anchor_corner": "sw",
    }, "minik")

    _nc_formula = {
        "type": "item_rect",
        "anchor": _nwp(_iw2_d, _c("NORTH_CTR_DEPTH")),
        "along": _W9W10_AL,
        "across": {"neg": _W9W10_IN},
        "width": _c("NORTH_CTR_LENGTH"),
        "depth": _c("NORTH_CTR_DEPTH"),
        "anchor_corner": "sw",
    }
    _f("north_counter", _nc_formula, "standard")
    _f("north_counter", _nc_formula, "daybed")

    # ===================================================================
    # COFFEE MAKER — all variants, different positions
    # ===================================================================
    # Standard: on north counter near east end
    _cm_d_std = {"sub": [{"add": [_iw2_d, _c("NORTH_CTR_LENGTH")]},
                          {"add": [_in(2), _c("COFFEE_W")]}]}
    _cm_i_std = _in(2)
    _cm_std = {
        "type": "item_rect",
        "anchor": _nwp(_cm_d_std, {"add": [_cm_i_std, _c("COFFEE_D")]}),
        "along": _W9W10_AL,
        "across": {"neg": _W9W10_IN},
        "width": _c("COFFEE_W"),
        "depth": _c("COFFEE_D"),
        "anchor_corner": "sw",
    }
    _f("coffee_maker", _cm_std, "standard")
    _f("coffee_maker", _cm_std, "daybed")

    # Minik: after microwave
    _cm_d_mk = {"add": [_iw2_d, _in(2), _c("MICROWAVE_W"), _in(3)]}
    _cm_i_mk = _in(3)
    _f("coffee_maker", {
        "type": "item_rect",
        "anchor": _nwp(_cm_d_mk, {"add": [_cm_i_mk, _c("COFFEE_D")]}),
        "along": _W9W10_AL,
        "across": {"neg": _W9W10_IN},
        "width": _c("COFFEE_W"),
        "depth": _c("COFFEE_D"),
        "anchor_corner": "sw",
    }, "minik")

    # ===================================================================
    # COOKTOP — minik only
    # ===================================================================
    _cp_d = {"add": [_iw2_d, _in(2), _c("MICROWAVE_W"), _in(3),
                      _c("COFFEE_W"), _in(3)]}
    _cp_i_far = {"sub": [_c("KITCHEN_CTR_DEPTH"), _in(2)]}
    _f("cooktop", {
        "type": "item_rect",
        "anchor": _nwp(_cp_d, _cp_i_far),
        "along": _W9W10_AL,
        "across": {"neg": _W9W10_IN},
        "width": _c("COOKTOP_W"),
        "depth": _c("COOKTOP_D"),
        "anchor_corner": "sw",
    }, "minik")

    # ===================================================================
    # TOASTER — minik only
    # ===================================================================
    _ts_d = {"add": [_cp_d, _c("COOKTOP_W"), _in(3)]}
    _ts_i = _in(3)
    _f("toaster", {
        "type": "item_rect",
        "anchor": _nwp(_ts_d, {"add": [_ts_i, _c("TOASTER_D")]}),
        "along": _W9W10_AL,
        "across": {"neg": _W9W10_IN},
        "width": _c("TOASTER_W"),
        "depth": _c("TOASTER_D"),
        "anchor_corner": "sw",
    }, "minik")

    # ===================================================================
    # DINING TABLE — all variants
    # ===================================================================
    # Table base center computed from kitchen reference point and IW12 corner
    # _space_n_ref for minik = _nwp(0, KITCHEN_CTR_DEPTH), else _nwp(0, 0) = W9
    # _tbl_ref for minik = _nwp(_st_d + STOVE_WIDTH + KITCHEN_APPL_GAP)
    # _tbl_ref for standard/db = _nwp(_ks_d + KITCHEN_SINK_WIDTH/2)
    # _tbl_d_al = proj(_tbl_ref, _iw12_corner, iw1_n_al)
    # _space_n_d_out = proj(_space_n_ref, _iw12_corner, iw1_n_out)
    # _tbl_n_offset = 30/12 + (28/12 if not minik else 0)
    # _tbl_n_d_out = _space_n_d_out - _tbl_n_offset
    # tbl_bc = _iw12_corner + _tbl_d_al * iw1_n_al + _tbl_n_d_out * iw1_n_out
    # toward_apex = -iw1_n_out (southward)
    # along_base = iw1_n_al (eastward)

    # Standard/daybed table
    _tbl_ref_std = _nwp({"add": [_ks_d, {"mul": [0.5, _c("KITCHEN_SINK_WIDTH")]}]})
    _tbl_d_al_std = {"proj": {"target": _tbl_ref_std,
                               "anchor": _IW12_CORNER, "dir": _IW1_N_AL}}
    _space_n_ref_std = "W9"
    _space_n_d_out_std = {"proj": {"target": _space_n_ref_std,
                                    "anchor": _IW12_CORNER, "dir": _IW1_N_OUT}}
    _tbl_n_d_out_std = {"sub": [_space_n_d_out_std, {"add": [_in(30), _in(28)]}]}
    _tbl_bc_std = _off(_off(_IW12_CORNER, _tbl_d_al_std, _IW1_N_AL),
                        _tbl_n_d_out_std, _IW1_N_OUT)
    _tbl_std = {
        "type": "dining_triangle",
        "base_center": _tbl_bc_std,
        "toward_apex": {"neg": _IW1_N_OUT},
        "along_base": _IW1_N_AL,
        "base_width": _c("DINING_TBL_BASE"),
        "height": _c("DINING_TBL_H"),
        "apex_radius": _in(12),
        "fillet_radius": _in(6),
    }
    _f("dining_table", _tbl_std, "standard")
    _f("dining_table", _tbl_std, "daybed")

    # Minik table
    _tbl_ref_mk = _nwp({"add": [_st_d, _c("STOVE_WIDTH"), _c("KITCHEN_APPL_GAP")]})
    _tbl_d_al_mk = {"proj": {"target": _tbl_ref_mk,
                               "anchor": _IW12_CORNER, "dir": _IW1_N_AL}}
    _space_n_ref_mk = _nwp(0, _c("KITCHEN_CTR_DEPTH"))
    _space_n_d_out_mk = {"proj": {"target": _space_n_ref_mk,
                                    "anchor": _IW12_CORNER, "dir": _IW1_N_OUT}}
    _tbl_n_d_out_mk = {"sub": [_space_n_d_out_mk, _in(30)]}
    _tbl_bc_mk = _off(_off(_IW12_CORNER, _tbl_d_al_mk, _IW1_N_AL),
                        _tbl_n_d_out_mk, _IW1_N_OUT)
    _f("dining_table", {
        "type": "dining_triangle",
        "base_center": _tbl_bc_mk,
        "toward_apex": {"neg": _IW1_N_OUT},
        "along_base": _IW1_N_AL,
        "base_width": _c("DINING_TBL_BASE"),
        "height": _c("DINING_TBL_H"),
        "apex_radius": _in(12),
        "fillet_radius": _in(6),
    }, "minik")

    # ===================================================================
    # DINING CHAIRS — all variants
    # ===================================================================
    _dc_common = {
        "type": "dining_chair",
        "table": "dining_table",
        "chair_width": _c("DINING_CHAIR_W"),
        "chair_depth": _c("DINING_CHAIR_D"),
        "gap": _in(2),
    }
    _dc1 = dict(_dc_common, side="ne_side")
    _dc2 = dict(_dc_common, side="nw_side")
    _f("dining_chair_1", _dc1)
    _f("dining_chair_2", _dc2)

    # ===================================================================
    # CHAIR + OTTOMAN — all variants
    # ===================================================================
    # chair orientation: W12-W13 direction rotated -45°
    # _ch_theta = atan2(w12w13_al[1], w12w13_al[0]) - pi/4
    # For the formula, we compute the direction from the segment and rotate
    # _ch_mid = midpoint(W11, W12)
    # _ch_base = offset(offset(_ch_mid, -1/12, w11w12_al), 8/12, w11w12_in)
    # ch_center = offset(_ch_base, 4/12, _ch_along)
    #
    # We use item_rect with computed along/across directions.
    # The along direction is W12→W13 rotated -45° around origin.
    # In formula spec: we need a "rotate" direction spec.
    # Instead, we'll compute the rotated directions directly.
    # _ch_along = (cos(theta), sin(theta)) where theta = angle(W12→W13) - pi/4
    # = cos(-pi/4) * w12w13_al - sin(-pi/4) * w12w13_perp
    # = (al + perp) / sqrt(2)   where perp = (-al[1], al[0])
    # = (al + left_perp) / sqrt(2)
    # segment_perp gives RIGHT perp (dy, -dx), left perp is (-dy, dx) = neg(segment_perp)
    # Actually: _ch_along = cos(-pi/4)*al + sin(-pi/4)*perp_right
    # = al/sqrt(2) - perp_right/sqrt(2)
    # where perp_right = segment_perp(W12, W13)
    # Hmm, this is complex. Let me use a different approach: store the
    # rotated direction as a point-based computation.
    #
    # Actually, item_rect with along/across works for any orientation.
    # The challenge is computing the center point and the direction vectors.
    # For the chair, the center is:
    #   _ch_mid = midpoint(W11, W12)
    #   _ch_base = offset(offset(_ch_mid, -1/12, w11w12_al), 8/12, w11w12_in)
    #   ch_center = offset(_ch_base, 4/12, _ch_along)
    # The along direction is rotated, which we can't easily express.
    #
    # Let's use four_corner with explicit corner computation.
    # ch_center = known, _ch_along and _ch_cross are known directions
    # corners at center ± half_w * cross ± half_d * along
    #
    # For the direction: _ch_theta = atan2(w12w13_al[1], w12w13_al[0]) - pi/4
    # This is the angle of W12→W13 minus 45°.
    # The unit vector is: (cos(theta), sin(theta))
    # = cos(angle(w12w13) - pi/4), sin(angle(w12w13) - pi/4)
    # = cos(a)*cos(pi/4) + sin(a)*sin(pi/4), sin(a)*cos(pi/4) - cos(a)*sin(pi/4)
    # = (al[0] + al[1])/sqrt(2), (al[1] - al[0])/sqrt(2)
    # where al = w12w13_al = (cos(a), sin(a))
    #
    # Actually w12w13_al = (cos(a), sin(a)), so:
    # cos(a - pi/4) = cos(a)*cos(pi/4) + sin(a)*sin(pi/4)
    #               = (cos(a) + sin(a)) / sqrt(2)
    #               = (al[0] + al[1]) / sqrt(2)
    # sin(a - pi/4) = sin(a)*cos(pi/4) - cos(a)*sin(pi/4)
    #               = (sin(a) - cos(a)) / sqrt(2)
    #               = (al[1] - al[0]) / sqrt(2)
    #
    # In formula spec, I can express this as:
    # _ch_along = normalize((al[0]+al[1], al[1]-al[0]))
    # Since al is already a unit vector, this vector has length 1.
    #
    # But there's no "normalize" or "add vectors" in the formula spec.
    # Let me add an "add_dirs" operation... or use a different approach.
    #
    # Simplest: compute the four corners explicitly using four_corner formula.
    # Each corner = center + ds * half_short * cross + dn * half_depth * along
    # We need to express this in terms of W12→W13 direction.
    #
    # Let me precompute the chair center and use literal offsets.
    # Actually, let me add a "rotate_dir" spec to the direction resolver.

    # For now, I'll use a simpler approach: compute the four corners
    # using the known relationship between W12-W13 direction and the
    # 45° rotation.

    _ch_mid = {"midpoint": ["W11", "W12"]}
    _ch_base = _off(_off(_ch_mid, _in(-1), _W11W12_AL), _in(8), _W11W12_IN)

    # The rotated direction can be expressed as a normalized sum of
    # w12w13_al and neg(w12w13_perp_right)
    # ch_along ≈ (al + left_perp) / sqrt(2)
    # where left_perp = neg(segment_perp)
    # In formula spec, I need a "rotate" or I can use four_corner
    # with explicit corner points computed from the base geometry.
    #
    # Let me add a "rotated" direction spec: rotate a direction by angle
    _f("chair", {
        "type": "four_corner",
        "sw": _off(_off(
            _off(_ch_base, _in(4),
                 {"rotated": _W12W13_AL, "angle": -0.7853981633974483}),
            {"neg": {"mul": [0.5, _c("CHAIR_WIDTH")]}},
            {"perp": {"rotated": _W12W13_AL, "angle": -0.7853981633974483}}),
            {"neg": {"mul": [0.5, _c("CHAIR_DEPTH")]}},
            {"rotated": _W12W13_AL, "angle": -0.7853981633974483}),
        "se": _off(_off(
            _off(_ch_base, _in(4),
                 {"rotated": _W12W13_AL, "angle": -0.7853981633974483}),
            {"mul": [0.5, _c("CHAIR_WIDTH")]},
            {"perp": {"rotated": _W12W13_AL, "angle": -0.7853981633974483}}),
            {"neg": {"mul": [0.5, _c("CHAIR_DEPTH")]}},
            {"rotated": _W12W13_AL, "angle": -0.7853981633974483}),
        "ne": _off(_off(
            _off(_ch_base, _in(4),
                 {"rotated": _W12W13_AL, "angle": -0.7853981633974483}),
            {"mul": [0.5, _c("CHAIR_WIDTH")]},
            {"perp": {"rotated": _W12W13_AL, "angle": -0.7853981633974483}}),
            {"mul": [0.5, _c("CHAIR_DEPTH")]},
            {"rotated": _W12W13_AL, "angle": -0.7853981633974483}),
        "nw": _off(_off(
            _off(_ch_base, _in(4),
                 {"rotated": _W12W13_AL, "angle": -0.7853981633974483}),
            {"neg": {"mul": [0.5, _c("CHAIR_WIDTH")]}},
            {"perp": {"rotated": _W12W13_AL, "angle": -0.7853981633974483}}),
            {"mul": [0.5, _c("CHAIR_DEPTH")]},
            {"rotated": _W12W13_AL, "angle": -0.7853981633974483}),
    })

    # ===================================================================
    # OTTOMAN — all variants
    # ===================================================================
    _ot_dist = _in(39)
    _ch_along_spec = {"rotated": _W12W13_AL, "angle": -0.7853981633974483}
    _ch_cross_spec = {"perp": _ch_along_spec}
    _ch_center = _off(_ch_base, _in(4), _ch_along_spec)
    _ot_center = _off(_ch_center, _ot_dist, _ch_along_spec)
    _f("ottoman", {
        "type": "four_corner",
        "sw": _off(_off(_ot_center,
                         {"neg": {"mul": [0.5, _c("OTTOMAN_SIZE")]}}, _ch_cross_spec),
                    {"neg": {"mul": [0.5, _c("OTTOMAN_SIZE")]}}, _ch_along_spec),
        "se": _off(_off(_ot_center,
                         {"mul": [0.5, _c("OTTOMAN_SIZE")]}, _ch_cross_spec),
                    {"neg": {"mul": [0.5, _c("OTTOMAN_SIZE")]}}, _ch_along_spec),
        "ne": _off(_off(_ot_center,
                         {"mul": [0.5, _c("OTTOMAN_SIZE")]}, _ch_cross_spec),
                    {"mul": [0.5, _c("OTTOMAN_SIZE")]}, _ch_along_spec),
        "nw": _off(_off(_ot_center,
                         {"neg": {"mul": [0.5, _c("OTTOMAN_SIZE")]}}, _ch_cross_spec),
                    {"mul": [0.5, _c("OTTOMAN_SIZE")]}, _ch_along_spec),
    })

    # ===================================================================
    # DESK + DESK CHAIR — all variants
    # ===================================================================
    _neg_w16w17_al = {"neg": _W16W17_AL}
    _f("desk", {
        "type": "item_rect",
        "anchor": "W17",
        "along": _neg_w16w17_al,
        "across": _W16W17_IN,
        "width": _c("DESK_WIDTH"),
        "depth": _c("DESK_DEPTH"),
        "anchor_corner": "sw",
    })

    _dc_sw = _off(
        _off("W17",
             {"sub": [{"mul": [0.5, _c("DESK_WIDTH")]},
                       {"mul": [0.5, _c("DESK_CHAIR_WIDTH")]}]},
             _neg_w16w17_al),
        {"add": [_c("DESK_DEPTH"), _c("DESK_CHAIR_GAP")]},
        _W16W17_IN,
    )
    _f("desk_chair", {
        "type": "item_rect",
        "anchor": _dc_sw,
        "along": _neg_w16w17_al,
        "across": _W16W17_IN,
        "width": _c("DESK_CHAIR_WIDTH"),
        "depth": _c("DESK_CHAIR_DEPTH"),
        "anchor_corner": "sw",
    })

    # ===================================================================
    # VARIANT-SPECIFIC SEATING
    # ===================================================================

    # --- Standard: loveseat, ET, loveseat2 ---
    _lv_perp = {"perp": _W12W13_AL}  # right perp = (-w12w13[1], w12w13[0])
    # But the procedural code uses (-w12w13_al[1], w12w13_al[0]) which is
    # the LEFT perpendicular. In evaluator, "perp" gives RIGHT perp (dy, -dx).
    # LEFT perp = neg(perp) = {"neg": {"perp": ...}}
    # Actually: for w12w13_al = (dx/L, dy/L):
    #   right perp (segment_perp) = (dy/L, -dx/L)
    #   left perp = (-dy/L, dx/L)
    # Procedural uses _lv_perp = (-w12w13_al[1], w12w13_al[0]) = (-dy/L, dx/L) = left perp
    # So _lv_perp = {"neg": {"segment_perp": ["W12", "W13"]}}
    _lv_perp_spec = {"neg": {"segment_perp": ["W12", "W13"]}}

    _lv_nw = _lwp(_c("LOVESEAT_OFFSET_IW4"), _c("LOVESEAT_OFFSET_IW1"))
    _lv_formula = {
        "type": "four_corner",
        "nw": _lv_nw,
        "sw": _off(_lv_nw, _c("LOVESEAT_LENGTH"), _W12W13_AL),
        "ne": _off(_lv_nw, _c("LOVESEAT_WIDTH"), _lv_perp_spec),
        "se": _off(_off(_lv_nw, _c("LOVESEAT_LENGTH"), _W12W13_AL),
                    _c("LOVESEAT_WIDTH"), _lv_perp_spec),
    }
    _f("loveseat", _lv_formula, "standard")

    # ET (standard): circle tangent to loveseat SE, on IW1 N face offset line
    _et_r_ft = {"mul": [{"mul": [_c("ET_RADIUS_CM"), 1.0 / 2.54]}, 1.0 / 12.0]}
    _et_from_iw1 = _off({"element": "IW1", "corner": "NW"},
                         {"add": [_c("STD_GAP"), _et_r_ft]}, _IW1_N_OUT)
    # Need quadratic solve for tangency to loveseat SE corner.
    # et_gap = et_r + STD_GAP
    # Solve: |_et_from_iw1 + t * iw1_n_al - loveseat.SE|² = et_gap²
    _lv_se_spec = _off(_off(_lv_nw, _c("LOVESEAT_LENGTH"), _W12W13_AL),
                        _c("LOVESEAT_WIDTH"), _lv_perp_spec)
    _et_gap = {"add": [_et_r_ft, _c("STD_GAP")]}
    _f("et", {
        "type": "item_circle",
        "center": {"ray_circle_isect": {
            "origin": _et_from_iw1,
            "dir": _IW1_N_AL,
            "target": _lv_se_spec,
            "distance": _et_gap,
            "select": "farthest",
        }},
        "radius": _et_r_ft,
    }, "standard")

    # LOVESEAT2 (standard): east of ET
    # _lv2_sw = offset(offset(et_center, et_r + STD_GAP, iw1_n_al), -et_r, iw1_n_out)
    _et_center = {"ray_circle_isect": {
        "origin": _et_from_iw1,
        "dir": _IW1_N_AL,
        "target": _lv_se_spec,
        "distance": _et_gap,
        "select": "farthest",
    }}
    _lv2_sw = _off(
        _off(_et_center, {"add": [_et_r_ft, _c("STD_GAP")]}, _IW1_N_AL),
        {"neg": _et_r_ft}, _IW1_N_OUT,
    )
    _f("loveseat2", {
        "type": "four_corner",
        "sw": _lv2_sw,
        "nw": _off(_lv2_sw, _c("LOVESEAT_WIDTH"), _IW1_N_OUT),
        "se": _off(_lv2_sw, _c("LOVESEAT_LENGTH"), _IW1_N_AL),
        "ne": _off(_off(_lv2_sw, _c("LOVESEAT_WIDTH"), _IW1_N_OUT),
                    _c("LOVESEAT_LENGTH"), _IW1_N_AL),
    }, "standard")

    # --- Minik: sofa, rocker ---
    # SOFA: near IW4/IW1 corner
    # _cx_d = -(6/12 + (SOFA_WIDTH - 24/12) / 2)
    _cx_d = {"neg": {"add": [_in(6), {"mul": [0.5, {"sub": [_c("SOFA_WIDTH"),
                                                              _in(24)]}]}]}}
    _sofa_sw = _lwp({"add": [_cx_d, {"mul": [0.5, _c("SOFA_FULL_W")]}]}, _in(2))
    _sofa_se = _lwp({"sub": [_cx_d, {"mul": [0.5, _c("SOFA_FULL_W")]}]}, _in(2))
    _f("sofa", {
        "type": "four_corner",
        "sw": _sofa_sw,
        "se": _sofa_se,
        "ne": _off(_sofa_se, _c("SOFA_FULL_D"), _IW1_N_OUT),
        "nw": _off(_sofa_sw, _c("SOFA_FULL_D"), _IW1_N_OUT),
    }, "minik")

    # ROCKER (minik): midpoint between ICE SE and SOFA NW, offset inward
    # _ice_se (minik) = _nwp(_ice_d_mk + ICE_WIDTH, 3/12 + ICE_DEPTH)
    _ice_se_mk = _nwp({"add": [_ice_d_mk, _c("ICE_WIDTH")]},
                       {"add": [_in(3), _c("ICE_DEPTH")]})
    _sofa_nw_mk = _off(_sofa_sw, _c("SOFA_FULL_D"), _IW1_N_OUT)
    _rk_mid = {"midpoint": [_ice_se_mk, _sofa_nw_mk]}
    _rk_center_mk = _off(_rk_mid, _in(18), _W9W10_IN)
    # Rocker oriented along W12-W13 direction
    _rk_cr = _lv_perp_spec  # left perp of W12→W13
    _rk_hw = {"mul": [0.5, _c("ROCKER_DEPTH")]}
    _rk_hh = {"mul": [0.5, _c("ROCKER_WIDTH")]}
    _f("rocker", {
        "type": "four_corner",
        "sw": _off(_off(_rk_center_mk, {"neg": _rk_hh}, _W12W13_AL),
                    {"neg": _rk_hw}, _rk_cr),
        "se": _off(_off(_rk_center_mk, {"neg": _rk_hh}, _W12W13_AL),
                    _rk_hw, _rk_cr),
        "ne": _off(_off(_rk_center_mk, _rk_hh, _W12W13_AL),
                    _rk_hw, _rk_cr),
        "nw": _off(_off(_rk_center_mk, _rk_hh, _W12W13_AL),
                    {"neg": _rk_hw}, _rk_cr),
    }, "minik")

    # --- Daybed variant: shelves2, et_east, daybed, et_west, rocker ---

    # SHELVES2 (daybed): KALLAX east of IW9
    _iw9_e_al = {"segment": [{"element": "IW9", "corner": "SE"},
                               {"element": "IW9", "corner": "NE"}]}
    _iw9_e_cw = {"segment_perp": [{"element": "IW9", "corner": "SE"},
                                    {"element": "IW9", "corner": "NE"}]}
    _sh2_nw = _li(
        _off(_off({"element": "IW1", "corner": "SW"}, _c("SHELVES_GAP_IW1"),
                   _W9W10_IN), 0, _W9W10_AL),
        _W9W10_AL,
        _off({"element": "IW9", "corner": "NE"}, _c("SHELVES_GAP_IW9"), _iw9_e_cw),
        _iw9_e_al,
    )
    _f("shelves2", {
        "type": "four_corner",
        "nw": _sh2_nw,
        "ne": _off(_sh2_nw, _c("SHELVES_LENGTH"), _iw9_e_cw),
        "sw": _off(_sh2_nw, _c("SHELVES_DEPTH"), _W9W10_IN),
        "se": _off(_off(_sh2_nw, _c("SHELVES_LENGTH"), _iw9_e_cw),
                    _c("SHELVES_DEPTH"), _W9W10_IN),
    }, "daybed")

    # ET_EAST (daybed): circle on IW1 N face offset, near east inner wall
    _et_from_iw1_db = _off({"element": "IW1", "corner": "NW"},
                            {"add": [_c("STD_GAP"), _et_r_ft]}, _IW1_N_OUT)
    # Find farthest intersection of the IW1-parallel line with inner_poly
    _f("et_east", {
        "type": "item_circle",
        "center": {"ray_circle_isect": {
            "origin": _et_from_iw1_db,
            "dir": _IW1_N_AL,
            "target": {"line_poly_isect": {
                "origin": _et_from_iw1_db,
                "dir": _IW1_N_AL,
                "poly": "inner_poly",
                "select": "farthest",
            }},
            "distance": {"add": [_c("STD_GAP"), _et_r_ft]},
            "select": "nearest",
        }},
        "radius": _et_r_ft,
    }, "daybed")

    # DAYBED (daybed): positioned relative to et_east
    # db_se = offset(et_east_center, et_r + 3/12, -iw1_n_al)
    # then project to IW1 N face offset line
    # db_sw = offset(db_se, -DAYBED_W, iw1_n_al)
    _neg_al = {"neg": _IW1_N_AL}
    _db_se_raw = _off({"element_centroid": "et_east"},
                       {"add": [_et_r_ft, _in(3)]}, _neg_al)
    # Project db_se onto the IW1 N face offset line (at STD_GAP from IW1.NW)
    _db_s_ref = _off({"element": "IW1", "corner": "NW"}, _c("STD_GAP"), _IW1_N_OUT)
    # _db_se_proj_out = proj(db_se_raw, _db_s_ref, iw1_n_out)
    _db_se_proj_out = {"proj": {"target": _db_se_raw, "anchor": _db_s_ref,
                                 "dir": _IW1_N_OUT}}
    _db_se = _off(_db_se_raw, {"neg": _db_se_proj_out}, _IW1_N_OUT)
    _db_sw = _off(_db_se, {"neg": _c("DAYBED_W")}, _IW1_N_AL)
    _f("daybed", {
        "type": "four_corner",
        "sw": _db_sw,
        "se": _db_se,
        "ne": _off(_db_se, _c("DAYBED_D"), _IW1_N_OUT),
        "nw": _off(_db_sw, _c("DAYBED_D"), _IW1_N_OUT),
    }, "daybed")

    # ET_WEST (daybed): circle west of daybed NW
    _neg_out = {"neg": _IW1_N_OUT}
    _et2_center = _off(
        _off({"element": "daybed", "corner": "NW"},
             {"neg": {"add": [_in(6), _et_r_ft]}}, _IW1_N_AL),
        _et_r_ft, _neg_out,
    )
    _f("et_west", {
        "type": "item_circle",
        "center": _et2_center,
        "radius": _et_r_ft,
    }, "daybed")

    # ROCKER (daybed): between daybed SW and RO1 east point
    _ro1_e_pt = _nwp({"add": [_iw2_d, _c("RO1_OFFSET_FROM_IW2"),
                                _c("IW1_RO_WIDTH")]}, 0)
    _ref_iw1 = {"element": "IW1", "corner": "NW"}
    _db_sw_d_al = {"proj": {"target": {"element": "daybed", "corner": "SW"},
                             "anchor": _ref_iw1, "dir": _IW1_N_AL}}
    _ro1_d_al = {"proj": {"target": _ro1_e_pt,
                            "anchor": _ref_iw1, "dir": _IW1_N_AL}}
    _rk_d_al_db = {"sub": [{"mul": [0.5, {"add": [_db_sw_d_al, _ro1_d_al]}]},
                            _in(8)]}
    # Rocker N offset: based on fridge south face projected
    _fr_s_pt = _nwp(0, {"add": [_in(3), _c("STD_FRIDGE_D")]})
    _fr_door_s_pt = _off(_fr_s_pt, _c("STD_FRIDGE_W"), _W9W10_IN)
    _fr_d_out = {"proj": {"target": _fr_door_s_pt,
                            "anchor": _ref_iw1, "dir": _IW1_N_OUT}}
    _rk_d_out_db = {"add": [{"mul": [0.5, _fr_d_out]}, _in(26)]}
    _rk_center_db = _off(_off(_ref_iw1, _rk_d_al_db, _IW1_N_AL),
                          _rk_d_out_db, _IW1_N_OUT)
    _f("rocker", {
        "type": "four_corner",
        "sw": _off(_off(_rk_center_db, {"neg": _rk_hh}, _W12W13_AL),
                    {"neg": _rk_hw}, _rk_cr),
        "se": _off(_off(_rk_center_db, {"neg": _rk_hh}, _W12W13_AL),
                    _rk_hw, _rk_cr),
        "ne": _off(_off(_rk_center_db, _rk_hh, _W12W13_AL),
                    _rk_hw, _rk_cr),
        "nw": _off(_off(_rk_center_db, _rk_hh, _W12W13_AL),
                    {"neg": _rk_hw}, _rk_cr),
    }, "daybed")

    return formulas
