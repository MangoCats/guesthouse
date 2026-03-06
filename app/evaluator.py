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

        return None

    def resolve_length(self, spec):
        """Resolve a length expression to a float value.

        Length specs:
          float/int          → literal value
          {"const": name}    → constant lookup
        """
        if spec is None:
            return None
        if isinstance(spec, (int, float)):
            return float(spec)
        if isinstance(spec, dict):
            if "const" in spec:
                val = self.constants.get(spec["const"])
                return float(val) if val is not None else None
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
                "depth", "radius", "origin"):
        if key in spec:
            _extract_deps_recursive(spec[key], deps)

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
