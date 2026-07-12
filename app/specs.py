"""Shared spec-resolution grammar — single source of truth.

Point, direction, and length specs are small JSON structures used in two
places: element formulas evaluated by app/evaluator.py (FormulaEvaluator) and
dimension/label anchors resolved by app/engine.py.  Both consumers previously
implemented overlapping grammars by hand; this module implements the union
once.  A *context* object supplies the named-entity lookups (points, chain
tangents, element polygons, radii, constants); the grammar itself lives here.

User-facing reference: app/FORMULAS.md (§ Point specs / Direction specs /
Length specs).  Add new spec forms here, then document them there.
"""
import math

_FACE_MAP = {"south": (0, 1), "east": (1, 2), "north": (2, 3), "west": (3, 0)}


def face_midpoint(poly, face):
    """Midpoint of a named face of a 4-vertex polygon.

    Convention: south=poly[0]→poly[1], east=poly[1]→poly[2],
    north=poly[2]→poly[3], west=poly[3]→poly[0].
    poly vertices are [E, N] lists or tuples.
    """
    indices = _FACE_MAP.get(face)
    if not indices or len(poly) < 4:
        return None
    i, j = indices
    a, b = poly[i], poly[j]
    return [(a[0] + b[0]) / 2, (a[1] + b[1]) / 2]


def face_vertices(poly, face):
    """Return (vertex_i, vertex_j) for a named face of a 4-vertex polygon."""
    indices = _FACE_MAP.get(face)
    if not indices or len(poly) < 4:
        return None, None
    i, j = indices
    return poly[i], poly[j]


def line_poly_intersections(origin, direction, poly):
    """Find intersections of a ray from origin along direction with polygon edges.

    Returns list of [E, N] intersection points sorted by distance from origin.
    """
    results = []
    n = len(poly)
    for i in range(n):
        a = poly[i]
        b = poly[(i + 1) % n]
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


class SpecContext:
    """Lookup interface for spec resolution.  All methods return None for
    unknown names.  Consumers subclass or duck-type this:

    - FormulaEvaluator backs every polygon lookup with its evaluated
      ``elements`` dict.
    - The engine's geometry-result adapter backs ``get_face_poly`` with
      interior walls and ``get_opening_poly`` with outer/rough openings
      (the historical anchor-spec namespaces).
    """

    def get_point(self, name):
        """Named geometry point → [E, N]."""
        return None

    def get_tangent(self, name):
        """Chain tangent at a named point → [dE, dN] (unit)."""
        return None

    def get_face_poly(self, name):
        """Polygon for face_mid / face_along / face_perp specs."""
        return None

    def get_element_poly(self, name):
        """Polygon for element / element_centroid specs."""
        return None

    def get_opening_poly(self, name):
        """Polygon for opening_face_mid / opening_centroid specs."""
        return None

    def get_radius(self, key):
        """Arc radius by key (e.g. 'R_a5')."""
        return None

    def get_constant(self, name):
        """Named constant value (feet)."""
        return None

    def get_boundary_poly(self, name):
        """Named boundary polygon (currently 'inner_poly')."""
        return None


def _unit(dx, dy):
    length = math.sqrt(dx * dx + dy * dy)
    if length < 1e-12:
        return None
    return [dx / length, dy / length]


def resolve_point(spec, ctx):
    """Resolve a point spec to [E, N] coordinates (None if unresolvable).

    Spec forms (see app/FORMULAS.md § Point specs):
      [E, N] literal; "NAME"; {ref, dAlong, dPerp}; {element, corner};
      {element_centroid}; {face_mid, face}; {opening_face_mid, face};
      {opening_centroid}; {midpoint: [A, B]}; {offset, dir, dist};
      {type: line_intersection, line1_point/dir, line2_point/dir};
      {line_poly_isect}; {arc_point}; {ray_circle_isect}.
    """
    if spec is None:
        return None

    # List/tuple passthrough: already [E, N]
    if isinstance(spec, (list, tuple)) and len(spec) == 2:
        return [float(spec[0]), float(spec[1])]

    # String: named point
    if isinstance(spec, str):
        return ctx.get_point(spec)

    if not isinstance(spec, dict):
        return None

    # Parametric reference with local-frame offset: {ref, dAlong, dPerp}
    # dAlong is along the chain tangent at the reference point; dPerp is
    # 90° CCW from that.
    if "ref" in spec:
        ref_pt = resolve_point(spec["ref"], ctx)
        if ref_pt is None:
            return None
        d_along = float(spec.get("dAlong", 0.0))
        d_perp = float(spec.get("dPerp", 0.0))
        if d_along == 0.0 and d_perp == 0.0:
            return list(ref_pt)
        ref_name = spec["ref"] if isinstance(spec["ref"], str) else None
        along = (ctx.get_tangent(ref_name) if ref_name else None) or [1.0, 0.0]
        perp = [-along[1], along[0]]
        return [
            ref_pt[0] + d_along * along[0] + d_perp * perp[0],
            ref_pt[1] + d_along * along[1] + d_perp * perp[1],
        ]

    # Element corner reference
    if "element" in spec:
        poly = ctx.get_element_poly(spec["element"])
        if not poly:
            return None
        corner = spec.get("corner")
        if isinstance(corner, int):
            if 0 <= corner < len(poly):
                return list(poly[corner])
            return None
        if isinstance(corner, str):
            corner_map = {"SW": 0, "SE": 1, "NE": 2, "NW": 3,
                          "sw": 0, "se": 1, "ne": 2, "nw": 3}
            idx = corner_map.get(corner)
            if idx is not None and idx < len(poly):
                return list(poly[idx])
            return None
        return None

    # Element centroid (vertex average)
    if "element_centroid" in spec:
        poly = ctx.get_element_poly(spec["element_centroid"])
        if not poly:
            return None
        n = len(poly)
        return [sum(p[0] for p in poly) / n, sum(p[1] for p in poly) / n]

    # Face midpoint
    if "face_mid" in spec:
        poly = ctx.get_face_poly(spec["face_mid"])
        if poly:
            return face_midpoint(poly, spec.get("face"))
        return None

    # Opening face midpoint
    if "opening_face_mid" in spec:
        poly = ctx.get_opening_poly(spec["opening_face_mid"])
        if poly:
            return face_midpoint(poly, spec.get("face"))
        return None

    # Opening centroid (vertex average)
    if "opening_centroid" in spec:
        poly = ctx.get_opening_poly(spec["opening_centroid"])
        if poly and len(poly) >= 3:
            n = len(poly)
            return [sum(p[0] for p in poly) / n, sum(p[1] for p in poly) / n]
        return None

    # Midpoint of two specs
    if "midpoint" in spec:
        pair = spec["midpoint"]
        if len(pair) != 2:
            return None
        a = resolve_point(pair[0], ctx)
        b = resolve_point(pair[1], ctx)
        if a and b:
            return [(a[0] + b[0]) / 2, (a[1] + b[1]) / 2]
        return None

    # Offset: base + dist * direction
    if "offset" in spec:
        base = resolve_point(spec["offset"], ctx)
        direction = resolve_dir(spec.get("dir"), ctx)
        dist = resolve_length(spec.get("dist", 0), ctx)
        if base and direction and dist is not None:
            return [base[0] + dist * direction[0],
                    base[1] + dist * direction[1]]
        return None

    # Line intersection
    if spec.get("type") == "line_intersection":
        p1 = resolve_point(spec.get("line1_point"), ctx)
        d1 = resolve_dir(spec.get("line1_dir"), ctx)
        p2 = resolve_point(spec.get("line2_point"), ctx)
        d2 = resolve_dir(spec.get("line2_dir"), ctx)
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
        origin = resolve_point(lpi.get("origin"), ctx)
        direction = resolve_dir(lpi.get("dir"), ctx)
        poly = ctx.get_boundary_poly(lpi.get("poly", "inner_poly")) or []
        select = lpi.get("select", "nearest")
        if origin and direction and poly:
            isects = line_poly_intersections(origin, direction, poly)
            if isects:
                return isects[0] if select == "nearest" else isects[-1]
        return None

    # Arc point: point on arc circle at the reference's northing
    if "arc_point" in spec:
        ap = spec["arc_point"]
        center = resolve_point(ap.get("center"), ctx)
        radius = ctx.get_radius(ap.get("radius_key"))
        ref = resolve_point(ap.get("reference"), ctx)
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

    # Ray-circle intersection: point on ray at given distance from target
    if "ray_circle_isect" in spec:
        rci = spec["ray_circle_isect"]
        origin = resolve_point(rci.get("origin"), ctx)
        direction = resolve_dir(rci.get("dir"), ctx)
        target = resolve_point(rci.get("target"), ctx)
        distance = resolve_length(rci.get("distance"), ctx)
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


def resolve_dir(spec, ctx):
    """Resolve a direction spec to a [dE, dN] unit vector (None if unresolvable).

    Spec forms (see app/FORMULAS.md § Direction specs):
      [dE, dN] literal; "east"/"west"/"north"/"south"; {ref_tangent,
      angle_offset_deg}; {face_along, face}; {face_perp, face};
      {segment: [A, B]}; {segment_perp: [A, B]}; {neg}; {rotated, angle};
      {perp}.
    """
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

    # Reference-tangent direction, optionally rotated by angle_offset_deg
    if "ref_tangent" in spec:
        along = list(ctx.get_tangent(spec["ref_tangent"]) or [1.0, 0.0])
        angle_offset = float(spec.get("angle_offset_deg", 0.0))
        if abs(angle_offset) > 1e-9:
            rad = math.radians(angle_offset)
            cos_a, sin_a = math.cos(rad), math.sin(rad)
            along = [along[0] * cos_a - along[1] * sin_a,
                     along[0] * sin_a + along[1] * cos_a]
        return along

    # Face along direction
    if "face_along" in spec:
        poly = ctx.get_face_poly(spec["face_along"])
        if poly:
            a, b = face_vertices(poly, spec.get("face"))
            if a is not None:
                return _unit(b[0] - a[0], b[1] - a[1])
        return None

    # Face perpendicular (right-hand)
    if "face_perp" in spec:
        poly = ctx.get_face_poly(spec["face_perp"])
        if poly:
            a, b = face_vertices(poly, spec.get("face"))
            if a is not None:
                d = _unit(b[0] - a[0], b[1] - a[1])
                return [d[1], -d[0]] if d else None
        return None

    # Segment direction (between two points)
    if "segment" in spec:
        pair = spec["segment"]
        if len(pair) != 2:
            return None
        a = resolve_point(pair[0], ctx)
        b = resolve_point(pair[1], ctx)
        if a and b:
            return _unit(b[0] - a[0], b[1] - a[1])
        return None

    # Segment perpendicular (right-hand)
    if "segment_perp" in spec:
        pair = spec["segment_perp"]
        if len(pair) != 2:
            return None
        a = resolve_point(pair[0], ctx)
        b = resolve_point(pair[1], ctx)
        if a and b:
            d = _unit(b[0] - a[0], b[1] - a[1])
            return [d[1], -d[0]] if d else None
        return None

    # Negate direction
    if "neg" in spec:
        d = resolve_dir(spec["neg"], ctx)
        if d:
            return [-d[0], -d[1]]
        return None

    # Rotate direction by angle (radians, CCW)
    if "rotated" in spec:
        d = resolve_dir(spec["rotated"], ctx)
        angle = spec.get("angle", 0)
        if d and isinstance(angle, (int, float)):
            cos_a = math.cos(angle)
            sin_a = math.sin(angle)
            return [d[0] * cos_a - d[1] * sin_a,
                    d[0] * sin_a + d[1] * cos_a]
        return None

    # Perpendicular of direction (rotate 90° CW: [dy, -dx])
    if "perp" in spec:
        d = resolve_dir(spec["perp"], ctx)
        if d:
            return [d[1], -d[0]]
        return None

    return None


def resolve_length(spec, ctx):
    """Resolve a length expression to a float value (None if unresolvable).

    Spec forms (see app/FORMULAS.md § Length specs):
      number literal; {const}; {radius_key}; {neg}; {add: [...]};
      {sub: [a, b]}; {mul: [...]}; {div: [a, b]}; {proj: {target, anchor,
      dir}}; {dist: [pt1, pt2]}.
    """
    if spec is None:
        return None
    if isinstance(spec, (int, float)):
        return float(spec)
    if isinstance(spec, dict):
        if "const" in spec:
            val = ctx.get_constant(spec["const"])
            return float(val) if val is not None else None
        if "radius_key" in spec:
            val = ctx.get_radius(spec["radius_key"])
            return float(val) if val is not None else None
        if "neg" in spec:
            val = resolve_length(spec["neg"], ctx)
            return -val if val is not None else None
        if "add" in spec:
            vals = [resolve_length(s, ctx) for s in spec["add"]]
            if any(v is None for v in vals):
                return None
            return sum(vals)
        if "sub" in spec:
            pair = spec["sub"]
            if len(pair) != 2:
                return None
            a = resolve_length(pair[0], ctx)
            b = resolve_length(pair[1], ctx)
            if a is None or b is None:
                return None
            return a - b
        if "mul" in spec:
            vals = [resolve_length(s, ctx) for s in spec["mul"]]
            if any(v is None for v in vals):
                return None
            result = 1.0
            for v in vals:
                result *= v
            return result
        if "div" in spec:
            pair = spec["div"]
            if len(pair) != 2:
                return None
            a = resolve_length(pair[0], ctx)
            b = resolve_length(pair[1], ctx)
            if a is None or b is None or abs(b) < 1e-15:
                return None
            return a / b
        if "proj" in spec:
            # Signed projection of (target - anchor) onto direction
            proj = spec["proj"]
            target = resolve_point(proj.get("target"), ctx)
            anchor = resolve_point(proj.get("anchor"), ctx)
            direction = resolve_dir(proj.get("dir"), ctx)
            if target and anchor and direction:
                dx = target[0] - anchor[0]
                dy = target[1] - anchor[1]
                return dx * direction[0] + dy * direction[1]
            return None
        if "dist" in spec:
            pair = spec["dist"]
            if len(pair) != 2:
                return None
            a = resolve_point(pair[0], ctx)
            b = resolve_point(pair[1], ctx)
            if a is None or b is None:
                return None
            return math.sqrt((b[0] - a[0])**2 + (b[1] - a[1])**2)
    return None
