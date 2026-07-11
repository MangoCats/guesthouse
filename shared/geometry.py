"""Pure geometry functions, path operations, polygon utilities, and formatting."""
import math
from .types import Point, BBox, LineSeg, ArcSeg, Segment

# Geometric tolerance for point/intersection coincidence checks
GEOM_EPS = 1e-9

# Default point count for arc-to-polygon conversions
ARC_POLY_PTS = 60

# ============================================================
# Error Type
# ============================================================
class GeometryError(ValueError):
    """Raised for impossible geometry operations."""


def require_pts(pts: dict[str, 'Point'], *names: str) -> None:
    """Raise KeyError with a descriptive message if any named point is missing."""
    missing = [n for n in names if n not in pts]
    if missing:
        raise KeyError(f"Missing geometry point(s): {', '.join(missing)}")


# ============================================================
# Geometry Utilities
# ============================================================
def left_norm(p1: Point, p2: Point) -> Point:
    """Unit normal vector to the left of the direction p1 → p2 (CCW perpendicular)."""
    dx = p2[0]-p1[0]; dy = p2[1]-p1[1]; Ln = math.hypot(dx, dy)
    if Ln < GEOM_EPS:
        raise GeometryError(f"Zero-length segment in left_norm: p1={p1}, p2={p2}")
    return (-dy/Ln, dx/Ln)

def off_pt(p: Point, n: Point, d: float) -> Point:
    """Offset point p by distance d along unit direction n.

    Equivalent to ``offset_pt(p, d, n)`` with different arg order.
    """
    return (p[0]+d*n[0], p[1]+d*n[1])

def seg_vec(p1: Point, p2: Point) -> tuple[float, float, float]:
    """Return (dE, dN, length) from p1 to p2."""
    dE = p2[0]-p1[0]; dN = p2[1]-p1[1]
    return dE, dN, math.hypot(dE, dN)

def seg_vecs(p1: Point, p2: Point) -> tuple[tuple[float, float], tuple[float, float]]:
    """Along-direction and CW-inward normal unit vectors for segment p1→p2.

    Returns (along, inward) where:
    - along: unit vector from p1 toward p2
    - inward: right perpendicular of along (CW-inward for CW outline traversal)
    """
    dx, dy = p2[0] - p1[0], p2[1] - p1[1]
    length = math.hypot(dx, dy)
    if length < GEOM_EPS:
        raise GeometryError(f"Zero-length segment in seg_vecs: p1={p1}, p2={p2}")
    along = (dx / length, dy / length)
    inward = (dy / length, -dx / length)
    return along, inward

def offset_pt(origin: Point, dist: float, direction: tuple[float, float]) -> Point:
    """Offset point by dist along direction vector.  Thin wrapper around off_pt."""
    return off_pt(origin, direction, dist)

def line_isect(p1: Point, d1: Point, p2: Point, d2: Point) -> Point:
    """Intersection of two lines (p1+t*d1) and (p2+s*d2). Raises GeometryError if parallel."""
    det = d1[0]*d2[1]-d1[1]*d2[0]
    if abs(det) < 1e-12:
        raise GeometryError(f"Parallel lines: det={det:.2e}")
    t = ((p2[0]-p1[0])*d2[1]-(p2[1]-p1[1])*d2[0])/det
    return (p1[0]+t*d1[0], p1[1]+t*d1[1])

def arc_poly(cx: float, cy: float, r: float, sa: float, ea: float, n: int = ARC_POLY_PTS) -> list[Point]:
    """Generate n+1 points along a circular arc from angle sa to ea (radians)."""
    return [(cx+r*math.cos(sa+(ea-sa)*i/n), cy+r*math.sin(sa+(ea-sa)*i/n))
            for i in range(n+1)]

def circle_circle_isect(c1: Point, r1: float, c2: Point, r2: float, near: Point) -> Point:
    """Intersection of two circles, returning the point nearest to *near*.

    Raises GeometryError if the circles don't intersect.
    """
    dx = c2[0]-c1[0]; dy = c2[1]-c1[1]; d = math.hypot(dx, dy)
    if d > r1 + r2 + GEOM_EPS:
        raise GeometryError(f"Circles too far apart: d={d:.6f}, r1+r2={r1+r2:.6f}")
    if d < abs(r1-r2) - GEOM_EPS:
        raise GeometryError(f"Circle contained: d={d:.6f}, |r1-r2|={abs(r1-r2):.6f}")
    a = (r1**2-r2**2+d**2)/(2*d)
    h_sq = r1**2-a**2
    if h_sq < -GEOM_EPS:
        raise GeometryError(f"No intersection: h^2={h_sq:.6f}")
    h = math.sqrt(max(0, h_sq))
    ux, uy = dx/d, dy/d; Mx, My = c1[0]+a*ux, c1[1]+a*uy
    I1 = (Mx+h*(-uy), My+h*ux); I2 = (Mx-h*(-uy), My-h*ux)
    d1 = (I1[0]-near[0])**2+(I1[1]-near[1])**2
    d2 = (I2[0]-near[0])**2+(I2[1]-near[1])**2
    return I1 if d1 < d2 else I2

def _line_circle_t_values(p: Point, d: Point, c: Point, r: float) -> tuple[float, float]:
    """Solve line-circle intersection for parameter values t1, t2.

    Line parametrised as p + t*d; circle centered at c with radius r.
    Returns (t1, t2) where t1 >= t2.
    """
    ax = p[0]-c[0]; ay = p[1]-c[1]
    A = d[0]**2+d[1]**2; B = 2*(ax*d[0]+ay*d[1]); C = ax**2+ay**2-r**2
    disc = B**2-4*A*C
    if disc < -GEOM_EPS:
        raise GeometryError(f"Line misses circle: disc={disc:.6f}")
    disc = max(0, disc)
    t1 = (-B+math.sqrt(disc))/(2*A); t2 = (-B-math.sqrt(disc))/(2*A)
    return t1, t2

def line_circle_isect_min_t_gt(p: Point, d: Point, c: Point, r: float, t_min: float) -> Point:
    """Line-circle intersection with smallest parameter t > t_min.

    Line parametrised as p + t*d; circle centered at c with radius r.
    """
    t1, t2 = _line_circle_t_values(p, d, c, r)
    candidates = [t for t in [t1, t2] if t > t_min]
    if not candidates:
        raise GeometryError(f"No intersection with t > {t_min}: t1={t1}, t2={t2}")
    t = min(candidates)
    return (p[0]+t*d[0], p[1]+t*d[1])

def line_circle_isect_min_abs_t(p: Point, d: Point, c: Point, r: float) -> Point:
    """Line-circle intersection with smallest |t|.

    Line parametrised as p + t*d; circle centered at c with radius r.
    """
    t1, t2 = _line_circle_t_values(p, d, c, r)
    t = min(t1, t2, key=lambda t: abs(t))
    return (p[0]+t*d[0], p[1]+t*d[1])

def bbox_from_points(pts: list[Point]) -> BBox:
    """Compute axis-aligned bounding box from a list of points."""
    return BBox(w=min(p[0] for p in pts), s=min(p[1] for p in pts),
                e=max(p[0] for p in pts), n=max(p[1] for p in pts))


def poly_centroid(poly: list[Point]) -> Point:
    """Area-weighted centroid of a simple polygon.

    Falls back to vertex average for degenerate (zero-area) polygons.
    """
    n = len(poly)
    area6 = 0.0
    cx = cy = 0.0
    for i in range(n):
        x0, y0 = poly[i]
        x1, y1 = poly[(i + 1) % n]
        cross = x0 * y1 - x1 * y0
        area6 += cross
        cx += (x0 + x1) * cross
        cy += (y0 + y1) * cross
    if abs(area6) < 1e-12:
        return (sum(p[0] for p in poly) / n, sum(p[1] for p in poly) / n)
    cx /= 3.0 * area6
    cy /= 3.0 * area6
    return (cx, cy)


def poly_area(verts: list[Point]) -> float:
    """Polygon area via the shoelace formula. Works for either winding order."""
    n = len(verts); a = 0
    for i in range(n):
        j = (i+1)%n; a += verts[i][0]*verts[j][1]-verts[j][0]*verts[i][1]
    return abs(a)/2

# ============================================================
# Path Operations
# ============================================================
def segment_polyline(seg: Segment, pts: dict[str, Point]) -> list[Point]:
    """Convert a LineSeg or ArcSeg to a polyline of coordinate points."""
    if isinstance(seg, LineSeg):
        return [pts[seg.start], pts[seg.end]]
    c = pts[seg.center]
    ang_s = math.atan2(pts[seg.start][1]-c[1], pts[seg.start][0]-c[0])
    ang_e = math.atan2(pts[seg.end][1]-c[1], pts[seg.end][0]-c[0])
    if seg.direction == "CW":
        sweep = (ang_s - ang_e) % (2*math.pi)
        return arc_poly(c[0], c[1], seg.radius, ang_s, ang_s - sweep, seg.n_pts)
    else:
        sweep = (ang_e - ang_s) % (2*math.pi)
        return arc_poly(c[0], c[1], seg.radius, ang_s, ang_s + sweep, seg.n_pts)

def path_polygon(segments: list[Segment], pts: dict[str, Point]) -> list[Point]:
    """Convert a closed segment path into a polygon vertex list."""
    polygon = []
    for i, seg in enumerate(segments):
        poly = segment_polyline(seg, pts)
        if i > 0:
            poly = poly[1:]
        polygon.extend(poly)
    polygon.pop()  # remove closing point (= polygon[0])
    return polygon

def arc_sweep_deg(seg: ArcSeg, pts: dict[str, Point]) -> float:
    """Sweep angle of an arc segment in degrees (always positive)."""
    c = pts[seg.center]
    ang_s = math.atan2(pts[seg.start][1]-c[1], pts[seg.start][0]-c[0])
    ang_e = math.atan2(pts[seg.end][1]-c[1], pts[seg.end][0]-c[0])
    if seg.direction == "CW":
        return math.degrees((ang_s - ang_e) % (2*math.pi))
    else:
        return math.degrees((ang_e - ang_s) % (2*math.pi))

# ============================================================
# Formatting Helpers
# ============================================================
def brg_dist(p1: Point, p2: Point) -> tuple[float, float]:
    """Bearing (degrees clockwise from North) and distance between two E/N points."""
    dE = p2[0]-p1[0]; dN = p2[1]-p1[1]
    d = math.hypot(dE, dN)
    b = math.degrees(math.atan2(dE, dN)) % 360
    return b, d

def fmt_brg(b: float) -> str:
    """Format bearing in degrees to DMS string, e.g. '257\u00b0 53' 45.0\"'."""
    d = int(b); m = int((b-d)*60); sc = (b-d-m/60)*3600
    return f"{d:d}\u00b0 {m:02d}' {sc:04.1f}\""

def fmt_dist(ft: float) -> str:
    """Format distance in feet to feet-inches string, e.g. \"2' 6\\\"\"."""
    total_in = round(ft * 12, 2)
    whole_ft = int(total_in // 12)
    remaining_in = total_in - whole_ft * 12
    if remaining_in >= 12.0:
        whole_ft += 1
        remaining_in = 0.0
    in_str = f"{remaining_in:.2f}".rstrip('0').rstrip('.')
    return f"{whole_ft}' {in_str}\""

# ============================================================
# Polygon Utilities
# ============================================================
def horiz_isects(poly: list[Point], n_val: float) -> list[float]:
    """Easting values where polygon boundary crosses a given northing."""
    r = []
    for i in range(len(poly)):
        j = (i+1)%len(poly); n1, n2 = poly[i][1], poly[j][1]
        if (n1 <= n_val < n2) or (n2 <= n_val < n1):
            t = (n_val-n1)/(n2-n1); r.append(poly[i][0]+t*(poly[j][0]-poly[i][0]))
    return r

def vert_isects(poly: list[Point], e_val: float) -> list[float]:
    """Northing values where polygon boundary crosses a given easting."""
    r = []
    for i in range(len(poly)):
        j = (i+1)%len(poly); e1, e2 = poly[i][0], poly[j][0]
        if (e1 <= e_val < e2) or (e2 <= e_val < e1):
            t = (e_val-e1)/(e2-e1); r.append(poly[i][1]+t*(poly[j][1]-poly[i][1]))
    return r

# ============================================================
# F8-F9 Corner Override
# ============================================================
def f8f9_corner_polyline(
    pts: dict[str, Point], inset: float, R_turn: float, n_arc: int = 20,
) -> list[Point]:
    """Straight-arc-straight polyline for the F8-F9 inner shell corner.

    At the F8-F9 concave corner, the inner shell goes straight in the
    F7→F8 exit direction, makes a tight CCW turn (radius *R_turn*),
    then goes straight along the F9→F10 bearing.

    Returns list of (E, N) points from the inset-F8 to inset-F9 position.
    Uses traversal directions from pts (rotation-safe).
    """
    F8, F9, F10, C7 = pts["F8"], pts["F9"], pts["F10"], pts["C7"]

    # CW traversal direction at F8 (tangent at exit of C7 arc)
    _r8x, _r8y = F8[0] - C7[0], F8[1] - C7[1]
    _r8_len = math.hypot(_r8x, _r8y)
    _dir_f8 = (_r8y / _r8_len, -_r8x / _r8_len)  # CW tangent = right normal of radius

    # CW traversal direction at F9 (F9→F10 line direction)
    _d9x, _d9y = F10[0] - F9[0], F10[1] - F9[1]
    _d9_len = math.hypot(_d9x, _d9y)
    _dir_f9 = (_d9x / _d9_len, _d9y / _d9_len)

    # Inset direction (right of CW direction = toward wall material)
    _ins_f8 = (_dir_f8[1], -_dir_f8[0])
    _ins_f9 = (_dir_f9[1], -_dir_f9[0])

    start = (F8[0] + inset * _ins_f8[0], F8[1] + inset * _ins_f8[1])
    end   = (F9[0] + inset * _ins_f9[0], F9[1] + inset * _ins_f9[1])

    # CCW turn arc center: offset R_turn to the LEFT of each direction
    _left_f8 = (-_dir_f8[1], _dir_f8[0])
    _left_f9 = (-_dir_f9[1], _dir_f9[0])

    # Line-line intersection to find arc center
    _P1 = (start[0] + R_turn * _left_f8[0], start[1] + R_turn * _left_f8[1])
    _P2 = (end[0]   + R_turn * _left_f9[0], end[1]   + R_turn * _left_f9[1])
    _dp = (_P2[0] - _P1[0], _P2[1] - _P1[1])
    _cross = _dir_f8[0] * _dir_f9[1] - _dir_f8[1] * _dir_f9[0]
    _t = (_dp[0] * _dir_f9[1] - _dp[1] * _dir_f9[0]) / _cross
    arc_cx = _P1[0] + _t * _dir_f8[0]
    arc_cy = _P1[1] + _t * _dir_f8[1]

    # Arc angles: from tangent-point-1 to tangent-point-2, CCW
    _entry = math.atan2(-_left_f8[1], -_left_f8[0])
    _exit  = math.atan2(-_left_f9[1], -_left_f9[0])
    _sweep = (_exit - _entry) % (2 * math.pi)

    polyline: list[Point] = [start]
    for i in range(n_arc + 1):
        angle = _entry + i * _sweep / n_arc
        polyline.append((arc_cx + R_turn * math.cos(angle),
                         arc_cy + R_turn * math.sin(angle)))
    polyline.append(end)
    return polyline


# ============================================================
# Inner Wall Computation
# ============================================================
def compute_inner_walls(
    outline_segs: list[Segment],
    pts: dict[str, Point],
    wall_t: float,
    radii: dict[str, float],
) -> list[Segment]:
    """Compute inner wall points and segments from outline.

    Mutates pts adding W-series points.
    Returns inner_segs list.

    radii keys: R_a1, R_a17, R_a15,
                R_a11, R_a8, R_a7, R_a5, R_a13, R_a10
    """
    def _inner_point(seg_b, seg_a):
        _wt = -wall_t  # negated: left_norm points exterior for CW traversal; negate to offset inward
        if isinstance(seg_b, LineSeg) and isinstance(seg_a, LineSeg):
            S1, E1 = pts[seg_b.start], pts[seg_b.end]
            S2, E2 = pts[seg_a.start], pts[seg_a.end]
            D1 = (E1[0]-S1[0], E1[1]-S1[1])
            D2 = (E2[0]-S2[0], E2[1]-S2[1])
            LN1 = left_norm(S1, E1); LN2 = left_norm(S2, E2)
            return line_isect(off_pt(S1, LN1, _wt), D1, off_pt(S2, LN2, _wt), D2)
        if not isinstance(seg_b, LineSeg) and not isinstance(seg_a, LineSeg):
            c1 = pts[seg_b.center]; c2 = pts[seg_a.center]
            r1 = (seg_b.radius + _wt) if seg_b.direction == "CW" else (seg_b.radius - _wt)
            dx = c2[0]-c1[0]; dy = c2[1]-c1[1]; d = math.hypot(dx, dy)
            return (c1[0]+r1*dx/d, c1[1]+r1*dy/d)
        ls = seg_b if isinstance(seg_b, LineSeg) else seg_a
        arc = seg_a if isinstance(seg_b, LineSeg) else seg_b
        c = pts[arc.center]; S = pts[ls.start]; E = pts[ls.end]
        D = (E[0]-S[0], E[1]-S[1]); LN = left_norm(S, E); P = off_pt(S, LN, _wt)
        t = ((c[0]-P[0])*D[0]+(c[1]-P[1])*D[1])/(D[0]**2+D[1]**2)
        return (P[0]+t*D[0], P[1]+t*D[1])

    n_segs = len(outline_segs)
    for i in range(n_segs):
        w_name = "W" + outline_segs[i].end[1:]
        pts[w_name] = _inner_point(outline_segs[i], outline_segs[(i+1)%n_segs])

    inner_segs = []
    for seg in outline_segs:
        w_start = "W" + seg.start[1:]
        w_end = "W" + seg.end[1:]
        if isinstance(seg, LineSeg):
            inner_segs.append(LineSeg(w_start, w_end))
        else:
            if seg.direction == "CW":
                r_inner = seg.radius - wall_t
            else:
                r_inner = seg.radius + wall_t
            inner_segs.append(ArcSeg(w_start, w_end, seg.center,
                                     r_inner, seg.direction, seg.n_pts))
    return inner_segs


def offset_outline_outward(
    outline_segs: list[Segment],
    pts: dict[str, Point],
    radii: dict[str, float],
    delta: float,
) -> list[Segment]:
    """Offset the exterior outline outward by ``delta`` feet (in place).

    The outline is traversed CW with the interior on the right, so the exterior
    is to the *left* of each segment.  Moving the whole outline out by ``delta``
    grows the outer-wall thickness on the exterior only: the F-series corner
    points move out, convex (CW) arc radii grow by ``delta`` and concave (CCW)
    arcs shrink by ``delta``, while arc *centers* stay put (arcs stay concentric).

    Combined with ``compute_inner_walls(offset_outline, wall_t)`` this keeps the
    interior face fixed when ``wall_t`` grows by the same ``delta`` — since
    offsetting out by ``delta`` then insetting by ``wall_t`` = insetting the
    original outline by ``wall_t - delta`` (the design-baseline wall).

    Mutates ``pts`` (F-series coordinates) and ``radii`` (``R_a<center>`` keys)
    and returns the new outline segment list.  ``delta == 0`` is a no-op.
    """
    if abs(delta) < GEOM_EPS:
        return outline_segs

    def _outer_point(seg_b, seg_a):
        _wt = delta  # +delta: left of CW traversal is the exterior
        if isinstance(seg_b, LineSeg) and isinstance(seg_a, LineSeg):
            S1, E1 = pts[seg_b.start], pts[seg_b.end]
            S2, E2 = pts[seg_a.start], pts[seg_a.end]
            D1 = (E1[0]-S1[0], E1[1]-S1[1])
            D2 = (E2[0]-S2[0], E2[1]-S2[1])
            LN1 = left_norm(S1, E1); LN2 = left_norm(S2, E2)
            return line_isect(off_pt(S1, LN1, _wt), D1, off_pt(S2, LN2, _wt), D2)
        if not isinstance(seg_b, LineSeg) and not isinstance(seg_a, LineSeg):
            c1 = pts[seg_b.center]; c2 = pts[seg_a.center]
            r1 = (seg_b.radius + _wt) if seg_b.direction == "CW" else (seg_b.radius - _wt)
            dx = c2[0]-c1[0]; dy = c2[1]-c1[1]; d = math.hypot(dx, dy)
            return (c1[0]+r1*dx/d, c1[1]+r1*dy/d)
        ls = seg_b if isinstance(seg_b, LineSeg) else seg_a
        arc = seg_a if isinstance(seg_b, LineSeg) else seg_b
        c = pts[arc.center]; S = pts[ls.start]; E = pts[ls.end]
        D = (E[0]-S[0], E[1]-S[1]); LN = left_norm(S, E); P = off_pt(S, LN, _wt)
        t = ((c[0]-P[0])*D[0]+(c[1]-P[1])*D[1])/(D[0]**2+D[1]**2)
        return (P[0]+t*D[0], P[1]+t*D[1])

    n_segs = len(outline_segs)
    # Compute all new corner points from the *original* coordinates first.
    new_corner = {}
    for i in range(n_segs):
        new_corner[outline_segs[i].end] = _outer_point(
            outline_segs[i], outline_segs[(i+1) % n_segs])
    pts.update(new_corner)

    new_segs = []
    for seg in outline_segs:
        if isinstance(seg, LineSeg):
            new_segs.append(seg)
        else:
            r_new = (seg.radius + delta) if seg.direction == "CW" else (seg.radius - delta)
            radii["R_a" + seg.center[1:]] = r_new
            new_segs.append(ArcSeg(seg.start, seg.end, seg.center,
                                   r_new, seg.direction, seg.n_pts))
    return new_segs
