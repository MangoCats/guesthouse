"""Pure geometry functions, path operations, polygon utilities, and formatting."""
import math
from .types import Point, LineSeg, ArcSeg, Segment

# ============================================================
# Error Type
# ============================================================
class GeometryError(ValueError):
    """Raised for impossible geometry operations."""

# ============================================================
# Geometry Utilities
# ============================================================
def left_norm(p1: Point, p2: Point) -> Point:
    """Unit normal vector to the left of the direction p1 → p2 (CCW perpendicular)."""
    dx = p2[0]-p1[0]; dy = p2[1]-p1[1]; Ln = math.sqrt(dx**2+dy**2)
    return (-dy/Ln, dx/Ln)

def off_pt(p: Point, n: Point, d: float) -> Point:
    """Offset point p by distance d along unit direction n."""
    return (p[0]+d*n[0], p[1]+d*n[1])

def line_isect(p1: Point, d1: Point, p2: Point, d2: Point) -> Point:
    """Intersection of two lines (p1+t*d1) and (p2+s*d2). Raises GeometryError if parallel."""
    det = d1[0]*d2[1]-d1[1]*d2[0]
    if abs(det) < 1e-12:
        raise GeometryError(f"Parallel lines: det={det:.2e}")
    t = ((p2[0]-p1[0])*d2[1]-(p2[1]-p1[1])*d2[0])/det
    return (p1[0]+t*d1[0], p1[1]+t*d1[1])

def arc_poly(cx: float, cy: float, r: float, sa: float, ea: float, n: int = 60) -> list[Point]:
    """Generate n+1 points along a circular arc from angle sa to ea (radians)."""
    return [(cx+r*math.cos(sa+(ea-sa)*i/n), cy+r*math.sin(sa+(ea-sa)*i/n))
            for i in range(n+1)]

def circle_circle_isect(c1: Point, r1: float, c2: Point, r2: float, near: Point) -> Point:
    """Intersection of two circles, returning the point nearest to *near*.

    Raises GeometryError if the circles don't intersect.
    """
    dx = c2[0]-c1[0]; dy = c2[1]-c1[1]; d = math.sqrt(dx**2+dy**2)
    if d > r1 + r2 + 1e-9:
        raise GeometryError(f"Circles too far apart: d={d:.6f}, r1+r2={r1+r2:.6f}")
    if d < abs(r1-r2) - 1e-9:
        raise GeometryError(f"Circle contained: d={d:.6f}, |r1-r2|={abs(r1-r2):.6f}")
    a = (r1**2-r2**2+d**2)/(2*d)
    h_sq = r1**2-a**2
    if h_sq < -1e-9:
        raise GeometryError(f"No intersection: h^2={h_sq:.6f}")
    h = math.sqrt(max(0, h_sq))
    ux, uy = dx/d, dy/d; Mx, My = c1[0]+a*ux, c1[1]+a*uy
    I1 = (Mx+h*(-uy), My+h*ux); I2 = (Mx-h*(-uy), My-h*ux)
    d1 = (I1[0]-near[0])**2+(I1[1]-near[1])**2
    d2 = (I2[0]-near[0])**2+(I2[1]-near[1])**2
    return I1 if d1 < d2 else I2

def line_circle_isect_min_t_gt(p: Point, d: Point, c: Point, r: float, t_min: float) -> Point:
    """Line-circle intersection with smallest parameter t > t_min.

    Line parametrised as p + t*d; circle centered at c with radius r.
    """
    ax = p[0]-c[0]; ay = p[1]-c[1]
    A = d[0]**2+d[1]**2; B = 2*(ax*d[0]+ay*d[1]); C = ax**2+ay**2-r**2
    disc = B**2-4*A*C
    if disc < -1e-9:
        raise GeometryError(f"Line misses circle: disc={disc:.6f}")
    disc = max(0, disc)
    t1 = (-B+math.sqrt(disc))/(2*A); t2 = (-B-math.sqrt(disc))/(2*A)
    candidates = [t for t in [t1, t2] if t > t_min]
    if not candidates:
        raise GeometryError(f"No intersection with t > {t_min}: t1={t1}, t2={t2}")
    t = min(candidates)
    return (p[0]+t*d[0], p[1]+t*d[1])

def line_circle_isect_min_abs_t(p: Point, d: Point, c: Point, r: float) -> Point:
    """Line-circle intersection with smallest |t|.

    Line parametrised as p + t*d; circle centered at c with radius r.
    """
    ax = p[0]-c[0]; ay = p[1]-c[1]
    A = d[0]**2+d[1]**2; B = 2*(ax*d[0]+ay*d[1]); C = ax**2+ay**2-r**2
    disc = B**2-4*A*C
    if disc < -1e-9:
        raise GeometryError(f"Line misses circle: disc={disc:.6f}")
    disc = max(0, disc)
    t1 = (-B+math.sqrt(disc))/(2*A); t2 = (-B-math.sqrt(disc))/(2*A)
    t = min(t1, t2, key=lambda t: abs(t))
    return (p[0]+t*d[0], p[1]+t*d[1])

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
    d = math.sqrt(dE**2+dN**2)
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

    At the F8-F9 concave corner, the inner shell goes straight south,
    makes a tight CCW turn (radius *R_turn*), then goes straight east
    along the F9-F10 bearing.

    Returns list of (E, N) points from the inset-F8 to inset-F9 position.
    """
    F8 = pts["F8"]
    C8 = pts["C8"]
    R_a8 = C8[0] - F8[0]  # SMALL_ARC_R

    start_E = F8[0] - inset
    start_N = F8[1]
    end_E   = C8[0]                    # F9 easting
    end_N   = F8[1] - R_a8 - inset     # F9 northing − inset

    d = R_a8 + inset - R_turn          # straight length (equal both sides)
    arc_cx = start_E + R_turn
    arc_cy = start_N - d

    polyline: list[Point] = [(start_E, start_N)]
    for i in range(n_arc + 1):
        angle = math.pi + i * (math.pi / 2) / n_arc
        polyline.append((arc_cx + R_turn * math.cos(angle),
                         arc_cy + R_turn * math.sin(angle)))
    polyline.append((end_E, end_N))
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

    radii keys: R_a0, R_a20, R_a19, R_a17, R_a15,
                R_a11, R_a8, R_a7, R_a5, R_a13, R_a10
    """
    def _inner_point(seg_b, seg_a):
        _wt = -wall_t  # negated: left_norm points exterior for CW traversal; negate to offset inward
        if not isinstance(seg_b, LineSeg) and not isinstance(seg_a, LineSeg):
            c1 = pts[seg_b.center]; c2 = pts[seg_a.center]
            r1 = (seg_b.radius + _wt) if seg_b.direction == "CW" else (seg_b.radius - _wt)
            dx = c2[0]-c1[0]; dy = c2[1]-c1[1]; d = math.sqrt(dx*dx+dy*dy)
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

    R = radii
    inner_segs = [
        ArcSeg("W0","W1","C0",R["R_a0"]-wall_t,"CW",20),
        LineSeg("W1","W2"),
        ArcSeg("W2","W3","C2",R["R_a2"]-wall_t,"CW",20),
        ArcSeg("W3","W4","C3",R["R_a3"]+wall_t,"CCW",20),
        LineSeg("W4","W5"),
        ArcSeg("W5","W6","C5",R["R_a5"]-wall_t,"CW",20),
        LineSeg("W6","W7"),
        ArcSeg("W7","W8","C7",R["R_a7"]-wall_t,"CW",20),
        ArcSeg("W8","W9","C8",R["R_a8"]+wall_t,"CCW",20),
        LineSeg("W9","W10"),
        ArcSeg("W10","W11","C10",R["R_a10"]+wall_t,"CCW",20),
        ArcSeg("W11","W12","C11",R["R_a11"]-wall_t,"CW",60),
        LineSeg("W12","W13"),
        ArcSeg("W13","W14","C13",R["R_a13"]-wall_t,"CW",60),
        LineSeg("W14","W15"),
        ArcSeg("W15","W16","C15",R["R_a15"]-wall_t,"CW",20),
        LineSeg("W16","W17"),
        ArcSeg("W17","W18","C17",R["R_a17"]-wall_t,"CW",20),
        LineSeg("W18","W19"),
        ArcSeg("W19","W20","C19",R["R_a19"]-wall_t,"CW",60),
        ArcSeg("W20","W21","C20",R["R_a20"]+wall_t,"CCW",60),
        LineSeg("W21","W0"),
    ]
    return inner_segs
