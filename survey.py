"""Shared geometry computation for the hut2 survey project.

Contains type definitions, geometry utilities, path operations,
traverse computation, three-arc system, inner wall computation,
and interior layout helpers.
"""
import math
from typing import NamedTuple

# ============================================================
# Types
# ============================================================
Point = tuple[float, float]

class LineSeg(NamedTuple):
    start: str; end: str

class ArcSeg(NamedTuple):
    start: str; end: str; center: str
    radius: float; direction: str  # "CW" or "CCW"
    n_pts: int

Segment = LineSeg | ArcSeg

# ============================================================
# Geometry Utilities
# ============================================================
def left_norm(p1: Point, p2: Point) -> Point:
    dx = p2[0]-p1[0]; dy = p2[1]-p1[1]; Ln = math.sqrt(dx**2+dy**2)
    return (-dy/Ln, dx/Ln)

def off_pt(p: Point, n: Point, d: float) -> Point:
    return (p[0]+d*n[0], p[1]+d*n[1])

def line_isect(p1: Point, d1: Point, p2: Point, d2: Point) -> Point | None:
    det = d1[0]*d2[1]-d1[1]*d2[0]
    if abs(det) < 1e-12: return None
    t = ((p2[0]-p1[0])*d2[1]-(p2[1]-p1[1])*d2[0])/det
    return (p1[0]+t*d1[0], p1[1]+t*d1[1])

def arc_poly(cx, cy, r, sa, ea, n=60):
    return [(cx+r*math.cos(sa+(ea-sa)*i/n), cy+r*math.sin(sa+(ea-sa)*i/n))
            for i in range(n+1)]

def circle_circle_isect(c1: Point, r1: float, c2: Point, r2: float, near: Point) -> Point:
    dx = c2[0]-c1[0]; dy = c2[1]-c1[1]; d = math.sqrt(dx**2+dy**2)
    a = (r1**2-r2**2+d**2)/(2*d); h = math.sqrt(r1**2-a**2)
    ux, uy = dx/d, dy/d; Mx, My = c1[0]+a*ux, c1[1]+a*uy
    I1 = (Mx+h*(-uy), My+h*ux); I2 = (Mx-h*(-uy), My-h*ux)
    d1 = (I1[0]-near[0])**2+(I1[1]-near[1])**2
    d2 = (I2[0]-near[0])**2+(I2[1]-near[1])**2
    return I1 if d1 < d2 else I2

def line_circle_isect_min_t_gt(p: Point, d: Point, c: Point, r: float, t_min: float) -> Point:
    ax = p[0]-c[0]; ay = p[1]-c[1]
    A = d[0]**2+d[1]**2; B = 2*(ax*d[0]+ay*d[1]); C = ax**2+ay**2-r**2
    disc = B**2-4*A*C
    t1 = (-B+math.sqrt(disc))/(2*A); t2 = (-B-math.sqrt(disc))/(2*A)
    t = min(t for t in [t1, t2] if t > t_min)
    return (p[0]+t*d[0], p[1]+t*d[1])

def line_circle_isect_min_abs_t(p: Point, d: Point, c: Point, r: float) -> Point:
    ax = p[0]-c[0]; ay = p[1]-c[1]
    A = d[0]**2+d[1]**2; B = 2*(ax*d[0]+ay*d[1]); C = ax**2+ay**2-r**2
    disc = B**2-4*A*C
    t1 = (-B+math.sqrt(disc))/(2*A); t2 = (-B-math.sqrt(disc))/(2*A)
    t = min(t1, t2, key=lambda t: abs(t))
    return (p[0]+t*d[0], p[1]+t*d[1])

def poly_area(verts):
    n = len(verts); a = 0
    for i in range(n):
        j = (i+1)%n; a += verts[i][0]*verts[j][1]-verts[j][0]*verts[i][1]
    return abs(a)/2

# ============================================================
# Path Operations
# ============================================================
def segment_polyline(seg: Segment, pts: dict) -> list[Point]:
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

def path_polygon(segments: list[Segment], pts: dict) -> list[Point]:
    polygon = []
    for i, seg in enumerate(segments):
        poly = segment_polyline(seg, pts)
        if i > 0:
            poly = poly[1:]
        polygon.extend(poly)
    polygon.pop()  # remove closing point (= polygon[0])
    return polygon

def arc_sweep_deg(seg: ArcSeg, pts: dict) -> float:
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
def brg_dist(p1, p2):
    dE = p2[0]-p1[0]; dN = p2[1]-p1[1]
    d = math.sqrt(dE**2+dN**2)
    b = math.degrees(math.atan2(dE, dN)) % 360
    return b, d

def fmt_brg(b):
    d = int(b); m = int((b-d)*60); sc = (b-d-m/60)*3600
    return f"{d:d}\u00b0 {m:02d}' {sc:04.1f}\""

def fmt_dist(ft):
    f = int(ft); i = (ft-f)*12
    return f"{f}' {i:.1f}\""

# ============================================================
# Traverse Computation
# ============================================================
def compute_traverse():
    """Compute traverse from raw survey legs.

    Returns (pts, p3_trav) where pts is a dict with keys
    P3/POB/P2/P4/P5 (P3-based coordinates) and p3_trav is
    the raw P3 position needed for SVG transform calibration.
    """
    legs = [(257,53,45,19,1.0),(180,54,31,26,11.0),(93,36,7,31,10.5),
            (56,36,31,13,2.5),(317,11,44,34,11.5)]
    _trav = [(0.0, 0.0)]
    for deg, mn, sec, ft, inch in legs:
        brg = deg + mn/60.0 + sec/3600.0
        dist_in = ft * 12 + inch
        brg_rad = math.radians(brg)
        dE = dist_in * math.sin(brg_rad); dN = dist_in * math.cos(brg_rad)
        last = _trav[-1]; _trav.append((last[0]+dE, last[1]+dN))
    _trav_ft = [(e/12, n/12) for e, n in _trav[:5]]
    _trav_ft[2] = (-19.1177, _trav_ft[3][1])
    _trav_ft[1] = (_trav_ft[2][0], _trav_ft[2][1] + 29.0)

    _p3_trav = _trav_ft[2]
    pts = {}
    pts["P3"]  = (0.0, 0.0)
    pts["POB"] = (_trav_ft[0][0] - _p3_trav[0], _trav_ft[0][1] - _p3_trav[1])
    pts["P2"]  = (_trav_ft[1][0] - _p3_trav[0], _trav_ft[1][1] - _p3_trav[1])
    pts["P4"]  = (_trav_ft[3][0] - _p3_trav[0], _trav_ft[3][1] - _p3_trav[1])
    pts["P5"]  = (_trav_ft[4][0] - _p3_trav[0], _trav_ft[4][1] - _p3_trav[1])
    return pts, _p3_trav

# ============================================================
# Three-Arc System
# ============================================================
def compute_three_arc(pts):
    """Compute three-arc boundary system.

    Mutates pts adding T1/C1/T2/C2/PA/T3/C3/PX.
    Returns dict with R1, R2, R3, uE, uN, nE, nN.
    """
    dE_l = pts["P5"][0]-pts["POB"][0]; dN_l = pts["P5"][1]-pts["POB"][1]
    L = math.sqrt(dE_l**2+dN_l**2)
    uE, uN = dE_l/L, dN_l/L
    nE, nN = -uN, uE

    R1, R2 = 10.0, 12.5
    T1_dist, T2_dist = 26.5, 5.75
    pts["T1"] = (pts["POB"][0]+T1_dist*uE, pts["POB"][1]+T1_dist*uN)
    pts["C1"] = (pts["T1"][0]+R1*nE, pts["T1"][1]+R1*nN)
    pts["T2"] = (pts["POB"][0]+T2_dist*uE, pts["POB"][1]+T2_dist*uN)
    pts["C2"] = (pts["T2"][0]+R2*nE, pts["T2"][1]+R2*nN)

    # PA: circle-circle intersection
    dx_cc = pts["C2"][0]-pts["C1"][0]; dy_cc = pts["C2"][1]-pts["C1"][1]
    d_cc = math.sqrt(dx_cc**2+dy_cc**2)
    a_cc = (R1**2-R2**2+d_cc**2)/(2*d_cc); h_cc = math.sqrt(R1**2-a_cc**2)
    ux_cc, uy_cc = dx_cc/d_cc, dy_cc/d_cc
    Mx, My = pts["C1"][0]+a_cc*ux_cc, pts["C1"][1]+a_cc*uy_cc
    I1 = (Mx+h_cc*(-uy_cc), My+h_cc*ux_cc); I2 = (Mx-h_cc*(-uy_cc), My-h_cc*ux_cc)
    ang_T2_C2 = math.atan2(pts["T2"][1]-pts["C2"][1], pts["T2"][0]-pts["C2"][0])
    def ccw_a(s, e): return (e-s)%(2*math.pi)
    s1 = ccw_a(ang_T2_C2, math.atan2(I1[1]-pts["C2"][1], I1[0]-pts["C2"][0]))
    s2 = ccw_a(ang_T2_C2, math.atan2(I2[1]-pts["C2"][1], I2[0]-pts["C2"][0]))
    pts["PA"] = I1 if s1 < s2 else I2

    R3 = 11.0
    T3_dist_from_P3 = 17.911244
    pts["T3"] = (pts["P3"][0]+T3_dist_from_P3, pts["P3"][1])
    pts["C3"] = (pts["T3"][0], pts["T3"][1]-R3)

    dxL = pts["P4"][0]-pts["P5"][0]; dyL = pts["P4"][1]-pts["P5"][1]
    pts["PX"] = line_circle_isect_min_t_gt(pts["P5"], (dxL, dyL), pts["C3"], R3, 1.0)

    return {"R1": R1, "R2": R2, "R3": R3, "uE": uE, "uN": uN, "nE": nE, "nN": nN}

# ============================================================
# Inner Wall Computation
# ============================================================
def compute_inner_walls(outline_segs, pts, wall_t, radii):
    """Compute inner wall points and segments from outline.

    Mutates pts adding W-series points.
    Returns inner_segs list.

    radii keys: R_fillet, R_w1, R_w2, R_wall, R_f_po5, R1i,
                R_turn3, R_turn2, R_turn1, R_fillet2, R_t4
    """
    def _inner_point(seg_b, seg_a):
        if not isinstance(seg_b, LineSeg) and not isinstance(seg_a, LineSeg):
            c1 = pts[seg_b.center]; c2 = pts[seg_a.center]
            r1 = (seg_b.radius + wall_t) if seg_b.direction == "CW" else (seg_b.radius - wall_t)
            dx = c2[0]-c1[0]; dy = c2[1]-c1[1]; d = math.sqrt(dx*dx+dy*dy)
            return (c1[0]+r1*dx/d, c1[1]+r1*dy/d)
        ls = seg_b if isinstance(seg_b, LineSeg) else seg_a
        arc = seg_a if isinstance(seg_b, LineSeg) else seg_b
        c = pts[arc.center]; S = pts[ls.start]; E = pts[ls.end]
        D = (E[0]-S[0], E[1]-S[1]); LN = left_norm(S, E); P = off_pt(S, LN, wall_t)
        t = ((c[0]-P[0])*D[0]+(c[1]-P[1])*D[1])/(D[0]**2+D[1]**2)
        return (P[0]+t*D[0], P[1]+t*D[1])

    n_segs = len(outline_segs)
    for i in range(n_segs):
        w_name = "W" + outline_segs[i].end[1:]
        pts[w_name] = _inner_point(outline_segs[i], outline_segs[(i+1)%n_segs])

    R = radii
    inner_segs = [
        LineSeg("W2","W1"), ArcSeg("W1","W0","Cf",R["R_fillet"]-wall_t,"CCW",20),
        LineSeg("W0","W15"), ArcSeg("W15","W14","Cw1",R["R_w1"]+wall_t,"CW",60),
        ArcSeg("W14","W13b","Cw2",R["R_w2"]-wall_t,"CCW",60), LineSeg("W13b","W13a"),
        ArcSeg("W13a","W13","Cw3",R["R_wall"]-wall_t,"CCW",20), LineSeg("W13","W12"),
        ArcSeg("W12","W11","Cf4",R["R_f_po5"]-wall_t,"CCW",20), LineSeg("W11","W10a"),
        ArcSeg("W10a","W10","Ct4",R["R_t4"]-wall_t,"CCW",20),
        ArcSeg("W10","W9","C1",R["R1i"]+wall_t,"CW",60), LineSeg("W9","W8"),
        ArcSeg("W8","W7","Ct3",R["R_turn3"]-wall_t,"CCW",20), LineSeg("W7","W6"),
        ArcSeg("W6","W5","Ct2",R["R_turn2"]+wall_t,"CW",20),
        ArcSeg("W5","W4","Ct1",R["R_turn1"]-wall_t,"CCW",20),
        LineSeg("W4","W3"), ArcSeg("W3","W2","Cf2",R["R_fillet2"]-wall_t,"CCW",20),
    ]
    return inner_segs

# ============================================================
# Interior Layout Helpers
# ============================================================
def horiz_isects(poly, n_val):
    """Easting values where polygon boundary crosses a given northing."""
    r = []
    for i in range(len(poly)):
        j = (i+1)%len(poly); n1, n2 = poly[i][1], poly[j][1]
        if (n1 <= n_val < n2) or (n2 <= n_val < n1):
            t = (n_val-n1)/(n2-n1); r.append(poly[i][0]+t*(poly[j][0]-poly[i][0]))
    return r

def vert_isects(poly, e_val):
    """Northing values where polygon boundary crosses a given easting."""
    r = []
    for i in range(len(poly)):
        j = (i+1)%len(poly); e1, e2 = poly[i][0], poly[j][0]
        if (e1 <= e_val < e2) or (e2 <= e_val < e1):
            t = (e_val-e1)/(e2-e1); r.append(poly[i][1]+t*(poly[j][1]-poly[i][1]))
    return r

def compute_interior_layout(pts, inner_poly):
    """Compute interior layout positions.

    Returns dict with iw1, iw2_w/e/s/n, dryer, washer, counter,
    wall8, iw3, iw4, wall5, bed, and intermediate values.
    """
    iwt = 6.0/12.0
    iw1_s = pts["W0"][1]+11.5; iw1_n = iw1_s+iwt
    si = horiz_isects(inner_poly, iw1_s)
    ni = horiz_isects(inner_poly, iw1_n)
    iw1 = [(min(si),iw1_s),(max(si),iw1_s),(max(ni),iw1_n),(min(ni),iw1_n)]

    iw2_w = pts["W1"][0]+6.5; iw2_e = iw2_w+iwt; iw2_s = iw1_n; iw2_n = pts["W3"][1]

    dryer_w = pts["W1"][0]+0.5; dryer_s = pts["W0"][1]+4.0/12
    dryer_e = dryer_w+35.0/12; dryer_n = dryer_s+30.0/12
    washer_w = dryer_w; washer_s = dryer_n+1.0/12
    washer_e = dryer_e; washer_n = washer_s+30.0/12

    ctr_w = dryer_e+3.0; ctr_e = ctr_w+2.0; ctr_s = pts["W0"][1]; ctr_n = ctr_s+6.0
    ctr_nw_r = 9.0/12.0

    iwt3 = 3.0/12; iwt4 = 4.0/12; cl2w = 30.0/12
    w8 = [(ctr_e,ctr_s),(ctr_e+iwt3,ctr_s),(ctr_e+iwt3,ctr_n),
          (ctr_e+iwt3+cl2w,ctr_n),(ctr_e+iwt3+cl2w,ctr_n+iwt3),(ctr_e,ctr_n+iwt3)]
    iw3_w = ctr_e+iwt3+cl2w; iw3_e = iw3_w+iwt4; iw3_s = ctr_s; iw3_n = iw1_s
    iw4_w = iw3_e+140.0/12; iw4_e = iw4_w+iwt4
    wall_south_n = 2.0/12
    cl1w = 30.0/12; cl1_top = ctr_n-1.0
    w5_w = iw4_e+cl1w; w5_e = w5_w+iwt3
    w5 = [(iw4_e,cl1_top+iwt3),(w5_e,cl1_top+iwt3),(w5_e,wall_south_n),
          (w5_w,wall_south_n),(w5_w,cl1_top),(iw4_e,cl1_top)]

    bed_cx = (iw3_e+iw4_w)/2
    bed_w = bed_cx-76.0/24; bed_e = bed_cx+76.0/24
    bed_s = ctr_s+2.0/12; bed_n = bed_s+94.0/12

    return {
        "iw1": iw1, "iw1_s": iw1_s, "iw1_n": iw1_n, "iwt": iwt,
        "iw2_w": iw2_w, "iw2_e": iw2_e, "iw2_s": iw2_s, "iw2_n": iw2_n,
        "dryer_w": dryer_w, "dryer_s": dryer_s, "dryer_e": dryer_e, "dryer_n": dryer_n,
        "washer_w": washer_w, "washer_s": washer_s, "washer_e": washer_e, "washer_n": washer_n,
        "ctr_w": ctr_w, "ctr_e": ctr_e, "ctr_s": ctr_s, "ctr_n": ctr_n, "ctr_nw_r": ctr_nw_r,
        "iwt3": iwt3, "iwt4": iwt4,
        "wall8": w8,
        "iw3_w": iw3_w, "iw3_e": iw3_e, "iw3_s": iw3_s, "iw3_n": iw3_n,
        "iw4_w": iw4_w, "iw4_e": iw4_e, "wall_south_n": wall_south_n,
        "wall5": w5, "w5_w": w5_w, "w5_e": w5_e,
        "cl1_top": cl1_top,
        "bed_w": bed_w, "bed_e": bed_e, "bed_s": bed_s, "bed_n": bed_n, "bed_cx": bed_cx,
    }
