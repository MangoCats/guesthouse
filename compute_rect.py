import math

# === Traverse (P3 adjusted so P3->P4 bearing = 90 deg exactly) ===
legs = [
    (257, 53, 45, 19, 1.0), (180, 54, 31, 26, 11.0),
    (93, 36, 7, 31, 10.5), (56, 36, 31, 13, 2.5),
    (317, 11, 44, 34, 11.5),
]
pts_in = [(0.0, 0.0)]
for deg, mn, sec, ft, inch in legs:
    brg = deg + mn/60.0 + sec/3600.0
    dist_in = ft * 12 + inch
    brg_rad = math.radians(brg)
    dE = dist_in * math.sin(brg_rad)
    dN = dist_in * math.cos(brg_rad)
    last = pts_in[-1]
    pts_in.append((last[0] + dE, last[1] + dN))
poly = [(e/12, n/12) for e, n in pts_in[:5]]

# Adjust P3: set P3->P4 bearing = 90 deg (same northing as P4)
P4_orig = poly[3]
P3_new = (-19.1177, P4_orig[1])
poly[2] = P3_new

# Adjust P2: P2->P3 bearing = 180 deg exactly, dist = 29' exactly
P2_new = (P3_new[0], P3_new[1] + 29.0)
poly[1] = P2_new
print(f"Adjusted P2: ({P2_new[0]:.4f}, {P2_new[1]:.4f})")
print(f"Adjusted P3: ({P3_new[0]:.4f}, {P3_new[1]:.4f})")
dE34 = P4_orig[0] - P3_new[0]; dN34 = P4_orig[1] - P3_new[1]
dist34 = math.sqrt(dE34**2 + dN34**2)
print(f"P3-P4: brg 90.0 deg, dist {dist34:.4f} ft")

# === Line geometry (P5-POB line) ===
P5, POB = poly[4], poly[0]
dE_l = P5[0] - POB[0]
dN_l = P5[1] - POB[1]
L = math.sqrt(dE_l**2 + dN_l**2)
uE, uN = dE_l/L, dN_l/L
nE, nN = -uN, uE

# === Arc intersection point ===
R1, R2 = 10.0, 12.5
T1_dist, T2_dist = 26.5, 5.75  # 26'6" and 5'9" from POB
T1 = (T1_dist*uE, T1_dist*uN)
C1 = (T1[0]+R1*nE, T1[1]+R1*nN)
T2 = (T2_dist*uE, T2_dist*uN)
C2 = (T2[0]+R2*nE, T2[1]+R2*nN)

dx = C2[0]-C1[0]; dy = C2[1]-C1[1]
d = math.sqrt(dx**2+dy**2)
a = (R1**2 - R2**2 + d**2)/(2*d)
h = math.sqrt(R1**2 - a**2)
ux, uy = dx/d, dy/d
Mx, My = C1[0]+a*ux, C1[1]+a*uy
I1 = (Mx + h*(-uy), My + h*ux)
I2 = (Mx - h*(-uy), My - h*ux)

ang_T2 = math.atan2(T2[1]-C2[1], T2[0]-C2[0])
def ccw(s,e): return (e-s)%(2*math.pi)
s1 = ccw(ang_T2, math.atan2(I1[1]-C2[1], I1[0]-C2[0]))
s2 = ccw(ang_T2, math.atan2(I2[1]-C2[1], I2[0]-C2[0]))
IX = I1 if s1 < s2 else I2
print(f"Intersection IX: ({IX[0]:.4f}, {IX[1]:.4f})")

# === Arc polyline generation ===
def arc_polyline(cx, cy, r, start_ang, end_ang, n=40):
    pts = []
    for i in range(n+1):
        t = start_ang + (end_ang - start_ang) * i / n
        pts.append((cx + r*math.cos(t), cy + r*math.sin(t)))
    return pts

# Arc 1 sweep (CW from T1 to IX)
ang_T1 = math.atan2(T1[1]-C1[1], T1[0]-C1[0])
ang_IX_C1 = math.atan2(IX[1]-C1[1], IX[0]-C1[0])
sweep1 = (ang_T1 - ang_IX_C1) % (2*math.pi)

# Arc 2 sweep (CCW from T2 to IX)
ang_T2_c = math.atan2(T2[1]-C2[1], T2[0]-C2[0])
ang_IX_C2 = math.atan2(IX[1]-C2[1], IX[0]-C2[0])
sweep2 = (ang_IX_C2 - ang_T2_c) % (2*math.pi)

# Arc 1: CW from T1 to IX
arc1_end_ang = ang_T1 - sweep1
arc1_pts = arc_polyline(C1[0], C1[1], R1, ang_T1, arc1_end_ang, 60)
# Arc 2: CCW from T2 to IX
arc2_end_ang = ang_T2_c + sweep2
arc2_pts = arc_polyline(C2[0], C2[1], R2, ang_T2_c, arc2_end_ang, 60)

print(f"Arc 1 sweep (CW): {math.degrees(sweep1):.1f} deg")
print(f"Arc 2 sweep (CCW): {math.degrees(sweep2):.1f} deg")

# === Build the "available region" polygon ===
# Traverse: POB -> P2 -> P3 -> P4 -> P5 -> (arcs) -> POB
# Replace the P5->POB edge with: P5 -> ... -> T1 -> (arc1 T1->IX) -> IX -> (arc2 IX->T2) -> T2 -> ... -> POB
# The P5->POB line from P5 to T1, then arc1 from T1 to IX, then arc2 from IX to T2, then T2 to POB

# Segment P5 to T1 (both on the P5-POB line)
# T1 is at T1_dist from POB along POB->P5 direction
# P5 is at distance L from POB
# So P5 to T1: we go from P5 back toward POB
# Points on P5-POB line: POB + t*(P5-POB)/L for t from 0 to L
# T1 is at t=T1_dist from POB: (T1_dist*uE, T1_dist*uN)
# T2 is at t=T2_dist from POB
# We want: P5 -> T1 (straight), T1 -> IX (arc1, reversed), IX -> T2 (arc2, reversed), T2 -> POB (straight)

# Wait: arc1 goes T1 to IX, arc2 goes T2 to IX.
# So the boundary is: P5 -> T1 (straight) -> arc1(T1->IX) -> arc2(IX->T2) [reversed] -> T2 -> POB (straight)
# Arc2 goes T2 to IX, so reversed it goes IX to T2.

# Build the available region polygon
avail_region = []
# Vertices: POB, P2, P3, P4, P5
avail_region.append(poly[0])  # POB
avail_region.append(poly[1])  # P2
avail_region.append(poly[2])  # P3
avail_region.append(poly[3])  # P4
avail_region.append(poly[4])  # P5
# P5 to T1 (straight along POB-P5 line) - T1 is between P5 and POB
avail_region.append(T1)
# Arc 1: T1 to IX (already in order)
for pt in arc1_pts[1:]:  # skip T1 (already added)
    avail_region.append(pt)
# Now at IX. Arc 2 goes T2 to IX, so reverse it: IX to T2
arc2_reversed = list(reversed(arc2_pts))
for pt in arc2_reversed[1:]:  # skip IX (already added)
    avail_region.append(pt)
# T2 to POB (straight along POB-P5 line)
# T2 is between POB and P5, close to POB. Just close the polygon back to POB.
# (T2 to POB is a short straight segment, polygon closes automatically)

print(f"Available region polygon: {len(avail_region)} vertices")

# === Point in polygon ===
def point_in_poly(pt, polygon):
    x, y = pt
    n = len(polygon)
    inside = False
    j = n - 1
    for i in range(n):
        xi, yi = polygon[i]
        xj, yj = polygon[j]
        if ((yi > y) != (yj > y)) and (x < (xj-xi)*(y-yi)/(yj-yi)+xi):
            inside = not inside
        j = i
    return inside

# === P4-P5 line clipping ===
P4 = poly[3]
dE45 = P5[0] - P4[0]; dN45 = P5[1] - P4[1]
n45E = -dN45; n45N = dE45
dot_pob = (POB[0]-P4[0])*n45E + (POB[1]-P4[1])*n45N
if dot_pob < 0:
    n45E, n45N = -n45E, -n45N

def side_of_p4p5(pt):
    return (pt[0]-P4[0])*n45E + (pt[1]-P4[1])*n45N

def clip_polygon_by_halfplane(vertices):
    output = []
    n = len(vertices)
    for i in range(n):
        cur = vertices[i]
        nxt = vertices[(i+1) % n]
        d_cur = side_of_p4p5(cur)
        d_nxt = side_of_p4p5(nxt)
        if d_cur >= 0:
            output.append(cur)
        if (d_cur >= 0) != (d_nxt >= 0):
            t = d_cur / (d_cur - d_nxt)
            ix = (cur[0] + t*(nxt[0]-cur[0]), cur[1] + t*(nxt[1]-cur[1]))
            output.append(ix)
    return output

def polygon_area(vertices):
    n = len(vertices)
    if n < 3:
        return 0
    area = 0
    for i in range(n):
        j = (i+1) % n
        area += vertices[i][0]*vertices[j][1]
        area -= vertices[j][0]*vertices[i][1]
    return abs(area) / 2

# === Rectangle helpers ===
def rect_corners(ne, ca, sa, w, hh):
    nw = (ne[0] - w*ca, ne[1] - w*sa)
    se = (ne[0] + hh*sa, ne[1] - hh*ca)
    sw = (nw[0] + hh*sa, nw[1] - hh*ca)
    return ne, nw, sw, se

def is_feasible(ne, ca, sa, w, hh):
    """Check constraints:
    - NE at IX (on boundary of available region) - OK
    - NW and SW must be inside available region (traverse with arc boundary)
    - SE can extend past P4-P5 line but must otherwise be in available region
    - No part of rectangle enters the arcs (enforced by using arc-bounded region)
    """
    ne_pt, nw, sw, se = rect_corners(ne, ca, sa, w, hh)

    # NW must be inside available region
    if not point_in_poly(nw, avail_region):
        return False

    # SW must be inside available region
    if not point_in_poly(sw, avail_region):
        return False

    # SE: either inside available region, or past P4-P5 line
    if not point_in_poly(se, avail_region):
        if side_of_p4p5(se) >= 0:
            return False  # Outside but not past P4-P5

    # Also check midpoints of edges to catch cases where corners are OK
    # but edges cross the arc boundary
    mid_ne_nw = ((ne_pt[0]+nw[0])/2, (ne_pt[1]+nw[1])/2)
    mid_nw_sw = ((nw[0]+sw[0])/2, (nw[1]+sw[1])/2)
    mid_ne_se = ((ne_pt[0]+se[0])/2, (ne_pt[1]+se[1])/2)
    mid_sw_se = ((sw[0]+se[0])/2, (sw[1]+se[1])/2)

    # Top edge (NE-NW) must be inside available region
    if not point_in_poly(mid_ne_nw, avail_region):
        return False

    # Left edge midpoint must be inside
    if not point_in_poly(mid_nw_sw, avail_region):
        return False

    # Right edge (NE-SE): allow past P4-P5
    if not point_in_poly(mid_ne_se, avail_region):
        if side_of_p4p5(mid_ne_se) >= 0:
            return False

    # Bottom edge (SW-SE): allow past P4-P5
    if not point_in_poly(mid_sw_se, avail_region):
        if side_of_p4p5(mid_sw_se) >= 0:
            return False

    # Check quarter points on top edge (NE-NW) for arc crossings
    for frac in [0.25, 0.75]:
        qpt = (ne_pt[0]+frac*(nw[0]-ne_pt[0]), ne_pt[1]+frac*(nw[1]-ne_pt[1]))
        if not point_in_poly(qpt, avail_region):
            return False

    return True

def effective_area(ne, ca, sa, w, hh):
    ne_pt, nw, sw, se = rect_corners(ne, ca, sa, w, hh)
    rect_verts = [ne_pt, nw, sw, se]
    any_outside = any(side_of_p4p5(v) < 0 for v in rect_verts)
    if not any_outside:
        return w * hh
    clipped = clip_polygon_by_halfplane(rect_verts)
    return polygon_area(clipped)

# === Search over angles ===
best_eff_area = 0
best_params = None

for a_deg10 in range(0, 3600):
    angle = math.radians(a_deg10 / 10.0)
    ca, sa = math.cos(angle), math.sin(angle)

    # Binary search for max width (NW inside available region)
    lo_w, hi_w = 0, 55
    for _ in range(50):
        mid_w = (lo_w + hi_w) / 2
        nw = (IX[0] - mid_w*ca, IX[1] - mid_w*sa)
        if point_in_poly(nw, avail_region):
            lo_w = mid_w
        else:
            hi_w = mid_w
    max_w = lo_w

    if max_w < 0.1:
        continue

    for wfrac in [1.0, 0.98, 0.95, 0.92, 0.9, 0.85, 0.8, 0.75, 0.7, 0.65, 0.6, 0.55, 0.5]:
        ww = max_w * wfrac
        lo_h, hi_h = 0, 55
        for _ in range(50):
            mid_h = (lo_h + hi_h) / 2
            if is_feasible(IX, ca, sa, ww, mid_h):
                lo_h = mid_h
            else:
                hi_h = mid_h
        hh = lo_h
        if hh < 0.1:
            continue
        ea = effective_area(IX, ca, sa, ww, hh)
        if ea > best_eff_area:
            best_eff_area = ea
            best_params = (angle, ww, hh, ea)

angle, w, h, eff_area = best_params
ca, sa = math.cos(angle), math.sin(angle)
full_area = w * h
print(f"\nBest rect: angle={math.degrees(angle):.1f} deg, w={w:.2f} ft, h={h:.2f} ft")
print(f"  Full area={full_area:.2f} sq ft, Effective area={eff_area:.2f} sq ft")

# Rectangle corners
NE = IX
NW = (IX[0]-w*ca, IX[1]-w*sa)
SE = (IX[0]+h*sa, IX[1]-h*ca)
SW = (NW[0]+h*sa, NW[1]-h*ca)

# SVG transform
s = (368.79 - 151.26) / 18.66
def to_svg(e, n):
    return (368.79 + e * s, 124.12 - n * s)

for name, pt in [("NE(IX)", NE), ("NW", NW), ("SW", SW), ("SE", SE)]:
    sx, sy = to_svg(*pt)
    ins_t = point_in_poly(pt, poly)
    ins_a = point_in_poly(pt, avail_region)
    past45 = side_of_p4p5(pt) < 0
    print(f"  {name}: ({pt[0]:.2f}, {pt[1]:.2f}) -> svg ({sx:.1f}, {sy:.1f}) in_trav={ins_t} in_avail={ins_a} past_P4P5={past45}")

rect_svg = " ".join(f"{to_svg(*p)[0]:.1f},{to_svg(*p)[1]:.1f}" for p in [NE, NW, SW, SE])
print(f"RECT_SVG: {rect_svg}")
print(f"Dimensions: {w:.2f} ft x {h:.2f} ft")
print(f"Effective area: {eff_area:.2f} sq ft")

# Clipped polygon for SVG
rect_verts = [NE, NW, SW, SE]
clipped = clip_polygon_by_halfplane(rect_verts)
clipped_svg = " ".join(f"{to_svg(*p)[0]:.1f},{to_svg(*p)[1]:.1f}" for p in clipped)
print(f"CLIPPED_SVG: {clipped_svg}")

# Dimension label positions
ne_s = to_svg(*NE); nw_s = to_svg(*NW); se_s = to_svg(*SE); sw_s = to_svg(*SW)
top_mid = ((ne_s[0]+nw_s[0])/2, (ne_s[1]+nw_s[1])/2)
bot_mid = ((se_s[0]+sw_s[0])/2, (se_s[1]+sw_s[1])/2)
left_mid = ((nw_s[0]+sw_s[0])/2, (nw_s[1]+sw_s[1])/2)
right_mid = ((ne_s[0]+se_s[0])/2, (ne_s[1]+se_s[1])/2)

top_ang = math.degrees(math.atan2(nw_s[1]-ne_s[1], nw_s[0]-ne_s[0]))
if top_ang > 90: top_ang -= 180
if top_ang < -90: top_ang += 180
left_ang = math.degrees(math.atan2(sw_s[1]-nw_s[1], sw_s[0]-nw_s[0]))
if left_ang > 90: left_ang -= 180
if left_ang < -90: left_ang += 180

print(f"Top label: ({top_mid[0]:.1f},{top_mid[1]:.1f}) ang={top_ang:.1f}")
print(f"Left label: ({left_mid[0]:.1f},{left_mid[1]:.1f}) ang={left_ang:.1f}")
print(f"Right label: ({right_mid[0]:.1f},{right_mid[1]:.1f})")
print(f"Bot label: ({bot_mid[0]:.1f},{bot_mid[1]:.1f})")

# Center of effective area for label
if len(clipped) >= 3:
    ccx = sum(p[0] for p in clipped) / len(clipped)
    ccy = sum(p[1] for p in clipped) / len(clipped)
    ccx_s, ccy_s = to_svg(ccx, ccy)
    print(f"Effective area center: ({ccx_s:.1f},{ccy_s:.1f})")

# === Arc geometry for SVG ===
print(f"\n--- Arc geometry ---")
print(f"C1: ({C1[0]:.4f}, {C1[1]:.4f})  R1={R1}")
print(f"C2: ({C2[0]:.4f}, {C2[1]:.4f})  R2={R2}")
print(f"T1: ({T1[0]:.4f}, {T1[1]:.4f})  ({T1_dist}' from POB)")
print(f"T2: ({T2[0]:.4f}, {T2[1]:.4f})  ({T2_dist}' from POB)")
print(f"Dist between centers: {d:.4f}")
print(f"Arc 1 sweep (CW): {math.degrees(sweep1):.1f} deg")
print(f"Arc 2 sweep (CCW): {math.degrees(sweep2):.1f} deg")

# Arc SVG polylines (use the 40-point versions for SVG)
arc1_svg_pts = arc_polyline(C1[0], C1[1], R1, ang_T1, arc1_end_ang, 40)
arc2_svg_pts = arc_polyline(C2[0], C2[1], R2, ang_T2_c, arc2_end_ang, 40)

print(f"\nArc 1 polyline (SVG):")
arc1_svg_str = " ".join(f"{to_svg(e,n)[0]:.1f},{to_svg(e,n)[1]:.1f}" for e,n in arc1_svg_pts)
print(arc1_svg_str)

print(f"\nArc 2 polyline (SVG):")
arc2_svg_str = " ".join(f"{to_svg(e,n)[0]:.1f},{to_svg(e,n)[1]:.1f}" for e,n in arc2_svg_pts)
print(arc2_svg_str)

# SVG positions for arc annotations
t1_svg = to_svg(*T1); t2_svg = to_svg(*T2); ix_svg = to_svg(*IX)
c1_svg = to_svg(*C1); c2_svg = to_svg(*C2)
print(f"\nT1 svg: ({t1_svg[0]:.1f},{t1_svg[1]:.1f})")
print(f"T2 svg: ({t2_svg[0]:.1f},{t2_svg[1]:.1f})")
print(f"IX svg: ({ix_svg[0]:.1f},{ix_svg[1]:.1f})")
print(f"C1 svg: ({c1_svg[0]:.1f},{c1_svg[1]:.1f})")
print(f"C2 svg: ({c2_svg[0]:.1f},{c2_svg[1]:.1f})")

# Traverse area
trav_area = polygon_area(poly)
print(f"\nTraverse area: {trav_area:.2f} sq ft")
