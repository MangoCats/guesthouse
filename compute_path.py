import math

# === Traverse ===
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

P4_orig = poly[3]
P3_new = (-19.1177, P4_orig[1])
poly[2] = P3_new
P2_new = (P3_new[0], P3_new[1] + 29.0)
poly[1] = P2_new

POB, P2, P3, P4, P5 = poly[0], poly[1], poly[2], poly[3], poly[4]
print(f"POB: ({POB[0]:.4f}, {POB[1]:.4f})")
print(f"P2:  ({P2[0]:.4f}, {P2[1]:.4f})")
print(f"P3:  ({P3[0]:.4f}, {P3[1]:.4f})")
print(f"P4:  ({P4[0]:.4f}, {P4[1]:.4f})")
print(f"P5:  ({P5[0]:.4f}, {P5[1]:.4f})")

# === P5-POB line geometry ===
dE_l = P5[0] - POB[0]; dN_l = P5[1] - POB[1]
L = math.sqrt(dE_l**2 + dN_l**2)
uE, uN = dE_l/L, dN_l/L
nE, nN = -uN, uE  # left normal of POB->P5

# === Arcs 1 & 2 (on P5-POB line) ===
R1, R2 = 10.0, 12.5
T1_dist, T2_dist = 26.5, 5.75
T1 = (T1_dist*uE, T1_dist*uN)
C1 = (T1[0]+R1*nE, T1[1]+R1*nN)
T2 = (T2_dist*uE, T2_dist*uN)
C2 = (T2[0]+R2*nE, T2[1]+R2*nN)

# Arc intersection PA
dx = C2[0]-C1[0]; dy = C2[1]-C1[1]
d_cc = math.sqrt(dx**2+dy**2)
a = (R1**2 - R2**2 + d_cc**2)/(2*d_cc)
h = math.sqrt(R1**2 - a**2)
ux, uy = dx/d_cc, dy/d_cc
Mx, My = C1[0]+a*ux, C1[1]+a*uy
I1 = (Mx + h*(-uy), My + h*ux)
I2 = (Mx - h*(-uy), My - h*ux)

ang_T2_from_C2 = math.atan2(T2[1]-C2[1], T2[0]-C2[0])
def ccw_angle(s,e): return (e-s)%(2*math.pi)
s1 = ccw_angle(ang_T2_from_C2, math.atan2(I1[1]-C2[1], I1[0]-C2[0]))
s2 = ccw_angle(ang_T2_from_C2, math.atan2(I2[1]-C2[1], I2[0]-C2[0]))
PA = I1 if s1 < s2 else I2

print(f"\nT1:  ({T1[0]:.4f}, {T1[1]:.4f})  [Arc 1 tangent, {T1_dist}' from POB]")
print(f"T2:  ({T2[0]:.4f}, {T2[1]:.4f})  [Arc 2 tangent, {T2_dist}' from POB]")
print(f"C1:  ({C1[0]:.4f}, {C1[1]:.4f})  R1={R1}")
print(f"C2:  ({C2[0]:.4f}, {C2[1]:.4f})  R2={R2}")
print(f"PA:  ({PA[0]:.4f}, {PA[1]:.4f})  [Arc 1/2 intersection]")

# === Arc 3 (on P3-P4 line) ===
R3 = 11.0
# T3: positioned east of P3 so that P5-PX = 20' exactly
# Solved: T3 offset = 17.911244' from P3
T3_dist_from_P3 = 17.911244
T3 = (P3[0] + T3_dist_from_P3, P3[1])
# C3: south of P3-P4 line by R3
C3 = (T3[0], T3[1] - R3)
print(f"\nT3:  ({T3[0]:.4f}, {T3[1]:.4f})  [Arc 3 tangent, 18' E of P3]")
print(f"C3:  ({C3[0]:.4f}, {C3[1]:.4f})  R3={R3}")

# === PX: intersection of extended P5-P4 line with Arc 3 ===
# Line from P5 through P4, parametric: P(t) = P5 + t*(P4-P5), t=1 at P4
dxL = P4[0]-P5[0]; dyL = P4[1]-P5[1]
# Substitute into circle equation: (x-C3x)^2 + (y-C3y)^2 = R3^2
# x = P5x + t*dxL, y = P5y + t*dyL
# (P5x + t*dxL - C3x)^2 + (P5y + t*dyL - C3y)^2 = R3^2
ax_coef = P5[0] - C3[0]
ay_coef = P5[1] - C3[1]
A = dxL**2 + dyL**2
B = 2*(ax_coef*dxL + ay_coef*dyL)
C = ax_coef**2 + ay_coef**2 - R3**2
disc = B**2 - 4*A*C
print(f"\nP5-P4 line extended, circle intersection:")
print(f"  Discriminant: {disc:.4f}")

t1_sol = (-B + math.sqrt(disc))/(2*A)
t2_sol = (-B - math.sqrt(disc))/(2*A)
print(f"  t1={t1_sol:.4f}, t2={t2_sol:.4f}  (t=1 is P4)")

# PX is the FIRST intersection past P4 (smallest t > 1)
candidates = [t for t in [t1_sol, t2_sol] if t > 1]
t_px = min(candidates)
PX = (P5[0] + t_px*dxL, P5[1] + t_px*dyL)
print(f"  PX: ({PX[0]:.4f}, {PX[1]:.4f})  [t={t_px:.4f}]")

# Verify PX is on circle
dist_px_c3 = math.sqrt((PX[0]-C3[0])**2 + (PX[1]-C3[1])**2)
print(f"  PX dist from C3: {dist_px_c3:.4f} (should be {R3})")

# === Arc sweep angles ===
def arc_polyline(cx, cy, r, start_ang, end_ang, n=60):
    pts = []
    for i in range(n+1):
        t = start_ang + (end_ang - start_ang) * i / n
        pts.append((cx + r*math.cos(t), cy + r*math.sin(t)))
    return pts

# Arc 1: CW from T1 to PA
ang_T1_C1 = math.atan2(T1[1]-C1[1], T1[0]-C1[0])
ang_PA_C1 = math.atan2(PA[1]-C1[1], PA[0]-C1[0])
sweep1 = (ang_T1_C1 - ang_PA_C1) % (2*math.pi)  # CW sweep
arc1_end = ang_T1_C1 - sweep1
arc1_pts = arc_polyline(C1[0], C1[1], R1, ang_T1_C1, arc1_end, 60)
print(f"\nArc 1: CW from T1 to PA, sweep={math.degrees(sweep1):.1f} deg")

# Arc 2: CCW from T2 to PA (for path, we go PA to T2, which is CW)
ang_T2_C2 = math.atan2(T2[1]-C2[1], T2[0]-C2[0])
ang_PA_C2 = math.atan2(PA[1]-C2[1], PA[0]-C2[0])
sweep2 = (ang_PA_C2 - ang_T2_C2) % (2*math.pi)  # CCW sweep T2->PA
# For path: PA to T2 (reverse of arc2), which is CW
arc2_fwd_pts = arc_polyline(C2[0], C2[1], R2, ang_T2_C2, ang_T2_C2 + sweep2, 60)
arc2_rev_pts = list(reversed(arc2_fwd_pts))  # PA to T2
print(f"Arc 2: CCW from T2 to PA, sweep={math.degrees(sweep2):.1f} deg")

# Arc 3: from T3 to PX, center south, going CW (east then south)
ang_T3_C3 = math.atan2(T3[1]-C3[1], T3[0]-C3[0])  # should be pi/2 (north)
ang_PX_C3 = math.atan2(PX[1]-C3[1], PX[0]-C3[0])
# CW sweep from T3 to PX
sweep3 = (ang_T3_C3 - ang_PX_C3) % (2*math.pi)
arc3_end = ang_T3_C3 - sweep3
arc3_pts = arc_polyline(C3[0], C3[1], R3, ang_T3_C3, arc3_end, 60)
print(f"Arc 3: CW from T3 to PX, sweep={math.degrees(sweep3):.1f} deg")
print(f"  T3 angle from C3: {math.degrees(ang_T3_C3):.1f} deg")
print(f"  PX angle from C3: {math.degrees(ang_PX_C3):.1f} deg")

# === Build full path polygon (dense, for area calculation) ===
path_pts = []
# 1. POB
path_pts.append(POB)
# 2. P2
path_pts.append(P2)
# 3. P3
path_pts.append(P3)
# 4. T3
path_pts.append(T3)
# 5. Arc 3: T3 to PX (CW)
for pt in arc3_pts[1:]:  # skip T3
    path_pts.append(pt)
# 6. PX to P4 (straight, back along P5-P4 extension) -- wait, check if needed
# Actually PX is past P4 on P5-P4 extension. But the path says PX - P4.
# So yes, PX straight to P4.

# But wait: does PX → P4 make sense? Let me check.
# PX is on the far side of P4 from P5. So PX is further south/west than P4.
# The path goes: Arc3 ends at PX, then straight to P4, then P4 to P5.
# PX → P4 is back along the P5-P4 line toward P4.
print(f"\nPX to P4 distance: {math.sqrt((PX[0]-P4[0])**2+(PX[1]-P4[1])**2):.4f}")

path_pts.append(P4)
# 7. P4 to P5
path_pts.append(P5)
# 8. P5 to T1 (straight along P5-POB line)
path_pts.append(T1)
# 9. Arc 1: T1 to PA (CW)
for pt in arc1_pts[1:]:  # skip T1
    path_pts.append(pt)
# 10. Arc 2: PA to T2 (reversed)
for pt in arc2_rev_pts[1:]:  # skip PA
    path_pts.append(pt)
# 11. T2 to POB (straight)
# Polygon closes back to POB

# === Compute area (shoelace) ===
def polygon_area(vertices):
    n = len(vertices)
    area = 0
    for i in range(n):
        j = (i+1) % n
        area += vertices[i][0]*vertices[j][1]
        area -= vertices[j][0]*vertices[i][1]
    return abs(area) / 2

area = polygon_area(path_pts)
print(f"\nPath area: {area:.2f} sq ft")
print(f"Path vertices: {len(path_pts)}")

# === Straight segment bearings and distances ===
def bearing_and_dist(p1, p2):
    dE = p2[0]-p1[0]; dN = p2[1]-p1[1]
    dist = math.sqrt(dE**2+dN**2)
    brg_rad = math.atan2(dE, dN)
    brg_deg = math.degrees(brg_rad) % 360
    return brg_deg, dist

def fmt_brg(brg_deg):
    d = int(brg_deg)
    m = int((brg_deg - d) * 60)
    s = (brg_deg - d - m/60) * 3600
    return f"{d:03d}\u00b0 {m:02d}' {s:04.1f}\""

def fmt_dist(ft):
    feet = int(ft)
    inches = (ft - feet) * 12
    return f"{feet}' {inches:.1f}\""

print(f"\n=== Straight segments ===")
straight_segs = [
    ("POB", "P2", POB, P2),
    ("P2", "P3", P2, P3),
    ("P3", "T3", P3, T3),
    ("PX", "P4", PX, P4),
    ("P4", "P5", P4, P5),
    ("P5", "T1", P5, T1),
    ("T2", "POB", T2, POB),
]
for name1, name2, p1, p2 in straight_segs:
    brg, dist = bearing_and_dist(p1, p2)
    print(f"  {name1} -> {name2}: brg {fmt_brg(brg)}, dist {fmt_dist(dist)} ({dist:.4f}')")

# === Arc segment info ===
print(f"\n=== Arc segments ===")
arc_len1 = R1 * sweep1
arc_len2 = R2 * sweep2
arc_len3 = R3 * sweep3
print(f"  Arc 1 (T1->PA): R={R1}', sweep={math.degrees(sweep1):.1f} deg, length={arc_len1:.2f}'")
print(f"  Arc 2 (PA->T2): R={R2}', sweep={math.degrees(sweep2):.1f} deg, length={arc_len2:.2f}'")
print(f"  Arc 3 (T3->PX): R={R3}', sweep={math.degrees(sweep3):.1f} deg, length={arc_len3:.2f}'")

# === SVG coordinates ===
s = (368.79 - 151.26) / 18.66
def to_svg(e, n):
    return (368.79 + e * s, 124.12 - n * s)

print(f"\n=== SVG coordinates ===")
for name, pt in [("POB", POB), ("P2", P2), ("P3", P3), ("P4", P4), ("P5", P5),
                  ("T1", T1), ("T2", T2), ("T3", T3), ("PA", PA), ("PX", PX),
                  ("C1", C1), ("C2", C2), ("C3", C3)]:
    sx, sy = to_svg(*pt)
    print(f"  {name}: ({sx:.1f}, {sy:.1f})")

# Traverse polygon SVG
trav_svg = " ".join(f"{to_svg(*p)[0]:.1f},{to_svg(*p)[1]:.1f}" for p in [POB, P2, P3, P4, P5])
print(f"\nTraverse SVG: {trav_svg}")

# Path polygon SVG (dense)
path_svg = " ".join(f"{to_svg(*p)[0]:.1f},{to_svg(*p)[1]:.1f}" for p in path_pts)
print(f"\nPath SVG points: {len(path_pts)}")

# Arc polylines for SVG
print(f"\nArc 1 SVG (T1->PA):")
a1_svg = " ".join(f"{to_svg(e,n)[0]:.1f},{to_svg(e,n)[1]:.1f}" for e,n in arc1_pts)
print(a1_svg)

print(f"\nArc 2 SVG (PA->T2, reversed):")
a2r_svg = " ".join(f"{to_svg(e,n)[0]:.1f},{to_svg(e,n)[1]:.1f}" for e,n in arc2_rev_pts)
print(a2r_svg)

print(f"\nArc 3 SVG (T3->PX):")
a3_svg = " ".join(f"{to_svg(e,n)[0]:.1f},{to_svg(e,n)[1]:.1f}" for e,n in arc3_pts)
print(a3_svg)

# Segment midpoints for labels (SVG)
print(f"\n=== Label positions (SVG midpoints) ===")
for name1, name2, p1, p2 in straight_segs:
    mid = ((p1[0]+p2[0])/2, (p1[1]+p2[1])/2)
    sx, sy = to_svg(*mid)
    brg, dist = bearing_and_dist(p1, p2)
    ang_svg = math.degrees(math.atan2(-(to_svg(*p2)[1]-to_svg(*p1)[1]), to_svg(*p2)[0]-to_svg(*p1)[0]))
    if ang_svg > 90: ang_svg -= 180
    if ang_svg < -90: ang_svg += 180
    print(f"  {name1}-{name2}: mid=({sx:.1f},{sy:.1f}) ang={ang_svg:.1f}")

# Arc label positions (at midpoint of arc)
for name, pts_list in [("Arc1", arc1_pts), ("Arc2_rev", arc2_rev_pts), ("Arc3", arc3_pts)]:
    mid_idx = len(pts_list) // 2
    mid = pts_list[mid_idx]
    sx, sy = to_svg(*mid)
    print(f"  {name} mid: ({sx:.1f},{sy:.1f})")

# Check the SVG bounds needed
all_pts = [POB, P2, P3, P4, P5, T1, T2, T3, PA, PX, C3]
for pt in arc3_pts:
    all_pts.append(pt)
min_e = min(p[0] for p in all_pts)
max_e = max(p[0] for p in all_pts)
min_n = min(p[1] for p in all_pts)
max_n = max(p[1] for p in all_pts)
print(f"\nSurvey bounds: E=[{min_e:.2f}, {max_e:.2f}], N=[{min_n:.2f}, {max_n:.2f}]")
min_sx, min_sy = to_svg(min_e, max_n)
max_sx, max_sy = to_svg(max_e, min_n)
print(f"SVG bounds: x=[{min_sx:.1f}, {max_sx:.1f}], y=[{min_sy:.1f}, {max_sy:.1f}]")
