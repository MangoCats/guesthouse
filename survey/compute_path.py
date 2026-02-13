import sys, os, math

_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _root not in sys.path:
    sys.path.insert(0, _root)

from shared.geometry import poly_area, arc_poly, brg_dist, fmt_brg, fmt_dist
from shared.survey import compute_traverse, compute_three_arc

# === Traverse (raw survey data, instrument at POB) ===
# All coordinates below are P3-based: P3 = (0, 0).
_pts, _p3_trav = compute_traverse()
P3  = _pts["P3"]
POB = _pts["POB"]
P2  = _pts["P2"]
P4  = _pts["P4"]
P5  = _pts["P5"]

print(f"P3:  ({P3[0]:.4f}, {P3[1]:.4f})")
print(f"POB: ({POB[0]:.4f}, {POB[1]:.4f})")
print(f"P2:  ({P2[0]:.4f}, {P2[1]:.4f})")
print(f"P4:  ({P4[0]:.4f}, {P4[1]:.4f})")
print(f"P5:  ({P5[0]:.4f}, {P5[1]:.4f})")

# === Three-arc system ===
_arc = compute_three_arc(_pts)
R1, R2, R3 = _arc["R1"], _arc["R2"], _arc["R3"]
uE, uN = _arc["uE"], _arc["uN"]
nE, nN = _arc["nE"], _arc["nN"]
T1, TC1 =_pts["T1"], _pts["TC1"]
T2, TC2 =_pts["T2"], _pts["TC2"]
T3, TC3 =_pts["T3"], _pts["TC3"]
PA = _pts["PA"]
PX = _pts["PX"]
T1_dist, T2_dist = 26.5, 5.75

print(f"\nT1:  ({T1[0]:.4f}, {T1[1]:.4f})  [Arc 1 tangent, {T1_dist}' from POB]")
print(f"T2:  ({T2[0]:.4f}, {T2[1]:.4f})  [Arc 2 tangent, {T2_dist}' from POB]")
print(f"C1:  ({TC1[0]:.4f}, {TC1[1]:.4f})  R1={R1}")
print(f"C2:  ({TC2[0]:.4f}, {TC2[1]:.4f})  R2={R2}")
print(f"PA:  ({PA[0]:.4f}, {PA[1]:.4f})  [Arc 1/2 intersection]")

print(f"\nT3:  ({T3[0]:.4f}, {T3[1]:.4f})  [Arc 3 tangent, 18' E of P3]")
print(f"C3:  ({TC3[0]:.4f}, {TC3[1]:.4f})  R3={R3}")

# === PX diagnostics ===
dxL = P4[0]-P5[0]; dyL = P4[1]-P5[1]
ax_coef = P5[0] - TC3[0]
ay_coef = P5[1] - TC3[1]
A = dxL**2 + dyL**2
B = 2*(ax_coef*dxL + ay_coef*dyL)
C = ax_coef**2 + ay_coef**2 - R3**2
disc = B**2 - 4*A*C
print(f"\nP5-P4 line extended, circle intersection:")
print(f"  Discriminant: {disc:.4f}")

t1_sol = (-B + math.sqrt(disc))/(2*A)
t2_sol = (-B - math.sqrt(disc))/(2*A)
print(f"  t1={t1_sol:.4f}, t2={t2_sol:.4f}  (t=1 is P4)")

t_px = min(t for t in [t1_sol, t2_sol] if t > 1)
print(f"  PX: ({PX[0]:.4f}, {PX[1]:.4f})  [t={t_px:.4f}]")

dist_px_c3 = math.sqrt((PX[0]-TC3[0])**2 + (PX[1]-TC3[1])**2)
print(f"  PX dist from C3: {dist_px_c3:.4f} (should be {R3})")

# === Arc sweep angles ===
# Arc 1: CW from T1 to PA
ang_T1_C1 = math.atan2(T1[1]-TC1[1], T1[0]-TC1[0])
ang_PA_C1 = math.atan2(PA[1]-TC1[1], PA[0]-TC1[0])
sweep1 = (ang_T1_C1 - ang_PA_C1) % (2*math.pi)  # CW sweep
arc1_end = ang_T1_C1 - sweep1
arc1_pts = arc_poly(TC1[0], TC1[1], R1, ang_T1_C1, arc1_end, 60)
print(f"\nArc 1: CW from T1 to PA, sweep={math.degrees(sweep1):.1f} deg")

# Arc 2: CCW from T2 to PA (for path, we go PA to T2, which is CW)
ang_T2_C2 = math.atan2(T2[1]-TC2[1], T2[0]-TC2[0])
ang_PA_C2 = math.atan2(PA[1]-TC2[1], PA[0]-TC2[0])
sweep2 = (ang_PA_C2 - ang_T2_C2) % (2*math.pi)  # CCW sweep T2->PA
# For path: PA to T2 (reverse of arc2), which is CW
arc2_fwd_pts = arc_poly(TC2[0], TC2[1], R2, ang_T2_C2, ang_T2_C2 + sweep2, 60)
arc2_rev_pts = list(reversed(arc2_fwd_pts))  # PA to T2
print(f"Arc 2: CCW from T2 to PA, sweep={math.degrees(sweep2):.1f} deg")

# Arc 3: from T3 to PX, center south, going CW (east then south)
ang_T3_C3 = math.atan2(T3[1]-TC3[1], T3[0]-TC3[0])  # should be pi/2 (north)
ang_PX_C3 = math.atan2(PX[1]-TC3[1], PX[0]-TC3[0])
# CW sweep from T3 to PX
sweep3 = (ang_T3_C3 - ang_PX_C3) % (2*math.pi)
arc3_end = ang_T3_C3 - sweep3
arc3_pts = arc_poly(TC3[0], TC3[1], R3, ang_T3_C3, arc3_end, 60)
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
area = poly_area(path_pts)
print(f"\nPath area: {area:.2f} sq ft")
print(f"Path vertices: {len(path_pts)}")

# === Straight segment bearings and distances ===
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
    brg, dist = brg_dist(p1, p2)
    print(f"  {name1} -> {name2}: brg {fmt_brg(brg)}, dist {fmt_dist(dist)} ({dist:.4f}')")

# === Arc segment info ===
print(f"\n=== Arc segments ===")
arc_len1 = R1 * sweep1
arc_len2 = R2 * sweep2
arc_len3 = R3 * sweep3
print(f"  Arc 1 (T1->PA): R={R1}', sweep={math.degrees(sweep1):.1f} deg, length={arc_len1:.2f}'")
print(f"  Arc 2 (PA->T2): R={R2}', sweep={math.degrees(sweep2):.1f} deg, length={arc_len2:.2f}'")
print(f"  Arc 3 (T3->PX): R={R3}', sweep={math.degrees(sweep3):.1f} deg, length={arc_len3:.2f}'")

# === SVG coordinates (P3-based) ===
s = (368.79 - 151.26) / 18.66
_p3_svg_x = 368.79 + _p3_trav[0] * s
_p3_svg_y = 124.12 - _p3_trav[1] * s
def to_svg(e, n):
    return (_p3_svg_x + e * s, _p3_svg_y - n * s)

print(f"\n=== SVG coordinates ===")
for name, pt in [("POB", POB), ("P2", P2), ("P3", P3), ("P4", P4), ("P5", P5),
                  ("T1", T1), ("T2", T2), ("T3", T3), ("PA", PA), ("PX", PX),
                  ("TC1", TC1), ("TC2", TC2), ("TC3", TC3)]:
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
    brg, dist = brg_dist(p1, p2)
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
all_pts = [POB, P2, P3, P4, P5, T1, T2, T3, PA, PX, TC3]
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
