import math

# All coordinates are P3-based: P3 = (0, 0).
legs = [
    (257, 53, 45,  19,  1.0),
    (180, 54, 31,  26, 11.0),
    ( 93, 36,  7,  31, 10.5),
    ( 56, 36, 31,  13,  2.5),
    (317, 11, 44,  34, 11.5),
]
_trav = [(0.0, 0.0)]
for deg, mn, sec, ft, inch in legs:
    brg = deg + mn/60.0 + sec/3600.0
    dist_in = ft * 12 + inch
    brg_rad = math.radians(brg)
    dE = dist_in * math.sin(brg_rad)
    dN = dist_in * math.cos(brg_rad)
    last = _trav[-1]
    _trav.append((last[0] + dE, last[1] + dN))

_trav_ft = [(e/12, n/12) for e, n in _trav[:5]]

# Adjust P3 and P2 (horizontal baseline, vertical west face)
_trav_ft[2] = (-19.1177, _trav_ft[3][1])
_trav_ft[1] = (_trav_ft[2][0], _trav_ft[2][1] + 29.0)

# Rebase to P3 = (0, 0)
_p3_trav = _trav_ft[2]
poly = [(_trav_ft[i][0] - _p3_trav[0], _trav_ft[i][1] - _p3_trav[1]) for i in range(5)]

def lerp_edge(a, b, n_val):
    if a[1] == b[1]:
        return None
    if not (min(a[1], b[1]) <= n_val <= max(a[1], b[1])):
        return None
    t = (n_val - a[1]) / (b[1] - a[1])
    return a[0] + t * (b[0] - a[0])

def west_at(n):
    for i, j in [(0,1), (1,2), (2,3)]:
        e = lerp_edge(poly[i], poly[j], n)
        if e is not None:
            return e
    return None

def east_at(n):
    for i, j in [(3,4), (4,0)]:
        e = lerp_edge(poly[i], poly[j], n)
        if e is not None:
            return e
    e = lerp_edge(poly[2], poly[3], n)
    if e is not None and e > 0:
        return e
    return None

# Key layout lines (P3-based, offsets from POB vertex)
H_OFF = poly[0][1] - 7.0          # 7' south of POB
H_WET_SPLIT = poly[0][1] - 13.0   # 13' south of POB
H_MID = poly[0][1] - 17.5         # 17.5' south of POB
H_BOT = poly[0][1] - 27.5         # 27.5' south of POB
V_WET = poly[0][0] - 8.5          # 8.5' west of POB
H_CLOSET = poly[0][1] - 10.0      # 10' south of POB

# Room polygons
office = [poly[0], poly[1], (west_at(H_OFF), H_OFF), (east_at(H_OFF), H_OFF)]
bath_util = [(west_at(H_OFF), H_OFF), (west_at(H_WET_SPLIT), H_WET_SPLIT),
             (V_WET, H_WET_SPLIT), (V_WET, H_OFF)]
bath = [(west_at(H_WET_SPLIT), H_WET_SPLIT), (west_at(H_MID), H_MID),
        (V_WET, H_MID), (V_WET, H_WET_SPLIT)]
bedroom = [(V_WET, H_OFF), (V_WET, H_MID),
           (east_at(H_MID), H_MID), (east_at(H_OFF), H_OFF)]
main_room = [(west_at(H_MID), H_MID), (west_at(H_BOT), H_BOT),
             (east_at(H_BOT), H_BOT), (east_at(H_MID), H_MID)]
storage = [(west_at(H_BOT), H_BOT), poly[2], poly[3], (east_at(H_BOT), H_BOT)]
closet = [(V_WET, H_OFF), (V_WET, H_CLOSET),
          (V_WET + 4.0, H_CLOSET), (V_WET + 4.0, H_OFF)]

def poly_area(pts):
    n = len(pts)
    a = 0
    for i in range(n):
        j = (i+1) % n
        a += pts[i][0] * pts[j][1] - pts[j][0] * pts[i][1]
    return abs(a) / 2.0

# SVG setup - Letter landscape
page_w, page_h = 792, 612
margin = 65
draw_w = page_w - 2 * margin
draw_h = page_h - 2 * margin
es = [p[0] for p in poly]
ns = [p[1] for p in poly]
e_min, e_max = min(es), max(es)
n_min, n_max = min(ns), max(ns)
e_range = e_max - e_min
n_range = n_max - n_min
scale = min(draw_w / e_range, draw_h / n_range) * 0.78
cx = margin + draw_w / 2 + 10
cy = margin + draw_h / 2 + 15
e_center = (e_min + e_max) / 2
n_center = (n_min + n_max) / 2

def to_svg(e, n):
    x = cx + (e - e_center) * scale
    y = cy - (n - n_center) * scale
    return x, y

def pts_str(pts):
    return " ".join(f"{to_svg(e,n)[0]:.1f},{to_svg(e,n)[1]:.1f}" for e, n in pts)

svg = []
svg.append(f'<svg xmlns="http://www.w3.org/2000/svg" width="{page_w}" height="{page_h}" viewBox="0 0 {page_w} {page_h}">')
svg.append(f'<rect width="{page_w}" height="{page_h}" fill="white"/>')
svg.append(f'<text x="{page_w/2}" y="30" text-anchor="middle" font-family="Arial" font-size="15" font-weight="bold">Option A Layout &#8212; ADU Floor Plan within Traverse</text>')

# Traverse outline (light background)
svg.append(f'<polygon points="{pts_str(poly)}" fill="#f5f5f5" stroke="#bbb" stroke-width="1" stroke-dasharray="4,2"/>')

# Room fills
rooms_data = [
    ("OFFICE", office, "#E3F2FD", "#1565C0"),
    ("BATH / UTILITY", bath_util, "#F3E5F5", "#7B1FA2"),
    ("BATH / SHOWER", bath, "#FCE4EC", "#C62828"),
    ("BEDROOM", bedroom, "#E8F5E9", "#2E7D32"),
    ("MAIN ROOM", main_room, "#FFF8E1", "#F57F17"),
    ("STORAGE / MECH", storage, "#EFEBE9", "#5D4037"),
]
for name, pts, fill, sc in rooms_data:
    svg.append(f'<polygon points="{pts_str(pts)}" fill="{fill}" stroke="{sc}" stroke-width="1.5" stroke-linejoin="round"/>')

# Closet overlay
svg.append(f'<polygon points="{pts_str(closet)}" fill="#C8E6C9" stroke="#2E7D32" stroke-width="1" stroke-dasharray="3,2"/>')

# --- Three-arc boundary system (matching gen_path_svg.py) ---
P5 = poly[4]
POB = poly[0]
P3_pt = poly[2]
P4_pt = poly[3]
dE_line = P5[0] - POB[0]
dN_line = P5[1] - POB[1]
L_line = math.sqrt(dE_line**2 + dN_line**2)
uE = dE_line / L_line
uN = dN_line / L_line
nE, nN = -uN, uE  # outward normal

R1, R2, R3 = 10.0, 12.5, 11.0
T1 = (POB[0] + 26.5 * uE, POB[1] + 26.5 * uN)
C1 = (T1[0] + R1 * nE, T1[1] + R1 * nN)
T2 = (POB[0] + 5.75 * uE, POB[1] + 5.75 * uN)
C2 = (T2[0] + R2 * nE, T2[1] + R2 * nN)

# PA: circle-circle intersection of Arc 1 and Arc 2
dx_cc = C2[0] - C1[0]; dy_cc = C2[1] - C1[1]
d_cc = math.sqrt(dx_cc**2 + dy_cc**2)
a_cc = (R1**2 - R2**2 + d_cc**2) / (2 * d_cc)
h_cc = math.sqrt(R1**2 - a_cc**2)
ux_cc, uy_cc = dx_cc / d_cc, dy_cc / d_cc
Mx, My = C1[0] + a_cc * ux_cc, C1[1] + a_cc * uy_cc
I1 = (Mx + h_cc * (-uy_cc), My + h_cc * ux_cc)
I2 = (Mx - h_cc * (-uy_cc), My - h_cc * ux_cc)
ang_T2_C2 = math.atan2(T2[1] - C2[1], T2[0] - C2[0])
_ccw = lambda s, e: (e - s) % (2 * math.pi)
s1 = _ccw(ang_T2_C2, math.atan2(I1[1] - C2[1], I1[0] - C2[0]))
s2 = _ccw(ang_T2_C2, math.atan2(I2[1] - C2[1], I2[0] - C2[0]))
PA = I1 if s1 < s2 else I2

T3 = (P3_pt[0] + 17.911244, P3_pt[1])
C3 = (T3[0], T3[1] - R3)

# PX: intersection of P5-P4 line (extended past P4) with Arc 3
dxL = P4_pt[0] - P5[0]; dyL = P4_pt[1] - P5[1]
ax_ = P5[0] - C3[0]; ay_ = P5[1] - C3[1]
A_ = dxL**2 + dyL**2; B_ = 2 * (ax_ * dxL + ay_ * dyL)
C_ = ax_**2 + ay_**2 - R3**2
disc_ = B_**2 - 4 * A_ * C_
t1_ = (-B_ + math.sqrt(disc_)) / (2 * A_)
t2_ = (-B_ - math.sqrt(disc_)) / (2 * A_)
t_ = min(t for t in [t1_, t2_] if t > 1.0)
PX = (P5[0] + t_ * dxL, P5[1] + t_ * dyL)

# Draw arcs (CW in survey coords -> sweep-flag=0 in SVG)
arc_color = "#0077B6"
for start, end, R in [(T1, PA, R1), (PA, T2, R2), (T3, PX, R3)]:
    sx_s, sy_s = to_svg(*start)
    sx_e, sy_e = to_svg(*end)
    svg_R = R * scale
    svg.append(f'<path d="M {sx_s:.1f},{sy_s:.1f} A {svg_R:.1f},{svg_R:.1f} 0 0 0 '
               f'{sx_e:.1f},{sy_e:.1f}" fill="none" stroke="{arc_color}" stroke-width="2.5"/>')

# Dots at arc junctions
for pt in [T1, PA, T2, T3, PX]:
    sx, sy = to_svg(*pt)
    svg.append(f'<circle cx="{sx:.1f}" cy="{sy:.1f}" r="3.5" fill="{arc_color}"/>')

# Center crosses and dashed radius lines
for ctr, tan_pt in [(C1, T1), (C2, T2), (C3, T3)]:
    sx_c, sy_c = to_svg(*ctr)
    sx_t, sy_t = to_svg(*tan_pt)
    svg.append(f'<line x1="{sx_c:.1f}" y1="{sy_c:.1f}" x2="{sx_t:.1f}" y2="{sy_t:.1f}" '
               f'stroke="{arc_color}" stroke-width="0.6" stroke-dasharray="4,3"/>')
    if 0 < sx_c < page_w and 0 < sy_c < page_h:
        svg.append(f'<line x1="{sx_c-4:.1f}" y1="{sy_c:.1f}" x2="{sx_c+4:.1f}" y2="{sy_c:.1f}" '
                   f'stroke="{arc_color}" stroke-width="0.8"/>')
        svg.append(f'<line x1="{sx_c:.1f}" y1="{sy_c-4:.1f}" x2="{sx_c:.1f}" y2="{sy_c+4:.1f}" '
                   f'stroke="{arc_color}" stroke-width="0.8"/>')

# Radius labels
for ctr, tan_pt, label in [(C1, T1, "R&#8321; = 10&#39;"),
                            (C2, T2, "R&#8322; = 12.5&#39;"),
                            (C3, T3, "R&#8323; = 11&#39;")]:
    sx_c, sy_c = to_svg(*ctr)
    sx_t, sy_t = to_svg(*tan_pt)
    rx_mid = (sx_c + sx_t) / 2
    ry_mid = (sy_c + sy_t) / 2
    r_ang = math.degrees(math.atan2(sy_t - sy_c, sx_t - sx_c))
    if r_ang > 90: r_ang -= 180
    if r_ang < -90: r_ang += 180
    svg.append(f'<text x="{rx_mid:.1f}" y="{ry_mid-4:.1f}" text-anchor="middle" font-family="Arial" '
               f'font-size="8" fill="{arc_color}" '
               f'transform="rotate({r_ang:.1f},{rx_mid:.1f},{ry_mid-4:.1f})">{label}</text>')

# Tangent point labels
sx_t1, sy_t1 = to_svg(*T1)
svg.append(f'<text x="{sx_t1+12:.1f}" y="{sy_t1+12:.1f}" text-anchor="start" font-family="Arial" '
           f'font-size="7.5" fill="{arc_color}">T1 (26.5&#39; from POB)</text>')
sx_pa, sy_pa = to_svg(*PA)
svg.append(f'<text x="{sx_pa+8:.1f}" y="{sy_pa:.1f}" text-anchor="start" font-family="Arial" '
           f'font-size="7.5" fill="{arc_color}">PA</text>')

# Traverse outline bold on top
svg.append(f'<polygon points="{pts_str(poly)}" fill="none" stroke="#333" stroke-width="2.5" stroke-linejoin="round"/>')

# Interior walls
def add_line(e1, n1, e2, n2, w=1.5):
    x1, y1 = to_svg(e1, n1)
    x2, y2 = to_svg(e2, n2)
    svg.append(f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" stroke="#333" stroke-width="{w}"/>')

add_line(west_at(H_OFF), H_OFF, east_at(H_OFF), H_OFF)
add_line(west_at(H_WET_SPLIT), H_WET_SPLIT, V_WET, H_WET_SPLIT)
add_line(west_at(H_MID), H_MID, east_at(H_MID), H_MID)
add_line(west_at(H_BOT), H_BOT, east_at(H_BOT), H_BOT)
add_line(V_WET, H_OFF, V_WET, H_MID)

# Room labels with area
for name, pts, fill, sc in rooms_data:
    area = poly_area(pts)
    ce = sum(p[0] for p in pts) / len(pts)
    cn = sum(p[1] for p in pts) / len(pts)
    x, y = to_svg(ce, cn)
    svg.append(f'<text x="{x:.1f}" y="{y-5:.1f}" text-anchor="middle" font-family="Arial" font-size="11" font-weight="bold" fill="{sc}">{name}</text>')
    svg.append(f'<text x="{x:.1f}" y="{y+9:.1f}" text-anchor="middle" font-family="Arial" font-size="9" fill="{sc}">{area:.0f} sq ft</text>')

# Closet label
x, y = to_svg(V_WET + 2.0, (H_OFF + H_CLOSET) / 2)
svg.append(f'<text x="{x:.1f}" y="{y+3:.1f}" text-anchor="middle" font-family="Arial" font-size="8" fill="#2E7D32">CLOSET</text>')

# Plumbing wall annotation
_plumb_n = H_WET_SPLIT + 1.0
wx, wy = to_svg(west_at(_plumb_n) - 1.8, _plumb_n)
svg.append(f'<text x="{wx:.1f}" y="{wy:.1f}" text-anchor="middle" font-family="Arial" font-size="7" fill="#7B1FA2" transform="rotate(-90,{wx:.1f},{wy:.1f})">plumbing wall</text>')

# Kitchen counter annotation
kce = (west_at(H_BOT + 1.5) + east_at(H_BOT + 1.5)) / 2
x, y = to_svg(kce, H_BOT + 1.5)
svg.append(f'<text x="{x:.1f}" y="{y:.1f}" text-anchor="middle" font-family="Arial" font-size="7.5" fill="#F57F17">kitchen counter along south wall</text>')

# Entry annotation
ex, ey = to_svg(poly[0][0], poly[0][1] + 1.0)
svg.append(f'<text x="{ex:.1f}" y="{ey:.1f}" text-anchor="middle" font-family="Arial" font-size="9" fill="#333">&#9650; ENTRY</text>')

# Fixture annotations
x, y = to_svg((west_at(H_CLOSET) + V_WET)/2, H_CLOSET)
svg.append(f'<text x="{x:.1f}" y="{y+3:.1f}" text-anchor="middle" font-family="Arial" font-size="7.5" fill="#7B1FA2">W/D</text>')

x, y = to_svg((V_WET + east_at(H_WET_SPLIT))/2, H_WET_SPLIT - 0.5)
svg.append(f'<text x="{x:.1f}" y="{y+3:.1f}" text-anchor="middle" font-family="Arial" font-size="7.5" fill="#2E7D32">king bed</text>')

x, y = to_svg((poly[2][0] + poly[3][0])/2, (poly[2][1] + poly[3][1])/2 + 1)
svg.append(f'<text x="{x:.1f}" y="{y+3:.1f}" text-anchor="middle" font-family="Arial" font-size="7.5" fill="#5D4037">W/H + mech</text>')

# Dimension annotations
def dim_label(e1, n1, e2, n2, text, side=1):
    x1, y1 = to_svg(e1, n1)
    x2, y2 = to_svg(e2, n2)
    mx, my = (x1+x2)/2, (y1+y2)/2
    dx, dy = x2-x1, y2-y1
    ln = math.sqrt(dx*dx + dy*dy) or 1
    nx, ny = -dy/ln * side * 12, dx/ln * side * 12
    ang = math.degrees(math.atan2(dy, dx))
    if ang > 90: ang -= 180
    if ang < -90: ang += 180
    svg.append(f'<text x="{mx+nx:.1f}" y="{my+ny:.1f}" text-anchor="middle" font-family="Arial" font-size="7.5" fill="#666" transform="rotate({ang:.1f},{mx+nx:.1f},{my+ny:.1f})">{text}</text>')

# Horizontal wall widths
w = east_at(H_OFF) - west_at(H_OFF)
dim_label(west_at(H_OFF), H_OFF, east_at(H_OFF), H_OFF, f"{w:.0f} ft", -1)

w = east_at(H_MID) - west_at(H_MID)
dim_label(west_at(H_MID), H_MID, east_at(H_MID), H_MID, f"{w:.0f} ft", -1)

# Vertical wall height
h = abs(H_MID - H_OFF)
dim_label(V_WET, H_OFF, V_WET, H_MID, f"{h:.0f} ft", -1)

# Main room depth
h = abs(H_BOT - H_MID)
_mid_n = (H_MID + H_BOT) / 2
dim_label(west_at(_mid_n) - 0.5, H_MID, west_at(_mid_n) - 0.5, H_BOT, f"{h:.0f} ft", -1)

# Vertex dots and labels
vlabels = ["POB", "P2", "P3", "P4", "P5"]
pcx = sum(p[0] for p in poly) / 5
pcy = sum(p[1] for p in poly) / 5
for i, (e, n) in enumerate(poly):
    x, y = to_svg(e, n)
    svg.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="3" fill="#d32f2f"/>')
    dx = e - pcx
    dy = n - pcy
    dist = math.sqrt(dx*dx + dy*dy) or 1
    ox = dx / dist * 2.2
    oy = dy / dist * 2.2
    sx, sy = to_svg(e + ox, n + oy)
    anc = "start" if ox > 0.5 else ("end" if ox < -0.5 else "middle")
    svg.append(f'<text x="{sx:.1f}" y="{sy-4:.1f}" text-anchor="{anc}" font-family="Arial" font-size="9" fill="#d32f2f">{vlabels[i]}</text>')

# North arrow
arx, ary = page_w - 50, page_h - 55
svg.append(f'<line x1="{arx}" y1="{ary+25}" x2="{arx}" y2="{ary-5}" stroke="#333" stroke-width="2" marker-end="url(#ah)"/>')
svg.append('<defs><marker id="ah" markerWidth="8" markerHeight="6" refX="8" refY="3" orient="auto"><polygon points="0 0, 8 3, 0 6" fill="#333"/></marker></defs>')
svg.append(f'<text x="{arx}" y="{ary-11}" text-anchor="middle" font-family="Arial" font-size="12" font-weight="bold">N</text>')

# Legend
lx, ly = 35, page_h - 100
legend = [
    ("OFFICE", "#E3F2FD", "#1565C0"),
    ("BATH / UTILITY", "#F3E5F5", "#7B1FA2"),
    ("BATH / SHOWER", "#FCE4EC", "#C62828"),
    ("BEDROOM", "#E8F5E9", "#2E7D32"),
    ("MAIN ROOM", "#FFF8E1", "#F57F17"),
    ("STORAGE / MECH", "#EFEBE9", "#5D4037"),
]
svg.append(f'<text x="{lx}" y="{ly-5}" font-family="Arial" font-size="9" font-weight="bold" fill="#333">Legend</text>')
for idx, (name, fill, sc) in enumerate(legend):
    ry = ly + idx * 14
    svg.append(f'<rect x="{lx}" y="{ry}" width="12" height="10" fill="{fill}" stroke="{sc}" stroke-width="1"/>')
    svg.append(f'<text x="{lx+16}" y="{ry+8.5}" font-family="Arial" font-size="8" fill="#333">{name}</text>')

# Footer
svg.append(f'<text x="{page_w/2}" y="{page_h-10}" text-anchor="middle" font-family="Arial" font-size="8" fill="#999">Traverse: {poly_area(poly):.0f} sq ft &#8226; Main room in wide south, wet rooms stacked on west plumbing wall, office in narrow north &#8226; Letter landscape</text>')
svg.append('</svg>')

outpath = r"c:\Users\Mango Cat\Dev\hut2\option_a_layout.svg"
with open(outpath, "w") as f:
    f.write("\n".join(svg))

print(f"Written to {outpath}")
for name, pts, _, _ in rooms_data:
    print(f"  {name:<18} {poly_area(pts):.0f} sq ft")
