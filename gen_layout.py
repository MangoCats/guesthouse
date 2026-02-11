import math

legs = [
    (257, 53, 45,  19,  1.0),
    (180, 54, 31,  26, 11.0),
    ( 93, 36,  7,  31, 10.5),
    ( 56, 36, 31,  13,  2.5),
    (317, 11, 44,  34, 11.5),
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

# Rebase origin to P4
_p4 = poly[3]
poly = [(e - _p4[0], n - _p4[1]) for e, n in poly]

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

# Key layout lines (shifted to P4 origin)
H_OFF = -7.0 - _p4[1]
H_WET_SPLIT = -13.0 - _p4[1]
H_MID = -17.5 - _p4[1]
H_BOT = -27.5 - _p4[1]
V_WET = -8.5 - _p4[0]
H_CLOSET = -10.0 - _p4[1]

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

# --- Circular arc: R=12', tangent to P5->POB line, 26' from POB, 45 deg CW ---
P5 = poly[4]
POB = poly[0]
dE_line = P5[0] - POB[0]
dN_line = P5[1] - POB[1]
L_line = math.sqrt(dE_line**2 + dN_line**2)
uE_line = dE_line / L_line  # unit POB -> P5
uN_line = dN_line / L_line

# Tangent point: 26' from POB toward P5
arc_T = (POB[0] + 26 * uE_line, POB[1] + 26 * uN_line)

# Outward normal (NE, away from polygon interior)
arc_nE = -uN_line  # outward perpendicular
arc_nN = uE_line

# Arc center: 12' outward from tangent point
arc_R = 12.0
arc_C = (arc_T[0] + arc_R * arc_nE, arc_T[1] + arc_R * arc_nN)

# Start angle (from center to tangent point)
arc_start_ang = math.atan2(arc_T[1] - arc_C[1], arc_T[0] - arc_C[0])

# End angle: 45 deg clockwise (decreasing math angle)
arc_end_ang = arc_start_ang - math.radians(45)
arc_end = (arc_C[0] + arc_R * math.cos(arc_end_ang),
           arc_C[1] + arc_R * math.sin(arc_end_ang))

# Convert to SVG coordinates
sx_start, sy_start = to_svg(arc_T[0], arc_T[1])
sx_end, sy_end = to_svg(arc_end[0], arc_end[1])
svg_R = arc_R * scale  # radius in SVG pixels

# SVG arc: CW in survey coords -> since Y is flipped in SVG, use sweep-flag=0
# 45 deg < 180 deg -> large-arc-flag=0
svg.append(f'<path d="M {sx_start:.1f},{sy_start:.1f} A {svg_R:.1f},{svg_R:.1f} 0 0 0 {sx_end:.1f},{sy_end:.1f}" '
           f'fill="none" stroke="#0077B6" stroke-width="2.5"/>')

# Dot at tangent point and end point
svg.append(f'<circle cx="{sx_start:.1f}" cy="{sy_start:.1f}" r="3.5" fill="#0077B6"/>')
svg.append(f'<circle cx="{sx_end:.1f}" cy="{sy_end:.1f}" r="3.5" fill="#0077B6"/>')

# Center cross
sx_c, sy_c = to_svg(arc_C[0], arc_C[1])
svg.append(f'<line x1="{sx_c-4:.1f}" y1="{sy_c:.1f}" x2="{sx_c+4:.1f}" y2="{sy_c:.1f}" stroke="#0077B6" stroke-width="0.8"/>')
svg.append(f'<line x1="{sx_c:.1f}" y1="{sy_c-4:.1f}" x2="{sx_c:.1f}" y2="{sy_c+4:.1f}" stroke="#0077B6" stroke-width="0.8"/>')

# Radius line from center to tangent point (dashed)
svg.append(f'<line x1="{sx_c:.1f}" y1="{sy_c:.1f}" x2="{sx_start:.1f}" y2="{sy_start:.1f}" '
           f'stroke="#0077B6" stroke-width="0.6" stroke-dasharray="4,3"/>')

# Labels
# "R=12'" near midpoint of radius line
rx_mid, ry_mid = (sx_c + sx_start) / 2, (sy_c + sy_start) / 2
r_ang = math.degrees(math.atan2(sy_start - sy_c, sx_start - sx_c))
if r_ang > 90: r_ang -= 180
if r_ang < -90: r_ang += 180
svg.append(f'<text x="{rx_mid+8:.1f}" y="{ry_mid-4:.1f}" text-anchor="middle" font-family="Arial" '
           f'font-size="8" fill="#0077B6" transform="rotate({r_ang:.1f},{rx_mid+8:.1f},{ry_mid-4:.1f})">R = 12&#39;</text>')

# "tangent pt" label
svg.append(f'<text x="{sx_start+12:.1f}" y="{sy_start+12:.1f}" text-anchor="start" font-family="Arial" '
           f'font-size="7.5" fill="#0077B6">tangent pt (26&#39; from POB)</text>')

# "45 deg" arc angle label near the arc midpoint
arc_mid_ang = arc_start_ang - math.radians(22.5)
arc_mid_pt = (arc_C[0] + (arc_R + 1.5) * math.cos(arc_mid_ang),
              arc_C[1] + (arc_R + 1.5) * math.sin(arc_mid_ang))
sx_m, sy_m = to_svg(arc_mid_pt[0], arc_mid_pt[1])
svg.append(f'<text x="{sx_m:.1f}" y="{sy_m:.1f}" text-anchor="middle" font-family="Arial" '
           f'font-size="8" fill="#0077B6">45&#176;</text>')

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
svg.append(f'<text x="{page_w/2}" y="{page_h-10}" text-anchor="middle" font-family="Arial" font-size="8" fill="#999">Traverse: 989 sq ft &#8226; Main room in wide south, wet rooms stacked on west plumbing wall, office in narrow north &#8226; Letter landscape</text>')
svg.append('</svg>')

outpath = r"c:\Users\Mango Cat\Dev\hut2\option_a_layout.svg"
with open(outpath, "w") as f:
    f.write("\n".join(svg))

print(f"Written to {outpath}")
for name, pts, _, _ in rooms_data:
    print(f"  {name:<18} {poly_area(pts):.0f} sq ft")
