"""Generate floorplan SVG with 8" wall inset from the outline path.

Runs gen_path_svg.py to obtain outline geometry, then computes an 8" inset.
Points C0-C15 correspond to outline points O0-O15.

Inset rules (CCW path, interior to left):
  - CW arcs  (center outside shape): increase radius by wall thickness
  - CCW arcs (center inside shape):  decrease radius by wall thickness
  - Line segments: offset inward along left normal
"""
import math, sys, os, io

# --- Load outline geometry by executing gen_path_svg.py ---
_parent = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_old_stdout = sys.stdout
sys.stdout = io.StringIO()
_ns = {}
exec(open(os.path.join(_parent, "gen_path_svg.py")).read(), _ns)
sys.stdout = _old_stdout

pts          = _ns["pts"]
outline_segs = _ns["outline_segs"]
LineSeg      = _ns["LineSeg"]
ArcSeg       = _ns["ArcSeg"]
left_norm    = _ns["left_norm"]
off_pt       = _ns["off_pt"]
to_svg       = _ns["to_svg"]
segment_polyline = _ns["segment_polyline"]
path_polygon     = _ns["path_polygon"]
poly_area        = _ns["poly_area"]
W, H         = _ns["W"], _ns["H"]

# Outline arc radii
R_fillet  = _ns["R_fillet"]     # Cf   CCW (inside)
R3i       = _ns["R3i"]          # C3   CW  (outside)
R_f_pox   = _ns["R_f_pox"]     # Cf3  CCW (inside)
R_f_po5   = _ns["R_f_po5"]     # Cf4  CCW (inside)
R1i       = _ns["R1i"]          # C1   CW  (outside)
R_turn3   = _ns["R_turn3"]     # Ct3  CCW (inside)
R_turn2   = _ns["R_turn2"]     # Ct2  CW  (outside)
R_turn1   = _ns["R_turn1"]     # Ct1  CCW (inside)
R_fillet2 = _ns["R_fillet2"]   # Cf2  CCW (inside)

# --- Wall thickness ---
wall_t = 8.0 / 12.0  # 8 inches in feet

# --- Compute inner wall points ---
# At each junction, tangency is preserved by the inset:
#   - Line offset and arc radius adjustment keep center-to-line distance = adjusted radius
#   - Arc-arc external tangency is preserved since +wall_t and -wall_t cancel

def _inner_point(seg_before, seg_after):
    """Inner wall point at junction of seg_before and seg_after."""
    b_line = isinstance(seg_before, LineSeg)
    a_line = isinstance(seg_after, LineSeg)

    if not b_line and not a_line:
        # Arc-Arc: new tangent point on the line between centers
        c1 = pts[seg_before.center]; c2 = pts[seg_after.center]
        r1 = (seg_before.radius + wall_t) if seg_before.direction == "CW" \
             else (seg_before.radius - wall_t)
        dx = c2[0] - c1[0]; dy = c2[1] - c1[1]
        d = math.sqrt(dx*dx + dy*dy)
        return (c1[0] + r1*dx/d, c1[1] + r1*dy/d)

    # Line-Arc or Arc-Line: perpendicular foot from arc center onto offset line
    line_seg = seg_before if b_line else seg_after
    arc_seg  = seg_after  if b_line else seg_before
    center = pts[arc_seg.center]
    S = pts[line_seg.start]; E = pts[line_seg.end]
    D = (E[0]-S[0], E[1]-S[1])
    LN = left_norm(S, E)
    P = off_pt(S, LN, wall_t)
    t = ((center[0]-P[0])*D[0] + (center[1]-P[1])*D[1]) / (D[0]**2 + D[1]**2)
    return (P[0]+t*D[0], P[1]+t*D[1])

for i in range(16):
    seg_b = outline_segs[i]
    seg_a = outline_segs[(i+1) % 16]
    n = int(seg_b.end[1:])   # "O7" -> 7
    pts[f"W{n}"] = _inner_point(seg_b, seg_a)

# --- Inner wall segments (same topology, adjusted radii) ---
inner_segs = [
    LineSeg("W2",  "W1"),
    ArcSeg("W1",  "W0",  "Cf",  R_fillet  - wall_t, "CCW", 20),
    LineSeg("W0",  "W15"),
    ArcSeg("W15", "W14", "C3",  R3i       + wall_t, "CW",  60),
    ArcSeg("W14", "W13", "Cf3", R_f_pox   - wall_t, "CCW", 20),
    LineSeg("W13", "W12"),
    ArcSeg("W12", "W11", "Cf4", R_f_po5   - wall_t, "CCW", 20),
    LineSeg("W11", "W10"),
    ArcSeg("W10", "W9",  "C1",  R1i       + wall_t, "CW",  60),
    LineSeg("W9",  "W8"),
    ArcSeg("W8",  "W7",  "Ct3", R_turn3   - wall_t, "CCW", 20),
    LineSeg("W7",  "W6"),
    ArcSeg("W6",  "W5",  "Ct2", R_turn2   + wall_t, "CW",  20),
    ArcSeg("W5",  "W4",  "Ct1", R_turn1   - wall_t, "CCW", 20),
    LineSeg("W4",  "W3"),
    ArcSeg("W3",  "W2",  "Cf2", R_fillet2 - wall_t, "CCW", 20),
]

outer_poly = path_polygon(outline_segs, pts)
inner_poly = path_polygon(inner_segs, pts)
outer_area = poly_area(outer_poly)
inner_area = poly_area(inner_poly)

# --- SVG ---
out = []
out.append(f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}">')
out.append(f'<rect width="{W}" height="{H}" fill="white"/>')
out.append('<defs>')
out.append('  <marker id="ah" markerWidth="8" markerHeight="6" refX="8" refY="3" orient="auto">'
           '<polygon points="0 0, 8 3, 0 6" fill="#333"/></marker>')
out.append('</defs>')
out.append(f'<text x="{W/2}" y="30" text-anchor="middle" font-family="Arial" font-size="14"'
           f' font-weight="bold">Floorplan &#8212; 8&#8243; Walls</text>')

# --- Interior wall: 6" thick, south face 11'6" north of inner C0-C15 wall ---
int_wall_t = 6.0 / 12.0   # 6 inches
int_wall_south = pts["W0"][1] + 11.5   # 11'6" north of inner face of C0-C15
int_wall_north = int_wall_south + int_wall_t

def _horiz_isects(polygon, n_val):
    """Easting values where polygon boundary crosses a given northing."""
    ints = []
    for i in range(len(polygon)):
        j = (i + 1) % len(polygon)
        n1, n2 = polygon[i][1], polygon[j][1]
        if (n1 <= n_val < n2) or (n2 <= n_val < n1):
            t = (n_val - n1) / (n2 - n1)
            ints.append(polygon[i][0] + t * (polygon[j][0] - polygon[i][0]))
    return ints

def _vert_isects(polygon, e_val):
    """Northing values where polygon boundary crosses a given easting."""
    ints = []
    for i in range(len(polygon)):
        j = (i + 1) % len(polygon)
        e1, e2 = polygon[i][0], polygon[j][0]
        if (e1 <= e_val < e2) or (e2 <= e_val < e1):
            t = (e_val - e1) / (e2 - e1)
            ints.append(polygon[i][1] + t * (polygon[j][1] - polygon[i][1]))
    return ints

_s_ints = _horiz_isects(inner_poly, int_wall_south)
_n_ints = _horiz_isects(inner_poly, int_wall_north)
iw_sw = (min(_s_ints), int_wall_south)
iw_se = (max(_s_ints), int_wall_south)
iw_nw = (min(_n_ints), int_wall_north)
iw_ne = (max(_n_ints), int_wall_north)

# Wall fill: outer gray, inner white cutout
outer_svg = " ".join(f"{to_svg(*p)[0]:.1f},{to_svg(*p)[1]:.1f}" for p in outer_poly)
inner_rev = list(reversed(inner_poly))
inner_svg = " ".join(f"{to_svg(*p)[0]:.1f},{to_svg(*p)[1]:.1f}" for p in inner_rev)
out.append(f'<polygon points="{outer_svg}" fill="rgba(160,160,160,0.35)" stroke="none"/>')
out.append(f'<polygon points="{inner_svg}" fill="white" stroke="none"/>')

# Outer outline strokes
for seg in outline_segs:
    if isinstance(seg, LineSeg):
        sx1, sy1 = to_svg(*pts[seg.start]); sx2, sy2 = to_svg(*pts[seg.end])
        out.append(f'<line x1="{sx1:.1f}" y1="{sy1:.1f}" x2="{sx2:.1f}" y2="{sy2:.1f}"'
                   f' stroke="#333" stroke-width="1.5"/>')
    else:
        poly = segment_polyline(seg, pts)
        svg_p = " ".join(f"{to_svg(e,n)[0]:.1f},{to_svg(e,n)[1]:.1f}" for e,n in poly)
        out.append(f'<polyline points="{svg_p}" fill="none" stroke="#333"'
                   f' stroke-width="1.5" stroke-linecap="round"/>')

# Inner outline strokes
for seg in inner_segs:
    if isinstance(seg, LineSeg):
        sx1, sy1 = to_svg(*pts[seg.start]); sx2, sy2 = to_svg(*pts[seg.end])
        out.append(f'<line x1="{sx1:.1f}" y1="{sy1:.1f}" x2="{sx2:.1f}" y2="{sy2:.1f}"'
                   f' stroke="#666" stroke-width="1.0"/>')
    else:
        poly = segment_polyline(seg, pts)
        svg_p = " ".join(f"{to_svg(e,n)[0]:.1f},{to_svg(e,n)[1]:.1f}" for e,n in poly)
        out.append(f'<polyline points="{svg_p}" fill="none" stroke="#666"'
                   f' stroke-width="1.0" stroke-linecap="round"/>')

# Interior wall IW1
iw_pts = [iw_sw, iw_se, iw_ne, iw_nw]
iw_svg = " ".join(f"{to_svg(*p)[0]:.1f},{to_svg(*p)[1]:.1f}" for p in iw_pts)
out.append(f'<polygon points="{iw_svg}" fill="rgba(160,160,160,0.35)" stroke="none"/>')
for a, b in [(iw_sw, iw_se), (iw_ne, iw_nw)]:
    sx1, sy1 = to_svg(*a); sx2, sy2 = to_svg(*b)
    out.append(f'<line x1="{sx1:.1f}" y1="{sy1:.1f}" x2="{sx2:.1f}" y2="{sy2:.1f}"'
               f' stroke="#666" stroke-width="1.0"/>')
# IW1 label
iw_mid_e = (iw_sw[0] + iw_se[0]) / 2
iw_mid_n = (int_wall_south + int_wall_north) / 2
iw_lx, iw_ly = to_svg(iw_mid_e, iw_mid_n)
out.append(f'<text x="{iw_lx:.1f}" y="{iw_ly+3.5:.1f}" text-anchor="middle" font-family="Arial"'
           f' font-size="8" fill="#666">IW1</text>')

# Interior wall IW2: 6" thick, N-S, west face 6'6" east of inner C1-C2 wall
iw2_w = pts["W1"][0] + 6.5      # west face: 6'6" east of inner C1-C2
iw2_e = iw2_w + int_wall_t      # east face: +6"
iw2_s = int_wall_north           # south end: IW1 north face
iw2_n = pts["W3"][1]             # north end: inner C3-C4 wall
iw2_pts = [(iw2_w, iw2_s), (iw2_e, iw2_s), (iw2_e, iw2_n), (iw2_w, iw2_n)]
iw2_svg = " ".join(f"{to_svg(*p)[0]:.1f},{to_svg(*p)[1]:.1f}" for p in iw2_pts)
out.append(f'<polygon points="{iw2_svg}" fill="rgba(160,160,160,0.35)" stroke="none"/>')
# IW2 stroke lines (west and east faces)
for e_val in [iw2_w, iw2_e]:
    sx1, sy1 = to_svg(e_val, iw2_s)
    sx2, sy2 = to_svg(e_val, iw2_n)
    out.append(f'<line x1="{sx1:.1f}" y1="{sy1:.1f}" x2="{sx2:.1f}" y2="{sy2:.1f}"'
               f' stroke="#666" stroke-width="1.0"/>')
# IW2 label
iw2_mid_e = (iw2_w + iw2_e) / 2
iw2_mid_n = (iw2_s + iw2_n) / 2
iw2_lx, iw2_ly = to_svg(iw2_mid_e, iw2_mid_n)
out.append(f'<text x="{iw2_lx:.1f}" y="{iw2_ly+3.5:.1f}" text-anchor="middle" font-family="Arial"'
           f' font-size="8" fill="#666" transform="rotate(-90,{iw2_lx:.1f},{iw2_ly+3.5:.1f})">IW2</text>')

# Dimension line: IW1 north face to C6-C7 south face (inner), mid-span
dim_e = (pts["O6"][0] + pts["O7"][0]) / 2
dim_top_n = pts["W6"][1]       # south face of C6-C7 wall (inner)
dim_bot_n = int_wall_north     # north face of IW1
dim_dist = dim_top_n - dim_bot_n
dim_ft = int(dim_dist)
dim_in = (dim_dist - dim_ft) * 12
dim_label = f"{dim_ft}' {dim_in:.1f}\""
dx1, dy1 = to_svg(dim_e, dim_bot_n)
dx2, dy2 = to_svg(dim_e, dim_top_n)
tick = 4  # half-width of tick marks in SVG px
out.append(f'<line x1="{dx1:.1f}" y1="{dy1:.1f}" x2="{dx2:.1f}" y2="{dy2:.1f}"'
           f' stroke="#999" stroke-width="0.8"/>')
out.append(f'<line x1="{dx1-tick:.1f}" y1="{dy1:.1f}" x2="{dx1+tick:.1f}" y2="{dy1:.1f}"'
           f' stroke="#999" stroke-width="0.8"/>')
out.append(f'<line x1="{dx2-tick:.1f}" y1="{dy2:.1f}" x2="{dx2+tick:.1f}" y2="{dy2:.1f}"'
           f' stroke="#999" stroke-width="0.8"/>')
dm_x, dm_y = (dx1 + dx2) / 2, (dy1 + dy2) / 2
out.append(f'<text x="{dm_x-6:.1f}" y="{dm_y+3:.1f}" text-anchor="middle" font-family="Arial"'
           f' font-size="8" fill="#999" transform="rotate(-90,{dm_x-6:.1f},{dm_y+3:.1f})">{dim_label}</text>')

# Dimension line: IW2 east face to inside C8-C9 wall, vertically centered in C8-C9 wall
dim2_n = (pts["O8"][1] + pts["O9"][1]) / 2
# West end: IW2 east face
dim2_west_e = iw2_e
# East end: interpolate along inner wall W9-W8 at dim2_n
_w9, _w8 = pts["W9"], pts["W8"]
_t_e = (dim2_n - _w9[1]) / (_w8[1] - _w9[1]) if _w8[1] != _w9[1] else 0.5
dim2_east_e = _w9[0] + _t_e * (_w8[0] - _w9[0])
dim2_dist = dim2_east_e - dim2_west_e
dim2_ft = int(dim2_dist)
dim2_in = (dim2_dist - dim2_ft) * 12
dim2_label = f"{dim2_ft}' {dim2_in:.1f}\""
d2x1, d2y1 = to_svg(dim2_west_e, dim2_n)
d2x2, d2y2 = to_svg(dim2_east_e, dim2_n)
out.append(f'<line x1="{d2x1:.1f}" y1="{d2y1:.1f}" x2="{d2x2:.1f}" y2="{d2y2:.1f}"'
           f' stroke="#999" stroke-width="0.8"/>')
out.append(f'<line x1="{d2x1:.1f}" y1="{d2y1-tick:.1f}" x2="{d2x1:.1f}" y2="{d2y1+tick:.1f}"'
           f' stroke="#999" stroke-width="0.8"/>')
out.append(f'<line x1="{d2x2:.1f}" y1="{d2y2-tick:.1f}" x2="{d2x2:.1f}" y2="{d2y2+tick:.1f}"'
           f' stroke="#999" stroke-width="0.8"/>')
d2m_x, d2m_y = (d2x1 + d2x2) / 2, (d2y1 + d2y2) / 2
out.append(f'<text x="{d2m_x:.1f}" y="{d2m_y-6:.1f}" text-anchor="middle" font-family="Arial"'
           f' font-size="8" fill="#999">{dim2_label}</text>')

# --- Appliances (dimensions from ../hut project) ---
# Washer & Dryer: 35" wide (E) x 30" deep (N), same offsets as ../hut
app_w = 35.0 / 12.0   # width in feet (easting)
app_d = 30.0 / 12.0   # depth in feet (northing)
app_offset_e = 6.0 / 12.0   # 6" from interior west wall
app_offset_n = 4.0 / 12.0   # 4" from interior south wall
app_gap = 1.0 / 12.0        # 1" gap between appliances

# Dryer in C0-C1 corner: offset from inner west wall (W1) and inner south wall (W0)
dryer_w = pts["W1"][0] + app_offset_e
dryer_s = pts["W0"][1] + app_offset_n
dryer_e = dryer_w + app_w
dryer_n = dryer_s + app_d

# Washer north of dryer with 1" gap
washer_w = dryer_w
washer_s = dryer_n + app_gap
washer_e = dryer_e
washer_n = washer_s + app_d

for label, sw_e, sw_n, ne_e, ne_n in [
    ("DRYER",  dryer_w,  dryer_s,  dryer_e,  dryer_n),
    ("WASHER", washer_w, washer_s, washer_e, washer_n),
]:
    sx1, sy1 = to_svg(sw_e, ne_n)  # SVG top-left = survey NW corner
    sx2, sy2 = to_svg(ne_e, sw_n)  # SVG bottom-right = survey SE corner
    sw = sx2 - sx1; sh = sy2 - sy1
    out.append(f'<rect x="{sx1:.1f}" y="{sy1:.1f}" width="{sw:.1f}" height="{sh:.1f}"'
               f' fill="rgba(100,150,200,0.2)" stroke="#4682B4" stroke-width="0.8"/>')
    cx, cy = (sx1 + sx2) / 2, (sy1 + sy2) / 2
    out.append(f'<text x="{cx:.1f}" y="{cy+3:.1f}" text-anchor="middle" font-family="Arial"'
               f' font-size="7" fill="#4682B4">{label}</text>')

# Counter from ../hut Utility room: 24" deep x 72" long, 9" NW corner radius
ctr_depth = 24.0 / 12.0    # 2' E-W
ctr_length = 72.0 / 12.0   # 6' N-S
ctr_nw_r = 9.0 / 12.0      # 0.75' NW corner radius
ctr_gap = 36.0 / 12.0      # 36" east of W/D (same offset as ../hut)

ctr_w = dryer_e + ctr_gap
ctr_e = ctr_w + ctr_depth
ctr_s = pts["W0"][1]        # south edge at inner south wall (same as ../hut)
ctr_n = ctr_s + ctr_length

# SVG path with rounded NW corner
_csw = to_svg(ctr_w, ctr_s)
_cse = to_svg(ctr_e, ctr_s)
_cne = to_svg(ctr_e, ctr_n)
_cnas = to_svg(ctr_w + ctr_nw_r, ctr_n)     # NW arc start (north edge)
_cnae = to_svg(ctr_w, ctr_n - ctr_nw_r)     # NW arc end (west edge)
_r_svg = abs(_cnas[0] - to_svg(ctr_w, ctr_n)[0])
ctr_path = (f'M {_csw[0]:.1f},{_csw[1]:.1f} '
            f'L {_cse[0]:.1f},{_cse[1]:.1f} '
            f'L {_cne[0]:.1f},{_cne[1]:.1f} '
            f'L {_cnas[0]:.1f},{_cnas[1]:.1f} '
            f'A {_r_svg:.1f} {_r_svg:.1f} 0 0 0 {_cnae[0]:.1f},{_cnae[1]:.1f} '
            f'Z')
out.append(f'<path d="{ctr_path}" fill="rgba(100,150,200,0.2)" stroke="#4682B4" stroke-width="0.8"/>')
_ccx = (_csw[0] + _cse[0]) / 2
_ccy = (_csw[1] + _cne[1]) / 2
out.append(f'<text x="{_ccx:.1f}" y="{_ccy:.1f}" text-anchor="middle" font-family="Arial"'
           f' font-size="7" fill="#4682B4" letter-spacing="0.5" transform="rotate(-90,{_ccx:.1f},{_ccy:.1f})">COUNTER</text>')

# Bedroom and closet walls from ../hut plan (placed relative to counter)
iw_thick_3 = 3.0 / 12.0    # 3" Wall 8 thickness
iw_thick_4 = 4.0 / 12.0    # 4" Wall 1 South thickness
closet2_width = 30.0 / 12.0 # 30" closet width

# Wall 8 L-shape (west/north walls of closet, east of counter)
#   Vertical section: ctr_e to ctr_e+3", from south wall to top of horizontal
#   Horizontal section: at counter north, connecting vertical to Wall 1 South
w8_poly = [
    (ctr_e, ctr_s),                                                  # SW
    (ctr_e + iw_thick_3, ctr_s),                                     # SE of vertical
    (ctr_e + iw_thick_3, ctr_n),                                     # inner corner
    (ctr_e + iw_thick_3 + closet2_width, ctr_n),                    # SE of horizontal
    (ctr_e + iw_thick_3 + closet2_width, ctr_n + iw_thick_3),       # NE of horizontal
    (ctr_e, ctr_n + iw_thick_3),                                     # NW
]
w8_svg = " ".join(f"{to_svg(*p)[0]:.1f},{to_svg(*p)[1]:.1f}" for p in w8_poly)
out.append(f'<polygon points="{w8_svg}" fill="rgba(160,160,160,0.35)" stroke="#666" stroke-width="0.8"/>')

# Wall 1 South (partition between closet 2 and bedroom, 4" thick)
w1s_w = ctr_e + iw_thick_3 + closet2_width
w1s_e = w1s_w + iw_thick_4
w1s_s = ctr_s
w1s_n = int_wall_south   # extends north to IW1 south face
w1s_poly = [(w1s_w, w1s_s), (w1s_e, w1s_s), (w1s_e, w1s_n), (w1s_w, w1s_n)]
w1s_svg = " ".join(f"{to_svg(*p)[0]:.1f},{to_svg(*p)[1]:.1f}" for p in w1s_poly)
out.append(f'<polygon points="{w1s_svg}" fill="rgba(160,160,160,0.35)" stroke="#666" stroke-width="0.8"/>')

# Wall 6 (bedroom east wall, 4" thick) — 11'8" east of Wall 1 South east face
bedroom_width = 140.0 / 12.0   # 11'8"
w6_w = w1s_e + bedroom_width
w6_e = w6_w + iw_thick_4
wall_south_n = -4.0 / 12.0  # south end of bedroom/closet walls: -4"
w6_south_w = wall_south_n
w6_south_e = wall_south_n
w6_poly = [(w6_w, w6_south_w), (w6_e, w6_south_e), (w6_e, int_wall_south), (w6_w, int_wall_south)]
w6_svg = " ".join(f"{to_svg(*p)[0]:.1f},{to_svg(*p)[1]:.1f}" for p in w6_poly)
out.append(f'<polygon points="{w6_svg}" fill="rgba(160,160,160,0.35)" stroke="#666" stroke-width="0.8"/>')

# Wall 5 (L-shaped, 3" thick — east/north walls of closet 1)
closet1_width = 30.0 / 12.0    # 30" closet
closet1_top = ctr_n - 1.0      # 12" shorter than closet 2 (matching ../hut)
w5_w = w6_e + closet1_width
w5_e = w5_w + iw_thick_3
w5_south_w = wall_south_n
w5_south_e = wall_south_n
# L-shape: horizontal bar at closet1_top, vertical bar on the east
w5_poly = [
    (w6_e, closet1_top + iw_thick_3),   # NW (top-left of horizontal)
    (w5_e, closet1_top + iw_thick_3),   # NE (top-right of vertical)
    (w5_e, w5_south_e),                  # SE (bottom of east face)
    (w5_w, w5_south_w),                  # SW of vertical (bottom of west face)
    (w5_w, closet1_top),                 # inner corner
    (w6_e, closet1_top),                 # SW of horizontal
]
w5_svg = " ".join(f"{to_svg(*p)[0]:.1f},{to_svg(*p)[1]:.1f}" for p in w5_poly)
out.append(f'<polygon points="{w5_svg}" fill="rgba(160,160,160,0.35)" stroke="#666" stroke-width="0.8"/>')

# Room labels
_cl_cx = (ctr_e + iw_thick_3 + w1s_w) / 2
_cl_cy = (ctr_s + ctr_n) / 2
_clx, _cly = to_svg(_cl_cx, _cl_cy)
out.append(f'<text x="{_clx:.1f}" y="{_cly+3:.1f}" text-anchor="middle" font-family="Arial"'
           f' font-size="7" fill="#666" transform="rotate(-90,{_clx:.1f},{_cly+3:.1f})">CLOSET</text>')

# King Bed (from ../hut: 76" wide x 94" long incl. frame, 2" from south wall, centered E-W)
bed_w_dim = 76.0 / 12.0    # E-W width
bed_l_dim = 94.0 / 12.0    # N-S length (incl. headboard/frame)
bed_offset_n = 2.0 / 12.0  # 2" from inner south wall
bed_cx = (w1s_e + w6_w) / 2
bed_w = bed_cx - bed_w_dim / 2
bed_e = bed_cx + bed_w_dim / 2
bed_s = ctr_s + bed_offset_n
bed_n = bed_s + bed_l_dim
_bed_sw = to_svg(bed_w, bed_n)   # SVG top-left = survey NW
_bed_se = to_svg(bed_e, bed_s)   # SVG bottom-right = survey SE
_bed_sw_x, _bed_sw_y = _bed_sw
_bed_se_x, _bed_se_y = _bed_se
_bed_w = _bed_se_x - _bed_sw_x
_bed_h = _bed_se_y - _bed_sw_y
out.append(f'<rect x="{_bed_sw_x:.1f}" y="{_bed_sw_y:.1f}" width="{_bed_w:.1f}" height="{_bed_h:.1f}"'
           f' fill="rgba(100,150,200,0.2)" stroke="#4682B4" stroke-width="0.8"/>')
_bed_cx_svg = (_bed_sw_x + _bed_se_x) / 2
# Label y: 76.5% down the bed rectangle
_bed_label_y = _bed_sw_y + 0.765 * _bed_h
out.append(f'<text x="{_bed_cx_svg:.1f}" y="{_bed_label_y+3:.1f}" text-anchor="middle" font-family="Arial"'
           f' font-size="7" fill="#4682B4">KING BED</text>')

_bd_cx = (w1s_e + w6_w) / 2   # centered in 11'8" bedroom width
_bd_cy = (ctr_s + int_wall_south) / 2
_bdx, _bdy = to_svg(_bd_cx, _bd_cy)
out.append(f'<text x="{_bdx:.1f}" y="{_bdy+3:.1f}" text-anchor="middle" font-family="Arial"'
           f' font-size="8" fill="#666">BEDROOM</text>')

# Bedroom interior dimension lines
# E-W dimension (horizontal): w1s_e to w6_w, placed at 25% from south wall
bd_ew_dist = w6_w - w1s_e
bd_ew_ft = int(bd_ew_dist)
bd_ew_in = (bd_ew_dist - bd_ew_ft) * 12
bd_ew_label = f"{bd_ew_ft}' {bd_ew_in:.0f}\""
bd_ew_n = ctr_s + 0.25 * (int_wall_south - ctr_s)
bd_ew_x1, bd_ew_y1 = to_svg(w1s_e, bd_ew_n)
bd_ew_x2, bd_ew_y2 = to_svg(w6_w, bd_ew_n)
out.append(f'<line x1="{bd_ew_x1:.1f}" y1="{bd_ew_y1:.1f}" x2="{bd_ew_x2:.1f}" y2="{bd_ew_y2:.1f}"'
           f' stroke="#999" stroke-width="0.8"/>')
out.append(f'<line x1="{bd_ew_x1:.1f}" y1="{bd_ew_y1-tick:.1f}" x2="{bd_ew_x1:.1f}" y2="{bd_ew_y1+tick:.1f}"'
           f' stroke="#999" stroke-width="0.8"/>')
out.append(f'<line x1="{bd_ew_x2:.1f}" y1="{bd_ew_y2-tick:.1f}" x2="{bd_ew_x2:.1f}" y2="{bd_ew_y2+tick:.1f}"'
           f' stroke="#999" stroke-width="0.8"/>')
bd_ew_mx = (bd_ew_x1 + bd_ew_x2) / 2
out.append(f'<text x="{bd_ew_mx:.1f}" y="{bd_ew_y1-6:.1f}" text-anchor="middle" font-family="Arial"'
           f' font-size="8" fill="#999">{bd_ew_label}</text>')

# N-S dimension (vertical): ctr_s to int_wall_south, centered between C15 and Wall 1 South west face
bd_ns_dist = int_wall_south - ctr_s
bd_ns_ft = int(bd_ns_dist)
bd_ns_in = (bd_ns_dist - bd_ns_ft) * 12
bd_ns_label = f"{bd_ns_ft}' {bd_ns_in:.0f}\""
bd_ns_e = (pts["O15"][0] + w1s_w) / 2
bd_ns_x1, bd_ns_y1 = to_svg(bd_ns_e, ctr_s)
bd_ns_x2, bd_ns_y2 = to_svg(bd_ns_e, int_wall_south)
out.append(f'<line x1="{bd_ns_x1:.1f}" y1="{bd_ns_y1:.1f}" x2="{bd_ns_x2:.1f}" y2="{bd_ns_y2:.1f}"'
           f' stroke="#999" stroke-width="0.8"/>')
out.append(f'<line x1="{bd_ns_x1-tick:.1f}" y1="{bd_ns_y1:.1f}" x2="{bd_ns_x1+tick:.1f}" y2="{bd_ns_y1:.1f}"'
           f' stroke="#999" stroke-width="0.8"/>')
out.append(f'<line x1="{bd_ns_x2-tick:.1f}" y1="{bd_ns_y2:.1f}" x2="{bd_ns_x2+tick:.1f}" y2="{bd_ns_y2:.1f}"'
           f' stroke="#999" stroke-width="0.8"/>')
bd_ns_mx, bd_ns_my = bd_ns_x1, (bd_ns_y1 + bd_ns_y2) / 2
out.append(f'<text x="{bd_ns_mx-6:.1f}" y="{bd_ns_my+3:.1f}" text-anchor="middle" font-family="Arial"'
           f' font-size="8" fill="#999" transform="rotate(-90,{bd_ns_mx-6:.1f},{bd_ns_my+3:.1f})">{bd_ns_label}</text>')

# East closet label (closet 1, between Wall 6 and Wall 5)
_cl1_cx = (w6_e + w5_w) / 2
_cl1_cy = (ctr_s + closet1_top) / 2
_cl1x, _cl1y = to_svg(_cl1_cx, _cl1_cy)
out.append(f'<text x="{_cl1x:.1f}" y="{_cl1y+3:.1f}" text-anchor="middle" font-family="Arial"'
           f' font-size="7" fill="#666" transform="rotate(-90,{_cl1x:.1f},{_cl1y+3:.1f})">CLOSET</text>')

# C-point labels (outer vertices displayed as C0-C15)
vs_map = _ns["outline_cfg"].vertex_styles
for i in range(16):
    o_name = f"O{i}"
    sx, sy = to_svg(*pts[o_name])
    out.append(f'<circle cx="{sx:.1f}" cy="{sy:.1f}" r="2.5" fill="#333"/>')
    if o_name in vs_map:
        vs = vs_map[o_name]
        out.append(f'<text x="{sx+vs.dx:.1f}" y="{sy+vs.dy:.1f}" text-anchor="{vs.anchor}"'
                   f' font-family="Arial" font-size="9" font-weight="bold"'
                   f' fill="#333">C{i}</text>')

# North arrow
out.append('<line x1="742" y1="560" x2="742" y2="524" stroke="#333" stroke-width="2"'
           ' marker-end="url(#ah)"/>')
out.append('<text x="742" y="518" text-anchor="middle" font-family="Arial"'
           ' font-size="13" font-weight="bold">N</text>')

# Title block (right edge aligned with N arrow, bottom with C4 label)
_c4_sx, _c4_sy = to_svg(*pts["O4"])
_c4_vs = vs_map["O4"]
tb_right = 752      # right edge, aligned with N arrow center + margin
tb_bottom = _c4_sy + _c4_vs.dy + 4  # bottom aligned with C4 label
tb_w = 130
tb_h = 58
tb_left = tb_right - tb_w
tb_top = tb_bottom - tb_h
tb_cx = (tb_left + tb_right) / 2
out.append(f'<rect x="{tb_left:.1f}" y="{tb_top:.1f}" width="{tb_w}" height="{tb_h}"'
           f' fill="white" stroke="#333" stroke-width="1"/>')
out.append(f'<text x="{tb_cx:.1f}" y="{tb_top+14:.1f}" text-anchor="middle"'
           f' font-family="Arial" font-size="11" font-weight="bold" fill="#333">'
           f'{inner_area:.2f} sq ft</text>')
out.append(f'<text x="{tb_cx:.1f}" y="{tb_top+26:.1f}" text-anchor="middle"'
           f' font-family="Arial" font-size="8" fill="#666">Interior area</text>')
out.append(f'<text x="{tb_cx:.1f}" y="{tb_top+40:.1f}" text-anchor="middle"'
           f' font-family="Arial" font-size="11" font-weight="bold" fill="#333">'
           f'{outer_area:.2f} sq ft</text>')
out.append(f'<text x="{tb_cx:.1f}" y="{tb_top+52:.1f}" text-anchor="middle"'
           f' font-family="Arial" font-size="8" fill="#666">Exterior area</text>')

# Footer
out.append(f'<text x="{W/2}" y="{H-2}" text-anchor="middle" font-family="Arial" font-size="7.5"'
           f' fill="#999">Wall thickness: 8&#8243; &#8226;'
           f' Outer: {outer_area:.2f} sq ft &#8226; Interior: {inner_area:.2f} sq ft</text>')
out.append('</svg>')

svg_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "floorplan.svg")
with open(svg_path, "w") as f:
    f.write("\n".join(out))

print(f"Floorplan written to {svg_path}")
print(f"Outer area:    {outer_area:.2f} sq ft")
print(f"Interior area: {inner_area:.2f} sq ft")
print(f"Wall area:     {outer_area - inner_area:.2f} sq ft")
print()
for i in range(16):
    o = pts[f"O{i}"]; w = pts[f"W{i}"]
    print(f"  C{i:<2d} ({o[0]:8.4f}, {o[1]:8.4f})  ->  inner ({w[0]:8.4f}, {w[1]:8.4f})")
