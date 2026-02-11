"""Generate floorplan SVG with 8" wall inset from the outline path.

Runs gen_path_svg.py to obtain outline geometry, then computes an 8" inset.
Points C0-C15 (plus C13a, C13b) correspond to outline points O0-O15 (plus O13a, O13b).

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
R_wall    = _ns["R_wall"]       # Cw3 wall arc (10")
R_w1      = _ns["R_w1"]        # Cw1 wall arc (28")
R_w2      = _ns["R_w2"]        # Cw2 wall arc (28")
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

for i in range(18):
    seg_b = outline_segs[i]
    seg_a = outline_segs[(i+1) % 18]
    w_name = "W" + seg_b.end[1:]   # "O7" -> "W7", "O13b" -> "W13b"
    pts[w_name] = _inner_point(seg_b, seg_a)

# --- Inner wall segments (same topology, adjusted radii) ---
inner_segs = [
    LineSeg("W2",  "W1"),
    ArcSeg("W1",  "W0",  "Cf",  R_fillet  - wall_t, "CCW", 20),
    LineSeg("W0",  "W15"),
    ArcSeg("W15", "W14", "Cw1", R_w1      + wall_t, "CW",  60),
    ArcSeg("W14", "W13b","Cw2", R_w2      - wall_t, "CCW", 60),
    LineSeg("W13b","W13a"),
    ArcSeg("W13a","W13", "Cw3", R_wall    - wall_t, "CCW", 20),
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
out.append(f'<text x="{iw_lx:.1f}" y="{iw_ly+3.5-8:.1f}" text-anchor="middle" font-family="Arial"'
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
out.append(f'<text x="{iw2_lx-4:.1f}" y="{iw2_ly+3.5:.1f}" text-anchor="middle" font-family="Arial"'
           f' font-size="8" fill="#666" transform="rotate(-90,{iw2_lx-4:.1f},{iw2_ly+3.5:.1f})">IW2</text>')

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
out.append(f'<text x="{dm_x-3:.1f}" y="{dm_y+3:.1f}" text-anchor="middle" font-family="Arial"'
           f' font-size="8" fill="#999" transform="rotate(-90,{dm_x-3:.1f},{dm_y+3:.1f})">{dim_label}</text>')

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
out.append(f'<text x="{d2m_x:.1f}" y="{d2m_y-3:.1f}" text-anchor="middle" font-family="Arial"'
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
iw_thick_4 = 4.0 / 12.0    # 4" IW3 thickness
closet2_width = 30.0 / 12.0 # 30" closet width

# Wall 8 L-shape (west/north walls of closet, east of counter)
#   Vertical section: ctr_e to ctr_e+3", from south wall to top of horizontal
#   Horizontal section: at counter north, connecting vertical to IW3
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

# IW3 (west bedroom wall, 4" thick)
iw3_w = ctr_e + iw_thick_3 + closet2_width
iw3_e = iw3_w + iw_thick_4
iw3_s = ctr_s
iw3_n = int_wall_south   # extends north to IW1 south face
iw3_poly = [(iw3_w, iw3_s), (iw3_e, iw3_s), (iw3_e, iw3_n), (iw3_w, iw3_n)]
iw3_svg = " ".join(f"{to_svg(*p)[0]:.1f},{to_svg(*p)[1]:.1f}" for p in iw3_poly)
out.append(f'<polygon points="{iw3_svg}" fill="rgba(160,160,160,0.35)" stroke="#666" stroke-width="0.8"/>')
# IW3 label
iw3_mid_e = (iw3_w + iw3_e) / 2
iw3_mid_n = (iw3_s + iw3_n) / 2
iw3_lx, iw3_ly = to_svg(iw3_mid_e, iw3_mid_n)
out.append(f'<text x="{iw3_lx-4:.1f}" y="{iw3_ly+3.5:.1f}" text-anchor="middle" font-family="Arial"'
           f' font-size="8" fill="#666" transform="rotate(-90,{iw3_lx-4:.1f},{iw3_ly+3.5:.1f})">IW3</text>')

# IW4 (bedroom east wall, 4" thick) — 11'8" east of IW3 east face
bedroom_width = 140.0 / 12.0   # 11'8"
iw4_w = iw3_e + bedroom_width
iw4_e = iw4_w + iw_thick_4
wall_south_n = -4.0 / 12.0  # south end of bedroom/closet walls: -4"
iw4_south_w = wall_south_n
iw4_south_e = wall_south_n
iw4_poly = [(iw4_w, iw4_south_w), (iw4_e, iw4_south_e), (iw4_e, int_wall_south), (iw4_w, int_wall_south)]
iw4_svg = " ".join(f"{to_svg(*p)[0]:.1f},{to_svg(*p)[1]:.1f}" for p in iw4_poly)
out.append(f'<polygon points="{iw4_svg}" fill="rgba(160,160,160,0.35)" stroke="#666" stroke-width="0.8"/>')
# IW4 label
iw4_mid_e = (iw4_w + iw4_e) / 2
iw4_mid_n = (iw4_south_w + int_wall_south) / 2
iw4_lx, iw4_ly = to_svg(iw4_mid_e, iw4_mid_n)
out.append(f'<text x="{iw4_lx-4:.1f}" y="{iw4_ly+3.5:.1f}" text-anchor="middle" font-family="Arial"'
           f' font-size="8" fill="#666" transform="rotate(-90,{iw4_lx-4:.1f},{iw4_ly+3.5:.1f})">IW4</text>')

# Wall 5 (L-shaped, 3" thick — east/north walls of closet 1)
closet1_width = 30.0 / 12.0    # 30" closet
closet1_top = wall_south_n + 6.0  # 6'0" inside N-S
w5_w = iw4_e + closet1_width
w5_e = w5_w + iw_thick_3
w5_south_w = wall_south_n
w5_south_e = wall_south_n
# L-shape: horizontal bar at closet1_top, vertical bar on the east
w5_poly = [
    (iw4_e, closet1_top + iw_thick_3),   # NW (top-left of horizontal)
    (w5_e, closet1_top + iw_thick_3),   # NE (top-right of vertical)
    (w5_e, w5_south_e),                  # SE (bottom of east face)
    (w5_w, w5_south_w),                  # SW of vertical (bottom of west face)
    (w5_w, closet1_top),                 # inner corner
    (iw4_e, closet1_top),                 # SW of horizontal
]
w5_svg = " ".join(f"{to_svg(*p)[0]:.1f},{to_svg(*p)[1]:.1f}" for p in w5_poly)
out.append(f'<polygon points="{w5_svg}" fill="rgba(160,160,160,0.35)" stroke="#666" stroke-width="0.8"/>')

# King Bed (from ../hut: 76" wide x 94" long incl. frame, 2" from south wall, centered E-W)
bed_w_dim = 76.0 / 12.0    # E-W width
bed_l_dim = 94.0 / 12.0    # N-S length (incl. headboard/frame)
bed_offset_n = 2.0 / 12.0  # 2" from inner south wall
bed_cx = (iw3_e + iw4_w) / 2
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

_bd_cx = (iw3_e + iw4_w) / 2   # centered in 11'8" bedroom width
_bd_cy = (ctr_s + int_wall_south) / 2
_bdx, _bdy = to_svg(_bd_cx, _bd_cy)
out.append(f'<text x="{_bdx:.1f}" y="{_bdy+3:.1f}" text-anchor="middle" font-family="Arial"'
           f' font-size="8" fill="#666">BEDROOM</text>')

# Bedroom interior dimension lines
# E-W dimension (horizontal): iw3_e to iw4_w, placed at 25% from south wall
bd_ew_dist = iw4_w - iw3_e
bd_ew_ft = int(bd_ew_dist)
bd_ew_in = (bd_ew_dist - bd_ew_ft) * 12
bd_ew_label = f"{bd_ew_ft}' {bd_ew_in:.0f}\""
bd_ew_n = ctr_s + 0.25 * (int_wall_south - ctr_s)
bd_ew_x1, bd_ew_y1 = to_svg(iw3_e, bd_ew_n)
bd_ew_x2, bd_ew_y2 = to_svg(iw4_w, bd_ew_n)
out.append(f'<line x1="{bd_ew_x1:.1f}" y1="{bd_ew_y1:.1f}" x2="{bd_ew_x2:.1f}" y2="{bd_ew_y2:.1f}"'
           f' stroke="#999" stroke-width="0.8"/>')
out.append(f'<line x1="{bd_ew_x1:.1f}" y1="{bd_ew_y1-tick:.1f}" x2="{bd_ew_x1:.1f}" y2="{bd_ew_y1+tick:.1f}"'
           f' stroke="#999" stroke-width="0.8"/>')
out.append(f'<line x1="{bd_ew_x2:.1f}" y1="{bd_ew_y2-tick:.1f}" x2="{bd_ew_x2:.1f}" y2="{bd_ew_y2+tick:.1f}"'
           f' stroke="#999" stroke-width="0.8"/>')
bd_ew_mx = (bd_ew_x1 + bd_ew_x2) / 2
out.append(f'<text x="{bd_ew_mx:.1f}" y="{bd_ew_y1-3:.1f}" text-anchor="middle" font-family="Arial"'
           f' font-size="8" fill="#999">{bd_ew_label}</text>')

# N-S dimension (vertical): ctr_s to int_wall_south, 2' east of IW3 east face
bd_ns_dist = int_wall_south - ctr_s
bd_ns_ft = int(bd_ns_dist)
bd_ns_in = (bd_ns_dist - bd_ns_ft) * 12
bd_ns_label = f"{bd_ns_ft}' {bd_ns_in:.0f}\""
bd_ns_e = iw3_e + 2.0
bd_ns_x1, bd_ns_y1 = to_svg(bd_ns_e, ctr_s)
bd_ns_x2, bd_ns_y2 = to_svg(bd_ns_e, int_wall_south)
out.append(f'<line x1="{bd_ns_x1:.1f}" y1="{bd_ns_y1:.1f}" x2="{bd_ns_x2:.1f}" y2="{bd_ns_y2:.1f}"'
           f' stroke="#999" stroke-width="0.8"/>')
out.append(f'<line x1="{bd_ns_x1-tick:.1f}" y1="{bd_ns_y1:.1f}" x2="{bd_ns_x1+tick:.1f}" y2="{bd_ns_y1:.1f}"'
           f' stroke="#999" stroke-width="0.8"/>')
out.append(f'<line x1="{bd_ns_x2-tick:.1f}" y1="{bd_ns_y2:.1f}" x2="{bd_ns_x2+tick:.1f}" y2="{bd_ns_y2:.1f}"'
           f' stroke="#999" stroke-width="0.8"/>')
bd_ns_mx, bd_ns_my = bd_ns_x1, (bd_ns_y1 + bd_ns_y2) / 2
out.append(f'<text x="{bd_ns_mx-3:.1f}" y="{bd_ns_my+3:.1f}" text-anchor="middle" font-family="Arial"'
           f' font-size="8" fill="#999" transform="rotate(-90,{bd_ns_mx-3:.1f},{bd_ns_my+3:.1f})">{bd_ns_label}</text>')

# N-S dimension line — west closet (closet 2): ctr_s to ctr_n
cl2_ns_dist = ctr_n - ctr_s
cl2_ns_ft = int(cl2_ns_dist)
cl2_ns_in = (cl2_ns_dist - cl2_ns_ft) * 12
cl2_ns_label = f"CLOSET {cl2_ns_ft}' {cl2_ns_in:.0f}\""
cl2_ns_e = (ctr_e + iw_thick_3 + iw3_w) / 2
cl2_x1, cl2_y1 = to_svg(cl2_ns_e, ctr_s)
cl2_x2, cl2_y2 = to_svg(cl2_ns_e, ctr_n)
out.append(f'<line x1="{cl2_x1:.1f}" y1="{cl2_y1:.1f}" x2="{cl2_x2:.1f}" y2="{cl2_y2:.1f}"'
           f' stroke="#999" stroke-width="0.8"/>')
out.append(f'<line x1="{cl2_x1-tick:.1f}" y1="{cl2_y1:.1f}" x2="{cl2_x1+tick:.1f}" y2="{cl2_y1:.1f}"'
           f' stroke="#999" stroke-width="0.8"/>')
out.append(f'<line x1="{cl2_x2-tick:.1f}" y1="{cl2_y2:.1f}" x2="{cl2_x2+tick:.1f}" y2="{cl2_y2:.1f}"'
           f' stroke="#999" stroke-width="0.8"/>')
cl2_mx, cl2_my = cl2_x1, (cl2_y1 + cl2_y2) / 2
out.append(f'<text x="{cl2_mx-3:.1f}" y="{cl2_my+3:.1f}" text-anchor="middle" font-family="Arial"'
           f' font-size="8" fill="#999" transform="rotate(-90,{cl2_mx-3:.1f},{cl2_my+3:.1f})">{cl2_ns_label}</text>')

# N-S dimension line — east closet (closet 1): wall_south_n to closet1_top
cl1_ns_dist = closet1_top - wall_south_n
cl1_ns_ft = int(cl1_ns_dist)
cl1_ns_in = (cl1_ns_dist - cl1_ns_ft) * 12
cl1_ns_label = f"CLOSET {cl1_ns_ft}' {cl1_ns_in:.0f}\""
cl1_ns_e = (iw4_e + w5_w) / 2
cl1_x1, cl1_y1 = to_svg(cl1_ns_e, wall_south_n)
cl1_x2, cl1_y2 = to_svg(cl1_ns_e, closet1_top)
out.append(f'<line x1="{cl1_x1:.1f}" y1="{cl1_y1:.1f}" x2="{cl1_x2:.1f}" y2="{cl1_y2:.1f}"'
           f' stroke="#999" stroke-width="0.8"/>')
out.append(f'<line x1="{cl1_x1-tick:.1f}" y1="{cl1_y1:.1f}" x2="{cl1_x1+tick:.1f}" y2="{cl1_y1:.1f}"'
           f' stroke="#999" stroke-width="0.8"/>')
out.append(f'<line x1="{cl1_x2-tick:.1f}" y1="{cl1_y2:.1f}" x2="{cl1_x2+tick:.1f}" y2="{cl1_y2:.1f}"'
           f' stroke="#999" stroke-width="0.8"/>')
cl1_mx, cl1_my = cl1_x1, (cl1_y1 + cl1_y2) / 2
out.append(f'<text x="{cl1_mx-3:.1f}" y="{cl1_my+3:.1f}" text-anchor="middle" font-family="Arial"'
           f' font-size="8" fill="#999" transform="rotate(-90,{cl1_mx-3:.1f},{cl1_my+3:.1f})">{cl1_ns_label}</text>')

# C-point labels (outer vertices displayed as C0-C15, C13a, C13b)
vs_map = _ns["outline_cfg"].vertex_styles
_o_names = []
for seg in outline_segs:
    if seg.start not in _o_names:
        _o_names.append(seg.start)
for o_name in _o_names:
    sx, sy = to_svg(*pts[o_name])
    out.append(f'<circle cx="{sx:.1f}" cy="{sy:.1f}" r="2.5" fill="#333"/>')
    if o_name in vs_map:
        vs = vs_map[o_name]
        c_label = "C" + o_name[1:]  # O13b -> C13b
        out.append(f'<text x="{sx+vs.dx:.1f}" y="{sy+vs.dy:.1f}" text-anchor="{vs.anchor}"'
                   f' font-family="Arial" font-size="9" font-weight="bold"'
                   f' fill="#333">{c_label}</text>')

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
_print_names = []
for seg in outline_segs:
    if seg.start not in _print_names:
        _print_names.append(seg.start)
for o_name in _print_names:
    w_name = "W" + o_name[1:]
    c_name = o_name.replace("O", "C")
    o = pts[o_name]; w = pts[w_name]
    print(f"  {c_name:<5s} ({o[0]:8.4f}, {o[1]:8.4f})  ->  inner ({w[0]:8.4f}, {w[1]:8.4f})")
