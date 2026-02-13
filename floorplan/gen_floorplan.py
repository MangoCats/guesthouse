"""Generate floorplan SVG with 8" wall inset from the outline path.

Computes geometry from shared/ and floorplan/ packages.
Outline points F0-F21, inner wall points W0-W21.
"""
import os, math

from shared.types import LineSeg, ArcSeg
from shared.geometry import (
    segment_polyline, path_polygon, poly_area,
    compute_inner_walls, horiz_isects, fmt_dist,
)
from shared.survey import compute_traverse, compute_three_arc, compute_inset
from shared.svg import make_svg_transform, W, H
from floorplan.geometry import compute_outline_geometry, OutlineAnchors
from floorplan.constants import WALL_OUTER

# --- Compute all geometry ---
pts, _p3_trav = compute_traverse()
to_svg = make_svg_transform(_p3_trav)
_arc_info = compute_three_arc(pts)
_inset = compute_inset(pts, _arc_info["R1"], _arc_info["R2"], _arc_info["R3"],
                       _arc_info["nE"], _arc_info["nN"])
pts.update(_inset.pts_update)
_anchors = OutlineAnchors(
    Pi2=pts["Pi2"], Pi3=pts["Pi3"], Ti3=pts["Ti3"],
    PiX=pts["PiX"], Pi5=pts["Pi5"],
    TC1=pts["TC1"], R1i=_inset.R1i,
)
_outline_geo = compute_outline_geometry(_anchors)
pts.update(_outline_geo.fp_pts)
outline_segs = _outline_geo.outline_segs
_radii = _outline_geo.radii

# --- SVG Helpers ---
def dim_line_h(out, e1, n, e2, label):
    """Horizontal (E-W) dimension line with vertical tick marks."""
    x1, y1 = to_svg(e1, n); x2, y2 = to_svg(e2, n)
    _t = 4
    out.append(f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" stroke="#999" stroke-width="0.8"/>')
    out.append(f'<line x1="{x1:.1f}" y1="{y1-_t:.1f}" x2="{x1:.1f}" y2="{y1+_t:.1f}" stroke="#999" stroke-width="0.8"/>')
    out.append(f'<line x1="{x2:.1f}" y1="{y2-_t:.1f}" x2="{x2:.1f}" y2="{y2+_t:.1f}" stroke="#999" stroke-width="0.8"/>')
    out.append(f'<text x="{(x1+x2)/2:.1f}" y="{y1-3:.1f}" text-anchor="middle" font-family="Arial" font-size="8" fill="#999">{label}</text>')

def dim_line_v(out, e, n1, n2, label):
    """Vertical (N-S) dimension line with horizontal tick marks."""
    x1, y1 = to_svg(e, n1); x2, y2 = to_svg(e, n2)
    _t = 4
    out.append(f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" stroke="#999" stroke-width="0.8"/>')
    out.append(f'<line x1="{x1-_t:.1f}" y1="{y1:.1f}" x2="{x1+_t:.1f}" y2="{y1:.1f}" stroke="#999" stroke-width="0.8"/>')
    out.append(f'<line x1="{x2-_t:.1f}" y1="{y2:.1f}" x2="{x2+_t:.1f}" y2="{y2:.1f}" stroke="#999" stroke-width="0.8"/>')
    lx, ly = x1 - 3, (y1 + y2) / 2 + 3
    out.append(f'<text x="{lx:.1f}" y="{ly:.1f}" text-anchor="middle" font-family="Arial" font-size="8" fill="#999" transform="rotate(-90,{lx:.1f},{ly:.1f})">{label}</text>')

def wall_poly(out, points, stroke=True):
    """Wall polygon with standard gray fill."""
    svg = " ".join(f"{to_svg(*p)[0]:.1f},{to_svg(*p)[1]:.1f}" for p in points)
    s = ' stroke="#666" stroke-width="0.8"' if stroke else ' stroke="none"'
    out.append(f'<polygon points="{svg}" fill="rgba(160,160,160,0.35)"{s}/>')

def wall_label(out, name, w, e, s, n, vertical=True):
    """Wall label text centered in bounding box."""
    mx, my = to_svg((w + e) / 2, (s + n) / 2)
    if vertical:
        lx, ly = mx - 4, my + 3.5
        out.append(f'<text x="{lx:.1f}" y="{ly:.1f}" text-anchor="middle" font-family="Arial"'
                   f' font-size="8" fill="#666" transform="rotate(-90,{lx:.1f},{ly:.1f})">{name}</text>')
    else:
        out.append(f'<text x="{mx:.1f}" y="{my-4.5:.1f}" text-anchor="middle" font-family="Arial"'
                   f' font-size="8" fill="#666">{name}</text>')

def stroke_segs(out, segs, color, width):
    """Render segment strokes (lines and arc polylines)."""
    for seg in segs:
        if isinstance(seg, LineSeg):
            sx1, sy1 = to_svg(*pts[seg.start]); sx2, sy2 = to_svg(*pts[seg.end])
            out.append(f'<line x1="{sx1:.1f}" y1="{sy1:.1f}" x2="{sx2:.1f}" y2="{sy2:.1f}"'
                       f' stroke="{color}" stroke-width="{width}"/>')
        else:
            poly = segment_polyline(seg, pts)
            svg_p = " ".join(f"{to_svg(e,n)[0]:.1f},{to_svg(e,n)[1]:.1f}" for e,n in poly)
            out.append(f'<polyline points="{svg_p}" fill="none" stroke="{color}"'
                       f' stroke-width="{width}" stroke-linecap="round"/>')

# --- Wall thickness ---
wall_t = WALL_OUTER

# --- Compute inner wall points and segments ---
inner_segs = compute_inner_walls(outline_segs, pts, wall_t, _radii)

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

# --- Interior wall: 6" thick, south face 11'6" north of inner F0-F21 wall ---
int_wall_t = 6.0 / 12.0   # 6 inches
int_wall_south = pts["W0"][1] + 11.5   # 11'6" north of inner face of F0-F21
int_wall_north = int_wall_south + int_wall_t

_s_ints = horiz_isects(inner_poly, int_wall_south)
_n_ints = horiz_isects(inner_poly, int_wall_north)
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

# Outer + inner outline strokes
stroke_segs(out, outline_segs, "#333", "1.5")
stroke_segs(out, inner_segs, "#666", "1.0")

# Interior wall IW1
iw_pts = [iw_sw, iw_se, iw_ne, iw_nw]
wall_poly(out, iw_pts, stroke=False)
for a, b in [(iw_sw, iw_se), (iw_ne, iw_nw)]:
    sx1, sy1 = to_svg(*a); sx2, sy2 = to_svg(*b)
    out.append(f'<line x1="{sx1:.1f}" y1="{sy1:.1f}" x2="{sx2:.1f}" y2="{sy2:.1f}"'
               f' stroke="#666" stroke-width="1.0"/>')
wall_label(out, "IW1", iw_sw[0], iw_se[0], int_wall_south, int_wall_north, vertical=False)

# Interior wall IW2: 6" thick, N-S, west face 6'6" east of inner W1-W2 wall
iw2_w = pts["W1"][0] + 6.5      # west face: 6'6" east of inner W1-W2
iw2_e = iw2_w + int_wall_t      # east face: +6"
iw2_s = int_wall_north           # south end: IW1 north face
iw2_n = pts["W6"][1]             # north end: inner W6-W7 wall
iw2_pts = [(iw2_w, iw2_s), (iw2_e, iw2_s), (iw2_e, iw2_n), (iw2_w, iw2_n)]
wall_poly(out, iw2_pts, stroke=False)
for e_val in [iw2_w, iw2_e]:
    sx1, sy1 = to_svg(e_val, iw2_s)
    sx2, sy2 = to_svg(e_val, iw2_n)
    out.append(f'<line x1="{sx1:.1f}" y1="{sy1:.1f}" x2="{sx2:.1f}" y2="{sy2:.1f}"'
               f' stroke="#666" stroke-width="1.0"/>')
wall_label(out, "IW2", iw2_w, iw2_e, iw2_s, iw2_n)

# Dimension line: IW1 north face to F9-F11 south face (inner), mid-span
dim_e = (pts["F9"][0] + pts["F11"][0]) / 2
dim_line_v(out, dim_e, int_wall_north, pts["W9"][1], fmt_dist(pts["W9"][1] - int_wall_north))

# Dimension line: IW2 east face to inside F12-F13 wall, vertically centered in F12-F13 wall
dim2_n = (pts["F12"][1] + pts["F13"][1]) / 2
_w9, _w8 = pts["W13"], pts["W12"]
_t_e = (dim2_n - _w9[1]) / (_w8[1] - _w9[1]) if _w8[1] != _w9[1] else 0.5
dim2_east_e = _w9[0] + _t_e * (_w8[0] - _w9[0])
dim_line_h(out, iw2_e, dim2_n, dim2_east_e, fmt_dist(dim2_east_e - iw2_e))

# --- Appliances (dimensions from ../hut project) ---
# Washer & Dryer: 35" wide (E) x 30" deep (N), same offsets as ../hut
app_w = 35.0 / 12.0   # width in feet (easting)
app_d = 30.0 / 12.0   # depth in feet (northing)
app_offset_e = 6.0 / 12.0   # 6" from interior west wall
app_offset_n = 4.0 / 12.0   # 4" from interior south wall
app_gap = 1.0 / 12.0        # 1" gap between appliances

# Dryer in F0-F1 corner: offset from inner west wall (W1) and inner south wall (W0)
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

# Water heater: 32" diameter circle, east of IW2, touching inner F7-F8 arc wall
wh_r = 14.0 / 12.0  # 14" radius in feet (28" diameter)
wh_e = iw2_e + wh_r  # center E: west side touches IW2 east face
# Inner arc W7-W8: center C7, inner radius R_a7 - wall_t
# Internal tangency: dist(C7, WH_center) = (R_a7 - wall_t) - wh_r
_wh_tangent_r = (_radii["R_a7"] - wall_t) - wh_r
_wh_dE = wh_e - pts["C7"][0]
wh_n = pts["C7"][1] + math.sqrt(_wh_tangent_r**2 - _wh_dE**2)
wh_sx, wh_sy = to_svg(wh_e, wh_n)
_wh_r_svg = (to_svg(wh_r, 0)[0] - to_svg(0, 0)[0])
out.append(f'<circle cx="{wh_sx:.1f}" cy="{wh_sy:.1f}" r="{_wh_r_svg:.1f}"'
           f' fill="rgba(100,150,200,0.2)" stroke="#4682B4" stroke-width="0.8"/>')
out.append(f'<text x="{wh_sx:.1f}" y="{wh_sy+3:.1f}" text-anchor="middle" font-family="Arial"'
           f' font-size="7" fill="#4682B4">WH</text>')

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
wall_poly(out, w8_poly)

# IW3 (west bedroom wall, 4" thick)
iw3_w = ctr_e + iw_thick_3 + closet2_width
iw3_e = iw3_w + iw_thick_4
iw3_s = ctr_s
iw3_n = int_wall_south   # extends north to IW1 south face
iw3_poly = [(iw3_w, iw3_s), (iw3_e, iw3_s), (iw3_e, iw3_n), (iw3_w, iw3_n)]
wall_poly(out, iw3_poly)
wall_label(out, "IW3", iw3_w, iw3_e, iw3_s, iw3_n)

# IW4 (bedroom east wall, 4" thick) — 11'8" east of IW3 east face
bedroom_width = 140.0 / 12.0   # 11'8"
iw4_w = iw3_e + bedroom_width
iw4_e = iw4_w + iw_thick_4
wall_south_n = 2.0 / 12.0   # south end of bedroom/closet walls: +2"
iw4_south_w = wall_south_n
iw4_south_e = wall_south_n
iw4_poly = [(iw4_w, iw4_south_w), (iw4_e, iw4_south_e), (iw4_e, int_wall_south), (iw4_w, int_wall_south)]
wall_poly(out, iw4_poly)
wall_label(out, "IW4", iw4_w, iw4_e, iw4_south_w, int_wall_south)

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
wall_poly(out, w5_poly)

# IW5: 3" thick, W-E in office, north face 30" south of IW1 south face
iw5_thick = 3.0 / 12.0
iw5_n = int_wall_south - 30.0 / 12.0   # north face
iw5_s = iw5_n - iw5_thick               # south face
iw5_w = iw4_e                            # west end at IW4 east face
iw5_e = pts["W15"][0]                    # east end at inner F14-F15 wall
iw5_poly = [(iw5_w, iw5_s), (iw5_e, iw5_s), (iw5_e, iw5_n), (iw5_w, iw5_n)]
wall_poly(out, iw5_poly, stroke=False)
for n_val in [iw5_s, iw5_n]:
    sx1, sy1 = to_svg(iw5_w, n_val)
    sx2, sy2 = to_svg(iw5_e, n_val)
    out.append(f'<line x1="{sx1:.1f}" y1="{sy1:.1f}" x2="{sx2:.1f}" y2="{sy2:.1f}"'
               f' stroke="#666" stroke-width="1.0"/>')
wall_label(out, "IW5", iw5_w, iw5_e, iw5_s, iw5_n, vertical=False)

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

# Office label (east of bedroom, above closet 1)
_of_cx = (iw4_e + pts["W15"][0]) / 2
_of_cy = (closet1_top + iw_thick_3 + int_wall_south) / 2 - 2.0
_ofx, _ofy = to_svg(_of_cx, _of_cy)
out.append(f'<text x="{_ofx:.1f}" y="{_ofy+3:.1f}" text-anchor="middle" font-family="Arial"'
           f' font-size="8" fill="#666">OFFICE</text>')

# Bedroom interior dimension lines
bd_ew_n = ctr_s + 0.25 * (int_wall_south - ctr_s)
dim_line_h(out, iw3_e, bd_ew_n, iw4_w, fmt_dist(iw4_w - iw3_e))
dim_line_v(out, iw3_e + 2.0, ctr_s, int_wall_south, fmt_dist(int_wall_south - ctr_s))

# Closet dimension lines
dim_line_v(out, (ctr_e + iw_thick_3 + iw3_w) / 2, ctr_s, ctr_n, f"CLOSET {fmt_dist(ctr_n - ctr_s)}")
dim_line_v(out, (iw4_e + w5_w) / 2, wall_south_n, closet1_top, f"CLOSET {fmt_dist(closet1_top - wall_south_n)}")

# Utility room E-W dimension
dim_line_h(out, pts["W1"][0], (ctr_s + ctr_n) / 2, ctr_e, fmt_dist(ctr_e - pts["W1"][0]))

# Office E-W dimension: Wall 5 east face to inner F14-F15 wall
dim_line_h(out, w5_e, 5.0, pts["W15"][0], fmt_dist(pts["W15"][0] - w5_e))

# Storage space between IW5 and IW1
dim_line_h(out, iw4_e, (iw5_n + int_wall_south) / 2, pts["W15"][0],
           f"STORAGE {fmt_dist(pts['W15'][0] - iw4_e)}")

# IW5 south face to F18-F19 wall north face
dim_line_v(out, pts["F18"][0], iw5_s, pts["W18"][1], fmt_dist(iw5_s - pts["W18"][1]))

# C3-F7 wall to IW1 north face
dim_line_v(out, pts["F6"][0] + 1.0, int_wall_north, pts["W6"][1],
           fmt_dist(pts["W6"][1] - int_wall_north))

# Unique ordered vertex names from outline segments (F-series)
_vert_names = []
for seg in outline_segs:
    if seg.start not in _vert_names:
        _vert_names.append(seg.start)

# Vertex label offsets (anchor, dx, dy) — local lookup, no gen_path_svg import needed
_vs_offsets = {
    "F1": ("end", -8, 0), "F2": ("end", -8, 0), "F3": ("end", -10, 0),
    "F4": ("end", -8, 0), "F5": ("end", -8, 0), "F8": ("end", -8, 0),
    "F11": ("start", 8, 0), "F12": ("start", 8, 0), "F13": ("start", 8, 0),
    "F14": ("start", 10, 0), "F15": ("start", 8, 0),
    "F0": ("middle", 0, 10), "F6": ("middle", 0, -6), "F7": ("middle", 0, -6),
    "F9": ("middle", 0, 17), "F10": ("middle", 0, 17),
    "F17": ("middle", 0, 13), "F18": ("middle", 0, 12), "F19": ("middle", 0, 12),
    "F20": ("middle", 0, 13), "F21": ("middle", 0, 10),
    "F16": ("start", 8, 4),
}
_vert_centered = {"F1","F2","F3","F4","F5","F8","F11","F12","F13","F14","F15"}
_horiz_centered = {"F0","F6","F7","F9","F10","F17","F18","F19","F20","F21"}

for f_name in _vert_names:
    sx, sy = to_svg(*pts[f_name])
    out.append(f'<circle cx="{sx:.1f}" cy="{sy:.1f}" r="1.25" fill="#333"/>')
    if f_name in _vs_offsets:
        anchor, dx, dy = _vs_offsets[f_name]
        if f_name in _vert_centered:
            out.append(f'<text x="{sx+dx:.1f}" y="{sy:.1f}" text-anchor="{anchor}"'
                       f' dominant-baseline="central"'
                       f' font-family="Arial" font-size="9" font-weight="bold"'
                       f' fill="#333">{f_name}</text>')
        elif f_name in _horiz_centered:
            out.append(f'<text x="{sx:.1f}" y="{sy+dy:.1f}" text-anchor="middle"'
                       f' font-family="Arial" font-size="9" font-weight="bold"'
                       f' fill="#333">{f_name}</text>')
        else:
            out.append(f'<text x="{sx+dx:.1f}" y="{sy+dy:.1f}" text-anchor="{anchor}"'
                       f' font-family="Arial" font-size="9" font-weight="bold"'
                       f' fill="#333">{f_name}</text>')

# North arrow
out.append('<line x1="742" y1="560" x2="742" y2="524" stroke="#333" stroke-width="2"'
           ' marker-end="url(#ah)"/>')
out.append('<text x="742" y="518" text-anchor="middle" font-family="Arial"'
           ' font-size="13" font-weight="bold">N</text>')

# Title block (right edge aligned with N arrow, bottom with F7 label)
_c4_sx, _c4_sy = to_svg(*pts["F7"])
_c4_dy = _vs_offsets["F7"][2]
tb_right = 752      # right edge, aligned with N arrow center + margin
tb_bottom = _c4_sy + _c4_dy + 4  # bottom aligned with F7 label
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
for f_name in _vert_names:
    w_name = "W" + f_name[1:]
    o = pts[f_name]; w = pts[w_name]
    print(f"  {f_name:<5s} ({o[0]:8.4f}, {o[1]:8.4f})  ->  inner ({w[0]:8.4f}, {w[1]:8.4f})")
