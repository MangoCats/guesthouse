"""Generate floorplan SVG with 8" wall inset from the outline path.

Computes geometry from shared/ and floorplan/ packages.
Outline points F0-F21, inner wall points W0-W21.
"""
import os, math, datetime

from shared.types import LineSeg, ArcSeg
from shared.geometry import (
    segment_polyline, path_polygon, poly_area,
    compute_inner_walls, horiz_isects, vert_isects, fmt_dist,
)
from shared.survey import compute_traverse, compute_three_arc, compute_inset
from shared.svg import make_svg_transform, W, H, git_describe
from floorplan.geometry import compute_outline_geometry, OutlineAnchors
from floorplan.constants import (
    WALL_OUTER, WALL_6IN, WALL_3IN, WALL_4IN,
    APPLIANCE_WIDTH, APPLIANCE_DEPTH, APPLIANCE_OFFSET_E,
    APPLIANCE_OFFSET_N, APPLIANCE_GAP,
    COUNTER_DEPTH, COUNTER_LENGTH, COUNTER_NW_RADIUS, COUNTER_GAP,
    BEDROOM_WIDTH, CLOSET_WIDTH, CLOSET1_HEIGHT,
    BED_WIDTH, BED_LENGTH, BED_OFFSET_N,
    IW1_OFFSET_N, IW2_OFFSET_E, WALL_SOUTH_N,
    WH_RADIUS, IW6_THICKNESS, IW6_OFFSET_N, IW5_OFFSET_N,
)

# ============================================================
# SVG Helpers
# ============================================================

def dim_line_h(out, e1, n, e2, label, to_svg):
    """Horizontal (E-W) dimension line with vertical tick marks."""
    x1, y1 = to_svg(e1, n); x2, y2 = to_svg(e2, n)
    _t = 4
    out.append(f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" stroke="#999" stroke-width="0.8"/>')
    out.append(f'<line x1="{x1:.1f}" y1="{y1-_t:.1f}" x2="{x1:.1f}" y2="{y1+_t:.1f}" stroke="#999" stroke-width="0.8"/>')
    out.append(f'<line x1="{x2:.1f}" y1="{y2-_t:.1f}" x2="{x2:.1f}" y2="{y2+_t:.1f}" stroke="#999" stroke-width="0.8"/>')
    out.append(f'<text x="{(x1+x2)/2:.1f}" y="{y1-3:.1f}" text-anchor="middle" font-family="Arial" font-size="8" fill="#999">{label}</text>')

def dim_line_v(out, e, n1, n2, label, to_svg):
    """Vertical (N-S) dimension line with horizontal tick marks."""
    x1, y1 = to_svg(e, n1); x2, y2 = to_svg(e, n2)
    _t = 4
    out.append(f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" stroke="#999" stroke-width="0.8"/>')
    out.append(f'<line x1="{x1-_t:.1f}" y1="{y1:.1f}" x2="{x1+_t:.1f}" y2="{y1:.1f}" stroke="#999" stroke-width="0.8"/>')
    out.append(f'<line x1="{x2-_t:.1f}" y1="{y2:.1f}" x2="{x2+_t:.1f}" y2="{y2:.1f}" stroke="#999" stroke-width="0.8"/>')
    lx, ly = x1 - 3, (y1 + y2) / 2 + 3
    out.append(f'<text x="{lx:.1f}" y="{ly:.1f}" text-anchor="middle" font-family="Arial" font-size="8" fill="#999" transform="rotate(-90,{lx:.1f},{ly:.1f})">{label}</text>')

def wall_poly(out, points, to_svg, stroke=True):
    """Wall polygon with standard gray fill."""
    svg = " ".join(f"{to_svg(*p)[0]:.1f},{to_svg(*p)[1]:.1f}" for p in points)
    s = ' stroke="#666" stroke-width="0.8"' if stroke else ' stroke="none"'
    out.append(f'<polygon points="{svg}" fill="rgba(160,160,160,0.35)"{s}/>')


def stroke_segs(out, segs, color, width, pts, to_svg):
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

# ============================================================
# Geometry computation
# ============================================================

def build_floorplan_data():
    """Compute all geometry needed for the floorplan SVG."""
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

    wall_t = WALL_OUTER
    inner_segs = compute_inner_walls(outline_segs, pts, wall_t, _radii)
    outer_poly = path_polygon(outline_segs, pts)
    inner_poly = path_polygon(inner_segs, pts)
    outer_area = poly_area(outer_poly)
    inner_area = poly_area(inner_poly)

    # --- Fit content on letter landscape (792x612) page ---
    _margin_top = 36   # 0.5" top margin
    _margin = 72       # 1" margins on left, right, bottom
    _f_svg = [to_svg(*pts[f"F{i}"]) for i in range(22)]
    _bldg_xmin = min(p[0] for p in _f_svg)
    _bldg_xmax = max(p[0] for p in _f_svg)
    _bldg_ymin = min(p[1] for p in _f_svg)
    _bldg_ymax = max(p[1] for p in _f_svg)
    _bldg_cx = (_bldg_xmin + _bldg_xmax) / 2

    _title_x = _bldg_cx
    _title_y = _bldg_ymin - 49

    _tb_w, _tb_h = 130, 92
    _tb_left = _bldg_xmax + 10
    _tb_right = _tb_left + _tb_w
    _tb_top = _title_y - 14 * 0.35
    _tb_bottom = _tb_top + _tb_h
    _tb_cx = (_tb_left + _tb_right) / 2

    _na_x = _tb_cx
    _na_text_y = _tb_bottom + 15
    _na_tip_y = _na_text_y + 6
    _na_base_y = _na_tip_y + 36

    _ext_w_x = to_svg(pts["F2"][0] - 2.7, 0)[0]
    _ext_s_y = to_svg(0, pts["F19"][1] - 3.0)[1]
    _cb_xmin = min(_bldg_xmin - 25, _ext_w_x - 10)
    _cb_xmax = _tb_right + 5
    _cb_ymin = _title_y - 14 - 5
    _cb_ymax = max(_bldg_ymax + 30, _ext_s_y + 10, _na_base_y + 5)

    # Drawing scale: 1:72 → 1 paper inch = 6 real feet
    _scale_ratio = 72
    _ft_per_inch = _scale_ratio / 12.0  # 6 feet per printed inch
    _svg_per_ft = to_svg(1, 0)[0] - to_svg(0, 0)[0]  # SVG units per survey foot
    _fit_scale = 72.0 / (_ft_per_inch * _svg_per_ft)   # CSS px per SVG unit

    _vb_w = W / _fit_scale
    _vb_h = H / _fit_scale
    _cb_cx = (_cb_xmin + _cb_xmax) / 2
    _vb_x = _cb_cx - _vb_w / 2
    _vb_y = _cb_ymin - _margin_top / _fit_scale

    return {
        "pts": pts, "to_svg": to_svg,
        "outline_segs": outline_segs, "inner_segs": inner_segs,
        "outer_poly": outer_poly, "inner_poly": inner_poly,
        "outer_area": outer_area, "inner_area": inner_area,
        "radii": _radii, "wall_t": wall_t,
        "vb_x": _vb_x, "vb_y": _vb_y, "vb_w": _vb_w, "vb_h": _vb_h,
        "title_x": _title_x, "title_y": _title_y,
        "tb_left": _tb_left, "tb_right": _tb_right, "tb_top": _tb_top,
        "tb_bottom": _tb_bottom, "tb_w": _tb_w, "tb_h": _tb_h, "tb_cx": _tb_cx,
        "na_x": _na_x, "na_text_y": _na_text_y, "na_tip_y": _na_tip_y, "na_base_y": _na_base_y,
        "ft_per_inch": _ft_per_inch,
    }

# ============================================================
# SVG rendering
# ============================================================

def render_floorplan_svg(data):
    """Render the complete floorplan SVG. Returns SVG string."""
    pts = data["pts"]
    to_svg = data["to_svg"]
    outline_segs = data["outline_segs"]
    inner_segs = data["inner_segs"]
    outer_poly = data["outer_poly"]
    inner_poly = data["inner_poly"]
    outer_area = data["outer_area"]
    inner_area = data["inner_area"]
    wall_t = data["wall_t"]

    out = []
    out.append(f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}"'
               f' viewBox="{data["vb_x"]:.2f} {data["vb_y"]:.2f} {data["vb_w"]:.2f} {data["vb_h"]:.2f}">')
    out.append(f'<rect x="{data["vb_x"]:.2f}" y="{data["vb_y"]:.2f}" width="{data["vb_w"]:.2f}" height="{data["vb_h"]:.2f}" fill="white"/>')
    out.append('<defs>')
    out.append('  <marker id="ah" markerWidth="8" markerHeight="6" refX="8" refY="3" orient="auto">'
               '<polygon points="0 0, 8 3, 0 6" fill="#333"/></marker>')
    out.append('</defs>')
    out.append(f'<text x="{data["title_x"]:.1f}" y="{data["title_y"]:.1f}" text-anchor="middle" font-family="Arial" font-size="14"'
               f' font-weight="bold">Parent Suite</text>')

    # --- Interior wall IW1: 6" thick, south face 11'6" north of inner F0-F21 ---
    int_wall_south = pts["W0"][1] + IW1_OFFSET_N
    int_wall_north = int_wall_south + WALL_6IN

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
    stroke_segs(out, outline_segs, "#333", "1.5", pts, to_svg)
    stroke_segs(out, inner_segs, "#666", "1.0", pts, to_svg)

    # Interior wall IW1
    iw_pts = [iw_sw, iw_se, iw_ne, iw_nw]
    wall_poly(out, iw_pts, to_svg, stroke=False)
    for a, b in [(iw_sw, iw_se), (iw_ne, iw_nw)]:
        sx1, sy1 = to_svg(*a); sx2, sy2 = to_svg(*b)
        out.append(f'<line x1="{sx1:.1f}" y1="{sy1:.1f}" x2="{sx2:.1f}" y2="{sy2:.1f}"'
                   f' stroke="#666" stroke-width="1.0"/>')

    # Interior wall IW2: 6" thick, N-S, west face 6'6" east of inner W1-W2 wall
    iw2_w = pts["W1"][0] + IW2_OFFSET_E
    iw2_e = iw2_w + WALL_6IN
    iw2_s = int_wall_north
    iw2_n = pts["W6"][1]
    iw2_pts = [(iw2_w, iw2_s), (iw2_e, iw2_s), (iw2_e, iw2_n), (iw2_w, iw2_n)]
    wall_poly(out, iw2_pts, to_svg, stroke=False)
    for e_val in [iw2_w, iw2_e]:
        sx1, sy1 = to_svg(e_val, iw2_s)
        sx2, sy2 = to_svg(e_val, iw2_n)
        out.append(f'<line x1="{sx1:.1f}" y1="{sy1:.1f}" x2="{sx2:.1f}" y2="{sy2:.1f}"'
                   f' stroke="#666" stroke-width="1.0"/>')

    # IW6: 1" thick, W-E, from inside of F4-F5 to IW2 west face
    iw6_n = pts["W6"][1] - IW6_OFFSET_N
    iw6_s = iw6_n - IW6_THICKNESS
    _iw6_n_ints = horiz_isects(inner_poly, iw6_n)
    _iw6_s_ints = horiz_isects(inner_poly, iw6_s)
    iw6_w_n = min(_iw6_n_ints)
    iw6_w_s = min(_iw6_s_ints)
    iw6_e = iw2_w
    iw6_poly = [(iw6_w_s, iw6_s), (iw6_e, iw6_s), (iw6_e, iw6_n), (iw6_w_n, iw6_n)]
    wall_poly(out, iw6_poly, to_svg, stroke=False)
    for w_e, n_val in [(iw6_w_s, iw6_s), (iw6_w_n, iw6_n)]:
        sx1, sy1 = to_svg(w_e, n_val)
        sx2, sy2 = to_svg(iw6_e, n_val)
        out.append(f'<line x1="{sx1:.1f}" y1="{sy1:.1f}" x2="{sx2:.1f}" y2="{sy2:.1f}"'
                   f' stroke="#666" stroke-width="1.0"/>')

    # Dimension line: IW1 north face to F9-F11 south face (inner), mid-span
    dim_e = (pts["F9"][0] + pts["F11"][0]) / 2
    dim_line_v(out, dim_e, int_wall_north, pts["W9"][1], fmt_dist(pts["W9"][1] - int_wall_north), to_svg)

    # Dimension line: IW2 east face to inside F12-F13 wall
    dim2_n = (pts["F12"][1] + pts["F13"][1]) / 2
    _w9, _w8 = pts["W13"], pts["W12"]
    _t_e = (dim2_n - _w9[1]) / (_w8[1] - _w9[1]) if _w8[1] != _w9[1] else 0.5
    dim2_east_e = _w9[0] + _t_e * (_w8[0] - _w9[0])
    dim_line_h(out, iw2_e, dim2_n, dim2_east_e, fmt_dist(dim2_east_e - iw2_e), to_svg)

    # --- Appliances ---
    dryer_w = pts["W1"][0] + APPLIANCE_OFFSET_E
    dryer_s = pts["W0"][1] + APPLIANCE_OFFSET_N
    dryer_e = dryer_w + APPLIANCE_WIDTH
    dryer_n = dryer_s + APPLIANCE_DEPTH

    washer_w = dryer_w
    washer_s = dryer_n + APPLIANCE_GAP
    washer_e = dryer_e
    washer_n = washer_s + APPLIANCE_DEPTH

    for label, sw_e, sw_n, ne_e, ne_n in [
        ("DRYER",  dryer_w,  dryer_s,  dryer_e,  dryer_n),
        ("WASHER", washer_w, washer_s, washer_e, washer_n),
    ]:
        sx1, sy1 = to_svg(sw_e, ne_n)
        sx2, sy2 = to_svg(ne_e, sw_n)
        sw = sx2 - sx1; sh = sy2 - sy1
        out.append(f'<rect x="{sx1:.1f}" y="{sy1:.1f}" width="{sw:.1f}" height="{sh:.1f}"'
                   f' fill="rgba(100,150,200,0.2)" stroke="#4682B4" stroke-width="0.8"/>')
        cx, cy = (sx1 + sx2) / 2, (sy1 + sy2) / 2
        out.append(f'<text x="{cx:.1f}" y="{cy+3:.1f}" text-anchor="middle" font-family="Arial"'
                   f' font-size="7" fill="#4682B4">{label}</text>')

    # Counter: 24" deep x 72" long, 9" NW corner radius
    ctr_w = dryer_e + COUNTER_GAP
    ctr_e = ctr_w + COUNTER_DEPTH
    ctr_s = pts["W0"][1]
    ctr_n = ctr_s + COUNTER_LENGTH

    _csw = to_svg(ctr_w, ctr_s)
    _cse = to_svg(ctr_e, ctr_s)
    _cne = to_svg(ctr_e, ctr_n)
    _cnas = to_svg(ctr_w + COUNTER_NW_RADIUS, ctr_n)
    _cnae = to_svg(ctr_w, ctr_n - COUNTER_NW_RADIUS)
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

    # Water heater: 28" diameter circle, east of IW2, touching inner F7-F8 arc wall
    wh_e = iw2_e + WH_RADIUS
    _wh_tangent_r = (data["radii"]["R_a7"] - wall_t) - WH_RADIUS
    _wh_dE = wh_e - pts["C7"][0]
    wh_n = pts["C7"][1] + math.sqrt(_wh_tangent_r**2 - _wh_dE**2)
    wh_sx, wh_sy = to_svg(wh_e, wh_n)
    _wh_r_svg = (to_svg(WH_RADIUS, 0)[0] - to_svg(0, 0)[0])
    out.append(f'<circle cx="{wh_sx:.1f}" cy="{wh_sy:.1f}" r="{_wh_r_svg:.1f}"'
               f' fill="rgba(100,150,200,0.2)" stroke="#4682B4" stroke-width="0.8"/>')
    out.append(f'<text x="{wh_sx:.1f}" y="{wh_sy+3:.1f}" text-anchor="middle" font-family="Arial"'
               f' font-size="7" fill="#4682B4">WH</text>')

    # --- Bedroom and closet walls ---
    # IW7 L-shape (west/north walls of closet, east of counter)
    iw7_poly = [
        (ctr_e, ctr_s),
        (ctr_e + WALL_3IN, ctr_s),
        (ctr_e + WALL_3IN, ctr_n),
        (ctr_e + WALL_3IN + CLOSET_WIDTH, ctr_n),
        (ctr_e + WALL_3IN + CLOSET_WIDTH, ctr_n + WALL_3IN),
        (ctr_e, ctr_n + WALL_3IN),
    ]
    wall_poly(out, iw7_poly, to_svg)

    # IW3 (west bedroom wall, 4" thick)
    iw3_w = ctr_e + WALL_3IN + CLOSET_WIDTH
    iw3_e = iw3_w + WALL_4IN
    iw3_s = ctr_s
    iw3_n = int_wall_south
    iw3_poly = [(iw3_w, iw3_s), (iw3_e, iw3_s), (iw3_e, iw3_n), (iw3_w, iw3_n)]
    wall_poly(out, iw3_poly, to_svg)

    # IW4 (bedroom east wall, 4" thick) — 11'8" east of IW3 east face
    iw4_w = iw3_e + BEDROOM_WIDTH
    iw4_e = iw4_w + WALL_4IN
    iw4_poly = [(iw4_w, WALL_SOUTH_N), (iw4_e, WALL_SOUTH_N), (iw4_e, int_wall_south), (iw4_w, int_wall_south)]
    wall_poly(out, iw4_poly, to_svg)

    # IW8 (L-shaped, 3" thick — east/north walls of closet 1)
    closet1_top = WALL_SOUTH_N + CLOSET1_HEIGHT
    w5_w = iw4_e + CLOSET_WIDTH
    w5_e = w5_w + WALL_3IN
    iw8_poly = [
        (iw4_e, closet1_top + WALL_3IN),
        (w5_e, closet1_top + WALL_3IN),
        (w5_e, WALL_SOUTH_N),
        (w5_w, WALL_SOUTH_N),
        (w5_w, closet1_top),
        (iw4_e, closet1_top),
    ]
    wall_poly(out, iw8_poly, to_svg)

    # IW5: 3" thick, W-E in office, north face 30" south of IW1 south face
    iw5_n = int_wall_south - IW5_OFFSET_N
    iw5_s = iw5_n - WALL_3IN
    iw5_w = iw4_e
    iw5_e = pts["W15"][0]
    iw5_poly = [(iw5_w, iw5_s), (iw5_e, iw5_s), (iw5_e, iw5_n), (iw5_w, iw5_n)]
    wall_poly(out, iw5_poly, to_svg, stroke=False)
    for n_val in [iw5_s, iw5_n]:
        sx1, sy1 = to_svg(iw5_w, n_val)
        sx2, sy2 = to_svg(iw5_e, n_val)
        out.append(f'<line x1="{sx1:.1f}" y1="{sy1:.1f}" x2="{sx2:.1f}" y2="{sy2:.1f}"'
                   f' stroke="#666" stroke-width="1.0"/>')

    # Subtract interior wall areas from interior area
    _iw_polys = [iw_pts, iw2_pts, iw6_poly, iw7_poly, iw3_poly, iw4_poly, iw8_poly, iw5_poly]
    inner_area -= sum(poly_area(p) for p in _iw_polys)

    # --- King Bed ---
    bed_cx = (iw3_e + iw4_w) / 2
    bed_w = bed_cx - BED_WIDTH / 2
    bed_e = bed_cx + BED_WIDTH / 2
    bed_s = ctr_s + BED_OFFSET_N
    bed_n = bed_s + BED_LENGTH
    _bed_sw = to_svg(bed_w, bed_n)
    _bed_se = to_svg(bed_e, bed_s)
    _bed_sw_x, _bed_sw_y = _bed_sw
    _bed_se_x, _bed_se_y = _bed_se
    _bed_w = _bed_se_x - _bed_sw_x
    _bed_h = _bed_se_y - _bed_sw_y
    out.append(f'<rect x="{_bed_sw_x:.1f}" y="{_bed_sw_y:.1f}" width="{_bed_w:.1f}" height="{_bed_h:.1f}"'
               f' fill="rgba(100,150,200,0.2)" stroke="#4682B4" stroke-width="0.8"/>')
    _bed_cx_svg = (_bed_sw_x + _bed_se_x) / 2
    _bed_label_y = _bed_sw_y + 0.765 * _bed_h
    out.append(f'<text x="{_bed_cx_svg:.1f}" y="{_bed_label_y+3:.1f}" text-anchor="middle" font-family="Arial"'
               f' font-size="7" fill="#4682B4">KING BED</text>')

    # Room labels
    _bd_cx = (iw3_e + iw4_w) / 2
    _bd_cy = (ctr_s + int_wall_south) / 2
    _bdx, _bdy = to_svg(_bd_cx, _bd_cy)
    out.append(f'<text x="{_bdx:.1f}" y="{_bdy+3:.1f}" text-anchor="middle" font-family="Arial"'
               f' font-size="8" fill="#666">BEDROOM</text>')

    _of_cx = (iw4_e + pts["W15"][0]) / 2
    _of_cy = (closet1_top + WALL_3IN + int_wall_south) / 2 - 2.0
    _ofx, _ofy = to_svg(_of_cx, _of_cy)
    out.append(f'<text x="{_ofx:.1f}" y="{_ofy+3:.1f}" text-anchor="middle" font-family="Arial"'
               f' font-size="8" fill="#666">OFFICE</text>')

    # --- Interior dimension lines ---
    bd_ew_n = ctr_s + 0.25 * (int_wall_south - ctr_s)
    dim_line_h(out, iw3_e, bd_ew_n, iw4_w, fmt_dist(iw4_w - iw3_e), to_svg)
    dim_line_v(out, iw3_e + 2.0, ctr_s, int_wall_south, fmt_dist(int_wall_south - ctr_s), to_svg)

    dim_line_v(out, (ctr_e + WALL_3IN + iw3_w) / 2, ctr_s, ctr_n, f"CLOSET {fmt_dist(ctr_n - ctr_s)}", to_svg)
    dim_line_v(out, (iw4_e + w5_w) / 2, WALL_SOUTH_N, closet1_top, f"CLOSET {fmt_dist(closet1_top - WALL_SOUTH_N)}", to_svg)

    dim_line_h(out, pts["W1"][0], (ctr_s + ctr_n) / 2, ctr_e, fmt_dist(ctr_e - pts["W1"][0]), to_svg)
    dim_line_h(out, w5_e, 5.0, pts["W15"][0], fmt_dist(pts["W15"][0] - w5_e), to_svg)

    dim_line_h(out, iw4_e, (iw5_n + int_wall_south) / 2, pts["W15"][0],
               f"STORAGE {fmt_dist(pts['W15'][0] - iw4_e)}", to_svg)

    _dim_f1f2_n = ctr_n + WALL_3IN + 1.0
    dim_line_h(out, pts["W2"][0], _dim_f1f2_n, iw3_w, fmt_dist(iw3_w - pts["W2"][0]), to_svg)
    dim_line_h(out, pts["W2"][0], pts["F2"][1], iw2_w, fmt_dist(iw2_w - pts["W2"][0]), to_svg)
    dim_line_h(out, pts["W5"][0], pts["F5"][1], iw2_w, fmt_dist(iw2_w - pts["W5"][0]), to_svg)

    dim_line_v(out, pts["F18"][0], iw5_s, pts["W18"][1], fmt_dist(iw5_s - pts["W18"][1]), to_svg)
    dim_line_v(out, pts["F6"][0] + 1.0, iw6_n, pts["W6"][1],
               fmt_dist(pts["W6"][1] - iw6_n), to_svg)
    dim_line_v(out, pts["F6"][0] + 1.0, int_wall_north, iw6_s,
               fmt_dist(iw6_s - int_wall_north), to_svg)

    # External dimensions
    _dim_ext_e = pts["F2"][0] - 2.7
    dim_line_v(out, _dim_ext_e, pts["F0"][1], pts["F6"][1],
               fmt_dist(pts["F6"][1] - pts["F0"][1]), to_svg)

    dim_line_h(out, pts["F8"][0], pts["F6"][1] + 1.0, pts["F11"][0],
               fmt_dist(pts["F11"][0] - pts["F8"][0]), to_svg)

    _dim_ext_n = pts["F19"][1] - 3.0
    dim_line_h(out, pts["F1"][0], _dim_ext_n, pts["F15"][0],
               fmt_dist(pts["F15"][0] - pts["F1"][0]), to_svg)

    _o9_dim_e = (bed_e + iw4_w) / 2
    dim_line_v(out, _o9_dim_e, pts["W18"][1], int_wall_south,
               fmt_dist(int_wall_south - pts["W18"][1]), to_svg)

    # --- Openings (numbered CW around the building outline) ---

    # O2: F1-F2, vertical, upper (near F2) — computed first, O1 depends on it
    _o2_n = pts["F2"][1] - 4.0 / 12.0
    _o2_s = pts["F2"][1] - 29.0 / 12.0
    _o2_poly = [
        (pts["F2"][0], _o2_s), (pts["F2"][0], _o2_n),
        (pts["W2"][0], _o2_n), (pts["W2"][0], _o2_s),
    ]

    # O1: F1-F2, vertical, lower (south of IW1)
    _o1_gap = _o2_s - int_wall_north
    _o1_n = int_wall_south - _o1_gap
    _o1_s = _o1_n - 25.0 / 12.0
    _o1_poly = [
        (pts["F2"][0], _o1_s), (pts["F2"][0], _o1_n),
        (pts["W2"][0], _o1_n), (pts["W2"][0], _o1_s),
    ]

    # O3: F4-F5, vertical
    _o3_cn = (pts["F4"][1] + pts["F5"][1]) / 2
    _o3_half = 16.0 / 12.0
    _o3_poly = [
        (pts["F4"][0], _o3_cn - _o3_half), (pts["F4"][0], _o3_cn + _o3_half),
        (pts["W4"][0], _o3_cn + _o3_half), (pts["W4"][0], _o3_cn - _o3_half),
    ]

    # O4: F6-F7, horizontal, centered on midpoint
    _o4_mid = (pts["F6"][0] + pts["F7"][0]) / 2
    _o4_w = _o4_mid - 4.5 / 12.0
    _o4_e = _o4_mid + 4.5 / 12.0
    _o4_poly = [
        (_o4_w, pts["W6"][1]), (_o4_e, pts["W6"][1]),
        (_o4_e, pts["F6"][1]), (_o4_w, pts["F6"][1]),
    ]

    # O5 & O6: F9-F10, horizontal
    # O6 edges (computed first so O5 can reference the 78" gap)
    _o6_e = pts["F10"][0] - 4.0 / 12.0
    _o6_w = pts["F10"][0] - 48.0 / 12.0
    # O5: 6' opening, 78" west of O6
    _o5_e = _o6_w - 78.0 / 12.0
    _o5_w = _o5_e - 6.0
    _o5_poly = [
        (_o5_w, pts["W9"][1]), (_o5_e, pts["W9"][1]),
        (_o5_e, pts["F9"][1]), (_o5_w, pts["F9"][1]),
    ]
    # O6
    _o6_poly = [
        (_o6_w, pts["W10"][1]), (_o6_e, pts["W10"][1]),
        (_o6_e, pts["F10"][1]), (_o6_w, pts["F10"][1]),
    ]

    # O7: F12-F13, diagonal
    _o7_dE = pts["F13"][0] - pts["F12"][0]
    _o7_dN = pts["F13"][1] - pts["F12"][1]
    _o7_len = math.sqrt(_o7_dE**2 + _o7_dN**2)
    _o7_half = 36.0 / 12.0
    _o7_ts = 0.5 - _o7_half / _o7_len
    _o7_te = 0.5 + _o7_half / _o7_len
    _o7_poly = [
        (pts["F12"][0] + _o7_ts * _o7_dE, pts["F12"][1] + _o7_ts * _o7_dN),
        (pts["F12"][0] + _o7_te * _o7_dE, pts["F12"][1] + _o7_te * _o7_dN),
        (pts["W12"][0] + _o7_te * (pts["W13"][0] - pts["W12"][0]),
         pts["W12"][1] + _o7_te * (pts["W13"][1] - pts["W12"][1])),
        (pts["W12"][0] + _o7_ts * (pts["W13"][0] - pts["W12"][0]),
         pts["W12"][1] + _o7_ts * (pts["W13"][1] - pts["W12"][1])),
    ]

    # O8: F14-F15, vertical
    _o8_cn = (iw5_s + pts["F15"][1]) / 2
    _o8_half = 12.5 / 12.0
    _o8_poly = [
        (pts["F15"][0], _o8_cn - _o8_half), (pts["F15"][0], _o8_cn + _o8_half),
        (pts["W15"][0], _o8_cn + _o8_half), (pts["W15"][0], _o8_cn - _o8_half),
    ]

    # O9: F18-F19, horizontal
    _o9_cn = (bed_e + iw4_w) / 2
    _o9_half = 12.5 / 12.0
    _o9_poly = [
        (_o9_cn - _o9_half, pts["F18"][1]), (_o9_cn + _o9_half, pts["F18"][1]),
        (_o9_cn + _o9_half, pts["W18"][1]), (_o9_cn - _o9_half, pts["W18"][1]),
    ]

    # O10: F21-F0, horizontal (bed area)
    _o10_cn = (bed_w + iw3_e) / 2
    _o10_half = 12.5 / 12.0
    _o10_poly = [
        (_o10_cn - _o10_half, pts["F0"][1]), (_o10_cn + _o10_half, pts["F0"][1]),
        (_o10_cn + _o10_half, pts["W0"][1]), (_o10_cn - _o10_half, pts["W0"][1]),
    ]

    # O11: F21-F0, horizontal (utility area)
    _o11_cn = (dryer_e + ctr_w) / 2
    _o11_half = 12.5 / 12.0
    _o11_poly = [
        (_o11_cn - _o11_half, pts["F0"][1]), (_o11_cn + _o11_half, pts["F0"][1]),
        (_o11_cn + _o11_half, pts["W0"][1]), (_o11_cn - _o11_half, pts["W0"][1]),
    ]

    _openings = [("O1", _o1_poly), ("O2", _o2_poly), ("O3", _o3_poly),
                 ("O4", _o4_poly), ("O5", _o5_poly), ("O6", _o6_poly),
                 ("O7", _o7_poly), ("O8", _o8_poly), ("O9", _o9_poly),
                 ("O10", _o10_poly), ("O11", _o11_poly)]
    _ns_openings = {"O1", "O2", "O3", "O8"}
    _o7_angle = math.degrees(math.atan2(-_o7_dN, _o7_dE))
    for _name, _poly in _openings:
        _svg = " ".join(f"{to_svg(*p)[0]:.1f},{to_svg(*p)[1]:.1f}" for p in _poly)
        out.append(f'<polygon points="{_svg}" fill="rgb(220,235,255)" stroke="#4682B4" stroke-width="1.0"/>')
        _cx = sum(p[0] for p in _poly) / len(_poly)
        _cy = sum(p[1] for p in _poly) / len(_poly)
        _sx, _sy = to_svg(_cx, _cy)
        if _name in _ns_openings:
            out.append(f'<text x="{_sx:.1f}" y="{_sy:.1f}" text-anchor="middle" dominant-baseline="central"'
                       f' font-family="Arial" font-size="7" fill="#4682B4"'
                       f' transform="rotate(-90,{_sx:.1f},{_sy:.1f})">{_name}</text>')
        elif _name == "O7":
            out.append(f'<text x="{_sx:.1f}" y="{_sy:.1f}" text-anchor="middle" dominant-baseline="central"'
                       f' font-family="Arial" font-size="7" fill="#4682B4"'
                       f' transform="rotate({_o7_angle:.1f},{_sx:.1f},{_sy:.1f})">{_name}</text>')
        else:
            out.append(f'<text x="{_sx:.1f}" y="{_sy+3:.1f}" text-anchor="middle" font-family="Arial"'
                       f' font-size="7" fill="#4682B4">{_name}</text>')

    # --- Vertex labels ---
    _vert_names = []
    for seg in outline_segs:
        if seg.start not in _vert_names:
            _vert_names.append(seg.start)

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
    out.append(f'<line x1="{data["na_x"]:.1f}" y1="{data["na_base_y"]:.1f}" x2="{data["na_x"]:.1f}" y2="{data["na_tip_y"]:.1f}" stroke="#333" stroke-width="2"'
               f' marker-end="url(#ah)"/>')
    out.append(f'<text x="{data["na_x"]:.1f}" y="{data["na_text_y"]:.1f}" text-anchor="middle" font-family="Arial"'
               f' font-size="13" font-weight="bold">N</text>')

    # Title block
    out.append(f'<rect x="{data["tb_left"]:.1f}" y="{data["tb_top"]:.1f}" width="{data["tb_w"]}" height="{data["tb_h"]}"'
               f' fill="white" stroke="#333" stroke-width="1"/>')
    out.append(f'<text x="{data["tb_cx"]:.1f}" y="{data["tb_top"]+14:.1f}" text-anchor="middle"'
               f' font-family="Arial" font-size="11" font-weight="bold" fill="#333">'
               f'{inner_area:.2f} sq ft</text>')
    out.append(f'<text x="{data["tb_cx"]:.1f}" y="{data["tb_top"]+26:.1f}" text-anchor="middle"'
               f' font-family="Arial" font-size="8" fill="#666">Interior area</text>')
    out.append(f'<text x="{data["tb_cx"]:.1f}" y="{data["tb_top"]+40:.1f}" text-anchor="middle"'
               f' font-family="Arial" font-size="11" font-weight="bold" fill="#333">'
               f'{outer_area:.2f} sq ft</text>')
    out.append(f'<text x="{data["tb_cx"]:.1f}" y="{data["tb_top"]+52:.1f}" text-anchor="middle"'
               f' font-family="Arial" font-size="8" fill="#666">Exterior area</text>')
    _ratio = data["ft_per_inch"] * 12  # paper inches to real inches
    _scale_label = f'Scale 1:{_ratio:.1f} 1&#8243; = {fmt_dist(data["ft_per_inch"])}'
    out.append(f'<text x="{data["tb_cx"]:.1f}" y="{data["tb_top"]+64:.1f}" text-anchor="middle"'
               f' font-family="Arial" font-size="8" fill="#666">{_scale_label}</text>')
    _now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    _git_desc = git_describe()
    out.append(f'<text x="{data["tb_cx"]:.1f}" y="{data["tb_top"]+76:.1f}" text-anchor="middle"'
               f' font-family="Arial" font-size="7.5" fill="#999">Generated {_now}</text>')
    out.append(f'<text x="{data["tb_cx"]:.1f}" y="{data["tb_top"]+86:.1f}" text-anchor="middle"'
               f' font-family="Arial" font-size="7.5" fill="#999">from {_git_desc}</text>')

    out.append('</svg>')

    return "\n".join(out), inner_area, outer_area, _vert_names


# ============================================================
# Main entry point
# ============================================================

if __name__ == "__main__":
    data = build_floorplan_data()
    svg_content, inner_area, outer_area, _vert_names = render_floorplan_svg(data)
    pts = data["pts"]

    svg_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "floorplan.svg")
    with open(svg_path, "w") as f:
        f.write(svg_content)

    print(f"Floorplan written to {svg_path}")
    print(f"Outer area:    {outer_area:.2f} sq ft")
    print(f"Interior area: {inner_area:.2f} sq ft")
    print(f"Wall area:     {outer_area - inner_area:.2f} sq ft")
    print()
    for f_name in _vert_names:
        w_name = "W" + f_name[1:]
        o = pts[f_name]; w = pts[w_name]
        print(f"  {f_name:<5s} ({o[0]:8.4f}, {o[1]:8.4f})  ->  inner ({w[0]:8.4f}, {w[1]:8.4f})")
