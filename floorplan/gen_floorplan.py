"""Generate floorplan SVG with 8" wall inset from the outline path.

Computes geometry from shared/ and floorplan/ packages.
Outline points F0-F21, inner wall points W0-W21.
"""
import os, math, datetime
from typing import NamedTuple, Any

from shared.types import LineSeg, ArcSeg
from shared.geometry import (
    segment_polyline, path_polygon, poly_area,
    compute_inner_walls, fmt_dist,
)
from shared.survey import compute_traverse, compute_three_arc, compute_inset
from shared.svg import make_svg_transform, W, H, git_describe
from floorplan.geometry import compute_outline_geometry, OutlineAnchors
from floorplan.constants import (
    WALL_OUTER, COUNTER_NW_RADIUS, WH_RADIUS,
    SINK_RX, SINK_RY,
    KITCHEN_SINK_WIDTH, KITCHEN_SINK_DEPTH,
    DW_WIDTH, DW_DEPTH, STOVE_WIDTH, STOVE_DEPTH,
    FRIDGE_SIZE,
    KITCHEN_CTR_LENGTH, KITCHEN_CTR_DEPTH,
    NORTH_CTR_LENGTH, NORTH_CTR_DEPTH,
    JAMB_WIDTH, STD_GAP, KITCHEN_APPL_GAP,
    WW_RADIUS,
    LOVESEAT_WIDTH, LOVESEAT_LENGTH, LOVESEAT_ANGLE_DEG,
    CHAIR_WIDTH, CHAIR_DEPTH, CHAIR_CORNER_R, CHAIR_ANGLE_DEG,
    OTTOMAN_SIZE, ET_RADIUS_CM,
    SHELVES_WIDTH, SHELVES_DEPTH,
)
from floorplan.layout import compute_interior_layout
from floorplan.openings import compute_outer_openings, compute_rough_openings

# ============================================================
# SVG Style Constants
# ============================================================
APPL_FILL = 'rgba(100,150,200,0.2)'
APPL_STROKE = '#4682B4'
APPL_SW = '0.8'
WALL_FILL = 'rgba(160,160,160,0.35)'
WALL_STROKE = '#666'
WALL_SW = '1.0'
OPENING_FILL = 'rgb(220,235,255)'
OPENING_STROKE = '#4682B4'
JAMB_COLOR = 'darkred'
DIM_COLOR = '#999'

# ============================================================
# SVG Helpers
# ============================================================

def dim_line_h(out, e1, n, e2, label, to_svg, label_offset_e=0.0):
    """Horizontal (E-W) dimension line with vertical tick marks."""
    x1, y1 = to_svg(e1, n); x2, y2 = to_svg(e2, n)
    _t = 4
    out.append(f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" stroke="{DIM_COLOR}" stroke-width="0.8"/>')
    out.append(f'<line x1="{x1:.1f}" y1="{y1-_t:.1f}" x2="{x1:.1f}" y2="{y1+_t:.1f}" stroke="{DIM_COLOR}" stroke-width="0.8"/>')
    out.append(f'<line x1="{x2:.1f}" y1="{y2-_t:.1f}" x2="{x2:.1f}" y2="{y2+_t:.1f}" stroke="{DIM_COLOR}" stroke-width="0.8"/>')
    lx, _ = to_svg((e1 + e2) / 2 + label_offset_e, n)
    out.append(f'<text x="{lx:.1f}" y="{y1-3:.1f}" text-anchor="middle" font-family="Arial" font-size="8" fill="{DIM_COLOR}">{label}</text>')

def dim_line_v(out, e, n1, n2, label, to_svg):
    """Vertical (N-S) dimension line with horizontal tick marks."""
    x1, y1 = to_svg(e, n1); x2, y2 = to_svg(e, n2)
    _t = 4
    out.append(f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" stroke="{DIM_COLOR}" stroke-width="0.8"/>')
    out.append(f'<line x1="{x1-_t:.1f}" y1="{y1:.1f}" x2="{x1+_t:.1f}" y2="{y1:.1f}" stroke="{DIM_COLOR}" stroke-width="0.8"/>')
    out.append(f'<line x1="{x2-_t:.1f}" y1="{y2:.1f}" x2="{x2+_t:.1f}" y2="{y2:.1f}" stroke="{DIM_COLOR}" stroke-width="0.8"/>')
    lx, ly = x1 - 3, (y1 + y2) / 2 + 3
    out.append(f'<text x="{lx:.1f}" y="{ly:.1f}" text-anchor="middle" font-family="Arial" font-size="8" fill="{DIM_COLOR}" transform="rotate(-90,{lx:.1f},{ly:.1f})">{label}</text>')

def wall_poly(out, points, to_svg, stroke=True):
    """Wall polygon with standard gray fill.  Stroke is inside-only via clip-path."""
    svg = " ".join(f"{to_svg(*p)[0]:.1f},{to_svg(*p)[1]:.1f}" for p in points)
    if stroke:
        cid = f"wc{len(out)}"
        out.append(f'<defs><clipPath id="{cid}"><polygon points="{svg}"/></clipPath></defs>')
        out.append(f'<polygon points="{svg}" fill="{WALL_FILL}"'
                   f' stroke="{WALL_STROKE}" stroke-width="1.6" clip-path="url(#{cid})"/>')
    else:
        out.append(f'<polygon points="{svg}" fill="{WALL_FILL}" stroke="none"/>')


# Toilet plan-view shape: (dx, dy) offsets in source SVG units from center of
# back edge, +dy = toward bowl.  Extracted from hut/floor_plan_2d.svg lines 264-269.
# 1 SVG unit = 10 cm; _SVG_TO_FT converts to feet.
_TOILET_SVG = [
    (-1.905, 0), (-1.905, 2.032), (-0.841, 2.032),
    (-1.078, 2.224), (-1.267, 2.455), (-1.408, 2.719),
    (-1.495, 3.005), (-1.524, 3.302),
    (-1.732, 5.461), (-1.699, 5.799), (-1.600, 6.124),
    (-1.440, 6.423), (-1.225, 6.686), (-0.962, 6.901),
    (-0.663, 7.061), (-0.338, 7.160), (0, 7.193),
    (0.338, 7.160), (0.663, 7.061), (0.962, 6.901),
    (1.225, 6.686), (1.440, 6.423), (1.600, 6.124),
    (1.699, 5.799), (1.732, 5.461),
    (1.524, 3.302), (1.495, 3.005), (1.408, 2.719),
    (1.267, 2.455), (1.078, 2.224), (0.847, 2.035),
    (0.841, 2.032), (1.905, 2.032), (1.905, 0),
]
_SVG_TO_FT = 10.0 / 30.48


def draw_toilet(out, center_e, back_n, face_north, to_svg):
    """Draw a toilet plan view against a wall.

    center_e: easting of toilet center
    back_n: northing of the wall face the tank sits against
    face_north: True = bowl faces north (+N), False = bowl faces south (-N)
    """
    sign = 1 if face_north else -1
    pts_survey = [(center_e + dx * _SVG_TO_FT, back_n + sign * dy * _SVG_TO_FT)
                  for dx, dy in _TOILET_SVG]
    svg_pts = " ".join(f"{to_svg(e, n)[0]:.1f},{to_svg(e, n)[1]:.1f}"
                       for e, n in pts_survey)
    out.append(f'<polygon points="{svg_pts}"'
               f' fill="{APPL_FILL}" stroke="{APPL_STROKE}" stroke-width="{APPL_SW}"/>')
    # Label at centroid
    cx = sum(p[0] for p in pts_survey) / len(pts_survey)
    cy = sum(p[1] for p in pts_survey) / len(pts_survey)
    sx, sy = to_svg(cx, cy)
    out.append(f'<text x="{sx:.1f}" y="{sy+3:.1f}" text-anchor="middle" font-family="Arial"'
               f' font-size="7" fill="{APPL_STROKE}">TOILET</text>')


def draw_sink(out, center_e, center_n, to_svg):
    """Draw a sink plan view as an ellipse."""
    sx, sy = to_svg(center_e, center_n)
    # Convert radii from feet to SVG pixel units
    rx_svg = abs(to_svg(SINK_RX, 0)[0] - to_svg(0, 0)[0])
    ry_svg = abs(to_svg(0, SINK_RY)[1] - to_svg(0, 0)[1])
    out.append(f'<ellipse cx="{sx:.1f}" cy="{sy:.1f}" rx="{rx_svg:.1f}" ry="{ry_svg:.1f}"'
               f' fill="{APPL_FILL}" stroke="{APPL_STROKE}" stroke-width="{APPL_SW}"/>')
    out.append(f'<text x="{sx:.1f}" y="{sy+3:.1f}" text-anchor="middle" font-family="Arial"'
               f' font-size="7" fill="{APPL_STROKE}">SINK</text>')


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

class FloorplanData(NamedTuple):
    """Typed container for all floorplan geometry and page layout."""
    pts: dict
    to_svg: Any  # Callable[[float, float], tuple[float, float]]
    outline_segs: list
    inner_segs: list
    outer_poly: list
    inner_poly: list
    outer_area: float
    inner_area: float
    radii: dict
    wall_t: float
    vb_x: float
    vb_y: float
    vb_w: float
    vb_h: float
    title_x: float
    title_y: float
    tb_left: float
    tb_right: float
    tb_top: float
    tb_bottom: float
    tb_w: float
    tb_h: float
    tb_cx: float
    na_x: float
    na_text_y: float
    na_tip_y: float
    na_base_y: float
    ft_per_inch: float


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

    return FloorplanData(
        pts=pts, to_svg=to_svg,
        outline_segs=outline_segs, inner_segs=inner_segs,
        outer_poly=outer_poly, inner_poly=inner_poly,
        outer_area=outer_area, inner_area=inner_area,
        radii=_radii, wall_t=wall_t,
        vb_x=_vb_x, vb_y=_vb_y, vb_w=_vb_w, vb_h=_vb_h,
        title_x=_title_x, title_y=_title_y,
        tb_left=_tb_left, tb_right=_tb_right, tb_top=_tb_top,
        tb_bottom=_tb_bottom, tb_w=_tb_w, tb_h=_tb_h, tb_cx=_tb_cx,
        na_x=_na_x, na_text_y=_na_text_y, na_tip_y=_na_tip_y, na_base_y=_na_base_y,
        ft_per_inch=_ft_per_inch,
    )

# ============================================================
# Render sub-functions
# ============================================================

def compute_iw_area(layout):
    """Compute total interior wall area from layout polygons."""
    iw2 = layout.iw2
    iw3 = layout.iw3
    iw5 = layout.iw5
    iw2_poly = [(iw2.w, iw2.s), (iw2.e, iw2.s), (iw2.e, iw2.n), (iw2.w, iw2.n)]
    iw3_poly = [(iw3.w, iw3.s), (iw3.e, iw3.s), (iw3.e, iw3.n), (iw3.w, iw3.n)]
    iw4_poly = [(layout.iw4_w, layout.wall_south_n), (layout.iw4_e, layout.wall_south_n),
                (layout.iw4_e, layout.iw1_s), (layout.iw4_w, layout.iw1_s)]
    iw5_poly = [(iw5.w, iw5.s), (iw5.e, iw5.s), (iw5.e, iw5.n), (iw5.w, iw5.n)]
    iw_polys = [layout.iw1, iw2_poly, layout.iw6_poly, layout.iw7,
                iw3_poly, iw4_poly, layout.iw8, iw5_poly]
    return sum(poly_area(p) for p in iw_polys)


def _render_walls(out, data, layout):
    """Render outer wall fill, outline strokes, and all interior walls with rough openings."""
    pts = data.pts
    to_svg = data.to_svg

    # Outer wall fill with inner cutout
    outer_svg = " ".join(f"{to_svg(*p)[0]:.1f},{to_svg(*p)[1]:.1f}" for p in data.outer_poly)
    inner_rev = list(reversed(data.inner_poly))
    inner_svg = " ".join(f"{to_svg(*p)[0]:.1f},{to_svg(*p)[1]:.1f}" for p in inner_rev)
    out.append(f'<polygon points="{outer_svg}" fill="{WALL_FILL}" stroke="none"/>')
    out.append(f'<polygon points="{inner_svg}" fill="white" stroke="none"/>')

    # Outline strokes
    stroke_segs(out, data.outline_segs, "#333", "1.5", pts, to_svg)
    stroke_segs(out, data.inner_segs, WALL_STROKE, WALL_SW, pts, to_svg)

    # Half stroke width in survey feet (for inside-only edge lines)
    svg_per_ft = abs(to_svg(1, 0)[0] - to_svg(0, 0)[0])
    half_sw = 0.5 / svg_per_ft

    # Rough openings
    rough_openings = compute_rough_openings(pts, layout)
    ro = {r.name: r.bbox for r in rough_openings}

    # ---- IW1 with RO1 ----
    iw_sw, iw_se, iw_ne, iw_nw = layout.iw1
    iw1_s, iw1_n = layout.iw1_s, layout.iw1_n
    ro1_w, ro1_e = ro["RO1"].w, ro["RO1"].e

    iw1_w_poly = [iw_sw, (ro1_w, iw1_s), (ro1_w, iw1_n), iw_nw]
    iw1_e_poly = [(ro1_e, iw1_s), iw_se, iw_ne, (ro1_e, iw1_n)]
    wall_poly(out, iw1_w_poly, to_svg, stroke=False)
    wall_poly(out, iw1_e_poly, to_svg, stroke=False)
    s_in = iw1_s + half_sw
    n_in = iw1_n - half_sw
    for a, b in [((iw_sw[0], s_in), (ro1_w, s_in)),
                 ((ro1_e, s_in), (iw_se[0], s_in)),
                 ((iw_nw[0], n_in), (ro1_w, n_in)),
                 ((ro1_e, n_in), (iw_ne[0], n_in))]:
        sx1, sy1 = to_svg(*a); sx2, sy2 = to_svg(*b)
        out.append(f'<line x1="{sx1:.1f}" y1="{sy1:.1f}" x2="{sx2:.1f}" y2="{sy2:.1f}"'
                   f' stroke="{WALL_STROKE}" stroke-width="{WALL_SW}"/>')
    for jamb_e in [ro1_w, ro1_e - JAMB_WIDTH]:
        jx1, jy1 = to_svg(jamb_e, iw1_n)
        jx2, jy2 = to_svg(jamb_e + JAMB_WIDTH, iw1_s)
        out.append(f'<rect x="{jx1:.1f}" y="{jy1:.1f}" width="{jx2 - jx1:.1f}" height="{jy2 - jy1:.1f}"'
                   f' fill="{JAMB_COLOR}" stroke="none"/>')

    # ---- IW2 with RO4 ----
    iw2 = layout.iw2
    ro4_n, ro4_s = ro["RO4"].n, ro["RO4"].s

    iw2_s_poly = [(iw2.w, iw2.s), (iw2.e, iw2.s), (iw2.e, ro4_s), (iw2.w, ro4_s)]
    iw2_n_poly = [(iw2.w, ro4_n), (iw2.e, ro4_n), (iw2.e, iw2.n), (iw2.w, iw2.n)]
    wall_poly(out, iw2_s_poly, to_svg, stroke=False)
    wall_poly(out, iw2_n_poly, to_svg, stroke=False)
    iw2_w_in = iw2.w + half_sw
    iw2_e_in = iw2.e - half_sw
    for a, b in [((iw2_w_in, iw2.s), (iw2_w_in, ro4_s)),
                 ((iw2_w_in, ro4_n), (iw2_w_in, iw2.n)),
                 ((iw2_e_in, iw2.s), (iw2_e_in, ro4_s)),
                 ((iw2_e_in, ro4_n), (iw2_e_in, iw2.n))]:
        sx1, sy1 = to_svg(*a); sx2, sy2 = to_svg(*b)
        out.append(f'<line x1="{sx1:.1f}" y1="{sy1:.1f}" x2="{sx2:.1f}" y2="{sy2:.1f}"'
                   f' stroke="{WALL_STROKE}" stroke-width="{WALL_SW}"/>')
    for jamb_n in [ro4_s, ro4_n - JAMB_WIDTH]:
        jx1, jy1 = to_svg(iw2.w, jamb_n + JAMB_WIDTH)
        jx2, jy2 = to_svg(iw2.e, jamb_n)
        out.append(f'<rect x="{jx1:.1f}" y="{jy1:.1f}" width="{jx2 - jx1:.1f}" height="{jy2 - jy1:.1f}"'
                   f' fill="{JAMB_COLOR}" stroke="none"/>')

    # ---- IW6 with RO5 ----
    iw6_s, iw6_n = layout.iw6_s, layout.iw6_n
    iw6_w_s = layout.iw6_poly[0][0]
    iw6_w_n = layout.iw6_poly[3][0]
    iw6_e = iw2.w
    ro5_e, ro5_w = ro["RO5"].e, ro["RO5"].w

    iw6_w_poly = [(iw6_w_s, iw6_s), (ro5_w, iw6_s), (ro5_w, iw6_n), (iw6_w_n, iw6_n)]
    iw6_e_poly = [(ro5_e, iw6_s), (iw6_e, iw6_s), (iw6_e, iw6_n), (ro5_e, iw6_n)]
    wall_poly(out, iw6_w_poly, to_svg, stroke=False)
    wall_poly(out, iw6_e_poly, to_svg, stroke=False)
    iw6_s_in = iw6_s + half_sw
    iw6_n_in = iw6_n - half_sw
    for a, b in [((iw6_w_s, iw6_s_in), (ro5_w, iw6_s_in)),
                 ((ro5_e, iw6_s_in), (iw6_e, iw6_s_in)),
                 ((iw6_w_n, iw6_n_in), (ro5_w, iw6_n_in)),
                 ((ro5_e, iw6_n_in), (iw6_e, iw6_n_in))]:
        sx1, sy1 = to_svg(*a); sx2, sy2 = to_svg(*b)
        out.append(f'<line x1="{sx1:.1f}" y1="{sy1:.1f}" x2="{sx2:.1f}" y2="{sy2:.1f}"'
                   f' stroke="{WALL_STROKE}" stroke-width="{WALL_SW}"/>')
    for jamb_e in [ro5_w, ro5_e - JAMB_WIDTH]:
        jx1, jy1 = to_svg(jamb_e, iw6_n)
        jx2, jy2 = to_svg(jamb_e + JAMB_WIDTH, iw6_s)
        out.append(f'<rect x="{jx1:.1f}" y="{jy1:.1f}" width="{jx2 - jx1:.1f}" height="{jy2 - jy1:.1f}"'
                   f' fill="{JAMB_COLOR}" stroke="none"/>')

    # IW1-north → F9-F11 south face dimension
    dim_e = (pts["F9"][0] + pts["F11"][0]) / 2
    dim_line_v(out, dim_e, iw1_n, pts["W9"][1], fmt_dist(pts["W9"][1] - iw1_n), to_svg)

    # IW2-east → inside F12-F13 wall dimension
    dim2_n = (pts["F12"][1] + pts["F13"][1]) / 2
    w13, w12 = pts["W13"], pts["W12"]
    t_e = (dim2_n - w13[1]) / (w12[1] - w13[1]) if w12[1] != w13[1] else 0.5
    dim2_east_e = w13[0] + t_e * (w12[0] - w13[0])
    dim_line_h(out, iw2.e, dim2_n, dim2_east_e, fmt_dist(dim2_east_e - iw2.e), to_svg,
               label_offset_e=-4.0)

    # ---- IW7 ----
    wall_poly(out, layout.iw7, to_svg)

    # ---- IW3 with RO3 ----
    iw3 = layout.iw3
    ro3_s, ro3_n = ro["RO3"].s, ro["RO3"].n

    iw3_s_poly = [(iw3.w, iw3.s), (iw3.e, iw3.s), (iw3.e, ro3_s), (iw3.w, ro3_s)]
    iw3_n_poly = [(iw3.w, ro3_n), (iw3.e, ro3_n), (iw3.e, iw3.n), (iw3.w, iw3.n)]
    wall_poly(out, iw3_s_poly, to_svg, stroke=False)
    wall_poly(out, iw3_n_poly, to_svg, stroke=False)
    iw3_w_in = iw3.w + half_sw
    iw3_e_in = iw3.e - half_sw
    for a, b in [((iw3_w_in, iw3.s), (iw3_w_in, ro3_s)),
                 ((iw3_w_in, ro3_n), (iw3_w_in, iw3.n)),
                 ((iw3_e_in, iw3.s), (iw3_e_in, ro3_s)),
                 ((iw3_e_in, ro3_n), (iw3_e_in, iw3.n))]:
        sx1, sy1 = to_svg(*a); sx2, sy2 = to_svg(*b)
        out.append(f'<line x1="{sx1:.1f}" y1="{sy1:.1f}" x2="{sx2:.1f}" y2="{sy2:.1f}"'
                   f' stroke="{WALL_STROKE}" stroke-width="{WALL_SW}"/>')
    for jamb_n in [ro3_s, ro3_n - JAMB_WIDTH]:
        jx1, jy1 = to_svg(iw3.w, jamb_n + JAMB_WIDTH)
        jx2, jy2 = to_svg(iw3.e, jamb_n)
        out.append(f'<rect x="{jx1:.1f}" y="{jy1:.1f}" width="{jx2 - jx1:.1f}" height="{jy2 - jy1:.1f}"'
                   f' fill="{JAMB_COLOR}" stroke="none"/>')

    # ---- IW4 with RO2 ----
    ro2_s, ro2_n = ro["RO2"].s, ro["RO2"].n
    iw4_poly = [(layout.iw4_w, layout.wall_south_n), (layout.iw4_e, layout.wall_south_n),
                (layout.iw4_e, layout.iw1_s), (layout.iw4_w, layout.iw1_s)]

    iw4_s_poly = [(layout.iw4_w, layout.wall_south_n), (layout.iw4_e, layout.wall_south_n),
                  (layout.iw4_e, ro2_s), (layout.iw4_w, ro2_s)]
    iw4_n_poly = [(layout.iw4_w, ro2_n), (layout.iw4_e, ro2_n),
                  (layout.iw4_e, layout.iw1_s), (layout.iw4_w, layout.iw1_s)]
    wall_poly(out, iw4_s_poly, to_svg, stroke=False)
    wall_poly(out, iw4_n_poly, to_svg, stroke=False)
    w_in = layout.iw4_w + half_sw
    e_in = layout.iw4_e - half_sw
    for a, b in [((w_in, layout.wall_south_n), (w_in, ro2_s)),
                 ((w_in, ro2_n), (w_in, layout.iw1_s)),
                 ((e_in, layout.wall_south_n), (e_in, ro2_s)),
                 ((e_in, ro2_n), (e_in, layout.iw1_s))]:
        sx1, sy1 = to_svg(*a); sx2, sy2 = to_svg(*b)
        out.append(f'<line x1="{sx1:.1f}" y1="{sy1:.1f}" x2="{sx2:.1f}" y2="{sy2:.1f}"'
                   f' stroke="{WALL_STROKE}" stroke-width="{WALL_SW}"/>')
    for jamb_n in [ro2_s, ro2_n - JAMB_WIDTH]:
        jx1, jy1 = to_svg(layout.iw4_w, jamb_n + JAMB_WIDTH)
        jx2, jy2 = to_svg(layout.iw4_e, jamb_n)
        out.append(f'<rect x="{jx1:.1f}" y="{jy1:.1f}" width="{jx2 - jx1:.1f}" height="{jy2 - jy1:.1f}"'
                   f' fill="{JAMB_COLOR}" stroke="none"/>')

    # ---- IW8 ----
    wall_poly(out, layout.iw8, to_svg)

    # ---- IW5 ----
    iw5 = layout.iw5
    iw5_poly = [(iw5.w, iw5.s), (iw5.e, iw5.s), (iw5.e, iw5.n), (iw5.w, iw5.n)]
    wall_poly(out, iw5_poly, to_svg, stroke=False)
    for n_val in [iw5.s + half_sw, iw5.n - half_sw]:
        sx1, sy1 = to_svg(iw5.w, n_val)
        sx2, sy2 = to_svg(iw5.e, n_val)
        out.append(f'<line x1="{sx1:.1f}" y1="{sy1:.1f}" x2="{sx2:.1f}" y2="{sy2:.1f}"'
                   f' stroke="{WALL_STROKE}" stroke-width="{WALL_SW}"/>')



def _render_appliances(out, data, layout):
    """Render utility room appliances: dryer, washer, counter, water heater, toilets, sinks."""
    pts = data.pts
    to_svg = data.to_svg

    # Dryer and washer
    for label, b in [("DRYER", layout.dryer), ("WASHER", layout.washer)]:
        sx1, sy1 = to_svg(b.w, b.n)
        sx2, sy2 = to_svg(b.e, b.s)
        sw = sx2 - sx1; sh = sy2 - sy1
        out.append(f'<rect x="{sx1:.1f}" y="{sy1:.1f}" width="{sw:.1f}" height="{sh:.1f}"'
                   f' fill="{APPL_FILL}" stroke="{APPL_STROKE}" stroke-width="{APPL_SW}"/>')
        cx, cy = (sx1 + sx2) / 2, (sy1 + sy2) / 2
        out.append(f'<text x="{cx:.1f}" y="{cy+3:.1f}" text-anchor="middle" font-family="Arial"'
                   f' font-size="7" fill="{APPL_STROKE}">{label}</text>')

    # Counter: 24" deep x 72" long, 9" NW corner radius
    ctr = layout.ctr
    csw = to_svg(ctr.w, ctr.s)
    cse = to_svg(ctr.e, ctr.s)
    cne = to_svg(ctr.e, ctr.n)
    cnas = to_svg(ctr.w + COUNTER_NW_RADIUS, ctr.n)
    cnae = to_svg(ctr.w, ctr.n - COUNTER_NW_RADIUS)
    r_svg = abs(cnas[0] - to_svg(ctr.w, ctr.n)[0])
    ctr_path = (f'M {csw[0]:.1f},{csw[1]:.1f} '
                f'L {cse[0]:.1f},{cse[1]:.1f} '
                f'L {cne[0]:.1f},{cne[1]:.1f} '
                f'L {cnas[0]:.1f},{cnas[1]:.1f} '
                f'A {r_svg:.1f} {r_svg:.1f} 0 0 0 {cnae[0]:.1f},{cnae[1]:.1f} '
                f'Z')
    out.append(f'<path d="{ctr_path}" fill="{APPL_FILL}" stroke="{APPL_STROKE}" stroke-width="{APPL_SW}"/>')
    ccx = (csw[0] + cse[0]) / 2
    ccy = (csw[1] + cne[1]) / 2
    out.append(f'<text x="{ccx:.1f}" y="{ccy:.1f}" text-anchor="middle" font-family="Arial"'
               f' font-size="7" fill="{APPL_STROKE}" letter-spacing="0.5" transform="rotate(-90,{ccx:.1f},{ccy:.1f})">COUNTER</text>')

    # Water heater: 28" diameter circle
    wh_e = layout.iw2.e + WH_RADIUS
    wh_tangent_r = (data.radii["R_a7"] - data.wall_t) - WH_RADIUS
    wh_dE = wh_e - pts["C7"][0]
    wh_n = pts["C7"][1] + math.sqrt(wh_tangent_r**2 - wh_dE**2)
    wh_sx, wh_sy = to_svg(wh_e, wh_n)
    wh_r_svg = (to_svg(WH_RADIUS, 0)[0] - to_svg(0, 0)[0])
    out.append(f'<circle cx="{wh_sx:.1f}" cy="{wh_sy:.1f}" r="{wh_r_svg:.1f}"'
               f' fill="{APPL_FILL}" stroke="{APPL_STROKE}" stroke-width="{APPL_SW}"/>')
    out.append(f'<text x="{wh_sx:.1f}" y="{wh_sy+3:.1f}" text-anchor="middle" font-family="Arial"'
               f' font-size="7" fill="{APPL_STROKE}">WH</text>')

    # Toilets and sinks
    toilet_e = (layout.dryer.w + layout.dryer.e) / 2
    sink_e = (layout.dryer.e + layout.ctr.w) / 2
    draw_toilet(out, toilet_e, layout.iw1_s, face_north=False, to_svg=to_svg)
    draw_sink(out, sink_e, layout.iw1_s - SINK_RY, to_svg=to_svg)
    draw_toilet(out, toilet_e, layout.iw1_n, face_north=True, to_svg=to_svg)
    draw_sink(out, sink_e, layout.iw1_n + SINK_RY, to_svg=to_svg)


def _render_kitchen(out, data, layout):
    """Render kitchen: D/W, sink, stove, shelves, fridge, counters.

    Returns ((ww1_cx, ww1_cy, ww1_r), (ww3_cx, ww3_cy, ww3_r)) work zone info.
    """
    pts = data.pts
    to_svg = data.to_svg
    back_n = pts["W9"][1]

    # Kitchen appliances
    dw_w = layout.iw2.e + NORTH_CTR_LENGTH + KITCHEN_APPL_GAP
    dw_e = dw_w + DW_WIDTH
    ks_w = dw_e + KITCHEN_APPL_GAP
    ks_e = ks_w + KITCHEN_SINK_WIDTH
    st_w = ks_e + KITCHEN_APPL_GAP
    st_e = st_w + STOVE_WIDTH
    appliances = [
        ("SINK",  ks_w, back_n - KITCHEN_SINK_DEPTH, ks_e, back_n,
         "https://www.webstaurantstore.com/advance-tabco-fs1181824l-45-fabricated-one-compartment-sink-with-24-left-drainboard-18-x-18-x-14-bowl/109FS1L241818.html"),
        ("D/W",   dw_w, back_n - DW_DEPTH,           dw_e, back_n, None),
        ("STOVE", st_w, back_n - KITCHEN_APPL_GAP - STOVE_DEPTH, st_e, back_n - KITCHEN_APPL_GAP, None),
    ]
    for label, sw_e, sw_n, ne_e, ne_n, href in appliances:
        sx1, sy1 = to_svg(sw_e, ne_n)
        sx2, sy2 = to_svg(ne_e, sw_n)
        sw = sx2 - sx1; sh = sy2 - sy1
        if href:
            out.append(f'<a href="{href}" target="_blank">')
        out.append(f'<rect x="{sx1:.1f}" y="{sy1:.1f}" width="{sw:.1f}" height="{sh:.1f}"'
                   f' fill="{APPL_FILL}" stroke="{APPL_STROKE}" stroke-width="{APPL_SW}"/>')
        cx, cy = (sx1 + sx2) / 2, (sy1 + sy2) / 2
        out.append(f'<text x="{cx:.1f}" y="{cy+3:.1f}" text-anchor="middle" font-family="Arial"'
                   f' font-size="7" fill="{APPL_STROKE}">{label}</text>')
        if href:
            out.append('</a>')

    # WW1: 30" radius constraint circle centered on SE corner of stove
    ww1_cx = st_e
    ww1_cy = back_n - KITCHEN_APPL_GAP - STOVE_DEPTH

    # SHELVES: 36" E-W x 15" N-S, against F9-F10 south face, 3" east of stove
    sh_w = st_e + KITCHEN_APPL_GAP
    sh_e = sh_w + SHELVES_WIDTH
    sh_n = back_n
    sh_s = sh_n - SHELVES_DEPTH
    sh_sx1, sh_sy1 = to_svg(sh_w, sh_n)
    sh_sx2, sh_sy2 = to_svg(sh_e, sh_s)
    sh_sw = sh_sx2 - sh_sx1; sh_sh = sh_sy2 - sh_sy1
    out.append('<a href="https://www.ikea.com/us/en/p/hemnes-bookcase-white-stain-light-brown-60413502/" target="_blank">')
    out.append(f'<rect x="{sh_sx1:.1f}" y="{sh_sy1:.1f}" width="{sh_sw:.1f}" height="{sh_sh:.1f}"'
               f' fill="{APPL_FILL}" stroke="{APPL_STROKE}" stroke-width="{APPL_SW}"/>')
    sh_cx = (sh_sx1 + sh_sx2) / 2
    sh_cy = (sh_sy1 + sh_sy2) / 2
    out.append(f'<text x="{sh_cx:.1f}" y="{sh_cy+3:.1f}" text-anchor="middle" font-family="Arial"'
               f' font-size="6" fill="{APPL_STROKE}">SHELVES</text>')
    out.append('</a>')

    # WW3: 30" radius constraint circle centered on SE corner of SHELVES
    ww3_cx = sh_e
    ww3_cy = sh_s

    # Fridge: 2" east of kitchen counter, 2" north of IW1 north face
    fr_w = layout.iw2.e + KITCHEN_CTR_LENGTH + STD_GAP
    fr_e = fr_w + FRIDGE_SIZE
    fr_s = layout.iw1_n + STD_GAP
    fr_n = fr_s + FRIDGE_SIZE
    sx1, sy1 = to_svg(fr_w, fr_n)
    sx2, sy2 = to_svg(fr_e, fr_s)
    sw = sx2 - sx1; sh = sy2 - sy1
    out.append(f'<rect x="{sx1:.1f}" y="{sy1:.1f}" width="{sw:.1f}" height="{sh:.1f}"'
               f' fill="{APPL_FILL}" stroke="{APPL_STROKE}" stroke-width="{APPL_SW}"/>')
    fr_cx = (sx1 + sx2) / 2
    fr_cy = (sy1 + sy2) / 2
    out.append(f'<text x="{fr_cx:.1f}" y="{fr_cy+3:.1f}" text-anchor="middle" font-family="Arial"'
               f' font-size="7" fill="{APPL_STROKE}">FRIDGE</text>')

    # Kitchen counter: along IW1 north face starting at IW2 east face
    kc_w = layout.iw2.e
    kc_e = kc_w + KITCHEN_CTR_LENGTH
    kc_s = layout.iw1_n
    kc_n = kc_s + KITCHEN_CTR_DEPTH
    kc_sx1, kc_sy1 = to_svg(kc_w, kc_n)
    kc_sx2, kc_sy2 = to_svg(kc_e, kc_s)
    kc_sw = kc_sx2 - kc_sx1; kc_sh = kc_sy2 - kc_sy1
    out.append('<a href="https://www.webstaurantstore.com/regency-spec-line-30-x-72-14-gauge-stainless-steel-commercial-work-table-with-4-backsplash-and-undershelf/600TSSB3072S.html" target="_blank">')
    out.append(f'<rect x="{kc_sx1:.1f}" y="{kc_sy1:.1f}" width="{kc_sw:.1f}" height="{kc_sh:.1f}"'
               f' fill="{APPL_FILL}" stroke="{APPL_STROKE}" stroke-width="{APPL_SW}"/>')
    kc_cx = (kc_sx1 + kc_sx2) / 2
    kc_cy = (kc_sy1 + kc_sy2) / 2
    out.append(f'<text x="{kc_cx:.1f}" y="{kc_cy+3:.1f}" text-anchor="middle" font-family="Arial"'
               f' font-size="7" fill="{APPL_STROKE}">COUNTER</text>')
    out.append('</a>')

    # North wall counter: south side against W9-W10, starting at IW2 east face
    nc_w = layout.iw2.e
    nc_e = nc_w + NORTH_CTR_LENGTH
    nc_n = pts["W9"][1]
    nc_s = nc_n - NORTH_CTR_DEPTH
    nc_sx1, nc_sy1 = to_svg(nc_w, nc_n)
    nc_sx2, nc_sy2 = to_svg(nc_e, nc_s)
    nc_sw = nc_sx2 - nc_sx1; nc_sh = nc_sy2 - nc_sy1
    out.append('<a href="https://www.webstaurantstore.com/regency-spec-line-30-x-36-14-gauge-stainless-steel-commercial-work-table-with-4-backsplash-and-undershelf/600TSSB3036S.html" target="_blank">')
    out.append(f'<rect x="{nc_sx1:.1f}" y="{nc_sy1:.1f}" width="{nc_sw:.1f}" height="{nc_sh:.1f}"'
               f' fill="{APPL_FILL}" stroke="{APPL_STROKE}" stroke-width="{APPL_SW}"/>')
    nc_cx = (nc_sx1 + nc_sx2) / 2
    nc_cy = (nc_sy1 + nc_sy2) / 2
    out.append(f'<text x="{nc_cx:.1f}" y="{nc_cy+3:.1f}" text-anchor="middle" font-family="Arial"'
               f' font-size="6" fill="{APPL_STROKE}">COUNTER</text>')
    out.append('</a>')

    return (ww1_cx, ww1_cy, WW_RADIUS), (ww3_cx, ww3_cy, WW_RADIUS)


def _render_furniture(out, data, layout, ww1, ww3):
    """Render furniture: bed, loveseat, ET, chair, ottoman, room labels."""
    pts = data.pts
    to_svg = data.to_svg
    bed = layout.bed
    ww1_cx, ww1_cy, ww1_r = ww1
    ww3_cx, ww3_cy, ww3_r = ww3

    # Bed
    bed_sw = to_svg(bed.w, bed.n)
    bed_se = to_svg(bed.e, bed.s)
    bed_sw_x, bed_sw_y = bed_sw
    bed_se_x, bed_se_y = bed_se
    bed_w = bed_se_x - bed_sw_x
    bed_h = bed_se_y - bed_sw_y
    out.append(f'<rect x="{bed_sw_x:.1f}" y="{bed_sw_y:.1f}" width="{bed_w:.1f}" height="{bed_h:.1f}"'
               f' fill="{APPL_FILL}" stroke="{APPL_STROKE}" stroke-width="{APPL_SW}"/>')
    bed_cx_svg = (bed_sw_x + bed_se_x) / 2
    bed_label_y = bed_sw_y + 0.765 * bed_h
    out.append(f'<text x="{bed_cx_svg:.1f}" y="{bed_label_y+3:.1f}" text-anchor="middle" font-family="Arial"'
               f' font-size="7" fill="{APPL_STROKE}">KING BED</text>')

    # Loveseat: 35" E-W x 65" N-S, rotated 15° CCW about SW corner
    lv_width = LOVESEAT_WIDTH
    lv_height = LOVESEAT_LENGTH
    lv_angle = math.radians(LOVESEAT_ANGLE_DEG)
    K = (math.cos(lv_angle) * (ww3_cy - ww1_cy)
         - math.sin(lv_angle) * (ww3_cx - ww1_cx)
         - ww3_r)
    theta = lv_angle + math.asin(K / ww1_r)
    lv_nw_e = ww1_cx + ww1_r * math.cos(theta)
    lv_nw_n = ww1_cy + ww1_r * math.sin(theta)
    lv_w = lv_nw_e + lv_height * math.sin(lv_angle)
    lv_s = lv_nw_n - lv_height * math.cos(lv_angle)

    # ET position: 2" N of IW1, 2" from loveseat SE corner
    et_r = (ET_RADIUS_CM / 2.54) / 12.0
    lv_se_e = lv_w + lv_width * math.cos(lv_angle)
    lv_se_n = lv_s + lv_width * math.sin(lv_angle)
    et_gap = et_r + STD_GAP
    et_cy = layout.iw1_n + STD_GAP + et_r
    et_cx = lv_se_e + math.sqrt(et_gap**2 - (et_cy - lv_se_n)**2)

    lv_e = lv_w + lv_width
    lv_n = lv_s + lv_height
    lv_sx1, lv_sy1 = to_svg(lv_w, lv_n)
    lv_sx2, lv_sy2 = to_svg(lv_e, lv_s)
    lv_sw = lv_sx2 - lv_sx1; lv_sh = lv_sy2 - lv_sy1
    lv_rot_x = lv_sx1
    lv_rot_y = lv_sy2
    out.append(f'<a href="https://www.ikea.com/us/en/p/saltsjoebaden-loveseat-tonerud-red-brown-s59579188/" target="_blank">')
    out.append(f'<g transform="rotate({int(-LOVESEAT_ANGLE_DEG)},{lv_rot_x:.1f},{lv_rot_y:.1f})">')
    out.append(f'<rect x="{lv_sx1:.1f}" y="{lv_sy1:.1f}" width="{lv_sw:.1f}" height="{lv_sh:.1f}"'
               f' fill="{APPL_FILL}" stroke="{APPL_STROKE}" stroke-width="{APPL_SW}"/>')
    lv_cx = (lv_sx1 + lv_sx2) / 2
    lv_cy = (lv_sy1 + lv_sy2) / 2
    out.append(f'<text x="{lv_cx:.1f}" y="{lv_cy+3:.1f}" text-anchor="middle" font-family="Arial"'
               f' font-size="6" fill="{APPL_STROKE}">LOVESEAT</text>')
    out.append('</g>')
    out.append('</a>')

    # ET: 50cm diameter endtable
    et_sx, et_sy = to_svg(et_cx, et_cy)
    et_r_svg = abs(to_svg(et_r, 0)[0] - to_svg(0, 0)[0])
    out.append('<a href="https://www.ikea.com/us/en/p/listerby-side-table-oak-veneer-30515314/" target="_blank">')
    out.append(f'<circle cx="{et_sx:.1f}" cy="{et_sy:.1f}" r="{et_r_svg:.1f}"'
               f' fill="{APPL_FILL}" stroke="{APPL_STROKE}" stroke-width="{APPL_SW}"/>')
    out.append(f'<text x="{et_sx:.1f}" y="{et_sy+3:.1f}" text-anchor="middle"'
               f' font-family="Arial" font-size="6" fill="{APPL_STROKE}">ET</text>')
    out.append('</a>')

    # LOVESEAT2: same as LOVESEAT but long side E-W (65" E-W x 35" N-S)
    lv2_w = et_cx + et_r + STD_GAP
    lv2_s = layout.iw1_n + STD_GAP
    lv2_e = lv2_w + lv_height  # 65" E-W
    lv2_n = lv2_s + lv_width   # 35" N-S
    lv2_sx1, lv2_sy1 = to_svg(lv2_w, lv2_n)
    lv2_sx2, lv2_sy2 = to_svg(lv2_e, lv2_s)
    lv2_sw = lv2_sx2 - lv2_sx1; lv2_sh = lv2_sy2 - lv2_sy1
    out.append('<a href="https://www.ikea.com/us/en/p/saltsjoebaden-loveseat-tonerud-red-brown-s59579188/" target="_blank">')
    out.append(f'<rect x="{lv2_sx1:.1f}" y="{lv2_sy1:.1f}" width="{lv2_sw:.1f}" height="{lv2_sh:.1f}"'
               f' fill="{APPL_FILL}" stroke="{APPL_STROKE}" stroke-width="{APPL_SW}"/>')
    lv2_cx = (lv2_sx1 + lv2_sx2) / 2
    lv2_cy = (lv2_sy1 + lv2_sy2) / 2
    out.append(f'<text x="{lv2_cx:.1f}" y="{lv2_cy+3:.1f}" text-anchor="middle" font-family="Arial"'
               f' font-size="6" fill="{APPL_STROKE}">LOVESEAT</text>')
    out.append('</a>')

    # CHAIR: 32" E-W x 37" N-S, rounded corners 3", centered between W11 and W12
    ch_angle = math.radians(CHAIR_ANGLE_DEG)
    ch_cx = ((pts["W11"][0] + pts["W12"][0]) / 2
             - 4.0 / 12.0 * math.sin(ch_angle)
             - 1.0 / 12.0)
    ch_cy = ((pts["W11"][1] + pts["W12"][1]) / 2 - 8.0 / 12.0
             - 4.0 / 12.0 * math.cos(ch_angle))
    ch_w = ch_cx - CHAIR_WIDTH / 2
    ch_e = ch_cx + CHAIR_WIDTH / 2
    ch_s = ch_cy - CHAIR_DEPTH / 2
    ch_n = ch_cy + CHAIR_DEPTH / 2
    ch_sx1, ch_sy1 = to_svg(ch_w, ch_n)
    ch_sx2, ch_sy2 = to_svg(ch_e, ch_s)
    ch_sw = ch_sx2 - ch_sx1; ch_sh = ch_sy2 - ch_sy1
    ch_r_svg = abs(to_svg(CHAIR_CORNER_R, 0)[0] - to_svg(0, 0)[0])
    ch_rot_x, ch_rot_y = to_svg(ch_cx, ch_cy)
    out.append(f'<g transform="rotate({int(CHAIR_ANGLE_DEG)},{ch_rot_x:.1f},{ch_rot_y:.1f})">')
    out.append('<a href="https://www.ikea.com/us/en/p/havberg-swivel-easy-chair-and-footstool-grann-bomstad-golden-brown-s59485321/" target="_blank">')
    out.append(f'<rect x="{ch_sx1:.1f}" y="{ch_sy1:.1f}" width="{ch_sw:.1f}" height="{ch_sh:.1f}"'
               f' rx="{ch_r_svg:.1f}" ry="{ch_r_svg:.1f}"'
               f' fill="{APPL_FILL}" stroke="{APPL_STROKE}" stroke-width="{APPL_SW}"/>')
    ch_label_x = (ch_sx1 + ch_sx2) / 2
    ch_label_y = (ch_sy1 + ch_sy2) / 2
    out.append(f'<text x="{ch_label_x:.1f}" y="{ch_label_y+3:.1f}" text-anchor="middle" font-family="Arial"'
               f' font-size="6" fill="{APPL_STROKE}">CHAIR</text>')
    out.append('</a>')
    out.append('</g>')

    # OTTO: 29" x 29", rounded corners 3", 30° CW, 6" SSW of CHAIR
    ot_dist = 39.0 / 12.0  # ch half-h 18.5" + 6" gap + ot half 14.5"
    ot_cx = ch_cx - ot_dist * math.sin(ch_angle)
    ot_cy = ch_cy - ot_dist * math.cos(ch_angle)
    ot_w = ot_cx - OTTOMAN_SIZE / 2
    ot_e = ot_cx + OTTOMAN_SIZE / 2
    ot_s = ot_cy - OTTOMAN_SIZE / 2
    ot_n = ot_cy + OTTOMAN_SIZE / 2
    ot_sx1, ot_sy1 = to_svg(ot_w, ot_n)
    ot_sx2, ot_sy2 = to_svg(ot_e, ot_s)
    ot_sw = ot_sx2 - ot_sx1; ot_sh = ot_sy2 - ot_sy1
    ot_r_svg = ch_r_svg  # same 3" corner radius
    ot_rot_x, ot_rot_y = to_svg(ot_cx, ot_cy)
    out.append(f'<g transform="rotate({int(CHAIR_ANGLE_DEG)},{ot_rot_x:.1f},{ot_rot_y:.1f})">')
    out.append('<a href="https://www.ikea.com/us/en/p/havberg-swivel-easy-chair-and-footstool-grann-bomstad-golden-brown-s59485321/" target="_blank">')
    out.append(f'<rect x="{ot_sx1:.1f}" y="{ot_sy1:.1f}" width="{ot_sw:.1f}" height="{ot_sh:.1f}"'
               f' rx="{ot_r_svg:.1f}" ry="{ot_r_svg:.1f}"'
               f' fill="{APPL_FILL}" stroke="{APPL_STROKE}" stroke-width="{APPL_SW}"/>')
    ot_label_x = (ot_sx1 + ot_sx2) / 2
    ot_label_y = (ot_sy1 + ot_sy2) / 2
    out.append(f'<text x="{ot_label_x:.1f}" y="{ot_label_y+3:.1f}" text-anchor="middle" font-family="Arial"'
               f' font-size="6" fill="{APPL_STROKE}">OTTO</text>')
    out.append('</a>')
    out.append('</g>')

    # Room labels
    bd_cx = layout.bed_cx
    bd_cy = (layout.ctr.s + layout.iw1_s) / 2
    bdx, bdy = to_svg(bd_cx, bd_cy)
    out.append(f'<text x="{bdx:.1f}" y="{bdy+3:.1f}" text-anchor="middle" font-family="Arial"'
               f' font-size="8" fill="#666">BEDROOM</text>')

    of_cx = (layout.iw4_e + pts["W15"][0]) / 2
    of_cy = (layout.cl1_top + layout.iwt3 + layout.iw1_s) / 2 - 2.0
    ofx, ofy = to_svg(of_cx, of_cy)
    out.append(f'<text x="{ofx:.1f}" y="{ofy+3:.1f}" text-anchor="middle" font-family="Arial"'
               f' font-size="8" fill="#666">OFFICE</text>')


def _render_dimensions(out, data, layout):
    """Render all dimension lines (interior and external)."""
    pts = data.pts
    to_svg = data.to_svg

    # Bedroom E-W and N-S
    bd_ew_n = layout.ctr.s + 0.25 * (layout.iw1_s - layout.ctr.s)
    dim_line_h(out, layout.iw3.e, bd_ew_n, layout.iw4_w,
               fmt_dist(layout.iw4_w - layout.iw3.e), to_svg)
    dim_line_v(out, layout.iw3.e + 2.0, layout.ctr.s, layout.iw1_s,
               fmt_dist(layout.iw1_s - layout.ctr.s), to_svg)

    # Closets
    dim_line_v(out, (layout.ctr.e + layout.iwt3 + layout.iw3.w) / 2,
               layout.ctr.s, layout.ctr.n,
               f"CLOSET {fmt_dist(layout.ctr.n - layout.ctr.s)}", to_svg)
    dim_line_v(out, (layout.iw4_e + layout.iw8_w) / 2,
               layout.wall_south_n, layout.cl1_top,
               f"CLOSET {fmt_dist(layout.cl1_top - layout.wall_south_n)}", to_svg)

    # Utility
    dim_line_h(out, pts["W1"][0], (layout.ctr.s + layout.ctr.n) / 2, layout.ctr.e,
               fmt_dist(layout.ctr.e - pts["W1"][0]), to_svg)
    dim_line_h(out, layout.iw8_e, 5.0, pts["W15"][0],
               fmt_dist(pts["W15"][0] - layout.iw8_e), to_svg)

    # Storage
    dim_line_h(out, layout.iw4_e, (layout.iw5.n + layout.iw1_s) / 2, pts["W15"][0],
               f"STORAGE {fmt_dist(pts['W15'][0] - layout.iw4_e)}", to_svg)

    # West wall interior widths
    dim_f1f2_n = layout.ctr.n + layout.iwt3 + 1.0
    dim_line_h(out, pts["W2"][0], dim_f1f2_n, layout.iw3.w,
               fmt_dist(layout.iw3.w - pts["W2"][0]), to_svg)
    dim_line_h(out, pts["W2"][0], pts["F2"][1], layout.iw2.w,
               fmt_dist(layout.iw2.w - pts["W2"][0]), to_svg)
    dim_line_h(out, pts["W5"][0], pts["F5"][1], layout.iw2.w,
               fmt_dist(layout.iw2.w - pts["W5"][0]), to_svg)

    # Office/bedroom verticals
    dim_line_v(out, pts["F18"][0], layout.iw5.s, pts["W18"][1],
               fmt_dist(layout.iw5.s - pts["W18"][1]), to_svg)
    dim_line_v(out, pts["F6"][0] + 1.0, layout.iw6_n, pts["W6"][1],
               fmt_dist(pts["W6"][1] - layout.iw6_n), to_svg)
    dim_line_v(out, pts["F6"][0] + 1.0, layout.iw1_n, layout.iw6_s,
               fmt_dist(layout.iw6_s - layout.iw1_n), to_svg)

    # External dimensions
    dim_ext_e = pts["F2"][0] - 2.7
    dim_line_v(out, dim_ext_e, pts["F0"][1], pts["F6"][1],
               fmt_dist(pts["F6"][1] - pts["F0"][1]), to_svg)

    dim_line_h(out, pts["F8"][0], pts["F6"][1] + 1.0, pts["F11"][0],
               fmt_dist(pts["F11"][0] - pts["F8"][0]), to_svg)

    dim_ext_n = pts["F19"][1] - 3.0
    dim_line_h(out, pts["F1"][0], dim_ext_n, pts["F15"][0],
               fmt_dist(pts["F15"][0] - pts["F1"][0]), to_svg)

    o9_dim_e = (layout.bed.e + layout.iw4_w) / 2
    dim_line_v(out, o9_dim_e, pts["W18"][1], layout.iw1_s,
               fmt_dist(layout.iw1_s - pts["W18"][1]), to_svg)


def _render_openings(out, data, layout):
    """Render O1-O11 opening polygons."""
    pts = data.pts
    to_svg = data.to_svg
    outer_openings = compute_outer_openings(pts, layout)
    for o in outer_openings:
        svg = " ".join(f"{to_svg(*p)[0]:.1f},{to_svg(*p)[1]:.1f}" for p in o.poly)
        out.append(f'<polygon points="{svg}" fill="{OPENING_FILL}" stroke="{OPENING_STROKE}" stroke-width="{WALL_SW}"/>')


def _render_title_block(out, data, inner_area):
    """Render north arrow and title block."""
    # North arrow
    out.append(f'<line x1="{data.na_x:.1f}" y1="{data.na_base_y:.1f}" x2="{data.na_x:.1f}" y2="{data.na_tip_y:.1f}" stroke="#333" stroke-width="2"'
               f' marker-end="url(#ah)"/>')
    out.append(f'<text x="{data.na_x:.1f}" y="{data.na_text_y:.1f}" text-anchor="middle" font-family="Arial"'
               f' font-size="13" font-weight="bold">N</text>')

    # Title block
    out.append(f'<rect x="{data.tb_left:.1f}" y="{data.tb_top:.1f}" width="{data.tb_w}" height="{data.tb_h}"'
               f' fill="white" stroke="#333" stroke-width="1"/>')
    out.append(f'<text x="{data.tb_cx:.1f}" y="{data.tb_top+14:.1f}" text-anchor="middle"'
               f' font-family="Arial" font-size="11" font-weight="bold" fill="#333">'
               f'{inner_area:.2f} sq ft</text>')
    out.append(f'<text x="{data.tb_cx:.1f}" y="{data.tb_top+26:.1f}" text-anchor="middle"'
               f' font-family="Arial" font-size="8" fill="#666">Interior area</text>')
    out.append(f'<text x="{data.tb_cx:.1f}" y="{data.tb_top+40:.1f}" text-anchor="middle"'
               f' font-family="Arial" font-size="11" font-weight="bold" fill="#333">'
               f'{data.outer_area:.2f} sq ft</text>')
    out.append(f'<text x="{data.tb_cx:.1f}" y="{data.tb_top+52:.1f}" text-anchor="middle"'
               f' font-family="Arial" font-size="8" fill="#666">Exterior area</text>')
    ratio = data.ft_per_inch * 12
    scale_label = f'Scale 1:{ratio:.1f} 1&#8243; = {fmt_dist(data.ft_per_inch)}'
    out.append(f'<text x="{data.tb_cx:.1f}" y="{data.tb_top+64:.1f}" text-anchor="middle"'
               f' font-family="Arial" font-size="8" fill="#666">{scale_label}</text>')
    now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    git_desc = git_describe()
    out.append(f'<text x="{data.tb_cx:.1f}" y="{data.tb_top+76:.1f}" text-anchor="middle"'
               f' font-family="Arial" font-size="7.5" fill="#999">Generated {now}</text>')
    out.append(f'<text x="{data.tb_cx:.1f}" y="{data.tb_top+86:.1f}" text-anchor="middle"'
               f' font-family="Arial" font-size="7.5" fill="#999">from {git_desc}</text>')


# ============================================================
# SVG rendering — orchestrator
# ============================================================

def render_floorplan_svg(data):
    """Render the complete floorplan SVG. Returns SVG string."""
    pts = data.pts
    to_svg = data.to_svg
    layout = compute_interior_layout(pts, data.inner_poly)

    out = []
    out.append(f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}"'
               f' viewBox="{data.vb_x:.2f} {data.vb_y:.2f} {data.vb_w:.2f} {data.vb_h:.2f}">')
    out.append(f'<rect x="{data.vb_x:.2f}" y="{data.vb_y:.2f}" width="{data.vb_w:.2f}" height="{data.vb_h:.2f}" fill="white"/>')
    out.append('<defs>')
    out.append('  <marker id="ah" markerWidth="8" markerHeight="6" refX="8" refY="3" orient="auto">'
               '<polygon points="0 0, 8 3, 0 6" fill="#333"/></marker>')
    out.append('</defs>')
    out.append(f'<text x="{data.title_x:.1f}" y="{data.title_y:.1f}" text-anchor="middle" font-family="Arial" font-size="14"'
               f' font-weight="bold">Parent Suite</text>')

    _render_walls(out, data, layout)
    _render_appliances(out, data, layout)
    ww1, ww3 = _render_kitchen(out, data, layout)
    _render_furniture(out, data, layout, ww1, ww3)
    _render_dimensions(out, data, layout)
    _render_openings(out, data, layout)

    inner_area = data.inner_area - compute_iw_area(layout)
    _render_title_block(out, data, inner_area)
    out.append('</svg>')

    return "\n".join(out)


# ============================================================
# Main entry point
# ============================================================

if __name__ == "__main__":
    data = build_floorplan_data()
    svg_content = render_floorplan_svg(data)
    pts = data.pts
    layout = compute_interior_layout(pts, data.inner_poly)
    inner_area = data.inner_area - compute_iw_area(layout)
    outer_area = data.outer_area

    svg_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "floorplan.svg")
    with open(svg_path, "w") as f:
        f.write(svg_content)

    print(f"Floorplan written to {svg_path}")
    print(f"Outer area:    {outer_area:.2f} sq ft")
    print(f"Interior area: {inner_area:.2f} sq ft")
    print(f"Wall area:     {outer_area - inner_area:.2f} sq ft")
    print()
    seen = set()
    for seg in data.outline_segs:
        if seg.start not in seen:
            seen.add(seg.start)
            f_name = seg.start
            w_name = "W" + f_name[1:]
            o = pts[f_name]; w = pts[w_name]
            print(f"  {f_name:<5s} ({o[0]:8.4f}, {o[1]:8.4f})  ->  inner ({w[0]:8.4f}, {w[1]:8.4f})")
