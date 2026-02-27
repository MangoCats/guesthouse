"""Generate floorplan SVG with 8" wall inset from the outline path.

Computes geometry from shared/ and floorplan/ packages.
Outline points F1-F18, inner wall points W1-W18.
"""
import os, sys, math, datetime
from typing import NamedTuple, Any

# Ensure project root is on sys.path for package imports
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

from shared.types import LineSeg, ArcSeg, BBox
from shared.geometry import (
    segment_polyline, path_polygon, poly_area,
    compute_inner_walls, fmt_dist, f8f9_corner_polyline,
    seg_vecs, offset_pt, line_isect, require_pts, GEOM_EPS,
)
from shared.survey import compute_traverse, compute_three_arc, compute_inset
from shared.svg import make_svg_transform, W, H, git_describe
from floorplan.geometry import compute_outline_geometry, align_pts_to_f_series
from floorplan.constants import (
    WALL_OUTER, WALL_3IN, SHELL_THICKNESS, AIR_GAP, OPENING_INSIDE_RADIUS,
    WH_RADIUS,
    SINK_RX, SINK_RY,
    KITCHEN_SINK_WIDTH, KITCHEN_SINK_DEPTH,
    DW_WIDTH, DW_DEPTH, STOVE_WIDTH, STOVE_DEPTH,
    FRIDGE_SIZE, MINIK_FRIDGE_W, MINIK_FRIDGE_D,
    KITCHEN_CTR_LENGTH, KITCHEN_CTR_DEPTH,
    NORTH_CTR_LENGTH, NORTH_CTR_DEPTH,
    JAMB_WIDTH, STD_GAP, KITCHEN_APPL_GAP,
    LOVESEAT_OFFSET_IW4, LOVESEAT_OFFSET_IW1,
    LOVESEAT_WIDTH, LOVESEAT_LENGTH,
    DESK_WIDTH, DESK_DEPTH, DESK_CHAIR_WIDTH, DESK_CHAIR_DEPTH, DESK_CHAIR_GAP,
    CHAIR_WIDTH, CHAIR_DEPTH, CHAIR_CORNER_R,
    OTTOMAN_SIZE, ET_RADIUS_CM,
    SOFA_WIDTH, SOFA_DEPTH,
    ICE_WIDTH, ICE_DEPTH,
    ROCKER_WIDTH, ROCKER_DEPTH, ROCKER_CORNER_R,
    RO1_OFFSET_FROM_IW2, IW1_RO_WIDTH,
    O3_WIDTH, O3_DOOR_WIDTH,
    O6_WIDTH, O6_DOOR_WIDTH, RO1_DOOR_WIDTH, RO2_DOOR_WIDTH,
    RO3_DOOR_WIDTH, RO4_DOOR_WIDTH, RO5_DOOR_WIDTH, RO6_DOOR_WIDTH,
    RO7_DOOR_WIDTH,
    IW4_RO_WIDTH, IW9_RO_WIDTH, IW11_RO_WIDTH, DOOR_FLAT_FACE, F8F9_INNER_TURN_R,
    SHELVES_LENGTH, SHELVES_DEPTH, SHELVES_GAP_IW1, SHELVES_GAP_IW9,
)
from floorplan.layout import compute_interior_layout
from floorplan.openings import (
    compute_outer_openings, compute_rough_openings, outer_to_wall_openings,
)
from shared.wall_shells import (
    compute_inset_path, lerp, openings_on_seg, solid_ranges,
    arc_strip_poly, line_strip_poly, partial_line_strip,
    uturn_polygon, enumerate_wall_sections, build_section_outlines,
)

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

_CASEMENT_NAMES = {"O8", "O9", "O10"}
_CASEMENT_URL = ("https://brogawindows.com/catalog/1894"
                 "?srsltid=AfmBOopybpD3TPjfatlYy33Uy4FHGXQocdwR0ZinpJycfhYyLZdzKXc-SSw")

# ============================================================
# SVG Helpers
# ============================================================

def _svg_angle(along):
    """SVG rotation angle (degrees, CW-positive) for a direction vector.

    Aligns the SVG width axis (positive-X) with the given survey direction.
    """
    return -math.degrees(math.atan2(along[1], along[0]))


def _rotated_dim(out, p1, p2, label, to_svg):
    """Rotated dimension line with tick marks and label.

    Label placement follows drafting convention: above for horizontal lines
    (reading left to right), left for vertical lines (reading bottom to top).
    """
    sx1, sy1 = to_svg(*p1)
    sx2, sy2 = to_svg(*p2)
    sdx = sx2 - sx1; sdy = sy2 - sy1
    slen = math.sqrt(sdx**2 + sdy**2)
    px = -sdy / slen; py = sdx / slen
    tk = 4
    out.append(f'<line x1="{sx1:.1f}" y1="{sy1:.1f}" x2="{sx2:.1f}" y2="{sy2:.1f}" stroke="{DIM_COLOR}" stroke-width="0.8"/>')
    for sx, sy in [(sx1, sy1), (sx2, sy2)]:
        out.append(f'<line x1="{sx - tk * px:.1f}" y1="{sy - tk * py:.1f}" '
                   f'x2="{sx + tk * px:.1f}" y2="{sy + tk * py:.1f}" '
                   f'stroke="{DIM_COLOR}" stroke-width="0.8"/>')
    lmx = (sx1 + sx2) / 2; lmy = (sy1 + sy2) / 2
    # Normalize text angle to [-90°, 90°) for readability
    ang = math.degrees(math.atan2(sdy, sdx))
    if ang >= 90:
        ang -= 180
    elif ang < -90:
        ang += 180
    # Label offset: above for horizontal, left for vertical (drafting convention)
    ang_rad = math.radians(ang)
    lx = lmx + 3 * math.sin(ang_rad)
    ly = lmy - 3 * math.cos(ang_rad)
    out.append(f'<text x="{lx:.1f}" y="{ly:.1f}" text-anchor="middle" font-family="Arial" '
               f'font-size="8" fill="{DIM_COLOR}" transform="rotate({ang:.1f},{lx:.1f},{ly:.1f})">'
               f'{label}</text>')

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


def _appl_poly(out, corners, to_svg, label=None, href=None,
               fill=APPL_FILL, stroke=APPL_STROKE, sw=APPL_SW,
               font_size="7", dash=False, close_href=True, text_rot=None,
               fill_color=None):
    """Render appliance/furniture polygon with optional label and link."""
    _fill = fill_color or fill
    if href:
        out.append(f'<a href="{href}" target="_blank">')
    pts_svg = " ".join(f"{to_svg(*p)[0]:.1f},{to_svg(*p)[1]:.1f}" for p in corners)
    attrs = f' stroke-dasharray="4,3"' if dash else ''
    out.append(f'<polygon points="{pts_svg}" fill="{_fill}" stroke="{stroke}" '
               f'stroke-width="{sw}"{attrs}/>')
    if label:
        cx = sum(p[0] for p in corners) / len(corners)
        cy = sum(p[1] for p in corners) / len(corners)
        scx, scy = to_svg(cx, cy)
        rot = f' transform="rotate({text_rot:.1f},{scx:.1f},{scy+3:.1f})"' if text_rot is not None else ''
        out.append(f'<text x="{scx:.1f}" y="{scy+3:.1f}" text-anchor="middle" '
                   f'font-family="Arial" font-size="{font_size}" fill="{stroke}"{rot}>{label}</text>')
    if href and close_href:
        out.append('</a>')


def _wall_stroke(out, p_start, p_end, half_sw, to_svg):
    """Render a wall face stroke line inset by half_sw from face endpoints."""
    _al, _out = seg_vecs(p_start, p_end)
    _inw = (-_out[0], -_out[1])
    _p1 = offset_pt(p_start, half_sw, _inw)
    _p2 = offset_pt(p_end, half_sw, _inw)
    sx1, sy1 = to_svg(*_p1)
    sx2, sy2 = to_svg(*_p2)
    out.append(f'<line x1="{sx1:.1f}" y1="{sy1:.1f}" x2="{sx2:.1f}" y2="{sy2:.1f}"'
               f' stroke="{WALL_STROKE}" stroke-width="{WALL_SW}"/>')


def _jamb_poly(out, j1, j2, along_dir, to_svg):
    """Render a jamb polygon from edge points j1, j2 + JAMB_WIDTH along along_dir."""
    _jn = (along_dir[0] * JAMB_WIDTH, along_dir[1] * JAMB_WIDTH)
    j_poly = [j1, j2, (j2[0] + _jn[0], j2[1] + _jn[1]),
              (j1[0] + _jn[0], j1[1] + _jn[1])]
    jp = " ".join(f"{to_svg(*p)[0]:.1f},{to_svg(*p)[1]:.1f}" for p in j_poly)
    out.append(f'<polygon points="{jp}" fill="{JAMB_COLOR}" stroke="none"/>')


def _swing_arc_svg(hinge, radius, dir_from, dir_to, to_svg, n_pts=20):
    """SVG polyline points string for a 90-degree door/window swing arc.

    Rotates ``dir_from`` toward ``dir_to`` (must be approximately
    perpendicular unit vectors).  Sweep direction (CW/CCW) is determined
    by the cross product of *dir_from* × *dir_to*.
    """
    cross = dir_from[0] * dir_to[1] - dir_from[1] * dir_to[0]
    sweep = math.pi / 2 if cross > 0 else -math.pi / 2
    pts = []
    for i in range(n_pts + 1):
        a = sweep * i / n_pts
        ca, sa = math.cos(a), math.sin(a)
        ae = hinge[0] + radius * (dir_from[0] * ca - dir_from[1] * sa)
        an = hinge[1] + radius * (dir_from[0] * sa + dir_from[1] * ca)
        sx, sy = to_svg(ae, an)
        pts.append(f"{sx:.1f},{sy:.1f}")
    return " ".join(pts)


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


def draw_toilet(out, center, facing, width, to_svg):
    """Draw a toilet plan view against a wall.

    center: (e, n) point on wall face at toilet center-line
    facing: unit vector from wall toward bowl
    width: unit vector along toilet width (perpendicular to facing)
    """
    pts_survey = [(center[0] + dx * _SVG_TO_FT * width[0] + dy * _SVG_TO_FT * facing[0],
                   center[1] + dx * _SVG_TO_FT * width[1] + dy * _SVG_TO_FT * facing[1])
                  for dx, dy in _TOILET_SVG]
    svg_pts = " ".join(f"{to_svg(e, n)[0]:.1f},{to_svg(e, n)[1]:.1f}"
                       for e, n in pts_survey)
    out.append(f'<polygon points="{svg_pts}"'
               f' fill="{APPL_FILL}" stroke="{APPL_STROKE}" stroke-width="{APPL_SW}"/>')


def draw_sink(out, center_e, center_n, to_svg):
    """Draw a sink plan view as a rectangle (24" Petten console, 23-5/8 x 16-1/2)."""
    sx, sy = to_svg(center_e, center_n)
    # Convert half-dimensions from feet to SVG pixel units
    rx_svg = abs(to_svg(SINK_RX, 0)[0] - to_svg(0, 0)[0])
    ry_svg = abs(to_svg(0, SINK_RY)[1] - to_svg(0, 0)[1])
    x0 = sx - rx_svg
    y0 = sy - ry_svg
    _url = ("https://www.magnushomeproducts.com/products/24-petten-matte-gray"
            "-vitreous-china-console-sink-with-black-powdercoat-steel-stand"
            "-and-shelves")
    out.append(f'<a href="{_url}" target="_blank">')
    out.append(f'<rect x="{x0:.1f}" y="{y0:.1f}" width="{2*rx_svg:.1f}" height="{2*ry_svg:.1f}"'
               f' fill="{APPL_FILL}" stroke="{APPL_STROKE}" stroke-width="{APPL_SW}"/>')
    out.append(f'<text x="{sx:.1f}" y="{sy+3:.1f}" text-anchor="middle" font-family="Arial"'
               f' font-size="7" fill="{APPL_STROKE}">SINK</text>')
    out.append('</a>')


def stroke_segs(out, segs, color, width, pts, to_svg, seg_overrides=None):
    """Render segment strokes (lines and arc polylines).

    seg_overrides: optional dict mapping seg index to replacement polyline
                   (list of (E, N) points) for non-standard segment paths.
    """
    for i, seg in enumerate(segs):
        if seg_overrides and i in seg_overrides:
            poly = seg_overrides[i]
            svg_p = " ".join(f"{to_svg(e,n)[0]:.1f},{to_svg(e,n)[1]:.1f}" for e,n in poly)
            out.append(f'<polyline points="{svg_p}" fill="none" stroke="{color}"'
                       f' stroke-width="{width}" stroke-linecap="round"/>')
        elif isinstance(seg, LineSeg):
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
    w_f8f9_poly: list      # W-series F8-F9 straight-arc-straight polyline
    s_segs: list            # S-series segments (2" inset, inner face of outer shell)
    g_segs: list            # G-series segments (6" inset, outer face of inner shell)
    openings: list          # WallOpening list (parametric outer wall openings)
    g_f8f9_poly: list       # G-series F8-F9 straight-arc-straight polyline
    layout: Any             # InteriorLayout


def build_floorplan_data():
    """Compute all geometry needed for the floorplan SVG."""
    pts = compute_traverse()
    to_svg = make_svg_transform()
    _arc_info = compute_three_arc(pts)
    _inset = compute_inset(pts, _arc_info["R1"], _arc_info["R2"], _arc_info["R3"],
                           _arc_info["nE"], _arc_info["nN"])
    pts.update(_inset.pts_update)
    # Align P/Pi with F-series coordinate space
    align_pts_to_f_series(pts)
    _outline_geo = compute_outline_geometry()
    pts.update(_outline_geo.fp_pts)
    outline_segs = _outline_geo.outline_segs
    _radii = _outline_geo.radii

    wall_t = WALL_OUTER
    inner_segs = compute_inner_walls(outline_segs, pts, wall_t, _radii)
    outer_poly = path_polygon(outline_segs, pts)
    inner_poly = path_polygon(inner_segs, pts)

    # Replace W8-W9 arc in inner_poly with straight-arc-straight path
    require_pts(pts, "W8", "W9")
    w_f8f9_poly = f8f9_corner_polyline(pts, WALL_OUTER, F8F9_INNER_TURN_R)
    w8 = pts["W8"]
    w9 = pts["W9"]
    w8_idx = next(i for i, p in enumerate(inner_poly)
                  if abs(p[0] - w8[0]) < GEOM_EPS and abs(p[1] - w8[1]) < GEOM_EPS)
    w9_idx = next(i for i, p in enumerate(inner_poly)
                  if i > w8_idx
                  and abs(p[0] - w9[0]) < GEOM_EPS and abs(p[1] - w9[1]) < GEOM_EPS)
    inner_poly[w8_idx:w9_idx + 1] = w_f8f9_poly

    # Compute S-series (2" inset = inner face of outer shell)
    s_pts, s_segs = compute_inset_path(outline_segs, pts, _radii,
                                        SHELL_THICKNESS, "S")
    pts.update(s_pts)

    # Compute G-series (6" inset = outer face of inner shell)
    g_pts, g_segs = compute_inset_path(outline_segs, pts, _radii,
                                        SHELL_THICKNESS + AIR_GAP, "G")
    pts.update(g_pts)

    # G-series F8-F9 straight-arc-straight polyline
    g_f8f9_poly = f8f9_corner_polyline(
        pts, SHELL_THICKNESS + AIR_GAP, OPENING_INSIDE_RADIUS)

    # Interior layout and wall openings
    layout = compute_interior_layout(pts, inner_poly)
    outer_openings = compute_outer_openings(pts, layout)
    openings = outer_to_wall_openings(outer_openings, outline_segs, pts)

    outer_area = poly_area(outer_poly)
    inner_area = poly_area(inner_poly)

    # --- Fit content on letter landscape (792x612) page ---
    _margin_top = 36   # 0.5" top margin
    _margin = 72       # 1" margins on left, right, bottom
    _f_names = [f"F{i}" for i in range(19) if i not in (0, 3, 4)]
    _f_svg = [to_svg(*pts[k]) for k in _f_names]
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
    _ext_s_y = to_svg(0, pts["F18"][1] - 3.0)[1]
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
        w_f8f9_poly=w_f8f9_poly,
        s_segs=s_segs, g_segs=g_segs,
        openings=openings,
        g_f8f9_poly=g_f8f9_poly,
        layout=layout,
    )

def compute_placement_points(pts, layout, radii):
    """Compute key appliance/furniture placement points for regression testing.

    Returns list of (name, (E, N)) tuples for the standard (non-minik, non-db)
    floorplan variant.
    """
    result = []
    wall_t = WALL_OUTER

    # ---- Utility room appliances ----
    w2w5_al, w2w5_in = seg_vecs(pts["W2"], pts["W5"])
    _iw2s_e_al, _iw2s_e_out = seg_vecs(layout.iw2s.poly[1], layout.iw2s.poly[2])

    # Water heater center (line-circle intersection)
    _wh_ref = offset_pt(layout.iw2s.poly[2], WH_RADIUS, _iw2s_e_out)
    wh_tangent_r = (radii["R_a7"] - wall_t) - WH_RADIUS
    _c7 = pts["C7"]
    _wh_d = (_wh_ref[0] - _c7[0], _wh_ref[1] - _c7[1])
    _wh_d_al = _wh_d[0] * _iw2s_e_al[0] + _wh_d[1] * _iw2s_e_al[1]
    _wh_d2 = _wh_d[0]**2 + _wh_d[1]**2
    _wh_t = -_wh_d_al + math.sqrt(wh_tangent_r**2 - _wh_d2 + _wh_d_al**2)
    wh_center = offset_pt(_wh_ref, _wh_t, _iw2s_e_al)
    result.append(("wh_center", wh_center))

    # ---- Kitchen appliances ----
    w9w10_al, w9w10_in = seg_vecs(pts["W9"], pts["W10"])
    _wall_anchor = pts["W9"]
    _iw2s_ne = layout.iw2s.poly[2]
    _iw2_d = ((_iw2s_ne[0] - _wall_anchor[0]) * w9w10_al[0] +
              (_iw2s_ne[1] - _wall_anchor[1]) * w9w10_al[1])

    def _nwp(d_along, d_inward=0):
        return offset_pt(offset_pt(_wall_anchor, d_along, w9w10_al),
                         d_inward, w9w10_in)

    _st_d = _iw2_d + NORTH_CTR_LENGTH + KITCHEN_APPL_GAP
    _ks_d = _st_d + STOVE_WIDTH + KITCHEN_APPL_GAP + 2.0 / 12.0
    _dw_d = _ks_d + KITCHEN_SINK_WIDTH + KITCHEN_APPL_GAP

    # SINK NW/SE
    result.append(("sink_NW", _nwp(_ks_d, 0)))
    result.append(("sink_SE", _nwp(_ks_d + KITCHEN_SINK_WIDTH, KITCHEN_SINK_DEPTH)))
    # STOVE NW/SE
    result.append(("stove_NW", _nwp(_st_d, KITCHEN_APPL_GAP)))
    result.append(("stove_SE", _nwp(_st_d + STOVE_WIDTH, KITCHEN_APPL_GAP + STOVE_DEPTH)))
    # D/W NW/SE
    result.append(("dw_NW", _nwp(_dw_d, 0)))
    result.append(("dw_SE", _nwp(_dw_d + DW_WIDTH, DW_DEPTH)))

    # Fridge NW/SE (standard variant: IW1/IW2 corner)
    _iw1_n_al, _iw1_n_cw = seg_vecs(layout.iw1.poly[3], layout.iw1.poly[2])
    _iw1_n_out = (-_iw1_n_cw[0], -_iw1_n_cw[1])
    _iw2_lo_e_al, _ = seg_vecs(layout.iw2.poly[1], layout.iw2.poly[2])
    _iw12_corner = line_isect(layout.iw2.poly[1], _iw2_lo_e_al,
                              layout.iw1.poly[3], _iw1_n_al)
    def _iwp(d_e, d_n=0):
        return offset_pt(offset_pt(_iw12_corner, d_e, _iw1_n_al),
                         d_n, _iw1_n_out)
    _fr_w2 = 32.75 / 12.0
    _fr_h1 = 35.0 / 12.0
    result.append(("fridge_NW", _iwp(KITCHEN_APPL_GAP, KITCHEN_APPL_GAP + _fr_h1)))
    result.append(("fridge_SE", _iwp(KITCHEN_APPL_GAP + _fr_w2, KITCHEN_APPL_GAP)))

    # ICE NW/SE (standard variant)
    _ice_d = _dw_d + DW_WIDTH + 6.0 / 12.0
    _ice_i = KITCHEN_APPL_GAP
    result.append(("ice_NW", _nwp(_ice_d, _ice_i)))
    result.append(("ice_SE", _nwp(_ice_d + ICE_WIDTH, _ice_i + ICE_DEPTH)))

    # North counter NW/SE (standard variant)
    result.append(("nctr_NW", _nwp(_iw2_d, 0)))
    result.append(("nctr_SE", _nwp(_iw2_d + NORTH_CTR_LENGTH, NORTH_CTR_DEPTH)))

    # Work counter NW/SE
    _wc_d2 = KITCHEN_APPL_GAP + _fr_w2 + KITCHEN_APPL_GAP
    result.append(("wctr_NW", _iwp(_wc_d2, 18.0 / 12.0)))
    result.append(("wctr_SE", _iwp(_wc_d2 + 60.0 / 12.0, 0)))

    # ---- Living room furniture ----
    _iw4_w_al, _ = seg_vecs(layout.iw4.poly[3], layout.iw4.poly[0])
    _iw41_corner = line_isect(layout.iw4.poly[3], _iw4_w_al,
                              layout.iw1.poly[3], _iw1_n_al)
    _iw1_w = (-_iw1_n_al[0], -_iw1_n_al[1])  # westward along IW1
    def _lwp(d_w=0, d_n=0):
        return offset_pt(offset_pt(_iw41_corner, d_w, _iw1_w),
                         d_n, _iw1_n_out)

    w12w13_al, _ = seg_vecs(pts["W12"], pts["W13"])
    w11w12_al, w11w12_in = seg_vecs(pts["W11"], pts["W12"])

    # Loveseat NW corner (rotation anchor)
    lv_nw = _lwp(LOVESEAT_OFFSET_IW4, LOVESEAT_OFFSET_IW1)
    result.append(("loveseat_NW", lv_nw))

    # ET center (from loveseat SE corner geometry)
    lv_angle = math.atan2(w12w13_al[0], -w12w13_al[1])
    lv_w = lv_nw[0] + LOVESEAT_LENGTH * math.sin(lv_angle)
    lv_s = lv_nw[1] - LOVESEAT_LENGTH * math.cos(lv_angle)
    lv_se_e = lv_w + LOVESEAT_WIDTH * math.cos(lv_angle)
    lv_se_n = lv_s + LOVESEAT_WIDTH * math.sin(lv_angle)
    et_r = (ET_RADIUS_CM / 2.54) / 12.0
    et_gap = et_r + STD_GAP
    # ET center: on IW1 offset line (STD_GAP + et_r out from IW1 north face),
    # at distance et_gap from loveseat SE corner.  Circle-line intersection
    # along IW1 direction (rotation-safe).
    _et_line_pt = offset_pt(layout.iw1.poly[3], STD_GAP + et_r, _iw1_n_out)
    _et_dx = _et_line_pt[0] - lv_se_e
    _et_dy = _et_line_pt[1] - lv_se_n
    _et_b = 2 * (_et_dx * _iw1_n_al[0] + _et_dy * _iw1_n_al[1])
    _et_c = _et_dx**2 + _et_dy**2 - et_gap**2
    _et_t = (-_et_b + math.sqrt(_et_b**2 - 4 * _et_c)) / 2
    et_center = offset_pt(_et_line_pt, _et_t, _iw1_n_al)
    result.append(("et_center", et_center))

    # Loveseat2 NW/SE: placed along IW1, starting at ET east edge + gap
    _lv2_sw = offset_pt(
        offset_pt(et_center, et_r + STD_GAP, _iw1_n_al),
        -et_r, _iw1_n_out)
    _lv2_nw = offset_pt(_lv2_sw, LOVESEAT_WIDTH, _iw1_n_out)
    _lv2_se = offset_pt(_lv2_sw, LOVESEAT_LENGTH, _iw1_n_al)
    result.append(("loveseat2_NW", _lv2_nw))
    result.append(("loveseat2_SE", _lv2_se))

    # Chair center
    _ch_svg_deg = _svg_angle(w12w13_al) - 45
    ch_angle = math.radians(_ch_svg_deg)
    _ch_mid = ((pts["W11"][0] + pts["W12"][0]) / 2,
               (pts["W11"][1] + pts["W12"][1]) / 2)
    _ch_base = offset_pt(offset_pt(_ch_mid, -1.0 / 12.0, w11w12_al),
                         8.0 / 12.0, w11w12_in)
    ch_cx = _ch_base[0] - 4.0 / 12.0 * math.sin(ch_angle)
    ch_cy = _ch_base[1] - 4.0 / 12.0 * math.cos(ch_angle)
    result.append(("chair_center", (ch_cx, ch_cy)))

    # Ottoman center
    ot_dist = 39.0 / 12.0
    ot_cx = ch_cx - ot_dist * math.sin(ch_angle)
    ot_cy = ch_cy - ot_dist * math.cos(ch_angle)
    result.append(("ottoman_center", (ot_cx, ot_cy)))

    # Desk SW corner (rotation anchor at W17)
    result.append(("desk_SW", pts["W17"]))

    return result


def compute_dimension_endpoints(pts, layout, radii, bare=False):
    """Compute dimension line endpoints using wall-relative operations.

    Returns list of (name, (E, N)) tuples — two endpoints per dimension line.
    All positions derived from wall polygon faces, segment directions, and
    line intersections (never from BBox fields or raw coordinate indexing).
    """
    result = []
    openings = compute_outer_openings(pts, layout)
    ro_by_name = {r.name: r for r in compute_rough_openings(pts, layout)}

    # Reference directions from axis-aligned wall segments
    _ew, _ = seg_vecs(pts["W9"], pts["W10"])      # east direction
    _ns, _ = seg_vecs(pts["W2"], pts["W5"])        # north direction
    _w18w1_al, _w18w1_in = seg_vecs(pts["W18"], pts["W1"])

    # Pre-compute commonly used IW face directions
    _iw1_n_al, _ = seg_vecs(layout.iw1.poly[3], layout.iw1.poly[2])
    _iw1_s_al, _ = seg_vecs(layout.iw1.poly[0], layout.iw1.poly[1])
    _iw2_e_al, _ = seg_vecs(layout.iw2.poly[1], layout.iw2.poly[2])
    _iw2_w_al, _ = seg_vecs(layout.iw2.poly[0], layout.iw2.poly[3])
    _w14w15_al, _w14w15_in = seg_vecs(pts["W14"], pts["W15"])

    # ---- dim01: IW1-north → W9 (vertical) ----
    _f9f11_mid = ((pts["F9"][0] + pts["F11"][0]) / 2,
                  (pts["F9"][1] + pts["F11"][1]) / 2)
    _dim01_ref = offset_pt(_f9f11_mid, 1.0, _ew)
    result.append(("dim01_A", line_isect(layout.iw1.poly[3], _iw1_n_al,
                                         _dim01_ref, _ns)))
    result.append(("dim01_B", line_isect(pts["W9"], _ew,
                                         _dim01_ref, _ns)))

    # ---- dim02: IW2-east → W12-W13 (horizontal) ----
    _f12f13_mid = ((pts["F12"][0] + pts["F13"][0]) / 2,
                   (pts["F12"][1] + pts["F13"][1]) / 2)
    _w13w12_al, _ = seg_vecs(pts["W13"], pts["W12"])
    result.append(("dim02_A", line_isect(layout.iw2.poly[1], _iw2_e_al,
                                         _f12f13_mid, _ew)))
    result.append(("dim02_B", line_isect(pts["W13"], _w13w12_al,
                                         _f12f13_mid, _ew)))

    # ---- dim03: East closet (rotated, IW12 south center → W18-W1) ----
    _iw12_s_mid = ((layout.iw12.poly[0][0] + layout.iw12.poly[1][0]) / 2,
                   (layout.iw12.poly[0][1] + layout.iw12.poly[1][1]) / 2)
    _iw11_e_al, _ = seg_vecs(layout.iw11.poly[1], layout.iw11.poly[2])
    result.append(("dim03_A", _iw12_s_mid))
    result.append(("dim03_B", line_isect(_iw12_s_mid, _iw11_e_al,
                                         pts["W18"], _w18w1_al)))

    # ---- dim04: West closet (rotated, IW7 south center → W18-W1) ----
    _iw7_s_mid = ((layout.iw7.poly[0][0] + layout.iw7.poly[1][0]) / 2,
                  (layout.iw7.poly[0][1] + layout.iw7.poly[1][1]) / 2)
    _iw7_w_al, _ = seg_vecs(layout.iw7.poly[0], layout.iw7.poly[3])
    result.append(("dim04_A", _iw7_s_mid))
    result.append(("dim04_B", line_isect(_iw7_s_mid, _iw7_w_al,
                                         pts["W18"], _w18w1_al)))

    # ---- dim05: W2-W5 → IW3 west face, parallel to W18-W1 at IW3 midpoint ----
    _iw3_w_mid = ((layout.iw3.poly[0][0] + layout.iw3.poly[3][0]) / 2,
                   (layout.iw3.poly[0][1] + layout.iw3.poly[3][1]) / 2)
    _dim05_B = _iw3_w_mid
    _dim05_A = line_isect(pts["W2"], _ns, _dim05_B, _w18w1_al)
    result.extend([("dim05_A", _dim05_A), ("dim05_B", _dim05_B)])

    # ---- dim06: IW4-east → W14-W15 (perp to W14-W15 at O8 center) ----
    _iw4_e_al, _ = seg_vecs(layout.iw4.poly[1], layout.iw4.poly[2])
    _o8 = [o for o in openings if o.name == "O8"][0]
    _o8_mid = ((_o8.poly[2][0] + _o8.poly[3][0]) / 2,
               (_o8.poly[2][1] + _o8.poly[3][1]) / 2)
    result.append(("dim06_A", line_isect(layout.iw4.poly[1], _iw4_e_al,
                                         _o8_mid, _w14w15_in)))
    result.append(("dim06_B", _o8_mid))

    # ---- dim07: Storage — IW11-east → W14-W15 (horizontal) ----
    _iw11_e_al, _ = seg_vecs(layout.iw11.poly[1], layout.iw11.poly[2])
    _iw5_n_mid = ((layout.iw5.poly[3][0] + layout.iw5.poly[2][0]) / 2,
                  (layout.iw5.poly[3][1] + layout.iw5.poly[2][1]) / 2)
    _iw1_s_mid = ((layout.iw1.poly[0][0] + layout.iw1.poly[1][0]) / 2,
                  (layout.iw1.poly[0][1] + layout.iw1.poly[1][1]) / 2)
    _stor_ref = ((_iw5_n_mid[0] + _iw1_s_mid[0]) / 2,
                 (_iw5_n_mid[1] + _iw1_s_mid[1]) / 2)
    result.append(("dim07_A", line_isect(layout.iw11.poly[1], _iw11_e_al,
                                         _stor_ref, _ew)))
    result.append(("dim07_B", line_isect(pts["W14"], _w14w15_al,
                                         _stor_ref, _ew)))

    # ---- dim08: O1 east center → RO3 west center (horizontal) ----
    _o1 = openings[0]
    _o1_e_ctr = ((_o1.poly[2][0] + _o1.poly[3][0]) / 2,
                 (_o1.poly[2][1] + _o1.poly[3][1]) / 2)
    _ro3 = ro_by_name["RO3"]
    # RO3 is now in IW9 (rotated): poly[1]+poly[2] = IW9 west face (toward IW3)
    _ro3_w_ctr = ((_ro3.poly[1][0] + _ro3.poly[2][0]) / 2,
                  (_ro3.poly[1][1] + _ro3.poly[2][1]) / 2)
    result.extend([("dim08_A", _o1_e_ctr), ("dim08_B", _ro3_w_ctr)])

    # ---- dim09: W2-W5 → IW2-west (horizontal at O2 center northing) ----
    _o2 = [o for o in openings if o.name == "O2"][0]
    _o2_ctr = (sum(p[0] for p in _o2.poly) / 4,
               sum(p[1] for p in _o2.poly) / 4)
    result.append(("dim09_A", line_isect(pts["W2"], _ns,
                                         _o2_ctr, _ew)))
    result.append(("dim09_B", line_isect(layout.iw2.poly[0], _iw2_w_al,
                                         _o2_ctr, _ew)))

    # ---- dim10: W2-W5 → IW2s-west (horizontal at F5 northing) ----
    _iw2s_w_al, _ = seg_vecs(layout.iw2s.poly[0], layout.iw2s.poly[3])
    result.append(("dim10_A", line_isect(pts["W2"], _ns,
                                         pts["F5"], _ew)))
    result.append(("dim10_B", line_isect(layout.iw2s.poly[0], _iw2s_w_al,
                                         pts["F5"], _ew)))

    # ---- dim11: IW5-south → W18 (vertical at F18 easting) ----
    _iw5_s_al, _ = seg_vecs(layout.iw5.poly[0], layout.iw5.poly[1])
    result.append(("dim11_A", line_isect(layout.iw5.poly[0], _iw5_s_al,
                                         pts["F18"], _ns)))
    result.append(("dim11_B", line_isect(pts["W18"], _ew,
                                         pts["F18"], _ns)))

    # ---- dim12: Office verticals ----
    if not bare:
        _dim12_ref = offset_pt(pts["F6"], 1.0, _ew)
        _iw6_n_al, _ = seg_vecs(layout.iw6.poly[3], layout.iw6.poly[2])
        _iw6_s_al, _ = seg_vecs(layout.iw6.poly[0], layout.iw6.poly[1])
        _iw8_n_al, _ = seg_vecs(layout.iw8.poly[3], layout.iw8.poly[2])
        # dim12a: IW6-north → W6
        result.append(("dim12a_A", line_isect(layout.iw6.poly[3], _iw6_n_al,
                                              _dim12_ref, _ns)))
        result.append(("dim12a_B", line_isect(pts["W6"], _ew,
                                              _dim12_ref, _ns)))
        # dim12b: IW8-north → IW6-south
        result.append(("dim12b_A", line_isect(layout.iw8.poly[3], _iw8_n_al,
                                              _dim12_ref, _ns)))
        result.append(("dim12b_B", line_isect(layout.iw6.poly[0], _iw6_s_al,
                                              _dim12_ref, _ns)))
    else:
        # dim12_bare: IW8-north → W6-W7, centered in O4
        _o4 = [o for o in openings if o.name == "O4"][0]
        _o4_ctr = (sum(p[0] for p in _o4.poly) / 4,
                   sum(p[1] for p in _o4.poly) / 4)
        _iw8_n_al, _ = seg_vecs(layout.iw8.poly[3], layout.iw8.poly[2])
        result.append(("dim12bare_A", line_isect(layout.iw8.poly[3], _iw8_n_al,
                                                 _o4_ctr, _ns)))
        result.append(("dim12bare_B", line_isect(pts["W6"], _ew,
                                                 _o4_ctr, _ns)))

    # ---- dim13: External F18 → F6 (vertical) ----
    _dim13_ref = offset_pt(pts["F2"], -2.7, _ew)
    result.append(("dim13_A", line_isect(pts["F18"], _ew, _dim13_ref, _ns)))
    result.append(("dim13_B", line_isect(pts["F6"], _ew, _dim13_ref, _ns)))

    # ---- dim14: Arc width F7-F8 to F11-F11a (horizontal) ----
    _dim14_arc_ref = offset_pt(pts["F7"], -4.0 / 12.0, _ns)
    _c7 = pts["C7"]; _r7 = radii["R_a7"]
    _c7_dn = ((_dim14_arc_ref[0] - _c7[0]) * _ns[0]
              + (_dim14_arc_ref[1] - _c7[1]) * _ns[1])
    _c7_de = math.sqrt(_r7**2 - _c7_dn**2)
    _arc7_pt = offset_pt(offset_pt(_c7, _c7_dn, _ns), _c7_de, _ew)
    _c11a = pts["C11a"]; _r11 = radii["R_a11"]
    _c11_dn = ((_dim14_arc_ref[0] - _c11a[0]) * _ns[0]
               + (_dim14_arc_ref[1] - _c11a[1]) * _ns[1])
    _c11_de = math.sqrt(_r11**2 - _c11_dn**2)
    _arc11_pt = offset_pt(offset_pt(_c11a, _c11_dn, _ns), -_c11_de, _ew)
    _dim14_render = offset_pt(pts["F6"], 1.0, _ns)
    result.append(("dim14_A", line_isect(_arc7_pt, _ns, _dim14_render, _ew)))
    result.append(("dim14_B", line_isect(_arc11_pt, _ns, _dim14_render, _ew)))

    # ---- dim15: External F2 → F15 (horizontal) ----
    _dim15_ref = offset_pt(pts["F18"], -3.0, _ns)
    result.append(("dim15_A", line_isect(pts["F2"], _ns, _dim15_ref, _ew)))
    result.append(("dim15_B", line_isect(pts["F15"], _ns, _dim15_ref, _ew)))

    # ---- dim17: O10 inner → IW1-south (rotated, perp to W18-W1) ----
    _o10 = openings[9]
    _o10_ic = ((_o10.poly[2][0] + _o10.poly[3][0]) / 2,
               (_o10.poly[2][1] + _o10.poly[3][1]) / 2)
    result.append(("dim17_A", _o10_ic))
    result.append(("dim17_B", line_isect(_o10_ic, _w18w1_in,
                                         layout.iw1.poly[0], _iw1_s_al)))

    # ---- dim18: IW9-east → IW11-west (rotated) ----
    _iw9_e_mid = ((layout.iw9.poly[1][0] + layout.iw9.poly[2][0]) / 2,
                  (layout.iw9.poly[1][1] + layout.iw9.poly[2][1]) / 2)
    _, _iw9_e_in = seg_vecs(layout.iw9.poly[1], layout.iw9.poly[2])
    _iw11_w_al, _ = seg_vecs(layout.iw11.poly[0], layout.iw11.poly[3])
    result.append(("dim18_A", _iw9_e_mid))
    result.append(("dim18_B", line_isect(_iw9_e_mid, _iw9_e_in,
                                         layout.iw11.poly[0], _iw11_w_al)))

    # ---- dim19: O11 inner → IW8-south (vertical) ----
    _o11 = openings[10]
    _o11_ic = ((_o11.poly[2][0] + _o11.poly[3][0]) / 2,
               (_o11.poly[2][1] + _o11.poly[3][1]) / 2)
    _iw8_s_al, _ = seg_vecs(layout.iw8.poly[0], layout.iw8.poly[1])
    result.append(("dim19_A", _o11_ic))
    result.append(("dim19_B", line_isect(layout.iw8.poly[0], _iw8_s_al,
                                         _o11_ic, _ns)))

    # ---- dim22: IW12-north-mid → IW5-south (rotated, perp to W18-W1) ----
    _iw12_n_mid = ((layout.iw12.poly[2][0] + layout.iw12.poly[3][0]) / 2,
                   (layout.iw12.poly[2][1] + layout.iw12.poly[3][1]) / 2)
    result.append(("dim22_A", _iw12_n_mid))
    result.append(("dim22_B", line_isect(_iw12_n_mid, _w18w1_in,
                                         layout.iw5.poly[0], _iw5_s_al)))

    # ---- Bare-only dimensions ----
    if bare:
        # dim20: IW2-east → W14 (horizontal at W14 northing)
        result.append(("dim20_A", line_isect(layout.iw2.poly[1], _iw2_e_al,
                                             pts["W14"], _ew)))
        result.append(("dim20_B", pts["W14"]))

        # dim21: W11a-W11b mid → IW1-north (vertical)
        _w11_mid = ((pts["W11a"][0] + pts["W11b"][0]) / 2,
                    (pts["W11a"][1] + pts["W11b"][1]) / 2)
        result.append(("dim21_A", line_isect(layout.iw1.poly[3], _iw1_n_al,
                                             _w11_mid, _ns)))
        result.append(("dim21_B", _w11_mid))

    return result


# ============================================================
# Render sub-functions
# ============================================================

def compute_iw_area(layout):
    """Compute total interior wall area from layout polygons."""
    iw_polys = [layout.iw1.poly, layout.iw8.poly,
                layout.iw2.poly, layout.iw2o.poly, layout.iw2s.poly,
                layout.iw3.poly, layout.iw7.poly, layout.iw9.poly, layout.iw6.poly,
                layout.iw4.poly, layout.iw11.poly, layout.iw12.poly,
                layout.iw5.poly]
    return sum(poly_area(p) for p in iw_polys)


def _svg_wall_poly(out, poly, to_svg):
    """Render a shell polygon with floorplan wall styling (no stroke, gray fill)."""
    svg = " ".join(f"{to_svg(*p)[0]:.1f},{to_svg(*p)[1]:.1f}" for p in poly)
    out.append(f'<polygon points="{svg}" fill="{WALL_FILL}" stroke="none"/>')


def _render_walls(out, data, layout, bare=False):
    """Render outer wall fill with double-shell detail, outline strokes, and interior walls."""
    pts = data.pts
    to_svg = data.to_svg
    outline_segs = data.outline_segs
    inner_segs = data.inner_segs
    s_segs = data.s_segs
    g_segs = data.g_segs
    openings = data.openings

    shell_t = SHELL_THICKNESS
    R_in = OPENING_INSIDE_RADIUS
    R_out = R_in + shell_t

    # --- Per-segment shell strips ---
    for seg_idx in range(len(outline_segs)):
        seg = outline_segs[seg_idx]
        inner_seg = inner_segs[seg_idx]
        s_seg = s_segs[seg_idx]
        g_seg = g_segs[seg_idx]

        seg_ops = openings_on_seg(openings, seg_idx)

        if isinstance(seg, ArcSeg):
            # Arc segments: full shell strips (no openings on arcs)
            outer_shell = arc_strip_poly(seg, pts, "F", s_seg)
            _svg_wall_poly(out, outer_shell, to_svg)

            if seg_idx == 5:
                inner_shell = (list(data.g_f8f9_poly)
                               + list(reversed(data.w_f8f9_poly)))
            else:
                inner_shell = arc_strip_poly(g_seg, pts, "G", inner_seg)
            _svg_wall_poly(out, inner_shell, to_svg)

        elif isinstance(seg, LineSeg):
            if not seg_ops:
                # Full rectangle strips
                outer_strip = line_strip_poly(pts, seg.start, seg.end,
                                              s_seg.start, s_seg.end)
                _svg_wall_poly(out, outer_strip, to_svg)

                inner_strip = line_strip_poly(pts, g_seg.start, g_seg.end,
                                              inner_seg.start, inner_seg.end)
                _svg_wall_poly(out, inner_strip, to_svg)
            else:
                # Segments with openings: partial strips + U-turns
                sr = solid_ranges(seg_ops)

                # Trim ranges where U-turn arcs will be
                F_A, F_B = pts[seg.start], pts[seg.end]
                seg_len = math.sqrt((F_B[0]-F_A[0])**2 +
                                    (F_B[1]-F_A[1])**2)
                delta_t = R_out / seg_len
                adjusted = []
                for t_s, t_e in sr:
                    if t_s > 1e-9:
                        t_s += delta_t
                    if t_e < 1.0 - 1e-9:
                        t_e -= delta_t
                    if t_e > t_s + 1e-9:
                        adjusted.append((t_s, t_e))

                for t_s, t_e in adjusted:
                    outer_strip = partial_line_strip(
                        pts, seg, s_seg, t_s, t_e)
                    _svg_wall_poly(out, outer_strip, to_svg)

                    inner_strip = partial_line_strip(
                        pts, g_seg, inner_seg, t_s, t_e)
                    _svg_wall_poly(out, inner_strip, to_svg)

                # U-turns at opening boundaries
                for op in seg_ops:
                    uturn_start = uturn_polygon(
                        pts, outline_segs, inner_segs, s_segs, g_segs,
                        seg_idx, op.t_start, "start", shell_t, R_in, WALL_OUTER)
                    _svg_wall_poly(out, uturn_start, to_svg)

                    uturn_end = uturn_polygon(
                        pts, outline_segs, inner_segs, s_segs, g_segs,
                        seg_idx, op.t_end, "end", shell_t, R_in, WALL_OUTER)
                    _svg_wall_poly(out, uturn_end, to_svg)

                # Opening void polygons (4" wide, centered on wall)
                # Skip doors (O3, O6) — only render window voids
                for op in seg_ops:
                    if op.name in ("O3", "O6"):
                        continue
                    S_A, S_B = pts[s_seg.start], pts[s_seg.end]
                    G_A, G_B = pts[g_seg.start], pts[g_seg.end]
                    o_poly = [
                        lerp(S_A, S_B, op.t_start),
                        lerp(S_A, S_B, op.t_end),
                        lerp(G_A, G_B, op.t_end),
                        lerp(G_A, G_B, op.t_start),
                    ]
                    svg = " ".join(f"{to_svg(*p)[0]:.1f},{to_svg(*p)[1]:.1f}"
                                   for p in o_poly)
                    _poly_el = (f'<polygon points="{svg}" fill="{OPENING_FILL}" '
                                f'stroke="{OPENING_STROKE}" stroke-width="{WALL_SW}"/>')
                    if op.name in _CASEMENT_NAMES:
                        out.append(f'<a href="{_CASEMENT_URL}">{_poly_el}</a>')
                    else:
                        out.append(_poly_el)

    # --- Continuous section outlines ---
    g_overrides = {5: data.g_f8f9_poly}
    w_overrides = {5: data.w_f8f9_poly}
    sections = enumerate_wall_sections(openings, outline_segs)
    for start_op, end_op in sections:
        outer_path, cavity_path = build_section_outlines(
            pts, outline_segs, inner_segs, s_segs, g_segs,
            start_op, end_op, shell_t, R_in, WALL_OUTER,
            g_seg_overrides=g_overrides, w_seg_overrides=w_overrides)
        for path in [outer_path, cavity_path]:
            svg_pts = " ".join(
                f"{to_svg(*p)[0]:.1f},{to_svg(*p)[1]:.1f}" for p in path)
            out.append(f'<polygon points="{svg_pts}" fill="none" '
                       f'stroke="#999" stroke-width="0.3"/>')

    # Half stroke width in survey feet (for inside-only edge lines)
    svg_per_ft = abs(to_svg(1, 0)[0] - to_svg(0, 0)[0])
    half_sw = 0.5 / svg_per_ft

    # Rough openings
    rough_openings = compute_rough_openings(pts, layout)

    # ---- IW1 with RO1 ----
    iw_sw, iw_se, iw_ne, iw_nw = layout.iw1.poly
    _ro1p = [r for r in rough_openings if r.name == "RO1"][0].poly
    iw1_w_poly = [iw_sw, _ro1p[0], _ro1p[3], iw_nw]
    iw1_e_poly = [_ro1p[1], iw_se, iw_ne, _ro1p[2]]
    wall_poly(out, iw1_w_poly, to_svg, stroke=False)
    wall_poly(out, iw1_e_poly, to_svg, stroke=False)
    _wall_stroke(out, iw_sw, iw_nw, half_sw, to_svg)        # west end
    _wall_stroke(out, iw_sw, _ro1p[0], half_sw, to_svg)     # south face west
    _wall_stroke(out, _ro1p[1], iw_se, half_sw, to_svg)     # south face east
    _wall_stroke(out, iw_nw, _ro1p[3], half_sw, to_svg)     # north face west
    _wall_stroke(out, _ro1p[2], iw_ne, half_sw, to_svg)     # north face east
    _iw1_al, _ = seg_vecs(iw_sw, iw_se)
    _jamb_poly(out, _ro1p[3], _ro1p[0], _iw1_al, to_svg)
    _neg_iw1_al = (-_iw1_al[0], -_iw1_al[1])
    _jamb_poly(out, _ro1p[1], _ro1p[2], _neg_iw1_al, to_svg)

    # ---- IW8 (no openings) ----
    wall_poly(out, layout.iw8.poly, to_svg, stroke=False)
    _wall_stroke(out, layout.iw8.poly[0], layout.iw8.poly[1], half_sw, to_svg)
    _wall_stroke(out, layout.iw8.poly[3], layout.iw8.poly[2], half_sw, to_svg)

    # ---- IW2 (lower, solid, no opening) ----
    wall_poly(out, layout.iw2.poly, to_svg, stroke=False)
    _wall_stroke(out, layout.iw2.poly[0], layout.iw2.poly[1], half_sw, to_svg)
    _wall_stroke(out, layout.iw2.poly[3], layout.iw2.poly[2], half_sw, to_svg)
    _wall_stroke(out, layout.iw2.poly[1], layout.iw2.poly[2], half_sw, to_svg)
    _wall_stroke(out, layout.iw2.poly[0], layout.iw2.poly[3], half_sw, to_svg)

    # ---- IW2s (shower, solid, no opening) ----
    wall_poly(out, layout.iw2s.poly, to_svg, stroke=False)
    _wall_stroke(out, layout.iw2s.poly[0], layout.iw2s.poly[1], half_sw, to_svg)
    _wall_stroke(out, layout.iw2s.poly[3], layout.iw2s.poly[2], half_sw, to_svg)
    _wall_stroke(out, layout.iw2s.poly[1], layout.iw2s.poly[2], half_sw, to_svg)
    _wall_stroke(out, layout.iw2s.poly[0], layout.iw2s.poly[3], half_sw, to_svg)

    # ---- IW2o (oblique) with RO4 ----
    _iw2o = layout.iw2o.poly
    _ro4p = [r for r in rough_openings if r.name == "RO4"][0].poly
    iw2o_s_poly = [_iw2o[0], _iw2o[1], _ro4p[1], _ro4p[0]]
    iw2o_n_poly = [_ro4p[3], _ro4p[2], _iw2o[2], _iw2o[3]]
    wall_poly(out, iw2o_s_poly, to_svg, stroke=False)
    wall_poly(out, iw2o_n_poly, to_svg, stroke=False)
    _wall_stroke(out, _iw2o[3], _ro4p[3], half_sw, to_svg)    # west face north
    _wall_stroke(out, _ro4p[0], _iw2o[0], half_sw, to_svg)    # west face south
    _wall_stroke(out, _iw2o[1], _ro4p[1], half_sw, to_svg)    # east face south
    _wall_stroke(out, _ro4p[2], _iw2o[2], half_sw, to_svg)    # east face north
    _iw2o_al, _ = seg_vecs(_iw2o[0], _iw2o[3])
    _jamb_poly(out, _ro4p[0], _ro4p[1], _iw2o_al, to_svg)
    _neg_iw2o_al = (-_iw2o_al[0], -_iw2o_al[1])
    _jamb_poly(out, _ro4p[2], _ro4p[3], _neg_iw2o_al, to_svg)

    # ---- IW3 (solid, no opening, rotated perpendicular to W18-W1) ----
    _iw3_sw, _iw3_se, _iw3_ne, _iw3_nw = layout.iw3.poly
    wall_poly(out, layout.iw3.poly, to_svg, stroke=False)
    _iw3_dx_t = _iw3_sw[0] - _iw3_se[0]
    _iw3_dy_t = _iw3_sw[1] - _iw3_se[1]
    _iw3_lt = math.sqrt(_iw3_dx_t**2 + _iw3_dy_t**2)
    _iw3_at = (_iw3_dx_t / _iw3_lt, _iw3_dy_t / _iw3_lt)
    for (p1, p2), (ox, oy) in [
        ((_iw3_se, _iw3_ne), _iw3_at),                            # east face
        ((_iw3_sw, _iw3_nw), (-_iw3_at[0], -_iw3_at[1])),        # west face
    ]:
        sx1, sy1 = to_svg(p1[0] + half_sw * ox, p1[1] + half_sw * oy)
        sx2, sy2 = to_svg(p2[0] + half_sw * ox, p2[1] + half_sw * oy)
        out.append(f'<line x1="{sx1:.1f}" y1="{sy1:.1f}" x2="{sx2:.1f}" y2="{sy2:.1f}"'
                   f' stroke="{WALL_STROKE}" stroke-width="{WALL_SW}"/>')

    # ---- IW7 (solid, no opening, rotated parallel to W18-W1) ----
    _iw7_sw, _iw7_se, _iw7_ne, _iw7_nw = layout.iw7.poly
    wall_poly(out, layout.iw7.poly, to_svg, stroke=False)
    # Thickness unit vector: NW - SW direction (norm direction)
    _iw7_dx_t = _iw7_nw[0] - _iw7_sw[0]
    _iw7_dy_t = _iw7_nw[1] - _iw7_sw[1]
    _iw7_lt = math.sqrt(_iw7_dx_t**2 + _iw7_dy_t**2)
    _iw7_at = (_iw7_dx_t / _iw7_lt, _iw7_dy_t / _iw7_lt)
    for (p1, p2), (ox, oy) in [
        ((_iw7_nw, _iw7_ne), (-_iw7_at[0], -_iw7_at[1])),        # north face
        ((_iw7_sw, _iw7_se), _iw7_at),                            # south face
    ]:
        sx1, sy1 = to_svg(p1[0] + half_sw * ox, p1[1] + half_sw * oy)
        sx2, sy2 = to_svg(p2[0] + half_sw * ox, p2[1] + half_sw * oy)
        out.append(f'<line x1="{sx1:.1f}" y1="{sy1:.1f}" x2="{sx2:.1f}" y2="{sy2:.1f}"'
                   f' stroke="{WALL_STROKE}" stroke-width="{WALL_SW}"/>')

    # ---- IW9 with RO7 and RO3 (rotated rectangle, split by two openings) ----
    _iw9_sw, _iw9_se, _iw9_ne, _iw9_nw = layout.iw9.poly
    # Thickness unit vector: SW - SE direction
    _iw9_dx_t = _iw9_sw[0] - _iw9_se[0]
    _iw9_dy_t = _iw9_sw[1] - _iw9_se[1]
    _iw9_lt = math.sqrt(_iw9_dx_t**2 + _iw9_dy_t**2)
    _iw9_at = (_iw9_dx_t / _iw9_lt, _iw9_dy_t / _iw9_lt)
    # Length unit vector: NE - SE direction
    _iw9_dx_n = _iw9_ne[0] - _iw9_se[0]
    _iw9_dy_n = _iw9_ne[1] - _iw9_se[1]
    _iw9_ln = math.sqrt(_iw9_dx_n**2 + _iw9_dy_n**2)
    _iw9_an = (_iw9_dx_n / _iw9_ln, _iw9_dy_n / _iw9_ln)

    # RO7 and RO3 polygon corners on IW9
    ro7_ro = [r for r in rough_openings if r.name == "RO7"][0]
    _ro7p = ro7_ro.poly  # [SW, SE, NE, NW]
    _ro3p = [r for r in rough_openings if r.name == "RO3"][0].poly

    # Three solid segments: S end→RO7, RO7→RO3, RO3→N end
    iw9_s_poly = [_iw9_sw, _iw9_se, _ro7p[0], _ro7p[1]]
    iw9_m_poly = [_ro7p[2], _ro7p[3], _ro3p[0], _ro3p[1]]
    iw9_n_poly = [_ro3p[2], _ro3p[3], _iw9_ne, _iw9_nw]
    wall_poly(out, iw9_s_poly, to_svg, stroke=False)
    wall_poly(out, iw9_m_poly, to_svg, stroke=False)
    wall_poly(out, iw9_n_poly, to_svg, stroke=False)

    _neg_iw9_at = (-_iw9_at[0], -_iw9_at[1])
    _neg_iw9_an = (-_iw9_an[0], -_iw9_an[1])
    # Stroke lines for south segment
    for (p1, p2), (ox, oy) in [
        ((_iw9_se, _ro7p[0]), _iw9_at),        # east face south
        ((_iw9_sw, _ro7p[1]), _neg_iw9_at),    # west face south
        ((_iw9_sw, _iw9_se), _iw9_an),         # south end
    ]:
        sx1, sy1 = to_svg(p1[0] + half_sw * ox, p1[1] + half_sw * oy)
        sx2, sy2 = to_svg(p2[0] + half_sw * ox, p2[1] + half_sw * oy)
        out.append(f'<line x1="{sx1:.1f}" y1="{sy1:.1f}" x2="{sx2:.1f}" y2="{sy2:.1f}"'
                   f' stroke="{WALL_STROKE}" stroke-width="{WALL_SW}"/>')
    # Stroke lines for middle segment (between RO7 and RO3)
    for (p1, p2), (ox, oy) in [
        ((_ro7p[3], _ro3p[0]), _iw9_at),       # east face middle
        ((_ro7p[2], _ro3p[1]), _neg_iw9_at),   # west face middle
    ]:
        sx1, sy1 = to_svg(p1[0] + half_sw * ox, p1[1] + half_sw * oy)
        sx2, sy2 = to_svg(p2[0] + half_sw * ox, p2[1] + half_sw * oy)
        out.append(f'<line x1="{sx1:.1f}" y1="{sy1:.1f}" x2="{sx2:.1f}" y2="{sy2:.1f}"'
                   f' stroke="{WALL_STROKE}" stroke-width="{WALL_SW}"/>')
    # Stroke lines for north segment
    for (p1, p2), (ox, oy) in [
        ((_ro3p[3], _iw9_ne), _iw9_at),        # east face north
        ((_ro3p[2], _iw9_nw), _neg_iw9_at),    # west face north
        ((_iw9_ne, _iw9_nw), _neg_iw9_an),     # north end
    ]:
        sx1, sy1 = to_svg(p1[0] + half_sw * ox, p1[1] + half_sw * oy)
        sx2, sy2 = to_svg(p2[0] + half_sw * ox, p2[1] + half_sw * oy)
        out.append(f'<line x1="{sx1:.1f}" y1="{sy1:.1f}" x2="{sx2:.1f}" y2="{sy2:.1f}"'
                   f' stroke="{WALL_STROKE}" stroke-width="{WALL_SW}"/>')
    # Jambs at RO7 and RO3 edges
    for (j1, j2) in [(_ro7p[0], _ro7p[1]), (_ro7p[3], _ro7p[2]),
                      (_ro3p[0], _ro3p[1]), (_ro3p[3], _ro3p[2])]:
        _jn = (_iw9_an[0] * JAMB_WIDTH, _iw9_an[1] * JAMB_WIDTH)
        j_poly = [j1, j2, (j2[0] + _jn[0], j2[1] + _jn[1]),
                  (j1[0] + _jn[0], j1[1] + _jn[1])]
        jp = " ".join(f"{to_svg(*p)[0]:.1f},{to_svg(*p)[1]:.1f}" for p in j_poly)
        out.append(f'<polygon points="{jp}" fill="{JAMB_COLOR}" stroke="none"/>')

    # ---- IW6 with RO5 ---- (omitted in bare)
    if not bare:
        _iw6 = layout.iw6.poly
        _ro5p = [r for r in rough_openings if r.name == "RO5"][0].poly
        iw6_w_poly = [_iw6[0], _ro5p[0], _ro5p[3], _iw6[3]]
        iw6_e_poly = [_ro5p[1], _iw6[1], _iw6[2], _ro5p[2]]
        wall_poly(out, iw6_w_poly, to_svg, stroke=False)
        wall_poly(out, iw6_e_poly, to_svg, stroke=False)
        _wall_stroke(out, _iw6[0], _ro5p[0], half_sw, to_svg)    # south face west
        _wall_stroke(out, _ro5p[1], _iw6[1], half_sw, to_svg)    # south face east
        _wall_stroke(out, _iw6[3], _ro5p[3], half_sw, to_svg)    # north face west
        _wall_stroke(out, _ro5p[2], _iw6[2], half_sw, to_svg)    # north face east
        _iw6_s_al, _ = seg_vecs(_iw6[0], _iw6[1])
        _jamb_poly(out, _ro5p[3], _ro5p[0], _iw6_s_al, to_svg)
        _neg_iw6_al = (-_iw6_s_al[0], -_iw6_s_al[1])
        _jamb_poly(out, _ro5p[1], _ro5p[2], _neg_iw6_al, to_svg)

    # ---- IW4 (solid, no opening) ----
    _iw4_sw, _iw4_se, _iw4_ne, _iw4_nw = layout.iw4.poly
    wall_poly(out, layout.iw4.poly, to_svg, stroke=False)
    _wall_stroke(out, _iw4_nw, _iw4_sw, half_sw, to_svg)   # west face
    _wall_stroke(out, _iw4_se, _iw4_ne, half_sw, to_svg)    # east face
    _wall_stroke(out, _iw4_ne, _iw4_nw, half_sw, to_svg)    # north face

    # ---- IW11 with RO6 and RO2 (rotated rectangle, split by two openings) ----
    _iw11_sw, _iw11_se, _iw11_ne, _iw11_nw = layout.iw11.poly
    # Unit vectors from polygon corners
    _iw11_dx_t = _iw11_sw[0] - _iw11_se[0]  # thickness direction (along)
    _iw11_dy_t = _iw11_sw[1] - _iw11_se[1]
    _iw11_lt = math.sqrt(_iw11_dx_t**2 + _iw11_dy_t**2)
    _iw11_at = (_iw11_dx_t / _iw11_lt, _iw11_dy_t / _iw11_lt)  # unit along
    _iw11_dx_n = _iw11_ne[0] - _iw11_se[0]  # length direction (normal)
    _iw11_dy_n = _iw11_ne[1] - _iw11_se[1]
    _iw11_ln = math.sqrt(_iw11_dx_n**2 + _iw11_dy_n**2)
    _iw11_an = (_iw11_dx_n / _iw11_ln, _iw11_dy_n / _iw11_ln)  # unit normal

    # RO6 and RO2 polygon corners on IW11
    ro6_ro = [r for r in rough_openings if r.name == "RO6"][0]
    _ro6p = ro6_ro.poly  # [SW, SE, NE, NW]
    ro2_ro = [r for r in rough_openings if r.name == "RO2"][0]
    _ro2p = ro2_ro.poly  # [SW, SE, NE, NW]

    # South segment: IW11 south end to RO6 south edge
    iw11_s_poly = [_iw11_sw, _iw11_se, _ro6p[0], _ro6p[1]]
    wall_poly(out, iw11_s_poly, to_svg, stroke=False)
    # Middle segment: RO6 north edge to RO2 south edge
    iw11_m_poly = [_ro6p[2], _ro6p[3], _ro2p[0], _ro2p[1]]
    wall_poly(out, iw11_m_poly, to_svg, stroke=False)
    # North segment: RO2 north edge to IW11 north end
    iw11_n_poly = [_ro2p[2], _ro2p[3], _iw11_ne, _iw11_nw]
    wall_poly(out, iw11_n_poly, to_svg, stroke=False)

    # Inset stroke lines for south segment
    for (p1, p2), (ox, oy) in [
        ((_iw11_se, _ro6p[0]), _iw11_at),     # east face south
        ((_iw11_sw, _ro6p[1]), (-_iw11_at[0], -_iw11_at[1])),  # west face south
        ((_iw11_sw, _iw11_se), _iw11_an),      # south end
    ]:
        sx1, sy1 = to_svg(p1[0] + half_sw * ox, p1[1] + half_sw * oy)
        sx2, sy2 = to_svg(p2[0] + half_sw * ox, p2[1] + half_sw * oy)
        out.append(f'<line x1="{sx1:.1f}" y1="{sy1:.1f}" x2="{sx2:.1f}" y2="{sy2:.1f}"'
                   f' stroke="{WALL_STROKE}" stroke-width="{WALL_SW}"/>')
    # Inset stroke lines for middle segment
    for (p1, p2), (ox, oy) in [
        ((_ro6p[3], _ro2p[0]), _iw11_at),     # east face middle
        ((_ro6p[2], _ro2p[1]), (-_iw11_at[0], -_iw11_at[1])),  # west face middle
    ]:
        sx1, sy1 = to_svg(p1[0] + half_sw * ox, p1[1] + half_sw * oy)
        sx2, sy2 = to_svg(p2[0] + half_sw * ox, p2[1] + half_sw * oy)
        out.append(f'<line x1="{sx1:.1f}" y1="{sy1:.1f}" x2="{sx2:.1f}" y2="{sy2:.1f}"'
                   f' stroke="{WALL_STROKE}" stroke-width="{WALL_SW}"/>')
    # Inset stroke lines for north segment
    for (p1, p2), (ox, oy) in [
        ((_ro2p[3], _iw11_ne), _iw11_at),     # east face north
        ((_ro2p[2], _iw11_nw), (-_iw11_at[0], -_iw11_at[1])),  # west face north
        ((_iw11_ne, _iw11_nw), (-_iw11_an[0], -_iw11_an[1])),  # north end
    ]:
        sx1, sy1 = to_svg(p1[0] + half_sw * ox, p1[1] + half_sw * oy)
        sx2, sy2 = to_svg(p2[0] + half_sw * ox, p2[1] + half_sw * oy)
        out.append(f'<line x1="{sx1:.1f}" y1="{sy1:.1f}" x2="{sx2:.1f}" y2="{sy2:.1f}"'
                   f' stroke="{WALL_STROKE}" stroke-width="{WALL_SW}"/>')
    # Jambs at RO6 edges
    for (j1, j2) in [(_ro6p[0], _ro6p[1]), (_ro6p[3], _ro6p[2])]:
        _jn = (_iw11_an[0] * JAMB_WIDTH, _iw11_an[1] * JAMB_WIDTH)
        j_poly = [j1, j2, (j2[0] + _jn[0], j2[1] + _jn[1]),
                  (j1[0] + _jn[0], j1[1] + _jn[1])]
        jp = " ".join(f"{to_svg(*p)[0]:.1f},{to_svg(*p)[1]:.1f}" for p in j_poly)
        out.append(f'<polygon points="{jp}" fill="{JAMB_COLOR}" stroke="none"/>')
    # Jambs at RO2 edges
    for (j1, j2) in [(_ro2p[0], _ro2p[1]), (_ro2p[3], _ro2p[2])]:
        _jn = (_iw11_an[0] * JAMB_WIDTH, _iw11_an[1] * JAMB_WIDTH)
        j_poly = [j1, j2, (j2[0] + _jn[0], j2[1] + _jn[1]),
                  (j1[0] + _jn[0], j1[1] + _jn[1])]
        jp = " ".join(f"{to_svg(*p)[0]:.1f},{to_svg(*p)[1]:.1f}" for p in j_poly)
        out.append(f'<polygon points="{jp}" fill="{JAMB_COLOR}" stroke="none"/>')

    # ---- IW12 (rotated rectangle) ----
    wall_poly(out, layout.iw12.poly, to_svg, stroke=False)
    _iw12_sw, _iw12_se, _iw12_ne, _iw12_nw = layout.iw12.poly
    _iw12_dx_t = _iw12_se[0] - _iw12_sw[0]  # length direction (-along)
    _iw12_dy_t = _iw12_se[1] - _iw12_sw[1]
    _iw12_lt = math.sqrt(_iw12_dx_t**2 + _iw12_dy_t**2)
    _iw12_al = (_iw12_dx_t / _iw12_lt, _iw12_dy_t / _iw12_lt)  # unit along length
    _iw12_dx_n = _iw12_nw[0] - _iw12_sw[0]  # thickness direction (norm)
    _iw12_dy_n = _iw12_nw[1] - _iw12_sw[1]
    _iw12_ln = math.sqrt(_iw12_dx_n**2 + _iw12_dy_n**2)
    _iw12_an = (_iw12_dx_n / _iw12_ln, _iw12_dy_n / _iw12_ln)  # unit normal
    for (p1, p2), (ox, oy) in [
        ((_iw12_sw, _iw12_se), _iw12_an),      # south face, inset toward north
        ((_iw12_nw, _iw12_ne), (-_iw12_an[0], -_iw12_an[1])),  # north face, inset toward south
        ((_iw12_sw, _iw12_nw), _iw12_al),      # west end, inset toward east
        ((_iw12_se, _iw12_ne), (-_iw12_al[0], -_iw12_al[1])),  # east end, inset toward west
    ]:
        sx1, sy1 = to_svg(p1[0] + half_sw * ox, p1[1] + half_sw * oy)
        sx2, sy2 = to_svg(p2[0] + half_sw * ox, p2[1] + half_sw * oy)
        out.append(f'<line x1="{sx1:.1f}" y1="{sy1:.1f}" x2="{sx2:.1f}" y2="{sy2:.1f}"'
                   f' stroke="{WALL_STROKE}" stroke-width="{WALL_SW}"/>')

    # ---- IW5 ----
    iw5 = layout.iw5
    wall_poly(out, iw5.poly, to_svg, stroke=False)
    # Stroke lines: inset half_sw from south and north faces (poly-based)
    for (p_start, p_end) in [(iw5.poly[0], iw5.poly[1]),   # south face SW→SE
                              (iw5.poly[3], iw5.poly[2])]:  # north face NW→NE
        _al, _out = seg_vecs(p_start, p_end)
        _inw = (-_out[0], -_out[1])
        _p1 = offset_pt(p_start, half_sw, _inw)
        _p2 = offset_pt(p_end, half_sw, _inw)
        sx1, sy1 = to_svg(*_p1)
        sx2, sy2 = to_svg(*_p2)
        out.append(f'<line x1="{sx1:.1f}" y1="{sy1:.1f}" x2="{sx2:.1f}" y2="{sy2:.1f}"'
                   f' stroke="{WALL_STROKE}" stroke-width="{WALL_SW}"/>')



def _render_appliances(out, data, layout, minik=False, db=False):
    """Render utility room appliances: dryer, washer, counter, water heater, toilets, sinks."""
    pts = data.pts
    to_svg = data.to_svg
    # East wall (W2-W5) direction vectors for wall-relative placement
    w2w5_al, w2w5_in = seg_vecs(pts["W2"], pts["W5"])
    # IW2s east face: along=north, outward=east (shower section)
    _iw2_e_al, _iw2_e_out = seg_vecs(layout.iw2s.poly[1], layout.iw2s.poly[2])

    # Dryer and washer
    minik_appl_w = 32.0 / 12.0   # 32" along w2w5_in in minik
    minik_appl_d = 27.0 / 12.0   # 27" along w2w5_al in minik
    minik_appl_links = {
        "DRYER": "https://www.lowes.com/pd/Electrolux-8-cu-ft-Stackable-Steam-Cycle-Electric-Dryer-Titanium-ENERGY-STAR/5015416377",
        "WASHER": "https://www.lowes.com/pd/Electrolux-Smartboost-Optic-Whites-and-Pure-Rinse-4-5-cu-ft-High-Efficiency-Stackable-Steam-Cycle-Front-Load-Washer-Titanium-ENERGY-STAR/5015416375",
    }
    # Standard shift: 4" inward from W2-W5, 2" opposite along-wall
    _shift = (4.0 / 12.0 * w2w5_in[0] + (-2.0 / 12.0) * w2w5_al[0],
              4.0 / 12.0 * w2w5_in[1] + (-2.0 / 12.0) * w2w5_al[1])
    _dryer_nw_mk = None
    _small_wd = minik or db
    for label, wall_obj in [("DRYER", layout.dryer), ("WASHER", layout.washer)]:
        poly = list(wall_obj.poly)  # [SW, SE, NE, NW]
        if not _small_wd:
            poly = [(p[0] + _shift[0], p[1] + _shift[1]) for p in poly]
        if _small_wd:
            if label == "DRYER":
                _sw = poly[0]
                _se = offset_pt(_sw, minik_appl_w, w2w5_in)
                _nw = offset_pt(_sw, minik_appl_d, w2w5_al)
                _ne = offset_pt(_se, minik_appl_d, w2w5_al)
                poly = [_sw, _se, _ne, _nw]
                _dryer_nw_mk = _nw
            else:  # WASHER: 1" along wall from dryer
                _sw = offset_pt(_dryer_nw_mk, 1.0 / 12.0, w2w5_al)
                _se = offset_pt(_sw, minik_appl_w, w2w5_in)
                _nw = offset_pt(_sw, minik_appl_d, w2w5_al)
                _ne = offset_pt(_se, minik_appl_d, w2w5_al)
                poly = [_sw, _se, _ne, _nw]
        link = minik_appl_links.get(label) if _small_wd else None
        _appl_poly(out, poly, to_svg, label=label, href=link, close_href=False)
        # Door: hinged at NE corner, swings from inward to -along
        _hinge = poly[2]  # NE corner
        _door_len = math.sqrt((poly[3][0] - poly[0][0])**2 + (poly[3][1] - poly[0][1])**2)
        _hx, _hy = to_svg(*_hinge)
        _tip = offset_pt(_hinge, _door_len, w2w5_in)
        _tx, _ty = to_svg(*_tip)
        out.append(f'<line x1="{_hx:.1f}" y1="{_hy:.1f}" x2="{_tx:.1f}" y2="{_ty:.1f}"'
                   f' stroke="{APPL_STROKE}" stroke-width="1.0"/>')
        _arc_pts = _swing_arc_svg(_hinge, _door_len, w2w5_in,
                                   (-w2w5_al[0], -w2w5_al[1]), to_svg)
        out.append(f'<polyline points="{_arc_pts}" fill="none"'
                   f' stroke="{APPL_STROKE}" stroke-width="0.5"/>')
        if link:
            out.append('</a>')
    washer_poly = poly  # save last poly for hamper positioning

    # Hamper: 31.5" x 19", 2" along wall from washer, 2" inward from W2-W5
    hm_ew = 31.5 / 12.0   # width along w2w5_in
    hm_ns = 19.0 / 12.0   # depth along w2w5_al
    # Anchor: 2" inward from W2, washer NW + 2" along wall
    _washer_nw_d = ((washer_poly[3][0] - pts["W2"][0]) * w2w5_al[0] +
                    (washer_poly[3][1] - pts["W2"][1]) * w2w5_al[1])
    _hm_sw = offset_pt(offset_pt(pts["W2"], _washer_nw_d + 2.0 / 12.0, w2w5_al),
                        2.0 / 12.0, w2w5_in)
    _hm_se = offset_pt(_hm_sw, hm_ew, w2w5_in)
    _hm_nw = offset_pt(_hm_sw, hm_ns, w2w5_al)
    _hm_ne = offset_pt(_hm_se, hm_ns, w2w5_al)
    hm_poly = [_hm_sw, _hm_se, _hm_ne, _hm_nw]
    hm_href = "https://www.homedepot.com/p/Casual-Home-Eco-Home-Laundry-Prep-Hamper-761-30/307595219"
    out.append(f'<a href="{hm_href}" target="_blank">')
    _hm_svg = " ".join(f"{to_svg(*p)[0]:.1f},{to_svg(*p)[1]:.1f}" for p in hm_poly)
    out.append(f'<polygon points="{_hm_svg}" fill="{APPL_FILL}" stroke="{APPL_STROKE}" stroke-width="{APPL_SW}"/>')
    # Dashed basket pull-out along wall
    _hm_bo_nw = offset_pt(_hm_nw, hm_ns, w2w5_al)
    _hm_bo_ne = offset_pt(_hm_ne, hm_ns, w2w5_al)
    hm_bo_poly = [_hm_nw, _hm_ne, _hm_bo_ne, _hm_bo_nw]
    _hm_bo_svg = " ".join(f"{to_svg(*p)[0]:.1f},{to_svg(*p)[1]:.1f}" for p in hm_bo_poly)
    out.append(f'<polygon points="{_hm_bo_svg}" fill="none" stroke="{APPL_STROKE}" stroke-width="{APPL_SW}"'
               f' stroke-dasharray="3,2"/>')
    _hm_cx = sum(p[0] for p in hm_poly) / 4
    _hm_cy = sum(p[1] for p in hm_poly) / 4
    _hm_scx, _hm_scy = to_svg(_hm_cx, _hm_cy)
    out.append(f'<text x="{_hm_scx:.1f}" y="{_hm_scy+3:.1f}" text-anchor="middle" font-family="Arial"'
               f' font-size="6" fill="{APPL_STROKE}">HAMPER</text>')
    out.append('</a>')

    # Counter: polygon clipped to W18-W1 south edge, IW3 west face east edge
    ctr_poly_svg = " ".join(f"{to_svg(*p)[0]:.1f},{to_svg(*p)[1]:.1f}"
                            for p in layout.ctr_clip)
    out.append(f'<polygon points="{ctr_poly_svg}"'
               f' fill="{APPL_FILL}" stroke="{APPL_STROKE}" stroke-width="{APPL_SW}"/>')
    _ctr_cx = sum(p[0] for p in layout.ctr_clip) / len(layout.ctr_clip)
    _ctr_cy = sum(p[1] for p in layout.ctr_clip) / len(layout.ctr_clip)
    ccx, ccy = to_svg(_ctr_cx, _ctr_cy)
    out.append(f'<text x="{ccx:.1f}" y="{ccy:.1f}" text-anchor="middle" font-family="Arial"'
               f' font-size="7" fill="{APPL_STROKE}" letter-spacing="0.5" transform="rotate(-90,{ccx:.1f},{ccy:.1f})">COUNTER</text>')

    # Water heater: 28" diameter circle, tangent to inner arc at C7
    # WH center lies on IW2s east face line at WH_RADIUS from face
    _wh_ref = offset_pt(layout.iw2s.poly[2], WH_RADIUS, _iw2_e_out)
    wh_tangent_r = (data.radii["R_a7"] - data.wall_t) - WH_RADIUS
    _c7 = pts["C7"]
    _wh_d = (_wh_ref[0] - _c7[0], _wh_ref[1] - _c7[1])
    _wh_d_al = _wh_d[0] * _iw2_e_al[0] + _wh_d[1] * _iw2_e_al[1]
    _wh_d2 = _wh_d[0]**2 + _wh_d[1]**2
    _wh_t = -_wh_d_al + math.sqrt(wh_tangent_r**2 - _wh_d2 + _wh_d_al**2)
    wh_center = offset_pt(_wh_ref, _wh_t, _iw2_e_al)
    wh_e, wh_n = wh_center
    wh_sx, wh_sy = to_svg(wh_e, wh_n)
    wh_r_svg = (to_svg(WH_RADIUS, 0)[0] - to_svg(0, 0)[0])
    out.append(f'<circle cx="{wh_sx:.1f}" cy="{wh_sy:.1f}" r="{wh_r_svg:.1f}"'
               f' fill="{APPL_FILL}" stroke="{APPL_STROKE}" stroke-width="{APPL_SW}"/>')
    out.append(f'<text x="{wh_sx:.1f}" y="{wh_sy+3:.1f}" text-anchor="middle" font-family="Arial"'
               f' font-size="7" fill="{APPL_STROKE}">WH</text>')

    # Toilets and sinks — oriented relative to IW8 wall direction
    _iw8_al, _iw8_in = seg_vecs(layout.iw8.poly[0], layout.iw8.poly[1])
    _iw8_out = (-_iw8_in[0], -_iw8_in[1])
    _iw8_s_ref = layout.iw8.poly[0]  # SW corner of IW8
    # Project dryer centroid onto IW8 south face for toilet position
    _dryer_cx = sum(p[0] for p in layout.dryer.poly) / 4
    _dryer_cy = sum(p[1] for p in layout.dryer.poly) / 4
    _d_dryer_al = ((_dryer_cx - _iw8_s_ref[0]) * _iw8_al[0] +
                   (_dryer_cy - _iw8_s_ref[1]) * _iw8_al[1])
    _toilet_s = offset_pt(_iw8_s_ref, _d_dryer_al, _iw8_al)
    # Project midpoint between dryer SE and counter SW for sink position
    _ctr_cx = sum(p[0] for p in layout.ctr.poly) / 4
    _ctr_cy = sum(p[1] for p in layout.ctr.poly) / 4
    _sink_mid = ((_dryer_cx + _ctr_cx) / 2, (_dryer_cy + _ctr_cy) / 2)
    _d_sink_al = ((_sink_mid[0] - _iw8_s_ref[0]) * _iw8_al[0] +
                  (_sink_mid[1] - _iw8_s_ref[1]) * _iw8_al[1])
    _sink_s = offset_pt(_iw8_s_ref, _d_sink_al, _iw8_al)
    # South toilet: faces along _iw8_in (away from IW8 south face)
    draw_toilet(out, _toilet_s, _iw8_in, _iw8_al, to_svg)
    _sk_s = offset_pt(_sink_s, SINK_RY, _iw8_in)
    draw_sink(out, _sk_s[0], _sk_s[1], to_svg=to_svg)
    # North toilet: faces opposite direction (away from IW8 north face)
    _toilet_n = offset_pt(_iw8_s_ref, _d_dryer_al, _iw8_al)
    # Shift to north face: IW8 poly[3]→poly[2] is north face
    _iw8_n_ref = layout.iw8.poly[3]
    _d_toilet_n_al = ((_dryer_cx - _iw8_n_ref[0]) * _iw8_al[0] +
                      (_dryer_cy - _iw8_n_ref[1]) * _iw8_al[1])
    _toilet_n = offset_pt(_iw8_n_ref, _d_toilet_n_al, _iw8_al)
    draw_toilet(out, _toilet_n, _iw8_out, _iw8_al, to_svg)
    _d_sink_n_al = ((_sink_mid[0] - _iw8_n_ref[0]) * _iw8_al[0] +
                    (_sink_mid[1] - _iw8_n_ref[1]) * _iw8_al[1])
    _sink_n_face = offset_pt(_iw8_n_ref, _d_sink_n_al, _iw8_al)
    _sk_n = offset_pt(_sink_n_face, SINK_RY, _iw8_out)
    draw_sink(out, _sk_n[0], _sk_n[1], to_svg=to_svg)


def _render_kitchen(out, data, layout, minik=False, db=False):
    """Render kitchen: D/W, sink, stove, shelves, fridge, counters."""
    pts = data.pts
    to_svg = data.to_svg
    # North wall (W9-W10) direction vectors for wall-relative placement
    w9w10_al, w9w10_in = seg_vecs(pts["W9"], pts["W10"])
    _wall_anchor = pts["W9"]
    # Distance from W9 to IW2s east face along wall
    _iw2s_ne = layout.iw2s.poly[2]
    _iw2_d = ((_iw2s_ne[0] - _wall_anchor[0]) * w9w10_al[0] +
              (_iw2s_ne[1] - _wall_anchor[1]) * w9w10_al[1])

    def _nwp(d_along, d_inward=0):
        """Point on north wall at d_along from W9, d_inward into room."""
        return offset_pt(offset_pt(_wall_anchor, d_along, w9w10_al),
                         d_inward, w9w10_in)

    # Kitchen appliance chain distances along wall from W9
    _st_d = _iw2_d + NORTH_CTR_LENGTH + KITCHEN_APPL_GAP
    _ks_d = _st_d + STOVE_WIDTH + KITCHEN_APPL_GAP + 2.0 / 12.0
    _dw_d = _ks_d + KITCHEN_SINK_WIDTH + KITCHEN_APPL_GAP
    # IW1/IW2 corner: for items in the fridge/counter/table area
    _iw2_lo_e_al, _ = seg_vecs(layout.iw2.poly[1], layout.iw2.poly[2])
    _iw1_n_al, _iw1_n_cw = seg_vecs(layout.iw1.poly[3], layout.iw1.poly[2])
    _iw1_n_out = (-_iw1_n_cw[0], -_iw1_n_cw[1])  # outward from IW1 north face
    _iw12_corner = line_isect(layout.iw2.poly[1], _iw2_lo_e_al,
                              layout.iw1.poly[3], _iw1_n_al)

    def _iwp(d_e, d_n=0):
        """Point from IW1/IW2 corner: d_e along IW1-along, d_n along IW1-outward."""
        return offset_pt(offset_pt(_iw12_corner, d_e, _iw1_n_al),
                         d_n, _iw1_n_out)


    # (label, along_start, along_width, inward_start, inward_depth, href)
    appliances = [
        ("SINK",  _ks_d, KITCHEN_SINK_WIDTH, 0, KITCHEN_SINK_DEPTH,
         "https://www.webstaurantstore.com/advance-tabco-fs1181824l-45-fabricated-one-compartment-sink-with-24-left-drainboard-18-x-18-x-14-bowl/109FS1L241818.html"),
        ("STOVE", _st_d, STOVE_WIDTH, KITCHEN_APPL_GAP, KITCHEN_APPL_GAP + STOVE_DEPTH, None),
        ("D/W",   _dw_d, DW_WIDTH, 0, DW_DEPTH, None),
    ]
    if minik:
        appliances = [(l, da, aw, di0, di1, h) for l, da, aw, di0, di1, h in appliances
                      if l not in ("STOVE", "D/W")]
    for label, da, aw, di0, di1, href in appliances:
        _c = [_nwp(da, di1), _nwp(da + aw, di1), _nwp(da + aw, di0), _nwp(da, di0)]
        _appl_poly(out, _c, to_svg, label=label, href=href)
        if label == "STOVE":
            _ext = 24.0 / 12.0
            _ec = [_nwp(da, di1), _nwp(da + aw, di1),
                   _nwp(da + aw, di1 + _ext), _nwp(da, di1 + _ext)]
            _appl_poly(out, _ec, to_svg, dash=True, fill_color="none")
        elif label == "D/W":
            _ext = 31.0 / 12.0
            _ec = [_nwp(da, di1), _nwp(da + aw, di1),
                   _nwp(da + aw, di1 + _ext), _nwp(da, di1 + _ext)]
            _appl_poly(out, _ec, to_svg, dash=True, fill_color="none")

    # Fridge — positioned relative to wall/IW vectors
    _fr_w2 = 32.75 / 12.0   # std/db fridge width along IW2-outward
    _fr_h1 = 35.0 / 12.0    # std/db fridge depth along IW1-outward
    if minik:
        # 3" along wall past SINK, 3" inward from north wall
        _fr_mk_d = _ks_d + KITCHEN_SINK_WIDTH + 3.0 / 12.0
        _fr_mk_i = 3.0 / 12.0
        fr_nw = _nwp(_fr_mk_d, _fr_mk_i)
        fr_se = _nwp(_fr_mk_d + MINIK_FRIDGE_W, _fr_mk_i + MINIK_FRIDGE_D)
        fr_ne = _nwp(_fr_mk_d + MINIK_FRIDGE_W, _fr_mk_i)
        fr_sw = _nwp(_fr_mk_d, _fr_mk_i + MINIK_FRIDGE_D)
    else:
        # IW1/IW2 corner, gap from each wall face
        fr_nw = _iwp(KITCHEN_APPL_GAP, KITCHEN_APPL_GAP + _fr_h1)
        fr_se = _iwp(KITCHEN_APPL_GAP + _fr_w2, KITCHEN_APPL_GAP)
        fr_ne = _iwp(KITCHEN_APPL_GAP + _fr_w2, KITCHEN_APPL_GAP + _fr_h1)
        fr_sw = _iwp(KITCHEN_APPL_GAP, KITCHEN_APPL_GAP)
    _fr_href = ("https://www.ikea.com/us/en/p/bergsnaes-bottom-freezer-refrigerator-stainless-steel-color-60607883/"
                if minik else
                "https://www.lowes.com/pd/LG-25-5-cu-ft-Bottom-Freezer-Refrigerator-with-Ice-Maker-Fingerprint-Resistant-Printproof-Stainless-Steel-ENERGY-STAR/1002543648")
    fr_fs = 6 if minik else 7
    _appl_poly(out, [fr_sw, fr_se, fr_ne, fr_nw], to_svg,
               label="FRIDGE", href=_fr_href, font_size=str(fr_fs), close_href=False)
    if minik:
        # Door arc: hinged at SE corner, sweeps from closed (back-along-wall) to open (inward)
        fr_door = MINIK_FRIDGE_W
        _close_dir = (-w9w10_al[0], -w9w10_al[1])  # back along wall
        _open_dir = w9w10_in                         # inward
        hx, hy = to_svg(*fr_se)
        tip = offset_pt(fr_se, fr_door, _open_dir)
        tip_x, tip_y = to_svg(*tip)
        out.append(f'<line x1="{hx:.1f}" y1="{hy:.1f}" x2="{tip_x:.1f}" y2="{tip_y:.1f}"'
                   f' stroke="{APPL_STROKE}" stroke-width="1.0"/>')
        arc_pts = _swing_arc_svg(fr_se, fr_door, _close_dir, _open_dir, to_svg)
        out.append(f'<polyline points="{arc_pts}" fill="none"'
                   f' stroke="{APPL_STROKE}" stroke-width="0.5"/>')
        out.append('</a>')
    if db or (not minik):
        # Door arc: hinged at NW corner, sweeps from open (northward) to closed (eastward)
        fr_door = _fr_w2
        _open_dir = _iw1_n_out   # outward from IW1 (northward) — door open
        _close_dir = _iw1_n_al   # along IW1 north face (eastward) — door closed
        hx, hy = to_svg(*fr_nw)
        tip = offset_pt(fr_nw, fr_door, _open_dir)
        tip_x, tip_y = to_svg(*tip)
        out.append(f'<line x1="{hx:.1f}" y1="{hy:.1f}" x2="{tip_x:.1f}" y2="{tip_y:.1f}"'
                   f' stroke="{APPL_STROKE}" stroke-width="1.0"/>')
        arc_pts = _swing_arc_svg(fr_nw, fr_door, _open_dir, _close_dir, to_svg)
        out.append(f'<polyline points="{arc_pts}" fill="none"'
                   f' stroke="{APPL_STROKE}" stroke-width="0.5"/>')
        out.append('</a>')

    # ICE — positioned along north wall
    if db:
        _ice_d = _dw_d + DW_WIDTH + 2.0 / 12.0    # 2" along wall past D/W
        _ice_i = KITCHEN_APPL_GAP                   # 3" inward from wall
    elif minik:
        _ice_d = _ks_d + KITCHEN_SINK_WIDTH + 3.0 / 12.0 + MINIK_FRIDGE_W + 3.0 / 12.0
        _ice_i = 3.0 / 12.0
    else:
        _ice_d = _dw_d + DW_WIDTH + 6.0 / 12.0    # 6" along wall past D/W
        _ice_i = KITCHEN_APPL_GAP
    _ice_c = [_nwp(_ice_d, _ice_i + ICE_DEPTH), _nwp(_ice_d + ICE_WIDTH, _ice_i + ICE_DEPTH),
              _nwp(_ice_d + ICE_WIDTH, _ice_i), _nwp(_ice_d, _ice_i)]
    _appl_poly(out, _ice_c, to_svg, label="ICE", font_size="6",
               href="https://www.homedepot.com/p/EUHOMY-17-3-in-100-lb-24H-Full-Ice-Sizes-Commercial-Ice-Maker-in-Black-33-lb-Storage-Bin-Ice-Full-Alert-and-Auto-Cleaning-CIM001-100BL-E/337185876")

    # Work counter: 60" along IW2-outward x 18" along IW1-outward, against IW1
    _wc_d2 = KITCHEN_APPL_GAP + _fr_w2 + KITCHEN_APPL_GAP  # gap + fridge + gap from IW2
    if not minik:
        _wc_c = [_iwp(_wc_d2, 0), _iwp(_wc_d2 + 60.0 / 12.0, 0),
                 _iwp(_wc_d2 + 60.0 / 12.0, 18.0 / 12.0), _iwp(_wc_d2, 18.0 / 12.0)]
        _appl_poly(out, _wc_c, to_svg,
                   href="https://www.webstaurantstore.com/table-s-s-18x60-s-s-under/600TS1860S.html")

    # Microwave on 18" counter (non-minik): 19.5" along IW2 x 16-5/8" along IW1
    if not minik:
        mw_ew = 19.5 / 12.0
        mw_ns = 16.625 / 12.0
        _mw_d2 = _wc_d2 + 2.0 / 12.0   # 2" from work counter start
        _mw_d1 = 2.0 / 12.0             # 2" from IW1 face
        _mw_c = [_iwp(_mw_d2, _mw_d1), _iwp(_mw_d2 + mw_ew, _mw_d1),
                 _iwp(_mw_d2 + mw_ew, _mw_d1 + mw_ns), _iwp(_mw_d2, _mw_d1 + mw_ns)]
        _appl_poly(out, _mw_c, to_svg, label="MICRO", font_size="5",
                   href="https://www.ikea.com/us/en/p/gatebo-microwave-oven-with-air-fryer-function-ikea-500-black-70603506/",
                   close_href=False)
        # Door: hinged at NE corner, sweeps from open (along IW2) to closed (back along IW2)
        _mw_door = mw_ew
        _mw_hinge = _iwp(_mw_d2 + mw_ew, _mw_d1 + mw_ns)  # NE corner
        _mw_open = _iw1_n_out                                 # northward from IW1
        _mw_close = (-_iw1_n_al[0], -_iw1_n_al[1])           # westward along IW1
        _mh_x, _mh_y = to_svg(*_mw_hinge)
        _mt = offset_pt(_mw_hinge, _mw_door, _mw_open)
        _mt_x, _mt_y = to_svg(*_mt)
        out.append(f'<line x1="{_mh_x:.1f}" y1="{_mh_y:.1f}" x2="{_mt_x:.1f}" y2="{_mt_y:.1f}"'
                   f' stroke="{APPL_STROKE}" stroke-width="1.0"/>')
        _mw_arc = []
        for _i in range(21):
            _t = _i / 20
            _ct = math.cos(_t * math.pi / 2)
            _st = math.sin(_t * math.pi / 2)
            _ae = _mw_hinge[0] + _mw_door * (_ct * _mw_open[0] + _st * _mw_close[0])
            _an = _mw_hinge[1] + _mw_door * (_ct * _mw_open[1] + _st * _mw_close[1])
            _ax, _ay = to_svg(_ae, _an)
            _mw_arc.append(f"{_ax:.1f},{_ay:.1f}")
        out.append(f'<polyline points="{" ".join(_mw_arc)}" fill="none"'
                   f' stroke="{APPL_STROKE}" stroke-width="0.5"/>')
        out.append('</a>')

    # Kitchen counter: along north wall from IW2 east face (minik only)
    if minik:
        _kc_c = [_nwp(_iw2_d, KITCHEN_CTR_DEPTH), _nwp(_iw2_d + KITCHEN_CTR_LENGTH, KITCHEN_CTR_DEPTH),
                 _nwp(_iw2_d + KITCHEN_CTR_LENGTH, 0), _nwp(_iw2_d, 0)]
        _appl_poly(out, _kc_c, to_svg,
                   href="https://www.webstaurantstore.com/regency-spec-line-30-x-72-14-gauge-stainless-steel-commercial-work-table-with-4-backsplash-and-undershelf/600TSSB3072S.html")

    # Minik: microwave on counter (19.5" along wall x 16-5/8" inward)
    if minik:
        mw_ew = 19.5 / 12.0
        mw_ns = 16.625 / 12.0
        _mw_mk_d = _iw2_d + 2.0 / 12.0   # 2" along wall from counter start
        _mw_mk_i = 3.0 / 12.0             # 3" inward from wall
        _mw_mk_c = [_nwp(_mw_mk_d, _mw_mk_i + mw_ns), _nwp(_mw_mk_d + mw_ew, _mw_mk_i + mw_ns),
                     _nwp(_mw_mk_d + mw_ew, _mw_mk_i), _nwp(_mw_mk_d, _mw_mk_i)]
        _appl_poly(out, _mw_mk_c, to_svg, label="MICRO", font_size="5",
                   href="https://www.ikea.com/us/en/p/gatebo-microwave-oven-with-air-fryer-function-ikea-500-black-70603506/",
                   close_href=False)
        # Door: hinged at wall-start/inward corner, sweeps from open (inward) to closed (along wall)
        _mw_door = mw_ew
        _mw_hinge = _nwp(_mw_mk_d, _mw_mk_i + mw_ns)
        _mw_open = w9w10_in                                   # inward
        _mw_close = w9w10_al                                   # along wall
        _mh_x, _mh_y = to_svg(*_mw_hinge)
        _mt = offset_pt(_mw_hinge, _mw_door, _mw_open)
        _mt_x, _mt_y = to_svg(*_mt)
        out.append(f'<line x1="{_mh_x:.1f}" y1="{_mh_y:.1f}" x2="{_mt_x:.1f}" y2="{_mt_y:.1f}"'
                   f' stroke="{APPL_STROKE}" stroke-width="1.0"/>')
        _mw_arc = []
        for _i in range(21):
            _t = _i / 20
            _ct = math.cos(_t * math.pi / 2)
            _st = math.sin(_t * math.pi / 2)
            _ae = _mw_hinge[0] + _mw_door * (_ct * _mw_open[0] + _st * _mw_close[0])
            _an = _mw_hinge[1] + _mw_door * (_ct * _mw_open[1] + _st * _mw_close[1])
            _ax, _ay = to_svg(_ae, _an)
            _mw_arc.append(f"{_ax:.1f},{_ay:.1f}")
        out.append(f'<polyline points="{" ".join(_mw_arc)}" fill="none"'
                   f' stroke="{APPL_STROKE}" stroke-width="0.5"/>')
        out.append('</a>')

    # Minik: coffee maker on counter (7.2" along wall x 9.2" inward)
    if minik:
        cm_ew = 7.2 / 12.0
        cm_ns = 9.2 / 12.0
        _cm_d = _mw_mk_d + mw_ew + 3.0 / 12.0  # 3" along wall past microwave
        _cm_i = 3.0 / 12.0                       # 3" inward from wall
        _cm_c = [_nwp(_cm_d, _cm_i + cm_ns), _nwp(_cm_d + cm_ew, _cm_i + cm_ns),
                 _nwp(_cm_d + cm_ew, _cm_i), _nwp(_cm_d, _cm_i)]
        _appl_poly(out, _cm_c, to_svg, label="C", font_size="5",
                   href="https://www.amazon.com/Holstein-Housewares-HH-0914701E-5-Cup-Coffee/dp/B08HSRCC4T/?th=1")

    # Minik: induction cooktop on counter (13.4" along wall x 20.5" inward)
    if minik:
        cp_ew = 13.4 / 12.0
        cp_ns = 20.5 / 12.0
        _cp_d = _cm_d + cm_ew + 3.0 / 12.0                 # 3" along wall past coffee maker
        _cp_i_far = KITCHEN_CTR_DEPTH - 2.0 / 12.0          # south edge: 2" from counter south
        _cp_i_near = _cp_i_far - cp_ns                      # north edge
        _cp_c = [_nwp(_cp_d, _cp_i_far), _nwp(_cp_d + cp_ew, _cp_i_far),
                 _nwp(_cp_d + cp_ew, _cp_i_near), _nwp(_cp_d, _cp_i_near)]
        out.append('<a href="https://www.homedepot.com/p/Empava-Portable-13-4-in-Induction-Electric-Cooktop-in-Black-with-2-Elements-EMPV-ID12/313815692" target="_blank">')
        _cp_svg = " ".join(f"{to_svg(*p)[0]:.1f},{to_svg(*p)[1]:.1f}" for p in _cp_c)
        out.append(f'<polygon points="{_cp_svg}" fill="#222" stroke="{APPL_STROKE}" stroke-width="{APPL_SW}"/>')
        # Two burner circles (6" diameter each), spaced along inward direction
        burner_r_ft = 3.0 / 12.0
        burner_r_svg = abs(to_svg(burner_r_ft, 0)[0] - to_svg(0, 0)[0])
        _cp_center = _nwp(_cp_d + cp_ew / 2, (_cp_i_near + _cp_i_far) / 2)
        _outward = (-w9w10_in[0], -w9w10_in[1])
        for sign in (-1, 1):
            bp = offset_pt(_cp_center, sign * cp_ns / 4, _outward)
            bsx, bsy = to_svg(*bp)
            out.append(f'<circle cx="{bsx:.1f}" cy="{bsy:.1f}" r="{burner_r_svg:.1f}"'
                       f' fill="none" stroke="#666" stroke-width="0.3"/>')
        out.append('</a>')

    # Minik: toaster 3" along wall past cooktop, 3" inward (13.7" x 12.5")
    if minik:
        ts_ew = 13.7 / 12.0
        ts_ns = 12.5 / 12.0
        _ts_d = _cp_d + cp_ew + 3.0 / 12.0  # 3" along wall past cooktop
        _ts_i = 3.0 / 12.0                   # 3" inward from wall
        _ts_c = [_nwp(_ts_d, _ts_i + ts_ns), _nwp(_ts_d + ts_ew, _ts_i + ts_ns),
                 _nwp(_ts_d + ts_ew, _ts_i), _nwp(_ts_d, _ts_i)]
        out.append('<a href="https://www.amazon.com/Roter-Mond-Stainless-Independent-Removable/dp/B0CGTQZTDZ?th=1" target="_blank">')
        _ts_svg = " ".join(f"{to_svg(*p)[0]:.1f},{to_svg(*p)[1]:.1f}" for p in _ts_c)
        out.append(f'<polygon points="{_ts_svg}" fill="{APPL_FILL}" stroke="{APPL_STROKE}" stroke-width="{APPL_SW}"/>')
        # Four toast slot lines, evenly spaced along wall
        slot_len_ft = 5.5 / 12.0
        slot_spacing = ts_ew / 5.0
        for i in range(4):
            _sl_c = _nwp(_ts_d + slot_spacing * (i + 1), _ts_i + ts_ns / 2)
            _sl_n = offset_pt(_sl_c, slot_len_ft / 2, _outward)
            _sl_s = offset_pt(_sl_c, slot_len_ft / 2, w9w10_in)
            sl_sx1, sl_sy1 = to_svg(*_sl_n)
            sl_sx2, sl_sy2 = to_svg(*_sl_s)
            out.append(f'<line x1="{sl_sx1:.1f}" y1="{sl_sy1:.1f}" x2="{sl_sx2:.1f}" y2="{sl_sy2:.1f}"'
                       f' stroke="#666" stroke-width="0.4"/>')
        out.append('</a>')

    # Oscar triangle dining set centered between north wall, IW1, IW2, RO1
    # Space bounds from wall-relative references (rotation-safe)
    _space_ne_ref = _iwp(RO1_OFFSET_FROM_IW2)
    _space_n_ref = _nwp(0, KITCHEN_CTR_DEPTH) if minik else _nwp(0, 0)

    # Table: base 31.5" (N), height 35.25", 24" arc at apex, 6" fillets
    tbl_base = 31.5 / 12.0
    tbl_h = 35.25 / 12.0
    apex_r = 12.0 / 12.0    # 24" diameter arc at apex
    fillet_r = 6.0 / 12.0   # 6" corner fillets

    # Position: base center along wall, 30"(+28") south of space north edge
    if minik:
        _tbl_ref = _nwp(_st_d + STOVE_WIDTH + KITCHEN_APPL_GAP)
    else:
        _tbl_ref = _nwp(_ks_d + KITCHEN_SINK_WIDTH / 2)
    _tbl_d_al = ((_tbl_ref[0] - _iw12_corner[0]) * _iw1_n_al[0] +
                 (_tbl_ref[1] - _iw12_corner[1]) * _iw1_n_al[1])
    _space_n_d_out = ((_space_n_ref[0] - _iw12_corner[0]) * _iw1_n_out[0] +
                      (_space_n_ref[1] - _iw12_corner[1]) * _iw1_n_out[1])
    _tbl_n_offset = 30.0 / 12.0 + (28.0 / 12.0 if not minik else 0)
    _tbl_n_d_out = _space_n_d_out - _tbl_n_offset
    tbl_bc = offset_pt(offset_pt(_iw12_corner, _tbl_d_al, _iw1_n_al),
                       _tbl_n_d_out, _iw1_n_out)
    _to_apex = (-_iw1_n_out[0], -_iw1_n_out[1])  # toward IW1

    # Base corners and arc center
    ne = offset_pt(tbl_bc, tbl_base / 2, _iw1_n_al)
    nw = offset_pt(tbl_bc, -tbl_base / 2, _iw1_n_al)
    _tbl_apex = offset_pt(tbl_bc, tbl_h, _to_apex)
    arc_c = offset_pt(_tbl_apex, apex_r, _iw1_n_out)

    # Mirror across symmetry axis (through tbl_bc, along _to_apex)
    def _sym(p):
        vx = p[0] - tbl_bc[0]; vy = p[1] - tbl_bc[1]
        v_dot = vx * _to_apex[0] + vy * _to_apex[1]
        return (tbl_bc[0] + 2 * v_dot * _to_apex[0] - vx,
                tbl_bc[1] + 2 * v_dot * _to_apex[1] - vy)

    # Right tangent from NE to apex arc
    dx_r = ne[0] - arc_c[0]
    dn_r = ne[1] - arc_c[1]
    dist_r = math.sqrt(dx_r**2 + dn_r**2)
    angle_cp = math.atan2(dn_r, dx_r)
    delta = math.acos(apex_r / dist_r)
    alpha_r = angle_cp - delta
    t_right = (arc_c[0] + apex_r * math.cos(alpha_r),
                arc_c[1] + apex_r * math.sin(alpha_r))
    t_left = _sym(t_right)

    # NE fillet between base direction (along IW1 face) and tangent line
    d_base_ne = (-_iw1_n_al[0], -_iw1_n_al[1])  # base runs NE→NW = reverse of IW1 along
    dtr = (t_right[0] - ne[0], t_right[1] - ne[1])
    dtr_len = math.sqrt(dtr[0]**2 + dtr[1]**2)
    d_tang_ne = (dtr[0] / dtr_len, dtr[1] / dtr_len)
    cos_th = d_base_ne[0] * d_tang_ne[0] + d_base_ne[1] * d_tang_ne[1]
    half_angle = math.acos(max(-1, min(1, cos_th))) / 2
    fillet_dist = fillet_r / math.sin(half_angle)
    bis_ne = (d_base_ne[0] + d_tang_ne[0], d_base_ne[1] + d_tang_ne[1])
    bis_ne_len = math.sqrt(bis_ne[0]**2 + bis_ne[1]**2)
    bis_ne = (bis_ne[0] / bis_ne_len, bis_ne[1] / bis_ne_len)
    fc_ne = (ne[0] + fillet_dist * bis_ne[0], ne[1] + fillet_dist * bis_ne[1])
    # f_ne_base: project fillet center onto base line (through tbl_bc along _iw1_n_al)
    _fc_ne_d = ((fc_ne[0] - tbl_bc[0]) * _iw1_n_al[0] +
                (fc_ne[1] - tbl_bc[1]) * _iw1_n_al[1])
    f_ne_base = offset_pt(tbl_bc, _fc_ne_d, _iw1_n_al)
    v_ne = (fc_ne[0] - ne[0], fc_ne[1] - ne[1])
    t_proj = v_ne[0] * d_tang_ne[0] + v_ne[1] * d_tang_ne[1]
    f_ne_tang = (ne[0] + t_proj * d_tang_ne[0], ne[1] + t_proj * d_tang_ne[1])

    # NW fillet by symmetry
    f_nw_base = _sym(f_ne_base)
    f_nw_tang = _sym(f_ne_tang)

    # SVG radii
    apex_r_svg = abs(to_svg(apex_r, 0)[0] - to_svg(0, 0)[0])
    fillet_r_svg = abs(to_svg(fillet_r, 0)[0] - to_svg(0, 0)[0])

    # Build path (sweep-flag=1 for convex outward arcs)
    s = lambda p: to_svg(*p)
    path_d = (
        f'M {s(f_nw_base)[0]:.2f},{s(f_nw_base)[1]:.2f} '
        f'L {s(f_ne_base)[0]:.2f},{s(f_ne_base)[1]:.2f} '
        f'A {fillet_r_svg:.2f},{fillet_r_svg:.2f} 0 0 1 '
        f'{s(f_ne_tang)[0]:.2f},{s(f_ne_tang)[1]:.2f} '
        f'L {s(t_right)[0]:.2f},{s(t_right)[1]:.2f} '
        f'A {apex_r_svg:.2f},{apex_r_svg:.2f} 0 0 1 '
        f'{s(t_left)[0]:.2f},{s(t_left)[1]:.2f} '
        f'L {s(f_nw_tang)[0]:.2f},{s(f_nw_tang)[1]:.2f} '
        f'A {fillet_r_svg:.2f},{fillet_r_svg:.2f} 0 0 1 '
        f'{s(f_nw_base)[0]:.2f},{s(f_nw_base)[1]:.2f} Z'
    )

    href_dining = 'https://www.homedepot.com/pep/NEW-CLASSIC-HOME-FURNISHINGS-New-Classic-Furniture-Oscar-3-Piece-Wood-Top-Triangle-Dining-Set-Walnut-40-1651-D2C/327836175'
    out.append(f'<a href="{href_dining}" target="_blank">')
    out.append(f'<path d="{path_d}" fill="{APPL_FILL}" '
               f'stroke="{APPL_STROKE}" stroke-width="{APPL_SW}"/>')
    out.append('</a>')

    # Two chairs: 18" x 21", one on each equal (tangent) side
    ch_short = 18.0 / 12.0   # width along side
    ch_long = 21.0 / 12.0    # depth perpendicular to side
    chair_gap = 2.0 / 12.0

    for side_start, side_end in [(f_ne_tang, t_right), (t_left, f_nw_tang)]:
        mid_e = (side_start[0] + side_end[0]) / 2
        mid_n = (side_start[1] + side_end[1]) / 2
        se_d = (side_end[0] - side_start[0], side_end[1] - side_start[1])
        sl = math.sqrt(se_d[0]**2 + se_d[1]**2)
        su = (se_d[0] / sl, se_d[1] / sl)
        # Outward normal (away from table center)
        sn = (-su[1], su[0])
        to_ctr = (tbl_bc[0] - mid_e, tbl_bc[1] - mid_n)
        if sn[0] * to_ctr[0] + sn[1] * to_ctr[1] > 0:
            sn = (-sn[0], -sn[1])
        cc_e = mid_e + sn[0] * (ch_long / 2 + chair_gap)
        cc_n = mid_n + sn[1] * (ch_long / 2 + chair_gap)
        corners = []
        for ds, dn in [(-1, -1), (1, -1), (1, 1), (-1, 1)]:
            ce = cc_e + su[0] * ds * ch_short / 2 + sn[0] * dn * ch_long / 2
            cn = cc_n + su[1] * ds * ch_short / 2 + sn[1] * dn * ch_long / 2
            corners.append(to_svg(ce, cn))
        ch_svg = " ".join(f'{p[0]:.1f},{p[1]:.1f}' for p in corners)
        out.append(f'<a href="{href_dining}" target="_blank">')
        out.append(f'<polygon points="{ch_svg}" fill="{APPL_FILL}" '
                   f'stroke="{APPL_STROKE}" stroke-width="{APPL_SW}"/>')
        out.append('</a>')

    # North wall counter: along north wall from IW2 east face
    if not minik:
        _nc_c = [_nwp(_iw2_d, NORTH_CTR_DEPTH), _nwp(_iw2_d + NORTH_CTR_LENGTH, NORTH_CTR_DEPTH),
                 _nwp(_iw2_d + NORTH_CTR_LENGTH, 0), _nwp(_iw2_d, 0)]
        _appl_poly(out, _nc_c, to_svg, label="COUNTER", font_size="6",
                   href="https://www.webstaurantstore.com/regency-spec-line-30-x-36-14-gauge-stainless-steel-commercial-work-table-with-4-backsplash-and-undershelf/600TSSB3036S.html")

        # Coffee maker on counter (7.2" along wall x 9.2" inward)
        cm_ew = 7.2 / 12.0
        cm_ns = 9.2 / 12.0
        _cm_d = _iw2_d + NORTH_CTR_LENGTH - 2.0 / 12.0 - cm_ew
        _cm_i = 2.0 / 12.0
        _cm_c = [_nwp(_cm_d, _cm_i + cm_ns), _nwp(_cm_d + cm_ew, _cm_i + cm_ns),
                 _nwp(_cm_d + cm_ew, _cm_i), _nwp(_cm_d, _cm_i)]
        _appl_poly(out, _cm_c, to_svg, label="C", font_size="5",
                   href="https://www.amazon.com/Holstein-Housewares-HH-0914701E-5-Cup-Coffee/dp/B08HSRCC4T/?th=1")

def _render_furniture(out, data, layout, minik=False, db=False):
    """Render furniture: bed, loveseat/sofa, ET, chair, ottoman, room labels."""
    pts = data.pts
    to_svg = data.to_svg
    # Reference wall direction vectors for rotation-invariant placement
    w12w13_al, _ = seg_vecs(pts["W12"], pts["W13"])
    w9w10_al, w9w10_in = seg_vecs(pts["W9"], pts["W10"])
    w11w12_al, w11w12_in = seg_vecs(pts["W11"], pts["W12"])
    # IW4 west face: along=south (needed for corner intersection)
    _iw4_w_al, _ = seg_vecs(layout.iw4.poly[3], layout.iw4.poly[0])
    # IW1 north face: along=east, CW-inward=south, outward=north
    _iw1_n_al, _iw1_n_cw = seg_vecs(layout.iw1.poly[3], layout.iw1.poly[2])
    _iw1_n_out = (-_iw1_n_cw[0], -_iw1_n_cw[1])
    _iw1_w = (-_iw1_n_al[0], -_iw1_n_al[1])  # westward along IW1
    # IW4/IW1 corner: intersection of IW4 west face and IW1 north face
    _iw41_corner = line_isect(layout.iw4.poly[3], _iw4_w_al,
                              layout.iw1.poly[3], _iw1_n_al)
    def _lwp(d_w=0, d_n=0):
        """Point offset from IW4/IW1 corner: d_w westward along IW1, d_n northward from IW1."""
        return offset_pt(offset_pt(_iw41_corner, d_w, _iw1_w),
                         d_n, _iw1_n_out)
    # North wall helper: from W9 anchor along W9-W10 wall
    _nw_anchor = pts["W9"]
    def _nwp(d_along, d_inward=0):
        """Point from W9 along north wall: d_along toward W10, d_inward from wall."""
        return offset_pt(offset_pt(_nw_anchor, d_along, w9w10_al),
                         d_inward, w9w10_in)
    # IW2s east face distance from W9 along north wall
    _iw2_d = ((layout.iw2s.poly[2][0] - _nw_anchor[0]) * w9w10_al[0] +
              (layout.iw2s.poly[2][1] - _nw_anchor[1]) * w9w10_al[1])

    # Bed (rotated polygon)
    _bp = layout.bed.poly
    _bp_svg = " ".join(f"{to_svg(*p)[0]:.1f},{to_svg(*p)[1]:.1f}" for p in _bp)
    out.append(f'<polygon points="{_bp_svg}" fill="{APPL_FILL}" stroke="{APPL_STROKE}" stroke-width="{APPL_SW}"/>')
    _bp_cx = sum(p[0] for p in _bp) / 4
    _bp_cy = sum(p[1] for p in _bp) / 4
    # Rotate label 90° CW from bed long axis (now parallel to W18-W1)
    _bed_dx = to_svg(*_bp[2])[0] - to_svg(*_bp[1])[0]
    _bed_dy = to_svg(*_bp[2])[1] - to_svg(*_bp[1])[1]
    _bed_ang = math.degrees(math.atan2(_bed_dy, _bed_dx)) + 90
    _bsx, _bsy = to_svg(_bp_cx, _bp_cy)
    out.append(f'<text x="{_bsx:.1f}" y="{_bsy:.1f}" text-anchor="middle" dominant-baseline="central" font-family="Arial"'
               f' font-size="7" fill="{APPL_STROKE}" transform="rotate({_bed_ang:.1f},{_bsx:.1f},{_bsy:.1f})">'
               f'KING BED</text>')

    # Dresser (34" E-W × 19" N-S) with 15" clearance zone on south side
    d = layout.dresser
    _appl_poly(out, d.poly, to_svg, label="DRESSER")
    # Dashed clearance polygon: 15" outward from south face
    _d_al, _d_out = seg_vecs(d.poly[0], d.poly[1])  # SW→SE, outward from CCW poly
    cl_sw = offset_pt(d.poly[0], 15.0 / 12.0, _d_out)
    cl_se = offset_pt(d.poly[1], 15.0 / 12.0, _d_out)
    _appl_poly(out, [d.poly[0], d.poly[1], cl_se, cl_sw], to_svg,
               dash=True, fill_color="none")

    # Shelves (KALLAX, linked to product page)
    _appl_poly(out, layout.shelves.poly, to_svg, label="SHELVES",
               href="https://www.ikea.com/us/en/p/kallax-shelving-unit-with-underframe-white-stained-oak-effect-black-s49442718/")

    if minik:
        # SOFA: 80.75" E-W x 34.625" N-S, centered on old sofa, 2" N of IW1
        _sofa_ew = 80.75 / 12.0
        _sofa_ns = 34.625 / 12.0
        # Old sofa center was 6" + half-overlap east of IW4 west face
        _cx_d = -(6.0 / 12.0 + (SOFA_WIDTH - 24.0 / 12.0) / 2)  # neg = east of IW4
        sofa_nw = _lwp(_cx_d + _sofa_ew / 2, 2.0 / 12.0 + _sofa_ns)
        sofa_ne = _lwp(_cx_d - _sofa_ew / 2, 2.0 / 12.0 + _sofa_ns)
        sofa_se = _lwp(_cx_d - _sofa_ew / 2, 2.0 / 12.0)
        sofa_sw = _lwp(_cx_d + _sofa_ew / 2, 2.0 / 12.0)
        _appl_poly(out, [sofa_sw, sofa_se, sofa_ne, sofa_nw], to_svg,
                   label="SOFA", href="https://www.ikea.com/us/en/p/saltsjoebaden-3-seat-sofa-gunnared-light-green-s89599953/")

        # ROCKER: midpoint between ICE SE corner and SOFA NW corner
        # Recompute ICE SE via north wall chain from IW2
        _st_d = _iw2_d + NORTH_CTR_LENGTH + KITCHEN_APPL_GAP
        _ks_d = _st_d + STOVE_WIDTH + KITCHEN_APPL_GAP
        _fr_mk_d = _ks_d + KITCHEN_SINK_WIDTH + 3.0 / 12.0
        _ice_d = _fr_mk_d + MINIK_FRIDGE_W + 3.0 / 12.0
        _ice_se = _nwp(_ice_d + ICE_WIDTH, 3.0 / 12.0 + ICE_DEPTH)
        _rk_mid = ((sofa_nw[0] + _ice_se[0]) / 2,
                    (sofa_nw[1] + _ice_se[1]) / 2)
        rk_center = offset_pt(_rk_mid, 18.0 / 12.0, w9w10_in)
        rk_cx, rk_cy = rk_center
        rk_hw = ROCKER_DEPTH / 2   # half short side
        rk_hh = ROCKER_WIDTH / 2   # half long side (along w12w13_al)
        _rk_cr = (-w12w13_al[1], w12w13_al[0])  # cross direction (90° CCW)
        rk_poly = [
            (rk_cx - rk_hh * w12w13_al[0] - rk_hw * _rk_cr[0],
             rk_cy - rk_hh * w12w13_al[1] - rk_hw * _rk_cr[1]),
            (rk_cx - rk_hh * w12w13_al[0] + rk_hw * _rk_cr[0],
             rk_cy - rk_hh * w12w13_al[1] + rk_hw * _rk_cr[1]),
            (rk_cx + rk_hh * w12w13_al[0] + rk_hw * _rk_cr[0],
             rk_cy + rk_hh * w12w13_al[1] + rk_hw * _rk_cr[1]),
            (rk_cx + rk_hh * w12w13_al[0] - rk_hw * _rk_cr[0],
             rk_cy + rk_hh * w12w13_al[1] - rk_hw * _rk_cr[1]),
        ]
        _appl_poly(out, rk_poly, to_svg, label="ROCKER", font_size="6",
                   href="https://www.ikea.com/us/en/p/poaeng-rocking-chair-brown-gunnared-beige-s39502048/",
                   text_rot=_svg_angle(w12w13_al) - 90)
    elif db:
        # Shelves2 (KALLAX copy, east side of IW9, same gaps from IW9 and IW1)
        _iw9_e_al, _iw9_e_cw = seg_vecs(layout.iw9.poly[1], layout.iw9.poly[2])
        _sh2_nw = line_isect(
            offset_pt(layout.iw1.poly[0], SHELVES_GAP_IW1, w9w10_in), w9w10_al,
            offset_pt(layout.iw9.poly[2], SHELVES_GAP_IW9, _iw9_e_cw), _iw9_e_al)
        _sh2_ne = offset_pt(_sh2_nw, SHELVES_LENGTH, _iw9_e_cw)
        _sh2_sw = offset_pt(_sh2_nw, SHELVES_DEPTH, w9w10_in)
        _sh2_se = offset_pt(_sh2_ne, SHELVES_DEPTH, w9w10_in)
        _appl_poly(out, [_sh2_sw, _sh2_se, _sh2_ne, _sh2_nw], to_svg,
                   label="SHELVES",
                   href="https://www.ikea.com/us/en/p/kallax-shelving-unit-with-underframe-white-stained-oak-effect-black-s49442718/")

        # DB variant: no loveseats; ET shifted east to 2" from W-series wall
        et_r = (ET_RADIUS_CM / 2.54) / 12.0
        _et_from_iw1 = offset_pt(layout.iw1.poly[3], STD_GAP + et_r, _iw1_n_out)
        # Ray-polygon intersection along _iw1_n_al to find east wall
        _et_t_max = 0
        for _i in range(len(data.inner_poly)):
            _j = (_i + 1) % len(data.inner_poly)
            _dx = data.inner_poly[_j][0] - data.inner_poly[_i][0]
            _dy = data.inner_poly[_j][1] - data.inner_poly[_i][1]
            _det = _iw1_n_al[0] * _dy - _iw1_n_al[1] * _dx
            if abs(_det) < 1e-12:
                continue
            _ox = data.inner_poly[_i][0] - _et_from_iw1[0]
            _oy = data.inner_poly[_i][1] - _et_from_iw1[1]
            _t = (_ox * _dy - _oy * _dx) / _det
            _s = (_ox * _iw1_n_al[1] - _oy * _iw1_n_al[0]) / _det
            if 0 <= _s <= 1 and _t > 0 and _t > _et_t_max:
                _et_t_max = _t
        et_cx, et_cy = offset_pt(_et_from_iw1, _et_t_max - STD_GAP - et_r, _iw1_n_al)

        # ET: 50cm diameter endtable
        et_sx, et_sy = to_svg(et_cx, et_cy)
        et_r_svg = abs(to_svg(et_r, 0)[0] - to_svg(0, 0)[0])
        out.append('<a href="https://www.ikea.com/us/en/p/listerby-side-table-oak-veneer-30515314/" target="_blank">')
        out.append(f'<circle cx="{et_sx:.1f}" cy="{et_sy:.1f}" r="{et_r_svg:.1f}"'
                   f' fill="{APPL_FILL}" stroke="{APPL_STROKE}" stroke-width="{APPL_SW}"/>')
        out.append(f'<text x="{et_sx:.1f}" y="{et_sy+3:.1f}" text-anchor="middle"'
                   f' font-family="Arial" font-size="6" fill="{APPL_STROKE}">ET</text>')
        out.append('</a>')

        # DAYBED: 86" along IW1 x 43" perpendicular, STD_GAP N of IW1, 3" W of ET
        db_ew = 86.0 / 12.0
        db_ns = 43.0 / 12.0
        _db_s_ref = offset_pt(layout.iw1.poly[3], STD_GAP, _iw1_n_out)
        # Project ET center onto _iw1_n_al from _db_s_ref to find east edge distance
        _neg_al = (-_iw1_n_al[0], -_iw1_n_al[1])
        db_se = offset_pt((et_cx, et_cy), et_r + 3.0 / 12.0, _neg_al)
        # Adjust to daybed's south edge height (STD_GAP from IW1, not STD_GAP+et_r)
        _db_se_proj_out = ((db_se[0] - _db_s_ref[0]) * _iw1_n_out[0] +
                           (db_se[1] - _db_s_ref[1]) * _iw1_n_out[1])
        db_se = offset_pt(db_se, -_db_se_proj_out, _iw1_n_out)
        db_sw = offset_pt(db_se, -db_ew, _iw1_n_al)
        db_ne = offset_pt(db_se, db_ns, _iw1_n_out)
        db_nw = offset_pt(db_sw, db_ns, _iw1_n_out)
        _appl_poly(out, [db_sw, db_se, db_ne, db_nw], to_svg, label="DAYBED")

        # ET west: 6" W of DAYBED along _iw1_n_al, et_r S of north edge along _iw1_n_out
        _neg_out = (-_iw1_n_out[0], -_iw1_n_out[1])
        et2_center = offset_pt(offset_pt(db_nw, -(6.0 / 12.0 + et_r), _iw1_n_al),
                               et_r, _neg_out)
        et2_cx, et2_cy = et2_center
        et2_sx, et2_sy = to_svg(et2_cx, et2_cy)
        out.append('<a href="https://www.ikea.com/us/en/p/listerby-side-table-oak-veneer-30515314/" target="_blank">')
        out.append(f'<circle cx="{et2_sx:.1f}" cy="{et2_sy:.1f}" r="{et_r_svg:.1f}"'
                   f' fill="{APPL_FILL}" stroke="{APPL_STROKE}" stroke-width="{APPL_SW}"/>')
        out.append(f'<text x="{et2_sx:.1f}" y="{et2_sy+3:.1f}" text-anchor="middle"'
                   f' font-family="Arial" font-size="6" fill="{APPL_STROKE}">ET</text>')
        out.append('</a>')

        # ROCKER: center between DAYBED and RO1 / IW1 and fridge door arc
        _ro1_e_pt = _nwp(_iw2_d + RO1_OFFSET_FROM_IW2 + IW1_RO_WIDTH, 0)
        _ref = layout.iw1.poly[3]
        _db_sw_d_al = (db_sw[0] - _ref[0]) * _iw1_n_al[0] + (db_sw[1] - _ref[1]) * _iw1_n_al[1]
        _ro1_d_al = (_ro1_e_pt[0] - _ref[0]) * _iw1_n_al[0] + (_ro1_e_pt[1] - _ref[1]) * _iw1_n_al[1]
        _rk_d_al = (_db_sw_d_al + _ro1_d_al) / 2 - 8.0 / 12.0
        _fr_s_pt = _nwp(0, 3.0 / 12.0 + 35.0 / 12.0)
        _fr_door_s_pt = offset_pt(_fr_s_pt, 32.75 / 12.0, w9w10_in)
        _fr_d_out = (_fr_door_s_pt[0] - _ref[0]) * _iw1_n_out[0] + (_fr_door_s_pt[1] - _ref[1]) * _iw1_n_out[1]
        _rk_d_out = _fr_d_out / 2 + 26.0 / 12.0
        rk_center = offset_pt(offset_pt(_ref, _rk_d_al, _iw1_n_al), _rk_d_out, _iw1_n_out)
        rk_cx, rk_cy = rk_center
        rk_hw = ROCKER_DEPTH / 2   # half short side
        rk_hh = ROCKER_WIDTH / 2   # half long side (along w12w13_al)
        _rk_cr = (-w12w13_al[1], w12w13_al[0])
        rk_poly = [
            (rk_cx - rk_hh * w12w13_al[0] - rk_hw * _rk_cr[0],
             rk_cy - rk_hh * w12w13_al[1] - rk_hw * _rk_cr[1]),
            (rk_cx - rk_hh * w12w13_al[0] + rk_hw * _rk_cr[0],
             rk_cy - rk_hh * w12w13_al[1] + rk_hw * _rk_cr[1]),
            (rk_cx + rk_hh * w12w13_al[0] + rk_hw * _rk_cr[0],
             rk_cy + rk_hh * w12w13_al[1] + rk_hw * _rk_cr[1]),
            (rk_cx + rk_hh * w12w13_al[0] - rk_hw * _rk_cr[0],
             rk_cy + rk_hh * w12w13_al[1] - rk_hw * _rk_cr[1]),
        ]
        _appl_poly(out, rk_poly, to_svg, label="ROCKER", font_size="6",
                   href="https://www.ikea.com/us/en/p/poaeng-rocking-chair-brown-gunnared-beige-s39502048/",
                   text_rot=_svg_angle(w12w13_al) - 90)
    else:
        # Loveseat: 35" wide x 65" long, long side along w12w13_al
        lv_width = LOVESEAT_WIDTH
        lv_height = LOVESEAT_LENGTH
        lv_nw = _lwp(LOVESEAT_OFFSET_IW4, LOVESEAT_OFFSET_IW1)
        _lv_perp = (-w12w13_al[1], w12w13_al[0])  # 90° CCW from w12w13_al
        lv_sw_pt = offset_pt(lv_nw, lv_height, w12w13_al)
        lv_ne_pt = offset_pt(lv_nw, lv_width, _lv_perp)
        lv_se_pt = offset_pt(lv_sw_pt, lv_width, _lv_perp)

        # ET position: STD_GAP N of IW1, STD_GAP from loveseat SE corner
        et_r = (ET_RADIUS_CM / 2.54) / 12.0
        et_gap = et_r + STD_GAP
        # Circle-line intersection along IW1 direction (rotation-safe)
        _et_from_iw1 = offset_pt(layout.iw1.poly[3], STD_GAP + et_r, _iw1_n_out)
        _et_dx = _et_from_iw1[0] - lv_se_pt[0]
        _et_dy = _et_from_iw1[1] - lv_se_pt[1]
        _et_b = 2 * (_et_dx * _iw1_n_al[0] + _et_dy * _iw1_n_al[1])
        _et_c = _et_dx**2 + _et_dy**2 - et_gap**2
        _et_t = (-_et_b + math.sqrt(_et_b**2 - 4 * _et_c)) / 2
        et_cx = _et_from_iw1[0] + _et_t * _iw1_n_al[0]
        et_cy = _et_from_iw1[1] + _et_t * _iw1_n_al[1]

        _appl_poly(out, [lv_sw_pt, lv_se_pt, lv_ne_pt, lv_nw], to_svg,
                   label="LOVESEAT", font_size="6",
                   href="https://www.ikea.com/us/en/p/saltsjoebaden-loveseat-tonerud-red-brown-s59579188/",
                   text_rot=_svg_angle(w12w13_al) - 90)

        # ET: 50cm diameter endtable
        et_sx, et_sy = to_svg(et_cx, et_cy)
        et_r_svg = abs(to_svg(et_r, 0)[0] - to_svg(0, 0)[0])
        out.append('<a href="https://www.ikea.com/us/en/p/listerby-side-table-oak-veneer-30515314/" target="_blank">')
        out.append(f'<circle cx="{et_sx:.1f}" cy="{et_sy:.1f}" r="{et_r_svg:.1f}"'
                   f' fill="{APPL_FILL}" stroke="{APPL_STROKE}" stroke-width="{APPL_SW}"/>')
        out.append(f'<text x="{et_sx:.1f}" y="{et_sy+3:.1f}" text-anchor="middle"'
                   f' font-family="Arial" font-size="6" fill="{APPL_STROKE}">ET</text>')
        out.append('</a>')

        # LOVESEAT2: same as LOVESEAT but long side along IW1, 35" across
        # SW corner: ET east edge + gap along IW1, offset to IW1 face
        _lv2_sw = offset_pt(
            offset_pt((et_cx, et_cy), et_r + STD_GAP, _iw1_n_al),
            -et_r, _iw1_n_out)
        _lv2_nw = offset_pt(_lv2_sw, lv_width, _iw1_n_out)
        _lv2_se = offset_pt(_lv2_sw, lv_height, _iw1_n_al)
        _lv2_ne = offset_pt(_lv2_nw, lv_height, _iw1_n_al)
        _lv2_pts = [_lv2_sw, _lv2_se, _lv2_ne, _lv2_nw]
        _lv2_svg = " ".join(f"{to_svg(*p)[0]:.1f},{to_svg(*p)[1]:.1f}" for p in _lv2_pts)
        _lv2_cx = sum(p[0] for p in _lv2_pts) / 4
        _lv2_cy = sum(p[1] for p in _lv2_pts) / 4
        _lv2_scx, _lv2_scy = to_svg(_lv2_cx, _lv2_cy)
        out.append('<a href="https://www.ikea.com/us/en/p/saltsjoebaden-loveseat-tonerud-red-brown-s59579188/" target="_blank">')
        out.append(f'<polygon points="{_lv2_svg}"'
                   f' fill="{APPL_FILL}" stroke="{APPL_STROKE}" stroke-width="{APPL_SW}"/>')
        out.append(f'<text x="{_lv2_scx:.1f}" y="{_lv2_scy+3:.1f}" text-anchor="middle" font-family="Arial"'
                   f' font-size="6" fill="{APPL_STROKE}">LOVESEAT</text>')
        out.append('</a>')

    # CHAIR: 32" x 37", 45° to W12-W13 wall, centered between W11 and W12
    _ch_theta = math.atan2(w12w13_al[1], w12w13_al[0]) - math.pi / 4
    _ch_along = (math.cos(_ch_theta), math.sin(_ch_theta))  # facing direction
    _ch_cross = (-math.sin(_ch_theta), math.cos(_ch_theta))  # width direction
    # Base position: midpoint of W11-W12 chord, offset 1" back along chord and 8" inward
    _ch_mid = ((pts["W11"][0] + pts["W12"][0]) / 2,
               (pts["W11"][1] + pts["W12"][1]) / 2)
    _ch_base = offset_pt(offset_pt(_ch_mid, -1.0 / 12.0, w11w12_al),
                         8.0 / 12.0, w11w12_in)
    # 4" offset in chair facing direction
    ch_cx, ch_cy = offset_pt(_ch_base, 4.0 / 12.0, _ch_along)
    _ch_hw = CHAIR_WIDTH / 2
    _ch_hd = CHAIR_DEPTH / 2
    ch_poly = [
        (ch_cx - _ch_hw * _ch_cross[0] - _ch_hd * _ch_along[0],
         ch_cy - _ch_hw * _ch_cross[1] - _ch_hd * _ch_along[1]),
        (ch_cx + _ch_hw * _ch_cross[0] - _ch_hd * _ch_along[0],
         ch_cy + _ch_hw * _ch_cross[1] - _ch_hd * _ch_along[1]),
        (ch_cx + _ch_hw * _ch_cross[0] + _ch_hd * _ch_along[0],
         ch_cy + _ch_hw * _ch_cross[1] + _ch_hd * _ch_along[1]),
        (ch_cx - _ch_hw * _ch_cross[0] + _ch_hd * _ch_along[0],
         ch_cy - _ch_hw * _ch_cross[1] + _ch_hd * _ch_along[1]),
    ]
    _ch_svg_deg = _svg_angle(w12w13_al) - 45
    _appl_poly(out, ch_poly, to_svg, label="CHAIR", font_size="6",
               href="https://www.ikea.com/us/en/p/havberg-swivel-easy-chair-and-footstool-grann-bomstad-golden-brown-s59485321/",
               text_rot=_ch_svg_deg)

    # OTTO: 29" x 29", same rotation, offset from CHAIR center in facing direction
    ot_dist = 39.0 / 12.0  # ch half-h 18.5" + 6" gap + ot half 14.5"
    ot_cx, ot_cy = offset_pt((ch_cx, ch_cy), ot_dist, _ch_along)
    _ot_hs = OTTOMAN_SIZE / 2
    ot_poly = [
        (ot_cx - _ot_hs * _ch_cross[0] - _ot_hs * _ch_along[0],
         ot_cy - _ot_hs * _ch_cross[1] - _ot_hs * _ch_along[1]),
        (ot_cx + _ot_hs * _ch_cross[0] - _ot_hs * _ch_along[0],
         ot_cy + _ot_hs * _ch_cross[1] - _ot_hs * _ch_along[1]),
        (ot_cx + _ot_hs * _ch_cross[0] + _ot_hs * _ch_along[0],
         ot_cy + _ot_hs * _ch_cross[1] + _ot_hs * _ch_along[1]),
        (ot_cx - _ot_hs * _ch_cross[0] + _ot_hs * _ch_along[0],
         ot_cy - _ot_hs * _ch_cross[1] + _ot_hs * _ch_along[1]),
    ]
    _appl_poly(out, ot_poly, to_svg, label="OTTO", font_size="6",
               href="https://www.ikea.com/us/en/p/havberg-swivel-easy-chair-and-footstool-grann-bomstad-golden-brown-s59485321/",
               text_rot=_ch_svg_deg)

    # DESK: 60" x 30", along W16-W17 wall (interior side)
    w16w17_al, w16w17_in = seg_vecs(pts["W16"], pts["W17"])
    _neg_w16w17_al = (-w16w17_al[0], -w16w17_al[1])
    dk_sw_pt = pts["W17"]
    dk_se_pt = offset_pt(dk_sw_pt, DESK_WIDTH, _neg_w16w17_al)
    dk_nw_pt = offset_pt(dk_sw_pt, DESK_DEPTH, w16w17_in)
    dk_ne_pt = offset_pt(dk_se_pt, DESK_DEPTH, w16w17_in)
    _appl_poly(out, [dk_sw_pt, dk_se_pt, dk_ne_pt, dk_nw_pt], to_svg,
               label="DESK", text_rot=_svg_angle(_neg_w16w17_al))

    # DESK CHAIR: 27" x 24", 12" in front of desk, centered
    dc_sw_pt = offset_pt(offset_pt(dk_sw_pt, DESK_WIDTH / 2 - DESK_CHAIR_WIDTH / 2, _neg_w16w17_al),
                         DESK_DEPTH + DESK_CHAIR_GAP, w16w17_in)
    dc_se_pt = offset_pt(dc_sw_pt, DESK_CHAIR_WIDTH, _neg_w16w17_al)
    dc_nw_pt = offset_pt(dc_sw_pt, DESK_CHAIR_DEPTH, w16w17_in)
    dc_ne_pt = offset_pt(dc_se_pt, DESK_CHAIR_DEPTH, w16w17_in)
    _appl_poly(out, [dc_sw_pt, dc_se_pt, dc_ne_pt, dc_nw_pt], to_svg,
               label="CHAIR",
               href="https://www.amazon.com/BESTFAIR-Ergonomic-Office-Chair-Adjustable/dp/B0FDQDMP2D?th=1",
               text_rot=_svg_angle(_neg_w16w17_al))

    # Room labels — use polygon face midpoints for rotation-safe positioning
    # BEDROOM: right edge at W end of RO1, top at N end of RO3
    _ro_list = compute_rough_openings(pts, layout)
    _ro1_bd = [r for r in _ro_list if r.name == "RO1"][0].poly
    _ro3_bd = [r for r in _ro_list if r.name == "RO3"][0].poly
    _ro1_w_mid = ((_ro1_bd[0][0] + _ro1_bd[3][0]) / 2,
                  (_ro1_bd[0][1] + _ro1_bd[3][1]) / 2)
    _ro3_n_mid = ((_ro3_bd[2][0] + _ro3_bd[3][0]) / 2,
                  (_ro3_bd[2][1] + _ro3_bd[3][1]) / 2)
    bdx, bdy = to_svg(_ro1_w_mid[0], _ro3_n_mid[1])
    out.append(f'<text x="{bdx:.1f}" y="{bdy:.1f}" text-anchor="end" dominant-baseline="hanging"'
               f' font-family="Arial" font-size="8" fill="#666">BEDROOM</text>')

    # UTIL: same northing as BEDROOM, centered in easting between south toilet and sink
    _iw8_al, _iw8_in = seg_vecs(layout.iw8.poly[0], layout.iw8.poly[1])
    _iw8_s_ref = layout.iw8.poly[0]
    _dryer_cx = sum(p[0] for p in layout.dryer.poly) / 4
    _dryer_cy = sum(p[1] for p in layout.dryer.poly) / 4
    _d_dryer_al = ((_dryer_cx - _iw8_s_ref[0]) * _iw8_al[0] +
                   (_dryer_cy - _iw8_s_ref[1]) * _iw8_al[1])
    _toilet_s = offset_pt(_iw8_s_ref, _d_dryer_al, _iw8_al)
    _ctr_cx = sum(p[0] for p in layout.ctr.poly) / 4
    _ctr_cy = sum(p[1] for p in layout.ctr.poly) / 4
    _sink_mid = ((_dryer_cx + _ctr_cx) / 2, (_dryer_cy + _ctr_cy) / 2)
    _d_sink_al = ((_sink_mid[0] - _iw8_s_ref[0]) * _iw8_al[0] +
                  (_sink_mid[1] - _iw8_s_ref[1]) * _iw8_al[1])
    _sink_s = offset_pt(_iw8_s_ref, _d_sink_al, _iw8_al)
    _sk_s = offset_pt(_sink_s, SINK_RY, _iw8_in)
    _util_e = (_toilet_s[0] + _sk_s[0]) / 2
    _util_n = _ro3_n_mid[1]
    utx, uty = to_svg(_util_e, _util_n)
    out.append(f'<text x="{utx:.1f}" y="{uty:.1f}" text-anchor="middle" dominant-baseline="hanging"'
               f' font-family="Arial" font-size="8" fill="#666">UTIL</text>')

    # KITCHEN: centered beneath kitchen sink, just above dim02 (25' 8.1") line
    _w9w10_al, _ = seg_vecs(pts["W9"], pts["W10"])
    _iw2s_ne_k = layout.iw2s.poly[2]
    _iw2_d_k = ((_iw2s_ne_k[0] - pts["W9"][0]) * _w9w10_al[0] +
                (_iw2s_ne_k[1] - pts["W9"][1]) * _w9w10_al[1])
    _st_d_k = _iw2_d_k + NORTH_CTR_LENGTH + KITCHEN_APPL_GAP
    _ks_d_k = _st_d_k + STOVE_WIDTH + KITCHEN_APPL_GAP + 2.0 / 12.0
    _sink_ctr = offset_pt(pts["W9"], _ks_d_k + KITCHEN_SINK_WIDTH / 2, _w9w10_al)
    _dim02_n = (pts["F12"][1] + pts["F13"][1]) / 2
    kx, ky = to_svg(_sink_ctr[0], _dim02_n + 3.0 / 12.0)
    out.append(f'<text x="{kx:.1f}" y="{ky:.1f}" text-anchor="middle"'
               f' font-family="Arial" font-size="8" fill="#666">KITCHEN</text>')

    # LIVING: centered under O6 at KITCHEN label's northing
    _o6_open = compute_outer_openings(pts, layout)
    _o6 = [o for o in _o6_open if o.name == "O6"][0]
    _o6_cx = sum(p[0] for p in _o6.poly) / 4
    lx, ly = to_svg(_o6_cx, _dim02_n + 3.0 / 12.0)
    out.append(f'<text x="{lx:.1f}" y="{ly:.1f}" text-anchor="middle"'
               f' font-family="Arial" font-size="8" fill="#666">LIVING</text>')

    # BATH: centered between IW2s west and W2-W5, at RO4 north end northing
    _bath_e = (layout.iw2s.poly[0][0] + pts["W2"][0]) / 2
    _ro4_bd = [r for r in compute_rough_openings(pts, layout) if r.name == "RO4"][0].poly
    _ro4_n_mid = ((_ro4_bd[2][0] + _ro4_bd[3][0]) / 2,
                  (_ro4_bd[2][1] + _ro4_bd[3][1]) / 2)
    bax, bay = to_svg(_bath_e, _ro4_n_mid[1])
    out.append(f'<text x="{bax:.1f}" y="{bay:.1f}" text-anchor="middle"'
               f' font-family="Arial" font-size="8" fill="#666">BATH</text>')

    # OFFICE: midpoint between IW4 east face and W15, vertically between ctr+5'+3" and IW1
    _iw4_e_mid = ((layout.iw4.poly[1][0] + layout.iw4.poly[2][0]) / 2,
                  (layout.iw4.poly[1][1] + layout.iw4.poly[2][1]) / 2)
    _of_ew = ((_iw4_e_mid[0] + pts["W15"][0]) / 2,
              (_iw4_e_mid[1] + pts["W15"][1]) / 2)
    # N-S: offset ctr south face by 5'+WALL_3IN, midpoint with IW1 south, adjust by -2'+8"
    _ctr_s_mid = ((layout.ctr.poly[0][0] + layout.ctr.poly[1][0]) / 2,
                  (layout.ctr.poly[0][1] + layout.ctr.poly[1][1]) / 2)
    _iw1_s_mid_r = ((layout.iw1.poly[0][0] + layout.iw1.poly[1][0]) / 2,
                    (layout.iw1.poly[0][1] + layout.iw1.poly[1][1]) / 2)
    _ctr_offset = offset_pt(_ctr_s_mid, 5.0 + WALL_3IN, _iw1_n_out)
    _of_ns = ((_ctr_offset[0] + _iw1_s_mid_r[0]) / 2,
              (_ctr_offset[1] + _iw1_s_mid_r[1]) / 2)
    _of_ns_adj = offset_pt(_of_ns, -2.0 + 26.0 / 12.0, _iw1_n_out)
    # Project onto _iw1_n_out from _of_ew
    _of_ns_d = ((_of_ns_adj[0] - _of_ew[0]) * _iw1_n_out[0] +
                (_of_ns_adj[1] - _of_ew[1]) * _iw1_n_out[1])
    of_cx = _of_ew[0] + _of_ns_d * _iw1_n_out[0]
    of_cy = _of_ew[1] + _of_ns_d * _iw1_n_out[1]
    ofx, ofy = to_svg(of_cx, of_cy)
    out.append(f'<text x="{ofx:.1f}" y="{ofy+3:.1f}" text-anchor="middle" font-family="Arial"'
               f' font-size="8" fill="#666">OFFICE</text>')


def _render_dimensions(out, data, layout, bare=False):
    """Render all dimension lines (interior and external).

    All endpoints computed wall-relative via compute_dimension_endpoints.
    """
    to_svg = data.to_svg
    ep = {name: pt for name, pt in compute_dimension_endpoints(
        data.pts, layout, data.radii, bare=bare)}

    def _edist(a, b):
        """Euclidean distance between two endpoint tuples."""
        return math.sqrt((b[0] - a[0])**2 + (b[1] - a[1])**2)

    # dim01: IW1-north → W9
    _rotated_dim(out, ep["dim01_A"], ep["dim01_B"],
                 fmt_dist(_edist(ep["dim01_A"], ep["dim01_B"])), to_svg)

    # dim02: IW2-east → W12-W13
    _rotated_dim(out, ep["dim02_A"], ep["dim02_B"],
                 fmt_dist(_edist(ep["dim02_A"], ep["dim02_B"])), to_svg)

    # dim03: East closet (rotated)
    _rotated_dim(out, ep["dim03_A"], ep["dim03_B"],
                 f"CLOSET {fmt_dist(_edist(ep['dim03_A'], ep['dim03_B']))}", to_svg)

    # dim04: West closet (rotated)
    _rotated_dim(out, ep["dim04_A"], ep["dim04_B"],
                 f"CLOSET {fmt_dist(_edist(ep['dim04_A'], ep['dim04_B']))}", to_svg)

    # dim05: W2-W5 → IW3 west face, parallel to W18-W1 at IW3 midpoint
    _rotated_dim(out, ep["dim05_A"], ep["dim05_B"],
                 fmt_dist(_edist(ep["dim05_A"], ep["dim05_B"])), to_svg)

    # dim06: IW4-east → W14-W15 (perp to W14-W15 at O8 center)
    _rotated_dim(out, ep["dim06_A"], ep["dim06_B"],
                 fmt_dist(_edist(ep["dim06_A"], ep["dim06_B"])), to_svg)

    # dim07: Storage
    _rotated_dim(out, ep["dim07_A"], ep["dim07_B"],
                 f"STORAGE {fmt_dist(_edist(ep['dim07_A'], ep['dim07_B']))}", to_svg)

    # dim08: O1 east → RO3 west
    _rotated_dim(out, ep["dim08_A"], ep["dim08_B"],
                 fmt_dist(_edist(ep["dim08_A"], ep["dim08_B"])), to_svg)

    # dim09: W2-W5 → IW2-west
    _rotated_dim(out, ep["dim09_A"], ep["dim09_B"],
                 fmt_dist(_edist(ep["dim09_A"], ep["dim09_B"])), to_svg)

    # dim10: W2 → IW2s-west
    _rotated_dim(out, ep["dim10_A"], ep["dim10_B"],
                 fmt_dist(_edist(ep["dim10_A"], ep["dim10_B"])), to_svg)

    # dim11: IW5-south → W18
    _rotated_dim(out, ep["dim11_A"], ep["dim11_B"],
                 fmt_dist(_edist(ep["dim11_A"], ep["dim11_B"])), to_svg)

    # dim12: Office verticals
    if bare:
        _rotated_dim(out, ep["dim12bare_A"], ep["dim12bare_B"],
                     fmt_dist(_edist(ep["dim12bare_A"], ep["dim12bare_B"])), to_svg)
    else:
        _rotated_dim(out, ep["dim12a_A"], ep["dim12a_B"],
                     fmt_dist(_edist(ep["dim12a_A"], ep["dim12a_B"])), to_svg)
        _rotated_dim(out, ep["dim12b_A"], ep["dim12b_B"],
                     fmt_dist(_edist(ep["dim12b_A"], ep["dim12b_B"])), to_svg)

    # dim13: External F18 → F6
    _rotated_dim(out, ep["dim13_A"], ep["dim13_B"],
                 fmt_dist(_edist(ep["dim13_A"], ep["dim13_B"])), to_svg)

    # dim14: Arc width
    _rotated_dim(out, ep["dim14_A"], ep["dim14_B"],
                 fmt_dist(_edist(ep["dim14_A"], ep["dim14_B"])), to_svg)

    # dim15: External F2 → F15
    _rotated_dim(out, ep["dim15_A"], ep["dim15_B"],
                 fmt_dist(_edist(ep["dim15_A"], ep["dim15_B"])), to_svg)

    # dim17: O10 → IW1-south (rotated)
    _rotated_dim(out, ep["dim17_A"], ep["dim17_B"],
                 fmt_dist(_edist(ep["dim17_A"], ep["dim17_B"])), to_svg)

    # dim18: IW9 → IW11 (rotated)
    _rotated_dim(out, ep["dim18_A"], ep["dim18_B"],
                 fmt_dist(_edist(ep["dim18_A"], ep["dim18_B"])), to_svg)

    # dim19: O11 → IW8-south
    _rotated_dim(out, ep["dim19_A"], ep["dim19_B"],
                 fmt_dist(_edist(ep["dim19_A"], ep["dim19_B"])), to_svg)

    # dim22: IW12-north-mid → IW5-south (rotated)
    _rotated_dim(out, ep["dim22_A"], ep["dim22_B"],
                 fmt_dist(_edist(ep["dim22_A"], ep["dim22_B"])), to_svg)

    # Bare-only dimensions
    if bare:
        # dim20: IW2-east → W14
        _rotated_dim(out, ep["dim20_A"], ep["dim20_B"],
                     fmt_dist(_edist(ep["dim20_A"], ep["dim20_B"])), to_svg)
        # dim21: W11a-W11b mid → IW1-north
        _rotated_dim(out, ep["dim21_A"], ep["dim21_B"],
                     fmt_dist(_edist(ep["dim21_A"], ep["dim21_B"])), to_svg)


def _render_openings(out, data, layout, bare=False):
    """Render door swings and jamb blocks for O3, O6, RO1-RO7.

    Opening fill polygons are rendered by _render_walls() as part of the
    double-shell wall section loop.
    """
    pts = data.pts
    to_svg = data.to_svg
    outer_openings = compute_outer_openings(pts, layout)

    # O3 door: 30" door on F2-F5 east wall, hinged F5-side, swings east
    o3 = [o for o in outer_openings if o.name == "O3"][0]
    # O3 poly: [outer_start, outer_end, inner_end, inner_start]
    # Wall direction (F2→F5) and cross direction (outer→inner)
    _o3_os, _o3_oe = o3.poly[0], o3.poly[1]
    _o3_is, _o3_ie = o3.poly[3], o3.poly[2]
    _o3_dE = _o3_oe[0] - _o3_os[0]
    _o3_dN = _o3_oe[1] - _o3_os[1]
    _o3_len = math.sqrt(_o3_dE**2 + _o3_dN**2)
    _o3_along = (_o3_dE / _o3_len, _o3_dN / _o3_len)  # unit along wall
    _o3_cross = (_o3_along[1], -_o3_along[0])  # unit outer→inner (right of along)
    # Midline endpoints (halfway between outer and inner faces)
    _o3_ms = ((_o3_os[0] + _o3_is[0]) / 2, (_o3_os[1] + _o3_is[1]) / 2)
    _o3_me = ((_o3_oe[0] + _o3_ie[0]) / 2, (_o3_oe[1] + _o3_ie[1]) / 2)
    gap = (_o3_len - O3_DOOR_WIDTH) / 2

    # Jamb blocks (rotated rectangles along wall direction)
    _hf = DOOR_FLAT_FACE / 2
    for _jc_t in [gap / 2, _o3_len - gap / 2]:  # center of each jamb along opening
        _jc = (_o3_ms[0] + _jc_t * _o3_along[0], _o3_ms[1] + _jc_t * _o3_along[1])
        _jl = gap / 2  # half-length of jamb block along wall
        _jp = [(_jc[0] + da * _o3_along[0] + dc * _o3_cross[0],
                _jc[1] + da * _o3_along[1] + dc * _o3_cross[1])
               for da, dc in [(-_jl, -_hf), (_jl, -_hf), (_jl, _hf), (-_jl, _hf)]]
        _svg_pts = " ".join(f"{to_svg(p[0], p[1])[0]:.1f},{to_svg(p[0], p[1])[1]:.1f}" for p in _jp)
        out.append(f'<polygon points="{_svg_pts}" fill="{JAMB_COLOR}" stroke="none"/>')

    # Door: hinge at wall midline, NE end (F5-side) - 1" toward SW
    hinge = (_o3_me[0] - gap * _o3_along[0],
             _o3_me[1] - gap * _o3_along[1])
    hx, hy = to_svg(hinge[0], hinge[1])
    # Open position: door extends perpendicular to wall (into interior)
    tip = (hinge[0] + O3_DOOR_WIDTH * _o3_cross[0],
           hinge[1] + O3_DOOR_WIDTH * _o3_cross[1])
    tx, ty = to_svg(tip[0], tip[1])
    out.append(f'<line x1="{hx:.1f}" y1="{hy:.1f}" x2="{tx:.1f}" y2="{ty:.1f}"'
               f' stroke="{JAMB_COLOR}" stroke-width="1.0"/>')
    # Arc from open (cross) sweeping 90° to closed (along wall toward F2)
    arc_pts_str = _swing_arc_svg(hinge, O3_DOOR_WIDTH, _o3_cross,
                                 (-_o3_along[0], -_o3_along[1]), to_svg)
    out.append(f'<polyline points="{arc_pts_str}" fill="none"'
               f' stroke="{JAMB_COLOR}" stroke-width="0.5"/>')

    # O6 door: 42" door, hinged east, swings inward (rotation-safe)
    o6 = [o for o in outer_openings if o.name == "O6"][0]
    # Along direction (start→end on inner face)
    _o6_dx = o6.poly[1][0] - o6.poly[0][0]
    _o6_dy = o6.poly[1][1] - o6.poly[0][1]
    _o6_len = math.sqrt(_o6_dx**2 + _o6_dy**2)
    _o6_al = (_o6_dx / _o6_len, _o6_dy / _o6_len)
    # Inward direction (outer→inner at end edge)
    _o6_in_dx = o6.poly[1][0] - o6.poly[2][0]
    _o6_in_dy = o6.poly[1][1] - o6.poly[2][1]
    _o6_in_len = math.sqrt(_o6_in_dx**2 + _o6_in_dy**2)
    _o6_inward = (_o6_in_dx / _o6_in_len, _o6_in_dy / _o6_in_len)
    # Edge midpoints (wall midline at start and end)
    _o6_s_mid = ((o6.poly[0][0] + o6.poly[3][0]) / 2,
                 (o6.poly[0][1] + o6.poly[3][1]) / 2)
    _o6_e_mid = ((o6.poly[1][0] + o6.poly[2][0]) / 2,
                 (o6.poly[1][1] + o6.poly[2][1]) / 2)
    gap = (O6_WIDTH - O6_DOOR_WIDTH) / 2
    _o6_door_s = (_o6_s_mid[0] + gap * _o6_al[0], _o6_s_mid[1] + gap * _o6_al[1])
    _o6_door_e = (_o6_e_mid[0] - gap * _o6_al[0], _o6_e_mid[1] - gap * _o6_al[1])
    _jht = DOOR_FLAT_FACE / 2

    # Jamb blocks as polygons (rotation-safe)
    for _jb_s, _jb_e in [(_o6_s_mid, _o6_door_s), (_o6_door_e, _o6_e_mid)]:
        _jp = [(_jb_s[0] - _jht * _o6_inward[0], _jb_s[1] - _jht * _o6_inward[1]),
               (_jb_e[0] - _jht * _o6_inward[0], _jb_e[1] - _jht * _o6_inward[1]),
               (_jb_e[0] + _jht * _o6_inward[0], _jb_e[1] + _jht * _o6_inward[1]),
               (_jb_s[0] + _jht * _o6_inward[0], _jb_s[1] + _jht * _o6_inward[1])]
        _jp_svg = " ".join(f"{to_svg(p[0], p[1])[0]:.1f},{to_svg(p[0], p[1])[1]:.1f}"
                           for p in _jp)
        out.append(f'<polygon points="{_jp_svg}" fill="{JAMB_COLOR}" stroke="none"/>')

    # Door: hinge at end side, swings inward
    hinge = _o6_door_e
    hx, hy = to_svg(*hinge)
    tip = (hinge[0] + O6_DOOR_WIDTH * _o6_inward[0],
           hinge[1] + O6_DOOR_WIDTH * _o6_inward[1])
    tx, ty = to_svg(*tip)
    out.append(f'<line x1="{hx:.1f}" y1="{hy:.1f}" x2="{tx:.1f}" y2="{ty:.1f}"'
               f' stroke="{JAMB_COLOR}" stroke-width="1.0"/>')
    # Arc from open (inward) sweeping 90° to closed (-along direction)
    arc_pts = _swing_arc_svg(hinge, O6_DOOR_WIDTH, _o6_inward,
                             (-_o6_al[0], -_o6_al[1]), to_svg)
    out.append(f'<polyline points="{arc_pts}" fill="none"'
               f' stroke="{JAMB_COLOR}" stroke-width="0.5"/>')

    # RO1 door: 36" door, hinged east, swings through wall (rotation-safe)
    rough_openings = compute_rough_openings(pts, layout)
    ro1 = [r for r in rough_openings if r.name == "RO1"][0]
    # Along direction: poly[0]→poly[1] (SW→SE)
    _ro1_dx = ro1.poly[1][0] - ro1.poly[0][0]
    _ro1_dy = ro1.poly[1][1] - ro1.poly[0][1]
    _ro1_len = math.sqrt(_ro1_dx**2 + _ro1_dy**2)
    _ro1_al = (_ro1_dx / _ro1_len, _ro1_dy / _ro1_len)
    # Through-wall direction (NE→SE = toward south face)
    _ro1_thru_dx = ro1.poly[1][0] - ro1.poly[2][0]
    _ro1_thru_dy = ro1.poly[1][1] - ro1.poly[2][1]
    _ro1_thru_len = math.sqrt(_ro1_thru_dx**2 + _ro1_thru_dy**2)
    _ro1_swing = (_ro1_thru_dx / _ro1_thru_len, _ro1_thru_dy / _ro1_thru_len)
    # End edge midpoint (east): midpoint of poly[1] and poly[2]
    _ro1_end = ((ro1.poly[1][0] + ro1.poly[2][0]) / 2,
                (ro1.poly[1][1] + ro1.poly[2][1]) / 2)
    ro1_gap = (_ro1_len - RO1_DOOR_WIDTH) / 2
    hinge = (_ro1_end[0] - ro1_gap * _ro1_al[0],
             _ro1_end[1] - ro1_gap * _ro1_al[1])
    hx, hy = to_svg(*hinge)
    tip = (hinge[0] + RO1_DOOR_WIDTH * _ro1_swing[0],
           hinge[1] + RO1_DOOR_WIDTH * _ro1_swing[1])
    tx, ty = to_svg(*tip)
    out.append(f'<line x1="{hx:.1f}" y1="{hy:.1f}" x2="{tx:.1f}" y2="{ty:.1f}"'
               f' stroke="{JAMB_COLOR}" stroke-width="1.0"/>')
    # Arc from open (swing) sweeping 90° to closed (-along direction)
    arc_pts = _swing_arc_svg(hinge, RO1_DOOR_WIDTH, _ro1_swing,
                             (-_ro1_al[0], -_ro1_al[1]), to_svg)
    out.append(f'<polyline points="{arc_pts}" fill="none"'
               f' stroke="{JAMB_COLOR}" stroke-width="0.5"/>')

    # RO2 door: 36" door in IW11 (rotated), hinged at NNE edge, swings into office
    ro2 = [r for r in rough_openings if r.name == "RO2"][0]
    _ro2p = ro2.poly  # [SW, SE, NE, NW]
    # IW11 unit vectors (recompute from polygon)
    _i11_sw, _i11_se, _i11_ne, _i11_nw = layout.iw11.poly
    _i11_dx_n = _i11_ne[0] - _i11_se[0]; _i11_dy_n = _i11_ne[1] - _i11_se[1]
    _i11_ln = math.sqrt(_i11_dx_n**2 + _i11_dy_n**2)
    _i11_an = (_i11_dx_n / _i11_ln, _i11_dy_n / _i11_ln)  # NNE (along length)
    _i11_dx_t = _i11_sw[0] - _i11_se[0]; _i11_dy_t = _i11_sw[1] - _i11_se[1]
    _i11_lt = math.sqrt(_i11_dx_t**2 + _i11_dy_t**2)
    _i11_at = (_i11_dx_t / _i11_lt, _i11_dy_t / _i11_lt)  # SE→SW (across thickness)
    # Hinge: NNE edge center, offset inward by centering gap
    _ro2_gap = (IW4_RO_WIDTH - RO2_DOOR_WIDTH) / 2
    _ro2_n_ctr = ((_ro2p[3][0] + _ro2p[2][0]) / 2,
                  (_ro2p[3][1] + _ro2p[2][1]) / 2)
    hinge_e = _ro2_n_ctr[0] - _ro2_gap * _i11_an[0]
    hinge_n = _ro2_n_ctr[1] - _ro2_gap * _i11_an[1]
    hx, hy = to_svg(hinge_e, hinge_n)
    # Open position: door swings into office (-_i11_at direction = ENE)
    tip_e = hinge_e - RO2_DOOR_WIDTH * _i11_at[0]
    tip_n = hinge_n - RO2_DOOR_WIDTH * _i11_at[1]
    tx, ty = to_svg(tip_e, tip_n)
    out.append(f'<line x1="{hx:.1f}" y1="{hy:.1f}" x2="{tx:.1f}" y2="{ty:.1f}"'
               f' stroke="{JAMB_COLOR}" stroke-width="1.0"/>')
    # Arc from closed (SSW, -_i11_an) sweeping CCW 90° to open (ENE, -_i11_at)
    arc_pts = _swing_arc_svg((hinge_e, hinge_n), RO2_DOOR_WIDTH,
                             (-_i11_an[0], -_i11_an[1]),
                             (-_i11_at[0], -_i11_at[1]), to_svg)
    out.append(f'<polyline points="{arc_pts}" fill="none"'
               f' stroke="{JAMB_COLOR}" stroke-width="0.5"/>')

    # RO4 door: centered in opening, hinged at SE end midpoint, swings 90° to NW end
    ro4 = [r for r in rough_openings if r.name == "RO4"][0]
    # Hinge: midpoint of SE end (poly[0] and poly[1])
    hinge = ((ro4.poly[0][0] + ro4.poly[1][0]) / 2,
             (ro4.poly[0][1] + ro4.poly[1][1]) / 2)
    # Closed tip: midpoint of NW end (poly[2] and poly[3])
    _ro4_closed = ((ro4.poly[2][0] + ro4.poly[3][0]) / 2,
                   (ro4.poly[2][1] + ro4.poly[3][1]) / 2)
    # Door length and closed direction (hinge → NW end mid)
    _ro4_dx = _ro4_closed[0] - hinge[0]
    _ro4_dy = _ro4_closed[1] - hinge[1]
    _ro4_door_len = math.sqrt(_ro4_dx**2 + _ro4_dy**2)
    _ro4_len_u = (_ro4_dx / _ro4_door_len, _ro4_dy / _ro4_door_len)
    # Swing direction (perpendicular, WSW): NW→NE = poly[3]→poly[2]
    _ro4_sw_dx = ro4.poly[2][0] - ro4.poly[3][0]
    _ro4_sw_dy = ro4.poly[2][1] - ro4.poly[3][1]
    _ro4_sw_len = math.sqrt(_ro4_sw_dx**2 + _ro4_sw_dy**2)
    _ro4_swing = (_ro4_sw_dx / _ro4_sw_len, _ro4_sw_dy / _ro4_sw_len)
    # Door line (open position, perpendicular to opening)
    tip = (hinge[0] + _ro4_door_len * _ro4_swing[0],
           hinge[1] + _ro4_door_len * _ro4_swing[1])
    hx, hy = to_svg(*hinge)
    tx, ty = to_svg(*tip)
    out.append(f'<line x1="{hx:.1f}" y1="{hy:.1f}" x2="{tx:.1f}" y2="{ty:.1f}"'
               f' stroke="{JAMB_COLOR}" stroke-width="1.0"/>')
    # Arc from open (swing) sweeping 90° to closed (length direction)
    arc_pts = _swing_arc_svg(hinge, _ro4_door_len, _ro4_swing,
                             _ro4_len_u, to_svg)
    out.append(f'<polyline points="{arc_pts}" fill="none"'
               f' stroke="{JAMB_COLOR}" stroke-width="0.5"/>')

    # RO5 door: 36" door, hinged east, swings toward north face (rotation-safe)
    if not bare:
        ro5 = [r for r in rough_openings if r.name == "RO5"][0]
        # Along direction: poly[0]→poly[1] (SW→SE)
        _ro5_dx = ro5.poly[1][0] - ro5.poly[0][0]
        _ro5_dy = ro5.poly[1][1] - ro5.poly[0][1]
        _ro5_len = math.sqrt(_ro5_dx**2 + _ro5_dy**2)
        _ro5_al = (_ro5_dx / _ro5_len, _ro5_dy / _ro5_len)
        # End edge midpoint (east): midpoint of poly[1] and poly[2]
        _ro5_end = ((ro5.poly[1][0] + ro5.poly[2][0]) / 2,
                    (ro5.poly[1][1] + ro5.poly[2][1]) / 2)
        # Swing direction (toward north face): SE→NE = poly[1]→poly[2]
        _ro5_sw_dx = ro5.poly[2][0] - ro5.poly[1][0]
        _ro5_sw_dy = ro5.poly[2][1] - ro5.poly[1][1]
        _ro5_sw_len = math.sqrt(_ro5_sw_dx**2 + _ro5_sw_dy**2)
        _ro5_swing = (_ro5_sw_dx / _ro5_sw_len, _ro5_sw_dy / _ro5_sw_len)
        ro5_gap = (_ro5_len - RO5_DOOR_WIDTH) / 2
        hinge = (_ro5_end[0] - ro5_gap * _ro5_al[0],
                 _ro5_end[1] - ro5_gap * _ro5_al[1])
        hx, hy = to_svg(*hinge)
        tip = (hinge[0] + RO5_DOOR_WIDTH * _ro5_swing[0],
               hinge[1] + RO5_DOOR_WIDTH * _ro5_swing[1])
        tx, ty = to_svg(*tip)
        out.append(f'<line x1="{hx:.1f}" y1="{hy:.1f}" x2="{tx:.1f}" y2="{ty:.1f}"'
                   f' stroke="{JAMB_COLOR}" stroke-width="1.0"/>')
        # Arc from open (swing) sweeping 90° to closed (-along direction)
        arc_pts = _swing_arc_svg(hinge, RO5_DOOR_WIDTH, _ro5_swing,
                                 (-_ro5_al[0], -_ro5_al[1]), to_svg)
        out.append(f'<polyline points="{arc_pts}" fill="none"'
                   f' stroke="{JAMB_COLOR}" stroke-width="0.5"/>')

    # IW9 unit vectors (shared by RO3 and RO7 doors)
    _i9_sw, _i9_se, _i9_ne, _i9_nw = layout.iw9.poly
    _i9_dx_n = _i9_ne[0] - _i9_se[0]; _i9_dy_n = _i9_ne[1] - _i9_se[1]
    _i9_ln = math.sqrt(_i9_dx_n**2 + _i9_dy_n**2)
    _i9_an = (_i9_dx_n / _i9_ln, _i9_dy_n / _i9_ln)  # NNE (along length)
    _i9_dx_t = _i9_sw[0] - _i9_se[0]; _i9_dy_t = _i9_sw[1] - _i9_se[1]
    _i9_lt = math.sqrt(_i9_dx_t**2 + _i9_dy_t**2)
    _i9_at = (_i9_dx_t / _i9_lt, _i9_dy_t / _i9_lt)  # SE→SW (across thickness)

    # RO3 door: 36" door in IW9 (rotated), hinged at south edge, swings west
    ro3 = [r for r in rough_openings if r.name == "RO3"][0]
    _ro3p_d = ro3.poly  # [SW, SE, NE, NW] in IW9 orientation
    _ro3_gap = (math.sqrt((_ro3p_d[3][0] - _ro3p_d[0][0])**2 +
                           (_ro3p_d[3][1] - _ro3p_d[0][1])**2) - RO3_DOOR_WIDTH) / 2
    # South edge center: midpoint of SW (poly[0]) and SE (poly[1])
    _ro3_s_ctr = ((_ro3p_d[0][0] + _ro3p_d[1][0]) / 2,
                  (_ro3p_d[0][1] + _ro3p_d[1][1]) / 2)
    h_ro3 = (_ro3_s_ctr[0] + _ro3_gap * _i9_an[0],
             _ro3_s_ctr[1] + _ro3_gap * _i9_an[1])
    hx, hy = to_svg(*h_ro3)
    # Open position: west (_i9_at direction)
    tip_ro3 = (h_ro3[0] + RO3_DOOR_WIDTH * _i9_at[0],
               h_ro3[1] + RO3_DOOR_WIDTH * _i9_at[1])
    tx, ty = to_svg(*tip_ro3)
    out.append(f'<line x1="{hx:.1f}" y1="{hy:.1f}" x2="{tx:.1f}" y2="{ty:.1f}"'
               f' stroke="{JAMB_COLOR}" stroke-width="1.0"/>')
    # Arc from west (open) sweeping CW to north (closed toward center)
    arc_pts = _swing_arc_svg(h_ro3, RO3_DOOR_WIDTH, _i9_at, _i9_an, to_svg)
    out.append(f'<polyline points="{arc_pts}" fill="none"'
               f' stroke="{JAMB_COLOR}" stroke-width="0.5"/>')

    # RO7 double door: 2×24" doors in IW9 (rotated), hinged at outer edges, open east
    ro7 = [r for r in rough_openings if r.name == "RO7"][0]
    _ro7p = ro7.poly  # [SW, SE, NE, NW]
    _ro7_gap = (IW9_RO_WIDTH - 2 * RO7_DOOR_WIDTH) / 2
    # South door: hinged at south edge (east face), swings east
    _ro7_s_ctr = ((_ro7p[0][0] + _ro7p[1][0]) / 2, (_ro7p[0][1] + _ro7p[1][1]) / 2)
    h_s = (_ro7_s_ctr[0] + _ro7_gap * _i9_an[0], _ro7_s_ctr[1] + _ro7_gap * _i9_an[1])
    hsx, hsy = to_svg(*h_s)
    # Open position: east
    tip_s = (h_s[0] - RO7_DOOR_WIDTH * _i9_at[0], h_s[1] - RO7_DOOR_WIDTH * _i9_at[1])
    tsx, tsy = to_svg(*tip_s)
    out.append(f'<line x1="{hsx:.1f}" y1="{hsy:.1f}" x2="{tsx:.1f}" y2="{tsy:.1f}"'
               f' stroke="{JAMB_COLOR}" stroke-width="1.0"/>')
    # Arc from east (open) sweeping CCW to NNE (closed toward center)
    arc_pts = _swing_arc_svg(h_s, RO7_DOOR_WIDTH,
                             (-_i9_at[0], -_i9_at[1]), _i9_an, to_svg)
    out.append(f'<polyline points="{arc_pts}" fill="none"'
               f' stroke="{JAMB_COLOR}" stroke-width="0.5"/>')
    # North door: hinged at north edge (east face), swings east
    _ro7_n_ctr = ((_ro7p[3][0] + _ro7p[2][0]) / 2, (_ro7p[3][1] + _ro7p[2][1]) / 2)
    h_n = (_ro7_n_ctr[0] - _ro7_gap * _i9_an[0], _ro7_n_ctr[1] - _ro7_gap * _i9_an[1])
    hnx, hny = to_svg(*h_n)
    # Open position: east
    tip_n = (h_n[0] - RO7_DOOR_WIDTH * _i9_at[0], h_n[1] - RO7_DOOR_WIDTH * _i9_at[1])
    tnx, tny = to_svg(*tip_n)
    out.append(f'<line x1="{hnx:.1f}" y1="{hny:.1f}" x2="{tnx:.1f}" y2="{tny:.1f}"'
               f' stroke="{JAMB_COLOR}" stroke-width="1.0"/>')
    # Arc from east (open) sweeping CW to SSW (closed toward center)
    arc_pts = _swing_arc_svg(h_n, RO7_DOOR_WIDTH,
                             (-_i9_at[0], -_i9_at[1]),
                             (-_i9_an[0], -_i9_an[1]), to_svg)
    out.append(f'<polyline points="{arc_pts}" fill="none"'
               f' stroke="{JAMB_COLOR}" stroke-width="0.5"/>')

    # RO6 double door: 2×24" doors in IW11 (rotated), hinged at outer edges, open west
    ro6 = [r for r in rough_openings if r.name == "RO6"][0]
    _ro6p = ro6.poly  # [SW, SE, NE, NW]
    # IW11 unit vectors already computed (_i11_an, _i11_at)
    _ro6_gap = (IW11_RO_WIDTH - 2 * RO6_DOOR_WIDTH) / 2
    # South door: hinged at south edge (west face), swings west
    _ro6_s_ctr = ((_ro6p[0][0] + _ro6p[1][0]) / 2, (_ro6p[0][1] + _ro6p[1][1]) / 2)
    h_s6 = (_ro6_s_ctr[0] + _ro6_gap * _i11_an[0], _ro6_s_ctr[1] + _ro6_gap * _i11_an[1])
    hs6x, hs6y = to_svg(*h_s6)
    # Open position: west
    tip_s6 = (h_s6[0] + RO6_DOOR_WIDTH * _i11_at[0], h_s6[1] + RO6_DOOR_WIDTH * _i11_at[1])
    ts6x, ts6y = to_svg(*tip_s6)
    out.append(f'<line x1="{hs6x:.1f}" y1="{hs6y:.1f}" x2="{ts6x:.1f}" y2="{ts6y:.1f}"'
               f' stroke="{JAMB_COLOR}" stroke-width="1.0"/>')
    # Arc from west (open) sweeping CW to NNE (closed toward center)
    arc_pts = _swing_arc_svg(h_s6, RO6_DOOR_WIDTH, _i11_at, _i11_an, to_svg)
    out.append(f'<polyline points="{arc_pts}" fill="none"'
               f' stroke="{JAMB_COLOR}" stroke-width="0.5"/>')
    # North door: hinged at north edge (west face), swings west
    _ro6_n_ctr = ((_ro6p[3][0] + _ro6p[2][0]) / 2, (_ro6p[3][1] + _ro6p[2][1]) / 2)
    h_n6 = (_ro6_n_ctr[0] - _ro6_gap * _i11_an[0], _ro6_n_ctr[1] - _ro6_gap * _i11_an[1])
    hn6x, hn6y = to_svg(*h_n6)
    # Open position: west
    tip_n6 = (h_n6[0] + RO6_DOOR_WIDTH * _i11_at[0], h_n6[1] + RO6_DOOR_WIDTH * _i11_at[1])
    tn6x, tn6y = to_svg(*tip_n6)
    out.append(f'<line x1="{hn6x:.1f}" y1="{hn6y:.1f}" x2="{tn6x:.1f}" y2="{tn6y:.1f}"'
               f' stroke="{JAMB_COLOR}" stroke-width="1.0"/>')
    # Arc from west (open) sweeping CCW to SSW (closed toward center)
    arc_pts = _swing_arc_svg(h_n6, RO6_DOOR_WIDTH, _i11_at,
                             (-_i11_an[0], -_i11_an[1]), to_svg)
    out.append(f'<polyline points="{arc_pts}" fill="none"'
               f' stroke="{JAMB_COLOR}" stroke-width="0.5"/>')

    # Casement windows: O8, O9, O10 (23" openings, 45° swing, hinged at S-series face)
    for oname, hinge_idx, close_idx in [("O8", 0, 1), ("O9", 1, 0), ("O10", 0, 1)]:
        o = [o for o in outer_openings if o.name == oname][0]
        # Inward direction (outer face midpoint toward inner face midpoint)
        _omid = ((o.poly[0][0] + o.poly[1][0]) / 2, (o.poly[0][1] + o.poly[1][1]) / 2)
        _imid = ((o.poly[2][0] + o.poly[3][0]) / 2, (o.poly[2][1] + o.poly[3][1]) / 2)
        _iE = _imid[0] - _omid[0]; _iN = _imid[1] - _omid[1]
        _ilen = math.sqrt(_iE**2 + _iN**2)
        _idir = (_iE / _ilen, _iN / _ilen)
        # Offset hinge/close from outer face inward by SHELL_THICKNESS to S-series
        _hinge = (o.poly[hinge_idx][0] + SHELL_THICKNESS * _idir[0],
                  o.poly[hinge_idx][1] + SHELL_THICKNESS * _idir[1])
        _close = (o.poly[close_idx][0] + SHELL_THICKNESS * _idir[0],
                  o.poly[close_idx][1] + SHELL_THICKNESS * _idir[1])
        _dE = _close[0] - _hinge[0]
        _dN = _close[1] - _hinge[1]
        _wlen = math.sqrt(_dE**2 + _dN**2)
        _cdir = (_dE / _wlen, _dN / _wlen)
        # Rotation sign: swing outward (away from interior)
        _cross = _cdir[0] * _idir[1] - _cdir[1] * _idir[0]
        _rsign = -1 if _cross > 0 else 1
        _oa = _rsign * math.pi / 4  # open angle
        # Open tip
        _cos_oa = math.cos(_oa); _sin_oa = math.sin(_oa)
        _odir = (_cdir[0] * _cos_oa - _cdir[1] * _sin_oa,
                 _cdir[0] * _sin_oa + _cdir[1] * _cos_oa)
        _tip = (_hinge[0] + _wlen * _odir[0], _hinge[1] + _wlen * _odir[1])
        hx, hy = to_svg(*_hinge)
        tx, ty = to_svg(*_tip)
        _line_el = (f'<line x1="{hx:.1f}" y1="{hy:.1f}" x2="{tx:.1f}" y2="{ty:.1f}"'
                    f' stroke="{OPENING_STROKE}" stroke-width="1.0"/>')
        # Arc from open to closed
        n_arc = 10
        _arc_pts = []
        for i in range(n_arc + 1):
            _a = _oa * (1 - i / n_arc)
            _ca = math.cos(_a); _sa = math.sin(_a)
            _d = (_cdir[0] * _ca - _cdir[1] * _sa,
                  _cdir[0] * _sa + _cdir[1] * _ca)
            _pt = (_hinge[0] + _wlen * _d[0], _hinge[1] + _wlen * _d[1])
            sx, sy = to_svg(*_pt)
            _arc_pts.append(f"{sx:.1f},{sy:.1f}")
        _arc_el = (f'<polyline points="{" ".join(_arc_pts)}" fill="none"'
                   f' stroke="{OPENING_STROKE}" stroke-width="0.5"/>')
        out.append(f'<a href="{_CASEMENT_URL}">{_line_el}{_arc_el}</a>')


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
# Room area computation
# ============================================================

def compute_room_areas(data, layout):
    """Compute floor area for each room.  Returns dict {name: area_sf}.

    The polygons trace along IW faces and W-surface boundaries so that
    interior-wall material is excluded from every room.  The sum of all
    returned areas should closely match  inner_area - IW_area  (the
    title-block "interior sq ft" value).
    """
    pts = data.pts

    # --- helpers reused by several rooms ---
    _ro_list = compute_rough_openings(pts, layout)
    _ro1_bd = [r for r in _ro_list if r.name == "RO1"][0].poly
    _ro1_w_nf = _ro1_bd[3]          # RO1 NW: IW1 north face at RO1 west end
    _w9w10_al, _ = seg_vecs(pts["W9"], pts["W10"])
    outer_openings = compute_outer_openings(pts, layout)
    o6 = [o for o in outer_openings if o.name == "O6"][0]
    o6_w = o6.poly[0]               # O6 SW: W-surface at O6 west edge

    # --- BEDROOM ---
    bedroom = poly_area([
        (layout.iw9.poly[2][0], layout.iw1.poly[0][1]),
        (layout.iw11.poly[3][0], layout.iw1.poly[0][1]),
        (layout.iw11.poly[3][0], pts["W1"][1]),
        (layout.iw9.poly[2][0], pts["W1"][1]),
    ])

    # --- UTIL (T-shape) ---
    _util_poly = [
        layout.iw8.poly[0],                               # IW8 SW
        layout.iw8.poly[1],                               # IW8 SE (= IW2 west at IW8 south)
        layout.iw1.poly[0],                               # IW1 SW (= IW2 west at IW1 south)
        layout.iw9.poly[3],                               # IW9 NW (on IW1 south)
        (layout.iw9.poly[0][0], layout.iw7.poly[2][1]),   # IW9 at IW7 north
        layout.iw7.poly[3],                               # IW7 NW (= IW3 NE)
        layout.iw3.poly[3],                               # IW3 NW
        (layout.iw3.poly[0][0], pts["W1"][1]),            # IW3 west at south wall
        pts["W1"],                                         # W1
    ]
    _util_poly.extend(segment_polyline(data.inner_segs[0], pts)[1:])
    util = poly_area(_util_poly)

    # --- KITCHEN ---
    _iw2s_e_al, _ = seg_vecs(layout.iw2s.poly[1], layout.iw2s.poly[2])
    _iw2s_e_at_w9 = line_isect(layout.iw2s.poly[1], _iw2s_e_al,
                                pts["W9"], _w9w10_al)
    kitchen = poly_area([
        o6_w,                       # NW: O6 SW on W-surface
        _iw2s_e_at_w9,             # NE: IW2s east face at W9 northing
        layout.iw2s.poly[1],       # IW2s SE
        layout.iw2o.poly[3],       # IW2o NW (kitchen-side, north end)
        layout.iw2o.poly[0],       # IW2o SW (kitchen-side, south end)
        layout.iw2.poly[2],        # IW2 NE
        layout.iw2.poly[1],        # IW2 SE (on IW1 north face)
        _ro1_w_nf,                 # SW: RO1 NW on IW1 north face
    ])

    # --- LIVING ---
    _iw1_ne = layout.iw1.poly[2]
    _living_poly = [o6_w]
    _living_poly.append(segment_polyline(data.inner_segs[6], pts)[-1])
    for _si in range(7, 13):
        _pl = segment_polyline(data.inner_segs[_si], pts)
        _living_poly.extend(_pl[1:])
    _living_poly.append(_iw1_ne)
    _living_poly.append(_ro1_w_nf)
    living = poly_area(_living_poly)

    # --- BATH (subtract IW6 which sits entirely inside) ---
    _seg2_pl = segment_polyline(data.inner_segs[2], pts)
    _seg3_pl = segment_polyline(data.inner_segs[3], pts)
    _bath_poly = [
        layout.iw8.poly[3],          # IW8 NW
        layout.iw8.poly[2],          # IW8 NE
        layout.iw2.poly[3],          # IW2 NW (bath side)
        layout.iw2o.poly[1],         # IW2o SE (bath side)
        layout.iw2o.poly[2],         # IW2o NE (bath side)
        layout.iw2s.poly[0],         # IW2s SW (bath side)
        layout.iw2s.poly[3],         # IW2s NW (on seg 3)
        _seg3_pl[0],                  # W6
    ]
    _bath_poly.extend(reversed(_seg2_pl[:-1]))
    bath = poly_area(_bath_poly) - poly_area(layout.iw6.poly)

    # --- OFFICE ---
    _office_poly = [layout.iw5.poly[0]]
    _office_poly.append((pts["W15"][0], layout.iw5.poly[0][1]))
    _office_poly.append(pts["W15"])
    for _si in [14, 15, 16]:
        _office_poly.extend(segment_polyline(data.inner_segs[_si], pts)[1:])
    _office_poly.append(layout.iw4.poly[1])
    _office_poly.append(layout.iw4.poly[2])
    _office_poly.append(layout.iw12.poly[2])
    _office_poly.append(layout.iw12.poly[3])
    office = poly_area(_office_poly)

    # --- E CLOSET ---
    e_closet = poly_area([
        (layout.iw11.poly[1][0], layout.iw12.poly[0][1]),
        (layout.iw4.poly[0][0], layout.iw12.poly[0][1]),
        layout.iw4.poly[0],
        (layout.iw11.poly[1][0], pts["W1"][1]),
    ])

    # --- W CLOSET ---
    w_closet = poly_area([
        (layout.iw3.poly[1][0], layout.iw7.poly[0][1]),
        (layout.iw9.poly[0][0], layout.iw7.poly[0][1]),
        (layout.iw9.poly[0][0], pts["W1"][1]),
        (layout.iw3.poly[1][0], pts["W1"][1]),
    ])

    # --- STORAGE ---
    storage = poly_area([
        (layout.iw11.poly[1][0], layout.iw5.poly[3][1]),
        (pts["W14"][0], layout.iw5.poly[3][1]),
        (pts["W14"][0], layout.iw1.poly[0][1]),
        (layout.iw11.poly[1][0], layout.iw1.poly[0][1]),
    ])

    # --- WH ---
    _seg3_pl = segment_polyline(data.inner_segs[3], pts)
    _seg4_pl = segment_polyline(data.inner_segs[4], pts)
    _seg5_pl = segment_polyline(data.inner_segs[5], pts)
    _wh_poly = [layout.iw2s.poly[2]]
    _wh_poly.append(_seg3_pl[-1])
    _wh_poly.extend(_seg4_pl[1:])
    _wh_poly.extend(_seg5_pl[1:])
    _wh_poly.append((layout.iw2s.poly[2][0], pts["W9"][1]))
    wh = poly_area(_wh_poly)

    return {
        "BEDROOM": bedroom, "UTIL": util, "KITCHEN": kitchen,
        "LIVING": living, "BATH": bath, "OFFICE": office,
        "E CLOSET": e_closet, "W CLOSET": w_closet,
        "STORAGE": storage, "WH": wh,
    }


# ============================================================
# SVG rendering — orchestrator
# ============================================================

def _render_sf_extras(out, data, layout):
    """Render SF-specific extras: BEDROOM/OFFICE labels and RO1–O6 dashed line."""
    pts = data.pts
    to_svg = data.to_svg

    # --- Room areas (single source of truth) ---
    _areas = compute_room_areas(data, layout)

    # --- Room labels (same positioning as _render_furniture) ---
    # SVG scale for sf-label gap computation
    _, _py0 = to_svg(0, 0)
    _, _py1 = to_svg(0, 1)
    _svg_per_ft = _py0 - _py1
    _half_gap = (3.0 / 12.0) * _svg_per_ft  # 50% of KITCHEN/LIVING gap

    # IW1 north face directions
    _iw1_n_al, _iw1_n_cw = seg_vecs(layout.iw1.poly[3], layout.iw1.poly[2])
    _iw1_n_out = (-_iw1_n_cw[0], -_iw1_n_cw[1])

    # BEDROOM: right edge at W end of RO1, top at N end of RO3
    _ro_list = compute_rough_openings(pts, layout)
    _ro1_bd = [r for r in _ro_list if r.name == "RO1"][0].poly
    _ro3_bd = [r for r in _ro_list if r.name == "RO3"][0].poly
    _ro1_w_nf = _ro1_bd[3]    # RO1 NW: IW1 north face at RO1 west end
    _ro3_n_mid = ((_ro3_bd[2][0] + _ro3_bd[3][0]) / 2,
                  (_ro3_bd[2][1] + _ro3_bd[3][1]) / 2)
    bdx, bdy = to_svg(_ro1_w_nf[0], _ro3_n_mid[1])
    out.append(f'<text x="{bdx:.1f}" y="{bdy:.1f}" text-anchor="end" dominant-baseline="hanging"'
               f' font-family="Arial" font-size="8" fill="#666">BEDROOM</text>')
    _bedroom_sf = _areas["BEDROOM"]
    _bd_sf_y = bdy + 8.0 + _half_gap
    out.append(f'<text x="{bdx:.1f}" y="{_bd_sf_y:.1f}" text-anchor="end" dominant-baseline="hanging"'
               f' font-family="Arial" font-size="8" fill="#666">{_bedroom_sf:.1f} sf</text>')

    # UTIL: same northing as BEDROOM, centered in easting between south toilet and sink
    _iw8_al, _iw8_in = seg_vecs(layout.iw8.poly[0], layout.iw8.poly[1])
    _iw8_s_ref = layout.iw8.poly[0]
    _dryer_cx = sum(p[0] for p in layout.dryer.poly) / 4
    _dryer_cy = sum(p[1] for p in layout.dryer.poly) / 4
    _d_dryer_al = ((_dryer_cx - _iw8_s_ref[0]) * _iw8_al[0] +
                   (_dryer_cy - _iw8_s_ref[1]) * _iw8_al[1])
    _toilet_s = offset_pt(_iw8_s_ref, _d_dryer_al, _iw8_al)
    _ctr_cx = sum(p[0] for p in layout.ctr.poly) / 4
    _ctr_cy = sum(p[1] for p in layout.ctr.poly) / 4
    _sink_mid = ((_dryer_cx + _ctr_cx) / 2, (_dryer_cy + _ctr_cy) / 2)
    _d_sink_al = ((_sink_mid[0] - _iw8_s_ref[0]) * _iw8_al[0] +
                  (_sink_mid[1] - _iw8_s_ref[1]) * _iw8_al[1])
    _sink_s = offset_pt(_iw8_s_ref, _d_sink_al, _iw8_al)
    _sk_s = offset_pt(_sink_s, SINK_RY, _iw8_in)
    _util_e = (_toilet_s[0] + _sk_s[0]) / 2
    _util_n = _ro3_n_mid[1]
    utx, uty = to_svg(_util_e, _util_n)
    out.append(f'<text x="{utx:.1f}" y="{uty:.1f}" text-anchor="middle" dominant-baseline="hanging"'
               f' font-family="Arial" font-size="8" fill="#666">UTIL</text>')
    _util_sf = _areas["UTIL"]
    _ut_sf_y = uty + 8.0 + _half_gap
    out.append(f'<text x="{utx:.1f}" y="{_ut_sf_y:.1f}" text-anchor="middle" dominant-baseline="hanging"'
               f' font-family="Arial" font-size="8" fill="#666">{_util_sf:.1f} sf</text>')

    # KITCHEN: centered beneath kitchen sink, just above dim02 (25' 8.1") line
    _w9w10_al, _ = seg_vecs(pts["W9"], pts["W10"])
    _iw2s_ne_k = layout.iw2s.poly[2]
    _iw2_d_k = ((_iw2s_ne_k[0] - pts["W9"][0]) * _w9w10_al[0] +
                (_iw2s_ne_k[1] - pts["W9"][1]) * _w9w10_al[1])
    _st_d_k = _iw2_d_k + NORTH_CTR_LENGTH + KITCHEN_APPL_GAP
    _ks_d_k = _st_d_k + STOVE_WIDTH + KITCHEN_APPL_GAP + 2.0 / 12.0
    _sink_ctr = offset_pt(pts["W9"], _ks_d_k + KITCHEN_SINK_WIDTH / 2, _w9w10_al)
    _dim02_n = (pts["F12"][1] + pts["F13"][1]) / 2
    kx, ky = to_svg(_sink_ctr[0], _dim02_n + 3.0 / 12.0)
    out.append(f'<text x="{kx:.1f}" y="{ky:.1f}" text-anchor="middle"'
               f' font-family="Arial" font-size="8" fill="#666">KITCHEN</text>')

    outer_openings = compute_outer_openings(pts, layout)
    o6 = [o for o in outer_openings if o.name == "O6"][0]
    o6_w = o6.poly[0]              # O6 SW: W-surface at O6 west edge
    _kitchen_sf = _areas["KITCHEN"]
    # SF label: centered under KITCHEN, equal distance below dim02 as KITCHEN is above
    sfx, sfy = to_svg(_sink_ctr[0], _dim02_n - 3.0 / 12.0)
    out.append(f'<text x="{sfx:.1f}" y="{sfy:.1f}" text-anchor="middle" dominant-baseline="hanging"'
               f' font-family="Arial" font-size="8" fill="#666">{_kitchen_sf:.1f} sf</text>')

    # LIVING: centered under O6 at KITCHEN label's northing
    _o6_cx = sum(p[0] for p in o6.poly) / 4
    _living_n = _dim02_n + 3.0 / 12.0
    lx, ly = to_svg(_o6_cx, _living_n)
    out.append(f'<text x="{lx:.1f}" y="{ly:.1f}" text-anchor="middle"'
               f' font-family="Arial" font-size="8" fill="#666">LIVING</text>')

    _living_sf = _areas["LIVING"]
    # SF label: centered under LIVING, same offset pattern
    lsfx, lsfy = to_svg(_o6_cx, _dim02_n - 3.0 / 12.0)
    out.append(f'<text x="{lsfx:.1f}" y="{lsfy:.1f}" text-anchor="middle" dominant-baseline="hanging"'
               f' font-family="Arial" font-size="8" fill="#666">{_living_sf:.1f} sf</text>')

    # BATH: centered between IW2s west and W2-W5, at RO4 north end northing
    _bath_e = (layout.iw2s.poly[0][0] + pts["W2"][0]) / 2
    _ro4_bd = [r for r in _ro_list if r.name == "RO4"][0].poly
    _ro4_n_mid = ((_ro4_bd[2][0] + _ro4_bd[3][0]) / 2,
                  (_ro4_bd[2][1] + _ro4_bd[3][1]) / 2)
    bax, bay = to_svg(_bath_e, _ro4_n_mid[1])
    out.append(f'<text x="{bax:.1f}" y="{bay:.1f}" text-anchor="middle"'
               f' font-family="Arial" font-size="8" fill="#666">BATH</text>')
    _bath_sf = _areas["BATH"]
    _ba_sf_y = bay + _half_gap
    out.append(f'<text x="{bax:.1f}" y="{_ba_sf_y:.1f}" text-anchor="middle" dominant-baseline="hanging"'
               f' font-family="Arial" font-size="8" fill="#666">{_bath_sf:.1f} sf</text>')

    # OFFICE: midpoint between IW4 east face and W15, vertically between ctr+5'+3" and IW1
    _iw4_e_mid = ((layout.iw4.poly[1][0] + layout.iw4.poly[2][0]) / 2,
                  (layout.iw4.poly[1][1] + layout.iw4.poly[2][1]) / 2)
    _of_ew = ((_iw4_e_mid[0] + pts["W15"][0]) / 2,
              (_iw4_e_mid[1] + pts["W15"][1]) / 2)
    _ctr_s_mid = ((layout.ctr.poly[0][0] + layout.ctr.poly[1][0]) / 2,
                  (layout.ctr.poly[0][1] + layout.ctr.poly[1][1]) / 2)
    _iw1_s_mid_r = ((layout.iw1.poly[0][0] + layout.iw1.poly[1][0]) / 2,
                    (layout.iw1.poly[0][1] + layout.iw1.poly[1][1]) / 2)
    _ctr_offset = offset_pt(_ctr_s_mid, 5.0 + WALL_3IN, _iw1_n_out)
    _of_ns = ((_ctr_offset[0] + _iw1_s_mid_r[0]) / 2,
              (_ctr_offset[1] + _iw1_s_mid_r[1]) / 2)
    _of_ns_adj = offset_pt(_of_ns, -2.0 + 26.0 / 12.0, _iw1_n_out)
    _of_ns_d = ((_of_ns_adj[0] - _of_ew[0]) * _iw1_n_out[0] +
                (_of_ns_adj[1] - _of_ew[1]) * _iw1_n_out[1])
    of_cx = _of_ew[0] + _of_ns_d * _iw1_n_out[0]
    of_cy = _of_ew[1] + _of_ns_d * _iw1_n_out[1]
    ofx, ofy = to_svg(of_cx, of_cy)
    out.append(f'<text x="{ofx:.1f}" y="{ofy+3:.1f}" text-anchor="middle" font-family="Arial"'
               f' font-size="8" fill="#666">OFFICE</text>')
    _office_sf = _areas["OFFICE"]
    _of_sf_y = ofy + 3 + _half_gap
    out.append(f'<text x="{ofx:.1f}" y="{_of_sf_y:.1f}" text-anchor="middle" dominant-baseline="hanging"'
               f' font-family="Arial" font-size="8" fill="#666">{_office_sf:.1f} sf</text>')

    # --- CLOSET and STORAGE sf labels (on opposite side of dim lines) ---
    _ep = {name: pt for name, pt in compute_dimension_endpoints(
        pts, layout, data.radii, bare=False)}

    _e_closet_sf = _areas["E CLOSET"]
    _w_closet_sf = _areas["W CLOSET"]
    _storage_sf = _areas["STORAGE"]

    # Place each sf label on the opposite side of its dim line from the dim label
    for _dim_key, _area, _rot in [
        ("dim03", _e_closet_sf, True),
        ("dim04", _w_closet_sf, True),
        ("dim07", _storage_sf, False),
    ]:
        _da, _db = _ep[f"{_dim_key}_A"], _ep[f"{_dim_key}_B"]
        _sx1, _sy1 = to_svg(*_da)
        _sx2, _sy2 = to_svg(*_db)
        _sdx = _sx2 - _sx1; _sdy = _sy2 - _sy1
        _slen = math.sqrt(_sdx**2 + _sdy**2)
        _ang = math.degrees(math.atan2(_sdy, _sdx))
        if _ang >= 90: _ang -= 180
        elif _ang < -90: _ang += 180
        _ang_rad = math.radians(_ang)
        _mx = (_sx1 + _sx2) / 2; _my = (_sy1 + _sy2) / 2
        # Opposite side: negate the label offset direction
        _sfx = _mx - 3 * math.sin(_ang_rad)
        _sfy = _my + 3 * math.cos(_ang_rad)
        _rot_attr = f' transform="rotate({_ang:.1f},{_sfx:.1f},{_sfy:.1f})"' if _rot else ""
        out.append(
            f'<text x="{_sfx:.1f}" y="{_sfy:.1f}" text-anchor="middle"'
            f' dominant-baseline="hanging" font-family="Arial" font-size="8"'
            f' fill="#666"{_rot_attr}>{_area:.1f} sf</text>')

    _wh_sf = _areas["WH"]
    _wh_e = (layout.iw2s.poly[2][0] + pts["W8"][0]) / 2
    _wh_n = (pts["W7"][1] + pts["W9"][1]) / 2
    _whx, _why = to_svg(_wh_e, _wh_n)
    # Vertically center the WH/sf pair: offset each by half the pair height
    _wh_pair_h = 6.0 + _half_gap / 2 + 6.0  # WH + gap + sf
    _why_wh = _why - _wh_pair_h / 2 + 6.0   # baseline of WH (bottom of text)
    _why_sf = _why_wh + _half_gap / 2         # top of sf text (hanging)
    out.append(f'<text x="{_whx:.1f}" y="{_why_wh:.1f}" text-anchor="middle"'
               f' font-family="Arial" font-size="6" fill="#666">WH</text>')
    out.append(f'<text x="{_whx:.1f}" y="{_why_sf:.1f}" text-anchor="middle"'
               f' dominant-baseline="hanging" font-family="Arial" font-size="6"'
               f' fill="#666">{_wh_sf:.1f} sf</text>')

    # --- Dashed line from IW1-north at RO1 west to W-surface at O6 west ---
    ro1_w = _ro1_w_nf  # RO1 NW on IW1 north face

    # o6_w = O6 SW on W-surface, already computed above
    r1x, r1y = to_svg(*ro1_w)
    o6x, o6y = to_svg(*o6_w)
    out.append(f'<line x1="{r1x:.1f}" y1="{r1y:.1f}" x2="{o6x:.1f}" y2="{o6y:.1f}"'
               f' stroke="#666" stroke-width="0.7" stroke-dasharray="4,3"/>')

    # --- Dashed line from W9 E-W to IW2s east face ---
    _iw2s_e_al, _ = seg_vecs(layout.iw2s.poly[1], layout.iw2s.poly[2])
    _iw2s_at_w9 = line_isect(layout.iw2s.poly[1], _iw2s_e_al,
                              pts["W9"], _w9w10_al)
    w9x, w9y = to_svg(*pts["W9"])
    i2x, i2y = to_svg(*_iw2s_at_w9)
    out.append(f'<line x1="{w9x:.1f}" y1="{w9y:.1f}" x2="{i2x:.1f}" y2="{i2y:.1f}"'
               f' stroke="#666" stroke-width="0.7" stroke-dasharray="4,3"/>')


def render_floorplan_svg(data, room_title="Parent Suite", minik=False, db=False, bare=False, sf=False):
    """Render the complete floorplan SVG. Returns SVG string.

    If bare=True, omit appliances, kitchen, and furniture (interior objects).
    If sf=True, render like bare but add BEDROOM/OFFICE labels and RO1–O6 dashed line.
    """
    pts = data.pts
    to_svg = data.to_svg
    layout = data.layout

    out = []
    out.append(f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}"'
               f' viewBox="{data.vb_x:.2f} {data.vb_y:.2f} {data.vb_w:.2f} {data.vb_h:.2f}">')
    out.append(f'<rect x="{data.vb_x:.2f}" y="{data.vb_y:.2f}" width="{data.vb_w:.2f}" height="{data.vb_h:.2f}" fill="white"/>')
    out.append('<defs>')
    out.append('  <marker id="ah" markerWidth="8" markerHeight="6" refX="8" refY="3" orient="auto">'
               '<polygon points="0 0, 8 3, 0 6" fill="#333"/></marker>')
    out.append('</defs>')
    out.append(f'<text x="{data.title_x:.1f}" y="{data.title_y:.1f}" text-anchor="middle" font-family="Arial" font-size="14"'
               f' font-weight="bold">{room_title}</text>')

    _render_walls(out, data, layout, bare=bare or sf)
    if not bare and not sf:
        _render_appliances(out, data, layout, minik=minik, db=db)
        _render_kitchen(out, data, layout, minik=minik, db=db)
        _render_furniture(out, data, layout, minik=minik, db=db)
    out.append('<g opacity="0.5">')
    _render_dimensions(out, data, layout, bare=bare or sf)
    out.append('</g>')
    _render_openings(out, data, layout, bare=bare or sf)
    if sf:
        _render_sf_extras(out, data, layout)

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

    base_dir = os.path.dirname(os.path.abspath(__file__))

    svg_path = os.path.join(base_dir, "floorplan.svg")
    with open(svg_path, "w") as f:
        f.write(svg_content)
    print(f"Floorplan written to {svg_path}")

    minik_content = render_floorplan_svg(data, room_title="Parent Suite w/Small Kitchen", minik=True)
    minik_path = os.path.join(base_dir, "floorplan_minik.svg")
    with open(minik_path, "w") as f:
        f.write(minik_content)
    print(f"Floorplan (minik) written to {minik_path}")

    db_content = render_floorplan_svg(data, room_title="Parent Suite with Daybed", db=True)
    db_path = os.path.join(base_dir, "floorplan_db.svg")
    with open(db_path, "w") as f:
        f.write(db_content)
    print(f"Floorplan (daybed) written to {db_path}")

    bare_content = render_floorplan_svg(data, room_title="Room Dimensions", bare=True)
    bare_path = os.path.join(base_dir, "floorplan_bare.svg")
    with open(bare_path, "w") as f:
        f.write(bare_content)
    print(f"Floorplan (bare) written to {bare_path}")

    sf_content = render_floorplan_svg(data, room_title="Room Dimensions", sf=True)
    sf_path = os.path.join(base_dir, "floorplan_sf.svg")
    with open(sf_path, "w") as f:
        f.write(sf_content)
    print(f"Floorplan (sf) written to {sf_path}")

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
