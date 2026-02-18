"""Generate floorplan SVG with 8" wall inset from the outline path.

Computes geometry from shared/ and floorplan/ packages.
Outline points F0-F20, inner wall points W0-W20.
"""
import os, math, datetime
from typing import NamedTuple, Any

from shared.types import LineSeg, ArcSeg, BBox
from shared.geometry import (
    segment_polyline, path_polygon, poly_area,
    compute_inner_walls, fmt_dist, f8f9_corner_polyline,
    horiz_isects,
)
from shared.survey import compute_traverse, compute_three_arc, compute_inset
from shared.svg import make_svg_transform, W, H, git_describe
from floorplan.geometry import compute_outline_geometry, OutlineAnchors
from floorplan.constants import (
    WALL_OUTER, SHELL_THICKNESS, AIR_GAP, OPENING_INSIDE_RADIUS,
    WH_RADIUS,
    SINK_RX, SINK_RY,
    KITCHEN_SINK_WIDTH, KITCHEN_SINK_DEPTH,
    DW_WIDTH, DW_DEPTH, STOVE_WIDTH, STOVE_DEPTH,
    FRIDGE_SIZE, MINIK_FRIDGE_W, MINIK_FRIDGE_D,
    KITCHEN_CTR_LENGTH, KITCHEN_CTR_DEPTH,
    NORTH_CTR_LENGTH, NORTH_CTR_DEPTH,
    JAMB_WIDTH, STD_GAP, KITCHEN_APPL_GAP,
    LOVESEAT_NW_E, LOVESEAT_NW_N,
    LOVESEAT_WIDTH, LOVESEAT_LENGTH, LOVESEAT_ANGLE_DEG,
    DESK_WIDTH, DESK_DEPTH, DESK_CHAIR_WIDTH, DESK_CHAIR_DEPTH, DESK_CHAIR_GAP,
    CHAIR_WIDTH, CHAIR_DEPTH, CHAIR_CORNER_R, CHAIR_ANGLE_DEG,
    OTTOMAN_SIZE, ET_RADIUS_CM,
    SOFA_WIDTH, SOFA_DEPTH,
    ICE_WIDTH, ICE_DEPTH,
    ROCKER_WIDTH, ROCKER_DEPTH, ROCKER_CORNER_R,
    RO1_OFFSET_E_IW2, IW1_RO_WIDTH,
    O3_HALF_WIDTH, O3_DOOR_WIDTH,
    O6_WIDTH, O6_DOOR_WIDTH, RO1_DOOR_WIDTH, RO2_DOOR_WIDTH,
    RO3_DOOR_WIDTH, RO4_DOOR_WIDTH, RO5_DOOR_WIDTH, IW4_RO_WIDTH, DOOR_FLAT_FACE, F8F9_INNER_TURN_R,
)
from floorplan.layout import compute_interior_layout
from floorplan.openings import (
    compute_outer_openings, compute_rough_openings, outer_to_wall_openings,
)
from shared.wall_shells import (
    compute_inset_path, lerp, openings_on_seg, solid_ranges,
    arc_strip_poly, line_strip_poly, partial_line_strip, partial_line_strip_2,
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

def dim_line_v(out, e, n1, n2, label, to_svg, label_n=None):
    """Vertical (N-S) dimension line with horizontal tick marks."""
    x1, y1 = to_svg(e, n1); x2, y2 = to_svg(e, n2)
    _t = 4
    out.append(f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" stroke="{DIM_COLOR}" stroke-width="0.8"/>')
    out.append(f'<line x1="{x1-_t:.1f}" y1="{y1:.1f}" x2="{x1+_t:.1f}" y2="{y1:.1f}" stroke="{DIM_COLOR}" stroke-width="0.8"/>')
    out.append(f'<line x1="{x2-_t:.1f}" y1="{y2:.1f}" x2="{x2+_t:.1f}" y2="{y2:.1f}" stroke="{DIM_COLOR}" stroke-width="0.8"/>')
    if label_n is not None:
        _, ly_base = to_svg(e, label_n)
        lx, ly = x1 - 3, ly_base + 3
    else:
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

    # Replace W8-W9 arc in inner_poly with straight-arc-straight path
    w_f8f9_poly = f8f9_corner_polyline(pts, WALL_OUTER, F8F9_INNER_TURN_R)
    w8 = pts["W8"]
    w9 = pts["W9"]
    w8_idx = next(i for i, p in enumerate(inner_poly)
                  if abs(p[0] - w8[0]) < 1e-9 and abs(p[1] - w8[1]) < 1e-9)
    w9_idx = next(i for i, p in enumerate(inner_poly)
                  if i > w8_idx
                  and abs(p[0] - w9[0]) < 1e-9 and abs(p[1] - w9[1]) < 1e-9)
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
    _f_names = [f"F{i}" for i in range(21)]
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
        w_f8f9_poly=w_f8f9_poly,
        s_segs=s_segs, g_segs=g_segs,
        openings=openings,
        g_f8f9_poly=g_f8f9_poly,
        layout=layout,
    )

# ============================================================
# Render sub-functions
# ============================================================

def compute_iw_area(layout):
    """Compute total interior wall area from layout polygons."""
    iw2 = layout.iw2
    iw5 = layout.iw5
    iw2_poly = [(iw2.w, iw2.s), (iw2.e, iw2.s), (iw2.e, iw2.n), (iw2.w, iw2.n)]
    iw3_poly = layout.iw3_poly
    iw7_poly = layout.iw7_poly
    iw9_poly = layout.iw9_poly
    _iw4_n_area = layout.iw12_poly[2][1]  # IW12 NE northing
    iw4_poly = [(layout.iw4_w, layout.iw4_s), (layout.iw4_e, layout.iw4_s),
                (layout.iw4_e, _iw4_n_area), (layout.iw4_w, _iw4_n_area)]
    iw5_poly = [(iw5.w, iw5.s), (iw5.e, iw5.s), (iw5.e, iw5.n), (iw5.w, iw5.n)]
    iw8 = layout.iw8
    iw8_poly = [(iw8.w, iw8.s), (iw8.e, iw8.s), (iw8.e, iw8.n), (iw8.w, iw8.n)]
    iw15 = layout.iw15
    iw15_poly = [(iw15.w, iw15.s), (iw15.e, iw15.s), (iw15.e, iw15.n), (iw15.w, iw15.n)]
    iw_polys = [layout.iw1, iw8_poly, iw2_poly, iw3_poly, iw7_poly, iw9_poly, layout.iw6_poly,
                iw4_poly, layout.iw11_poly, layout.iw12_poly,
                layout.iw14_poly, iw5_poly, layout.iw16_poly, iw15_poly]
    return sum(poly_area(p) for p in iw_polys)


def _svg_wall_poly(out, poly, to_svg):
    """Render a shell polygon with floorplan wall styling (no stroke, gray fill)."""
    svg = " ".join(f"{to_svg(*p)[0]:.1f},{to_svg(*p)[1]:.1f}" for p in poly)
    out.append(f'<polygon points="{svg}" fill="{WALL_FILL}" stroke="none"/>')


def _render_walls(out, data, layout):
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

            if seg_idx == 8:
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
                        pts, seg, s_seg.start, s_seg.end, t_s, t_e)
                    _svg_wall_poly(out, outer_strip, to_svg)

                    inner_strip = partial_line_strip_2(
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
                    out.append(f'<polygon points="{svg}" fill="{OPENING_FILL}" '
                               f'stroke="{OPENING_STROKE}" stroke-width="{WALL_SW}"/>')

    # --- Continuous section outlines ---
    g_overrides = {8: data.g_f8f9_poly}
    w_overrides = {8: data.w_f8f9_poly}
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

    # ---- IW8 (no openings) ----
    iw8 = layout.iw8
    iw8_poly = [(iw8.w, iw8.s), (iw8.e, iw8.s), (iw8.e, iw8.n), (iw8.w, iw8.n)]
    wall_poly(out, iw8_poly, to_svg, stroke=False)
    for n_val in [iw8.s + half_sw, iw8.n - half_sw]:
        sx1, sy1 = to_svg(iw8.w, n_val)
        sx2, sy2 = to_svg(iw8.e, n_val)
        out.append(f'<line x1="{sx1:.1f}" y1="{sy1:.1f}" x2="{sx2:.1f}" y2="{sy2:.1f}"'
                   f' stroke="{WALL_STROKE}" stroke-width="{WALL_SW}"/>')

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

    # ---- IW3 (solid, no opening, rotated perpendicular to W20-W0) ----
    _iw3_sw, _iw3_se, _iw3_ne, _iw3_nw = layout.iw3_poly
    wall_poly(out, layout.iw3_poly, to_svg, stroke=False)
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

    # ---- IW7 (solid, no opening, rotated parallel to W20-W0) ----
    _iw7_sw, _iw7_se, _iw7_ne, _iw7_nw = layout.iw7_poly
    wall_poly(out, layout.iw7_poly, to_svg, stroke=False)
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

    # ---- IW9 (solid, no opening, rotated perpendicular to W20-W0) ----
    _iw9_sw, _iw9_se, _iw9_ne, _iw9_nw = layout.iw9_poly
    wall_poly(out, layout.iw9_poly, to_svg, stroke=False)
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
    for (p1, p2), (ox, oy) in [
        ((_iw9_se, _iw9_ne), _iw9_at),                            # east face
        ((_iw9_sw, _iw9_nw), (-_iw9_at[0], -_iw9_at[1])),        # west face
    ]:
        sx1, sy1 = to_svg(p1[0] + half_sw * ox, p1[1] + half_sw * oy)
        sx2, sy2 = to_svg(p2[0] + half_sw * ox, p2[1] + half_sw * oy)
        out.append(f'<line x1="{sx1:.1f}" y1="{sy1:.1f}" x2="{sx2:.1f}" y2="{sy2:.1f}"'
                   f' stroke="{WALL_STROKE}" stroke-width="{WALL_SW}"/>')

    # ---- IW16 with RO3 ----
    _iw16 = layout.iw16_poly
    ro3_n, ro3_s = ro["RO3"].n, ro["RO3"].s
    _iw16_w = _iw16[0][0]
    _iw16_e = _iw16[1][0]
    _iw16_s = _iw16[0][1]
    _iw16_n = _iw16[2][1]

    iw16_s_poly = [(_iw16_w, _iw16_s), (_iw16_e, _iw16_s), (_iw16_e, ro3_s), (_iw16_w, ro3_s)]
    iw16_n_poly = [(_iw16_w, ro3_n), (_iw16_e, ro3_n), (_iw16_e, _iw16_n), (_iw16_w, _iw16_n)]
    wall_poly(out, iw16_s_poly, to_svg, stroke=False)
    wall_poly(out, iw16_n_poly, to_svg, stroke=False)
    _iw16_w_in = _iw16_w + half_sw
    _iw16_e_in = _iw16_e - half_sw
    for a, b in [((_iw16_w_in, _iw16_s), (_iw16_w_in, ro3_s)),
                 ((_iw16_w_in, ro3_n), (_iw16_w_in, _iw16_n)),
                 ((_iw16_e_in, _iw16_s), (_iw16_e_in, ro3_s)),
                 ((_iw16_e_in, ro3_n), (_iw16_e_in, _iw16_n))]:
        sx1, sy1 = to_svg(*a); sx2, sy2 = to_svg(*b)
        out.append(f'<line x1="{sx1:.1f}" y1="{sy1:.1f}" x2="{sx2:.1f}" y2="{sy2:.1f}"'
                   f' stroke="{WALL_STROKE}" stroke-width="{WALL_SW}"/>')
    for jamb_n in [ro3_s, ro3_n - JAMB_WIDTH]:
        jx1, jy1 = to_svg(_iw16_w, jamb_n + JAMB_WIDTH)
        jx2, jy2 = to_svg(_iw16_e, jamb_n)
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

    # ---- IW4 (solid, no opening) — north end at IW12 north face ----
    iw4_n = layout.iw12_poly[2][1]  # IW12 NE northing
    iw4_poly = [(layout.iw4_w, layout.iw4_s), (layout.iw4_e, layout.iw4_s),
                (layout.iw4_e, iw4_n), (layout.iw4_w, iw4_n)]
    wall_poly(out, iw4_poly, to_svg, stroke=False)
    w_in = layout.iw4_w + half_sw
    e_in = layout.iw4_e - half_sw
    for a, b in [((w_in, layout.iw4_s), (w_in, iw4_n)),
                 ((e_in, layout.iw4_s), (e_in, iw4_n))]:
        sx1, sy1 = to_svg(*a); sx2, sy2 = to_svg(*b)
        out.append(f'<line x1="{sx1:.1f}" y1="{sy1:.1f}" x2="{sx2:.1f}" y2="{sy2:.1f}"'
                   f' stroke="{WALL_STROKE}" stroke-width="{WALL_SW}"/>')

    # ---- IW11 with RO2 (rotated rectangle, split by opening) ----
    _iw11_sw, _iw11_se, _iw11_ne, _iw11_nw = layout.iw11_poly
    # Unit vectors from polygon corners
    _iw11_dx_t = _iw11_sw[0] - _iw11_se[0]  # thickness direction (along)
    _iw11_dy_t = _iw11_sw[1] - _iw11_se[1]
    _iw11_lt = math.sqrt(_iw11_dx_t**2 + _iw11_dy_t**2)
    _iw11_at = (_iw11_dx_t / _iw11_lt, _iw11_dy_t / _iw11_lt)  # unit along
    _iw11_dx_n = _iw11_ne[0] - _iw11_se[0]  # length direction (normal)
    _iw11_dy_n = _iw11_ne[1] - _iw11_se[1]
    _iw11_ln = math.sqrt(_iw11_dx_n**2 + _iw11_dy_n**2)
    _iw11_an = (_iw11_dx_n / _iw11_ln, _iw11_dy_n / _iw11_ln)  # unit normal

    # RO2 polygon corners on IW11
    ro2_ro = [r for r in rough_openings if r.name == "RO2"][0]
    _ro2p = ro2_ro.poly  # [SW, SE, NE, NW]

    # South segment: IW11 south end to RO2 south edge
    iw11_s_poly = [_iw11_sw, _iw11_se, _ro2p[0], _ro2p[1]]
    wall_poly(out, iw11_s_poly, to_svg, stroke=False)
    # North segment: RO2 north edge to IW11 north end
    iw11_n_poly = [_ro2p[2], _ro2p[3], _iw11_ne, _iw11_nw]
    wall_poly(out, iw11_n_poly, to_svg, stroke=False)

    # Inset stroke lines for south segment
    for (p1, p2), (ox, oy) in [
        ((_iw11_se, _ro2p[0]), _iw11_at),     # east face south
        ((_iw11_sw, _ro2p[1]), (-_iw11_at[0], -_iw11_at[1])),  # west face south
        ((_iw11_sw, _iw11_se), _iw11_an),      # south end
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
    # Jambs at RO2 edges
    for (j1, j2) in [(_ro2p[0], _ro2p[1]), (_ro2p[3], _ro2p[2])]:
        # Jamb is JAMB_WIDTH wide along the opening direction
        _jn = (_iw11_an[0] * JAMB_WIDTH, _iw11_an[1] * JAMB_WIDTH)
        j_poly = [j1, j2, (j2[0] + _jn[0], j2[1] + _jn[1]),
                  (j1[0] + _jn[0], j1[1] + _jn[1])]
        jp = " ".join(f"{to_svg(*p)[0]:.1f},{to_svg(*p)[1]:.1f}" for p in j_poly)
        out.append(f'<polygon points="{jp}" fill="{JAMB_COLOR}" stroke="none"/>')

    # ---- IW12 (rotated rectangle) ----
    wall_poly(out, layout.iw12_poly, to_svg, stroke=False)
    _iw12_sw, _iw12_se, _iw12_ne, _iw12_nw = layout.iw12_poly
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

    # ---- IW14 (rotated rectangle, 3" thick, parallel to IW12) ----
    wall_poly(out, layout.iw14_poly, to_svg, stroke=False)
    _iw14_sw, _iw14_se, _iw14_ne, _iw14_nw = layout.iw14_poly
    _iw14_dx_t = _iw14_se[0] - _iw14_sw[0]  # length direction (-along)
    _iw14_dy_t = _iw14_se[1] - _iw14_sw[1]
    _iw14_lt = math.sqrt(_iw14_dx_t**2 + _iw14_dy_t**2)
    _iw14_al = (_iw14_dx_t / _iw14_lt, _iw14_dy_t / _iw14_lt)  # unit along length
    _iw14_dx_n = _iw14_nw[0] - _iw14_sw[0]  # thickness direction (norm)
    _iw14_dy_n = _iw14_nw[1] - _iw14_sw[1]
    _iw14_ln = math.sqrt(_iw14_dx_n**2 + _iw14_dy_n**2)
    _iw14_an = (_iw14_dx_n / _iw14_ln, _iw14_dy_n / _iw14_ln)  # unit normal
    for (p1, p2), (ox, oy) in [
        ((_iw14_sw, _iw14_se), _iw14_an),      # south face, inset toward north
        ((_iw14_nw, _iw14_ne), (-_iw14_an[0], -_iw14_an[1])),  # north face
        ((_iw14_sw, _iw14_nw), _iw14_al),      # west end, inset toward east
        ((_iw14_se, _iw14_ne), (-_iw14_al[0], -_iw14_al[1])),  # east end
    ]:
        sx1, sy1 = to_svg(p1[0] + half_sw * ox, p1[1] + half_sw * oy)
        sx2, sy2 = to_svg(p2[0] + half_sw * ox, p2[1] + half_sw * oy)
        out.append(f'<line x1="{sx1:.1f}" y1="{sy1:.1f}" x2="{sx2:.1f}" y2="{sy2:.1f}"'
                   f' stroke="{WALL_STROKE}" stroke-width="{WALL_SW}"/>')

    # ---- IW15 (no openings, N-S from IW11 north to IW1 south) ----
    iw15 = layout.iw15
    iw15_poly = [(iw15.w, iw15.s), (iw15.e, iw15.s), (iw15.e, iw15.n), (iw15.w, iw15.n)]
    wall_poly(out, iw15_poly, to_svg, stroke=False)
    for e_val in [iw15.w + half_sw, iw15.e - half_sw]:
        sx1, sy1 = to_svg(e_val, iw15.s)
        sx2, sy2 = to_svg(e_val, iw15.n)
        out.append(f'<line x1="{sx1:.1f}" y1="{sy1:.1f}" x2="{sx2:.1f}" y2="{sy2:.1f}"'
                   f' stroke="{WALL_STROKE}" stroke-width="{WALL_SW}"/>')

    # ---- IW5 ----
    iw5 = layout.iw5
    iw5_poly = [(iw5.w, iw5.s), (iw5.e, iw5.s), (iw5.e, iw5.n), (iw5.w, iw5.n)]
    wall_poly(out, iw5_poly, to_svg, stroke=False)
    for n_val in [iw5.s + half_sw, iw5.n - half_sw]:
        sx1, sy1 = to_svg(iw5.w, n_val)
        sx2, sy2 = to_svg(iw5.e, n_val)
        out.append(f'<line x1="{sx1:.1f}" y1="{sy1:.1f}" x2="{sx2:.1f}" y2="{sy2:.1f}"'
                   f' stroke="{WALL_STROKE}" stroke-width="{WALL_SW}"/>')



def _render_appliances(out, data, layout, minik=False, db=False):
    """Render utility room appliances: dryer, washer, counter, water heater, toilets, sinks."""
    pts = data.pts
    to_svg = data.to_svg

    # Dryer and washer
    minik_appl_w = 32.0 / 12.0   # 32" E-W in minik
    minik_appl_d = 27.0 / 12.0   # 27" N-S in minik
    minik_appl_links = {
        "DRYER": "https://www.lowes.com/pd/Electrolux-8-cu-ft-Stackable-Steam-Cycle-Electric-Dryer-Titanium-ENERGY-STAR/5015416377",
        "WASHER": "https://www.lowes.com/pd/Electrolux-Smartboost-Optic-Whites-and-Pure-Rinse-4-5-cu-ft-High-Efficiency-Stackable-Steam-Cycle-Front-Load-Washer-Titanium-ENERGY-STAR/5015416375",
    }
    _appl_shift_e = 4.0 / 12.0   # 4" east of layout position
    _appl_shift_s = -2.0 / 12.0  # 2" south of layout position
    minik_dryer_n = None
    _small_wd = minik or db
    for label, b in [("DRYER", layout.dryer), ("WASHER", layout.washer)]:
        if not _small_wd:
            b = BBox(w=b.w + _appl_shift_e, s=b.s + _appl_shift_s,
                     e=b.e + _appl_shift_e, n=b.n + _appl_shift_s)
        if _small_wd:
            if label == "DRYER":
                b = BBox(w=b.w, s=b.s, e=b.w + minik_appl_w, n=b.s + minik_appl_d)
                minik_dryer_n = b.n
            else:  # WASHER: 1" north of dryer
                ws = minik_dryer_n + 1.0 / 12.0
                b = BBox(w=b.w, s=ws, e=b.w + minik_appl_w, n=ws + minik_appl_d)
        link = minik_appl_links.get(label) if _small_wd else None
        if link:
            out.append(f'<a href="{link}" target="_blank">')
        sx1, sy1 = to_svg(b.w, b.n)
        sx2, sy2 = to_svg(b.e, b.s)
        sw = sx2 - sx1; sh = sy2 - sy1
        out.append(f'<rect x="{sx1:.1f}" y="{sy1:.1f}" width="{sw:.1f}" height="{sh:.1f}"'
                   f' fill="{APPL_FILL}" stroke="{APPL_STROKE}" stroke-width="{APPL_SW}"/>')
        cx, cy = (sx1 + sx2) / 2, (sy1 + sy2) / 2
        out.append(f'<text x="{cx:.1f}" y="{cy+3:.1f}" text-anchor="middle" font-family="Arial"'
                   f' font-size="7" fill="{APPL_STROKE}">{label}</text>')
        if link:
            out.append('</a>')

    # Counter: polygon clipped to W20-W0 south edge, IW3/IW16 west face east edge
    ctr_poly_svg = " ".join(f"{to_svg(*p)[0]:.1f},{to_svg(*p)[1]:.1f}"
                            for p in layout.ctr_poly)
    out.append(f'<polygon points="{ctr_poly_svg}"'
               f' fill="{APPL_FILL}" stroke="{APPL_STROKE}" stroke-width="{APPL_SW}"/>')
    _ctr_cx = sum(p[0] for p in layout.ctr_poly) / len(layout.ctr_poly)
    _ctr_cy = sum(p[1] for p in layout.ctr_poly) / len(layout.ctr_poly)
    ccx, ccy = to_svg(_ctr_cx, _ctr_cy)
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
    draw_toilet(out, toilet_e, layout.iw8.s, face_north=False, to_svg=to_svg)
    draw_sink(out, sink_e, layout.iw8.s - SINK_RY, to_svg=to_svg)
    draw_toilet(out, toilet_e, layout.iw8.n, face_north=True, to_svg=to_svg)
    draw_sink(out, sink_e, layout.iw8.n + SINK_RY, to_svg=to_svg)


def _render_kitchen(out, data, layout, minik=False, db=False):
    """Render kitchen: D/W, sink, stove, shelves, fridge, counters."""
    pts = data.pts
    to_svg = data.to_svg
    back_n = pts["W9"][1]

    # Kitchen appliances
    st_w = layout.iw2.e + NORTH_CTR_LENGTH + KITCHEN_APPL_GAP
    st_e = st_w + STOVE_WIDTH
    ks_w = st_e + KITCHEN_APPL_GAP
    ks_e = ks_w + KITCHEN_SINK_WIDTH
    dw_w = ks_e + KITCHEN_APPL_GAP
    dw_e = dw_w + DW_WIDTH
    appliances = [
        ("SINK",  ks_w, back_n - KITCHEN_SINK_DEPTH, ks_e, back_n,
         "https://www.webstaurantstore.com/advance-tabco-fs1181824l-45-fabricated-one-compartment-sink-with-24-left-drainboard-18-x-18-x-14-bowl/109FS1L241818.html"),
        ("STOVE", st_w, back_n - KITCHEN_APPL_GAP - STOVE_DEPTH, st_e, back_n - KITCHEN_APPL_GAP, None),
        ("D/W",   dw_w, back_n - DW_DEPTH,           dw_e, back_n, None),
    ]
    if minik:
        appliances = [(l, w, s, e, n, h) for l, w, s, e, n, h in appliances
                      if l not in ("STOVE", "D/W")]
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

    # Fridge
    if minik:
        # 3" east of SINK, 3" south of W9-W10 inner face
        fr_w = ks_e + 3.0 / 12.0
        fr_e = fr_w + MINIK_FRIDGE_W
        fr_n = back_n - 3.0 / 12.0
        fr_s = fr_n - MINIK_FRIDGE_D
    elif db:
        # 6" east of D/W, back 3" south of W9-W10 wall
        fr_w = dw_e + 6.0 / 12.0
        fr_e = fr_w + 32.75 / 12.0
        fr_n = back_n - 3.0 / 12.0
        fr_s = fr_n - 35.0 / 12.0
    else:
        # East edge 6" west of RO1, 2" north of IW1 north face
        ro1_w = layout.iw2.e + RO1_OFFSET_E_IW2
        fr_e = ro1_w - 6.0 / 12.0
        fr_s = layout.iw1_n + STD_GAP
        fr_w = fr_e - 32.75 / 12.0
        fr_n = fr_s + 35.0 / 12.0
    sx1, sy1 = to_svg(fr_w, fr_n)
    sx2, sy2 = to_svg(fr_e, fr_s)
    sw = sx2 - sx1; sh = sy2 - sy1
    if minik:
        out.append('<a href="https://www.ikea.com/us/en/p/bergsnaes-bottom-freezer-refrigerator-stainless-steel-color-60607883/" target="_blank">')
    else:
        out.append('<a href="https://www.lowes.com/pd/LG-25-5-cu-ft-Bottom-Freezer-Refrigerator-with-Ice-Maker-Fingerprint-Resistant-Printproof-Stainless-Steel-ENERGY-STAR/1002543648" target="_blank">')
    out.append(f'<rect x="{sx1:.1f}" y="{sy1:.1f}" width="{sw:.1f}" height="{sh:.1f}"'
               f' fill="{APPL_FILL}" stroke="{APPL_STROKE}" stroke-width="{APPL_SW}"/>')
    fr_cx = (sx1 + sx2) / 2
    fr_cy = (sy1 + sy2) / 2
    fr_fs = 6 if minik else 7
    out.append(f'<text x="{fr_cx:.1f}" y="{fr_cy+3:.1f}" text-anchor="middle" font-family="Arial"'
               f' font-size="{fr_fs}" fill="{APPL_STROKE}">FRIDGE</text>')
    if minik:
        # Door arc: hinged at SE corner, 23.375" door, sweeps from south (open) to west (closed)
        fr_door = MINIK_FRIDGE_W
        hx, hy = to_svg(fr_e, fr_s)
        tip_x, tip_y = to_svg(fr_e, fr_s - fr_door)
        out.append(f'<line x1="{hx:.1f}" y1="{hy:.1f}" x2="{tip_x:.1f}" y2="{tip_y:.1f}"'
                   f' stroke="{APPL_STROKE}" stroke-width="1.0"/>')
        n_arc = 20
        arc_pts = []
        for i in range(n_arc + 1):
            angle = math.pi + i * (math.pi / 2) / n_arc  # 180° to 270°
            ae = fr_e + fr_door * math.cos(angle)
            an = fr_s + fr_door * math.sin(angle)
            ax, ay = to_svg(ae, an)
            arc_pts.append(f"{ax:.1f},{ay:.1f}")
        out.append(f'<polyline points="{" ".join(arc_pts)}" fill="none"'
                   f' stroke="{APPL_STROKE}" stroke-width="0.5"/>')
        out.append('</a>')
    if db:
        # Door arc: hinged at SE corner, 32.75" door, sweeps from west to south
        fr_door = 32.75 / 12.0
        hx, hy = to_svg(fr_e, fr_s)
        tip_x, tip_y = to_svg(fr_e, fr_s - fr_door)
        out.append(f'<line x1="{hx:.1f}" y1="{hy:.1f}" x2="{tip_x:.1f}" y2="{tip_y:.1f}"'
                   f' stroke="{APPL_STROKE}" stroke-width="1.0"/>')
        n_arc = 20
        arc_pts = []
        for i in range(n_arc + 1):
            angle = math.pi + i * (math.pi / 2) / n_arc  # 180° to 270°
            ae = fr_e + fr_door * math.cos(angle)
            an = fr_s + fr_door * math.sin(angle)
            ax, ay = to_svg(ae, an)
            arc_pts.append(f"{ax:.1f},{ay:.1f}")
        out.append(f'<polyline points="{" ".join(arc_pts)}" fill="none"'
                   f' stroke="{APPL_STROKE}" stroke-width="0.5"/>')
        out.append('</a>')
    elif not minik:
        # Door arc: hinged at NE corner, 32.75" door, sweeps from west to north
        fr_door = 32.75 / 12.0
        hx, hy = to_svg(fr_e, fr_n)
        tip_x, tip_y = to_svg(fr_e, fr_n + fr_door)
        out.append(f'<line x1="{hx:.1f}" y1="{hy:.1f}" x2="{tip_x:.1f}" y2="{tip_y:.1f}"'
                   f' stroke="{APPL_STROKE}" stroke-width="1.0"/>')
        n_arc = 20
        arc_pts = []
        for i in range(n_arc + 1):
            angle = math.pi / 2 + i * (math.pi / 2) / n_arc  # 90° to 180°
            ae = fr_e + fr_door * math.cos(angle)
            an = fr_n + fr_door * math.sin(angle)
            ax, ay = to_svg(ae, an)
            arc_pts.append(f"{ax:.1f},{ay:.1f}")
        out.append(f'<polyline points="{" ".join(arc_pts)}" fill="none"'
                   f' stroke="{APPL_STROKE}" stroke-width="0.5"/>')
        out.append('</a>')

    if db:
        # ICE: in corner of IW1 and IW2, 2" from each wall
        ice_w = layout.iw2.e + 2.0 / 12.0
        ice_s = layout.iw1_n + 2.0 / 12.0
    elif minik:
        # ICE: 3" east of fridge, against W9-W10 wall (3" south)
        ice_w = fr_e + 3.0 / 12.0
        ice_s = back_n - 3.0 / 12.0 - ICE_DEPTH
    else:
        # ICE: 6" east of D/W, against W9-W10 wall
        ice_w = dw_e + 6.0 / 12.0
        ice_s = back_n - ICE_DEPTH
    ice_e = ice_w + ICE_WIDTH
    ice_n = ice_s + ICE_DEPTH
    ix1, iy1 = to_svg(ice_w, ice_n)
    ix2, iy2 = to_svg(ice_e, ice_s)
    isw = ix2 - ix1; ish = iy2 - iy1
    out.append('<a href="https://www.homedepot.com/p/EUHOMY-17-3-in-100-lb-24H-Full-Ice-Sizes-Commercial-Ice-Maker-in-Black-33-lb-Storage-Bin-Ice-Full-Alert-and-Auto-Cleaning-CIM001-100BL-E/337185876" target="_blank">')
    out.append(f'<rect x="{ix1:.1f}" y="{iy1:.1f}" width="{isw:.1f}" height="{ish:.1f}"'
               f' fill="{APPL_FILL}" stroke="{APPL_STROKE}" stroke-width="{APPL_SW}"/>')
    ice_cx = (ix1 + ix2) / 2
    ice_cy = (iy1 + iy2) / 2
    out.append(f'<text x="{ice_cx:.1f}" y="{ice_cy+3:.1f}" text-anchor="middle" font-family="Arial"'
               f' font-size="6" fill="{APPL_STROKE}">ICE</text>')
    out.append('</a>')

    # Kitchen counter: starting at IW2 east face (minik only)
    if minik:
        kc_w = layout.iw2.e
        kc_e = kc_w + KITCHEN_CTR_LENGTH
        # Against W9-W10 (north wall) and IW2
        kc_n = pts["W9"][1]
        kc_s = kc_n - KITCHEN_CTR_DEPTH
        kc_sx1, kc_sy1 = to_svg(kc_w, kc_n)
        kc_sx2, kc_sy2 = to_svg(kc_e, kc_s)
        kc_sw = kc_sx2 - kc_sx1; kc_sh = kc_sy2 - kc_sy1
        out.append('<a href="https://www.webstaurantstore.com/regency-spec-line-30-x-72-14-gauge-stainless-steel-commercial-work-table-with-4-backsplash-and-undershelf/600TSSB3072S.html" target="_blank">')
        out.append(f'<rect x="{kc_sx1:.1f}" y="{kc_sy1:.1f}" width="{kc_sw:.1f}" height="{kc_sh:.1f}"'
                   f' fill="{APPL_FILL}" stroke="{APPL_STROKE}" stroke-width="{APPL_SW}"/>')
        out.append('</a>')

    # Minik: microwave on counter (19.5" E-W x 16-5/8" N-S)
    if minik:
        mw_ew = 19.5 / 12.0
        mw_ns = 16.625 / 12.0
        mw_w = kc_w + 2.0 / 12.0
        mw_e = mw_w + mw_ew
        mw_n = kc_n - 3.0 / 12.0
        mw_s = mw_n - mw_ns
        mw_sx1, mw_sy1 = to_svg(mw_w, mw_n)
        mw_sx2, mw_sy2 = to_svg(mw_e, mw_s)
        mw_sw = mw_sx2 - mw_sx1
        mw_sh = mw_sy2 - mw_sy1
        out.append('<a href="https://www.ikea.com/us/en/p/gatebo-microwave-oven-with-air-fryer-function-ikea-500-black-70603506/" target="_blank">')
        out.append(f'<rect x="{mw_sx1:.1f}" y="{mw_sy1:.1f}" width="{mw_sw:.1f}" height="{mw_sh:.1f}"'
                   f' fill="{APPL_FILL}" stroke="{APPL_STROKE}" stroke-width="{APPL_SW}"/>')
        mw_cx = (mw_sx1 + mw_sx2) / 2
        mw_cy = (mw_sy1 + mw_sy2) / 2
        out.append(f'<text x="{mw_cx:.1f}" y="{mw_cy+3:.1f}" text-anchor="middle" font-family="Arial"'
                   f' font-size="5" fill="{APPL_STROKE}">MICRO</text>')
        out.append('</a>')

    # Minik: coffee maker on counter (7.2" E-W x 9.2" N-S)
    if minik:
        cm_ew = 7.2 / 12.0
        cm_ns = 9.2 / 12.0
        cm_w = mw_e + 3.0 / 12.0
        cm_e = cm_w + cm_ew
        cm_n = kc_n - 3.0 / 12.0
        cm_s = cm_n - cm_ns
        cm_sx1, cm_sy1 = to_svg(cm_w, cm_n)
        cm_sx2, cm_sy2 = to_svg(cm_e, cm_s)
        cm_sw = cm_sx2 - cm_sx1
        cm_sh = cm_sy2 - cm_sy1
        out.append('<a href="https://www.amazon.com/Holstein-Housewares-HH-0914701E-5-Cup-Coffee/dp/B08HSRCC4T/?th=1" target="_blank">')
        out.append(f'<rect x="{cm_sx1:.1f}" y="{cm_sy1:.1f}" width="{cm_sw:.1f}" height="{cm_sh:.1f}"'
                   f' fill="{APPL_FILL}" stroke="{APPL_STROKE}" stroke-width="{APPL_SW}"/>')
        cm_cx = (cm_sx1 + cm_sx2) / 2
        cm_cy = (cm_sy1 + cm_sy2) / 2
        out.append(f'<text x="{cm_cx:.1f}" y="{cm_cy+3:.1f}" text-anchor="middle" font-family="Arial"'
                   f' font-size="5" fill="{APPL_STROKE}">C</text>')
        out.append('</a>')

    # Minik: induction cooktop on counter (rotated: 13.4" E-W x 20.5" N-S)
    if minik:
        cp_ew = 13.4 / 12.0
        cp_ns = 20.5 / 12.0
        cp_w = cm_e + 3.0 / 12.0
        cp_e = cp_w + cp_ew
        cp_cx_e = (cp_w + cp_e) / 2
        cp_s = kc_s + 2.0 / 12.0
        cp_n = cp_s + cp_ns
        cp_cy_n = (cp_s + cp_n) / 2
        cp_sx1, cp_sy1 = to_svg(cp_w, cp_n)
        cp_sx2, cp_sy2 = to_svg(cp_e, cp_s)
        cp_sw = cp_sx2 - cp_sx1
        cp_sh = cp_sy2 - cp_sy1
        cp_r = abs(to_svg(1.0 / 12.0, 0)[0] - to_svg(0, 0)[0])  # 1" corner radius
        out.append('<a href="https://www.homedepot.com/p/Empava-Portable-13-4-in-Induction-Electric-Cooktop-in-Black-with-2-Elements-EMPV-ID12/313815692" target="_blank">')
        out.append(f'<rect x="{cp_sx1:.1f}" y="{cp_sy1:.1f}" width="{cp_sw:.1f}" height="{cp_sh:.1f}"'
                   f' rx="{cp_r:.1f}" ry="{cp_r:.1f}"'
                   f' fill="#222" stroke="{APPL_STROKE}" stroke-width="{APPL_SW}"/>')
        # Two burner circles (6" diameter each), spaced evenly N-S
        burner_r_ft = 3.0 / 12.0
        burner_r_svg = abs(to_svg(burner_r_ft, 0)[0] - to_svg(0, 0)[0])
        for sign in (-1, 1):
            bx = cp_cx_e
            by = cp_cy_n + sign * cp_ns / 4
            bsx, bsy = to_svg(bx, by)
            out.append(f'<circle cx="{bsx:.1f}" cy="{bsy:.1f}" r="{burner_r_svg:.1f}"'
                       f' fill="none" stroke="#666" stroke-width="0.3"/>')
        out.append('</a>')

    # Minik: toaster 3" east of cooktop, 3" south of W9-W10 (13.7" E-W x 12.5" N-S)
    if minik:
        ts_ew = 13.7 / 12.0
        ts_ns = 12.5 / 12.0
        ts_w = cp_e + 3.0 / 12.0
        ts_e = ts_w + ts_ew
        ts_n = kc_n - 3.0 / 12.0
        ts_s = ts_n - ts_ns
        ts_sx1, ts_sy1 = to_svg(ts_w, ts_n)
        ts_sx2, ts_sy2 = to_svg(ts_e, ts_s)
        ts_sw = ts_sx2 - ts_sx1
        ts_sh = ts_sy2 - ts_sy1
        out.append('<a href="https://www.amazon.com/Roter-Mond-Stainless-Independent-Removable/dp/B0CGTQZTDZ?th=1" target="_blank">')
        out.append(f'<rect x="{ts_sx1:.1f}" y="{ts_sy1:.1f}" width="{ts_sw:.1f}" height="{ts_sh:.1f}"'
                   f' fill="{APPL_FILL}" stroke="{APPL_STROKE}" stroke-width="{APPL_SW}"/>')
        # Four toast slot lines (1.5" wide slots, evenly spaced E-W)
        slot_len_ft = 5.5 / 12.0  # 5.5" slot length
        ts_cx_e = (ts_w + ts_e) / 2
        ts_cy_n = (ts_s + ts_n) / 2
        slot_spacing = ts_ew / 5.0  # 5 gaps for 4 slots
        for i in range(4):
            sl_e = ts_w + slot_spacing * (i + 1)
            sl_n = ts_cy_n + slot_len_ft / 2
            sl_s = ts_cy_n - slot_len_ft / 2
            sl_sx, sl_sy1 = to_svg(sl_e, sl_n)
            _, sl_sy2 = to_svg(sl_e, sl_s)
            out.append(f'<line x1="{sl_sx:.1f}" y1="{sl_sy1:.1f}" x2="{sl_sx:.1f}" y2="{sl_sy2:.1f}"'
                       f' stroke="#666" stroke-width="0.4"/>')
        out.append('</a>')

    # Oscar triangle dining set centered between north wall, IW1, IW2, RO1
    ro1_w_pos = layout.iw2.e + RO1_OFFSET_E_IW2
    space_w = layout.iw2.e
    space_s = layout.iw1_n
    space_e = ro1_w_pos
    space_n = kc_s if minik else back_n
    space_cx = (space_w + space_e) / 2
    space_cy = (space_s + space_n) / 2

    # Table: base 31.5" (N), height 35.25", 24" arc at apex, 6" fillets
    tbl_base = 31.5 / 12.0
    tbl_h = 35.25 / 12.0
    apex_r = 12.0 / 12.0    # 24" diameter arc at apex
    fillet_r = 6.0 / 12.0   # 6" corner fillets

    # Position: north side 30" south of space north, centered E-W
    if db:
        # Center under SINK west end
        tbl_cx = ks_w
    elif not minik:
        # Center between fridge west and IW2 east + 1.125" east
        _ro1_w = layout.iw2.e + RO1_OFFSET_E_IW2
        _fr_w = _ro1_w - 6.0 / 12.0 - 32.75 / 12.0
        tbl_cx = (_fr_w + layout.iw2.e) / 2 + 1.125 / 12.0
    else:
        # Center on SINK west end
        _st_w = layout.iw2.e + NORTH_CTR_LENGTH + KITCHEN_APPL_GAP
        tbl_cx = _st_w + STOVE_WIDTH + KITCHEN_APPL_GAP
    tbl_n = space_n - 30.0 / 12.0 - (28.0 / 12.0 if not minik else 0)
    tbl_s_y = tbl_n - tbl_h

    # Base corners and arc center
    ne = (tbl_cx + tbl_base / 2, tbl_n)
    nw = (tbl_cx - tbl_base / 2, tbl_n)
    arc_c = (tbl_cx, tbl_s_y + apex_r)  # 12" north of south edge

    # Right tangent from NE to apex arc
    dx_r = ne[0] - arc_c[0]
    dn_r = ne[1] - arc_c[1]
    dist_r = math.sqrt(dx_r**2 + dn_r**2)
    angle_cp = math.atan2(dn_r, dx_r)
    delta = math.acos(apex_r / dist_r)
    alpha_r = angle_cp - delta
    t_right = (arc_c[0] + apex_r * math.cos(alpha_r),
                arc_c[1] + apex_r * math.sin(alpha_r))
    t_left = (2 * tbl_cx - t_right[0], t_right[1])

    # NE fillet between base (west) and tangent line (toward t_right)
    d_base_ne = (-1.0, 0.0)
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
    f_ne_base = (fc_ne[0], tbl_n)  # tangent to base
    v_ne = (fc_ne[0] - ne[0], fc_ne[1] - ne[1])
    t_proj = v_ne[0] * d_tang_ne[0] + v_ne[1] * d_tang_ne[1]
    f_ne_tang = (ne[0] + t_proj * d_tang_ne[0], ne[1] + t_proj * d_tang_ne[1])

    # NW fillet by symmetry
    f_nw_base = (2 * tbl_cx - f_ne_base[0], tbl_n)
    f_nw_tang = (2 * tbl_cx - f_ne_tang[0], f_ne_tang[1])

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
        to_ctr = (tbl_cx - mid_e, space_cy - mid_n)
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

    # North wall counter: south side against W9-W10, starting at IW2 east face
    if not minik:
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

def _render_furniture(out, data, layout, minik=False, db=False):
    """Render furniture: bed, loveseat/sofa, ET, chair, ottoman, room labels."""
    pts = data.pts
    to_svg = data.to_svg
    # Bed (rotated polygon)
    _bp = layout.bed_poly
    _bp_svg = " ".join(f"{to_svg(*p)[0]:.1f},{to_svg(*p)[1]:.1f}" for p in _bp)
    out.append(f'<polygon points="{_bp_svg}" fill="{APPL_FILL}" stroke="{APPL_STROKE}" stroke-width="{APPL_SW}"/>')
    _bp_cx = sum(p[0] for p in _bp) / 4
    _bp_cy = sum(p[1] for p in _bp) / 4
    # Rotate label 90° CW from bed long axis (now parallel to W20-W0)
    _bed_dx = to_svg(*_bp[2])[0] - to_svg(*_bp[1])[0]
    _bed_dy = to_svg(*_bp[2])[1] - to_svg(*_bp[1])[1]
    _bed_ang = math.degrees(math.atan2(_bed_dy, _bed_dx)) + 90
    # Position at 1/2 the perpendicular distance from W20-W0 to IW1
    _w20b = pts["W20"]; _w0b = pts["W0"]
    _dEwb = _w0b[0] - _w20b[0]; _dNwb = _w0b[1] - _w20b[1]
    _wlb = math.sqrt(_dEwb**2 + _dNwb**2)
    _nEb = _dNwb / _wlb; _nNb = -_dEwb / _wlb  # inward normal
    _ucx = _bp_cx - _w20b[0]; _ucy = _bp_cy - _w20b[1]
    _tw = (_ucx * _dEwb + _ucy * _dNwb) / (_dEwb**2 + _dNwb**2)
    _wb = (_w20b[0] + _tw * _dEwb, _w20b[1] + _tw * _dNwb)
    _d_iw1 = (layout.iw1_s - _wb[1]) / _nNb
    _lbl_e = _wb[0] + _d_iw1 / 2 * _nEb
    _lbl_n = _wb[1] + _d_iw1 / 2 * _nNb
    _bsx, _bsy = to_svg(_lbl_e, _lbl_n)
    out.append(f'<text x="{_bsx:.1f}" y="{_bsy+3:.1f}" text-anchor="middle" font-family="Arial"'
               f' font-size="7" fill="{APPL_STROKE}" transform="rotate({_bed_ang:.1f},{_bsx:.1f},{_bsy+3:.1f})">'
               f'KING BED</text>')

    if minik:
        # SOFA: 73.2" E-W x 24.6" N-S, 6" east of IW4 west, 2" north of IW1
        sofa_w = layout.iw4_w + 6.0 / 12.0
        sofa_e = sofa_w + SOFA_WIDTH - 24.0 / 12.0
        sofa_s = layout.iw1_n + 2.0 / 12.0
        sofa_n = sofa_s + SOFA_DEPTH
        sf_sx1, sf_sy1 = to_svg(sofa_w, sofa_n)
        sf_sx2, sf_sy2 = to_svg(sofa_e, sofa_s)
        sf_sw = sf_sx2 - sf_sx1; sf_sh = sf_sy2 - sf_sy1
        out.append('<a href="https://www.homedepot.com/p/AURA-OUTDOOR-4-Piece-Metal-Outdoor-Sectional-Sofa-Set-Patio-Furniture-Set-with-6-in-Olefin-Cushion-Gray-SIS006-GY/335858535" target="_blank">')
        out.append(f'<rect x="{sf_sx1:.1f}" y="{sf_sy1:.1f}" width="{sf_sw:.1f}" height="{sf_sh:.1f}"'
                   f' fill="{APPL_FILL}" stroke="{APPL_STROKE}" stroke-width="{APPL_SW}"/>')
        sf_cx = (sf_sx1 + sf_sx2) / 2
        sf_cy = (sf_sy1 + sf_sy2) / 2
        out.append(f'<text x="{sf_cx:.1f}" y="{sf_cy+3:.1f}" text-anchor="middle" font-family="Arial"'
                   f' font-size="6" fill="{APPL_STROKE}">SOFA</text>')
        out.append('</a>')

        # ROCKER: midpoint between ICE SE corner and SOFA NW corner
        # Recompute ICE SE: fridge_e = ks_e + 3" + MINIK_FRIDGE_W
        _st_e = layout.iw2.e + NORTH_CTR_LENGTH + KITCHEN_APPL_GAP + STOVE_WIDTH
        _ks_e = _st_e + KITCHEN_APPL_GAP + KITCHEN_SINK_WIDTH
        _fr_e = _ks_e + 3.0 / 12.0 + MINIK_FRIDGE_W
        _ice_e = _fr_e + 3.0 / 12.0 + ICE_WIDTH
        _back_n = pts["W9"][1]
        _ice_s = _back_n - 3.0 / 12.0 - ICE_DEPTH
        rk_cx = (_ice_e + sofa_w) / 2
        rk_cy = (_ice_s + sofa_n) / 2 - 18.0 / 12.0
        rk_hw = ROCKER_DEPTH / 2   # half E-W (rotated 90°)
        rk_hh = ROCKER_WIDTH / 2   # half N-S (rotated 90°)
        rk_r = ROCKER_CORNER_R
        rk_scx, rk_scy = to_svg(rk_cx, rk_cy)
        rk_sx1, rk_sy1 = to_svg(rk_cx - rk_hw, rk_cy + rk_hh)
        rk_sx2, rk_sy2 = to_svg(rk_cx + rk_hw, rk_cy - rk_hh)
        rk_sw = rk_sx2 - rk_sx1; rk_sh = rk_sy2 - rk_sy1
        rk_sr = abs(to_svg(rk_r, 0)[0] - to_svg(0, 0)[0])
        rk_angle = -15.0  # 15° CCW in plan = -15° in SVG (Y-axis flipped)
        out.append(f'<a href="https://www.ikea.com/us/en/p/poaeng-rocking-chair-brown-gunnared-beige-s39502048/" target="_blank">')
        out.append(f'<g transform="rotate({rk_angle:.1f},{rk_scx:.1f},{rk_scy:.1f})">')
        out.append(f'<rect x="{rk_sx1:.1f}" y="{rk_sy1:.1f}" width="{rk_sw:.1f}" height="{rk_sh:.1f}"'
                   f' rx="{rk_sr:.1f}" ry="{rk_sr:.1f}"'
                   f' fill="{APPL_FILL}" stroke="{APPL_STROKE}" stroke-width="{APPL_SW}"/>')
        out.append(f'<text x="{rk_scx:.1f}" y="{rk_scy+3:.1f}" text-anchor="middle" font-family="Arial"'
                   f' font-size="6" fill="{APPL_STROKE}">ROCKER</text>')
        out.append('</g>')
        out.append('</a>')
    elif db:
        # DB variant: no loveseats; ET shifted east to 2" from W-series wall
        et_r = (ET_RADIUS_CM / 2.54) / 12.0
        et_cy = layout.iw1_n + STD_GAP + et_r
        # Find easternmost W-series wall easting at ET's northing
        wall_e = max(horiz_isects(data.inner_poly, et_cy))
        et_cx = wall_e - STD_GAP - et_r

        # ET: 50cm diameter endtable
        et_sx, et_sy = to_svg(et_cx, et_cy)
        et_r_svg = abs(to_svg(et_r, 0)[0] - to_svg(0, 0)[0])
        out.append('<a href="https://www.ikea.com/us/en/p/listerby-side-table-oak-veneer-30515314/" target="_blank">')
        out.append(f'<circle cx="{et_sx:.1f}" cy="{et_sy:.1f}" r="{et_r_svg:.1f}"'
                   f' fill="{APPL_FILL}" stroke="{APPL_STROKE}" stroke-width="{APPL_SW}"/>')
        out.append(f'<text x="{et_sx:.1f}" y="{et_sy+3:.1f}" text-anchor="middle"'
                   f' font-family="Arial" font-size="6" fill="{APPL_STROKE}">ET</text>')
        out.append('</a>')
    else:
        # Loveseat: 35" E-W x 65" N-S, rotated 15° CCW about SW corner
        lv_width = LOVESEAT_WIDTH
        lv_height = LOVESEAT_LENGTH
        lv_angle = math.radians(LOVESEAT_ANGLE_DEG)
        lv_nw_e = LOVESEAT_NW_E
        lv_nw_n = LOVESEAT_NW_N
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

    # DESK: 60" x 30", along W16-W17 wall (30° from horizontal)
    dk_w17 = pts["W17"]
    dk_sw_e = dk_w17[0]
    dk_sw_n = dk_w17[1]
    dk_ne_e = dk_sw_e + DESK_WIDTH
    dk_ne_n = dk_sw_n + DESK_DEPTH
    dk_sx1, dk_sy1 = to_svg(dk_sw_e, dk_ne_n)
    dk_sx2, dk_sy2 = to_svg(dk_ne_e, dk_sw_n)
    dk_sw = dk_sx2 - dk_sx1; dk_sh = dk_sy2 - dk_sy1
    dk_rot_x = dk_sx1
    dk_rot_y = dk_sy2
    out.append(f'<g transform="rotate(-30,{dk_rot_x:.1f},{dk_rot_y:.1f})">')
    out.append(f'<rect x="{dk_sx1:.1f}" y="{dk_sy1:.1f}" width="{dk_sw:.1f}" height="{dk_sh:.1f}"'
               f' fill="{APPL_FILL}" stroke="{APPL_STROKE}" stroke-width="{APPL_SW}"/>')
    dk_cx = (dk_sx1 + dk_sx2) / 2
    dk_cy = (dk_sy1 + dk_sy2) / 2
    out.append(f'<text x="{dk_cx:.1f}" y="{dk_cy+3:.1f}" text-anchor="middle" font-family="Arial"'
               f' font-size="7" fill="{APPL_STROKE}">DESK</text>')
    out.append('</g>')

    # DESK CHAIR: 27" x 24", 3" rounded corners, 12" in front of desk, centered
    dc_w = dk_sw_e + DESK_WIDTH / 2 - DESK_CHAIR_WIDTH / 2
    dc_e = dc_w + DESK_CHAIR_WIDTH
    dc_s = dk_sw_n + DESK_DEPTH + DESK_CHAIR_GAP
    dc_n = dc_s + DESK_CHAIR_DEPTH
    dc_sx1, dc_sy1 = to_svg(dc_w, dc_n)
    dc_sx2, dc_sy2 = to_svg(dc_e, dc_s)
    dc_sw = dc_sx2 - dc_sx1; dc_sh = dc_sy2 - dc_sy1
    dc_r_svg = abs(to_svg(CHAIR_CORNER_R, 0)[0] - to_svg(0, 0)[0])
    out.append(f'<g transform="rotate(-30,{dk_rot_x:.1f},{dk_rot_y:.1f})">')
    out.append('<a href="https://www.amazon.com/BESTFAIR-Ergonomic-Office-Chair-Adjustable/dp/B0FDQDMP2D?th=1" target="_blank">')
    out.append(f'<rect x="{dc_sx1:.1f}" y="{dc_sy1:.1f}" width="{dc_sw:.1f}" height="{dc_sh:.1f}"'
               f' rx="{dc_r_svg:.1f}" ry="{dc_r_svg:.1f}"'
               f' fill="{APPL_FILL}" stroke="{APPL_STROKE}" stroke-width="{APPL_SW}"/>')
    dc_cx = (dc_sx1 + dc_sx2) / 2
    dc_cy = (dc_sy1 + dc_sy2) / 2
    out.append(f'<text x="{dc_cx:.1f}" y="{dc_cy+3:.1f}" text-anchor="middle" font-family="Arial"'
               f' font-size="7" fill="{APPL_STROKE}">CHAIR</text>')
    out.append('</a>')
    out.append('</g>')

    # Room labels
    bd_cx = layout.iw2.e + 139.0 / 12.0  # 11'7" east of IW2 east face
    bd_cy = (layout.ctr.s + layout.iw1_s) / 2
    bdx, bdy = to_svg(bd_cx, bd_cy)
    out.append(f'<text x="{bdx:.1f}" y="{bdy+3:.1f}" text-anchor="middle" font-family="Arial"'
               f' font-size="8" fill="#666">BEDROOM</text>')

    of_cx = (layout.iw4_e + pts["W15"][0]) / 2
    of_cy = (layout.ctr.s + 5.0 + layout.iwt3 + layout.iw1_s) / 2 - 2.0 + 8.0 / 12.0
    ofx, ofy = to_svg(of_cx, of_cy)
    out.append(f'<text x="{ofx:.1f}" y="{ofy+3:.1f}" text-anchor="middle" font-family="Arial"'
               f' font-size="8" fill="#666">OFFICE</text>')


def _render_dimensions(out, data, layout):
    """Render all dimension lines (interior and external)."""
    pts = data.pts
    to_svg = data.to_svg

    # IW1-north → F9-F11 south face dimension
    iw1_n = layout.iw1_n
    dim_e = (pts["F9"][0] + pts["F11"][0]) / 2
    dim_line_v(out, dim_e, iw1_n, pts["W9"][1], fmt_dist(pts["W9"][1] - iw1_n), to_svg)

    # IW2-east → inside F12-F13 wall dimension
    dim2_n = (pts["F12"][1] + pts["F13"][1]) / 2
    w13, w12 = pts["W13"], pts["W12"]
    t_e = (dim2_n - w13[1]) / (w12[1] - w13[1]) if w12[1] != w13[1] else 0.5
    dim2_east_e = w13[0] + t_e * (w12[0] - w13[0])
    dim_line_h(out, layout.iw2.e, dim2_n, dim2_east_e, fmt_dist(dim2_east_e - layout.iw2.e), to_svg,
               label_offset_e=-4.0)

    # East closet (rotated dimension, parallel to IW11)
    _iw12_sw = layout.iw12_poly[0]
    _iw12_se = layout.iw12_poly[1]
    _dim_s = ((_iw12_sw[0] + _iw12_se[0]) / 2,
              (_iw12_sw[1] + _iw12_se[1]) / 2)
    _dn = (layout.iw11_poly[2][0] - layout.iw11_poly[1][0],
           layout.iw11_poly[2][1] - layout.iw11_poly[1][1])
    _dl = math.sqrt(_dn[0]**2 + _dn[1]**2)
    _nrm = (_dn[0] / _dl, _dn[1] / _dl)
    # Ray-line intersection: ray from _dim_s in direction -_nrm hits W20-W0
    _w20 = pts["W20"]; _w0 = pts["W0"]
    _dw = (_w0[0] - _w20[0], _w0[1] - _w20[1])
    _u = (_dim_s[0] - _w20[0], _dim_s[1] - _w20[1])
    _det = _nrm[0] * _dw[1] - _nrm[1] * _dw[0]
    _t_line = (_u[0] * _dw[1] - _u[1] * _dw[0]) / _det
    _dim_e = (_dim_s[0] - _t_line * _nrm[0], _dim_s[1] - _t_line * _nrm[1])
    _dsx1, _dsy1 = to_svg(*_dim_s)
    _dsx2, _dsy2 = to_svg(*_dim_e)
    _sdx = _dsx2 - _dsx1; _sdy = _dsy2 - _dsy1
    _slen = math.sqrt(_sdx**2 + _sdy**2)
    _px = -_sdy / _slen; _py = _sdx / _slen
    _tk = 4
    out.append(f'<line x1="{_dsx1:.1f}" y1="{_dsy1:.1f}" x2="{_dsx2:.1f}" y2="{_dsy2:.1f}" stroke="{DIM_COLOR}" stroke-width="0.8"/>')
    for _sx, _sy in [(_dsx1, _dsy1), (_dsx2, _dsy2)]:
        out.append(f'<line x1="{_sx - _tk * _px:.1f}" y1="{_sy - _tk * _py:.1f}" '
                   f'x2="{_sx + _tk * _px:.1f}" y2="{_sy + _tk * _py:.1f}" '
                   f'stroke="{DIM_COLOR}" stroke-width="0.8"/>')
    _lmx = (_dsx1 + _dsx2) / 2; _lmy = (_dsy1 + _dsy2) / 2
    _up_dx = _dsx1 - _dsx2; _up_dy = _dsy1 - _dsy2
    _up_ang = math.degrees(math.atan2(_up_dy, _up_dx))
    _lx = _lmx + 3 * _up_dy / _slen; _ly = _lmy - 3 * _up_dx / _slen
    out.append(f'<text x="{_lx:.1f}" y="{_ly:.1f}" text-anchor="middle" font-family="Arial" '
               f'font-size="8" fill="{DIM_COLOR}" transform="rotate({_up_ang:.1f},{_lx:.1f},{_ly:.1f})">'
               f'CLOSET {fmt_dist(_t_line)}</text>')

    # West closet (rotated dimension, perpendicular to IW7 south face)
    _iw7_sw7 = layout.iw7_poly[0]
    _iw7_se7 = layout.iw7_poly[1]
    _dim7_s = ((_iw7_sw7[0] + _iw7_se7[0]) / 2,
               (_iw7_sw7[1] + _iw7_se7[1]) / 2)
    _dn7 = (layout.iw7_poly[3][0] - _iw7_sw7[0],
            layout.iw7_poly[3][1] - _iw7_sw7[1])
    _dl7 = math.sqrt(_dn7[0]**2 + _dn7[1]**2)
    _nrm7 = (_dn7[0] / _dl7, _dn7[1] / _dl7)
    _w20 = pts["W20"]; _w0 = pts["W0"]
    _dw7 = (_w0[0] - _w20[0], _w0[1] - _w20[1])
    _u7 = (_dim7_s[0] - _w20[0], _dim7_s[1] - _w20[1])
    _det7 = _nrm7[0] * _dw7[1] - _nrm7[1] * _dw7[0]
    _t7 = (_u7[0] * _dw7[1] - _u7[1] * _dw7[0]) / _det7
    _dim7_e = (_dim7_s[0] - _t7 * _nrm7[0], _dim7_s[1] - _t7 * _nrm7[1])
    _dsx1, _dsy1 = to_svg(*_dim7_s)
    _dsx2, _dsy2 = to_svg(*_dim7_e)
    _sdx = _dsx2 - _dsx1; _sdy = _dsy2 - _dsy1
    _slen = math.sqrt(_sdx**2 + _sdy**2)
    _px = -_sdy / _slen; _py = _sdx / _slen
    _tk = 4
    out.append(f'<line x1="{_dsx1:.1f}" y1="{_dsy1:.1f}" x2="{_dsx2:.1f}" y2="{_dsy2:.1f}" stroke="{DIM_COLOR}" stroke-width="0.8"/>')
    for _sx, _sy in [(_dsx1, _dsy1), (_dsx2, _dsy2)]:
        out.append(f'<line x1="{_sx - _tk * _px:.1f}" y1="{_sy - _tk * _py:.1f}" '
                   f'x2="{_sx + _tk * _px:.1f}" y2="{_sy + _tk * _py:.1f}" '
                   f'stroke="{DIM_COLOR}" stroke-width="0.8"/>')
    _lmx = (_dsx1 + _dsx2) / 2; _lmy = (_dsy1 + _dsy2) / 2
    _up_dx = _dsx1 - _dsx2; _up_dy = _dsy1 - _dsy2
    _up_ang = math.degrees(math.atan2(_up_dy, _up_dx))
    _lx = _lmx + 3 * _up_dy / _slen; _ly = _lmy - 3 * _up_dx / _slen
    out.append(f'<text x="{_lx:.1f}" y="{_ly:.1f}" text-anchor="middle" font-family="Arial" '
               f'font-size="8" fill="{DIM_COLOR}" transform="rotate({_up_ang:.1f},{_lx:.1f},{_ly:.1f})">'
               f'CLOSET {fmt_dist(_t7)}</text>')

    # Utility: W1 to IW3 west face, at northing where distance = 8'
    _iw3_sw, _iw3_nw = layout.iw3_poly[0], layout.iw3_poly[3]
    _target_dist = 8.0
    _t_iw3 = (_target_dist - (_iw3_sw[0] - pts["W1"][0])) / (_iw3_nw[0] - _iw3_sw[0])
    _dim_n = _iw3_sw[1] + _t_iw3 * (_iw3_nw[1] - _iw3_sw[1])
    _dim_e = _iw3_sw[0] + _t_iw3 * (_iw3_nw[0] - _iw3_sw[0])
    dim_line_h(out, pts["W1"][0], _dim_n, _dim_e,
               fmt_dist(_dim_e - pts["W1"][0]), to_svg)
    dim_line_h(out, layout.iw4_e, 5.0, pts["W15"][0],
               fmt_dist(pts["W15"][0] - layout.iw4_e), to_svg)

    # Storage — west end at IW15 east face
    _stor_n = (layout.iw5.n + layout.iw1_s) / 2
    _stor_w = layout.iw15.e
    dim_line_h(out, _stor_w, _stor_n, pts["W15"][0],
               f"STORAGE {fmt_dist(pts['W15'][0] - _stor_w)}", to_svg)

    # O1 east face center to RO3 west face center
    _o1 = compute_outer_openings(pts, layout)[0]  # O1
    _o1_e_ctr = ((_o1.poly[2][0] + _o1.poly[3][0]) / 2,
                 (_o1.poly[2][1] + _o1.poly[3][1]) / 2)
    _ro3 = [r for r in compute_rough_openings(pts, layout) if r.name == "RO3"][0]
    _ro3_w_ctr = (_ro3.bbox.w, (_ro3.bbox.s + _ro3.bbox.n) / 2)
    dim_line_h(out, _o1_e_ctr[0], _o1_e_ctr[1], _ro3_w_ctr[0],
               fmt_dist(_ro3_w_ctr[0] - _o1_e_ctr[0]), to_svg)

    # West wall interior widths
    dim_line_h(out, pts["W2"][0], pts["F2"][1], layout.iw2.w,
               fmt_dist(layout.iw2.w - pts["W2"][0]), to_svg)
    dim_line_h(out, pts["W5"][0], pts["F5"][1], layout.iw2.w,
               fmt_dist(layout.iw2.w - pts["W5"][0]), to_svg)

    # Office/bedroom verticals
    dim_line_v(out, pts["F18"][0], layout.iw1_s, pts["W18"][1],
               fmt_dist(layout.iw1_s - pts["W18"][1]), to_svg,
               label_n=(layout.iw1_s + pts["W18"][1]) / 2 + 2.5)
    dim_line_v(out, pts["F6"][0] + 1.0, layout.iw6_n, pts["W6"][1],
               fmt_dist(pts["W6"][1] - layout.iw6_n), to_svg)
    dim_line_v(out, pts["F6"][0] + 1.0, layout.iw8.n, layout.iw6_s,
               fmt_dist(layout.iw6_s - layout.iw8.n), to_svg)

    # External dimensions
    dim_ext_e = pts["F2"][0] - 2.7
    dim_line_v(out, dim_ext_e, pts["F0"][1], pts["F6"][1],
               fmt_dist(pts["F6"][1] - pts["F0"][1]), to_svg)

    # Top exterior dim: endpoints on arcs F7-F8 and F11-F12, 4" south of F7
    _dim_n_arc = pts["F7"][1] - 4.0 / 12.0
    # West end: point on arc F7-F8 (center C7, radius R_a7) at that northing
    _c7 = pts["C7"]; _r7 = data.radii["R_a7"]
    _sin7 = (_dim_n_arc - _c7[1]) / _r7
    _dim_w_e = _c7[0] + _r7 * math.sqrt(1.0 - _sin7 ** 2)
    # East end: westernmost point on arc F11-F12 (center C11, radius R_a11) at that northing
    _c11 = pts["C11"]; _r11 = data.radii["R_a11"]
    _sin11 = (_dim_n_arc - _c11[1]) / _r11
    _dim_e_e = _c11[0] - _r11 * math.sqrt(1.0 - _sin11 ** 2)
    dim_line_h(out, _dim_w_e, pts["F6"][1] + 1.0, _dim_e_e,
               fmt_dist(_dim_e_e - _dim_w_e), to_svg)

    dim_ext_n = pts["F19"][1] - 3.0
    dim_line_h(out, pts["F1"][0], dim_ext_n, pts["F15"][0],
               fmt_dist(pts["F15"][0] - pts["F1"][0]), to_svg)

    # Bedroom: O9 inner center perpendicular to W20-W0 up to IW1 south face
    _o9_open = compute_outer_openings(pts, layout)[8]  # O9
    _o9_ic = ((_o9_open.poly[2][0] + _o9_open.poly[3][0]) / 2,
              (_o9_open.poly[2][1] + _o9_open.poly[3][1]) / 2)
    _dEw = pts["W0"][0] - pts["W20"][0]
    _dNw = pts["W0"][1] - pts["W20"][1]
    _wlen = math.sqrt(_dEw**2 + _dNw**2)
    _nrmE = _dNw / _wlen; _nrmN = -_dEw / _wlen  # inward normal (NNE)
    _t_iw1 = (layout.iw1_s - _o9_ic[1]) / _nrmN
    _dim_end = (_o9_ic[0] + _t_iw1 * _nrmE, layout.iw1_s)
    _dim_len = abs(_t_iw1)
    _dsx1, _dsy1 = to_svg(*_o9_ic)
    _dsx2, _dsy2 = to_svg(*_dim_end)
    _sdx = _dsx2 - _dsx1; _sdy = _dsy2 - _dsy1
    _slen = math.sqrt(_sdx**2 + _sdy**2)
    _px = -_sdy / _slen; _py = _sdx / _slen
    _tk = 4
    out.append(f'<line x1="{_dsx1:.1f}" y1="{_dsy1:.1f}" x2="{_dsx2:.1f}" y2="{_dsy2:.1f}" stroke="{DIM_COLOR}" stroke-width="0.8"/>')
    for _sx, _sy in [(_dsx1, _dsy1), (_dsx2, _dsy2)]:
        out.append(f'<line x1="{_sx - _tk * _px:.1f}" y1="{_sy - _tk * _py:.1f}" '
                   f'x2="{_sx + _tk * _px:.1f}" y2="{_sy + _tk * _py:.1f}" '
                   f'stroke="{DIM_COLOR}" stroke-width="0.8"/>')
    _lmx = (_dsx1 + _dsx2) / 2; _lmy = (_dsy1 + _dsy2) / 2
    _up_dx = _dsx2 - _dsx1; _up_dy = _dsy2 - _dsy1
    _up_ang = math.degrees(math.atan2(_up_dy, _up_dx))
    _lx = _lmx - 3 * _px; _ly = _lmy - 3 * _py
    out.append(f'<text x="{_lx:.1f}" y="{_ly:.1f}" text-anchor="middle" font-family="Arial" '
               f'font-size="8" fill="{DIM_COLOR}" transform="rotate({_up_ang:.1f},{_lx:.1f},{_ly:.1f})">'
               f'{fmt_dist(_dim_len)}</text>')

    # O10 inner center perpendicular to W20-W0 up to IW1 south face
    _o10_open = compute_outer_openings(pts, layout)[9]  # O10
    _o10_ic = ((_o10_open.poly[2][0] + _o10_open.poly[3][0]) / 2,
               (_o10_open.poly[2][1] + _o10_open.poly[3][1]) / 2)
    _t_iw1_10 = (layout.iw1_s - _o10_ic[1]) / _nrmN
    _dim_end10 = (_o10_ic[0] + _t_iw1_10 * _nrmE, layout.iw1_s)
    _dim_len10 = abs(_t_iw1_10)
    _dsx1, _dsy1 = to_svg(*_o10_ic)
    _dsx2, _dsy2 = to_svg(*_dim_end10)
    _sdx = _dsx2 - _dsx1; _sdy = _dsy2 - _dsy1
    _slen = math.sqrt(_sdx**2 + _sdy**2)
    _px = -_sdy / _slen; _py = _sdx / _slen
    out.append(f'<line x1="{_dsx1:.1f}" y1="{_dsy1:.1f}" x2="{_dsx2:.1f}" y2="{_dsy2:.1f}" stroke="{DIM_COLOR}" stroke-width="0.8"/>')
    for _sx, _sy in [(_dsx1, _dsy1), (_dsx2, _dsy2)]:
        out.append(f'<line x1="{_sx - _tk * _px:.1f}" y1="{_sy - _tk * _py:.1f}" '
                   f'x2="{_sx + _tk * _px:.1f}" y2="{_sy + _tk * _py:.1f}" '
                   f'stroke="{DIM_COLOR}" stroke-width="0.8"/>')
    _lmx = (_dsx1 + _dsx2) / 2; _lmy = (_dsy1 + _dsy2) / 2
    _up_dx = _dsx2 - _dsx1; _up_dy = _dsy2 - _dsy1
    _up_ang = math.degrees(math.atan2(_up_dy, _up_dx))
    _lx = _lmx - 3 * _px; _ly = _lmy - 3 * _py
    out.append(f'<text x="{_lx:.1f}" y="{_ly:.1f}" text-anchor="middle" font-family="Arial" '
               f'font-size="8" fill="{DIM_COLOR}" transform="rotate({_up_ang:.1f},{_lx:.1f},{_ly:.1f})">'
               f'{fmt_dist(_dim_len10)}</text>')

    # IW9 east face to IW11 west face, perpendicular to IW9 east face
    _iw9_se9 = layout.iw9_poly[1]
    _iw9_ne9 = layout.iw9_poly[2]
    _dim3_s = ((_iw9_se9[0] + _iw9_ne9[0]) / 2,
               (_iw9_se9[1] + _iw9_ne9[1]) / 2)
    _fd3 = (_iw9_ne9[0] - _iw9_se9[0], _iw9_ne9[1] - _iw9_se9[1])
    _fd3_len = math.sqrt(_fd3[0]**2 + _fd3[1]**2)
    _perp3 = (_fd3[1] / _fd3_len, -_fd3[0] / _fd3_len)  # CW perp, toward IW11
    _iw11_sw3 = layout.iw11_poly[0]
    _iw11_nw3 = layout.iw11_poly[3]
    _dw11f = (_iw11_nw3[0] - _iw11_sw3[0], _iw11_nw3[1] - _iw11_sw3[1])
    _dx3 = _iw11_sw3[0] - _dim3_s[0]
    _dy3 = _iw11_sw3[1] - _dim3_s[1]
    _det3 = _dw11f[0] * _perp3[1] - _dw11f[1] * _perp3[0]
    _t3 = (_dw11f[0] * _dy3 - _dw11f[1] * _dx3) / _det3
    _dim3_e = (_dim3_s[0] + _t3 * _perp3[0], _dim3_s[1] + _t3 * _perp3[1])
    _dsx1, _dsy1 = to_svg(*_dim3_s)
    _dsx2, _dsy2 = to_svg(*_dim3_e)
    _sdx = _dsx2 - _dsx1; _sdy = _dsy2 - _dsy1
    _slen = math.sqrt(_sdx**2 + _sdy**2)
    _px = -_sdy / _slen; _py = _sdx / _slen
    _tk = 4
    out.append(f'<line x1="{_dsx1:.1f}" y1="{_dsy1:.1f}" x2="{_dsx2:.1f}" y2="{_dsy2:.1f}" stroke="{DIM_COLOR}" stroke-width="0.8"/>')
    for _sx, _sy in [(_dsx1, _dsy1), (_dsx2, _dsy2)]:
        out.append(f'<line x1="{_sx - _tk * _px:.1f}" y1="{_sy - _tk * _py:.1f}" '
                   f'x2="{_sx + _tk * _px:.1f}" y2="{_sy + _tk * _py:.1f}" '
                   f'stroke="{DIM_COLOR}" stroke-width="0.8"/>')
    _lmx = (_dsx1 + _dsx2) / 2; _lmy = (_dsy1 + _dsy2) / 2
    _up_dx = _dsx2 - _dsx1; _up_dy = _dsy2 - _dsy1
    _up_ang = math.degrees(math.atan2(_up_dy, _up_dx))
    _lx = _lmx - 3 * _px; _ly = _lmy - 3 * _py
    out.append(f'<text x="{_lx:.1f}" y="{_ly:.1f}" text-anchor="middle" font-family="Arial" '
               f'font-size="8" fill="{DIM_COLOR}" transform="rotate({_up_ang:.1f},{_lx:.1f},{_ly:.1f})">'
               f'{fmt_dist(_t3)}</text>')

    # Utility area N-S: IW8 south face to W20-W0 at O11 center
    _o11 = compute_outer_openings(pts, layout)[10]  # O11
    _o11_inner_mid = ((_o11.poly[2][0] + _o11.poly[3][0]) / 2,
                      (_o11.poly[2][1] + _o11.poly[3][1]) / 2)
    dim_line_v(out, _o11_inner_mid[0], _o11_inner_mid[1], layout.iw8.s,
               fmt_dist(layout.iw8.s - _o11_inner_mid[1]), to_svg)


def _render_openings(out, data, layout):
    """Render door swings and jamb blocks for O3, O6, RO1-RO5.

    Opening fill polygons are rendered by _render_walls() as part of the
    double-shell wall section loop.
    """
    pts = data.pts
    to_svg = data.to_svg
    outer_openings = compute_outer_openings(pts, layout)

    # O3 door: 30" door, hinged north, swings east
    o3 = [o for o in outer_openings if o.name == "O3"][0]
    # O3 poly: [(F4_e, south), (F4_e, north), (W4_e, north), (W4_e, south)]
    wall_w = pts["F4"][0]  # outer face
    wall_e = pts["W4"][0]  # inner face
    wall_mid = (wall_w + wall_e) / 2
    o3_s = o3.poly[0][1]
    o3_n = o3.poly[1][1]
    o3_opening = o3_n - o3_s
    gap = (o3_opening - O3_DOOR_WIDTH) / 2
    door_s = o3_s + gap
    door_n = o3_n - gap

    # Jamb blocks
    block_w = wall_mid - DOOR_FLAT_FACE / 2
    block_e = wall_mid + DOOR_FLAT_FACE / 2
    # South block
    bx1, by1 = to_svg(block_w, o3_s + gap)
    bx2, by2 = to_svg(block_e, o3_s)
    out.append(f'<rect x="{bx1:.1f}" y="{by1:.1f}" width="{bx2-bx1:.1f}" height="{by2-by1:.1f}"'
               f' fill="{JAMB_COLOR}" stroke="none"/>')
    # North block
    bx1, by1 = to_svg(block_w, o3_n)
    bx2, by2 = to_svg(block_e, o3_n - gap)
    out.append(f'<rect x="{bx1:.1f}" y="{by1:.1f}" width="{bx2-bx1:.1f}" height="{by2-by1:.1f}"'
               f' fill="{JAMB_COLOR}" stroke="none"/>')

    # Door: hinge at north side, swings east
    hinge_e, hinge_n = wall_mid, door_n
    hx, hy = to_svg(hinge_e, hinge_n)
    # Straight line from hinge eastward (door in open position)
    tip_e, tip_n = hinge_e + O3_DOOR_WIDTH, hinge_n
    tx, ty = to_svg(tip_e, tip_n)
    out.append(f'<line x1="{hx:.1f}" y1="{hy:.1f}" x2="{tx:.1f}" y2="{ty:.1f}"'
               f' stroke="{JAMB_COLOR}" stroke-width="1.0"/>')
    # Arc from open (east) sweeping CW 90° to closed (south)
    n_arc = 20
    arc_pts = []
    for i in range(n_arc + 1):
        angle = -i * (math.pi / 2) / n_arc  # 0° to -90°
        ae = hinge_e + O3_DOOR_WIDTH * math.cos(angle)
        an = hinge_n + O3_DOOR_WIDTH * math.sin(angle)
        sx, sy = to_svg(ae, an)
        arc_pts.append(f"{sx:.1f},{sy:.1f}")
    out.append(f'<polyline points="{" ".join(arc_pts)}" fill="none"'
               f' stroke="{JAMB_COLOR}" stroke-width="0.5"/>')

    # O6 door: 42" door, hinged east, swings south
    o6 = [o for o in outer_openings if o.name == "O6"][0]
    o6_w = o6.poly[0][0]
    o6_e = o6.poly[1][0]
    wall_n = pts["F9"][1]  # outer face
    wall_s = pts["W9"][1]  # inner face
    wall_mid = (wall_n + wall_s) / 2
    gap = (O6_WIDTH - O6_DOOR_WIDTH) / 2
    door_w = o6_w + gap
    door_e = o6_e - gap

    # Jamb blocks: fill gap between opening edge and door edge
    # Block thickness = flat face of U-turn
    block_s = wall_mid - DOOR_FLAT_FACE / 2
    block_n = wall_mid + DOOR_FLAT_FACE / 2
    # West block
    bx1, by1 = to_svg(o6_w, block_n)
    bx2, by2 = to_svg(door_w, block_s)
    out.append(f'<rect x="{bx1:.1f}" y="{by1:.1f}" width="{bx2-bx1:.1f}" height="{by2-by1:.1f}"'
               f' fill="{JAMB_COLOR}" stroke="none"/>')
    # East block
    bx1, by1 = to_svg(door_e, block_n)
    bx2, by2 = to_svg(o6_e, block_s)
    out.append(f'<rect x="{bx1:.1f}" y="{by1:.1f}" width="{bx2-bx1:.1f}" height="{by2-by1:.1f}"'
               f' fill="{JAMB_COLOR}" stroke="none"/>')

    # Door: hinge at east side, swings south
    hinge_e, hinge_n = door_e, wall_mid
    hx, hy = to_svg(hinge_e, hinge_n)
    # Straight line from hinge southward (door in open position)
    tip_e, tip_n = hinge_e, hinge_n - O6_DOOR_WIDTH
    tx, ty = to_svg(tip_e, tip_n)
    out.append(f'<line x1="{hx:.1f}" y1="{hy:.1f}" x2="{tx:.1f}" y2="{ty:.1f}"'
               f' stroke="{JAMB_COLOR}" stroke-width="1.0"/>')
    # Arc from open (south) sweeping CW 90° to closed (west)
    # In survey coords: center=hinge, from (hinge_e, hinge_n - door_w) to (hinge_e - door_w, hinge_n)
    # In SVG coords: y is flipped, so CW in survey becomes CCW in SVG
    n_arc = 20
    arc_pts = []
    for i in range(n_arc + 1):
        angle = -math.pi / 2 - i * (math.pi / 2) / n_arc  # -90° to -180° in survey
        ae = hinge_e + O6_DOOR_WIDTH * math.cos(angle)
        an = hinge_n + O6_DOOR_WIDTH * math.sin(angle)
        sx, sy = to_svg(ae, an)
        arc_pts.append(f"{sx:.1f},{sy:.1f}")
    out.append(f'<polyline points="{" ".join(arc_pts)}" fill="none"'
               f' stroke="{JAMB_COLOR}" stroke-width="0.5"/>')

    # RO1 door: 36" door, hinged east, swings south
    rough_openings = compute_rough_openings(pts, layout)
    ro1 = [r for r in rough_openings if r.name == "RO1"][0]
    ro1_mid = (ro1.bbox.s + ro1.bbox.n) / 2
    ro1_gap = (ro1.bbox.e - ro1.bbox.w - RO1_DOOR_WIDTH) / 2
    hinge_e, hinge_n = ro1.bbox.e - ro1_gap, ro1_mid
    hx, hy = to_svg(hinge_e, hinge_n)
    # Straight line from hinge southward (door in open position)
    tip_e, tip_n = hinge_e, hinge_n - RO1_DOOR_WIDTH
    tx, ty = to_svg(tip_e, tip_n)
    out.append(f'<line x1="{hx:.1f}" y1="{hy:.1f}" x2="{tx:.1f}" y2="{ty:.1f}"'
               f' stroke="{JAMB_COLOR}" stroke-width="1.0"/>')
    # Arc from open (south) sweeping CW 90° to closed (west)
    n_arc = 20
    arc_pts = []
    for i in range(n_arc + 1):
        angle = -math.pi / 2 - i * (math.pi / 2) / n_arc
        ae = hinge_e + RO1_DOOR_WIDTH * math.cos(angle)
        an = hinge_n + RO1_DOOR_WIDTH * math.sin(angle)
        sx, sy = to_svg(ae, an)
        arc_pts.append(f"{sx:.1f},{sy:.1f}")
    out.append(f'<polyline points="{" ".join(arc_pts)}" fill="none"'
               f' stroke="{JAMB_COLOR}" stroke-width="0.5"/>')

    # RO2 door: 36" door in IW11 (rotated), hinged at NNE edge, swings into office
    ro2 = [r for r in rough_openings if r.name == "RO2"][0]
    _ro2p = ro2.poly  # [SW, SE, NE, NW]
    # IW11 unit vectors (recompute from polygon)
    _i11_sw, _i11_se, _i11_ne, _i11_nw = layout.iw11_poly
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
    _closed_ang = math.atan2(-_i11_an[1], -_i11_an[0])
    n_arc = 20
    arc_pts = []
    for i in range(n_arc + 1):
        angle = _closed_ang + i * (math.pi / 2) / n_arc
        ae = hinge_e + RO2_DOOR_WIDTH * math.cos(angle)
        an = hinge_n + RO2_DOOR_WIDTH * math.sin(angle)
        sx, sy = to_svg(ae, an)
        arc_pts.append(f"{sx:.1f},{sy:.1f}")
    out.append(f'<polyline points="{" ".join(arc_pts)}" fill="none"'
               f' stroke="{JAMB_COLOR}" stroke-width="0.5"/>')

    # RO4 door: 36" door, hinged north, swings west
    ro4 = [r for r in rough_openings if r.name == "RO4"][0]
    ro4_mid = (ro4.bbox.w + ro4.bbox.e) / 2
    ro4_gap = (ro4.bbox.n - ro4.bbox.s - RO4_DOOR_WIDTH) / 2
    hinge_e, hinge_n = ro4_mid, ro4.bbox.n - ro4_gap
    hx, hy = to_svg(hinge_e, hinge_n)
    # Straight line from hinge westward (door in open position)
    tip_e, tip_n = hinge_e - RO4_DOOR_WIDTH, hinge_n
    tx, ty = to_svg(tip_e, tip_n)
    out.append(f'<line x1="{hx:.1f}" y1="{hy:.1f}" x2="{tx:.1f}" y2="{ty:.1f}"'
               f' stroke="{JAMB_COLOR}" stroke-width="1.0"/>')
    # Arc from open (west) sweeping 90° to closed (south)
    n_arc = 20
    arc_pts = []
    for i in range(n_arc + 1):
        angle = math.pi + i * (math.pi / 2) / n_arc  # 180° to 270°
        ae = hinge_e + RO4_DOOR_WIDTH * math.cos(angle)
        an = hinge_n + RO4_DOOR_WIDTH * math.sin(angle)
        sx, sy = to_svg(ae, an)
        arc_pts.append(f"{sx:.1f},{sy:.1f}")
    out.append(f'<polyline points="{" ".join(arc_pts)}" fill="none"'
               f' stroke="{JAMB_COLOR}" stroke-width="0.5"/>')

    # RO5 door: 36" door, hinged east, swings north
    ro5 = [r for r in rough_openings if r.name == "RO5"][0]
    ro5_mid = (ro5.bbox.s + ro5.bbox.n) / 2
    ro5_gap = (ro5.bbox.e - ro5.bbox.w - RO5_DOOR_WIDTH) / 2
    hinge_e, hinge_n = ro5.bbox.e - ro5_gap, ro5_mid
    hx, hy = to_svg(hinge_e, hinge_n)
    # Straight line from hinge northward (door in open position)
    tip_e, tip_n = hinge_e, hinge_n + RO5_DOOR_WIDTH
    tx, ty = to_svg(tip_e, tip_n)
    out.append(f'<line x1="{hx:.1f}" y1="{hy:.1f}" x2="{tx:.1f}" y2="{ty:.1f}"'
               f' stroke="{JAMB_COLOR}" stroke-width="1.0"/>')
    # Arc from open (north) sweeping 90° to closed (west)
    n_arc = 20
    arc_pts = []
    for i in range(n_arc + 1):
        angle = math.pi / 2 + i * (math.pi / 2) / n_arc  # 90° to 180°
        ae = hinge_e + RO5_DOOR_WIDTH * math.cos(angle)
        an = hinge_n + RO5_DOOR_WIDTH * math.sin(angle)
        sx, sy = to_svg(ae, an)
        arc_pts.append(f"{sx:.1f},{sy:.1f}")
    out.append(f'<polyline points="{" ".join(arc_pts)}" fill="none"'
               f' stroke="{JAMB_COLOR}" stroke-width="0.5"/>')

    # RO3 door: 36" door, hinged north, swings east
    ro3 = [r for r in rough_openings if r.name == "RO3"][0]
    ro3_mid = (ro3.bbox.w + ro3.bbox.e) / 2
    ro3_gap = (ro3.bbox.n - ro3.bbox.s - RO3_DOOR_WIDTH) / 2
    hinge_e, hinge_n = ro3_mid, ro3.bbox.n - ro3_gap
    hx, hy = to_svg(hinge_e, hinge_n)
    # Straight line from hinge eastward (door in open position)
    tip_e, tip_n = hinge_e + RO3_DOOR_WIDTH, hinge_n
    tx, ty = to_svg(tip_e, tip_n)
    out.append(f'<line x1="{hx:.1f}" y1="{hy:.1f}" x2="{tx:.1f}" y2="{ty:.1f}"'
               f' stroke="{JAMB_COLOR}" stroke-width="1.0"/>')
    # Arc from open (east) sweeping 90° to closed (south)
    n_arc = 20
    arc_pts = []
    for i in range(n_arc + 1):
        angle = 0 + i * (-math.pi / 2) / n_arc  # 0° to -90° (east to south)
        ae = hinge_e + RO3_DOOR_WIDTH * math.cos(angle)
        an = hinge_n + RO3_DOOR_WIDTH * math.sin(angle)
        sx, sy = to_svg(ae, an)
        arc_pts.append(f"{sx:.1f},{sy:.1f}")
    out.append(f'<polyline points="{" ".join(arc_pts)}" fill="none"'
               f' stroke="{JAMB_COLOR}" stroke-width="0.5"/>')


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

def render_floorplan_svg(data, room_title="Parent Suite", minik=False, db=False):
    """Render the complete floorplan SVG. Returns SVG string."""
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

    _render_walls(out, data, layout)
    _render_appliances(out, data, layout, minik=minik, db=db)
    _render_kitchen(out, data, layout, minik=minik, db=db)
    _render_furniture(out, data, layout, minik=minik, db=db)
    out.append('<g opacity="0.5">')
    _render_dimensions(out, data, layout)
    out.append('</g>')
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

    db_content = render_floorplan_svg(data, db=True)
    db_path = os.path.join(base_dir, "floorplan_db.svg")
    with open(db_path, "w") as f:
        f.write(db_content)
    print(f"Floorplan (daybed) written to {db_path}")

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
