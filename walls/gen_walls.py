"""Generate outer wall detail SVG showing double-shell concrete construction.

Outer walls are double-shell 3D-printed concrete: two 2" shells separated
by a 4" air gap (8" total = WALL_OUTER). At openings, shells connect via
90-degree corner turns with configurable radii.

Outputs walls/walls.svg at 1:72 scale.
"""
import os, sys, math, datetime
from typing import NamedTuple

# Ensure project root is on sys.path for package imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from shared.types import LineSeg, ArcSeg
from shared.geometry import fmt_dist, left_norm
from shared.svg import make_svg_transform, W, H, git_describe
from floorplan.gen_floorplan import build_floorplan_data
from floorplan.constants import WALL_OUTER
from floorplan.openings import compute_rough_openings
from walls.constants import SHELL_THICKNESS, AIR_GAP, OPENING_INSIDE_RADIUS
from shared.wall_shells import (
    lerp, openings_on_seg, solid_ranges,
    arc_strip_poly, line_strip_poly, partial_line_strip, partial_line_strip_2,
    uturn_arc_data, uturn_polygon,
    trace_boundary_path, enumerate_wall_sections, build_section_outlines,
)




# ============================================================
# Length computation
# ============================================================

def _seg_arc_sweep(seg, pts):
    """Compute sweep angle (radians) of an arc segment."""
    c = pts[seg.center]
    s = pts[seg.start]
    e = pts[seg.end]
    ang_s = math.atan2(s[1] - c[1], s[0] - c[0])
    ang_e = math.atan2(e[1] - c[1], e[0] - c[0])
    if seg.direction == "CW":
        return (ang_s - ang_e) % (2 * math.pi)
    else:
        return (ang_e - ang_s) % (2 * math.pi)


def _path_length_between(pts, outline_segs, start_seg_idx, start_t,
                         end_seg_idx, end_t, inset):
    """Compute path length between two parametric positions at a given inset.

    For line segments, inset does not change length (perpendicular offset).
    For arc segments, length = (R ± inset) * sweep_angle.
      CW arcs (convex): R_adj = R - inset
      CCW arcs (concave): R_adj = R + inset
    Returns length in feet.
    """
    n_segs = len(outline_segs)
    total = 0.0

    if start_seg_idx == end_seg_idx:
        # Same segment (must be a LineSeg — openings only on lines)
        seg = outline_segs[start_seg_idx]
        A, B = pts[seg.start], pts[seg.end]
        seg_len = math.sqrt((B[0] - A[0])**2 + (B[1] - A[1])**2)
        return seg_len * (end_t - start_t)

    # Start segment: partial from start_t to 1.0
    seg = outline_segs[start_seg_idx]
    A, B = pts[seg.start], pts[seg.end]
    seg_len = math.sqrt((B[0] - A[0])**2 + (B[1] - A[1])**2)
    total += seg_len * (1.0 - start_t)

    # Intermediate full segments
    idx = (start_seg_idx + 1) % n_segs
    while idx != end_seg_idx:
        seg = outline_segs[idx]
        if isinstance(seg, ArcSeg):
            if idx == 8 and inset >= SHELL_THICKNESS + AIR_GAP:
                # F8-F9 inner shell: straight-arc-straight path length
                R_a8 = seg.radius
                R_turn = OPENING_INSIDE_RADIUS + (inset - (SHELL_THICKNESS + AIR_GAP))
                d_straight = R_a8 + inset - R_turn
                total += 2 * d_straight + R_turn * (math.pi / 2)
            else:
                sweep = _seg_arc_sweep(seg, pts)
                R = seg.radius
                R_adj = (R - inset) if seg.direction == "CW" else (R + inset)
                total += R_adj * sweep
        else:
            A, B = pts[seg.start], pts[seg.end]
            total += math.sqrt((B[0] - A[0])**2 + (B[1] - A[1])**2)
        idx = (idx + 1) % n_segs

    # End segment: partial from 0.0 to end_t
    seg = outline_segs[end_seg_idx]
    A, B = pts[seg.start], pts[seg.end]
    seg_len = math.sqrt((B[0] - A[0])**2 + (B[1] - A[1])**2)
    total += seg_len * end_t

    return total


# ============================================================
# Data computation
# ============================================================

class WallData(NamedTuple):
    """Typed container for all wall detail geometry and page layout."""
    pts: dict
    to_svg: object  # Callable[[float, float], tuple[float, float]]
    outline_segs: list
    inner_segs: list
    s_segs: list
    g_segs: list
    radii: dict
    openings: list
    layout: object  # InteriorLayout
    inner_poly: list
    outer_area: float
    inner_area: float
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
    ft_per_inch: float
    g_f8f9_poly: list      # G-series F8-F9 straight-arc-straight polyline
    w_f8f9_poly: list      # W-series F8-F9 straight-arc-straight polyline


def build_wall_data():
    """Compute all geometry needed for the wall detail SVG."""
    fp_data = build_floorplan_data()
    pts = fp_data.pts
    to_svg = fp_data.to_svg
    outline_segs = fp_data.outline_segs
    inner_segs = fp_data.inner_segs
    radii = fp_data.radii
    inner_poly = fp_data.inner_poly

    # S/G series, openings, layout, F8-F9 polylines from FloorplanData
    s_segs = fp_data.s_segs
    g_segs = fp_data.g_segs
    layout = fp_data.layout
    openings = fp_data.openings
    g_f8f9_poly = fp_data.g_f8f9_poly
    w_f8f9_poly = fp_data.w_f8f9_poly

    # --- Page layout: 1:72 scale ---
    _f_svg = [to_svg(*pts[f"F{i}"]) for i in range(22)]
    _bldg_xmin = min(p[0] for p in _f_svg)
    _bldg_xmax = max(p[0] for p in _f_svg)
    _bldg_ymin = min(p[1] for p in _f_svg)
    _bldg_ymax = max(p[1] for p in _f_svg)
    _bldg_cx = (_bldg_xmin + _bldg_xmax) / 2

    _title_y = _bldg_ymin - 49

    _tb_w, _tb_h = 130, 80
    _tb_left = _bldg_xmax + 10
    _tb_right = _tb_left + _tb_w
    _tb_top = _title_y - 14 * 0.35
    _tb_bottom = _tb_top + _tb_h
    _tb_cx = (_tb_left + _tb_right) / 2

    _scale_ratio = 72
    _ft_per_inch = _scale_ratio / 12.0
    _svg_per_ft = to_svg(1, 0)[0] - to_svg(0, 0)[0]
    _fit_scale = 72.0 / (_ft_per_inch * _svg_per_ft)

    _margin_top = 36
    _vb_w = W / _fit_scale
    _vb_h = H / _fit_scale
    _cb_xmin = _bldg_xmin - 25
    _cb_xmax = _tb_right + 5
    _cb_cx = (_cb_xmin + _cb_xmax) / 2
    _cb_ymin = _title_y - 14 - 5
    _vb_x = _cb_cx - _vb_w / 2
    _vb_y = _cb_ymin - _margin_top / _fit_scale

    return WallData(
        pts=pts, to_svg=to_svg,
        outline_segs=outline_segs, inner_segs=inner_segs,
        s_segs=s_segs, g_segs=g_segs,
        radii=radii, openings=openings,
        layout=layout, inner_poly=inner_poly,
        outer_area=fp_data.outer_area,
        inner_area=fp_data.inner_area,
        vb_x=_vb_x, vb_y=_vb_y, vb_w=_vb_w, vb_h=_vb_h,
        title_x=_bldg_cx, title_y=_title_y,
        tb_left=_tb_left, tb_right=_tb_right, tb_top=_tb_top,
        tb_bottom=_tb_bottom, tb_w=_tb_w, tb_h=_tb_h, tb_cx=_tb_cx,
        ft_per_inch=_ft_per_inch,
        g_f8f9_poly=g_f8f9_poly,
        w_f8f9_poly=w_f8f9_poly,
    )


# ============================================================
# SVG rendering
# ============================================================

def _svg_polygon(out, poly, to_svg, fill, stroke="#666", stroke_width="0.5"):
    """Render a polygon as an SVG element."""
    svg = " ".join(f"{to_svg(*p)[0]:.2f},{to_svg(*p)[1]:.2f}" for p in poly)
    out.append(f'<polygon points="{svg}" fill="{fill}" '
               f'stroke="{stroke}" stroke-width="{stroke_width}"/>')


def _render_interior_walls(out, data):
    """Render interior wall polygons and labels into the SVG output list."""
    pts = data.pts
    to_svg = data.to_svg
    layout = data.layout
    inner_poly = data.inner_poly

    IW_FILL = "rgba(160,160,160,0.35)"
    IW_STROKE = "#666"
    IW_SW = "0.5"
    LABEL_SIZE = "6"
    LABEL_GAP = 3.0  # SVG px from wall face to label center

    def iw_poly(poly):
        svg = " ".join(f"{to_svg(*p)[0]:.2f},{to_svg(*p)[1]:.2f}" for p in poly)
        out.append(f'<polygon points="{svg}" fill="{IW_FILL}" '
                   f'stroke="{IW_STROKE}" stroke-width="{IW_SW}"/>')

    def iw_rect(w, e, s, n):
        iw_poly([(w, s), (e, s), (e, n), (w, n)])

    def iw_label(name, w, e, s, n, vertical=True, n_shift=0.0):
        """Label just outside the wall: north for horizontal, west for vertical.

        n_shift: survey-feet offset applied to the label's northing
        (negative = south).
        """
        if vertical:
            # N-S wall: label west of wall, centered vertically
            lx, ly = to_svg(w, (s + n) / 2 + n_shift)
            lx -= LABEL_GAP
            rot = f' transform="rotate(-90 {lx:.1f} {ly:.1f})"'
        else:
            # W-E wall: label north of wall, centered horizontally
            lx, ly = to_svg((w + e) / 2, n + n_shift)
            ly -= LABEL_GAP
            rot = ""
        out.append(f'<text x="{lx:.1f}" y="{ly:.1f}" text-anchor="middle"'
                   f' dominant-baseline="central" font-family="Arial"'
                   f' font-size="{LABEL_SIZE}" fill="#666"{rot}>{name}</text>')

    rough_openings = compute_rough_openings(pts, layout)
    ro_map = {ro.name: ro.bbox for ro in rough_openings}

    # IW1 (horizontal, 6")
    iw_poly(layout.iw1)
    iw_label("IW1", layout.iw1[0][0], layout.iw1[1][0],
             layout.iw1_s, layout.iw1_n, vertical=False)

    # IW2 (vertical, 6")
    iw_rect(layout.iw2.w, layout.iw2.e, layout.iw2.s, layout.iw2.n)
    iw_label("IW2", layout.iw2.w, layout.iw2.e, layout.iw2.s, layout.iw2.n)

    # IW6 (horizontal, 1" partition)
    iw_poly(layout.iw6_poly)
    iw_label("IW6", min(layout.iw6_poly[0][0], layout.iw6_poly[3][0]),
             layout.iw2.w, layout.iw6_s, layout.iw6_n, vertical=False)

    # IW7 (straight N-S, 3")
    iw_poly(layout.iw7)
    iw_label("IW7", layout.iw7[0][0], layout.iw7[1][0],
             layout.iw7[0][1], layout.iw7[3][1])

    # IW3 (vertical, 4") — west face aligned with IW7
    iw_rect(layout.iw3.w, layout.iw3.e, layout.iw3.s, layout.iw3.n)
    # Label centered between RO3 north end and IW3 north end
    ro3_n = ro_map["RO3"].n
    iw3_label_shift = (ro3_n - layout.iw3.s) / 2
    iw_label("IW3", layout.iw3.w, layout.iw3.e, layout.iw3.s, layout.iw3.n,
             n_shift=iw3_label_shift)

    # IW9 (vertical, 4") — old IW3 position, south of IW7 L north face
    iw_rect(layout.iw9.w, layout.iw9.e, layout.iw9.s, layout.iw9.n)
    iw_label("IW9", layout.iw9.w, layout.iw9.e, layout.iw9.s, layout.iw9.n)

    # IW10 (horizontal, 4") — closet north wall
    iw_rect(layout.iw10.w, layout.iw10.e, layout.iw10.s, layout.iw10.n)
    iw_label("IW10", layout.iw10.w, layout.iw10.e, layout.iw10.s, layout.iw10.n,
             vertical=False)

    # IW4 (vertical, 4")
    iw_rect(layout.iw4_w, layout.iw4_e, layout.iw4_s, layout.iw3.n)
    iw_label("IW4", layout.iw4_w, layout.iw4_e, layout.iw4_s, layout.iw3.n)

    # IW8 (N-S, 3")
    iw_poly(layout.iw8)
    iw_label("IW8", layout.iw8[3][0], layout.iw8[2][0],
             layout.iw8[2][1], layout.iw8[1][1])

    # IW11 (vertical, 4")
    iw_rect(layout.iw11.w, layout.iw11.e, layout.iw11.s, layout.iw11.n)
    iw_label("IW11", layout.iw11.w, layout.iw11.e, layout.iw11.s, layout.iw11.n)

    # IW12 (horizontal, 4")
    iw_rect(layout.iw12.w, layout.iw12.e, layout.iw12.s, layout.iw12.n)
    iw_label("IW12", layout.iw12.w, layout.iw12.e, layout.iw12.s, layout.iw12.n,
             vertical=False)

    # IW5 (horizontal, 3")
    iw_rect(layout.iw5.w, layout.iw5.e, layout.iw5.s, layout.iw5.n)
    iw_label("IW5", layout.iw5.w, layout.iw5.e, layout.iw5.s, layout.iw5.n,
             vertical=False)

    # Rough openings (RO1-RO6) — dark red outline box with X
    _RO_COLOR = "darkred"
    _RO_SW = "0.5"
    for ro in rough_openings:
        b = ro.bbox
        x1, y1 = to_svg(b.w, b.n)  # NW corner (SVG top-left)
        x2, y2 = to_svg(b.e, b.s)  # SE corner (SVG bottom-right)
        out.append(f'<rect x="{x1:.1f}" y="{y1:.1f}" width="{x2 - x1:.1f}" height="{y2 - y1:.1f}"'
                   f' fill="none" stroke="{_RO_COLOR}" stroke-width="{_RO_SW}"/>')
        out.append(f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}"'
                   f' stroke="{_RO_COLOR}" stroke-width="{_RO_SW}"/>')
        out.append(f'<line x1="{x2:.1f}" y1="{y1:.1f}" x2="{x1:.1f}" y2="{y2:.1f}"'
                   f' stroke="{_RO_COLOR}" stroke-width="{_RO_SW}"/>')

    # RO labels — same placement convention as IW labels
    for ro in rough_openings:
        b = ro.bbox
        if ro.orientation == "H":
            # Horizontal opening: label centered above (north)
            lx, ly = to_svg((b.w + b.e) / 2, b.n)
            ly -= LABEL_GAP
            rot = ""
        else:
            # Vertical opening: label centered left (west)
            lx, ly = to_svg(b.w, (b.s + b.n) / 2)
            lx -= LABEL_GAP
            rot = f' transform="rotate(-90 {lx:.1f} {ly:.1f})"'
        out.append(f'<text x="{lx:.1f}" y="{ly:.1f}" text-anchor="middle"'
                   f' dominant-baseline="central" font-family="Arial"'
                   f' font-size="{LABEL_SIZE}" fill="{_RO_COLOR}"{rot}>{ro.name}</text>')


def _render_opening_dims(out, data):
    """Render opening width dimension lines on the exterior face."""
    pts = data.pts
    to_svg = data.to_svg
    outline_segs = data.outline_segs
    openings = data.openings

    DIM_OFFSET = 0.75       # feet outside F-face
    TICK = 3.0              # SVG pts, tick half-length
    EXT_OVERSHOOT = 1.5     # SVG pts, extension line past dim line
    LABEL_OFFSET = 4.0      # SVG pts, label offset from dim line toward exterior
    DIM_COLOR = "#4682B4"
    DIM_SW = "0.4"
    FONT_SIZE = "5"

    for op in openings:
        seg = outline_segs[op.seg_idx]
        F_A, F_B = pts[seg.start], pts[seg.end]

        # Boundary points on F-face
        p1 = lerp(F_A, F_B, op.t_start)
        p2 = lerp(F_A, F_B, op.t_end)

        # Exterior normal (left of CW traversal)
        n = left_norm(F_A, F_B)

        # Dim line position (offset outward from F-face)
        q1 = (p1[0] + n[0] * DIM_OFFSET, p1[1] + n[1] * DIM_OFFSET)
        q2 = (p2[0] + n[0] * DIM_OFFSET, p2[1] + n[1] * DIM_OFFSET)

        # SVG coords
        sx1, sy1 = to_svg(*q1)
        sx2, sy2 = to_svg(*q2)
        fx1, fy1 = to_svg(*p1)
        fx2, fy2 = to_svg(*p2)

        # Dim line direction in SVG (unit vector)
        ddx, ddy = sx2 - sx1, sy2 - sy1
        dim_len = math.sqrt(ddx**2 + ddy**2)
        if dim_len < 1e-6:
            continue
        udx, udy = ddx / dim_len, ddy / dim_len

        # Tick direction (perpendicular to dim line in SVG)
        tx, ty = -udy, udx

        # Extension line normal in SVG (from F-face toward dim line)
        enx, eny = sx1 - fx1, sy1 - fy1
        en_len = math.sqrt(enx**2 + eny**2)
        if en_len > 1e-6:
            enx, eny = enx / en_len, eny / en_len

        # Extension lines (from F-face to dim line + overshoot)
        out.append(f'<line x1="{fx1:.2f}" y1="{fy1:.2f}"'
                   f' x2="{sx1 + enx * EXT_OVERSHOOT:.2f}"'
                   f' y2="{sy1 + eny * EXT_OVERSHOOT:.2f}"'
                   f' stroke="{DIM_COLOR}" stroke-width="{DIM_SW}"/>')
        out.append(f'<line x1="{fx2:.2f}" y1="{fy2:.2f}"'
                   f' x2="{sx2 + enx * EXT_OVERSHOOT:.2f}"'
                   f' y2="{sy2 + eny * EXT_OVERSHOOT:.2f}"'
                   f' stroke="{DIM_COLOR}" stroke-width="{DIM_SW}"/>')

        # Tick marks at dim line endpoints
        out.append(f'<line x1="{sx1 - tx * TICK:.2f}" y1="{sy1 - ty * TICK:.2f}"'
                   f' x2="{sx1 + tx * TICK:.2f}" y2="{sy1 + ty * TICK:.2f}"'
                   f' stroke="{DIM_COLOR}" stroke-width="{DIM_SW}"/>')
        out.append(f'<line x1="{sx2 - tx * TICK:.2f}" y1="{sy2 - ty * TICK:.2f}"'
                   f' x2="{sx2 + tx * TICK:.2f}" y2="{sy2 + ty * TICK:.2f}"'
                   f' stroke="{DIM_COLOR}" stroke-width="{DIM_SW}"/>')

        # Opening width in inches
        seg_len = math.sqrt((F_B[0] - F_A[0])**2 + (F_B[1] - F_A[1])**2)
        width_ft = seg_len * (op.t_end - op.t_start)
        inches = width_ft * 12
        label = f"{inches:.2f}".rstrip('0').rstrip('.') + '&#8243;'

        # Connecting line between endpoints
        out.append(f'<line x1="{sx1:.2f}" y1="{sy1:.2f}"'
                   f' x2="{sx2:.2f}" y2="{sy2:.2f}"'
                   f' stroke="{DIM_COLOR}" stroke-width="{DIM_SW}"/>')

        # Label: centered, offset toward exterior
        mx = (sx1 + sx2) / 2 + enx * LABEL_OFFSET
        my = (sy1 + sy2) / 2 + eny * LABEL_OFFSET

        # Rotation: text parallel to dim line, kept readable
        svg_angle = math.degrees(math.atan2(udy, udx))
        if svg_angle > 90:
            svg_angle -= 180
        elif svg_angle < -90:
            svg_angle += 180
        rot = (f' transform="rotate({svg_angle:.1f},{mx:.1f},{my:.1f})"'
               if abs(svg_angle) > 0.1 else "")

        out.append(f'<text x="{mx:.1f}" y="{my:.1f}" text-anchor="middle"'
                   f' dominant-baseline="central" font-family="Arial"'
                   f' font-size="{FONT_SIZE}" fill="{DIM_COLOR}"{rot}>'
                   f'{label}</text>')


def render_walls_svg(data, *, title="Outer Walls", include_interior=False):
    """Render the wall detail SVG. Returns SVG string."""
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

    out = []
    out.append(f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}"'
               f' viewBox="{data.vb_x:.2f} {data.vb_y:.2f}'
               f' {data.vb_w:.2f} {data.vb_h:.2f}">')
    out.append(f'<rect x="{data.vb_x:.2f}" y="{data.vb_y:.2f}"'
               f' width="{data.vb_w:.2f}" height="{data.vb_h:.2f}" fill="white"/>')

    # Title
    out.append(f'<text x="{data.title_x:.1f}" y="{data.title_y:.1f}"'
               f' text-anchor="middle" font-family="Arial" font-size="14"'
               f' font-weight="bold">{title}</text>')

    WALL_FILL = "rgba(180,180,180,0.5)"
    OPENING_FILL = "rgb(220,235,255)"

    # --- Draw wall sections ---
    for seg_idx in range(22):
        seg = outline_segs[seg_idx]
        inner_seg = inner_segs[seg_idx]
        s_seg = s_segs[seg_idx]
        g_seg = g_segs[seg_idx]

        seg_openings = openings_on_seg(openings, seg_idx)

        if isinstance(seg, ArcSeg):
            # Arc segments have no openings — draw full strips
            # Outer shell: F-arc to S-arc
            outer_shell = arc_strip_poly(seg, pts, "F", s_seg)
            _svg_polygon(out, outer_shell, to_svg, WALL_FILL, stroke="none")

            # Inner shell: G-arc to W-arc
            if seg_idx == 8:
                # F8-F9: straight-arc-straight path for inner shell
                inner_shell = (list(data.g_f8f9_poly)
                               + list(reversed(data.w_f8f9_poly)))
            else:
                inner_shell = arc_strip_poly(g_seg, pts, "G", inner_seg)
            _svg_polygon(out, inner_shell, to_svg, WALL_FILL, stroke="none")

        elif isinstance(seg, LineSeg):
            if not seg_openings:
                # No openings — draw full rectangle strips
                outer_strip = line_strip_poly(pts, seg.start, seg.end,
                                               s_seg.start, s_seg.end)
                _svg_polygon(out, outer_strip, to_svg, WALL_FILL, stroke="none")

                inner_strip = line_strip_poly(pts, g_seg.start, g_seg.end,
                                               inner_seg.start, inner_seg.end)
                _svg_polygon(out, inner_strip, to_svg, WALL_FILL, stroke="none")
            else:
                # Has openings — draw solid sections and U-turns
                sr = solid_ranges(seg_openings)

                # Shrink ranges so shells end where U-turn arcs begin
                F_A, F_B = pts[seg.start], pts[seg.end]
                seg_len = math.sqrt((F_B[0]-F_A[0])**2 +
                                    (F_B[1]-F_A[1])**2)
                delta_t = R_out / seg_len
                adjusted = []
                for t_s, t_e in sr:
                    if t_s > 1e-9:   # borders an opening end
                        t_s += delta_t
                    if t_e < 1.0 - 1e-9:  # borders an opening start
                        t_e -= delta_t
                    if t_e > t_s + 1e-9:
                        adjusted.append((t_s, t_e))

                for t_s, t_e in adjusted:
                    # Outer shell partial strip
                    outer_strip = partial_line_strip(
                        pts, seg, s_seg.start, s_seg.end, t_s, t_e)
                    _svg_polygon(out, outer_strip, to_svg, WALL_FILL, stroke="none")

                    # Inner shell partial strip
                    inner_strip = partial_line_strip_2(
                        pts, g_seg, inner_seg, t_s, t_e)
                    _svg_polygon(out, inner_strip, to_svg, WALL_FILL, stroke="none")

                # Draw U-turns at each opening boundary
                for op in seg_openings:
                    # U-turn at opening start (wall→opening transition)
                    uturn_start = uturn_polygon(
                        pts, outline_segs, inner_segs, s_segs, g_segs,
                        seg_idx, op.t_start, "start", shell_t, R_in, WALL_OUTER)
                    _svg_polygon(out, uturn_start, to_svg, WALL_FILL,
                                 stroke="none")

                    # U-turn at opening end (opening→wall transition)
                    uturn_end = uturn_polygon(
                        pts, outline_segs, inner_segs, s_segs, g_segs,
                        seg_idx, op.t_end, "end", shell_t, R_in, WALL_OUTER)
                    _svg_polygon(out, uturn_end, to_svg, WALL_FILL,
                                 stroke="none")

                # Draw opening regions
                for op in seg_openings:
                    F_A, F_B = pts[seg.start], pts[seg.end]
                    W_A = pts[inner_seg.start]
                    W_B = pts[inner_seg.end]
                    o_poly = [
                        lerp(F_A, F_B, op.t_start),
                        lerp(F_A, F_B, op.t_end),
                        lerp(W_A, W_B, op.t_end),
                        lerp(W_A, W_B, op.t_start),
                    ]
                    _svg_polygon(out, o_poly, to_svg, OPENING_FILL,
                                 stroke="#4682B4", stroke_width="0.5")

    # --- Continuous outlines per wall section ---
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
                f"{to_svg(*p)[0]:.2f},{to_svg(*p)[1]:.2f}" for p in path)
            out.append(f'<polygon points="{svg_pts}" fill="none" '
                       f'stroke="#999" stroke-width="0.3"/>')

    # --- Interior walls (optional) ---
    if include_interior:
        _render_interior_walls(out, data)
        _render_opening_dims(out, data)

    # --- Opening labels ---
    for op in openings:
        seg = outline_segs[op.seg_idx]
        inner_seg = inner_segs[op.seg_idx]
        t_mid = (op.t_start + op.t_end) / 2
        F_A, F_B = pts[seg.start], pts[seg.end]
        W_A, W_B = pts[inner_seg.start], pts[inner_seg.end]
        f_mid = lerp(F_A, F_B, t_mid)
        w_mid = lerp(W_A, W_B, t_mid)
        cx, cn = (f_mid[0] + w_mid[0]) / 2, (f_mid[1] + w_mid[1]) / 2
        sx, sy = to_svg(cx, cn)
        dE, dN = F_B[0] - F_A[0], F_B[1] - F_A[1]
        svg_angle = -math.degrees(math.atan2(dN, dE))
        if svg_angle > 90:
            svg_angle -= 180
        elif svg_angle < -90:
            svg_angle += 180
        rot = (f' transform="rotate({svg_angle:.1f},{sx:.1f},{sy:.1f})"'
               if abs(svg_angle) > 0.1 else "")
        out.append(f'<text x="{sx:.1f}" y="{sy:.1f}" text-anchor="middle"'
                   f' dominant-baseline="central" font-family="Arial"'
                   f' font-size="5" fill="#4682B4" font-weight="bold"'
                   f'{rot}>{op.name}</text>')

    # --- Title block ---
    out.append(f'<rect x="{data.tb_left:.1f}" y="{data.tb_top:.1f}"'
               f' width="{data.tb_w}" height="{data.tb_h}"'
               f' fill="white" stroke="#333" stroke-width="1"/>')
    out.append(f'<text x="{data.tb_cx:.1f}" y="{data.tb_top+14:.1f}"'
               f' text-anchor="middle" font-family="Arial" font-size="11"'
               f' font-weight="bold" fill="#333">'
               f'{data.outer_area:.2f} sq ft</text>')
    out.append(f'<text x="{data.tb_cx:.1f}" y="{data.tb_top+26:.1f}"'
               f' text-anchor="middle" font-family="Arial" font-size="8"'
               f' fill="#666">Exterior area</text>')

    _ratio = data.ft_per_inch * 12
    _scale_label = f'Scale 1:{_ratio:.1f} 1&#8243; = {fmt_dist(data.ft_per_inch)}'
    out.append(f'<text x="{data.tb_cx:.1f}" y="{data.tb_top+40:.1f}"'
               f' text-anchor="middle" font-family="Arial" font-size="8"'
               f' fill="#666">{_scale_label}</text>')

    _now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    _git_desc = git_describe()
    out.append(f'<text x="{data.tb_cx:.1f}" y="{data.tb_top+54:.1f}"'
               f' text-anchor="middle" font-family="Arial" font-size="7.5"'
               f' fill="#999">Generated {_now}</text>')
    out.append(f'<text x="{data.tb_cx:.1f}" y="{data.tb_top+64:.1f}"'
               f' text-anchor="middle" font-family="Arial" font-size="7.5"'
               f' fill="#999">from {_git_desc}</text>')

    # Wall construction note
    out.append(f'<text x="{data.tb_cx:.1f}" y="{data.tb_top+76:.1f}"'
               f' text-anchor="middle" font-family="Arial" font-size="7"'
               f' fill="#999">2&#8243; shell / 4&#8243; gap / 2&#8243; shell</text>')

    # --- Wall segment table ---
    sections = enumerate_wall_sections(openings, outline_segs)
    # Rotate so O11-O1 (last section) comes first
    sections = sections[-1:] + sections[:-1]

    tbl_left = data.tb_left
    tbl_top = data.tb_bottom + 12
    row_h = 7.5
    # Column right-edges (From-To is left-aligned, others right-aligned)
    col_r = [tbl_left + 32, tbl_left + 62, tbl_left + 92, tbl_left + 128]

    # U-turn centerline length (same for every section)
    R_mid = R_in + shell_t / 2          # centerline radius through shell
    uturn_straight = WALL_OUTER - 2 * (shell_t + R_in)
    uturn_cl = 2 * (math.pi / 2) * R_mid + uturn_straight  # feet

    # Compute row data
    table_rows = []
    for start_op, end_op in sections:
        label = f"{start_op.name}&#8211;{end_op.name}"
        s_seg = start_op.seg_idx
        s_t = start_op.t_end
        e_seg = end_op.seg_idx
        e_t = end_op.t_start

        outer_ft = _path_length_between(
            pts, outline_segs, s_seg, s_t, e_seg, e_t, 0.0)
        inner_ft = _path_length_between(
            pts, outline_segs, s_seg, s_t, e_seg, e_t, WALL_OUTER)
        outer_cl_ft = _path_length_between(
            pts, outline_segs, s_seg, s_t, e_seg, e_t, shell_t / 2)
        inner_cl_ft = _path_length_between(
            pts, outline_segs, s_seg, s_t, e_seg, e_t,
            WALL_OUTER - shell_t / 2)
        shell_ft = (outer_cl_ft - 2 * R_out) + (inner_cl_ft - 2 * R_out) + 2 * uturn_cl

        table_rows.append((label,
                           outer_ft * 12, inner_ft * 12, shell_ft * 12))

    # Table title
    out.append(f'<text x="{(tbl_left + col_r[-1]) / 2:.1f}" y="{tbl_top:.1f}"'
               f' text-anchor="middle" font-family="Arial" font-size="7"'
               f' font-weight="bold" fill="#333">Wall Segments</text>')

    # Column headers
    hdr_y = tbl_top + 10
    hdrs = ["From&#8211;To", "Outer (in)", "Inner (in)", "Shell (in)"]
    hdr_x = [tbl_left + 2, col_r[1] - 2, col_r[2] - 2, col_r[3] - 2]
    hdr_anchor = ["start", "end", "end", "end"]
    for hx, ha, hd in zip(hdr_x, hdr_anchor, hdrs):
        out.append(f'<text x="{hx:.1f}" y="{hdr_y:.1f}"'
                   f' text-anchor="{ha}" font-family="Arial" font-size="6"'
                   f' font-weight="bold" fill="#333">{hd}</text>')

    # Header underline
    line_y = hdr_y + 2.5
    out.append(f'<line x1="{tbl_left:.1f}" y1="{line_y:.1f}"'
               f' x2="{col_r[-1]:.1f}" y2="{line_y:.1f}"'
               f' stroke="#999" stroke-width="0.5"/>')

    # Data rows
    for ri, (label, o_in, i_in, s_in) in enumerate(table_rows):
        y = line_y + (ri + 1) * row_h
        vals = [label, f"{o_in:.2f}", f"{i_in:.2f}", f"{s_in:.2f}"]
        for vx, va, vv in zip(hdr_x, hdr_anchor, vals):
            out.append(f'<text x="{vx:.1f}" y="{y:.1f}"'
                       f' text-anchor="{va}" font-family="Arial"'
                       f' font-size="6" fill="#333">{vv}</text>')

    # Total row (separated by a line)
    total_line_y = line_y + len(table_rows) * row_h + 2
    out.append(f'<line x1="{tbl_left:.1f}" y1="{total_line_y:.1f}"'
               f' x2="{col_r[-1]:.1f}" y2="{total_line_y:.1f}"'
               f' stroke="#999" stroke-width="0.5"/>')
    tot_o = sum(r[1] for r in table_rows)
    tot_i = sum(r[2] for r in table_rows)
    tot_s = sum(r[3] for r in table_rows)
    tot_y = total_line_y + row_h
    tot_vals = ["Total", f"{tot_o:.1f}", f"{tot_i:.1f}", f"{tot_s:.1f}"]
    for vx, va, vv in zip(hdr_x, hdr_anchor, tot_vals):
        out.append(f'<text x="{vx:.1f}" y="{tot_y:.1f}"'
                   f' text-anchor="{va}" font-family="Arial"'
                   f' font-size="6" font-weight="bold" fill="#333">{vv}</text>')

    # "in feet" row
    ft_y = tot_y + row_h
    ft_vals = ["in feet", f"{tot_o / 12:.1f}", f"{tot_i / 12:.1f}", f"{tot_s / 12:.1f}"]
    for vx, va, vv in zip(hdr_x, hdr_anchor, ft_vals):
        out.append(f'<text x="{vx:.1f}" y="{ft_y:.1f}"'
                   f' text-anchor="{va}" font-family="Arial"'
                   f' font-size="6" fill="#333">{vv}</text>')

    # Table border
    tbl_border_top = tbl_top - 8.5
    tbl_border_bottom = ft_y + 3
    out.append(f'<rect x="{tbl_left:.1f}" y="{tbl_border_top:.1f}"'
               f' width="{col_r[-1] - tbl_left:.1f}"'
               f' height="{tbl_border_bottom - tbl_border_top:.1f}"'
               f' fill="none" stroke="#999" stroke-width="0.5"/>')

    out.append('</svg>')
    return "\n".join(out)



# ============================================================
# Main entry point
# ============================================================

if __name__ == "__main__":
    data = build_wall_data()
    _dir = os.path.dirname(os.path.abspath(__file__))

    svg_content = render_walls_svg(data)
    svg_path = os.path.join(_dir, "walls.svg")
    with open(svg_path, "w", encoding="utf-8") as f:
        f.write(svg_content)
    print(f"Wall detail written to {svg_path}")

    all_svg = render_walls_svg(data, title="Walls", include_interior=True)
    all_path = os.path.join(_dir, "all_walls.svg")
    with open(all_path, "w", encoding="utf-8") as f:
        f.write(all_svg)
    print(f"All walls written to {all_path}")

    print(f"Shell: {SHELL_THICKNESS * 12:.0f}\" / "
          f"Gap: {AIR_GAP * 12:.0f}\" / "
          f"Shell: {SHELL_THICKNESS * 12:.0f}\"")
    print(f"Opening corner inside radius: {OPENING_INSIDE_RADIUS * 12:.1f}\"")
