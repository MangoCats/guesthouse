"""Generate outer wall detail SVG showing double-shell concrete construction.

Outer walls are double-shell 3D-printed concrete: two 2" shells separated
by a 6" air gap (10" total = WALL_OUTER). At openings, shells connect via
90-degree corner turns with configurable radii.

Outputs walls/walls.svg at 1:72 scale.
"""
import os, sys, math, datetime
from typing import NamedTuple

# Ensure project root is on sys.path for package imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from shared.types import LineSeg, ArcSeg
from shared.geometry import fmt_dist, left_norm, segment_polyline
from shared.svg import make_svg_transform, W, H, git_describe, normalize_svg_angle, svg_polygon_pts
from floorplan.constants import WALL_OUTER, SHELL_THICKNESS, AIR_GAP, OPENING_INSIDE_RADIUS
from floorplan.openings import compute_rough_openings
from shared.wall_shells import (
    lerp, openings_on_seg, solid_ranges,
    arc_strip_poly, line_strip_poly, partial_line_strip,
    uturn_arc_data, uturn_polygon,
    trace_boundary_path, enumerate_wall_sections, build_section_outlines,
)


# ── Color palette ─────────────────────────────────────────────
CLR_IW_FILL = "rgba(160,160,160,0.35)"
CLR_IW_STROKE = "#666"
CLR_DIM = "#4682B4"
CLR_SHELL_FILL = "rgba(180,180,180,0.5)"
CLR_OPENING_FILL = "rgb(220,235,255)"
CLR_LABEL = "#666"
CLR_TITLE = "#333"


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
            if seg.start == "F8" and seg.end == "F9" and inset >= SHELL_THICKNESS + AIR_GAP:
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
    roof: object           # RoofGeometry
    constants: dict        # DB constants (SHELL_THICKNESS, AIR_GAP, etc.)
    iw_polys: dict         # DB interior wall polygons, or None for seed path
    ro_polys: list         # DB rough opening dicts, or None for seed path
    roof_poly: list        # Roof outline polygon (precomputed)
    glazing_rows: list     # Room glazing table rows, or None


def _get_roof_poly(gd):
    """Return roof polygon from a GeneratorData, handling both old and DB roof types."""
    roof = getattr(gd, 'roof', None)
    if roof is None:
        return []
    # DbRoofResult has .poly directly
    if hasattr(roof, 'poly') and roof.poly:
        return roof.poly
    # Old RoofGeometry requires roof_polyline()
    try:
        from floorplan.roof import roof_polyline
        return roof_polyline(roof)
    except Exception:
        return []


def build_wall_data(gd=None):
    """Compute all geometry needed for the wall detail SVG.

    If gd (GeneratorData) is provided, uses it as the geometry source.
    Otherwise constructs one from the hardcoded procedural modules.
    """
    if gd is None:
        from floorplan.gen_floorplan import build_floorplan_data, _default_constants_dict
        from app.gen_provider import GeneratorData, compute_native_geometry
        constants_dict = _default_constants_dict()
        pts, outline_segs, inner_segs, radii = compute_native_geometry(
            constants_dict)
        gd = GeneratorData(pts, outline_segs, inner_segs, radii, constants_dict)

    to_svg = make_svg_transform()

    # --- Page layout: 1:72 scale ---
    # Use actual outline segment point names (works with any naming scheme)
    _f_names = list(dict.fromkeys(
        n for s in gd.outline_segs for n in (s.start, s.end) if n in gd.pts))
    _f_svg = [to_svg(*gd.pts[k]) for k in _f_names]
    _bldg_xmin = min(p[0] for p in _f_svg)
    _bldg_xmax = max(p[0] for p in _f_svg)
    _bldg_ymin = min(p[1] for p in _f_svg)
    _bldg_cx = (_bldg_xmin + _bldg_xmax) / 2

    _title_y = _bldg_ymin - 49

    _tb_w, _tb_h = 130, 80
    _tb_left = _bldg_xmax + 30
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
        pts=gd.pts, to_svg=to_svg,
        outline_segs=gd.outline_segs, inner_segs=gd.inner_segs,
        s_segs=gd.s_segs, g_segs=gd.g_segs,
        radii=gd.radii, openings=gd.openings,
        layout=gd.layout, inner_poly=gd.inner_poly,
        outer_area=gd.outer_area,
        inner_area=gd.inner_area,
        vb_x=_vb_x, vb_y=_vb_y, vb_w=_vb_w, vb_h=_vb_h,
        title_x=_bldg_cx, title_y=_title_y,
        tb_left=_tb_left, tb_right=_tb_right, tb_top=_tb_top,
        tb_bottom=_tb_bottom, tb_w=_tb_w, tb_h=_tb_h, tb_cx=_tb_cx,
        ft_per_inch=_ft_per_inch,
        g_f8f9_poly=gd.g_f8f9_poly,
        w_f8f9_poly=gd.w_f8f9_poly,
        roof=gd.roof,
        constants=gd.constants,
        iw_polys=getattr(gd, 'iw_polys', None),
        ro_polys=getattr(gd, 'ro_polys', None),
        roof_poly=_get_roof_poly(gd),
        glazing_rows=getattr(gd, 'glazing_rows', None),
    )


# ============================================================
# SVG rendering
# ============================================================

def _svg_polygon(out, poly, to_svg, fill, stroke=CLR_IW_STROKE, stroke_width="0.5"):
    """Render a polygon as an SVG element."""
    svg = svg_polygon_pts(poly, to_svg, prec=2)
    out.append(f'<polygon points="{svg}" fill="{fill}" '
               f'stroke="{stroke}" stroke-width="{stroke_width}"/>')


def _render_interior_walls(out, data):
    """Render interior wall polygons and labels into the SVG output list."""
    pts = data.pts
    to_svg = data.to_svg
    layout = data.layout

    if layout is None:
        # DB path — render iw_polys and rough openings from DB
        _RO_COLOR = "darkred"
        _RO_SW = "0.5"
        _LABEL_GAP = 3.0

        if data.iw_polys:
            for name, poly in sorted(data.iw_polys.items()):
                svg = svg_polygon_pts(poly, to_svg, prec=2)
                out.append(f'<polygon points="{svg}" fill="{CLR_IW_FILL}"'
                           f' stroke="{CLR_IW_STROKE}" stroke-width="0.5"/>')
                if len(poly) >= 3:
                    cx = sum(p[0] for p in poly) / len(poly)
                    cy = sum(p[1] for p in poly) / len(poly)
                    sx, sy = to_svg(cx, cy)
                    out.append(f'<text x="{sx:.1f}" y="{sy:.1f}"'
                               f' text-anchor="middle" dominant-baseline="central"'
                               f' font-family="Arial" font-size="6"'
                               f' fill="{CLR_LABEL}">{name}</text>')

        for ro in (data.ro_polys or []):
            b = ro.get("bbox") or {}
            poly = ro.get("poly")
            orient = ro.get("orientation", "H")
            name = ro.get("name", "")

            if poly:
                _rp_svg = [to_svg(*p) for p in poly]
                pts_str = " ".join(f"{x:.1f},{y:.1f}" for x, y in _rp_svg)
                out.append(f'<polygon points="{pts_str}" fill="none"'
                           f' stroke="{_RO_COLOR}" stroke-width="{_RO_SW}"/>')
                out.append(f'<line x1="{_rp_svg[0][0]:.1f}" y1="{_rp_svg[0][1]:.1f}"'
                           f' x2="{_rp_svg[2][0]:.1f}" y2="{_rp_svg[2][1]:.1f}"'
                           f' stroke="{_RO_COLOR}" stroke-width="{_RO_SW}"/>')
                out.append(f'<line x1="{_rp_svg[1][0]:.1f}" y1="{_rp_svg[1][1]:.1f}"'
                           f' x2="{_rp_svg[3][0]:.1f}" y2="{_rp_svg[3][1]:.1f}"'
                           f' stroke="{_RO_COLOR}" stroke-width="{_RO_SW}"/>')
                cx = sum(p[0] for p in poly) / len(poly)
                cy = sum(p[1] for p in poly) / len(poly)
                lx, ly = to_svg(cx, cy)
                rot = ""
            elif b:
                x1, y1 = to_svg(b["w"], b["n"])
                x2, y2 = to_svg(b["e"], b["s"])
                out.append(f'<rect x="{x1:.1f}" y="{y1:.1f}"'
                           f' width="{x2 - x1:.1f}" height="{y2 - y1:.1f}"'
                           f' fill="none" stroke="{_RO_COLOR}" stroke-width="{_RO_SW}"/>')
                out.append(f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}"'
                           f' stroke="{_RO_COLOR}" stroke-width="{_RO_SW}"/>')
                out.append(f'<line x1="{x2:.1f}" y1="{y1:.1f}" x2="{x1:.1f}" y2="{y2:.1f}"'
                           f' stroke="{_RO_COLOR}" stroke-width="{_RO_SW}"/>')
                if orient == "H":
                    lx, ly = to_svg((b["w"] + b["e"]) / 2, b["n"])
                    ly -= _LABEL_GAP
                    rot = ""
                else:
                    lx, ly = to_svg(b["w"], (b["s"] + b["n"]) / 2)
                    lx -= _LABEL_GAP
                    rot = f' transform="rotate(-90 {lx:.1f} {ly:.1f})"'
            else:
                continue

            out.append(f'<text x="{lx:.1f}" y="{ly:.1f}" text-anchor="middle"'
                       f' dominant-baseline="central" font-family="Arial"'
                       f' font-size="6" fill="{_RO_COLOR}"{rot}>{name}</text>')

        return

    inner_poly = data.inner_poly

    IW_FILL = CLR_IW_FILL
    IW_STROKE = CLR_IW_STROKE
    IW_SW = "0.5"
    LABEL_SIZE = "6"
    LABEL_GAP = 3.0  # SVG px from wall face to label center

    def iw_poly(poly):
        svg = svg_polygon_pts(poly, to_svg, prec=2)
        out.append(f'<polygon points="{svg}" fill="{IW_FILL}" '
                   f'stroke="{IW_STROKE}" stroke-width="{IW_SW}"/>')

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
                   f' font-size="{LABEL_SIZE}" fill="{CLR_LABEL}"{rot}>{name}</text>')

    rough_openings = compute_rough_openings(pts, layout)
    ro_map = {ro.name: ro.bbox for ro in rough_openings}

    # IW1 (horizontal, 6")
    iw_poly(layout.iw1.poly)
    iw_label("IW1", layout.iw1.poly[0][0], layout.iw1.poly[1][0],
             layout.iw1.s, layout.iw1.n, vertical=False)

    # IW8 (horizontal, 6" — west extension of IW1)
    iw_poly(layout.iw8.poly)
    iw_label("IW8", layout.iw8.w, layout.iw8.e, layout.iw8.s, layout.iw8.n,
             vertical=False)

    # IW2 (vertical, 6" — lower segment)
    iw_poly(layout.iw2.poly)
    iw_label("IW2", layout.iw2.w, layout.iw2.e, layout.iw2.s, layout.iw2.n)

    # IW2o (oblique, 6" — connector with RO4)
    iw_poly(layout.iw2o.poly)
    iw_label("IW2o", layout.iw2o.w, layout.iw2o.e, layout.iw2o.s, layout.iw2o.n)

    # IW2s (vertical, 6" — shower segment)
    iw_poly(layout.iw2s.poly)
    iw_label("IW2s", layout.iw2s.w, layout.iw2s.e, layout.iw2s.s, layout.iw2s.n)

    # IW3 (rotated, 4" thick, perpendicular to W18-W1)
    iw_poly(layout.iw3.poly)
    iw_label("IW3", layout.iw3.w, layout.iw3.e, layout.iw3.s, layout.iw3.n)

    # IW7 (rotated, 4" thick, parallel to W18-W1)
    iw_poly(layout.iw7.poly)
    iw_label("IW7", layout.iw7.w, layout.iw7.e, layout.iw7.s, layout.iw7.n,
             vertical=False)

    # IW9 (rotated, 4" thick, perpendicular to W18-W1)
    iw_poly(layout.iw9.poly)
    iw_label("IW9", layout.iw9.w, layout.iw9.e,
             ro_map["RO7"].n, ro_map["RO3"].s)

    # IW6 (horizontal, 1" partition)
    iw_poly(layout.iw6.poly)
    _ro5_w = ro_map["RO5"].w
    iw_label("IW6", _ro5_w, _ro5_w, layout.iw6.s, layout.iw6.n, vertical=False)

    # IW4 (vertical, 4")
    iw_poly(layout.iw4.poly)
    iw_label("IW4", layout.iw4.w, layout.iw4.e, layout.iw4.s, layout.iw4.n)

    # IW11 (vertical, 4")
    iw_poly(layout.iw11.poly)
    iw_label("IW11", layout.iw11.w, layout.iw11.e, layout.iw11.s, layout.iw11.n)

    # IW12 (rotated, 4")
    iw_poly(layout.iw12.poly)
    iw_label("IW12", layout.iw12.w, layout.iw12.e, layout.iw12.s, layout.iw12.n,
             vertical=False)

    # IW5 (horizontal, 3")
    iw_poly(layout.iw5.poly)
    iw_label("IW5", layout.iw5.w, layout.iw5.e, layout.iw5.s, layout.iw5.n,
             vertical=False)

    # Rough openings (RO1-RO7) — dark red outline box with X (7 openings)
    _RO_COLOR = "darkred"
    _RO_SW = "0.5"
    for ro in rough_openings:
        if ro.poly is not None:
            # Rotated opening: draw polygon + diagonals
            _rp = ro.poly  # [SW, SE, NE, NW]
            _rp_svg = [(to_svg(*p)) for p in _rp]
            pts_str = " ".join(f"{x:.1f},{y:.1f}" for x, y in _rp_svg)
            out.append(f'<polygon points="{pts_str}" fill="none"'
                       f' stroke="{_RO_COLOR}" stroke-width="{_RO_SW}"/>')
            out.append(f'<line x1="{_rp_svg[0][0]:.1f}" y1="{_rp_svg[0][1]:.1f}"'
                       f' x2="{_rp_svg[2][0]:.1f}" y2="{_rp_svg[2][1]:.1f}"'
                       f' stroke="{_RO_COLOR}" stroke-width="{_RO_SW}"/>')
            out.append(f'<line x1="{_rp_svg[1][0]:.1f}" y1="{_rp_svg[1][1]:.1f}"'
                       f' x2="{_rp_svg[3][0]:.1f}" y2="{_rp_svg[3][1]:.1f}"'
                       f' stroke="{_RO_COLOR}" stroke-width="{_RO_SW}"/>')
        else:
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
        if ro.orientation == "R" and ro.poly is not None:
            # Rotated opening: label on WNW face of IW11
            import math as _m
            _iw11 = layout.iw11.poly  # [SW, SE, NE, NW]
            # WNW face direction: SW→NW
            _face_dx = _iw11[3][0] - _iw11[0][0]
            _face_dy = _iw11[3][1] - _iw11[0][1]
            _face_len = _m.sqrt(_face_dx**2 + _face_dy**2)
            _face_ux = _face_dx / _face_len
            _face_uy = _face_dy / _face_len
            # Label at center of RO2 poly, offset toward WNW face
            _rc = (sum(p[0] for p in ro.poly) / 4,
                   sum(p[1] for p in ro.poly) / 4)
            # Offset toward west face (SE→SW direction = _iw11_at)
            _thick_ux = _iw11[0][0] - _iw11[1][0]
            _thick_uy = _iw11[0][1] - _iw11[1][1]
            _thick_l = _m.sqrt(_thick_ux**2 + _thick_uy**2)
            _lpt = (_rc[0] + (_thick_ux / _thick_l) * LABEL_GAP / abs(to_svg(1, 0)[0] - to_svg(0, 0)[0]),
                    _rc[1] + (_thick_uy / _thick_l) * LABEL_GAP / abs(to_svg(1, 0)[0] - to_svg(0, 0)[0]))
            lx, ly = to_svg(*_lpt)
            # SVG rotation: angle of face direction in SVG coords
            _sx1, _sy1 = to_svg(*_iw11[0])
            _sx2, _sy2 = to_svg(*_iw11[3])
            _svg_ang = _m.degrees(_m.atan2(_sy2 - _sy1, _sx2 - _sx1))
            # Shift 1/3 font height "up" (outward normal to WNW face)
            _fdx = _sx2 - _sx1; _fdy = _sy2 - _sy1
            _fl = _m.sqrt(_fdx**2 + _fdy**2)
            _nudge = float(LABEL_SIZE) / 3.0
            lx += (_fdy / _fl) * _nudge
            ly += (-_fdx / _fl) * _nudge
            rot = f' transform="rotate({_svg_ang:.1f} {lx:.1f} {ly:.1f})"'
        elif ro.orientation == "H":
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
    DIM_COLOR = CLR_DIM
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
        svg_angle = normalize_svg_angle(math.degrees(math.atan2(udy, udx)))
        rot = (f' transform="rotate({svg_angle:.1f},{mx:.1f},{my:.1f})"'
               if abs(svg_angle) > 0.1 else "")

        out.append(f'<text x="{mx:.1f}" y="{my:.1f}" text-anchor="middle"'
                   f' dominant-baseline="central" font-family="Arial"'
                   f' font-size="{FONT_SIZE}" fill="{DIM_COLOR}"{rot}>'
                   f'{label}</text>')


def _render_wall_segments(out, data):
    """Render wall section fills (arcs, lines, openings, U-turns)."""
    pts = data.pts
    to_svg = data.to_svg
    outline_segs = data.outline_segs
    inner_segs = data.inner_segs
    s_segs = data.s_segs
    g_segs = data.g_segs
    openings = data.openings

    _c = data.constants
    shell_t = _c.get("SHELL_THICKNESS", SHELL_THICKNESS)
    R_in = _c.get("OPENING_INSIDE_RADIUS", OPENING_INSIDE_RADIUS)
    R_out = R_in + shell_t
    _wall_outer = _c.get("WALL_OUTER", WALL_OUTER)

    WALL_FILL = CLR_SHELL_FILL
    OPENING_FILL = CLR_OPENING_FILL

    for seg_idx in range(len(outline_segs)):
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
            if seg.start == "F8" and seg.end == "F9":
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
                        pts, seg, s_seg, t_s, t_e)
                    _svg_polygon(out, outer_strip, to_svg, WALL_FILL, stroke="none")

                    # Inner shell partial strip
                    inner_strip = partial_line_strip(
                        pts, g_seg, inner_seg, t_s, t_e)
                    _svg_polygon(out, inner_strip, to_svg, WALL_FILL, stroke="none")

                # Draw U-turns at each opening boundary
                for op in seg_openings:
                    # U-turn at opening start (wall→opening transition)
                    uturn_start = uturn_polygon(
                        pts, outline_segs, inner_segs, s_segs, g_segs,
                        seg_idx, op.t_start, "start", shell_t, R_in, _wall_outer)
                    _svg_polygon(out, uturn_start, to_svg, WALL_FILL,
                                 stroke="none")

                    # U-turn at opening end (opening→wall transition)
                    uturn_end = uturn_polygon(
                        pts, outline_segs, inner_segs, s_segs, g_segs,
                        seg_idx, op.t_end, "end", shell_t, R_in, _wall_outer)
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
                                 stroke=CLR_DIM, stroke_width="0.5")


def _render_section_outlines(out, data):
    """Render continuous outlines per wall section."""
    pts = data.pts
    to_svg = data.to_svg
    _c = data.constants
    shell_t = _c.get("SHELL_THICKNESS", SHELL_THICKNESS)
    R_in = _c.get("OPENING_INSIDE_RADIUS", OPENING_INSIDE_RADIUS)
    _wall_outer = _c.get("WALL_OUTER", WALL_OUTER)

    _f8f9_idx = next((i for i, s in enumerate(data.outline_segs)
                      if s.start == "F8" and s.end == "F9"), None)
    g_overrides = {_f8f9_idx: data.g_f8f9_poly} if _f8f9_idx is not None else {}
    w_overrides = {_f8f9_idx: data.w_f8f9_poly} if _f8f9_idx is not None else {}
    sections = enumerate_wall_sections(data.openings, data.outline_segs)
    for start_op, end_op in sections:
        outer_path, cavity_path = build_section_outlines(
            pts, data.outline_segs, data.inner_segs, data.s_segs, data.g_segs,
            start_op, end_op, shell_t, R_in, _wall_outer,
            g_seg_overrides=g_overrides, w_seg_overrides=w_overrides)
        for path in [outer_path, cavity_path]:
            svg_pts = svg_polygon_pts(path, to_svg, prec=2)
            out.append(f'<polygon points="{svg_pts}" fill="none" '
                       f'stroke="#999" stroke-width="0.3"/>')


def _render_opening_labels(out, data):
    """Render opening name labels on the wall face."""
    pts = data.pts
    to_svg = data.to_svg
    outline_segs = data.outline_segs
    inner_segs = data.inner_segs

    for op in data.openings:
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
        svg_angle = normalize_svg_angle(-math.degrees(math.atan2(dN, dE)))
        rot = (f' transform="rotate({svg_angle:.1f},{sx:.1f},{sy:.1f})"'
               if abs(svg_angle) > 0.1 else "")
        out.append(f'<text x="{sx:.1f}" y="{sy:.1f}" text-anchor="middle"'
                   f' dominant-baseline="central" font-family="Arial"'
                   f' font-size="5" fill="{CLR_DIM}" font-weight="bold"'
                   f'{rot}>{op.name}</text>')


def _render_title_block(out, data):
    """Render the title block with area, scale, timestamp."""
    out.append(f'<rect x="{data.tb_left:.1f}" y="{data.tb_top:.1f}"'
               f' width="{data.tb_w}" height="{data.tb_h}"'
               f' fill="white" stroke="{CLR_TITLE}" stroke-width="1"/>')
    out.append(f'<text x="{data.tb_cx:.1f}" y="{data.tb_top+14:.1f}"'
               f' text-anchor="middle" font-family="Arial" font-size="11"'
               f' font-weight="bold" fill="{CLR_TITLE}">'
               f'{data.outer_area:.2f} sq ft</text>')
    out.append(f'<text x="{data.tb_cx:.1f}" y="{data.tb_top+26:.1f}"'
               f' text-anchor="middle" font-family="Arial" font-size="8"'
               f' fill="{CLR_LABEL}">Exterior area</text>')

    _ratio = data.ft_per_inch * 12
    _scale_label = f'Scale 1:{_ratio:.1f} 1&#8243; = {fmt_dist(data.ft_per_inch)}'
    out.append(f'<text x="{data.tb_cx:.1f}" y="{data.tb_top+40:.1f}"'
               f' text-anchor="middle" font-family="Arial" font-size="8"'
               f' fill="{CLR_LABEL}">{_scale_label}</text>')

    _now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    _git_desc = git_describe()
    out.append(f'<text x="{data.tb_cx:.1f}" y="{data.tb_top+54:.1f}"'
               f' text-anchor="middle" font-family="Arial" font-size="7.5"'
               f' fill="#999">Generated {_now}</text>')
    out.append(f'<text x="{data.tb_cx:.1f}" y="{data.tb_top+64:.1f}"'
               f' text-anchor="middle" font-family="Arial" font-size="7.5"'
               f' fill="#999">from {_git_desc}</text>')

    # Wall construction note
    _c = data.constants
    _shell_in = _c.get("SHELL_THICKNESS", SHELL_THICKNESS) * 12
    _gap_in = _c.get("AIR_GAP", AIR_GAP) * 12
    out.append(f'<text x="{data.tb_cx:.1f}" y="{data.tb_top+76:.1f}"'
               f' text-anchor="middle" font-family="Arial" font-size="7"'
               f' fill="#999">{_shell_in:.0f}&#8243; shell / {_gap_in:.0f}&#8243; gap / {_shell_in:.0f}&#8243; shell</text>')


def _compute_slope_wall_area(pts, outline_segs, inner_segs):
    """Compute total wall face area from 80" up to the 2:12 roof bottom.

    The 2:12 roof rises 2" per foot of northing from a reference elevation
    of 7'6" (90") at F18 northing (-13.5 ft).  Height above 80" at northing y:
        h(y) = 2*y + 37 inches

    Returns total slope area in sq ft (outer F-series + inner W-series faces).
    """
    # Roof slope: 2" rise per foot of northing
    # At F18 northing (-13.5 ft): roof bottom = 90" (7'6")
    # h_above_80(y) = 2*y + 37 inches
    def _h(y):
        return 2.0 * y + 37.0  # inches above 80" at northing y (feet)

    def _seg_slope_area(segs):
        """Sum trapezoidal slope area for a list of segments."""
        total = 0.0
        for seg in segs:
            poly = segment_polyline(seg, pts)
            for i in range(len(poly) - 1):
                p1, p2 = poly[i], poly[i + 1]
                dx = p2[0] - p1[0]
                dy = p2[1] - p1[1]
                seg_len_ft = math.sqrt(dx * dx + dy * dy)
                h1 = _h(p1[1])
                h2 = _h(p2[1])
                # Trapezoidal area: length_ft * avg_height_in / 12 → sq ft
                total += seg_len_ft * (h1 + h2) / 2.0 / 12.0
        return total

    outer_area = _seg_slope_area(outline_segs)
    inner_area = _seg_slope_area(inner_segs)
    return outer_area + inner_area


def compute_wall_table_rows(pts, outline_segs, openings, constants=None):
    """Compute wall segment measurements for the wall table.

    Returns a list of (label, outer_inches, inner_inches, shell_inches,
    area_2_12_sqft) tuples, one per wall section, rotated so O11-O1 first.
    """
    _c = constants or {}
    shell_t = _c.get("SHELL_THICKNESS", SHELL_THICKNESS)
    R_in = _c.get("OPENING_INSIDE_RADIUS", OPENING_INSIDE_RADIUS)
    _wall_outer = _c.get("WALL_OUTER", WALL_OUTER)
    R_out = R_in + shell_t

    sections = enumerate_wall_sections(openings, outline_segs)
    # Rotate so O11-O1 (last section) comes first
    sections = sections[-1:] + sections[:-1]

    # U-turn centerline length (same for every section)
    R_mid = R_in + shell_t / 2          # centerline radius through shell
    uturn_straight = _wall_outer - 2 * (shell_t + R_in)
    uturn_cl = 2 * (math.pi / 2) * R_mid + uturn_straight  # feet

    rows = []
    for start_op, end_op in sections:
        label = f"{start_op.name}&#8211;{end_op.name}"
        s_seg = start_op.seg_idx
        s_t = start_op.t_end
        e_seg = end_op.seg_idx
        e_t = end_op.t_start

        outer_ft = _path_length_between(
            pts, outline_segs, s_seg, s_t, e_seg, e_t, 0.0)
        inner_ft = _path_length_between(
            pts, outline_segs, s_seg, s_t, e_seg, e_t, _wall_outer)
        outer_cl_ft = _path_length_between(
            pts, outline_segs, s_seg, s_t, e_seg, e_t, shell_t / 2)
        inner_cl_ft = _path_length_between(
            pts, outline_segs, s_seg, s_t, e_seg, e_t,
            _wall_outer - shell_t / 2)
        shell_ft = (outer_cl_ft - 2 * R_out) + (inner_cl_ft - 2 * R_out) + 2 * uturn_cl

        outer_in = outer_ft * 12
        inner_in = inner_ft * 12
        # 2:12 area: both faces × 80" height, converted to sq ft
        area_2_12 = (outer_in + inner_in) * 80.0 / 144.0
        rows.append((label, outer_in, inner_in, shell_ft * 12, area_2_12))

    return rows


def _render_wall_table(out, data):
    """Render the wall segment measurement table. Returns table bottom y."""
    table_rows = compute_wall_table_rows(
        data.pts, data.outline_segs, data.openings, data.constants)
    slope_area = _compute_slope_wall_area(
        data.pts, data.outline_segs, data.inner_segs)

    tbl_left = data.tb_left
    tbl_top = data.tb_bottom + 12
    row_h = 7.5
    # Column right-edges (From-To is left-aligned, others right-aligned)
    col_r = [tbl_left + 32, tbl_left + 62, tbl_left + 92,
             tbl_left + 128, tbl_left + 168]

    # Table title
    out.append(f'<text x="{(tbl_left + col_r[-1]) / 2:.1f}" y="{tbl_top:.1f}"'
               f' text-anchor="middle" font-family="Arial" font-size="7"'
               f' font-weight="bold" fill="{CLR_TITLE}">Wall Segments</text>')

    # Column headers
    hdr_y = tbl_top + 10
    hdrs = ["From&#8211;To", "Outer (in)", "Inner (in)", "Shell (in)"]
    hdr_x = [tbl_left + 2, col_r[1] - 2, col_r[2] - 2, col_r[3] - 2,
             col_r[4] - 2]
    hdr_anchor = ["start", "end", "end", "end", "end"]
    for hx, ha, hd in zip(hdr_x, hdr_anchor, hdrs):
        out.append(f'<text x="{hx:.1f}" y="{hdr_y:.1f}"'
                   f' text-anchor="{ha}" font-family="Arial" font-size="6"'
                   f' font-weight="bold" fill="{CLR_TITLE}">{hd}</text>')
    # 2:12 area header — two lines: "2:12 area" at title level, "(sq ft)" at hdr level
    out.append(f'<text x="{hdr_x[4]:.1f}" y="{tbl_top:.1f}"'
               f' text-anchor="end" font-family="Arial" font-size="6"'
               f' font-weight="bold" fill="{CLR_TITLE}">2:12 area</text>')
    out.append(f'<text x="{hdr_x[4]:.1f}" y="{hdr_y:.1f}"'
               f' text-anchor="end" font-family="Arial" font-size="6"'
               f' font-weight="bold" fill="{CLR_TITLE}">(sq ft)</text>')

    # Header underline
    line_y = hdr_y + 2.5
    out.append(f'<line x1="{tbl_left:.1f}" y1="{line_y:.1f}"'
               f' x2="{col_r[-1]:.1f}" y2="{line_y:.1f}"'
               f' stroke="#999" stroke-width="0.5"/>')

    # Data rows
    for ri, (label, o_in, i_in, s_in, a_sqft) in enumerate(table_rows):
        y = line_y + (ri + 1) * row_h
        vals = [label, f"{o_in:.2f}", f"{i_in:.2f}", f"{s_in:.2f}",
                f"{a_sqft:.1f}"]
        for vx, va, vv in zip(hdr_x, hdr_anchor, vals):
            out.append(f'<text x="{vx:.1f}" y="{y:.1f}"'
                       f' text-anchor="{va}" font-family="Arial"'
                       f' font-size="6" fill="{CLR_TITLE}">{vv}</text>')

    # Slope wall row
    slope_y = line_y + (len(table_rows) + 1) * row_h
    slope_vals = ["Slope wall", "", "", "", f"{slope_area:.1f}"]
    for vx, va, vv in zip(hdr_x, hdr_anchor, slope_vals):
        if vv:
            out.append(f'<text x="{vx:.1f}" y="{slope_y:.1f}"'
                       f' text-anchor="{va}" font-family="Arial"'
                       f' font-size="6" fill="{CLR_TITLE}">{vv}</text>')

    # Total row (separated by a line)
    total_line_y = slope_y + 2
    out.append(f'<line x1="{tbl_left:.1f}" y1="{total_line_y:.1f}"'
               f' x2="{col_r[-1]:.1f}" y2="{total_line_y:.1f}"'
               f' stroke="#999" stroke-width="0.5"/>')
    tot_o = sum(r[1] for r in table_rows)
    tot_i = sum(r[2] for r in table_rows)
    tot_s = sum(r[3] for r in table_rows)
    tot_a = sum(r[4] for r in table_rows) + slope_area
    tot_y = total_line_y + row_h
    tot_vals = ["Total", f"{tot_o:.1f}", f"{tot_i:.1f}", f"{tot_s:.1f}",
                ""]
    for vx, va, vv in zip(hdr_x, hdr_anchor, tot_vals):
        if vv:
            out.append(f'<text x="{vx:.1f}" y="{tot_y:.1f}"'
                       f' text-anchor="{va}" font-family="Arial"'
                       f' font-size="6" font-weight="bold" fill="{CLR_TITLE}">{vv}</text>')

    # "in feet" row
    ft_y = tot_y + row_h
    ft_vals = ["in feet", f"{tot_o / 12:.1f}", f"{tot_i / 12:.1f}",
               f"{tot_s / 12:.1f}", f"{tot_a:.1f}"]
    for vx, va, vv in zip(hdr_x, hdr_anchor, ft_vals):
        out.append(f'<text x="{vx:.1f}" y="{ft_y:.1f}"'
                   f' text-anchor="{va}" font-family="Arial"'
                   f' font-size="6" fill="{CLR_TITLE}">{vv}</text>')

    # Table border
    tbl_border_top = tbl_top - 8.5
    tbl_border_bottom = ft_y + 3
    out.append(f'<rect x="{tbl_left:.1f}" y="{tbl_border_top:.1f}"'
               f' width="{col_r[-1] - tbl_left:.1f}"'
               f' height="{tbl_border_bottom - tbl_border_top:.1f}"'
               f' fill="none" stroke="#999" stroke-width="0.5"/>')

    return tbl_border_bottom


def _render_interior_walls_table_db(out, data, tbl_border_bottom):
    """Render interior walls table from DB iw_polys + ro_polys."""
    if not data.iw_polys:
        return
    tbl_left = data.tb_left

    def _poly_edges(poly):
        edges = []
        for i in range(len(poly)):
            j = (i + 1) % len(poly)
            dx = poly[j][0] - poly[i][0]
            dy = poly[j][1] - poly[i][1]
            edges.append(math.sqrt(dx**2 + dy**2))
        return edges

    # Build RO lookup: wall_name → list of "ROx width" strings
    ro_by_wall: dict[str, list[str]] = {}
    for ro in (data.ro_polys or []):
        poly = ro.get("poly")
        if poly and len(poly) >= 4:
            edges = [math.hypot(poly[(i+1)%4][0]-poly[i][0],
                                poly[(i+1)%4][1]-poly[i][1]) for i in range(4)]
            w_in = max(edges) * 12
        else:
            b = ro.get("bbox") or {}
            orient = ro.get("orientation", "H")
            w_in = (b.get("e", 0) - b.get("w", 0) if orient == "H"
                    else b.get("n", 0) - b.get("s", 0)) * 12
        w_str = f"{w_in:.2f}".rstrip("0").rstrip(".")
        ro_by_wall.setdefault(ro.get("wall_name", ""), []).append(
            f"{ro.get('name','')} {w_str}&#8243;")

    # Build rows: (id, thickness_in, length_ft)
    iw_rows = []
    for name, poly in sorted(data.iw_polys.items()):
        if len(poly) < 3:
            continue
        edges = _poly_edges(poly)
        thick_in = min(edges) * 12
        length_ft = max(edges)
        iw_rows.append((name, thick_in, length_ft))

    iw_tbl_top = tbl_border_bottom + 14
    iw_row_h = 7.5
    iw_col = [tbl_left + 20, tbl_left + 48, tbl_left + 82, tbl_left + 168]

    out.append(f'<text x="{(tbl_left + iw_col[-1]) / 2:.1f}" y="{iw_tbl_top:.1f}"'
               f' text-anchor="middle" font-family="Arial" font-size="7"'
               f' font-weight="bold" fill="{CLR_TITLE}">Interior Walls</text>')

    iw_hdr_y = iw_tbl_top + 10
    iw_hdrs = ["ID", "Thk", "Length", "Openings"]
    iw_hdr_x = [tbl_left + 2, iw_col[1] - 2, iw_col[2] - 2, iw_col[2] + 2]
    iw_hdr_a = ["start", "end", "end", "start"]
    for hx, ha, hd in zip(iw_hdr_x, iw_hdr_a, iw_hdrs):
        out.append(f'<text x="{hx:.1f}" y="{iw_hdr_y:.1f}"'
                   f' text-anchor="{ha}" font-family="Arial" font-size="6"'
                   f' font-weight="bold" fill="{CLR_TITLE}">{hd}</text>')

    iw_line_y = iw_hdr_y + 2.5
    out.append(f'<line x1="{tbl_left:.1f}" y1="{iw_line_y:.1f}"'
               f' x2="{iw_col[-1]:.1f}" y2="{iw_line_y:.1f}"'
               f' stroke="#999" stroke-width="0.5"/>')

    for ri, (iw_id, thick_in, length_ft) in enumerate(iw_rows):
        y = iw_line_y + (ri + 1) * iw_row_h
        total_in = round(length_ft * 12, 1)
        ft = int(total_in) // 12
        remain_in = total_in - ft * 12
        if remain_in >= 12.0:
            ft += 1
            remain_in = 0.0
        length_str = f"{ft}&#8242; {remain_in:.1f}&#8243;"
        ros = ", ".join(ro_by_wall.get(iw_id, []))
        thk_str = f'{thick_in:.1f}&#8243;'
        vals = [iw_id, thk_str, length_str, ros]
        for vx, va, vv in zip(iw_hdr_x, iw_hdr_a, vals):
            out.append(f'<text x="{vx:.1f}" y="{y:.1f}"'
                       f' text-anchor="{va}" font-family="Arial"'
                       f' font-size="6" fill="{CLR_TITLE}">{vv}</text>')

    iw_border_top = iw_tbl_top - 8.5
    iw_border_bottom = iw_line_y + len(iw_rows) * iw_row_h + 3
    out.append(f'<rect x="{tbl_left:.1f}" y="{iw_border_top:.1f}"'
               f' width="{iw_col[-1] - tbl_left:.1f}"'
               f' height="{iw_border_bottom - iw_border_top:.1f}"'
               f' fill="none" stroke="#999" stroke-width="0.5"/>')
    return iw_border_bottom


def _render_interior_walls_table(out, data, tbl_border_bottom):
    """Render the interior walls summary table. Returns table bottom y."""
    if data.layout is None:
        # DB path — build table from iw_polys + ro_polys
        return _render_interior_walls_table_db(out, data, tbl_border_bottom)
    pts = data.pts
    layout = data.layout
    tbl_left = data.tb_left

    rough_openings = compute_rough_openings(pts, layout)
    # Map IW name → list of "RO# width" strings
    ro_by_wall: dict[str, list[str]] = {}
    for ro in rough_openings:
        w_in = ro.width * 12
        w_str = f"{w_in:.2f}".rstrip("0").rstrip(".")
        ro_by_wall.setdefault(ro.wall_name, []).append(f"{ro.name} {w_str}&#8243;")

    def _poly_len(poly):
        """Length of longest edge of a 4-point polygon (the wall length)."""
        edges = []
        for i in range(4):
            j = (i + 1) % 4
            dx = poly[j][0] - poly[i][0]
            dy = poly[j][1] - poly[i][1]
            edges.append(math.sqrt(dx**2 + dy**2))
        return max(edges)

    # (id, thickness_inches, length_ft, vertical)
    iw_rows = [
        ("IW1",  6, _poly_len(layout.iw1.poly), True),
        ("IW2",  6, _poly_len(layout.iw2.poly), True),
        ("IW2o", 6, _poly_len(layout.iw2o.poly), True),
        ("IW2s", 6, _poly_len(layout.iw2s.poly), True),
        ("IW3",  4, _poly_len(layout.iw3.poly), True),
        ("IW4",  4, _poly_len(layout.iw4.poly), True),
        ("IW5",  3, _poly_len(layout.iw5.poly), False),
        ("IW6",  1, _poly_len(layout.iw6.poly), False),
        ("IW7",  4, _poly_len(layout.iw7.poly), False),
        ("IW8",  6, _poly_len(layout.iw8.poly), False),
        ("IW9",  4, _poly_len(layout.iw9.poly), True),
        ("IW11", 4, _poly_len(layout.iw11.poly), True),
        ("IW12", 4, _poly_len(layout.iw12.poly), False),
    ]

    iw_tbl_top = tbl_border_bottom + 14
    iw_row_h = 7.5
    iw_col = [tbl_left + 20, tbl_left + 48, tbl_left + 82, tbl_left + 168]

    # Table title
    out.append(f'<text x="{(tbl_left + iw_col[-1]) / 2:.1f}" y="{iw_tbl_top:.1f}"'
               f' text-anchor="middle" font-family="Arial" font-size="7"'
               f' font-weight="bold" fill="{CLR_TITLE}">Interior Walls</text>')

    # Column headers
    iw_hdr_y = iw_tbl_top + 10
    iw_hdrs = ["ID", "Thk", "Length", "Openings"]
    iw_hdr_x = [tbl_left + 2, iw_col[1] - 2, iw_col[2] - 2, iw_col[2] + 2]
    iw_hdr_a = ["start", "end", "end", "start"]
    for hx, ha, hd in zip(iw_hdr_x, iw_hdr_a, iw_hdrs):
        out.append(f'<text x="{hx:.1f}" y="{iw_hdr_y:.1f}"'
                   f' text-anchor="{ha}" font-family="Arial" font-size="6"'
                   f' font-weight="bold" fill="{CLR_TITLE}">{hd}</text>')

    # Header underline
    iw_line_y = iw_hdr_y + 2.5
    out.append(f'<line x1="{tbl_left:.1f}" y1="{iw_line_y:.1f}"'
               f' x2="{iw_col[-1]:.1f}" y2="{iw_line_y:.1f}"'
               f' stroke="#999" stroke-width="0.5"/>')

    # Data rows
    for ri, (iw_id, thick, length_ft, _vert) in enumerate(iw_rows):
        y = iw_line_y + (ri + 1) * iw_row_h
        total_in = round(length_ft * 12, 1)
        ft = int(total_in) // 12
        remain_in = total_in - ft * 12
        if remain_in >= 12.0:
            ft += 1
            remain_in = 0.0
        length_str = f"{ft}&#8242; {remain_in:.1f}&#8243;"
        ros = ", ".join(ro_by_wall.get(iw_id, []))
        vals = [iw_id, f'{thick}&#8243;', length_str, ros]
        for vx, va, vv in zip(iw_hdr_x, iw_hdr_a, vals):
            out.append(f'<text x="{vx:.1f}" y="{y:.1f}"'
                       f' text-anchor="{va}" font-family="Arial"'
                       f' font-size="6" fill="{CLR_TITLE}">{vv}</text>')

    # Table border
    iw_border_top = iw_tbl_top - 8.5
    iw_border_bottom = iw_line_y + len(iw_rows) * iw_row_h + 3
    out.append(f'<rect x="{tbl_left:.1f}" y="{iw_border_top:.1f}"'
               f' width="{iw_col[-1] - tbl_left:.1f}"'
               f' height="{iw_border_bottom - iw_border_top:.1f}"'
               f' fill="none" stroke="#999" stroke-width="0.5"/>')
    return iw_border_bottom


def _render_room_glazing_table(out, data, tbl_border_bottom):
    """Render room window-to-floor glazing summary table below Interior Walls."""
    rows = data.glazing_rows
    if not rows:
        return

    tbl_left  = data.tb_left
    tbl_right = tbl_left + 168   # match Interior Walls table right edge

    tbl_top = tbl_border_bottom + 14
    row_h   = 7.5

    # Columns: Room | Floor | Glass (with opening IDs) | Ratio
    # Interior Walls right edge = tbl_left + 168
    col_x = [tbl_left + 2,  tbl_left + 83, tbl_left + 148, tbl_left + 166]
    col_a = ["start",        "end",          "end",           "end"]
    hdrs  = ["Room",         "Floor",        "Glass",         "Ratio"]

    out.append(
        f'<text x="{(tbl_left + tbl_right) / 2:.1f}" y="{tbl_top:.1f}"'
        f' text-anchor="middle" font-family="Arial" font-size="7"'
        f' font-weight="bold" fill="{CLR_TITLE}">Room Glazing</text>'
    )

    hdr_y = tbl_top + 10
    for hx, ha, hd in zip(col_x, col_a, hdrs):
        out.append(
            f'<text x="{hx:.1f}" y="{hdr_y:.1f}"'
            f' text-anchor="{ha}" font-family="Arial" font-size="6"'
            f' font-weight="bold" fill="{CLR_TITLE}">{hd}</text>'
        )

    line_y = hdr_y + 2.5
    out.append(
        f'<line x1="{tbl_left:.1f}" y1="{line_y:.1f}"'
        f' x2="{tbl_right:.1f}" y2="{line_y:.1f}"'
        f' stroke="#999" stroke-width="0.5"/>'
    )

    for ri, row in enumerate(rows):
        y = line_y + (ri + 1) * row_h
        label      = row['label']
        floor_sqft = row['floor_sqft']
        glass_sqft = row['glass_sqft']
        pct        = row['pct']

        # Opening IDs with widths (no type label — just IDs)
        op_ids = sorted({
            o['name']
            for w in row['walls'] for o in w['openings']
        })
        glass_str = (
            f"({', '.join(op_ids)})  {glass_sqft:.1f} sf"
            if op_ids else "&#8212;"
        )

        vals = [
            label,
            f"{floor_sqft:.1f} sf",
            glass_str,
            f"{pct:.1f}%",
        ]
        for vx, va, vv in zip(col_x, col_a, vals):
            out.append(
                f'<text x="{vx:.1f}" y="{y:.1f}"'
                f' text-anchor="{va}" font-family="Arial"'
                f' font-size="6" fill="{CLR_TITLE}">{vv}</text>'
            )

    border_top    = tbl_top - 8.5
    border_bottom = line_y + len(rows) * row_h + 3
    out.append(
        f'<rect x="{tbl_left:.1f}" y="{border_top:.1f}"'
        f' width="{tbl_right - tbl_left:.1f}"'
        f' height="{border_bottom - border_top:.1f}"'
        f' fill="none" stroke="#999" stroke-width="0.5"/>'
    )


def _render_f_labels(out, data):
    """Render outline point labels on the exterior side of the outline."""
    pts = data.pts
    to_svg = data.to_svg
    outline_segs = data.outline_segs

    FONT_SIZE = "5"
    COLOR = "#d32f2f"
    DOT_R = "1.0"
    LABEL_DIST = 5.0  # SVG px outward from point

    all_names = list(dict.fromkeys(
        n for s in outline_segs for n in (s.start, s.end) if n in pts))
    if not all_names:
        return

    # Centroid in SVG coords for outward label placement
    _csx, _csy = to_svg(
        sum(pts[n][0] for n in all_names) / len(all_names),
        sum(pts[n][1] for n in all_names) / len(all_names),
    )

    for name in all_names:
        sx, sy = to_svg(*pts[name])
        dx, dy = sx - _csx, sy - _csy
        d = math.sqrt(dx * dx + dy * dy) or 1.0
        lx = sx + dx / d * LABEL_DIST
        ly = sy + dy / d * LABEL_DIST
        anchor = "end" if dx < -1 else ("start" if dx > 1 else "middle")
        out.append(f'<circle cx="{sx:.1f}" cy="{sy:.1f}" r="{DOT_R}"'
                   f' fill="{COLOR}"/>')
        out.append(f'<text x="{lx:.1f}" y="{ly:.1f}"'
                   f' text-anchor="{anchor}" dominant-baseline="central"'
                   f' font-family="Arial" font-size="{FONT_SIZE}"'
                   f' fill="{COLOR}">{name}</text>')


def render_walls_svg(data, *, title="Outer Walls", include_interior=False):
    """Render the wall detail SVG. Returns SVG string."""
    to_svg = data.to_svg

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

    _render_wall_segments(out, data)
    _render_section_outlines(out, data)

    if include_interior:
        _render_interior_walls(out, data)
        _render_opening_dims(out, data)

    _render_f_labels(out, data)
    _render_opening_labels(out, data)
    _render_title_block(out, data)
    tbl_bottom = _render_wall_table(out, data)

    if include_interior:
        iw_bottom = _render_interior_walls_table(out, data, tbl_bottom)
        if iw_bottom is not None:
            _render_room_glazing_table(out, data, iw_bottom)

    # --- Roof outline (dotted) ---
    roof_poly = data.roof_poly
    roof_svg = svg_polygon_pts(roof_poly, to_svg, prec=2)
    out.append(f'<polygon points="{roof_svg}" fill="none"'
               f' stroke="#333" stroke-width="0.6" stroke-dasharray="3,2"/>')

    out.append('</svg>')
    return "\n".join(out)



# ============================================================
# Main entry point
# ============================================================

if __name__ == "__main__":
    _dir = os.path.dirname(os.path.abspath(__file__))
    _db = os.path.join(_dir, "..", "app", "adu.db")
    if os.path.exists(_db):
        from app.database import get_constants_dict, get_outline_chain
        from app.gen_provider import build_generator_data
        _constants = get_constants_dict(_db)
        _chain = get_outline_chain(_db)
        _gd = build_generator_data(_constants, chain_rows=_chain, db_path=_db)
        data = build_wall_data(_gd)
    else:
        data = build_wall_data()

    svg_content = render_walls_svg(data)
    svg_path = os.path.join(_dir, "walls.svg")
    with open(svg_path, "w", encoding="utf-8") as f:
        f.write(svg_content)
    print(f"Wall detail written to {svg_path}")

    all_svg = render_walls_svg(data, title="All Walls", include_interior=True)
    all_path = os.path.join(_dir, "all_walls.svg")
    with open(all_path, "w", encoding="utf-8") as f:
        f.write(all_svg)
    print(f"All walls written to {all_path}")

    print(f"Shell: {SHELL_THICKNESS * 12:.0f}\" / "
          f"Gap: {AIR_GAP * 12:.0f}\" / "
          f"Shell: {SHELL_THICKNESS * 12:.0f}\"")
    print(f"Opening corner inside radius: {OPENING_INSIDE_RADIUS * 12:.1f}\"")
