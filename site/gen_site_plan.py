"""Generate site_plan.pdf: building outline overlaid on site survey."""

import datetime
import math
import sys
import os
from collections import namedtuple

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import fitz  # pymupdf
from floorplan.gen_floorplan import build_floorplan_data
from shared.geometry import vert_isects, path_polygon, segment_polyline, left_norm, off_pt, seg_vecs
from floorplan.constants import ARC_180_R
from shared.types import LineSeg, ArcSeg
from shared.svg import git_describe

# Survey coordinate calibration constants (PDF coords from least-squares fitting)
LINE_TOP = (698.9, 55.2)     # 251.53' meets 216.73' (upper-right corner)
LINE_BOT = (817.9, 557.8)    # 216.73' meets 275.08' (lower-right corner)
BOT_LEFT = (160.0, 561.9)    # 275.08' left endpoint
TL_251 = (108.0, 174.5)      # 251.53' upper-left corner of parcel

# Named corner aliases (intersections of the four labeled property lines)
CORNER_NW = TL_251     # 251.53' meets 163.69'
CORNER_NE = LINE_TOP   # 251.53' meets 216.73'
CORNER_SE = LINE_BOT   # 216.73' meets 275.08'
CORNER_SW = BOT_LEFT   # 275.08' meets 163.69'

# FC (building center) distances from property lines (feet)
FC_DIST_216 = 29.097863567855153  # FC distance from 216.73' line
FC_DIST_275 = 45.786428974476834  # FC distance from 275.08' line

# Existing FRAME & STONE RESIDENCE lower-right corner (14.4'/28.2' wall midline isect)
RESIDENCE_LR = (661.35, 380.80)

COLOR_PROPOSED = (0, 0, 0.6)  # dark blue for proposed/new construction

# Property line vectors (precomputed from corner coordinates)
_LINE_216_DX = LINE_BOT[0] - LINE_TOP[0]
_LINE_216_DY = LINE_BOT[1] - LINE_TOP[1]
_LINE_216_LEN = math.hypot(_LINE_216_DX, _LINE_216_DY)
_LINE_275_DX = LINE_BOT[0] - BOT_LEFT[0]
_LINE_275_DY = LINE_BOT[1] - BOT_LEFT[1]
_LINE_275_LEN = math.hypot(_LINE_275_DX, _LINE_275_DY)

# Rendering constants
DIM_LABEL_FS = 6.0
DIM_LABEL_OFFSET = -3.5
SETBACK_LABEL_FS = 9.0
OUTLINE_STROKE_W = 1.6
BLDG_LABEL_FS = 8.0

SitePlanData = namedtuple("SitePlanData", [
    "pts",              # building points dict
    "building_to_pdf",  # transform function (E,N) → (pdf_x, pdf_y)
    "rotation_deg",     # rotation angle in degrees
    "f15_pdf",          # F15 position in PDF coords
    "ew_dim_ft",        # E-W external dimension (F2→F15 along building E-W)
    "ns_dim_ft",        # N-S external dimension (F18→F6 along building N-S)
    "min_setback_216",  # min perpendicular dist of any F point from 216.73' line (ft)
    "min_setback_275",  # min perpendicular dist of any F point from 275.08' line (ft)
    "draw_poly",        # interpolated building outline (building coords)
    "inner_poly",       # inner wall polygon
    "span_s_pdf",       # N-S span south endpoint in PDF coords
    "span_n_pdf",       # N-S span north endpoint in PDF coords
    "f2_pdf",           # F2 position in PDF coords
    "SCALE",            # PDF pts per foot (2.4)
    "f_series_pdf",     # dict of F-series + FC points in PDF coords
    "residence_dist_ft", # distance from RESIDENCE_LR to closest F point (ft)
    "residence_closest",  # name of closest F point to RESIDENCE_LR
    "p_series_pdf",     # dict of angularly-corrected P-series points in PDF coords
])


def build_site_plan_data():
    """Compute all site plan geometry — no PDF I/O."""
    data = build_floorplan_data()
    pts = data.pts
    outer_poly = data.outer_poly  # F-series outline (350 pts, arcs sampled)
    inner_poly = data.inner_poly  # W-series inner wall (352 pts)

    # --- Survey coordinate calibration ---
    # Scale: 1 inch = 30 ft on the survey; 1 inch = 72 PDF pts
    SCALE = 72.0 / 30.0  # 2.4 PDF pts per foot

    # Direction of 216.73' line in PDF coords (x-right, y-down)
    ldx, ldy, llen = _LINE_216_DX, _LINE_216_DY, _LINE_216_LEN

    # In real-world coords (E=x-right, N=y-up), the line direction is:
    # dE = ldx = +119, dN = -ldy = -502.6
    # Angle of property line from E-axis (real-world):
    prop_angle = math.atan2(-ldy, ldx)  # atan2(-502.6, 119) ≈ -76.7°

    # F16→F17 direction in building coords
    f16 = pts["F16"]
    f17 = pts["F17"]
    f16f17_angle = math.atan2(f17[1] - f16[1], f17[0] - f16[0])

    # Rotation needed: rotate building so F16→F17 is parallel to property line
    rotation = prop_angle - f16f17_angle  # ≈ 73° CCW

    # --- Rotation and placement ---
    cos_r = math.cos(rotation)
    sin_r = math.sin(rotation)

    def rotate_pt(e, n, ce, cn):
        """Rotate point (e,n) by `rotation` around center (ce,cn)."""
        de = e - ce
        dn = n - cn
        return (ce + de * cos_r - dn * sin_r,
                cn + de * sin_r + dn * cos_r)

    # F15 is the reference point for placement
    f15 = pts["F15"]
    fc = pts["FC"]

    # --- Constraint-based placement (2×2 linear system) ---
    # After rotation, each building point has a fixed PDF offset from F15.
    # Unknowns: f15_pdf_x, f15_pdf_y (F15's position on the PDF page).
    # Constraints use FC (building center) distances from property lines.

    # PDF offset of FC from F15 (fixed after rotation)
    off_fc_x = ((fc[0] - f15[0]) * cos_r - (fc[1] - f15[1]) * sin_r) * SCALE
    off_fc_y = -((fc[0] - f15[0]) * sin_r + (fc[1] - f15[1]) * cos_r) * SCALE

    # Constraint A: FC is FC_DIST_216 inside the 216.73' line.
    # Signed distance to left of LINE_TOP→LINE_BOT = property interior.
    a1 = -ldy / llen
    b1 = ldx / llen
    c1 = (FC_DIST_216 * SCALE
          + a1 * (LINE_TOP[0] - off_fc_x) + b1 * (LINE_TOP[1] - off_fc_y))

    # Constraint B: FC is FC_DIST_275 inside the 275.08' line.
    # Direction BOT_LEFT→LINE_BOT; interior is to the right.
    bdx, bdy, blen = _LINE_275_DX, _LINE_275_DY, _LINE_275_LEN
    a2 = bdy / blen
    b2 = -bdx / blen
    c2 = (FC_DIST_275 * SCALE
          + a2 * (BOT_LEFT[0] - off_fc_x) + b2 * (BOT_LEFT[1] - off_fc_y))

    # Solve: a1*fx + b1*fy = c1,  a2*fx + b2*fy = c2
    det = a1 * b2 - a2 * b1
    f15_pdf_x = (c1 * b2 - c2 * b1) / det
    f15_pdf_y = (a1 * c2 - a2 * c1) / det

    def building_to_pdf(e, n):
        """Transform building coords (E,N) → PDF coords (x,y)."""
        re, rn = rotate_pt(e, n, f15[0], f15[1])
        pdf_x = f15_pdf_x + (re - f15[0]) * SCALE
        pdf_y = f15_pdf_y - (rn - f15[1]) * SCALE  # N-up → y-down
        return pdf_x, pdf_y

    # --- Build drawing polygon (outline offset inward by half stroke width) ---
    WALL_T = 8.0 / 12.0
    half_stroke_ft = (OUTLINE_STROKE_W / 2.0) / SCALE
    frac = half_stroke_ft / WALL_T

    n_out = len(outer_poly)
    n_inn = len(inner_poly)
    draw_poly = []
    for i, (oe, on) in enumerate(outer_poly):
        j = round(i * (n_inn - 1) / (n_out - 1))
        ie, inn = inner_poly[j]
        draw_poly.append((oe + frac * (ie - oe), on + frac * (inn - on)))

    # --- External dimensions (building-aligned, matching floorplan dim15/dim13) ---
    _bld_ew, _ = seg_vecs(pts["W9"], pts["W10"])   # building E-W direction
    _bld_ns, _ = seg_vecs(pts["W2"], pts["W5"])     # building N-S direction

    _df_ew = (pts["F15"][0] - pts["F2"][0], pts["F15"][1] - pts["F2"][1])
    ew_dim_ft = abs(_df_ew[0] * _bld_ew[0] + _df_ew[1] * _bld_ew[1])

    _df_ns = (pts["F6"][0] - pts["F18"][0], pts["F6"][1] - pts["F18"][1])
    ns_dim_ft = abs(_df_ns[0] * _bld_ns[0] + _df_ns[1] * _bld_ns[1])

    # --- N-S Interior Max Span (dimension line position only) ---
    _inch = 1.0 / 12.0
    _e_min = min(p[0] for p in inner_poly)
    _e_max = max(p[0] for p in inner_poly)
    _best_span, _best_e, _best_s, _best_n = 0, 0, 0, 0
    _e = _e_min
    while _e <= _e_max + 1e-9:
        _ns = vert_isects(inner_poly, _e)
        if len(_ns) >= 2:
            _s, _n = min(_ns), max(_ns)
            if _n - _s > _best_span:
                _best_span, _best_e, _best_s, _best_n = _n - _s, _e, _s, _n
        _e += _inch

    span_s_pdf = building_to_pdf(_best_e, _best_s)
    span_n_pdf = building_to_pdf(_best_e, _best_n)

    f2_pdf = building_to_pdf(*pts["F2"])
    f15_pdf = (f15_pdf_x, f15_pdf_y)

    # --- F-series PDF coordinates ---
    _f_names = [f"F{i}" for i in range(1, 21) if i not in (3, 4)] + ["F11a", "F11b", "FC"]
    f_series_pdf = {name: building_to_pdf(*pts[name]) for name in _f_names}

    # --- Min setback distances (F-points only, excluding FC) ---
    _f_struct = [f"F{i}" for i in range(1, 21) if i not in (3, 4)] + ["F11a", "F11b"]
    min_setback_216 = min(
        ((pt[0] - LINE_TOP[0]) * (-ldy) + (pt[1] - LINE_TOP[1]) * ldx)
        / (llen * SCALE)
        for pt in (f_series_pdf[n] for n in _f_struct))
    min_setback_275 = min(
        ((pt[0] - BOT_LEFT[0]) * bdy - (pt[1] - BOT_LEFT[1]) * bdx)
        / (blen * SCALE)
        for pt in (f_series_pdf[n] for n in _f_struct))

    # --- Distance from existing residence corner to closest F point ---
    _res_best_name, _res_best_dist = None, float("inf")
    for n in _f_struct:
        pt = f_series_pdf[n]
        d = math.hypot(pt[0] - RESIDENCE_LR[0], pt[1] - RESIDENCE_LR[1])
        if d < _res_best_dist:
            _res_best_dist, _res_best_name = d, n
    residence_dist_ft = _res_best_dist / SCALE

    # --- P-series PDF coordinates (angularly corrected, path_area.svg pivot) ---
    _corrected_bld = _correct_p_series(pts)
    p_series_pdf = {n: building_to_pdf(*p) for n, p in _corrected_bld.items()}

    return SitePlanData(
        pts=pts,
        building_to_pdf=building_to_pdf,
        rotation_deg=math.degrees(rotation),
        f15_pdf=f15_pdf,
        ew_dim_ft=ew_dim_ft,
        ns_dim_ft=ns_dim_ft,
        min_setback_216=min_setback_216,
        min_setback_275=min_setback_275,
        draw_poly=draw_poly,
        inner_poly=inner_poly,
        span_s_pdf=span_s_pdf,
        span_n_pdf=span_n_pdf,
        f2_pdf=f2_pdf,
        SCALE=SCALE,
        f_series_pdf=f_series_pdf,
        residence_dist_ft=residence_dist_ft,
        residence_closest=_res_best_name,
        p_series_pdf=p_series_pdf,
    )


def _draw_dim_line(shape, page, pt1, pt2, label, color,
                   fs=DIM_LABEL_FS, offset=DIM_LABEL_OFFSET):
    """Draw a dimension line between pt1 and pt2 with centered rotated label."""
    shape.draw_line(fitz.Point(*pt1), fitz.Point(*pt2))
    shape.finish(color=color, width=0.3)
    dx = pt2[0] - pt1[0]
    dy = pt2[1] - pt1[1]
    length = math.hypot(dx, dy)
    deg = math.degrees(math.atan2(dy, dx))
    mid_x = (pt1[0] + pt2[0]) / 2.0 + offset * dy / length
    mid_y = (pt1[1] + pt2[1]) / 2.0 - offset * dx / length
    tw = fitz.get_text_length(label, fontname="helv", fontsize=fs)
    page.insert_text(
        fitz.Point(mid_x - tw / 2.0, mid_y + fs / 3.0),
        label, fontname="helv", fontsize=fs, color=color,
        morph=(fitz.Point(mid_x, mid_y), fitz.Matrix(-deg - 180)))


def _draw_setback_label(page, pt_pdf, line_p1, line_p2, value_ft, color,
                        fs=SETBACK_LABEL_FS):
    """Draw a label midway between pt and its perpendicular projection onto a line."""
    ldx = line_p2[0] - line_p1[0]
    ldy = line_p2[1] - line_p1[1]
    llen_sq = ldx * ldx + ldy * ldy
    t_proj = ((pt_pdf[0] - line_p1[0]) * ldx +
              (pt_pdf[1] - line_p1[1]) * ldy) / llen_sq
    proj_x = line_p1[0] + t_proj * ldx
    proj_y = line_p1[1] + t_proj * ldy
    cap_x = (pt_pdf[0] + proj_x) / 2.0
    cap_y = (pt_pdf[1] + proj_y) / 2.0
    perp_deg = math.degrees(math.atan2(proj_y - pt_pdf[1], proj_x - pt_pdf[0]))
    text = f"{value_ft:.1f}'"
    tw = fitz.get_text_length(text, fontname="helv", fontsize=fs)
    page.insert_text(
        fitz.Point(cap_x - tw / 2.0, cap_y + fs / 3.0),
        text, fontname="helv", fontsize=fs, color=color,
        morph=(fitz.Point(cap_x, cap_y), fitz.Matrix(-perp_deg)))


def render_site_plan(sp, corners=True):
    """Render base site plan PDF overlay. Returns unsaved fitz.Document."""
    src = fitz.open("site/site_survey.pdf")
    doc = fitz.open()
    doc.insert_pdf(src)
    src.close()
    page = doc[0]

    SCALE = sp.SCALE
    pts = sp.pts
    building_to_pdf = sp.building_to_pdf

    # --- Draw building outline ---
    shape = page.new_shape()
    x0, y0 = building_to_pdf(*sp.draw_poly[0])
    shape.draw_line(fitz.Point(x0, y0), fitz.Point(x0, y0))
    for pt in sp.draw_poly[1:]:
        x1, y1 = building_to_pdf(*pt)
        shape.draw_line(fitz.Point(x0, y0), fitz.Point(x1, y1))
        x0, y0 = x1, y1
    x1, y1 = building_to_pdf(*sp.draw_poly[0])
    shape.draw_line(fitz.Point(x0, y0), fitz.Point(x1, y1))
    shape.finish(color=COLOR_PROPOSED, width=OUTLINE_STROKE_W)

    # --- Parcel corner circles (2' radius) ---
    if corners:
        _corner_r = 2.0 * SCALE
        for cx, cy in (CORNER_NW, CORNER_NE, CORNER_SE, CORNER_SW):
            shape.draw_circle(fitz.Point(cx, cy), _corner_r)
        shape.finish(color=(1, 0, 0), width=0.5, fill=None,
                     stroke_opacity=0.4)

    # --- F15 to F2-F5 dimension line ---
    f15 = pts["F15"]
    f15_pdf = building_to_pdf(*f15)
    foot_pdf = building_to_pdf(pts["F2"][0], f15[1])
    _draw_dim_line(shape, page, f15_pdf, foot_pdf,
                   f"{sp.ew_dim_ft:.1f}'", COLOR_PROPOSED)

    # --- N-S Interior Max Span dimension line ---
    _draw_dim_line(shape, page, sp.span_s_pdf, sp.span_n_pdf,
                   f"{sp.ns_dim_ft:.1f}'", COLOR_PROPOSED)

    # --- "PROPOSED CONC. GUEST HOUSE" label ---
    _cx = sum(p[0] for p in sp.inner_poly) / len(sp.inner_poly)
    _cy = sum(p[1] for p in sp.inner_poly) / len(sp.inner_poly)
    label_pdf_raw = building_to_pdf(_cx, _cy + 2.0)
    label_pdf = (label_pdf_raw[0], label_pdf_raw[1] + 5.0 * SCALE)
    label_lines = ["     PROPOSED", "CONC.", "GUEST", "HOUSE"]
    label_lh = BLDG_LABEL_FS * 1.15
    block_h = label_lh * len(label_lines)
    start_y = label_pdf[1] - block_h / 2.0 + BLDG_LABEL_FS
    for i, line in enumerate(label_lines):
        lw = fitz.get_text_length(line, fontname="helv", fontsize=BLDG_LABEL_FS)
        page.insert_text(
            fitz.Point(label_pdf[0] - lw / 2.0, start_y + i * label_lh),
            line, fontname="helv", fontsize=BLDG_LABEL_FS, color=COLOR_PROPOSED)

    # --- 11.5' setback caption (from 216.73' line) ---
    f16_pdf = building_to_pdf(*pts["F16"])
    f17_pdf = building_to_pdf(*pts["F17"])
    mid_f16f17 = ((f16_pdf[0] + f17_pdf[0]) / 2.0,
                  (f16_pdf[1] + f17_pdf[1]) / 2.0)
    _draw_setback_label(page, mid_f16f17, LINE_TOP, LINE_BOT,
                        sp.min_setback_216, COLOR_PROPOSED)

    # --- Min setback from 275.08' line caption ---
    _draw_setback_label(page, sp.f2_pdf, BOT_LEFT, LINE_BOT,
                        sp.min_setback_275, COLOR_PROPOSED)

    # --- Distance from residence corner to closest F point ---
    _res_pt = sp.f_series_pdf[sp.residence_closest]
    _res_mid_x = (RESIDENCE_LR[0] + _res_pt[0]) / 2.0
    _res_mid_y = (RESIDENCE_LR[1] + _res_pt[1]) / 2.0
    _res_deg = math.degrees(math.atan2(
        _res_pt[1] - RESIDENCE_LR[1], _res_pt[0] - RESIDENCE_LR[0]))
    _res_text = f"{sp.residence_dist_ft:.1f}'"
    _res_fs = 9.0
    _res_tw = fitz.get_text_length(_res_text, fontname="helv", fontsize=_res_fs)
    page.insert_text(
        fitz.Point(_res_mid_x - _res_tw / 2.0, _res_mid_y + _res_fs / 3.0),
        _res_text, fontname="helv", fontsize=_res_fs, color=COLOR_PROPOSED,
        morph=(fitz.Point(_res_mid_x, _res_mid_y), fitz.Matrix(-_res_deg)))

    # --- "FRONT ↑" annotation above 251.53' line ---
    _mid_251_x = (TL_251[0] + LINE_TOP[0]) / 2.0
    _mid_251_y = (TL_251[1] + LINE_TOP[1]) / 2.0

    front_fs = 16.8
    front_gap = 4.0 + 4.0 * SCALE

    front_text = "FRONT"
    front_tw = fitz.get_text_length(front_text, fontname="helv", fontsize=front_fs)

    _arr_space = front_fs * 0.3
    _arr_h = front_fs * 0.65
    _arr_hw = front_fs * 0.15
    _arr_ah = front_fs * 0.2

    _total_w = front_tw + _arr_space

    _front_x = _mid_251_x - _total_w / 2.0
    _front_y = _mid_251_y - front_gap

    page.insert_text(
        fitz.Point(_front_x, _front_y),
        front_text, fontname="helv", fontsize=front_fs, color=(0, 0, 0))

    # Draw up-arrow next to "FRONT"
    _arr_cx = _front_x + front_tw + _arr_space
    _arr_bot = _front_y
    _arr_top = _front_y - _arr_h

    shape.draw_line(fitz.Point(_arr_cx, _arr_bot), fitz.Point(_arr_cx, _arr_top))
    shape.draw_line(fitz.Point(_arr_cx, _arr_top),
                    fitz.Point(_arr_cx - _arr_hw, _arr_top + _arr_ah))
    shape.draw_line(fitz.Point(_arr_cx, _arr_top),
                    fitz.Point(_arr_cx + _arr_hw, _arr_top + _arr_ah))
    shape.finish(color=(0, 0, 0), width=0.8)

    # --- Git describe / timestamp caption (right of GRAPHIC SCALE) ---
    _cap_x = 680.0
    _cap_y = 600.0
    _cap_fs = 5.5
    _cap_now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    _cap_git = git_describe()
    _cap_text = f"Generated {_cap_now}  {_cap_git}"
    _cap_tw = fitz.get_text_length(_cap_text, fontname="helv", fontsize=_cap_fs)
    page.insert_text(
        fitz.Point(_cap_x - _cap_tw / 2.0, _cap_y),
        _cap_text, fontname="helv", fontsize=_cap_fs, color=(0.4, 0.4, 0.4))

    shape.commit()
    return doc


def render_site_plan_df(doc, sp):
    """Add drainfield annotations to an existing site plan document."""
    page = doc[0]
    SCALE = sp.SCALE
    building_to_pdf = sp.building_to_pdf

    # Drainfield: dashed rounded rectangle, 30'W x 8'H, 2' corner radius
    _df_line_x = 661.5
    _df_line_y = 316.0

    _df_w = 25.0 * SCALE
    _df_h = 10.0 * SCALE
    _df_r = 2.0 * SCALE
    _df_left = _df_line_x + 8.0 * SCALE
    _df_top = _df_line_y - _df_h / 2.0
    _df_right = _df_left + _df_w
    _df_bot = _df_line_y + _df_h / 2.0

    _df_fs = 7.2
    _df_lh = _df_fs * 1.15
    _df_lines_existing = ["EXISTING", "DRAINFIELD"]
    _n_arc = 8

    def _draw_drainfield(left, top, right, bot, r, lines=None, angle_deg=0,
                         color=(0, 0, 0)):
        if lines is None:
            lines = _df_lines_existing
        cx = (left + right) / 2.0
        cy = (top + bot) / 2.0
        hw = (right - left) / 2.0
        hh = (bot - top) / 2.0

        rel = []
        rel.append((-hw + r, -hh))
        rel.append((hw - r, -hh))
        for i in range(1, _n_arc + 1):
            a = -math.pi / 2 + math.pi / 2 * i / _n_arc
            rel.append((hw - r + r * math.cos(a), -hh + r + r * math.sin(a)))
        rel.append((hw, hh - r))
        for i in range(1, _n_arc + 1):
            a = math.pi / 2 * i / _n_arc
            rel.append((hw - r + r * math.cos(a), hh - r + r * math.sin(a)))
        rel.append((-hw + r, hh))
        for i in range(1, _n_arc + 1):
            a = math.pi / 2 + math.pi / 2 * i / _n_arc
            rel.append((-hw + r + r * math.cos(a), hh - r + r * math.sin(a)))
        rel.append((-hw, -hh + r))
        for i in range(1, _n_arc + 1):
            a = math.pi + math.pi / 2 * i / _n_arc
            rel.append((-hw + r + r * math.cos(a), -hh + r + r * math.sin(a)))

        rad = math.radians(angle_deg)
        cos_r = math.cos(rad)
        sin_r = math.sin(rad)
        pts = []
        for rx, ry in rel:
            pts.append(fitz.Point(cx + rx * cos_r - ry * sin_r,
                                  cy + rx * sin_r + ry * cos_r))

        page.draw_polyline(pts, color=color, width=0.8,
                           dashes="[4 3] 0", closePath=True)

        block_h = _df_lh * len(lines)
        start_y = cy - block_h / 2.0 + _df_fs
        morph = None
        if angle_deg != 0:
            morph = (fitz.Point(cx, cy), fitz.Matrix(-angle_deg))
        for li, lt in enumerate(lines):
            tw = fitz.get_text_length(lt, fontname="helv", fontsize=_df_fs)
            page.insert_text(
                fitz.Point(cx - tw / 2.0, start_y + li * _df_lh),
                lt, fontname="helv", fontsize=_df_fs, color=color,
                morph=morph)

    # Right drainfield
    _draw_drainfield(_df_left, _df_top, _df_right, _df_bot, _df_r)

    # Left drainfield
    _res_left = 534.0
    _res_bl_y = 348.0
    _df2_right = _res_left - 8.0 * SCALE
    _df2_left = _df2_right - _df_w
    _df2_cy = _res_bl_y - 1.0 * SCALE
    _df2_top = _df2_cy - _df_h / 2.0
    _df2_bot = _df2_cy + _df_h / 2.0
    _draw_drainfield(_df2_left, _df2_top, _df2_right, _df2_bot, _df_r)

    # New drainfield: midway between F2 (PDF) and 275.08'/216.73' corner
    f2_pdf = sp.f2_pdf
    _ndf_dx = LINE_BOT[0] - f2_pdf[0]
    _ndf_dy = LINE_BOT[1] - f2_pdf[1]
    _ndf_len = math.hypot(_ndf_dx, _ndf_dy)
    _ndf_cx = (f2_pdf[0] + LINE_BOT[0]) / 2.0 - 4.0 * SCALE * _ndf_dx / _ndf_len
    _ndf_cy = (f2_pdf[1] + LINE_BOT[1]) / 2.0 - 4.0 * SCALE * _ndf_dy / _ndf_len
    _ndf_angle = math.degrees(math.atan2(_ndf_dy, _ndf_dx))
    _ndf_left = _ndf_cx - _df_w / 2.0
    _ndf_top = _ndf_cy - _df_h / 2.0
    _ndf_right = _ndf_cx + _df_w / 2.0
    _ndf_bot = _ndf_cy + _df_h / 2.0
    _draw_drainfield(_ndf_left, _ndf_top, _ndf_right, _ndf_bot, _df_r,
                     lines=["NEW", "DRAINFIELD"], angle_deg=_ndf_angle,
                     color=COLOR_PROPOSED)

    return doc


COLOR_OUTER_PATH = (0, 0.4, 0)  # dark green for survey outer path


def _correct_p_series(pts):
    """Correct P-series points by rotation about the path_area.svg pivot.

    The rotation angle aligns PX→P5 with F17→F16 (both parallel to the
    216.73' property line, field-verified).  The pivot matches
    gen_path_svg.py: projection of the C15 arc center onto the PiX→Pi5
    line, using F15 easting to locate C15.

    Returns dict of corrected building-coord points for the 13 P-series names.
    """
    pxp5_angle = math.atan2(pts["P5"][1] - pts["PX"][1],
                            pts["P5"][0] - pts["PX"][0])
    f17f16_angle = math.atan2(pts["F16"][1] - pts["F17"][1],
                              pts["F16"][0] - pts["F17"][0])
    corr = f17f16_angle - pxp5_angle
    cos_c = math.cos(corr)
    sin_c = math.sin(corr)

    # Pivot: same as gen_path_svg.py — project C15 arc-center estimate onto
    # the PiX→Pi5 line, where C15.E = F15.E - ARC_180_R.
    d_pip = (pts["Pi5"][0] - pts["PiX"][0], pts["Pi5"][1] - pts["PiX"][1])
    L_pip = math.hypot(*d_pip)
    d_pip_u = (d_pip[0] / L_pip, d_pip[1] / L_pip)
    ln_pip = left_norm(pts["PiX"], pts["Pi5"])
    o_pip = off_pt(pts["PiX"], ln_pip, ARC_180_R)
    f15_e = pts["F15"][0]
    t_cf = (f15_e - ARC_180_R - o_pip[0]) / d_pip_u[0]
    cf = (f15_e - ARC_180_R, o_pip[1] + t_cf * d_pip_u[1])
    t_piv = ((cf[0] - pts["PiX"][0]) * d_pip[0]
             + (cf[1] - pts["PiX"][1]) * d_pip[1]) / (d_pip[0] ** 2 + d_pip[1] ** 2)
    a = pts["PiX"][0] + t_piv * d_pip[0]
    b = pts["PiX"][1] + t_piv * d_pip[1]

    names = ["POB", "P2", "P3", "P4", "P5",
             "T1", "T2", "T3", "PA", "PX", "TC1", "TC2", "TC3"]
    return {name: (a + (pts[name][0] - a) * cos_c - (pts[name][1] - b) * sin_c,
                   b + (pts[name][0] - a) * sin_c + (pts[name][1] - b) * cos_c)
            for name in names}


def render_site_plan_fs(doc, sp):
    """Add outer survey path (P-series traverse) overlay to an existing site plan."""
    page = doc[0]
    building_to_pdf = sp.building_to_pdf

    # Outer path segments (same definition as survey/gen_path_svg.py)
    R1, R2, R3 = 10.0, 12.5, 11.0
    outer_segs = [
        LineSeg("POB", "P2"), LineSeg("P2", "P3"), LineSeg("P3", "T3"),
        ArcSeg("T3", "PX", "TC3", R3, "CW", 60),
        LineSeg("PX", "P4"), LineSeg("P4", "P5"), LineSeg("P5", "T1"),
        ArcSeg("T1", "PA", "TC1", R1, "CW", 60),
        ArcSeg("PA", "T2", "TC2", R2, "CW", 60),
        LineSeg("T2", "POB"),
    ]

    # Corrected building coords → polygon → PDF coords
    corrected_pts = dict(sp.pts)
    corrected_pts.update(_correct_p_series(sp.pts))
    poly_bld = path_polygon(outer_segs, corrected_pts)
    poly_pdf = [fitz.Point(*building_to_pdf(*p)) for p in poly_bld]

    # Draw closed polyline: 40% opaque dark green
    page.draw_polyline(poly_pdf, color=COLOR_OUTER_PATH, width=1.2,
                       closePath=True, stroke_opacity=0.4)

    return doc


def main():
    sp = build_site_plan_data()

    # --- Base site_plan.pdf (includes corner circles) ---
    doc = render_site_plan(sp)
    out_path = os.path.join(os.path.dirname(__file__), "site_plan.pdf")
    doc.save(out_path)
    doc.close()
    print(f"Written to {out_path}")

    # --- site_plan_df.pdf (fresh render without corner circles) ---
    doc_df = render_site_plan(sp, corners=False)
    render_site_plan_df(doc_df, sp)
    df_path = os.path.join(os.path.dirname(__file__), "site_plan_df.pdf")
    doc_df.save(df_path)
    doc_df.close()
    print(f"Written to {df_path}")

    # --- site_plan_fs.pdf (df content + outer survey path) ---
    doc_fs = render_site_plan(sp, corners=False)
    render_site_plan_df(doc_fs, sp)
    render_site_plan_fs(doc_fs, sp)
    fs_path = os.path.join(os.path.dirname(__file__), "site_plan_fs.pdf")
    doc_fs.save(fs_path)
    doc_fs.close()
    print(f"Written to {fs_path}")

    print(f"Rotation: {sp.rotation_deg:.1f}\u00b0 CCW")
    print(f"F15 PDF position: ({sp.f15_pdf[0]:.1f}, {sp.f15_pdf[1]:.1f})")


if __name__ == "__main__":
    main()
