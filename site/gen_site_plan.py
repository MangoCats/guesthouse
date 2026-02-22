"""Generate site_plan.pdf: building outline overlaid on site survey."""

import datetime
import math
import sys
import os
from collections import namedtuple

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import fitz  # pymupdf
from floorplan.gen_floorplan import build_floorplan_data
from shared.geometry import vert_isects
from shared.svg import git_describe

# Survey coordinate calibration constants (PDF coords from least-squares fitting)
LINE_TOP = (698.9, 55.2)     # 251.53' meets 216.73' (upper-right corner)
LINE_BOT = (817.9, 557.8)    # 216.73' meets 275.08' (lower-right corner)
BOT_LEFT = (160.0, 561.9)    # 275.08' left endpoint
TL_251 = (108.0, 174.5)      # 251.53' upper-left corner of parcel

SitePlanData = namedtuple("SitePlanData", [
    "pts",              # building points dict
    "building_to_pdf",  # transform function (E,N) → (pdf_x, pdf_y)
    "rotation_deg",     # rotation angle in degrees
    "f15_pdf",          # F15 position in PDF coords
    "ew_dim_ft",        # E-W dimension (F15.E - F2.E)
    "ns_dim_ft",        # N-S dimension (F6.N - F18.N)
    "f2_275_dist_ft",   # perpendicular distance from F2 to 275.08' line
    "draw_poly",        # interpolated building outline (building coords)
    "inner_poly",       # inner wall polygon
    "span_s_pdf",       # N-S span south endpoint in PDF coords
    "span_n_pdf",       # N-S span north endpoint in PDF coords
    "f2_pdf",           # F2 position in PDF coords
    "SCALE",            # PDF pts per foot (2.4)
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
    ldx = LINE_BOT[0] - LINE_TOP[0]  # +119 (slants left going up)
    ldy = LINE_BOT[1] - LINE_TOP[1]  # +502.6
    llen = math.hypot(ldx, ldy)

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

    # Inward perpendicular from 216.73' line (into property = roughly west)
    line_de = ldx / llen
    line_dn = -ldy / llen
    # Right perp of (dE, dN) = (dN, -dE) → points into property (west)
    inward_e = line_dn
    inward_n = -line_de

    # Target: F15 placed near the bottom PATIO level on 216.73' line, 11' inside
    patio_y = 435.0 - 12.15 * SCALE * (ldy / llen)
    t_patio = (patio_y - LINE_TOP[1]) / (LINE_BOT[1] - LINE_TOP[1])
    line_at_patio_x = LINE_TOP[0] + t_patio * (LINE_BOT[0] - LINE_TOP[0])

    # F15 target in PDF coords: 11.5' inward from property line
    setback = 11.5  # feet
    f15_pdf_x = line_at_patio_x + setback * SCALE * inward_e
    f15_pdf_y = patio_y + setback * SCALE * (-inward_n)  # y-flip for PDF

    def building_to_pdf(e, n):
        """Transform building coords (E,N) → PDF coords (x,y)."""
        re, rn = rotate_pt(e, n, f15[0], f15[1])
        pdf_x = f15_pdf_x + (re - f15[0]) * SCALE
        pdf_y = f15_pdf_y - (rn - f15[1]) * SCALE  # N-up → y-down
        return pdf_x, pdf_y

    # --- Build drawing polygon (outline offset inward by half stroke width) ---
    STROKE_W = 1.6
    WALL_T = 8.0 / 12.0
    half_stroke_ft = (STROKE_W / 2.0) / SCALE
    frac = half_stroke_ft / WALL_T

    n_out = len(outer_poly)
    n_inn = len(inner_poly)
    draw_poly = []
    for i, (oe, on) in enumerate(outer_poly):
        j = round(i * (n_inn - 1) / (n_out - 1))
        ie, inn = inner_poly[j]
        draw_poly.append((oe + frac * (ie - oe), on + frac * (inn - on)))

    # --- Dimension values ---
    ew_dim_ft = pts["F15"][0] - pts["F2"][0]
    ns_dim_ft = pts["F6"][1] - pts["F18"][1]

    # --- N-S Interior Max Span ---
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

    # --- F2 distance to 275.08' line ---
    bdx = LINE_BOT[0] - BOT_LEFT[0]
    bdy = LINE_BOT[1] - BOT_LEFT[1]
    blen = math.hypot(bdx, bdy)

    f2_pdf = building_to_pdf(*pts["F2"])

    t_f2 = ((f2_pdf[0] - BOT_LEFT[0]) * bdx +
            (f2_pdf[1] - BOT_LEFT[1]) * bdy) / (blen * blen)
    proj_f2_x = BOT_LEFT[0] + t_f2 * bdx
    proj_f2_y = BOT_LEFT[1] + t_f2 * bdy

    dist_pts = math.hypot(f2_pdf[0] - proj_f2_x, f2_pdf[1] - proj_f2_y)
    f2_275_dist_ft = dist_pts / SCALE

    f15_pdf = (f15_pdf_x, f15_pdf_y)

    return SitePlanData(
        pts=pts,
        building_to_pdf=building_to_pdf,
        rotation_deg=math.degrees(rotation),
        f15_pdf=f15_pdf,
        ew_dim_ft=ew_dim_ft,
        ns_dim_ft=ns_dim_ft,
        f2_275_dist_ft=f2_275_dist_ft,
        draw_poly=draw_poly,
        inner_poly=inner_poly,
        span_s_pdf=span_s_pdf,
        span_n_pdf=span_n_pdf,
        f2_pdf=f2_pdf,
        SCALE=SCALE,
    )


def render_site_plan(sp):
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
    STROKE_W = 1.6
    shape = page.new_shape()
    x0, y0 = building_to_pdf(*sp.draw_poly[0])
    shape.draw_line(fitz.Point(x0, y0), fitz.Point(x0, y0))
    for pt in sp.draw_poly[1:]:
        x1, y1 = building_to_pdf(*pt)
        shape.draw_line(fitz.Point(x0, y0), fitz.Point(x1, y1))
        x0, y0 = x1, y1
    x1, y1 = building_to_pdf(*sp.draw_poly[0])
    shape.draw_line(fitz.Point(x0, y0), fitz.Point(x1, y1))
    shape.finish(color=(0, 0, 0), width=STROKE_W)

    # --- F15 to F2-F3 dimension line ---
    f15 = pts["F15"]
    f15_pdf = building_to_pdf(*f15)
    foot_pdf = building_to_pdf(pts["F2"][0], f15[1])
    shape.draw_line(fitz.Point(*f15_pdf), fitz.Point(*foot_pdf))
    shape.finish(color=(0, 0, 0), width=0.3)

    # Caption "36.0'"
    dim_dx = foot_pdf[0] - f15_pdf[0]
    dim_dy = foot_pdf[1] - f15_pdf[1]
    dim_len = math.hypot(dim_dx, dim_dy)
    dim_deg = math.degrees(math.atan2(dim_dy, dim_dx))
    dim_shift = -3.5
    dim_mid_x = (f15_pdf[0] + foot_pdf[0]) / 2.0 + dim_shift * dim_dy / dim_len
    dim_mid_y = (f15_pdf[1] + foot_pdf[1]) / 2.0 - dim_shift * dim_dx / dim_len
    dim_text = f"{sp.ew_dim_ft:.1f}'"
    dim_fs = 6.0
    dim_tw = fitz.get_text_length(dim_text, fontname="helv", fontsize=dim_fs)
    page.insert_text(
        fitz.Point(dim_mid_x - dim_tw / 2.0, dim_mid_y + dim_fs / 3.0),
        dim_text, fontname="helv", fontsize=dim_fs, color=(0, 0, 0),
        morph=(fitz.Point(dim_mid_x, dim_mid_y), fitz.Matrix(-dim_deg - 180)))

    # --- N-S Interior Max Span dimension line ---
    shape.draw_line(fitz.Point(*sp.span_s_pdf), fitz.Point(*sp.span_n_pdf))
    shape.finish(color=(0, 0, 0), width=0.3)

    # Caption
    ns_dx = sp.span_n_pdf[0] - sp.span_s_pdf[0]
    ns_dy = sp.span_n_pdf[1] - sp.span_s_pdf[1]
    ns_len = math.hypot(ns_dx, ns_dy)
    ns_deg = math.degrees(math.atan2(ns_dy, ns_dx))
    ns_shift = -3.5
    ns_mid_x = (sp.span_s_pdf[0] + sp.span_n_pdf[0]) / 2.0 + ns_shift * ns_dy / ns_len
    ns_mid_y = (sp.span_s_pdf[1] + sp.span_n_pdf[1]) / 2.0 - ns_shift * ns_dx / ns_len
    ns_text = f"{sp.ns_dim_ft:.1f}'"
    ns_fs = 6.0
    ns_tw = fitz.get_text_length(ns_text, fontname="helv", fontsize=ns_fs)
    page.insert_text(
        fitz.Point(ns_mid_x - ns_tw / 2.0, ns_mid_y + ns_fs / 3.0),
        ns_text, fontname="helv", fontsize=ns_fs, color=(0, 0, 0),
        morph=(fitz.Point(ns_mid_x, ns_mid_y), fitz.Matrix(-ns_deg - 180)))

    # --- "PROPOSED CONC. GUEST HOUSE" label ---
    _cx = sum(p[0] for p in sp.inner_poly) / len(sp.inner_poly)
    _cy = sum(p[1] for p in sp.inner_poly) / len(sp.inner_poly)
    label_pdf_raw = building_to_pdf(_cx, _cy + 2.0)
    label_pdf = (label_pdf_raw[0], label_pdf_raw[1] + 5.0 * SCALE)
    label_fs = 8.0
    label_lines = ["     PROPOSED", "CONC.", "GUEST", "HOUSE"]
    label_lh = label_fs * 1.15
    block_h = label_lh * len(label_lines)
    start_y = label_pdf[1] - block_h / 2.0 + label_fs
    for i, line in enumerate(label_lines):
        lw = fitz.get_text_length(line, fontname="helv", fontsize=label_fs)
        page.insert_text(
            fitz.Point(label_pdf[0] - lw / 2.0, start_y + i * label_lh),
            line, fontname="helv", fontsize=label_fs, color=(0, 0, 0))

    # --- 11.5' setback caption ---
    f16 = pts["F16"]
    f17 = pts["F17"]
    f16_pdf = building_to_pdf(*f16)
    f17_pdf = building_to_pdf(*f17)
    mid_f16f17_x = (f16_pdf[0] + f17_pdf[0]) / 2.0
    mid_f16f17_y = (f16_pdf[1] + f17_pdf[1]) / 2.0

    ldx = LINE_BOT[0] - LINE_TOP[0]
    ldy = LINE_BOT[1] - LINE_TOP[1]
    llen = math.hypot(ldx, ldy)

    t_proj = ((mid_f16f17_x - LINE_TOP[0]) * ldx +
              (mid_f16f17_y - LINE_TOP[1]) * ldy) / (llen * llen)
    proj_x = LINE_TOP[0] + t_proj * ldx
    proj_y = LINE_TOP[1] + t_proj * ldy

    cap_x = (mid_f16f17_x + proj_x) / 2.0
    cap_y = (mid_f16f17_y + proj_y) / 2.0

    perp_deg = math.degrees(math.atan2(proj_y - mid_f16f17_y,
                                       proj_x - mid_f16f17_x))

    text = "11.5'"
    fs = 9.0
    tw = fitz.get_text_length(text, fontname="helv", fontsize=fs)
    page.insert_text(
        fitz.Point(cap_x - tw / 2.0, cap_y + fs / 3.0),
        text, fontname="helv", fontsize=fs, color=(0, 0, 0),
        morph=(fitz.Point(cap_x, cap_y), fitz.Matrix(-perp_deg)))

    # --- F2 distance to 275.08' line caption ---
    bdx = LINE_BOT[0] - BOT_LEFT[0]
    bdy = LINE_BOT[1] - BOT_LEFT[1]
    blen = math.hypot(bdx, bdy)

    t_f2 = ((sp.f2_pdf[0] - BOT_LEFT[0]) * bdx +
            (sp.f2_pdf[1] - BOT_LEFT[1]) * bdy) / (blen * blen)
    proj_f2_x = BOT_LEFT[0] + t_f2 * bdx
    proj_f2_y = BOT_LEFT[1] + t_f2 * bdy

    cap2_x = (sp.f2_pdf[0] + proj_f2_x) / 2.0
    cap2_y = (sp.f2_pdf[1] + proj_f2_y) / 2.0

    perp2_deg = math.degrees(math.atan2(proj_f2_y - sp.f2_pdf[1],
                                        proj_f2_x - sp.f2_pdf[0]))

    text2 = f"{sp.f2_275_dist_ft:.1f}'"
    fs2 = 9.0
    tw2 = fitz.get_text_length(text2, fontname="helv", fontsize=fs2)
    page.insert_text(
        fitz.Point(cap2_x - tw2 / 2.0, cap2_y + fs2 / 3.0),
        text2, fontname="helv", fontsize=fs2, color=(0, 0, 0),
        morph=(fitz.Point(cap2_x, cap2_y), fitz.Matrix(-perp2_deg)))

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

    def _draw_drainfield(left, top, right, bot, r, lines=None, angle_deg=0):
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

        page.draw_polyline(pts, color=(0, 0, 0), width=0.8,
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
                lt, fontname="helv", fontsize=_df_fs, color=(0, 0, 0),
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
                     lines=["NEW", "DRAINFIELD"], angle_deg=_ndf_angle)

    return doc


def main():
    sp = build_site_plan_data()

    # --- Base site_plan.pdf ---
    doc = render_site_plan(sp)
    out_path = os.path.join(os.path.dirname(__file__), "site_plan.pdf")
    doc.save(out_path)
    print(f"Written to {out_path}")

    # --- site_plan_df.pdf (with drainfield annotations) ---
    render_site_plan_df(doc, sp)
    df_path = os.path.join(os.path.dirname(__file__), "site_plan_df.pdf")
    doc.save(df_path)
    doc.close()
    print(f"Written to {df_path}")
    print(f"Rotation: {sp.rotation_deg:.1f}\u00b0 CCW")
    print(f"F15 PDF position: ({sp.f15_pdf[0]:.1f}, {sp.f15_pdf[1]:.1f})")


if __name__ == "__main__":
    main()
