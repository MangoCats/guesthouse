"""Generate site_plan.pdf: building outline overlaid on site survey."""

import math
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import fitz  # pymupdf
from floorplan.gen_floorplan import build_floorplan_data
from shared.geometry import vert_isects


def main():
    data = build_floorplan_data()
    pts = data.pts
    outer_poly = data.outer_poly  # F-series outline (350 pts, arcs sampled)
    inner_poly = data.inner_poly  # W-series inner wall (352 pts)

    # --- Survey coordinate calibration ---
    # Scale: 1 inch = 30 ft on the survey; 1 inch = 72 PDF pts
    SCALE = 72.0 / 30.0  # 2.4 PDF pts per foot

    # 216.73' line endpoints (PDF coordinates, from least-squares line fitting)
    # Upper-right corner (251.53' meets 216.73')
    line_top = (698.9, 55.2)
    # Lower-right corner (216.73' meets 275.08')
    line_bot = (817.9, 557.8)

    # Direction of 216.73' line in PDF coords (x-right, y-down)
    ldx = line_bot[0] - line_top[0]  # +119 (slants left going up)
    ldy = line_bot[1] - line_top[1]  # +502.6
    llen = math.hypot(ldx, ldy)

    # In real-world coords (E=x-right, N=y-up), the line direction is:
    # dE = ldx = +119, dN = -ldy = -502.6
    # Angle of property line from E-axis (real-world):
    prop_angle = math.atan2(-ldy, ldx)  # atan2(-502.6, 119) ≈ -76.7°

    # F16→F17 direction in building coords
    f16 = pts["F16"]
    f17 = pts["F17"]
    f16f17_angle = math.atan2(f17[1] - f16[1], f17[0] - f16[0])  # ≈ -150° = 210°

    # Rotation needed: rotate building so F16→F17 is parallel to property line
    rotation = prop_angle - f16f17_angle  # ≈ 73° CCW

    # --- Rotation and placement ---
    cos_r = math.cos(rotation)
    sin_r = math.sin(rotation)

    def rotate_pt(e, n, ce, cn):
        """Rotate point (e,n) by `rotation` around center (ce,cn) in real-world coords."""
        de = e - ce
        dn = n - cn
        return (ce + de * cos_r - dn * sin_r,
                cn + de * sin_r + dn * cos_r)

    # F15 is the reference point for placement
    f15 = pts["F15"]  # (36.5, 5.0)

    # Inward perpendicular from 216.73' line (into property = roughly west)
    # Line direction real-world: (dE, dN) = (ldx, -ldy) = (+119, -502.6)
    # Line goes roughly south with slight east; right perp points west with slight south
    line_de = ldx / llen
    line_dn = -ldy / llen
    # Right perp of (dE, dN) = (dN, -dE) → points into property (west)
    inward_e = line_dn
    inward_n = -line_de
    # line_de ≈ +0.230, line_dn ≈ -0.973
    # inward_e = -0.973, inward_n = -0.230 → points roughly west with slight south ✓

    # Target: F15 placed near the bottom PATIO level on 216.73' line, 11' inside
    # Bottom PATIO is at approximately PDF y ≈ 435, shifted 20' up along the line
    patio_y = 435.0 - 10.0 * SCALE * (ldy / llen)  # move 10' toward upper corner
    # Find point on 216.73' line at this y
    t_patio = (patio_y - line_top[1]) / (line_bot[1] - line_top[1])
    line_at_patio_x = line_top[0] + t_patio * (line_bot[0] - line_top[0])

    # F15 target in PDF coords: 11' inward from property line
    setback = 11.0  # feet
    f15_pdf_x = line_at_patio_x + setback * SCALE * inward_e
    f15_pdf_y = patio_y + setback * SCALE * (-inward_n)  # y-flip for PDF

    def building_to_pdf(e, n):
        """Transform building coords (E,N) → PDF coords (x,y)."""
        # 1. Rotate around F15 in real-world coords
        re, rn = rotate_pt(e, n, f15[0], f15[1])
        # 2. Translate and scale: F15 maps to (f15_pdf_x, f15_pdf_y)
        pdf_x = f15_pdf_x + (re - f15[0]) * SCALE
        pdf_y = f15_pdf_y - (rn - f15[1]) * SCALE  # N-up → y-down
        return pdf_x, pdf_y

    # --- Generate PDF ---
    src = fitz.open("site/site_survey.pdf")
    doc = fitz.open()  # new PDF
    doc.insert_pdf(src)
    page = doc[0]

    # Draw building outline — single line, black pixels inside F boundary
    # Stroke width ≈ 80% of the survey property lines (~1pt → 0.8pt)
    STROKE_W = 1.6
    # Offset path inward by half stroke width so outer edge aligns with F boundary.
    # outer_poly→inner_poly span 8" of wall; we need half_stroke / SCALE feet inward.
    WALL_T = 8.0 / 12.0  # 8" F-to-W gap in feet
    half_stroke_ft = (STROKE_W / 2.0) / SCALE
    frac = half_stroke_ft / WALL_T  # fraction of the way from outer to inner

    # Interpolate: outer_poly has 350 pts, inner_poly has 352 pts (extra arc samples).
    # Resample inner_poly to match outer_poly length via nearest-index mapping.
    n_out = len(outer_poly)
    n_inn = len(inner_poly)
    draw_poly = []
    for i, (oe, on) in enumerate(outer_poly):
        j = round(i * (n_inn - 1) / (n_out - 1))
        ie, inn = inner_poly[j]
        draw_poly.append((oe + frac * (ie - oe), on + frac * (inn - on)))

    shape = page.new_shape()
    x0, y0 = building_to_pdf(*draw_poly[0])
    shape.draw_line(fitz.Point(x0, y0), fitz.Point(x0, y0))
    for pt in draw_poly[1:]:
        x1, y1 = building_to_pdf(*pt)
        shape.draw_line(fitz.Point(x0, y0), fitz.Point(x1, y1))
        x0, y0 = x1, y1
    x1, y1 = building_to_pdf(*draw_poly[0])
    shape.draw_line(fitz.Point(x0, y0), fitz.Point(x1, y1))
    shape.finish(color=(0, 0, 0), width=STROKE_W)

    # --- F15 to F2-F3 dimension line (perpendicular to F2-F3) ---
    # F2-F3 is vertical at E=0.5; perpendicular foot from F15 is (0.5, F15.N)
    f15_pdf = building_to_pdf(*f15)
    foot_pdf = building_to_pdf(pts["F2"][0], f15[1])  # (0.5, 5.0)
    shape.draw_line(fitz.Point(*f15_pdf), fitz.Point(*foot_pdf))
    shape.finish(color=(0, 0, 0), width=0.3)

    # Caption "36.0'" at midpoint, rotated 180° and shifted above the line
    dim_dx = foot_pdf[0] - f15_pdf[0]
    dim_dy = foot_pdf[1] - f15_pdf[1]
    dim_len = math.hypot(dim_dx, dim_dy)
    dim_deg = math.degrees(math.atan2(dim_dy, dim_dx))
    # Shift "above" the line from the flipped text's perspective:
    # reading direction after 180° flip is (-dim_dx, -dim_dy);
    # "above" = left perp of reading dir = (dim_dy, -dim_dx) normalised
    dim_shift = -3.5  # PDF pts; negative = other side of line
    dim_mid_x = (f15_pdf[0] + foot_pdf[0]) / 2.0 + dim_shift * dim_dy / dim_len
    dim_mid_y = (f15_pdf[1] + foot_pdf[1]) / 2.0 - dim_shift * dim_dx / dim_len
    dim_text = "36.0'"
    dim_fs = 6.0
    dim_tw = fitz.get_text_length(dim_text, fontname="helv", fontsize=dim_fs)
    page.insert_text(
        fitz.Point(dim_mid_x - dim_tw / 2.0, dim_mid_y + dim_fs / 3.0),
        dim_text, fontname="helv", fontsize=dim_fs, color=(0, 0, 0),
        morph=(fitz.Point(dim_mid_x, dim_mid_y), fitz.Matrix(-dim_deg - 180)))

    # --- N-S Interior Max Span dimension line (perpendicular to 36.0' line) ---
    # Find the easting where N-S interior span is maximum
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

    # Draw thin vertical line (in building coords) from south to north wall
    span_s_pdf = building_to_pdf(_best_e, _best_s)
    span_n_pdf = building_to_pdf(_best_e, _best_n)
    shape.draw_line(fitz.Point(*span_s_pdf), fitz.Point(*span_n_pdf))
    shape.finish(color=(0, 0, 0), width=0.3)

    # Caption "27.1'" at midpoint, same style as 36.0' line
    ns_dx = span_n_pdf[0] - span_s_pdf[0]
    ns_dy = span_n_pdf[1] - span_s_pdf[1]
    ns_len = math.hypot(ns_dx, ns_dy)
    ns_deg = math.degrees(math.atan2(ns_dy, ns_dx))
    ns_shift = -3.5  # same side convention as 36.0'
    ns_mid_x = (span_s_pdf[0] + span_n_pdf[0]) / 2.0 + ns_shift * ns_dy / ns_len
    ns_mid_y = (span_s_pdf[1] + span_n_pdf[1]) / 2.0 - ns_shift * ns_dx / ns_len
    ns_text = "27.1'"
    ns_fs = 6.0
    ns_tw = fitz.get_text_length(ns_text, fontname="helv", fontsize=ns_fs)
    page.insert_text(
        fitz.Point(ns_mid_x - ns_tw / 2.0, ns_mid_y + ns_fs / 3.0),
        ns_text, fontname="helv", fontsize=ns_fs, color=(0, 0, 0),
        morph=(fitz.Point(ns_mid_x, ns_mid_y), fitz.Matrix(-ns_deg - 180)))

    # --- "PROPOSED CONC. GUEST HOUSE" label inside outline ---
    # Place at inner polygon centroid, horizontal (parallel to page bottom)
    _cx = sum(p[0] for p in inner_poly) / len(inner_poly)
    _cy = sum(p[1] for p in inner_poly) / len(inner_poly)
    label_pdf_raw = building_to_pdf(_cx, _cy + 2.0)  # 2' north of centroid
    label_pdf = (label_pdf_raw[0], label_pdf_raw[1] + 5.0 * SCALE)  # then 5' down on page
    label_fs = 8.0
    label_lines = ["PROPOSED", "CONC.", "GUEST", "HOUSE"]
    label_lh = label_fs * 1.15  # line height
    # Total block height; start y so block is vertically centered
    block_h = label_lh * len(label_lines)
    start_y = label_pdf[1] - block_h / 2.0 + label_fs  # baseline of first line
    for i, line in enumerate(label_lines):
        lw = fitz.get_text_length(line, fontname="helv", fontsize=label_fs)
        page.insert_text(
            fitz.Point(label_pdf[0] - lw / 2.0, start_y + i * label_lh),
            line, fontname="helv", fontsize=label_fs, color=(0, 0, 0))

    # --- 11.0' setback caption ---
    # Midpoint of F16-F17 in PDF coords
    f16_pdf = building_to_pdf(*f16)
    f17_pdf = building_to_pdf(*f17)
    mid_f16f17_x = (f16_pdf[0] + f17_pdf[0]) / 2.0
    mid_f16f17_y = (f16_pdf[1] + f17_pdf[1]) / 2.0

    # Project midpoint perpendicularly onto 216.73' line.
    # Line parameterised as P = line_top + t*(line_bot - line_top).
    # t = dot(mid - line_top, line_dir) / llen^2
    t_proj = ((mid_f16f17_x - line_top[0]) * ldx +
              (mid_f16f17_y - line_top[1]) * ldy) / (llen * llen)
    proj_x = line_top[0] + t_proj * ldx
    proj_y = line_top[1] + t_proj * ldy

    # Caption at midpoint between F16-F17 midpoint and its projection on the line
    cap_x = (mid_f16f17_x + proj_x) / 2.0
    cap_y = (mid_f16f17_y + proj_y) / 2.0

    # Rotate text to sit on a perpendicular to the 216.73' line,
    # matching the survey captions (46.7', 39.5', etc.)
    perp_deg = math.degrees(math.atan2(proj_y - mid_f16f17_y,
                                       proj_x - mid_f16f17_x))  # ≈ -13.3°

    text = "11.0'"
    fs = 9.0
    tw = fitz.get_text_length(text, fontname="helv", fontsize=fs)
    # Place centered at caption point, rotated to match perpendicular
    page.insert_text(
        fitz.Point(cap_x - tw / 2.0, cap_y + fs / 3.0),
        text, fontname="helv", fontsize=fs, color=(0, 0, 0),
        morph=(fitz.Point(cap_x, cap_y), fitz.Matrix(-perp_deg)))

    # --- F2 distance to 275.08' line caption ---
    # 275.08' line endpoints (from least-squares fit on rasterized survey)
    # Lower-right corner shared with 216.73' line
    bot_right = line_bot  # (817.9, 557.8)
    bot_left = (160.0, 561.9)

    # Direction and length of 275.08' line in PDF coords
    bdx = bot_right[0] - bot_left[0]
    bdy = bot_right[1] - bot_left[1]
    blen = math.hypot(bdx, bdy)

    # F2 in PDF coords
    f2 = pts["F2"]
    f2_pdf = building_to_pdf(*f2)

    # Project F2 perpendicularly onto 275.08' line
    t_f2 = ((f2_pdf[0] - bot_left[0]) * bdx +
            (f2_pdf[1] - bot_left[1]) * bdy) / (blen * blen)
    proj_f2_x = bot_left[0] + t_f2 * bdx
    proj_f2_y = bot_left[1] + t_f2 * bdy

    # Distance in feet
    dist_pts = math.hypot(f2_pdf[0] - proj_f2_x, f2_pdf[1] - proj_f2_y)
    dist_ft = dist_pts / SCALE

    # Caption at midpoint of perpendicular
    cap2_x = (f2_pdf[0] + proj_f2_x) / 2.0
    cap2_y = (f2_pdf[1] + proj_f2_y) / 2.0

    # Rotation: direction from F2 toward 275.08' line, matching 80.6'/74.9' style
    perp2_deg = math.degrees(math.atan2(proj_f2_y - f2_pdf[1],
                                        proj_f2_x - f2_pdf[0]))

    text2 = f"{dist_ft:.1f}'"
    fs2 = 9.0
    tw2 = fitz.get_text_length(text2, fontname="helv", fontsize=fs2)
    page.insert_text(
        fitz.Point(cap2_x - tw2 / 2.0, cap2_y + fs2 / 3.0),
        text2, fontname="helv", fontsize=fs2, color=(0, 0, 0),
        morph=(fitz.Point(cap2_x, cap2_y), fitz.Matrix(-perp2_deg)))

    shape.commit()

    out_path = os.path.join(os.path.dirname(__file__), "site_plan.pdf")
    doc.save(out_path)
    doc.close()
    src.close()
    print(f"Written to {out_path}")
    print(f"Rotation: {math.degrees(rotation):.1f}° CCW")
    print(f"F15 PDF position: ({f15_pdf_x:.1f}, {f15_pdf_y:.1f})")


if __name__ == "__main__":
    main()
