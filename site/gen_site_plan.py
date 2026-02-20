"""Generate site_plan.pdf: building outline overlaid on site survey."""

import math
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import fitz  # pymupdf
from floorplan.gen_floorplan import build_floorplan_data


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

    text = "11.0\u2032"
    fs = 8
    tw = fitz.get_text_length(text, fontname="helv", fontsize=fs)
    # Place centered at caption point, rotated to match perpendicular
    page.insert_text(
        fitz.Point(cap_x - tw / 2.0, cap_y + fs / 3.0),
        text, fontname="helv", fontsize=fs, color=(0, 0, 0),
        morph=(fitz.Point(cap_x, cap_y), fitz.Matrix(-perp_deg)))

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
