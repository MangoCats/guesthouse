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

    # Target: F15 placed at the bottom PATIO level on 216.73' line, 11' inside
    # Bottom PATIO is at approximately PDF y ≈ 435
    patio_y = 435.0
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

    # Draw F-series outline (outer wall)
    shape = page.new_shape()
    x0, y0 = building_to_pdf(*outer_poly[0])
    shape.draw_line(fitz.Point(x0, y0), fitz.Point(x0, y0))
    for pt in outer_poly[1:]:
        x1, y1 = building_to_pdf(*pt)
        shape.draw_line(fitz.Point(x0, y0), fitz.Point(x1, y1))
        x0, y0 = x1, y1
    # Close
    x1, y1 = building_to_pdf(*outer_poly[0])
    shape.draw_line(fitz.Point(x0, y0), fitz.Point(x1, y1))
    shape.finish(color=(0, 0, 0), width=1.5)

    # Draw W-series outline (inner wall)
    x0, y0 = building_to_pdf(*inner_poly[0])
    for pt in inner_poly[1:]:
        x1, y1 = building_to_pdf(*pt)
        shape.draw_line(fitz.Point(x0, y0), fitz.Point(x1, y1))
        x0, y0 = x1, y1
    x1, y1 = building_to_pdf(*inner_poly[0])
    shape.draw_line(fitz.Point(x0, y0), fitz.Point(x1, y1))
    shape.finish(color=(0, 0, 0), width=0.75)

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
