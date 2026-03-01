"""Generate plumbing plan SVG.

Based on the daybed floorplan variant, showing only plumbing-relevant
fixtures: washer, toilets, sinks, water heater, fridge, and dishwasher.
"""
import math
import os, sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

from floorplan.gen_floorplan import build_floorplan_data, render_floorplan_svg

# Property corner coordinates in site survey PDF coords
# (from site/gen_site_plan.py calibration constants)
_PDF_NE = (698.9, 55.2)     # 251.53' meets 216.73'
_PDF_SE = (817.9, 557.8)    # 216.73' meets 275.08'
_PDF_SW = (160.0, 561.9)    # 275.08' meets 163.69'
_PDF_NW = (108.0, 174.5)    # 251.53' meets 163.69'

# Placement constraints (same as site/gen_site_plan.py)
_P45_DIST_216 = 11.0    # feet
_P3_DIST_275 = 25.5     # feet
_SCALE = 72.0 / 30.0    # 2.4 PDF pts per foot


def _compute_boundary_corners(pts):
    """Compute property boundary corners in building (E,N) coordinates.

    Replicates the building_to_pdf transform from site/gen_site_plan.py
    and inverts it to map PDF corner coords → building coords.
    """
    # Direction of 216.73' property line (NE→SE) in PDF coords
    ldx = _PDF_SE[0] - _PDF_NE[0]
    ldy = _PDF_SE[1] - _PDF_NE[1]
    llen = math.hypot(ldx, ldy)

    # Property line angle in real-world coords (E-right, N-up)
    prop_angle = math.atan2(-ldy, ldx)

    # F16→F17 angle in building coords
    f16, f17 = pts["F16"], pts["F17"]
    bld_angle = math.atan2(f17[1] - f16[1], f17[0] - f16[0])

    # Rotation: align building F16→F17 with property line
    rotation = prop_angle - bld_angle
    cos_r = math.cos(rotation)
    sin_r = math.sin(rotation)

    # Solve for F15 PDF position using placement constraints
    f15 = pts["F15"]
    p4, p3 = pts["P4"], pts["P3"]

    off_p4_x = ((p4[0] - f15[0]) * cos_r - (p4[1] - f15[1]) * sin_r) * _SCALE
    off_p4_y = -((p4[0] - f15[0]) * sin_r + (p4[1] - f15[1]) * cos_r) * _SCALE
    off_p3_x = ((p3[0] - f15[0]) * cos_r - (p3[1] - f15[1]) * sin_r) * _SCALE
    off_p3_y = -((p3[0] - f15[0]) * sin_r + (p3[1] - f15[1]) * cos_r) * _SCALE

    # Constraint A: P4 is P45_DIST_216 inside the 216.73' line
    a1 = -ldy / llen
    b1 = ldx / llen
    c1 = (_P45_DIST_216 * _SCALE
          + a1 * (_PDF_NE[0] - off_p4_x) + b1 * (_PDF_NE[1] - off_p4_y))

    # Constraint B: P3 is P3_DIST_275 inside the 275.08' line
    bdx = _PDF_SE[0] - _PDF_SW[0]
    bdy = _PDF_SE[1] - _PDF_SW[1]
    blen = math.hypot(bdx, bdy)
    a2 = bdy / blen
    b2 = -bdx / blen
    c2 = (_P3_DIST_275 * _SCALE
          + a2 * (_PDF_SW[0] - off_p3_x) + b2 * (_PDF_SW[1] - off_p3_y))

    # Solve 2×2 linear system for F15 PDF position
    det = a1 * b2 - a2 * b1
    f15_px = (c1 * b2 - c2 * b1) / det
    f15_py = (a1 * c2 - a2 * c1) / det

    # Inverse transform: PDF coords → building coords
    def pdf_to_bld(px, py):
        dre = (px - f15_px) / _SCALE
        drn = -(py - f15_py) / _SCALE
        e = f15[0] + dre * cos_r + drn * sin_r
        n = f15[1] - dre * sin_r + drn * cos_r
        return (e, n)

    sw = pdf_to_bld(*_PDF_SW)
    nw = pdf_to_bld(*_PDF_NW)

    # Clip west boundary start to 40' from SW toward NW (avoids
    # extending the SVG all the way to the distant NW corner).
    dx_w = nw[0] - sw[0]
    dy_w = nw[1] - sw[1]
    west_len = math.hypot(dx_w, dy_w)
    west_start = (sw[0] + 40 * dx_w / west_len,
                  sw[1] + 40 * dy_w / west_len)

    return {
        'sw': sw,
        'se': pdf_to_bld(*_PDF_SE),
        'nw': west_start,
    }


if __name__ == "__main__":
    data = build_floorplan_data()
    boundary = _compute_boundary_corners(data.pts)
    svg = render_floorplan_svg(data, room_title="Plumbing Plan",
                               db=True, plumbing=True, boundary=boundary)
    base_dir = os.path.dirname(os.path.abspath(__file__))
    path = os.path.join(base_dir, "plumbing.svg")
    with open(path, "w") as f:
        f.write(svg)
    print(f"Plumbing plan written to {path}")
