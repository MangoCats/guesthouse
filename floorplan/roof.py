"""Compute R-series roof corner geometry from F-series outline."""
import math
from typing import NamedTuple

from shared.types import Point
from shared.geometry import left_norm, line_isect, poly_area
from floorplan.constants import ROOF_OVERHANG


class RoofGeometry(NamedTuple):
    """Roof corner points and arc parameters."""
    pts: dict[str, Point]      # R1-R7, R3s, R3e, R4s, R4e
    r3_radius: float           # R3 fillet radius (R_a5 + overhang)
    r3_center: Point           # R3 fillet center (= C5)
    r4_radius: float           # R4 fillet radius (R_a11 + overhang)
    r4_center: Point           # R4 fillet center (= C11)
    area: float                # Roof area in sq ft


def _offset_line(p1: Point, p2: Point, dist: float) -> tuple[Point, Point]:
    """Two points on the line through p1->p2 offset left by dist."""
    ln = left_norm(p1, p2)
    return ((p1[0] + dist * ln[0], p1[1] + dist * ln[1]),
            (p2[0] + dist * ln[0], p2[1] + dist * ln[1]))


def _isect(oL_a: tuple[Point, Point], oL_b: tuple[Point, Point]) -> Point:
    """Intersection of two offset lines (each given as two points)."""
    a1, a2 = oL_a
    d1 = (a2[0] - a1[0], a2[1] - a1[1])
    b1, b2 = oL_b
    d2 = (b2[0] - b1[0], b2[1] - b1[1])
    return line_isect(a1, d1, b1, d2)


def _arc_tangent_pt(center: Point, tangent_pt: Point,
                    orig_r: float, new_r: float) -> Point:
    """Point on circle(center, new_r) in the direction of tangent_pt."""
    scale = new_r / orig_r
    return (center[0] + (tangent_pt[0] - center[0]) * scale,
            center[1] + (tangent_pt[1] - center[1]) * scale)


def compute_roof_geometry(fp_pts: dict[str, Point],
                          radii: dict[str, float]) -> RoofGeometry:
    """Compute R-series roof corner points from F-series outline."""
    oh = ROOF_OVERHANG
    R_a5 = radii["R_a5"]
    R_a11 = radii["R_a11"]

    # Offset straight wall segments outward (left of CW traversal = exterior)
    oL_F20_F1 = _offset_line(fp_pts["F20"], fp_pts["F1"], oh)
    oL_F2_F3 = _offset_line(fp_pts["F2"], fp_pts["F3"], oh)
    oL_F4_F5 = _offset_line(fp_pts["F4"], fp_pts["F5"], oh)
    oL_F6_F7 = _offset_line(fp_pts["F6"], fp_pts["F7"], oh)
    oL_F11a_F11b = _offset_line(fp_pts["F11a"], fp_pts["F11b"], oh)
    oL_F12_F13 = _offset_line(fp_pts["F12"], fp_pts["F13"], oh)
    oL_F14_F15 = _offset_line(fp_pts["F14"], fp_pts["F15"], oh)
    oL_F16_F17 = _offset_line(fp_pts["F16"], fp_pts["F17"], oh)

    pts: dict[str, Point] = {}

    # Sharp corners
    pts["R1"] = _isect(oL_F20_F1, oL_F2_F3)
    pts["R2"] = _isect(oL_F2_F3, oL_F4_F5)
    pts["R5"] = _isect(oL_F12_F13, oL_F14_F15)
    pts["R6"] = _isect(oL_F14_F15, oL_F16_F17)
    pts["R7"] = _isect(oL_F16_F17, oL_F20_F1)

    # R3: filleted corner (F4-F5 / F6-F7), radius R_a5 + oh, center C5
    r3_r = R_a5 + oh
    r3_c = fp_pts["C5"]
    pts["R3s"] = _arc_tangent_pt(r3_c, fp_pts["F5"], R_a5, r3_r)
    pts["R3e"] = _arc_tangent_pt(r3_c, fp_pts["F6"], R_a5, r3_r)

    # R4: filleted corner (F11a-F11b / F12-F13), radius R_a11 + oh, center C11
    r4_r = R_a11 + oh
    r4_c = fp_pts["C11"]
    pts["R4s"] = _arc_tangent_pt(r4_c, fp_pts["F11b"], R_a11, r4_r)
    pts["R4e"] = _arc_tangent_pt(r4_c, fp_pts["F12"], R_a11, r4_r)

    # Compute area: polygon with arc segments sampled into points
    poly: list[Point] = [pts["R1"], pts["R2"], pts["R3s"]]

    # R3 arc: CW from R3s to R3e around r3_c
    a3_start = math.atan2(pts["R3s"][1] - r3_c[1], pts["R3s"][0] - r3_c[0])
    a3_end = math.atan2(pts["R3e"][1] - r3_c[1], pts["R3e"][0] - r3_c[0])
    if a3_end > a3_start:
        a3_end -= 2 * math.pi
    n_arc = 30
    for i in range(1, n_arc):
        a = a3_start + (a3_end - a3_start) * i / n_arc
        poly.append((r3_c[0] + r3_r * math.cos(a), r3_c[1] + r3_r * math.sin(a)))
    poly.append(pts["R3e"])

    # Straight segment across the top (R3e to R4s)
    poly.append(pts["R4s"])

    # R4 arc: CW from R4s to R4e around r4_c
    a4_start = math.atan2(pts["R4s"][1] - r4_c[1], pts["R4s"][0] - r4_c[0])
    a4_end = math.atan2(pts["R4e"][1] - r4_c[1], pts["R4e"][0] - r4_c[0])
    if a4_end > a4_start:
        a4_end -= 2 * math.pi
    for i in range(1, n_arc):
        a = a4_start + (a4_end - a4_start) * i / n_arc
        poly.append((r4_c[0] + r4_r * math.cos(a), r4_c[1] + r4_r * math.sin(a)))
    poly.append(pts["R4e"])

    poly.extend([pts["R5"], pts["R6"], pts["R7"]])

    area = poly_area(poly)

    return RoofGeometry(
        pts=pts,
        r3_radius=r3_r, r3_center=r3_c,
        r4_radius=r4_r, r4_center=r4_c,
        area=area,
    )


def roof_segments(roof: RoofGeometry) -> list:
    """Build the roof outline as line/arc segment elements.

    Returns list of ("line", x1, y1, x2, y2) and
    ("arc", cx, cy, r, a1_deg, a2_deg) tuples — same format as T-path
    elements.  The outline is CW (R1→R2→R3 arc→...→R7→R1).
    """
    pts = roof.pts
    r3_c, r3_r = roof.r3_center, roof.r3_radius
    r4_c, r4_r = roof.r4_center, roof.r4_radius

    def _line(a, b):
        return ("line", a[0], a[1], b[0], b[1])

    def _arc_cw(center, radius, sp, ep):
        a1 = math.degrees(math.atan2(sp[1] - center[1], sp[0] - center[0]))
        a2 = math.degrees(math.atan2(ep[1] - center[1], ep[0] - center[0]))
        sweep = (a1 - a2) % 360
        a2 = a1 - sweep
        return ("arc", center[0], center[1], radius, a1, a2)

    return [
        _line(pts["R1"], pts["R2"]),
        _line(pts["R2"], pts["R3s"]),
        _arc_cw(r3_c, r3_r, pts["R3s"], pts["R3e"]),
        _line(pts["R3e"], pts["R4s"]),
        _arc_cw(r4_c, r4_r, pts["R4s"], pts["R4e"]),
        _line(pts["R4e"], pts["R5"]),
        _line(pts["R5"], pts["R6"]),
        _line(pts["R6"], pts["R7"]),
        _line(pts["R7"], pts["R1"]),
    ]


def roof_polyline(roof: RoofGeometry, n_arc: int = 30) -> list[Point]:
    """Build the roof outline as a list of (E,N) points, sampling arcs."""
    pts = roof.pts
    poly: list[Point] = [pts["R1"], pts["R2"], pts["R3s"]]

    # R3 arc: CW from R3s to R3e around r3_center
    r3_c = roof.r3_center
    r3_r = roof.r3_radius
    a_start = math.atan2(pts["R3s"][1] - r3_c[1], pts["R3s"][0] - r3_c[0])
    a_end = math.atan2(pts["R3e"][1] - r3_c[1], pts["R3e"][0] - r3_c[0])
    if a_end > a_start:
        a_end -= 2 * math.pi
    for i in range(1, n_arc):
        a = a_start + (a_end - a_start) * i / n_arc
        poly.append((r3_c[0] + r3_r * math.cos(a), r3_c[1] + r3_r * math.sin(a)))
    poly.append(pts["R3e"])

    poly.append(pts["R4s"])

    # R4 arc: CW from R4s to R4e around r4_center
    r4_c = roof.r4_center
    r4_r = roof.r4_radius
    a_start = math.atan2(pts["R4s"][1] - r4_c[1], pts["R4s"][0] - r4_c[0])
    a_end = math.atan2(pts["R4e"][1] - r4_c[1], pts["R4e"][0] - r4_c[0])
    if a_end > a_start:
        a_end -= 2 * math.pi
    for i in range(1, n_arc):
        a = a_start + (a_end - a_start) * i / n_arc
        poly.append((r4_c[0] + r4_r * math.cos(a), r4_c[1] + r4_r * math.sin(a)))
    poly.append(pts["R4e"])

    poly.extend([pts["R5"], pts["R6"], pts["R7"]])
    return poly
