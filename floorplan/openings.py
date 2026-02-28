"""Opening position computation — outer wall openings and interior rough openings.

Single source of truth for all opening positions, consumed by both
gen_floorplan.py (polygon rendering) and gen_walls.py (parametric wall openings).
"""
import math
from typing import NamedTuple

from shared.types import Point, BBox, LineSeg
from shared.geometry import seg_vec, bbox_from_points
from floorplan.constants import (
    O1_WIDTH, O1_GAP_O2, O2_WIDTH, O2_GAP_O3,
    O3_GAP_F5, O3_WIDTH, O4_HALF_WIDTH,
    O5_OFFSET_FROM_IW2, O5_WIDTH, O6_WIDTH, O6_GAP_F10,
    O7_NW_GAP, O7_HALF_WIDTH,
    O8_HALF_WIDTH,
    STD_GAP,
    RO1_OFFSET_FROM_IW2, IW1_RO_WIDTH,
    IW2_RO_WIDTH,
    IW4_RO_WIDTH, IW9_RO_WIDTH, IW11_RO_WIDTH,
    RO3_WIDTH, RO3_IW7_GAP,
    IW6_RO_OFFSET_W, IW6_RO_WIDTH,
    WALL_4IN,
)


class OuterOpening(NamedTuple):
    """Opening in the outer wall, positioned as a 4-point polygon."""
    name: str
    seg_start: str       # e.g., "F2" — outline segment start point
    seg_end: str         # e.g., "F5" — outline segment end point
    poly: list[Point]    # 4 vertices: [outer_start, outer_end, inner_end, inner_start]


class RoughOpening(NamedTuple):
    """Rough opening in an interior wall."""
    name: str
    bbox: BBox           # w, s, e, n in survey coords (axis-aligned BB)
    wall_name: str       # "IW1", "IW2", etc.
    orientation: str     # "H", "V", or "R" (rotated)
    width: float         # opening width in feet
    poly: list[Point] | None = None  # [SW, SE, NE, NW] for rotated openings


class WallOpening(NamedTuple):
    """Opening on an outline segment, parameterized along the segment."""
    name: str
    seg_idx: int    # index in outline_segs (0-based)
    t_start: float  # parametric position [0, 1] along the segment
    t_end: float    # parametric position [0, 1] along the segment


def compute_outer_openings(pts: dict[str, Point], layout) -> list[OuterOpening]:
    """Compute all 11 outer-wall opening polygons.

    Each polygon has 4 vertices spanning from the F-face (outer) to the W-face (inner).
    Returns openings in order: O1, O2, ..., O11.
    """
    openings = []

    # O3: F2-F5, 8" from F5 along F5-F2 line (compute first; O2 and O1 depend on it)
    _dE2, _dN2, _seg2_len = seg_vec(pts["F2"], pts["F5"])
    _t3_end = 1 - O3_GAP_F5 / _seg2_len       # closer to F5
    _t3_start = 1 - (O3_GAP_F5 + O3_WIDTH) / _seg2_len  # farther from F5
    openings.append(OuterOpening("O3", "F2", "F5", [
        (pts["F2"][0] + _t3_start * _dE2, pts["F2"][1] + _t3_start * _dN2),
        (pts["F2"][0] + _t3_end * _dE2, pts["F2"][1] + _t3_end * _dN2),
        (pts["W2"][0] + _t3_end * (pts["W5"][0] - pts["W2"][0]),
         pts["W2"][1] + _t3_end * (pts["W5"][1] - pts["W2"][1])),
        (pts["W2"][0] + _t3_start * (pts["W5"][0] - pts["W2"][0]),
         pts["W2"][1] + _t3_start * (pts["W5"][1] - pts["W2"][1])),
    ]))

    # O2: F2-F5, 48" south of O3
    _t2_end = _t3_start - O2_GAP_O3 / _seg2_len    # north edge, 48" south of O3
    _t2_start = _t2_end - O2_WIDTH / _seg2_len      # south edge
    openings.append(OuterOpening("O2", "F2", "F5", [
        (pts["F2"][0] + _t2_start * _dE2, pts["F2"][1] + _t2_start * _dN2),
        (pts["F2"][0] + _t2_end * _dE2, pts["F2"][1] + _t2_end * _dN2),
        (pts["W2"][0] + _t2_end * (pts["W5"][0] - pts["W2"][0]),
         pts["W2"][1] + _t2_end * (pts["W5"][1] - pts["W2"][1])),
        (pts["W2"][0] + _t2_start * (pts["W5"][0] - pts["W2"][0]),
         pts["W2"][1] + _t2_start * (pts["W5"][1] - pts["W2"][1])),
    ]))

    # O1: F2-F5, 72" south of O2
    _t1_end = _t2_start - O1_GAP_O2 / _seg2_len    # north edge, 72" south of O2
    _t1_start = _t1_end - O1_WIDTH / _seg2_len      # south edge
    openings.append(OuterOpening("O1", "F2", "F5", [
        (pts["F2"][0] + _t1_start * _dE2, pts["F2"][1] + _t1_start * _dN2),
        (pts["F2"][0] + _t1_end * _dE2, pts["F2"][1] + _t1_end * _dN2),
        (pts["W2"][0] + _t1_end * (pts["W5"][0] - pts["W2"][0]),
         pts["W2"][1] + _t1_end * (pts["W5"][1] - pts["W2"][1])),
        (pts["W2"][0] + _t1_start * (pts["W5"][0] - pts["W2"][0]),
         pts["W2"][1] + _t1_start * (pts["W5"][1] - pts["W2"][1])),
    ]))

    # O4: F6-F7, centered on segment at t=0.5
    _dE4, _dN4, _seg4_len = seg_vec(pts["F6"], pts["F7"])
    _t4_half = O4_HALF_WIDTH / _seg4_len
    _t4_start = 0.5 - _t4_half
    _t4_end = 0.5 + _t4_half
    openings.append(OuterOpening("O4", "F6", "F7", [
        (pts["W6"][0] + _t4_start * (pts["W7"][0] - pts["W6"][0]),
         pts["W6"][1] + _t4_start * (pts["W7"][1] - pts["W6"][1])),
        (pts["W6"][0] + _t4_end * (pts["W7"][0] - pts["W6"][0]),
         pts["W6"][1] + _t4_end * (pts["W7"][1] - pts["W6"][1])),
        (pts["F6"][0] + _t4_end * _dE4, pts["F6"][1] + _t4_end * _dN4),
        (pts["F6"][0] + _t4_start * _dE4, pts["F6"][1] + _t4_start * _dN4),
    ]))

    # O5 and O6 share segment F9-F10
    _dE56, _dN56, _seg56_len = seg_vec(pts["F9"], pts["F10"])

    # O5: F9-F10, anchored at IW2s east face projection + offset
    _iw2_e_mid = ((layout.iw2s.poly[1][0] + layout.iw2s.poly[2][0]) / 2,
                  (layout.iw2s.poly[1][1] + layout.iw2s.poly[2][1]) / 2)
    _t5_ref = ((_iw2_e_mid[0] - pts["F9"][0]) * _dE56
               + (_iw2_e_mid[1] - pts["F9"][1]) * _dN56) / (_dE56**2 + _dN56**2)
    _t5_end = _t5_ref + O5_OFFSET_FROM_IW2 / _seg56_len
    _t5_start = _t5_end - O5_WIDTH / _seg56_len
    openings.append(OuterOpening("O5", "F9", "F10", [
        (pts["W9"][0] + _t5_start * (pts["W10"][0] - pts["W9"][0]),
         pts["W9"][1] + _t5_start * (pts["W10"][1] - pts["W9"][1])),
        (pts["W9"][0] + _t5_end * (pts["W10"][0] - pts["W9"][0]),
         pts["W9"][1] + _t5_end * (pts["W10"][1] - pts["W9"][1])),
        (pts["F9"][0] + _t5_end * _dE56, pts["F9"][1] + _t5_end * _dN56),
        (pts["F9"][0] + _t5_start * _dE56, pts["F9"][1] + _t5_start * _dN56),
    ]))

    # O6: F9-F10, offset from F10 end by O6_GAP_F10
    _t6_end = 1.0 - O6_GAP_F10 / _seg56_len
    _t6_start = _t6_end - O6_WIDTH / _seg56_len
    openings.append(OuterOpening("O6", "F9", "F10", [
        (pts["W9"][0] + _t6_start * (pts["W10"][0] - pts["W9"][0]),
         pts["W9"][1] + _t6_start * (pts["W10"][1] - pts["W9"][1])),
        (pts["W9"][0] + _t6_end * (pts["W10"][0] - pts["W9"][0]),
         pts["W9"][1] + _t6_end * (pts["W10"][1] - pts["W9"][1])),
        (pts["F9"][0] + _t6_end * _dE56, pts["F9"][1] + _t6_end * _dN56),
        (pts["F9"][0] + _t6_start * _dE56, pts["F9"][1] + _t6_start * _dN56),
    ]))

    # O7: F12-F13, diagonal — NW end 2' from F12, 6' opening
    dE, dN, seg_len = seg_vec(pts["F12"], pts["F13"])
    ts = O7_NW_GAP / seg_len
    te = ts + 2 * O7_HALF_WIDTH / seg_len
    openings.append(OuterOpening("O7", "F12", "F13", [
        (pts["F12"][0] + ts * dE, pts["F12"][1] + ts * dN),
        (pts["F12"][0] + te * dE, pts["F12"][1] + te * dN),
        (pts["W12"][0] + te * (pts["W13"][0] - pts["W12"][0]),
         pts["W12"][1] + te * (pts["W13"][1] - pts["W12"][1])),
        (pts["W12"][0] + ts * (pts["W13"][0] - pts["W12"][0]),
         pts["W12"][1] + ts * (pts["W13"][1] - pts["W12"][1])),
    ]))

    # O8: F14-F15, centered between IW5 south face projection and F15
    _dE8, _dN8, _seg8_len = seg_vec(pts["F14"], pts["F15"])
    _iw5_s_mid = ((layout.iw5.poly[0][0] + layout.iw5.poly[1][0]) / 2,
                  (layout.iw5.poly[0][1] + layout.iw5.poly[1][1]) / 2)
    _t8_iw5 = ((_iw5_s_mid[0] - pts["F14"][0]) * _dE8
               + (_iw5_s_mid[1] - pts["F14"][1]) * _dN8) / (_dE8**2 + _dN8**2)
    _t8_ctr = (_t8_iw5 + 1.0) / 2
    _t8_half = O8_HALF_WIDTH / _seg8_len
    _t8_start = _t8_ctr - _t8_half   # toward F14
    _t8_end = _t8_ctr + _t8_half     # toward F15
    openings.append(OuterOpening("O8", "F14", "F15", [
        (pts["F14"][0] + _t8_end * _dE8, pts["F14"][1] + _t8_end * _dN8),
        (pts["F14"][0] + _t8_start * _dE8, pts["F14"][1] + _t8_start * _dN8),
        (pts["W14"][0] + _t8_start * (pts["W15"][0] - pts["W14"][0]),
         pts["W14"][1] + _t8_start * (pts["W15"][1] - pts["W14"][1])),
        (pts["W14"][0] + _t8_end * (pts["W15"][0] - pts["W14"][0]),
         pts["W14"][1] + _t8_end * (pts["W15"][1] - pts["W14"][1])),
    ]))

    # O8a, O9, O10, O11: F18-F1 — parametric positions from layout (single source)
    _dE9, _dN9, _ = seg_vec(pts["F18"], pts["F1"])
    for _name, _ts, _te in [("O8a", layout.sw_t_o8a_start, layout.sw_t_o8a_end),
                             ("O9",  layout.sw_t_o9_start,  layout.sw_t_o9_end),
                             ("O10", layout.sw_t_o10_start, layout.sw_t_o10_end),
                             ("O11", layout.sw_t_o11_start, layout.sw_t_o11_end)]:
        openings.append(OuterOpening(_name, "F18", "F1", [
            (pts["F18"][0] + _ts * _dE9, pts["F18"][1] + _ts * _dN9),
            (pts["F18"][0] + _te * _dE9, pts["F18"][1] + _te * _dN9),
            (pts["W18"][0] + _te * (pts["W1"][0] - pts["W18"][0]),
             pts["W18"][1] + _te * (pts["W1"][1] - pts["W18"][1])),
            (pts["W18"][0] + _ts * (pts["W1"][0] - pts["W18"][0]),
             pts["W18"][1] + _ts * (pts["W1"][1] - pts["W18"][1])),
        ]))

    return openings


def _ro_poly_bbox(face_a: Point, face_b: Point,
                   unit: tuple[float, float],
                   start_d: float, end_d: float
                   ) -> tuple[list[Point], BBox]:
    """Build a rough-opening polygon and bbox on a wall.

    face_a, face_b: two reference corners on opposite faces of the wall.
    unit: unit vector along the wall's length axis.
    start_d, end_d: distances along the unit vector from the reference corners.

    Returns 4-point polygon [a+start, b+start, b+end, a+end] and BBox.
    """
    p0 = (face_a[0] + start_d * unit[0], face_a[1] + start_d * unit[1])
    p1 = (face_b[0] + start_d * unit[0], face_b[1] + start_d * unit[1])
    p2 = (face_b[0] + end_d * unit[0],   face_b[1] + end_d * unit[1])
    p3 = (face_a[0] + end_d * unit[0],   face_a[1] + end_d * unit[1])
    poly = [p0, p1, p2, p3]
    return poly, bbox_from_points(poly)


def _wall_unit(p1: Point, p2: Point) -> tuple[tuple[float, float], float]:
    """Unit vector from p1 toward p2 and their distance."""
    dx, dy, length = seg_vec(p1, p2)
    return (dx / length, dy / length), length


def _project(pt: Point, origin: Point, unit: tuple[float, float]) -> float:
    """Project pt onto axis defined by origin + t*unit; return t."""
    return (pt[0] - origin[0]) * unit[0] + (pt[1] - origin[1]) * unit[1]


def compute_rough_openings(pts: dict[str, Point], layout) -> list[RoughOpening]:
    """Compute all 7 interior rough-opening polygons and bounding boxes."""

    # RO1: in IW1, positioned relative to IW2 east face along IW1 length
    # RO1 uses a unique winding (along-face-a first) so we build its polygon
    # directly rather than using _ro_poly_bbox.
    _iw1_sw, _iw1_se, _iw1_nw = layout.iw1.poly[0], layout.iw1.poly[1], layout.iw1.poly[3]
    _un1_al, _ = _wall_unit(_iw1_sw, _iw1_se)
    _iw2_e_mid_r = ((layout.iw2.poly[1][0] + layout.iw2.poly[2][0]) / 2,
                    (layout.iw2.poly[1][1] + layout.iw2.poly[2][1]) / 2)
    _ro1_ref_d = _project(_iw2_e_mid_r, _iw1_sw, _un1_al)
    _ro1_start_d = _ro1_ref_d + RO1_OFFSET_FROM_IW2
    _ro1_end_d = _ro1_start_d + IW1_RO_WIDTH
    _ro1_poly = [
        (_iw1_sw[0] + _ro1_start_d * _un1_al[0], _iw1_sw[1] + _ro1_start_d * _un1_al[1]),
        (_iw1_sw[0] + _ro1_end_d * _un1_al[0],   _iw1_sw[1] + _ro1_end_d * _un1_al[1]),
        (_iw1_nw[0] + _ro1_end_d * _un1_al[0],   _iw1_nw[1] + _ro1_end_d * _un1_al[1]),
        (_iw1_nw[0] + _ro1_start_d * _un1_al[0], _iw1_nw[1] + _ro1_start_d * _un1_al[1]),
    ]
    _ro1_bb = bbox_from_points(_ro1_poly)

    # RO2: in IW11 (rotated), centered between IW12 N face and IW5 S face
    _iw11_se, _iw11_ne, _iw11_sw = layout.iw11.poly[1], layout.iw11.poly[2], layout.iw11.poly[0]
    _un11, _ = _wall_unit(_iw11_se, _iw11_ne)
    _ro2_iw12_n_d = _project(layout.iw12.poly[3], _iw11_se, _un11)
    _ro2_iw5_s_d = _project(layout.iw5.poly[0], _iw11_se, _un11)
    _ro2_center_d = (_ro2_iw12_n_d + _ro2_iw5_s_d) / 2
    _ro2_half = IW4_RO_WIDTH / 2
    _ro2_poly, _ro2_bb = _ro_poly_bbox(_iw11_se, _iw11_sw, _un11,
                                       _ro2_center_d - _ro2_half,
                                       _ro2_center_d + _ro2_half)

    # RO6: in IW11 (rotated), 50" centered between IW12 S face and W18-W1
    _ro6_iw12_s_d = _project(layout.iw12.poly[0], _iw11_se, _un11)
    _ro6_w18_d = _project(pts["W18"], _iw11_se, _un11)
    _ro6_center_d = (_ro6_iw12_s_d + _ro6_w18_d) / 2
    _ro6_half = IW11_RO_WIDTH / 2
    _ro6_poly, _ro6_bb = _ro_poly_bbox(_iw11_se, _iw11_sw, _un11,
                                       _ro6_center_d - _ro6_half,
                                       _ro6_center_d + _ro6_half)

    # RO7: in IW9 (rotated), 62" centered between IW7 S face and IW9 S end
    _iw9_se, _iw9_ne, _iw9_sw = layout.iw9.poly[1], layout.iw9.poly[2], layout.iw9.poly[0]
    _un9, _ = _wall_unit(_iw9_se, _iw9_ne)
    _ro7_iw7_s_d = _project(layout.iw7.poly[0], _iw9_se, _un9)
    _ro7_center_d = _ro7_iw7_s_d / 2
    _ro7_half = IW9_RO_WIDTH / 2
    _ro7_poly, _ro7_bb = _ro_poly_bbox(_iw9_se, _iw9_sw, _un9,
                                       _ro7_center_d - _ro7_half,
                                       _ro7_center_d + _ro7_half)

    # RO3: in IW9 (rotated), south edge 5" N of IW7 N face
    _ro3_iw7_n_d = _project(layout.iw7.poly[2], _iw9_se, _un9)
    _ro3_start_d = _ro3_iw7_n_d + RO3_IW7_GAP
    _ro3_poly, _ro3_bb = _ro_poly_bbox(_iw9_se, _iw9_sw, _un9,
                                       _ro3_start_d,
                                       _ro3_start_d + RO3_WIDTH)

    # RO4: in IW2o (oblique), centered along IW2o length
    _iw2o_sw, _iw2o_se, _iw2o_nw = layout.iw2o.poly[0], layout.iw2o.poly[1], layout.iw2o.poly[3]
    _un2_al, _len2r = _wall_unit(_iw2o_sw, _iw2o_nw)
    _ro4_half = IW2_RO_WIDTH / 2
    _ro4_poly, _ro4_bb = _ro_poly_bbox(_iw2o_sw, _iw2o_se, _un2_al,
                                       _len2r / 2 - _ro4_half,
                                       _len2r / 2 + _ro4_half)

    # RO5: in IW6, positioned relative to IW2s west face
    # IW6 is trapezoidal — compute each face independently.
    _iw6_sw, _iw6_se = layout.iw6.poly[0], layout.iw6.poly[1]
    _iw6_ne, _iw6_nw = layout.iw6.poly[2], layout.iw6.poly[3]
    _iw2_w_mid = ((layout.iw2s.poly[0][0] + layout.iw2s.poly[3][0]) / 2,
                  (layout.iw2s.poly[0][1] + layout.iw2s.poly[3][1]) / 2)
    _un6s, _ = _wall_unit(_iw6_sw, _iw6_se)
    _ref6s = _project(_iw2_w_mid, _iw6_sw, _un6s)
    _end6s = _ref6s - IW6_RO_OFFSET_W
    _start6s = _end6s - IW6_RO_WIDTH
    _un6n, _ = _wall_unit(_iw6_nw, _iw6_ne)
    _ref6n = _project(_iw2_w_mid, _iw6_nw, _un6n)
    _end6n = _ref6n - IW6_RO_OFFSET_W
    _start6n = _end6n - IW6_RO_WIDTH
    _ro5_sw = (_iw6_sw[0] + _start6s * _un6s[0], _iw6_sw[1] + _start6s * _un6s[1])
    _ro5_se = (_iw6_sw[0] + _end6s * _un6s[0],   _iw6_sw[1] + _end6s * _un6s[1])
    _ro5_ne = (_iw6_nw[0] + _end6n * _un6n[0],   _iw6_nw[1] + _end6n * _un6n[1])
    _ro5_nw = (_iw6_nw[0] + _start6n * _un6n[0], _iw6_nw[1] + _start6n * _un6n[1])
    _ro5_poly = [_ro5_sw, _ro5_se, _ro5_ne, _ro5_nw]
    _ro5_bb = bbox_from_points(_ro5_poly)

    return [
        RoughOpening("RO1", _ro1_bb, "IW1", "H", IW1_RO_WIDTH, _ro1_poly),
        RoughOpening("RO2", _ro2_bb, "IW11", "R", IW4_RO_WIDTH, _ro2_poly),
        RoughOpening("RO3", _ro3_bb, "IW9", "R", RO3_WIDTH, _ro3_poly),
        RoughOpening("RO4", _ro4_bb, "IW2o", "R", IW2_RO_WIDTH, _ro4_poly),
        RoughOpening("RO5", _ro5_bb, "IW6", "H", IW6_RO_WIDTH, _ro5_poly),
        RoughOpening("RO6", _ro6_bb, "IW11", "R", IW11_RO_WIDTH, _ro6_poly),
        RoughOpening("RO7", _ro7_bb, "IW9", "R", IW9_RO_WIDTH, _ro7_poly),
    ]


def _seg_param(pts, seg, point):
    """Compute parametric position t of a point along a LineSeg."""
    A = pts[seg.start]
    B = pts[seg.end]
    dx = B[0] - A[0]
    dy = B[1] - A[1]
    if abs(dx) < 1e-12 and abs(dy) < 1e-12:
        return 0.0  # degenerate zero-length segment
    if abs(dx) > abs(dy):
        return (point[0] - A[0]) / dx
    else:
        return (point[1] - A[1]) / dy


def outer_to_wall_openings(openings, outline_segs, pts):
    """Convert OuterOpenings to WallOpenings (parametric on outline segments).

    Uses the first two polygon vertices (on the outer face) to compute
    parametric positions along the outline segment.
    """
    seg_map = {(seg.start, seg.end): i for i, seg in enumerate(outline_segs)}
    result = []
    for o in openings:
        idx = seg_map[(o.seg_start, o.seg_end)]
        seg = outline_segs[idx]
        t1 = _seg_param(pts, seg, o.poly[0])
        t2 = _seg_param(pts, seg, o.poly[1])
        result.append(WallOpening(o.name, idx, min(t1, t2), max(t1, t2)))
    return result
