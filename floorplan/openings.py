"""Opening position computation — outer wall openings and interior rough openings.

Single source of truth for all opening positions, consumed by both
gen_floorplan.py (polygon rendering) and gen_walls.py (parametric wall openings).

All openings are defined as polygons using building axis vectors derived from
W20→W1 (south wall direction).  No axis-aligned assumptions.
"""
import math
from typing import NamedTuple

from shared.types import Point
from floorplan.constants import (
    O1_WIDTH, O2_WIDTH,
    IW2_RO_OFFSET_S, IW2_RO_WIDTH,
    O3_GAP_F5, O3_WIDTH, O4_HALF_WIDTH,
    O5_E_FROM_IW2, O5_WIDTH, O6_WIDTH, O6_GAP_F10,
    O7_NW_GAP, O7_HALF_WIDTH,
    O8_HALF_WIDTH,
    RO1_OFFSET_E_IW2, IW1_RO_WIDTH,
    IW4_RO_WIDTH, IW9_RO_WIDTH, IW11_RO_WIDTH, IW16_RO_WIDTH,
    IW6_RO_OFFSET_W, IW6_RO_WIDTH,
)


class OuterOpening(NamedTuple):
    """Opening in the outer wall, positioned as a 4-point polygon."""
    name: str
    seg_start: str       # e.g., "F2" — outline segment start point
    seg_end: str         # e.g., "F3" — outline segment end point
    poly: list[Point]    # 4 vertices: [outer_start, outer_end, inner_end, inner_start]


class RoughOpening(NamedTuple):
    """Rough opening in an interior wall."""
    name: str
    poly: list[Point]    # [SW, SE, NE, NW] — always present
    wall_name: str       # "IW1", "IW2", etc.
    orientation: str     # "H", "V", or "R" (rotated)


class WallOpening(NamedTuple):
    """Opening on an outline segment, parameterized along the segment."""
    name: str
    seg_idx: int    # index in outline_segs (0-based)
    t_start: float  # parametric position [0, 1] along the segment
    t_end: float    # parametric position [0, 1] along the segment


def compute_outer_openings(pts, layout) -> list[OuterOpening]:
    """Compute all 11 outer-wall opening polygons.

    Each polygon has 4 vertices spanning from the F-face (outer) to the W-face (inner).
    All positions use building axis vectors for rotation independence.
    Returns openings in order: O1, O2, ..., O11.
    """
    # Building axis vectors
    _w20, _w1 = pts["W20"], pts["W1"]
    _dE = _w1[0] - _w20[0]; _dN = _w1[1] - _w20[1]
    _seg_len = math.sqrt(_dE**2 + _dN**2)
    _along_E = _dE / _seg_len; _along_N = _dN / _seg_len
    _norm_E = _along_N; _norm_N = -_along_E
    _eE = -_along_E; _eN = -_along_N

    def _bn(p):
        """Building northing (perpendicular distance from south wall)."""
        return (p[0] - _w1[0]) * _norm_E + (p[1] - _w1[1]) * _norm_N

    def _be(p):
        """Building easting (distance east along south wall from W1)."""
        return (p[0] - _w1[0]) * _eE + (p[1] - _w1[1]) * _eN

    def _bp(be, bn):
        """World-coordinate point from building easting/northing."""
        return (_w1[0] + be * _eE + bn * _norm_E,
                _w1[1] + be * _eN + bn * _norm_N)

    def _seg_opening(seg_start, seg_end, t_start, t_end):
        """Build outer opening polygon from parametric t on F-face and W-face."""
        f_a, f_b = pts[seg_start], pts[seg_end]
        w_a, w_b = pts["W" + seg_start[1:]], pts["W" + seg_end[1:]]
        dF = (f_b[0] - f_a[0], f_b[1] - f_a[1])
        dW = (w_b[0] - w_a[0], w_b[1] - w_a[1])
        return [
            (f_a[0] + t_start * dF[0], f_a[1] + t_start * dF[1]),
            (f_a[0] + t_end * dF[0],   f_a[1] + t_end * dF[1]),
            (w_a[0] + t_end * dW[0],   w_a[1] + t_end * dW[1]),
            (w_a[0] + t_start * dW[0], w_a[1] + t_start * dW[1]),
        ]

    def _t_from_bn(seg_start, seg_end, bn):
        """Parametric t on F-segment at given building northing."""
        return (bn - _bn(pts[seg_start])) / (_bn(pts[seg_end]) - _bn(pts[seg_start]))

    openings = []

    # O1: F2→F3, building-axis rectangle centered at IW16 midpoint northing
    _iw16_ctr_bn = (_bn(layout.iw16_poly[0]) + _bn(layout.iw16_poly[2])) / 2
    _o1_s_bn = _iw16_ctr_bn - O1_WIDTH / 2
    _o1_n_bn = _iw16_ctr_bn + O1_WIDTH / 2
    _f3_be = _be(pts["F3"])
    _w3_be = _be(pts["W3"])
    openings.append(OuterOpening("O1", "F2", "F3", [
        _bp(_f3_be, _o1_s_bn), _bp(_f3_be, _o1_n_bn),
        _bp(_w3_be, _o1_n_bn), _bp(_w3_be, _o1_s_bn),
    ]))

    # O2: F4→F5, centered at RO4 center building-northing
    _iw6_s_bn = _bn(layout.iw6_poly[0])
    _ro4_ctr_bn = _iw6_s_bn - IW2_RO_OFFSET_S - IW2_RO_WIDTH / 2
    _f4f5_len = math.sqrt((pts["F5"][0] - pts["F4"][0])**2
                          + (pts["F5"][1] - pts["F4"][1])**2)
    _t2_ctr = _t_from_bn("F4", "F5", _ro4_ctr_bn)
    _t2_half = (O2_WIDTH / 2) / _f4f5_len
    openings.append(OuterOpening("O2", "F4", "F5",
        _seg_opening("F4", "F5", _t2_ctr - _t2_half, _t2_ctr + _t2_half)))

    # O3: F4→F5, O3_GAP_F5 from F5
    _t3_end = 1 - O3_GAP_F5 / _f4f5_len
    _t3_start = 1 - (O3_GAP_F5 + O3_WIDTH) / _f4f5_len
    openings.append(OuterOpening("O3", "F4", "F5",
        _seg_opening("F4", "F5", _t3_start, _t3_end)))

    # O4: F6→F7, building-axis rectangle centered on F6-F7 midpoint easting
    _f6_be = _be(pts["F6"])
    _f7_be = _be(pts["F7"])
    _o4_mid_be = (_f6_be + _f7_be) / 2
    _o4_w_be = _o4_mid_be - O4_HALF_WIDTH
    _o4_e_be = _o4_mid_be + O4_HALF_WIDTH
    _w6_bn = _bn(pts["W6"])
    _f6_bn = _bn(pts["F6"])
    openings.append(OuterOpening("O4", "F6", "F7", [
        _bp(_o4_w_be, _w6_bn), _bp(_o4_e_be, _w6_bn),
        _bp(_o4_e_be, _f6_bn), _bp(_o4_w_be, _f6_bn),
    ]))

    # O5: F9→F10, building-axis rectangle, east edge at IW2 east + O5_E_FROM_IW2
    _iw2_e_be = _be(layout.iw2_poly[1])
    _o5_e_be = _iw2_e_be + O5_E_FROM_IW2
    _o5_w_be = _o5_e_be - O5_WIDTH
    _w9_bn = _bn(pts["W9"])
    _f9_bn = _bn(pts["F9"])
    openings.append(OuterOpening("O5", "F9", "F10", [
        _bp(_o5_w_be, _w9_bn), _bp(_o5_e_be, _w9_bn),
        _bp(_o5_e_be, _f9_bn), _bp(_o5_w_be, _f9_bn),
    ]))

    # O6: F9→F10, building-axis rectangle, east edge O6_GAP_F10 west of F10
    _f10_be = _be(pts["F10"])
    _o6_e_be = _f10_be - O6_GAP_F10
    _o6_w_be = _o6_e_be - O6_WIDTH
    openings.append(OuterOpening("O6", "F9", "F10", [
        _bp(_o6_w_be, _w9_bn), _bp(_o6_e_be, _w9_bn),
        _bp(_o6_e_be, _f9_bn), _bp(_o6_w_be, _f9_bn),
    ]))

    # O7: F12→F13, O7_NW_GAP from F12
    _f12f13_len = math.sqrt((pts["F13"][0] - pts["F12"][0])**2
                            + (pts["F13"][1] - pts["F12"][1])**2)
    _t7_s = O7_NW_GAP / _f12f13_len
    _t7_e = _t7_s + 2 * O7_HALF_WIDTH / _f12f13_len
    openings.append(OuterOpening("O7", "F12", "F13",
        _seg_opening("F12", "F13", _t7_s, _t7_e)))

    # O8: F14→F15, building-axis rectangle centered between IW5 south and F15
    _iw5_s_bn = _bn(layout.iw5_poly[0])
    _f15_bn = _bn(pts["F15"])
    _o8_ctr_bn = (_iw5_s_bn + _f15_bn) / 2
    _f15_be = _be(pts["F15"])
    _w15_be = _be(pts["W15"])
    openings.append(OuterOpening("O8", "F14", "F15", [
        _bp(_f15_be, _o8_ctr_bn - O8_HALF_WIDTH),
        _bp(_f15_be, _o8_ctr_bn + O8_HALF_WIDTH),
        _bp(_w15_be, _o8_ctr_bn + O8_HALF_WIDTH),
        _bp(_w15_be, _o8_ctr_bn - O8_HALF_WIDTH),
    ]))

    # O9, O10, O11: parametric positions from layout
    for _name, _ts, _te in [("O9",  layout.sw_t_o9_start,  layout.sw_t_o9_end),
                             ("O10", layout.sw_t_o10_start, layout.sw_t_o10_end),
                             ("O11", layout.sw_t_o11_start, layout.sw_t_o11_end)]:
        openings.append(OuterOpening(_name, "F20", "F1",
            _seg_opening("F20", "F1", _ts, _te)))

    return openings


def compute_rough_openings(pts, layout) -> list[RoughOpening]:
    """Compute all 7 interior rough-opening polygons.

    All openings are [SW, SE, NE, NW] polygons constructed using building axis
    vectors.  No axis-aligned BBox assumptions.
    """
    # Building axis vectors
    _w20, _w1 = pts["W20"], pts["W1"]
    _dE = _w1[0] - _w20[0]; _dN = _w1[1] - _w20[1]
    _seg_len = math.sqrt(_dE**2 + _dN**2)
    _along_E = _dE / _seg_len; _along_N = _dN / _seg_len
    _norm_E = _along_N; _norm_N = -_along_E
    _eE = -_along_E; _eN = -_along_N

    def _bn(p):
        """Building northing (perpendicular distance from south wall)."""
        return (p[0] - _w1[0]) * _norm_E + (p[1] - _w1[1]) * _norm_N

    def _be(p):
        """Building easting (distance east along south wall from W1)."""
        return (p[0] - _w1[0]) * _eE + (p[1] - _w1[1]) * _eN

    def _bp(be, bn):
        """World-coordinate point from building easting/northing."""
        return (_w1[0] + be * _eE + bn * _norm_E,
                _w1[1] + be * _eN + bn * _norm_N)

    # ── RO1: in IW1, E-W wall ───────────────────────────────────
    _iw1_s_bn = _bn(layout.iw1_poly[0])
    _iw1_n_bn = _bn(layout.iw1_poly[3])
    _iw2_e_be = _be(layout.iw2_poly[1])
    _ro1_w_be = _iw2_e_be + RO1_OFFSET_E_IW2
    _ro1_e_be = _ro1_w_be + IW1_RO_WIDTH
    ro1_poly = [
        _bp(_ro1_w_be, _iw1_s_bn),
        _bp(_ro1_e_be, _iw1_s_bn),
        _bp(_ro1_e_be, _iw1_n_bn),
        _bp(_ro1_w_be, _iw1_n_bn),
    ]

    # ── RO2: in IW11 (rotated), 3" NNE of IW12 north face along IW11 ─
    _iw11_se, _iw11_ne = layout.iw11_poly[1], layout.iw11_poly[2]
    _iw11_sw = layout.iw11_poly[0]
    _dx11 = _iw11_ne[0] - _iw11_se[0]
    _dy11 = _iw11_ne[1] - _iw11_se[1]
    _len11 = math.sqrt(_dx11**2 + _dy11**2)
    _un11 = (_dx11 / _len11, _dy11 / _len11)
    _iw12_nw = layout.iw12_poly[3]
    _ro2_start_d = ((_iw12_nw[0] - _iw11_se[0]) * _un11[0]
                    + (_iw12_nw[1] - _iw11_se[1]) * _un11[1]) + 3.0 / 12.0
    _ro2_end_d = _ro2_start_d + IW4_RO_WIDTH
    _ro2_sw = (_iw11_se[0] + _ro2_start_d * _un11[0],
               _iw11_se[1] + _ro2_start_d * _un11[1])
    _ro2_se = (_iw11_sw[0] + _ro2_start_d * _un11[0],
               _iw11_sw[1] + _ro2_start_d * _un11[1])
    _ro2_ne = (_iw11_sw[0] + _ro2_end_d * _un11[0],
               _iw11_sw[1] + _ro2_end_d * _un11[1])
    _ro2_nw = (_iw11_se[0] + _ro2_end_d * _un11[0],
               _iw11_se[1] + _ro2_end_d * _un11[1])
    ro2_poly = [_ro2_sw, _ro2_se, _ro2_ne, _ro2_nw]

    # ── RO3: in IW16, N-S wall, centered ─────────────────────────
    _iw16 = layout.iw16_poly
    _iw16_w_be = _be(_iw16[0])
    _iw16_e_be = _be(_iw16[1])
    _iw16_sw_bn = _bn(_iw16[0])
    _iw16_nw_bn = _bn(_iw16[3])
    _iw16_mid_bn = (_iw16_sw_bn + _iw16_nw_bn) / 2
    _ro3_s_bn = _iw16_mid_bn - IW16_RO_WIDTH / 2
    _ro3_n_bn = _iw16_mid_bn + IW16_RO_WIDTH / 2
    ro3_poly = [
        _bp(_iw16_w_be, _ro3_s_bn),
        _bp(_iw16_e_be, _ro3_s_bn),
        _bp(_iw16_e_be, _ro3_n_bn),
        _bp(_iw16_w_be, _ro3_n_bn),
    ]

    # ── RO4: in IW2, N-S wall ────────────────────────────────────
    _iw6_s_bn = _bn(layout.iw6_poly[0])
    _ro4_n_bn = _iw6_s_bn - IW2_RO_OFFSET_S
    _ro4_s_bn = _ro4_n_bn - IW2_RO_WIDTH
    _iw2_w_be = _be(layout.iw2_poly[0])
    ro4_poly = [
        _bp(_iw2_w_be, _ro4_s_bn),
        _bp(_iw2_e_be, _ro4_s_bn),
        _bp(_iw2_e_be, _ro4_n_bn),
        _bp(_iw2_w_be, _ro4_n_bn),
    ]

    # ── RO5: in IW6, E-W wall ────────────────────────────────────
    _ro5_e_be = _iw2_w_be - IW6_RO_OFFSET_W
    _ro5_w_be = _ro5_e_be - IW6_RO_WIDTH
    _iw6_n_bn = _bn(layout.iw6_poly[3])
    ro5_poly = [
        _bp(_ro5_w_be, _iw6_s_bn),
        _bp(_ro5_e_be, _iw6_s_bn),
        _bp(_ro5_e_be, _iw6_n_bn),
        _bp(_ro5_w_be, _iw6_n_bn),
    ]

    # ── RO6: in IW11 (rotated), centered between IW12 S face and IW11 S end ─
    _iw12_sw = layout.iw12_poly[0]
    _ro6_iw12_s_d = ((_iw12_sw[0] - _iw11_se[0]) * _un11[0]
                     + (_iw12_sw[1] - _iw11_se[1]) * _un11[1])
    _ro6_center_d = _ro6_iw12_s_d / 2
    _ro6_half = IW11_RO_WIDTH / 2
    _ro6_start_d = _ro6_center_d - _ro6_half
    _ro6_end_d = _ro6_center_d + _ro6_half
    _ro6_sw = (_iw11_se[0] + _ro6_start_d * _un11[0],
               _iw11_se[1] + _ro6_start_d * _un11[1])
    _ro6_se = (_iw11_sw[0] + _ro6_start_d * _un11[0],
               _iw11_sw[1] + _ro6_start_d * _un11[1])
    _ro6_ne = (_iw11_sw[0] + _ro6_end_d * _un11[0],
               _iw11_sw[1] + _ro6_end_d * _un11[1])
    _ro6_nw = (_iw11_se[0] + _ro6_end_d * _un11[0],
               _iw11_se[1] + _ro6_end_d * _un11[1])
    ro6_poly = [_ro6_sw, _ro6_se, _ro6_ne, _ro6_nw]

    # ── RO7: in IW9 (rotated), centered between IW7 S face and IW9 S end ─
    _iw9_se, _iw9_ne = layout.iw9_poly[1], layout.iw9_poly[2]
    _iw9_sw = layout.iw9_poly[0]
    _dx9 = _iw9_ne[0] - _iw9_se[0]
    _dy9 = _iw9_ne[1] - _iw9_se[1]
    _len9 = math.sqrt(_dx9**2 + _dy9**2)
    _un9 = (_dx9 / _len9, _dy9 / _len9)
    _iw7_sw = layout.iw7_poly[0]
    _ro7_iw7_s_d = ((_iw7_sw[0] - _iw9_se[0]) * _un9[0]
                    + (_iw7_sw[1] - _iw9_se[1]) * _un9[1])
    _ro7_center_d = _ro7_iw7_s_d / 2
    _ro7_half = IW9_RO_WIDTH / 2
    _ro7_start_d = _ro7_center_d - _ro7_half
    _ro7_end_d = _ro7_center_d + _ro7_half
    _ro7_sw = (_iw9_se[0] + _ro7_start_d * _un9[0],
               _iw9_se[1] + _ro7_start_d * _un9[1])
    _ro7_se = (_iw9_sw[0] + _ro7_start_d * _un9[0],
               _iw9_sw[1] + _ro7_start_d * _un9[1])
    _ro7_ne = (_iw9_sw[0] + _ro7_end_d * _un9[0],
               _iw9_sw[1] + _ro7_end_d * _un9[1])
    _ro7_nw = (_iw9_se[0] + _ro7_end_d * _un9[0],
               _iw9_se[1] + _ro7_end_d * _un9[1])
    ro7_poly = [_ro7_sw, _ro7_se, _ro7_ne, _ro7_nw]

    return [
        RoughOpening("RO1", ro1_poly, "IW1", "H"),
        RoughOpening("RO2", ro2_poly, "IW11", "R"),
        RoughOpening("RO3", ro3_poly, "IW16", "V"),
        RoughOpening("RO4", ro4_poly, "IW2", "V"),
        RoughOpening("RO5", ro5_poly, "IW6", "H"),
        RoughOpening("RO6", ro6_poly, "IW11", "R"),
        RoughOpening("RO7", ro7_poly, "IW9", "R"),
    ]


def _seg_param(pts, seg, point):
    """Compute parametric position t of a point along a LineSeg."""
    A = pts[seg.start]
    B = pts[seg.end]
    dx = B[0] - A[0]
    dy = B[1] - A[1]
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
