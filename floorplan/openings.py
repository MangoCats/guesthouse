"""Opening position computation — outer wall openings and interior rough openings.

Single source of truth for all opening positions, consumed by both
gen_floorplan.py (polygon rendering) and gen_walls.py (parametric wall openings).
"""
import math
from typing import NamedTuple

from shared.types import Point, BBox, LineSeg
from shared.geometry import seg_vec
from floorplan.constants import (
    O1_WIDTH, O2_WIDTH,
    O3_GAP_F5, O3_WIDTH, O4_HALF_WIDTH,
    O5_OFFSET_FROM_IW2, O5_WIDTH, O6_WIDTH, O6_GAP_F10,
    O7_NW_GAP, O7_HALF_WIDTH,
    O8_HALF_WIDTH,
    STD_GAP,
    RO1_OFFSET_FROM_IW2, IW1_RO_WIDTH,
    IW2_RO_OFFSET_FROM_IW6, IW2_RO_WIDTH,
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
    bbox: BBox           # w, s, e, n in survey coords (axis-aligned BB)
    wall_name: str       # "IW1", "IW2", etc.
    orientation: str     # "H", "V", or "R" (rotated)
    poly: list[Point] | None = None  # [SW, SE, NE, NW] for rotated openings


class WallOpening(NamedTuple):
    """Opening on an outline segment, parameterized along the segment."""
    name: str
    seg_idx: int    # index in outline_segs (0-based)
    t_start: float  # parametric position [0, 1] along the segment
    t_end: float    # parametric position [0, 1] along the segment


def compute_outer_openings(pts, layout) -> list[OuterOpening]:
    """Compute all 11 outer-wall opening polygons.

    Each polygon has 4 vertices spanning from the F-face (outer) to the W-face (inner).
    Returns openings in order: O1, O2, ..., O11.
    """
    openings = []

    # O1: F2-F3, centered at IW16 midpoint projected onto segment
    _dE1, _dN1, _seg1_len = seg_vec(pts["F2"], pts["F3"])
    _iw16_ctr = ((layout.iw16.poly[0][0] + layout.iw16.poly[2][0]) / 2,
                 (layout.iw16.poly[0][1] + layout.iw16.poly[2][1]) / 2)
    _t1_ctr = ((_iw16_ctr[0] - pts["F2"][0]) * _dE1
               + (_iw16_ctr[1] - pts["F2"][1]) * _dN1) / (_dE1**2 + _dN1**2)
    _t1_half = (O1_WIDTH / 2) / _seg1_len
    _t1_start = _t1_ctr - _t1_half
    _t1_end = _t1_ctr + _t1_half
    openings.append(OuterOpening("O1", "F2", "F3", [
        (pts["F2"][0] + _t1_start * _dE1, pts["F2"][1] + _t1_start * _dN1),
        (pts["F2"][0] + _t1_end * _dE1, pts["F2"][1] + _t1_end * _dN1),
        (pts["W2"][0] + _t1_end * (pts["W3"][0] - pts["W2"][0]),
         pts["W2"][1] + _t1_end * (pts["W3"][1] - pts["W2"][1])),
        (pts["W2"][0] + _t1_start * (pts["W3"][0] - pts["W2"][0]),
         pts["W2"][1] + _t1_start * (pts["W3"][1] - pts["W2"][1])),
    ]))

    # O2: F4-F5, diagonal, centered at RO4 position
    # Normal from RO4 center (perpendicular to IW2) intersects F4-F5 line.
    _dE2, _dN2, _seg2_len = seg_vec(pts["F4"], pts["F5"])
    # RO4 center: offset from IW6 south face midpoint in the IW6 outward direction
    _iw6_s_mid = ((layout.iw6.poly[0][0] + layout.iw6.poly[1][0]) / 2,
                  (layout.iw6.poly[0][1] + layout.iw6.poly[1][1]) / 2)
    _iw6_sdx = layout.iw6.poly[1][0] - layout.iw6.poly[0][0]
    _iw6_sdy = layout.iw6.poly[1][1] - layout.iw6.poly[0][1]
    _iw6_slen = math.sqrt(_iw6_sdx**2 + _iw6_sdy**2)
    _iw6_s_out = (_iw6_sdy / _iw6_slen, -_iw6_sdx / _iw6_slen)
    _ro4_off = IW2_RO_OFFSET_FROM_IW6 + IW2_RO_WIDTH / 2
    _ro4_ctr = (_iw6_s_mid[0] + _ro4_off * _iw6_s_out[0],
                _iw6_s_mid[1] + _ro4_off * _iw6_s_out[1])
    # IW2 normal: left normal of IW2 length direction (SW→NW)
    _iw2_ldx = layout.iw2.poly[3][0] - layout.iw2.poly[0][0]
    _iw2_ldy = layout.iw2.poly[3][1] - layout.iw2.poly[0][1]
    _iw2_llen = math.sqrt(_iw2_ldx**2 + _iw2_ldy**2)
    _iw2_norm = (-_iw2_ldy / _iw2_llen, _iw2_ldx / _iw2_llen)
    # Line-line intersection: RO4_ctr + s*IW2_norm ∩ F4 + t*(F5-F4)
    _cross_den = _dE2 * _iw2_norm[1] - _dN2 * _iw2_norm[0]
    _t2_ctr = ((_ro4_ctr[0] - pts["F4"][0]) * _iw2_norm[1]
               - (_ro4_ctr[1] - pts["F4"][1]) * _iw2_norm[0]) / _cross_den
    _t2_half = (O2_WIDTH / 2) / _seg2_len
    _t2_start = _t2_ctr - _t2_half
    _t2_end = _t2_ctr + _t2_half
    openings.append(OuterOpening("O2", "F4", "F5", [
        (pts["F4"][0] + _t2_start * _dE2, pts["F4"][1] + _t2_start * _dN2),
        (pts["F4"][0] + _t2_end * _dE2, pts["F4"][1] + _t2_end * _dN2),
        (pts["W4"][0] + _t2_end * (pts["W5"][0] - pts["W4"][0]),
         pts["W4"][1] + _t2_end * (pts["W5"][1] - pts["W4"][1])),
        (pts["W4"][0] + _t2_start * (pts["W5"][0] - pts["W4"][0]),
         pts["W4"][1] + _t2_start * (pts["W5"][1] - pts["W4"][1])),
    ]))

    # O3: F4-F5, diagonal, 4" from F5 along F5-F4 line
    # reuse _dE2, _dN2, _seg2_len from O2 (same segment F4-F5)
    _t3_end = 1 - O3_GAP_F5 / _seg2_len       # closer to F5
    _t3_start = 1 - (O3_GAP_F5 + O3_WIDTH) / _seg2_len  # farther from F5
    openings.append(OuterOpening("O3", "F4", "F5", [
        (pts["F4"][0] + _t3_start * _dE2, pts["F4"][1] + _t3_start * _dN2),
        (pts["F4"][0] + _t3_end * _dE2, pts["F4"][1] + _t3_end * _dN2),
        (pts["W4"][0] + _t3_end * (pts["W5"][0] - pts["W4"][0]),
         pts["W4"][1] + _t3_end * (pts["W5"][1] - pts["W4"][1])),
        (pts["W4"][0] + _t3_start * (pts["W5"][0] - pts["W4"][0]),
         pts["W4"][1] + _t3_start * (pts["W5"][1] - pts["W4"][1])),
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

    # O5: F9-F10, anchored at IW2 east face projection + offset
    _iw2_e_mid = ((layout.iw2.poly[1][0] + layout.iw2.poly[2][0]) / 2,
                  (layout.iw2.poly[1][1] + layout.iw2.poly[2][1]) / 2)
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

    # O9, O10, O11: F20-F1 — parametric positions from layout (single source)
    _dE9, _dN9, _ = seg_vec(pts["F20"], pts["F1"])
    for _name, _ts, _te in [("O9",  layout.sw_t_o9_start,  layout.sw_t_o9_end),
                             ("O10", layout.sw_t_o10_start, layout.sw_t_o10_end),
                             ("O11", layout.sw_t_o11_start, layout.sw_t_o11_end)]:
        openings.append(OuterOpening(_name, "F20", "F1", [
            (pts["F20"][0] + _ts * _dE9, pts["F20"][1] + _ts * _dN9),
            (pts["F20"][0] + _te * _dE9, pts["F20"][1] + _te * _dN9),
            (pts["W20"][0] + _te * (pts["W1"][0] - pts["W20"][0]),
             pts["W20"][1] + _te * (pts["W1"][1] - pts["W20"][1])),
            (pts["W20"][0] + _ts * (pts["W1"][0] - pts["W20"][0]),
             pts["W20"][1] + _ts * (pts["W1"][1] - pts["W20"][1])),
        ]))

    return openings


def compute_rough_openings(pts, layout) -> list[RoughOpening]:
    """Compute all 7 interior rough-opening polygons and bounding boxes."""

    # RO1: in IW1, positioned relative to IW2 east face along IW1 length
    _iw1_sw, _iw1_se = layout.iw1.poly[0], layout.iw1.poly[1]
    _iw1_nw = layout.iw1.poly[3]
    _dx1r, _dy1r, _len1r = seg_vec(_iw1_sw, _iw1_se)
    _un1_al = (_dx1r / _len1r, _dy1r / _len1r)  # unit along IW1 length
    # IW2 east face midpoint projected onto IW1 along-axis
    _iw2_e_mid_r = ((layout.iw2.poly[1][0] + layout.iw2.poly[2][0]) / 2,
                    (layout.iw2.poly[1][1] + layout.iw2.poly[2][1]) / 2)
    _ro1_ref_d = ((_iw2_e_mid_r[0] - _iw1_sw[0]) * _un1_al[0]
                  + (_iw2_e_mid_r[1] - _iw1_sw[1]) * _un1_al[1])
    _ro1_start_d = _ro1_ref_d + RO1_OFFSET_FROM_IW2
    _ro1_end_d = _ro1_start_d + IW1_RO_WIDTH
    _ro1_sw = (_iw1_sw[0] + _ro1_start_d * _un1_al[0],
               _iw1_sw[1] + _ro1_start_d * _un1_al[1])
    _ro1_se = (_iw1_sw[0] + _ro1_end_d * _un1_al[0],
               _iw1_sw[1] + _ro1_end_d * _un1_al[1])
    _ro1_ne = (_iw1_nw[0] + _ro1_end_d * _un1_al[0],
               _iw1_nw[1] + _ro1_end_d * _un1_al[1])
    _ro1_nw = (_iw1_nw[0] + _ro1_start_d * _un1_al[0],
               _iw1_nw[1] + _ro1_start_d * _un1_al[1])
    _ro1_poly = [_ro1_sw, _ro1_se, _ro1_ne, _ro1_nw]
    _ro1_bb = BBox(w=min(p[0] for p in _ro1_poly), s=min(p[1] for p in _ro1_poly),
                   e=max(p[0] for p in _ro1_poly), n=max(p[1] for p in _ro1_poly))

    # RO2: in IW11 (rotated), 3" NNE of IW12 north face along IW11
    _iw11_se, _iw11_ne = layout.iw11.poly[1], layout.iw11.poly[2]
    _iw11_sw = layout.iw11.poly[0]
    _dx11, _dy11, _len11 = seg_vec(_iw11_se, _iw11_ne)
    _un11 = (_dx11 / _len11, _dy11 / _len11)  # unit along IW11 length (NNE)
    _dx11t, _dy11t, _lt11 = seg_vec(_iw11_se, _iw11_sw)
    _ut11 = (_dx11t / _lt11, _dy11t / _lt11)  # unit along IW11 thickness
    # IW12 NW corner projected onto IW11 length axis
    _iw12_nw = layout.iw12.poly[3]
    _ro2_start_d = ((_iw12_nw[0] - _iw11_se[0]) * _un11[0]
                    + (_iw12_nw[1] - _iw11_se[1]) * _un11[1]) + 3.0 / 12.0
    _ro2_end_d = _ro2_start_d + IW4_RO_WIDTH
    # RO2 polygon [SW, SE, NE, NW] in IW11 coords
    _ro2_sw = (_iw11_se[0] + _ro2_start_d * _un11[0],
               _iw11_se[1] + _ro2_start_d * _un11[1])
    _ro2_se = (_iw11_sw[0] + _ro2_start_d * _un11[0],
               _iw11_sw[1] + _ro2_start_d * _un11[1])
    _ro2_ne = (_iw11_sw[0] + _ro2_end_d * _un11[0],
               _iw11_sw[1] + _ro2_end_d * _un11[1])
    _ro2_nw = (_iw11_se[0] + _ro2_end_d * _un11[0],
               _iw11_se[1] + _ro2_end_d * _un11[1])
    _ro2_poly = [_ro2_sw, _ro2_se, _ro2_ne, _ro2_nw]
    _ro2_bb = BBox(w=min(p[0] for p in _ro2_poly), s=min(p[1] for p in _ro2_poly),
                   e=max(p[0] for p in _ro2_poly), n=max(p[1] for p in _ro2_poly))

    # RO6: in IW11 (rotated), 62" centered between IW12 S face and IW11 S end
    # IW12 SW corner projected onto IW11 length axis = distance of IW12 south face
    _iw12_sw = layout.iw12.poly[0]
    _ro6_iw12_s_d = ((_iw12_sw[0] - _iw11_se[0]) * _un11[0]
                     + (_iw12_sw[1] - _iw11_se[1]) * _un11[1])
    _ro6_center_d = _ro6_iw12_s_d / 2  # centered between 0 (IW11 S end) and IW12 S face
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
    _ro6_poly = [_ro6_sw, _ro6_se, _ro6_ne, _ro6_nw]
    _ro6_bb = BBox(w=min(p[0] for p in _ro6_poly), s=min(p[1] for p in _ro6_poly),
                   e=max(p[0] for p in _ro6_poly), n=max(p[1] for p in _ro6_poly))

    # RO7: in IW9 (rotated), 62" centered between IW7 S face and IW9 S end
    _iw9_se, _iw9_ne = layout.iw9.poly[1], layout.iw9.poly[2]
    _iw9_sw = layout.iw9.poly[0]
    _dx9, _dy9, _len9 = seg_vec(_iw9_se, _iw9_ne)
    _un9 = (_dx9 / _len9, _dy9 / _len9)  # unit along IW9 length (NNE)
    # IW7 SW corner projected onto IW9 length axis = distance of IW7 south face
    _iw7_sw = layout.iw7.poly[0]
    _ro7_iw7_s_d = ((_iw7_sw[0] - _iw9_se[0]) * _un9[0]
                    + (_iw7_sw[1] - _iw9_se[1]) * _un9[1])
    _ro7_center_d = _ro7_iw7_s_d / 2  # centered between 0 (IW9 S end) and IW7 S face
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
    _ro7_poly = [_ro7_sw, _ro7_se, _ro7_ne, _ro7_nw]
    _ro7_bb = BBox(w=min(p[0] for p in _ro7_poly), s=min(p[1] for p in _ro7_poly),
                   e=max(p[0] for p in _ro7_poly), n=max(p[1] for p in _ro7_poly))

    # RO3: in IW16, centered along IW16 length
    _iw16_sw, _iw16_se = layout.iw16.poly[0], layout.iw16.poly[1]
    _iw16_nw = layout.iw16.poly[3]
    _dx16, _dy16, _len16 = seg_vec(_iw16_sw, _iw16_nw)
    _un16_al = (_dx16 / _len16, _dy16 / _len16)  # unit along IW16 length
    _ro3_center_d = _len16 / 2
    _ro3_half = IW16_RO_WIDTH / 2
    _ro3_start_d = _ro3_center_d - _ro3_half
    _ro3_end_d = _ro3_center_d + _ro3_half
    _ro3_sw = (_iw16_sw[0] + _ro3_start_d * _un16_al[0],
               _iw16_sw[1] + _ro3_start_d * _un16_al[1])
    _ro3_se = (_iw16_se[0] + _ro3_start_d * _un16_al[0],
               _iw16_se[1] + _ro3_start_d * _un16_al[1])
    _ro3_ne = (_iw16_se[0] + _ro3_end_d * _un16_al[0],
               _iw16_se[1] + _ro3_end_d * _un16_al[1])
    _ro3_nw = (_iw16_sw[0] + _ro3_end_d * _un16_al[0],
               _iw16_sw[1] + _ro3_end_d * _un16_al[1])
    _ro3_poly = [_ro3_sw, _ro3_se, _ro3_ne, _ro3_nw]
    _ro3_bb = BBox(w=min(p[0] for p in _ro3_poly), s=min(p[1] for p in _ro3_poly),
                   e=max(p[0] for p in _ro3_poly), n=max(p[1] for p in _ro3_poly))

    # RO4: in IW2, positioned relative to IW6 south face along IW2 length
    _iw2_sw, _iw2_se = layout.iw2.poly[0], layout.iw2.poly[1]
    _iw2_nw = layout.iw2.poly[3]
    _dx2r, _dy2r, _len2r = seg_vec(_iw2_sw, _iw2_nw)
    _un2_al = (_dx2r / _len2r, _dy2r / _len2r)  # unit along IW2 length
    # IW6 south face midpoint projected onto IW2 along-axis
    _iw6_s_mid = ((layout.iw6.poly[0][0] + layout.iw6.poly[1][0]) / 2,
                  (layout.iw6.poly[0][1] + layout.iw6.poly[1][1]) / 2)
    _ro4_ref_d = ((_iw6_s_mid[0] - _iw2_sw[0]) * _un2_al[0]
                  + (_iw6_s_mid[1] - _iw2_sw[1]) * _un2_al[1])
    _ro4_end_d = _ro4_ref_d - IW2_RO_OFFSET_FROM_IW6
    _ro4_start_d = _ro4_end_d - IW2_RO_WIDTH
    _ro4_sw = (_iw2_sw[0] + _ro4_start_d * _un2_al[0],
               _iw2_sw[1] + _ro4_start_d * _un2_al[1])
    _ro4_se = (_iw2_se[0] + _ro4_start_d * _un2_al[0],
               _iw2_se[1] + _ro4_start_d * _un2_al[1])
    _ro4_ne = (_iw2_se[0] + _ro4_end_d * _un2_al[0],
               _iw2_se[1] + _ro4_end_d * _un2_al[1])
    _ro4_nw = (_iw2_sw[0] + _ro4_end_d * _un2_al[0],
               _iw2_sw[1] + _ro4_end_d * _un2_al[1])
    _ro4_poly = [_ro4_sw, _ro4_se, _ro4_ne, _ro4_nw]
    _ro4_bb = BBox(w=min(p[0] for p in _ro4_poly), s=min(p[1] for p in _ro4_poly),
                   e=max(p[0] for p in _ro4_poly), n=max(p[1] for p in _ro4_poly))

    # RO5: in IW6, positioned relative to IW2 west face
    # IW6 is trapezoidal (west end meets inner polygon at slightly different
    # points on each face), so compute each face independently.
    _iw6_sw, _iw6_se = layout.iw6.poly[0], layout.iw6.poly[1]
    _iw6_ne, _iw6_nw = layout.iw6.poly[2], layout.iw6.poly[3]
    _iw2_w_mid = ((layout.iw2.poly[0][0] + layout.iw2.poly[3][0]) / 2,
                  (layout.iw2.poly[0][1] + layout.iw2.poly[3][1]) / 2)
    # South face (SW→SE)
    _dx6s, _dy6s, _len6s = seg_vec(_iw6_sw, _iw6_se)
    _un6s = (_dx6s / _len6s, _dy6s / _len6s)
    _ref6s = ((_iw2_w_mid[0] - _iw6_sw[0]) * _un6s[0]
              + (_iw2_w_mid[1] - _iw6_sw[1]) * _un6s[1])
    _end6s = _ref6s - IW6_RO_OFFSET_W
    _start6s = _end6s - IW6_RO_WIDTH
    # North face (NW→NE)
    _dx6n, _dy6n, _len6n = seg_vec(_iw6_nw, _iw6_ne)
    _un6n = (_dx6n / _len6n, _dy6n / _len6n)
    _ref6n = ((_iw2_w_mid[0] - _iw6_nw[0]) * _un6n[0]
              + (_iw2_w_mid[1] - _iw6_nw[1]) * _un6n[1])
    _end6n = _ref6n - IW6_RO_OFFSET_W
    _start6n = _end6n - IW6_RO_WIDTH
    _ro5_sw = (_iw6_sw[0] + _start6s * _un6s[0],
               _iw6_sw[1] + _start6s * _un6s[1])
    _ro5_se = (_iw6_sw[0] + _end6s * _un6s[0],
               _iw6_sw[1] + _end6s * _un6s[1])
    _ro5_ne = (_iw6_nw[0] + _end6n * _un6n[0],
               _iw6_nw[1] + _end6n * _un6n[1])
    _ro5_nw = (_iw6_nw[0] + _start6n * _un6n[0],
               _iw6_nw[1] + _start6n * _un6n[1])
    _ro5_poly = [_ro5_sw, _ro5_se, _ro5_ne, _ro5_nw]
    _ro5_bb = BBox(w=min(p[0] for p in _ro5_poly), s=min(p[1] for p in _ro5_poly),
                   e=max(p[0] for p in _ro5_poly), n=max(p[1] for p in _ro5_poly))

    return [
        RoughOpening("RO1", _ro1_bb, "IW1", "H", _ro1_poly),
        RoughOpening("RO2", _ro2_bb, "IW11", "R", _ro2_poly),
        RoughOpening("RO3", _ro3_bb, "IW16", "V", _ro3_poly),
        RoughOpening("RO4", _ro4_bb, "IW2", "V", _ro4_poly),
        RoughOpening("RO5", _ro5_bb, "IW6", "H", _ro5_poly),
        RoughOpening("RO6", _ro6_bb, "IW11", "R", _ro6_poly),
        RoughOpening("RO7", _ro7_bb, "IW9", "R", _ro7_poly),
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
