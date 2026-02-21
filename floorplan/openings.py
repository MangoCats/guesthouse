"""Opening position computation — outer wall openings and interior rough openings.

Single source of truth for all opening positions, consumed by both
gen_floorplan.py (polygon rendering) and gen_walls.py (parametric wall openings).
"""
import math
from typing import NamedTuple

from shared.types import Point, BBox, LineSeg
from floorplan.constants import (
    O1_WIDTH, O2_WIDTH,
    IW2_RO_OFFSET_S, IW2_RO_WIDTH,
    O3_GAP_F5, O3_WIDTH, O4_HALF_WIDTH, O4_OFFSET_W_IW2,
    O5_E_FROM_IW2, O5_WIDTH, O6_WIDTH, O6_GAP_F10,
    O7_NW_GAP, O7_HALF_WIDTH,
    O8_HALF_WIDTH,
    IW5_OFFSET_N, WALL_3IN,
    STD_GAP,
    RO1_OFFSET_E_IW2, IW1_RO_WIDTH,
    IW2_RO_OFFSET_S, IW2_RO_WIDTH,
    IW4_RO_WIDTH, IW16_RO_WIDTH,
    IW6_THICKNESS, IW6_OFFSET_N, IW6_RO_OFFSET_W, IW6_RO_WIDTH,
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

    # O1: F2-F3, vertical, centered at RO3/IW16 center northing
    _iw16_ctr_n = (layout.iw16_poly[0][1] + layout.iw16_poly[2][1]) / 2
    o1_n = _iw16_ctr_n + O1_WIDTH / 2
    o1_s = _iw16_ctr_n - O1_WIDTH / 2
    openings.append(OuterOpening("O1", "F2", "F3", [
        (pts["F3"][0], o1_s), (pts["F3"][0], o1_n),
        (pts["W3"][0], o1_n), (pts["W3"][0], o1_s),
    ]))

    # O2: F4-F5, diagonal, centered at RO4 northing center
    _dE2 = pts["F5"][0] - pts["F4"][0]
    _dN2 = pts["F5"][1] - pts["F4"][1]
    _seg2_len = math.sqrt(_dE2**2 + _dN2**2)
    _ro4_ctr_n = layout.iw6_s - IW2_RO_OFFSET_S - IW2_RO_WIDTH / 2
    _t2_ctr = (_ro4_ctr_n - pts["F4"][1]) / _dN2
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
    _dE3 = pts["F5"][0] - pts["F4"][0]
    _dN3 = pts["F5"][1] - pts["F4"][1]
    _seg3_len = math.sqrt(_dE3**2 + _dN3**2)
    _t3_end = 1 - O3_GAP_F5 / _seg3_len       # closer to F5
    _t3_start = 1 - (O3_GAP_F5 + O3_WIDTH) / _seg3_len  # farther from F5
    openings.append(OuterOpening("O3", "F4", "F5", [
        (pts["F4"][0] + _t3_start * _dE3, pts["F4"][1] + _t3_start * _dN3),
        (pts["F4"][0] + _t3_end * _dE3, pts["F4"][1] + _t3_end * _dN3),
        (pts["W4"][0] + _t3_end * (pts["W5"][0] - pts["W4"][0]),
         pts["W4"][1] + _t3_end * (pts["W5"][1] - pts["W4"][1])),
        (pts["W4"][0] + _t3_start * (pts["W5"][0] - pts["W4"][0]),
         pts["W4"][1] + _t3_start * (pts["W5"][1] - pts["W4"][1])),
    ]))

    # O4: F6-F7, horizontal, relative to IW2 west face
    o4_mid = layout.iw2.w - O4_OFFSET_W_IW2
    o4_w = o4_mid - O4_HALF_WIDTH
    o4_e = o4_mid + O4_HALF_WIDTH
    openings.append(OuterOpening("O4", "F6", "F7", [
        (o4_w, pts["W6"][1]), (o4_e, pts["W6"][1]),
        (o4_e, pts["F6"][1]), (o4_w, pts["F6"][1]),
    ]))

    # O5: F9-F10, horizontal
    o5_e = layout.iw2.e + O5_E_FROM_IW2
    o5_w = o5_e - O5_WIDTH
    openings.append(OuterOpening("O5", "F9", "F10", [
        (o5_w, pts["W9"][1]), (o5_e, pts["W9"][1]),
        (o5_e, pts["F9"][1]), (o5_w, pts["F9"][1]),
    ]))

    # O6: F9-F10, horizontal, 6" west of F10
    o6_e = pts["F10"][0] - O6_GAP_F10
    o6_w = o6_e - O6_WIDTH
    openings.append(OuterOpening("O6", "F9", "F10", [
        (o6_w, pts["W9"][1]), (o6_e, pts["W9"][1]),
        (o6_e, pts["F9"][1]), (o6_w, pts["F9"][1]),
    ]))

    # O7: F12-F13, diagonal — NW end 2' from F12, 6' opening
    dE = pts["F13"][0] - pts["F12"][0]
    dN = pts["F13"][1] - pts["F12"][1]
    seg_len = math.sqrt(dE**2 + dN**2)
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

    # O8: F14-F15, vertical — centered between IW5 south face and F15
    o8_cn = (layout.iw5.s + pts["F15"][1]) / 2
    openings.append(OuterOpening("O8", "F14", "F15", [
        (pts["F15"][0], o8_cn - O8_HALF_WIDTH), (pts["F15"][0], o8_cn + O8_HALF_WIDTH),
        (pts["W15"][0], o8_cn + O8_HALF_WIDTH), (pts["W15"][0], o8_cn - O8_HALF_WIDTH),
    ]))

    # O9, O10, O11: F20-F1 — parametric positions from layout (single source)
    _dE9 = pts["F1"][0] - pts["F20"][0]
    _dN9 = pts["F1"][1] - pts["F20"][1]
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
    """Compute all 5 interior rough-opening bounding boxes."""
    iw1_s = layout.iw1_s
    iw1_n = layout.iw1_n
    iw6_n = pts["W6"][1] - IW6_OFFSET_N
    iw6_s = iw6_n - IW6_THICKNESS
    iw5_n = iw1_s - IW5_OFFSET_N
    iw5_s = iw5_n - WALL_3IN

    # RO1: in IW1, horizontal
    ro1_w = layout.iw2.e + RO1_OFFSET_E_IW2
    ro1_e = ro1_w + IW1_RO_WIDTH

    # RO2: in IW11 (rotated), 3" NNE of IW12 north face along IW11
    _iw11_se, _iw11_ne = layout.iw11_poly[1], layout.iw11_poly[2]
    _iw11_sw = layout.iw11_poly[0]
    _dx11 = _iw11_ne[0] - _iw11_se[0]
    _dy11 = _iw11_ne[1] - _iw11_se[1]
    _len11 = math.sqrt(_dx11**2 + _dy11**2)
    _un11 = (_dx11 / _len11, _dy11 / _len11)  # unit along IW11 length (NNE)
    _dx11t = _iw11_sw[0] - _iw11_se[0]
    _dy11t = _iw11_sw[1] - _iw11_se[1]
    _lt11 = math.sqrt(_dx11t**2 + _dy11t**2)
    _ut11 = (_dx11t / _lt11, _dy11t / _lt11)  # unit along IW11 thickness
    # IW12 NW corner projected onto IW11 length axis
    _iw12_nw = layout.iw12_poly[3]
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

    # RO3: in IW16 (axis-aligned N-S), 38" centered
    _iw16 = layout.iw16_poly  # [(w,s), (e,s), (e,n), (w,n)]
    _iw16_w = _iw16[0][0]
    _iw16_e = _iw16[1][0]
    _iw16_s = _iw16[0][1]
    _iw16_n = _iw16[2][1]
    _iw16_mid_n = (_iw16_s + _iw16_n) / 2
    ro3_s = _iw16_mid_n - IW16_RO_WIDTH / 2
    ro3_n = _iw16_mid_n + IW16_RO_WIDTH / 2

    # RO4: in IW2, vertical
    ro4_n = iw6_s - IW2_RO_OFFSET_S
    ro4_s = ro4_n - IW2_RO_WIDTH

    # RO5: in IW6, horizontal
    ro5_e = layout.iw2.w - IW6_RO_OFFSET_W
    ro5_w = ro5_e - IW6_RO_WIDTH

    return [
        RoughOpening("RO1", BBox(w=ro1_w, s=iw1_s, e=ro1_e, n=iw1_n), "IW1", "H"),
        RoughOpening("RO2", _ro2_bb, "IW11", "R", _ro2_poly),
        RoughOpening("RO3", BBox(w=_iw16_w, s=ro3_s, e=_iw16_e, n=ro3_n), "IW16", "V"),
        RoughOpening("RO4", BBox(w=layout.iw2.w, s=ro4_s, e=layout.iw2.e, n=ro4_n), "IW2", "V"),
        RoughOpening("RO5", BBox(w=ro5_w, s=iw6_s, e=ro5_e, n=iw6_n), "IW6", "H"),
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
