"""Opening position computation — outer wall openings and interior rough openings.

Single source of truth for all opening positions, consumed by both
gen_floorplan.py (polygon rendering) and gen_walls.py (parametric wall openings).
"""
import math
from typing import NamedTuple

from shared.types import Point, BBox, LineSeg
from floorplan.constants import (
    O1_OFFSET_S, O1_WIDTH, O2_OFFSET_S, O2_WIDTH,
    O3_HALF_WIDTH, O4_HALF_WIDTH,
    O5_E_FROM_F7, O5_WIDTH, O6_E_FROM_F9, O6_WIDTH,
    O7_NW_GAP, O7_HALF_WIDTH,
    O8_HALF_WIDTH, O9_HALF_WIDTH, O10_HALF_WIDTH, O11_HALF_WIDTH,
    IW1_OFFSET_N, IW5_OFFSET_N, WALL_3IN,
    KITCHEN_CTR_LENGTH, STD_GAP, FRIDGE_SIZE,
    IW1_RO_OFFSET_E, IW1_RO_WIDTH,
    IW2_RO_OFFSET_S, IW2_RO_WIDTH,
    IW3_RO_OFFSET_N, IW3_RO_WIDTH,
    IW4_RO_WIDTH, CLOSET1_HEIGHT, WALL_SOUTH_N,
    IW6_THICKNESS, IW6_OFFSET_N, IW6_RO_OFFSET_W, IW6_RO_WIDTH,
)


class OuterOpening(NamedTuple):
    """Opening in the outer wall, positioned as a 4-point polygon."""
    name: str
    seg_start: str       # e.g., "F1" — outline segment start point
    seg_end: str         # e.g., "F2" — outline segment end point
    poly: list[Point]    # 4 vertices: [outer_start, outer_end, inner_end, inner_start]


class RoughOpening(NamedTuple):
    """Rough opening in an interior wall."""
    name: str
    bbox: BBox           # w, s, e, n in survey coords
    wall_name: str       # "IW1", "IW2", etc.
    orientation: str     # "H" or "V"


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

    # O1: F1-F2, vertical, lower (south of IW1)
    o1_n = pts["F2"][1] - O1_OFFSET_S
    o1_s = o1_n - O1_WIDTH
    openings.append(OuterOpening("O1", "F1", "F2", [
        (pts["F2"][0], o1_s), (pts["F2"][0], o1_n),
        (pts["W2"][0], o1_n), (pts["W2"][0], o1_s),
    ]))

    # O2: F1-F2, vertical, upper (near F2)
    o2_n = pts["F2"][1] - O2_OFFSET_S
    o2_s = o2_n - O2_WIDTH
    openings.append(OuterOpening("O2", "F1", "F2", [
        (pts["F2"][0], o2_s), (pts["F2"][0], o2_n),
        (pts["W2"][0], o2_n), (pts["W2"][0], o2_s),
    ]))

    # O3: F4-F5, vertical, centered
    o3_cn = (pts["F4"][1] + pts["F5"][1]) / 2
    openings.append(OuterOpening("O3", "F4", "F5", [
        (pts["F4"][0], o3_cn - O3_HALF_WIDTH), (pts["F4"][0], o3_cn + O3_HALF_WIDTH),
        (pts["W4"][0], o3_cn + O3_HALF_WIDTH), (pts["W4"][0], o3_cn - O3_HALF_WIDTH),
    ]))

    # O4: F6-F7, horizontal, centered on midpoint
    o4_mid = (pts["F6"][0] + pts["F7"][0]) / 2
    o4_w = o4_mid - O4_HALF_WIDTH
    o4_e = o4_mid + O4_HALF_WIDTH
    openings.append(OuterOpening("O4", "F6", "F7", [
        (o4_w, pts["W6"][1]), (o4_e, pts["W6"][1]),
        (o4_e, pts["F6"][1]), (o4_w, pts["F6"][1]),
    ]))

    # O5: F9-F10, horizontal
    o5_e = pts["F7"][0] + O5_E_FROM_F7
    o5_w = o5_e - O5_WIDTH
    openings.append(OuterOpening("O5", "F9", "F10", [
        (o5_w, pts["W9"][1]), (o5_e, pts["W9"][1]),
        (o5_e, pts["F9"][1]), (o5_w, pts["F9"][1]),
    ]))

    # O6: F9-F10, horizontal
    o6_e = pts["F9"][0] + O6_E_FROM_F9
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
    iw5_n = pts["W0"][1] + IW1_OFFSET_N - IW5_OFFSET_N
    iw5_s = iw5_n - WALL_3IN
    o8_cn = (iw5_s + pts["F15"][1]) / 2
    openings.append(OuterOpening("O8", "F14", "F15", [
        (pts["F15"][0], o8_cn - O8_HALF_WIDTH), (pts["F15"][0], o8_cn + O8_HALF_WIDTH),
        (pts["W15"][0], o8_cn + O8_HALF_WIDTH), (pts["W15"][0], o8_cn - O8_HALF_WIDTH),
    ]))

    # O9: F18-F19, horizontal — centered between bed east and IW4 west
    o9_cn = (layout.bed.e + layout.iw4_w) / 2
    openings.append(OuterOpening("O9", "F18", "F19", [
        (o9_cn - O9_HALF_WIDTH, pts["F18"][1]), (o9_cn + O9_HALF_WIDTH, pts["F18"][1]),
        (o9_cn + O9_HALF_WIDTH, pts["W18"][1]), (o9_cn - O9_HALF_WIDTH, pts["W18"][1]),
    ]))

    # O10: F21-F0, horizontal (bed area) — centered between bed west and IW3 east
    o10_cn = (layout.bed.w + layout.iw3.e) / 2
    openings.append(OuterOpening("O10", "F21", "F0", [
        (o10_cn - O10_HALF_WIDTH, pts["F0"][1]), (o10_cn + O10_HALF_WIDTH, pts["F0"][1]),
        (o10_cn + O10_HALF_WIDTH, pts["W0"][1]), (o10_cn - O10_HALF_WIDTH, pts["W0"][1]),
    ]))

    # O11: F21-F0, horizontal (utility area) — centered between dryer and counter
    o11_cn = (layout.dryer.e + layout.ctr.w) / 2
    openings.append(OuterOpening("O11", "F21", "F0", [
        (o11_cn - O11_HALF_WIDTH, pts["F0"][1]), (o11_cn + O11_HALF_WIDTH, pts["F0"][1]),
        (o11_cn + O11_HALF_WIDTH, pts["W0"][1]), (o11_cn - O11_HALF_WIDTH, pts["W0"][1]),
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
    closet1_top = WALL_SOUTH_N + CLOSET1_HEIGHT
    iw8_n_face = closet1_top + WALL_3IN

    # RO1: in IW1, horizontal
    ro1_w = layout.iw2.e + KITCHEN_CTR_LENGTH + STD_GAP + FRIDGE_SIZE + IW1_RO_OFFSET_E
    ro1_e = ro1_w + IW1_RO_WIDTH

    # RO2: in IW4, vertical, centered between IW5 south and IW8 north
    ro2_center = (iw5_s + iw8_n_face) / 2
    ro2_s = ro2_center - IW4_RO_WIDTH / 2
    ro2_n = ro2_center + IW4_RO_WIDTH / 2

    # RO3: in IW3, vertical
    ro3_s = layout.ctr.n + WALL_3IN + IW3_RO_OFFSET_N
    ro3_n = ro3_s + IW3_RO_WIDTH

    # RO4: in IW2, vertical
    ro4_n = iw6_s - IW2_RO_OFFSET_S
    ro4_s = ro4_n - IW2_RO_WIDTH

    # RO5: in IW6, horizontal
    ro5_e = layout.iw2.w - IW6_RO_OFFSET_W
    ro5_w = ro5_e - IW6_RO_WIDTH

    return [
        RoughOpening("RO1", BBox(w=ro1_w, s=iw1_s, e=ro1_e, n=iw1_n), "IW1", "H"),
        RoughOpening("RO2", BBox(w=layout.iw4_w, s=ro2_s, e=layout.iw4_e, n=ro2_n), "IW4", "V"),
        RoughOpening("RO3", BBox(w=layout.iw3.w, s=ro3_s, e=layout.iw3.e, n=ro3_n), "IW3", "V"),
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
