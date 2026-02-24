"""Interior layout computation — rooms, walls, appliances, furniture."""
import math
from typing import NamedTuple

from shared.types import Point, Wall
from shared.geometry import line_isect, seg_vec, seg_vecs, offset_pt
from floorplan.constants import (
    WALL_6IN, WALL_4IN, WALL_3IN,
    APPLIANCE_WIDTH, APPLIANCE_DEPTH, APPLIANCE_OFFSET_FROM_W2,
    APPLIANCE_OFFSET_FROM_W1, APPLIANCE_GAP,
    COUNTER_DEPTH, COUNTER_GAP, COUNTER_LENGTH,
    BED_WIDTH, BED_LENGTH, BED_WALL_GAP,
    O9_HALF_WIDTH, O10_HALF_WIDTH, O11_HALF_WIDTH,
    O9_OFFSET_IW11, O9_O10_WALL, O10_O11_WALL, BED_GAP_O9,
    IW1_OFFSET_FROM_W9, IW1_OFFSET_FROM_W2, IW2_OFFSET_FROM_W2,
    IW3_LENGTH, IW3_OFFSET_IW9, IW3_DIST_W2W5,
    IW9_LENGTH,
    IW4_OFFSET_FROM_IW2,
    IW5_OFFSET_FROM_IW1, IW6_THICKNESS, IW6_OFFSET_FROM_W6,
    IW8_OFFSET_FROM_W20W1,
    IW11_RADIUS_FROM_IW4, IW9_IW11_GAP, IW12_OFFSET_IW11, IW12_SHORTEN,
    DRESSER_WIDTH, DRESSER_DEPTH, DRESSER_GAP_IW15, DRESSER_GAP_IW1,
)


_seg_vecs = seg_vecs
_offset = offset_pt


def _dot(a, b):
    """Dot product of two 2D vectors."""
    return a[0] * b[0] + a[1] * b[1]


def _proj(target, anchor, direction):
    """Signed projection of (target - anchor) onto direction vector."""
    return _dot((target[0] - anchor[0], target[1] - anchor[1]), direction)


def _make_wall(poly):
    """Create Wall from polygon corners, computing axis-aligned bounding box."""
    xs = [p[0] for p in poly]
    ys = [p[1] for p in poly]
    return Wall(poly=poly, w=min(xs), s=min(ys), e=max(xs), n=max(ys))


def _line_poly_isects(poly, origin, direction):
    """Parametric t values where ray origin+t*direction crosses polygon edges."""
    results = []
    n = len(poly)
    for i in range(n):
        j = (i + 1) % n
        ex, ey = poly[j][0] - poly[i][0], poly[j][1] - poly[i][1]
        det = direction[0] * ey - direction[1] * ex
        if abs(det) < 1e-12:
            continue
        ox, oy = poly[i][0] - origin[0], poly[i][1] - origin[1]
        t = (ox * ey - oy * ex) / det
        s = (ox * direction[1] - oy * direction[0]) / det
        if 0 <= s <= 1:
            results.append(t)
    return results


class InteriorLayout(NamedTuple):
    """Interior layout positions for walls, appliances, and furniture."""
    # Interior walls
    iw1: Wall
    iw2: Wall
    iw3: Wall
    iw4: Wall
    iw5: Wall
    iw6: Wall
    iw7: Wall
    iw8: Wall
    iw9: Wall
    iw11: Wall
    iw12: Wall
    iw16: Wall
    # Appliances & furniture
    dryer: Wall
    washer: Wall
    ctr: Wall              # rectangular counter bounds
    ctr_clip: list[Point]  # clipped rendering polygon (5 vertices)
    ctr_nw_r: float        # NW corner radius (currently 0)
    dresser: Wall
    bed: Wall
    # South wall opening parametric positions
    sw_t_o9_start: float
    sw_t_o9_end: float
    sw_t_o10_start: float
    sw_t_o10_end: float
    sw_t_o11_start: float
    sw_t_o11_end: float


def _compute_segment_vectors(pts):
    """Compute along/inward unit vectors for the 4 reference wall segments.

    Returns ((w2w5_al, w2w5_in), (w6w7_al, w6w7_in),
             (w9w10_al, w9w10_in), (w20w1_al, w20w1_in)).
    """
    return (
        _seg_vecs(pts["W2"], pts["W5"]),    # along≈(0,1), inward≈(1,0)
        _seg_vecs(pts["W6"], pts["W7"]),    # along≈(1,0), inward≈(0,-1)
        _seg_vecs(pts["W9"], pts["W10"]),   # along≈(1,0), inward≈(0,-1)
        _seg_vecs(pts["W20"], pts["W1"]),
    )


def _compute_south_walls(pts, w20w1, w2w5, iw4_sw, iw11_sw, iw11_se):
    """Compute walls along the south (W20-W1) wall: bed, IW9, IW3, IW7.

    Returns (bed, iw9, iw3, iw7, opening_params) where opening_params is
    (ts9, te9, ts10, te10, ts11, te11) — parametric positions on F20-F1.
    """
    w20w1_al, w20w1_in = w20w1
    w2w5_al, _ = w2w5
    _dE, _dN, _seg_len = seg_vec(pts["W20"], pts["W1"])

    # --- Bed ---
    _dE9, _dN9, _seg9_len = seg_vec(pts["F20"], pts["F1"])
    _t_sw9 = ((iw11_sw[0] - pts["F20"][0]) * _dE9
              + (iw11_sw[1] - pts["F20"][1]) * _dN9) / (_dE9**2 + _dN9**2)
    _ts9 = _t_sw9 + O9_OFFSET_IW11 / _seg9_len
    _te9 = _ts9 + 2 * O9_HALF_WIDTH / _seg9_len
    _bed_t = _te9 + BED_GAP_O9 / _seg_len
    _bed_se_wall = _offset(pts["W20"], _bed_t, (_dE, _dN))
    bed_se = _offset(_bed_se_wall, BED_WALL_GAP, w20w1_in)
    bed_sw = _offset(bed_se, BED_WIDTH, w20w1_al)
    bed_ne = _offset(bed_se, BED_LENGTH, w20w1_in)
    bed_nw = _offset(bed_sw, BED_LENGTH, w20w1_in)

    # --- Openings O10, O11 parametric positions ---
    _ts10 = _te9 + O9_O10_WALL / _seg9_len
    _te10 = _ts10 + 2 * O10_HALF_WIDTH / _seg9_len
    _ts11 = _te10 + O10_O11_WALL / _seg9_len
    _te11 = _ts11 + 2 * O11_HALF_WIDTH / _seg9_len

    # --- IW3 (anchored at IW3_DIST_W2W5 from W2-W5) ---
    _w2w5_ref = line_isect(pts["W2"], w2w5_al, pts["W1"], w20w1_al)
    iw3_sw = _offset(_w2w5_ref, -IW3_DIST_W2W5, w20w1_al)
    iw3_se = _offset(iw3_sw, -WALL_4IN, w20w1_al)
    iw3_ne = _offset(iw3_se, IW3_LENGTH, w20w1_in)
    iw3_nw = _offset(iw3_sw, IW3_LENGTH, w20w1_in)

    # --- IW9 (IW3_OFFSET_IW9 from IW3 east face, toward W20) ---
    iw9_sw = _offset(iw3_se, -IW3_OFFSET_IW9, w20w1_al)
    iw9_se = _offset(iw9_sw, -WALL_4IN, w20w1_al)
    iw9_ne = _offset(iw9_se, IW9_LENGTH, w20w1_in)
    iw9_nw = _offset(iw9_sw, IW9_LENGTH, w20w1_in)

    # --- IW7 ---
    iw7_nw = iw3_ne
    iw7_ne = _offset(iw7_nw, -IW3_OFFSET_IW9, w20w1_al)
    iw7_sw = _offset(iw7_nw, -WALL_4IN, w20w1_in)
    iw7_se = _offset(iw7_ne, -WALL_4IN, w20w1_in)

    return (
        _make_wall([bed_sw, bed_se, bed_ne, bed_nw]),
        _make_wall([iw9_sw, iw9_se, iw9_ne, iw9_nw]),
        _make_wall([iw3_sw, iw3_se, iw3_ne, iw3_nw]),
        _make_wall([iw7_sw, iw7_se, iw7_ne, iw7_nw]),
        (_ts9, _te9, _ts10, _te10, _ts11, _te11),
    )


def compute_interior_layout(pts, inner_poly) -> InteriorLayout:
    """Compute interior layout positions.

    pts must contain W-series (W1-W20) and F-series (F1-F20).
    """
    # Segment vectors for wall segments
    ((_w2w5_al, _w2w5_in), (_w6w7_al, _w6w7_in),
     (_w9w10_al, _w9w10_in), (_w20w1_al, _w20w1_in)) = _compute_segment_vectors(pts)
    _neg_w20w1_al = (-_w20w1_al[0], -_w20w1_al[1])
    _neg_w6w7_al = (-_w6w7_al[0], -_w6w7_al[1])

    # Old east-wall direction (perpendicular to W6-W7, unchanged by F2 move).
    # All IW positioning uses these to preserve interior layout positions.
    _iw_al = (-_w6w7_al[1], _w6w7_al[0])    # ≈ old W2-W3 along (roughly north)
    _iw_in = _w6w7_al                         # ≈ old W2-W3 inward (roughly east)
    # Virtual W2 reference: 8' from W7 toward W6, preserving old IW offsets
    _iw_w2_ref = _offset(pts["W7"], 8.0, _neg_w6w7_al)

    # --- IW1 ---
    _iw1_n_anchor = _offset(pts["W9"], IW1_OFFSET_FROM_W9, _w9w10_in)
    _iw1_w_anchor = _offset(_iw_w2_ref, IW1_OFFSET_FROM_W2, _iw_in)
    iw1_nw = line_isect(_iw1_n_anchor, _w9w10_al, _iw1_w_anchor, _iw_al)
    iw1_sw = _offset(iw1_nw, WALL_6IN, _w9w10_in)
    _ne_ts = _line_poly_isects(inner_poly, iw1_nw, _w9w10_al)
    iw1_ne = _offset(iw1_nw, max(_ne_ts), _w9w10_al)
    _se_ts = _line_poly_isects(inner_poly, iw1_sw, _w9w10_al)
    iw1_se = _offset(iw1_sw, max(_se_ts), _w9w10_al)

    # --- IW2 ---
    _iw2_w_anchor = _offset(_iw_w2_ref, IW2_OFFSET_FROM_W2, _iw_in)
    _iw2_e_anchor = _offset(_iw2_w_anchor, WALL_6IN, _iw_in)
    iw2_sw = line_isect(_iw2_w_anchor, _iw_al, iw1_nw, _w9w10_al)
    iw2_se = line_isect(_iw2_e_anchor, _iw_al, iw1_nw, _w9w10_al)
    iw2_nw = line_isect(_iw2_w_anchor, _iw_al, pts["W6"], _w6w7_al)
    iw2_ne = line_isect(_iw2_e_anchor, _iw_al, pts["W6"], _w6w7_al)

    # --- Dryer (aligned to new W2-W5 east wall) ---
    _dryer_sw = line_isect(
        _offset(pts["W2"], APPLIANCE_OFFSET_FROM_W2, _w2w5_in), _w2w5_al,
        _offset(pts["W1"], APPLIANCE_OFFSET_FROM_W1, _w2w5_al), _w9w10_al)
    _dryer_se = _offset(_dryer_sw, APPLIANCE_WIDTH, _w2w5_in)
    _dryer_nw = _offset(_dryer_sw, APPLIANCE_DEPTH, _w2w5_al)
    _dryer_ne = _offset(_dryer_se, APPLIANCE_DEPTH, _w2w5_al)

    # --- Washer (aligned to new W2-W5 east wall) ---
    _washer_sw = _offset(_dryer_nw, APPLIANCE_GAP, _w2w5_al)
    _washer_se = _offset(_washer_sw, APPLIANCE_WIDTH, _w2w5_in)
    _washer_nw = _offset(_washer_sw, APPLIANCE_DEPTH, _w2w5_al)
    _washer_ne = _offset(_washer_se, APPLIANCE_DEPTH, _w2w5_al)

    # --- Counter (aligned to new W2-W5 east wall) ---
    _ctr_sw_anchor = _offset(_dryer_se, COUNTER_GAP, _w2w5_in)
    _ctr_se_anchor = _offset(_ctr_sw_anchor, COUNTER_DEPTH, _w2w5_in)
    ctr_sw = line_isect(_ctr_sw_anchor, _w2w5_al, pts["W1"], _w9w10_al)
    ctr_nw = _offset(ctr_sw, COUNTER_LENGTH, _w2w5_al)
    ctr_se = line_isect(_ctr_se_anchor, _w2w5_al, pts["W1"], _w9w10_al)
    ctr_ne = _offset(ctr_se, COUNTER_LENGTH, _w2w5_al)
    ctr_nw_r = 0

    # --- IW4 ---
    _iw4_w_anchor = _offset(iw2_se, IW4_OFFSET_FROM_IW2, _iw_in)
    _iw4_e_anchor = _offset(_iw4_w_anchor, WALL_4IN, _iw_in)
    _iw4_t_s = _proj(pts["W19"], _iw4_w_anchor, _iw_al)
    iw4_sw = _offset(_iw4_w_anchor, _iw4_t_s, _iw_al)
    iw4_se = _offset(_iw4_e_anchor, _iw4_t_s, _iw_al)

    # Shift IW4/IW11 group along W20-W1 for IW9-IW11 gap constraint
    _w20 = pts["W20"]
    _w2w5_ref_s = line_isect(pts["W2"], _w2w5_al, pts["W1"], _w20w1_al)
    _iw9_se_pre = _offset(_w2w5_ref_s,
                           -(IW3_DIST_W2W5 + 2 * WALL_4IN + IW3_OFFSET_IW9),
                           _w20w1_al)
    _h = iw4_sw[1] - _w20[1]
    _circle_dx = math.sqrt(IW11_RADIUS_FROM_IW4**2 - _h**2)
    _iw_shift = (iw4_sw[0] - _circle_dx - WALL_4IN) - (_iw9_se_pre[0] + IW9_IW11_GAP)
    _iw4_w_anchor = _offset(_iw4_w_anchor, _iw_shift, _w20w1_al)
    _iw4_e_anchor = _offset(_iw4_e_anchor, _iw_shift, _w20w1_al)
    iw4_sw = _offset(iw4_sw, _iw_shift, _w20w1_al)
    iw4_se = _offset(iw4_se, _iw_shift, _w20w1_al)
    # IW4 north end computed after IW12 (depends on IW12 north face)

    # --- IW11 ---
    # SE corner: circle(IW4_SW, IW11_RADIUS_FROM_IW4) intersect W20-W1
    _w1 = pts["W1"]
    _dE, _dN, _seg_len = seg_vec(_w20, _w1)
    _uE = _w20[0] - iw4_sw[0]
    _uN = _w20[1] - iw4_sw[1]
    _qa = _dE**2 + _dN**2
    _qb = 2 * (_uE * _dE + _uN * _dN)
    _qc = _uE**2 + _uN**2 - IW11_RADIUS_FROM_IW4**2
    _disc = _qb**2 - 4 * _qa * _qc
    _t = (-_qb + math.sqrt(_disc)) / (2 * _qa)
    iw11_se = (_w20[0] + _t * _dE, _w20[1] + _t * _dN)
    iw11_sw = _offset(iw11_se, WALL_4IN, _w20w1_al)
    iw11_ne = line_isect(iw11_se, _w20w1_in, iw1_sw, _w9w10_al)
    iw11_nw = line_isect(iw11_sw, _w20w1_in, iw1_sw, _w9w10_al)

    # --- IW12 ---
    _iw12_base = _offset(iw11_sw, IW12_OFFSET_IW11, _w20w1_in)
    iw12_sw = _offset(_iw12_base, -IW12_SHORTEN, _w20w1_al)
    iw12_nw = _offset(iw12_sw, WALL_4IN, _w20w1_in)
    iw12_se = line_isect(iw12_sw, _neg_w20w1_al, iw4_sw, _iw_al)
    iw12_ne = line_isect(iw12_nw, _neg_w20w1_al, iw4_sw, _iw_al)

    # IW4 north end: truncate at IW12 north face
    _iw12_n_dir = (iw12_ne[0] - iw12_nw[0], iw12_ne[1] - iw12_nw[1])
    iw4_nw = line_isect(_iw4_w_anchor, _iw_al, iw12_nw, _iw12_n_dir)
    iw4_ne = line_isect(_iw4_e_anchor, _iw_al, iw12_nw, _iw12_n_dir)

    # --- Bed, IW9, IW3, IW7 (south wall chain) ---
    (bed, iw9, iw3, iw7,
     (_ts9, _te9, _ts10, _te10, _ts11, _te11)) = _compute_south_walls(
        pts, (_w20w1_al, _w20w1_in), (_w2w5_al, _w2w5_in),
        iw4_sw, iw11_sw, iw11_se)
    iw3_sw, iw3_nw = iw3.poly[0], iw3.poly[3]

    # --- IW16 ---
    _iw16_w_anchor = iw3_nw
    _iw16_e_anchor = _offset(_iw16_w_anchor, WALL_4IN, _iw_in)
    _iw16_t_n = _proj(iw1_sw, _iw16_w_anchor, _iw_al)
    iw16_nw = _offset(_iw16_w_anchor, _iw16_t_n, _iw_al)
    iw16_ne = _offset(_iw16_e_anchor, _iw16_t_n, _iw_al)

    # Counter clip polygon (south edge follows W20-W1, east clipped at IW3/IW16)
    _ctr_poly_sw = line_isect(_ctr_sw_anchor, _w2w5_al, pts["W20"], _w20w1_al)
    _iw3_w_dir = _w20w1_in
    _ctr_vs_iw3 = _proj(ctr_nw, iw3_nw, _w2w5_al)
    if _ctr_vs_iw3 > 0:
        _ctr_ne_clip = line_isect(ctr_nw, _w2w5_in, _iw16_w_anchor, _iw_al)
        ctr_clip = [ctr_nw, _ctr_ne_clip, iw3_nw, iw3_sw, _ctr_poly_sw]
    else:
        _ctr_ne_clip = line_isect(ctr_nw, _w2w5_in, iw3_sw, _iw3_w_dir)
        ctr_clip = [ctr_nw, _ctr_ne_clip, iw3_sw, _ctr_poly_sw]

    # --- IW8 (perpendicular to W2-W5, 12' from W20-W1; west at W2-W5, east at IW1) ---
    _iw8_s_anchor = _offset(pts["W20"], IW8_OFFSET_FROM_W20W1, _w20w1_in)
    _iw8_n_anchor = _offset(_iw8_s_anchor, WALL_6IN, _w2w5_al)
    iw8_sw = line_isect(_iw8_s_anchor, _w2w5_in, pts["W2"], _w2w5_al)
    iw8_nw = line_isect(_iw8_n_anchor, _w2w5_in, pts["W2"], _w2w5_al)
    iw8_se = line_isect(_iw8_s_anchor, _w2w5_in, iw1_nw, _iw_al)
    iw8_ne = line_isect(_iw8_n_anchor, _w2w5_in, iw1_nw, _iw_al)

    # --- IW5 ---
    _iw5_n_anchor = _offset(iw1_sw, IW5_OFFSET_FROM_IW1, _w9w10_in)
    _iw5_s_anchor = _offset(_iw5_n_anchor, WALL_3IN, _w9w10_in)
    _iw5_t_e = _proj(pts["W15"], _iw5_n_anchor, _w9w10_al)
    iw5_ne = _offset(_iw5_n_anchor, _iw5_t_e, _w9w10_al)
    iw5_se = _offset(_iw5_s_anchor, _iw5_t_e, _w9w10_al)

    # IW5 west end (at IW11 east face)
    iw5_nw = line_isect(_iw5_n_anchor, _w9w10_al, iw11_se, _w20w1_in)
    iw5_sw = line_isect(_iw5_s_anchor, _w9w10_al, iw11_se, _w20w1_in)

    # --- Dresser ---
    dresser_ne = line_isect(
        _offset(iw11_nw, -DRESSER_GAP_IW15, _iw_in), _iw_al,
        _offset(iw1_sw, DRESSER_GAP_IW1, _w9w10_in), _w9w10_al)
    dresser_nw = _offset(dresser_ne, -DRESSER_WIDTH, _iw_in)
    dresser_se = _offset(dresser_ne, DRESSER_DEPTH, _w9w10_in)
    dresser_sw = _offset(dresser_nw, DRESSER_DEPTH, _w9w10_in)

    # --- IW6 ---
    _iw6_n_anchor = _offset(pts["W6"], IW6_OFFSET_FROM_W6, _w6w7_in)
    _iw6_s_anchor = _offset(_iw6_n_anchor, IW6_THICKNESS, _w6w7_in)
    iw6_ne = line_isect(_iw6_n_anchor, _w6w7_al, _iw2_w_anchor, _iw_al)
    iw6_se = line_isect(_iw6_s_anchor, _w6w7_al, _iw2_w_anchor, _iw_al)
    _iw6_n_ts = _line_poly_isects(inner_poly, iw6_ne, _neg_w6w7_al)
    iw6_nw = _offset(iw6_ne, max(_iw6_n_ts), _neg_w6w7_al)
    _iw6_s_ts = _line_poly_isects(inner_poly, iw6_se, _neg_w6w7_al)
    iw6_sw = _offset(iw6_se, max(_iw6_s_ts), _neg_w6w7_al)

    return InteriorLayout(
        iw1=_make_wall([iw1_sw, iw1_se, iw1_ne, iw1_nw]),
        iw2=_make_wall([iw2_sw, iw2_se, iw2_ne, iw2_nw]),
        iw3=iw3,
        iw4=_make_wall([iw4_sw, iw4_se, iw4_ne, iw4_nw]),
        iw5=_make_wall([iw5_sw, iw5_se, iw5_ne, iw5_nw]),
        iw6=_make_wall([iw6_sw, iw6_se, iw6_ne, iw6_nw]),
        iw7=iw7,
        iw8=_make_wall([iw8_sw, iw8_se, iw8_ne, iw8_nw]),
        iw9=iw9,
        iw11=_make_wall([iw11_sw, iw11_se, iw11_ne, iw11_nw]),
        iw12=_make_wall([iw12_sw, iw12_se, iw12_ne, iw12_nw]),
        iw16=_make_wall([_iw16_w_anchor, _iw16_e_anchor, iw16_ne, iw16_nw]),
        dryer=_make_wall([_dryer_sw, _dryer_se, _dryer_ne, _dryer_nw]),
        washer=_make_wall([_washer_sw, _washer_se, _washer_ne, _washer_nw]),
        ctr=_make_wall([ctr_sw, ctr_se, ctr_ne, ctr_nw]),
        ctr_clip=ctr_clip,
        ctr_nw_r=ctr_nw_r,
        dresser=_make_wall([dresser_sw, dresser_se, dresser_ne, dresser_nw]),
        bed=bed,
        sw_t_o9_start=_ts9, sw_t_o9_end=_te9,
        sw_t_o10_start=_ts10, sw_t_o10_end=_te10,
        sw_t_o11_start=_ts11, sw_t_o11_end=_te11,
    )
