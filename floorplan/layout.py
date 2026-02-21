"""Interior layout computation — rooms, walls, appliances, furniture.

All walls/appliances/furniture use polygon [SW, SE, NE, NW] representations
constructed from building axis vectors (derived from W20→W1 south wall direction).
No axis-aligned BBox or scalar coordinate fields.
"""
import math
from typing import NamedTuple

from shared.types import Point
from shared.geometry import directed_poly_isects, line_isect
from floorplan.constants import (
    WALL_6IN, WALL_4IN, WALL_3IN,
    APPLIANCE_WIDTH, APPLIANCE_DEPTH, APPLIANCE_OFFSET_E,
    APPLIANCE_OFFSET_N, APPLIANCE_GAP,
    COUNTER_DEPTH, COUNTER_GAP,
    BED_WIDTH, BED_LENGTH,
    O9_HALF_WIDTH, O10_HALF_WIDTH, O11_HALF_WIDTH,
    O9_OFFSET_IW11, O9_O10_WALL, O10_O11_WALL, BED_GAP_O9,
    IW1_DIST_FROM_NORTH, IW1_WEST_OFFSET_E, IW2_OFFSET_E,
    IW3_LENGTH, IW3_OFFSET_IW9,
    IW9_LENGTH, IW9_OFFSET_O10,
    IW4_OFFSET_E_IW2,
    IW5_OFFSET_N, IW6_THICKNESS, IW6_OFFSET_N,
    IW4_RO_WIDTH, IW8_OFFSET_N_IW1,
)


class InteriorLayout(NamedTuple):
    """Interior layout positions — all items as [SW, SE, NE, NW] polygons."""
    # Interior walls
    iw1_poly:  list[Point]
    iw2_poly:  list[Point]
    iw3_poly:  list[Point]
    iw4_poly:  list[Point]
    iw5_poly:  list[Point]
    iw6_poly:  list[Point]
    iw7_poly:  list[Point]
    iw8_poly:  list[Point]
    iw9_poly:  list[Point]
    iw11_poly: list[Point]
    iw12_poly: list[Point]
    iw14_poly: list[Point]
    iw15_poly: list[Point]
    iw16_poly: list[Point]
    # Wall thicknesses (scalar dimensions, not positions)
    iwt:  float   # 6" (IW1, IW2, IW8)
    iwt3: float   # 3" (IW5, IW14)
    iwt4: float   # 4" (IW3, IW4, IW7, IW9, IW11, IW12, IW15, IW16)
    # Furniture / appliances
    bed_poly:     list[Point]   # [SW, SE, NE, NW]
    dresser_poly: list[Point]
    dryer_poly:   list[Point]
    washer_poly:  list[Point]
    ctr_poly:     list[Point]   # counter (irregular polygon)
    ctr_nw_r: float
    # South wall opening anchors (parametric t along F20→F1)
    sw_t_o9_start:  float
    sw_t_o9_end:    float
    sw_t_o10_start: float
    sw_t_o10_end:   float
    sw_t_o11_start: float
    sw_t_o11_end:   float


def compute_interior_layout(pts, inner_poly) -> InteriorLayout:
    """Compute interior layout positions.

    pts must contain W-series (W1-W20) and F-series (F1-F20).
    All geometry uses building axis vectors from W20→W1 (south wall).
    """
    # ── Building axis vectors ──────────────────────────────────
    _w20, _w1 = pts["W20"], pts["W1"]
    _dE = _w1[0] - _w20[0]
    _dN = _w1[1] - _w20[1]
    _seg_len = math.sqrt(_dE**2 + _dN**2)
    _along_E = _dE / _seg_len    # unit along south wall (W20→W1 ≈ westward)
    _along_N = _dN / _seg_len
    _norm_E = _along_N            # right normal = into building (≈ northward)
    _norm_N = -_along_E
    _eE = -_along_E               # building east unit vector
    _eN = -_along_N

    # Helpers: building northing/easting relative to W1 on the south wall line
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

    # ── IW1: E-W wall, 6" thick, separating utility/bedroom zones ──
    _w9_bn = _bn(pts["W9"])
    _iw1_n_bn = _w9_bn - IW1_DIST_FROM_NORTH
    _iw1_s_bn = _iw1_n_bn - WALL_6IN
    _w2_be = _be(pts["W2"])
    _iw1_w_be = _w2_be + IW1_WEST_OFFSET_E
    iw1_nw = _bp(_iw1_w_be, _iw1_n_bn)
    iw1_sw = _bp(_iw1_w_be, _iw1_s_bn)
    # East ends: ray east from corners, intersect inner_poly
    _ne_ts = directed_poly_isects(inner_poly, iw1_nw, (_eE, _eN))
    _se_ts = directed_poly_isects(inner_poly, iw1_sw, (_eE, _eN))
    _ne_t = max(t for t in _ne_ts if t > 0)
    _se_t = max(t for t in _se_ts if t > 0)
    iw1_ne = (iw1_nw[0] + _ne_t * _eE, iw1_nw[1] + _ne_t * _eN)
    iw1_se = (iw1_sw[0] + _se_t * _eE, iw1_sw[1] + _se_t * _eN)
    iw1_poly = [iw1_sw, iw1_se, iw1_ne, iw1_nw]

    # ── IW2: N-S wall, 6" thick, west of utility area ─────────
    _iw2_w_be = _w2_be + IW2_OFFSET_E
    _iw2_e_be = _iw2_w_be + WALL_6IN
    _w6_bn = _bn(pts["W6"])
    iw2_poly = [
        _bp(_iw2_w_be, _iw1_n_bn),   # SW
        _bp(_iw2_e_be, _iw1_n_bn),   # SE
        _bp(_iw2_e_be, _w6_bn),      # NE
        _bp(_iw2_w_be, _w6_bn),      # NW
    ]

    # ── IW4: N-S wall, 4" thick, east bedroom wall (partial) ──
    _iw4_w_be = _iw2_e_be + IW4_OFFSET_E_IW2
    _iw4_e_be = _iw4_w_be + WALL_4IN
    _w19_bn = _bn(pts["W19"])
    _iw4_sw = _bp(_iw4_w_be, _w19_bn)

    # ── IW11: 4" thick, perpendicular to W20-W0 ───────────────
    # SE corner: circle(IW4_SW, 32") ∩ W20-W0
    _uE = _w20[0] - _iw4_sw[0]
    _uN = _w20[1] - _iw4_sw[1]
    _r = 32.0 / 12.0
    _qa = _dE**2 + _dN**2
    _qb = 2 * (_uE * _dE + _uN * _dN)
    _qc = _uE**2 + _uN**2 - _r**2
    _disc = _qb**2 - 4 * _qa * _qc
    _t = (-_qb + math.sqrt(_disc)) / (2 * _qa)
    iw11_se = (_w20[0] + _t * _dE, _w20[1] + _t * _dN)
    _iw14_d = 6.0 + WALL_4IN + 3.0 / 12.0 + IW4_RO_WIDTH + 3.0 / 12.0
    _iw11_len = _iw14_d + WALL_3IN
    iw11_sw = (iw11_se[0] + WALL_4IN * _along_E,
               iw11_se[1] + WALL_4IN * _along_N)
    iw11_ne = (iw11_se[0] + _iw11_len * _norm_E,
               iw11_se[1] + _iw11_len * _norm_N)
    iw11_nw = (iw11_sw[0] + _iw11_len * _norm_E,
               iw11_sw[1] + _iw11_len * _norm_N)
    iw11_poly = [iw11_sw, iw11_se, iw11_ne, iw11_nw]

    # ── IW12: 4" thick, perpendicular to IW11 ─────────────────
    _iw12_shorten = 4.0 / 12.0
    _iw12_base = (iw11_sw[0] + 6.0 * _norm_E, iw11_sw[1] + 6.0 * _norm_N)
    iw12_sw = (_iw12_base[0] - _iw12_shorten * _along_E,
               _iw12_base[1] - _iw12_shorten * _along_N)
    # SE: south face line (-_along) meets IW4 west face line (_norm through _iw4_sw)
    iw12_se = line_isect(iw12_sw, (-_along_E, -_along_N),
                         _iw4_sw, (_norm_E, _norm_N))
    iw12_nw = (iw12_sw[0] + WALL_4IN * _norm_E,
               iw12_sw[1] + WALL_4IN * _norm_N)
    iw12_ne = (iw12_se[0] + WALL_4IN * _norm_E,
               iw12_se[1] + WALL_4IN * _norm_N)
    iw12_poly = [iw12_sw, iw12_se, iw12_ne, iw12_nw]

    # ── IW4 complete: north end at IW12 NE building-northing ───
    _iw12_ne_bn = _bn(iw12_ne)
    iw4_poly = [
        _bp(_iw4_w_be, _w19_bn),      # SW
        _bp(_iw4_e_be, _w19_bn),      # SE
        _bp(_iw4_e_be, _iw12_ne_bn),  # NE
        _bp(_iw4_w_be, _iw12_ne_bn),  # NW
    ]

    # ── Bed: rotated, long sides perpendicular to W20-W0 ──────
    _dE9 = pts["F1"][0] - pts["F20"][0]
    _dN9 = pts["F1"][1] - pts["F20"][1]
    _seg9_len = math.sqrt(_dE9**2 + _dN9**2)
    _t_sw9 = ((iw11_sw[0] - pts["F20"][0]) * _dE9
              + (iw11_sw[1] - pts["F20"][1]) * _dN9) / (_dE9**2 + _dN9**2)
    _ts9 = _t_sw9 + O9_OFFSET_IW11 / _seg9_len
    _te9 = _ts9 + 2 * O9_HALF_WIDTH / _seg9_len
    _bed_t = _te9 + BED_GAP_O9 / _seg_len
    _bed_se_wall = (_w20[0] + _bed_t * (_w1[0] - _w20[0]),
                    _w20[1] + _bed_t * (_w1[1] - _w20[1]))
    bed_se = (_bed_se_wall[0] + 2.0 / 12.0 * _norm_E,
              _bed_se_wall[1] + 2.0 / 12.0 * _norm_N)
    bed_sw = (bed_se[0] + BED_WIDTH * _along_E,
              bed_se[1] + BED_WIDTH * _along_N)
    bed_ne = (bed_se[0] + BED_LENGTH * _norm_E,
              bed_se[1] + BED_LENGTH * _norm_N)
    bed_nw = (bed_sw[0] + BED_LENGTH * _norm_E,
              bed_sw[1] + BED_LENGTH * _norm_N)
    bed_poly = [bed_sw, bed_se, bed_ne, bed_nw]

    # ── Opening t-values along F20→F1 ─────────────────────────
    _ts10 = _te9 + O9_O10_WALL / _seg9_len
    _te10 = _ts10 + 2 * O10_HALF_WIDTH / _seg9_len
    _ts11 = _te10 + O10_O11_WALL / _seg9_len
    _te11 = _ts11 + 2 * O11_HALF_WIDTH / _seg9_len

    # ── IW9: 4" thick, perpendicular to W20-W0 ────────────────
    _o10_end = (_w20[0] + _te10 * _dE, _w20[1] + _te10 * _dN)
    iw9_se = (_o10_end[0] + IW9_OFFSET_O10 * _along_E,
              _o10_end[1] + IW9_OFFSET_O10 * _along_N)
    iw9_sw = (iw9_se[0] + WALL_4IN * _along_E,
              iw9_se[1] + WALL_4IN * _along_N)
    iw9_ne = (iw9_se[0] + IW9_LENGTH * _norm_E,
              iw9_se[1] + IW9_LENGTH * _norm_N)
    iw9_nw = (iw9_sw[0] + IW9_LENGTH * _norm_E,
              iw9_sw[1] + IW9_LENGTH * _norm_N)
    iw9_poly = [iw9_sw, iw9_se, iw9_ne, iw9_nw]

    # ── IW3: 4" thick, perpendicular to W20-W0 ────────────────
    iw3_se = (iw9_sw[0] + IW3_OFFSET_IW9 * _along_E,
              iw9_sw[1] + IW3_OFFSET_IW9 * _along_N)
    iw3_sw = (iw3_se[0] + WALL_4IN * _along_E,
              iw3_se[1] + WALL_4IN * _along_N)
    iw3_ne = (iw3_se[0] + IW3_LENGTH * _norm_E,
              iw3_se[1] + IW3_LENGTH * _norm_N)
    iw3_nw = (iw3_sw[0] + IW3_LENGTH * _norm_E,
              iw3_sw[1] + IW3_LENGTH * _norm_N)
    iw3_poly = [iw3_sw, iw3_se, iw3_ne, iw3_nw]

    # ── IW7: 4" thick, parallel to W20-W0 ─────────────────────
    iw7_nw = iw3_ne
    iw7_ne = (iw7_nw[0] - IW3_OFFSET_IW9 * _along_E,
              iw7_nw[1] - IW3_OFFSET_IW9 * _along_N)
    iw7_sw = (iw7_nw[0] - WALL_4IN * _norm_E,
              iw7_nw[1] - WALL_4IN * _norm_N)
    iw7_se = (iw7_ne[0] - WALL_4IN * _norm_E,
              iw7_ne[1] - WALL_4IN * _norm_N)
    iw7_poly = [iw7_sw, iw7_se, iw7_ne, iw7_nw]

    # ── IW16: N-S wall, 4" thick, west face at IW3 NW ─────────
    _iw3_nw_be = _be(iw3_nw)
    _iw3_nw_bn = _bn(iw3_nw)
    iw16_poly = [
        iw3_nw,                                            # SW (= IW3 NW)
        (iw3_nw[0] - WALL_4IN * _along_E,                 # SE
         iw3_nw[1] - WALL_4IN * _along_N),
        _bp(_iw3_nw_be + WALL_4IN, _iw1_s_bn),            # NE
        _bp(_iw3_nw_be, _iw1_s_bn),                       # NW
    ]

    # ── Counter polygon ────────────────────────────────────────
    _dryer_e_be = _w2_be + APPLIANCE_OFFSET_E + APPLIANCE_WIDTH
    _ctr_w_be = _dryer_e_be + COUNTER_GAP
    _ctr_n_bn = 6.0  # 6' north of south wall
    ctr_nw = _bp(_ctr_w_be, _ctr_n_bn)
    ctr_sw = _bp(_ctr_w_be, 0)
    if _ctr_n_bn > _iw3_nw_bn:
        # Counter top above IW3 NW: clip at IW16 west then IW3 west
        ctr_poly = [
            ctr_nw,
            _bp(_iw3_nw_be, _ctr_n_bn),     # at IW16 west face, counter top
            iw3_nw,
            iw3_sw,
            ctr_sw,
        ]
    else:
        # Counter top below IW3 NW: clip at IW3 west face
        _iw3_sw_bn = _bn(iw3_sw)
        _t_face = (_ctr_n_bn - _iw3_sw_bn) / (_iw3_nw_bn - _iw3_sw_bn)
        _ctr_ne = (iw3_sw[0] + _t_face * (iw3_nw[0] - iw3_sw[0]),
                   iw3_sw[1] + _t_face * (iw3_nw[1] - iw3_sw[1]))
        ctr_poly = [
            ctr_nw,
            _ctr_ne,
            iw3_sw,
            ctr_sw,
        ]

    # ── IW14: 3" thick, parallel to IW12 ──────────────────────
    iw14_sw = (iw11_se[0] + _iw14_d * _norm_E,
               iw11_se[1] + _iw14_d * _norm_N)
    iw14_nw = (iw14_sw[0] + WALL_3IN * _norm_E,
               iw14_sw[1] + WALL_3IN * _norm_N)
    # SE: south face line (-_along) meets IW4 west face line
    iw14_se = line_isect(iw14_sw, (-_along_E, -_along_N),
                         _iw4_sw, (_norm_E, _norm_N))
    iw14_ne = (iw14_se[0] + WALL_3IN * _norm_E,
               iw14_se[1] + WALL_3IN * _norm_N)
    iw14_poly = [iw14_sw, iw14_se, iw14_ne, iw14_nw]

    # ── IW5: E-W wall, 3" thick ───────────────────────────────
    _iw5_n_bn = _iw1_s_bn - IW5_OFFSET_N
    _iw5_s_bn = _iw5_n_bn - WALL_3IN
    _w15_be = _be(pts["W15"])
    _iw14_se_be = _be(iw14_se)
    iw5_poly = [
        _bp(_iw14_se_be, _iw5_s_bn),  # SW
        _bp(_w15_be, _iw5_s_bn),      # SE
        _bp(_w15_be, _iw5_n_bn),      # NE
        _bp(_iw14_se_be, _iw5_n_bn),  # NW
    ]

    # ── IW15: N-S wall, 4" thick ──────────────────────────────
    _iw11_nw_be = _be(iw11_nw)
    _iw11_nw_bn = _bn(iw11_nw)
    iw15_poly = [
        _bp(_iw11_nw_be, _iw11_nw_bn),              # SW
        _bp(_iw11_nw_be + WALL_4IN, _iw11_nw_bn),   # SE
        _bp(_iw11_nw_be + WALL_4IN, _iw1_s_bn),     # NE
        _bp(_iw11_nw_be, _iw1_s_bn),                # NW
    ]

    # ── IW8: E-W wall, 6" thick ───────────────────────────────
    _iw8_n_bn = _iw1_n_bn + IW8_OFFSET_N_IW1
    _iw8_s_bn = _iw1_s_bn + IW8_OFFSET_N_IW1
    iw8_poly = [
        _bp(_w2_be, _iw8_s_bn),      # SW
        _bp(_iw1_w_be, _iw8_s_bn),   # SE
        _bp(_iw1_w_be, _iw8_n_bn),   # NE
        _bp(_w2_be, _iw8_n_bn),      # NW
    ]

    # ── IW6: E-W wall, thin ───────────────────────────────────
    _iw6_n_bn = _w6_bn - IW6_OFFSET_N
    _iw6_s_bn = _iw6_n_bn - IW6_THICKNESS
    _iw6_ne = _bp(_iw2_w_be, _iw6_n_bn)
    _iw6_se = _bp(_iw2_w_be, _iw6_s_bn)
    # West ends: ray westward, intersect inner_poly
    _iw6_nw_ts = directed_poly_isects(inner_poly, _iw6_ne, (_along_E, _along_N))
    _iw6_sw_ts = directed_poly_isects(inner_poly, _iw6_se, (_along_E, _along_N))
    _iw6_nw_t = min(t for t in _iw6_nw_ts if t > 0)
    _iw6_sw_t = min(t for t in _iw6_sw_ts if t > 0)
    iw6_nw = (_iw6_ne[0] + _iw6_nw_t * _along_E,
              _iw6_ne[1] + _iw6_nw_t * _along_N)
    iw6_sw = (_iw6_se[0] + _iw6_sw_t * _along_E,
              _iw6_se[1] + _iw6_sw_t * _along_N)
    iw6_poly = [iw6_sw, _iw6_se, _iw6_ne, iw6_nw]

    # ── Dryer ──────────────────────────────────────────────────
    _dryer_w_be = _w2_be + APPLIANCE_OFFSET_E
    _dryer_s_bn = APPLIANCE_OFFSET_N  # offset from south wall (W1 bn = 0)
    _dryer_n_bn = _dryer_s_bn + APPLIANCE_DEPTH
    dryer_poly = [
        _bp(_dryer_w_be, _dryer_s_bn),   # SW
        _bp(_dryer_e_be, _dryer_s_bn),   # SE
        _bp(_dryer_e_be, _dryer_n_bn),   # NE
        _bp(_dryer_w_be, _dryer_n_bn),   # NW
    ]

    # ── Washer ─────────────────────────────────────────────────
    _washer_s_bn = _dryer_n_bn + APPLIANCE_GAP
    _washer_n_bn = _washer_s_bn + APPLIANCE_DEPTH
    washer_poly = [
        _bp(_dryer_w_be, _washer_s_bn),  # SW
        _bp(_dryer_e_be, _washer_s_bn),  # SE
        _bp(_dryer_e_be, _washer_n_bn),  # NE
        _bp(_dryer_w_be, _washer_n_bn),  # NW
    ]

    # ── Dresser ────────────────────────────────────────────────
    _dresser_e_be = _iw11_nw_be - 2.0 / 12.0  # 2" west of IW15 west face
    _dresser_w_be = _dresser_e_be - 34.0 / 12.0
    _dresser_n_bn = _iw1_s_bn - 1.0 / 12.0
    _dresser_s_bn = _dresser_n_bn - 19.0 / 12.0
    dresser_poly = [
        _bp(_dresser_w_be, _dresser_s_bn),  # SW
        _bp(_dresser_e_be, _dresser_s_bn),  # SE
        _bp(_dresser_e_be, _dresser_n_bn),  # NE
        _bp(_dresser_w_be, _dresser_n_bn),  # NW
    ]

    return InteriorLayout(
        iw1_poly=iw1_poly, iw2_poly=iw2_poly, iw3_poly=iw3_poly,
        iw4_poly=iw4_poly, iw5_poly=iw5_poly, iw6_poly=iw6_poly,
        iw7_poly=iw7_poly, iw8_poly=iw8_poly, iw9_poly=iw9_poly,
        iw11_poly=iw11_poly, iw12_poly=iw12_poly, iw14_poly=iw14_poly,
        iw15_poly=iw15_poly, iw16_poly=iw16_poly,
        iwt=WALL_6IN, iwt3=WALL_3IN, iwt4=WALL_4IN,
        bed_poly=bed_poly, dresser_poly=dresser_poly,
        dryer_poly=dryer_poly, washer_poly=washer_poly,
        ctr_poly=ctr_poly, ctr_nw_r=0,
        sw_t_o9_start=_ts9, sw_t_o9_end=_te9,
        sw_t_o10_start=_ts10, sw_t_o10_end=_te10,
        sw_t_o11_start=_ts11, sw_t_o11_end=_te11,
    )
