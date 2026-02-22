"""Interior layout computation — rooms, walls, appliances, furniture."""
import math
from typing import NamedTuple

from shared.types import Point, BBox
from shared.geometry import horiz_isects
from floorplan.constants import (
    WALL_6IN, WALL_4IN, WALL_3IN,
    APPLIANCE_WIDTH, APPLIANCE_DEPTH, APPLIANCE_OFFSET_FROM_W2,
    APPLIANCE_OFFSET_FROM_W1, APPLIANCE_GAP,
    COUNTER_DEPTH, COUNTER_GAP,
    BED_WIDTH, BED_LENGTH,
    O9_HALF_WIDTH, O10_HALF_WIDTH, O11_HALF_WIDTH,
    O9_OFFSET_IW11, O9_O10_WALL, O10_O11_WALL, BED_GAP_O9,
    IW1_OFFSET_FROM_W9, IW1_OFFSET_FROM_W2, IW2_OFFSET_FROM_W2,
    IW3_LENGTH, IW3_OFFSET_IW9,
    IW9_LENGTH, IW9_OFFSET_O10,
    IW4_OFFSET_FROM_IW2,
    IW5_OFFSET_FROM_IW1, IW6_THICKNESS, IW6_OFFSET_FROM_W6,
    IW4_RO_WIDTH, IW8_OFFSET_FROM_IW1,
)


def _seg_vecs(p1, p2):
    """Along-direction and CW-inward normal for segment p1->p2."""
    dx, dy = p2[0] - p1[0], p2[1] - p1[1]
    length = math.sqrt(dx * dx + dy * dy)
    along = (dx / length, dy / length)
    inward = (dy / length, -dx / length)  # right perp = CW inward
    return along, inward


def _offset(origin, dist, direction):
    """Offset point by dist along direction vector."""
    return (origin[0] + dist * direction[0], origin[1] + dist * direction[1])


class InteriorLayout(NamedTuple):
    """Interior layout positions for walls, appliances, and furniture."""
    # Interior wall 1 (IW1) — horizontal wall separating utility/bedroom zones
    iw1: list[Point]
    iw1_s: float
    iw1_n: float
    iwt: float
    # Interior wall 2 (IW2) — vertical wall west of utility area
    iw2: BBox
    # Dryer
    dryer: BBox
    # Washer
    washer: BBox
    # Counter
    ctr: BBox
    ctr_poly: list[Point]  # polygon clipped to W20-W0 and IW3/IW16 west faces
    ctr_nw_r: float
    # Interior wall 3 (IW3) — perpendicular to W20-W0, 30" from IW9 W face
    iw3: BBox
    iw3_poly: list[Point]  # [SW, SE, NE, NW]
    # Interior wall 7 (IW7) — parallel to W20-W0, between IW3 and IW9
    iw7: BBox
    iw7_poly: list[Point]  # [SW, SE, NE, NW]
    # Interior wall 9 (IW9) — perpendicular to W20-W0, 8" past O10
    iw9: BBox
    iw9_poly: list[Point]  # [SW, SE, NE, NW]
    # Wall thicknesses
    iwt3: float
    iwt4: float
    # Interior wall 4 (IW4) — east bedroom wall
    iw4_w: float
    iw4_e: float
    iw4_s: float
    wall_south_n: float
    # Bed (rotated polygon, long sides perpendicular to W20-W0)
    bed_poly: list[Point]  # [SW, SE, NE, NW]
    # IW11 (4" thick, N-S — south extension below IW4)
    iw11: BBox         # bounding box (for labels/dimensions)
    iw11_poly: list[Point]  # actual polygon (SE corner on W20-W0)
    # IW12 (4" thick, perpendicular to IW11 — connects IW11 NW to IW4 west)
    iw12: BBox
    iw12_poly: list[Point]  # actual polygon [SW, SE, NE, NW]
    # IW5 (3" thick, horizontal in office)
    iw5: BBox
    # IW8 (6" thick, horizontal — west extension of IW1, W2-W3 to IW1 west end)
    iw8: BBox
    # IW14 (3" thick, parallel to IW12, north of RO2)
    iw14: BBox
    iw14_poly: list[Point]  # [SW, SE, NE, NW]
    # IW15 (4" thick, N-S, west face at IW11 NW easting, IW11 north to IW1 south)
    iw15: BBox
    # Dresser (34" E-W × 19" N-S, 1" south of IW1, 6" west of IW15)
    dresser: BBox
    # IW16 (4" thick, N-S, west face at IW3 NW easting, IW1 to IW3 end)
    iw16_poly: list[Point]  # [SW, SE, NE, NW]
    # IW6 (1" thick, horizontal above kitchen)
    iw6_poly: list[Point]
    iw6_n: float
    iw6_s: float
    # South wall opening anchors (parametric t along F20→F1)
    sw_t_o9_start: float
    sw_t_o9_end: float
    sw_t_o10_start: float
    sw_t_o10_end: float
    sw_t_o11_start: float
    sw_t_o11_end: float


def compute_interior_layout(pts, inner_poly) -> InteriorLayout:
    """Compute interior layout positions.

    pts must contain W-series (W1-W20) and F-series (F1-F20).
    """
    # Segment vectors for axis-aligned wall segments
    _w2w3_al, _w2w3_in = _seg_vecs(pts["W2"], pts["W3"])   # along=(0,1), inward=(1,0)
    _w6w7_al, _w6w7_in = _seg_vecs(pts["W6"], pts["W7"])   # along=(1,0), inward=(0,-1)
    _w9w10_al, _w9w10_in = _seg_vecs(pts["W9"], pts["W10"])  # along=(1,0), inward=(0,-1)

    iw1_n = _offset(pts["W9"], IW1_OFFSET_FROM_W9, _w9w10_in)[1]
    iw1_s = iw1_n - WALL_6IN
    si = horiz_isects(inner_poly, iw1_s)
    ni = horiz_isects(inner_poly, iw1_n)
    iw1_w = _offset(pts["W2"], IW1_OFFSET_FROM_W2, _w2w3_in)[0]
    iw1 = [(iw1_w, iw1_s), (max(si), iw1_s), (max(ni), iw1_n), (iw1_w, iw1_n)]

    iw2_w = _offset(pts["W2"], IW2_OFFSET_FROM_W2, _w2w3_in)[0]
    iw2_e = iw2_w + WALL_6IN
    iw2_s = iw1_n
    iw2_n = pts["W6"][1]

    dryer_w = _offset(pts["W2"], APPLIANCE_OFFSET_FROM_W2, _w2w3_in)[0]
    dryer_s = _offset(pts["W1"], APPLIANCE_OFFSET_FROM_W1, _w2w3_al)[1]
    dryer_e = dryer_w + APPLIANCE_WIDTH
    dryer_n = dryer_s + APPLIANCE_DEPTH
    washer_w = dryer_w
    washer_s = dryer_n + APPLIANCE_GAP
    washer_e = dryer_e
    washer_n = washer_s + APPLIANCE_DEPTH

    ctr_w = dryer_e + COUNTER_GAP
    ctr_e = ctr_w + COUNTER_DEPTH
    ctr_s = pts["W1"][1]
    ctr_n = ctr_s + 6.0  # 6' along W2-W3 from W20-W0 south face
    ctr_nw_r = 0

    iw2_e = iw2_w + WALL_6IN
    iw4_w = iw2_e + IW4_OFFSET_FROM_IW2
    iw4_e = iw4_w + WALL_4IN
    iw4_s = pts["W19"][1]
    wall_south_n = pts["W19"][1]

    # IW11: 4" thick, normal to W20-W0, 6' long
    # SE corner: circle(IW4_SW, 32") ∩ W20-W0
    _iw4_sw = (iw4_w, iw4_s)
    _w20 = pts["W20"]
    _w1 = pts["W1"]
    _dE = _w1[0] - _w20[0]
    _dN = _w1[1] - _w20[1]
    _seg_len = math.sqrt(_dE**2 + _dN**2)
    _uE = _w20[0] - _iw4_sw[0]
    _uN = _w20[1] - _iw4_sw[1]
    _r = 32.0 / 12.0
    _qa = _dE**2 + _dN**2
    _qb = 2 * (_uE * _dE + _uN * _dN)
    _qc = _uE**2 + _uN**2 - _r**2
    _disc = _qb**2 - 4 * _qa * _qc
    _t = (-_qb + math.sqrt(_disc)) / (2 * _qa)  # westward along W20-W0
    iw11_se = (_w20[0] + _t * _dE, _w20[1] + _t * _dN)
    # Unit vectors: along W20-W0 and inward normal
    _along_E = _dE / _seg_len
    _along_N = _dN / _seg_len
    _norm_E = _along_N    # right normal = inward
    _norm_N = -_along_E
    _iw11_thick = WALL_4IN
    _iw14_d = 6.0 + WALL_4IN + 3.0 / 12.0 + IW4_RO_WIDTH + 3.0 / 12.0
    _iw11_len = _iw14_d + WALL_3IN  # extend to IW14 north face
    iw11_sw = (iw11_se[0] + _iw11_thick * _along_E,
               iw11_se[1] + _iw11_thick * _along_N)
    iw11_ne = (iw11_se[0] + _iw11_len * _norm_E,
               iw11_se[1] + _iw11_len * _norm_N)
    iw11_nw = (iw11_sw[0] + _iw11_len * _norm_E,
               iw11_sw[1] + _iw11_len * _norm_N)
    iw11_poly = [iw11_sw, iw11_se, iw11_ne, iw11_nw]
    # Bounding box (for IW12 connection and labels)
    iw11_w = min(p[0] for p in iw11_poly)
    iw11_e = max(p[0] for p in iw11_poly)
    iw11_s = min(p[1] for p in iw11_poly)
    iw11_n = max(p[1] for p in iw11_poly)

    # IW12: 4" thick, perpendicular to IW11, from IW11 NW corner to IW4 west face
    _iw12_shorten = 4.0 / 12.0
    _iw12_base = (iw11_sw[0] + 6.0 * _norm_E, iw11_sw[1] + 6.0 * _norm_N)
    iw12_sw = (_iw12_base[0] - _iw12_shorten * _along_E,
               _iw12_base[1] - _iw12_shorten * _along_N)
    # IW12 east end: SE and NE at iw4_w easting
    # South edge: line from iw12_sw in -_along direction, solve for easting = iw4_w
    _t_se = (iw4_w - iw12_sw[0]) / (-_along_E)
    iw12_se = (iw4_w, iw12_sw[1] - _t_se * _along_N)
    iw12_nw = (iw12_sw[0] + WALL_4IN * _norm_E,
               iw12_sw[1] + WALL_4IN * _norm_N)
    iw12_ne = (iw4_w, iw12_nw[1] + (iw12_se[1] - iw12_sw[1]))
    iw12_poly = [iw12_sw, iw12_se, iw12_ne, iw12_nw]
    # Bounding box (for RO2 center, dimension lines, labels)
    iw12_w = min(p[0] for p in iw12_poly)
    iw12_e = max(p[0] for p in iw12_poly)
    iw12_s = min(p[1] for p in iw12_poly)
    iw12_n = max(p[1] for p in iw12_poly)

    # Bed: rotated, long sides perpendicular to W20-W0
    # SE corner = 1" past O9 NW corner along W20-W0, 2" from wall
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

    # IW9: 4" thick, perpendicular to W20-W0, past O10 along inner wall
    _ts10 = _te9 + O9_O10_WALL / _seg9_len
    _te10 = _ts10 + 2 * O10_HALF_WIDTH / _seg9_len
    _ts11 = _te10 + O10_O11_WALL / _seg9_len
    _te11 = _ts11 + 2 * O11_HALF_WIDTH / _seg9_len
    _o10_end = (_w20[0] + _te10 * _dE, _w20[1] + _te10 * _dN)
    iw9_base = (_o10_end[0] + IW9_OFFSET_O10 * _along_E,
                _o10_end[1] + IW9_OFFSET_O10 * _along_N)
    iw9_se = iw9_base
    iw9_sw = (iw9_se[0] + WALL_4IN * _along_E,
              iw9_se[1] + WALL_4IN * _along_N)
    iw9_ne = (iw9_se[0] + IW9_LENGTH * _norm_E,
              iw9_se[1] + IW9_LENGTH * _norm_N)
    iw9_nw = (iw9_sw[0] + IW9_LENGTH * _norm_E,
              iw9_sw[1] + IW9_LENGTH * _norm_N)
    iw9_poly = [iw9_sw, iw9_se, iw9_ne, iw9_nw]
    iw9_w = min(p[0] for p in iw9_poly)
    iw9_e = max(p[0] for p in iw9_poly)
    iw9_s = min(p[1] for p in iw9_poly)
    iw9_n = max(p[1] for p in iw9_poly)

    # IW3: 4" thick, perpendicular to W20-W0, E face 30" from IW9 W face
    iw3_se = (iw9_sw[0] + IW3_OFFSET_IW9 * _along_E,
              iw9_sw[1] + IW3_OFFSET_IW9 * _along_N)
    iw3_sw = (iw3_se[0] + WALL_4IN * _along_E,
              iw3_se[1] + WALL_4IN * _along_N)
    iw3_ne = (iw3_se[0] + IW3_LENGTH * _norm_E,
              iw3_se[1] + IW3_LENGTH * _norm_N)
    iw3_nw = (iw3_sw[0] + IW3_LENGTH * _norm_E,
              iw3_sw[1] + IW3_LENGTH * _norm_N)
    iw3_poly = [iw3_sw, iw3_se, iw3_ne, iw3_nw]
    iw3_w = min(p[0] for p in iw3_poly)
    iw3_e = max(p[0] for p in iw3_poly)
    iw3_s = min(p[1] for p in iw3_poly)
    iw3_n = max(p[1] for p in iw3_poly)

    # IW7: 4" thick, parallel to W20-W0, NW corner at IW3 NE, spans to IW9
    iw7_nw = iw3_ne
    iw7_ne = (iw7_nw[0] - IW3_OFFSET_IW9 * _along_E,
              iw7_nw[1] - IW3_OFFSET_IW9 * _along_N)
    iw7_sw = (iw7_nw[0] - WALL_4IN * _norm_E,
              iw7_nw[1] - WALL_4IN * _norm_N)
    iw7_se = (iw7_ne[0] - WALL_4IN * _norm_E,
              iw7_ne[1] - WALL_4IN * _norm_N)
    iw7_poly = [iw7_sw, iw7_se, iw7_ne, iw7_nw]
    iw7_w = min(p[0] for p in iw7_poly)
    iw7_e = max(p[0] for p in iw7_poly)
    iw7_s = min(p[1] for p in iw7_poly)
    iw7_n = max(p[1] for p in iw7_poly)

    # IW16: 4" thick, N-S, west face at IW3 NW easting, IW1 south to IW3 NW northing
    iw16_w = iw3_nw[0]
    iw16_e = iw16_w + WALL_4IN
    iw16_n = iw1_s
    iw16_s = iw3_nw[1]
    iw16_poly = [(iw16_w, iw16_s), (iw16_e, iw16_s),
                 (iw16_e, iw16_n), (iw16_w, iw16_n)]

    # Counter polygon: south edge follows W20-W0, east edge clipped at IW3/IW16
    _t_ctr_w = (ctr_w - _w20[0]) / (_w1[0] - _w20[0])
    _ctr_sw_n = _w20[1] + _t_ctr_w * (_w1[1] - _w20[1])
    if ctr_n > iw3_nw[1]:
        # Counter top above IW3 NW: clip at IW16 west face then IW3 west face
        ctr_poly = [
            (ctr_w, ctr_n),
            (iw16_w, ctr_n),
            iw3_nw,
            iw3_sw,
            (ctr_w, _ctr_sw_n),
        ]
    else:
        # Counter top below IW3 NW: clip at IW3 west face only
        _t_face = (ctr_n - iw3_sw[1]) / (iw3_nw[1] - iw3_sw[1])
        _ctr_ne_e = iw3_sw[0] + _t_face * (iw3_nw[0] - iw3_sw[0])
        ctr_poly = [
            (ctr_w, ctr_n),
            (_ctr_ne_e, ctr_n),
            iw3_sw,
            (ctr_w, _ctr_sw_n),
        ]

    # IW8: 6" thick, horizontal, from W2-W3 face to IW1 west end
    iw8_w = pts["W2"][0]
    iw8_e = iw1_w
    iw8 = BBox(w=iw8_w, s=iw1_s + IW8_OFFSET_FROM_IW1, e=iw8_e, n=iw1_n + IW8_OFFSET_FROM_IW1)

    # IW5: 3" thick, north face IW5_OFFSET_FROM_IW1 from IW1, CW-normal to W9-W10
    iw5_n = iw1_s - IW5_OFFSET_FROM_IW1
    iw5_s = iw5_n - WALL_3IN
    iw5_e = pts["W15"][0]  # east at W15

    # IW14: 3" thick, parallel to IW12, 3" past RO2 north edge along IW11
    # _iw14_d computed earlier (used for IW11 length)
    iw14_sw = (iw11_se[0] + _iw14_d * _norm_E,
               iw11_se[1] + _iw14_d * _norm_N)
    iw14_nw = (iw14_sw[0] + WALL_3IN * _norm_E,
               iw14_sw[1] + WALL_3IN * _norm_N)
    # SE corner: south face line (-_along direction) hits iw5_s northing
    _t14 = (iw14_sw[1] - iw5_s) / _along_N
    iw14_se = (iw14_sw[0] - _t14 * _along_E, iw5_s)
    iw14_ne = (iw14_se[0] + WALL_3IN * _norm_E,
               iw14_se[1] + WALL_3IN * _norm_N)
    iw14_poly = [iw14_sw, iw14_se, iw14_ne, iw14_nw]
    iw14_w = min(p[0] for p in iw14_poly)
    iw14_e = max(p[0] for p in iw14_poly)
    iw14_s = min(p[1] for p in iw14_poly)
    iw14_n = max(p[1] for p in iw14_poly)

    # IW15: 4" thick, N-S, west face at IW11 NW easting, IW11 north to IW1 south
    iw15_w = iw11_nw[0]
    iw15 = BBox(w=iw15_w, s=iw11_nw[1], e=iw15_w + WALL_4IN, n=iw1_s)

    # Dresser: 34" E-W × 19" N-S, 1" south of IW1, 2" west of IW15
    dresser_e = iw15_w - 2.0 / 12.0
    dresser_w = dresser_e - 34.0 / 12.0
    dresser_n = iw1_s - 1.0 / 12.0
    dresser_s = dresser_n - 19.0 / 12.0
    dresser = BBox(w=dresser_w, s=dresser_s, e=dresser_e, n=dresser_n)

    # IW5 west end: meets IW14 SE corner
    iw5_w = iw14_se[0]

    # IW6: IW6_THICKNESS thick, IW6_OFFSET_FROM_W6 from W6, CW-normal to W6-W7
    iw6_n = _offset(pts["W6"], IW6_OFFSET_FROM_W6, _w6w7_in)[1]
    iw6_s = iw6_n - IW6_THICKNESS
    _iw6_n_ints = horiz_isects(inner_poly, iw6_n)
    _iw6_s_ints = horiz_isects(inner_poly, iw6_s)
    iw6_w_n = min(_iw6_n_ints)
    iw6_w_s = min(_iw6_s_ints)
    iw6_e = iw2_w
    iw6_poly = [(iw6_w_s, iw6_s), (iw6_e, iw6_s), (iw6_e, iw6_n), (iw6_w_n, iw6_n)]

    return InteriorLayout(
        iw1=iw1, iw1_s=iw1_s, iw1_n=iw1_n, iwt=WALL_6IN,
        iw2=BBox(w=iw2_w, s=iw2_s, e=iw2_e, n=iw2_n),
        iw3=BBox(w=iw3_w, s=iw3_s, e=iw3_e, n=iw3_n), iw3_poly=iw3_poly,
        iw7=BBox(w=iw7_w, s=iw7_s, e=iw7_e, n=iw7_n), iw7_poly=iw7_poly,
        iw9=BBox(w=iw9_w, s=iw9_s, e=iw9_e, n=iw9_n), iw9_poly=iw9_poly,
        dryer=BBox(w=dryer_w, s=dryer_s, e=dryer_e, n=dryer_n),
        washer=BBox(w=washer_w, s=washer_s, e=washer_e, n=washer_n),
        ctr=BBox(w=ctr_w, s=ctr_s, e=ctr_e, n=ctr_n), ctr_poly=ctr_poly, ctr_nw_r=ctr_nw_r,
        iwt3=WALL_3IN, iwt4=WALL_4IN,
        iw4_w=iw4_w, iw4_e=iw4_e, iw4_s=iw4_s, wall_south_n=wall_south_n,
        iw11=BBox(w=iw11_w, s=iw11_s, e=iw11_e, n=iw11_n), iw11_poly=iw11_poly,
        iw12=BBox(w=iw12_w, s=iw12_s, e=iw12_e, n=iw12_n), iw12_poly=iw12_poly,
        bed_poly=bed_poly,
        iw8=iw8,
        iw5=BBox(w=iw5_w, s=iw5_s, e=iw5_e, n=iw5_n),
        iw14=BBox(w=iw14_w, s=iw14_s, e=iw14_e, n=iw14_n), iw14_poly=iw14_poly,
        iw15=iw15, dresser=dresser,
        iw16_poly=iw16_poly,
        iw6_poly=iw6_poly, iw6_n=iw6_n, iw6_s=iw6_s,
        sw_t_o9_start=_ts9, sw_t_o9_end=_te9,
        sw_t_o10_start=_ts10, sw_t_o10_end=_te10,
        sw_t_o11_start=_ts11, sw_t_o11_end=_te11,
    )
