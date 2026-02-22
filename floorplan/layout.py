"""Interior layout computation — rooms, walls, appliances, furniture."""
import math
from typing import NamedTuple

from shared.types import Point, BBox
from shared.geometry import line_isect
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


def _dot(a, b):
    """Dot product of two 2D vectors."""
    return a[0] * b[0] + a[1] * b[1]


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
    # Segment vectors for wall segments
    _w2w3_al, _w2w3_in = _seg_vecs(pts["W2"], pts["W3"])   # along=(0,1), inward=(1,0)
    _w6w7_al, _w6w7_in = _seg_vecs(pts["W6"], pts["W7"])   # along=(1,0), inward=(0,-1)
    _w9w10_al, _w9w10_in = _seg_vecs(pts["W9"], pts["W10"])  # along=(1,0), inward=(0,-1)
    _w20w1_al, _w20w1_in = _seg_vecs(pts["W20"], pts["W1"])

    _iw1_n_anchor = _offset(pts["W9"], IW1_OFFSET_FROM_W9, _w9w10_in)
    _iw1_w_anchor = _offset(pts["W2"], IW1_OFFSET_FROM_W2, _w2w3_in)
    iw1_nw = line_isect(_iw1_n_anchor, _w9w10_al, _iw1_w_anchor, _w2w3_al)
    iw1_sw = _offset(iw1_nw, WALL_6IN, _w9w10_in)
    _ne_ts = _line_poly_isects(inner_poly, iw1_nw, _w9w10_al)
    iw1_ne = _offset(iw1_nw, max(_ne_ts), _w9w10_al)
    _se_ts = _line_poly_isects(inner_poly, iw1_sw, _w9w10_al)
    iw1_se = _offset(iw1_sw, max(_se_ts), _w9w10_al)
    iw1 = [iw1_sw, iw1_se, iw1_ne, iw1_nw]
    iw1_s, iw1_n = iw1_sw[1], iw1_nw[1]

    _iw2_w_anchor = _offset(pts["W2"], IW2_OFFSET_FROM_W2, _w2w3_in)
    _iw2_e_anchor = _offset(_iw2_w_anchor, WALL_6IN, _w2w3_in)
    iw2_sw = line_isect(_iw2_w_anchor, _w2w3_al, iw1_nw, _w9w10_al)
    iw2_se = line_isect(_iw2_e_anchor, _w2w3_al, iw1_nw, _w9w10_al)
    iw2_nw = line_isect(_iw2_w_anchor, _w2w3_al, pts["W6"], _w6w7_al)
    iw2_ne = line_isect(_iw2_e_anchor, _w2w3_al, pts["W6"], _w6w7_al)
    iw2_w, iw2_e = iw2_sw[0], iw2_se[0]
    iw2_s, iw2_n = iw2_sw[1], iw2_nw[1]

    _dryer_sw = line_isect(
        _offset(pts["W2"], APPLIANCE_OFFSET_FROM_W2, _w2w3_in), _w2w3_al,
        _offset(pts["W1"], APPLIANCE_OFFSET_FROM_W1, _w2w3_al), _w9w10_al)
    _dryer_se = _offset(_dryer_sw, APPLIANCE_WIDTH, _w2w3_in)
    _dryer_nw = _offset(_dryer_sw, APPLIANCE_DEPTH, _w2w3_al)
    _washer_sw = _offset(_dryer_nw, APPLIANCE_GAP, _w2w3_al)
    _washer_se = _offset(_washer_sw, APPLIANCE_WIDTH, _w2w3_in)
    _washer_nw = _offset(_washer_sw, APPLIANCE_DEPTH, _w2w3_al)
    dryer_w, dryer_s = _dryer_sw[0], _dryer_sw[1]
    dryer_e, dryer_n = _dryer_se[0], _dryer_nw[1]
    washer_w, washer_s = _washer_sw[0], _washer_sw[1]
    washer_e, washer_n = _washer_se[0], _washer_nw[1]

    _ctr_sw_anchor = _offset(_dryer_se, COUNTER_GAP, _w2w3_in)
    _ctr_se_anchor = _offset(_ctr_sw_anchor, COUNTER_DEPTH, _w2w3_in)
    # BBox: south at W1 projection onto counter west face
    _ctr_s_pt = line_isect(_ctr_sw_anchor, _w2w3_al, pts["W1"], _w9w10_al)
    _ctr_n_pt = _offset(_ctr_s_pt, 6.0, _w2w3_al)
    ctr_w, ctr_s = _ctr_sw_anchor[0], _ctr_s_pt[1]
    ctr_e = _ctr_se_anchor[0]
    ctr_n = _ctr_n_pt[1]
    ctr_nw_r = 0

    _iw4_w_anchor = _offset(iw2_se, IW4_OFFSET_FROM_IW2, _w2w3_in)
    _iw4_e_anchor = _offset(_iw4_w_anchor, WALL_4IN, _w2w3_in)
    _iw4_t_s = _dot((pts["W19"][0] - _iw4_w_anchor[0],
                      pts["W19"][1] - _iw4_w_anchor[1]), _w2w3_al)
    iw4_sw = _offset(_iw4_w_anchor, _iw4_t_s, _w2w3_al)
    iw4_se = _offset(_iw4_e_anchor, _iw4_t_s, _w2w3_al)
    iw4_w, iw4_e = iw4_sw[0], iw4_se[0]
    iw4_s = iw4_sw[1]
    wall_south_n = iw4_s

    # IW11: 4" thick, normal to W20-W0, 6' long
    # SE corner: circle(IW4_SW, 32") ∩ W20-W0
    _iw4_sw = iw4_sw
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

    # IW16: 4" thick, N-S, west face at IW3 NW, south at IW3 NW, north at IW1 south
    _iw16_w_anchor = iw3_nw
    _iw16_e_anchor = _offset(_iw16_w_anchor, WALL_4IN, _w2w3_in)
    _iw16_t_n = _dot((iw1_sw[0] - _iw16_w_anchor[0],
                       iw1_sw[1] - _iw16_w_anchor[1]), _w2w3_al)
    iw16_nw = _offset(_iw16_w_anchor, _iw16_t_n, _w2w3_al)
    iw16_ne = _offset(_iw16_e_anchor, _iw16_t_n, _w2w3_al)
    iw16_w = _iw16_w_anchor[0]
    iw16_poly = [_iw16_w_anchor, _iw16_e_anchor, iw16_ne, iw16_nw]

    # Counter polygon: south edge follows W20-W1, east edge clipped at IW3/IW16
    _ctr_poly_sw = line_isect(_ctr_sw_anchor, _w2w3_al, pts["W20"], _w20w1_al)
    _iw3_w_dir = (_norm_E, _norm_N)  # IW3 west face direction (along IW3 length)
    _ctr_vs_iw3 = _dot((_ctr_n_pt[0] - iw3_nw[0],
                         _ctr_n_pt[1] - iw3_nw[1]), _w2w3_al)
    if _ctr_vs_iw3 > 0:
        # Counter top past IW3 NW: clip at IW16 west face then IW3 west face
        _ctr_ne_clip = line_isect(_ctr_n_pt, _w2w3_in, _iw16_w_anchor, _w2w3_al)
        ctr_poly = [_ctr_n_pt, _ctr_ne_clip, iw3_nw, iw3_sw, _ctr_poly_sw]
    else:
        # Counter top below IW3 NW: clip at IW3 west face only
        _ctr_ne_clip = line_isect(_ctr_n_pt, _w2w3_in, iw3_sw, _iw3_w_dir)
        ctr_poly = [_ctr_n_pt, _ctr_ne_clip, iw3_sw, _ctr_poly_sw]

    # IW8: 6" thick, parallel to W9-W10, from W2-W3 face to IW1 west end
    _iw8_n_anchor = _offset(iw1_nw, -IW8_OFFSET_FROM_IW1, _w9w10_in)
    _iw8_s_anchor = _offset(iw1_sw, -IW8_OFFSET_FROM_IW1, _w9w10_in)
    iw8_nw = line_isect(_iw8_n_anchor, _w9w10_al, pts["W2"], _w2w3_al)
    iw8_sw = line_isect(_iw8_s_anchor, _w9w10_al, pts["W2"], _w2w3_al)
    iw8 = BBox(w=iw8_nw[0], s=iw8_sw[1], e=iw1_nw[0], n=iw8_nw[1])

    # IW5: 3" thick, parallel to W9-W10, offset from IW1 south face
    _iw5_n_anchor = _offset(iw1_sw, IW5_OFFSET_FROM_IW1, _w9w10_in)
    _iw5_s_anchor = _offset(_iw5_n_anchor, WALL_3IN, _w9w10_in)
    _iw5_t_e = _dot((pts["W15"][0] - _iw5_n_anchor[0],
                      pts["W15"][1] - _iw5_n_anchor[1]), _w9w10_al)
    iw5_ne = _offset(_iw5_n_anchor, _iw5_t_e, _w9w10_al)
    iw5_n, iw5_s, iw5_e = _iw5_n_anchor[1], _iw5_s_anchor[1], iw5_ne[0]

    # IW14: 3" thick, parallel to IW12, 3" past RO2 north edge along IW11
    # _iw14_d computed earlier (used for IW11 length)
    iw14_sw = (iw11_se[0] + _iw14_d * _norm_E,
               iw11_se[1] + _iw14_d * _norm_N)
    iw14_nw = (iw14_sw[0] + WALL_3IN * _norm_E,
               iw14_sw[1] + WALL_3IN * _norm_N)
    # SE corner: IW14 south face meets IW5 south face
    iw14_se = line_isect(iw14_sw, (-_along_E, -_along_N),
                         _iw5_s_anchor, _w9w10_al)
    iw14_ne = (iw14_se[0] + WALL_3IN * _norm_E,
               iw14_se[1] + WALL_3IN * _norm_N)
    iw14_poly = [iw14_sw, iw14_se, iw14_ne, iw14_nw]
    iw14_w = min(p[0] for p in iw14_poly)
    iw14_e = max(p[0] for p in iw14_poly)
    iw14_s = min(p[1] for p in iw14_poly)
    iw14_n = max(p[1] for p in iw14_poly)

    # IW15: 4" thick, parallel to W2-W3, west face at IW11 NW, north at IW1 south
    _iw15_w_anchor = iw11_nw
    _iw15_e_anchor = _offset(_iw15_w_anchor, WALL_4IN, _w2w3_in)
    _iw15_t_n = _dot((iw1_sw[0] - _iw15_w_anchor[0],
                       iw1_sw[1] - _iw15_w_anchor[1]), _w2w3_al)
    iw15_nw = _offset(_iw15_w_anchor, _iw15_t_n, _w2w3_al)
    iw15 = BBox(w=_iw15_w_anchor[0], s=_iw15_w_anchor[1],
                e=_iw15_e_anchor[0], n=iw15_nw[1])

    # Dresser: 34" wide x 19" deep, 1" from IW1 south, 2" from IW15 west
    dresser_ne = line_isect(
        _offset(_iw15_w_anchor, -2.0 / 12.0, _w2w3_in), _w2w3_al,
        _offset(iw1_sw, 1.0 / 12.0, _w9w10_in), _w9w10_al)
    dresser_nw = _offset(dresser_ne, -34.0 / 12.0, _w2w3_in)
    dresser_se = _offset(dresser_ne, 19.0 / 12.0, _w9w10_in)
    dresser = BBox(w=dresser_nw[0], s=dresser_se[1],
                   e=dresser_ne[0], n=dresser_ne[1])

    # IW5 west end: meets IW14 SE corner
    iw5_w = iw14_se[0]

    # IW6: parallel to W6-W7, offset from W6
    _iw6_n_anchor = _offset(pts["W6"], IW6_OFFSET_FROM_W6, _w6w7_in)
    _iw6_s_anchor = _offset(_iw6_n_anchor, IW6_THICKNESS, _w6w7_in)
    iw6_ne = line_isect(_iw6_n_anchor, _w6w7_al, _iw2_w_anchor, _w2w3_al)
    iw6_se = line_isect(_iw6_s_anchor, _w6w7_al, _iw2_w_anchor, _w2w3_al)
    _neg_w6w7_al = (-_w6w7_al[0], -_w6w7_al[1])
    _iw6_n_ts = _line_poly_isects(inner_poly, iw6_ne, _neg_w6w7_al)
    iw6_nw = _offset(iw6_ne, max(_iw6_n_ts), _neg_w6w7_al)
    _iw6_s_ts = _line_poly_isects(inner_poly, iw6_se, _neg_w6w7_al)
    iw6_sw = _offset(iw6_se, max(_iw6_s_ts), _neg_w6w7_al)
    iw6_poly = [iw6_sw, iw6_se, iw6_ne, iw6_nw]
    iw6_n, iw6_s = iw6_ne[1], iw6_se[1]

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
