"""Interior layout computation — rooms, walls, appliances, furniture."""
import math
from typing import NamedTuple

from shared.types import Point, BBox
from shared.geometry import horiz_isects
from floorplan.constants import (
    WALL_6IN, WALL_4IN, WALL_3IN,
    APPLIANCE_WIDTH, APPLIANCE_DEPTH, APPLIANCE_OFFSET_E,
    APPLIANCE_OFFSET_N, APPLIANCE_GAP,
    COUNTER_DEPTH, COUNTER_GAP,
    BEDROOM_WIDTH, CLOSET_WIDTH, CLOSET2_WIDTH,
    BED_WIDTH, BED_LENGTH,
    O9_HALF_WIDTH,
    IW1_DIST_FROM_NORTH, IW1_WEST_OFFSET_E, IW2_OFFSET_E, WALL_SOUTH_N,
    IW5_OFFSET_N, IW6_THICKNESS, IW6_OFFSET_N,
)


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
    ctr_nw_r: float
    # Wall thicknesses
    iwt3: float
    iwt4: float
    # IW7 (L-shaped, west/north walls of closet)
    iw7: list[Point]
    # Interior wall 3 (IW3) — west bedroom wall
    iw3: BBox
    # Interior wall 4 (IW4) — east bedroom wall
    iw4_w: float
    iw4_e: float
    iw4_s: float
    wall_south_n: float
    # Closet 1
    cl1_top: float
    # Bed (rotated polygon, long sides perpendicular to W20-W20a)
    bed_poly: list[Point]  # [SW, SE, NE, NW]
    # IW9 (vertical, 4" — old IW3 position, south of IW7 L north face)
    iw9: BBox
    # IW10 (4" thick, horizontal — closet north wall, IW3.e to IW9.e)
    iw10: BBox
    # IW11 (4" thick, N-S — south extension below IW4)
    iw11: BBox         # bounding box (for labels/dimensions)
    iw11_poly: list[Point]  # actual polygon (SE corner on W20-W20a)
    # IW12 (4" thick, perpendicular to IW11 — connects IW11 NW to IW4 west)
    iw12: BBox
    iw12_poly: list[Point]  # actual polygon [SW, SE, NE, NW]
    # IW5 (3" thick, horizontal in office)
    iw5: BBox
    # IW8 (6" thick, horizontal — west extension of IW1, W1-W2 to IW1 west end)
    iw8: BBox
    # IW6 (1" thick, horizontal above kitchen)
    iw6_poly: list[Point]
    iw6_n: float
    iw6_s: float


def compute_interior_layout(pts, inner_poly) -> InteriorLayout:
    """Compute interior layout positions.

    pts must contain W-series (W0-W21) and F-series (F0-F21).
    """
    iw1_n = pts["W9"][1] - IW1_DIST_FROM_NORTH
    iw1_s = iw1_n - WALL_6IN
    si = horiz_isects(inner_poly, iw1_s)
    ni = horiz_isects(inner_poly, iw1_n)
    iw1_w = pts["W1"][0] + IW1_WEST_OFFSET_E
    iw1 = [(iw1_w, iw1_s), (max(si), iw1_s), (max(ni), iw1_n), (iw1_w, iw1_n)]

    iw2_w = pts["W1"][0] + IW2_OFFSET_E
    iw2_e = iw2_w + WALL_6IN
    iw2_s = iw1_n
    iw2_n = pts["W6"][1]

    dryer_w = pts["W1"][0] + APPLIANCE_OFFSET_E
    dryer_s = pts["W0"][1] + APPLIANCE_OFFSET_N
    dryer_e = dryer_w + APPLIANCE_WIDTH
    dryer_n = dryer_s + APPLIANCE_DEPTH
    washer_w = dryer_w
    washer_s = dryer_n + APPLIANCE_GAP
    washer_e = dryer_e
    washer_n = washer_s + APPLIANCE_DEPTH

    ctr_w = dryer_e + COUNTER_GAP
    ctr_e = ctr_w + COUNTER_DEPTH
    ctr_s = pts["W0"][1]
    iw7_n = ctr_s + 6.0  # 6' north of W21-W0 face
    ctr_n = iw7_n         # counter north = IW10 south face
    ctr_nw_r = 0
    iw7_poly = [(ctr_e, ctr_s), (ctr_e + WALL_3IN, ctr_s),
                (ctr_e + WALL_3IN, iw7_n), (ctr_e, iw7_n)]
    iw3_e = ctr_e + WALL_3IN          # east face = IW7 east face
    iw3_w = iw3_e - WALL_4IN
    iw3_s = iw7_n
    iw3_n = iw1_s
    iw9_w = ctr_e + WALL_3IN + CLOSET_WIDTH
    iw9_e = iw9_w + WALL_4IN
    iw9_s = ctr_s
    iw9_n = iw7_n
    # IW10: 4" thick E-W wall from IW3 east face to IW9 east face
    iw10_w = iw3_e
    iw10_e = iw9_e
    iw10_s = iw7_n
    iw10_n = iw7_n + WALL_4IN

    iw4_w = iw9_e + BEDROOM_WIDTH + WALL_4IN + CLOSET2_WIDTH
    iw4_e = iw4_w + WALL_4IN
    iw4_s = WALL_SOUTH_N
    wall_south_n = WALL_SOUTH_N
    cl1_top = iw7_n - 1.0

    # IW11: 4" thick, normal to W20-W20a, 6' long
    # SE corner: circle(IW4_SW, 32") ∩ W20-W20a
    _iw4_sw = (iw4_w, iw4_s)
    _w20 = pts["W20"]
    _w20a = pts["W20a"]
    _dE = _w20a[0] - _w20[0]
    _dN = _w20a[1] - _w20[1]
    _seg_len = math.sqrt(_dE**2 + _dN**2)
    _uE = _w20[0] - _iw4_sw[0]
    _uN = _w20[1] - _iw4_sw[1]
    _r = 32.0 / 12.0
    _qa = _dE**2 + _dN**2
    _qb = 2 * (_uE * _dE + _uN * _dN)
    _qc = _uE**2 + _uN**2 - _r**2
    _disc = _qb**2 - 4 * _qa * _qc
    _t = (-_qb + math.sqrt(_disc)) / (2 * _qa)  # westward along W20-W20a
    iw11_se = (_w20[0] + _t * _dE, _w20[1] + _t * _dN)
    # Unit vectors: along W20-W20a and inward normal
    _along_E = _dE / _seg_len
    _along_N = _dN / _seg_len
    _norm_E = _along_N    # right normal = inward
    _norm_N = -_along_E
    _iw11_thick = WALL_4IN
    _iw11_len = 6.0  # 6'
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
    iw12_sw = (iw11_nw[0] - _iw12_shorten * _along_E,
               iw11_nw[1] - _iw12_shorten * _along_N)
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

    # Bed: rotated, long sides perpendicular to W20-W20a
    # SE corner = 1" past O9 NW corner along W20-W20a, 2" from wall
    _dE9 = pts["F20a"][0] - pts["F20"][0]
    _dN9 = pts["F20a"][1] - pts["F20"][1]
    _seg9_len = math.sqrt(_dE9**2 + _dN9**2)
    _t_sw9 = ((iw11_sw[0] - pts["F20"][0]) * _dE9
              + (iw11_sw[1] - pts["F20"][1]) * _dN9) / (_dE9**2 + _dN9**2)
    _ts9 = _t_sw9 + 4.0 / 12.0 / _seg9_len
    _te9 = _ts9 + 2 * O9_HALF_WIDTH / _seg9_len
    _bed_t = _te9 + 1.0 / 12.0 / _seg_len  # 1" past O9 NW along W20-W20a
    _bed_se_wall = (_w20[0] + _bed_t * (_w20a[0] - _w20[0]),
                    _w20[1] + _bed_t * (_w20a[1] - _w20[1]))
    bed_se = (_bed_se_wall[0] + 2.0 / 12.0 * _norm_E,
              _bed_se_wall[1] + 2.0 / 12.0 * _norm_N)
    bed_sw = (bed_se[0] + BED_WIDTH * _along_E,
              bed_se[1] + BED_WIDTH * _along_N)
    bed_ne = (bed_se[0] + BED_LENGTH * _norm_E,
              bed_se[1] + BED_LENGTH * _norm_N)
    bed_nw = (bed_sw[0] + BED_LENGTH * _norm_E,
              bed_sw[1] + BED_LENGTH * _norm_N)
    bed_poly = [bed_sw, bed_se, bed_ne, bed_nw]

    # IW8: 6" thick, horizontal, from W1-W2 face to IW1 west end
    iw8_w = pts["W1"][0]
    iw8_e = iw1_w
    iw8 = BBox(w=iw8_w, s=iw1_s + 18.0 / 12.0, e=iw8_e, n=iw1_n + 18.0 / 12.0)

    # IW5: 3" thick, north face IW5_OFFSET_N south of IW1 south face
    iw5_n = iw1_s - IW5_OFFSET_N
    iw5_s = iw5_n - WALL_3IN
    iw5_w = iw4_e
    iw5_e = pts["W15"][0]

    # IW6: IW6_THICKNESS thick, south face IW6_OFFSET_N south of W6
    iw6_n = pts["W6"][1] - IW6_OFFSET_N
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
        dryer=BBox(w=dryer_w, s=dryer_s, e=dryer_e, n=dryer_n),
        washer=BBox(w=washer_w, s=washer_s, e=washer_e, n=washer_n),
        ctr=BBox(w=ctr_w, s=ctr_s, e=ctr_e, n=ctr_n), ctr_nw_r=ctr_nw_r,
        iwt3=WALL_3IN, iwt4=WALL_4IN,
        iw7=iw7_poly,
        iw3=BBox(w=iw3_w, s=iw3_s, e=iw3_e, n=iw3_n),
        iw9=BBox(w=iw9_w, s=iw9_s, e=iw9_e, n=iw9_n),
        iw10=BBox(w=iw10_w, s=iw10_s, e=iw10_e, n=iw10_n),
        iw4_w=iw4_w, iw4_e=iw4_e, iw4_s=iw4_s, wall_south_n=wall_south_n,
        iw11=BBox(w=iw11_w, s=iw11_s, e=iw11_e, n=iw11_n), iw11_poly=iw11_poly,
        iw12=BBox(w=iw12_w, s=iw12_s, e=iw12_e, n=iw12_n), iw12_poly=iw12_poly,
        cl1_top=cl1_top,
        bed_poly=bed_poly,
        iw8=iw8,
        iw5=BBox(w=iw5_w, s=iw5_s, e=iw5_e, n=iw5_n),
        iw6_poly=iw6_poly, iw6_n=iw6_n, iw6_s=iw6_s,
    )
