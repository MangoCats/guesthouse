"""Interior layout computation — rooms, walls, appliances, furniture."""
from typing import NamedTuple

from shared.types import Point, BBox
from shared.geometry import horiz_isects
from floorplan.constants import (
    WALL_6IN, WALL_4IN, WALL_3IN,
    APPLIANCE_WIDTH, APPLIANCE_DEPTH, APPLIANCE_OFFSET_E,
    APPLIANCE_OFFSET_N, APPLIANCE_GAP,
    COUNTER_DEPTH, COUNTER_LENGTH, COUNTER_NW_RADIUS, COUNTER_GAP,
    BEDROOM_WIDTH, CLOSET_WIDTH,
    BED_WIDTH, BED_LENGTH, BED_OFFSET_N,
    IW1_OFFSET_N, IW2_OFFSET_E, WALL_SOUTH_N,
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
    wall_south_n: float
    # IW8 (L-shaped, east closet wall)
    iw8: list[Point]
    iw8_w: float
    iw8_e: float
    # Closet 1
    cl1_top: float
    # Bed
    bed: BBox
    bed_cx: float
    # IW5 (3" thick, horizontal in office)
    iw5: BBox
    # IW6 (1" thick, horizontal above kitchen)
    iw6_poly: list[Point]
    iw6_n: float
    iw6_s: float


def compute_interior_layout(pts, inner_poly) -> InteriorLayout:
    """Compute interior layout positions.

    pts must contain W-series (W0-W21) and F-series (F0-F21).
    """
    iw1_s = pts["W0"][1] + IW1_OFFSET_N
    iw1_n = iw1_s + WALL_6IN
    si = horiz_isects(inner_poly, iw1_s)
    ni = horiz_isects(inner_poly, iw1_n)
    iw1 = [(min(si), iw1_s), (max(si), iw1_s), (max(ni), iw1_n), (min(ni), iw1_n)]

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
    ctr_n = ctr_s + COUNTER_LENGTH
    ctr_nw_r = COUNTER_NW_RADIUS

    iw7_poly = [(ctr_e, ctr_s), (ctr_e + WALL_3IN, ctr_s), (ctr_e + WALL_3IN, ctr_n),
                (ctr_e + WALL_3IN + CLOSET_WIDTH, ctr_n),
                (ctr_e + WALL_3IN + CLOSET_WIDTH, ctr_n + WALL_3IN),
                (ctr_e, ctr_n + WALL_3IN)]
    iw3_w = ctr_e + WALL_3IN + CLOSET_WIDTH
    iw3_e = iw3_w + WALL_4IN
    iw3_s = ctr_s
    iw3_n = iw1_s
    iw4_w = iw3_e + BEDROOM_WIDTH
    iw4_e = iw4_w + WALL_4IN
    wall_south_n = WALL_SOUTH_N
    cl1_top = ctr_n - 1.0
    iw8_w = iw4_e + CLOSET_WIDTH
    iw8_e = iw8_w + WALL_3IN
    iw8_poly = [(iw4_e, cl1_top + WALL_3IN), (iw8_e, cl1_top + WALL_3IN), (iw8_e, wall_south_n),
                (iw8_w, wall_south_n), (iw8_w, cl1_top), (iw4_e, cl1_top)]

    bed_cx = (iw3_e + iw4_w) / 2
    bed_w = bed_cx - BED_WIDTH / 2
    bed_e = bed_cx + BED_WIDTH / 2
    bed_s = ctr_s + BED_OFFSET_N
    bed_n = bed_s + BED_LENGTH

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
        iw4_w=iw4_w, iw4_e=iw4_e, wall_south_n=wall_south_n,
        iw8=iw8_poly, iw8_w=iw8_w, iw8_e=iw8_e,
        cl1_top=cl1_top,
        bed=BBox(w=bed_w, s=bed_s, e=bed_e, n=bed_n), bed_cx=bed_cx,
        iw5=BBox(w=iw5_w, s=iw5_s, e=iw5_e, n=iw5_n),
        iw6_poly=iw6_poly, iw6_n=iw6_n, iw6_s=iw6_s,
    )
