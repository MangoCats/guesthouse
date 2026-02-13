"""Interior layout computation — rooms, walls, appliances, furniture."""
from typing import NamedTuple

from shared.types import Point
from shared.geometry import horiz_isects
from floorplan.constants import (
    WALL_6IN, WALL_4IN, WALL_3IN,
    APPLIANCE_WIDTH, APPLIANCE_DEPTH, APPLIANCE_OFFSET_E,
    APPLIANCE_OFFSET_N, APPLIANCE_GAP,
    COUNTER_DEPTH, COUNTER_LENGTH, COUNTER_NW_RADIUS, COUNTER_GAP,
    BEDROOM_WIDTH, CLOSET_WIDTH,
    BED_WIDTH, BED_LENGTH, BED_OFFSET_N,
    IW1_OFFSET_N, IW2_OFFSET_E, WALL_SOUTH_N,
)


class InteriorLayout(NamedTuple):
    """Interior layout positions for walls, appliances, and furniture."""
    # Interior wall 1 (IW1) — horizontal wall separating utility/bedroom zones
    iw1: list[Point]
    iw1_s: float
    iw1_n: float
    iwt: float
    # Interior wall 2 (IW2) — vertical wall west of utility area
    iw2_w: float
    iw2_e: float
    iw2_s: float
    iw2_n: float
    # Dryer
    dryer_w: float
    dryer_s: float
    dryer_e: float
    dryer_n: float
    # Washer
    washer_w: float
    washer_s: float
    washer_e: float
    washer_n: float
    # Counter
    ctr_w: float
    ctr_e: float
    ctr_s: float
    ctr_n: float
    ctr_nw_r: float
    # Wall thicknesses
    iwt3: float
    iwt4: float
    # IW7 (L-shaped, west/north walls of closet)
    iw7: list[Point]
    # Interior wall 3 (IW3) — west bedroom wall
    iw3_w: float
    iw3_e: float
    iw3_s: float
    iw3_n: float
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
    bed_w: float
    bed_e: float
    bed_s: float
    bed_n: float
    bed_cx: float


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

    return InteriorLayout(
        iw1=iw1, iw1_s=iw1_s, iw1_n=iw1_n, iwt=WALL_6IN,
        iw2_w=iw2_w, iw2_e=iw2_e, iw2_s=iw2_s, iw2_n=iw2_n,
        dryer_w=dryer_w, dryer_s=dryer_s, dryer_e=dryer_e, dryer_n=dryer_n,
        washer_w=washer_w, washer_s=washer_s, washer_e=washer_e, washer_n=washer_n,
        ctr_w=ctr_w, ctr_e=ctr_e, ctr_s=ctr_s, ctr_n=ctr_n, ctr_nw_r=ctr_nw_r,
        iwt3=WALL_3IN, iwt4=WALL_4IN,
        iw7=iw7_poly,
        iw3_w=iw3_w, iw3_e=iw3_e, iw3_s=iw3_s, iw3_n=iw3_n,
        iw4_w=iw4_w, iw4_e=iw4_e, wall_south_n=wall_south_n,
        iw8=iw8_poly, iw8_w=iw8_w, iw8_e=iw8_e,
        cl1_top=cl1_top,
        bed_w=bed_w, bed_e=bed_e, bed_s=bed_s, bed_n=bed_n, bed_cx=bed_cx,
    )
