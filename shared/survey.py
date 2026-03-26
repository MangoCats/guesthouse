"""Survey computation: traverse, three-arc system, and inset path."""
import math
from typing import NamedTuple

from .types import Point, LineSeg, ArcSeg, Segment
from .geometry import (
    left_norm, off_pt, line_isect,
    line_circle_isect_min_t_gt, line_circle_isect_min_abs_t,
    circle_circle_isect,
)

# FC (building center) position in P3-based coordinates.
# Used to shift survey coordinates from P3 origin to FC origin.
FC_IN_P3 = (18.5141152720, 13.3968094375)

# Coordinate rotation: CCW angle from FC-based pre-rotation to the primary
# working frame.  Best-fit alignment at C17 = arctan(7/12): T3 on F18-F1,
# F12-F13 tangent to TC1 arc, F16-F17 on Pi5-PiX line (RMS ≈ 0.085 ft).
COORD_ROTATION = 0.0015153784


def rotate_pts(pts: dict[str, 'Point'], angle: float) -> None:
    """Rotate all points in pts CCW by angle (radians) around the origin."""
    c, s = math.cos(angle), math.sin(angle)
    for k in list(pts):
        e, n = pts[k]
        pts[k] = (e * c - n * s, e * s + n * c)

# ============================================================
# Traverse Computation
# ============================================================
# Raw traverse leg accumulation (all inputs are hardcoded survey constants).
def _accumulate_legs():
    legs = [(257,53,45,19,1.0),(180,54,31,26,11.0),(93,36,7,31,10.5),
            (56,36,31,13,2.5),(317,11,44,34,11.5)]
    trav = [(0.0, 0.0)]
    for deg, mn, sec, ft, inch in legs:
        brg = deg + mn/60.0 + sec/3600.0
        dist_in = ft * 12 + inch
        brg_rad = math.radians(brg)
        dE = dist_in * math.sin(brg_rad); dN = dist_in * math.cos(brg_rad)
        last = trav[-1]; trav.append((last[0]+dE, last[1]+dN))
    trav_ft = [(e/12, n/12) for e, n in trav[:5]]
    trav_ft[2] = (-19.1177, trav_ft[3][1])
    trav_ft[1] = (trav_ft[2][0], trav_ft[2][1] + 29.0)
    return trav_ft

_TRAV_FT = _accumulate_legs()

# Raw P3 traverse position — a constant needed only for SVG calibration.
_P3_TRAV = _TRAV_FT[2]


def compute_traverse() -> dict[str, Point]:
    """Compute traverse producing FC-based coordinates (building center = origin)."""
    p3 = _P3_TRAV
    pts = {}
    pts["P3"]  = (0.0, 0.0)
    pts["POB"] = (_TRAV_FT[0][0] - p3[0], _TRAV_FT[0][1] - p3[1])
    pts["P2"]  = (_TRAV_FT[1][0] - p3[0], _TRAV_FT[1][1] - p3[1])
    pts["P4"]  = (_TRAV_FT[3][0] - p3[0], _TRAV_FT[3][1] - p3[1])
    pts["P5"]  = (_TRAV_FT[4][0] - p3[0], _TRAV_FT[4][1] - p3[1])

    # Shift all points from P3 origin to FC (building center) origin
    for k in list(pts):
        pts[k] = (pts[k][0] - FC_IN_P3[0], pts[k][1] - FC_IN_P3[1])

    return pts

# ============================================================
# Three-Arc System
# ============================================================
def compute_three_arc(pts: dict[str, Point]) -> dict[str, float]:
    """Compute three-arc boundary system.

    Mutates pts adding T1/TC1/T2/TC2/PA/T3/TC3/PX.
    Returns dict with R1, R2, R3, uE, uN, nE, nN.
    """
    dE_l = pts["P5"][0]-pts["POB"][0]; dN_l = pts["P5"][1]-pts["POB"][1]
    L = math.hypot(dE_l, dN_l)
    uE, uN = dE_l/L, dN_l/L
    nE, nN = -uN, uE

    R1, R2 = 10.0, 12.5
    T1_dist, T2_dist = 26.5, 5.75
    pts["T1"] = (pts["POB"][0]+T1_dist*uE, pts["POB"][1]+T1_dist*uN)
    pts["TC1"] = (pts["T1"][0]+R1*nE, pts["T1"][1]+R1*nN)
    pts["T2"] = (pts["POB"][0]+T2_dist*uE, pts["POB"][1]+T2_dist*uN)
    pts["TC2"] = (pts["T2"][0]+R2*nE, pts["T2"][1]+R2*nN)

    # PA: circle-circle intersection (nearest to T2 on arc 2)
    pts["PA"] = circle_circle_isect(pts["TC1"], R1, pts["TC2"], R2, near=pts["T2"])

    R3 = 11.0
    T3_dist_from_P3 = 17.911244
    pts["T3"] = (pts["P3"][0]+T3_dist_from_P3, pts["P3"][1])
    pts["TC3"] = (pts["T3"][0], pts["T3"][1]-R3)

    dxL = pts["P4"][0]-pts["P5"][0]; dyL = pts["P4"][1]-pts["P5"][1]
    pts["PX"] = line_circle_isect_min_t_gt(pts["P5"], (dxL, dyL), pts["TC3"], R3, 1.0)

    return {"R1": R1, "R2": R2, "R3": R3, "uE": uE, "uN": uN, "nE": nE, "nN": nN}

def compute_pt1(pts: dict[str, Point], R1: float) -> Point:
    """Compute PT1: tangency point where the TC1 arc meets the P4-P5 line extension.

    PT1 is the foot of the perpendicular from TC1 to the P4-P5 line, projected
    outward onto the arc (radius R1).  The P4-P5 line is tangent to the arc at
    PT1 within survey measurement precision.
    """
    p4, p5, tc1 = pts["P4"], pts["P5"], pts["TC1"]
    dx, dy = p5[0] - p4[0], p5[1] - p4[1]
    t = ((tc1[0] - p4[0]) * dx + (tc1[1] - p4[1]) * dy) / (dx * dx + dy * dy)
    foot_e = p4[0] + t * dx
    foot_n = p4[1] + t * dy
    dir_e, dir_n = foot_e - tc1[0], foot_n - tc1[1]
    dir_len = math.hypot(dir_e, dir_n)
    return (tc1[0] + R1 * dir_e / dir_len, tc1[1] + R1 * dir_n / dir_len)


# ============================================================
# Inset Path Computation
# ============================================================
class InsetResult(NamedTuple):
    pts_update: dict[str, Point]   # PiOB, Pi2, Pi3, Pi4, Pi5, PTi1, Ti1-3, PiX, Ai2
    inset_segs: list[Segment]
    R1i: float; R2i: float; R3i: float

def compute_inset(
    pts: dict[str, Point], R1: float, R2: float, R3: float,
    nE: float, nN: float, delta: float = 0.5,
) -> InsetResult:
    """Compute inset path (6" inside outer path).

    Pure function — does not mutate pts. Returns InsetResult with
    pts_update dict, inset_segs list, and inset radii.
    """
    R1i, R2i, R3i = R1+delta, R2+delta, R3+delta

    d_e1 = (pts["P2"][0]-pts["POB"][0], pts["P2"][1]-pts["POB"][1])
    d_e2 = (pts["P3"][0]-pts["P2"][0], pts["P3"][1]-pts["P2"][1])
    d_e3 = (pts["T3"][0]-pts["P3"][0], pts["T3"][1]-pts["P3"][1])
    d_e5 = (pts["P4"][0]-pts["PX"][0], pts["P4"][1]-pts["PX"][1])
    d_e10 = (pts["POB"][0]-pts["T2"][0], pts["POB"][1]-pts["T2"][1])

    ln1 = left_norm(pts["POB"], pts["P2"])
    ln2 = left_norm(pts["P2"], pts["P3"])
    ln3 = left_norm(pts["P3"], pts["T3"])
    ln5 = left_norm(pts["PX"], pts["P4"])
    ln10 = left_norm(pts["T2"], pts["POB"])

    o1 = off_pt(pts["POB"], ln1, delta)
    o2 = off_pt(pts["P2"], ln2, delta)
    o3 = off_pt(pts["P3"], ln3, delta)
    o5 = off_pt(pts["PX"], ln5, delta)
    o10 = off_pt(pts["T2"], ln10, delta)

    update = {}
    update["PiOB"] = line_isect(o10, d_e10, o1, d_e1)
    update["Pi2"] = line_isect(o1, d_e1, o2, d_e2)
    update["Pi3"] = line_isect(o2, d_e2, o3, d_e3)
    update["Pi4"] = off_pt(pts["P4"], ln5, delta)  # collinear PX-P4-P5-PT1
    update["Pi5"] = off_pt(pts["P5"], ln5, delta)  # collinear with Pi4 and PTi1

    # PTi1: intersection of Pi4-Pi5 line (extended beyond Pi5) with TC1 arc (radius R1i)
    # Line: Pi4 + t*(Pi5-Pi4); t=1 is Pi5; t>1 is the extension toward the arc
    _pi4 = update["Pi4"]
    _pi5 = update["Pi5"]
    _d56 = (_pi5[0] - _pi4[0], _pi5[1] - _pi4[1])
    update["PTi1"] = line_circle_isect_min_t_gt(
        _pi4, _d56, pts["TC1"], R1i, 1.0
    )

    update["Ti3"] = (pts["TC3"][0], pts["P3"][1] + delta)
    update["Ti1"] = (pts["T1"][0] - delta*nE, pts["T1"][1] - delta*nN)
    update["Ti2"] = (pts["T2"][0] - delta*nE, pts["T2"][1] - delta*nN)

    update["PiX"] = line_circle_isect_min_abs_t(o5, d_e5, pts["TC3"], R3i)
    update["Ai2"] = circle_circle_isect(pts["TC1"], R1i, pts["TC2"], R2i, near=pts["PA"])

    inset_segs = [
        LineSeg("PiOB", "Pi2"), LineSeg("Pi2", "Pi3"), LineSeg("Pi3", "Ti3"),
        ArcSeg("Ti3", "PiX", "TC3", R3i, "CW", 60),
        LineSeg("PiX", "Pi4"), LineSeg("Pi4", "Pi5"), LineSeg("Pi5", "PTi1"),
        ArcSeg("PTi1", "Ti1", "TC1", R1i, "CW", 60),
        ArcSeg("Ti1", "Ai2", "TC1", R1i, "CW", 60),
        ArcSeg("Ai2", "Ti2", "TC2", R2i, "CW", 60),
        LineSeg("Ti2", "PiOB"),
    ]

    return InsetResult(pts_update=update, inset_segs=inset_segs,
                       R1i=R1i, R2i=R2i, R3i=R3i)
