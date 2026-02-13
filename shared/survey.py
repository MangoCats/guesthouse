"""Survey computation: traverse, three-arc system, and inset path."""
import math
from typing import NamedTuple

from .types import Point, LineSeg, ArcSeg, Segment
from .geometry import (
    left_norm, off_pt, line_isect,
    line_circle_isect_min_t_gt, line_circle_isect_min_abs_t,
    circle_circle_isect,
)

# ============================================================
# Traverse Computation
# ============================================================
def compute_traverse() -> tuple[dict[str, Point], Point]:
    """Compute traverse from raw survey legs.

    Returns (pts, p3_trav) where pts is a dict with keys
    P3/POB/P2/P4/P5 (P3-based coordinates) and p3_trav is
    the raw P3 position needed for SVG transform calibration.
    """
    legs = [(257,53,45,19,1.0),(180,54,31,26,11.0),(93,36,7,31,10.5),
            (56,36,31,13,2.5),(317,11,44,34,11.5)]
    _trav = [(0.0, 0.0)]
    for deg, mn, sec, ft, inch in legs:
        brg = deg + mn/60.0 + sec/3600.0
        dist_in = ft * 12 + inch
        brg_rad = math.radians(brg)
        dE = dist_in * math.sin(brg_rad); dN = dist_in * math.cos(brg_rad)
        last = _trav[-1]; _trav.append((last[0]+dE, last[1]+dN))
    _trav_ft = [(e/12, n/12) for e, n in _trav[:5]]
    _trav_ft[2] = (-19.1177, _trav_ft[3][1])
    _trav_ft[1] = (_trav_ft[2][0], _trav_ft[2][1] + 29.0)

    _p3_trav = _trav_ft[2]
    pts = {}
    pts["P3"]  = (0.0, 0.0)
    pts["POB"] = (_trav_ft[0][0] - _p3_trav[0], _trav_ft[0][1] - _p3_trav[1])
    pts["P2"]  = (_trav_ft[1][0] - _p3_trav[0], _trav_ft[1][1] - _p3_trav[1])
    pts["P4"]  = (_trav_ft[3][0] - _p3_trav[0], _trav_ft[3][1] - _p3_trav[1])
    pts["P5"]  = (_trav_ft[4][0] - _p3_trav[0], _trav_ft[4][1] - _p3_trav[1])
    return pts, _p3_trav

# ============================================================
# Three-Arc System
# ============================================================
def compute_three_arc(pts: dict[str, Point]) -> dict[str, float]:
    """Compute three-arc boundary system.

    Mutates pts adding T1/TC1/T2/TC2/PA/T3/TC3/PX.
    Returns dict with R1, R2, R3, uE, uN, nE, nN.
    """
    dE_l = pts["P5"][0]-pts["POB"][0]; dN_l = pts["P5"][1]-pts["POB"][1]
    L = math.sqrt(dE_l**2+dN_l**2)
    uE, uN = dE_l/L, dN_l/L
    nE, nN = -uN, uE

    R1, R2 = 10.0, 12.5
    T1_dist, T2_dist = 26.5, 5.75
    pts["T1"] = (pts["POB"][0]+T1_dist*uE, pts["POB"][1]+T1_dist*uN)
    pts["TC1"] = (pts["T1"][0]+R1*nE, pts["T1"][1]+R1*nN)
    pts["T2"] = (pts["POB"][0]+T2_dist*uE, pts["POB"][1]+T2_dist*uN)
    pts["TC2"] = (pts["T2"][0]+R2*nE, pts["T2"][1]+R2*nN)

    # PA: circle-circle intersection
    dx_cc = pts["TC2"][0]-pts["TC1"][0]; dy_cc = pts["TC2"][1]-pts["TC1"][1]
    d_cc = math.sqrt(dx_cc**2+dy_cc**2)
    a_cc = (R1**2-R2**2+d_cc**2)/(2*d_cc); h_cc = math.sqrt(R1**2-a_cc**2)
    ux_cc, uy_cc = dx_cc/d_cc, dy_cc/d_cc
    Mx, My = pts["TC1"][0]+a_cc*ux_cc, pts["TC1"][1]+a_cc*uy_cc
    I1 = (Mx+h_cc*(-uy_cc), My+h_cc*ux_cc); I2 = (Mx-h_cc*(-uy_cc), My-h_cc*ux_cc)
    ang_T2_C2 = math.atan2(pts["T2"][1]-pts["TC2"][1], pts["T2"][0]-pts["TC2"][0])
    def ccw_a(s, e): return (e-s)%(2*math.pi)
    s1 = ccw_a(ang_T2_C2, math.atan2(I1[1]-pts["TC2"][1], I1[0]-pts["TC2"][0]))
    s2 = ccw_a(ang_T2_C2, math.atan2(I2[1]-pts["TC2"][1], I2[0]-pts["TC2"][0]))
    pts["PA"] = I1 if s1 < s2 else I2

    R3 = 11.0
    T3_dist_from_P3 = 17.911244
    pts["T3"] = (pts["P3"][0]+T3_dist_from_P3, pts["P3"][1])
    pts["TC3"] = (pts["T3"][0], pts["T3"][1]-R3)

    dxL = pts["P4"][0]-pts["P5"][0]; dyL = pts["P4"][1]-pts["P5"][1]
    pts["PX"] = line_circle_isect_min_t_gt(pts["P5"], (dxL, dyL), pts["TC3"], R3, 1.0)

    return {"R1": R1, "R2": R2, "R3": R3, "uE": uE, "uN": uN, "nE": nE, "nN": nN}

# ============================================================
# Inset Path Computation
# ============================================================
class InsetResult(NamedTuple):
    pts_update: dict[str, Point]   # PiOB, Pi2, Pi3, Pi4, Pi5, Ti1-3, PiX, Ai2
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
    d_e6 = (pts["P5"][0]-pts["P4"][0], pts["P5"][1]-pts["P4"][1])
    d_e7 = (pts["T1"][0]-pts["P5"][0], pts["T1"][1]-pts["P5"][1])
    d_e10 = (pts["POB"][0]-pts["T2"][0], pts["POB"][1]-pts["T2"][1])

    ln1 = left_norm(pts["POB"], pts["P2"])
    ln2 = left_norm(pts["P2"], pts["P3"])
    ln3 = left_norm(pts["P3"], pts["T3"])
    ln5 = left_norm(pts["PX"], pts["P4"])
    ln6 = left_norm(pts["P4"], pts["P5"])
    ln7 = left_norm(pts["P5"], pts["T1"])
    ln10 = left_norm(pts["T2"], pts["POB"])

    o1 = off_pt(pts["POB"], ln1, delta)
    o2 = off_pt(pts["P2"], ln2, delta)
    o3 = off_pt(pts["P3"], ln3, delta)
    o5 = off_pt(pts["PX"], ln5, delta)
    o6 = off_pt(pts["P4"], ln6, delta)
    o7 = off_pt(pts["P5"], ln7, delta)
    o10 = off_pt(pts["T2"], ln10, delta)

    update = {}
    update["PiOB"] = line_isect(o10, d_e10, o1, d_e1)
    update["Pi2"] = line_isect(o1, d_e1, o2, d_e2)
    update["Pi3"] = line_isect(o2, d_e2, o3, d_e3)
    update["Pi4"] = off_pt(pts["P4"], ln5, delta)  # collinear PX-P4-P5
    update["Pi5"] = line_isect(o6, d_e6, o7, d_e7)

    update["Ti3"] = (pts["TC3"][0], pts["P3"][1] + delta)
    update["Ti1"] = (pts["T1"][0] - delta*nE, pts["T1"][1] - delta*nN)
    update["Ti2"] = (pts["T2"][0] - delta*nE, pts["T2"][1] - delta*nN)

    update["PiX"] = line_circle_isect_min_abs_t(o5, d_e5, pts["TC3"], R3i)
    update["Ai2"] = circle_circle_isect(pts["TC1"], R1i, pts["TC2"], R2i, near=pts["PA"])

    inset_segs = [
        LineSeg("PiOB", "Pi2"), LineSeg("Pi2", "Pi3"), LineSeg("Pi3", "Ti3"),
        ArcSeg("Ti3", "PiX", "TC3", R3i, "CW", 60),
        LineSeg("PiX", "Pi4"), LineSeg("Pi4", "Pi5"), LineSeg("Pi5", "Ti1"),
        ArcSeg("Ti1", "Ai2", "TC1", R1i, "CW", 60),
        ArcSeg("Ai2", "Ti2", "TC2", R2i, "CW", 60),
        LineSeg("Ti2", "PiOB"),
    ]

    return InsetResult(pts_update=update, inset_segs=inset_segs,
                       R1i=R1i, R2i=R2i, R3i=R3i)
