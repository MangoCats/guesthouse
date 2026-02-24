"""Distance-squared regression tests.

For every F-series outline point, every P-series survey/inset point,
every corner of every IW wall, every corner/center of every appliance
or piece of furniture, every corner of every rough opening, every corner
of every outer opening, and every door tip point, capture distance squared
to each of: F1, F6, F12, F15.

Values must match within +/- 0.000001 (in feet-squared).
"""
import math
import pytest
from shared.types import BBox
from shared.geometry import f8f9_corner_polyline
from floorplan.openings import compute_outer_openings, compute_rough_openings
from floorplan.gen_floorplan import compute_placement_points, compute_dimension_endpoints
from floorplan.geometry import walk_outline_chain, F2_E, F2_N
from floorplan.constants import (
    WALL_OUTER, SHELL_THICKNESS, AIR_GAP,
    F8F9_INNER_TURN_R, OPENING_INSIDE_RADIUS,
    O3_DOOR_WIDTH, O6_DOOR_WIDTH, O6_WIDTH,
    RO1_DOOR_WIDTH, RO2_DOOR_WIDTH, RO3_DOOR_WIDTH,
    RO4_DOOR_WIDTH, RO5_DOOR_WIDTH, RO6_DOOR_WIDTH, RO7_DOOR_WIDTH,
    IW4_RO_WIDTH, IW9_RO_WIDTH, IW11_RO_WIDTH,
)

TOL = 1e-6  # tolerance in ft^2

from conftest import dist_sq as _dist_sq


def _bbox_corners(bbox):
    return [(bbox.w, bbox.s), (bbox.e, bbox.s), (bbox.e, bbox.n), (bbox.w, bbox.n)]


def _collect_all_points(pts, layout, radii):
    """Collect all named points from layout, openings, and door tips.

    Returns list of (name, point) tuples.
    """
    outer_openings = compute_outer_openings(pts, layout)
    rough_openings = compute_rough_openings(pts, layout)
    result = []

    # ---- F-series outline points (from chain walk) ----
    chain_pts = walk_outline_chain()
    f_names = [f"F{i}" for i in range(1, 21) if i not in (3, 4)] + ["F11a", "F11b"]
    for name in f_names:
        result.append((name, chain_pts[name]))

    # ---- FC (building center) ----
    result.append(("FC", pts["FC"]))

    # ---- P-series survey points ----
    # Traverse stations
    for name in ["P3", "POB", "P2", "P4", "P5"]:
        result.append((name, pts[name]))
    # Three-arc points (centers, tangent pts, auxiliary)
    for name in ["T1", "TC1", "T2", "TC2", "PA", "T3", "TC3", "PX"]:
        result.append((name, pts[name]))
    # Inset points
    for name in ["PiOB", "Pi2", "Pi3", "Pi4", "Pi5",
                  "Ti1", "Ti2", "Ti3", "PiX", "Ai2"]:
        result.append((name, pts[name]))

    # ---- IW walls, appliances, furniture (all Wall type with .poly) ----
    for prefix, wall in [
        ("IW1", layout.iw1), ("IW2", layout.iw2),
        ("IW3", layout.iw3), ("IW4", layout.iw4),
        ("IW5", layout.iw5), ("IW6", layout.iw6),
        ("IW7", layout.iw7), ("IW8", layout.iw8),
        ("IW9", layout.iw9), ("IW11", layout.iw11),
        ("IW12", layout.iw12), ("IW16", layout.iw16),
        ("dryer", layout.dryer), ("washer", layout.washer),
    ]:
        for i, label in enumerate(["SW", "SE", "NE", "NW"]):
            result.append((f"{prefix}_{label}", wall.poly[i]))

    # ---- Counter clip polygon (4 or 5 pts) ----
    for i, p in enumerate(layout.ctr_clip):
        result.append((f"ctr_poly_{i}", p))

    # ---- Furniture (bed, dresser) ----
    for prefix, wall in [("bed", layout.bed), ("dresser", layout.dresser)]:
        for i, label in enumerate(["SW", "SE", "NE", "NW"]):
            result.append((f"{prefix}_{label}", wall.poly[i]))

    # ---- Rough openings ----
    for ro in rough_openings:
        if ro.poly:
            for i, label in enumerate(["SW", "SE", "NE", "NW"]):
                result.append((f"{ro.name}_{label}", ro.poly[i]))
        else:
            for label, pt in zip(["SW", "SE", "NE", "NW"], _bbox_corners(ro.bbox)):
                result.append((f"{ro.name}_{label}", pt))

    # ---- Outer openings (O1-O11) ----
    oo_labels = ["outer_start", "outer_end", "inner_end", "inner_start"]
    for oo in outer_openings:
        for i, label in enumerate(oo_labels):
            result.append((f"{oo.name}_{label}", oo.poly[i]))

    # ---- Door tips ----
    result.extend(_compute_door_tips(pts, layout, outer_openings, rough_openings))

    # ---- W8-W9 shell polyline special points ----
    result.extend(_collect_f8f9_shell_points(pts))

    # ---- Appliance/furniture placement points (standard variant) ----
    result.extend(compute_placement_points(pts, layout, radii))

    # ---- Dimension line endpoints (standard variant) ----
    result.extend(compute_dimension_endpoints(pts, layout, radii, bare=False))

    return result


def _compute_door_tips(pts, layout, outer_openings, rough_openings):
    """Compute all door tip points, matching gen_floorplan.py logic."""
    tips = []

    # O3 door tip
    o3 = [o for o in outer_openings if o.name == "O3"][0]
    _o3_os, _o3_oe = o3.poly[0], o3.poly[1]
    _o3_is, _o3_ie = o3.poly[3], o3.poly[2]
    _o3_dE = _o3_oe[0] - _o3_os[0]
    _o3_dN = _o3_oe[1] - _o3_os[1]
    _o3_len = math.sqrt(_o3_dE**2 + _o3_dN**2)
    _o3_along = (_o3_dE / _o3_len, _o3_dN / _o3_len)
    _o3_cross = (_o3_along[1], -_o3_along[0])
    _o3_me = ((_o3_oe[0] + _o3_ie[0]) / 2, (_o3_oe[1] + _o3_ie[1]) / 2)
    gap = (_o3_len - O3_DOOR_WIDTH) / 2
    hinge = (_o3_me[0] - gap * _o3_along[0], _o3_me[1] - gap * _o3_along[1])
    tips.append(("O3_door_tip",
                 (hinge[0] + O3_DOOR_WIDTH * _o3_cross[0],
                  hinge[1] + O3_DOOR_WIDTH * _o3_cross[1])))

    # O6 door tip (poly-based, rotation-safe)
    o6 = [o for o in outer_openings if o.name == "O6"][0]
    # End edge midpoint (wall midpoint at hinge side)
    _o6_end = ((o6.poly[1][0] + o6.poly[2][0]) / 2,
               (o6.poly[1][1] + o6.poly[2][1]) / 2)
    # Along direction (start→end on inner face)
    _o6_dx = o6.poly[1][0] - o6.poly[0][0]
    _o6_dy = o6.poly[1][1] - o6.poly[0][1]
    _o6_len = math.sqrt(_o6_dx**2 + _o6_dy**2)
    _o6_al = (_o6_dx / _o6_len, _o6_dy / _o6_len)
    # Inward direction (outer→inner at end edge)
    _o6_in_dx = o6.poly[1][0] - o6.poly[2][0]
    _o6_in_dy = o6.poly[1][1] - o6.poly[2][1]
    _o6_in_len = math.sqrt(_o6_in_dx**2 + _o6_in_dy**2)
    _o6_inward = (_o6_in_dx / _o6_in_len, _o6_in_dy / _o6_in_len)
    gap = (O6_WIDTH - O6_DOOR_WIDTH) / 2
    _o6_hinge = (_o6_end[0] - gap * _o6_al[0], _o6_end[1] - gap * _o6_al[1])
    tips.append(("O6_door_tip",
                 (_o6_hinge[0] + O6_DOOR_WIDTH * _o6_inward[0],
                  _o6_hinge[1] + O6_DOOR_WIDTH * _o6_inward[1])))

    # RO1 door tip (poly-based, rotation-safe)
    ro1 = [r for r in rough_openings if r.name == "RO1"][0]
    # Along direction: poly[0]→poly[1] (SW→SE)
    _ro1_dx = ro1.poly[1][0] - ro1.poly[0][0]
    _ro1_dy = ro1.poly[1][1] - ro1.poly[0][1]
    _ro1_len = math.sqrt(_ro1_dx**2 + _ro1_dy**2)
    _ro1_al = (_ro1_dx / _ro1_len, _ro1_dy / _ro1_len)
    # End edge midpoint (east): midpoint of poly[1] (SE) and poly[2] (NE)
    _ro1_end = ((ro1.poly[1][0] + ro1.poly[2][0]) / 2,
                (ro1.poly[1][1] + ro1.poly[2][1]) / 2)
    # Through-wall direction (NE→SE = toward south face)
    _ro1_sw_dx = ro1.poly[1][0] - ro1.poly[2][0]
    _ro1_sw_dy = ro1.poly[1][1] - ro1.poly[2][1]
    _ro1_sw_len = math.sqrt(_ro1_sw_dx**2 + _ro1_sw_dy**2)
    _ro1_swing = (_ro1_sw_dx / _ro1_sw_len, _ro1_sw_dy / _ro1_sw_len)
    ro1_gap = (_ro1_len - RO1_DOOR_WIDTH) / 2
    _ro1_hinge = (_ro1_end[0] - ro1_gap * _ro1_al[0],
                  _ro1_end[1] - ro1_gap * _ro1_al[1])
    tips.append(("RO1_door_tip",
                 (_ro1_hinge[0] + RO1_DOOR_WIDTH * _ro1_swing[0],
                  _ro1_hinge[1] + RO1_DOOR_WIDTH * _ro1_swing[1])))

    # RO2 door tip
    ro2 = [r for r in rough_openings if r.name == "RO2"][0]
    _i11_sw, _i11_se, _i11_ne, _i11_nw = layout.iw11.poly
    _i11_dx_n = _i11_ne[0] - _i11_se[0]; _i11_dy_n = _i11_ne[1] - _i11_se[1]
    _i11_ln = math.sqrt(_i11_dx_n**2 + _i11_dy_n**2)
    _i11_an = (_i11_dx_n / _i11_ln, _i11_dy_n / _i11_ln)
    _i11_dx_t = _i11_sw[0] - _i11_se[0]; _i11_dy_t = _i11_sw[1] - _i11_se[1]
    _i11_lt = math.sqrt(_i11_dx_t**2 + _i11_dy_t**2)
    _i11_at = (_i11_dx_t / _i11_lt, _i11_dy_t / _i11_lt)
    _ro2_gap = (IW4_RO_WIDTH - RO2_DOOR_WIDTH) / 2
    _ro2_n_ctr = ((ro2.poly[3][0] + ro2.poly[2][0]) / 2,
                  (ro2.poly[3][1] + ro2.poly[2][1]) / 2)
    hinge_e = _ro2_n_ctr[0] - _ro2_gap * _i11_an[0]
    hinge_n = _ro2_n_ctr[1] - _ro2_gap * _i11_an[1]
    tips.append(("RO2_door_tip",
                 (hinge_e - RO2_DOOR_WIDTH * _i11_at[0],
                  hinge_n - RO2_DOOR_WIDTH * _i11_at[1])))

    # RO3 door tip (poly-based, rotation-safe)
    ro3 = [r for r in rough_openings if r.name == "RO3"][0]
    # Length direction: poly[0]→poly[3] (SW→NW)
    _ro3_dx = ro3.poly[3][0] - ro3.poly[0][0]
    _ro3_dy = ro3.poly[3][1] - ro3.poly[0][1]
    _ro3_len = math.sqrt(_ro3_dx**2 + _ro3_dy**2)
    _ro3_len_u = (_ro3_dx / _ro3_len, _ro3_dy / _ro3_len)
    # End edge midpoint (north): midpoint of poly[2] (NE) and poly[3] (NW)
    _ro3_end = ((ro3.poly[2][0] + ro3.poly[3][0]) / 2,
                (ro3.poly[2][1] + ro3.poly[3][1]) / 2)
    # Door swing direction (toward west face): NE→NW = poly[2]→poly[3]
    _ro3_sw_dx = ro3.poly[3][0] - ro3.poly[2][0]
    _ro3_sw_dy = ro3.poly[3][1] - ro3.poly[2][1]
    _ro3_sw_len = math.sqrt(_ro3_sw_dx**2 + _ro3_sw_dy**2)
    _ro3_swing = (_ro3_sw_dx / _ro3_sw_len, _ro3_sw_dy / _ro3_sw_len)
    ro3_gap = (_ro3_len - RO3_DOOR_WIDTH) / 2
    _ro3_hinge = (_ro3_end[0] - ro3_gap * _ro3_len_u[0],
                  _ro3_end[1] - ro3_gap * _ro3_len_u[1])
    tips.append(("RO3_door_tip",
                 (_ro3_hinge[0] + RO3_DOOR_WIDTH * _ro3_swing[0],
                  _ro3_hinge[1] + RO3_DOOR_WIDTH * _ro3_swing[1])))

    # RO4 door tip (poly-based, rotation-safe)
    ro4 = [r for r in rough_openings if r.name == "RO4"][0]
    _ro4_dx = ro4.poly[3][0] - ro4.poly[0][0]
    _ro4_dy = ro4.poly[3][1] - ro4.poly[0][1]
    _ro4_len = math.sqrt(_ro4_dx**2 + _ro4_dy**2)
    _ro4_len_u = (_ro4_dx / _ro4_len, _ro4_dy / _ro4_len)
    _ro4_end = ((ro4.poly[2][0] + ro4.poly[3][0]) / 2,
                (ro4.poly[2][1] + ro4.poly[3][1]) / 2)
    _ro4_sw_dx = ro4.poly[3][0] - ro4.poly[2][0]
    _ro4_sw_dy = ro4.poly[3][1] - ro4.poly[2][1]
    _ro4_sw_len = math.sqrt(_ro4_sw_dx**2 + _ro4_sw_dy**2)
    _ro4_swing = (_ro4_sw_dx / _ro4_sw_len, _ro4_sw_dy / _ro4_sw_len)
    ro4_gap = (_ro4_len - RO4_DOOR_WIDTH) / 2
    _ro4_hinge = (_ro4_end[0] - ro4_gap * _ro4_len_u[0],
                  _ro4_end[1] - ro4_gap * _ro4_len_u[1])
    tips.append(("RO4_door_tip",
                 (_ro4_hinge[0] + RO4_DOOR_WIDTH * _ro4_swing[0],
                  _ro4_hinge[1] + RO4_DOOR_WIDTH * _ro4_swing[1])))

    # RO5 door tip (poly-based, rotation-safe)
    ro5 = [r for r in rough_openings if r.name == "RO5"][0]
    # Along direction: poly[0]→poly[1] (SW→SE)
    _ro5_dx = ro5.poly[1][0] - ro5.poly[0][0]
    _ro5_dy = ro5.poly[1][1] - ro5.poly[0][1]
    _ro5_len = math.sqrt(_ro5_dx**2 + _ro5_dy**2)
    _ro5_al = (_ro5_dx / _ro5_len, _ro5_dy / _ro5_len)
    # End edge midpoint (east): midpoint of poly[1] (SE) and poly[2] (NE)
    _ro5_end = ((ro5.poly[1][0] + ro5.poly[2][0]) / 2,
                (ro5.poly[1][1] + ro5.poly[2][1]) / 2)
    # Door swing direction (toward north face): SE→NE = poly[1]→poly[2]
    _ro5_sw_dx = ro5.poly[2][0] - ro5.poly[1][0]
    _ro5_sw_dy = ro5.poly[2][1] - ro5.poly[1][1]
    _ro5_sw_len = math.sqrt(_ro5_sw_dx**2 + _ro5_sw_dy**2)
    _ro5_swing = (_ro5_sw_dx / _ro5_sw_len, _ro5_sw_dy / _ro5_sw_len)
    ro5_gap = (_ro5_len - RO5_DOOR_WIDTH) / 2
    _ro5_hinge = (_ro5_end[0] - ro5_gap * _ro5_al[0],
                  _ro5_end[1] - ro5_gap * _ro5_al[1])
    tips.append(("RO5_door_tip",
                 (_ro5_hinge[0] + RO5_DOOR_WIDTH * _ro5_swing[0],
                  _ro5_hinge[1] + RO5_DOOR_WIDTH * _ro5_swing[1])))

    # RO6 double door tips (south and north leaves)
    ro6 = [r for r in rough_openings if r.name == "RO6"][0]
    _ro6_gap = (IW11_RO_WIDTH - 2 * RO6_DOOR_WIDTH) / 2
    _ro6_s_ctr = ((ro6.poly[0][0] + ro6.poly[1][0]) / 2,
                  (ro6.poly[0][1] + ro6.poly[1][1]) / 2)
    h_s6 = (_ro6_s_ctr[0] + _ro6_gap * _i11_an[0],
            _ro6_s_ctr[1] + _ro6_gap * _i11_an[1])
    tips.append(("RO6_door_tip_S",
                 (h_s6[0] + RO6_DOOR_WIDTH * _i11_at[0],
                  h_s6[1] + RO6_DOOR_WIDTH * _i11_at[1])))
    _ro6_n_ctr = ((ro6.poly[3][0] + ro6.poly[2][0]) / 2,
                  (ro6.poly[3][1] + ro6.poly[2][1]) / 2)
    h_n6 = (_ro6_n_ctr[0] - _ro6_gap * _i11_an[0],
            _ro6_n_ctr[1] - _ro6_gap * _i11_an[1])
    tips.append(("RO6_door_tip_N",
                 (h_n6[0] + RO6_DOOR_WIDTH * _i11_at[0],
                  h_n6[1] + RO6_DOOR_WIDTH * _i11_at[1])))

    # RO7 double door tips (south and north leaves)
    ro7 = [r for r in rough_openings if r.name == "RO7"][0]
    _i9_sw, _i9_se, _i9_ne, _i9_nw = layout.iw9.poly
    _i9_dx_n = _i9_ne[0] - _i9_se[0]; _i9_dy_n = _i9_ne[1] - _i9_se[1]
    _i9_ln = math.sqrt(_i9_dx_n**2 + _i9_dy_n**2)
    _i9_an = (_i9_dx_n / _i9_ln, _i9_dy_n / _i9_ln)
    _i9_dx_t = _i9_sw[0] - _i9_se[0]; _i9_dy_t = _i9_sw[1] - _i9_se[1]
    _i9_lt = math.sqrt(_i9_dx_t**2 + _i9_dy_t**2)
    _i9_at = (_i9_dx_t / _i9_lt, _i9_dy_t / _i9_lt)
    _ro7_gap = (IW9_RO_WIDTH - 2 * RO7_DOOR_WIDTH) / 2
    _ro7_s_ctr = ((ro7.poly[0][0] + ro7.poly[1][0]) / 2,
                  (ro7.poly[0][1] + ro7.poly[1][1]) / 2)
    h_s7 = (_ro7_s_ctr[0] + _ro7_gap * _i9_an[0],
            _ro7_s_ctr[1] + _ro7_gap * _i9_an[1])
    tips.append(("RO7_door_tip_S",
                 (h_s7[0] - RO7_DOOR_WIDTH * _i9_at[0],
                  h_s7[1] - RO7_DOOR_WIDTH * _i9_at[1])))
    _ro7_n_ctr = ((ro7.poly[3][0] + ro7.poly[2][0]) / 2,
                  (ro7.poly[3][1] + ro7.poly[2][1]) / 2)
    h_n7 = (_ro7_n_ctr[0] - _ro7_gap * _i9_an[0],
            _ro7_n_ctr[1] - _ro7_gap * _i9_an[1])
    tips.append(("RO7_door_tip_N",
                 (h_n7[0] - RO7_DOOR_WIDTH * _i9_at[0],
                  h_n7[1] - RO7_DOOR_WIDTH * _i9_at[1])))

    return tips


def _collect_f8f9_shell_points(pts):
    """Collect special points of the W8-W9 interior shell polylines.

    For both the W-series (WALL_OUTER inset) and G-series (SHELL+GAP inset)
    polylines: start, arc tangent entry, arc center, arc tangent exit, end.
    Arc center computed via traversal directions (rotation-safe).
    """
    F8, F9, F10, C7 = pts["F8"], pts["F9"], pts["F10"], pts["C7"]

    # CW traversal direction at F8 (tangent at exit of C7 arc)
    _r8x, _r8y = F8[0] - C7[0], F8[1] - C7[1]
    _r8_len = math.sqrt(_r8x**2 + _r8y**2)
    _dir_f8 = (_r8y / _r8_len, -_r8x / _r8_len)

    # CW traversal direction at F9 (F9→F10 line direction)
    _d9x, _d9y = F10[0] - F9[0], F10[1] - F9[1]
    _d9_len = math.sqrt(_d9x**2 + _d9y**2)
    _dir_f9 = (_d9x / _d9_len, _d9y / _d9_len)

    # Inset direction (right of CW direction = toward interior)
    _ins_f8 = (_dir_f8[1], -_dir_f8[0])
    _ins_f9 = (_dir_f9[1], -_dir_f9[0])

    # Left of CW direction (toward arc center for CCW turn)
    _left_f8 = (-_dir_f8[1], _dir_f8[0])
    _left_f9 = (-_dir_f9[1], _dir_f9[0])

    result = []

    for prefix, inset, R_turn in [
        ("w_f8f9", WALL_OUTER, F8F9_INNER_TURN_R),
        ("g_f8f9", SHELL_THICKNESS + AIR_GAP, OPENING_INSIDE_RADIUS),
    ]:
        poly = f8f9_corner_polyline(pts, inset, R_turn)
        start = (F8[0] + inset * _ins_f8[0], F8[1] + inset * _ins_f8[1])
        end   = (F9[0] + inset * _ins_f9[0], F9[1] + inset * _ins_f9[1])

        # Arc center via line-line intersection (matches f8f9_corner_polyline)
        _P1 = (start[0] + R_turn * _left_f8[0], start[1] + R_turn * _left_f8[1])
        _P2 = (end[0]   + R_turn * _left_f9[0], end[1]   + R_turn * _left_f9[1])
        _dp = (_P2[0] - _P1[0], _P2[1] - _P1[1])
        _cross = _dir_f8[0] * _dir_f9[1] - _dir_f8[1] * _dir_f9[0]
        _t = (_dp[0] * _dir_f9[1] - _dp[1] * _dir_f9[0]) / _cross
        arc_cx = _P1[0] + _t * _dir_f8[0]
        arc_cy = _P1[1] + _t * _dir_f8[1]

        result.append((f"{prefix}_start", poly[0]))
        result.append((f"{prefix}_arc_entry", poly[1]))
        result.append((f"{prefix}_arc_center", (arc_cx, arc_cy)))
        result.append((f"{prefix}_arc_exit", poly[-2]))
        result.append((f"{prefix}_end", poly[-1]))

    return result


# Expected distance-squared values: (name, d2_F1, d2_F6, d2_F12, d2_F15)
EXPECTED = [
("F1", 0.000000, 749.361111, 1610.150678, 1269.687726),
    ("F2", 1.388889, 707.694444, 1620.360460, 1320.146256),
    ("F5", 625.694444, 10.888889, 1007.143645, 1667.619769),
    ("F6", 749.361111, 0.000000, 867.141976, 1600.434218),
    ("F7", 783.111111, 20.250000, 622.823956, 1317.660118),
    ("F8", 694.444444, 52.138889, 509.458790, 1091.598455),
    ("F9", 688.944444, 55.250000, 502.194612, 1076.284526),
    ("F10", 1211.301029, 529.952916, 42.998506, 481.342250),
    ("F11", 1300.151063, 573.665322, 30.334251, 490.728544),
    ("F12", 1610.150678, 867.141976, 0.000000, 413.327009),
    ("F13", 1330.403050, 1429.151431, 261.762749, 18.057535),
    ("F14", 1323.877839, 1457.848697, 283.191787, 12.956791),
    ("F15", 1269.687726, 1600.434218, 413.327009, 0.000000),
    ("F16", 1180.972108, 1623.652437, 488.898341, 4.984551),
    ("F17", 968.827936, 1544.210817, 582.178645, 34.706201),
    ("F18", 784.239564, 1434.839254, 650.186178, 81.143273),
    ("F19", 708.316722, 1371.520371, 667.997262, 104.689332),
    ("F20", 601.625705, 1277.402638, 696.116618, 146.041882),
    ("F11a", 1511.225535, 683.436565, 13.578291, 523.461048),
    ("F11b", 1567.510775, 736.721804, 8.070637, 509.407599),
    ("FC", 494.361111, 452.722222, 321.531742, 366.754560),
    ("P3", 0.694444, 752.555556, 1662.339201, 1328.997744),
    ("POB", 1417.872460, 312.896625, 212.551497, 1024.389218),
    ("P2", 841.694444, 8.222222, 1018.312357, 1837.799293),
    ("P4", 961.650783, 1617.980410, 655.573271, 50.100969),
    ("P5", 1820.080898, 2045.939326, 460.305014, 49.562685),
    ("T1", 1498.618268, 1402.474612, 176.213533, 61.280634),
    ("TC1", 2314.261138, 1824.647647, 190.539879, 283.027793),
    ("T2", 1316.080277, 430.001990, 85.354345, 696.100327),
    ("TC2", 2366.883865, 988.968283, 134.512276, 1004.534276),
    ("PA", 1856.242937, 1122.428701, 16.466506, 379.778367),
    ("T3", 291.655033, 989.782412, 846.512200, 360.104913),
    ("TC3", 412.655033, 1712.115745, 1530.798244, 607.111222),
    ("PX", 656.083155, 1533.759161, 891.812666, 186.211076),
    ("PiOB", 1375.815362, 301.233108, 209.316065, 1000.934661),
    ("Pi2", 817.617319, 4.945518, 984.228223, 1783.225762),
    ("Pi3", 0.361111, 723.388889, 1605.505510, 1287.767507),
    ("Pi4", 944.833994, 1579.167718, 634.383015, 47.857699),
    ("Pi5", 1757.496093, 1980.578004, 440.391554, 39.902039),
    ("Ti1", 1463.086125, 1386.615960, 180.747216, 55.443276),
    ("Ti2", 1280.548133, 414.143338, 89.888027, 690.262969),
    ("Ti3", 291.905033, 962.699079, 821.158289, 354.627354),
    ("PiX", 655.754857, 1497.522986, 853.057463, 172.696145),
    ("Ai2", 1754.548362, 1064.584740, 10.942580, 352.387187),
    ("IW1_SW", 180.694444, 224.111111, 864.144689, 988.759896),
    ("IW1_SE", 1294.301009, 1250.199316, 175.080358, 50.390893),
    ("IW1_NE", 1298.213029, 1227.179926, 161.667145, 57.998951),
    ("IW1_NW", 193.611111, 209.694444, 851.457444, 995.949003),
    ("IW2_SW", 193.611111, 209.694444, 851.457444, 995.949003),
    ("IW2_SE", 198.361111, 212.944444, 825.310997, 965.529659),
    ("IW2_NE", 736.111111, 12.694444, 671.755398, 1348.635553),
    ("IW2_NW", 731.361111, 9.444444, 697.901844, 1379.054897),
    ("IW3_SW", 69.888889, 757.805556, 1130.960242, 745.784128),
    ("IW3_SE", 75.555556, 762.472222, 1116.029278, 728.004565),
    ("IW3_NE", 128.888889, 451.361111, 827.977130, 704.970438),
    ("IW3_NW", 123.222222, 446.694444, 842.908094, 722.750001),
    ("IW4_SW", 702.731982, 1330.062989, 635.669261, 99.694913),
    ("IW4_SE", 720.508206, 1346.839213, 632.847854, 94.024908),
    ("IW4_NE", 768.902739, 1055.104174, 362.633338, 71.147037),
    ("IW4_NW", 751.126515, 1038.327949, 365.454744, 76.817043),
    ("IW5_SW", 666.368056, 802.118056, 295.977208, 146.053767),
    ("IW5_SE", 1288.775051, 1392.517017, 259.094541, 17.993063),
    ("IW5_NE", 1293.795884, 1383.871184, 251.313419, 20.150117),
    ("IW5_NW", 671.388889, 793.472222, 288.196085, 148.210820),
    ("IW6_SW", 444.534722, 41.840278, 985.352631, 1484.433471),
    ("IW6_SE", 464.756944, 48.062500, 717.208018, 1176.408478),
    ("IW6_NE", 468.277778, 47.027778, 716.461533, 1178.974718),
    ("IW6_NW", 448.055556, 40.805556, 984.606146, 1486.999711),
    ("IW7_SW", 124.111111, 464.805556, 840.268626, 704.011034),
    ("IW7_SE", 173.694444, 506.888889, 735.369726, 577.747645),
    ("IW7_NE", 178.472222, 493.444444, 723.078230, 578.707050),
    ("IW7_NW", 128.888889, 451.361111, 827.977130, 704.970438),
    ("IW8_SW", 160.472222, 217.888889, 1132.289302, 1296.784889),
    ("IW8_SE", 180.694444, 224.111111, 864.144689, 988.759896),
    ("IW8_NE", 193.611111, 209.694444, 851.457444, 995.949003),
    ("IW8_NW", 173.388889, 203.472222, 1119.602057, 1303.973996),
    ("IW9_SW", 125.138889, 804.555556, 1011.130378, 601.741176),
    ("IW9_SE", 132.694444, 811.111111, 998.088302, 585.850502),
    ("IW9_NE", 186.027778, 500.000000, 710.036154, 562.816375),
    ("IW9_NW", 178.472222, 493.444444, 723.078230, 578.707050),
    ("IW11_SW", 552.694444, 1195.111111, 676.573583, 161.786236),
    ("IW11_SE", 568.472222, 1209.888889, 671.753730, 154.117784),
    ("IW11_NE", 728.472222, 713.888889, 217.259863, 176.656356),
    ("IW11_NW", 712.694444, 699.111111, 222.079716, 184.324808),
    ("IW12_SW", 612.472222, 925.888889, 408.506796, 129.387070),
    ("IW12_SE", 746.570959, 1051.994616, 377.968463, 76.079860),
    ("IW12_NE", 751.126515, 1038.327949, 365.454744, 76.817043),
    ("IW12_NW", 617.027778, 912.222222, 395.993078, 130.124252),
    ("IW16_SW", 123.222222, 446.694444, 842.908094, 722.750001),
    ("IW16_SE", 128.888889, 451.361111, 827.977130, 704.970438),
    ("IW16_NE", 235.555556, 266.472222, 661.535411, 750.543137),
    ("IW16_NW", 229.888889, 261.805556, 676.466375, 768.322700),
    ("dryer_SW", 1.111111, 694.805556, 1539.456336, 1235.897488),
    ("dryer_SE", 11.562500, 696.506944, 1369.678453, 1041.194368),
    ("dryer_NE", 22.812500, 571.090278, 1252.908898, 1023.806570),
    ("dryer_NW", 12.361111, 569.388889, 1422.686781, 1218.509691),
    ("washer_SW", 12.951389, 565.423611, 1419.009740, 1218.145375),
    ("washer_SE", 23.402778, 567.125000, 1249.231857, 1023.442255),
    ("washer_NE", 47.569444, 454.625000, 1145.378968, 1018.971124),
    ("washer_NW", 37.118056, 452.923611, 1315.156851, 1213.674244),
    ("ctr_poly_0", 83.506944, 449.673611, 966.066559, 837.210405),
    ("ctr_poly_1", 113.888889, 473.805556, 867.713309, 721.053414),
    ("ctr_poly_2", 69.888889, 757.805556, 1130.960242, 745.784128),
    ("ctr_poly_3", 39.506944, 733.673611, 1229.313492, 861.941119),
    ("bed_SW", 208.534722, 869.090278, 885.164671, 454.627084),
    ("bed_SE", 431.256944, 1072.812500, 716.531902, 231.870944),
    ("bed_NE", 505.673611, 719.006944, 389.820628, 216.555845),
    ("bed_NW", 282.951389, 515.284722, 558.453397, 439.311985),
    ("dresser_SW", 541.250000, 627.777778, 321.360322, 242.988295),
    ("dresser_SE", 665.444444, 743.472222, 270.474902, 167.889788),
    ("dresser_NE", 702.784722, 694.256944, 226.736128, 187.092794),
    ("dresser_NW", 578.590278, 578.562500, 277.621548, 262.191301),
    ("RO1_SW", 375.555556, 388.472222, 430.778051, 468.511004),
    ("RO1_SE", 478.472222, 481.888889, 338.017223, 348.688489),
    ("RO1_NE", 491.388889, 467.472222, 325.329978, 355.877596),
    ("RO1_NW", 388.472222, 374.055556, 418.090807, 475.700111),
    ("RO2_SW", 620.590278, 902.118056, 386.753622, 130.822973),
    ("RO2_SE", 604.812500, 887.340278, 391.573476, 138.491425),
    ("RO2_NE", 660.756944, 770.173611, 285.359816, 158.161326),
    ("RO2_NW", 676.534722, 784.951389, 280.539963, 150.492874),
    ("RO3_SW", 140.284722, 404.534722, 804.495453, 727.402789),
    ("RO3_SE", 145.951389, 409.201389, 789.564489, 709.623226),
    ("RO3_NE", 209.284722, 299.423611, 690.739718, 736.682016),
    ("RO3_NW", 203.618056, 294.756944, 705.670683, 754.461579),
    ("RO4_SW", 314.944444, 112.361111, 767.959489, 1071.461861),
    ("RO4_SE", 319.694444, 115.611111, 741.813042, 1041.042516),
    ("RO4_NE", 438.444444, 61.250000, 698.404938, 1123.517973),
    ("RO4_NW", 433.694444, 58.000000, 724.551385, 1153.937317),
    ("RO5_SW", 445.680556, 39.236111, 909.257348, 1397.655943),
    ("RO5_SE", 462.569444, 46.625000, 730.468742, 1191.805650),
    ("RO5_NE", 466.090278, 45.590278, 729.722256, 1194.371890),
    ("RO5_NW", 449.201389, 38.201389, 908.510863, 1400.222183),
    ("RO6_SW", 569.201389, 1187.840278, 651.146304, 150.073984),
    ("RO6_SE", 553.423611, 1173.062500, 655.966157, 157.742436),
    ("RO6_NE", 591.312500, 928.506944, 429.281298, 136.446544),
    ("RO6_NW", 607.090278, 943.284722, 424.461445, 128.778092),
    ("RO7_SW", 133.812500, 780.340278, 969.335128, 580.286405),
    ("RO7_SE", 126.256944, 773.784722, 982.377204, 596.177079),
    ("RO7_NE", 165.868056, 530.951389, 757.414567, 576.603409),
    ("RO7_NW", 173.423611, 537.506944, 744.372491, 560.712735),
    ("O1_outer_start", 85.487847, 333.960069, 1275.593907, 1308.308593),
    ("O1_outer_end", 117.154514, 279.071181, 1226.181522, 1321.837988),
    ("O1_inner_end", 116.487847, 276.404514, 1184.319593, 1274.278862),
    ("O1_inner_start", 84.821181, 331.293403, 1233.731978, 1260.749467),
    ("O2_outer_start", 323.196181, 93.335069, 1065.233796, 1445.784635),
    ("O2_outer_end", 382.571181, 66.154514, 1043.529744, 1487.022363),
    ("O2_inner_end", 381.904514, 63.487847, 1001.667816, 1439.463237),
    ("O2_inner_start", 322.529514, 90.668403, 1023.371867, 1398.225509),
    ("O3_outer_start", 484.694444, 33.888889, 1019.767112, 1560.985126),
    ("O3_outer_end", 609.138889, 12.555556, 1007.657364, 1654.882587),
    ("O3_inner_end", 608.472222, 9.888889, 965.795435, 1607.323461),
    ("O3_inner_start", 484.027778, 31.222222, 977.905183, 1513.426000),
    ("O4_outer_start", 722.501736, 3.960069, 758.559474, 1449.326547),
    ("O4_outer_end", 728.126736, 7.335069, 717.839804, 1402.197530),
    ("O4_inner_end", 764.126736, 6.890625, 719.701256, 1430.560784),
    ("O4_inner_start", 758.501736, 3.515625, 760.420926, 1477.689801),
    ("O5_outer_start", 671.138889, 71.388889, 467.033527, 1007.500143),
    ("O5_outer_end", 809.027778, 192.277778, 254.762687, 746.803128),
    ("O5_inner_end", 841.694444, 188.500000, 253.290806, 771.833049),
    ("O5_inner_start", 703.805556, 67.611111, 465.561645, 1032.530063),
    ("O6_outer_start", 992.790838, 360.387169, 116.097252, 563.546675),
    ("O6_outer_end", 1154.499806, 511.096137, 51.232278, 467.347117),
    ("O6_inner_end", 1187.166473, 507.318360, 49.760397, 492.377038),
    ("O6_inner_start", 1025.457505, 356.609391, 114.625370, 588.576595),
    ("O7_outer_start", 1547.211118, 908.257391, 4.000000, 336.107054),
    ("O7_outer_end", 1406.392436, 1079.603637, 64.000000, 152.447190),
    ("O7_inner_end", 1358.209606, 1042.785176, 64.444444, 151.535533),
    ("O7_inner_start", 1499.028288, 871.438930, 4.444444, 335.195397),
    ("O8_outer_start", 1283.994260, 1552.627378, 369.450179, 1.290996),
    ("O8_outer_end", 1313.979027, 1477.834368, 301.286620, 9.320119),
    ("O8_inner_end", 1267.531012, 1433.386353, 296.033867, 9.764563),
    ("O8_inner_start", 1237.546245, 1508.179364, 364.197426, 1.735441),
    ("O9_outer_start", 529.000000, 1209.361111, 717.914133, 180.897882),
    ("O9_outer_end", 444.506944, 1130.618056, 751.857457, 231.220646),
    ("O9_inner_end", 444.951389, 1094.618056, 718.163353, 224.028345),
    ("O9_inner_start", 529.444444, 1173.361111, 684.220030, 173.705580),
    ("O10_outer_start", 193.673611, 901.284722, 943.873192, 484.481250),
    ("O10_outer_end", 144.000000, 857.361111, 1012.635959, 569.623459),
    ("O10_inner_end", 144.444444, 821.361111, 978.941856, 562.431158),
    ("O10_inner_start", 194.118056, 865.284722, 910.179088, 477.288948),
    ("O11_outer_start", 36.000000, 767.361111, 1275.393319, 883.655592),
    ("O11_outer_end", 19.506944, 755.618056, 1356.739011, 978.532128),
    ("O11_inner_end", 19.951389, 719.618056, 1323.044907, 971.339826),
    ("O11_inner_start", 36.444444, 731.361111, 1241.699215, 876.463291),
    ("O3_door_tip", 608.340278, 7.812500, 836.046207, 1455.745617),
    ("O6_door_tip", 1007.498213, 537.455656, 71.538197, 362.481625),
    ("RO1_door_tip", 413.402778, 567.402778, 418.918697, 320.982179),
    ("RO2_door_tip", 817.888889, 922.361111, 251.081584, 93.508682),
    ("RO3_door_tip", 162.500000, 266.694444, 843.899810, 913.591068),
    ("RO4_door_tip", 413.125000, 50.236111, 878.179716, 1327.753860),
    ("RO5_door_tip", 599.376736, 17.404514, 716.668315, 1299.631846),
    ("RO6_door_tip_S", 449.388889, 1071.472222, 691.847534, 216.876728),
    ("RO6_door_tip_N", 486.055556, 834.805556, 472.475090, 196.267800),
    ("RO7_door_tip_S", 193.138889, 828.111111, 880.210781, 474.534593),
    ("RO7_door_tip_N", 231.472222, 593.111111, 662.505003, 455.592332),
    ("w_f8f9_start", 683.777778, 43.472222, 539.987386, 1127.824247),
    ("w_f8f9_arc_entry", 652.486641, 46.832004, 541.154750, 1103.794032),
    ("w_f8f9_arc_center", 655.585049, 49.331987, 531.927017, 1092.861625),
    ("w_f8f9_arc_exit", 645.903964, 50.555538, 532.460615, 1085.465551),
    ("w_f8f9_end", 656.277778, 59.027778, 503.666493, 1051.254606),
    ("g_f8f9_start", 686.361111, 45.555556, 532.271904, 1118.684466),
    ("g_f8f9_arc_entry", 655.069974, 48.915337, 533.439268, 1094.654251),
    ("g_f8f9_arc_center", 655.585049, 49.331987, 531.927017, 1092.861625),
    ("g_f8f9_arc_exit", 653.987297, 49.527761, 532.009312, 1091.639698),
    ("g_f8f9_end", 664.361111, 58.000000, 503.215190, 1057.428753),
    ("wh_center", 686.820226, 25.244557, 611.579575, 1230.974734),
    ("sink_NW", 708.722222, 103.472222, 391.329889, 916.129213),
    ("sink_SE", 713.868056, 206.701389, 267.168017, 680.415201),
    ("stove_NW", 640.069444, 57.236111, 515.708339, 1055.494724),
    ("stove_SE", 585.006944, 117.673611, 421.381206, 850.359520),
    ("dw_NW", 814.055556, 196.805556, 249.491649, 740.107791),
    ("dw_SE", 786.590278, 285.340278, 193.060554, 574.292924),
    ("fridge_NW", 294.340278, 135.062500, 743.683837, 1007.649610),
    ("fridge_SE", 243.674045, 235.653212, 673.730648, 798.428136),
    ("ice_NW", 896.006944, 283.923611, 169.165130, 625.630008),
    ("ice_SE", 890.035625, 345.505069, 139.472925, 530.991696),
    ("nctr_NW", 609.027778, 22.277778, 672.691620, 1250.190017),
    ("nctr_SE", 533.444444, 74.361111, 539.749162, 992.228414),
    ("wctr_NW", 282.830295, 205.726128, 633.428897, 805.680815),
    ("wctr_SE", 348.371962, 338.267795, 466.817831, 536.711716),
    ("loveseat_NW", 835.440201, 492.448567, 126.212204, 344.357581),
    ("et_center", 938.545808, 832.669105, 144.995391, 134.982558),
    ("loveseat2_NW", 1056.886834, 833.443203, 94.999473, 159.892815),
    ("loveseat2_SE", 1304.978035, 1224.728848, 157.735317, 60.393548),
    ("chair_center", 1377.073949, 704.734596, 10.295107, 404.455617),
    ("ottoman_center", 1159.205859, 645.605352, 36.969112, 333.129274),
    ("desk_SW", 946.292855, 1493.504426, 556.555947, 33.807080),
    ("dim01_A", 495.362473, 471.112617, 322.440655, 352.039151),
    ("dim01_B", 906.029140, 280.445950, 169.821277, 636.699508),
    ("dim02_A", 341.484588, 103.324949, 731.682941, 1055.692145),
    ("dim02_B", 1359.773942, 1039.714334, 63.026927, 153.735996),
    ("dim03_A", 677.746918, 987.167080, 391.462957, 100.958792),
    ("dim03_B", 633.746918, 1271.167080, 654.709890, 125.689506),
    ("dim04_A", 147.340278, 484.284722, 786.256676, 639.316839),
    ("dim04_B", 98.784722, 781.951389, 1062.017328, 663.310371),
    ("dim05_A", 16.027778, 547.222222, 1431.645983, 1251.618142),
    ("dim05_B", 85.444444, 591.138889, 975.823057, 723.155954),
    ("dim06_A", 781.088192, 1022.347445, 332.720170, 73.915087),
    ("dim06_B", 1251.620226, 1469.864455, 329.197244, 4.831599),
    ("dim07_A", 698.368056, 752.118056, 251.165474, 160.871088),
    ("dim07_B", 1320.775051, 1342.517017, 214.282808, 32.810385),
    ("dim08_A", 100.027778, 303.222222, 1208.399050, 1266.887428),
    ("dim08_B", 169.444444, 347.138889, 752.576124, 738.425240),
    ("dim10_A", 625.027778, 8.222222, 965.281717, 1620.060643),
    ("dim10_B", 645.250000, 14.444444, 697.137104, 1312.035651),
    ("dim11_A", 944.611182, 917.518917, 175.744458, 99.508026),
    ("dim11_B", 785.036204, 1399.414354, 617.053071, 74.349876),
    ("dim12a_A", 454.277778, 39.027778, 826.047320, 1305.652096),
    ("dim12a_B", 717.361111, 1.444444, 807.487631, 1505.732275),
    ("dim12b_A", 179.611111, 201.694444, 961.043231, 1122.626381),
    ("dim12b_B", 450.756944, 40.062500, 826.793805, 1303.085856),
    ("dim13_A", 12.557270, 757.765892, 1827.223762, 1527.683722),
    ("dim13_B", 759.595556, 25.334444, 1188.400650, 1964.706729),
    ("dim14_A", 854.644427, 33.511098, 568.883635, 1293.201732),
    ("dim14_B", 1501.892058, 623.041972, 27.302806, 587.209035),
    ("dim15_A", 9.014621, 918.560179, 1818.366885, 1370.359886),
    ("dim15_B", 1245.202964, 2046.740489, 829.845604, 74.167079),
    ("dim16_A", 486.279514, 1133.071181, 700.273288, 197.948560),
    ("dim16_B", 646.279514, 637.071181, 245.779422, 220.487132),
    ("dim17_A", 168.362847, 842.404514, 943.642069, 518.941650),
    ("dim17_B", 328.362847, 346.404514, 489.148202, 541.480222),
    ("dim18_A", 148.250000, 644.444444, 842.951117, 563.222328),
    ("dim18_B", 568.250000, 1028.444444, 521.436398, 139.158061),
    ("dim19_A", 27.571181, 724.862847, 1281.745325, 923.274822),
    ("dim19_B", 187.571181, 228.862847, 827.251459, 945.813394),
]

REF_NAMES = ["F1", "F6", "F12", "F15"]


class TestDistanceSquaredRegression:
    """Verify that every interior point stays at the expected distance
    from F1, F6, F12, F15 (within +/- 0.000001 ft^2)."""

    @pytest.fixture(scope="class")
    def all_pts(self, pts_with_outline, layout, outline_geo):
        pts, _ = pts_with_outline
        return _collect_all_points(pts, layout, outline_geo.radii)

    @pytest.fixture(scope="class")
    def ref_pts(self, pts_with_outline):
        pts, _ = pts_with_outline
        return {n: pts[n] for n in REF_NAMES}

    def test_point_count(self, all_pts):
        """Sanity check: total collected points matches expected."""
        assert len(all_pts) == len(EXPECTED), (
            f"collected {len(all_pts)} points, expected {len(EXPECTED)}")

    def test_point_names_match(self, all_pts):
        """Verify collected point names match expected names in order."""
        actual_names = [name for name, _ in all_pts]
        expected_names = [name for name, *_ in EXPECTED]
        assert actual_names == expected_names

    _F_NAMES = [f"F{i}" for i in range(1, 21) if i not in (3, 4)] + ["F11a", "F11b"]

    def test_outline_closure(self, outline_geo):
        """Chain walk must close: endpoint F2 matches starting F2."""
        fp_pts = outline_geo.fp_pts
        assert abs(fp_pts["F2"][0] - F2_E) < 1e-9, "closure failed (E)"
        assert abs(fp_pts["F2"][1] - F2_N) < 1e-9, "closure failed (N)"

    def test_fc_is_origin(self, outline_geo):
        """FC (building center) must be at origin."""
        assert outline_geo.fp_pts["FC"] == (0.0, 0.0)

    @pytest.mark.parametrize("idx", range(len(EXPECTED)),
                             ids=[e[0] for e in EXPECTED])
    def test_d2_regression(self, all_pts, ref_pts, idx):
        exp_name, exp_f1, exp_f6, exp_f12, exp_f15 = EXPECTED[idx]
        act_name, pt = all_pts[idx]
        assert act_name == exp_name, f"name mismatch at {idx}"
        expected = [exp_f1, exp_f6, exp_f12, exp_f15]
        for i, rn in enumerate(REF_NAMES):
            actual = _dist_sq(pt, ref_pts[rn])
            assert abs(actual - expected[i]) < TOL, (
                f"{act_name} -> {rn}: d^2 = {actual:.12f}, "
                f"expected {expected[i]:.12f}, delta {actual - expected[i]:.2e}")
