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
        ("IW2o", layout.iw2o), ("IW2s", layout.iw2s),
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
    ("F7", 791.555556, 26.694444, 590.073138, 1279.212103),
    ("F8", 706.000000, 61.694444, 479.819084, 1056.261551),
    ("F9", 700.722222, 65.027778, 472.777128, 1041.169845),
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
    ("IW1_SW", 187.138889, 228.555556, 829.393871, 948.311881),
    ("IW1_SE", 1294.301009, 1250.199316, 175.080358, 50.390893),
    ("IW1_NE", 1298.213029, 1227.179926, 161.667145, 57.998951),
    ("IW1_NW", 200.055556, 214.138889, 816.706626, 955.500989),
    ("IW2_SW", 213.472222, 224.055556, 758.031584, 886.855852),
    ("IW2_SE", 220.055556, 229.138889, 733.718471, 858.269840),
    ("IW2_NE", 324.472222, 142.222222, 658.907760, 922.593591),
    ("IW2_NW", 317.888889, 137.138889, 683.220873, 951.179602),
    ("IW2o_SW", 326.673889, 140.627222, 658.142808, 924.696786),
    ("IW2o_SE", 315.687222, 138.733889, 683.985825, 949.076407),
    ("IW2o_NE", 451.020556, 58.900556, 686.579578, 1120.651068),
    ("IW2o_NW", 462.007222, 60.793889, 660.736561, 1096.271447),
    ("IW2s_SW", 453.805556, 57.888889, 686.397960, 1123.337596),
    ("IW2s_SE", 459.222222, 61.805556, 660.918180, 1093.584918),
    ("IW2s_NE", 743.222222, 17.805556, 637.671247, 1308.854204),
    ("IW2s_NW", 737.805556, 13.888889, 663.151027, 1338.606882),
    ("IW3_SW", 69.888889, 757.805556, 1130.960242, 745.784128),
    ("IW3_SE", 75.555556, 762.472222, 1116.029278, 728.004565),
    ("IW3_NE", 128.888889, 451.361111, 827.977130, 704.970438),
    ("IW3_NW", 123.222222, 446.694444, 842.908094, 722.750001),
    ("IW4_SW", 693.888889, 1327.805556, 642.688163, 103.687728),
    ("IW4_SE", 711.555556, 1344.472222, 639.757199, 97.908165),
    ("IW4_NE", 760.111111, 1046.805556, 363.996547, 73.914634),
    ("IW4_NW", 742.444444, 1030.138889, 366.927511, 79.694197),
    ("IW5_SW", 666.368056, 802.118056, 295.977208, 146.053767),
    ("IW5_SE", 1288.775051, 1392.517017, 259.094541, 17.993063),
    ("IW5_NE", 1293.795884, 1383.871184, 251.313419, 20.150117),
    ("IW5_NW", 671.388889, 793.472222, 288.196085, 148.210820),
    ("IW6_SW", 444.534722, 41.840278, 985.352631, 1484.433471),
    ("IW6_SE", 471.201389, 52.506944, 682.457201, 1135.960463),
    ("IW6_NE", 474.722222, 51.472222, 681.710715, 1138.526703),
    ("IW6_NW", 448.055556, 40.805556, 984.606146, 1486.999711),
    ("IW7_SW", 124.111111, 464.805556, 840.268626, 704.011034),
    ("IW7_SE", 173.694444, 506.888889, 735.369726, 577.747645),
    ("IW7_NE", 178.472222, 493.444444, 723.078230, 578.707050),
    ("IW7_NW", 128.888889, 451.361111, 827.977130, 704.970438),
    ("IW8_SW", 160.472222, 217.888889, 1132.289302, 1296.784889),
    ("IW8_SE", 187.138889, 228.555556, 829.393871, 948.311881),
    ("IW8_NE", 200.055556, 214.138889, 816.706626, 955.500989),
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
    ("IW12_SE", 737.888889, 1043.805556, 379.441230, 78.957014),
    ("IW12_NE", 742.444444, 1030.138889, 366.927511, 79.694197),
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
    ("RO1_SW", 432.694444, 440.111111, 374.629969, 396.695630),
    ("RO1_SE", 547.222222, 545.138889, 293.480252, 288.484226),
    ("RO1_NE", 560.138889, 530.722222, 280.793007, 295.673333),
    ("RO1_NW", 445.611111, 425.694444, 361.942725, 403.884737),
    ("RO2_SW", 620.590278, 902.118056, 386.753622, 130.822973),
    ("RO2_SE", 604.812500, 887.340278, 391.573476, 138.491425),
    ("RO2_NE", 660.756944, 770.173611, 285.359816, 158.161326),
    ("RO2_NW", 676.534722, 784.951389, 280.539963, 150.492874),
    ("RO3_SW", 140.284722, 404.534722, 804.495453, 727.402789),
    ("RO3_SE", 145.951389, 409.201389, 789.564489, 709.623226),
    ("RO3_NE", 209.284722, 299.423611, 690.739718, 736.682016),
    ("RO3_NW", 203.618056, 294.756944, 705.670683, 754.461579),
    ("RO4_SW", 341.080556, 129.213889, 656.620725, 943.452412),
    ("RO4_SE", 330.093889, 127.320556, 682.463742, 967.832033),
    ("RO4_NE", 432.947222, 66.647222, 684.434995, 1098.228775),
    ("RO4_NW", 443.933889, 68.540556, 658.591977, 1073.849154),
    ("RO5_SW", 447.569444, 39.125000, 869.950975, 1352.652372),
    ("RO5_SE", 468.680556, 50.736111, 695.384591, 1151.024302),
    ("RO5_NE", 472.201389, 49.701389, 694.638105, 1153.590542),
    ("RO5_NW", 451.090278, 38.090278, 869.204490, 1355.218613),
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
    ("O2_outer_start", 252.930604, 136.578752, 1101.295573, 1399.304340),
    ("O2_outer_end", 305.730372, 102.822965, 1073.016289, 1433.966837),
    ("O2_inner_end", 305.063706, 100.156298, 1031.154360, 1386.407711),
    ("O2_inner_start", 252.263937, 133.912085, 1059.433644, 1351.745214),
    ("O3_outer_start", 484.694444, 33.888889, 1019.767112, 1560.985126),
    ("O3_outer_end", 609.138889, 12.555556, 1007.657364, 1654.882587),
    ("O3_inner_end", 608.472222, 9.888889, 965.795435, 1607.323461),
    ("O3_inner_start", 484.027778, 31.222222, 977.905183, 1513.426000),
    ("O4_outer_start", 724.862847, 5.321181, 740.322954, 1428.241428),
    ("O4_outer_end", 730.987847, 9.196181, 700.103284, 1381.612412),
    ("O4_inner_end", 766.987847, 8.751736, 701.964736, 1409.975666),
    ("O4_inner_start", 760.862847, 4.876736, 742.184406, 1456.604682),
    ("O5_outer_start", 684.027778, 82.277778, 438.727153, 973.496572),
    ("O5_outer_end", 829.472222, 210.722222, 234.011869, 720.355113),
    ("O5_inner_end", 862.138889, 206.944444, 232.539988, 745.385034),
    ("O5_inner_start", 716.694444, 78.500000, 437.255272, 998.526493),
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
    ("RO1_door_tip", 481.847222, 630.347222, 374.076171, 260.472361),
    ("RO2_door_tip", 817.888889, 922.361111, 251.081584, 93.508682),
    ("RO3_door_tip", 162.500000, 266.694444, 843.899810, 913.591068),
    ("RO4_door_tip", 510.334444, 89.231111, 525.084063, 945.010303),
    ("RO5_door_tip", 605.376736, 21.404514, 681.473053, 1258.739386),
    ("RO6_door_tip_S", 449.388889, 1071.472222, 691.847534, 216.876728),
    ("RO6_door_tip_N", 486.055556, 834.805556, 472.475090, 196.267800),
    ("RO7_door_tip_S", 193.138889, 828.111111, 880.210781, 474.534593),
    ("RO7_door_tip_N", 231.472222, 593.111111, 662.505003, 455.592332),
    ("w_f8f9_start", 694.444444, 52.138889, 509.458790, 1091.598455),
    ("w_f8f9_arc_entry", 663.153307, 55.498670, 510.626155, 1067.568239),
    ("w_f8f9_arc_center", 666.517682, 58.264620, 501.664388, 1056.901800),
    ("w_f8f9_arc_exit", 656.836597, 59.488172, 502.197987, 1049.505725),
    ("w_f8f9_end", 668.055556, 68.805556, 474.249009, 1016.139924),
    ("g_f8f9_start", 697.250000, 54.444444, 501.965530, 1082.680895),
    ("g_f8f9_arc_entry", 665.958863, 57.804226, 503.132895, 1058.650680),
    ("g_f8f9_arc_center", 666.517682, 58.264620, 501.664388, 1056.901800),
    ("g_f8f9_arc_exit", 664.919930, 58.460394, 501.746683, 1055.679872),
    ("g_f8f9_end", 676.138889, 67.777778, 473.797705, 1022.314071),
    ("wh_center", 695.486893, 31.911223, 579.050980, 1192.748941),
    ("sink_NW", 724.055556, 116.805556, 365.467960, 884.570087),
    ("sink_SE", 734.201389, 225.034722, 246.306088, 653.856075),
    ("stove_NW", 651.513889, 66.680556, 485.957521, 1020.046710),
    ("stove_SE", 599.784722, 130.451389, 394.963721, 818.244838),
    ("dw_NW", 834.722222, 215.472222, 228.963054, 713.881998),
    ("dw_SE", 810.368056, 307.118056, 175.643070, 551.178242),
    ("fridge_NW", 316.951389, 152.173611, 653.007977, 901.306458),
    ("fridge_SE", 276.292101, 262.771267, 593.061733, 702.091929),
    ("ice_NW", 920.451389, 306.368056, 152.414312, 603.181993),
    ("ice_SE", 916.446736, 369.916181, 124.688774, 510.510348),
    ("nctr_NW", 616.138889, 27.388889, 638.607469, 1210.408669),
    ("nctr_SE", 544.555556, 83.472222, 509.665011, 956.447066),
    ("wctr_NW", 316.365017, 233.760851, 553.676648, 710.261274),
    ("wctr_SE", 400.240017, 384.635851, 405.398915, 459.625509),
    ("loveseat_NW", 828.244635, 485.746011, 129.171476, 348.721240),
    ("et_center", 929.642675, 824.258982, 146.247095, 137.638650),
    ("loveseat2_NW", 1047.659342, 824.708720, 95.926818, 162.224547),
    ("loveseat2_SE", 1293.970229, 1214.214051, 156.882348, 60.944966),
    ("chair_center", 1377.073949, 704.734596, 10.295107, 404.455617),
    ("ottoman_center", 1159.205859, 645.605352, 36.969112, 333.129274),
    ("desk_SW", 946.292855, 1493.504426, 556.555947, 33.807080),
    ("dim01_A", 507.436515, 482.186659, 313.917066, 340.666963),
    ("dim01_B", 918.103182, 291.519992, 161.297688, 625.327321),
    ("dim02_A", 363.179032, 119.519394, 640.090414, 948.432327),
    ("dim02_B", 1359.773942, 1039.714334, 63.026927, 153.735996),
    ("dim03_A", 673.618056, 983.284722, 392.411513, 102.609542),
    ("dim03_B", 629.618056, 1267.284722, 655.658446, 127.340256),
    ("dim04_A", 147.340278, 484.284722, 786.256676, 639.316839),
    ("dim04_B", 98.784722, 781.951389, 1062.017328, 663.310371),
    ("dim05_A", 16.027778, 547.222222, 1431.645983, 1251.618142),
    ("dim05_B", 85.444444, 591.138889, 975.823057, 723.155954),
    ("dim06_A", 772.296564, 1014.048827, 334.083380, 76.682684),
    ("dim06_B", 1251.620226, 1469.864455, 329.197244, 4.831599),
    ("dim07_A", 698.368056, 752.118056, 251.165474, 160.871088),
    ("dim07_B", 1320.775051, 1342.517017, 214.282808, 32.810385),
    ("dim08_A", 100.027778, 303.222222, 1208.399050, 1266.887428),
    ("dim08_B", 169.444444, 347.138889, 752.576124, 738.425240),
    ("dim10_A", 625.027778, 8.222222, 965.281717, 1620.060643),
    ("dim10_B", 665.111111, 28.805556, 603.711244, 1202.942499),
    ("dim11_A", 944.611182, 917.518917, 175.744458, 99.508026),
    ("dim11_B", 785.036204, 1399.414354, 617.053071, 74.349876),
    ("dim12a_A", 454.277778, 39.027778, 826.047320, 1305.652096),
    ("dim12a_B", 717.361111, 1.444444, 807.487631, 1505.732275),
    ("dim12b_A", 179.611111, 201.694444, 961.043231, 1122.626381),
    ("dim12b_B", 450.756944, 40.062500, 826.793805, 1303.085856),
    ("dim13_A", 12.557270, 757.765892, 1827.223762, 1527.683722),
    ("dim13_B", 759.595556, 25.334444, 1188.400650, 1964.706729),
    ("dim14_A", 864.691339, 41.558010, 537.735285, 1256.356184),
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
