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
    ("F1", 0.000000, 749.361111, 1691.773596, 1269.687726),
    ("F2", 1.388889, 707.694444, 1700.862856, 1320.146256),
    ("F5", 625.694444, 10.888889, 1035.532028, 1667.619769),
    ("F6", 749.361111, 0.000000, 888.604422, 1600.434218),
    ("F7", 792.241152, 27.226372, 604.819998, 1276.294745),
    ("F8", 706.924637, 62.465412, 497.942443, 1053.583233),
    ("F9", 701.663933, 65.815820, 491.141665, 1038.508601),
    ("F10", 1211.301029, 529.952916, 51.273995, 481.342250),
    ("F11", 1346.954556, 593.456569, 30.334251, 501.227599),
    ("F12", 1691.773596, 888.604422, 0.000000, 454.047720),
    ("F13", 1521.499024, 1224.267556, 103.102251, 132.206814),
    ("F14", 1502.960407, 1255.097059, 122.495615, 112.027964),
    ("F15", 1269.687726, 1600.434218, 454.047720, 0.000000),
    ("F16", 1180.972108, 1623.652437, 534.733316, 4.984551),
    ("F17", 968.827936, 1544.210817, 635.377806, 34.706201),
    ("F18", 784.239564, 1434.839254, 708.493872, 81.143273),
    ("F19", 708.316722, 1371.520371, 727.765388, 104.689332),
    ("F20", 601.625705, 1277.402638, 757.827250, 146.041882),
    ("F11a", 1563.486251, 696.271400, 13.578291, 551.090483),
    ("F11b", 1626.093497, 749.860071, 8.070637, 542.524244),
    ("FC", 494.361111, 452.722222, 359.700600, 366.754560),
    ("P3", 0.694444, 752.555556, 1744.638632, 1328.997744),
    ("POB", 1417.872460, 312.896625, 208.346427, 1024.389218),
    ("P2", 841.694444, 8.222222, 1038.074971, 1837.799293),
    ("P4", 961.650783, 1617.980410, 712.021355, 50.100969),
    ("P5", 1820.080898, 2045.939326, 492.124580, 49.562685),
    ("T1", 1498.618268, 1402.474612, 199.316757, 61.280634),
    ("TC1", 2314.261138, 1824.647647, 193.033991, 283.027793),
    ("T2", 1316.080277, 430.001990, 87.074659, 696.100327),
    ("TC2", 2366.883865, 988.968283, 110.471202, 1004.534276),
    ("PA", 1856.242937, 1122.428701, 15.355491, 379.778367),
    ("T3", 291.655033, 989.782412, 914.270987, 360.104913),
    ("TC3", 412.655033, 1712.115745, 1622.277893, 607.111222),
    ("PX", 656.083155, 1533.759161, 960.924625, 186.211076),
    ("PiOB", 1375.815362, 301.233108, 206.431636, 1000.934661),
    ("Pi2", 817.617319, 4.945518, 1004.464617, 1783.225762),
    ("Pi3", 0.361111, 723.388889, 1686.320812, 1287.767507),
    ("Pi4", 944.833994, 1579.167718, 690.154254, 47.857699),
    ("Pi5", 1757.496093, 1980.578004, 472.633871, 39.902039),
    ("Ti1", 1463.086125, 1386.615960, 204.880895, 55.443276),
    ("Ti2", 1280.548133, 414.143338, 92.638798, 690.262969),
    ("Ti3", 291.905033, 962.699079, 887.838855, 354.627354),
    ("PiX", 655.754857, 1497.522986, 920.696631, 172.696145),
    ("Ai2", 1754.548362, 1064.584740, 12.209105, 352.387187),
    ("IW1_SW", 187.670816, 228.933814, 876.832602, 945.240854),
    ("IW1_SE", 1350.879217, 1304.287851, 206.675337, 48.595653),
    ("IW1_NE", 1363.795884, 1289.871184, 192.909872, 55.784760),
    ("IW1_NW", 200.587483, 214.517148, 863.067136, 952.429962),
    ("IW2_SW", 200.587483, 214.517148, 863.067136, 952.429962),
    ("IW2_SE", 206.055373, 218.485037, 837.232671, 922.728507),
    ("IW2_NE", 743.805373, 18.235037, 654.565105, 1305.834400),
    ("IW2_NW", 738.337483, 14.267148, 680.399570, 1335.535855),
    ("IW3_SW", 69.888889, 757.805556, 1204.380394, 745.784128),
    ("IW3_SE", 75.555556, 762.472222, 1189.178824, 728.004565),
    ("IW3_NE", 128.888889, 451.361111, 886.750396, 704.970438),
    ("IW3_NW", 123.222222, 446.694444, 901.951966, 722.750001),
    ("IW4_SW", 702.731982, 1330.062989, 694.102965, 99.694913),
    ("IW4_SE", 720.508206, 1346.839213, 691.010953, 94.024908),
    ("IW4_NE", 768.902739, 1055.104174, 407.379309, 71.147037),
    ("IW4_NW", 751.126515, 1038.327949, 410.471321, 76.817043),
    ("IW5_SW", 666.368056, 802.118056, 336.867114, 146.053767),
    ("IW5_SE", 1288.775051, 1392.517017, 291.322897, 17.993063),
    ("IW5_NE", 1293.795884, 1383.871184, 283.002664, 20.150117),
    ("IW5_NW", 671.388889, 793.472222, 328.546881, 148.210820),
    ("IW6_SW", 444.534722, 41.840278, 1021.645867, 1484.433471),
    ("IW6_SE", 471.733316, 52.885203, 711.745878, 1132.889436),
    ("IW6_NE", 475.254150, 51.850481, 710.819689, 1135.455676),
    ("IW6_NW", 448.055556, 40.805556, 1020.719678, 1486.999711),
    ("IW7_SW", 124.111111, 464.805556, 899.760706, 704.011034),
    ("IW7_SE", 173.694444, 506.888889, 792.832265, 577.747645),
    ("IW7_NE", 178.472222, 493.444444, 779.821955, 578.707050),
    ("IW7_NW", 128.888889, 451.361111, 886.750396, 704.970438),
    ("IW8_SW", 160.472222, 217.888889, 1186.732591, 1296.784889),
    ("IW8_SE", 187.670816, 228.933814, 876.832602, 945.240854),
    ("IW8_NE", 200.587483, 214.517148, 863.067136, 952.429962),
    ("IW8_NW", 173.388889, 203.472222, 1172.967125, 1303.973996),
    ("IW9_SW", 125.138889, 804.555556, 1082.250382, 601.741176),
    ("IW9_SE", 132.694444, 811.111111, 1068.937701, 585.850502),
    ("IW9_NE", 186.027778, 500.000000, 766.509274, 562.816375),
    ("IW9_NW", 178.472222, 493.444444, 779.821955, 578.707050),
    ("IW11_SW", 552.694444, 1195.111111, 737.681183, 161.786236),
    ("IW11_SE", 568.472222, 1209.888889, 732.590724, 154.117784),
    ("IW11_NE", 728.472222, 713.888889, 252.219554, 176.656356),
    ("IW11_NW", 712.694444, 699.111111, 257.310013, 184.324808),
    ("IW12_SW", 612.472222, 925.888889, 456.405139, 129.387070),
    ("IW12_SE", 746.570959, 1051.994616, 423.703853, 76.079860),
    ("IW12_NE", 751.126515, 1038.327949, 410.471321, 76.817043),
    ("IW12_NW", 617.027778, 912.222222, 443.172607, 130.124252),
    ("IW16_SW", 123.222222, 446.694444, 901.951966, 722.750001),
    ("IW16_SE", 128.888889, 451.361111, 886.750396, 704.970438),
    ("IW16_NE", 235.555556, 266.472222, 708.807654, 750.543137),
    ("IW16_NW", 229.888889, 261.805556, 724.009224, 768.322700),
    ("dryer_SW", 1.111111, 694.805556, 1618.652207, 1235.897488),
    ("dryer_SE", 11.562500, 696.506944, 1446.506525, 1041.194368),
    ("dryer_NE", 22.812500, 571.090278, 1324.345865, 1023.806570),
    ("dryer_NW", 12.361111, 569.388889, 1496.491546, 1218.509691),
    ("washer_SW", 12.951389, 565.423611, 1492.634802, 1218.145375),
    ("washer_SE", 23.402778, 567.125000, 1320.489120, 1023.442255),
    ("washer_NE", 47.569444, 454.625000, 1211.245127, 1018.971124),
    ("washer_NW", 37.118056, 452.923611, 1383.390808, 1213.674244),
    ("ctr_poly_0", 83.506944, 449.673611, 1028.239343, 837.210405),
    ("ctr_poly_1", 113.888889, 473.805556, 928.194809, 721.053414),
    ("ctr_poly_2", 69.888889, 757.805556, 1204.380394, 745.784128),
    ("ctr_poly_3", 39.506944, 733.673611, 1304.424928, 861.941119),
    ("bed_SW", 208.534722, 869.090278, 953.286865, 454.627084),
    ("bed_SE", 431.256944, 1072.812500, 779.512591, 231.870944),
    ("bed_NE", 505.673611, 719.006944, 435.909189, 216.555845),
    ("bed_NW", 282.951389, 515.284722, 609.683462, 439.311985),
    ("dresser_SW", 541.250000, 627.777778, 362.620139, 242.988295),
    ("dresser_SE", 665.444444, 743.472222, 309.434572, 167.889788),
    ("dresser_NE", 702.784722, 694.256944, 262.281431, 187.092794),
    ("dresser_NW", 578.590278, 578.562500, 315.466998, 262.191301),
    ("RO1_SW", 397.129016, 407.892014, 449.809585, 439.589051),
    ("RO1_SE", 504.592317, 505.855315, 359.024638, 324.313170),
    ("RO1_NE", 517.508984, 491.438648, 345.259172, 331.502277),
    ("RO1_NW", 410.045683, 393.475347, 436.044119, 446.778158),
    ("RO2_SW", 620.590278, 902.118056, 433.394041, 130.822973),
    ("RO2_SE", 604.812500, 887.340278, 438.484500, 138.491425),
    ("RO2_NE", 660.756944, 770.173611, 325.442108, 158.161326),
    ("RO2_NW", 676.534722, 784.951389, 320.351649, 150.492874),
    ("RO3_SW", 140.284722, 404.534722, 861.203180, 727.402789),
    ("RO3_SE", 145.951389, 409.201389, 846.001610, 709.623226),
    ("RO3_NE", 209.284722, 299.423611, 740.348106, 736.682016),
    ("RO3_NW", 203.618056, 294.756944, 755.549676, 754.461579),
    ("RO4_SW", 321.920816, 117.183814, 770.943413, 1027.942819),
    ("RO4_SE", 327.388706, 121.151704, 745.108948, 998.241364),
    ("RO4_NE", 446.138706, 66.790593, 694.872111, 1080.716821),
    ("RO4_NW", 440.670816, 62.822703, 720.706576, 1110.418276),
    ("RO5_SW", 447.751349, 39.153235, 901.663336, 1349.231322),
    ("RO5_SE", 469.186872, 51.088758, 724.850611, 1147.927664),
    ("RO5_NE", 472.707705, 50.054036, 723.924422, 1150.493904),
    ("RO5_NW", 451.272182, 38.118513, 900.737147, 1351.797562),
    ("RO6_SW", 569.201389, 1187.840278, 711.084781, 150.073984),
    ("RO6_SE", 553.423611, 1173.062500, 716.175240, 157.742436),
    ("RO6_NE", 591.312500, 928.506944, 478.348764, 136.446544),
    ("RO6_NW", 607.090278, 943.284722, 473.258305, 128.778092),
    ("RO7_SW", 133.812500, 780.340278, 1038.926603, 580.286405),
    ("RO7_SE", 126.256944, 773.784722, 1052.239284, 596.177079),
    ("RO7_NE", 165.868056, 530.951389, 816.135030, 576.603409),
    ("RO7_NW", 173.423611, 537.506944, 802.822349, 560.712735),
    ("O1_outer_start", 85.487847, 333.960069, 1338.036102, 1308.308593),
    ("O1_outer_end", 117.154514, 279.071181, 1285.209351, 1321.837988),
    ("O1_inner_end", 116.487847, 276.404514, 1242.806211, 1274.278862),
    ("O1_inner_start", 84.821181, 331.293403, 1295.632962, 1260.749467),
    ("O2_outer_start", 323.196181, 93.335069, 1108.807124, 1445.784635),
    ("O2_outer_end", 382.571181, 66.154514, 1083.688706, 1487.022363),
    ("O2_inner_end", 381.904514, 63.487847, 1041.285566, 1439.463237),
    ("O2_inner_start", 322.529514, 90.668403, 1066.403984, 1398.225509),
    ("O3_outer_start", 484.694444, 33.888889, 1054.624820, 1560.985126),
    ("O3_outer_end", 609.138889, 12.555556, 1036.764560, 1654.882587),
    ("O3_inner_end", 608.472222, 9.888889, 994.361420, 1607.323461),
    ("O3_inner_start", 484.027778, 31.222222, 1012.221681, 1513.426000),
    ("O4_outer_start", 725.053455, 5.434954, 760.017476, 1426.630559),
    ("O4_outer_end", 731.216872, 9.348371, 719.227360, 1380.039959),
    ("O4_inner_end", 767.216872, 8.903927, 719.651184, 1408.403213),
    ("O4_inner_start", 761.053455, 4.990509, 760.441299, 1454.993813),
    ("O5_outer_start", 685.054860, 83.151192, 457.938177, 970.920701),
    ("O5_outer_end", 831.079831, 212.176163, 249.203125, 718.359768),
    ("O5_inner_end", 863.746498, 208.398385, 246.293616, 743.389689),
    ("O5_inner_start", 717.721527, 79.373414, 455.028667, 995.950621),
    ("O6_outer_start", 992.790838, 360.387169, 129.192938, 563.546675),
    ("O6_outer_end", 1154.499806, 511.096137, 61.351303, 467.347117),
    ("O6_inner_end", 1187.166473, 507.318360, 58.441794, 492.377038),
    ("O6_inner_start", 1025.457505, 356.609391, 126.283428, 588.576595),
    ("O7_outer_start", 1641.927079, 938.411500, 4.000000, 374.347469),
    ("O7_outer_end", 1540.387529, 1135.832734, 64.000000, 183.246715),
    ("O7_inner_end", 1489.010774, 1099.581156, 64.444444, 178.326341),
    ("O7_inner_start", 1590.550325, 902.159922, 4.444444, 369.427095),
    ("O8_outer_start", 1283.994260, 1552.627378, 407.720697, 1.290996),
    ("O8_outer_end", 1313.979027, 1477.834368, 335.423958, 9.320119),
    ("O8_inner_end", 1267.531012, 1433.386353, 330.712416, 9.764563),
    ("O8_inner_start", 1237.546245, 1508.179364, 403.009156, 1.735441),
    ("O9_outer_start", 529.000000, 1209.361111, 780.865270, 180.897882),
    ("O9_outer_end", 444.506944, 1130.618056, 816.364575, 231.220646),
    ("O9_inner_end", 444.951389, 1094.618056, 781.232843, 224.028345),
    ("O9_inner_start", 529.444444, 1173.361111, 745.733538, 173.705580),
    ("O10_outer_start", 193.673611, 901.284722, 1014.198329, 484.481250),
    ("O10_outer_end", 144.000000, 857.361111, 1084.517078, 569.623459),
    ("O10_inner_end", 144.444444, 821.361111, 1049.385346, 562.431158),
    ("O10_inner_start", 194.118056, 865.284722, 979.066597, 477.288948),
    ("O11_outer_start", 36.000000, 767.361111, 1352.145337, 883.655592),
    ("O11_outer_end", 19.506944, 755.618056, 1434.776406, 978.532128),
    ("O11_inner_end", 19.951389, 719.618056, 1399.644674, 971.339826),
    ("O11_inner_start", 36.444444, 731.361111, 1317.013605, 876.463291),
    ("O3_door_tip", 608.340278, 7.812500, 863.032960, 1455.745617),
    ("O6_door_tip", 1007.498213, 537.455656, 88.553606, 362.481625),
    ("RO1_door_tip", 439.403224, 591.249555, 445.804331, 296.487212),
    ("RO2_door_tip", 817.888889, 922.361111, 288.772826, 93.508682),
    ("RO3_door_tip", 162.500000, 266.694444, 896.258654, 913.591068),
    ("RO4_door_tip", 416.152979, 51.110421, 872.798713, 1280.286426),
    ("RO5_door_tip", 605.874515, 21.748624, 704.439010, 1255.634211),
    ("RO6_door_tip_S", 449.388889, 1071.472222, 753.771152, 216.876728),
    ("RO6_door_tip_N", 486.055556, 834.805556, 523.616498, 196.267800),
    ("RO7_door_tip_S", 193.138889, 828.111111, 947.728313, 474.534593),
    ("RO7_door_tip_N", 231.472222, 593.111111, 719.240326, 455.592332),
    ("w_f8f9_start", 695.300784, 52.841560, 528.055063, 1088.851840),
    ("w_f8f9_arc_entry", 664.009647, 56.201341, 530.589306, 1064.821624),
    ("w_f8f9_arc_center", 667.394457, 58.987726, 521.486038, 1054.175620),
    ("w_f8f9_arc_exit", 657.713372, 60.211278, 522.449793, 1046.779545),
    ("w_f8f9_end", 668.997267, 69.593598, 494.051174, 1013.478681),
    ("g_f8f9_start", 698.123414, 55.164189, 520.443575, 1079.951355),
    ("g_f8f9_arc_entry", 666.832277, 58.523971, 522.977818, 1055.921140),
    ("g_f8f9_arc_center", 667.394457, 58.987726, 521.486038, 1054.175620),
    ("g_f8f9_arc_exit", 665.796705, 59.183500, 521.639082, 1052.953692),
    ("g_f8f9_end", 677.080600, 68.565820, 493.240464, 1019.652828),
    ("wh_center", 696.189563, 32.460225, 597.694752, 1189.848658),
    ("sink_NW", 725.270456, 117.866787, 383.378471, 882.182033),
    ("sink_SE", 735.800461, 226.480125, 265.869342, 651.852193),
    ("stove_NW", 652.429988, 67.442986, 506.476140, 1017.359855),
    ("stove_SE", 600.956936, 131.469934, 418.560908, 815.814098),
    ("dw_NW", 836.346906, 216.943237, 244.036081, 711.903728),
    ("dw_SE", 812.231780, 308.828111, 193.912893, 549.439012),
    ("fridge_NW", 302.393484, 140.962038, 748.932768, 965.207403),
    ("fridge_SE", 255.645732, 245.471230, 686.972100, 759.904409),
    ("ice_NW", 922.366336, 308.129334, 166.016566, 601.493986),
    ("ice_SE", 918.512791, 371.828566, 140.084022, 508.973448),
    ("nctr_NW", 616.722039, 27.818371, 660.892432, 1207.388865),
    ("nctr_SE", 545.446044, 84.209041, 535.212967, 953.734600),
    ("wctr_NW", 295.160927, 215.903092, 644.130787, 767.516033),
    ("wctr_SE", 367.881490, 355.623654, 483.874197, 505.725830),
    ("loveseat_NW", 835.440201, 492.448567, 149.505895, 344.357581),
    ("et_center", 951.910497, 845.299274, 172.076288, 131.125510),
    ("loveseat2_NW", 1070.734778, 846.556625, 117.241574, 156.519021),
    ("loveseat2_SE", 1321.478414, 1240.494706, 184.522135, 59.672189),
    ("chair_center", 1445.440250, 730.268402, 10.295107, 425.587281),
    ("ottoman_center", 1231.606778, 683.080432, 36.969112, 337.592326),
    ("desk_SW", 946.292855, 1493.504426, 608.917791, 33.807080),
    ("dim01_A", 516.467770, 490.481688, 345.976179, 332.436488),
    ("dim01_B", 927.134437, 299.815021, 169.635939, 617.096846),
    ("dim02_A", 514.640069, 46.735737, 676.927057, 1130.661683),
    ("dim02_B", 1533.919595, 985.866496, 23.643142, 270.754806),
    ("dim03_A", 677.746918, 987.167080, 438.279824, 100.958792),
    ("dim03_B", 633.746918, 1271.167080, 714.465409, 125.689506),
    ("dim04_A", 147.340278, 484.284722, 844.733986, 639.316839),
    ("dim04_B", 98.784722, 781.951389, 1134.152103, 663.310371),
    ("dim05_A", 16.027778, 547.222222, 1504.778436, 1251.618142),
    ("dim05_B", 85.444444, 591.138889, 1042.055069, 723.155954),
    ("dim06_A", 781.088192, 1022.347445, 375.693302, 73.915087),
    ("dim06_B", 1251.620226, 1469.864455, 365.942383, 4.831599),
    ("dim07_A", 698.368056, 752.118056, 288.820718, 160.871088),
    ("dim07_B", 1320.775051, 1342.517017, 243.276501, 32.810385),
    ("dim08_A", 100.027778, 303.222222, 1268.592851, 1266.887428),
    ("dim08_B", 169.444444, 347.138889, 805.869484, 738.425240),
    ("dim10_A", 625.027778, 8.222222, 993.128888, 1620.060643),
    ("dim10_B", 652.226372, 19.267148, 683.228899, 1268.516609),
    ("dim11_A", 944.611182, 917.518917, 207.319164, 99.508026),
    ("dim11_B", 785.036204, 1399.414354, 673.931929, 74.349876),
    ("dim12a_A", 454.277778, 39.027778, 859.996008, 1305.652096),
    ("dim12a_B", 717.361111, 1.444444, 829.575888, 1505.732275),
    ("dim12b_A", 179.611111, 201.694444, 1012.243454, 1122.626381),
    ("dim12b_B", 450.756944, 40.062500, 860.922196, 1303.085856),
    ("dim13_A", 12.557270, 757.765892, 1911.133155, 1527.683722),
    ("dim13_B", 759.595556, 25.334444, 1213.949239, 1964.706729),
    ("dim14_A", 865.500060, 42.213062, 549.473145, 1253.561951),
    ("dim14_B", 1488.301863, 610.226541, 28.939228, 591.784116),
    ("dim15_A", 9.014621, 918.560179, 1906.886512, 1370.359886),
    ("dim15_B", 1245.202964, 2046.740489, 889.137661, 74.167079),
    ("dim16_A", 486.279514, 1133.071181, 762.564788, 197.948560),
    ("dim16_B", 646.279514, 637.071181, 282.193618, 220.487132),
    ("dim17_A", 168.362847, 842.404514, 1013.307569, 518.941650),
    ("dim17_B", 328.362847, 346.404514, 532.936399, 541.480222),
    ("dim18_A", 148.250000, 644.444444, 906.612376, 563.222328),
    ("dim18_B", 568.250000, 1028.444444, 575.355858, 139.158061),
    ("dim19_A", 27.571181, 724.862847, 1357.702404, 923.274822),
    ("dim19_B", 187.571181, 228.862847, 877.331233, 945.813394),
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
