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
    f_names = [f"F{i}" for i in range(1, 19) if i not in (3, 4)] + ["F11a", "F11b"]
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
    ("F1", 0.000000, 749.361111, 1691.790045, 1257.780120),
    ("F2", 1.388889, 707.694444, 1704.156621, 1310.126930),
    ("F5", 625.694444, 10.888889, 1090.939806, 1712.489968),
    ("F6", 749.361111, 0.000000, 944.899113, 1650.616591),
    ("F7", 791.555556, 26.694444, 654.458151, 1329.422146),
    ("F8", 706.000000, 61.694444, 538.165073, 1101.184412),
    ("F9", 700.722222, 65.027778, 530.691758, 1085.715050),
    ("F10", 1276.086627, 590.856284, 42.998506, 499.737837),
    ("F11", 1367.535778, 637.167807, 30.334251, 513.478834),
    ("F12", 1691.790045, 944.899113, 0.000000, 450.361446),
    ("F13", 1450.196321, 1284.390512, 124.742314, 104.532468),
    ("F14", 1437.388889, 1306.805556, 139.667991, 91.675886),
    ("F15", 1257.780120, 1650.616591, 450.361446, 0.000000),
    ("F16", 1162.611009, 1669.413831, 533.780530, 5.444444),
    ("F17", 905.386867, 1547.394705, 644.070549, 44.105085),
    ("F18", 836.444868, 1499.041909, 666.247714, 60.089573),
    ("F11a", 1584.443500, 752.772300, 13.578291, 555.984678),
    ("F11b", 1643.316892, 808.645692, 8.070637, 544.524737),
    ("FC", 494.361111, 452.722222, 357.447071, 385.603993),
    ("P3", 0.694444, 752.555556, 1746.135361, 1317.085676),
    ("POB", 1417.872460, 312.896625, 246.868124, 1087.344490),
    ("P2", 841.694444, 8.222222, 1102.108517, 1891.754655),
    ("P4", 961.650783, 1617.980410, 656.952731, 38.359441),
    ("P5", 1820.080898, 2045.939326, 433.142215, 54.390856),
    ("T1", 1498.618268, 1402.474612, 163.925993, 80.172951),
    ("TC1", 2314.261138, 1824.647647, 159.263673, 317.392804),
    ("T2", 1316.080277, 430.001990, 109.558747, 749.494769),
    ("TC2", 2366.883865, 988.968283, 134.980847, 1077.269586),
    ("PA", 1856.242937, 1122.428701, 7.687269, 425.326321),
    ("T3", 291.655033, 989.782412, 883.951320, 348.288769),
    ("TC3", 412.655033, 1712.115745, 1568.237364, 570.310881),
    ("PX", 656.083155, 1533.759161, 907.868429, 165.949483),
    ("PiOB", 1375.815362, 301.233108, 244.086899, 1062.648076),
    ("Pi2", 817.617319, 4.945518, 1066.730307, 1836.257264),
    ("Pi3", 0.361111, 723.388889, 1688.007594, 1276.993761),
    ("Pi4", 944.833994, 1579.167718, 636.474677, 37.062882),
    ("Pi5", 1757.496093, 1980.578004, 415.215292, 44.937129),
    ("Ti1", 1463.086125, 1386.615960, 169.409109, 73.561959),
    ("Ti2", 1280.548133, 414.143338, 115.041863, 742.883777),
    ("Ti3", 291.905033, 962.699079, 858.597409, 343.946855),
    ("PiX", 655.754857, 1497.522986, 868.902999, 153.916762),
    ("Ai2", 1754.548362, 1064.584740, 4.580883, 396.224660),
    ("IW1_SW", 187.138889, 228.555556, 897.661113, 965.201628),
    ("IW1_SE", 1350.694444, 1304.111111, 172.705977, 65.646073),
    ("IW1_NE", 1363.611111, 1289.694444, 160.018733, 73.970825),
    ("IW1_NW", 200.055556, 214.138889, 884.973869, 973.526381),
    ("IW2_SW", 213.472222, 224.055556, 823.279315, 904.887492),
    ("IW2_SE", 220.055556, 229.138889, 797.672125, 876.304158),
    ("IW2_NE", 324.472222, 142.222222, 722.861414, 948.577426),
    ("IW2_NW", 317.888889, 137.138889, 748.468604, 977.160759),
    ("IW2o_SW", 326.673889, 140.627222, 722.122343, 950.839558),
    ("IW2o_SE", 315.687222, 138.733889, 749.207674, 974.898627),
    ("IW2o_NE", 451.020556, 58.900556, 754.820939, 1155.552203),
    ("IW2o_NW", 462.007222, 60.793889, 727.735609, 1131.493134),
    ("IW2s_SW", 453.805556, 57.888889, 754.665202, 1158.397668),
    ("IW2s_SE", 459.222222, 61.805556, 727.891346, 1128.647668),
    ("IW2s_NE", 743.222222, 17.805556, 704.644412, 1357.544698),
    ("IW2s_NW", 737.805556, 13.888889, 731.418269, 1387.294698),
    ("IW3_SW", 69.888889, 757.805556, 1191.031666, 735.435346),
    ("IW3_SE", 75.555556, 762.472222, 1175.237984, 717.657568),
    ("IW3_NE", 128.888889, 451.361111, 887.185836, 709.765379),
    ("IW3_NW", 123.222222, 446.694444, 902.979518, 727.543157),
    ("IW4_SW", 693.888889, 1327.805556, 656.172833, 93.435346),
    ("IW4_SE", 711.555556, 1344.472222, 652.379151, 87.657568),
    ("IW4_NE", 760.111111, 1046.805556, 376.618499, 78.048877),
    ("IW4_NW", 742.444444, 1030.138889, 380.412181, 83.826655),
    ("IW5_SW", 666.368056, 802.118056, 315.932260, 156.797434),
    ("IW5_SE", 1288.590278, 1392.340278, 251.423321, 28.797434),
    ("IW5_NE", 1293.611111, 1383.694444, 243.642199, 31.522310),
    ("IW5_NW", 671.388889, 793.472222, 308.151138, 159.522310),
    ("IW6_SW", 444.534722, 41.840278, 1067.423357, 1520.411351),
    ("IW6_SE", 471.201389, 52.506944, 750.724443, 1171.966907),
    ("IW6_NE", 474.722222, 51.472222, 749.977958, 1174.722421),
    ("IW6_NW", 448.055556, 40.805556, 1066.676871, 1523.166865),
    ("IW7_SW", 124.111111, 464.805556, 899.477332, 708.048877),
    ("IW7_SE", 173.694444, 506.888889, 788.108050, 581.798877),
    ("IW7_NE", 178.472222, 493.444444, 775.816554, 583.515379),
    ("IW7_NW", 128.888889, 451.361111, 887.185836, 709.765379),
    ("IW8_SW", 160.472222, 217.888889, 1214.360027, 1313.646073),
    ("IW8_SE", 187.138889, 228.555556, 897.661113, 965.201628),
    ("IW8_NE", 200.055556, 214.138889, 884.973869, 973.526381),
    ("IW8_NW", 173.388889, 203.472222, 1201.672782, 1321.970825),
    ("IW9_SW", 125.138889, 804.555556, 1063.868702, 591.407568),
    ("IW9_SE", 132.694444, 811.111111, 1049.963909, 575.518679),
    ("IW9_NE", 186.027778, 500.000000, 761.911761, 567.626490),
    ("IW9_NW", 178.472222, 493.444444, 775.816554, 583.515379),
    ("IW11_SW", 552.694444, 1195.111111, 697.391353, 151.518679),
    ("IW11_SE", 568.472222, 1209.888889, 691.708782, 143.852012),
    ("IW11_NE", 728.472222, 713.888889, 237.214915, 193.646073),
    ("IW11_NW", 712.694444, 699.111111, 242.897486, 201.312739),
    ("IW12_SW", 612.472222, 925.888889, 428.461849, 132.749042),
    ("IW12_SE", 737.888889, 1043.805556, 392.925900, 82.332376),
    ("IW12_NE", 742.444444, 1030.138889, 380.412181, 83.826655),
    ("IW12_NW", 617.027778, 912.222222, 415.948130, 134.243322),
    ("IW16_SW", 123.222222, 446.694444, 902.979518, 727.543157),
    ("IW16_SE", 128.888889, 451.361111, 887.185836, 709.765379),
    ("IW16_NE", 235.555556, 266.472222, 720.744118, 767.451628),
    ("IW16_NW", 229.888889, 261.805556, 736.537800, 785.229406),
    ("dryer_SW", 1.111111, 694.805556, 1620.232985, 1226.262958),
    ("dryer_SE", 11.562500, 696.506944, 1442.906322, 1031.575458),
    ("dryer_NE", 22.812500, 571.090278, 1326.136767, 1019.865888),
    ("dryer_NW", 12.361111, 569.388889, 1503.463429, 1214.553388),
    ("washer_SW", 12.951389, 565.423611, 1499.786389, 1214.378346),
    ("washer_SE", 23.402778, 567.125000, 1322.459726, 1019.690846),
    ("washer_NE", 47.569444, 454.625000, 1218.606837, 1020.897942),
    ("washer_NW", 37.118056, 452.923611, 1395.933500, 1215.585442),
    ("ctr_poly_0", 83.506944, 449.673611, 1031.529968, 840.478209),
    ("ctr_poly_1", 113.888889, 473.805556, 927.784733, 724.332376),
    ("ctr_poly_2", 69.888889, 757.805556, 1191.031666, 735.435346),
    ("ctr_poly_3", 39.506944, 733.673611, 1294.776902, 851.581179),
    ("bed_SW", 208.534722, 869.090278, 929.491498, 444.689430),
    ("bed_SE", 431.256944, 1072.812500, 744.467093, 221.967208),
    ("bed_NE", 505.673611, 719.006944, 417.755819, 224.443886),
    ("bed_NW", 282.951389, 515.284722, 602.780224, 447.166108),
    ("dresser_SW", 541.250000, 627.777778, 349.942551, 256.174675),
    ("dresser_SE", 665.444444, 743.472222, 291.724031, 181.091342),
    ("dresser_NE", 702.784722, 694.256944, 247.985257, 203.890558),
    ("dresser_NW", 578.590278, 578.562500, 306.203777, 278.973892),
    ("RO1_SW", 432.694444, 440.111111, 413.564811, 413.646073),
    ("RO1_SE", 547.222222, 545.138889, 324.219275, 305.451628),
    ("RO1_NE", 560.138889, 530.722222, 311.532031, 313.776381),
    ("RO1_NW", 445.611111, 425.694444, 400.877566, 421.970825),
    ("RO2_SW", 620.590278, 902.118056, 406.708675, 135.509865),
    ("RO2_SE", 604.812500, 887.340278, 412.391246, 143.176531),
    ("RO2_NE", 660.756944, 770.173611, 306.177586, 170.038853),
    ("RO2_NW", 676.534722, 784.951389, 300.495015, 162.372186),
    ("RO3_SW", 140.284722, 404.534722, 864.566878, 734.656510),
    ("RO3_SE", 145.951389, 409.201389, 848.773195, 716.878732),
    ("RO3_NE", 209.284722, 299.423611, 749.948425, 751.129942),
    ("RO3_NW", 203.618056, 294.756944, 765.742107, 768.907720),
    ("RO4_SW", 341.080556, 129.213889, 720.962602, 970.684654),
    ("RO4_SE", 330.093889, 127.320556, 748.047933, 994.743723),
    ("RO4_NE", 432.947222, 66.647222, 752.314014, 1132.040441),
    ("RO4_NW", 443.933889, 68.540556, 725.228683, 1107.981371),
    ("RO5_SW", 447.569444, 39.125000, 947.061073, 1388.640518),
    ("RO5_SE", 468.680556, 50.736111, 764.298871, 1187.029407),
    ("RO5_NE", 472.201389, 49.701389, 763.552386, 1189.784921),
    ("RO5_NW", 451.090278, 38.090278, 946.314588, 1391.396032),
    ("RO6_SW", 569.201389, 1187.840278, 671.101356, 140.754584),
    ("RO6_SE", 553.423611, 1173.062500, 676.783927, 148.421250),
    ("RO6_NE", 591.312500, 928.506944, 450.099068, 138.860360),
    ("RO6_NW", 607.090278, 943.284722, 444.416497, 131.193693),
    ("RO7_SW", 133.812500, 780.340278, 1021.210735, 571.279501),
    ("RO7_SE", 126.256944, 773.784722, 1035.115528, 587.168390),
    ("RO7_NE", 165.868056, 530.951389, 810.152891, 579.329722),
    ("RO7_NW", 173.423611, 537.506944, 796.248098, 563.440833),
    ("O1_outer_start", 85.487847, 333.960069, 1359.390068, 1317.311326),
    ("O1_outer_end", 117.154514, 279.071181, 1309.977683, 1334.436931),
    ("O1_inner_end", 116.487847, 276.404514, 1266.390318, 1286.881376),
    ("O1_inner_start", 84.821181, 331.293403, 1315.802704, 1269.755770),
    ("O2_outer_start", 252.930604, 136.578752, 1185.091733, 1423.464784),
    ("O2_outer_end", 305.730372, 102.822965, 1156.812450, 1461.723491),
    ("O2_inner_end", 305.063706, 100.156298, 1113.225086, 1414.167935),
    ("O2_inner_start", 252.263937, 133.912085, 1141.504369, 1375.909228),
    ("O3_outer_start", 484.694444, 33.888889, 1103.563273, 1599.041453),
    ("O3_outer_end", 609.138889, 12.555556, 1091.453525, 1698.995688),
    ("O3_inner_end", 608.472222, 9.888889, 1047.866160, 1651.440133),
    ("O3_inner_start", 484.027778, 31.222222, 1059.975908, 1551.485897),
    ("O4_outer_start", 724.862847, 5.321181, 812.364586, 1476.921435),
    ("O4_outer_end", 730.987847, 9.196181, 770.203802, 1430.296435),
    ("O4_inner_end", 766.987847, 8.751736, 772.065254, 1460.173882),
    ("O4_inner_start", 760.862847, 4.876736, 814.226038, 1506.798882),
    ("O5_outer_start", 684.027778, 82.277778, 494.484989, 1016.532047),
    ("O5_outer_end", 829.472222, 210.722222, 275.103505, 763.420936),
    ("O5_inner_end", 862.138889, 206.944444, 273.631623, 789.965050),
    ("O5_inner_start", 716.694444, 78.500000, 493.013108, 1043.076161),
    ("O6_outer_start", 1046.792464, 410.506566, 116.097252, 569.621783),
    ("O6_outer_end", 1217.991327, 570.705429, 51.232278, 482.931757),
    ("O6_inner_end", 1250.657994, 566.927651, 49.760397, 509.475871),
    ("O6_inner_start", 1079.459131, 406.728788, 114.625370, 596.165897),
    ("O7_outer_start", 1630.190211, 987.354255, 4.000000, 370.096194),
    ("O7_outer_end", 1493.390708, 1162.719681, 64.000000, 177.300436),
    ("O7_inner_end", 1443.541236, 1124.234576, 64.444444, 174.326785),
    ("O7_inner_start", 1580.340739, 948.869150, 4.444444, 367.122543),
    ("O8_outer_start", 1276.333520, 1576.015644, 381.654519, 2.903762),
    ("O8_outer_end", 1304.141633, 1499.045980, 311.314306, 13.109538),
    ("O8_inner_end", 1257.697189, 1454.601535, 307.790559, 13.553982),
    ("O8_inner_start", 1229.889076, 1531.571200, 378.130772, 3.348207),
    ("O9_outer_start", 529.000000, 1209.361111, 740.025980, 169.113453),
    ("O9_outer_end", 444.506944, 1130.618056, 778.929930, 219.425953),
    ("O9_inner_end", 444.951389, 1094.618056, 745.235826, 213.747846),
    ("O9_inner_start", 529.444444, 1173.361111, 706.331876, 163.435346),
    ("O10_outer_start", 193.673611, 901.284722, 989.494095, 472.648176),
    ("O10_outer_end", 144.000000, 857.361111, 1063.217489, 557.780120),
    ("O10_inner_end", 144.444444, 821.361111, 1029.523385, 552.102012),
    ("O10_inner_start", 194.118056, 865.284722, 955.799991, 466.970068),
    ("O11_outer_start", 36.000000, 767.361111, 1341.503767, 871.780120),
    ("O11_outer_end", 19.506944, 755.618056, 1426.947368, 966.648176),
    ("O11_inner_end", 19.951389, 719.618056, 1393.253264, 960.970068),
    ("O11_inner_start", 36.444444, 731.361111, 1307.809663, 866.102012),
    ("O3_door_tip", 608.340278, 7.812500, 912.509267, 1499.684619),
    ("O6_door_tip", 1070.774055, 596.849268, 71.538197, 370.657719),
    ("RO1_door_tip", 481.847222, 630.347222, 405.030873, 271.193267),
    ("RO2_door_tip", 817.888889, 922.361111, 263.703536, 105.213894),
    ("RO3_door_tip", 162.500000, 266.694444, 911.304335, 927.832761),
    ("RO4_door_tip", 510.334444, 89.231111, 584.827655, 980.723974),
    ("RO5_door_tip", 605.376736, 21.404514, 750.603013, 1301.652554),
    ("RO6_door_tip_S", 449.388889, 1071.472222, 718.704328, 207.732320),
    ("RO6_door_tip_N", 486.055556, 834.805556, 499.331883, 198.479845),
    ("RO7_door_tip_S", 193.138889, 828.111111, 926.047364, 465.729460),
    ("RO7_door_tip_N", 231.472222, 593.111111, 708.341586, 458.143652),
    ("w_f8f9_start", 694.444444, 52.138889, 569.530215, 1136.517746),
    ("w_f8f9_arc_entry", 663.153307, 55.498670, 570.697579, 1111.047854),
    ("w_f8f9_arc_center", 666.517682, 58.264620, 561.219541, 1100.382483),
    ("w_f8f9_arc_exit", 656.836597, 59.488172, 561.753139, 1092.533342),
    ("w_f8f9_end", 668.055556, 68.805556, 532.163639, 1059.170936),
    ("g_f8f9_start", 697.250000, 54.444444, 561.605596, 1127.601079),
    ("g_f8f9_arc_entry", 665.958863, 57.804226, 562.772960, 1102.131187),
    ("g_f8f9_arc_center", 666.517682, 58.264620, 561.219541, 1100.382483),
    ("g_f8f9_arc_exit", 664.919930, 58.460394, 561.301835, 1099.086037),
    ("g_f8f9_end", 676.138889, 67.777778, 531.712335, 1065.723631),
    ("wh_center", 695.486893, 31.911223, 643.004634, 1238.730896),
    ("sink_NW", 724.055556, 116.805556, 416.480849, 927.615380),
    ("sink_SE", 734.201389, 225.034722, 287.613403, 692.378870),
    ("stove_NW", 651.513889, 66.680556, 544.519190, 1062.508560),
    ("stove_SE", 599.784722, 130.451389, 447.055007, 855.609673),
    ("dw_NW", 834.722222, 215.472222, 269.623330, 756.948714),
    ("dw_SE", 810.368056, 307.118056, 210.264322, 589.147050),
    ("fridge_NW", 316.951389, 152.173611, 716.314593, 926.534535),
    ("fridge_SE", 276.292101, 262.771267, 649.304847, 720.710024),
    ("ice_NW", 920.451389, 306.368056, 185.741488, 645.696060),
    ("ice_SE", 916.446736, 369.916181, 154.198425, 550.041781),
    ("nctr_NW", 616.138889, 27.388889, 705.580635, 1253.420936),
    ("nctr_SE", 544.555556, 83.472222, 568.873718, 993.797173),
    ("wctr_NW", 316.365017, 233.760851, 609.272725, 731.719822),
    ("wctr_SE", 400.240017, 384.635851, 448.054226, 477.703898),
    ("loveseat_NW", 828.244635, 485.746011, 154.361693, 379.576471),
    ("et_center", 929.642675, 824.258982, 157.991000, 158.022486),
    ("loveseat2_NW", 1047.659342, 824.708720, 105.116534, 187.375331),
    ("loveseat2_SE", 1293.970229, 1214.214051, 152.052902, 79.500162),
    ("chair_center", 1450.938952, 774.717369, 10.295107, 431.137842),
    ("ottoman_center", 1228.865113, 711.382376, 36.969112, 349.204312),
    ("desk_SW", 886.133598, 1497.579620, 616.712778, 42.993974),
    ("dim01_A", 531.507954, 504.316984, 330.665310, 337.325654),
    ("dim01_B", 942.174621, 313.650317, 178.045932, 646.970209),
    ("dim02_A", 455.130062, 79.189312, 672.083617, 1046.695828),
    ("dim02_B", 1493.854067, 1039.576249, 29.666817, 248.850791),
    ("dim03_A", 673.618056, 983.284722, 409.131374, 105.978209),
    ("dim03_B", 629.618056, 1267.284722, 672.378307, 117.081179),
    ("dim04_A", 147.340278, 484.284722, 842.230191, 643.361377),
    ("dim04_B", 98.784722, 781.951389, 1117.990843, 652.970068),
    ("dim05_A", 16.027778, 547.222222, 1513.716708, 1248.794807),
    ("dim05_B", 85.444444, 591.138889, 1035.894481, 720.378140),
    ("dim06_A", 763.735840, 1036.529076, 367.221570, 79.338247),
    ("dim06_B", 1242.874729, 1492.167965, 342.042263, 7.532692),
    ("dim07_A", 698.368056, 752.118056, 271.120526, 175.021691),
    ("dim07_B", 1320.590278, 1342.340278, 206.611588, 47.021691),
    ("dim08_A", 100.027778, 303.222222, 1290.469775, 1277.691837),
    ("dim08_B", 169.444444, 347.138889, 812.647548, 749.275170),
    ("dim10_A", 625.027778, 8.222222, 1047.352442, 1664.934412),
    ("dim10_B", 665.111111, 28.805556, 668.958974, 1247.851079),
    ("dim11_A", 996.889312, 967.041909, 178.059744, 104.205526),
    ("dim11_B", 836.889312, 1463.041909, 632.553610, 54.411466),
    ("dim12a_A", 454.277778, 39.027778, 901.216303, 1341.833532),
    ("dim12a_B", 717.361111, 1.444444, 882.656615, 1554.405810),
    ("dim12b_A", 179.611111, 201.694444, 1036.212214, 1140.637492),
    ("dim12b_B", 450.756944, 40.062500, 901.962789, 1339.078018),
    ("dim13_A", 12.484444, 772.445556, 1931.754186, 1518.775676),
    ("dim13_B", 759.595556, 25.334444, 1279.184824, 2014.862146),
    ("dim14_A", 864.691339, 41.558010, 599.009725, 1308.843955),
    ("dim14_B", 1571.999449, 689.267134, 27.302806, 618.886946),
    ("dim15_A", 9.694444, 925.555556, 1908.758828, 1353.637161),
    ("dim15_B", 1245.694444, 2053.555556, 827.041161, 57.637161),
    ("dim16_A", 486.279514, 1133.071181, 724.865448, 187.673193),
    ("dim16_B", 646.279514, 637.071181, 270.371582, 237.467253),
    ("dim17_A", 168.362847, 842.404514, 991.743286, 508.617637),
    ("dim17_B", 328.362847, 346.404514, 537.249419, 558.411698),
    ("dim18_A", 148.250000, 644.444444, 894.826723, 560.461473),
    ("dim18_B", 568.250000, 1028.444444, 542.254168, 136.461473),
    ("dim19_A", 27.571181, 724.862847, 1349.904728, 912.909304),
    ("dim19_B", 187.571181, 228.862847, 895.410861, 962.703364),
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

    _F_NAMES = [f"F{i}" for i in range(1, 19) if i not in (3, 4)] + ["F11a", "F11b"]

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
