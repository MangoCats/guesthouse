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
    RO5_DOOR_WIDTH, RO6_DOOR_WIDTH, RO7_DOOR_WIDTH,
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
        ("IW12", layout.iw12),
        ("dryer", layout.dryer), ("washer", layout.washer),
    ]:
        for i, label in enumerate(["SW", "SE", "NE", "NW"]):
            result.append((f"{prefix}_{label}", wall.poly[i]))

    # ---- Counter clip polygon (4 or 5 pts) ----
    for i, p in enumerate(layout.ctr_clip):
        result.append((f"ctr_poly_{i}", p))

    # ---- Furniture (bed, dresser) ----
    for prefix, wall in [("bed", layout.bed), ("dresser", layout.dresser),
                          ("shelves", layout.shelves)]:
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

    # IW9 unit vectors (shared by RO3 and RO7)
    _i9_sw, _i9_se, _i9_ne, _i9_nw = layout.iw9.poly
    _i9_dx_n = _i9_ne[0] - _i9_se[0]; _i9_dy_n = _i9_ne[1] - _i9_se[1]
    _i9_ln = math.sqrt(_i9_dx_n**2 + _i9_dy_n**2)
    _i9_an = (_i9_dx_n / _i9_ln, _i9_dy_n / _i9_ln)
    _i9_dx_t = _i9_sw[0] - _i9_se[0]; _i9_dy_t = _i9_sw[1] - _i9_se[1]
    _i9_lt = math.sqrt(_i9_dx_t**2 + _i9_dy_t**2)
    _i9_at = (_i9_dx_t / _i9_lt, _i9_dy_t / _i9_lt)

    # RO3 door tip: in IW9, hinged at south edge, swings west
    ro3 = [r for r in rough_openings if r.name == "RO3"][0]
    _ro3_dx_l = ro3.poly[3][0] - ro3.poly[0][0]
    _ro3_dy_l = ro3.poly[3][1] - ro3.poly[0][1]
    _ro3_len = math.sqrt(_ro3_dx_l**2 + _ro3_dy_l**2)
    _ro3_gap = (_ro3_len - RO3_DOOR_WIDTH) / 2
    _ro3_s_ctr = ((ro3.poly[0][0] + ro3.poly[1][0]) / 2,
                  (ro3.poly[0][1] + ro3.poly[1][1]) / 2)
    _ro3_hinge = (_ro3_s_ctr[0] + _ro3_gap * _i9_an[0],
                  _ro3_s_ctr[1] + _ro3_gap * _i9_an[1])
    tips.append(("RO3_door_tip",
                 (_ro3_hinge[0] + RO3_DOOR_WIDTH * _i9_at[0],
                  _ro3_hinge[1] + RO3_DOOR_WIDTH * _i9_at[1])))

    # RO4 door tip (hinged at SE end midpoint, swings WSW)
    ro4 = [r for r in rough_openings if r.name == "RO4"][0]
    _ro4_hinge = ((ro4.poly[0][0] + ro4.poly[1][0]) / 2,
                  (ro4.poly[0][1] + ro4.poly[1][1]) / 2)
    _ro4_closed = ((ro4.poly[2][0] + ro4.poly[3][0]) / 2,
                   (ro4.poly[2][1] + ro4.poly[3][1]) / 2)
    _ro4_dx = _ro4_closed[0] - _ro4_hinge[0]
    _ro4_dy = _ro4_closed[1] - _ro4_hinge[1]
    _ro4_door_len = math.sqrt(_ro4_dx**2 + _ro4_dy**2)
    _ro4_sw_dx = ro4.poly[2][0] - ro4.poly[3][0]
    _ro4_sw_dy = ro4.poly[2][1] - ro4.poly[3][1]
    _ro4_sw_len = math.sqrt(_ro4_sw_dx**2 + _ro4_sw_dy**2)
    _ro4_swing = (_ro4_sw_dx / _ro4_sw_len, _ro4_sw_dy / _ro4_sw_len)
    tips.append(("RO4_door_tip",
                 (_ro4_hinge[0] + _ro4_door_len * _ro4_swing[0],
                  _ro4_hinge[1] + _ro4_door_len * _ro4_swing[1])))

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
    ("F1", 0.000000, 749.361111, 1645.186559, 1259.288658),
    ("F2", 1.388889, 707.694444, 1655.927050, 1311.366429),
    ("F5", 625.694444, 10.888889, 1036.235807, 1705.927361),
    ("F6", 749.361111, 0.000000, 893.497918, 1643.300677),
    ("F7", 791.555556, 26.694444, 611.754497, 1322.106232),
    ("F8", 706.000000, 61.694444, 500.014457, 1094.621805),
    ("F9", 700.722222, 65.027778, 492.866359, 1079.206251),
    ("F10", 1225.138889, 542.944444, 45.469427, 513.456251),
    ("F11", 1334.368950, 598.961875, 29.454377, 525.820978),
    ("F12", 1645.186559, 893.497918, 0.000000, 454.926314),
    ("F13", 1442.249887, 1290.289948, 138.757963, 96.234096),
    ("F14", 1432.437656, 1311.466947, 152.568830, 85.331202),
    ("F15", 1259.288658, 1643.300677, 454.926314, 0.000000),
    ("F16", 1195.786713, 1661.940790, 513.288814, 2.689679),
    ("F17", 923.197764, 1562.446685, 645.633848, 42.307632),
    ("F18", 859.485883, 1520.896025, 666.566042, 56.812934),
    ("F11a", 1538.702356, 709.435551, 12.872709, 559.297314),
    ("F11b", 1595.972819, 763.706014, 7.445520, 546.234444),
    ("FC", 494.361111, 452.722222, 336.966812, 382.754113),
    ("P3", 0.730711, 796.500679, 1698.993662, 1279.022934),
    ("POB", 1387.503253, 284.080887, 246.529512, 1110.169917),
    ("P2", 793.470345, 9.734599, 1085.887680, 1892.068392),
    ("P4", 1002.286178, 1620.956317, 629.308333, 29.158984),
    ("P5", 1862.679270, 2024.090528, 420.000821, 64.805511),
    ("T1", 1523.561765, 1378.940128, 153.882297, 93.590455),
    ("TC1", 2337.250715, 1781.731507, 160.957191, 346.452123),
    ("T2", 1297.712883, 402.332184, 107.114315, 770.278855),
    ("TC2", 2346.074070, 937.071408, 147.207932, 1117.605940),
    ("PA", 1859.922595, 1083.461569, 9.330267, 453.049877),
    ("T3", 314.527150, 1010.683521, 847.776296, 326.460239),
    ("TC3", 453.832806, 1749.112034, 1520.333737, 533.925755),
    ("PX", 695.709196, 1549.499821, 872.766651, 146.663022),
    ("PiOB", 1346.131617, 273.442572, 243.058298, 1084.591441),
    ("Pi2", 770.709556, 6.411500, 1050.380655, 1836.484355),
    ("Pi3", 0.202775, 765.959130, 1641.705152, 1240.045868),
    ("Pi4", 984.423828, 1581.886825, 609.106911, 28.165456),
    ("Pi5", 1798.961262, 1959.580763, 401.703009, 54.779052),
    ("Ti1", 1488.127317, 1364.050559, 158.778552, 86.197372),
    ("Ti2", 1262.278436, 387.442615, 112.010570, 762.885771),
    ("Ti3", 313.945074, 982.868588, 822.955503, 322.779988),
    ("PiX", 694.398776, 1512.204560, 834.546560, 135.567254),
    ("Ai2", 1758.286706, 1027.918057, 4.851343, 422.107929),
    ("IW1_SW", 200.555556, 238.472222, 796.631081, 893.981897),
    ("IW1_SE", 1350.694444, 1304.111111, 180.786119, 63.065230),
    ("IW1_NE", 1363.611111, 1289.694444, 167.964921, 71.228560),
    ("IW1_NW", 213.472222, 224.055556, 783.809883, 902.145226),
    ("IW2_SW", 213.472222, 224.055556, 783.809883, 902.145226),
    ("IW2_SE", 220.055556, 229.138889, 759.044390, 873.561893),
    ("IW2_NE", 324.472222, 142.222222, 683.296004, 944.705200),
    ("IW2_NW", 317.888889, 137.138889, 708.061496, 973.288534),
    ("IW2o_SW", 326.673889, 140.627222, 682.521346, 946.944733),
    ("IW2o_SE", 315.687222, 138.733889, 708.836154, 971.049001),
    ("IW2o_NE", 451.020556, 58.900556, 711.413829, 1150.411193),
    ("IW2o_NW", 462.007222, 60.793889, 685.099021, 1126.306926),
    ("IW2s_SW", 453.805556, 57.888889, 711.222504, 1153.234060),
    ("IW2s_SE", 459.222222, 61.805556, 685.290345, 1123.484060),
    ("IW2s_NE", 743.222222, 17.805556, 660.435968, 1350.444015),
    ("IW2s_NW", 737.805556, 13.888889, 686.368127, 1380.194015),
    ("IW3_SW", 69.888889, 757.805556, 1158.277866, 736.728653),
    ("IW3_SE", 75.555556, 762.472222, 1143.045316, 718.950875),
    ("IW3_NE", 124.111111, 464.805556, 865.587917, 707.297494),
    ("IW3_NW", 118.444444, 460.138889, 880.820468, 725.075272),
    ("IW4_SW", 693.888889, 1327.805556, 653.720139, 94.728653),
    ("IW4_SE", 711.555556, 1344.472222, 650.487589, 88.950875),
    ("IW4_NE", 755.555556, 1060.472222, 385.633211, 75.910830),
    ("IW4_NW", 737.888889, 1043.805556, 388.865762, 81.688608),
    ("IW5_SW", 671.388889, 793.472222, 298.944837, 157.748582),
    ("IW5_SE", 1293.611111, 1383.694444, 252.392110, 29.748582),
    ("IW5_NE", 1298.756944, 1375.173611, 244.669011, 32.517747),
    ("IW5_NW", 676.534722, 784.951389, 291.221738, 160.517747),
    ("IW6_SW", 444.534722, 41.840278, 1014.890925, 1515.113223),
    ("IW6_SE", 471.201389, 52.506944, 707.170117, 1166.668779),
    ("IW6_NE", 474.722222, 51.472222, 706.401306, 1169.397389),
    ("IW6_NW", 448.055556, 40.805556, 1014.122114, 1517.841834),
    ("IW7_SW", 119.555556, 478.472222, 878.190938, 705.910830),
    ("IW7_SE", 169.138889, 520.555556, 771.030143, 579.660830),
    ("IW7_NE", 173.694444, 506.888889, 758.427122, 581.047494),
    ("IW7_NW", 124.111111, 464.805556, 865.587917, 707.297494),
    ("IW8_SW", 180.034722, 196.451389, 1145.038185, 1323.497725),
    ("IW8_SE", 220.118056, 217.034722, 777.586784, 906.414391),
    ("IW8_NE", 233.784722, 203.368056, 765.515586, 915.327721),
    ("IW8_NW", 193.701389, 182.784722, 1132.966987, 1332.411054),
    ("IW9_SW", 125.138889, 804.555556, 1035.884520, 592.700875),
    ("IW9_SE", 132.694444, 811.111111, 1022.540859, 576.811986),
    ("IW9_NE", 292.694444, 315.111111, 564.832104, 622.731897),
    ("IW9_NW", 285.138889, 308.555556, 578.175765, 638.620786),
    ("IW11_SW", 552.694444, 1195.111111, 690.169041, 152.811986),
    ("IW11_SE", 568.472222, 1209.888889, 685.047601, 145.145319),
    ("IW11_NE", 728.472222, 713.888889, 227.338847, 191.065230),
    ("IW11_NW", 712.694444, 699.111111, 232.460286, 198.731897),
    ("IW12_SW", 608.138889, 939.777778, 433.018467, 130.940833),
    ("IW12_SE", 733.555556, 1057.694444, 401.691005, 80.524166),
    ("IW12_NE", 737.888889, 1043.805556, 388.865762, 81.688608),
    ("IW12_NW", 612.472222, 925.888889, 420.193224, 132.105275),
    ("dryer_SW", 1.111111, 694.805556, 1573.922724, 1227.448650),
    ("dryer_SE", 11.562500, 696.506944, 1401.505963, 1032.761150),
    ("dryer_NE", 22.812500, 571.090278, 1284.066639, 1020.244465),
    ("dryer_NW", 12.361111, 569.388889, 1456.483400, 1214.931965),
    ("washer_SW", 12.951389, 565.423611, 1452.784034, 1214.730020),
    ("washer_SE", 23.402778, 567.125000, 1280.367272, 1020.042520),
    ("washer_NE", 47.569444, 454.625000, 1175.844615, 1020.442501),
    ("washer_NW", 37.118056, 452.923611, 1348.261376, 1215.130001),
    ("ctr_poly_0", 83.506944, 449.673611, 993.661652, 839.834441),
    ("ctr_poly_1", 113.888889, 473.805556, 893.423489, 723.688608),
    ("ctr_poly_2", 69.888889, 757.805556, 1158.277866, 736.728653),
    ("ctr_poly_3", 39.506944, 733.673611, 1258.516029, 852.874486),
    ("bed_SW", 208.534722, 869.090278, 906.933698, 445.928929),
    ("bed_SE", 431.256944, 1072.812500, 732.570794, 223.206707),
    ("bed_NE", 505.673611, 719.006944, 403.760912, 223.154426),
    ("bed_NW", 282.951389, 515.284722, 578.123816, 445.876649),
    ("dresser_SW", 541.250000, 627.777778, 334.901679, 254.131909),
    ("dresser_SE", 665.444444, 743.472222, 281.452777, 179.048576),
    ("dresser_NE", 702.784722, 694.256944, 237.289816, 201.336620),
    ("dresser_NW", 578.590278, 578.562500, 290.738718, 276.419953),
    ("shelves_SW", 169.133030, 278.261502, 832.887345, 875.670491),
    ("shelves_SE", 253.383355, 348.043077, 614.752122, 620.709011),
    ("shelves_NE", 284.146701, 308.764757, 579.562190, 639.291932),
    ("shelves_NW", 199.896376, 238.983181, 797.697413, 894.253412),
    ("RO1_SW", 432.694444, 440.111111, 391.343846, 411.065230),
    ("RO1_SE", 547.222222, 545.138889, 307.329061, 302.870786),
    ("RO1_NE", 560.138889, 530.722222, 294.507863, 311.034115),
    ("RO1_NW", 445.611111, 425.694444, 378.522648, 419.228560),
    ("RO2_SW", 614.722222, 919.027778, 413.863936, 132.770829),
    ("RO2_SE", 598.944444, 904.250000, 418.985375, 140.437496),
    ("RO2_NE", 652.250000, 784.444444, 309.284454, 163.638583),
    ("RO2_NW", 668.027778, 799.222222, 304.163014, 155.971917),
    ("RO3_SW", 187.256944, 496.673611, 729.642184, 567.204436),
    ("RO3_SE", 179.701389, 490.118056, 742.985846, 583.093325),
    ("RO3_NE", 236.701389, 374.006944, 636.979369, 609.988857),
    ("RO3_NW", 244.256944, 380.562500, 623.635707, 594.099968),
    ("RO4_SW", 341.080556, 129.213889, 680.997333, 966.634863),
    ("RO4_SE", 330.093889, 127.320556, 707.312142, 990.739131),
    ("RO4_NE", 432.947222, 66.647222, 709.271175, 1127.054397),
    ("RO4_NW", 443.933889, 68.540556, 682.956366, 1102.950129),
    ("RO5_SW", 447.569444, 39.125000, 897.755149, 1383.342390),
    ("RO5_SE", 468.680556, 50.736111, 720.323697, 1181.731279),
    ("RO5_NE", 472.201389, 49.701389, 719.554886, 1184.459889),
    ("RO5_NW", 451.090278, 38.090278, 896.986338, 1386.071000),
    ("RO6_SW", 570.034722, 1170.451389, 648.003304, 139.577814),
    ("RO6_SE", 554.256944, 1155.673611, 653.124744, 147.244480),
    ("RO6_NE", 583.423611, 957.062500, 467.809204, 136.800005),
    ("RO6_NW", 599.201389, 971.840278, 462.687764, 129.133338),
    ("RO7_SW", 134.756944, 763.062500, 977.417273, 570.160035),
    ("RO7_SE", 127.201389, 756.506944, 990.760935, 586.048924),
    ("RO7_NE", 157.756944, 559.284722, 806.834284, 576.993337),
    ("RO7_NW", 165.312500, 565.840278, 793.490622, 561.104448),
    ("O1_outer_start", 68.071181, 371.210069, 1342.975836, 1307.936999),
    ("O1_outer_end", 96.571181, 313.154514, 1289.972597, 1321.384765),
    ("O1_inner_end", 95.904514, 310.487847, 1247.507496, 1273.829210),
    ("O1_inner_start", 67.404514, 368.543403, 1300.510735, 1260.381444),
    ("O2_outer_start", 252.930604, 136.578752, 1132.830529, 1419.845903),
    ("O2_outer_end", 305.730372, 102.822965, 1104.127059, 1457.593437),
    ("O2_inner_end", 305.063706, 100.156298, 1061.661958, 1410.037882),
    ("O2_inner_start", 252.263937, 133.912085, 1090.365428, 1372.290347),
    ("O3_outer_start", 484.694444, 33.888889, 1049.662996, 1593.447383),
    ("O3_outer_end", 609.138889, 12.555556, 1036.838828, 1692.540696),
    ("O3_inner_end", 608.472222, 9.888889, 994.373727, 1644.985141),
    ("O3_inner_start", 484.027778, 31.222222, 1007.197895, 1545.891827),
    ("O4_outer_start", 724.862847, 5.321181, 764.859494, 1469.820751),
    ("O4_outer_end", 730.987847, 9.196181, 723.961255, 1423.195751),
    ("O4_inner_end", 766.987847, 8.751736, 725.644102, 1452.857968),
    ("O4_inner_start", 760.862847, 4.876736, 766.542341, 1499.482968),
    ("O5_outer_start", 684.027778, 82.277778, 458.241024, 1010.238478),
    ("O5_outer_end", 829.472222, 210.722222, 248.398777, 757.127367),
    ("O5_inner_end", 862.138889, 206.944444, 246.748291, 783.456251),
    ("O5_inner_start", 716.694444, 78.500000, 456.590538, 1036.567362),
    ("O6_outer_start", 1004.277778, 371.027778, 120.165684, 591.988478),
    ("O6_outer_end", 1168.055556, 523.805556, 54.052073, 497.877367),
    ("O6_inner_end", 1200.722222, 520.027778, 52.401586, 524.206251),
    ("O6_inner_start", 1036.944444, 367.250000, 118.515197, 618.317362),
    ("O7_outer_start", 1591.171706, 941.308406, 4.000000, 374.466401),
    ("O7_outer_end", 1477.127146, 1132.739868, 64.000000, 181.086661),
    ("O7_inner_end", 1427.066078, 1096.100911, 64.444444, 177.512481),
    ("O7_inner_start", 1541.110637, 904.669449, 4.444444, 370.892220),
    ("O8_outer_start", 1278.966138, 1567.402703, 384.607156, 3.056663),
    ("O8_outer_end", 1307.562812, 1491.221599, 314.542015, 13.432211),
    ("O8_inner_end", 1261.118367, 1446.777155, 309.896005, 13.876656),
    ("O8_inner_start", 1232.521694, 1522.958259, 379.961146, 3.501107),
    ("O9_outer_start", 529.000000, 1209.361111, 732.140575, 170.621991),
    ("O9_outer_end", 444.506944, 1130.618056, 767.818018, 220.934491),
    ("O9_inner_end", 444.951389, 1094.618056, 733.945310, 215.041153),
    ("O9_inner_start", 529.444444, 1173.361111, 698.267867, 164.728653),
    ("O10_outer_start", 193.673611, 901.284722, 966.317854, 474.156713),
    ("O10_outer_end", 144.000000, 857.361111, 1036.814742, 559.288658),
    ("O10_inner_end", 144.444444, 821.361111, 1002.942033, 553.395319),
    ("O10_inner_start", 194.118056, 865.284722, 932.445145, 468.263375),
    ("O11_outer_start", 36.000000, 767.361111, 1305.000650, 873.288658),
    ("O11_outer_end", 19.506944, 755.618056, 1387.778876, 968.156713),
    ("O11_inner_end", 19.951389, 719.618056, 1353.906168, 962.263375),
    ("O11_inner_start", 36.444444, 731.361111, 1271.127942, 867.395319),
    ("O3_door_tip", 608.340278, 7.812500, 862.686516, 1493.256530),
    ("O6_door_tip", 1021.006944, 550.118056, 75.234743, 386.794335),
    ("RO1_door_tip", 481.847222, 630.347222, 388.737122, 269.500251),
    ("RO2_door_tip", 809.451389, 936.701389, 272.232924, 98.909973),
    ("RO3_door_tip", 125.694444, 440.055556, 862.332603, 727.571935),
    ("RO4_door_tip", 275.970278, 126.241389, 870.780467, 1141.312637),
    ("RO5_door_tip", 605.376736, 21.404514, 705.672671, 1295.372437),
    ("RO6_door_tip_S", 471.694444, 1074.055556, 681.218294, 192.834202),
    ("RO6_door_tip_N", 499.694444, 883.388889, 503.315376, 182.807506),
    ("RO7_door_tip_S", 180.555556, 798.805556, 903.980546, 486.221978),
    ("RO7_door_tip_N", 209.888889, 609.472222, 727.410961, 477.528615),
    ("w_f8f9_start", 694.444444, 52.138889, 530.257335, 1129.955138),
    ("w_f8f9_arc_entry", 663.153307, 55.498670, 531.594515, 1104.689885),
    ("w_f8f9_arc_center", 666.517682, 58.264620, 522.452272, 1094.024514),
    ("w_f8f9_arc_exit", 656.836597, 59.488172, 523.039311, 1086.239773),
    ("w_f8f9_end", 668.055556, 68.805556, 494.516845, 1052.877367),
    ("g_f8f9_start", 697.250000, 54.444444, 522.613282, 1121.038472),
    ("g_f8f9_arc_entry", 665.958863, 57.804226, 523.950462, 1095.773219),
    ("g_f8f9_arc_center", 666.517682, 58.264620, 522.452272, 1094.024514),
    ("g_f8f9_arc_exit", 664.919930, 58.460394, 522.543356, 1092.738661),
    ("g_f8f9_end", 676.138889, 67.777778, 494.020890, 1059.376255),
    ("wh_center", 695.486893, 31.911223, 601.080369, 1232.016098),
    ("sink_NW", 724.055556, 116.805556, 383.323108, 921.321811),
    ("sink_SE", 734.201389, 225.034722, 261.304207, 686.730993),
    ("stove_NW", 651.513889, 66.680556, 506.518524, 1056.295702),
    ("stove_SE", 599.784722, 130.451389, 413.865620, 850.123219),
    ("dw_NW", 834.722222, 215.472222, 243.199168, 750.655145),
    ("dw_SE", 810.368056, 307.118056, 188.370873, 583.579884),
    ("fridge_NW", 316.951389, 152.173611, 677.259334, 922.769925),
    ("fridge_SE", 276.292101, 262.771267, 615.625250, 717.887047),
    ("ice_NW", 920.451389, 306.368056, 164.153921, 639.483202),
    ("ice_SE", 916.446736, 369.916181, 135.446610, 544.254003),
    ("nctr_NW", 616.138889, 27.388889, 662.041958, 1247.127367),
    ("nctr_SE", 544.555556, 83.472222, 531.054995, 988.310719),
    ("wctr_NW", 316.365017, 233.760851, 575.679092, 728.493288),
    ("wctr_SE", 400.240017, 384.635851, 423.279428, 474.961633),
    ("loveseat_NW", 828.244635, 485.746011, 139.443784, 375.023217),
    ("et_center", 939.612646, 833.677067, 151.493291, 151.992920),
    ("loveseat2_NW", 1057.992408, 834.489901, 100.081568, 181.032028),
    ("loveseat2_SE", 1306.296217, 1225.988154, 158.910642, 76.091415),
    ("chair_center", 1410.024070, 731.441066, 10.061010, 438.142055),
    ("ottoman_center", 1195.223360, 679.910586, 36.022036, 353.201759),
    ("desk_SW", 900.159074, 1511.131807, 619.476396, 41.760705),
    ("dim01_A", 517.246679, 491.197557, 321.994696, 347.092114),
    ("dim01_B", 927.913345, 300.530890, 166.428337, 653.185366),
    ("dim02_A", 452.730255, 80.039977, 632.229910, 1039.845179),
    ("dim02_B", 1462.845462, 1013.129357, 32.565639, 244.132672),
    ("dim03_A", 669.284722, 997.173611, 415.792236, 104.169999),
    ("dim03_B", 629.618056, 1267.284722, 667.821370, 118.374486),
    ("dim04_A", 142.784722, 497.951389, 823.048041, 641.223330),
    ("dim04_B", 98.784722, 781.951389, 1087.902418, 654.263375),
    ("dim05_A", 14.722222, 555.027778, 1473.034760, 1249.290851),
    ("dim05_B", 84.138889, 598.944444, 1009.521389, 720.874185),
    ("dim06_A", 766.762739, 1028.310415, 356.002887, 79.576034),
    ("dim06_B", 1245.901628, 1483.949304, 344.010172, 7.770479),
    ("dim07_A", 701.237847, 748.154514, 258.014667, 174.525863),
    ("dim07_B", 1323.460069, 1338.376736, 211.461940, 46.525863),
    ("dim08_A", 81.027778, 338.888889, 1273.382380, 1266.478591),
    ("dim08_B", 205.694444, 429.555556, 687.475663, 594.034146),
    ("dim09_A", 278.037085, 116.407456, 1075.386957, 1390.537378),
    ("dim09_B", 318.120419, 136.990789, 707.935555, 973.454045),
    ("dim10_A", 625.027778, 8.222222, 993.770706, 1658.371805),
    ("dim10_B", 651.694444, 18.888889, 686.049898, 1309.927361),
    ("dim11_A", 962.846994, 1068.479359, 246.590569, 63.522859),
    ("dim11_B", 859.930328, 1484.896025, 632.693333, 50.919596),
    ("dim12a_A", 454.277778, 39.027778, 853.150599, 1336.508500),
    ("dim12a_B", 717.361111, 1.444444, 833.117420, 1547.305126),
    ("dim12b_A", 199.923611, 181.006944, 971.995472, 1151.077721),
    ("dim12b_B", 450.756944, 40.062500, 853.919410, 1333.779890),
    ("dim13_A", 12.484444, 772.445556, 1879.202706, 1520.284213),
    ("dim13_B", 759.595556, 25.334444, 1219.310542, 2007.546232),
    ("dim14_A", 864.691339, 41.558010, 558.061351, 1301.205195),
    ("dim14_B", 1528.184788, 647.856867, 26.232610, 623.803218),
    ("dim15_A", 9.694444, 925.555556, 1861.556236, 1356.114235),
    ("dim15_B", 1245.694444, 2053.555556, 840.440782, 60.114235),
    ("dim17_A", 168.362847, 842.404514, 966.775186, 509.910944),
    ("dim17_B", 328.362847, 346.404514, 509.066431, 555.830855),
    ("dim18_A", 176.694444, 527.111111, 757.686481, 563.771941),
    ("dim18_B", 596.694444, 911.111111, 425.314663, 139.771941),
    ("dim19_A", 27.571181, 724.862847, 1311.890319, 914.202611),
    ("dim19_B", 207.133681, 207.425347, 835.137267, 972.555016),
    ("dim22_A", 673.618056, 983.284722, 402.966993, 105.334441),
    ("dim22_B", 732.534722, 850.868056, 281.718606, 130.977749),
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
