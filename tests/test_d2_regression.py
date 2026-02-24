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
    ("F1", 0.000000, 749.361111, 1638.423559, 1255.259576),
    ("F2", 1.388889, 707.694444, 1648.983961, 1308.078364),
    ("F5", 625.694444, 10.888889, 1029.292719, 1724.128763),
    ("F6", 749.361111, 0.000000, 887.059077, 1663.576924),
    ("F7", 791.555556, 26.694444, 606.432204, 1342.382480),
    ("F8", 706.000000, 61.694444, 495.196412, 1112.823207),
    ("F9", 700.722222, 65.027778, 488.084332, 1097.259450),
    ("F10", 1230.203495, 547.701706, 42.700532, 529.368534),
    ("F11", 1319.825904, 592.186487, 29.910639, 541.719223),
    ("F12", 1638.423559, 887.059077, 0.000000, 474.620574),
    ("F13", 1432.913635, 1298.674405, 146.924724, 98.611710),
    ("F14", 1423.472222, 1320.222222, 161.126468, 87.571028),
    ("F15", 1255.259576, 1663.576924, 474.620574, 0.000000),
    ("F16", 1181.425744, 1678.086260, 537.924564, 3.271495),
    ("F17", 895.993878, 1542.311443, 651.735050, 43.935533),
    ("F18", 842.613885, 1504.891560, 667.334877, 56.250665),
    ("F11a", 1532.633833, 703.691187, 12.872709, 581.104774),
    ("F11b", 1589.688190, 757.745544, 7.445520, 567.825797),
    ("FC", 494.361111, 452.722222, 334.021688, 390.729492),
    ("P3", 0.694444, 752.555556, 1691.185958, 1314.565132),
    ("POB", 1417.872460, 312.896625, 219.414001, 1103.467484),
    ("P2", 841.694444, 8.222222, 1039.389800, 1905.658944),
    ("P4", 961.650783, 1617.980410, 662.490712, 35.838897),
    ("P5", 1820.080898, 2045.939326, 457.680479, 55.987438),
    ("T1", 1498.618268, 1402.474612, 175.884436, 85.284272),
    ("TC1", 2314.261138, 1824.647647, 183.337876, 326.352629),
    ("T2", 1316.080277, 430.001990, 90.656407, 763.228438),
    ("TC2", 2366.883865, 988.968283, 131.223207, 1095.813884),
    ("PA", 1856.242937, 1122.428701, 13.925903, 437.086354),
    ("T3", 291.655033, 989.782412, 863.024343, 345.768225),
    ("TC3", 412.655033, 1712.115745, 1550.257368, 561.560227),
    ("PX", 656.083155, 1533.759161, 903.636548, 161.311931),
    ("PiOB", 1375.815362, 301.233108, 216.445795, 1078.461632),
    ("Pi2", 817.617319, 4.945518, 1005.070630, 1849.930509),
    ("Pi3", 0.361111, 723.388889, 1633.873987, 1274.756404),
    ("Pi4", 944.833994, 1579.167718, 641.378116, 34.778780),
    ("Pi5", 1757.496093, 1980.578004, 438.270703, 46.586333),
    ("Ti1", 1463.086125, 1386.615960, 180.761764, 78.480854),
    ("Ti2", 1280.548133, 414.143338, 95.533735, 756.425020),
    ("Ti3", 291.905033, 962.699079, 837.536478, 341.709498),
    ("PiX", 655.754857, 1497.522986, 864.650627, 149.648709),
    ("Ai2", 1754.548362, 1064.584740, 9.246400, 407.559412),
    ("IW1_SW", 187.138889, 228.555556, 850.715223, 969.855149),
    ("IW1_SE", 1350.694444, 1304.111111, 181.478783, 70.299594),
    ("IW1_NE", 1363.611111, 1289.694444, 168.657585, 78.907533),
    ("IW1_NW", 200.055556, 214.138889, 837.894025, 978.463089),
    ("IW2_SW", 213.472222, 224.055556, 778.415555, 909.824200),
    ("IW2_SE", 220.055556, 229.138889, 753.758116, 881.240866),
    ("IW2_NE", 324.472222, 142.222222, 678.009729, 955.496441),
    ("IW2_NW", 317.888889, 137.138889, 702.667168, 984.079775),
    ("IW2o_SW", 326.673889, 140.627222, 677.232910, 957.798219),
    ("IW2o_SE", 315.687222, 138.733889, 703.443987, 981.777996),
    ("IW2o_NE", 451.020556, 58.900556, 705.769538, 1164.697066),
    ("IW2o_NW", 462.007222, 60.793889, 679.558461, 1140.717289),
    ("IW2s_SW", 453.805556, 57.888889, 705.576053, 1167.582178),
    ("IW2s_SE", 459.222222, 61.805556, 679.751947, 1137.832178),
    ("IW2s_NE", 743.222222, 17.805556, 654.897569, 1370.127449),
    ("IW2s_NW", 737.805556, 13.888889, 680.721675, 1399.877449),
    ("IW3_SW", 69.888889, 757.805556, 1153.315751, 733.292384),
    ("IW3_SE", 75.555556, 762.472222, 1138.155236, 715.514606),
    ("IW3_NE", 128.888889, 451.361111, 848.317039, 711.398241),
    ("IW3_NW", 123.222222, 446.694444, 863.477554, 729.176019),
    ("IW4_SW", 693.888889, 1327.805556, 652.647936, 91.292384),
    ("IW4_SE", 711.555556, 1344.472222, 649.487421, 85.514606),
    ("IW4_NE", 760.111111, 1046.805556, 372.030022, 79.492948),
    ("IW4_NW", 742.444444, 1030.138889, 375.190537, 85.270726),
    ("IW5_SW", 666.368056, 802.118056, 305.180467, 159.893428),
    ("IW5_SE", 1288.590278, 1392.340278, 260.932873, 31.893428),
    ("IW5_NE", 1293.611111, 1383.694444, 253.084774, 34.759897),
    ("IW5_NW", 671.388889, 793.472222, 297.332368, 162.759897),
    ("IW6_SW", 444.534722, 41.840278, 1008.091907, 1529.831850),
    ("IW6_SE", 471.201389, 52.506944, 701.523665, 1181.387405),
    ("IW6_NE", 474.722222, 51.472222, 700.754855, 1184.190117),
    ("IW6_NW", 448.055556, 40.805556, 1007.323096, 1532.634562),
    ("IW7_SW", 124.111111, 464.805556, 860.697837, 709.492948),
    ("IW7_SE", 173.694444, 506.888889, 754.077307, 583.242948),
    ("IW7_NE", 178.472222, 493.444444, 741.696509, 585.148241),
    ("IW7_NW", 128.888889, 451.361111, 848.317039, 711.398241),
    ("IW8_SW", 160.472222, 217.888889, 1157.283464, 1318.299594),
    ("IW8_SE", 187.138889, 228.555556, 850.715223, 969.855149),
    ("IW8_NE", 200.055556, 214.138889, 837.894025, 978.463089),
    ("IW8_NW", 173.388889, 203.472222, 1144.462266, 1326.907533),
    ("IW9_SW", 125.138889, 804.555556, 1031.534706, 589.264606),
    ("IW9_SE", 132.694444, 811.111111, 1018.263080, 573.375717),
    ("IW9_NE", 186.027778, 500.000000, 728.424883, 569.259352),
    ("IW9_NW", 178.472222, 493.444444, 741.696509, 585.148241),
    ("IW11_SW", 552.694444, 1195.111111, 688.484536, 149.375717),
    ("IW11_SE", 568.472222, 1209.888889, 683.435132, 141.709051),
    ("IW11_NE", 728.472222, 713.888889, 225.726377, 198.299594),
    ("IW11_NW", 712.694444, 699.111111, 230.775781, 205.966260),
    ("IW12_SW", 612.472222, 925.888889, 418.580755, 134.004322),
    ("IW12_SE", 737.888889, 1043.805556, 387.793558, 83.587655),
    ("IW12_NE", 742.444444, 1030.138889, 375.190537, 85.270726),
    ("IW12_NW", 617.027778, 912.222222, 405.977734, 135.687393),
    ("IW16_SW", 123.222222, 446.694444, 863.477554, 729.176019),
    ("IW16_SE", 128.888889, 451.361111, 848.317039, 711.398241),
    ("IW16_NE", 235.555556, 266.472222, 680.446481, 772.105149),
    ("IW16_NW", 229.888889, 261.805556, 695.606996, 789.882927),
    ("dryer_SW", 1.111111, 694.805556, 1567.231759, 1224.308788),
    ("dryer_SE", 11.562500, 696.506944, 1395.445307, 1029.621288),
    ("dryer_NE", 22.812500, 571.090278, 1278.005983, 1019.327651),
    ("dryer_NW", 12.361111, 569.388889, 1449.792435, 1214.015151),
    ("washer_SW", 12.951389, 565.423611, 1446.093069, 1213.887308),
    ("washer_SE", 23.402778, 567.125000, 1274.306617, 1019.199808),
    ("washer_NE", 47.569444, 454.625000, 1169.783960, 1021.822837),
    ("washer_NW", 37.118056, 452.923611, 1341.570411, 1216.510337),
    ("ctr_poly_0", 83.506944, 449.673611, 988.249315, 841.733489),
    ("ctr_poly_1", 113.888889, 473.805556, 888.461373, 725.587655),
    ("ctr_poly_2", 69.888889, 757.805556, 1153.315751, 733.292384),
    ("ctr_poly_3", 39.506944, 733.673611, 1253.103692, 849.438217),
    ("bed_SW", 208.534722, 869.090278, 903.286229, 442.640864),
    ("bed_SE", 431.256944, 1072.812500, 730.291998, 219.918642),
    ("bed_NE", 505.673611, 719.006944, 401.482116, 226.831913),
    ("bed_NW", 282.951389, 515.284722, 574.476347, 449.554135),
    ("dresser_SW", 541.250000, 627.777778, 332.568855, 259.884240),
    ("dresser_SE", 665.444444, 743.472222, 279.732255, 184.800907),
    ("dresser_NE", 702.784722, 694.256944, 235.569294, 208.496882),
    ("dresser_NW", 578.590278, 578.562500, 288.405895, 283.580215),
    ("RO1_SW", 432.694444, 440.111111, 388.146598, 418.299594),
    ("RO1_SE", 547.222222, 545.138889, 304.816149, 310.105149),
    ("RO1_NE", 560.138889, 530.722222, 291.994951, 318.713089),
    ("RO1_NW", 445.611111, 425.694444, 375.325400, 426.907533),
    ("RO2_SW", 620.590278, 902.118056, 396.671302, 137.095529),
    ("RO2_SE", 604.812500, 887.340278, 401.720705, 144.762196),
    ("RO2_NE", 660.756944, 770.173611, 294.658673, 173.418034),
    ("RO2_NW", 676.534722, 784.951389, 289.609269, 165.751367),
    ("RO3_SW", 140.284722, 404.534722, 824.774680, 736.902943),
    ("RO3_SE", 145.951389, 409.201389, 809.614165, 719.125165),
    ("RO3_NE", 209.284722, 299.423611, 709.941021, 755.169892),
    ("RO3_NW", 203.618056, 294.756944, 725.101536, 772.947670),
    ("RO4_SW", 341.080556, 129.213889, 675.678643, 977.915175),
    ("RO4_SE", 330.093889, 127.320556, 701.889720, 1001.894952),
    ("RO4_NE", 432.947222, 66.647222, 703.657139, 1140.913445),
    ("RO4_NW", 443.933889, 68.540556, 677.446062, 1116.933668),
    ("RO5_SW", 447.569444, 39.125000, 891.370334, 1398.061016),
    ("RO5_SE", 468.680556, 50.736111, 714.623218, 1196.449905),
    ("RO5_NE", 472.201389, 49.701389, 713.854408, 1199.252617),
    ("RO5_NW", 451.090278, 38.090278, 890.601523, 1400.863728),
    ("RO6_SW", 569.201389, 1187.840278, 662.716078, 138.847611),
    ("RO6_SE", 553.423611, 1173.062500, 667.765482, 146.514278),
    ("RO6_NE", 591.312500, 928.506944, 439.696435, 139.879651),
    ("RO6_NW", 607.090278, 943.284722, 434.647031, 132.212984),
    ("RO7_SW", 133.812500, 780.340278, 989.353626, 569.466924),
    ("RO7_SE", 126.256944, 773.784722, 1002.625253, 585.355813),
    ("RO7_NE", 165.868056, 530.951389, 776.278428, 580.443408),
    ("RO7_NW", 173.423611, 537.506944, 763.006801, 564.554519),
    ("O1_outer_start", 85.487847, 333.960069, 1301.973684, 1320.006139),
    ("O1_outer_end", 117.154514, 279.071181, 1252.137112, 1338.028502),
    ("O1_inner_end", 116.487847, 276.404514, 1209.816082, 1290.472946),
    ("O1_inner_start", 84.821181, 331.293403, 1259.652654, 1272.450583),
    ("O2_outer_start", 252.930604, 136.578752, 1125.887440, 1429.939353),
    ("O2_outer_end", 305.730372, 102.822965, 1097.183970, 1469.094819),
    ("O2_inner_end", 305.063706, 100.156298, 1054.862940, 1421.539263),
    ("O2_inner_start", 252.263937, 133.912085, 1083.566410, 1382.383798),
    ("O3_outer_start", 484.694444, 33.888889, 1042.719907, 1608.981127),
    ("O3_outer_end", 609.138889, 12.555556, 1029.895740, 1710.445692),
    ("O3_inner_end", 608.472222, 9.888889, 987.574709, 1662.890137),
    ("O3_inner_start", 484.027778, 31.222222, 1000.398877, 1561.425572),
    ("O4_outer_start", 724.862847, 5.321181, 758.897887, 1489.504186),
    ("O4_outer_end", 730.987847, 9.196181, 718.161728, 1442.879186),
    ("O4_inner_end", 766.987847, 8.751736, 719.844575, 1473.134216),
    ("O4_inner_start", 760.862847, 4.876736, 760.580734, 1519.759216),
    ("O5_outer_start", 684.027778, 82.277778, 453.639086, 1027.698864),
    ("O5_outer_end", 829.472222, 210.722222, 245.021440, 774.587753),
    ("O5_inner_end", 862.138889, 206.944444, 243.370954, 801.509450),
    ("O5_inner_start", 716.694444, 78.500000, 451.988600, 1054.620561),
    ("O6_outer_start", 1008.488649, 374.931304, 115.642611, 606.454213),
    ("O6_outer_end", 1173.017713, 528.460369, 51.072676, 513.094389),
    ("O6_inner_end", 1205.684380, 524.682591, 49.422190, 540.016085),
    ("O6_inner_start", 1041.155315, 371.153527, 113.992124, 633.375910),
    ("O7_outer_start", 1584.272028, 934.732887, 4.000000, 392.336808),
    ("O7_outer_end", 1469.817436, 1125.754316, 64.000000, 193.485511),
    ("O7_inner_end", 1419.893045, 1089.252037, 64.444444, 189.860544),
    ("O7_inner_start", 1534.347637, 898.230608, 4.444444, 388.711842),
    ("O8_outer_start", 1274.570637, 1581.993200, 398.926159, 3.406373),
    ("O8_outer_end", 1301.835975, 1504.480761, 327.529683, 14.154924),
    ("O8_inner_end", 1255.391531, 1460.036316, 322.739602, 14.599368),
    ("O8_inner_start", 1228.126192, 1537.548755, 394.136078, 3.850818),
    ("O9_outer_start", 529.000000, 1209.361111, 730.348018, 166.592909),
    ("O9_outer_end", 444.506944, 1130.618056, 765.611257, 216.905409),
    ("O9_inner_end", 444.951389, 1094.618056, 731.738549, 211.604884),
    ("O9_inner_start", 529.444444, 1173.361111, 696.475309, 161.292384),
    ("O10_outer_start", 193.673611, 901.284722, 962.562332, 470.127632),
    ("O10_outer_end", 144.000000, 857.361111, 1032.645016, 555.259576),
    ("O10_inner_end", 144.444444, 821.361111, 998.772307, 549.959051),
    ("O10_inner_start", 194.118056, 865.284722, 928.689623, 464.827106),
    ("O11_outer_start", 36.000000, 767.361111, 1299.534287, 869.259576),
    ("O11_outer_end", 19.506944, 755.618056, 1381.970345, 964.127632),
    ("O11_inner_end", 19.951389, 719.618056, 1348.097637, 958.827106),
    ("O11_inner_start", 36.444444, 731.361111, 1265.661579, 863.959051),
    ("O3_door_tip", 608.340278, 7.812500, 856.355728, 1511.087425),
    ("O6_door_tip", 1025.952027, 554.755794, 72.220263, 399.178421),
    ("RO1_door_tip", 481.847222, 630.347222, 386.206201, 274.289261),
    ("RO2_door_tip", 817.888889, 922.361111, 258.222035, 108.545877),
    ("RO3_door_tip", 162.500000, 266.694444, 865.304170, 931.825513),
    ("RO4_door_tip", 510.334444, 89.231111, 541.919182, 990.067068),
    ("RO5_door_tip", 605.376736, 21.404514, 699.954184, 1312.795772),
    ("RO6_door_tip_S", 449.388889, 1071.472222, 705.231388, 205.872545),
    ("RO6_door_tip_N", 486.055556, 834.805556, 484.519407, 199.451938),
    ("RO7_door_tip_S", 193.138889, 828.111111, 898.600099, 463.964081),
    ("RO7_door_tip_N", 231.472222, 593.111111, 679.554784, 459.210140),
    ("w_f8f9_start", 694.444444, 52.138889, 525.295220, 1148.156541),
    ("w_f8f9_arc_entry", 663.153307, 55.498670, 526.632400, 1122.327648),
    ("w_f8f9_arc_center", 666.517682, 58.264620, 517.533264, 1111.662277),
    ("w_f8f9_arc_exit", 656.836597, 59.488172, 518.120304, 1103.700159),
    ("w_f8f9_end", 668.055556, 68.805556, 489.734818, 1070.337753),
    ("g_f8f9_start", 697.250000, 54.444444, 517.687185, 1139.239874),
    ("g_f8f9_arc_entry", 665.958863, 57.804226, 519.024364, 1113.410982),
    ("g_f8f9_arc_center", 666.517682, 58.264620, 517.533264, 1111.662277),
    ("g_f8f9_arc_exit", 664.919930, 58.460394, 517.624349, 1110.347250),
    ("g_f8f9_end", 676.138889, 67.777778, 489.238863, 1076.984844),
    ("wh_center", 695.486893, 31.911223, 595.794095, 1250.636682),
    ("sink_NW", 724.055556, 116.805556, 379.117364, 938.782197),
    ("sink_SE", 734.201389, 225.034722, 257.908862, 702.412940),
    ("stove_NW", 651.513889, 66.680556, 501.682470, 1073.533783),
    ("stove_SE", 599.784722, 130.451389, 409.569832, 865.360556),
    ("dw_NW", 834.722222, 215.472222, 239.857850, 768.115531),
    ("dw_SE", 810.368056, 307.118056, 185.533802, 599.039526),
    ("fridge_NW", 316.951389, 152.173611, 672.027086, 933.264760),
    ("fridge_SE", 276.292101, 262.771267, 610.982792, 725.788326),
    ("ice_NW", 920.451389, 306.368056, 161.424904, 656.721283),
    ("ice_SE", 916.446736, 369.916181, 133.036349, 560.321279),
    ("nctr_NW", 616.138889, 27.388889, 656.503560, 1264.587753),
    ("nctr_SE", 544.555556, 83.472222, 526.164915, 1003.548056),
    ("wctr_NW", 316.365017, 233.760851, 571.090660, 737.506090),
    ("wctr_SE", 400.240017, 384.635851, 419.771528, 482.640606),
    ("loveseat_NW", 828.244635, 485.746011, 137.394187, 387.690226),
    ("et_center", 939.612646, 833.677067, 150.606193, 160.549443),
    ("loveseat2_NW", 1057.992408, 834.489901, 99.407740, 191.452761),
    ("loveseat2_SE", 1306.296217, 1225.988154, 159.407390, 83.918591),
    ("chair_center", 1401.514306, 725.962623, 10.174896, 455.410197),
    ("ottoman_center", 1187.313350, 675.031897, 36.421579, 368.484795),
    ("desk_SW", 876.763350, 1492.519099, 623.612177, 43.174159),
    ("dim01_A", 514.502344, 488.675650, 321.127426, 357.241465),
    ("dim01_B", 925.169011, 298.008984, 165.561067, 673.116130),
    ("dim02_A", 446.224608, 82.394684, 628.781156, 1048.625227),
    ("dim02_B", 1452.830077, 1012.136407, 34.531397, 253.202851),
    ("dim03_A", 673.618056, 983.284722, 401.624657, 107.233489),
    ("dim03_B", 629.618056, 1267.284722, 666.479034, 114.938217),
    ("dim04_A", 147.340278, 484.284722, 805.825072, 644.805448),
    ("dim04_B", 98.784722, 781.951389, 1083.282471, 650.827106),
    ("dim05_A", 16.027778, 547.222222, 1458.962009, 1248.539757),
    ("dim05_B", 85.444444, 591.138889, 997.285541, 720.123090),
    ("dim06_A", 761.701570, 1042.235244, 367.819320, 80.112246),
    ("dim06_B", 1240.840459, 1497.874133, 357.519437, 8.306690),
    ("dim07_A", 698.368056, 752.118056, 259.966873, 178.967245),
    ("dim07_B", 1320.590278, 1342.340278, 215.719279, 50.967245),
    ("dim08_A", 100.027778, 303.222222, 1234.107632, 1280.835029),
    ("dim08_B", 169.444444, 347.138889, 772.431164, 752.418362),
    ("dim10_A", 625.027778, 8.222222, 986.971688, 1676.573207),
    ("dim10_B", 665.111111, 28.805556, 620.924977, 1259.489874),
    ("dim11_A", 1003.058330, 972.891560, 175.753414, 107.540682),
    ("dim11_B", 843.058330, 1468.891560, 633.462169, 50.950139),
    ("dim12a_A", 454.277778, 39.027778, 846.927864, 1351.301228),
    ("dim12a_B", 717.361111, 1.444444, 826.894685, 1566.988561),
    ("dim12b_A", 179.611111, 201.694444, 984.067034, 1145.574200),
    ("dim12b_B", 450.756944, 40.062500, 847.696675, 1348.498516),
    ("dim13_A", 12.484444, 772.445556, 1871.676130, 1516.255132),
    ("dim13_B", 759.595556, 25.334444, 1211.783966, 2027.822480),
    ("dim14_A", 864.691339, 41.558010, 552.998786, 1322.370662),
    ("dim14_B", 1522.375992, 642.372231, 26.232610, 646.759625),
    ("dim15_A", 9.694444, 925.555556, 1854.613147, 1349.417496),
    ("dim15_B", 1245.694444, 2053.555556, 841.277516, 53.417496),
    ("dim16_A", 486.279514, 1133.071181, 713.188526, 185.530231),
    ("dim16_B", 646.279514, 637.071181, 255.479771, 242.120774),
    ("dim17_A", 168.362847, 842.404514, 962.812562, 506.474676),
    ("dim17_B", 328.362847, 346.404514, 505.103807, 563.065219),
    ("dim18_A", 148.250000, 644.444444, 862.232870, 560.206424),
    ("dim18_B", 568.250000, 1028.444444, 532.454327, 136.206424),
    ("dim19_A", 27.571181, 724.862847, 1306.252872, 910.766342),
    ("dim19_B", 187.571181, 228.862847, 848.544117, 967.356885),
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
