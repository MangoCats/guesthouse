"""Distance-squared regression tests.

For every corner of every IW wall, every corner/center of every appliance
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
from floorplan.constants import (
    WALL_OUTER, SHELL_THICKNESS, AIR_GAP,
    F8F9_INNER_TURN_R, OPENING_INSIDE_RADIUS,
    O3_DOOR_WIDTH, O6_DOOR_WIDTH, O6_WIDTH,
    RO1_DOOR_WIDTH, RO2_DOOR_WIDTH, RO3_DOOR_WIDTH,
    RO4_DOOR_WIDTH, RO5_DOOR_WIDTH, RO6_DOOR_WIDTH, RO7_DOOR_WIDTH,
    IW4_RO_WIDTH, IW9_RO_WIDTH, IW11_RO_WIDTH,
)

TOL = 1e-6  # tolerance in ft^2

# ============================================================
# F-series outline chain definition
# ============================================================
# Starting point and initial bearing
_F2_E = 0.500000000000
_F2_N = 3.000000000000
_INITIAL_BRG = 0.0  # bearing 0 = north

# Sweep angle constants
_A19 = math.atan(1.0 / 9.0)   # arctan(1/9) for F3-F4, F19-F20
_A9 = math.atan(9.0)           # arctan(9) for F5-F6, F1-F2
_PI_2 = math.pi / 2            # 90 deg
_5PI_12 = 5 * math.pi / 12    # 75 deg
_PI_12 = math.pi / 12          # 15 deg
_PI_3 = math.pi / 3            # 60 deg
_PI_6 = math.pi / 6            # 30 deg

# Chain: ("L", distance) | ("CW", radius, sweep) | ("CCW", radius, sweep)
# Starting at F2, bearing north, traversing CW to produce:
# F3, F4, F5, F6, F7, F8, F9, F10, F11, F11a, F11b, F12, F13, F14, F15,
# F16, F17, F18, F19, F20, F1, and closing back to F2.
_CHAIN = [
    ("L",   12.083333333333),                       # F2->F3
    ("CW",   8.351795046046, _A19),                 # F3->F4
    ("L",    9.476667232982),                        # F4->F5
    ("CW",   2.333333333333, _A9),                  # F5->F6
    ("L",    5.250000000000),                        # F6->F7
    ("CW",   2.333333333333, _PI_2),                # F7->F8
    ("CCW",  0.166666666667, _PI_2),                 # F8->F9
    ("L",   15.166666666667),                        # F9->F10
    ("CCW",  1.039662132188, _5PI_12),               # F10->F11
    ("CW",   2.333333333333, _5PI_12),               # F11->F11a
    ("L",    1.000000000000),                        # F11a->F11b
    ("CW",   2.333333333333, _5PI_12),               # F11b->F12
    ("L",   11.858994000010),                        # F12->F13
    ("CW",   2.507553207938, _PI_12),                # F13->F14
    ("L",    8.666666666667),                        # F14->F15
    ("CW",   2.473295271375, _PI_3),                 # F15->F16
    ("L",    5.000000000000),                        # F16->F17
    ("CW",   6.404672887007, _PI_6),                 # F17->F18
    ("L",    1.397555568554),                        # F18->F19
    ("CW",  18.888718471469, _A19),                  # F19->F20
    ("L",   23.147693701700),                        # F20->F1
    ("CW",   0.833333333333, _A9),                   # F1->F2
]

# Point names produced by the chain (one per segment, in order)
_CHAIN_NAMES = [
    "F3", "F4", "F5", "F6", "F7", "F8", "F9", "F10", "F11",
    "F11a", "F11b", "F12", "F13", "F14", "F15", "F16", "F17",
    "F18", "F19", "F20", "F1", "F2",
]


def _reconstruct_f_points():
    """Walk the chain to reconstruct all F-series coordinates.

    Returns dict mapping point name to (E, N).
    """
    E, N = _F2_E, _F2_N
    brg = _INITIAL_BRG
    pts = {"F2": (E, N)}

    for seg, name in zip(_CHAIN, _CHAIN_NAMES):
        if seg[0] == "L":
            d = seg[1]
            E += d * math.sin(brg)
            N += d * math.cos(brg)
        elif seg[0] == "CW":
            R, sweep = seg[1], seg[2]
            cx = E + R * math.cos(brg)
            cy = N - R * math.sin(brg)
            alpha = math.atan2(N - cy, E - cx)
            alpha -= sweep
            E = cx + R * math.cos(alpha)
            N = cy + R * math.sin(alpha)
            brg += sweep
        elif seg[0] == "CCW":
            R, sweep = seg[1], seg[2]
            cx = E - R * math.cos(brg)
            cy = N + R * math.sin(brg)
            alpha = math.atan2(N - cy, E - cx)
            alpha += sweep
            E = cx + R * math.cos(alpha)
            N = cy + R * math.sin(alpha)
            brg -= sweep
        pts[name] = (E, N)

    return pts


from conftest import dist_sq as _dist_sq


def _bbox_corners(bbox):
    return [(bbox.w, bbox.s), (bbox.e, bbox.s), (bbox.e, bbox.n), (bbox.w, bbox.n)]


def _collect_all_points(pts, layout):
    """Collect all named points from layout, openings, and door tips.

    Returns list of (name, point) tuples.
    """
    outer_openings = compute_outer_openings(pts, layout)
    rough_openings = compute_rough_openings(pts, layout)
    result = []

    # ---- F-series outline points (from chain definition) ----
    chain_pts = _reconstruct_f_points()
    f_names = [f"F{i}" for i in range(1, 21)] + ["F11a", "F11b"]
    for name in f_names:
        result.append((name, chain_pts[name]))

    # ---- IW walls, appliances, furniture (all Wall type with .poly) ----
    for prefix, wall in [
        ("IW1", layout.iw1), ("IW2", layout.iw2),
        ("IW3", layout.iw3), ("IW4", layout.iw4),
        ("IW5", layout.iw5), ("IW6", layout.iw6),
        ("IW7", layout.iw7), ("IW8", layout.iw8),
        ("IW9", layout.iw9), ("IW11", layout.iw11),
        ("IW12", layout.iw12), ("IW14", layout.iw14),
        ("IW15", layout.iw15), ("IW16", layout.iw16),
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

    # O6 door tip
    o6 = [o for o in outer_openings if o.name == "O6"][0]
    wall_mid = (pts["F9"][1] + pts["W9"][1]) / 2
    gap = (O6_WIDTH - O6_DOOR_WIDTH) / 2
    door_e = o6.poly[1][0] - gap
    tips.append(("O6_door_tip", (door_e, wall_mid - O6_DOOR_WIDTH)))

    # RO1 door tip
    ro1 = [r for r in rough_openings if r.name == "RO1"][0]
    ro1_mid = (ro1.bbox.s + ro1.bbox.n) / 2
    ro1_gap = (ro1.bbox.e - ro1.bbox.w - RO1_DOOR_WIDTH) / 2
    tips.append(("RO1_door_tip", (ro1.bbox.e - ro1_gap, ro1_mid - RO1_DOOR_WIDTH)))

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

    # RO3 door tip
    ro3 = [r for r in rough_openings if r.name == "RO3"][0]
    ro3_mid = (ro3.bbox.w + ro3.bbox.e) / 2
    ro3_gap = (ro3.bbox.n - ro3.bbox.s - RO3_DOOR_WIDTH) / 2
    tips.append(("RO3_door_tip", (ro3_mid - RO3_DOOR_WIDTH, ro3.bbox.n - ro3_gap)))

    # RO4 door tip
    ro4 = [r for r in rough_openings if r.name == "RO4"][0]
    ro4_mid = (ro4.bbox.w + ro4.bbox.e) / 2
    ro4_gap = (ro4.bbox.n - ro4.bbox.s - RO4_DOOR_WIDTH) / 2
    tips.append(("RO4_door_tip", (ro4_mid - RO4_DOOR_WIDTH, ro4.bbox.n - ro4_gap)))

    # RO5 door tip
    ro5 = [r for r in rough_openings if r.name == "RO5"][0]
    ro5_mid = (ro5.bbox.s + ro5.bbox.n) / 2
    ro5_gap = (ro5.bbox.e - ro5.bbox.w - RO5_DOOR_WIDTH) / 2
    tips.append(("RO5_door_tip", (ro5.bbox.e - ro5_gap, ro5_mid + RO5_DOOR_WIDTH)))

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
    """
    F8 = pts["F8"]
    C8 = pts["C8"]
    R_a8 = C8[0] - F8[0]
    result = []

    for prefix, inset, R_turn in [
        ("w_f8f9", WALL_OUTER, F8F9_INNER_TURN_R),
        ("g_f8f9", SHELL_THICKNESS + AIR_GAP, OPENING_INSIDE_RADIUS),
    ]:
        poly = f8f9_corner_polyline(pts, inset, R_turn)
        d = R_a8 + inset - R_turn
        arc_cx = F8[0] - inset + R_turn
        arc_cy = F8[1] - d
        result.append((f"{prefix}_start", poly[0]))
        result.append((f"{prefix}_arc_entry", poly[1]))
        result.append((f"{prefix}_arc_center", (arc_cx, arc_cy)))
        result.append((f"{prefix}_arc_exit", poly[-2]))
        result.append((f"{prefix}_end", poly[-1]))

    return result


# Expected distance-squared values: (name, d2_F1, d2_F6, d2_F12, d2_F15)
EXPECTED = [
    # F-series outline points
    ("F1", 0.000000000000, 648.677110458661, 1587.565495574659, 1251.174349116076),
    ("F2", 1.235511769340, 611.923611111111, 1597.307931538801, 1300.000000000000),
    ("F3", 167.258170314851, 165.847222222222, 1193.025913162186, 1397.673611111111),
    ("F4", 191.852404716788, 143.447592751261, 1170.809854023777, 1413.448717626247),
    ("F5", 540.809269886441, 9.686412271629, 1008.031197327137, 1635.330877054027),
    ("F6", 648.677110458661, 0.000000000000, 869.030869921322, 1567.923611111111),
    ("F7", 704.330886311816, 27.562500000000, 587.593521140158, 1253.361111111111),
    ("F8", 634.006349945043, 62.951388888889, 477.578506839646, 1031.694444444445),
    ("F9", 629.816502109322, 66.312500000000, 470.553624865800, 1016.694444444445),
    ("F10", 1176.080187907325, 531.423611111111, 42.998506164660, 493.444444444445),
    ("F11", 1264.265226111227, 575.200509747991, 30.334250997684, 506.454597300219),
    ("F12", 1587.565495574659, 869.030869921322, 0.000000000000, 441.369514523197),
    ("F13", 1384.631969215298, 1229.939467698260, 140.635738692275, 86.788992995208),
    ("F14", 1375.308225283584, 1253.034722222222, 156.457277587939, 75.111111111111),
    ("F15", 1251.174349116076, 1567.923611111111, 441.369514523197, 0.000000000000),
    ("F16", 1157.970574145304, 1589.839740474726, 528.662529511229, 6.117189499405),
    ("F17", 884.899519389345, 1466.579161769001, 651.613619484118, 52.536554860110),
    ("F18", 708.836266563858, 1351.117180668459, 721.664011621445, 107.147309451050),
    ("F19", 636.748102911967, 1286.506944444444, 739.309862891971, 133.611111111111),
    ("F20", 535.815723707719, 1190.883422017860, 766.856578431436, 179.118169391529),
    ("F11a", 1473.847237361418, 685.116495743439, 13.578290920455, 547.320195418060),
    ("F11b", 1532.547416769082, 738.465955940979, 8.070637064439, 535.502988948933),
    # IW wall corners
    ("IW1_SW", 154.954316204608, 229.173611111111, 826.760584640841, 892.722222222222),
    ("IW1_SE", 1310.277351734232, 1233.784722222222, 173.562983243486, 61.805555555556),
    ("IW1_NE", 1321.188921513127, 1219.368055555556, 160.875738804913, 69.888888888889),
    ("IW1_NW", 165.865885983502, 214.756944444444, 814.073340202269, 900.805555555556),
    ("IW2_SW", 165.865885983502, 214.756944444444, 814.073340202269, 900.805555555556),
    ("IW2_SE", 172.541245588565, 218.756944444444, 788.644783175491, 872.222222222222),
    ("IW2_NE", 656.153629618721, 18.506944444444, 635.089183334032, 1279.472222222223),
    ("IW2_NW", 649.478270013659, 14.506944444444, 660.517740360810, 1308.055555555556),
    ("IW3_SW", 56.626534767540, 673.414631034161, 1168.935129387146, 778.215248855451),
    ("IW3_SE", 61.734626759281, 678.614749297631, 1154.508883050812, 760.169530880165),
    ("IW3_NE", 115.067960092614, 392.364532295004, 847.845388184693, 724.108071284885),
    ("IW3_NW", 109.959868100874, 387.164414031534, 862.271634521027, 742.153789260171),
    ("IW4_SW", 658.972222222222, 1272.284722222222, 697.964979947096, 116.805555555555),
    ("IW4_SE", 676.144684181153, 1287.673611111111, 693.734830818134, 110.472222222222),
    ("IW4_NE", 785.793340802040, 755.673611111111, 205.546860596518, 148.472222222222),
    ("IW4_NW", 768.620878843110, 740.284722222222, 209.777009725481, 154.805555555556),
    ("IW5_SW", 707.074771949109, 819.142843680960, 291.211352112973, 123.297572800617),
    ("IW5_SE", 1259.201217950312, 1322.013888888889, 252.280327655635, 26.284722222222),
    ("IW5_NE", 1263.219502839759, 1313.368055555556, 244.499205436349, 28.888888888889),
    ("IW5_NW", 712.504197709980, 811.759520129881, 283.057819841120, 125.354064064353),
    ("IW6_SW", 364.296357839032, 43.523501719796, 1015.249196364029, 1467.836174220158),
    ("IW6_SE", 405.264351927112, 53.125000000000, 679.823914369314, 1095.423611111111),
    ("IW6_NE", 408.451002445817, 52.090277777778, 679.077429185108, 1098.138888888889),
    ("IW6_NW", 367.493524670252, 42.449752113698, 1013.918710517674, 1469.909029552332),
    ("IW7_SW", 110.290182314836, 404.565932034024, 861.067451816888, 723.800033153538),
    ("IW7_SE", 155.684205586224, 450.650152343379, 759.953937627716, 595.540481672229),
    ("IW7_NE", 160.461983364002, 438.448752604359, 746.731873995521, 595.848519803576),
    ("IW7_NW", 115.067960092614, 392.364532295004, 847.845388184693, 724.108071284885),
    ("IW8_SW", 149.943223416408, 178.736111111111, 1164.370829711247, 1337.118055555556),
    ("IW8_SE", 191.222898282218, 185.236111111111, 788.299588363139, 920.034722222222),
    ("IW8_NE", 203.717801394446, 172.402777777778, 777.195677257900, 929.701388888889),
    ("IW8_NW", 162.438126528636, 165.902777777778, 1153.266918606007, 1346.784722222222),
    ("IW9_SW", 107.128650030668, 724.698969606986, 1053.395368861640, 631.909979398856),
    ("IW9_SE", 114.125630911298, 731.787976759344, 1040.858011414195, 615.753150312460),
    ("IW9_NE", 167.458964244631, 445.537759756717, 734.194516548076, 579.691690717180),
    ("IW9_NW", 160.461983364002, 438.448752604359, 746.731873995521, 595.848519803577),
    ("IW11_SW", 514.016942613957, 1134.992234244248, 737.513143306170, 182.107303202177),
    ("IW11_SE", 529.236145716808, 1150.303463618829, 733.198008080947, 174.172696338002),
    ("IW11_NE", 647.965312383475, 746.922921643957, 298.432051390955, 155.457368876926),
    ("IW11_NW", 632.746109280623, 731.611692269376, 302.747186616178, 163.391975741101),
    ("IW12_SW", 573.236145716808, 888.678268316465, 453.200862701440, 137.717382702251),
    ("IW12_SE", 669.883341965209, 985.881295448437, 431.862310653484, 94.517528644112),
    ("IW12_NE", 672.587767306587, 971.596318131714, 418.737376621343, 95.324859890951),
    ("IW12_NW", 577.791701272364, 876.254646355222, 439.756576847022, 137.803198611375),
    ("IW14_SW", 642.569479050142, 754.261471448222, 306.536099115101, 153.413840278416),
    ("IW14_SE", 707.074771949109, 819.142843680960, 291.211352112973, 123.297572800617),
    ("IW14_NE", 712.470605282442, 811.804293876695, 283.107304388827, 125.341101399127),
    ("IW14_NW", 647.965312383475, 746.922921643957, 298.432051390955, 155.457368876926),
    ("IW15_SW", 632.746109280623, 731.611692269376, 302.747186616178, 163.391975741101),
    ("IW15_SE", 648.676603390712, 745.758613309423, 297.275069638373, 155.816674558925),
    ("IW15_NE", 692.668540814074, 672.516942963880, 232.038157444507, 186.717911299305),
    ("IW15_NW", 676.738046703986, 658.370021923833, 237.510274422312, 194.293212481480),
    ("IW16_SW", 109.959868100874, 387.164414031534, 862.271634521027, 742.153789260171),
    ("IW16_SE", 115.577284692528, 390.998257553147, 846.486440024787, 724.265410559561),
    ("IW16_NE", 187.505136843275, 250.127111983168, 720.145296082256, 772.431407745544),
    ("IW16_NW", 181.887720251620, 246.293268461556, 735.930490578496, 790.319786446154),
    ("dryer_SW", 1.172792399273, 597.124002442594, 1514.421689752237, 1216.718485543654),
    ("dryer_SE", 12.161001206581, 592.505946887039, 1338.137051540480, 1022.030985543653),
    ("dryer_NE", 23.390613655494, 477.094377108144, 1231.372592902056, 1019.119415764759),
    ("dryer_NW", 12.402404848185, 481.712432663699, 1407.657231113813, 1213.806915764759),
    ("washer_SW", 12.992003040927, 478.080658115514, 1404.313693603644, 1213.925141216574),
    ("washer_SE", 23.980211848235, 473.462602559958, 1228.029055391886, 1019.237641216573),
    ("washer_NE", 48.126490963814, 370.967699447731, 1134.181263420129, 1029.242738104345),
    ("washer_NW", 37.138282156506, 375.585755003286, 1310.465901631886, 1223.930238104346),
    ("ctr_poly_0", 84.611391498934, 361.850833165988, 950.510714243831, 850.875871822603),
    ("ctr_poly_1", 112.608744369013, 379.588546071989, 855.393916232578, 743.660936046535),
    ("ctr_poly_2", 109.959868100874, 387.164414031534, 862.271634521027, 742.153789260171),
    ("ctr_poly_3", 56.626534767540, 673.414631034161, 1168.935129387146, 778.215248855451),
    ("ctr_poly_4", 40.222446462744, 656.682425569572, 1222.180191158826, 844.365437954458),
    ("bed_SW", 185.078380283472, 795.056367250748, 931.885379710730, 481.875192649939),
    ("bed_SE", 397.187683682098, 1008.914169812226, 772.842254875939, 254.062106675068),
    ("bed_NE", 471.604350348764, 684.320164834139, 424.262648408249, 223.439891650614),
    ("bed_NW", 259.495046950139, 470.462362272662, 583.305773243040, 451.252977625486),
    ("dresser_SW", 510.602811561319, 592.714399230079, 342.661253128909, 249.137589787726),
    ("dresser_SE", 636.095344830404, 703.046561403810, 286.231592150901, 174.830863072568),
    ("dresser_NE", 667.086149130237, 653.831283626032, 242.492818095421, 196.865585294790),
    ("dresser_NW", 541.593615861152, 543.499121452301, 298.922479073429, 271.172312009948),
    ("RO1_SW", 388.964405951987, 408.784722222222, 407.991036207476, 409.805555555556),
    ("RO1_SE", 504.075016784049, 506.951388888889, 319.776841704552, 301.611111111111),
    ("RO1_NE", 514.986586562943, 492.534722222222, 307.089597265979, 309.694444444444),
    ("RO1_NW", 399.875975730882, 394.368055555556, 395.303791768904, 417.888888888889),
    ("RO2_SW", 581.354201272364, 867.082763217624, 429.819195789543, 138.013393876552),
    ("RO2_SE", 566.134998169512, 851.771533843043, 434.134331014766, 145.948000740727),
    ("RO2_NE", 622.079442613957, 746.413791877907, 319.080282064470, 159.429918544080),
    ("RO2_NW", 637.298645716808, 761.725021252487, 314.765146839248, 151.495311679906),
    ("RO3_SW", 116.908674363341, 367.919803511250, 844.815521092168, 746.177750117867),
    ("RO3_SE", 122.526090954995, 371.753647032863, 829.030326595929, 728.289371417258),
    ("RO3_NE", 176.746752471835, 265.562144394480, 733.791831402143, 764.597868778875),
    ("RO3_NW", 171.129335880181, 261.728300872868, 749.577025898382, 782.486247479485),
    ("RO4_SW", 271.158444214660, 117.423611111111, 730.575384693688, 983.472222222222),
    ("RO4_SE", 277.833803819722, 121.423611111111, 705.146827666911, 954.888888888889),
    ("RO4_NE", 383.884856863833, 67.062500000000, 661.738724000396, 1043.027777777778),
    ("RO4_NW", 377.209497258770, 63.062500000000, 687.167281027173, 1071.611111111111),
    ("RO5_SW", 373.031339070297, 39.173611111111, 866.967665163405, 1304.125000000000),
    ("RO5_SE", 402.114172124581, 51.312500000000, 692.725692882703, 1109.902777777778),
    ("RO5_NE", 405.300822643286, 50.277777777778, 691.979207698497, 1112.618055555556),
    ("RO5_NW", 376.217989589002, 38.138888888889, 866.221179979198, 1306.840277777778),
    ("RO6_SW", 529.965312383475, 1129.808658389498, 711.427372985148, 169.314688446631),
    ("RO6_SE", 514.746109280623, 1114.497429014918, 715.742508210371, 177.249295310805),
    ("RO6_NE", 552.634998169512, 889.209066393437, 474.633855244684, 145.857219680019),
    ("RO6_NW", 567.854201272364, 904.520295768018, 470.318720019461, 137.922612815844),
    ("RO7_SW", 115.243686466853, 703.192471660503, 1010.476344502298, 609.049161486761),
    ("RO7_SE", 108.246705586224, 696.103464508144, 1023.013701949744, 625.205990573158),
    ("RO7_NE", 147.857816697335, 472.537324108887, 783.627271206279, 595.536137164594),
    ("RO7_NW", 154.854797577964, 479.626331261245, 771.089913758834, 579.379308078197),
    ("O1_outer_start", 60.914517467422, 319.983173278988, 1329.376211517656, 1320.415990789896),
    ("O1_outer_end", 88.024848225842, 266.887421959796, 1281.756963920763, 1338.570239470705),
    ("O1_inner_end", 87.480883254814, 262.776310848685, 1238.407776773948, 1291.014683915149),
    ("O1_inner_start", 60.370552496394, 315.872062167876, 1286.027024370841, 1272.860435234340),
    ("O2_outer_start", 254.866568734819, 97.594017478682, 1118.411946994519, 1447.803634160523),
    ("O2_outer_end", 307.439515113277, 69.515870652350, 1085.485646675607, 1479.145317218436),
    ("O2_inner_end", 304.932360650941, 67.192768733472, 1043.909815557122, 1430.330542822047),
    ("O2_inner_start", 252.359414272482, 95.270915559804, 1076.836115876033, 1398.988859764135),
    ("O3_outer_start", 410.911044520690, 32.600784557060, 1040.131544651140, 1545.660308506152),
    ("O3_outer_end", 525.487244845802, 11.343564747788, 1010.709013696471, 1624.478591659819),
    ("O3_inner_end", 522.980090383466, 9.020462828910, 969.133182577985, 1575.663817263430),
    ("O3_inner_start", 408.403890058354, 30.277682638182, 998.555713532654, 1496.845534109763),
    ("O4_outer_start", 632.452191198471, 5.506944444444, 739.803411441142, 1396.805555555556),
    ("O4_outer_end", 640.402730606065, 9.444444444444, 699.598075900976, 1351.868055555556),
    ("O4_inner_end", 673.729268089036, 9.000000000000, 701.459527760657, 1381.423611111111),
    ("O4_inner_start", 665.778728681442, 5.062500000000, 741.664863300823, 1426.361111111111),
    ("O5_outer_start", 617.893341745899, 83.701388888889, 436.589022405934, 949.777777777778),
    ("O5_outer_end", 777.602972825495, 213.090277777778, 232.454264991345, 709.888888888889),
    ("O5_inner_end", 807.596176975132, 209.312500000000, 230.982383517692, 736.111111111111),
    ("O5_inner_start", 647.886545895536, 79.923611111111, 435.117140932282, 976.000000000000),
    ("O6_outer_start", 950.181209271056, 361.590277777778, 116.097251750347, 565.138888888889),
    ("O6_outer_end", 1120.744957485958, 512.534722222222, 51.232277998423, 477.138888888889),
    ("O6_inner_end", 1150.738161635596, 508.756944444445, 49.760396524771, 503.361111111111),
    ("O6_inner_start", 980.174413420694, 357.812500000000, 114.625370276695, 591.361111111111),
    ("O7_outer_start", 1533.623099271676, 910.179528413682, 4.000000000000, 361.852099771179),
    ("O7_outer_end", 1419.795910362730, 1081.625503890764, 64.000000000000, 171.299855515127),
    ("O7_inner_end", 1370.749921792452, 1044.765687197968, 64.444444444444, 168.639578366924),
    ("O7_inner_start", 1484.577110701398, 873.319711720886, 4.444444444444, 359.191822622977),
    ("O8_outer_start", 1262.637375638132, 1499.180555555555, 378.102962689939, 2.506944444444),
    ("O8_outer_end", 1283.222004235006, 1422.673611111111, 308.225470119856, 12.250000000000),
    ("O8_inner_end", 1236.654858094923, 1379.673611111111, 304.463546155559, 12.694444444444),
    ("O8_inner_start", 1216.070229498049, 1456.180555555556, 374.341038725643, 2.951388888889),
    ("O9_outer_start", 491.160360181902, 1145.955967438196, 779.957751186172, 202.920915013522),
    ("O9_outer_end", 409.879109007171, 1064.145565201023, 810.998945397870, 254.774071149192),
    ("O9_inner_end", 410.323553451616, 1030.631654611872, 775.443707022369, 246.279036300775),
    ("O9_inner_start", 491.604804626346, 1112.442056849044, 744.402512810671, 194.425880165105),
    ("O10_outer_start", 171.055131184750, 823.343022536428, 992.163241629052, 513.757007617833),
    ("O10_outer_end", 124.593324454464, 776.352064743700, 1058.023880285195, 600.429608197948),
    ("O10_inner_end", 125.037768898908, 742.838154154548, 1022.468641909694, 591.934573349531),
    ("O10_inner_start", 171.499575629194, 789.829111947276, 956.608003253551, 505.261972769416),
    ("O11_outer_start", 26.647668603134, 676.749936001248, 1311.696314339207, 919.252531753089),
    ("O11_outer_end", 12.807842753478, 662.472985360878, 1390.644595547905, 1015.393303246807),
    ("O11_inner_end", 13.252287197922, 628.959074771727, 1355.089357172405, 1006.898268398391),
    ("O11_inner_start", 27.092113047579, 643.236025412096, 1276.141075963706, 910.757496904673),
    ("O3_door_tip", 517.174943231824, 8.058280822707, 840.854796646352, 1420.476340237726),
    ("O6_door_tip", 986.241066729893, 538.888888888889, 71.538196724957, 366.423611111111),
    ("RO1_door_tip", 449.712711954840, 592.340277777778, 400.558667843386, 268.680555555556),
    ("RO2_door_tip", 773.904927646603, 903.358320871178, 289.829513332901, 92.411643801585),
    ("RO3_door_tip", 130.669713115298, 240.650421944265, 894.972733083503, 942.297257439771),
    ("RO4_door_tip", 346.384757801112, 51.173611111111, 836.847218747545, 1235.180555555556),
    ("RO5_door_tip", 526.402897789791, 21.967013888889, 678.805618163627, 1222.862847222222),
    ("RO6_door_tip_S", 414.621409782884, 1009.483751235875, 747.858216656660, 238.083830670642),
    ("RO6_door_tip_N", 451.288076449551, 791.462755150572, 514.527262173737, 207.704402640848),
    ("RO7_door_tip_S", 170.660052631259, 755.008227347813, 924.652382017689, 501.271422970377),
    ("RO7_door_tip_N", 208.993385964592, 538.653897929176, 692.988094201433, 472.558661607250),
    # W8-W9 shell polyline special points (W-series: inner wall face)
    ("w_f8f9_start", 620.772537138294, 53.284722222222, 507.149916208683, 1065.472222222222),
    ("w_f8f9_arc_entry", 592.023294535794, 56.644503775119, 508.317280606279, 1040.308381727875),
    ("w_f8f9_arc_center", 595.889824977499, 59.443699151287, 499.375949294599, 1030.108451994682),
    ("w_f8f9_arc_exit", 587.008673504063, 60.667250931724, 499.909547850179, 1022.355625822363),
    ("w_f8f9_end", 599.823297959684, 70.090277777778, 472.025506339452, 990.472222222222),
    # W8-W9 shell polyline special points (G-series: outer face of inner shell)
    ("g_f8f9_start", 623.997657006648, 55.618055555556, 499.673730533090, 1056.944444444445),
    ("g_f8f9_arc_entry", 595.248414404148, 58.977837108452, 500.841094930686, 1031.780603950097),
    ("g_f8f9_arc_center", 595.889824977499, 59.443699151287, 499.375949294599, 1030.108451994682),
    ("g_f8f9_arc_exit", 594.423641208139, 59.639473153946, 499.458244148433, 1028.827848044585),
    ("g_f8f9_end", 607.238265663760, 69.062500000000, 471.574202637706, 996.944444444445),
]

REF_NAMES = ["F1", "F6", "F12", "F15"]


class TestDistanceSquaredRegression:
    """Verify that every interior point stays at the expected distance
    from F1, F6, F12, F15 (within +/- 0.000001 ft^2)."""

    @pytest.fixture(scope="class")
    def all_pts(self, pts_with_outline, layout):
        pts, _ = pts_with_outline
        return _collect_all_points(pts, layout)

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

    _F_NAMES = [f"F{i}" for i in range(1, 21)] + ["F11a", "F11b"]

    @pytest.mark.parametrize("fname", _F_NAMES)
    def test_chain_matches_geometry(self, pts_with_outline, fname):
        """Chain-reconstructed F points must match geometry-computed ones."""
        pts, _ = pts_with_outline
        chain_pts = _reconstruct_f_points()
        geo_pt = pts[fname]
        chain_pt = chain_pts[fname]
        # 1e-9 ft positional tolerance (~0.000012 inches)
        assert abs(chain_pt[0] - geo_pt[0]) < 1e-9, (
            f"{fname} E: chain={chain_pt[0]:.12f}, geo={geo_pt[0]:.12f}")
        assert abs(chain_pt[1] - geo_pt[1]) < 1e-9, (
            f"{fname} N: chain={chain_pt[1]:.12f}, geo={geo_pt[1]:.12f}")

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
