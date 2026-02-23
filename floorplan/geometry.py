"""Compute F-series outline geometry from chain walk and design constants."""
import math
from typing import NamedTuple


from shared.types import Point, LineSeg, ArcSeg, Segment
from shared.geometry import left_norm, off_pt, poly_area
from floorplan.constants import (
    CORNER_NE_R, CORNER_NW_R, UPPER_E_R, SMALL_ARC_R, ARC_180_R,
    F6_F7_LENGTH, F3_OFFSET_FROM_IW8, F6_OFFSET_ADJ,
    F6_HEIGHT, NW_SHIFT,
    F14_F15_SEG, F14_F15_DIST, ARC_F13_R_BASELINE, F13_EXIT_BRG,
    SOUTH_WALL_FACE, PIX_PI5_TARGET_BRG, F15_OFFSET_FROM_IW8, ARC_F17_SWEEP, F16_F17_MIN,
    F18_OFFSET_FROM_IW4, F19_OFFSET_FROM_IW4, ARC_F19_R,
    WALL_OUTER, WALL_EXTRA, WALL_6IN, WALL_3IN, WALL_4IN,
    APPLIANCE_WIDTH, COUNTER_GAP, COUNTER_DEPTH,
    CLOSET_WIDTH, CLOSET2_WIDTH, BEDROOM_WIDTH, APPLIANCE_OFFSET_FROM_W2,
    IW1_OFFSET_FROM_W1, IW1_OFFSET_FROM_W9, IW8_OFFSET_FROM_IW1, F14_OFFSET_FROM_IW1,
    F10_OFFSET_FROM_F9, FLAT_SEG_11, F11AB_TARGET,
)


class OutlineAnchors(NamedTuple):
    """Anchor points from inset path computation."""
    Pi2: Point       # NW anchor (Po2 alias)
    Pi3: Point       # NE anchor (sets F2 easting)
    Ti3: Point       # F21 northing anchor
    PiX: Point       # south wall line start
    Pi5: Point       # south wall line end
    TC1: Point       # Arc 1 center (for west_E)
    R1i: float       # Arc 1 inset radius


class OutlineGeometry(NamedTuple):
    """Complete outline geometry result."""
    fp_pts: dict[str, Point]     # F1-F20, F11a, F11b + C1-C19, C11a + FC
    outline_segs: list[Segment]  # 22 segments with F-series names
    radii: dict[str, float]      # R_a1 through R_a19


# ============================================================
# F-series outline chain: single source of truth
# ============================================================

# FC (building center) = origin, by definition.
# F2 position and initial bearing define the chain starting point.
F2_E = -18.0           # F2 easting relative to FC
F2_N = -10.5           # F2 northing relative to FC
F2_BRG = 0.0           # initial bearing: due north (radians)

# Sweep angle constants (radians)
_A19 = math.atan(1.0 / 9.0)   # arctan(1/9) for F3-F4, F19-F20
_A9 = math.atan(9.0)           # arctan(9) for F5-F6, F1-F2
_PI_2 = math.pi / 2            # 90 deg
_5PI_12 = 5 * math.pi / 12    # 75 deg
_PI_12 = math.pi / 12          # 15 deg
_PI_3 = math.pi / 3            # 60 deg
_PI_6 = math.pi / 6            # 30 deg

# Chain: ("L", distance) for lines
#        ("CW"/"CCW", radius, sweep, center_name, n_pts) for arcs
# Starting at F2, bearing north (0 rad), CW traversal.
OUTLINE_CHAIN = [
    ("L",   12.083333333333),                                  # F2->F3
    ("CW",   8.351795046046, _A19, "C3", 20),                 # F3->F4
    ("L",    9.476667232982),                                  # F4->F5
    ("CW",   2.333333333333, _A9, "C5", 20),                  # F5->F6
    ("L",    5.250000000000),                                  # F6->F7
    ("CW",   2.333333333333, _PI_2, "C7", 20),                # F7->F8
    ("CCW",  0.166666666667, _PI_2, "C8", 20),                # F8->F9
    ("L",   15.166666666667),                                  # F9->F10
    ("CCW",  1.039662132188, _5PI_12, "C10", 20),             # F10->F11
    ("CW",   2.333333333333, _5PI_12, "C11a", 30),            # F11->F11a
    ("L",    1.000000000000),                                  # F11a->F11b
    ("CW",   2.333333333333, _5PI_12, "C11", 30),             # F11b->F12
    ("L",   11.858994000010),                                  # F12->F13
    ("CW",   2.507553207938, _PI_12, "C13", 60),              # F13->F14
    ("L",    8.666666666667),                                  # F14->F15
    ("CW",   2.473295271375, _PI_3, "C15", 20),               # F15->F16
    ("L",    5.000000000000),                                  # F16->F17
    ("CW",   6.404672887007, _PI_6, "C17", 20),               # F17->F18
    ("L",    1.397555568554),                                  # F18->F19
    ("CW",  18.888718471469, _A19, "C19", 60),                # F19->F20
    ("L",   23.147693701700),                                  # F20->F1
    ("CW",   0.833333333333, _A9, "C1", 20),                  # F1->F2
]

# Point names produced by each chain segment (one per segment, in order)
CHAIN_POINT_NAMES = [
    "F3", "F4", "F5", "F6", "F7", "F8", "F9", "F10", "F11",
    "F11a", "F11b", "F12", "F13", "F14", "F15", "F16", "F17",
    "F18", "F19", "F20", "F1", "F2",
]


def walk_outline_chain() -> dict[str, Point]:
    """Walk F-series chain from F2 bearing north. Returns all F/C/FC points.

    FC = (0, 0) by definition. All F-series and arc center points are
    computed from the chain walk starting at (F2_E, F2_N) with bearing F2_BRG.
    """
    E, N = F2_E, F2_N
    brg = F2_BRG
    fp_pts: dict[str, Point] = {"FC": (0.0, 0.0)}

    for seg, name in zip(OUTLINE_CHAIN, CHAIN_POINT_NAMES):
        if seg[0] == "L":
            d = seg[1]
            E += d * math.sin(brg)
            N += d * math.cos(brg)
        else:
            direction, R, sweep, center_name = seg[0], seg[1], seg[2], seg[3]
            if direction == "CW":
                cx = E + R * math.cos(brg)
                cy = N - R * math.sin(brg)
                alpha = math.atan2(N - cy, E - cx) - sweep
                E, N = cx + R * math.cos(alpha), cy + R * math.sin(alpha)
                brg += sweep
            else:  # CCW
                cx = E - R * math.cos(brg)
                cy = N + R * math.sin(brg)
                alpha = math.atan2(N - cy, E - cx) + sweep
                E, N = cx + R * math.cos(alpha), cy + R * math.sin(alpha)
                brg -= sweep
            fp_pts[center_name] = (cx, cy)
        fp_pts[name] = (E, N)

    return fp_pts


def _build_outline_segs() -> list[Segment]:
    """Build outline segment list from chain definition.

    Outline convention starts at F1, so we rotate the chain (which starts
    at F2) so the F1->F2 arc comes first.
    """
    start_names = ["F2"] + CHAIN_POINT_NAMES[:-1]  # F2, F3, ..., F20, F1
    end_names = CHAIN_POINT_NAMES                    # F3, F4, ..., F1, F2

    segs: list[Segment] = []
    for entry, start, end in zip(OUTLINE_CHAIN, start_names, end_names):
        if entry[0] == "L":
            segs.append(LineSeg(start, end))
        else:
            segs.append(ArcSeg(start, end, entry[3], entry[1], entry[0], entry[4]))

    # Rotate so F1->F2 comes first (matches outline convention)
    return segs[-1:] + segs[:-1]


def _build_radii() -> dict[str, float]:
    """Extract radii dict from chain arc entries."""
    radii: dict[str, float] = {}
    for entry in OUTLINE_CHAIN:
        if entry[0] != "L":
            center_name = entry[3]
            ra_name = "R_a" + center_name[1:]  # "C3" -> "R_a3"
            if ra_name == "R_a11a":
                ra_name = "R_a11"  # C11a and C11 share the same radius
            radii[ra_name] = entry[1]
    return radii


# ============================================================
# Per-section helpers (mutate fp_pts, return computed radii)
# ============================================================

def _compute_ne_corner(fp_pts: dict[str, Point], anchors: OutlineAnchors) -> float:
    """NE corner arc: F1, F2, C1. Returns R_a1."""
    R_a1 = CORNER_NE_R
    fp_pts["C1"] = (anchors.Pi3[0] + R_a1, anchors.Pi3[1] + R_a1)
    fp_pts["F2"] = (anchors.Pi3[0], fp_pts["C1"][1])
    fp_pts["F1"] = (fp_pts["C1"][0], anchors.Pi3[1])
    return R_a1


def _compute_nw_corner(fp_pts: dict[str, Point], anchors: OutlineAnchors) -> float:
    """NW corner: F6, C5. Returns R_a5.

    F6 exits at bearing 90° (east). F5 computed later by tangent constraint.
    """
    R_a5 = CORNER_NW_R
    corner2_N = fp_pts["F1"][1] + F6_HEIGHT
    fp_pts["C5"] = (anchors.Pi2[0] + NW_SHIFT + R_a5 + F6_OFFSET_ADJ, corner2_N - R_a5)
    fp_pts["F6"] = (fp_pts["C5"][0], corner2_N)
    return R_a5


def _compute_west_wall(
    fp_pts: dict[str, Point], anchors: OutlineAnchors,
) -> float:
    """West wall: F3, F4, F5, C3; update C5, F6. Returns R_a3.

    F3 northing = 2" north of IW8 north face.
    F3-F4 arc sweep = arctan(1/9).  R_a3 solved from F6-F7 target length,
    then C5 easting (and F6) solved from tangency constraint.
    Depends on F2, F6, C5 already in fp_pts.
    """
    R_a5 = CORNER_NW_R

    F3_E = fp_pts["F2"][0]
    _, C5_N = fp_pts["C5"]

    # F3_N: 2" north of IW8 north face
    # F9_N = F6_N - UPPER_E_R - SMALL_ARC_R + WALL_EXTRA (arcs cancel WALL_EXTRA)
    _F9_N = fp_pts["F6"][1] - (UPPER_E_R + SMALL_ARC_R) + WALL_EXTRA
    _W9_N = _F9_N - WALL_OUTER
    _IW1_n = _W9_N - IW1_OFFSET_FROM_W9
    _IW8_n = _IW1_n + IW8_OFFSET_FROM_IW1
    F3_N = _IW8_n + F3_OFFSET_FROM_IW8

    fp_pts["F3"] = (F3_E, F3_N)

    # F3-F4 arc sweep: arctan(1/9)
    _gamma = math.atan(1.0 / 9.0)
    _cos_g = math.cos(_gamma)
    _tan_g = math.tan(_gamma)

    # Solve R_a3 from F6-F7 target length:
    # F7_E (independent of F6 easting) - C5_E = F6_F7_LENGTH
    # C5_E = F3_E + R_a3 + ((C5_N - F3_N)*sin(γ) - (R_a3 - R_a5)) / cos(γ)
    # Combining: R_a3 = _K / (1 - 1/cos(γ))
    F7_E = anchors.Pi2[0] + CORNER_NW_R + 5.5 + 10.0 / 12.0
    _K = (F7_E - F6_F7_LENGTH - F3_E
          - (C5_N - F3_N) * _tan_g - R_a5 / _cos_g)
    R_a3 = _K / (1.0 - 1.0 / _cos_g)

    fp_pts["C3"] = (F3_E + R_a3, F3_N)

    # Solve C5_E from tangency with F3 fixed:
    # (C5_N - F3_N)·sin(γ) - (C5_E - F3_E - R_a3)·cos(γ) = R_a3 - R_a5
    C5_E = F3_E + R_a3 + ((C5_N - F3_N) * math.sin(_gamma) - (R_a3 - R_a5)) / _cos_g
    fp_pts["C5"] = (C5_E, C5_N)
    fp_pts["F6"] = (C5_E, C5_N + R_a5)

    # F4 and F5 at angle (π - γ) from their respective centers
    _angle = math.pi - _gamma
    C3_E, C3_N = fp_pts["C3"]
    fp_pts["F4"] = (C3_E + R_a3 * math.cos(_angle),
                    C3_N + R_a3 * math.sin(_angle))
    fp_pts["F5"] = (C5_E + R_a5 * math.cos(_angle),
                    C5_N + R_a5 * math.sin(_angle))

    return R_a3


def _compute_upper_east_arcs(
    fp_pts: dict[str, Point], anchors: OutlineAnchors,
) -> tuple[float, float]:
    """Upper east arcs: F7, F8, C7, F9, C8. Returns (R_a7, R_a8).

    Depends on F6 already in fp_pts.
    C8 placed for tangency with C7: |C7-C8| = R_a7 + R_a8.
    """
    # F7: computed from anchors (independent of F6 easting shift).
    # Original: F6_E + 5.5 + 6/12 - NW_SHIFT + 4/12 - F6_OFFSET_ADJ
    # where F6_E = Pi2[0] + NW_SHIFT + R_a5 + F6_OFFSET_ADJ;
    # NW_SHIFT and F6_OFFSET_ADJ cancel, leaving:
    F7_E = anchors.Pi2[0] + CORNER_NW_R + 5.5 + 10.0/12.0
    fp_pts["F7"] = (F7_E, fp_pts["F6"][1])
    R_a7 = UPPER_E_R
    fp_pts["C7"] = (fp_pts["F7"][0], fp_pts["F7"][1] - R_a7)

    # F8-F9: small arc (R=2"), C8 via tangency with C7
    R_a8 = SMALL_ARC_R
    _sum_R = R_a7 + R_a8
    _dE = math.sqrt(_sum_R**2 - WALL_EXTRA**2)
    fp_pts["C8"] = (fp_pts["C7"][0] + _dE, fp_pts["C7"][1] + WALL_EXTRA)
    # F8: tangent point on C7→C8 line
    fp_pts["F8"] = (fp_pts["C7"][0] + R_a7 * _dE / _sum_R,
                    fp_pts["C7"][1] + R_a7 * WALL_EXTRA / _sum_R)
    # F9: bottom of C8 circle
    fp_pts["F9"] = (fp_pts["C8"][0], fp_pts["C8"][1] - R_a8)

    return R_a7, R_a8


def _compute_central_region(
    fp_pts: dict[str, Point], anchors: OutlineAnchors,
) -> tuple[float, float, float, float]:
    """Central region: F10-F16, C10, C11, C13, C15. Returns (R_a10, R_a11, R_a13, R_a15).

    Depends on F0, F2, F9, F16 already in fp_pts.
    """
    R_a11 = ARC_180_R

    # R_a15 baseline: use existing F14_N formula to keep F15/F16 fixed
    _F14_N_R15 = fp_pts["F9"][1] - WALL_OUTER - IW1_OFFSET_FROM_W9 + F14_OFFSET_FROM_IW1

    # Arc at Po5 corner (exits North)
    d_in_po5 = (anchors.Pi5[0] - anchors.PiX[0], anchors.Pi5[1] - anchors.PiX[1])
    L_in = math.sqrt(d_in_po5[0]**2 + d_in_po5[1]**2)
    d_in_u = (d_in_po5[0]/L_in, d_in_po5[1]/L_in)
    # F15 E-coordinate
    _iw8_e = (fp_pts["F2"][0] + WALL_OUTER + APPLIANCE_OFFSET_FROM_W2 + APPLIANCE_WIDTH + COUNTER_GAP
              + COUNTER_DEPTH + (WALL_3IN + CLOSET_WIDTH + WALL_4IN
              + BEDROOM_WIDTH + WALL_4IN + CLOSET_WIDTH + WALL_3IN))
    F15_E = _iw8_e + F15_OFFSET_FROM_IW8
    ln_in_po5 = left_norm(anchors.PiX, anchors.Pi5)

    # R_a15 from baseline F14_N (keeps F15/F16 fixed)
    _A_15 = anchors.PiX[1] + (F15_E - anchors.PiX[0]) * d_in_u[1] / d_in_u[0]
    _B_15 = ln_in_po5[1] - (1.0 + ln_in_po5[0]) * d_in_u[1] / d_in_u[0]
    R_a15 = (_F14_N_R15 - F14_F15_SEG - _A_15) / _B_15
    o_in_po5 = off_pt(anchors.PiX, ln_in_po5, R_a15)
    t_cf4 = (F15_E - R_a15 - o_in_po5[0]) / d_in_u[0]
    fp_pts["C15"] = (F15_E - R_a15, o_in_po5[1] + t_cf4 * d_in_u[1])
    fp_pts["F15"] = (F15_E, fp_pts["C15"][1])
    # F16: tangent point for 60° incoming bearing
    _brg_f4 = math.radians(PIX_PI5_TARGET_BRG)
    fp_pts["F16"] = (fp_pts["C15"][0] + R_a15 * math.cos(_brg_f4),
                     fp_pts["C15"][1] - R_a15 * math.sin(_brg_f4))

    # F13-F14 arc direction: bearing F13->F12 = 345°
    _brg_off = math.radians(360.0 - F13_EXIT_BRG)
    _nx_t = math.cos(_brg_off)
    _ny_t = math.sin(_brg_off)
    # C11a: keep F11a fixed using baseline R_a13 and F14_N (original design)
    _C11_N = fp_pts["C7"][1]
    _F14_N_ref = fp_pts["F1"][1] + WALL_OUTER + IW1_OFFSET_FROM_W1 + WALL_6IN + (SOUTH_WALL_FACE + WALL_OUTER)
    _R_a13_base = ARC_F13_R_BASELINE
    _C13_E_base = F15_E - _R_a13_base
    _C11_E_ref = _C13_E_base + (_R_a13_base - R_a11 - (_C11_N - _F14_N_ref) * _ny_t) / _nx_t
    _C11a_E = _C11_E_ref - FLAT_SEG_11
    fp_pts["C11a"] = (_C11a_E, _C11_N)
    # F10: 15'2" east of nominal F9 easting
    _nominal_F9_E = fp_pts["C7"][0] + UPPER_E_R + SMALL_ARC_R - WALL_EXTRA
    _F10_E = _nominal_F9_E + F10_OFFSET_FROM_F9
    _F10_N = fp_pts["F9"][1]
    fp_pts["F10"] = (_F10_E, _F10_N)
    # R_a10: tangent to F11-F11a arc. C10 = (F10_E, F10_N + R_a10).
    # |C10 - C11a| = R_a10 + R_a11  =>  solve for R_a10
    _dE = _C11a_E - _F10_E
    _dN_base = _C11_N - _F10_N
    R_a10 = (_dE**2 + _dN_base**2 - R_a11**2) / (2 * (_dN_base + R_a11))
    fp_pts["C10"] = (_F10_E, _F10_N + R_a10)
    # F11: tangent point on line C10 -> C11a
    _dist_cc = R_a10 + R_a11
    _C10_N = _F10_N + R_a10
    fp_pts["F11"] = (_F10_E + R_a10 * (_C11a_E - _F10_E) / _dist_cc,
                     _C10_N + R_a10 * (_C11_N - _C10_N) / _dist_cc)
    fp_pts["F11a"] = (_C11a_E, _C11_N + R_a11)
    # Solve R_a13 from two constraints: F11a-F11b = F11AB_TARGET, F14-F15 = F14_F15_DIST
    _target_E = fp_pts["F11a"][0] + F11AB_TARGET
    _F14_N = fp_pts["F15"][1] + F14_F15_DIST
    R_a13 = (R_a11 + (_C11_N - _F14_N) * _ny_t - (F15_E - _target_E) * _nx_t) / (1 - _nx_t)
    # Place F14, C13, C11, F11b, F12, F13
    _C13_E = F15_E - R_a13
    fp_pts["C13"] = (_C13_E, _F14_N)
    fp_pts["F14"] = (F15_E, _F14_N)
    _C11_E = _C13_E + (R_a13 - R_a11 - (_C11_N - _F14_N) * _ny_t) / _nx_t
    fp_pts["C11"] = (_C11_E, _C11_N)
    fp_pts["F11b"] = (_C11_E, _C11_N + R_a11)
    fp_pts["F13"] = (_C13_E + R_a13 * _nx_t, _F14_N + R_a13 * _ny_t)
    fp_pts["F12"] = (_C11_E + R_a11 * _nx_t, _C11_N + R_a11 * _ny_t)

    return R_a10, R_a11, R_a13, R_a15


def _compute_south_wall(
    fp_pts: dict[str, Point],
) -> tuple[float, float]:
    """South wall: F17-F20, C17, C19. Returns (R_a17, R_a19).

    Depends on F16 already in fp_pts.

    F19, F20 are fixed (layout-determined).  The F16→F17 ray and the
    horizontal through F19 meet at point P at 150°.  The F17-F18 arc
    is the largest fillet of that corner (30° CW sweep, tangent to both
    lines) subject to:
      - F16-F17 ≥ F16_F17_MIN  (5')
      - F18 ≥ F18_OFFSET_FROM_IW4 east of IW4 east face
    """
    _sweep19 = math.atan(1.0 / 9.0)  # arctan(1/9) for F19-F20 arc
    R_a19 = ARC_F19_R

    # IW4 east face
    _iw4_e = (fp_pts["F2"][0] + WALL_OUTER
              + APPLIANCE_OFFSET_FROM_W2 + APPLIANCE_WIDTH + COUNTER_GAP + COUNTER_DEPTH
              + WALL_3IN + CLOSET_WIDTH + WALL_4IN      # closet 1 (IW7-IW3-IW9)
              + BEDROOM_WIDTH + WALL_4IN + CLOSET2_WIDTH  # bedroom + closet 2
              + WALL_4IN)                                 # IW4 thickness

    # ── F19, F20 (fixed, independent of F17-F18 fillet) ──
    F19_E = _iw4_e + F19_OFFSET_FROM_IW4
    fp_pts["F19"] = (F19_E, SOUTH_WALL_FACE)
    fp_pts["C19"] = (F19_E, SOUTH_WALL_FACE + R_a19)
    _theta_f20 = -math.pi / 2 - _sweep19  # F19 at -π/2, sweep CW
    fp_pts["F20"] = (fp_pts["C19"][0] + R_a19 * math.cos(_theta_f20),
                     fp_pts["C19"][1] + R_a19 * math.sin(_theta_f20))

    # ── F17-F18 fillet ──
    # Intersection P of F16→F17 ray with horizontal at SOUTH_WALL_FACE
    _brg = math.radians(PIX_PI5_TARGET_BRG)
    _sin_b, _cos_b = math.sin(_brg), math.cos(_brg)
    _t_P = (fp_pts["F16"][1] - SOUTH_WALL_FACE) / _cos_b   # parametric dist F16→P
    _P_E = fp_pts["F16"][0] - _t_P * _sin_b

    # Tangency requires sweep + bearing = 90°
    assert abs(ARC_F17_SWEEP + PIX_PI5_TARGET_BRG - 90.0) < 1e-9, \
        f"sweep {ARC_F17_SWEEP}° + bearing {PIX_PI5_TARGET_BRG}° != 90°"

    # Fillet tangent length d = R / tan(half_angle), half_angle = (180° - sweep)/2
    _sw_rad = math.radians(ARC_F17_SWEEP)
    _tan_ha = math.tan((math.pi - _sw_rad) / 2)

    # Max d from each constraint (larger d → larger R → shorter F16-F17, F18 more west)
    F18_min_E = _iw4_e + F18_OFFSET_FROM_IW4
    _d_max_seg = _t_P - F16_F17_MIN          # F16-F17 ≥ 5'
    _d_max_F18 = _P_E - F18_min_E            # F18_E ≥ min
    _d = min(_d_max_seg, _d_max_F18)
    assert _d > 0, f"Cannot fit fillet: d_seg={_d_max_seg:.4f}, d_F18={_d_max_F18:.4f}"
    R_a17 = _d * _tan_ha

    # F17: at tangent length d from P, back toward F16 on the ray
    fp_pts["F17"] = (_P_E + _d * _sin_b, SOUTH_WALL_FACE + _d * _cos_b)
    # F18: at tangent length d from P, west on horizontal
    F18_E = _P_E - _d
    fp_pts["F18"] = (F18_E, SOUTH_WALL_FACE)
    # C17: directly above F18 (perpendicular to horizontal → tangent at F18)
    fp_pts["C17"] = (F18_E, SOUTH_WALL_FACE + R_a17)

    return R_a17, R_a19


# ============================================================
# Main entry point
# ============================================================

def compute_outline_geometry(anchors: OutlineAnchors) -> OutlineGeometry:
    """Compute F-series outline. Chain walk is single source of truth."""
    fp_pts = walk_outline_chain()
    outline_segs = _build_outline_segs()
    radii = _build_radii()
    return OutlineGeometry(fp_pts=fp_pts, outline_segs=outline_segs, radii=radii)
