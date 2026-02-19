"""Compute F-series outline geometry from inset anchor points and design constants."""
import math
from typing import NamedTuple


from shared.types import Point, LineSeg, ArcSeg, Segment
from shared.geometry import left_norm, off_pt, poly_area
from floorplan.constants import (
    CORNER_NE_R, CORNER_NW_R, UPPER_E_R, SMALL_ARC_R, ARC_180_R,
    ARC_F3_R, ARC_F3_SWEEP, F6_EAST_ADJ,
    F6_HEIGHT, NW_SHIFT,
    F14_F15_SEG, ARC_F13_R, F13_EXIT_BRG,
    SOUTH_WALL_N, PIX_PI5_TARGET_BRG, F15_OFFSET_E, F16_F17_SEG,
    F18_OFFSET_E, F18_F19_GAP, ARC_F19_R,
    WALL_OUTER, WALL_6IN, WALL_3IN, WALL_4IN,
    APPLIANCE_WIDTH, COUNTER_GAP, COUNTER_DEPTH,
    CLOSET_WIDTH, CLOSET2_WIDTH, BEDROOM_WIDTH, APPLIANCE_OFFSET_E,
    IW1_OFFSET_N, WALL_SOUTH_N,
    O6_E_FROM_F9, F10_O6_CLEARANCE,
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
    fp_pts: dict[str, Point]     # F1-F20 + C1-C19
    outline_segs: list[Segment]  # 20 segments with F-series names
    radii: dict[str, float]      # R_a1 through R_a19


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
    fp_pts["C5"] = (anchors.Pi2[0] + NW_SHIFT + R_a5 + F6_EAST_ADJ, corner2_N - R_a5)
    fp_pts["F6"] = (fp_pts["C5"][0], corner2_N)
    return R_a5


def _compute_west_wall(fp_pts: dict[str, Point]) -> float:
    """West wall: F3, F4, F5, C3. Returns R_a3.

    F3 northing solved from fixed sweep angle ARC_F3_SWEEP and tangent
    constraint (F4-F5 tangent to both arcs F3-F4 and F5-F6).
    Depends on F2, F6, C5 already in fp_pts.
    """
    R_a3 = ARC_F3_R
    R_a5 = CORNER_NW_R

    F3_E = fp_pts["F2"][0]
    C5_E, C5_N = fp_pts["C5"]

    # Fixed sweep angle γ = ARC_F3_SWEEP
    _gamma = math.radians(ARC_F3_SWEEP)

    # Solve F3_N from tangent constraint:
    # (C5_N - F3_N)·sin(γ) - (C5_E - F3_E - R_a3)·cos(γ) = R_a3 - R_a5
    _dE_fixed = C5_E - F3_E - R_a3
    F3_N = C5_N - (R_a3 - R_a5 + _dE_fixed * math.cos(_gamma)) / math.sin(_gamma)

    fp_pts["F3"] = (F3_E, F3_N)
    fp_pts["C3"] = (F3_E + R_a3, F3_N)

    # F4 and F5 at angle (π - γ) from their respective centers
    _angle = math.pi - _gamma
    C3_E, C3_N = fp_pts["C3"]
    fp_pts["F4"] = (C3_E + R_a3 * math.cos(_angle),
                    C3_N + R_a3 * math.sin(_angle))
    fp_pts["F5"] = (C5_E + R_a5 * math.cos(_angle),
                    C5_N + R_a5 * math.sin(_angle))

    return R_a3


def _compute_upper_east_arcs(fp_pts: dict[str, Point]) -> tuple[float, float]:
    """Upper east arcs: F7, F8, C7, F9, C8. Returns (R_a7, R_a8).

    Depends on F6 already in fp_pts.
    """
    # F7: east of F6, adjusted for F6 east shift to keep F7 fixed
    fp_pts["F7"] = (fp_pts["F6"][0] + 5.5 + 6.0/12 - NW_SHIFT + 4.0/12 - F6_EAST_ADJ,
                    fp_pts["F6"][1])
    R_a7 = UPPER_E_R
    fp_pts["C7"] = (fp_pts["F7"][0], fp_pts["F7"][1] - R_a7)
    fp_pts["F8"] = (fp_pts["C7"][0] + R_a7, fp_pts["C7"][1])

    # F8-F9: small arc (R=2")
    R_a8 = SMALL_ARC_R
    fp_pts["C8"] = (fp_pts["F8"][0] + R_a8, fp_pts["F8"][1])
    fp_pts["F9"] = (fp_pts["C8"][0], fp_pts["C8"][1] - R_a8)

    return R_a7, R_a8


def _compute_central_region(
    fp_pts: dict[str, Point], anchors: OutlineAnchors,
) -> tuple[float, float, float, float]:
    """Central region: F10-F16, C10, C11, C13, C15. Returns (R_a10, R_a11, R_a13, R_a15).

    Depends on F0, F2, F9, F16 already in fp_pts.
    """
    R_a11 = ARC_180_R
    R_a13 = ARC_F13_R

    # F14 Northing from IW1 constraint
    _iw1_n_face = fp_pts["F1"][1] + WALL_OUTER + IW1_OFFSET_N + WALL_6IN
    _F14_N = _iw1_n_face + WALL_SOUTH_N

    # Arc at Po5 corner (exits North)
    d_in_po5 = (anchors.Pi5[0] - anchors.PiX[0], anchors.Pi5[1] - anchors.PiX[1])
    L_in = math.sqrt(d_in_po5[0]**2 + d_in_po5[1]**2)
    d_in_u = (d_in_po5[0]/L_in, d_in_po5[1]/L_in)
    # F15 E-coordinate
    _iw8_e = (fp_pts["F2"][0] + WALL_OUTER + APPLIANCE_OFFSET_E + APPLIANCE_WIDTH + COUNTER_GAP
              + COUNTER_DEPTH + (WALL_3IN + CLOSET_WIDTH + WALL_4IN
              + BEDROOM_WIDTH + WALL_4IN + CLOSET_WIDTH + WALL_3IN))
    F15_E = _iw8_e + F15_OFFSET_E
    ln_in_po5 = left_norm(anchors.PiX, anchors.Pi5)

    # R_a15 from constraint: F14-F15 segment length
    _A_15 = anchors.PiX[1] + (F15_E - anchors.PiX[0]) * d_in_u[1] / d_in_u[0]
    _B_15 = ln_in_po5[1] - (1.0 + ln_in_po5[0]) * d_in_u[1] / d_in_u[0]
    R_a15 = (_F14_N - F14_F15_SEG - _A_15) / _B_15
    o_in_po5 = off_pt(anchors.PiX, ln_in_po5, R_a15)
    t_cf4 = (F15_E - R_a15 - o_in_po5[0]) / d_in_u[0]
    fp_pts["C15"] = (F15_E - R_a15, o_in_po5[1] + t_cf4 * d_in_u[1])
    fp_pts["F15"] = (F15_E, fp_pts["C15"][1])
    # F16: tangent point for 60° incoming bearing
    _brg_f4 = math.radians(PIX_PI5_TARGET_BRG)
    fp_pts["F16"] = (fp_pts["C15"][0] + R_a15 * math.cos(_brg_f4),
                     fp_pts["C15"][1] - R_a15 * math.sin(_brg_f4))

    # F13-F14 arc: R_a13 is a fixed constant, bearing F13->F12 = 345°
    _brg_off = math.radians(360.0 - F13_EXIT_BRG)
    _nx_t = math.cos(_brg_off)
    _ny_t = math.sin(_brg_off)
    fp_pts["C13"] = (F15_E - R_a13, _F14_N)
    fp_pts["F14"] = (F15_E, _F14_N)
    # C11 from tangent constraint with C13 (position unchanged)
    _C11_N = fp_pts["F9"][1] + SMALL_ARC_R
    _C13_E, _C13_N = fp_pts["C13"]
    _C11_E = _C13_E + (R_a13 - R_a11 - (_C11_N - _C13_N) * _ny_t) / _nx_t
    fp_pts["C11"] = (_C11_E, _C11_N)
    # F10: 4" east of O6 east edge (O6 east = F9 + O6_E_FROM_F9)
    _F10_E = fp_pts["F9"][0] + O6_E_FROM_F9 + F10_O6_CLEARANCE
    _F10_N = fp_pts["F9"][1]
    fp_pts["F10"] = (_F10_E, _F10_N)
    # R_a10: tangent to F11-F12 arc. C10 = (F10_E, F10_N + R_a10).
    # |C10 - C11| = R_a10 + R_a11  =>  solve for R_a10
    _dE = _C11_E - _F10_E
    _dN_base = _C11_N - _F10_N
    R_a10 = (_dE**2 + _dN_base**2 - R_a11**2) / (2 * (_dN_base + R_a11))
    fp_pts["C10"] = (_F10_E, _F10_N + R_a10)
    # F11: tangent point on line C10 -> C11
    _dist_cc = R_a10 + R_a11
    _C10_N = _F10_N + R_a10
    fp_pts["F11"] = (_F10_E + R_a10 * (_C11_E - _F10_E) / _dist_cc,
                     _C10_N + R_a10 * (_C11_N - _C10_N) / _dist_cc)
    # F13, F12: tangent points (unchanged — C11, C13, R_a11, R_a13 unchanged)
    fp_pts["F13"] = (fp_pts["C13"][0] + R_a13 * _nx_t, fp_pts["C13"][1] + R_a13 * _ny_t)
    fp_pts["F12"] = (fp_pts["C11"][0] + R_a11 * _nx_t, fp_pts["C11"][1] + R_a11 * _ny_t)

    return R_a10, R_a11, R_a13, R_a15


def _compute_south_wall(
    fp_pts: dict[str, Point],
) -> tuple[float, float]:
    """South wall: F17-F20, C17, C19. Returns (R_a17, R_a19).

    Depends on F16 already in fp_pts.
    F19-F20: CW arc (C19). F20 exit bearing used by caller for tangency to F0.
    """
    _sweep = math.asin(1.0 / 9.0)  # arcsin(1/9)
    R_a19 = ARC_F19_R

    # F18: 4" east of IW4 east face, at SOUTH_WALL_N
    _iw4_e = (fp_pts["F2"][0] + WALL_OUTER
              + APPLIANCE_OFFSET_E + APPLIANCE_WIDTH + COUNTER_GAP + COUNTER_DEPTH
              + WALL_3IN + CLOSET_WIDTH + WALL_4IN      # closet 1 (IW7-IW3-IW9)
              + BEDROOM_WIDTH + WALL_4IN + CLOSET2_WIDTH  # bedroom + closet 2
              + WALL_4IN)                                 # IW4 thickness
    F18_E = _iw4_e + F18_OFFSET_E
    fp_pts["F18"] = (F18_E, SOUTH_WALL_N)

    # F17: fixed 5' from F16 along PiX-Pi5 bearing
    _brg = math.radians(PIX_PI5_TARGET_BRG)
    _sin_b = math.sin(_brg)
    _cos_b = math.cos(_brg)
    fp_pts["F17"] = (fp_pts["F16"][0] - F16_F17_SEG * _sin_b,
                     fp_pts["F16"][1] - F16_F17_SEG * _cos_b)

    # R_a17: arc from F17 to F18, center C17 directly above F18
    _dx17 = F18_E - fp_pts["F17"][0]
    _dy17 = SOUTH_WALL_N - fp_pts["F17"][1]
    R_a17 = -(_dx17**2 + _dy17**2) / (2 * _dy17)
    fp_pts["C17"] = (F18_E, SOUTH_WALL_N + R_a17)
    # F19: 12" west of F18
    F19_E = fp_pts["F18"][0] - F18_F19_GAP
    fp_pts["F19"] = (F19_E, SOUTH_WALL_N)
    fp_pts["C19"] = (F19_E, SOUTH_WALL_N + R_a19)
    # F20: CW arc from F19 with sweep = arcsin(1/9)
    _theta_f20 = -math.pi / 2 - _sweep  # F19 at -π/2, sweep CW
    fp_pts["F20"] = (fp_pts["C19"][0] + R_a19 * math.cos(_theta_f20),
                     fp_pts["C19"][1] + R_a19 * math.sin(_theta_f20))

    return R_a17, R_a19


# ============================================================
# Main entry point
# ============================================================

def compute_outline_geometry(anchors: OutlineAnchors) -> OutlineGeometry:
    """Compute F-series outline from inset anchor points + design constants."""
    fp_pts: dict[str, Point] = {}

    R_a1 = _compute_ne_corner(fp_pts, anchors)
    R_a5 = _compute_nw_corner(fp_pts, anchors)
    R_a3 = _compute_west_wall(fp_pts)
    R_a7, R_a8 = _compute_upper_east_arcs(fp_pts)
    R_a10, R_a11, R_a13, R_a15 = _compute_central_region(fp_pts, anchors)
    R_a17, R_a19 = _compute_south_wall(fp_pts)

    # F1: tangent to line from F20 on arc C1 (CW, radius R_a1, center at F2 northing)
    # Exit bearing from F20 (same sweep as C19 arc)
    _sweep = math.asin(1.0 / 9.0)
    _exit_angle = -math.pi - _sweep  # CW tangent at F20: radius_angle - π/2
    _ex = math.cos(_exit_angle)
    _ey = math.sin(_exit_angle)
    # Solve for F2_N from tangency: F1 on line from F20 and on circle C1
    _F2_E = fp_pts["F2"][0]
    _F20_E, _F20_N = fp_pts["F20"]
    _dE = _F2_E - _F20_E
    _F2_N = _F20_N + (_dE + R_a1 * (1 - _ey)) * _ey / _ex - R_a1 * _ex
    fp_pts["F2"] = (_F2_E, _F2_N)
    fp_pts["C1"] = (_F2_E + R_a1, _F2_N)
    fp_pts["F1"] = (_F2_E + R_a1 * (1 - _ey), _F2_N + R_a1 * _ex)

    # --- Build outline segments (F-series) ---
    outline_segs: list[Segment] = [
        ArcSeg("F1", "F2", "C1", R_a1, "CW", 20),
        LineSeg("F2", "F3"),
        ArcSeg("F3", "F4", "C3", R_a3, "CW", 20),
        LineSeg("F4", "F5"),
        ArcSeg("F5", "F6", "C5", R_a5, "CW", 20),
        LineSeg("F6", "F7"),
        ArcSeg("F7", "F8", "C7", R_a7, "CW", 20),
        ArcSeg("F8", "F9", "C8", R_a8, "CCW", 20),
        LineSeg("F9", "F10"),
        ArcSeg("F10", "F11", "C10", R_a10, "CCW", 20),
        ArcSeg("F11", "F12", "C11", R_a11, "CW", 60),
        LineSeg("F12", "F13"),
        ArcSeg("F13", "F14", "C13", R_a13, "CW", 60),
        LineSeg("F14", "F15"),
        ArcSeg("F15", "F16", "C15", R_a15, "CW", 20),
        LineSeg("F16", "F17"),
        ArcSeg("F17", "F18", "C17", R_a17, "CW", 20),
        LineSeg("F18", "F19"),
        ArcSeg("F19", "F20", "C19", R_a19, "CW", 60),
        LineSeg("F20", "F1"),
    ]

    radii = {
        "R_a1": R_a1, "R_a3": R_a3, "R_a5": R_a5,
        "R_a7": R_a7, "R_a8": R_a8, "R_a10": R_a10, "R_a11": R_a11,
        "R_a13": R_a13, "R_a15": R_a15, "R_a17": R_a17,
        "R_a19": R_a19,
    }

    return OutlineGeometry(fp_pts=fp_pts, outline_segs=outline_segs, radii=radii)
