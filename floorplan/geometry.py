"""Compute F-series outline geometry from inset anchor points and design constants."""
import math
from typing import NamedTuple


from shared.types import Point, LineSeg, ArcSeg, Segment
from shared.geometry import left_norm, off_pt, poly_area
from floorplan.constants import (
    CORNER_NE_R, CORNER_NW_R, UPPER_E_R, SMALL_ARC_R, ARC_180_R,
    R_a2_a3_DELTA, F6_HEIGHT, NW_SHIFT, F1_F2_TARGET, F4_F5_DROP,
    F16_F17_SEG, F14_F15_SEG, F9_F10_DIST, F13_EXIT_BRG,
    SOUTH_WALL_N, PIX_PI5_TARGET_BRG, F15_OFFSET_E,
    WALL_OUTER, WALL_6IN, WALL_3IN, WALL_4IN,
    APPLIANCE_WIDTH, COUNTER_GAP, COUNTER_DEPTH,
    CLOSET_WIDTH, BEDROOM_WIDTH, APPLIANCE_OFFSET_E,
    IW1_OFFSET_N, WALL_SOUTH_N,
)


class OutlineAnchors(NamedTuple):
    """Anchor points from inset path computation."""
    Pi2: Point       # NW anchor (Po2 alias)
    Pi3: Point       # NE anchor (sets F1 easting)
    Ti3: Point       # F21 northing anchor
    PiX: Point       # south wall line start
    Pi5: Point       # south wall line end
    TC1: Point       # Arc 1 center (for west_E)
    R1i: float       # Arc 1 inset radius


class OutlineGeometry(NamedTuple):
    """Complete outline geometry result."""
    fp_pts: dict[str, Point]     # F0-F21 + C0-C20
    outline_segs: list[Segment]  # 22 segments with F-series names
    radii: dict[str, float]      # R_a0 through R_a20


# ============================================================
# Per-section helpers (mutate fp_pts, return computed radii)
# ============================================================

def _compute_ne_corner(fp_pts: dict[str, Point], anchors: OutlineAnchors) -> float:
    """NE corner arc: F0, F1, C0. Returns R_a0."""
    R_a0 = CORNER_NE_R
    fp_pts["C0"] = (anchors.Pi3[0] + R_a0, anchors.Pi3[1] + R_a0)
    fp_pts["F1"] = (anchors.Pi3[0], fp_pts["C0"][1])
    fp_pts["F0"] = (fp_pts["C0"][0], anchors.Pi3[1])
    return R_a0


def _compute_nw_corner(fp_pts: dict[str, Point], anchors: OutlineAnchors) -> float:
    """NW corner arc: F4, F5, F6, C5. Returns R_a5."""
    R_a5 = CORNER_NW_R
    corner2_N = fp_pts["F0"][1] + F6_HEIGHT
    fp_pts["C5"] = (anchors.Pi2[0] + NW_SHIFT + R_a5, corner2_N - R_a5)
    fp_pts["F5"] = (anchors.Pi2[0] + NW_SHIFT, fp_pts["C5"][1])
    fp_pts["F6"] = (fp_pts["C5"][0], corner2_N)
    fp_pts["F4"] = (fp_pts["F5"][0], fp_pts["F5"][1] - F4_F5_DROP + CORNER_NW_R)
    return R_a5


def _compute_east_wall_arcs(fp_pts: dict[str, Point]) -> tuple[float, float, float, float]:
    """East wall arcs: F2, F3, C2, C3, F7, F8, C7, F9, C8. Returns (R_a2, R_a3, R_a7, R_a8).

    Depends on F1, F4, F6 already in fp_pts.
    """
    # Arcs F2-F3-F4: R_a2 = R_a3 + 8/12, F1-F2 = 16'8"
    _a23 = fp_pts["F1"][0] - fp_pts["F4"][0]
    _K23 = (fp_pts["F4"][1] - fp_pts["F1"][1]) - F1_F2_TARGET
    _S23 = -(_K23**2 + _a23**2) / (2 * _a23)
    R_a3 = (_S23 - R_a2_a3_DELTA) / 2
    R_a2 = R_a3 + R_a2_a3_DELTA
    fp_pts["C3"] = (fp_pts["F4"][0] - R_a3, fp_pts["F4"][1])
    _cc2_E = fp_pts["F1"][0] + R_a2
    _dE_cc = _cc2_E - fp_pts["C3"][0]
    _dN_cc = -math.sqrt((R_a3 + R_a2)**2 - _dE_cc**2)
    fp_pts["F2"] = (fp_pts["F1"][0], fp_pts["C3"][1] + _dN_cc)
    fp_pts["C2"] = (_cc2_E, fp_pts["F2"][1])
    # F3: tangent point on C3->C2 line
    _f1b = R_a3 / (R_a3 + R_a2)
    fp_pts["F3"] = (fp_pts["C3"][0] + _f1b * (fp_pts["C2"][0] - fp_pts["C3"][0]),
                    fp_pts["C3"][1] + _f1b * (fp_pts["C2"][1] - fp_pts["C3"][1]))

    # F7: east of F6, arc C7 (R=28")
    fp_pts["F7"] = (fp_pts["F6"][0] + 5.5 + 6.0/12 - NW_SHIFT + 4.0/12, fp_pts["F6"][1])
    R_a7 = UPPER_E_R
    fp_pts["C7"] = (fp_pts["F7"][0], fp_pts["F7"][1] - R_a7)
    fp_pts["F8"] = (fp_pts["C7"][0] + R_a7, fp_pts["C7"][1])

    # F8-F9: small arc (R=2")
    R_a8 = SMALL_ARC_R
    fp_pts["C8"] = (fp_pts["F8"][0] + R_a8, fp_pts["F8"][1])
    fp_pts["F9"] = (fp_pts["C8"][0], fp_pts["C8"][1] - R_a8)

    return R_a2, R_a3, R_a7, R_a8


def _compute_central_region(
    fp_pts: dict[str, Point], anchors: OutlineAnchors,
) -> tuple[float, float, float, float]:
    """Central region: F10-F16, C10, C11, C13, C15. Returns (R_a10, R_a11, R_a13, R_a15).

    Depends on F0, F1, F9, F16 already in fp_pts.
    """
    R_a11 = ARC_180_R
    R_a10 = SMALL_ARC_R

    # F14 Northing from IW1 constraint
    _iw1_n_face = fp_pts["F0"][1] + WALL_OUTER + IW1_OFFSET_N + WALL_6IN
    _F14_N = _iw1_n_face + WALL_SOUTH_N

    # Arc at Po5 corner (exits North)
    d_in_po5 = (anchors.Pi5[0] - anchors.PiX[0], anchors.Pi5[1] - anchors.PiX[1])
    L_in = math.sqrt(d_in_po5[0]**2 + d_in_po5[1]**2)
    d_in_u = (d_in_po5[0]/L_in, d_in_po5[1]/L_in)
    # F15 E-coordinate
    _iw8_e = (fp_pts["F1"][0] + WALL_OUTER + APPLIANCE_OFFSET_E + APPLIANCE_WIDTH + COUNTER_GAP
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

    # F13-F14 arc: bearing F13->F12 = 345°, R_a13 derived from F9-F10 distance
    _brg_off = math.radians(360.0 - F13_EXIT_BRG)
    _nx_t = math.cos(_brg_off)
    _ny_t = math.sin(_brg_off)
    _inv_nx = 1.0 / _nx_t
    _dN = (fp_pts["F9"][1] + SMALL_ARC_R) - _F14_N
    _target_F10_E = fp_pts["F9"][0] + F9_F10_DIST
    R_a13 = (_target_F10_E - F15_E + R_a11 * (_inv_nx + 1)
             + _dN * _ny_t * _inv_nx + R_a10) / (_inv_nx - 1)
    fp_pts["C13"] = (F15_E - R_a13, _F14_N)
    fp_pts["F14"] = (F15_E, _F14_N)
    # C11 easting from tangent constraint
    _C11_N = fp_pts["F9"][1] + R_a10
    _C13_E, _C13_N = fp_pts["C13"]
    _C11_E = _C13_E + (R_a13 - R_a11 - (_C11_N - _C13_N) * _ny_t) / _nx_t
    _corner_E = _C11_E - R_a11
    fp_pts["F10"] = (_corner_E - R_a10, fp_pts["F9"][1])
    fp_pts["C10"] = (_corner_E - R_a10, fp_pts["F9"][1] + R_a10)
    fp_pts["F11"] = (_corner_E, fp_pts["F9"][1] + R_a10)
    fp_pts["C11"] = (_C11_E, _C11_N)
    # F13, F12: tangent points
    fp_pts["F13"] = (fp_pts["C13"][0] + R_a13 * _nx_t, fp_pts["C13"][1] + R_a13 * _ny_t)
    fp_pts["F12"] = (fp_pts["C11"][0] + R_a11 * _nx_t, fp_pts["C11"][1] + R_a11 * _ny_t)

    return R_a10, R_a11, R_a13, R_a15


def _compute_south_wall(
    fp_pts: dict[str, Point], anchors: OutlineAnchors, R_a2: float, R_a3: float,
) -> tuple[float, float, float]:
    """South wall: F17-F21, C17, C19, C20. Returns (R_a17, R_a19, R_a20).

    Depends on F1, F16 already in fp_pts.
    """
    R_a20 = R_a3
    R_a19 = R_a2
    dN_c = (SOUTH_WALL_N + R_a19) - (anchors.Ti3[1] - R_a20)
    dE_c = math.sqrt((R_a20 + R_a19)**2 - dN_c**2)
    # Align F19 with east side of king bed
    _bed_e_align = fp_pts["F1"][0] + WALL_OUTER + 20.5
    fp_pts["F21"] = (_bed_e_align - dE_c - 2.0/12, anchors.Ti3[1])

    # F17 on line from F16 at bearing 60°
    _brg_13 = math.radians(PIX_PI5_TARGET_BRG)
    _sin_b = math.sin(_brg_13)
    _cos_b = math.cos(_brg_13)
    R_a17 = (fp_pts["F16"][1] - SOUTH_WALL_N - F16_F17_SEG * _cos_b) / (1.0 - _sin_b)
    F17_N = SOUTH_WALL_N + R_a17 * (1.0 - _sin_b)
    _t_13 = (fp_pts["F16"][1] - F17_N) / _cos_b
    fp_pts["F17"] = (fp_pts["F16"][0] - _t_13 * _sin_b, F17_N)
    fp_pts["C17"] = (fp_pts["F17"][0] - R_a17 * _cos_b, SOUTH_WALL_N + R_a17)
    fp_pts["F18"] = (fp_pts["C17"][0], SOUTH_WALL_N)
    # C20, C19, F19, F20
    fp_pts["C20"] = (fp_pts["F21"][0], fp_pts["F21"][1] - R_a20)
    F19_E = fp_pts["F21"][0] + dE_c
    fp_pts["F19"] = (F19_E, SOUTH_WALL_N)
    fp_pts["C19"] = (F19_E, SOUTH_WALL_N + R_a19)
    _f_w = R_a20 / (R_a20 + R_a19)
    fp_pts["F20"] = (fp_pts["C20"][0] + _f_w * (fp_pts["C19"][0] - fp_pts["C20"][0]),
                     fp_pts["C20"][1] + _f_w * (fp_pts["C19"][1] - fp_pts["C20"][1]))

    return R_a17, R_a19, R_a20


# ============================================================
# Main entry point
# ============================================================

def compute_outline_geometry(anchors: OutlineAnchors) -> OutlineGeometry:
    """Compute F-series outline from inset anchor points + design constants."""
    fp_pts: dict[str, Point] = {}

    R_a0 = _compute_ne_corner(fp_pts, anchors)
    R_a5 = _compute_nw_corner(fp_pts, anchors)
    R_a2, R_a3, R_a7, R_a8 = _compute_east_wall_arcs(fp_pts)
    R_a10, R_a11, R_a13, R_a15 = _compute_central_region(fp_pts, anchors)
    R_a17, R_a19, R_a20 = _compute_south_wall(fp_pts, anchors, R_a2, R_a3)

    # --- Build outline segments (F-series) ---
    outline_segs: list[Segment] = [
        ArcSeg("F0", "F1", "C0", R_a0, "CW", 20),
        LineSeg("F1", "F2"),
        ArcSeg("F2", "F3", "C2", R_a2, "CW", 20),
        ArcSeg("F3", "F4", "C3", R_a3, "CCW", 20),
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
        ArcSeg("F20", "F21", "C20", R_a20, "CCW", 60),
        LineSeg("F21", "F0"),
    ]

    radii = {
        "R_a0": R_a0, "R_a2": R_a2, "R_a3": R_a3, "R_a5": R_a5,
        "R_a7": R_a7, "R_a8": R_a8, "R_a10": R_a10, "R_a11": R_a11,
        "R_a13": R_a13, "R_a15": R_a15, "R_a17": R_a17,
        "R_a19": R_a19, "R_a20": R_a20,
    }

    return OutlineGeometry(fp_pts=fp_pts, outline_segs=outline_segs, radii=radii)
