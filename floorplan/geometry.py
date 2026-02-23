"""Compute F-series outline geometry from chain walk and design constants."""
import math
from typing import NamedTuple


from shared.types import Point, LineSeg, ArcSeg, Segment
from shared.survey import COORD_ROTATION, rotate_pts


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
# Pre-rotation values (axis-aligned system where F2_BRG would be 0):
_F2_E0 = -18.0
_F2_N0 = -10.5
# Rotate F2 CCW by COORD_ROTATION so F4-F5 aligns to bearing 0.
_cos_R = math.cos(COORD_ROTATION)
_sin_R = math.sin(COORD_ROTATION)
F2_E = _F2_E0 * _cos_R - _F2_N0 * _sin_R
F2_N = _F2_E0 * _sin_R + _F2_N0 * _cos_R
F2_BRG = -COORD_ROTATION  # initial bearing rotated by -COORD_ROTATION

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
# Main entry point
# ============================================================

def compute_outline_geometry() -> OutlineGeometry:
    """Compute F-series outline. Chain walk is single source of truth."""
    fp_pts = walk_outline_chain()
    outline_segs = _build_outline_segs()
    radii = _build_radii()
    return OutlineGeometry(fp_pts=fp_pts, outline_segs=outline_segs, radii=radii)


def align_pts_to_f_series(pts: dict[str, Point]) -> None:
    """Rigid transform P/Pi survey points into F-series coordinate space.

    Rotation: PiX->Pi5 parallel to F17->F16
    Translation: F16 on PiX-Pi5 line AND F3 on P2-P3 line
    Modifies pts in place.
    """
    fp = walk_outline_chain()
    # Rotation: align PiX->Pi5 direction with F17->F16 direction
    pip = (pts["Pi5"][0] - pts["PiX"][0], pts["Pi5"][1] - pts["PiX"][1])
    f16 = (fp["F16"][0] - fp["F17"][0], fp["F16"][1] - fp["F17"][1])
    rot = math.atan2(f16[1], f16[0]) - math.atan2(pip[1], pip[0])
    rotate_pts(pts, rot)
    # Translation: 2x2 solve from two line-containment constraints
    p23 = (pts["P3"][0] - pts["P2"][0], pts["P3"][1] - pts["P2"][1])
    n1 = (-f16[1], f16[0])          # normal to PiX-Pi5
    n2 = (-p23[1], p23[0])          # normal to P2-P3
    d1 = (fp["F16"][0] - pts["PiX"][0]) * n1[0] + (fp["F16"][1] - pts["PiX"][1]) * n1[1]
    d2 = (fp["F3"][0] - pts["P2"][0]) * n2[0] + (fp["F3"][1] - pts["P2"][1]) * n2[1]
    det = n1[0] * n2[1] - n1[1] * n2[0]
    tx = (d1 * n2[1] - d2 * n1[1]) / det
    ty = (n1[0] * d2 - n2[0] * d1) / det
    for k in list(pts):
        pts[k] = (pts[k][0] + tx, pts[k][1] + ty)
