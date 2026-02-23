"""Compute F-series outline geometry from chain walk and design constants."""
import math
from typing import NamedTuple


from shared.types import Point, LineSeg, ArcSeg, Segment
from shared.survey import COORD_ROTATION, rotate_pts
from floorplan.constants import CORNER_NE_R


class OutlineGeometry(NamedTuple):
    """Complete outline geometry result."""
    fp_pts: dict[str, Point]     # F1-F20, F11a, F11b + C1-C19, C11a + FC
    outline_segs: list[Segment]  # 20 segments with F-series names
    radii: dict[str, float]      # R_a1 through R_a19


# ============================================================
# F-series outline chain: single source of truth
# ============================================================

# Sweep angle constants (radians)
_A19 = math.atan(1.0 / 9.0)   # arctan(1/9) for F19-F20
_A9 = math.atan(9.0)           # arctan(9) for F5-F6
_PI_2 = math.pi / 2            # 90 deg
_5PI_12 = 5 * math.pi / 12    # 75 deg
_PI_12 = math.pi / 12          # 15 deg
_PI_3 = math.pi / 3            # 60 deg
_PI_6 = math.pi / 6            # 30 deg

# Chain: ("L", distance) for lines
#        ("CW"/"CCW", radius, sweep, center_name, n_pts) for arcs


def _chain_offset(chain, start_brg=0.0):
    """Walk chain entries from (0,0) and return (delta_E, delta_N, exit_brg)."""
    E, N, brg = 0.0, 0.0, start_brg
    for seg in chain:
        if seg[0] == "L":
            d = seg[1]
            E += d * math.sin(brg)
            N += d * math.cos(brg)
        else:
            direction, R, sweep = seg[0], seg[1], seg[2]
            if direction == "CW":
                cx = E + R * math.cos(brg)
                cy = N - R * math.sin(brg)
                alpha = math.atan2(N - cy, E - cx) - sweep
                E = cx + R * math.cos(alpha)
                N = cy + R * math.sin(alpha)
                brg += sweep
            else:  # CCW
                cx = E - R * math.cos(brg)
                cy = N + R * math.sin(brg)
                alpha = math.atan2(N - cy, E - cx) + sweep
                E = cx + R * math.cos(alpha)
                N = cy + R * math.sin(alpha)
                brg -= sweep
    return E, N, brg


# Fixed chain entries: F5→F6→...→F19→F20 (17 entries, unchanged from original)
_CHAIN_F5_TO_F20 = [
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
]

# Derive NE corner closure constraints from CORNER_NE_R.
# Walk fixed F5→F20 chain to get F20 position relative to F5, then compute
# the F20→F1 distance that places F1 at the correct easting for the 90° arc.
_dE_20, _dN_20, _brg_20 = _chain_offset(_CHAIN_F5_TO_F20, start_brg=0.0)
_R_a1 = CORNER_NE_R
_d_F20_F1 = (_R_a1 - _dE_20) / math.sin(_brg_20)
_F1_N_rel = _dN_20 + _d_F20_F1 * math.cos(_brg_20)
_d_F2_F5 = -(_F1_N_rel + _R_a1)

# Anchor: compute old F5 position to keep F5-F20 at their current locations.
# Old F2 was at pre-rotation (-18.0, -10.5), rotated by arctan(1/9).
# Frozen at arctan(1/9) so F-series position is independent of COORD_ROTATION.
_F_ANCHOR_ROTATION = math.atan(1.0 / 9.0)
_cos_R = math.cos(_F_ANCHOR_ROTATION)
_sin_R = math.sin(_F_ANCHOR_ROTATION)
_F2_E_old = -18.0 * _cos_R - (-10.5) * _sin_R
_F2_N_old = -18.0 * _sin_R + (-10.5) * _cos_R
_F2_BRG_old = -_F_ANCHOR_ROTATION
# Walk old F2→F3→F4→F5 chain entries (now removed) to find F5 anchor position.
_OLD_F2_TO_F5 = [
    ("L",   12.083333333333),                                  # F2->F3 (removed)
    ("CW",   8.351795046046, _A19, "C3", 20),                 # F3->F4 (removed)
    ("L",    9.476667232982),                                  # F4->F5 (removed)
]
_F5_off_E, _F5_off_N, _ = _chain_offset(_OLD_F2_TO_F5, start_brg=_F2_BRG_old)
_F5_E = _F2_E_old + _F5_off_E
_F5_N = _F2_N_old + _F5_off_N

# New F2: directly south of F5, at the base of the east wall.
F2_E = _F5_E                     # same easting as F5
F2_N = _F5_N - _d_F2_F5          # F5.N - line_distance
F2_BRG = 0.0                     # bearing 0 (due north) in rotated frame

# Full outline chain: F2→F5 + fixed F5→F20 + derived F20→F1 + F1→F2 arc = 20 entries.
# Starting at F2, bearing 0 (due north), CW traversal.
OUTLINE_CHAIN = [
    ("L",   _d_F2_F5),                                        # F2->F5
] + _CHAIN_F5_TO_F20 + [
    ("L",   _d_F20_F1),                                       # F20->F1
    ("CW",  _R_a1, _PI_2, "C1", 20),                          # F1->F2
]

# Point names produced by each chain segment (one per segment, in order)
CHAIN_POINT_NAMES = [
    "F5", "F6", "F7", "F8", "F9", "F10", "F11",
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
    start_names = ["F2"] + CHAIN_POINT_NAMES[:-1]  # F2, F5, ..., F20, F1
    end_names = CHAIN_POINT_NAMES                    # F5, F6, ..., F1, F2

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
            ra_name = "R_a" + center_name[1:]  # "C5" -> "R_a5"
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
    """Rotate P/Pi survey points into FC-based rotated coordinate space.

    Pure CCW rotation by COORD_ROTATION around the origin (FC).
    Modifies pts in place.
    """
    rotate_pts(pts, COORD_ROTATION)
