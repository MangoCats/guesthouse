"""Compute F-series outline geometry from chain walk and design constants."""
import math
from typing import NamedTuple


from shared.types import Point, LineSeg, ArcSeg, Segment
from shared.survey import COORD_ROTATION, rotate_pts
from floorplan.constants import CORNER_SW_R, F11AB_TARGET


class OutlineGeometry(NamedTuple):
    """Complete outline geometry result."""
    fp_pts: dict[str, Point]     # F1-F18, F11a, F11b + C1-C17, C11a + FC
    outline_segs: list[Segment]  # 18 segments with F-series names
    radii: dict[str, float]      # R_a1 through R_a17


# ============================================================
# F-series outline chain: single source of truth
# ============================================================

# Sweep angle constants (radians)
_PI_2 = math.pi / 2            # 90 deg
_5PI_12 = 5 * math.pi / 12    # 75 deg
_PI_12 = math.pi / 12          # 15 deg
_C15_SWEEP = math.pi / 2 - math.atan(7.0 / 12.0)  # F15→F16: π/2 − arctan(7/12)
_C17_SWEEP = math.atan(7.0 / 12.0)                 # F17→F18: arctan(7/12)
_C10_SWEEP = math.pi / 2 - math.atan(1.0 / 3.0)  # F10→F11: 90° - arctan(1/3)
_C11_SWEEP = _C10_SWEEP                           # F11b→F12: same as C10
_C13_SWEEP = math.atan(1.0 / 3.0)                 # F13→F14: arctan(1/3) → exit bearing = π

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


# Segment definitions — all parameters are explicit constants.
# Only d_F2_F5 and d_F18_F1 are computed (closure solver below).
_CHAIN_F5_TO_F18 = [
    ("CW",   2.333333333333, _PI_2, "C5", 20),                # F5->F6
    ("L",    5.166666666667),                                  # F6->F7
    ("CW",   2.333333333333, _PI_2, "C7", 20),                # F7->F8
    ("CCW",  0.166666666667, _PI_2, "C8", 20),                # F8->F9
    ("L",   15.500000000000),                                  # F9->F10
    ("CCW",  1.322854905602, _C10_SWEEP, "C10", 20),          # F10->F11
    ("CW",   2.333333333333, _C10_SWEEP, "C11a", 30),         # F11->F11a
    ("L",    F11AB_TARGET),                                    # F11a->F11b
    ("CW",   2.333333333333, _C11_SWEEP, "C11", 30),          # F11b->F12
    ("L",   11.779557008578),                                  # F12->F13
    ("CW",   1.808727374505, _C13_SWEEP, "C13", 60),          # F13->F14
    ("L",    9.662743475808),                                  # F14->F15
    ("CW",   1.808727374505, _C15_SWEEP, "C15", 20),          # F15->F16
    ("L",    5.000000000000),                                  # F16->F17
    ("CW",   1.808727374505, _C17_SWEEP, "C17", 20),          # F17->F18
]

# Closure solver: compute d_F2_F5 and d_F18_F1 so the chain closes.
_R_a1 = CORNER_SW_R
_dE_18, _dN_18, _brg_18 = _chain_offset(_CHAIN_F5_TO_F18, start_brg=0.0)
_d_F18_F1 = (_R_a1 - _dE_18) / math.sin(_brg_18)
_F1_N_rel = _dN_18 + _d_F18_F1 * math.cos(_brg_18)
_d_F2_F5 = -(_F1_N_rel + _R_a1)

# F2 position: east face at E = -18'6", F1.N = -13'6" exactly.
# F2.N derived from F1.N + CORNER_SW_R (F1→F2 is 90° CW arc of radius R_a1).
F2_E = -18.5                     # -18'6" east face easting
F2_N = -13.5 + CORNER_SW_R      # F1.N + R_a1
F2_BRG = 0.0                     # bearing 0 (due north)

# Full outline chain: F2→F5 + F5→F18 + F18→F1 + F1→F2 arc = 18 entries.
# Starting at F2, bearing 0 (due north), CW traversal.
OUTLINE_CHAIN = [
    ("L",   _d_F2_F5),                                        # F2->F5
] + _CHAIN_F5_TO_F18 + [
    ("L",   _d_F18_F1),                                       # F18->F1
    ("CW",  _R_a1, _PI_2, "C1", 20),                          # F1->F2
]

# Point names produced by each chain segment (one per segment, in order)
CHAIN_POINT_NAMES = [
    "F5", "F6", "F7", "F8", "F9", "F10", "F11",
    "F11a", "F11b", "F12", "F13", "F14", "F15", "F16", "F17",
    "F18", "F1", "F2",
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
    start_names = ["F2"] + CHAIN_POINT_NAMES[:-1]  # F2, F5, ..., F18, F1
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
