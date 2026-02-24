"""Compute F-series outline geometry from chain walk and design constants."""
import math
from typing import NamedTuple


from shared.types import Point, LineSeg, ArcSeg, Segment
from shared.survey import COORD_ROTATION, rotate_pts
from floorplan.constants import CORNER_NE_R, F11AB_TARGET


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
_C7_SWEEP = _PI_2 + _A19       # F7→F8: 90° + arctan(1/9)
_C10_SWEEP = _5PI_12 + _A19    # F10→F11: 75° + arctan(1/9)
_C13_SWEEP = _PI_12 + _A19     # F13→F14: 15° + arctan(1/9) → exit bearing = π
_C15_SWEEP = _PI_3 - _A19      # F15→F16: 60° - arctan(1/9)

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


# Unchanged bookend sub-chains
_CHAIN_F5_TO_F7 = [
    ("CW",   2.333333333333, _A9, "C5", 20),                  # F5->F6
    ("L",    5.250000000000),                                  # F6->F7
]
_CHAIN_F17_TO_F20 = [
    ("CW",   6.404672887007, _PI_6, "C17", 20),               # F17->F18
    ("L",    1.397555568554),                                  # F18->F19
    ("CW",  18.888718471469, _A19, "C19", 60),                # F19->F20
]

# Original F7→F14 entries (target displacement, before sweep changes)
_CHAIN_F7_TO_F14_ORIG = [
    ("CW",   2.333333333333, _PI_2, "C7", 20),
    ("CCW",  0.166666666667, _PI_2, "C8", 20),
    ("L",   15.166666666667),
    ("CCW",  1.039662132188, _5PI_12, "C10", 20),
    ("CW",   2.333333333333, _5PI_12, "C11a", 30),
    ("L",    1.000000000000),
    ("CW",   2.333333333333, _5PI_12, "C11", 30),
    ("L",   11.858994000010),
    ("CW",   2.507553207938, _PI_12, "C13", 60),
]

# Solve for R_C10 and L_F12F13 to maintain F7/F14 positions with new sweeps.
_, _, _brg_F7 = _chain_offset(_CHAIN_F5_TO_F7, start_brg=0.0)
_tgt_E, _tgt_N, _ = _chain_offset(_CHAIN_F7_TO_F14_ORIG, start_brg=_brg_F7)

_dE_A, _dN_A, _brg_F10 = _chain_offset([
    ("CW",  2.333333333333, _C7_SWEEP, "C7", 20),
    ("CCW", 0.166666666667, _PI_2, "C8", 20),
    ("L",  15.166666666667),
], start_brg=_brg_F7)

_a_E, _a_N, _brg_post_arc = _chain_offset(
    [("CCW", 1.0, _C10_SWEEP, "C10", 20)], start_brg=_brg_F10)

_b_E, _b_N, _brg_F12 = _chain_offset([
    ("CW",  2.333333333333, _5PI_12, "C11a", 30),
    ("L",   F11AB_TARGET),
    ("CW",  2.333333333333, _5PI_12, "C11", 30),
], start_brg=_brg_post_arc)

_c_E, _c_N, _ = _chain_offset(
    [("CW", 2.507553207938, _PI_12, "C13", 60)], start_brg=_brg_F12)

_rhs_E = _tgt_E - _dE_A - _b_E - _c_E
_rhs_N = _tgt_N - _dN_A - _b_N - _c_N
_s12, _c12 = math.sin(_brg_F12), math.cos(_brg_F12)
_det = _a_E * _c12 - _a_N * _s12
_R_C10 = (_rhs_E * _c12 - _rhs_N * _s12) / _det
_L_F12F13 = (_a_E * _rhs_N - _a_N * _rhs_E) / _det

# Solve for F14→F15 distance and F16→F17 distance to keep F17 fixed.
# Target: original F13→F17 displacement (old C13/C15 sweeps and line distances)
_tgt2_E, _tgt2_N, _ = _chain_offset([
    ("CW",  2.507553207938, _PI_12, "C13", 60),
    ("L",   8.666666666667),
    ("CW",  2.473295271375, _PI_3, "C15", 20),
    ("L",   5.000000000000),
], start_brg=_brg_F12)

# Arc offsets with new sweeps
_c13_dE, _c13_dN, _ = _chain_offset(
    [("CW", 2.507553207938, _C13_SWEEP, "C13", 60)], start_brg=_brg_F12)
_c15_dE, _c15_dN, _ = _chain_offset(
    [("CW", 2.473295271375, _C15_SWEEP, "C15", 20)], start_brg=math.pi)

# 2×2 linear solve: d_F14F15 (line at bearing π) and L_F16F17 (line at bearing π+S)
_rhs2_E = _tgt2_E - _c13_dE - _c15_dE
_rhs2_N = _tgt2_N - _c13_dN - _c15_dN
_sinS = math.sin(_C15_SWEEP)
_cosS = math.cos(_C15_SWEEP)
_d_F14F15 = (_rhs2_E * _cosS - _rhs2_N * _sinS) / _sinS
_L_F16F17 = -_rhs2_E / _sinS

_CHAIN_F14_TO_F20 = [
    ("L",    _d_F14F15),                                       # F14->F15
    ("CW",   2.473295271375, _C15_SWEEP, "C15", 20),          # F15->F16
    ("L",    _L_F16F17),                                       # F16->F17
] + _CHAIN_F17_TO_F20

# Full chain: F5→F20 with computed R_C10 and L_F12F13
_CHAIN_F5_TO_F20 = _CHAIN_F5_TO_F7 + [
    ("CW",   2.333333333333, _C7_SWEEP, "C7", 20),            # F7->F8
    ("CCW",  0.166666666667, _PI_2, "C8", 20),                # F8->F9
    ("L",   15.166666666667),                                  # F9->F10
    ("CCW",  _R_C10, _C10_SWEEP, "C10", 20),                  # F10->F11
    ("CW",   2.333333333333, _5PI_12, "C11a", 30),            # F11->F11a
    ("L",    F11AB_TARGET),                                    # F11a->F11b
    ("CW",   2.333333333333, _5PI_12, "C11", 30),             # F11b->F12
    ("L",    _L_F12F13),                                       # F12->F13
    ("CW",   2.507553207938, _C13_SWEEP, "C13", 60),          # F13->F14
] + _CHAIN_F14_TO_F20

# Derive NE corner closure constraints from CORNER_NE_R.
# Walk F5→F20 chain to get F20 position relative to F5, then compute
# the F20→F1 distance that places F1 at the correct easting for the 90° arc.
_dE_20, _dN_20, _brg_20 = _chain_offset(_CHAIN_F5_TO_F20, start_brg=0.0)
_R_a1 = CORNER_NE_R
_d_F20_F1 = (_R_a1 - _dE_20) / math.sin(_brg_20)
_F1_N_rel = _dN_20 + _d_F20_F1 * math.cos(_brg_20)
_d_F2_F5 = -(_F1_N_rel + _R_a1)

# F2 position: east face at E = -18'0", F1.N = -13'0" exactly.
# F2.N derived from F1.N + CORNER_NE_R (F1→F2 is 90° CW arc of radius R_a1).
F2_E = -18.0                     # -18'0" east face easting
F2_N = -13.0 + CORNER_NE_R      # F1.N + R_a1
F2_BRG = 0.0                     # bearing 0 (due north)

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
