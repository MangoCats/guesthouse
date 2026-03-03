"""Outline chain walk and closure solver for the ADU Editor.

Replicates floorplan/geometry.py's chain walk and closure logic as
pure math.  No imports from floorplan/ (NF-4 compliance).  After
Phase 5 the outline_chain database table is authoritative.
"""
import math
from typing import NamedTuple


class ChainEntry(NamedTuple):
    """One segment in the outline chain."""
    seg_type: str               # 'L', 'CW', 'CCW'
    distance: float | None      # for lines
    radius: float | None        # for arcs
    sweep: float | None         # radians, for arcs
    center_name: str | None     # for arcs
    n_pts: int                  # arc discretisation
    end_name: str


class SolverResult(NamedTuple):
    """Result of the closure solver."""
    valid: bool
    d_F2_F5: float              # solved distance for first line segment
    d_F18_F1: float             # solved distance for second-to-last line segment
    closure_error: float        # residual (should be ~0 for valid)
    exit_bearing: float         # bearing at exit of inner chain


class WalkResult(NamedTuple):
    """Full chain walk result including all points."""
    points: dict[str, tuple[float, float]]  # F/C/FC series
    radii: dict[str, float]


# ---------------------------------------------------------------------------
# Core chain walk — bit-identical replication of geometry.py::_chain_offset
# ---------------------------------------------------------------------------

def chain_offset(chain, start_brg=0.0):
    """Walk chain entries from (0,0). Returns (delta_E, delta_N, exit_brg).

    Bit-identical replication of floorplan/geometry.py::_chain_offset.
    """
    E, N, brg = 0.0, 0.0, start_brg
    for seg in chain:
        if seg.seg_type == "L":
            d = seg.distance
            E += d * math.sin(brg)
            N += d * math.cos(brg)
        else:
            direction = seg.seg_type  # 'CW' or 'CCW'
            R, sweep = seg.radius, seg.sweep
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


# ---------------------------------------------------------------------------
# Closure solver
# ---------------------------------------------------------------------------

def solve_closure(chain, R_a1):
    """Solve for d_F2_F5 (first line) and d_F18_F1 (second-to-last line).

    Uses the same 2-variable linear system as floorplan/geometry.py:
    - Chain starts at F2 heading due north (bearing 0)
    - First segment is F2→F5 line due north (solved)
    - Second-to-last segment is F18→F1 line at exit bearing (solved)
    - Last segment is F1→F2: 90° CW arc of radius R_a1 (fixed)
    - Two unknowns, two equations (easting closure, northing closure)

    The "inner chain" is everything between the first and second-to-last
    segments (seq 1 through N-3 inclusive).
    """
    n = len(chain)
    # Inner chain: skip first (F2→F5) and last two (F18→F1, F1→F2)
    inner_chain = chain[1:n - 2]

    dE_18, dN_18, brg_18 = chain_offset(inner_chain, start_brg=0.0)

    # Solve: d_F18_F1 such that easting closes
    sin_brg = math.sin(brg_18)
    if abs(sin_brg) < 1e-15:
        # Degenerate: F18→F1 line is purely N-S, can't solve easting
        return SolverResult(
            valid=False, d_F2_F5=0.0, d_F18_F1=0.0,
            closure_error=abs(R_a1 - dE_18), exit_bearing=brg_18)

    d_F18_F1 = (R_a1 - dE_18) / sin_brg
    F1_N_rel = dN_18 + d_F18_F1 * math.cos(brg_18)
    d_F2_F5 = -(F1_N_rel + R_a1)

    # Validate: both distances should be positive
    if d_F2_F5 <= 0 or d_F18_F1 <= 0:
        error = max(0, -d_F2_F5) + max(0, -d_F18_F1)
        return SolverResult(
            valid=False, d_F2_F5=d_F2_F5, d_F18_F1=d_F18_F1,
            closure_error=error, exit_bearing=brg_18)

    return SolverResult(
        valid=True,
        d_F2_F5=d_F2_F5,
        d_F18_F1=d_F18_F1,
        closure_error=0.0,
        exit_bearing=brg_18,
    )


# ---------------------------------------------------------------------------
# DB row conversion
# ---------------------------------------------------------------------------

def db_rows_to_chain(rows):
    """Convert outline_chain DB rows (list of dicts) to ChainEntry list."""
    chain = []
    for row in rows:
        chain.append(ChainEntry(
            seg_type=row["seg_type"],
            distance=row.get("distance"),
            radius=row.get("radius"),
            sweep=row.get("sweep"),
            center_name=row.get("center_name"),
            n_pts=row.get("n_pts") or 60,
            end_name=row["end_name"],
        ))
    return chain


# ---------------------------------------------------------------------------
# Full chain walk (produces all F/C/FC points)
# ---------------------------------------------------------------------------

def walk_chain(chain, F2_E, F2_N, F2_BRG=0.0):
    """Walk the full outline chain from F2 and produce all points.

    Replicates floorplan/geometry.py::walk_outline_chain.
    Returns WalkResult with points and radii dicts.
    """
    E, N = F2_E, F2_N
    brg = F2_BRG
    pts = {"FC": (0.0, 0.0)}
    radii = {}

    for seg in chain:
        if seg.seg_type == "L":
            d = seg.distance
            E += d * math.sin(brg)
            N += d * math.cos(brg)
        else:
            direction = seg.seg_type
            R, sweep = seg.radius, seg.sweep
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
            pts[seg.center_name] = (cx, cy)
            # Build radii dict
            ra_name = "R_a" + seg.center_name[1:]
            if ra_name == "R_a11a":
                ra_name = "R_a11"  # C11a and C11 share radius
            radii[ra_name] = R
        pts[seg.end_name] = (E, N)

    return WalkResult(points=pts, radii=radii)


# ---------------------------------------------------------------------------
# Validation (dry-run for API-17)
# ---------------------------------------------------------------------------

def validate_chain(chain, R_a1):
    """Validate a chain without committing.  Returns status dict."""
    result = solve_closure(chain, R_a1)
    return {
        "valid": result.valid,
        "closure_error": result.closure_error,
        "d_F2_F5": result.d_F2_F5,
        "d_F18_F1": result.d_F18_F1,
    }


# ---------------------------------------------------------------------------
# Constraint solver (OE-3) — secant method
# ---------------------------------------------------------------------------

def solve_for_constraint(chain, R_a1, F2_E, F2_N,
                         from_point, to_point,
                         target_distance, adjust_seq, adjust_param,
                         max_iter=50, tol=1e-10):
    """Adjust one chain parameter to achieve a target distance between points.

    Uses the secant method on the residual:
      f(x) = distance(from_point, to_point) - target_distance
    where x is the value of chain[adjust_seq].adjust_param.

    Returns (modified_chain, solver_result) or None if unsolvable.
    """
    def _measure(ch):
        """Solve closure, walk chain, measure distance between points."""
        sr = solve_closure(ch, R_a1)
        if not sr.valid:
            return None
        # Inject solved distances
        ch2 = list(ch)
        ch2[0] = ch2[0]._replace(distance=sr.d_F2_F5)
        ch2[-2] = ch2[-2]._replace(distance=sr.d_F18_F1)
        wr = walk_chain(ch2, F2_E, F2_N)
        p1 = wr.points.get(from_point)
        p2 = wr.points.get(to_point)
        if p1 is None or p2 is None:
            return None
        return math.sqrt((p2[0] - p1[0])**2 + (p2[1] - p1[1])**2)

    def _set_param(ch, val):
        """Return chain with adjusted parameter at adjust_seq."""
        ch2 = list(ch)
        seg = ch2[adjust_seq]
        if adjust_param == "distance":
            ch2[adjust_seq] = seg._replace(distance=val)
        elif adjust_param == "radius":
            ch2[adjust_seq] = seg._replace(radius=val)
        elif adjust_param == "sweep":
            ch2[adjust_seq] = seg._replace(sweep=val)
        else:
            return None
        return ch2

    def _get_param(ch):
        seg = ch[adjust_seq]
        if adjust_param == "distance":
            return seg.distance
        elif adjust_param == "radius":
            return seg.radius
        elif adjust_param == "sweep":
            return seg.sweep
        return None

    x0 = _get_param(chain)
    if x0 is None:
        return None

    d0 = _measure(chain)
    if d0 is None:
        return None

    f0 = d0 - target_distance
    if abs(f0) < tol:
        return (chain, solve_closure(chain, R_a1))

    # Second guess: small perturbation
    delta = max(abs(x0) * 0.01, 0.001)
    x1 = x0 + delta
    ch1 = _set_param(chain, x1)
    if ch1 is None:
        return None

    d1 = _measure(ch1)
    if d1 is None:
        return None

    f1 = d1 - target_distance

    for _ in range(max_iter):
        if abs(f1) < tol:
            final_chain = _set_param(chain, x1)
            sr = solve_closure(final_chain, R_a1)
            if sr.valid:
                return (final_chain, sr)
            return None

        if abs(f1 - f0) < 1e-30:
            return None  # secant denominator too small

        x2 = x1 - f1 * (x1 - x0) / (f1 - f0)
        x0, f0 = x1, f1
        x1 = x2

        ch_new = _set_param(chain, x1)
        if ch_new is None:
            return None
        d_new = _measure(ch_new)
        if d_new is None:
            return None
        f1 = d_new - target_distance

    return None  # did not converge
