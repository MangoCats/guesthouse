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


class FlexSpec(NamedTuple):
    """One flex (solver-controlled) segment parameter."""
    seq: int
    param: str   # "distance", "radius", or "sweep"


class SolverResult(NamedTuple):
    """Result of the closure solver."""
    valid: bool
    solved_values: dict         # {seq: (param, value)}
    closure_error: float        # residual (should be ~0 for valid)
    exit_bearing: float         # bearing at exit of chain


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

def _inject_value(chain, seq, param, value):
    """Return a new chain list with chain[seq].param replaced by value."""
    ch = list(chain)
    ch[seq] = ch[seq]._replace(**{param: value})
    return ch


def _walk_residual(chain, start_brg=0.0):
    """Walk the full chain and return (dE, dN, exit_brg) closure residual."""
    return chain_offset(chain, start_brg)


def solve_closure_general(chain, flex_specs):
    """Solve for closure given user-designated flex segments.

    flex_specs: list of 3 FlexSpec — exactly 1 with param='sweep',
                exactly 2 with param in ('distance', 'radius').

    Algorithm:
    1. Angular closure: flex-sweep = 2π − Σ(other signed sweeps)
    2. Positional closure: solve for the 2 positional flex values so
       that the chain's total (dE, dN) = (0, 0).
       - Two line distances → 2×2 linear system (closed-form)
       - Any arc radius involved → Newton iteration with finite differences
    """
    TWO_PI = 2.0 * math.pi

    # Separate sweep flex from positional flex
    sweep_flex = [f for f in flex_specs if f.param == "sweep"]
    pos_flex = [f for f in flex_specs if f.param != "sweep"]

    if len(sweep_flex) != 1 or len(pos_flex) != 2:
        return SolverResult(valid=False, solved_values={},
                            closure_error=float("inf"), exit_bearing=0.0)

    sf = sweep_flex[0]
    pf_a, pf_b = pos_flex[0], pos_flex[1]

    # --- Step 1: Angular closure ---
    # Total bearing change for CW outline must be 2π.
    # Bearing change = sum of CW sweeps - sum of CCW sweeps.
    # The flex-sweep arc absorbs the remainder.
    total_sweep = 0.0
    for i, seg in enumerate(chain):
        if seg.seg_type != "L" and i != sf.seq:
            if seg.seg_type == "CW":
                total_sweep += seg.sweep
            else:  # CCW
                total_sweep -= seg.sweep

    flex_seg = chain[sf.seq]
    if flex_seg.seg_type == "CW":
        solved_sweep = (TWO_PI - total_sweep) % TWO_PI
    else:  # CCW
        solved_sweep = (total_sweep - TWO_PI) % TWO_PI
    if solved_sweep < 1e-6:
        solved_sweep = TWO_PI

    # Inject solved sweep
    ch = _inject_value(chain, sf.seq, "sweep", solved_sweep)

    # --- Step 2: Positional closure ---
    # Both positional flex are line distances → linear solve
    both_linear = (pf_a.param == "distance" and pf_b.param == "distance")

    if both_linear:
        result = _solve_positional_linear(ch, pf_a, pf_b)
    else:
        result = _solve_positional_newton(ch, pf_a, pf_b)

    if result is None:
        return SolverResult(valid=False,
                            solved_values={sf.seq: ("sweep", solved_sweep)},
                            closure_error=float("inf"), exit_bearing=0.0)

    val_a, val_b = result
    solved = {
        sf.seq: ("sweep", solved_sweep),
        pf_a.seq: (pf_a.param, val_a),
        pf_b.seq: (pf_b.param, val_b),
    }

    # Validate: distances and radii must be positive
    for seq, (param, val) in solved.items():
        if val <= 0:
            return SolverResult(valid=False, solved_values=solved,
                                closure_error=abs(val),
                                exit_bearing=0.0)

    # Verify closure
    ch_final = ch
    ch_final = _inject_value(ch_final, pf_a.seq, pf_a.param, val_a)
    ch_final = _inject_value(ch_final, pf_b.seq, pf_b.param, val_b)
    dE, dN, brg = _walk_residual(ch_final)
    err = math.hypot(dE, dN)
    if err > 1e-6:
        return SolverResult(valid=False, solved_values=solved,
                            closure_error=err, exit_bearing=brg)

    return SolverResult(valid=True, solved_values=solved,
                        closure_error=err, exit_bearing=brg)


def _solve_positional_linear(chain, pf_a, pf_b):
    """Solve for two flex line distances using a 2×2 linear system.

    Walk the chain in pieces to determine the bearing at each flex line,
    then solve: [sin(brg_a) sin(brg_b)] [d_a]   [-E_rest]
                [cos(brg_a) cos(brg_b)] [d_b] = [-N_rest]
    """
    # Walk chain with both flex distances set to 0
    ch = _inject_value(chain, pf_a.seq, "distance", 0.0)
    ch = _inject_value(ch, pf_b.seq, "distance", 0.0)
    E_rest, N_rest, _ = chain_offset(ch)

    # Find bearing at each flex line by walking up to that point
    brg_a = _bearing_at_seq(chain, pf_a.seq)
    brg_b = _bearing_at_seq(chain, pf_b.seq)

    # 2×2 linear system
    sa, ca = math.sin(brg_a), math.cos(brg_a)
    sb, cb = math.sin(brg_b), math.cos(brg_b)
    det = sa * cb - sb * ca
    if abs(det) < 1e-15:
        return None  # parallel flex lines — singular

    d_a = (-E_rest * cb + sb * N_rest) / det
    d_b = (E_rest * ca - sa * N_rest) / det
    return (d_a, d_b)


def _solve_positional_newton(chain, pf_a, pf_b, max_iter=30, tol=1e-10):
    """Solve for two positional flex values using Newton iteration.

    Used when at least one flex is an arc radius (nonlinear in position).
    """
    # Initial guesses from current chain values
    val_a = _get_flex_value(chain, pf_a)
    val_b = _get_flex_value(chain, pf_b)
    if val_a is None or val_b is None:
        return None

    eps = 1e-8  # finite difference step

    for _ in range(max_iter):
        ch = _inject_value(chain, pf_a.seq, pf_a.param, val_a)
        ch = _inject_value(ch, pf_b.seq, pf_b.param, val_b)
        dE, dN, _ = chain_offset(ch)

        if math.hypot(dE, dN) < tol:
            return (val_a, val_b)

        # Jacobian by finite differences
        ch_a = _inject_value(ch, pf_a.seq, pf_a.param, val_a + eps)
        dEa, dNa, _ = chain_offset(ch_a)
        dEdva = (dEa - dE) / eps
        dNdva = (dNa - dN) / eps

        ch_b = _inject_value(ch, pf_b.seq, pf_b.param, val_b + eps)
        dEb, dNb, _ = chain_offset(ch_b)
        dEdvb = (dEb - dE) / eps
        dNdvb = (dNb - dN) / eps

        det = dEdva * dNdvb - dEdvb * dNdva
        if abs(det) < 1e-30:
            return None  # singular Jacobian

        # Newton step: [va, vb] -= J^{-1} @ [dE, dN]
        delta_a = (dNdvb * dE - dEdvb * dN) / det
        delta_b = (dEdva * dN - dNdva * dE) / det
        val_a -= delta_a
        val_b -= delta_b

    return None  # did not converge


def _bearing_at_seq(chain, seq):
    """Compute the bearing at the start of chain[seq]."""
    brg = 0.0
    for i in range(seq):
        seg = chain[i]
        if seg.seg_type == "CW":
            brg += seg.sweep
        elif seg.seg_type == "CCW":
            brg -= seg.sweep
        # Lines don't change bearing
    return brg


def _get_flex_value(chain, fs):
    """Get current value of a flex parameter from the chain."""
    seg = chain[fs.seq]
    return getattr(seg, fs.param)


def solve_closure(chain, R_a1):
    """Legacy wrapper: solve with default flex (seq 0, n-2, n-1)."""
    n = len(chain)
    flex_specs = [
        FlexSpec(0, "distance"),
        FlexSpec(n - 2, "distance"),
        FlexSpec(n - 1, "sweep"),
    ]
    return solve_closure_general(chain, flex_specs)


def flex_specs_from_chain_rows(chain_rows):
    """Extract FlexSpec list from DB outline_chain rows.

    Falls back to legacy defaults (seq 0, n-2, n-1) if no flex column.
    """
    specs = []
    for row in chain_rows:
        flex = row.get("flex")
        if flex:
            specs.append(FlexSpec(row["seq"], flex))
    if len(specs) == 3:
        return specs
    # Fallback to defaults
    n = len(chain_rows)
    return [
        FlexSpec(0, "distance"),
        FlexSpec(n - 2, "distance"),
        FlexSpec(n - 1, "sweep"),
    ]


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

def validate_chain(chain, flex_specs=None):
    """Validate a chain without committing.  Returns status dict."""
    if flex_specs is None:
        n = len(chain)
        flex_specs = [FlexSpec(0, "distance"), FlexSpec(n - 2, "distance"),
                      FlexSpec(n - 1, "sweep")]
    result = solve_closure_general(chain, flex_specs)
    return {
        "valid": result.valid,
        "closure_error": result.closure_error,
        "solved_values": {str(k): {"param": v[0], "value": v[1]}
                          for k, v in result.solved_values.items()},
    }


# ---------------------------------------------------------------------------
# Constraint solver (OE-3) — secant method
# ---------------------------------------------------------------------------

def solve_for_constraint(chain, R_a1, F2_E, F2_N,
                         from_point, to_point,
                         target_distance, adjust_seq, adjust_param,
                         flex_specs=None,
                         max_iter=50, tol=1e-10):
    """Adjust one chain parameter to achieve a target distance between points.

    Uses the secant method on the residual:
      f(x) = distance(from_point, to_point) - target_distance
    where x is the value of chain[adjust_seq].adjust_param.

    Returns (modified_chain, solver_result) or None if unsolvable.
    """
    if flex_specs is None:
        n = len(chain)
        flex_specs = [FlexSpec(0, "distance"), FlexSpec(n - 2, "distance"),
                      FlexSpec(n - 1, "sweep")]

    def _measure(ch):
        """Solve closure, walk chain, measure distance between points."""
        sr = solve_closure_general(ch, flex_specs)
        if not sr.valid:
            return None
        # Inject solved values
        ch2 = list(ch)
        for seq, (param, value) in sr.solved_values.items():
            ch2[seq] = ch2[seq]._replace(**{param: value})
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
        return (chain, solve_closure_general(chain, flex_specs))

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
            sr = solve_closure_general(final_chain, flex_specs)
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
