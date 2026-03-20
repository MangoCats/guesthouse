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
    bearing_flex: int = 0       # 1 = line bearing may rotate (opt-in); 0 = fixed


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


def _solve_positional_linear(chain, pf_a, pf_b, target=(0.0, 0.0),
                              start_brg=0.0):
    """Solve for two flex line distances using a 2×2 linear system.

    Walk the chain in pieces to determine the bearing at each flex line,
    then solve: [sin(brg_a) sin(brg_b)] [d_a]   [tE - E_rest]
                [cos(brg_a) cos(brg_b)] [d_b] = [tN - N_rest]

    target: (tE, tN) endpoint to hit — (0,0) for closure, or a specific point
            for sub-chain solving with a pivot.
    start_brg: bearing at the start of the chain (0.0 for full-chain).
    """
    tE, tN = target
    # Walk chain with both flex distances set to 0
    ch = _inject_value(chain, pf_a.seq, "distance", 0.0)
    ch = _inject_value(ch, pf_b.seq, "distance", 0.0)
    E_rest, N_rest, _ = chain_offset(ch, start_brg)

    # Find bearing at each flex line by walking up to that point
    brg_a = _bearing_at_seq(chain, pf_a.seq, start_brg)
    brg_b = _bearing_at_seq(chain, pf_b.seq, start_brg)

    # 2×2 linear system
    sa, ca = math.sin(brg_a), math.cos(brg_a)
    sb, cb = math.sin(brg_b), math.cos(brg_b)
    det = sa * cb - sb * ca
    if abs(det) < 1e-15:
        return None  # parallel flex lines — singular

    rhs_E = tE - E_rest
    rhs_N = tN - N_rest
    d_a = (rhs_E * cb - sb * rhs_N) / det
    d_b = (-rhs_E * ca + sa * rhs_N) / det
    return (d_a, d_b)


def _solve_positional_newton(chain, pf_a, pf_b, target=(0.0, 0.0),
                             max_iter=30, tol=1e-10, start_brg=0.0):
    """Solve for two positional flex values using Newton iteration.

    Used when at least one flex is an arc radius (nonlinear in position).
    target: (tE, tN) endpoint to hit — (0,0) for closure.
    start_brg: bearing at the start of the chain (0.0 for full-chain).
    """
    tE, tN = target
    # Initial guesses from current chain values
    val_a = _get_flex_value(chain, pf_a)
    val_b = _get_flex_value(chain, pf_b)
    if val_a is None or val_b is None:
        return None

    eps = 1e-8  # finite difference step

    for _ in range(max_iter):
        ch = _inject_value(chain, pf_a.seq, pf_a.param, val_a)
        ch = _inject_value(ch, pf_b.seq, pf_b.param, val_b)
        dE, dN, _ = chain_offset(ch, start_brg)
        rE, rN = dE - tE, dN - tN

        if math.hypot(rE, rN) < tol:
            return (val_a, val_b)

        # Jacobian by finite differences
        ch_a = _inject_value(ch, pf_a.seq, pf_a.param, val_a + eps)
        dEa, dNa, _ = chain_offset(ch_a, start_brg)
        dEdva = (dEa - dE) / eps
        dNdva = (dNa - dN) / eps

        ch_b = _inject_value(ch, pf_b.seq, pf_b.param, val_b + eps)
        dEb, dNb, _ = chain_offset(ch_b, start_brg)
        dEdvb = (dEb - dE) / eps
        dNdvb = (dNb - dN) / eps

        det = dEdva * dNdvb - dEdvb * dNdva
        if abs(det) < 1e-30:
            return None  # singular Jacobian

        # Newton step: [va, vb] -= J^{-1} @ [rE, rN]
        delta_a = (dNdvb * rE - dEdvb * rN) / det
        delta_b = (dEdva * rN - dNdva * rE) / det
        val_a -= delta_a
        val_b -= delta_b

    return None  # did not converge


def _bearing_at_seq(chain, seq, start_brg=0.0):
    """Compute the bearing at the start of chain[seq]."""
    brg = start_brg
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

    Returns the 3-spec list for whole-chain solvers (solve_closure_general).
    Falls back to legacy defaults (seq 0, n-2, n-1) if no flex column or
    if the DB has more than 3 (pivot mode — caller should use
    all_flex_specs_from_chain_rows instead).
    """
    specs = []
    for row in chain_rows:
        flex = row.get("flex")
        if flex:
            specs.append(FlexSpec(row["seq"], flex))
    if len(specs) == 3:
        return specs
    # Fallback to defaults (covers 0 specs and pivot-mode 6-spec DB)
    n = len(chain_rows)
    return [
        FlexSpec(0, "distance"),
        FlexSpec(n - 2, "distance"),
        FlexSpec(n - 1, "sweep"),
    ]


def all_flex_specs_from_chain_rows(chain_rows):
    """Extract all FlexSpec entries from DB rows without count gating.

    Used by pivot-aware code that needs the full 6-spec set or whatever
    was explicitly written to the DB.  Falls back to defaults only when
    no specs exist at all.
    """
    specs = []
    for row in chain_rows:
        flex = row.get("flex")
        if flex:
            specs.append(FlexSpec(row["seq"], flex))
    if specs:
        return specs
    n = len(chain_rows)
    return [
        FlexSpec(0, "distance"),
        FlexSpec(n - 2, "distance"),
        FlexSpec(n - 1, "sweep"),
    ]


# ---------------------------------------------------------------------------
# Pivot-aware solver
# ---------------------------------------------------------------------------

def point_name_to_seq(chain, point_name):
    """Find the seq index of the segment whose end_name matches point_name.

    Returns the seq index, or None if not found.
    The anchor/pivot seq is the index *after* this point: i.e., seq+1
    (wrapping) is where that section's sub-chain starts.
    """
    for i, seg in enumerate(chain):
        if seg.end_name == point_name:
            return i
    return None


def identify_section(seq, anchor_start, pivot_start, n):
    """Determine which section a seq index belongs to.

    anchor_start: first seq of section A (seg after anchor point)
    pivot_start:  first seq of section B (seg after pivot point)
    n: total chain length

    Returns 'A' or 'B'.
    """
    # Build section A seq set: anchor_start forward (wrapping) to pivot_start-1
    a_seqs = set()
    i = anchor_start
    while i != pivot_start:
        a_seqs.add(i)
        i = (i + 1) % n
    return "A" if seq in a_seqs else "B"


def section_seqs(anchor_start, pivot_start, n):
    """Return (section_a_seqs, section_b_seqs) as lists in chain-walk order."""
    a = []
    i = anchor_start
    while i != pivot_start:
        a.append(i)
        i = (i + 1) % n
    b = []
    i = pivot_start
    while i != anchor_start:
        b.append(i)
        i = (i + 1) % n
    return a, b


def validate_pivot_placement(chain, anchor_start, pivot_start):
    """Check that both sections have enough segments for 3 flex vars each.

    Each section needs >= 3 segments, >= 1 arc, >= 2 adjustable params.
    Returns (valid, error_message).
    """
    n = len(chain)
    a_seqs, b_seqs = section_seqs(anchor_start, pivot_start, n)

    for label, seqs in [("A", a_seqs), ("B", b_seqs)]:
        if len(seqs) < 3:
            return False, f"Section {label} has only {len(seqs)} segments (need >= 3)"
        arcs = [s for s in seqs if chain[s].seg_type != "L"]
        if len(arcs) < 1:
            return False, f"Section {label} has no arcs (need >= 1 for sweep flex)"
        # Need 2 segments with adjustable positional params (distance or radius)
        pos_count = 0
        for s in seqs:
            if chain[s].seg_type == "L":
                pos_count += 1  # can flex distance
            else:
                pos_count += 1  # can flex radius
        if pos_count < 2:
            return False, f"Section {label} needs >= 2 adjustable segments"
    return True, ""


def auto_assign_section_flex(chain, seqs, edited_seq=None):
    """Auto-pick 3 flex vars for a section: 1 sweep + 2 positional.

    Sweep placement logic (bearing_flex aware):

    When edited_seq is an arc in this section, we know the user's edit
    will rotate everything downstream of that arc by some Δ.  The sweep
    flex arc (at position m) compensates with −Δ, so lines *after* m
    have net-zero bearing change.  Lines *between* edited_seq and m
    absorb the full Δ — so those must be bearing_flex=1 (opt-in).

    Algorithm: starting from the position after edited_seq, scan forward
    until a run of bearing_flex=1 lines ends at an arc.  Place the sweep
    flex at that arc, so only the bearing_flex=1 lines between the edit
    and the sweep arc change bearing.

    When edited_seq is not an arc, or not in this section, fall back to:
    placing the sweep flex after the last fixed-bearing (bearing_flex=0)
    line in the section (section-wide safe placement).

    Last resort: last arc in the section (minimises downstream exposure).

    Strategy: positional → first and last lines, fall back to arcs.
    Returns [FlexSpec, FlexSpec, FlexSpec].
    """
    arcs = [s for s in seqs if chain[s].seg_type != "L"]
    seq_pos = {s: i for i, s in enumerate(seqs)}

    # --- Sweep flex ---
    sweep_seq = None

    # Dynamic placement: only when the edited segment is an arc in this section.
    if (edited_seq is not None
            and edited_seq in seq_pos
            and chain[edited_seq].seg_type != "L"):
        k = seq_pos[edited_seq]
        # Scan forward from k+1; find first arc m such that every line
        # between k+1 and m-1 (exclusive) has bearing_flex=1.
        for i in range(k + 1, len(seqs)):
            s = seqs[i]
            if chain[s].seg_type != "L":
                # All lines between k+1 and i are bearing_flex=1?
                lines_between = [
                    seqs[j] for j in range(k + 1, i)
                    if chain[seqs[j]].seg_type == "L"
                ]
                if all(chain[l].bearing_flex for l in lines_between):
                    sweep_seq = s
                    break
        # If nothing found forward, try arcs before the edit (wrap-around
        # rarely needed, but handle gracefully by falling through to below).

    # Section-wide fallback: sweep arc after all fixed-bearing lines.
    if sweep_seq is None:
        fixed_bearing_positions = [
            seq_pos[s] for s in seqs
            if chain[s].seg_type == "L" and not chain[s].bearing_flex
        ]
        last_fixed_pos = max(fixed_bearing_positions) if fixed_bearing_positions else -1
        valid_sweep = [s for s in arcs if seq_pos[s] > last_fixed_pos]
        sweep_seq = valid_sweep[0] if valid_sweep else arcs[-1]

    # --- Positional flex: prefer first and last lines, fall back to arcs ---
    pos_candidates = [s for s in seqs if s != sweep_seq]
    line_candidates = [s for s in pos_candidates if chain[s].seg_type == "L"]
    arc_candidates = [s for s in pos_candidates if chain[s].seg_type != "L"]

    pos_specs = []
    if len(line_candidates) >= 2:
        pos_specs.append(FlexSpec(line_candidates[0], "distance"))
        pos_specs.append(FlexSpec(line_candidates[-1], "distance"))
    elif len(line_candidates) == 1:
        pos_specs.append(FlexSpec(line_candidates[0], "distance"))
        pos_specs.append(FlexSpec(arc_candidates[0], "radius"))
    else:
        pos_specs.append(FlexSpec(arc_candidates[0], "radius"))
        pos_specs.append(FlexSpec(arc_candidates[-1], "radius"))

    return [FlexSpec(sweep_seq, "sweep")] + pos_specs


def solve_subchain(subchain, flex_specs, target_dE, target_dN,
                   target_brg_change, start_brg=0.0):
    """Solve a sub-chain to hit a specific endpoint and bearing change.

    Like solve_closure_general but targets (target_dE, target_dN) instead
    of (0, 0), and target_brg_change instead of 2π.

    subchain: list of ChainEntry (re-indexed starting from 0)
    flex_specs: 3 FlexSpec with seq indices relative to subchain (0-based)
    target_dE, target_dN: required displacement from sub-chain start
    target_brg_change: required total bearing change across the sub-chain
    start_brg: bearing at the start of the sub-chain

    Returns SolverResult with seq indices relative to subchain.
    """
    # Separate sweep flex from positional flex
    sweep_flex = [f for f in flex_specs if f.param == "sweep"]
    pos_flex = [f for f in flex_specs if f.param != "sweep"]

    if len(sweep_flex) != 1 or len(pos_flex) != 2:
        return SolverResult(valid=False, solved_values={},
                            closure_error=float("inf"), exit_bearing=0.0)

    sf = sweep_flex[0]
    pf_a, pf_b = pos_flex[0], pos_flex[1]

    # --- Step 1: Angular closure ---
    # Total bearing change must equal target_brg_change
    total_sweep = 0.0
    for i, seg in enumerate(subchain):
        if seg.seg_type != "L" and i != sf.seq:
            if seg.seg_type == "CW":
                total_sweep += seg.sweep
            else:
                total_sweep -= seg.sweep

    flex_seg = subchain[sf.seq]
    if flex_seg.seg_type == "CW":
        solved_sweep = target_brg_change - total_sweep
    else:  # CCW
        solved_sweep = -(target_brg_change - total_sweep)

    # Normalise to positive
    TWO_PI = 2.0 * math.pi
    if solved_sweep < 0:
        solved_sweep = solved_sweep % TWO_PI
    if solved_sweep < 1e-6:
        solved_sweep = TWO_PI

    # Inject solved sweep
    ch = _inject_value(subchain, sf.seq, "sweep", solved_sweep)

    # --- Step 2: Positional — solve to hit target endpoint ---
    both_linear = (pf_a.param == "distance" and pf_b.param == "distance")
    target = (target_dE, target_dN)

    if both_linear:
        result = _solve_positional_linear(ch, pf_a, pf_b, target=target,
                                          start_brg=start_brg)
    else:
        result = _solve_positional_newton(ch, pf_a, pf_b, target=target,
                                          start_brg=start_brg)

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
                                closure_error=abs(val), exit_bearing=0.0)

    # Verify endpoint
    ch_final = ch
    ch_final = _inject_value(ch_final, pf_a.seq, pf_a.param, val_a)
    ch_final = _inject_value(ch_final, pf_b.seq, pf_b.param, val_b)
    dE, dN, brg = chain_offset(ch_final, start_brg)
    err = math.hypot(dE - target_dE, dN - target_dN)
    if err > 1e-6:
        return SolverResult(valid=False, solved_values=solved,
                            closure_error=err, exit_bearing=brg)

    return SolverResult(valid=True, solved_values=solved,
                        closure_error=err, exit_bearing=brg)


def solve_with_pivot(chain, anchor_start, pivot_start,
                     section_a_flex, section_b_flex,
                     edited_seq, start_E, start_N, start_brg=0.0,
                     ref_chain=None):
    """Solve the section containing the edited segment, keeping the other fixed.

    anchor_start: seq index where section A begins (seg after anchor point)
    pivot_start: seq index where section B begins (seg after pivot point)
    section_a_flex: 3 FlexSpec with original (full-chain) seq indices
    section_b_flex: 3 FlexSpec with original (full-chain) seq indices
    edited_seq: which segment was edited (original seq index)
    ref_chain: pre-edit chain used to compute target positions (if None, uses chain)

    Returns SolverResult with solved_values keyed by original seq indices.
    """
    n = len(chain)
    a_seqs, b_seqs = section_seqs(anchor_start, pivot_start, n)

    # Walk the REFERENCE chain rotated to start at anchor_start.
    #
    # The rotation ensures that the global walk and the sub-chain B walk
    # share identical entry bearings at every section-boundary segment.
    # Without rotation, the bearing at seg 0 in the global walk (= start_brg)
    # differs from the bearing at seg 0 in the sub-chain B walk (accumulated
    # through the pivot arc and subsequent arcs), causing solved flex distances
    # to be geometrically inconsistent with the renderer — manifesting as all
    # chain points drifting after a Section B edit.
    #
    # start_E, start_N: absolute position of the anchor point.
    # start_brg: bearing entering chain[anchor_start] (first Section A seg).
    walk_src = ref_chain if ref_chain is not None else chain
    rotated_src = [walk_src[(anchor_start + i) % n] for i in range(n)]
    wr = walk_chain(rotated_src, start_E, start_N, start_brg)

    # Find anchor and pivot point names
    anchor_point_seq = (anchor_start - 1) % n
    pivot_point_seq = (pivot_start - 1) % n
    anchor_name = chain[anchor_point_seq].end_name
    pivot_name = chain[pivot_point_seq].end_name

    anchor_pos = wr.points.get(anchor_name)
    pivot_pos = wr.points.get(pivot_name)
    if anchor_pos is None or pivot_pos is None:
        return SolverResult(valid=False, solved_values={},
                            closure_error=float("inf"), exit_bearing=0.0)

    # Compute bearings in rotated order.
    # bearings[i] = bearing entering rotated_src[i] = chain[(anchor_start+i)%n]
    # bearings[0] = start_brg (= brg_at_anchor_start, since rotation starts here)
    # bearings[n] = exit bearing (≈ start_brg + 2π for a closed chain)
    bearings = [start_brg]
    brg = start_brg
    for seg in rotated_src:
        if seg.seg_type == "CW":
            brg += seg.sweep
        elif seg.seg_type == "CCW":
            brg -= seg.sweep
        bearings.append(brg)

    brg_at_anchor_start = start_brg  # bearings[0]
    pivot_rotated_idx = (pivot_start - anchor_start) % n
    brg_at_pivot_start = bearings[pivot_rotated_idx]

    # Determine which section the edit is in
    side = identify_section(edited_seq, anchor_start, pivot_start, n)

    if side == "A":
        sub_seqs = a_seqs
        flex_specs = section_a_flex
        # Sub-chain goes from anchor to pivot
        sub_start_pos = anchor_pos
        sub_start_brg = brg_at_anchor_start
        target_pos = pivot_pos
        target_brg = brg_at_pivot_start
    else:
        sub_seqs = b_seqs
        flex_specs = section_b_flex
        # Sub-chain goes from pivot to anchor
        sub_start_pos = pivot_pos
        sub_start_brg = brg_at_pivot_start
        target_pos = anchor_pos
        target_brg = brg_at_anchor_start + 2.0 * math.pi  # wrap

    # Target displacement relative to sub-chain start
    target_dE = target_pos[0] - sub_start_pos[0]
    target_dN = target_pos[1] - sub_start_pos[1]
    target_brg_change = target_brg - sub_start_brg

    # Extract sub-chain and remap flex specs to 0-based indices
    subchain = [chain[s] for s in sub_seqs]
    seq_map = {orig: local for local, orig in enumerate(sub_seqs)}
    local_flex = [FlexSpec(seq_map[f.seq], f.param) for f in flex_specs]

    result = solve_subchain(subchain, local_flex,
                            target_dE, target_dN, target_brg_change,
                            start_brg=sub_start_brg)

    # Remap solved values back to original seq indices
    remapped = {}
    for local_seq, (param, value) in result.solved_values.items():
        orig_seq = sub_seqs[local_seq]
        remapped[orig_seq] = (param, value)

    return SolverResult(valid=result.valid, solved_values=remapped,
                        closure_error=result.closure_error,
                        exit_bearing=result.exit_bearing)


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
            bearing_flex=row.get("bearing_flex") or 0,
        ))
    return chain


# ---------------------------------------------------------------------------
# Full chain walk (produces all F/C/FC points)
# ---------------------------------------------------------------------------

def walk_chain(chain, start_E, start_N, start_brg=0.0):
    """Walk the full outline chain from the anchor and produce all points.

    start_E, start_N: absolute position of the anchor (chain walk origin).
    start_brg: bearing at the start of the first segment.
    Returns WalkResult with points and radii dicts.
    """
    E, N = start_E, start_N
    brg = start_brg
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

def solve_for_constraint(chain, R_a1, start_E, start_N,
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
        wr = walk_chain(ch2, start_E, start_N)
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
