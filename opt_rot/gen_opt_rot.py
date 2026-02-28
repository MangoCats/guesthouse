"""Generate path_area SVGs for C17 sweep angles with alignment constraints.

For each C17 sweep angle (30-42 deg), finds the rigid-body placement
(rotation + translation) that best satisfies the 4 alignment constraints:
  1. T3 on F18-F1 line
  2. F12-F13 tangent to TC1 arc (signed dist = +R1)
  3. F16 on Pi5-PiX line
  4. F17 on Pi5-PiX line

With 4 constraints and 4 unknowns (rot, dx, dy, C17), the system is
fully constrained — a unique C17 exists with zero residual.
"""
import os, sys, math, datetime

# Ensure project root is on sys.path for package imports
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

from shared.types import Point, LineSeg, ArcSeg, Segment
from shared.geometry import poly_area, path_polygon
from shared.svg import W, H, git_describe
from floorplan.geometry import (
    OUTLINE_CHAIN, CHAIN_POINT_NAMES,
    F2_E, F2_N, F2_BRG,
)
from floorplan.constants import CORNER_SW_R, F11AB_TARGET
from survey.gen_path_svg import (
    compute_all, render_layer, build_outline_cfg,
    outer_cfg, inset_cfg,
)

try:
    from scipy.optimize import least_squares, fsolve
except ImportError:
    sys.exit("scipy required: pip install -e '.[adjust]'")

# ============================================================
# Constants
# ============================================================
_PI_2 = math.pi / 2
_C10_SWEEP = math.pi / 2 - math.atan(1.0 / 3.0)
_C11_SWEEP = math.pi / 2 - math.atan(1.0 / 3.0)
_C13_SWEEP = math.atan(1.0 / 3.0)

_R_a1 = CORNER_SW_R
_R_C15C17 = 1.808727374505       # C15/C17 arc radius (feet)
_L_F16F17 = 5.000000000000       # F16-F17 line length (feet)

# Locked chain: F5 -> F14 (first 11 entries of _CHAIN_F5_TO_F18)
_CHAIN_F5_TO_F14 = [
    ("CW",   2.333333333333, _PI_2, "C5", 20),
    ("L",    5.166666666667),
    ("CW",   2.333333333333, _PI_2, "C7", 20),
    ("CCW",  0.166666666667, _PI_2, "C8", 20),
    ("L",   15.500000000000),
    ("CCW",  1.322854905602, _C10_SWEEP, "C10", 20),
    ("CW",   2.333333333333, _C10_SWEEP, "C11a", 30),
    ("L",    F11AB_TARGET),
    ("CW",   2.333333333333, _C11_SWEEP, "C11", 30),
    ("L",   11.779557008578),
    ("CW",   1.808727374505, _C13_SWEEP, "C13", 60),
]

# d_F2_F5: extracted from current OUTLINE_CHAIN (locked)
_d_F2_F5 = OUTLINE_CHAIN[0][1]


# ============================================================
# Chain walk
# ============================================================
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
            else:
                cx = E - R * math.cos(brg)
                cy = N + R * math.sin(brg)
                alpha = math.atan2(N - cy, E - cx) + sweep
                E = cx + R * math.cos(alpha)
                N = cy + R * math.sin(alpha)
                brg -= sweep
    return E, N, brg


# Walk locked chain once: F5 -> F14 offset from (0,0) bearing 0
_dE14, _dN14, _brg14 = _chain_offset(_CHAIN_F5_TO_F14, start_brg=0.0)


# ============================================================
# Parameterized geometry
# ============================================================
def solve_for_angle(c17_deg):
    """Solve closure for a given C17 sweep angle (degrees).

    Returns (d_F14_F15, d_F18_F1, c15_sweep, c17_sweep, outline_chain).
    """
    s2 = math.radians(c17_deg)
    s1 = _PI_2 - s2   # C15 sweep

    # Derived line lengths from closure equations
    d_F14_F15 = _d_F2_F5 + _dN14 + _R_a1 - _R_C15C17 - _L_F16F17 * math.cos(s1)
    d_F18_F1 = _dE14 - _R_C15C17 - _L_F16F17 * math.sin(s1) - _R_a1

    # Build full outline chain
    chain_F14_to_F18 = [
        ("L",   d_F14_F15),
        ("CW",  _R_C15C17, s1, "C15", 20),
        ("L",   _L_F16F17),
        ("CW",  _R_C15C17, s2, "C17", 20),
    ]
    outline_chain = [
        ("L", _d_F2_F5),
    ] + _CHAIN_F5_TO_F14 + chain_F14_to_F18 + [
        ("L", d_F18_F1),
        ("CW", _R_a1, _PI_2, "C1", 20),
    ]

    return d_F14_F15, d_F18_F1, s1, s2, outline_chain


def walk_chain(outline_chain):
    """Walk outline chain from F2, return pts dict with all F/C/FC points."""
    E, N = F2_E, F2_N
    brg = F2_BRG
    fp_pts: dict[str, Point] = {"FC": (0.0, 0.0)}

    for seg, name in zip(outline_chain, CHAIN_POINT_NAMES):
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
            else:
                cx = E - R * math.cos(brg)
                cy = N + R * math.sin(brg)
                alpha = math.atan2(N - cy, E - cx) + sweep
                E, N = cx + R * math.cos(alpha), cy + R * math.sin(alpha)
                brg -= sweep
            fp_pts[center_name] = (cx, cy)
        fp_pts[name] = (E, N)

    return fp_pts


def build_outline_segs(outline_chain):
    """Build outline segment list from chain definition (rotated: F1->F2 first)."""
    start_names = ["F2"] + CHAIN_POINT_NAMES[:-1]
    end_names = CHAIN_POINT_NAMES

    segs: list[Segment] = []
    for entry, start, end in zip(outline_chain, start_names, end_names):
        if entry[0] == "L":
            segs.append(LineSeg(start, end))
        else:
            segs.append(ArcSeg(start, end, entry[3], entry[1], entry[0], entry[4]))

    # Rotate so F1->F2 comes first (matches outline convention)
    return segs[-1:] + segs[:-1]


def build_radii(outline_chain):
    """Extract radii dict from chain arc entries."""
    radii: dict[str, float] = {}
    for entry in outline_chain:
        if entry[0] != "L":
            center_name = entry[3]
            ra_name = "R_a" + center_name[1:]
            if ra_name == "R_a11a":
                ra_name = "R_a11"
            radii[ra_name] = entry[1]
    return radii


# ============================================================
# Alignment constraints
# ============================================================
def signed_dist(P, A, B):
    """Signed perpendicular distance from P to line A->B. Left-positive."""
    dx, dy = B[0] - A[0], B[1] - A[1]
    L = math.hypot(dx, dy)
    if L < 1e-15:
        return 0.0
    return (dx * (P[1] - A[1]) - dy * (P[0] - A[0])) / L


def xform_pts(std_pts, rot, dx, dy):
    """Rotate CW by rot around origin, then translate by (dx, dy)."""
    c, s = math.cos(rot), math.sin(rot)
    return {name: (e * c + n * s + dx, -e * s + n * c + dy)
            for name, (e, n) in std_pts.items()}


def alignment_residuals(placement, c17_deg, survey_info):
    """4 alignment residuals for 3 placement unknowns at fixed C17.

    survey_info: (T3, TC1, R1, Pi5, PiX)
    """
    rot, dx, dy = placement
    T3, TC1, R1, Pi5, PiX = survey_info

    d1, d2, _, _, chain = solve_for_angle(c17_deg)
    if d1 <= 0 or d2 <= 0:
        return [1e6, 1e6, 1e6, 1e6]

    xf = xform_pts(walk_chain(chain), rot, dx, dy)

    return [
        signed_dist(T3, xf["F18"], xf["F1"]),             # T3 on F18-F1
        signed_dist(TC1, xf["F12"], xf["F13"]) - R1,      # F12-F13 tangent to TC1
        signed_dist(xf["F16"], Pi5, PiX),                  # F16 on Pi5-PiX
        signed_dist(xf["F17"], Pi5, PiX),                  # F17 on Pi5-PiX
    ]


def full_residuals(x, survey_info):
    """4 residuals for full 4-variable system (rot, dx, dy, C17_deg)."""
    return alignment_residuals(x[:3], x[3], survey_info)


# ============================================================
# Main
# ============================================================
if __name__ == "__main__":
    # Compute survey/inset data once (shared across all angles)
    data = compute_all()
    base_pts = dict(data["pts"])
    to_svg = data["to_svg"]
    outer_segs = data["outer_segs"]
    inset_segs = data["inset_segs"]

    # Survey alignment targets (fixed in primary frame)
    T3 = base_pts["T3"]
    TC1 = base_pts["TC1"]
    Pi5 = base_pts["Pi5"]
    PiX = base_pts["PiX"]

    # R1 (TC1 arc radius) from outer segments
    R1 = next(s.radius for s in outer_segs
              if isinstance(s, ArcSeg) and s.center == "TC1")

    # Verify tangent sign: TC1 is on the LEFT of F12->F13 (positive signed dist)
    tang_dist = signed_dist(TC1, base_pts["F12"], base_pts["F13"])
    assert tang_dist > 0, f"Expected TC1 left of F12->F13, got {tang_dist}"
    print(f"Tangent check: signed_dist(TC1, F12-F13) = {tang_dist:.6f}, R1 = {R1:.6f}")

    survey_info = (T3, TC1, R1, Pi5, PiX)

    # Verify current design satisfies all constraints
    c17_current = math.degrees(0.629724265938)
    r0 = alignment_residuals([0, 0, 0], c17_current, survey_info)
    print(f"\nCurrent design (C17={c17_current:.4f}\u00b0) residuals:")
    for name, val in zip(
        ["T3 on F18-F1", "F12-F13 tangent", "F16 on line", "F17 on line"], r0
    ):
        print(f"  {name}: {val:.2e} ft")

    # Solve for exact C17 + placement
    sol, info, ier, msg = fsolve(
        full_residuals, [0.0, 0.0, 0.0, c17_current],
        args=(survey_info,), full_output=True,
    )
    rot_sol, dx_sol, dy_sol, c17_sol = sol
    d1_sol, d2_sol, s1_sol, _, _ = solve_for_angle(c17_sol)
    r_sol = full_residuals(sol, survey_info)
    print(f"\n{'='*60}")
    print(f"EXACT SOLUTION (fsolve ier={ier}):")
    print(f"  C17  = {c17_sol:.6f}\u00b0")
    print(f"  C15  = {math.degrees(s1_sol):.6f}\u00b0")
    print(f"  rot  = {rot_sol:.2e} rad ({math.degrees(rot_sol):.4e}\u00b0)")
    print(f"  dx   = {dx_sol:.2e} ft ({dx_sol*12:.2e}\")")
    print(f"  dy   = {dy_sol:.2e} ft ({dy_sol*12:.2e}\")")
    print(f"  d(F14-F15) = {d1_sol:.6f}'")
    print(f"  d(F18-F1)  = {d2_sol:.6f}'")
    print(f"  Residuals: {[f'{r:.2e}' for r in r_sol]}")
    print(f"{'='*60}")

    # Sweep: for each C17, find best-fit placement
    out_dir = os.path.dirname(os.path.abspath(__file__))
    count = 0
    prev_placement = [0.0, 0.0, 0.0]

    # Build list of (c17_deg, filename_label) pairs
    angles = [(x / 10.0, f"{x:03d}") for x in range(300, 421, 10)]
    angles.append((math.degrees(math.atan(7.0 / 12.0)), "712"))
    angles.append((math.degrees(math.atan(8.0 / 12.0)), "812"))

    print(f"\n{'C17':>8s} {'C15':>6s} {'rot':>9s} {'dx':>9s} {'dy':>9s} "
          f"{'RMS':>10s} {'d14-15':>8s} {'d18-1':>8s} {'area':>8s}")
    print("-" * 95)

    for c17_deg, label in angles:
        d1, d2, s1, s2, chain = solve_for_angle(c17_deg)

        if d1 <= 0 or d2 <= 0:
            print(f"{c17_deg:8.4f}\u00b0  SKIPPED (non-positive length)")
            continue

        # Find best-fit placement (3 unknowns, 4 constraints → least squares)
        result = least_squares(
            alignment_residuals, prev_placement,
            args=(c17_deg, survey_info),
        )
        rot, dx, dy = result.x
        rms = math.sqrt(result.cost / 4)
        prev_placement = list(result.x)

        # Walk chain and apply transform
        std_pts = walk_chain(chain)
        xf = xform_pts(std_pts, rot, dx, dy)

        # Merge with survey points
        pts = dict(base_pts)
        pts.update(xf)

        # Build outline geometry
        outline_segs = build_outline_segs(chain)
        radii = build_radii(chain)
        outer_poly = path_polygon(outline_segs, pts)
        outline_area = poly_area(outer_poly)
        outline_cfg = build_outline_cfg(outline_segs, pts, radii)

        print(f"{c17_deg:8.4f}\u00b0 {math.degrees(s1):5.1f}\u00b0 "
              f"{math.degrees(rot):+9.4f}\u00b0 {dx:+9.4f}' {dy:+9.4f}' "
              f"{rms:10.4e} {d1:8.4f}' {d2:8.4f}' {outline_area:8.2f}")

        # Individual constraint residuals
        resid = alignment_residuals([rot, dx, dy], c17_deg, survey_info)

        # --- Render SVG ---
        lines: list[str] = []
        lines.append(f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}"'
                     f' viewBox="0 0 {W} {H}">')
        lines.append(f'<rect width="{W}" height="{H}" fill="white"/>')
        lines.append('<defs>')
        lines.append('  <marker id="ah" markerWidth="8" markerHeight="6" refX="8" refY="3"'
                     ' orient="auto"><polygon points="0 0, 8 3, 0 6" fill="#333"/></marker>')
        lines.append(f'  <clipPath id="page"><rect width="{W}" height="{H}"/></clipPath>')
        lines.append('</defs>')

        # Title
        _atan_labels = {"712": "atan(7/12)", "812": "atan(8/12)"}
        if label in _atan_labels:
            angle_str = f"{_atan_labels[label]} = {c17_deg:.4f}"
        else:
            angle_str = f"{c17_deg:.1f}"
        lines.append(f'<text x="{W/2}" y="22" text-anchor="middle" font-family="Arial"'
                     f' font-size="13" font-weight="bold">'
                     f'C17 = {angle_str}\u00b0'
                     f'   C15 = {math.degrees(s1):.1f}\u00b0'
                     f'   RMS = {rms:.4e} ft</text>')
        sub2 = (f"d(F14-F15)={d1:.3f}\u2032  d(F18-F1)={d2:.3f}\u2032"
                f"  Area={outline_area:.2f} sq ft"
                f"  rot={math.degrees(rot):+.3f}\u00b0"
                f"  \u0394=({dx:+.3f}\u2032, {dy:+.3f}\u2032)")
        lines.append(f'<text x="{W/2}" y="36" text-anchor="middle" font-family="Arial"'
                     f' font-size="9" fill="#666">{sub2}</text>')

        # Constraint residual subtitle
        rms_color = "#2E7D32" if rms < 0.01 else "#BF360C"
        lines.append(f'<text x="{W/2}" y="48" text-anchor="middle" font-family="Arial"'
                     f' font-size="8" fill="{rms_color}">'
                     f'T3\u2190F18F1: {resid[0]:.3e}'
                     f'   TC1\u2190F12F13: {resid[1]:.3e}'
                     f'   F16\u2192line: {resid[2]:.3e}'
                     f'   F17\u2192line: {resid[3]:.3e}</text>')

        # Three layers
        render_layer(lines, outer_segs, pts, outer_cfg, to_svg)
        render_layer(lines, inset_segs, pts, inset_cfg, to_svg)
        render_layer(lines, outline_segs, pts, outline_cfg, to_svg)

        # Area label: halfway between FC and midpoint of F9-F10
        _f9f10_mid = ((pts["F9"][0] + pts["F10"][0]) / 2,
                      (pts["F9"][1] + pts["F10"][1]) / 2)
        cx_a = (pts["FC"][0] + _f9f10_mid[0]) / 2
        cy_a = (pts["FC"][1] + _f9f10_mid[1]) / 2
        ax, ay = to_svg(cx_a, cy_a)
        lines.append(f'<text x="{ax:.1f}" y="{ay:.1f}" text-anchor="middle"'
                     f' font-family="Arial" font-size="12" fill="#333"'
                     f' font-weight="bold">{outline_area:.2f} sq ft</text>')
        lines.append(f'<text x="{ax:.1f}" y="{ay+14:.1f}" text-anchor="middle"'
                     f' font-family="Arial" font-size="9" fill="#666">'
                     f'(Outline enclosed area)</text>')

        # North arrow
        lines.append('<line x1="742" y1="560" x2="742" y2="524" stroke="#333"'
                     ' stroke-width="2" marker-end="url(#ah)"/>')
        lines.append('<text x="742" y="518" text-anchor="middle" font-family="Arial"'
                     ' font-size="13" font-weight="bold">N</text>')

        # Legend
        ly = 554
        lines.append(f'<rect x="40" y="{ly}" width="14" height="8" fill="#e8edf5"'
                     f' stroke="#333" stroke-width="1" opacity="0.3"/>')
        lines.append(f'<text x="60" y="{ly+7}" font-family="Arial" font-size="8"'
                     f' fill="#999">Outer path at 20%</text>')
        ly += 12
        lines.append(f'<rect x="40" y="{ly}" width="14" height="8"'
                     f' fill="rgba(255,152,0,0.35)" stroke="#BF360C" stroke-width="1"'
                     f' opacity="0.3"/>')
        lines.append(f'<text x="60" y="{ly+7}" font-family="Arial" font-size="8"'
                     f' fill="#999">Inset path at 20%</text>')
        ly += 12
        lines.append(f'<line x1="40" y1="{ly+4}" x2="54" y2="{ly+4}" stroke="#333"'
                     f' stroke-width="2.0"/>')
        lines.append(f'<text x="60" y="{ly+7}" font-family="Arial" font-size="8"'
                     f' fill="#333">Outline ({outline_area:.2f} sq ft)</text>')

        # Footer
        _now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        _git_desc = git_describe()
        lines.append(f'<text x="{W/2}" y="{H-2}" text-anchor="middle" font-family="Arial"'
                     f' font-size="7.5" fill="#999">Generated {_now} from {_git_desc}</text>')
        lines.append('</svg>')

        svg_content = "\n".join(lines)
        svg_path = os.path.join(out_dir, f"path_area_{label}.svg")
        with open(svg_path, "w", encoding="utf-8") as f:
            f.write(svg_content)
        count += 1

    print(f"\nGenerated {count} SVGs in {out_dir}")
