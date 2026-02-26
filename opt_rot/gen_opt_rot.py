"""Generate path_area SVGs for C17 sweep angles from 30.0 to 42.0 degrees."""
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
from floorplan.constants import CORNER_NE_R, F11AB_TARGET
from survey.gen_path_svg import (
    compute_all, render_layer, build_outline_cfg,
    outer_cfg, inset_cfg,
)

# ============================================================
# Constants
# ============================================================
_PI_2 = math.pi / 2
_C10_SWEEP = math.pi / 2 - math.atan(1.0 / 3.0)
_C11_SWEEP = math.pi / 2 - math.atan(1.0 / 3.0)
_C13_SWEEP = math.atan(1.0 / 3.0)

_R_a1 = CORNER_NE_R
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
# Main
# ============================================================
if __name__ == "__main__":
    # Compute survey/inset data once (shared across all angles)
    data = compute_all()
    base_pts = dict(data["pts"])
    to_svg = data["to_svg"]
    outer_segs = data["outer_segs"]
    inset_segs = data["inset_segs"]
    outer_area = data["outer_area"]
    inset_area = data["inset_area"]

    out_dir = os.path.dirname(os.path.abspath(__file__))
    count = 0

    for c17_deg_x10 in range(300, 421, 10):
        c17_deg = c17_deg_x10 / 10.0
        d1, d2, s1, s2, chain = solve_for_angle(c17_deg)

        print(f"C17={c17_deg:5.1f}\u00b0  C15={math.degrees(s1):5.1f}\u00b0  "
              f"d_F14_F15={d1:.4f}'  d_F18_F1={d2:.4f}'")

        if d1 <= 0 or d2 <= 0:
            print(f"  WARNING: non-positive length, skipping")
            continue

        # Walk chain to get F/C points
        fp_pts = walk_chain(chain)

        # Merge with survey points (copy base, override with new outline)
        pts = dict(base_pts)
        pts.update(fp_pts)

        # Build outline segments and radii
        outline_segs = build_outline_segs(chain)
        radii = build_radii(chain)

        # Compute outline polygon and area
        outer_poly = path_polygon(outline_segs, pts)
        outline_area = poly_area(outer_poly)

        # Build outline layer config
        outline_cfg = build_outline_cfg(outline_segs, pts, radii)

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
        lines.append(f'<text x="{W/2}" y="26" text-anchor="middle" font-family="Arial"'
                     f' font-size="14" font-weight="bold">'
                     f'C17 Sweep = {c17_deg:.1f}\u00b0'
                     f'  (C15 = {math.degrees(s1):.1f}\u00b0)</text>')
        lines.append(f'<text x="{W/2}" y="42" text-anchor="middle" font-family="Arial"'
                     f' font-size="10" fill="#666">'
                     f"d(F14-F15) = {d1:.3f}'   d(F18-F1) = {d2:.3f}'"
                     f'   Area = {outline_area:.2f} sq ft</text>')

        # Three layers
        render_layer(lines, outer_segs, pts, outer_cfg, to_svg)
        render_layer(lines, inset_segs, pts, inset_cfg, to_svg)
        render_layer(lines, outline_segs, pts, outline_cfg, to_svg)

        # North arrow
        lines.append('<line x1="742" y1="560" x2="742" y2="524" stroke="#333"'
                     ' stroke-width="2" marker-end="url(#ah)"/>')
        lines.append('<text x="742" y="518" text-anchor="middle" font-family="Arial"'
                     ' font-size="13" font-weight="bold">N</text>')

        # Legend
        ly = 550
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
        svg_path = os.path.join(out_dir, f"path_area_{c17_deg_x10:03d}.svg")
        with open(svg_path, "w", encoding="utf-8") as f:
            f.write(svg_content)
        count += 1

    print(f"\nGenerated {count} SVGs in {out_dir}")
