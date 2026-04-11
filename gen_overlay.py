"""One-off overlay SVG: exterior wall outlines for Base.db vs Mark9.db."""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))

from app.database import get_constants_dict, get_outline_chain
from app.gen_provider import compute_native_geometry
from shared.geometry import path_polygon

CONFIGS_DIR = os.path.join(os.path.dirname(__file__), "app", "configs")
OUTPUT = os.path.join(os.path.dirname(__file__), "overlay.svg")

DBS = [
    ("Base",  os.path.join(CONFIGS_DIR, "Base.db"),  "#4a90d9", "6"),
    ("Mark9", os.path.join(CONFIGS_DIR, "Mark9.db"), "#e05c3a", "3"),
]

SCALE = 30          # pixels per foot
MARGIN = 40         # px


def outline_polygon(db_path):
    constants = get_constants_dict(db_path)
    chain_rows = get_outline_chain(db_path)
    pts, outline_segs, _inner, _radii = compute_native_geometry(
        constants, chain_rows=chain_rows, db_path=db_path)
    return path_polygon(outline_segs, pts)   # list of (E, N)


def to_svg_pt(E, N, min_E, max_N):
    """Convert (E, N) world coords to SVG pixel coords (flip N axis)."""
    x = (E - min_E) * SCALE + MARGIN
    y = (max_N - N) * SCALE + MARGIN
    return x, y


def poly_to_path(poly, min_E, max_N):
    pts = [to_svg_pt(E, N, min_E, max_N) for E, N in poly]
    d = "M {:.2f},{:.2f}".format(*pts[0])
    d += " ".join(" L {:.2f},{:.2f}".format(x, y) for x, y in pts[1:])
    d += " Z"
    return d


def main():
    polys = []
    for label, db_path, color, stroke_w in DBS:
        poly = outline_polygon(db_path)
        polys.append((label, poly, color, stroke_w))
        print(f"  {label}: {len(poly)} outline points")

    all_pts = [pt for _, poly, _, _ in polys for pt in poly]
    min_E = min(p[0] for p in all_pts)
    max_E = max(p[0] for p in all_pts)
    min_N = min(p[1] for p in all_pts)
    max_N = max(p[1] for p in all_pts)

    W = (max_E - min_E) * SCALE + 2 * MARGIN
    H = (max_N - min_N) * SCALE + 2 * MARGIN

    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{W:.0f}" height="{H:.0f}" '
        f'viewBox="0 0 {W:.0f} {H:.0f}">',
        '  <rect width="100%" height="100%" fill="#1e1e2e"/>',
        '  <!-- Overlay: Base (blue) vs Mark9 (orange) -->',
    ]

    for label, poly, color, stroke_w in polys:
        d = poly_to_path(poly, min_E, max_N)
        lines.append(
            f'  <path d="{d}" fill="{color}" fill-opacity="0.12" '
            f'stroke="{color}" stroke-width="{stroke_w}" stroke-linejoin="round"/>'
        )

    # Legend
    legend_x, legend_y = MARGIN, MARGIN - 10
    for i, (label, _, color, _) in enumerate(DBS):
        lx = legend_x + i * 110
        lines.append(
            f'  <rect x="{lx}" y="{legend_y - 14}" width="20" height="14" '
            f'fill="{color}" fill-opacity="0.5" stroke="{color}" stroke-width="1.5"/>'
        )
        lines.append(
            f'  <text x="{lx + 26}" y="{legend_y}" font-family="monospace" '
            f'font-size="13" fill="#cdd6f4">{label}</text>'
        )

    lines.append('</svg>')

    svg_text = "\n".join(lines)
    with open(OUTPUT, "w") as f:
        f.write(svg_text)
    print(f"Written: {OUTPUT}  ({W:.0f}×{H:.0f} px)")


if __name__ == "__main__":
    main()
