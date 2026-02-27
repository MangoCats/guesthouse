"""Generate path_area_wo.svg — path_area.svg plus O-series openings on the F-path."""
import os, sys, datetime

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

from shared.svg import make_svg_transform, W, H, git_describe, svg_polygon_pts
from survey.gen_path_svg import (
    compute_all, render_layer, render_floorplan,
    build_outline_cfg, outer_cfg, inset_cfg,
)
from floorplan.openings import compute_outer_openings

OPENING_FILL = 'rgb(220,235,255)'
OPENING_STROKE = '#4682B4'


def render_openings(lines, openings, to_svg):
    """Draw O-series opening polygons and labels on the outline."""
    for o in openings:
        svg_pts = svg_polygon_pts(o.poly, to_svg)
        lines.append(f'<polygon points="{svg_pts}" fill="{OPENING_FILL}"'
                     f' stroke="{OPENING_STROKE}" stroke-width="1.0"/>')
        # Label at polygon centroid
        cx = sum(p[0] for p in o.poly) / 4
        cy = sum(p[1] for p in o.poly) / 4
        sx, sy = to_svg(cx, cy)
        lines.append(f'<text x="{sx:.1f}" y="{sy + 3:.1f}" text-anchor="middle"'
                     f' font-family="Arial" font-size="7" fill="{OPENING_STROKE}"'
                     f' font-weight="bold">{o.name}</text>')


if __name__ == "__main__":
    data = compute_all()
    pts = data["pts"]; to_svg = data["to_svg"]
    outer_segs = data["outer_segs"]; inset_segs = data["inset_segs"]
    outline_segs = data["outline_segs"]
    outer_area = data["outer_area"]; inset_area = data["inset_area"]
    outline_area = data["outline_area"]
    radii = data["radii"]; layout = data["layout"]

    outline_cfg = build_outline_cfg(outline_segs, pts, radii)

    lines = []
    lines.append(f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}">')
    lines.append(f'<rect width="{W}" height="{H}" fill="white"/>')
    lines.append('<defs>')
    lines.append('  <marker id="ah" markerWidth="8" markerHeight="6" refX="8" refY="3" orient="auto">'
                 '<polygon points="0 0, 8 3, 0 6" fill="#333"/></marker>')
    lines.append(f'  <clipPath id="page"><rect width="{W}" height="{H}"/></clipPath>')
    lines.append('</defs>')
    lines.append(f'<text x="{W/2}" y="30" text-anchor="middle" font-family="Arial" font-size="14"'
                 f' font-weight="bold">Site Path \u2014 Outline + Openings</text>')

    render_layer(lines, outer_segs, pts, outer_cfg, to_svg)
    render_layer(lines, inset_segs, pts, inset_cfg, to_svg)
    render_layer(lines, outline_segs, pts, outline_cfg, to_svg)
    render_floorplan(lines, to_svg, pts, data["outer_poly"], data["inner_poly"],
                     data["inner_segs"], layout)

    # O-series openings
    outer_openings = compute_outer_openings(pts, layout)
    render_openings(lines, outer_openings, to_svg)

    # Area label
    _w9w10_mid = ((pts["W9"][0] + pts["W10"][0]) / 2, (pts["W9"][1] + pts["W10"][1]) / 2)
    cx_o = (pts["FC"][0] + _w9w10_mid[0]) / 2
    cy_o = (pts["FC"][1] + _w9w10_mid[1]) / 2
    sx, sy = to_svg(cx_o, cy_o)
    lines.append(f'<text x="{sx:.1f}" y="{sy:.1f}" text-anchor="middle" font-family="Arial"'
                 f' font-size="12" fill="#333" font-weight="bold">{outline_area:.2f} sq ft</text>')
    lines.append(f'<text x="{sx:.1f}" y="{sy+14:.1f}" text-anchor="middle" font-family="Arial"'
                 f' font-size="9" fill="#666">(Outline enclosed area)</text>')

    # North arrow
    lines.append('<line x1="742" y1="560" x2="742" y2="524" stroke="#333" stroke-width="2" marker-end="url(#ah)"/>')
    lines.append('<text x="742" y="518" text-anchor="middle" font-family="Arial" font-size="13" font-weight="bold">N</text>')

    # Legend
    ly = 550
    lines.append(f'<rect x="40" y="{ly}" width="14" height="8" fill="#e8edf5" stroke="#333" stroke-width="1" opacity="0.3"/>')
    lines.append(f'<text x="60" y="{ly+7}" font-family="Arial" font-size="8" fill="#999">Outer path at 20% ({outer_area:.2f} sq ft)</text>')
    ly += 12
    lines.append(f'<rect x="40" y="{ly}" width="14" height="8" fill="rgba(255,152,0,0.35)" stroke="#BF360C" stroke-width="1" opacity="0.3"/>')
    lines.append(f'<text x="60" y="{ly+7}" font-family="Arial" font-size="8" fill="#999">Inset path at 20% ({inset_area:.2f} sq ft)</text>')
    ly += 12
    lines.append(f'<line x1="40" y1="{ly+4}" x2="54" y2="{ly+4}" stroke="#333" stroke-width="2.0"/>')
    lines.append(f'<text x="60" y="{ly+7}" font-family="Arial" font-size="8" fill="#333">Outline path ({outline_area:.2f} sq ft)</text>')
    ly += 12
    lines.append(f'<rect x="40" y="{ly}" width="14" height="8" fill="{OPENING_FILL}" stroke="{OPENING_STROKE}" stroke-width="1"/>')
    lines.append(f'<text x="60" y="{ly+7}" font-family="Arial" font-size="8" fill="{OPENING_STROKE}">Openings (O1\u2013O11)</text>')

    # Footer
    _now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    _git_desc = git_describe()
    lines.append(f'<text x="{W/2}" y="{H-2}" text-anchor="middle" font-family="Arial" font-size="7.5"'
                 f' fill="#999">Generated {_now} from {_git_desc}</text>')
    lines.append('</svg>')

    svg_content = "\n".join(lines)
    svg_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "path_area_wo.svg")
    with open(svg_path, "w", encoding="utf-8") as f:
        f.write(svg_content)

    print(f"\nSVG written to path_area_wo.svg")
    print(f"Openings rendered: {', '.join(o.name for o in outer_openings)}")
