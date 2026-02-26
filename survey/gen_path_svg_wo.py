"""Generate path_area_wo.svg — path_area.svg plus O-series openings on the F-path."""
import os, sys, math, datetime

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

    # Campfire circle (same as path_area.svg)
    _circ_r_ft = 48.0 / 12.0 / 2.0
    _pob_r = pts["POB"]; _p2_r = pts["P2"]
    _dx_pp = _p2_r[0] - _pob_r[0]; _dy_pp = _p2_r[1] - _pob_r[1]
    _L_pp = math.sqrt(_dx_pp**2 + _dy_pp**2)
    _ux_pp, _uy_pp = _dx_pp / _L_pp, _dy_pp / _L_pp
    _tan_e = _pob_r[0] + 3.0 * _ux_pp
    _tan_n = _pob_r[1] + 3.0 * _uy_pp
    _ln_e, _ln_n = -_dy_pp / _L_pp, _dx_pp / _L_pp
    _circ_ce = _tan_e + _circ_r_ft * _ln_e
    _circ_cn = _tan_n + _circ_r_ft * _ln_n
    _ccx, _ccy = to_svg(_circ_ce, _circ_cn)
    _cr_svg = to_svg(_circ_ce + _circ_r_ft, _circ_cn)[0] - _ccx
    _fireplace_url = ("https://www.wayfair.com/outdoor/pdp/big-horn-outdoors-wellington"
                      "-819-h-concrete-outdoor-fireplace-famu1012.html")
    lines.append(f'<a href="{_fireplace_url}" target="_blank">')
    lines.append(f'<circle cx="{_ccx:.1f}" cy="{_ccy:.1f}" r="{_cr_svg:.1f}"'
                 f' fill="transparent" stroke="#333" stroke-width="1.5" cursor="pointer"/>')
    _icon_s = 3.0 * _cr_svg / _circ_r_ft
    _fire_dy = 0.15 * _icon_s
    lines.append(f'<g transform="translate({_ccx:.1f},{_ccy + _fire_dy:.1f}) scale({_icon_s:.2f})">')
    lines.append('<path d="M 0,0.05 C -0.22,-0.15 -0.32,-0.45 -0.18,-0.7'
                 ' C -0.1,-0.82 -0.04,-0.92 0,-0.95'
                 ' C 0.04,-0.92 0.1,-0.82 0.18,-0.7'
                 ' C 0.32,-0.45 0.22,-0.15 0,0.05 Z"'
                 ' fill="#B22222" opacity="0.85"/>')
    lines.append('<path d="M -0.06,0.03 C -0.2,-0.08 -0.28,-0.35 -0.15,-0.58'
                 ' C -0.08,-0.72 -0.03,-0.62 -0.02,-0.45'
                 ' C -0.01,-0.3 0.0,-0.1 -0.06,0.03 Z" fill="#D44000"/>')
    lines.append('<path d="M 0.06,0.03 C 0.2,-0.08 0.28,-0.35 0.15,-0.58'
                 ' C 0.08,-0.72 0.03,-0.62 0.02,-0.45'
                 ' C 0.01,-0.3 0.0,-0.1 0.06,0.03 Z" fill="#D44000"/>')
    lines.append('<path d="M 0,0.04 C -0.14,-0.12 -0.2,-0.4 -0.08,-0.62'
                 ' C -0.03,-0.72 0,-0.82 0,-0.82'
                 ' C 0,-0.82 0.03,-0.72 0.08,-0.62'
                 ' C 0.2,-0.4 0.14,-0.12 0,0.04 Z" fill="#E86100"/>')
    lines.append('<path d="M -0.04,0.0 C -0.15,-0.1 -0.18,-0.32 -0.08,-0.48'
                 ' C -0.02,-0.38 0.0,-0.15 -0.04,0.0 Z" fill="#FF7518"/>')
    lines.append('<path d="M 0.04,0.0 C 0.15,-0.1 0.18,-0.32 0.08,-0.48'
                 ' C 0.02,-0.38 0.0,-0.15 0.04,0.0 Z" fill="#FF7518"/>')
    lines.append('<path d="M 0,0.02 C -0.09,-0.08 -0.12,-0.3 -0.04,-0.5'
                 ' C -0.01,-0.58 0,-0.65 0,-0.65'
                 ' C 0,-0.65 0.01,-0.58 0.04,-0.5'
                 ' C 0.12,-0.3 0.09,-0.08 0,0.02 Z" fill="#FFA500"/>')
    lines.append('<path d="M 0,0.0 C -0.05,-0.06 -0.07,-0.2 -0.02,-0.38'
                 ' C 0,-0.44 0,-0.44'
                 ' C 0.02,-0.38 0.07,-0.2 0.05,-0.06 Z" fill="#FFD700"/>')
    lines.append('<path d="M 0,-0.02 C -0.025,-0.08 -0.03,-0.15 0,-0.28'
                 ' C 0.03,-0.15 0.025,-0.08 0,-0.02 Z" fill="#FFF176"/>')
    _lw = 0.065
    lines.append(f'<line x1="-0.4" y1="0.32" x2="0.12" y2="0.06"'
                 f' stroke="#6B3410" stroke-width="{_lw}" stroke-linecap="round"/>')
    lines.append(f'<line x1="0.4" y1="0.32" x2="-0.12" y2="0.06"'
                 f' stroke="#8B4513" stroke-width="{_lw}" stroke-linecap="round"/>')
    lines.append(f'<line x1="-0.3" y1="0.25" x2="0.3" y2="0.25"'
                 f' stroke="#A0522D" stroke-width="{_lw * 0.85:.4f}" stroke-linecap="round"/>')
    for _le, _ln in [(-0.4, 0.32), (0.4, 0.32), (-0.3, 0.25), (0.3, 0.25)]:
        lines.append(f'<circle cx="{_le}" cy="{_ln}" r="0.028" fill="#A0522D" stroke="#5C2D0E" stroke-width="0.008"/>')
    lines.append('<ellipse cx="0" cy="0.18" rx="0.18" ry="0.06" fill="#FF4500" opacity="0.3"/>')
    lines.append('</g>')
    lines.append('</a>')

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
