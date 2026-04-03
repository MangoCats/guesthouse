"""Generate minimum-span rotation SVG.

Finds the rotation angle that minimises the maximum interior span via a
three-pass search (0-175 deg coarse, 0.5 deg fine, 0.1 deg superfine).
If the result exceeds 90 deg, 180 deg is subtracted to normalise.

Outputs a single portrait-Letter SVG with:
  - Three-curve graph (blue total, green south-to-IW, cyan IW-to-north)
  - Dashed line at max total span (red)
  - Dashed line at max unsupported span (orange) = max of green / cyan
  - Rotated structure outline beneath the graph

Output: span/span_min.svg
"""
import os, sys, math

# Ensure project root is on sys.path for package imports
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

from shared.geometry import fmt_dist
from shared.svg import git_describe
from span._common import (
    build_geometry, extract_iw_centerlines,
    compute_rotation_data, find_min_span_angle, fmt_angle,
    render_y_grid, render_x_grid, render_span_curves,
)


# ── SVG generation ─────────────────────────────────────────────

def _generate_svg(pts, outer_poly, inner_poly, layout, roof_poly):
    iw_cls = extract_iw_centerlines(layout)
    cx = sum(p[0] for p in inner_poly) / len(inner_poly)
    cy = sum(p[1] for p in inner_poly) / len(inner_poly)

    angle, _ = find_min_span_angle(inner_poly, iw_cls, cx, cy)
    d = compute_rotation_data(angle, outer_poly, inner_poly, iw_cls, cx, cy,
                              roof_poly)

    max_span = d['max_span']
    max_e = d['max_e']
    max_unsup = max(max(d['s_spans']), max(d['n_spans']))
    max_unsup_e_s = d['eastings'][d['s_spans'].index(max(d['s_spans']))]
    max_unsup_e_n = d['eastings'][d['n_spans'].index(max(d['n_spans']))]
    max_roof_span = d['max_roof_span']

    ang_str = fmt_angle(angle)

    # bounding box
    E0 = d['e_min']
    N_MIN, N_MAX = d['n_min'], d['n_max']

    # page layout (portrait US Letter @ 72 dpi)
    PW, PH = 612, 792
    ML, MR, MT, MB = 62, 20, 40, 25
    GAP = 32

    plot_w = PW - ML - MR
    geo_w = d['e_max'] - d['e_min']
    xs = plot_w / geo_w                             # pt per survey-foot (horiz)

    outline_h = (N_MAX - N_MIN) * xs                # 1:1 aspect outline
    graph_h = PH - MT - MB - GAP - outline_h
    if graph_h < 160:
        graph_h = 200
        outline_h = PH - MT - MB - GAP - graph_h

    y_max_ft = math.ceil(max(max_span, max_roof_span))  # graph Y cap (feet)
    ys = graph_h / y_max_ft                         # pt per span-foot

    g_top = MT;            g_bot = MT + graph_h
    o_top = g_bot + GAP;   o_bot = o_top + outline_h

    # coordinate mappers
    def ex(e):  return ML + (e - E0) * xs
    def sy(s):  return g_bot - s * ys
    def ny(n):  return o_bot - (n - N_MIN) * xs

    o: list[str] = []
    o.append(f'<svg xmlns="http://www.w3.org/2000/svg"'
             f' width="{PW}" height="{PH}" viewBox="0 0 {PW} {PH}">')
    o.append('<rect width="100%" height="100%" fill="white"/>')

    # title
    o.append(f'<text x="{PW / 2}" y="{MT - 14}" text-anchor="middle"'
             f' font-family="Arial" font-size="13" font-weight="bold"'
             f' fill="#222">Minimum-Span Rotation: {ang_str}\u00b0</text>')

    # graph frame
    o.append(f'<rect x="{ML}" y="{g_top}" width="{plot_w}" height="{graph_h}"'
             f' fill="none" stroke="#bbb" stroke-width="0.5"/>')

    # Y grid + labels
    render_y_grid(o, ML, plot_w, g_bot, y_max_ft, ys)

    # X grid + ticks + labels (every 5 ft)
    render_x_grid(o, ML, plot_w, g_top, g_bot, xs, geo_w)

    # span curves: roof (grey), south-to-IW (green), IW-to-north (cyan), total (blue)
    render_span_curves(o, d['eastings'],
                       (d['roof_spans'], d['s_spans'], d['n_spans'], d['spans']),
                       ex, sy)

    # max roof span dashed line + label (grey)
    if max_roof_span > 0:
        ry = sy(max_roof_span)
        o.append(f'<line x1="{ML}" y1="{ry:.1f}" x2="{ML + plot_w}" y2="{ry:.1f}"'
                 f' stroke="#999" stroke-width="0.7" stroke-dasharray="6,3"/>')
        ry_label_y = ry + 10 if ry < sy(max_span) - 6 else ry - 3
        o.append(f'<text x="{ML + plot_w - 3}" y="{ry_label_y:.1f}" text-anchor="end"'
                 f' font-family="Arial" font-size="8" fill="#999"'
                 f' font-weight="bold">max roof: {fmt_dist(max_roof_span)}</text>')

    # max total span dashed line + label (red)
    my = sy(max_span)
    o.append(f'<line x1="{ML}" y1="{my:.1f}" x2="{ML + plot_w}" y2="{my:.1f}"'
             f' stroke="#C62828" stroke-width="0.7" stroke-dasharray="6,3"/>')
    o.append(f'<text x="{ML + plot_w - 3}" y="{my - 3:.1f}" text-anchor="end"'
             f' font-family="Arial" font-size="8" fill="#C62828"'
             f' font-weight="bold">max span: {fmt_dist(max_span)}</text>')

    # max unsupported span dashed line + label (orange)
    uy = sy(max_unsup)
    o.append(f'<line x1="{ML}" y1="{uy:.1f}" x2="{ML + plot_w}" y2="{uy:.1f}"'
             f' stroke="#E65100" stroke-width="0.7" stroke-dasharray="6,3"/>')
    label_y = uy + 11 if uy < my - 6 else uy - 3
    o.append(f'<text x="{ML + plot_w - 3}" y="{label_y:.1f}" text-anchor="end"'
             f' font-family="Arial" font-size="8" fill="#E65100"'
             f' font-weight="bold">max unsupported: {fmt_dist(max_unsup)}</text>')

    # Y-axis title (rotated)
    lx, ly = 10, (g_top + g_bot) / 2
    o.append(f'<text x="{lx}" y="{ly:.1f}" text-anchor="middle"'
             f' font-family="Arial" font-size="9" fill="#333"'
             f' transform="rotate(-90,{lx},{ly:.1f})">Span (ft)</text>')

    # legend
    leg_x = ML + 8
    leg_y = g_top + 10
    for i, (color, label) in enumerate([
        ("#1565C0", "Total span"),
        ("#2E7D32", "Bottom \u2192 central IW mid"),
        ("#00ACC1", "Central IW mid \u2192 top"),
        ("#999", "Roof span"),
        ("#C62828", "Max total span"),
        ("#E65100", "Max unsupported span"),
    ]):
        ly_i = leg_y + i * 11
        dash = "4,2" if i >= 4 else "none"
        o.append(f'<line x1="{leg_x}" y1="{ly_i}" x2="{leg_x + 14}" y2="{ly_i}"'
                 f' stroke="{color}" stroke-width="1.2"'
                 f' stroke-dasharray="{dash}"/>')
        o.append(f'<text x="{leg_x + 17}" y="{ly_i + 3}" font-family="Arial"'
                 f' font-size="6" fill="#444">{label}</text>')

    # ── outline panel ─────────────────────────────────────────
    fp = " ".join(f"{ex(p[0]):.1f},{ny(p[1]):.1f}" for p in d['r_outer'])
    o.append(f'<polygon points="{fp}" fill="rgba(180,180,180,0.30)"'
             f' stroke="#555" stroke-width="0.6"/>')

    wp = " ".join(f"{ex(p[0]):.1f},{ny(p[1]):.1f}" for p in d['r_inner'])
    o.append(f'<polygon points="{wp}" fill="white"'
             f' stroke="#1565C0" stroke-width="0.6"/>')

    # IW centerlines on outline
    for cl in d['r_cls']:
        x1, y1 = ex(cl[0][0]), ny(cl[0][1])
        x2, y2 = ex(cl[1][0]), ny(cl[1][1])
        o.append(f'<line x1="{x1:.1f}" y1="{y1:.1f}"'
                 f' x2="{x2:.1f}" y2="{y2:.1f}"'
                 f' stroke="#666" stroke-width="0.5"'
                 f' stroke-dasharray="2,1"/>')

    # vertical reference lines at max-span and max-unsupported eastings
    mx = ex(max_e)
    o.append(f'<line x1="{mx:.1f}" y1="{g_top}" x2="{mx:.1f}" y2="{o_bot}"'
             f' stroke="#C62828" stroke-width="0.4" stroke-dasharray="2,3" opacity="0.35"/>')

    # version stamp
    ver = git_describe()
    o.append(f'<text x="{PW - 5}" y="{PH - 5}" text-anchor="end"'
             f' font-family="Arial" font-size="5" fill="#ccc">{ver}</text>')

    o.append('</svg>')
    return '\n'.join(o)


# ── entry point ────────────────────────────────────────────────

def main():
    pts, _, outer_poly, inner_poly, layout, roof_poly = build_geometry()
    svg = _generate_svg(pts, outer_poly, inner_poly, layout, roof_poly)
    outpath = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           "span_min.svg")
    with open(outpath, 'w', encoding='utf-8') as f:
        f.write(svg)
    print(f"Wrote {outpath}")


if __name__ == "__main__":
    main()
