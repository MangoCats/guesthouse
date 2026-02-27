"""Generate span-vs-rotation analysis SVG.

Rotates the structure from 0 to 175 deg in 5 deg steps (36 panels).
For each rotation, measures the span from bottom to top at every inch
of the rotated x-axis.  Three curves per panel:
  blue  — total span
  green — bottom surface to nearest 6" IW centerline (IW1, IW2, IW8)
  cyan  — uppermost 6" IW centerline to top surface
When no IW centerline is intersected, green and cyan equal the full span.

The rotation with the lowest maximum total span is highlighted.
A fine search (0.5 deg steps, +/-4.5 deg around the coarse minimum) and a
superfine search (0.1 deg steps, +/-0.9 deg around the fine minimum) refine
the exact minimum angle.

Output: span/span_minmax.svg
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

    # coarse rotations (5 deg steps)
    all_data = [compute_rotation_data(a, outer_poly, inner_poly, iw_cls, cx, cy,
                                      roof_poly)
                for a in range(0, 176, 5)]

    # coarse best angle (for highlighting in panel grid)
    mm_span = min(d['max_span'] for d in all_data)
    mm_angle = next(d['angle'] for d in all_data if d['max_span'] == mm_span)

    # refined minimum via three-pass search (unnormalised for panel rendering)
    sf_best_angle, sf_best_span = find_min_span_angle(
        inner_poly, iw_cls, cx, cy, normalize=False)

    # full data for the refined minimum angle
    refined = compute_rotation_data(
        sf_best_angle, outer_poly, inner_poly, iw_cls, cx, cy, roof_poly)

    ref_ang_str = fmt_angle(sf_best_angle)

    # global bounds (include refined panel)
    all_panels = all_data + [refined]
    g_e_min = min(d['e_min'] for d in all_panels)
    g_e_max = max(d['e_max'] for d in all_panels)
    g_max_span = max(d['max_span'] for d in all_panels)
    g_max_roof = max(d['max_roof_span'] for d in all_panels)

    # page layout (portrait Letter width)
    PW = 612
    ML, MR, MT = 62, 20, 40
    PGAP = 20
    pw = PW - ML - MR
    gw = g_e_max - g_e_min
    xs = pw / gw                                  # pt per survey-foot
    GH = 140                                      # graph height (pt)
    y_cap = math.ceil(max(g_max_span, g_max_roof))  # Y-axis cap (ft)
    ys = GH / y_cap                               # pt per span-foot

    # total SVG height (coarse panels + refined panel)
    th = MT + 16                                  # title + legend
    for d in all_panels:
        oh = (d['n_max'] - d['n_min']) * xs
        th += 18 + GH + 14 + 8 + oh + PGAP
    th += 20                                      # bottom margin

    o: list[str] = []
    o.append(f'<svg xmlns="http://www.w3.org/2000/svg"'
             f' width="{PW}" height="{th:.0f}" viewBox="0 0 {PW} {th:.0f}">')
    o.append('<rect width="100%" height="100%" fill="white"/>')

    # title
    o.append(f'<text x="{PW / 2}" y="{MT - 12}" text-anchor="middle"'
             f' font-family="Arial" font-size="13" font-weight="bold"'
             f' fill="#222">Span vs. Rotation (0\u00b0\u2013175\u00b0)'
             f' \u2014 min max: {fmt_dist(sf_best_span)}'
             f' at {ref_ang_str}\u00b0</text>')

    # legend
    ly = MT + 2
    for i, (col, lab) in enumerate([
        ("#1565C0", "Total span"),
        ("#2E7D32", "Bottom \u2192 nearest IW mid"),
        ("#00ACC1", "Uppermost IW mid \u2192 top"),
        ("#999", "Roof span"),
    ]):
        lx = ML + i * 120
        o.append(f'<line x1="{lx}" y1="{ly}" x2="{lx + 14}" y2="{ly}"'
                 f' stroke="{col}" stroke-width="1.2"/>')
        o.append(f'<text x="{lx + 17}" y="{ly + 3}" font-family="Arial"'
                 f' font-size="6.5" fill="#444">{lab}</text>')

    yc = MT + 16

    def ex(e):
        return ML + (e - g_e_min) * xs

    # ── panel renderer (nested to share layout vars) ──────────
    def _panel(d, title_text, title_color, title_weight, highlight):
        nonlocal yc
        oh = (d['n_max'] - d['n_min']) * xs

        yc += 14
        o.append(f'<text x="{ML}" y="{yc}" font-family="Arial"'
                 f' font-size="9" fill="{title_color}"'
                 f' font-weight="{title_weight}">{title_text}</text>')
        yc += 4

        gt = yc                                   # graph top
        gb = gt + GH                              # graph bottom

        if highlight:
            ph = GH + 14 + 8 + oh
            o.append(f'<rect x="{ML - 4}" y="{gt - 2}" width="{pw + 8}"'
                     f' height="{ph + 4:.0f}" fill="rgba(198,40,40,0.04)"'
                     f' stroke="#C62828" stroke-width="0.5" rx="3"/>')

        o.append(f'<rect x="{ML}" y="{gt}" width="{pw}" height="{GH}"'
                 f' fill="none" stroke="#bbb" stroke-width="0.5"/>')

        render_y_grid(o, ML, pw, gb, y_cap, ys, font_size="6")

        render_x_grid(o, ML, pw, gt, gb, xs, gw,
                      font_size="5", label_y_offset=10)

        def _sy(s): return gb - s * ys
        render_span_curves(o, d['eastings'],
                           (d['roof_spans'], d['s_spans'], d['n_spans'], d['spans']),
                           ex, _sy, stroke_widths=(0.6, 0.6, 0.6, 0.8))

        if d['max_roof_span'] > 0:
            ry = gb - d['max_roof_span'] * ys
            o.append(f'<line x1="{ML}" y1="{ry:.1f}"'
                     f' x2="{ML + pw}" y2="{ry:.1f}"'
                     f' stroke="#999" stroke-width="0.5"'
                     f' stroke-dasharray="4,2"/>')

        my = gb - d['max_span'] * ys
        o.append(f'<line x1="{ML}" y1="{my:.1f}"'
                 f' x2="{ML + pw}" y2="{my:.1f}"'
                 f' stroke="#C62828" stroke-width="0.5"'
                 f' stroke-dasharray="4,2"/>')

        yc = gb + 14

        ot = yc + 8
        ob = ot + oh
        nm = d['n_min']

        fp = " ".join(f"{ex(p[0]):.1f},{ob - (p[1] - nm) * xs:.1f}"
                      for p in d['r_outer'])
        o.append(f'<polygon points="{fp}" fill="rgba(180,180,180,0.25)"'
                 f' stroke="#555" stroke-width="0.4"/>')

        wp = " ".join(f"{ex(p[0]):.1f},{ob - (p[1] - nm) * xs:.1f}"
                      for p in d['r_inner'])
        o.append(f'<polygon points="{wp}" fill="white"'
                 f' stroke="#1565C0" stroke-width="0.4"/>')

        for cl in d['r_cls']:
            x1, y1 = ex(cl[0][0]), ob - (cl[0][1] - nm) * xs
            x2, y2 = ex(cl[1][0]), ob - (cl[1][1] - nm) * xs
            o.append(f'<line x1="{x1:.1f}" y1="{y1:.1f}"'
                     f' x2="{x2:.1f}" y2="{y2:.1f}"'
                     f' stroke="#666" stroke-width="0.5"'
                     f' stroke-dasharray="2,1"/>')

        mx = ex(d['max_e'])
        o.append(f'<line x1="{mx:.1f}" y1="{gt}"'
                 f' x2="{mx:.1f}" y2="{ob}"'
                 f' stroke="#C62828" stroke-width="0.3"'
                 f' stroke-dasharray="2,3" opacity="0.3"/>')

        yc = ob + PGAP

    # ── render coarse panels ──────────────────────────────────
    for d in all_data:
        ang = d['angle']
        is_min = ang == mm_angle
        if is_min:
            tag = (f" \u2190 MIN (refined: {fmt_dist(sf_best_span)}"
                   f" at {ref_ang_str}\u00b0)")
        else:
            tag = ""
        _panel(d,
               f'{ang}\u00b0  max: {fmt_dist(d["max_span"])}{tag}',
               "#C62828" if is_min else "#333",
               "bold" if is_min else "normal",
               is_min)

    # ── render refined minimum panel ──────────────────────────
    _panel(refined,
           f'{ref_ang_str}\u00b0  max: {fmt_dist(refined["max_span"])}'
           f' \u2014 REFINED MIN',
           "#C62828", "bold", True)

    # version stamp
    ver = git_describe()
    o.append(f'<text x="{PW - 5}" y="{th - 5}" text-anchor="end"'
             f' font-family="Arial" font-size="5" fill="#ccc">{ver}</text>')

    o.append('</svg>')
    return '\n'.join(o)


# ── entry point ────────────────────────────────────────────────

def main():
    pts, _, outer_poly, inner_poly, layout, roof_poly = build_geometry()
    svg = _generate_svg(pts, outer_poly, inner_poly, layout, roof_poly)
    outpath = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           "span_minmax.svg")
    with open(outpath, 'w', encoding='utf-8') as f:
        f.write(svg)
    print(f"Wrote {outpath}")


if __name__ == "__main__":
    main()
