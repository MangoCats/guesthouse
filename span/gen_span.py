"""Generate N-S span measurement SVG.

Measures the north-south interior span at every inch of easting across
the W-series (inner wall) polygon.  Outputs a portrait-Letter SVG with
the span graph directly above a plan-view structure outline.

Output: span/span.svg
"""
import os, math

from shared.geometry import (
    path_polygon, vert_isects, compute_inner_walls,
    f8f9_corner_polyline, fmt_dist,
)
from shared.survey import compute_traverse, compute_three_arc, compute_inset
from shared.svg import git_describe
from floorplan.geometry import compute_outline_geometry, OutlineAnchors
from floorplan.constants import WALL_OUTER, F8F9_INNER_TURN_R


# ── geometry bootstrap ─────────────────────────────────────────

def _build_geometry():
    """Return (pts, outline_segs, outer_poly, inner_poly)."""
    pts, _ = compute_traverse()
    ai = compute_three_arc(pts)
    ins = compute_inset(pts, ai["R1"], ai["R2"], ai["R3"], ai["nE"], ai["nN"])
    pts.update(ins.pts_update)
    anch = OutlineAnchors(
        Pi2=pts["Pi2"], Pi3=pts["Pi3"], Ti3=pts["Ti3"],
        PiX=pts["PiX"], Pi5=pts["Pi5"], TC1=pts["TC1"], R1i=ins.R1i,
    )
    geo = compute_outline_geometry(anch)
    pts.update(geo.fp_pts)
    inner_segs = compute_inner_walls(geo.outline_segs, pts, WALL_OUTER, geo.radii)
    outer_poly = path_polygon(geo.outline_segs, pts)
    inner_poly = path_polygon(inner_segs, pts)
    # patch concave W8-W9 corner
    w_f8f9 = f8f9_corner_polyline(pts, WALL_OUTER, F8F9_INNER_TURN_R)
    w8, w9 = pts["W8"], pts["W9"]
    i8 = next(i for i, p in enumerate(inner_poly)
              if abs(p[0] - w8[0]) < 1e-9 and abs(p[1] - w8[1]) < 1e-9)
    i9 = next(i for i, p in enumerate(inner_poly)
              if i > i8 and abs(p[0] - w9[0]) < 1e-9 and abs(p[1] - w9[1]) < 1e-9)
    inner_poly[i8:i9 + 1] = w_f8f9
    return pts, geo.outline_segs, outer_poly, inner_poly


# ── span computation ───────────────────────────────────────────

def _compute_spans(inner_poly):
    """N-S span at every inch of easting.  Returns (eastings, spans) in feet."""
    e_min = min(p[0] for p in inner_poly)
    e_max = max(p[0] for p in inner_poly)
    inch = 1.0 / 12.0
    eastings, spans = [], []
    e = e_min
    while e <= e_max + 1e-9:
        ns = vert_isects(inner_poly, e)
        spans.append(max(ns) - min(ns) if len(ns) >= 2 else 0.0)
        eastings.append(e)
        e += inch
    return eastings, spans


# ── SVG generation ─────────────────────────────────────────────

def _generate_svg(pts, outer_poly, inner_poly):
    eastings, spans = _compute_spans(inner_poly)
    max_span = max(spans)
    max_span_e = eastings[spans.index(max_span)]

    # bounding box of all visible geometry
    all_pts = outer_poly + inner_poly
    E0 = min(p[0] for p in all_pts)           # westmost easting (X origin)
    N_MIN = min(p[1] for p in all_pts)
    N_MAX = max(p[1] for p in all_pts)

    # ── page layout (portrait US Letter @ 72 dpi) ──
    PW, PH = 612, 792
    ML, MR, MT, MB = 62, 20, 40, 25           # margins
    GAP = 32                                    # graph ↔ outline gap

    plot_w = PW - ML - MR                      # usable width
    AXIS_FT = 36.0                              # X-axis range (feet)
    xs = plot_w / AXIS_FT                       # pt per survey-foot (horiz)

    outline_h = (N_MAX - N_MIN) * xs            # 1:1 aspect outline
    graph_h = PH - MT - MB - GAP - outline_h
    if graph_h < 160:                           # safety: compress outline
        graph_h = 200
        outline_h = PH - MT - MB - GAP - graph_h

    y_max_ft = math.ceil(max_span)              # graph Y cap (feet)
    ys = graph_h / y_max_ft                     # pt per span-foot

    g_top = MT;            g_bot = MT + graph_h
    o_top = g_bot + GAP;   o_bot = o_top + outline_h

    # coordinate mappers
    def ex(e):  return ML + (e - E0) * xs       # easting → SVG x
    def sy(s):  return g_bot - s * ys            # span → graph SVG y
    def ny(n):  return o_bot - (n - N_MIN) * xs  # northing → outline SVG y

    o: list[str] = []
    o.append(f'<svg xmlns="http://www.w3.org/2000/svg"'
             f' width="{PW}" height="{PH}" viewBox="0 0 {PW} {PH}">')
    o.append('<rect width="100%" height="100%" fill="white"/>')

    # ── title ─────────────────────────────────────────────────
    o.append(f'<text x="{PW / 2}" y="{MT - 14}" text-anchor="middle"'
             f' font-family="Arial" font-size="13" font-weight="bold"'
             f' fill="#222">N\u2013S Interior Span</text>')

    # ── graph frame ───────────────────────────────────────────
    o.append(f'<rect x="{ML}" y="{g_top}" width="{plot_w}" height="{graph_h}"'
             f' fill="none" stroke="#bbb" stroke-width="0.5"/>')

    # Y grid + labels (every 2 ft if tall, else every 1 ft)
    y_step = 2 if y_max_ft > 14 else 1
    for ft in range(0, y_max_ft + 1, y_step):
        y = sy(ft)
        if ft > 0:
            o.append(f'<line x1="{ML}" y1="{y:.1f}" x2="{ML + plot_w}" y2="{y:.1f}"'
                     f' stroke="#eaeaea" stroke-width="0.3"/>')
        o.append(f'<line x1="{ML - 3}" y1="{y:.1f}" x2="{ML}" y2="{y:.1f}"'
                 f' stroke="#666" stroke-width="0.5"/>')
        o.append(f'<text x="{ML - 5}" y="{y + 2.5:.1f}" text-anchor="end"'
                 f' font-family="Arial" font-size="7" fill="#555">{ft}\'</text>')

    # X grid + ticks + labels (every 5 ft)
    for ft in range(0, int(AXIS_FT) + 1, 5):
        x = ML + ft * xs
        if 0 < ft < AXIS_FT:
            o.append(f'<line x1="{x:.1f}" y1="{g_top}" x2="{x:.1f}" y2="{g_bot}"'
                     f' stroke="#f0f0f0" stroke-width="0.3"/>')
        o.append(f'<line x1="{x:.1f}" y1="{g_bot}" x2="{x:.1f}" y2="{g_bot + 3}"'
                 f' stroke="#666" stroke-width="0.5"/>')
        o.append(f'<text x="{x:.1f}" y="{g_bot + 12}" text-anchor="middle"'
                 f' font-family="Arial" font-size="7" fill="#555">{ft}\'</text>')

    # span curve
    crv = [f"{ex(e):.1f},{sy(s):.1f}" for e, s in zip(eastings, spans) if s > 0]
    if crv:
        o.append(f'<polyline points="{" ".join(crv)}" fill="none"'
                 f' stroke="#1565C0" stroke-width="1.0" stroke-linejoin="round"/>')

    # max-span dashed line + label
    my = sy(max_span)
    o.append(f'<line x1="{ML}" y1="{my:.1f}" x2="{ML + plot_w}" y2="{my:.1f}"'
             f' stroke="#C62828" stroke-width="0.7" stroke-dasharray="6,3"/>')
    o.append(f'<text x="{ML + plot_w - 3}" y="{my - 3:.1f}" text-anchor="end"'
             f' font-family="Arial" font-size="8" fill="#C62828"'
             f' font-weight="bold">{fmt_dist(max_span)}</text>')

    # Y-axis title (rotated)
    lx, ly = 10, (g_top + g_bot) / 2
    o.append(f'<text x="{lx}" y="{ly:.1f}" text-anchor="middle"'
             f' font-family="Arial" font-size="9" fill="#333"'
             f' transform="rotate(-90,{lx},{ly:.1f})">Span (ft)</text>')

    # ── outline panel ─────────────────────────────────────────
    # F polygon (outer) — gray fill shows walls
    fp = " ".join(f"{ex(p[0]):.1f},{ny(p[1]):.1f}" for p in outer_poly)
    o.append(f'<polygon points="{fp}" fill="rgba(180,180,180,0.30)"'
             f' stroke="#555" stroke-width="0.6"/>')
    # W polygon (inner) — white fill cuts out interior
    wp = " ".join(f"{ex(p[0]):.1f},{ny(p[1]):.1f}" for p in inner_poly)
    o.append(f'<polygon points="{wp}" fill="white"'
             f' stroke="#1565C0" stroke-width="0.6"/>')

    # F-series dots + labels
    for i in range(21):
        if i == 4:
            continue
        nm = f"F{i}"
        if nm in pts:
            x, y = ex(pts[nm][0]), ny(pts[nm][1])
            o.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="1.2" fill="#555"/>')
            o.append(f'<text x="{x + 3:.1f}" y="{y - 2:.1f}" font-family="Arial"'
                     f' font-size="4.5" fill="#555">{nm}</text>')

    # W-series dots
    for i in range(21):
        nm = f"W{i}"
        if nm in pts:
            x, y = ex(pts[nm][0]), ny(pts[nm][1])
            o.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="1.0" fill="#1565C0"/>')

    # vertical reference line at max-span easting (across both panels)
    mx = ex(max_span_e)
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
    pts, _, outer_poly, inner_poly = _build_geometry()
    svg = _generate_svg(pts, outer_poly, inner_poly)
    outpath = os.path.join(os.path.dirname(os.path.abspath(__file__)), "span.svg")
    with open(outpath, 'w', encoding='utf-8') as f:
        f.write(svg)
    print(f"Wrote {outpath}")


if __name__ == "__main__":
    main()
