"""Generate N-S span measurement SVG.

Measures the north-south interior span at every inch of easting across
the W-series (inner wall) polygon.  Outputs a portrait-Letter SVG with
the span graph directly above a plan-view structure outline.

Three curves:
  blue  — total N-S span (south W surface to north W surface)
  green — south W surface to midline of IW8 or IW1 (whichever is hit first)
  cyan  — midline of IW8/IW1 to north W surface
When no 6" IW is intersected, green and cyan both equal the full span.

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
from floorplan.layout import compute_interior_layout
from floorplan.roof import compute_roof_geometry, roof_polyline


# ── geometry bootstrap ─────────────────────────────────────────

def _build_geometry():
    """Return (pts, outline_segs, outer_poly, inner_poly, layout)."""
    pts, _ = compute_traverse()
    ai = compute_three_arc(pts)
    ins = compute_inset(pts, ai["R1"], ai["R2"], ai["R3"], ai["nE"], ai["nN"])
    pts.update(ins.pts_update)
    anch = OutlineAnchors(
        Pi2=pts["Pi2"], Pi3=pts["Pi3"], Ti3=pts["Ti3"],
        PiX=pts["PiX"], Pi5=pts["Pi5"], TC1=pts["TC1"], R1i=ins.R1i,
    )
    geo = compute_outline_geometry(anch)
    _off = geo.origin_offset
    for _k in list(pts):
        pts[_k] = (pts[_k][0] - _off[0], pts[_k][1] - _off[1])
    _cos_r = math.cos(geo.rotation_angle)
    _sin_r = math.sin(geo.rotation_angle)
    for _k in list(pts):
        _e, _n = pts[_k]
        pts[_k] = (_e * _cos_r - _n * _sin_r, _e * _sin_r + _n * _cos_r)
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
    layout = compute_interior_layout(pts, inner_poly)
    roof = compute_roof_geometry(pts, geo.radii)
    roof_poly = roof_polyline(roof)
    return pts, geo.outline_segs, outer_poly, inner_poly, layout, roof_poly


# ── span computation ───────────────────────────────────────────

def _compute_spans(inner_poly, layout):
    """N-S span at every inch of easting.

    Returns (eastings, spans, south_spans, north_spans) all in feet.
      spans       — total N-S span (blue curve)
      south_spans — south W surface to IW midline (green curve)
      north_spans — IW midline to north W surface (cyan curve)
    When no IW is intersected, south_spans = north_spans = spans.
    """
    iw1_poly = layout.iw1                          # 4-vertex polygon
    iw8_poly = [(layout.iw8.w, layout.iw8.s),      # convert BBox → polygon
                (layout.iw8.e, layout.iw8.s),
                (layout.iw8.e, layout.iw8.n),
                (layout.iw8.w, layout.iw8.n)]

    e_min = min(p[0] for p in inner_poly)
    e_max = max(p[0] for p in inner_poly)
    inch = 1.0 / 12.0
    eastings, spans, south_spans, north_spans = [], [], [], []
    e = e_min
    while e <= e_max + 1e-9:
        ns = vert_isects(inner_poly, e)
        if len(ns) >= 2:
            south_n = min(ns)
            north_n = max(ns)
            span = north_n - south_n
        else:
            span = 0.0
            south_n = north_n = 0.0

        spans.append(span)

        # find IW midline at this easting (IW8 or IW1, whichever first from south)
        mid_n = None
        iw1_ns = vert_isects(iw1_poly, e)
        if len(iw1_ns) >= 2:
            mid_n = (min(iw1_ns) + max(iw1_ns)) / 2
        iw8_ns = vert_isects(iw8_poly, e)
        if len(iw8_ns) >= 2:
            iw8_mid = (min(iw8_ns) + max(iw8_ns)) / 2
            if mid_n is None or iw8_mid < mid_n:   # closer to south
                mid_n = iw8_mid

        if mid_n is not None and span > 0:
            south_spans.append(mid_n - south_n)
            north_spans.append(north_n - mid_n)
        else:
            south_spans.append(span)
            north_spans.append(0.0)

        eastings.append(e)
        e += inch
    return eastings, spans, south_spans, north_spans


# ── SVG generation ─────────────────────────────────────────────

def _generate_svg(pts, outer_poly, inner_poly, layout, roof_poly):
    eastings, spans, south_spans, north_spans = _compute_spans(inner_poly, layout)

    # roof spans at same eastings
    roof_spans = []
    for e in eastings:
        rns = vert_isects(roof_poly, e)
        if len(rns) >= 2:
            roof_spans.append(max(rns) - min(rns))
        else:
            roof_spans.append(0.0)
    max_span = max(spans)
    max_span_e = eastings[spans.index(max_span)]
    max_roof_span = max(roof_spans) if roof_spans else 0.0

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

    y_max_ft = math.ceil(max(max_span, max_roof_span))  # graph Y cap (feet)
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

    # grey curve — roof span
    rsp = [f"{ex(e):.1f},{sy(s):.1f}" for e, s in zip(eastings, roof_spans) if s > 0]
    if rsp:
        o.append(f'<polyline points="{" ".join(rsp)}" fill="none"'
                 f' stroke="#999" stroke-width="0.8" stroke-linejoin="round"/>')

    # green curve — south W surface to IW midline
    grn = [f"{ex(e):.1f},{sy(s):.1f}" for e, s in zip(eastings, south_spans) if s > 0]
    if grn:
        o.append(f'<polyline points="{" ".join(grn)}" fill="none"'
                 f' stroke="#2E7D32" stroke-width="0.8" stroke-linejoin="round"/>')

    # cyan curve — IW midline to north W surface
    cyn = [f"{ex(e):.1f},{sy(s):.1f}" for e, s in zip(eastings, north_spans) if s > 0]
    if cyn:
        o.append(f'<polyline points="{" ".join(cyn)}" fill="none"'
                 f' stroke="#00ACC1" stroke-width="0.8" stroke-linejoin="round"/>')

    # blue curve — total span (on top so it's visible where it overlaps)
    crv = [f"{ex(e):.1f},{sy(s):.1f}" for e, s in zip(eastings, spans) if s > 0]
    if crv:
        o.append(f'<polyline points="{" ".join(crv)}" fill="none"'
                 f' stroke="#1565C0" stroke-width="1.0" stroke-linejoin="round"/>')

    # max roof span dashed line + label (grey)
    if max_roof_span > 0:
        ry = sy(max_roof_span)
        o.append(f'<line x1="{ML}" y1="{ry:.1f}" x2="{ML + plot_w}" y2="{ry:.1f}"'
                 f' stroke="#999" stroke-width="0.7" stroke-dasharray="6,3"/>')
        o.append(f'<text x="{ML + plot_w - 3}" y="{ry + 10:.1f}" text-anchor="end"'
                 f' font-family="Arial" font-size="8" fill="#999"'
                 f' font-weight="bold">max roof: {fmt_dist(max_roof_span)}</text>')

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

    # legend
    leg_x = ML + 8
    leg_y = g_top + 10
    for i, (color, label) in enumerate([
        ("#1565C0", "Total span"),
        ("#2E7D32", "S wall \u2192 IW mid"),
        ("#00ACC1", "IW mid \u2192 N wall"),
        ("#999", "Roof span"),
    ]):
        ly_i = leg_y + i * 11
        o.append(f'<line x1="{leg_x}" y1="{ly_i}" x2="{leg_x + 14}" y2="{ly_i}"'
                 f' stroke="{color}" stroke-width="1.2"/>')
        o.append(f'<text x="{leg_x + 17}" y="{ly_i + 3}" font-family="Arial"'
                 f' font-size="6" fill="#444">{label}</text>')

    # ── outline panel ─────────────────────────────────────────
    # F polygon (outer) — gray fill shows walls
    fp = " ".join(f"{ex(p[0]):.1f},{ny(p[1]):.1f}" for p in outer_poly)
    o.append(f'<polygon points="{fp}" fill="rgba(180,180,180,0.30)"'
             f' stroke="#555" stroke-width="0.6"/>')
    # W polygon (inner) — white fill cuts out interior
    wp = " ".join(f"{ex(p[0]):.1f},{ny(p[1]):.1f}" for p in inner_poly)
    o.append(f'<polygon points="{wp}" fill="white"'
             f' stroke="#1565C0" stroke-width="0.6"/>')

    # IW1 and IW8 in the outline (show the dividing walls)
    for wall_poly in [layout.iw1,
                      [(layout.iw8.w, layout.iw8.s), (layout.iw8.e, layout.iw8.s),
                       (layout.iw8.e, layout.iw8.n), (layout.iw8.w, layout.iw8.n)]]:
        wpts = " ".join(f"{ex(p[0]):.1f},{ny(p[1]):.1f}" for p in wall_poly)
        o.append(f'<polygon points="{wpts}" fill="rgba(100,100,100,0.25)"'
                 f' stroke="#666" stroke-width="0.4"/>')

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
    pts, _, outer_poly, inner_poly, layout, roof_poly = _build_geometry()
    svg = _generate_svg(pts, outer_poly, inner_poly, layout, roof_poly)
    outpath = os.path.join(os.path.dirname(os.path.abspath(__file__)), "span.svg")
    with open(outpath, 'w', encoding='utf-8') as f:
        f.write(svg)
    print(f"Wrote {outpath}")


if __name__ == "__main__":
    main()
