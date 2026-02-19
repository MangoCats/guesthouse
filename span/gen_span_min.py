"""Generate minimum-span rotation SVG.

Finds the rotation angle that minimises the maximum interior span via a
three-pass search (5 deg coarse, 0.5 deg fine, 0.1 deg superfine).
If the result exceeds 90 deg, 180 deg is subtracted to normalise.

Outputs a single portrait-Letter SVG with:
  - Three-curve graph (blue total, green south-to-IW, cyan IW-to-north)
  - Dashed line at max total span (red)
  - Dashed line at max unsupported span (orange) = max of green / cyan
  - Rotated structure outline beneath the graph

Output: span/span_min.svg
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
    return pts, geo.outline_segs, outer_poly, inner_poly, layout


# ── IW centerlines ─────────────────────────────────────────────

def _extract_iw_centerlines(layout):
    """Midlines of IW1, IW2, IW8 as line segments [(p1, p2), ...]."""
    iw1 = layout.iw1
    mid_n1 = (layout.iw1_s + layout.iw1_n) / 2
    cl1 = ((iw1[0][0], mid_n1), ((iw1[1][0] + iw1[2][0]) / 2, mid_n1))
    iw2 = layout.iw2
    cl2 = (((iw2.w + iw2.e) / 2, iw2.s), ((iw2.w + iw2.e) / 2, iw2.n))
    iw8 = layout.iw8
    mid_n8 = (iw8.s + iw8.n) / 2
    cl8 = ((iw8.w, mid_n8), (iw8.e, mid_n8))
    return [cl1, cl2, cl8]


# ── rotation helpers ───────────────────────────────────────────

def _rot_pt(p, cx, cy, cos_a, sin_a):
    dx, dy = p[0] - cx, p[1] - cy
    return (cx + dx * cos_a - dy * sin_a, cy + dx * sin_a + dy * cos_a)


def _rot_poly(poly, cx, cy, cos_a, sin_a):
    return [_rot_pt(p, cx, cy, cos_a, sin_a) for p in poly]


def _seg_vert_isect(p1, p2, e):
    """Northing where vertical x=e crosses segment p1-p2, or None."""
    de = p2[0] - p1[0]
    if abs(de) < 1e-12:
        return None
    t = (e - p1[0]) / de
    if -1e-9 <= t <= 1 + 1e-9:
        return p1[1] + t * (p2[1] - p1[1])
    return None


# ── max-span-at-angle helper ───────────────────────────────────

def _max_span_at_angle(inner_poly, iw_cls, angle, cx, cy):
    """Return just the max total span for a given rotation angle."""
    rad = math.radians(angle)
    ca, sa = math.cos(rad), math.sin(rad)
    r_inner = _rot_poly(inner_poly, cx, cy, ca, sa)
    ie_min = min(p[0] for p in r_inner)
    ie_max = max(p[0] for p in r_inner)
    inch = 1.0 / 12.0
    best = 0.0
    e = ie_min
    while e <= ie_max + 1e-9:
        ns = vert_isects(r_inner, e)
        if len(ns) >= 2:
            span = max(ns) - min(ns)
            if span > best:
                best = span
        e += inch
    return best


# ── three-pass refinement ─────────────────────────────────────

def _find_min_span_angle(inner_poly, iw_cls, cx, cy):
    """Three-pass search returning the optimal angle (normalised)."""
    # coarse: 5 deg steps
    best_span = float('inf')
    best_angle = 0.0
    for a in range(5, 176, 5):
        ms = _max_span_at_angle(inner_poly, iw_cls, a, cx, cy)
        if ms < best_span:
            best_span = ms
            best_angle = float(a)

    # fine: 0.5 deg steps, +/-4.5 deg around coarse minimum
    fine_a = best_angle - 4.5
    while fine_a <= best_angle + 4.5 + 1e-9:
        ms = _max_span_at_angle(inner_poly, iw_cls, fine_a, cx, cy)
        if ms < best_span:
            best_span = ms
            best_angle = fine_a
        fine_a += 0.5

    # superfine: 0.1 deg steps, +/-0.9 deg around fine minimum
    center = best_angle
    sf_a = center - 0.9
    while sf_a <= center + 0.9 + 1e-9:
        ms = _max_span_at_angle(inner_poly, iw_cls, sf_a, cx, cy)
        if ms < best_span:
            best_span = ms
            best_angle = sf_a
        sf_a += 0.1

    # normalise: if > 90, subtract 180
    if best_angle > 90:
        best_angle -= 180

    return best_angle


# ── full rotation data ─────────────────────────────────────────

def _compute_rotation_data(angle, outer_poly, inner_poly, iw_cls, cx, cy):
    """Compute full span data for one rotation angle."""
    rad = math.radians(angle)
    ca, sa = math.cos(rad), math.sin(rad)
    r_outer = _rot_poly(outer_poly, cx, cy, ca, sa)
    r_inner = _rot_poly(inner_poly, cx, cy, ca, sa)
    r_cls = [(_rot_pt(c[0], cx, cy, ca, sa),
              _rot_pt(c[1], cx, cy, ca, sa)) for c in iw_cls]

    all_vis = r_outer + r_inner
    e_min = min(p[0] for p in all_vis)
    e_max = max(p[0] for p in all_vis)
    n_min = min(p[1] for p in all_vis)
    n_max = max(p[1] for p in all_vis)

    ie_min = min(p[0] for p in r_inner)
    ie_max = max(p[0] for p in r_inner)
    inch = 1.0 / 12.0
    eastings, spans, s_spans, n_spans = [], [], [], []
    e = ie_min
    while e <= ie_max + 1e-9:
        ns = vert_isects(r_inner, e)
        if len(ns) >= 2:
            bot, top = min(ns), max(ns)
            span = top - bot
        else:
            span = bot = top = 0.0
        spans.append(span)

        iw_ns = []
        for cl in r_cls:
            n = _seg_vert_isect(cl[0], cl[1], e)
            if n is not None and span > 0 and bot < n < top:
                iw_ns.append(n)
        if iw_ns and span > 0:
            s_spans.append(min(iw_ns) - bot)
            n_spans.append(top - max(iw_ns))
        else:
            s_spans.append(span)
            n_spans.append(span)

        eastings.append(e)
        e += inch

    max_span = max(spans) if spans else 0
    max_e = eastings[spans.index(max_span)] if spans else 0
    return dict(
        angle=angle, eastings=eastings, spans=spans,
        s_spans=s_spans, n_spans=n_spans,
        r_outer=r_outer, r_inner=r_inner, r_cls=r_cls,
        max_span=max_span, max_e=max_e,
        e_min=e_min, e_max=e_max, n_min=n_min, n_max=n_max,
    )


# ── SVG generation ─────────────────────────────────────────────

def _generate_svg(pts, outer_poly, inner_poly, layout):
    iw_cls = _extract_iw_centerlines(layout)
    cx = sum(p[0] for p in inner_poly) / len(inner_poly)
    cy = sum(p[1] for p in inner_poly) / len(inner_poly)

    angle = _find_min_span_angle(inner_poly, iw_cls, cx, cy)
    d = _compute_rotation_data(angle, outer_poly, inner_poly, iw_cls, cx, cy)

    max_span = d['max_span']
    max_e = d['max_e']
    max_unsup = max(max(d['s_spans']), max(d['n_spans']))
    max_unsup_e_s = d['eastings'][d['s_spans'].index(max(d['s_spans']))]
    max_unsup_e_n = d['eastings'][d['n_spans'].index(max(d['n_spans']))]

    def _fmt_angle(a):
        if a == int(a):
            return f"{int(a)}"
        s = f"{a:.1f}"
        return s.rstrip('0').rstrip('.')

    ang_str = _fmt_angle(angle)

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

    y_max_ft = math.ceil(max_span)                  # graph Y cap (feet)
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
    for ft in range(0, int(math.ceil(geo_w)) + 1, 5):
        x = ML + ft * xs
        if x > ML + plot_w + 0.5:
            break
        if 0 < ft < geo_w:
            o.append(f'<line x1="{x:.1f}" y1="{g_top}" x2="{x:.1f}" y2="{g_bot}"'
                     f' stroke="#f0f0f0" stroke-width="0.3"/>')
        o.append(f'<line x1="{x:.1f}" y1="{g_bot}" x2="{x:.1f}" y2="{g_bot + 3}"'
                 f' stroke="#666" stroke-width="0.5"/>')
        o.append(f'<text x="{x:.1f}" y="{g_bot + 12}" text-anchor="middle"'
                 f' font-family="Arial" font-size="7" fill="#555">{ft}\'</text>')

    # green curve -- south W surface to IW midline
    grn = [f"{ex(e):.1f},{sy(s):.1f}"
           for e, s in zip(d['eastings'], d['s_spans']) if s > 0]
    if grn:
        o.append(f'<polyline points="{" ".join(grn)}" fill="none"'
                 f' stroke="#2E7D32" stroke-width="0.8" stroke-linejoin="round"/>')

    # cyan curve -- IW midline to north W surface
    cyn = [f"{ex(e):.1f},{sy(s):.1f}"
           for e, s in zip(d['eastings'], d['n_spans']) if s > 0]
    if cyn:
        o.append(f'<polyline points="{" ".join(cyn)}" fill="none"'
                 f' stroke="#00ACC1" stroke-width="0.8" stroke-linejoin="round"/>')

    # blue curve -- total span (on top)
    crv = [f"{ex(e):.1f},{sy(s):.1f}"
           for e, s in zip(d['eastings'], d['spans']) if s > 0]
    if crv:
        o.append(f'<polyline points="{" ".join(crv)}" fill="none"'
                 f' stroke="#1565C0" stroke-width="1.0" stroke-linejoin="round"/>')

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
        ("#2E7D32", "Bottom \u2192 nearest IW mid"),
        ("#00ACC1", "Uppermost IW mid \u2192 top"),
        ("#C62828", "Max total span"),
        ("#E65100", "Max unsupported span"),
    ]):
        ly_i = leg_y + i * 11
        dash = "4,2" if i >= 3 else "none"
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
    pts, _, outer_poly, inner_poly, layout = _build_geometry()
    svg = _generate_svg(pts, outer_poly, inner_poly, layout)
    outpath = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           "span_min.svg")
    with open(outpath, 'w', encoding='utf-8') as f:
        f.write(svg)
    print(f"Wrote {outpath}")


if __name__ == "__main__":
    main()
