"""Generate span-vs-rotation analysis SVG.

Rotates the structure from 5 to 175 deg in 5 deg steps (35 panels).
For each rotation, measures the span from bottom to top at every inch
of the rotated x-axis.  Three curves per panel:
  blue  — total span
  green — bottom surface to nearest 6" IW centerline (IW1, IW2, IW8)
  cyan  — uppermost 6" IW centerline to top surface
When no IW centerline is intersected, green and cyan equal the full span.

The rotation with the lowest maximum total span is highlighted.
A fine search (0.5 deg steps, +/-4.5 deg around the coarse minimum) refines
the exact minimum angle.

Output: span/span_minmax.svg
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
    # IW1: polygon [SW, SE, NE, NW]; midline across the middle
    iw1 = layout.iw1
    mid_n1 = (layout.iw1_s + layout.iw1_n) / 2
    cl1 = ((iw1[0][0], mid_n1), ((iw1[1][0] + iw1[2][0]) / 2, mid_n1))
    # IW2: vertical BBox; midline is a vertical segment
    iw2 = layout.iw2
    cl2 = (((iw2.w + iw2.e) / 2, iw2.s), ((iw2.w + iw2.e) / 2, iw2.n))
    # IW8: horizontal BBox; midline is a horizontal segment
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
    """Northing where vertical x=e crosses segment p1–p2, or None."""
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


# ── SVG generation ─────────────────────────────────────────────

def _generate_svg(pts, outer_poly, inner_poly, layout):
    iw_cls = _extract_iw_centerlines(layout)
    cx = sum(p[0] for p in inner_poly) / len(inner_poly)
    cy = sum(p[1] for p in inner_poly) / len(inner_poly)
    angles = list(range(5, 176, 5))

    # pre-compute all rotations
    all_data = []
    for angle in angles:
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
        all_data.append(dict(
            angle=angle, eastings=eastings, spans=spans,
            s_spans=s_spans, n_spans=n_spans,
            r_outer=r_outer, r_inner=r_inner, r_cls=r_cls,
            max_span=max_span, max_e=max_e,
            e_min=e_min, e_max=e_max, n_min=n_min, n_max=n_max,
        ))

    # global bounds and coarse best angle
    g_e_min = min(d['e_min'] for d in all_data)
    g_e_max = max(d['e_max'] for d in all_data)
    g_max_span = max(d['max_span'] for d in all_data)
    mm_span = min(d['max_span'] for d in all_data)
    mm_angle = next(d['angle'] for d in all_data if d['max_span'] == mm_span)

    # fine search: 0.5 deg steps, +/-4.5 deg around coarse minimum
    fine_best_span = mm_span
    fine_best_angle = float(mm_angle)
    fine_a = mm_angle - 4.5
    while fine_a <= mm_angle + 4.5 + 1e-9:
        ms = _max_span_at_angle(inner_poly, iw_cls, fine_a, cx, cy)
        if ms < fine_best_span:
            fine_best_span = ms
            fine_best_angle = fine_a
        fine_a += 0.5

    # page layout (portrait Letter width)
    PW = 612
    ML, MR, MT = 62, 20, 40
    PGAP = 20
    pw = PW - ML - MR
    gw = g_e_max - g_e_min
    xs = pw / gw                                  # pt per survey-foot
    GH = 140                                      # graph height (pt)
    y_cap = math.ceil(g_max_span)                 # Y-axis cap (ft)
    ys = GH / y_cap                               # pt per span-foot

    # total SVG height
    th = MT + 16                                  # title + legend
    for d in all_data:
        oh = (d['n_max'] - d['n_min']) * xs
        th += 18 + GH + 14 + 8 + oh + PGAP
    th += 20                                      # bottom margin

    o: list[str] = []
    o.append(f'<svg xmlns="http://www.w3.org/2000/svg"'
             f' width="{PW}" height="{th:.0f}" viewBox="0 0 {PW} {th:.0f}">')
    o.append('<rect width="100%" height="100%" fill="white"/>')

    # title
    fine_ang_str = (f"{fine_best_angle:.1f}".rstrip('0').rstrip('.')
                    if fine_best_angle != int(fine_best_angle)
                    else f"{int(fine_best_angle)}")
    o.append(f'<text x="{PW / 2}" y="{MT - 12}" text-anchor="middle"'
             f' font-family="Arial" font-size="13" font-weight="bold"'
             f' fill="#222">Span vs. Rotation (5\u00b0\u2013175\u00b0)'
             f' \u2014 min max: {fmt_dist(fine_best_span)}'
             f' at {fine_ang_str}\u00b0</text>')

    # legend
    ly = MT + 2
    for i, (col, lab) in enumerate([
        ("#1565C0", "Total span"),
        ("#2E7D32", "Bottom \u2192 nearest IW mid"),
        ("#00ACC1", "Uppermost IW mid \u2192 top"),
    ]):
        lx = ML + i * 135
        o.append(f'<line x1="{lx}" y1="{ly}" x2="{lx + 14}" y2="{ly}"'
                 f' stroke="{col}" stroke-width="1.2"/>')
        o.append(f'<text x="{lx + 17}" y="{ly + 3}" font-family="Arial"'
                 f' font-size="6.5" fill="#444">{lab}</text>')

    yc = MT + 16

    def ex(e):
        return ML + (e - g_e_min) * xs

    for d in all_data:
        ang = d['angle']
        is_min = ang == mm_angle
        oh = (d['n_max'] - d['n_min']) * xs

        # panel title
        yc += 14
        tc = "#C62828" if is_min else "#333"
        tw = "bold" if is_min else "normal"
        if is_min:
            tag = (f" \u2190 MIN (refined: {fmt_dist(fine_best_span)}"
                   f" at {fine_ang_str}\u00b0)")
        else:
            tag = ""
        o.append(f'<text x="{ML}" y="{yc}" font-family="Arial"'
                 f' font-size="9" fill="{tc}" font-weight="{tw}">'
                 f'{ang}\u00b0  max: {fmt_dist(d["max_span"])}{tag}</text>')
        yc += 4

        gt = yc                                   # graph top
        gb = gt + GH                              # graph bottom

        # highlight min panel
        if is_min:
            ph = GH + 14 + 8 + oh
            o.append(f'<rect x="{ML - 4}" y="{gt - 2}" width="{pw + 8}"'
                     f' height="{ph + 4:.0f}" fill="rgba(198,40,40,0.04)"'
                     f' stroke="#C62828" stroke-width="0.5" rx="3"/>')

        # graph frame
        o.append(f'<rect x="{ML}" y="{gt}" width="{pw}" height="{GH}"'
                 f' fill="none" stroke="#bbb" stroke-width="0.5"/>')

        # Y grid + labels
        ystep = 2 if y_cap > 14 else 1
        for ft in range(0, y_cap + 1, ystep):
            yy = gb - ft * ys
            if ft > 0:
                o.append(f'<line x1="{ML}" y1="{yy:.1f}" x2="{ML + pw}" y2="{yy:.1f}"'
                         f' stroke="#eaeaea" stroke-width="0.3"/>')
            o.append(f'<line x1="{ML - 3}" y1="{yy:.1f}" x2="{ML}" y2="{yy:.1f}"'
                     f' stroke="#666" stroke-width="0.5"/>')
            o.append(f'<text x="{ML - 5}" y="{yy + 2.5:.1f}" text-anchor="end"'
                     f' font-family="Arial" font-size="6" fill="#555">{ft}\'</text>')

        # X ticks (every 5 ft)
        for ft in range(0, int(math.ceil(gw)) + 1, 5):
            xx = ML + ft * xs
            if xx > ML + pw + 0.5:
                break
            if 0 < ft < gw:
                o.append(f'<line x1="{xx:.1f}" y1="{gt}" x2="{xx:.1f}" y2="{gb}"'
                         f' stroke="#f0f0f0" stroke-width="0.3"/>')
            o.append(f'<line x1="{xx:.1f}" y1="{gb}" x2="{xx:.1f}" y2="{gb + 3}"'
                     f' stroke="#666" stroke-width="0.5"/>')
            o.append(f'<text x="{xx:.1f}" y="{gb + 10}" text-anchor="middle"'
                     f' font-family="Arial" font-size="5" fill="#555">{ft}\'</text>')

        # green curve — bottom to nearest IW centerline
        grn = [f"{ex(e):.1f},{gb - s * ys:.1f}"
               for e, s in zip(d['eastings'], d['s_spans']) if s > 0]
        if grn:
            o.append(f'<polyline points="{" ".join(grn)}" fill="none"'
                     f' stroke="#2E7D32" stroke-width="0.6" stroke-linejoin="round"/>')

        # cyan curve — uppermost IW centerline to top
        cyn = [f"{ex(e):.1f},{gb - s * ys:.1f}"
               for e, s in zip(d['eastings'], d['n_spans']) if s > 0]
        if cyn:
            o.append(f'<polyline points="{" ".join(cyn)}" fill="none"'
                     f' stroke="#00ACC1" stroke-width="0.6" stroke-linejoin="round"/>')

        # blue curve — total span (on top)
        crv = [f"{ex(e):.1f},{gb - s * ys:.1f}"
               for e, s in zip(d['eastings'], d['spans']) if s > 0]
        if crv:
            o.append(f'<polyline points="{" ".join(crv)}" fill="none"'
                     f' stroke="#1565C0" stroke-width="0.8" stroke-linejoin="round"/>')

        # max-span dashed line
        my = gb - d['max_span'] * ys
        o.append(f'<line x1="{ML}" y1="{my:.1f}" x2="{ML + pw}" y2="{my:.1f}"'
                 f' stroke="#C62828" stroke-width="0.5" stroke-dasharray="4,2"/>')

        yc = gb + 14                              # past x-axis labels

        # ── outline panel ─────────────────────────────────────
        ot = yc + 8                               # outline top
        ob = ot + oh                              # outline bottom
        nm = d['n_min']

        # outer polygon
        fp = " ".join(f"{ex(p[0]):.1f},{ob - (p[1] - nm) * xs:.1f}"
                      for p in d['r_outer'])
        o.append(f'<polygon points="{fp}" fill="rgba(180,180,180,0.25)"'
                 f' stroke="#555" stroke-width="0.4"/>')

        # inner polygon
        wp = " ".join(f"{ex(p[0]):.1f},{ob - (p[1] - nm) * xs:.1f}"
                      for p in d['r_inner'])
        o.append(f'<polygon points="{wp}" fill="white"'
                 f' stroke="#1565C0" stroke-width="0.4"/>')

        # IW centerlines
        for cl in d['r_cls']:
            x1, y1 = ex(cl[0][0]), ob - (cl[0][1] - nm) * xs
            x2, y2 = ex(cl[1][0]), ob - (cl[1][1] - nm) * xs
            o.append(f'<line x1="{x1:.1f}" y1="{y1:.1f}"'
                     f' x2="{x2:.1f}" y2="{y2:.1f}"'
                     f' stroke="#666" stroke-width="0.5" stroke-dasharray="2,1"/>')

        # max-span easting reference line
        mx = ex(d['max_e'])
        o.append(f'<line x1="{mx:.1f}" y1="{gt}" x2="{mx:.1f}" y2="{ob}"'
                 f' stroke="#C62828" stroke-width="0.3"'
                 f' stroke-dasharray="2,3" opacity="0.3"/>')

        yc = ob + PGAP

    # version stamp
    ver = git_describe()
    o.append(f'<text x="{PW - 5}" y="{th - 5}" text-anchor="end"'
             f' font-family="Arial" font-size="5" fill="#ccc">{ver}</text>')

    o.append('</svg>')
    return '\n'.join(o)


# ── entry point ────────────────────────────────────────────────

def main():
    pts, _, outer_poly, inner_poly, layout = _build_geometry()
    svg = _generate_svg(pts, outer_poly, inner_poly, layout)
    outpath = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           "span_minmax.svg")
    with open(outpath, 'w', encoding='utf-8') as f:
        f.write(svg)
    print(f"Wrote {outpath}")


if __name__ == "__main__":
    main()
