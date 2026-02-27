"""Common geometry bootstrap and rotation helpers for span generators."""
import math

from shared.geometry import (
    path_polygon, vert_isects, compute_inner_walls,
    f8f9_corner_polyline, GEOM_EPS, require_pts,
)
from shared.survey import compute_traverse, compute_three_arc, compute_inset
from floorplan.geometry import compute_outline_geometry, align_pts_to_f_series
from floorplan.constants import WALL_OUTER, F8F9_INNER_TURN_R
from floorplan.layout import compute_interior_layout
from floorplan.roof import compute_roof_geometry, roof_polyline


# ── geometry bootstrap ─────────────────────────────────────────

def build_geometry():
    """Return (pts, outline_segs, outer_poly, inner_poly, layout, roof_poly)."""
    pts = compute_traverse()
    ai = compute_three_arc(pts)
    ins = compute_inset(pts, ai["R1"], ai["R2"], ai["R3"], ai["nE"], ai["nN"])
    pts.update(ins.pts_update)
    align_pts_to_f_series(pts)
    geo = compute_outline_geometry()
    pts.update(geo.fp_pts)
    inner_segs = compute_inner_walls(geo.outline_segs, pts, WALL_OUTER, geo.radii)
    outer_poly = path_polygon(geo.outline_segs, pts)
    inner_poly = path_polygon(inner_segs, pts)
    # patch concave W8-W9 corner
    require_pts(pts, "W8", "W9")
    w_f8f9 = f8f9_corner_polyline(pts, WALL_OUTER, F8F9_INNER_TURN_R)
    w8, w9 = pts["W8"], pts["W9"]
    i8 = next(i for i, p in enumerate(inner_poly)
              if abs(p[0] - w8[0]) < GEOM_EPS and abs(p[1] - w8[1]) < GEOM_EPS)
    i9 = next(i for i, p in enumerate(inner_poly)
              if i > i8 and abs(p[0] - w9[0]) < GEOM_EPS and abs(p[1] - w9[1]) < GEOM_EPS)
    inner_poly[i8:i9 + 1] = w_f8f9
    layout = compute_interior_layout(pts, inner_poly)
    roof = compute_roof_geometry(pts, geo.radii)
    roof_poly = roof_polyline(roof)
    return pts, geo.outline_segs, outer_poly, inner_poly, layout, roof_poly


# ── IW centerlines ─────────────────────────────────────────────

def extract_iw_centerlines(layout):
    """Midlines of IW1, IW2, IW8 as line segments [(p1, p2), ...]."""
    iw1 = layout.iw1
    mid_n1 = (iw1.s + iw1.n) / 2
    cl1 = ((iw1.poly[0][0], mid_n1), ((iw1.poly[1][0] + iw1.poly[2][0]) / 2, mid_n1))
    iw2 = layout.iw2
    cl2 = (((iw2.w + iw2.e) / 2, iw2.s), ((iw2.w + iw2.e) / 2, iw2.n))
    iw8 = layout.iw8
    mid_n8 = (iw8.s + iw8.n) / 2
    cl8 = ((iw8.w, mid_n8), (iw8.e, mid_n8))
    return [cl1, cl2, cl8]


# ── rotation helpers ───────────────────────────────────────────

def rot_pt(p, cx, cy, cos_a, sin_a):
    """Rotate point p around (cx, cy) by pre-computed cos/sin."""
    dx, dy = p[0] - cx, p[1] - cy
    return (cx + dx * cos_a - dy * sin_a, cy + dx * sin_a + dy * cos_a)


def rot_poly(poly, cx, cy, cos_a, sin_a):
    """Rotate all points in poly around (cx, cy)."""
    return [rot_pt(p, cx, cy, cos_a, sin_a) for p in poly]


def seg_vert_isect(p1, p2, e):
    """Northing where vertical x=e crosses segment p1-p2, or None."""
    de = p2[0] - p1[0]
    if abs(de) < 1e-12:
        return None
    t = (e - p1[0]) / de
    if -GEOM_EPS <= t <= 1 + GEOM_EPS:
        return p1[1] + t * (p2[1] - p1[1])
    return None


# ── max-span-at-angle helper ───────────────────────────────────

def max_span_at_angle(inner_poly, iw_cls, angle, cx, cy):
    """Return just the max total span for a given rotation angle."""
    rad = math.radians(angle)
    ca, sa = math.cos(rad), math.sin(rad)
    r_inner = rot_poly(inner_poly, cx, cy, ca, sa)
    ie_min = min(p[0] for p in r_inner)
    ie_max = max(p[0] for p in r_inner)
    inch = 1.0 / 12.0
    n_steps = int((ie_max - ie_min) / inch) + 1
    best = 0.0
    for i in range(n_steps):
        e = ie_min + i * inch
        ns = vert_isects(r_inner, e)
        if len(ns) >= 2:
            span = max(ns) - min(ns)
            if span > best:
                best = span
    return best


# ── three-pass refinement ─────────────────────────────────────

def find_min_span_angle(inner_poly, iw_cls, cx, cy, normalize=True):
    """Three-pass search returning (angle, span) for the minimum-span rotation.

    If normalize is True (default), angles > 90 are shifted by -180.
    """
    # coarse: 5 deg steps
    best_span = float('inf')
    best_angle = 0.0
    for a in range(0, 176, 5):
        ms = max_span_at_angle(inner_poly, iw_cls, a, cx, cy)
        if ms < best_span:
            best_span = ms
            best_angle = float(a)

    # fine: 0.5 deg steps, +/-4.5 deg around coarse minimum
    fine_base = best_angle - 4.5
    for i in range(19):  # 0..18 → -4.5 to +4.5 in 0.5 steps
        a = fine_base + i * 0.5
        ms = max_span_at_angle(inner_poly, iw_cls, a, cx, cy)
        if ms < best_span:
            best_span = ms
            best_angle = a

    # superfine: 0.1 deg steps, +/-0.9 deg around fine minimum
    sf_base = best_angle - 0.9
    for i in range(19):  # 0..18 → -0.9 to +0.9 in 0.1 steps
        a = sf_base + i * 0.1
        ms = max_span_at_angle(inner_poly, iw_cls, a, cx, cy)
        if ms < best_span:
            best_span = ms
            best_angle = a

    if normalize and best_angle > 90:
        best_angle -= 180

    return best_angle, best_span


# ── full rotation data ─────────────────────────────────────────

def compute_rotation_data(angle, outer_poly, inner_poly, iw_cls, cx, cy,
                          roof_poly=None):
    """Compute full span data for one rotation angle."""
    rad = math.radians(angle)
    ca, sa = math.cos(rad), math.sin(rad)
    r_outer = rot_poly(outer_poly, cx, cy, ca, sa)
    r_inner = rot_poly(inner_poly, cx, cy, ca, sa)
    r_cls = [(rot_pt(c[0], cx, cy, ca, sa),
              rot_pt(c[1], cx, cy, ca, sa)) for c in iw_cls]
    r_roof = rot_poly(roof_poly, cx, cy, ca, sa) if roof_poly else None

    all_vis = r_outer + r_inner
    e_min = min(p[0] for p in all_vis)
    e_max = max(p[0] for p in all_vis)
    n_min = min(p[1] for p in all_vis)
    n_max = max(p[1] for p in all_vis)

    ie_min = min(p[0] for p in r_inner)
    ie_max = max(p[0] for p in r_inner)
    inch = 1.0 / 12.0
    n_steps = int((ie_max - ie_min) / inch) + 1
    eastings, spans, s_spans, n_spans, roof_spans = [], [], [], [], []
    for i in range(n_steps):
        e = ie_min + i * inch
        ns = vert_isects(r_inner, e)
        if len(ns) >= 2:
            bot, top = min(ns), max(ns)
            span = top - bot
        else:
            span = bot = top = 0.0
        spans.append(span)

        iw_ns = []
        for cl in r_cls:
            n = seg_vert_isect(cl[0], cl[1], e)
            if n is not None and span > 0 and bot < n < top:
                iw_ns.append(n)
        if iw_ns and span > 0:
            s_spans.append(min(iw_ns) - bot)
            n_spans.append(top - max(iw_ns))
        else:
            s_spans.append(span)
            n_spans.append(span)

        if r_roof:
            rns = vert_isects(r_roof, e)
            roof_spans.append(max(rns) - min(rns) if len(rns) >= 2 else 0.0)

        eastings.append(e)

    max_span = max(spans) if spans else 0
    max_e = eastings[spans.index(max_span)] if spans else 0
    max_roof_span = max(roof_spans) if roof_spans else 0
    return dict(
        angle=angle, eastings=eastings, spans=spans,
        s_spans=s_spans, n_spans=n_spans, roof_spans=roof_spans,
        r_outer=r_outer, r_inner=r_inner, r_cls=r_cls,
        max_span=max_span, max_e=max_e, max_roof_span=max_roof_span,
        e_min=e_min, e_max=e_max, n_min=n_min, n_max=n_max,
    )


# ── formatting ─────────────────────────────────────────────────

def fmt_angle(a):
    """Format angle for display: integer if whole, else 1 decimal trimmed."""
    if a == int(a):
        return f"{int(a)}"
    s = f"{a:.1f}"
    return s.rstrip('0').rstrip('.')


# ── shared SVG chart helpers ──────────────────────────────────

# Span curve color palette (order: roof, south-to-IW, IW-to-north, total)
CURVE_COLORS = ["#999", "#2E7D32", "#00ACC1", "#1565C0"]


def render_y_grid(o, ML, plot_w, g_bot, y_max_ft, ys, font_size="7"):
    """Render Y-axis grid lines, tick marks, and labels."""
    y_step = 2 if y_max_ft > 14 else 1
    for ft in range(0, y_max_ft + 1, y_step):
        y = g_bot - ft * ys
        if ft > 0:
            o.append(f'<line x1="{ML}" y1="{y:.1f}" x2="{ML + plot_w}" y2="{y:.1f}"'
                     f' stroke="#eaeaea" stroke-width="0.3"/>')
        o.append(f'<line x1="{ML - 3}" y1="{y:.1f}" x2="{ML}" y2="{y:.1f}"'
                 f' stroke="#666" stroke-width="0.5"/>')
        o.append(f'<text x="{ML - 5}" y="{y + 2.5:.1f}" text-anchor="end"'
                 f' font-family="Arial" font-size="{font_size}" fill="#555">{ft}\'</text>')


def render_x_grid(o, ML, plot_w, g_top, g_bot, xs, axis_ft,
                  font_size="7", label_y_offset=12):
    """Render X-axis grid lines, tick marks, and labels."""
    for ft in range(0, int(math.ceil(axis_ft)) + 1, 5):
        x = ML + ft * xs
        if x > ML + plot_w + 0.5:
            break
        if 0 < ft < axis_ft:
            o.append(f'<line x1="{x:.1f}" y1="{g_top}" x2="{x:.1f}" y2="{g_bot}"'
                     f' stroke="#f0f0f0" stroke-width="0.3"/>')
        o.append(f'<line x1="{x:.1f}" y1="{g_bot}" x2="{x:.1f}" y2="{g_bot + 3}"'
                 f' stroke="#666" stroke-width="0.5"/>')
        o.append(f'<text x="{x:.1f}" y="{g_bot + label_y_offset}" text-anchor="middle"'
                 f' font-family="Arial" font-size="{font_size}" fill="#555">{ft}\'</text>')


def render_span_curves(o, eastings, data_arrays, ex_fn, sy_fn,
                       stroke_widths=(0.8, 0.8, 0.8, 1.0)):
    """Render the 4 standard span polyline curves.

    data_arrays: (roof_spans, s_spans, n_spans, spans) — in rendering order.
    stroke_widths: per-curve stroke widths (default: 0.8 for first 3, 1.0 for total).
    """
    for arr, color, sw in zip(data_arrays, CURVE_COLORS, stroke_widths):
        if not arr:
            continue
        pts_str = [f"{ex_fn(e):.1f},{sy_fn(s):.1f}"
                   for e, s in zip(eastings, arr) if s > 0]
        if pts_str:
            o.append(f'<polyline points="{" ".join(pts_str)}" fill="none"'
                     f' stroke="{color}" stroke-width="{sw}"'
                     f' stroke-linejoin="round"/>')
