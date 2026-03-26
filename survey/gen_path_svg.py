import os, sys, math, datetime
from typing import NamedTuple

# Ensure project root is on sys.path for package imports
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

from shared.types import Point, LineSeg, ArcSeg, Segment
from shared.geometry import (
    poly_area, segment_polyline, path_polygon, arc_sweep_deg,
    brg_dist, fmt_brg, fmt_dist,
    compute_inner_walls,
)
from shared.survey import compute_traverse, compute_three_arc, compute_inset, compute_pt1
from shared.svg import make_svg_transform, W, H, git_describe, normalize_svg_angle, svg_polygon_pts
from floorplan.geometry import compute_outline_geometry, align_pts_to_f_series
from floorplan.layout import compute_interior_layout
from floorplan.constants import WALL_OUTER

# ============================================================
# Section 1: Rendering Types
# ============================================================
class VertexStyle(NamedTuple):
    display_name: str; anchor: str; dx: float; dy: float
    color: str; dot_radius: float; font_size: float

class BrgDistLabel(NamedTuple):
    offset: float

class ArcLabel(NamedTuple):
    text1: str; text2: str | None; anchor: str
    dx: float; dy: float; dy2: float; color: str

class CenterMark(NamedTuple):
    center: str; tangent_to: str; color: str

class LayerConfig(NamedTuple):
    opacity: float; fill_color: str
    line_stroke: str; line_width: float
    arc_styles: dict
    vertex_styles: dict
    brg_dist_labels: dict | None
    arc_labels: dict | None
    center_marks: list | None
    traverse_pts: list | None
    traverse_stroke: str | None
    brg_decimal: bool = False

# W, H imported from shared.svg

# ============================================================
# Section 4: Generic SVG Renderer
# ============================================================
def render_layer(lines: list, segments: list[Segment], pts: dict, cfg: LayerConfig, to_svg=None):
    if cfg.opacity < 1.0:
        lines.append(f'<g opacity="{cfg.opacity}">')

    # Dashed traverse overlay
    if cfg.traverse_pts:
        trav = svg_polygon_pts([pts[n] for n in cfg.traverse_pts], to_svg)
        lines.append(f'<polygon points="{trav}" fill="none" stroke="{cfg.traverse_stroke}"'
                     f' stroke-width="0.8" stroke-dasharray="4,4"/>')

    # Filled polygon
    dense = path_polygon(segments, pts)
    poly_svg = svg_polygon_pts(dense, to_svg)
    lines.append(f'<polygon points="{poly_svg}" fill="{cfg.fill_color}" stroke="none"/>')

    # Straight segments
    for seg in segments:
        if isinstance(seg, LineSeg):
            sx1, sy1 = to_svg(*pts[seg.start]); sx2, sy2 = to_svg(*pts[seg.end])
            lines.append(f'<line x1="{sx1:.1f}" y1="{sy1:.1f}" x2="{sx2:.1f}" y2="{sy2:.1f}"'
                         f' stroke="{cfg.line_stroke}" stroke-width="{cfg.line_width}"/>')

    # Arc polylines — emit in arc_styles order (not path order)
    arc_seg_map = {(s.start, s.end): s for s in segments if isinstance(s, ArcSeg)}
    for key, (stroke, width) in cfg.arc_styles.items():
        seg = arc_seg_map[key]
        poly = segment_polyline(seg, pts)
        svg_pts = svg_polygon_pts(poly, to_svg)
        lines.append(f'<polyline points="{svg_pts}" fill="none" stroke="{stroke}"'
                     f' stroke-width="{width}" stroke-linecap="round"/>')

    # Center marks and dashed radius lines
    if cfg.center_marks:
        lines.append('<g clip-path="url(#page)">')
        for cm in cfg.center_marks:
            cx_s, cy_s = to_svg(*pts[cm.center]); tx, ty = to_svg(*pts[cm.tangent_to])
            lines.append(f'<line x1="{cx_s:.1f}" y1="{cy_s:.1f}" x2="{tx:.1f}" y2="{ty:.1f}"'
                         f' stroke="{cm.color}" stroke-width="0.6" stroke-dasharray="4,3"/>')
            if 0 < cx_s < W and 0 < cy_s < H:
                lines.append(f'<line x1="{cx_s-4:.1f}" y1="{cy_s:.1f}" x2="{cx_s+4:.1f}" y2="{cy_s:.1f}"'
                             f' stroke="{cm.color}" stroke-width="0.7"/>')
                lines.append(f'<line x1="{cx_s:.1f}" y1="{cy_s-4:.1f}" x2="{cx_s:.1f}" y2="{cy_s+4:.1f}"'
                             f' stroke="{cm.color}" stroke-width="0.7"/>')
        lines.append('</g>')

    # Vertex dots and labels
    seen = set(); vertex_names = []
    for seg in segments:
        if seg.start not in seen: vertex_names.append(seg.start); seen.add(seg.start)
        if seg.end not in seen: vertex_names.append(seg.end); seen.add(seg.end)
    # Include styled vertices not on any segment (e.g. FC)
    for vname in cfg.vertex_styles:
        if vname not in seen and vname in pts:
            vertex_names.append(vname)
    for vname in vertex_names:
        if vname not in cfg.vertex_styles: continue
        vs = cfg.vertex_styles[vname]
        sx, sy = to_svg(*pts[vname])
        lines.append(f'<circle cx="{sx:.1f}" cy="{sy:.1f}" r="{vs.dot_radius}" fill="{vs.color}"/>')
        lines.append(f'<text x="{sx+vs.dx:.1f}" y="{sy+vs.dy:.1f}" text-anchor="{vs.anchor}"'
                     f' font-family="Arial" font-size="{vs.font_size}" font-weight="bold"'
                     f' fill="{vs.color}">{vs.display_name}</text>')

    # Bearing/distance labels
    if cfg.brg_dist_labels:
        for seg in segments:
            if not isinstance(seg, LineSeg): continue
            key = (seg.start, seg.end)
            if key not in cfg.brg_dist_labels: continue
            bdl = cfg.brg_dist_labels[key]
            b, d = brg_dist(pts[seg.start], pts[seg.end])
            sx1, sy1 = to_svg(*pts[seg.start]); sx2, sy2 = to_svg(*pts[seg.end])
            mx, my = (sx1+sx2)/2, (sy1+sy2)/2
            ang = normalize_svg_angle(math.degrees(math.atan2(-(sy2-sy1), sx2-sx1)))
            ddx = sx2-sx1; ddy = sy2-sy1; ll = math.sqrt(ddx**2+ddy**2)
            if ll > 0: nx, ny = -ddy/ll, ddx/ll; mx += nx*bdl.offset; my += ny*bdl.offset
            lines.append(f'<g transform="translate({mx:.1f},{my:.1f}) rotate({ang:.1f})">')
            brg_text = f"{b:.2f}\u00b0" if cfg.brg_decimal else f"Brg {fmt_brg(b)}"
            dist_text = fmt_dist(d) if cfg.brg_decimal else f"Dist {fmt_dist(d)}"
            lines.append(f'  <text x="0" y="-3" text-anchor="middle" font-family="Arial"'
                         f' font-size="8.5" fill="#1a237e">{brg_text}</text>')
            lines.append(f'  <text x="0" y="8" text-anchor="middle" font-family="Arial"'
                         f' font-size="8.5" fill="#555">{dist_text}</text>')
            lines.append('</g>')

    # Arc labels — emit in arc_labels order (not path order)
    if cfg.arc_labels:
        for key, al in cfg.arc_labels.items():
            seg = arc_seg_map[key]
            poly = segment_polyline(seg, pts)
            mid = poly[len(poly)//2]
            sx, sy = to_svg(*mid)
            lines.append(f'<text x="{sx+al.dx:.1f}" y="{sy+al.dy:.1f}" text-anchor="{al.anchor}"'
                         f' font-family="Arial" font-size="8.5" fill="{al.color}">{al.text1}</text>')
            if al.text2 is not None:
                lines.append(f'<text x="{sx+al.dx:.1f}" y="{sy+al.dy+al.dy2:.1f}" text-anchor="{al.anchor}"'
                             f' font-family="Arial" font-size="8.5" fill="{al.color}">{al.text2}</text>')

    if cfg.opacity < 1.0:
        lines.append('</g>')

# ============================================================
# Section 5: Geometry Computation
# ============================================================
def compute_all():
    """Compute all geometry. Returns dict with everything needed for rendering."""
    pts = compute_traverse()
    to_svg = make_svg_transform()
    arc_info = compute_three_arc(pts)
    R1, R2, R3 = arc_info["R1"], arc_info["R2"], arc_info["R3"]
    nE, nN = arc_info["nE"], arc_info["nN"]

    # PT1: tangency point where TC1 arc meets the P4-P5 line extension
    pts["PT1"] = compute_pt1(pts, R1)

    # Outer path (corner at PT1 replaces P5; PT1→T1 is CW arc on TC1 circle)
    outer_segs = [
        LineSeg("POB", "P2"), LineSeg("P2", "P3"), LineSeg("P3", "T3"),
        ArcSeg("T3", "PX", "TC3", R3, "CW", 60),
        LineSeg("PX", "P4"), LineSeg("P4", "P5"), LineSeg("P5", "PT1"),
        ArcSeg("PT1", "T1", "TC1", R1, "CW", 60),
        ArcSeg("T1", "PA", "TC1", R1, "CW", 60),
        ArcSeg("PA", "T2", "TC2", R2, "CW", 60),
        LineSeg("T2", "POB"),
    ]
    outer_area = poly_area(path_polygon(outer_segs, pts))

    # Inset path (6" inside)
    inset = compute_inset(pts, R1, R2, R3, nE, nN)
    pts.update(inset.pts_update)
    inset_segs = inset.inset_segs
    inset_area = poly_area(path_polygon(inset_segs, pts))

    # F-series outline (chain walk — independent of survey data)
    outline_geo = compute_outline_geometry()

    # Align P/Pi with F-series: rigid transform (rotate + translate)
    align_pts_to_f_series(pts)
    pts.update(outline_geo.fp_pts)
    outline_segs = outline_geo.outline_segs
    radii = outline_geo.radii

    # Derive U-series as aliases (downstream from F-series)
    for i in [j for j in range(19) if j not in (0, 3, 4)]:
        pts[f"U{i}"] = outline_geo.fp_pts[f"F{i}"]

    # Inner walls + layout
    inner_segs = compute_inner_walls(outline_segs, pts, WALL_OUTER, radii)
    outer_poly = path_polygon(outline_segs, pts)
    inner_poly = path_polygon(inner_segs, pts)
    outline_area = poly_area(outer_poly)
    layout = compute_interior_layout(pts, inner_poly)

    return {
        "pts": pts, "to_svg": to_svg,
        "outer_segs": outer_segs, "inset_segs": inset_segs, "outline_segs": outline_segs,
        "outer_area": outer_area, "inset_area": inset_area, "outline_area": outline_area,
        "R1i": inset.R1i, "R2i": inset.R2i, "R3i": inset.R3i,
        "radii": radii, "outer_poly": outer_poly, "inner_poly": inner_poly,
        "inner_segs": inner_segs, "layout": layout,
    }

def render_floorplan(lines, to_svg, pts, outer_poly, inner_poly, inner_segs,
                     layout=None, iw_polys=None):
    lines.append('<g opacity="0.5">')
    # Wall band (outer - inner hole, evenodd fill)
    od = "M "+" L ".join(f"{to_svg(*p)[0]:.1f},{to_svg(*p)[1]:.1f}" for p in outer_poly)+" Z"
    id_ = "M "+" L ".join(f"{to_svg(*p)[0]:.1f},{to_svg(*p)[1]:.1f}" for p in reversed(inner_poly))+" Z"
    lines.append(f'<path d="{od} {id_}" fill="rgba(160,160,160,0.5)" fill-rule="evenodd" stroke="none"/>')
    # Inner wall strokes
    for seg in inner_segs:
        if isinstance(seg, LineSeg):
            s1,s2 = to_svg(*pts[seg.start]),to_svg(*pts[seg.end])
            lines.append(f'<line x1="{s1[0]:.1f}" y1="{s1[1]:.1f}" x2="{s2[0]:.1f}" y2="{s2[1]:.1f}" stroke="#666" stroke-width="1.0"/>')
        else:
            pl = segment_polyline(seg, pts)
            lines.append(f'<polyline points="{svg_polygon_pts(pl, to_svg)}" fill="none" stroke="#666" stroke-width="1.0" stroke-linecap="round"/>')
    # Interior wall polygons — DB-sourced if provided, else fall back to seed layout
    if iw_polys is not None:
        _polys = list(iw_polys.values())
    elif layout is not None:
        L = layout
        _polys = [iw.poly for iw in [L.iw1, L.iw2, L.iw2s, L.iw2o, L.iw3, L.iw4, L.iw5, L.iw6,
                                     L.iw7, L.iw8, L.iw9, L.iw11, L.iw12]]
    else:
        _polys = []
    for poly in _polys:
        lines.append(f'<polygon points="{svg_polygon_pts(poly, to_svg)}" fill="rgba(160,160,160,0.5)" stroke="#666" stroke-width="0.8"/>')
    lines.append('</g>')

# ============================================================
# Section 8: Style Configurations (static — no computed values)
# ============================================================
outer_cfg = LayerConfig(
    opacity=0.2, fill_color="#e8edf5",
    line_stroke="#333", line_width=2.0,
    arc_styles={
        ("PT1","T1"): ("#0077B6", 2.5), ("T1","PA"): ("#0077B6", 2.5),
        ("PA","T2"): ("#0077B6", 2.5),
        ("T3","PX"): ("#2E7D32", 2.5),
    },
    vertex_styles={
        "POB": VertexStyle("P.O.B.", "middle", 0, -16, "#d32f2f", 3.5, 10),
        "P2":  VertexStyle("P2",  "end",  -10, -8, "#d32f2f", 3.5, 10),
        "P3":  VertexStyle("P3",  "end",  -10,  7, "#d32f2f", 3.5, 10),
        "T3":  VertexStyle("T3",  "middle", 0, -12, "#2E7D32", 3.5, 10),
        "PX":  VertexStyle("PX",  "start",  8, -6, "#2E7D32", 3.5, 10),
        "P4":  VertexStyle("P4",  "start", 10, 10, "#d32f2f", 3.5, 10),
        "P5":  VertexStyle("P5",  "start", 10,  4, "#d32f2f", 3.5, 10),
        "PT1": VertexStyle("PT1", "start", 10,  3, "#0077B6", 3.5, 10),
        "T1":  VertexStyle("T1",  "start", 10, -4, "#0077B6", 3.5, 10),
        "PA":  VertexStyle("PA",  "start",  8, -8, "#0077B6", 3.5, 10),
        "T2":  VertexStyle("T2",  "end",  -10, -4, "#0077B6", 3.5, 10),
    },
    brg_dist_labels={
        ("POB","P2"): BrgDistLabel(-18), ("P2","P3"): BrgDistLabel(-18),
        ("P3","T3"): BrgDistLabel(18), ("PX","P4"): BrgDistLabel(-16),
        ("P4","P5"): BrgDistLabel(-16), ("P5","PT1"): BrgDistLabel(-16),
        ("T2","POB"): BrgDistLabel(16),
    },
    arc_labels={
        ("PT1","T1"): ArcLabel("Arc 1: R=10'", None, "start", 12, 4, 0, "#0077B6"),
        ("PA","T2"): ArcLabel("Arc 2: R=12' 6\"", None, "middle", 0, 16, 0, "#0077B6"),
        ("T3","PX"): ArcLabel("Arc 3: R=11'", None, "start", 12, 4, 0, "#2E7D32"),
    },
    center_marks=None,
    traverse_pts=["POB","P2","P3","P4","P5"],
    traverse_stroke="#bbb",
)

inset_cfg = LayerConfig(
    opacity=0.2, fill_color="rgba(255,152,0,0.3)",
    line_stroke="#BF360C", line_width=1.5,
    arc_styles={
        ("PTi1","Ti1"): ("#BF360C", 1.5),
        ("Ti1","Ai2"): ("#BF360C", 1.5), ("Ai2","Ti2"): ("#BF360C", 1.5),
        ("Ti3","PiX"): ("#BF360C", 1.5),
    },
    vertex_styles={
        "PiOB": VertexStyle("PiOB", "middle", 0, 14, "#BF360C", 2.5, 8.5),
        "Pi2":  VertexStyle("Pi2",  "start",  6,  4, "#BF360C", 2.5, 8.5),
        "Pi3":  VertexStyle("Pi3",  "start",  6, -4, "#BF360C", 2.5, 8.5),
        "Ti3":  VertexStyle("Ti3",  "middle", 0, 14, "#BF360C", 2.5, 8.5),
        "PiX":  VertexStyle("PiX",  "end",   -8, -6, "#BF360C", 2.5, 8.5),
        "Pi4":  VertexStyle("Pi4",  "end",   -8, 12, "#BF360C", 2.5, 8.5),
        "Pi5":  VertexStyle("Pi5",  "end",   -8,  4, "#BF360C", 2.5, 8.5),
        "PTi1": VertexStyle("PTi1", "end",  -12, -8, "#BF360C", 2.5, 8.5),
        "Ti1":  VertexStyle("Ti1",  "end",   -8, 12, "#BF360C", 2.5, 8.5),
        "Ai2":  VertexStyle("Ai2",  "end",   -8,  4, "#BF360C", 2.5, 8.5),
        "Ti2":  VertexStyle("Ti2",  "start",  8, 12, "#BF360C", 2.5, 8.5),
    },
    brg_dist_labels=None, arc_labels=None, center_marks=None,
    traverse_pts=None, traverse_stroke=None,
)

def build_outline_cfg(outline_segs, pts, radii):
    """Build outline layer config (needs computed sweep angles and radii)."""
    R = radii
    sw = {i: arc_sweep_deg(outline_segs[i], pts)
          for i in [0,2,4,5,7,8,10,12,14,16]}
    return LayerConfig(
        opacity=1.0, fill_color="rgba(200,230,255,0.25)",
        line_stroke="#333", line_width=2.0,
        arc_styles={
            ("F1","F2"):    ("#333", 2.0),
            ("F5","F6"):    ("#333", 2.0),
            ("F7","F8"):    ("#333", 2.0),
            ("F8","F9"):    ("#333", 2.0),
            ("F10","F11"):  ("#333", 2.0),
            ("F11","F11a"): ("#333", 2.0),
            ("F11b","F12"): ("#333", 2.0),
            ("F13","F14"):  ("#333", 2.0),
            ("F15","F16"):  ("#333", 2.0),
            ("F17","F18"):  ("#333", 2.0),
        },
        vertex_styles={
            "F5":   VertexStyle("F5",   "end",    -8,  4,  "#d32f2f", 1.75, 10),
            "F2":   VertexStyle("F2",   "end",    -8,  4,  "#d32f2f", 1.75, 10),
            "F8":   VertexStyle("F8",   "end",    -8,  4,  "#d32f2f", 1.75, 10),
            "F11":  VertexStyle("F11",  "start",   8,  4,  "#d32f2f", 1.75, 10),
            "F11a": VertexStyle("F11a", "end",    -2, -6,  "#d32f2f", 1.75, 10),
            "F11b": VertexStyle("F11b", "start",   2, -6,  "#d32f2f", 1.75, 10),
            "F12":  VertexStyle("F12",  "start",   8,  4,  "#d32f2f", 1.75, 10),
            "F13":  VertexStyle("F13",  "start",   8,  4,  "#d32f2f", 1.75, 10),
            "F14":  VertexStyle("F14",  "start",  10,  4,  "#d32f2f", 1.75, 10),
            "F15":  VertexStyle("F15",  "start",   8,  4,  "#d32f2f", 1.75, 10),
            "F1":   VertexStyle("F1",   "middle",  0, 10,  "#d32f2f", 1.75, 10),
            "F6":   VertexStyle("F6",   "middle",  0, -6,  "#d32f2f", 1.75, 10),
            "F7":   VertexStyle("F7",   "middle",  0, -6,  "#d32f2f", 1.75, 10),
            "F9":   VertexStyle("F9",   "middle",  0, 17,  "#d32f2f", 1.75, 10),
            "F10":  VertexStyle("F10",  "middle",  0, 17,  "#d32f2f", 1.75, 10),
            "F17":  VertexStyle("F17",  "middle",  0, 13,  "#d32f2f", 1.75, 10),
            "F18":  VertexStyle("F18",  "middle",  0, 12,  "#d32f2f", 1.75, 10),
            "F16":  VertexStyle("F16",  "start",   8,  4,  "#d32f2f", 1.75, 10),
            "FC":   VertexStyle("FC",   "start",   4,  -2, "#d32f2f", 1.75, 10),
        },
        brg_dist_labels={
            ("F2","F5"): BrgDistLabel(18),
            ("F6","F7"): BrgDistLabel(-16),
            ("F9","F10"): BrgDistLabel(-16),
            ("F11a","F11b"): BrgDistLabel(-28),
            ("F12","F13"): BrgDistLabel(-16),
            ("F14","F15"): BrgDistLabel(-16),
            ("F16","F17"): BrgDistLabel(16),
            ("F18","F1"): BrgDistLabel(16),
        },
        arc_labels={
            ("F1","F2"): ArcLabel(f"Arc R={R['R_a1']*12:.0f}\"",
                f"{sw[0]:.1f}\u00b0", "end", -10, 14, 11, "#333"),
            ("F5","F6"): ArcLabel(f"Arc R={R['R_a5']*12:.0f}\"",
                f"{sw[2]:.1f}\u00b0", "end", -10, -14, 11, "#333"),
            ("F7","F8"): ArcLabel(f"Arc R={R['R_a7']*12:.0f}\"",
                f"{sw[4]:.1f}\u00b0", "start", 12, 0, 11, "#333"),
            ("F8","F9"): ArcLabel(f"Arc R={R['R_a8']*12:.0f}\"",
                f"{sw[5]:.1f}\u00b0", "end", -10, 14, 11, "#333"),
            ("F10","F11"): ArcLabel(f"Arc R={R['R_a10']*12:.1f}\"",
                f"{sw[7]:.1f}\u00b0", "end", -10, -10, 11, "#333"),
            ("F11","F11a"): ArcLabel(f"Arc R={R['R_a11']*12:.0f}\"",
                f"{sw[8]:.1f}\u00b0", "end", -10, -20, 11, "#333"),
            ("F11b","F12"): ArcLabel(f"Arc R={R['R_a11']*12:.0f}\"",
                f"{sw[10]:.1f}\u00b0", "start", 12, 0, 11, "#333"),
            ("F13","F14"): ArcLabel(f"Arc R={R['R_a13']*12:.0f}\"",
                f"{sw[12]:.1f}\u00b0", "start", 12, 0, 11, "#333"),
            ("F15","F16"): ArcLabel(f"Arc R={R['R_a15']*12:.0f}\"",
                f"{sw[14]:.1f}\u00b0", "start", 10, -10, 11, "#333"),
            ("F17","F18"): ArcLabel(f"Arc R={R['R_a17']*12:.0f}\"",
                f"{sw[16]:.1f}\u00b0", "end", -10, -10, 11, "#333"),
        },
        center_marks=[
            CenterMark("C1", "F1", "#333"),
            CenterMark("C5", "F5", "#333"),
            CenterMark("C7", "F7", "#333"), CenterMark("C8", "F8", "#333"),
            CenterMark("C10", "F10", "#333"),
            CenterMark("C11a", "F11", "#333"), CenterMark("C11", "F11b", "#333"),
            CenterMark("C13", "F13", "#333"), CenterMark("C15", "F15", "#333"),
            CenterMark("C17", "F17", "#333"),
        ],
        traverse_pts=None, traverse_stroke=None,
        brg_decimal=True,
    )

def build_outline_cfg_db(outline_segs, pts):
    """Build outline layer config dynamically from DB chain naming convention.

    Works with any point naming (PA/PB, F-series, or mixed).  Labels are
    placed radially outward from the outline centroid so they don't pile up.
    """
    all_names = list(dict.fromkeys(
        n for s in outline_segs for n in (s.start, s.end) if n in pts
    ))
    if all_names:
        cx = sum(pts[n][0] for n in all_names) / len(all_names)
        cy = sum(pts[n][1] for n in all_names) / len(all_names)
    else:
        cx, cy = 0.0, 0.0

    def _outward(name, label_r=12):
        de = pts[name][0] - cx
        dn = pts[name][1] - cy
        dist = math.hypot(de, dn) or 1.0
        dx = de / dist * label_r
        dy = -dn / dist * label_r  # SVG y-down inverts N
        return dx, dy, "start" if dx >= 0 else "end"

    arc_styles = {}; arc_labels = {}; brg_dist_labels = {}; vertex_styles = {}
    center_marks = []
    seen_verts = set()
    for seg in outline_segs:
        for name in (seg.start, seg.end):
            if name not in seen_verts:
                seen_verts.add(name)
                if name in pts:
                    dx, dy, anchor = _outward(name)
                    vertex_styles[name] = VertexStyle(name, anchor, dx, dy, "#333", 2.5, 10)
        if isinstance(seg, ArcSeg):
            arc_styles[(seg.start, seg.end)] = ("#333", 2.0)
            sw = arc_sweep_deg(seg, pts)
            r_in = seg.radius * 12
            r_str = f"{r_in:.0f}\"" if abs(r_in - round(r_in)) < 0.05 else f"{r_in:.1f}\""
            arc_labels[(seg.start, seg.end)] = ArcLabel(
                f"R={r_str}", f"{sw:.1f}\u00b0", "start", 10, 4, 11, "#333")
            if seg.center in pts:
                center_marks.append(CenterMark(seg.center, seg.start, "#333"))
        else:
            brg_dist_labels[(seg.start, seg.end)] = BrgDistLabel(16)

    return LayerConfig(
        opacity=1.0, fill_color="rgba(200,230,255,0.25)",
        line_stroke="#333", line_width=2.0,
        arc_styles=arc_styles, vertex_styles=vertex_styles,
        brg_dist_labels=brg_dist_labels, arc_labels=arc_labels,
        center_marks=center_marks if center_marks else None,
        traverse_pts=None, traverse_stroke=None,
        brg_decimal=True,
    )

# ============================================================
# Section 9: SVG Assembly
# ============================================================
def generate_svg(data, iw_polys=None, out_path=None):
    """Assemble and write path_area.svg.

    iw_polys — dict {name: poly} of DB-sourced interior wall polygons.
               If None, falls back to the seed layout in data["layout"].
    out_path  — output file path; defaults to survey/path_area.svg next to this file.
    """
    pts = data["pts"]; to_svg = data["to_svg"]
    outer_segs = data["outer_segs"]; inset_segs = data["inset_segs"]
    outline_segs = data["outline_segs"]
    outer_area = data["outer_area"]; inset_area = data["inset_area"]
    outline_area = data["outline_area"]

    _seg_starts = {s.start for s in outline_segs}
    if any(n.startswith("F") and n[1:].lstrip("0123456789") in ("", "a", "b") for n in _seg_starts):
        outline_cfg = build_outline_cfg(outline_segs, pts, data["radii"])
    else:
        outline_cfg = build_outline_cfg_db(outline_segs, pts)

    lines = []
    lines.append(f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}">')
    lines.append(f'<rect width="{W}" height="{H}" fill="white"/>')
    lines.append('<defs>')
    lines.append('  <marker id="ah" markerWidth="8" markerHeight="6" refX="8" refY="3" orient="auto"><polygon points="0 0, 8 3, 0 6" fill="#333"/></marker>')
    lines.append(f'  <clipPath id="page"><rect width="{W}" height="{H}"/></clipPath>')
    lines.append('</defs>')
    lines.append(f'<text x="{W/2}" y="30" text-anchor="middle" font-family="Arial" font-size="14"'
                 f' font-weight="bold">Site Path \u2014 Outline</text>')

    render_layer(lines, outer_segs, pts, outer_cfg, to_svg)
    render_layer(lines, inset_segs, pts, inset_cfg, to_svg)
    render_layer(lines, outline_segs, pts, outline_cfg, to_svg)
    render_floorplan(lines, to_svg, pts, data["outer_poly"], data["inner_poly"],
                     data["inner_segs"], layout=data.get("layout"), iw_polys=iw_polys)

    # Area label: halfway between FC and midpoint of W9-W10 (fall back to outline centroid)
    if "W9" in pts and "W10" in pts and "FC" in pts:
        _w9w10_mid = ((pts["W9"][0] + pts["W10"][0]) / 2, (pts["W9"][1] + pts["W10"][1]) / 2)
        cx_o = (pts["FC"][0] + _w9w10_mid[0]) / 2
        cy_o = (pts["FC"][1] + _w9w10_mid[1]) / 2
    elif "FC" in pts:
        cx_o, cy_o = pts["FC"]
    else:
        _op = data["outer_poly"]
        cx_o = sum(p[0] for p in _op) / len(_op)
        cy_o = sum(p[1] for p in _op) / len(_op)
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

    # Footer
    _now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    _git_desc = git_describe()
    lines.append(f'<text x="{W/2}" y="{H-2}" text-anchor="middle" font-family="Arial" font-size="7.5"'
                 f' fill="#999">Generated {_now} from {_git_desc}</text>')
    lines.append('</svg>')

    svg_content = "\n".join(lines)
    if out_path is None:
        out_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "path_area.svg")
    with open(out_path, "w", encoding="utf-8") as f:
        f.write(svg_content)
    return svg_content


if __name__ == "__main__":
    data = compute_all()
    pts = data["pts"]; radii = data["radii"]

    print(f'=== INSET PATH (6" inside) ===')
    print(f"  delta=0.5' R1i={data['R1i']}' R2i={data['R2i']}' R3i={data['R3i']}'")
    print(f"  Inset area: {data['inset_area']:.2f} sq ft")
    print(f'=== OUTLINE PATH ===')
    _pt_notes = [
        ("F1", "arc tangent"), ("F2", "arc tangent"),
        ("F5", "arc tangent"),
        ("F6", "arc tangent"), ("F7", "east of F6"),
        ("F8", "C7/C8 arc junction"), ("F9", "arc tangent"),
        ("F10", "arc E-W tangent"),
        ("F11", "arc tangent C10/C11a"),
        ("F11a", "top of C11a / flat seg west"),
        ("F11b", "top of C11 / flat seg east"),
        ("F12", "line / arc tangent"),
        ("F13", f"{radii['R_a13']*12:.1f}\" arc / line tangent"),
        ("F14", f"{radii['R_a13']*12:.1f}\" arc tangent to N-S line"),
        ("F15", "arc C15, exits North"),
        ("F16", "arc C15, incoming tangent"),
        ("F17", "on PiX-Pi5 line"), ("F18", "arc C17 tangent"),
    ]
    for name, note in _pt_notes:
        print(f"  {name:<5s} ({pts[name][0]:.4f}, {pts[name][1]:.4f})  ({note})")
    print(f"  C17:  ({pts['C17'][0]:.4f}, {pts['C17'][1]:.4f})  (F17->F18 arc center)")
    print(f"  C15:  ({pts['C15'][0]:.4f}, {pts['C15'][1]:.4f})  (arc C15 center, R={radii['R_a15']:.4f}')")
    print(f"  C13:  ({pts['C13'][0]:.4f}, {pts['C13'][1]:.4f})  ({radii['R_a13']*12:.1f}\" arc center, R={radii['R_a13']:.4f}')")
    print(f"  Outline area: {data['outline_area']:.2f} sq ft")

    generate_svg(data)
    print(f"\nSVG written to path_area.svg")
    print(f"Outer path area: {data['outer_area']:.2f} sq ft (rendered at 20%)")
    print(f"Inset path area: {inset_area:.2f} sq ft (rendered at 20%)")
    print(f"Outline path area: {outline_area:.2f} sq ft (rendered at 100%)")
    print(f"Outline: F1->ArcC1->F2->F5->ArcC5->F6->F7->ArcC7->F8->ArcC8->F9->F10->ArcC10->F11->ArcC11a->F11a->F11b->ArcC11->F12->F13->ArcC13->F14->F15->ArcC15->F16->F17->ArcC17->F18->F1")
