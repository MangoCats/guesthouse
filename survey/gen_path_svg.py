import os, math, datetime
from typing import NamedTuple

from shared.types import Point, LineSeg, ArcSeg, Segment
from shared.geometry import (
    poly_area, left_norm, off_pt, line_isect, arc_poly,
    circle_circle_isect, line_circle_isect_min_t_gt, line_circle_isect_min_abs_t,
    segment_polyline, path_polygon, arc_sweep_deg,
    brg_dist, fmt_brg, fmt_dist,
    compute_inner_walls,
)
from shared.survey import compute_traverse, compute_three_arc, compute_inset
from shared.svg import make_svg_transform, W, H, git_describe
from floorplan.geometry import compute_outline_geometry, OutlineAnchors
from floorplan.layout import compute_interior_layout
from floorplan.constants import (
    WALL_OUTER, WALL_3IN, WALL_4IN,
    APPLIANCE_OFFSET_E, APPLIANCE_WIDTH,
    COUNTER_GAP, COUNTER_DEPTH,
    CLOSET_WIDTH, BEDROOM_WIDTH,
    ARC_180_R, F15_OFFSET_E, PIX_PI5_TARGET_BRG,
)

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
        trav = " ".join(f"{to_svg(*pts[n])[0]:.1f},{to_svg(*pts[n])[1]:.1f}"
                        for n in cfg.traverse_pts)
        lines.append(f'<polygon points="{trav}" fill="none" stroke="{cfg.traverse_stroke}"'
                     f' stroke-width="0.8" stroke-dasharray="4,4"/>')

    # Filled polygon
    dense = path_polygon(segments, pts)
    poly_svg = " ".join(f"{to_svg(*p)[0]:.1f},{to_svg(*p)[1]:.1f}" for p in dense)
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
        svg_pts = " ".join(f"{to_svg(e,n)[0]:.1f},{to_svg(e,n)[1]:.1f}" for e, n in poly)
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
            ang = math.degrees(math.atan2(-(sy2-sy1), sx2-sx1))
            if ang > 90: ang -= 180
            if ang < -90: ang += 180
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
    pts, p3_trav = compute_traverse()
    to_svg = make_svg_transform(p3_trav)
    arc_info = compute_three_arc(pts)
    R1, R2, R3 = arc_info["R1"], arc_info["R2"], arc_info["R3"]
    nE, nN = arc_info["nE"], arc_info["nN"]

    # Outer path
    outer_segs = [
        LineSeg("POB", "P2"), LineSeg("P2", "P3"), LineSeg("P3", "T3"),
        ArcSeg("T3", "PX", "TC3", R3, "CW", 60),
        LineSeg("PX", "P4"), LineSeg("P4", "P5"), LineSeg("P5", "T1"),
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

    # Rotate outer/inset points about pivot so PiX-Pi5 bearing = 60°.
    # The pivot is chosen so that F15 (after outline geometry) lands at the
    # correct easting — computed here from the pre-rotation dimension chain.
    _R_fp = ARC_180_R
    _d_pip = (pts["Pi5"][0] - pts["PiX"][0], pts["Pi5"][1] - pts["PiX"][1])
    _L_pip = math.sqrt(_d_pip[0]**2 + _d_pip[1]**2)
    _d_pip_u = (_d_pip[0]/_L_pip, _d_pip[1]/_L_pip)
    _ln_pip = left_norm(pts["PiX"], pts["Pi5"])
    _o_pip = off_pt(pts["PiX"], _ln_pip, _R_fp)
    _pre_U1_E = pts["Pi3"][0]
    # Dimension chain west→east: outer wall + appliance offset + appliance width
    # + counter gap + counter depth + closet/bedroom/closet walls
    _pre_iw8_e = (_pre_U1_E + WALL_OUTER + APPLIANCE_OFFSET_E + APPLIANCE_WIDTH
                 + COUNTER_GAP + COUNTER_DEPTH
                 + WALL_3IN + CLOSET_WIDTH + WALL_4IN + BEDROOM_WIDTH
                 + WALL_4IN + CLOSET_WIDTH + WALL_3IN)
    _pre_F15_E = _pre_iw8_e + F15_OFFSET_E
    _t_cf4 = (_pre_F15_E - _R_fp - _o_pip[0]) / _d_pip_u[0]
    _cf4 = (_pre_F15_E - _R_fp, _o_pip[1] + _t_cf4 * _d_pip_u[1])
    _t_o12 = ((_cf4[0]-pts["PiX"][0])*_d_pip[0] + (_cf4[1]-pts["PiX"][1])*_d_pip[1]) \
             / (_d_pip[0]**2 + _d_pip[1]**2)
    _pivot = (pts["PiX"][0] + _t_o12*_d_pip[0], pts["PiX"][1] + _t_o12*_d_pip[1])
    _brg_pip = math.degrees(math.atan2(_d_pip[0], _d_pip[1])) % 360
    _rot_deg = PIX_PI5_TARGET_BRG - _brg_pip
    _rot_rad = math.radians(_rot_deg)
    _cos_r = math.cos(_rot_rad)
    _sin_r = math.sin(_rot_rad)
    def _rot_cw(p):
        dE, dN = p[0] - _pivot[0], p[1] - _pivot[1]
        return (_pivot[0] + dE*_cos_r + dN*_sin_r, _pivot[1] - dE*_sin_r + dN*_cos_r)
    pts_rot = dict(pts)
    for _n in ["POB","P2","P3","P4","P5","T1","T2","T3","PA","PX","TC1","TC2","TC3",
               "PiOB","Pi2","Pi3","Pi4","Pi5","Ti1","Ti2","Ti3","Ai2","PiX"]:
        pts_rot[_n] = _rot_cw(pts[_n])

    # F-series outline geometry (from floorplan)
    anchors = OutlineAnchors(
        Pi2=pts["Pi2"], Pi3=pts["Pi3"], Ti3=pts["Ti3"],
        PiX=pts["PiX"], Pi5=pts["Pi5"],
        TC1=pts["TC1"], R1i=inset.R1i,
    )
    outline_geo = compute_outline_geometry(anchors)
    pts.update(outline_geo.fp_pts)
    outline_segs = outline_geo.outline_segs
    radii = outline_geo.radii

    # Derive U-series as aliases (downstream from F-series)
    for i in range(22):
        pts[f"U{i}"] = outline_geo.fp_pts[f"F{i}"]

    # Inner walls + layout
    inner_segs = compute_inner_walls(outline_segs, pts, WALL_OUTER, radii)
    outer_poly = path_polygon(outline_segs, pts)
    inner_poly = path_polygon(inner_segs, pts)
    outline_area = poly_area(outer_poly)
    layout = compute_interior_layout(pts, inner_poly)

    return {
        "pts": pts, "pts_rot": pts_rot, "to_svg": to_svg,
        "outer_segs": outer_segs, "inset_segs": inset_segs, "outline_segs": outline_segs,
        "outer_area": outer_area, "inset_area": inset_area, "outline_area": outline_area,
        "R1i": inset.R1i, "R2i": inset.R2i, "R3i": inset.R3i,
        "radii": radii, "outer_poly": outer_poly, "inner_poly": inner_poly,
        "inner_segs": inner_segs, "layout": layout,
    }

def render_floorplan(lines, to_svg, pts, outer_poly, inner_poly, inner_segs, layout):
    L = layout
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
            sp = " ".join(f"{to_svg(e,n)[0]:.1f},{to_svg(e,n)[1]:.1f}" for e,n in pl)
            lines.append(f'<polyline points="{sp}" fill="none" stroke="#666" stroke-width="1.0" stroke-linecap="round"/>')
    # IW1
    iw1 = L.iw1
    svg = " ".join(f"{to_svg(*p)[0]:.1f},{to_svg(*p)[1]:.1f}" for p in iw1)
    lines.append(f'<polygon points="{svg}" fill="rgba(160,160,160,0.5)" stroke="none"/>')
    for a,b in [(iw1[0],iw1[1]),(iw1[3],iw1[2])]:
        s1,s2 = to_svg(*a),to_svg(*b)
        lines.append(f'<line x1="{s1[0]:.1f}" y1="{s1[1]:.1f}" x2="{s2[0]:.1f}" y2="{s2[1]:.1f}" stroke="#666" stroke-width="1.0"/>')
    # IW2
    iw2 = [(L.iw2_w,L.iw2_s),(L.iw2_e,L.iw2_s),(L.iw2_e,L.iw2_n),(L.iw2_w,L.iw2_n)]
    svg = " ".join(f"{to_svg(*p)[0]:.1f},{to_svg(*p)[1]:.1f}" for p in iw2)
    lines.append(f'<polygon points="{svg}" fill="rgba(160,160,160,0.5)" stroke="none"/>')
    for ev in [L.iw2_w,L.iw2_e]:
        s1,s2 = to_svg(ev,L.iw2_s),to_svg(ev,L.iw2_n)
        lines.append(f'<line x1="{s1[0]:.1f}" y1="{s1[1]:.1f}" x2="{s2[0]:.1f}" y2="{s2[1]:.1f}" stroke="#666" stroke-width="1.0"/>')
    # IW7 L-shape
    svg = " ".join(f"{to_svg(*p)[0]:.1f},{to_svg(*p)[1]:.1f}" for p in L.iw7)
    lines.append(f'<polygon points="{svg}" fill="rgba(160,160,160,0.5)" stroke="#666" stroke-width="0.8"/>')
    # IW3 (west bedroom wall)
    iw3_poly = [(L.iw3_w,L.iw3_s),(L.iw3_e,L.iw3_s),(L.iw3_e,L.iw3_n),(L.iw3_w,L.iw3_n)]
    svg = " ".join(f"{to_svg(*p)[0]:.1f},{to_svg(*p)[1]:.1f}" for p in iw3_poly)
    lines.append(f'<polygon points="{svg}" fill="rgba(160,160,160,0.5)" stroke="#666" stroke-width="0.8"/>')
    # IW4 (east bedroom wall)
    wsn = L.wall_south_n
    iw4_poly = [(L.iw4_w,wsn),(L.iw4_e,wsn),(L.iw4_e,L.iw1_s),(L.iw4_w,L.iw1_s)]
    svg = " ".join(f"{to_svg(*p)[0]:.1f},{to_svg(*p)[1]:.1f}" for p in iw4_poly)
    lines.append(f'<polygon points="{svg}" fill="rgba(160,160,160,0.5)" stroke="#666" stroke-width="0.8"/>')
    # IW8 L-shape
    cl1_top = L.cl1_top; iwt3 = L.iwt3; iw8_w = L.iw8_w; iw8_e = L.iw8_e
    iw8_poly = [(L.iw4_e,cl1_top+iwt3),(iw8_e,cl1_top+iwt3),(iw8_e,wsn),
                (iw8_w,wsn),(iw8_w,cl1_top),(L.iw4_e,cl1_top)]
    svg = " ".join(f"{to_svg(*p)[0]:.1f},{to_svg(*p)[1]:.1f}" for p in iw8_poly)
    lines.append(f'<polygon points="{svg}" fill="rgba(160,160,160,0.5)" stroke="#666" stroke-width="0.8"/>')
    # Appliances
    for lbl,sw_e,sw_n,ne_e,ne_n in [("DRYER",L.dryer_w,L.dryer_s,L.dryer_e,L.dryer_n),
                                      ("WASHER",L.washer_w,L.washer_s,L.washer_e,L.washer_n)]:
        s1,s2 = to_svg(sw_e,ne_n),to_svg(ne_e,sw_n); w=s2[0]-s1[0]; h=s2[1]-s1[1]
        lines.append(f'<rect x="{s1[0]:.1f}" y="{s1[1]:.1f}" width="{w:.1f}" height="{h:.1f}" fill="rgba(100,150,200,0.3)" stroke="#4682B4" stroke-width="0.8"/>')
        cx,cy = (s1[0]+s2[0])/2,(s1[1]+s2[1])/2
        lines.append(f'<text x="{cx:.1f}" y="{cy+3:.1f}" text-anchor="middle" font-family="Arial" font-size="7" fill="#4682B4">{lbl}</text>')
    # Counter with rounded NW corner
    ctr_w, ctr_e, ctr_s, ctr_n = L.ctr_w, L.ctr_e, L.ctr_s, L.ctr_n
    ctr_nw_r = L.ctr_nw_r
    csw,cse,cne = to_svg(ctr_w,ctr_s),to_svg(ctr_e,ctr_s),to_svg(ctr_e,ctr_n)
    cnas,cnae = to_svg(ctr_w+ctr_nw_r,ctr_n),to_svg(ctr_w,ctr_n-ctr_nw_r)
    rsv = abs(cnas[0]-to_svg(ctr_w,ctr_n)[0])
    cd = (f'M {csw[0]:.1f},{csw[1]:.1f} L {cse[0]:.1f},{cse[1]:.1f} L {cne[0]:.1f},{cne[1]:.1f} '
          f'L {cnas[0]:.1f},{cnas[1]:.1f} A {rsv:.1f} {rsv:.1f} 0 0 0 {cnae[0]:.1f},{cnae[1]:.1f} Z')
    lines.append(f'<path d="{cd}" fill="rgba(100,150,200,0.3)" stroke="#4682B4" stroke-width="0.8"/>')
    ccx,ccy = (csw[0]+cse[0])/2,(csw[1]+cne[1])/2
    lines.append(f'<text x="{ccx:.1f}" y="{ccy:.1f}" text-anchor="middle" font-family="Arial" font-size="7" fill="#4682B4" letter-spacing="0.5" transform="rotate(-90,{ccx:.1f},{ccy:.1f})">COUNTER</text>')
    # King Bed
    bs,be = to_svg(L.bed_w,L.bed_n),to_svg(L.bed_e,L.bed_s); bw,bh = be[0]-bs[0],be[1]-bs[1]
    lines.append(f'<rect x="{bs[0]:.1f}" y="{bs[1]:.1f}" width="{bw:.1f}" height="{bh:.1f}" fill="rgba(100,150,200,0.3)" stroke="#4682B4" stroke-width="0.8"/>')
    bcx,bly = (bs[0]+be[0])/2,bs[1]+0.765*bh
    lines.append(f'<text x="{bcx:.1f}" y="{bly+3:.1f}" text-anchor="middle" font-family="Arial" font-size="7" fill="#4682B4">KING BED</text>')
    # Room labels
    cx,cy = to_svg((ctr_e+iwt3+L.iw3_w)/2,(ctr_s+ctr_n)/2)
    lines.append(f'<text x="{cx:.1f}" y="{cy+3:.1f}" text-anchor="middle" font-family="Arial" font-size="7" fill="#666" transform="rotate(-90,{cx:.1f},{cy+3:.1f})">CLOSET</text>')
    bx,by = to_svg((L.iw3_e+L.iw4_w)/2,(ctr_s+L.iw1_s)/2)
    lines.append(f'<text x="{bx:.1f}" y="{by+3:.1f}" text-anchor="middle" font-family="Arial" font-size="8" fill="#666">BEDROOM</text>')
    cx,cy = to_svg((L.iw4_e+iw8_w)/2,(ctr_s+cl1_top)/2)
    lines.append(f'<text x="{cx:.1f}" y="{cy+3:.1f}" text-anchor="middle" font-family="Arial" font-size="7" fill="#666" transform="rotate(-90,{cx:.1f},{cy+3:.1f})">CLOSET</text>')
    lines.append('</g>')

# ============================================================
# Section 8: Style Configurations (static — no computed values)
# ============================================================
outer_cfg = LayerConfig(
    opacity=0.2, fill_color="#e8edf5",
    line_stroke="#333", line_width=2.0,
    arc_styles={
        ("T1","PA"): ("#0077B6", 2.5), ("PA","T2"): ("#0077B6", 2.5),
        ("T3","PX"): ("#2E7D32", 2.5),
    },
    vertex_styles={
        "POB": VertexStyle("P.O.B.", "middle", 0, -16, "#d32f2f", 3.5, 10),
        "P2":  VertexStyle("P2",  "end",  -10, -8, "#d32f2f", 3.5, 10),
        "P3":  VertexStyle("P3",  "end",  -10,  7, "#d32f2f", 3.5, 10),
        "T3":  VertexStyle("T3",  "middle", 0, -12, "#2E7D32", 3.5, 10),
        "PX":  VertexStyle("PX",  "start",  8, -6, "#2E7D32", 3.5, 10),
        "P4":  VertexStyle("P4",  "start", 10, 10, "#d32f2f", 3.5, 10),
        "P5":  VertexStyle("P5",  "start", 14,  3, "#d32f2f", 3.5, 10),
        "T1":  VertexStyle("T1",  "start", 10, -4, "#0077B6", 3.5, 10),
        "PA":  VertexStyle("PA",  "start",  8, -8, "#0077B6", 3.5, 10),
        "T2":  VertexStyle("T2",  "end",  -10, -4, "#0077B6", 3.5, 10),
    },
    brg_dist_labels={
        ("POB","P2"): BrgDistLabel(-18), ("P2","P3"): BrgDistLabel(-18),
        ("P3","T3"): BrgDistLabel(18), ("PX","P4"): BrgDistLabel(-16),
        ("P4","P5"): BrgDistLabel(-16), ("P5","T1"): BrgDistLabel(16),
        ("T2","POB"): BrgDistLabel(16),
    },
    arc_labels={
        ("T1","PA"): ArcLabel("Arc 1: R=10'", None, "end", -20, 0, 0, "#0077B6"),
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
        "Pi5":  VertexStyle("Pi5",  "end",   -8, -8, "#BF360C", 2.5, 8.5),
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
          for i in [0,2,3,5,7,8,10,11,13,15,17,19,20]}
    return LayerConfig(
        opacity=1.0, fill_color="rgba(200,230,255,0.25)",
        line_stroke="#333", line_width=2.0,
        arc_styles={
            ("F0","F1"):    ("#333", 2.0),
            ("F2","F3"):    ("#333", 2.0),
            ("F3","F4"):    ("#333", 2.0),
            ("F5","F6"):    ("#333", 2.0),
            ("F7","F8"):    ("#333", 2.0),
            ("F8","F9"):    ("#333", 2.0),
            ("F10","F11"):  ("#333", 2.0),
            ("F11","F12"):  ("#333", 2.0),
            ("F13","F14"):  ("#333", 2.0),
            ("F15","F16"):  ("#333", 2.0),
            ("F17","F18"):  ("#333", 2.0),
            ("F19","F20"):  ("#333", 2.0),
            ("F20","F21"):  ("#333", 2.0),
        },
        vertex_styles={
            "F5":   VertexStyle("F5",   "end",    -8,  4,  "#d32f2f", 1.75, 10),
            "F4":   VertexStyle("F4",   "end",    -8,  4,  "#d32f2f", 1.75, 10),
            "F3":   VertexStyle("F3",   "end",   -10,  4,  "#d32f2f", 1.75, 10),
            "F2":   VertexStyle("F2",   "end",    -8,  4,  "#d32f2f", 1.75, 10),
            "F1":   VertexStyle("F1",   "end",    -8,  4,  "#d32f2f", 1.75, 10),
            "F8":   VertexStyle("F8",   "end",    -8,  4,  "#d32f2f", 1.75, 10),
            "F11":  VertexStyle("F11",  "start",   8,  4,  "#d32f2f", 1.75, 10),
            "F12":  VertexStyle("F12",  "start",   8,  4,  "#d32f2f", 1.75, 10),
            "F13":  VertexStyle("F13",  "start",   8,  4,  "#d32f2f", 1.75, 10),
            "F14":  VertexStyle("F14",  "start",  10,  4,  "#d32f2f", 1.75, 10),
            "F15":  VertexStyle("F15",  "start",   8,  4,  "#d32f2f", 1.75, 10),
            "F0":   VertexStyle("F0",   "middle",  0, 10,  "#d32f2f", 1.75, 10),
            "F6":   VertexStyle("F6",   "middle",  0, -6,  "#d32f2f", 1.75, 10),
            "F7":   VertexStyle("F7",   "middle",  0, -6,  "#d32f2f", 1.75, 10),
            "F9":   VertexStyle("F9",   "middle",  0, 17,  "#d32f2f", 1.75, 10),
            "F10":  VertexStyle("F10",  "middle",  0, 17,  "#d32f2f", 1.75, 10),
            "F17":  VertexStyle("F17",  "middle",  0, 13,  "#d32f2f", 1.75, 10),
            "F18":  VertexStyle("F18",  "middle",  0, 12,  "#d32f2f", 1.75, 10),
            "F19":  VertexStyle("F19",  "middle",  0, 12,  "#d32f2f", 1.75, 10),
            "F20":  VertexStyle("F20",  "middle",  0, 13,  "#d32f2f", 1.75, 10),
            "F21":  VertexStyle("F21",  "middle",  0, 10,  "#d32f2f", 1.75, 10),
            "F16":  VertexStyle("F16",  "start",   8,  4,  "#d32f2f", 1.75, 10),
        },
        brg_dist_labels={
            ("F1","F2"): BrgDistLabel(18),
            ("F4","F5"): BrgDistLabel(18),
            ("F6","F7"): BrgDistLabel(-16),
            ("F9","F10"): BrgDistLabel(-16),
            ("F12","F13"): BrgDistLabel(-16),
            ("F14","F15"): BrgDistLabel(-16),
            ("F16","F17"): BrgDistLabel(16),
            ("F18","F19"): BrgDistLabel(-16),
            ("F21","F0"): BrgDistLabel(-18),
        },
        arc_labels={
            ("F0","F1"): ArcLabel("Arc R=10\"",
                f"{sw[0]:.1f}\u00b0", "end", -10, 14, 11, "#333"),
            ("F2","F3"): ArcLabel(f"Arc R={R['R_a2']*12:.1f}\"",
                f"{sw[2]:.1f}\u00b0 CCW", "start", 12, 0, 11, "#333"),
            ("F3","F4"): ArcLabel(f"Arc R={R['R_a3']*12:.1f}\"",
                f"{sw[3]:.1f}\u00b0 CW", "end", -10, -10, 11, "#333"),
            ("F5","F6"): ArcLabel("Arc R=28\"",
                f"{sw[5]:.1f}\u00b0", "end", -10, -14, 11, "#333"),
            ("F7","F8"): ArcLabel("Arc R=28\"",
                f"{sw[7]:.1f}\u00b0", "start", 12, 0, 11, "#333"),
            ("F8","F9"): ArcLabel("Arc R=2\"",
                f"{sw[8]:.1f}\u00b0", "end", -10, 14, 11, "#333"),
            ("F10","F11"): ArcLabel("Arc R=2\"",
                f"{sw[10]:.1f}\u00b0", "end", -10, -10, 11, "#333"),
            ("F11","F12"): ArcLabel("Arc R=28\"",
                f"{sw[11]:.1f}\u00b0 CCW", "start", 12, 0, 11, "#333"),
            ("F13","F14"): ArcLabel(f"Arc R={R['R_a13']*12:.0f}\"",
                f"{sw[13]:.1f}\u00b0", "start", 12, 0, 11, "#333"),
            ("F15","F16"): ArcLabel(f"Arc R={R['R_a15']*12:.0f}\"",
                f"{sw[15]:.1f}\u00b0", "start", 10, -10, 11, "#333"),
            ("F17","F18"): ArcLabel(f"Arc R={R['R_a17']*12:.0f}\"",
                f"{sw[17]:.1f}\u00b0", "end", -10, -10, 11, "#333"),
            ("F19","F20"): ArcLabel(f"Arc R={R['R_a19']*12:.1f}\"",
                f"{sw[19]:.1f}\u00b0", "start", 12, -10, 11, "#333"),
            ("F20","F21"): ArcLabel(f"Arc R={R['R_a20']*12:.1f}\"",
                f"{sw[20]:.1f}\u00b0 CW", "start", 12, 4, 11, "#333"),
        },
        center_marks=[
            CenterMark("C0", "F0", "#333"), CenterMark("C2", "F2", "#333"),
            CenterMark("C3", "F3", "#333"), CenterMark("C5", "F5", "#333"),
            CenterMark("C7", "F7", "#333"), CenterMark("C8", "F8", "#333"),
            CenterMark("C10", "F10", "#333"), CenterMark("C11", "F11", "#333"),
            CenterMark("C13", "F13", "#333"), CenterMark("C15", "F15", "#333"),
            CenterMark("C17", "F17", "#333"), CenterMark("C19", "F19", "#333"),
            CenterMark("C20", "F20", "#333"),
        ],
        traverse_pts=None, traverse_stroke=None,
        brg_decimal=True,
    )

# ============================================================
# Section 9: SVG Assembly + Output (only when run directly)
# ============================================================
if __name__ == "__main__":
    data = compute_all()
    pts = data["pts"]; pts_rot = data["pts_rot"]; to_svg = data["to_svg"]
    outer_segs = data["outer_segs"]; inset_segs = data["inset_segs"]
    outline_segs = data["outline_segs"]
    outer_area = data["outer_area"]; inset_area = data["inset_area"]
    outline_area = data["outline_area"]
    radii = data["radii"]

    print(f'=== INSET PATH (6" inside) ===')
    print(f"  delta=0.5' R1i={data['R1i']}' R2i={data['R2i']}' R3i={data['R3i']}'")
    print(f"  Inset area: {inset_area:.2f} sq ft")
    print(f'=== OUTLINE PATH ===')
    _pt_notes = [
        ("F0", "arc tangent"), ("F1", "arc tangent"),
        ("F2", "same E as F1"), ("F3", "arc tangent point"),
        ("F4", "5'8\" south of F5"), ("F5", "arc tangent"),
        ("F6", "arc tangent"), ("F7", "6.0' east of F6"),
        ("F8", "C7/C8 arc junction"), ("F9", "arc tangent"),
        ("F10", "arc E-W tangent"),
        ("F11", "180\u00b0 arc west end / arc N-S tangent"),
        ("F12", "line / 180\u00b0 arc tangent"),
        ("F13", f"{radii['R_a13']*12:.1f}\" arc / line tangent"),
        ("F14", f"{radii['R_a13']*12:.1f}\" arc tangent to N-S line"),
        ("F15", "arc C15, exits North"),
        ("F16", "arc C15, incoming tangent"),
        ("F17", "on PiX-Pi5 line"), ("F18", "arc C17 tangent"),
        ("F19", "arc C19 tangent"), ("F20", "arc junction"),
        ("F21", "south wall anchor"),
    ]
    for name, note in _pt_notes:
        print(f"  {name:<5s} ({pts[name][0]:.4f}, {pts[name][1]:.4f})  ({note})")
    print(f"  C20:  ({pts['C20'][0]:.4f}, {pts['C20'][1]:.4f})  (F20->F21 arc center)")
    print(f"  C19:  ({pts['C19'][0]:.4f}, {pts['C19'][1]:.4f})  (F19->F20 arc center)")
    print(f"  C17:  ({pts['C17'][0]:.4f}, {pts['C17'][1]:.4f})  (F17->F18 arc center)")
    print(f"  F18-F19 segment length = {abs(pts['F18'][0]-pts['F19'][0])*12:.1f}\"")
    print(f"  C15:  ({pts['C15'][0]:.4f}, {pts['C15'][1]:.4f})  (arc C15 center, R={radii['R_a15']:.4f}')")
    print(f"  C13:  ({pts['C13'][0]:.4f}, {pts['C13'][1]:.4f})  ({radii['R_a13']*12:.1f}\" arc center, R={radii['R_a13']:.4f}')")
    print(f"  Outline area: {outline_area:.2f} sq ft")

    outline_cfg = build_outline_cfg(outline_segs, pts, radii)

    lines = []
    lines.append(f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}">')
    lines.append(f'<rect width="{W}" height="{H}" fill="white"/>')
    lines.append('<defs>')
    lines.append('  <marker id="ah" markerWidth="8" markerHeight="6" refX="8" refY="3" orient="auto"><polygon points="0 0, 8 3, 0 6" fill="#333"/></marker>')
    lines.append(f'  <clipPath id="page"><rect width="{W}" height="{H}"/></clipPath>')
    lines.append('</defs>')
    lines.append(f'<text x="{W/2}" y="30" text-anchor="middle" font-family="Arial" font-size="14"'
                 f' font-weight="bold">Site Path \u2014 Outline</text>')

    render_layer(lines, outer_segs, pts_rot, outer_cfg, to_svg)
    render_layer(lines, inset_segs, pts_rot, inset_cfg, to_svg)
    render_layer(lines, outline_segs, pts, outline_cfg, to_svg)
    render_floorplan(lines, to_svg, pts, data["outer_poly"], data["inner_poly"],
                     data["inner_segs"], data["layout"])

    # Area label centered in outline
    centroid_names = [f"F{i}" for i in range(22)]
    cx_o = sum(pts[n][0] for n in centroid_names) / len(centroid_names)
    cy_o = sum(pts[n][1] for n in centroid_names) / len(centroid_names)
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
    svg_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "path_area.svg")
    with open(svg_path, "w", encoding="utf-8") as f:
        f.write(svg_content)

    print(f"\nSVG written to path_area.svg")
    print(f"Outer path area: {outer_area:.2f} sq ft (rendered at 20%)")
    print(f"Inset path area: {inset_area:.2f} sq ft (rendered at 20%)")
    print(f"Outline path area: {outline_area:.2f} sq ft (rendered at 100%)")
    print(f"Outline: F0->ArcC0->F1->F2->ArcC2->F3->ArcC3->F4->F5->ArcC5->F6->F7->ArcC7->F8->ArcC8->F9->F10->ArcC10->F11->ArcC11_180->F12->F13->ArcC13->F14->F15->ArcC15->F16->F17->ArcC17->F18->F19->ArcC19->F20->ArcC20->F21->F0")
