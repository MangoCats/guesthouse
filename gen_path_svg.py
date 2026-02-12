import math
from typing import NamedTuple
from survey import (
    Point, LineSeg, ArcSeg, Segment,
    poly_area, left_norm, off_pt, line_isect, arc_poly,
    circle_circle_isect, line_circle_isect_min_t_gt, line_circle_isect_min_abs_t,
    segment_polyline, path_polygon, arc_sweep_deg,
    brg_dist, fmt_brg, fmt_dist,
    compute_traverse, compute_three_arc,
    compute_inner_walls, compute_interior_layout,
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

# ============================================================
# SVG Constants
# ============================================================
W, H = 792, 612
_s = (368.79 - 151.26) / 18.66

# ============================================================
# Section 4: Generic SVG Renderer
# ============================================================
def render_layer(lines: list, segments: list[Segment], pts: dict, cfg: LayerConfig):
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
# Section 5: Base Geometry + Outer Path
# ============================================================
# All coordinates are P3-based: P3 = (0, 0).
pts: dict[str, Point] = {}

pts, _p3_trav = compute_traverse()

# SVG transform (P3-based coordinates)
_p3_svg_x = 368.79 + _p3_trav[0] * _s
_p3_svg_y = 124.12 - _p3_trav[1] * _s
def to_svg(e, n):
    return (_p3_svg_x + e*_s, _p3_svg_y - n*_s)

_arc_info = compute_three_arc(pts)
R1, R2, R3 = _arc_info["R1"], _arc_info["R2"], _arc_info["R3"]
nE, nN = _arc_info["nE"], _arc_info["nN"]

# Outer path
outer_segs = [
    LineSeg("POB", "P2"), LineSeg("P2", "P3"), LineSeg("P3", "T3"),
    ArcSeg("T3", "PX", "C3", R3, "CW", 60),
    LineSeg("PX", "P4"), LineSeg("P4", "P5"), LineSeg("P5", "T1"),
    ArcSeg("T1", "PA", "C1", R1, "CW", 60),
    ArcSeg("PA", "T2", "C2", R2, "CW", 60),
    LineSeg("T2", "POB"),
]
outer_area = poly_area(path_polygon(outer_segs, pts))

# ============================================================
# Section 6: Inset Path (6" inside)
# ============================================================
delta = 0.5
R1i, R2i, R3i = R1+delta, R2+delta, R3+delta

d_e1 = (pts["P2"][0]-pts["POB"][0], pts["P2"][1]-pts["POB"][1])
d_e2 = (pts["P3"][0]-pts["P2"][0], pts["P3"][1]-pts["P2"][1])
d_e3 = (pts["T3"][0]-pts["P3"][0], pts["T3"][1]-pts["P3"][1])
d_e5 = (pts["P4"][0]-pts["PX"][0], pts["P4"][1]-pts["PX"][1])
d_e6 = (pts["P5"][0]-pts["P4"][0], pts["P5"][1]-pts["P4"][1])
d_e7 = (pts["T1"][0]-pts["P5"][0], pts["T1"][1]-pts["P5"][1])
d_e10 = (pts["POB"][0]-pts["T2"][0], pts["POB"][1]-pts["T2"][1])

ln1 = left_norm(pts["POB"], pts["P2"])
ln2 = left_norm(pts["P2"], pts["P3"])
ln3 = left_norm(pts["P3"], pts["T3"])
ln5 = left_norm(pts["PX"], pts["P4"])
ln6 = left_norm(pts["P4"], pts["P5"])
ln7 = left_norm(pts["P5"], pts["T1"])
ln10 = left_norm(pts["T2"], pts["POB"])

o1 = off_pt(pts["POB"], ln1, delta)
o2 = off_pt(pts["P2"], ln2, delta)
o3 = off_pt(pts["P3"], ln3, delta)
o5 = off_pt(pts["PX"], ln5, delta)
o6 = off_pt(pts["P4"], ln6, delta)
o7 = off_pt(pts["P5"], ln7, delta)
o10 = off_pt(pts["T2"], ln10, delta)

pts["PiOB"] = line_isect(o10, d_e10, o1, d_e1)
pts["Pi2"] = line_isect(o1, d_e1, o2, d_e2)
pts["Pi3"] = line_isect(o2, d_e2, o3, d_e3)
pts["Pi4"] = off_pt(pts["P4"], ln5, delta)  # collinear PX-P4-P5
pts["Pi5"] = line_isect(o6, d_e6, o7, d_e7)

pts["Ti3"] = (pts["C3"][0], pts["P3"][1] + delta)
pts["Ti1"] = (pts["T1"][0] - delta*nE, pts["T1"][1] - delta*nN)
pts["Ti2"] = (pts["T2"][0] - delta*nE, pts["T2"][1] - delta*nN)

pts["PiX"] = line_circle_isect_min_abs_t(o5, d_e5, pts["C3"], R3i)
pts["Ai2"] = circle_circle_isect(pts["C1"], R1i, pts["C2"], R2i, near=pts["PA"])

inset_segs = [
    LineSeg("PiOB", "Pi2"), LineSeg("Pi2", "Pi3"), LineSeg("Pi3", "Ti3"),
    ArcSeg("Ti3", "PiX", "C3", R3i, "CW", 60),
    LineSeg("PiX", "Pi4"), LineSeg("Pi4", "Pi5"), LineSeg("Pi5", "Ti1"),
    ArcSeg("Ti1", "Ai2", "C1", R1i, "CW", 60),
    ArcSeg("Ai2", "Ti2", "C2", R2i, "CW", 60),
    LineSeg("Ti2", "PiOB"),
]
inset_area = poly_area(path_polygon(inset_segs, pts))

# ============================================================
# Section 7: Outline Path (inset + fillets at Po3 and Po2)
# ============================================================
pts["Po2"]  = pts["Pi2"]
pts["O15"]  = pts["Ti3"]   # was To3
pts["To2"]  = pts["Ti2"]

R_fillet = 10.0 / 12.0
pts["Cf"] = (pts["Pi3"][0] + R_fillet, pts["Pi3"][1] + R_fillet)
pts["O1"] = (pts["Pi3"][0], pts["Cf"][1])
pts["O0"] = (pts["Cf"][0], pts["Pi3"][1])

# --- Fillet at Po2 corner (R = 28", exact 90°) ---
R_fillet2 = 28.0 / 12.0
corner2_N = pts["O0"][1] + 26 - 2.0/12  # O3-O4 line is 25'10" north of O0
# C2 and C3 shift 1' east
_nw_shift = 1.0
pts["Cf2"] = (pts["Po2"][0] + _nw_shift + R_fillet2, corner2_N - R_fillet2)
pts["O2"] = (pts["Po2"][0] + _nw_shift, pts["Cf2"][1])
pts["O3"] = (pts["Cf2"][0], corner2_N)
# C1c: 5'8" south of C2, same easting
pts["O1c"] = (pts["O2"][0], pts["O2"][1] - (5.0 + 8.0/12) + 28.0/12)
# Arc C1c->C1b: R=20", center due West of C1c (CW)
R_1c = 20.0 / 12.0
pts["Cc1"] = (pts["O1c"][0] - R_1c, pts["O1c"][1])
# Arc C1b->C1a: R=28", center due East of C1a (CCW), tangent to C1c-C1b arc
R_1a = 28.0 / 12.0
# External tangency: |Cc1 - Cc2| = R_1c + R_1a, Cc2 = (O1_E + R_1a, O1a_N)
_cc2_E = pts["O1"][0] + R_1a
_dE_cc = _cc2_E - pts["Cc1"][0]
_dN_cc = -math.sqrt((R_1c + R_1a)**2 - _dE_cc**2)  # south
pts["O1a"] = (pts["O1"][0], pts["Cc1"][1] + _dN_cc)
pts["Cc2"] = (_cc2_E, pts["O1a"][1])
# O1b: tangent point on Cc1->Cc2 line
_f1b = R_1c / (R_1c + R_1a)
pts["O1b"] = (pts["Cc1"][0] + _f1b * (pts["Cc2"][0] - pts["Cc1"][0]),
              pts["Cc1"][1] + _f1b * (pts["Cc2"][1] - pts["Cc1"][1]))
# O4: 6.0' east of O3, adjusted 1' west to undo C2/C3 shift
pts["O4"] = (pts["O3"][0] + 5.5 + 6.0/12 - _nw_shift, pts["O3"][1])
# Turn 1 at O4: 90° CW east->south, R=28"
R_turn1 = 28.0 / 12.0
pts["Ct1"] = (pts["O4"][0], pts["O4"][1] - R_turn1)
pts["O5"] = (pts["Ct1"][0] + R_turn1, pts["Ct1"][1])
# Turn 2 at O5: 90° CCW south->east, R=2"
R_turn2 = 2.0 / 12.0
pts["Ct2"] = (pts["O5"][0] + R_turn2, pts["O5"][1])
pts["O6"] = (pts["Ct2"][0], pts["Ct2"][1] - R_turn2)
# Turn 3: 180° CCW arc + 90° CW fillet at C6a (R=2")
R_turn3 = 28.0 / 12.0
R_6a = 2.0 / 12.0  # 90° fillet at C6a
west_E = pts["C1"][0] - R1i  # westernmost E on Arc 1o
# O6a, Ct6a, O7, Ct3 eastings determined below from bearing constraint

# --- Fillet at Po5 corner (R = 10", exits North) ---
R_f_po5 = 28.0 / 12.0
d_in_po5 = (pts["Pi5"][0] - pts["PiX"][0], pts["Pi5"][1] - pts["PiX"][1])
L_in = math.sqrt(d_in_po5[0]**2 + d_in_po5[1]**2)
d_in_u = (d_in_po5[0]/L_in, d_in_po5[1]/L_in)
# C11 E-coordinate: 9' 2" east of east closet wall (w5_e)
_w5_e = pts["O1"][0] + 8.0/12 + 0.5 + 35.0/12 + 3.0 + 2.0 + (3+30+4+140+4+30+3)/12.0
C11_E = _w5_e + 9.0 + 1.0/12
# Outgoing direction is North; Cf4 is R_f_po5 west of C11
ln_in_po5 = left_norm(pts["PiX"], pts["Pi5"])
o_in_po5 = off_pt(pts["PiX"], ln_in_po5, R_f_po5)
t_cf4 = (C11_E - R_f_po5 - o_in_po5[0]) / d_in_u[0]
pts["Cf4"] = (C11_E - R_f_po5, o_in_po5[1] + t_cf4 * d_in_u[1])
pts["O11"] = (C11_E, pts["Cf4"][1])  # directly east of Cf4
# O12: perpendicular foot from Cf4 onto incoming line
t_in = ((pts["Cf4"][0]-pts["PiX"][0])*d_in_po5[0] + (pts["Cf4"][1]-pts["PiX"][1])*d_in_po5[1]) \
       / (d_in_po5[0]**2 + d_in_po5[1]**2)
pts["O12"] = (pts["PiX"][0] + t_in*d_in_po5[0], pts["PiX"][1] + t_in*d_in_po5[1])

# --- Arc O10a-O9: bearing O9->O8 = 345°, O9_N = IW1_north + 24" ---
# C10a: centered 2" north of IW1 north face
# IW1 north = inner south wall (O0 + 8" wall) + 11'6" + 6" IW1 thickness
_iw1_n_face = pts["O0"][1] + 8.0/12 + 11.5 + 6.0/12
_O10a_N = _iw1_n_face + 2.0/12
_target_O9_N = _iw1_n_face + 24.0/12
_dN_t4 = _target_O9_N - _O10a_N
# Bearing 345° => tangent direction (-sin15°, cos15°), normal n = (cos15°, sin15°)
_brg_off = math.radians(360.0 - 345.0)  # 15°
_nx_t = math.cos(_brg_off)
_ny_t = math.sin(_brg_off)
# R_t4 from O9_N constraint: O9_N = _O10a_N + R_t4 * ny
R_t4 = _dN_t4 / _ny_t
pts["Ct4"] = (C11_E - R_t4, _O10a_N)
pts["O10a"] = (C11_E, _O10a_N)
# Ct3 easting from tangent constraint: (Ct3 - Ct4) · n = R_t4 - R_turn3
_Ct3_N = pts["O6"][1] + R_6a  # O7_N = Ct3_N (independent of _corner_E)
_Ct4_E, _Ct4_N = pts["Ct4"]
_Ct3_E = _Ct4_E + (R_t4 - R_turn3 - (_Ct3_N - _Ct4_N) * _ny_t) / _nx_t
# Derive corner easting and dependent points (C6a, C7, Ct3)
_corner_E = _Ct3_E - R_turn3
pts["O6a"] = (_corner_E - R_6a, pts["O6"][1])
pts["Ct6a"] = (_corner_E - R_6a, pts["O6"][1] + R_6a)
pts["O7"] = (_corner_E, pts["O6"][1] + R_6a)
pts["Ct3"] = (_Ct3_E, _Ct3_N)
# O9, O8: tangent points (shared outward normal for external tangent)
pts["O9"] = (pts["Ct4"][0] + R_t4 * _nx_t, pts["Ct4"][1] + R_t4 * _ny_t)
pts["O8"] = (pts["Ct3"][0] + R_turn3 * _nx_t, pts["Ct3"][1] + R_turn3 * _ny_t)

# --- C13-C14-C15 wall geometry (wall at -6" Northing) ---
R_wall = 28.0 / 12.0          # 28" connecting arc radius (Cw3)
R_w1   = 20.0 / 12.0          # 20" arc radius (Cw1: C15->C14)
R_w2   = 28.0 / 12.0          # 28" arc radius (Cw2: C14->C13b)
wall_south_N = -6.0 / 12.0    # south face of wall at -6" Northing
# Tangency distance between Cw1 and Cw2 arc centers
dN_c = (wall_south_N + R_w2) - (pts["O15"][1] - R_w1)
dE_c = math.sqrt((R_w1 + R_w2)**2 - dN_c**2)
# Align C13b with east side of king bed
# bed_e = inner_W1_E + 20.5' (dryer+counter+closet2+W1S+half_bedroom+half_bed)
_bed_e_align = pts["O1"][0] + 8.0/12 + 20.5
pts["O15"] = (_bed_e_align - dE_c, pts["O15"][1])
# O13 on line from O12 at bearing 60° (O13->O12), tangent to Cw3 arc
_brg_13 = math.radians(60.0)
C13_N = wall_south_N + R_wall * (1.0 - math.sin(_brg_13))
_t_13 = (pts["O12"][1] - C13_N) / math.cos(_brg_13)
pts["O13"] = (pts["O12"][0] - _t_13 * math.sin(_brg_13), C13_N)
# C13a->C13 arc center (Cw3) and tangent point C13a
pts["Cw3"] = (pts["O13"][0] - R_wall * math.cos(_brg_13), wall_south_N + R_wall)
pts["O13a"] = (pts["Cw3"][0], wall_south_N)
# C15->C14 arc center (Cw1): center south of O15
pts["Cw1"] = (pts["O15"][0], pts["O15"][1] - R_w1)
# C14->C13b arc center (Cw2): from external tangency |Cw1-Cw2| = R_w1+R_w2
C13b_E = pts["O15"][0] + dE_c
pts["O13b"] = (C13b_E, wall_south_N)
pts["Cw2"] = (C13b_E, wall_south_N + R_w2)
# C14 = tangent point between the two arcs (R_w1 from Cw1 along Cw1->Cw2)
_f_w = R_w1 / (R_w1 + R_w2)
pts["O14"] = (pts["Cw1"][0] + _f_w * (pts["Cw2"][0] - pts["Cw1"][0]),
              pts["Cw1"][1] + _f_w * (pts["Cw2"][1] - pts["Cw1"][1]))

outline_segs = [
    LineSeg("O2", "O1c"),                                    # 0
    ArcSeg("O1c", "O1b", "Cc1", R_1c, "CW", 20),          # 1: 20" arc
    ArcSeg("O1b", "O1a", "Cc2", R_1a, "CCW", 20),         # 2: 28" arc
    LineSeg("O1a", "O1"),                                    # 3
    ArcSeg("O1", "O0", "Cf", R_fillet, "CCW", 20),         # 4
    LineSeg("O0", "O15"),                                    # 5
    ArcSeg("O15", "O14", "Cw1", R_w1, "CW", 60),           # 6: C15->C14
    ArcSeg("O14", "O13b", "Cw2", R_w2, "CCW", 60),         # 7: C14->C13b
    LineSeg("O13b", "O13a"),                                 # 8: wall segment
    ArcSeg("O13a", "O13", "Cw3", R_wall, "CCW", 20),       # 9: C13a->C13
    LineSeg("O13", "O12"),                                   # 10
    ArcSeg("O12", "O11", "Cf4", R_f_po5, "CCW", 20),       # 11
    LineSeg("O11", "O10a"),                                  # 12
    ArcSeg("O10a", "O9", "Ct4", R_t4, "CCW", 60),          # 13: 28" arc
    LineSeg("O9", "O8"),                                     # 14
    ArcSeg("O8", "O7", "Ct3", R_turn3, "CCW", 60),         # 15: 180° arc
    ArcSeg("O7", "O6a", "Ct6a", R_6a, "CW", 20),          # 16: 90° fillet at C6a
    LineSeg("O6a", "O6"),                                    # 17
    ArcSeg("O6", "O5", "Ct2", R_turn2, "CW", 20),          # 18
    ArcSeg("O5", "O4", "Ct1", R_turn1, "CCW", 20),         # 19
    LineSeg("O4", "O3"),                                     # 20
    ArcSeg("O3", "O2", "Cf2", R_fillet2, "CCW", 20),       # 21
]
outline_area = poly_area(path_polygon(outline_segs, pts))

# ============================================================
# Section 7b: Floorplan Inner Walls + Interior Elements
# ============================================================
_wall_t = 8.0 / 12.0
_radii = {
    "R_fillet": R_fillet, "R_w1": R_w1, "R_w2": R_w2, "R_wall": R_wall,
    "R_f_po5": R_f_po5, "R_turn3": R_turn3, "R_turn2": R_turn2,
    "R_turn1": R_turn1, "R_fillet2": R_fillet2, "R_t4": R_t4, "R_6a": R_6a,
    "R_1c": R_1c, "R_1a": R_1a,
}
_fp_inner_segs = compute_inner_walls(outline_segs, pts, _wall_t, _radii)
_fp_outer_poly = path_polygon(outline_segs, pts)
_fp_inner_poly = path_polygon(_fp_inner_segs, pts)
_fp_inner_area = poly_area(_fp_inner_poly)

_layout = compute_interior_layout(pts, _fp_inner_poly)
_iw1 = _layout["iw1"]
_iw2_w = _layout["iw2_w"]; _iw2_e = _layout["iw2_e"]
_iw2_s = _layout["iw2_s"]; _iw2_n = _layout["iw2_n"]
_dryer_w = _layout["dryer_w"]; _dryer_s = _layout["dryer_s"]
_dryer_e = _layout["dryer_e"]; _dryer_n = _layout["dryer_n"]
_washer_w = _layout["washer_w"]; _washer_s = _layout["washer_s"]
_washer_e = _layout["washer_e"]; _washer_n = _layout["washer_n"]
_ctr_w = _layout["ctr_w"]; _ctr_e = _layout["ctr_e"]
_ctr_s = _layout["ctr_s"]; _ctr_n = _layout["ctr_n"]
_ctr_nw_r = _layout["ctr_nw_r"]
_iwt3 = _layout["iwt3"]
_w8 = _layout["wall8"]
_iw3_w = _layout["iw3_w"]; _iw3_e = _layout["iw3_e"]
_iw3_s = _layout["iw3_s"]; _iw3_n = _layout["iw3_n"]
_iw4_w = _layout["iw4_w"]; _iw4_e = _layout["iw4_e"]
_wall_south_n = _layout["wall_south_n"]
_iw1_s = _layout["iw1_s"]
_iw4_sw = _wall_south_n; _iw4_se = _wall_south_n
_cl1_top = _layout["cl1_top"]
_w5_w = _layout["w5_w"]; _w5_e_inner = _layout["w5_e"]
_w5_sw = _wall_south_n; _w5_se = _wall_south_n
_bed_w = _layout["bed_w"]; _bed_e = _layout["bed_e"]
_bed_s = _layout["bed_s"]; _bed_n = _layout["bed_n"]

def render_floorplan(lines):
    lines.append('<g opacity="0.5">')
    # Wall band (outer - inner hole, evenodd fill)
    od = "M "+" L ".join(f"{to_svg(*p)[0]:.1f},{to_svg(*p)[1]:.1f}" for p in _fp_outer_poly)+" Z"
    id_ = "M "+" L ".join(f"{to_svg(*p)[0]:.1f},{to_svg(*p)[1]:.1f}" for p in reversed(_fp_inner_poly))+" Z"
    lines.append(f'<path d="{od} {id_}" fill="rgba(160,160,160,0.5)" fill-rule="evenodd" stroke="none"/>')
    # Inner wall strokes
    for seg in _fp_inner_segs:
        if isinstance(seg, LineSeg):
            s1,s2 = to_svg(*pts[seg.start]),to_svg(*pts[seg.end])
            lines.append(f'<line x1="{s1[0]:.1f}" y1="{s1[1]:.1f}" x2="{s2[0]:.1f}" y2="{s2[1]:.1f}" stroke="#666" stroke-width="1.0"/>')
        else:
            pl = segment_polyline(seg, pts)
            sp = " ".join(f"{to_svg(e,n)[0]:.1f},{to_svg(e,n)[1]:.1f}" for e,n in pl)
            lines.append(f'<polyline points="{sp}" fill="none" stroke="#666" stroke-width="1.0" stroke-linecap="round"/>')
    # IW1
    svg = " ".join(f"{to_svg(*p)[0]:.1f},{to_svg(*p)[1]:.1f}" for p in _iw1)
    lines.append(f'<polygon points="{svg}" fill="rgba(160,160,160,0.5)" stroke="none"/>')
    for a,b in [(_iw1[0],_iw1[1]),(_iw1[3],_iw1[2])]:
        s1,s2 = to_svg(*a),to_svg(*b)
        lines.append(f'<line x1="{s1[0]:.1f}" y1="{s1[1]:.1f}" x2="{s2[0]:.1f}" y2="{s2[1]:.1f}" stroke="#666" stroke-width="1.0"/>')
    # IW2
    iw2 = [(_iw2_w,_iw2_s),(_iw2_e,_iw2_s),(_iw2_e,_iw2_n),(_iw2_w,_iw2_n)]
    svg = " ".join(f"{to_svg(*p)[0]:.1f},{to_svg(*p)[1]:.1f}" for p in iw2)
    lines.append(f'<polygon points="{svg}" fill="rgba(160,160,160,0.5)" stroke="none"/>')
    for ev in [_iw2_w,_iw2_e]:
        s1,s2 = to_svg(ev,_iw2_s),to_svg(ev,_iw2_n)
        lines.append(f'<line x1="{s1[0]:.1f}" y1="{s1[1]:.1f}" x2="{s2[0]:.1f}" y2="{s2[1]:.1f}" stroke="#666" stroke-width="1.0"/>')
    # Wall 8 L-shape
    svg = " ".join(f"{to_svg(*p)[0]:.1f},{to_svg(*p)[1]:.1f}" for p in _w8)
    lines.append(f'<polygon points="{svg}" fill="rgba(160,160,160,0.5)" stroke="#666" stroke-width="0.8"/>')
    # IW3 (west bedroom wall)
    w1sp = [(_iw3_w,_iw3_s),(_iw3_e,_iw3_s),(_iw3_e,_iw3_n),(_iw3_w,_iw3_n)]
    svg = " ".join(f"{to_svg(*p)[0]:.1f},{to_svg(*p)[1]:.1f}" for p in w1sp)
    lines.append(f'<polygon points="{svg}" fill="rgba(160,160,160,0.5)" stroke="#666" stroke-width="0.8"/>')
    # IW4 (east bedroom wall)
    w6p = [(_iw4_w,_iw4_sw),(_iw4_e,_iw4_se),(_iw4_e,_iw1_s),(_iw4_w,_iw1_s)]
    svg = " ".join(f"{to_svg(*p)[0]:.1f},{to_svg(*p)[1]:.1f}" for p in w6p)
    lines.append(f'<polygon points="{svg}" fill="rgba(160,160,160,0.5)" stroke="#666" stroke-width="0.8"/>')
    # Wall 5 L-shape
    w5p = [(_iw4_e,_cl1_top+_iwt3),(_w5_e_inner,_cl1_top+_iwt3),(_w5_e_inner,_w5_se),
           (_w5_w,_w5_sw),(_w5_w,_cl1_top),(_iw4_e,_cl1_top)]
    svg = " ".join(f"{to_svg(*p)[0]:.1f},{to_svg(*p)[1]:.1f}" for p in w5p)
    lines.append(f'<polygon points="{svg}" fill="rgba(160,160,160,0.5)" stroke="#666" stroke-width="0.8"/>')
    # Appliances
    for lbl,sw_e,sw_n,ne_e,ne_n in [("DRYER",_dryer_w,_dryer_s,_dryer_e,_dryer_n),
                                      ("WASHER",_washer_w,_washer_s,_washer_e,_washer_n)]:
        s1,s2 = to_svg(sw_e,ne_n),to_svg(ne_e,sw_n); w=s2[0]-s1[0]; h=s2[1]-s1[1]
        lines.append(f'<rect x="{s1[0]:.1f}" y="{s1[1]:.1f}" width="{w:.1f}" height="{h:.1f}" fill="rgba(100,150,200,0.3)" stroke="#4682B4" stroke-width="0.8"/>')
        cx,cy = (s1[0]+s2[0])/2,(s1[1]+s2[1])/2
        lines.append(f'<text x="{cx:.1f}" y="{cy+3:.1f}" text-anchor="middle" font-family="Arial" font-size="7" fill="#4682B4">{lbl}</text>')
    # Counter with rounded NW corner
    csw,cse,cne = to_svg(_ctr_w,_ctr_s),to_svg(_ctr_e,_ctr_s),to_svg(_ctr_e,_ctr_n)
    cnas,cnae = to_svg(_ctr_w+_ctr_nw_r,_ctr_n),to_svg(_ctr_w,_ctr_n-_ctr_nw_r)
    rsv = abs(cnas[0]-to_svg(_ctr_w,_ctr_n)[0])
    cd = (f'M {csw[0]:.1f},{csw[1]:.1f} L {cse[0]:.1f},{cse[1]:.1f} L {cne[0]:.1f},{cne[1]:.1f} '
          f'L {cnas[0]:.1f},{cnas[1]:.1f} A {rsv:.1f} {rsv:.1f} 0 0 0 {cnae[0]:.1f},{cnae[1]:.1f} Z')
    lines.append(f'<path d="{cd}" fill="rgba(100,150,200,0.3)" stroke="#4682B4" stroke-width="0.8"/>')
    ccx,ccy = (csw[0]+cse[0])/2,(csw[1]+cne[1])/2
    lines.append(f'<text x="{ccx:.1f}" y="{ccy:.1f}" text-anchor="middle" font-family="Arial" font-size="7" fill="#4682B4" letter-spacing="0.5" transform="rotate(-90,{ccx:.1f},{ccy:.1f})">COUNTER</text>')
    # King Bed
    bs,be = to_svg(_bed_w,_bed_n),to_svg(_bed_e,_bed_s); bw,bh = be[0]-bs[0],be[1]-bs[1]
    lines.append(f'<rect x="{bs[0]:.1f}" y="{bs[1]:.1f}" width="{bw:.1f}" height="{bh:.1f}" fill="rgba(100,150,200,0.3)" stroke="#4682B4" stroke-width="0.8"/>')
    bcx,bly = (bs[0]+be[0])/2,bs[1]+0.765*bh
    lines.append(f'<text x="{bcx:.1f}" y="{bly+3:.1f}" text-anchor="middle" font-family="Arial" font-size="7" fill="#4682B4">KING BED</text>')
    # Room labels
    cx,cy = to_svg((_ctr_e+_iwt3+_iw3_w)/2,(_ctr_s+_ctr_n)/2)
    lines.append(f'<text x="{cx:.1f}" y="{cy+3:.1f}" text-anchor="middle" font-family="Arial" font-size="7" fill="#666" transform="rotate(-90,{cx:.1f},{cy+3:.1f})">CLOSET</text>')
    bx,by = to_svg((_iw3_e+_iw4_w)/2,(_ctr_s+_iw1_s)/2)
    lines.append(f'<text x="{bx:.1f}" y="{by+3:.1f}" text-anchor="middle" font-family="Arial" font-size="8" fill="#666">BEDROOM</text>')
    cx,cy = to_svg((_iw4_e+_w5_w)/2,(_ctr_s+_cl1_top)/2)
    lines.append(f'<text x="{cx:.1f}" y="{cy+3:.1f}" text-anchor="middle" font-family="Arial" font-size="7" fill="#666" transform="rotate(-90,{cx:.1f},{cy+3:.1f})">CLOSET</text>')
    lines.append('</g>')

# ============================================================
# Section 8: Style Configurations
# ============================================================

# Compute sweep angles for outline arc labels
sweep_1c = arc_sweep_deg(outline_segs[1], pts)   # O1c->O1b (Cc1)
sweep_1a = arc_sweep_deg(outline_segs[2], pts)   # O1b->O1a (Cc2)
sweep_f = arc_sweep_deg(outline_segs[4], pts)    # O1->O0
sweep_w1 = arc_sweep_deg(outline_segs[6], pts)   # O15->O14 (Cw1)
sweep_w2 = arc_sweep_deg(outline_segs[7], pts)   # O14->O13b (Cw2)
sweep_w3 = arc_sweep_deg(outline_segs[9], pts)   # O13a->O13 (Cw3)
sweep_f4 = arc_sweep_deg(outline_segs[11], pts)  # O12->O11 (fillet Po5)
sweep_t4 = arc_sweep_deg(outline_segs[13], pts)  # O10a->O9 (28" arc)
sweep_t3 = arc_sweep_deg(outline_segs[15], pts)  # O8->O7
sweep_6a = arc_sweep_deg(outline_segs[16], pts)  # O7->O6a (fillet)
sweep_t2 = arc_sweep_deg(outline_segs[18], pts)  # O6->O5
sweep_t1 = arc_sweep_deg(outline_segs[19], pts)  # O5->O4
sweep_f2 = arc_sweep_deg(outline_segs[21], pts)  # O3->O2

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

outline_cfg = LayerConfig(
    opacity=1.0, fill_color="rgba(200,230,255,0.25)",
    line_stroke="#333", line_width=2.0,
    arc_styles={
        ("O1c","O1b"):  ("#333", 2.0),
        ("O1b","O1a"):  ("#333", 2.0),
        ("O10a","O9"):  ("#333", 2.0),
        ("O15","O14"):  ("#2E7D32", 2.5),
        ("O14","O13b"): ("#2E7D32", 2.5),
        ("O13a","O13"): ("#2E7D32", 2.5),
        ("O12","O11"):  ("#333", 2.0),
        ("O1","O0"):    ("#333", 2.0),
        ("O3","O2"):    ("#333", 2.0),
        ("O6","O5"):    ("#333", 2.0),
        ("O5","O4"):    ("#333", 2.0),
        ("O8","O7"):    ("#333", 2.0),
        ("O7","O6a"):   ("#333", 2.0),
    },
    vertex_styles={
        "O2":   VertexStyle("O2",   "end",   -8, -6,  "#d32f2f", 3.5, 10),
        "O1c":  VertexStyle("O1c",  "end",   -8, -4,  "#d32f2f", 3.5, 10),
        "O1b":  VertexStyle("O1b",  "end",  -10, 10,  "#d32f2f", 3.5, 10),
        "O1a":  VertexStyle("O1a",  "end",   -8, -4,  "#d32f2f", 3.5, 10),
        "O1":   VertexStyle("O1",   "end",   -8, -4,  "#d32f2f", 3.5, 10),
        "O0":   VertexStyle("O0",   "end",   -8, 10,  "#d32f2f", 3.5, 10),
        "O15":  VertexStyle("O15",  "middle", 0, -12, "#2E7D32", 3.5, 10),
        "O14":  VertexStyle("O14",  "start",  8, -6,  "#2E7D32", 3.5, 10),
        "O13b": VertexStyle("O13b", "start",  8, 12,  "#2E7D32", 3.5, 10),
        "O13a": VertexStyle("O13a", "start",  8, 12,  "#2E7D32", 3.5, 10),
        "O13":  VertexStyle("O13",  "start", 10, 10,  "#2E7D32", 3.5, 10),
        "O12":  VertexStyle("O12",  "start",  8,  4,  "#d32f2f", 3.5, 10),
        "O11":  VertexStyle("O11",  "start",  8,  4,  "#d32f2f", 3.5, 10),
        "O10a": VertexStyle("O10a", "start", 10, -4,  "#d32f2f", 3.5, 10),
        "O9":   VertexStyle("O9",   "start",  8,  8,  "#d32f2f", 3.5, 10),
        "O8":   VertexStyle("O8",   "start",  8, 12,  "#d32f2f", 3.5, 10),
        "O7":   VertexStyle("O7",   "start",  8, 12,  "#d32f2f", 3.5, 10),
        "O6a":  VertexStyle("O6a",  "end",   -8, -6,  "#d32f2f", 3.5, 10),
        "O6":   VertexStyle("O6",   "middle", 0, 14,  "#d32f2f", 3.5, 10),
        "O5":   VertexStyle("O5",   "end",   -8, -4,  "#d32f2f", 3.5, 10),
        "O4":   VertexStyle("O4",   "middle", 0, -12, "#d32f2f", 3.5, 10),
        "O3":   VertexStyle("O3",   "end",  -10, -6,  "#d32f2f", 3.5, 10),
    },
    brg_dist_labels={
        ("O2","O1c"): BrgDistLabel(-18),
        ("O1a","O1"): BrgDistLabel(-18),
        ("O0","O15"): BrgDistLabel(18),
        ("O13b","O13a"): BrgDistLabel(16),
        ("O13","O12"): BrgDistLabel(-16),
        ("O11","O10a"): BrgDistLabel(16),
        ("O9","O8"): BrgDistLabel(16),
        ("O6a","O6"): BrgDistLabel(16),
        ("O4","O3"): BrgDistLabel(16),
    },
    arc_labels={
        ("O8","O7"): ArcLabel("Turn R=28\"",
            f"{sweep_t3:.1f}\u00b0 CCW", "start", 12, 0, 11, "#333"),
        ("O7","O6a"): ArcLabel("Fillet R=2\"",
            f"{sweep_6a:.1f}\u00b0", "end", -10, -10, 11, "#333"),
        ("O15","O14"): ArcLabel("Wall R=20\"",
            f"{sweep_w1:.1f}\u00b0 CW", "start", 12, 4, 11, "#2E7D32"),
        ("O14","O13b"): ArcLabel("Wall R=28\"",
            f"{sweep_w2:.1f}\u00b0", "start", 12, -10, 11, "#2E7D32"),
        ("O13a","O13"): ArcLabel("Wall R=28\"",
            f"{sweep_w3:.1f}\u00b0", "end", -10, -10, 11, "#2E7D32"),
        ("O10a","O9"): ArcLabel(f"Turn R={R_t4*12:.0f}\"",
            f"{sweep_t4:.1f}\u00b0", "start", 12, 0, 11, "#333"),
        ("O12","O11"): ArcLabel("Fillet R=28\"",
            f"{sweep_f4:.1f}\u00b0", "start", 10, -10, 11, "#333"),
        ("O1c","O1b"): ArcLabel("Arc R=20\"",
            f"{sweep_1c:.1f}\u00b0 CW", "end", -10, -10, 11, "#333"),
        ("O1b","O1a"): ArcLabel("Arc R=28\"",
            f"{sweep_1a:.1f}\u00b0 CCW", "start", 12, 0, 11, "#333"),
        ("O1","O0"): ArcLabel("Fillet R=10\"",
            f"{sweep_f:.1f}\u00b0", "end", -10, 14, 11, "#333"),
        ("O3","O2"): ArcLabel("Fillet R=28\"",
            f"{sweep_f2:.1f}\u00b0", "end", -10, -14, 11, "#333"),
        ("O6","O5"): ArcLabel("Turn R=2\"",
            f"{sweep_t2:.1f}\u00b0", "end", -10, 14, 11, "#333"),
        ("O5","O4"): ArcLabel("Turn R=28\"",
            f"{sweep_t1:.1f}\u00b0", "start", 12, 0, 11, "#333"),
    },
    center_marks=[
        CenterMark("Cc1", "O1c", "#333"), CenterMark("Cc2", "O1a", "#333"),
        CenterMark("Cw1", "O15", "#2E7D32"), CenterMark("Cf", "O1", "#333"),
        CenterMark("Cf2", "O3", "#333"), CenterMark("Cw2", "O14", "#2E7D32"),
        CenterMark("Cw3", "O13a", "#2E7D32"), CenterMark("Cf4", "O12", "#333"),
        CenterMark("Ct1", "O4", "#333"), CenterMark("Ct2", "O5", "#333"),
        CenterMark("Ct3", "O7", "#333"), CenterMark("Ct6a", "O6a", "#333"),
        CenterMark("Ct4", "O10a", "#333"),
    ],
    traverse_pts=None, traverse_stroke=None,
    brg_decimal=True,
)

# ============================================================
# Section 9: SVG Assembly + Output (only when run directly)
# ============================================================
if __name__ == "__main__":
    print(f'=== INSET PATH (6" inside) ===')
    print(f"  delta={delta}' R1i={R1i}' R2i={R2i}' R3i={R3i}'")
    print(f"  Inset area: {inset_area:.2f} sq ft")
    print(f'=== OUTLINE PATH (wall arcs R=28", fillets R=10" & R=28") ===')
    print(f"  O2:   ({pts['O2'][0]:.4f}, {pts['O2'][1]:.4f})  (fillet2 tangent)")
    print(f"  O1c:  ({pts['O1c'][0]:.4f}, {pts['O1c'][1]:.4f})  (5'8\" south of O2)")
    print(f"  O1b:  ({pts['O1b'][0]:.4f}, {pts['O1b'][1]:.4f})  (arc tangent point)")
    print(f"  O1a:  ({pts['O1a'][0]:.4f}, {pts['O1a'][1]:.4f})  (same E as O1)")
    print(f"  O1:   ({pts['O1'][0]:.4f}, {pts['O1'][1]:.4f})  (fillet tangent)")
    print(f"  O0:   ({pts['O0'][0]:.4f}, {pts['O0'][1]:.4f})  (fillet tangent)")
    print(f"  O15:  ({pts['O15'][0]:.4f}, {pts['O15'][1]:.4f})  (C15, was To3)")
    print(f"  O14:  ({pts['O14'][0]:.4f}, {pts['O14'][1]:.4f})  (C14, arc junction)")
    print(f"  O13b: ({pts['O13b'][0]:.4f}, {pts['O13b'][1]:.4f})  (C13b, wall west end)")
    print(f"  O13a: ({pts['O13a'][0]:.4f}, {pts['O13a'][1]:.4f})  (C13a, wall east end)")
    print(f"  O13:  ({pts['O13'][0]:.4f}, {pts['O13'][1]:.4f})  (C13, on PiX-Pi5 line)")
    print(f"  O12:  ({pts['O12'][0]:.4f}, {pts['O12'][1]:.4f})  (fillet Po5, incoming tangent)")
    print(f"  O11:  ({pts['O11'][0]:.4f}, {pts['O11'][1]:.4f})  (fillet Po5, exits North)")
    print(f"  O10a: ({pts['O10a'][0]:.4f}, {pts['O10a'][1]:.4f})  ({R_t4*12:.1f}\" arc tangent to N-S line)")
    print(f"  O9:   ({pts['O9'][0]:.4f}, {pts['O9'][1]:.4f})  ({R_t4*12:.1f}\" arc / line tangent)")
    print(f"  O8:   ({pts['O8'][0]:.4f}, {pts['O8'][1]:.4f})  (line / 180° arc tangent)")
    print(f"  O7:   ({pts['O7'][0]:.4f}, {pts['O7'][1]:.4f})  (180° arc west end / fillet N-S tangent)")
    print(f"  O6a:  ({pts['O6a'][0]:.4f}, {pts['O6a'][1]:.4f})  (fillet E-W tangent)")
    print(f"  O6:   ({pts['O6'][0]:.4f}, {pts['O6'][1]:.4f})  (turn2 tangent)")
    print(f"  O5:   ({pts['O5'][0]:.4f}, {pts['O5'][1]:.4f})  (turn1/turn2 junction)")
    print(f"  O4:   ({pts['O4'][0]:.4f}, {pts['O4'][1]:.4f})  (6.0' east of O3)")
    print(f"  O3:   ({pts['O3'][0]:.4f}, {pts['O3'][1]:.4f})  (fillet2 tangent)")
    print(f"  Cw1:  ({pts['Cw1'][0]:.4f}, {pts['Cw1'][1]:.4f})  (C15->C14 arc center)")
    print(f"  Cw2:  ({pts['Cw2'][0]:.4f}, {pts['Cw2'][1]:.4f})  (C14->C13b arc center)")
    print(f"  Cw3:  ({pts['Cw3'][0]:.4f}, {pts['Cw3'][1]:.4f})  (C13a->C13 arc center)")
    print(f"  Wall segment: O13b to O13a, length = {abs(pts['O13a'][0]-pts['O13b'][0])*12:.1f}\"")
    print(f"  Cf4:  ({pts['Cf4'][0]:.4f}, {pts['Cf4'][1]:.4f})  (fillet Po5 center, R={R_f_po5:.4f}')")
    print(f"  Ct4:  ({pts['Ct4'][0]:.4f}, {pts['Ct4'][1]:.4f})  ({R_t4*12:.1f}\" arc center, R={R_t4:.4f}')")
    print(f"  Outline area: {outline_area:.2f} sq ft")

    lines = []
    lines.append(f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}">')
    lines.append(f'<rect width="{W}" height="{H}" fill="white"/>')
    lines.append('<defs>')
    lines.append('  <marker id="ah" markerWidth="8" markerHeight="6" refX="8" refY="3" orient="auto"><polygon points="0 0, 8 3, 0 6" fill="#333"/></marker>')
    lines.append(f'  <clipPath id="page"><rect width="{W}" height="{H}"/></clipPath>')
    lines.append('</defs>')
    lines.append(f'<text x="{W/2}" y="30" text-anchor="middle" font-family="Arial" font-size="14"'
                 f' font-weight="bold">Site Path \u2014 Outline (wall arcs R=28\u2033, fillets R=10\u2033 &amp; R=28\u2033)</text>')

    render_layer(lines, outer_segs, pts, outer_cfg)
    render_layer(lines, inset_segs, pts, inset_cfg)
    render_layer(lines, outline_segs, pts, outline_cfg)
    render_floorplan(lines)

    # Area label centered in outline
    centroid_names = ["O2","O1c","O1b","O1a","O1","O0","O15","O14","O13b","O13a","O13","O12","O11","O10a","O9","O8","O7","O6a","O6","O5","O4","O3"]
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
    lines.append(f'<text x="60" y="{ly+7}" font-family="Arial" font-size="8" fill="#333">Outline path ({outline_area:.2f} sq ft) \u2014 wall arcs R=28", fillets R=10" &amp; R=28"</text>')
    ly += 12
    lines.append(f'<line x1="40" y1="{ly+4}" x2="54" y2="{ly+4}" stroke="#2E7D32" stroke-width="2.5"/>')
    lines.append(f'<text x="60" y="{ly+7}" font-family="Arial" font-size="8" fill="#333">Wall arcs (R=28")</text>')

    # Footer
    lines.append(f'<text x="{W/2}" y="{H-2}" text-anchor="middle" font-family="Arial" font-size="7.5"'
                 f' fill="#999">Bearings as adjusted \u2022 Distances in feet and inches \u2022'
                 f' Outline arcs: wall R=28", fillets R=10" &amp; R=28"</text>')
    lines.append('</svg>')

    svg_content = "\n".join(lines)
    with open(r"c:\Users\Mango Cat\Dev\hut2\path_area.svg", "w", encoding="utf-8") as f:
        f.write(svg_content)

    print(f"\nSVG written to path_area.svg")
    print(f"Outer path area: {outer_area:.2f} sq ft (rendered at 20%)")
    print(f"Inset path area: {inset_area:.2f} sq ft (rendered at 20%)")
    print(f"Outline path area: {outline_area:.2f} sq ft (rendered at 100%)")
    print(f"Outline: O2->O1c->Arc20->O1b->Arc28->O1a->O1->Fillet->O0->O15->Cw1->O14->Cw2->O13b->Wall->O13a->Cw3->O13->O12->FilletPo5->O11->O10a->Turn4->O9->O8->Turn3_180->O7->Fillet6a->O6a->O6->Turn2->O5->Turn1->O4->O3->Fillet2->O2")
