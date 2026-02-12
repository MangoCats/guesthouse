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
    ArcSeg("T3", "PX", "TC3", R3, "CW", 60),
    LineSeg("PX", "P4"), LineSeg("P4", "P5"), LineSeg("P5", "T1"),
    ArcSeg("T1", "PA", "TC1", R1, "CW", 60),
    ArcSeg("PA", "T2", "TC2", R2, "CW", 60),
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

pts["Ti3"] = (pts["TC3"][0], pts["P3"][1] + delta)
pts["Ti1"] = (pts["T1"][0] - delta*nE, pts["T1"][1] - delta*nN)
pts["Ti2"] = (pts["T2"][0] - delta*nE, pts["T2"][1] - delta*nN)

pts["PiX"] = line_circle_isect_min_abs_t(o5, d_e5, pts["TC3"], R3i)
pts["Ai2"] = circle_circle_isect(pts["TC1"], R1i, pts["TC2"], R2i, near=pts["PA"])

inset_segs = [
    LineSeg("PiOB", "Pi2"), LineSeg("Pi2", "Pi3"), LineSeg("Pi3", "Ti3"),
    ArcSeg("Ti3", "PiX", "TC3", R3i, "CW", 60),
    LineSeg("PiX", "Pi4"), LineSeg("Pi4", "Pi5"), LineSeg("Pi5", "Ti1"),
    ArcSeg("Ti1", "Ai2", "TC1", R1i, "CW", 60),
    ArcSeg("Ai2", "Ti2", "TC2", R2i, "CW", 60),
    LineSeg("Ti2", "PiOB"),
]
inset_area = poly_area(path_polygon(inset_segs, pts))

# --- Rotate outer/inset points about U16 so PiX-Pi5 bearing = 60° ---
# Pre-compute U16 position as pivot (perpendicular foot from C15 onto PiX-Pi5)
_R_fp = 28.0 / 12.0
_d_pip = (pts["Pi5"][0] - pts["PiX"][0], pts["Pi5"][1] - pts["PiX"][1])
_L_pip = math.sqrt(_d_pip[0]**2 + _d_pip[1]**2)
_d_pip_u = (_d_pip[0]/_L_pip, _d_pip[1]/_L_pip)
_ln_pip = left_norm(pts["PiX"], pts["Pi5"])
_o_pip = off_pt(pts["PiX"], _ln_pip, _R_fp)
_pre_U1_E = pts["Pi3"][0]
_pre_w5_e = _pre_U1_E + 8.0/12 + 0.5 + 35.0/12 + 3.0 + 2.0 + (3+30+4+140+4+30+3)/12.0
_pre_F15_E = _pre_w5_e + 9.0 + 1.0/12
_t_cf4 = (_pre_F15_E - _R_fp - _o_pip[0]) / _d_pip_u[0]
_cf4 = (_pre_F15_E - _R_fp, _o_pip[1] + _t_cf4 * _d_pip_u[1])
_t_o12 = ((_cf4[0]-pts["PiX"][0])*_d_pip[0] + (_cf4[1]-pts["PiX"][1])*_d_pip[1]) \
         / (_d_pip[0]**2 + _d_pip[1]**2)
_pivot = (pts["PiX"][0] + _t_o12*_d_pip[0], pts["PiX"][1] + _t_o12*_d_pip[1])
# CW rotation angle: target bearing (60°) minus current PiX->Pi5 bearing
_brg_pip = math.degrees(math.atan2(_d_pip[0], _d_pip[1])) % 360
_rot_deg = 60.0 - _brg_pip
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

# ============================================================
# Section 7: Outline Path (inset + arcs at Po3 and Po2)
# ============================================================
pts["Po2"]  = pts["Pi2"]
pts["U21"]  = pts["Ti3"]   # was To3
pts["To2"]  = pts["Ti2"]

R_a0 = 10.0 / 12.0
pts["C0"] = (pts["Pi3"][0] + R_a0, pts["Pi3"][1] + R_a0)
pts["U1"] = (pts["Pi3"][0], pts["C0"][1])
pts["U0"] = (pts["C0"][0], pts["Pi3"][1])

# --- Arc at Po2 corner (R = 28", exact 90°) ---
R_a5 = 28.0 / 12.0
corner2_N = pts["U0"][1] + 26 - 2.0/12  # U6-U7 line is 25'10" north of U0
# C5 and C3 shift 1' east
_nw_shift = 1.0
pts["C5"] = (pts["Po2"][0] + _nw_shift + R_a5, corner2_N - R_a5)
pts["U5"] = (pts["Po2"][0] + _nw_shift, pts["C5"][1])
pts["U6"] = (pts["C5"][0], corner2_N)
# U4: 5'8" south of C5, same easting
pts["U4"] = (pts["U5"][0], pts["U5"][1] - (5.0 + 8.0/12) + 28.0/12)
# Arc U3->U4 and U2->U3: radii from constraints
# Constraint: R_a2 = R_a3 + 8/12 (8" larger), U1-U2 segment = 16'8"
_dR_23 = 8.0 / 12.0
_tgt_12 = 16.0 + 8.0 / 12.0  # target U1-U2 distance
_a23 = pts["U1"][0] - pts["U4"][0]
_K23 = (pts["U4"][1] - pts["U1"][1]) - _tgt_12
_S23 = -(_K23**2 + _a23**2) / (2 * _a23)
R_a3 = (_S23 - _dR_23) / 2
R_a2 = R_a3 + _dR_23
pts["C3"] = (pts["U4"][0] - R_a3, pts["U4"][1])
# External tangency: |C3 - C2| = R_a3 + R_a2, C2 = (U1_E + R_a2, U2_N)
_cc2_E = pts["U1"][0] + R_a2
_dE_cc = _cc2_E - pts["C3"][0]
_dN_cc = -math.sqrt((R_a3 + R_a2)**2 - _dE_cc**2)  # south
pts["U2"] = (pts["U1"][0], pts["C3"][1] + _dN_cc)
pts["C2"] = (_cc2_E, pts["U2"][1])
# U3: tangent point on C3->C2 line
_f1b = R_a3 / (R_a3 + R_a2)
pts["U3"] = (pts["C3"][0] + _f1b * (pts["C2"][0] - pts["C3"][0]),
              pts["C3"][1] + _f1b * (pts["C2"][1] - pts["C3"][1]))
# U7: 6.0' east of U6, adjusted 1' west to undo C5/C3 shift
pts["U7"] = (pts["U6"][0] + 5.5 + 6.0/12 - _nw_shift, pts["U6"][1])
# Arc at U7: 90° CW east->south, R=28"
R_a7 = 28.0 / 12.0
pts["C7"] = (pts["U7"][0], pts["U7"][1] - R_a7)
pts["U8"] = (pts["C7"][0] + R_a7, pts["C7"][1])
# Arc at U8: 90° CCW south->east, R=2"
R_a8 = 2.0 / 12.0
pts["C8"] = (pts["U8"][0] + R_a8, pts["U8"][1])
pts["U9"] = (pts["C8"][0], pts["C8"][1] - R_a8)
# Arc 3: 180° CCW arc + 90° CW arc at F10 (R=2")
R_a11 = 28.0 / 12.0
R_a10 = 2.0 / 12.0  # 90° arc at F10
west_E = pts["TC1"][0] - R1i  # westernmost E on Arc 1o
# U10, C10, U11, C11 eastings determined below from bearing constraint

# --- Arc at Po5 corner (exits North) ---
# Pre-compute U14 Northing (needed for R_a15 constraint below)
_iw1_n_face = pts["U0"][1] + 8.0/12 + 11.5 + 6.0/12
_U14_N = _iw1_n_face + 2.0/12
d_in_po5 = (pts["Pi5"][0] - pts["PiX"][0], pts["Pi5"][1] - pts["PiX"][1])
L_in = math.sqrt(d_in_po5[0]**2 + d_in_po5[1]**2)
d_in_u = (d_in_po5[0]/L_in, d_in_po5[1]/L_in)
# F15 E-coordinate: 9' 2" east of east closet wall (w5_e)
_w5_e = pts["U1"][0] + 8.0/12 + 0.5 + 35.0/12 + 3.0 + 2.0 + (3+30+4+140+4+30+3)/12.0
F15_E = _w5_e + 9.0 + 1.0/12
ln_in_po5 = left_norm(pts["PiX"], pts["Pi5"])
# R_a15 from constraint: U14-U15 segment length = 8'4"
# C15_N = PiX_N + (F15_E-PiX_E)*dN/dE + R_a15*(lnN - (1+lnE)*dN/dE)
_seg_14_15 = 8.0 + 4.0/12
_A_15 = pts["PiX"][1] + (F15_E - pts["PiX"][0]) * d_in_u[1] / d_in_u[0]
_B_15 = ln_in_po5[1] - (1.0 + ln_in_po5[0]) * d_in_u[1] / d_in_u[0]
R_a15 = (_U14_N - _seg_14_15 - _A_15) / _B_15
# Outgoing direction is North; C15 is R_a15 west of F15
o_in_po5 = off_pt(pts["PiX"], ln_in_po5, R_a15)
t_cf4 = (F15_E - R_a15 - o_in_po5[0]) / d_in_u[0]
pts["C15"] = (F15_E - R_a15, o_in_po5[1] + t_cf4 * d_in_u[1])
pts["U15"] = (F15_E, pts["C15"][1])  # directly east of C15
# U16: tangent point on arc circle for 60° incoming bearing (U17->U16)
_brg_f4 = math.radians(60.0)
pts["U16"] = (pts["C15"][0] + R_a15 * math.cos(_brg_f4),
              pts["C15"][1] - R_a15 * math.sin(_brg_f4))

# --- Arc U14-U13: bearing U13->U12 = 345°, U12-U13 length = 9'8" ---
# F14: centered 2" north of IW1 north face (_U14_N computed above)
# Bearing 345° => tangent direction (-sin15°, cos15°), normal n = (cos15°, sin15°)
_brg_off = math.radians(360.0 - 345.0)  # 15°
_nx_t = math.cos(_brg_off)
_ny_t = math.sin(_brg_off)
# R_a13 from U12-U13 distance constraint (9'8"), C13_N = _U14_N fixed
_seg_12_13 = 9.0 + 8.0/12  # target U12-U13 length
_dN_centers = (pts["U9"][1] + R_a10) - _U14_N  # C11_N - C13_N
R_a13 = (_dN_centers - _seg_12_13 * _nx_t) / _ny_t + R_a11
pts["C13"] = (F15_E - R_a13, _U14_N)
pts["U14"] = (F15_E, _U14_N)
# C11 easting from tangent constraint: (C11 - C13) · n = R_a13 - R_a11
_C11_N = pts["U9"][1] + R_a10  # U11_N = C11_N (independent of _corner_E)
_C13_E, _C13_N = pts["C13"]
_C11_E = _C13_E + (R_a13 - R_a11 - (_C11_N - _C13_N) * _ny_t) / _nx_t
# Derive corner easting and dependent points (U10, U11, C11)
_corner_E = _C11_E - R_a11
pts["U10"] = (_corner_E - R_a10, pts["U9"][1])
pts["C10"] = (_corner_E - R_a10, pts["U9"][1] + R_a10)
pts["U11"] = (_corner_E, pts["U9"][1] + R_a10)
pts["C11"] = (_C11_E, _C11_N)
# U13, U12: tangent points (shared outward normal for external tangent)
pts["U13"] = (pts["C13"][0] + R_a13 * _nx_t, pts["C13"][1] + R_a13 * _ny_t)
pts["U12"] = (pts["C11"][0] + R_a11 * _nx_t, pts["C11"][1] + R_a11 * _ny_t)

# --- F17-F20-F21 wall geometry (wall at -6" Northing) ---
R_a20   = R_a3                  # arc C20 radius = arc C3 radius
R_a19   = R_a2                  # arc C19 radius = arc C2 radius
wall_south_N = -6.0 / 12.0    # south face of wall at -6" Northing
# Tangency distance between C20 and C19 arc centers
dN_c = (wall_south_N + R_a19) - (pts["U21"][1] - R_a20)
dE_c = math.sqrt((R_a20 + R_a19)**2 - dN_c**2)
# Align F19 with east side of king bed
# bed_e = inner_W1_E + 20.5' (dryer+counter+closet2+W1S+half_bedroom+half_bed)
_bed_e_align = pts["U1"][0] + 8.0/12 + 20.5
pts["U21"] = (_bed_e_align - dE_c - 2.0/12, pts["U21"][1])
# U17 on line from U16 at bearing 60° (U17->U16), tangent to C17 arc
_brg_13 = math.radians(60.0)
_sin_b = math.sin(_brg_13)
_cos_b = math.cos(_brg_13)
# R_a17 from constraint: U16-U17 segment length = 5'
_seg_16_17 = 5.0
R_a17 = (pts["U16"][1] - wall_south_N - _seg_16_17 * _cos_b) / (1.0 - _sin_b)
F17_N = wall_south_N + R_a17 * (1.0 - _sin_b)
_t_13 = (pts["U16"][1] - F17_N) / _cos_b
pts["U17"] = (pts["U16"][0] - _t_13 * _sin_b, F17_N)
# U17->U18 arc center (C17) and tangent point U18
pts["C17"] = (pts["U17"][0] - R_a17 * _cos_b, wall_south_N + R_a17)
pts["U18"] = (pts["C17"][0], wall_south_N)
# U20->U21 arc center (C20): center south of U21
pts["C20"] = (pts["U21"][0], pts["U21"][1] - R_a20)
# U19->U20 arc center (C19): from external tangency |C20-C19| = R_a20+R_a19
F19_E = pts["U21"][0] + dE_c
pts["U19"] = (F19_E, wall_south_N)
pts["C19"] = (F19_E, wall_south_N + R_a19)
# U20 = tangent point between the two arcs (R_a20 from C20 along C20->C19)
_f_w = R_a20 / (R_a20 + R_a19)
pts["U20"] = (pts["C20"][0] + _f_w * (pts["C19"][0] - pts["C20"][0]),
              pts["C20"][1] + _f_w * (pts["C19"][1] - pts["C20"][1]))

outline_segs = [
    ArcSeg("U0", "U1", "C0", R_a0, "CW", 20),           # 0: arc C0
    LineSeg("U1", "U2"),                                      # 1
    ArcSeg("U2", "U3", "C2", R_a2, "CW", 20),             # 2: arc C2
    ArcSeg("U3", "U4", "C3", R_a3, "CCW", 20),            # 3: arc C3
    LineSeg("U4", "U5"),                                      # 4
    ArcSeg("U5", "U6", "C5", R_a5, "CW", 20),        # 5: arc C5
    LineSeg("U6", "U7"),                                      # 6
    ArcSeg("U7", "U8", "C7", R_a7, "CW", 20),          # 7: arc C7
    ArcSeg("U8", "U9", "C8", R_a8, "CCW", 20),         # 8: arc C8
    LineSeg("U9", "U10"),                                     # 9
    ArcSeg("U10", "U11", "C10", R_a10, "CCW", 20),         # 10: arc C10
    ArcSeg("U11", "U12", "C11", R_a11, "CW", 60),        # 11: 180° arc
    LineSeg("U12", "U13"),                                    # 12
    ArcSeg("U13", "U14", "C13", R_a13, "CW", 60),           # 13: 28" arc
    LineSeg("U14", "U15"),                                    # 14
    ArcSeg("U15", "U16", "C15", R_a15, "CW", 20),        # 15: arc C15
    LineSeg("U16", "U17"),                                    # 16
    ArcSeg("U17", "U18", "C17", R_a17, "CW", 20),         # 17: arc C17
    LineSeg("U18", "U19"),                                    # 18
    ArcSeg("U19", "U20", "C19", R_a19, "CW", 60),           # 19: arc C19
    ArcSeg("U20", "U21", "C20", R_a20, "CCW", 60),          # 20: arc C20
    LineSeg("U21", "U0"),                                     # 21
]
outline_area = poly_area(path_polygon(outline_segs, pts))

# ============================================================
# Section 7b: Floorplan Inner Walls + Interior Elements
# ============================================================
_wall_t = 8.0 / 12.0
_radii = {
    "R_a0": R_a0, "R_a20": R_a20, "R_a19": R_a19, "R_a17": R_a17,
    "R_a15": R_a15, "R_a11": R_a11, "R_a8": R_a8,
    "R_a7": R_a7, "R_a5": R_a5, "R_a13": R_a13, "R_a10": R_a10,
    "R_a3": R_a3, "R_a2": R_a2,
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
sweep_a0 = arc_sweep_deg(outline_segs[0], pts)    # U0->U1 (C0)
sweep_a2 = arc_sweep_deg(outline_segs[2], pts)   # U2->U3 (C2)
sweep_a3 = arc_sweep_deg(outline_segs[3], pts)   # U3->U4 (C3)
sweep_a5 = arc_sweep_deg(outline_segs[5], pts)   # U5->U6 (C5)
sweep_a7 = arc_sweep_deg(outline_segs[7], pts)   # U7->U8 (C7)
sweep_a8 = arc_sweep_deg(outline_segs[8], pts)   # U8->U9 (C8)
sweep_a10 = arc_sweep_deg(outline_segs[10], pts)  # U10->U11 (C10)
sweep_a11 = arc_sweep_deg(outline_segs[11], pts)  # U11->U12 (C11)
sweep_a13 = arc_sweep_deg(outline_segs[13], pts)  # U13->U14 (C13)
sweep_a15 = arc_sweep_deg(outline_segs[15], pts)  # U15->U16 (C15)
sweep_a17 = arc_sweep_deg(outline_segs[17], pts)  # U17->U18 (C17)
sweep_a19 = arc_sweep_deg(outline_segs[19], pts)  # U19->U20 (C19)
sweep_a20 = arc_sweep_deg(outline_segs[20], pts)  # U20->U21 (C20)

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
        ("U0","U1"):    ("#333", 2.0),
        ("U2","U3"):    ("#333", 2.0),
        ("U3","U4"):    ("#333", 2.0),
        ("U5","U6"):    ("#333", 2.0),
        ("U7","U8"):    ("#333", 2.0),
        ("U8","U9"):    ("#333", 2.0),
        ("U10","U11"):  ("#333", 2.0),
        ("U11","U12"):  ("#333", 2.0),
        ("U13","U14"):  ("#333", 2.0),
        ("U15","U16"):  ("#333", 2.0),
        ("U17","U18"):  ("#333", 2.0),
        ("U19","U20"):  ("#333", 2.0),
        ("U20","U21"):  ("#333", 2.0),
    },
    vertex_styles={
        # Vertically centered at Northing (label offset E/W, dy≈+4 centers font-10 baseline)
        "U5":   VertexStyle("U5",   "end",    -8,  4,  "#d32f2f", 1.75, 10),
        "U4":   VertexStyle("U4",   "end",    -8,  4,  "#d32f2f", 1.75, 10),
        "U3":   VertexStyle("U3",   "end",   -10,  4,  "#d32f2f", 1.75, 10),
        "U2":   VertexStyle("U2",   "end",    -8,  4,  "#d32f2f", 1.75, 10),
        "U1":   VertexStyle("U1",   "end",    -8,  4,  "#d32f2f", 1.75, 10),
        "U8":   VertexStyle("U8",   "end",    -8,  4,  "#d32f2f", 1.75, 10),
        "U11":  VertexStyle("U11",  "start",   8,  4,  "#d32f2f", 1.75, 10),
        "U12":  VertexStyle("U12",  "start",   8,  4,  "#d32f2f", 1.75, 10),
        "U13":  VertexStyle("U13",  "start",   8,  4,  "#d32f2f", 1.75, 10),
        "U14":  VertexStyle("U14",  "start",  10,  4,  "#d32f2f", 1.75, 10),
        "U15":  VertexStyle("U15",  "start",   8,  4,  "#d32f2f", 1.75, 10),
        # Horizontally centered at Easting (label offset N/S)
        "U0":   VertexStyle("U0",   "middle",  0, 10,  "#d32f2f", 1.75, 10),
        "U6":   VertexStyle("U6",   "middle",  0, -6,  "#d32f2f", 1.75, 10),
        "U7":   VertexStyle("U7",   "middle",  0, -6,  "#d32f2f", 1.75, 10),
        "U9":   VertexStyle("U9",   "middle",  0, 17,  "#d32f2f", 1.75, 10),
        "U10":  VertexStyle("U10",  "middle",  0, 17,  "#d32f2f", 1.75, 10),
        "U17":  VertexStyle("U17",  "middle",  0, 13,  "#d32f2f", 1.75, 10),
        "U18":  VertexStyle("U18",  "middle",  0, 12,  "#d32f2f", 1.75, 10),
        "U19":  VertexStyle("U19",  "middle",  0, 12,  "#d32f2f", 1.75, 10),
        "U20":  VertexStyle("U20",  "middle",  0, 13,  "#d32f2f", 1.75, 10),
        "U21":  VertexStyle("U21",  "middle",  0, 10,  "#d32f2f", 1.75, 10),
        # Default (diagonal segment)
        "U16":  VertexStyle("U16",  "start",   8,  4,  "#d32f2f", 1.75, 10),
    },
    brg_dist_labels={
        ("U1","U2"): BrgDistLabel(18),
        ("U4","U5"): BrgDistLabel(18),
        ("U6","U7"): BrgDistLabel(-16),
        ("U9","U10"): BrgDistLabel(-16),
        ("U12","U13"): BrgDistLabel(-16),
        ("U14","U15"): BrgDistLabel(-16),
        ("U16","U17"): BrgDistLabel(16),
        ("U18","U19"): BrgDistLabel(-16),
        ("U21","U0"): BrgDistLabel(-18),
    },
    arc_labels={
        ("U0","U1"): ArcLabel("Arc R=10\"",
            f"{sweep_a0:.1f}\u00b0", "end", -10, 14, 11, "#333"),
        ("U2","U3"): ArcLabel(f"Arc R={R_a2*12:.1f}\"",
            f"{sweep_a2:.1f}\u00b0 CCW", "start", 12, 0, 11, "#333"),
        ("U3","U4"): ArcLabel(f"Arc R={R_a3*12:.1f}\"",
            f"{sweep_a3:.1f}\u00b0 CW", "end", -10, -10, 11, "#333"),
        ("U5","U6"): ArcLabel("Arc R=28\"",
            f"{sweep_a5:.1f}\u00b0", "end", -10, -14, 11, "#333"),
        ("U7","U8"): ArcLabel("Arc R=28\"",
            f"{sweep_a7:.1f}\u00b0", "start", 12, 0, 11, "#333"),
        ("U8","U9"): ArcLabel("Arc R=2\"",
            f"{sweep_a8:.1f}\u00b0", "end", -10, 14, 11, "#333"),
        ("U10","U11"): ArcLabel("Arc R=2\"",
            f"{sweep_a10:.1f}\u00b0", "end", -10, -10, 11, "#333"),
        ("U11","U12"): ArcLabel("Arc R=28\"",
            f"{sweep_a11:.1f}\u00b0 CCW", "start", 12, 0, 11, "#333"),
        ("U13","U14"): ArcLabel(f"Arc R={R_a13*12:.0f}\"",
            f"{sweep_a13:.1f}\u00b0", "start", 12, 0, 11, "#333"),
        ("U15","U16"): ArcLabel(f"Arc R={R_a15*12:.0f}\"",
            f"{sweep_a15:.1f}\u00b0", "start", 10, -10, 11, "#333"),
        ("U17","U18"): ArcLabel(f"Arc R={R_a17*12:.0f}\"",
            f"{sweep_a17:.1f}\u00b0", "end", -10, -10, 11, "#333"),
        ("U19","U20"): ArcLabel(f"Arc R={R_a19*12:.1f}\"",
            f"{sweep_a19:.1f}\u00b0", "start", 12, -10, 11, "#333"),
        ("U20","U21"): ArcLabel(f"Arc R={R_a20*12:.1f}\"",
            f"{sweep_a20:.1f}\u00b0 CW", "start", 12, 4, 11, "#333"),
    },
    center_marks=[
        CenterMark("C0", "U0", "#333"), CenterMark("C2", "U2", "#333"),
        CenterMark("C3", "U3", "#333"), CenterMark("C5", "U5", "#333"),
        CenterMark("C7", "U7", "#333"), CenterMark("C8", "U8", "#333"),
        CenterMark("C10", "U10", "#333"), CenterMark("C11", "U11", "#333"),
        CenterMark("C13", "U13", "#333"), CenterMark("C15", "U15", "#333"),
        CenterMark("C17", "U17", "#333"), CenterMark("C19", "U19", "#333"),
        CenterMark("C20", "U20", "#333"),
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
    print(f'=== OUTLINE PATH ===')
    print(f"  U0:   ({pts['U0'][0]:.4f}, {pts['U0'][1]:.4f})  (arc tangent)")
    print(f"  U1:   ({pts['U1'][0]:.4f}, {pts['U1'][1]:.4f})  (arc tangent)")
    print(f"  U2:   ({pts['U2'][0]:.4f}, {pts['U2'][1]:.4f})  (same E as U1)")
    print(f"  U3:   ({pts['U3'][0]:.4f}, {pts['U3'][1]:.4f})  (arc tangent point)")
    print(f"  U4:   ({pts['U4'][0]:.4f}, {pts['U4'][1]:.4f})  (5'8\" south of U5)")
    print(f"  U5:   ({pts['U5'][0]:.4f}, {pts['U5'][1]:.4f})  (arc tangent)")
    print(f"  U6:   ({pts['U6'][0]:.4f}, {pts['U6'][1]:.4f})  (arc tangent)")
    print(f"  U7:   ({pts['U7'][0]:.4f}, {pts['U7'][1]:.4f})  (6.0' east of U6)")
    print(f"  U8:   ({pts['U8'][0]:.4f}, {pts['U8'][1]:.4f})  (C7/C8 arc junction)")
    print(f"  U9:   ({pts['U9'][0]:.4f}, {pts['U9'][1]:.4f})  (arc tangent)")
    print(f"  U10:  ({pts['U10'][0]:.4f}, {pts['U10'][1]:.4f})  (arc E-W tangent)")
    print(f"  U11:  ({pts['U11'][0]:.4f}, {pts['U11'][1]:.4f})  (180° arc west end / arc N-S tangent)")
    print(f"  U12:  ({pts['U12'][0]:.4f}, {pts['U12'][1]:.4f})  (line / 180° arc tangent)")
    print(f"  U13:  ({pts['U13'][0]:.4f}, {pts['U13'][1]:.4f})  ({R_a13*12:.1f}\" arc / line tangent)")
    print(f"  U14:  ({pts['U14'][0]:.4f}, {pts['U14'][1]:.4f})  ({R_a13*12:.1f}\" arc tangent to N-S line)")
    print(f"  U15:  ({pts['U15'][0]:.4f}, {pts['U15'][1]:.4f})  (arc C15, exits North)")
    print(f"  U16:  ({pts['U16'][0]:.4f}, {pts['U16'][1]:.4f})  (arc C15, incoming tangent)")
    print(f"  U17:  ({pts['U17'][0]:.4f}, {pts['U17'][1]:.4f})  (F17, on PiX-Pi5 line)")
    print(f"  U18:  ({pts['U18'][0]:.4f}, {pts['U18'][1]:.4f})  (arc C17 tangent)")
    print(f"  U19:  ({pts['U19'][0]:.4f}, {pts['U19'][1]:.4f})  (arc C19 tangent)")
    print(f"  U20:  ({pts['U20'][0]:.4f}, {pts['U20'][1]:.4f})  (F20, arc junction)")
    print(f"  U21:  ({pts['U21'][0]:.4f}, {pts['U21'][1]:.4f})  (F21, was To3)")
    print(f"  C20:  ({pts['C20'][0]:.4f}, {pts['C20'][1]:.4f})  (U20->U21 arc center)")
    print(f"  C19:  ({pts['C19'][0]:.4f}, {pts['C19'][1]:.4f})  (U19->U20 arc center)")
    print(f"  C17:  ({pts['C17'][0]:.4f}, {pts['C17'][1]:.4f})  (U17->U18 arc center)")
    print(f"  U18-U19 segment length = {abs(pts['U18'][0]-pts['U19'][0])*12:.1f}\"")
    print(f"  C15:  ({pts['C15'][0]:.4f}, {pts['C15'][1]:.4f})  (arc C15 center, R={R_a15:.4f}')")
    print(f"  C13:  ({pts['C13'][0]:.4f}, {pts['C13'][1]:.4f})  ({R_a13*12:.1f}\" arc center, R={R_a13:.4f}')")
    print(f"  Outline area: {outline_area:.2f} sq ft")

    lines = []
    lines.append(f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}">')
    lines.append(f'<rect width="{W}" height="{H}" fill="white"/>')
    lines.append('<defs>')
    lines.append('  <marker id="ah" markerWidth="8" markerHeight="6" refX="8" refY="3" orient="auto"><polygon points="0 0, 8 3, 0 6" fill="#333"/></marker>')
    lines.append(f'  <clipPath id="page"><rect width="{W}" height="{H}"/></clipPath>')
    lines.append('</defs>')
    lines.append(f'<text x="{W/2}" y="30" text-anchor="middle" font-family="Arial" font-size="14"'
                 f' font-weight="bold">Site Path \u2014 Outline</text>')

    render_layer(lines, outer_segs, pts_rot, outer_cfg)
    render_layer(lines, inset_segs, pts_rot, inset_cfg)
    render_layer(lines, outline_segs, pts, outline_cfg)
    render_floorplan(lines)

    # Area label centered in outline
    centroid_names = ["U0","U1","U2","U3","U4","U5","U6","U7","U8","U9","U10","U11","U12","U13","U14","U15","U16","U17","U18","U19","U20","U21"]
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
    lines.append(f'<text x="{W/2}" y="{H-2}" text-anchor="middle" font-family="Arial" font-size="7.5"'
                 f' fill="#999">Bearings as adjusted \u2022 Distances in feet and inches</text>')
    lines.append('</svg>')

    svg_content = "\n".join(lines)
    with open(r"c:\Users\Mango Cat\Dev\hut2\path_area.svg", "w", encoding="utf-8") as f:
        f.write(svg_content)

    print(f"\nSVG written to path_area.svg")
    print(f"Outer path area: {outer_area:.2f} sq ft (rendered at 20%)")
    print(f"Inset path area: {inset_area:.2f} sq ft (rendered at 20%)")
    print(f"Outline path area: {outline_area:.2f} sq ft (rendered at 100%)")
    print(f"Outline: U0->ArcC0->U1->U2->ArcC2->U3->ArcC3->U4->U5->ArcC5->U6->U7->ArcC7->U8->ArcC8->U9->U10->ArcC10->U11->ArcC11_180->U12->U13->ArcC13->U14->U15->ArcC15->U16->U17->ArcC17->U18->U19->ArcC19->U20->ArcC20->U21->U0")
