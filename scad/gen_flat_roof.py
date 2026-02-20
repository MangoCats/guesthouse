"""Generate flat_roof.scad — OpenSCAD 3D model of double-shell concrete walls.

Models wall sections between openings, each as an outer boundary minus
its air-gap cavity, using circular arcs and line segments.
Matches the wall section geometry from walls.svg.
All coordinates in feet.
"""
import math
import os

from floorplan.gen_floorplan import build_floorplan_data
from floorplan.constants import (WALL_OUTER, SHELL_THICKNESS, AIR_GAP,
                                  OPENING_INSIDE_RADIUS, F8F9_INNER_TURN_R)
from shared.wall_shells import enumerate_wall_sections, lerp
from shared.types import ArcSeg

_DIR = os.path.dirname(os.path.abspath(__file__))
_OUT = os.path.join(_DIR, "flat_roof.scad")

WALL_HEIGHT_IN = 80.0  # inches
WALL_HEIGHT_FT = WALL_HEIGHT_IN / 12.0


# ── element types ────────────────────────────────────────────
# ("pts", [(x,y), ...])            — polyline vertices
# ("arc", cx, cy, r, a1_deg, a2_deg) — circular arc


def _arc_n(sweep_deg):
    """Segment count for an arc (~3 deg per segment, min 4)."""
    return max(4, round(abs(sweep_deg) / 3))


def _seg_arc_elem(seg, pts):
    """Convert an ArcSeg to an arc element."""
    c = pts[seg.center]
    a1 = math.degrees(math.atan2(pts[seg.start][1] - c[1],
                                  pts[seg.start][0] - c[0]))
    a2 = math.degrees(math.atan2(pts[seg.end][1] - c[1],
                                  pts[seg.end][0] - c[0]))
    if seg.direction == "CW":
        sweep = (a1 - a2) % 360
        a2 = a1 - sweep
    else:
        sweep = (a2 - a1) % 360
        a2 = a1 + sweep
    return ("arc", c[0], c[1], seg.radius, a1, a2)


def _qarc_elem(cx, cy, R, u0, u1):
    """Quarter-circle arc from direction u0 toward u1.

    Returns (element, start_pt, end_pt).
    """
    a1 = math.degrees(math.atan2(u0[1], u0[0]))
    a2 = math.degrees(math.atan2(u1[1], u1[0]))
    cross = u0[0] * u1[1] - u0[1] * u1[0]
    if cross > 0:          # CCW
        while a2 < a1:
            a2 += 360
        if a2 - a1 > 180:
            a2 -= 360
    elif cross < 0:        # CW
        while a2 > a1:
            a2 -= 360
        if a1 - a2 > 180:
            a2 += 360
    sp = (cx + R * u0[0], cy + R * u0[1])
    ep = (cx + R * u1[0], cy + R * u1[1])
    return ("arc", cx, cy, R, a1, a2), sp, ep


def _f8f9_elems(pts, inset, R_turn):
    """F8-F9 corner: straight south, arc 180-270 CCW, straight east."""
    F8, C8 = pts["F8"], pts["C8"]
    R_a8 = C8[0] - F8[0]
    se, sn = F8[0] - inset, F8[1]
    ee, en_ = C8[0], F8[1] - R_a8 - inset
    d = R_a8 + inset - R_turn
    acx, acy = se + R_turn, sn - d
    return [
        ("pts", [(se, sn), (se, acy)]),
        ("arc", acx, acy, R_turn, 180, 270),
        ("pts", [(acx, acy - R_turn), (ee, en_)]),
    ]


# ── U-turn arcs ──────────────────────────────────────────────

def _uturn_arcs(pts, outline_segs, inner_segs, seg_idx, t_param, side,
                shell_t, R_in, wall_total):
    """U-turn arc elements at an opening boundary.

    Returns dict: key -> (element, start_pt, end_pt).
    """
    R_out = R_in + shell_t
    seg = outline_segs[seg_idx]
    F_pt = lerp(pts[seg.start], pts[seg.end], t_param)
    iseg = inner_segs[seg_idx]
    W_pt = lerp(pts[iseg.start], pts[iseg.end], t_param)

    dx = pts[seg.end][0] - pts[seg.start][0]
    dy = pts[seg.end][1] - pts[seg.start][1]
    ln = math.sqrt(dx * dx + dy * dy)
    th = (dx / ln, dy / ln)
    ne = (-th[1], th[0])
    od = th if side == "start" else (-th[0], -th[1])
    wd = (-od[0], -od[1])
    ni = (-ne[0], -ne[1])

    oc = (F_pt[0] - R_out * ne[0] + R_out * wd[0],
          F_pt[1] - R_out * ne[1] + R_out * wd[1])
    ic = (W_pt[0] + R_out * ne[0] + R_out * wd[0],
          W_pt[1] + R_out * ne[1] + R_out * wd[1])

    return {
        'oc_F': _qarc_elem(oc[0], oc[1], R_out, ne, od),
        'oc_S': _qarc_elem(oc[0], oc[1], R_in, ne, od),
        'ic_W': _qarc_elem(ic[0], ic[1], R_out, od, ni),
        'ic_G': _qarc_elem(ic[0], ic[1], R_in, od, ni),
    }


# ── boundary tracing ─────────────────────────────────────────

def _trace_elems(pts, segs, start_seg_idx, start_t, end_seg_idx, end_t,
                 R_out, seg_elem_overrides=None):
    """Trace boundary path as element list."""
    n = len(segs)
    elems = []

    def slen(seg):
        a, b = pts[seg.start], pts[seg.end]
        return math.sqrt((b[0] - a[0]) ** 2 + (b[1] - a[1]) ** 2)

    if start_seg_idx == end_seg_idx:
        seg = segs[start_seg_idx]
        a, b = pts[seg.start], pts[seg.end]
        dt = R_out / slen(seg)
        elems.append(("pts", [lerp(a, b, start_t + dt),
                               lerp(a, b, end_t - dt)]))
        return elems

    # Start segment: partial line
    seg = segs[start_seg_idx]
    a, b = pts[seg.start], pts[seg.end]
    dt = R_out / slen(seg)
    elems.append(("pts", [lerp(a, b, start_t + dt), b]))

    # Intermediate full segments
    idx = (start_seg_idx + 1) % n
    while idx != end_seg_idx:
        seg = segs[idx]
        if seg_elem_overrides and idx in seg_elem_overrides:
            elems.extend(seg_elem_overrides[idx])
        elif isinstance(seg, ArcSeg):
            elems.append(_seg_arc_elem(seg, pts))
        else:
            elems.append(("pts", [pts[seg.start], pts[seg.end]]))
        idx = (idx + 1) % n

    # End segment: partial line
    seg = segs[end_seg_idx]
    a, b = pts[seg.start], pts[seg.end]
    dt = R_out / slen(seg)
    elems.append(("pts", [a, lerp(a, b, end_t - dt)]))

    return elems


def _rev_elem(e):
    """Reverse a single element."""
    if e[0] == "pts":
        return ("pts", list(reversed(e[1])))
    _, cx, cy, r, a1, a2 = e
    return ("arc", cx, cy, r, a2, a1)


def _rev_elems(elems):
    """Reverse element list and each element's direction."""
    return [_rev_elem(e) for e in reversed(elems)]


# ── section assembly ──────────────────────────────────────────

def _build_section(pts, outline_segs, inner_segs, s_segs, g_segs,
                   start_op, end_op, shell_t, R_in, wall_total,
                   g_ov, w_ov):
    """Build (outer_elements, cavity_elements) for one wall section."""
    R_out = R_in + shell_t

    su = _uturn_arcs(pts, outline_segs, inner_segs,
                     start_op.seg_idx, start_op.t_end, "end",
                     shell_t, R_in, wall_total)
    eu = _uturn_arcs(pts, outline_segs, inner_segs,
                     end_op.seg_idx, end_op.t_start, "start",
                     shell_t, R_in, wall_total)

    f = _trace_elems(pts, outline_segs, start_op.seg_idx, start_op.t_end,
                     end_op.seg_idx, end_op.t_start, R_out)
    s = _trace_elems(pts, s_segs, start_op.seg_idx, start_op.t_end,
                     end_op.seg_idx, end_op.t_start, R_out)
    g = _trace_elems(pts, g_segs, start_op.seg_idx, start_op.t_end,
                     end_op.seg_idx, end_op.t_start, R_out,
                     seg_elem_overrides=g_ov)
    w = _trace_elems(pts, inner_segs, start_op.seg_idx, start_op.t_end,
                     end_op.seg_idx, end_op.t_start, R_out,
                     seg_elem_overrides=w_ov)

    # Outer: F fwd → end oc_F → xwall → end ic_W → W bwd
    #      → start ic_W rev → xwall → start oc_F rev
    e_ocF, _, e_ocF_ep = eu['oc_F']
    e_icW, e_icW_sp, _ = eu['ic_W']
    s_icW, s_icW_sp, s_icW_ep = su['ic_W']
    s_ocF, s_ocF_sp, s_ocF_ep = su['oc_F']

    outer = list(f)
    outer.append(e_ocF)
    outer.append(("pts", [e_ocF_ep, e_icW_sp]))
    outer.append(e_icW)
    outer.extend(_rev_elems(w))
    outer.append(_rev_elem(s_icW))
    outer.append(("pts", [s_icW_sp, s_ocF_ep]))
    outer.append(_rev_elem(s_ocF))

    # Cavity: S fwd → end oc_S → xwall → end ic_G → G bwd
    #       → start ic_G rev → xwall → start oc_S rev
    e_ocS, _, e_ocS_ep = eu['oc_S']
    e_icG, e_icG_sp, _ = eu['ic_G']
    s_icG, s_icG_sp, s_icG_ep = su['ic_G']
    s_ocS, s_ocS_sp, s_ocS_ep = su['oc_S']

    cavity = list(s)
    cavity.append(e_ocS)
    cavity.append(("pts", [e_ocS_ep, e_icG_sp]))
    cavity.append(e_icG)
    cavity.extend(_rev_elems(g))
    cavity.append(_rev_elem(s_icG))
    cavity.append(("pts", [s_icG_sp, s_ocS_ep]))
    cavity.append(_rev_elem(s_ocS))

    return outer, cavity


# ── OpenSCAD output ───────────────────────────────────────────

def _scad_elem(elem):
    """Render one element as an OpenSCAD expression."""
    if elem[0] == "pts":
        body = ", ".join(f"[{p[0]:.6f}, {p[1]:.6f}]" for p in elem[1])
        return f"[{body}]"
    _, cx, cy, r, a1, a2 = elem
    n = _arc_n(a2 - a1)
    return f"arc_pts({cx:.6f}, {cy:.6f}, {r:.6f}, {a1:.3f}, {a2:.3f}, {n})"


def _scad_polygon(name, elements):
    """Emit OpenSCAD module for a polygon built from arc/line elements."""
    lines = [f"module {name}() {{"]
    lines.append("  polygon(points = concat(")
    for i, elem in enumerate(elements):
        expr = _scad_elem(elem)
        if i > 0:
            expr = f"tail({expr})"
        comma = "," if i < len(elements) - 1 else ""
        lines.append(f"    {expr}{comma}")
    lines.append("  ));")
    lines.append("}")
    return "\n".join(lines)


# ── main ──────────────────────────────────────────────────────

def generate():
    fp = build_floorplan_data()
    pts = fp.pts
    outline_segs = fp.outline_segs
    inner_segs = fp.inner_segs
    s_segs = fp.s_segs
    g_segs = fp.g_segs
    openings = fp.openings

    shell_t = SHELL_THICKNESS
    R_in = OPENING_INSIDE_RADIUS

    g_ov = {7: _f8f9_elems(pts, SHELL_THICKNESS + AIR_GAP,
                            OPENING_INSIDE_RADIUS)}
    w_ov = {7: _f8f9_elems(pts, WALL_OUTER, F8F9_INNER_TURN_R)}

    sections = enumerate_wall_sections(openings, outline_segs)
    sections = sections[-1:] + sections[:-1]

    section_data = []
    for start_op, end_op in sections:
        outer, cavity = _build_section(
            pts, outline_segs, inner_segs, s_segs, g_segs,
            start_op, end_op, shell_t, R_in, WALL_OUTER, g_ov, w_ov)
        label = f"{start_op.name}_{end_op.name}"
        section_data.append((label, outer, cavity))

    out = []
    out.append("// flat_roof.scad — Double-shell concrete wall extrusion")
    out.append(f"// Wall height: {WALL_HEIGHT_IN:.0f}\" ({WALL_HEIGHT_FT:.4f} ft)")
    out.append("// Construction: 2\" outer shell / 4\" air gap / 2\" inner shell")
    out.append(f"// {len(section_data)} wall sections between openings")
    out.append("// Units: feet")
    out.append("// Generated by scad/gen_flat_roof.py")
    out.append("")
    out.append("// --- Helper functions ---")
    out.append("function arc_pts(cx, cy, r, a1, a2, n) =")
    out.append("  [for (i=[0:n]) [cx + r*cos(a1 + i*(a2-a1)/n),")
    out.append("                   cy + r*sin(a1 + i*(a2-a1)/n)]];")
    out.append("")
    out.append("function tail(v) = len(v) > 1")
    out.append("  ? [for (i=[1:len(v)-1]) v[i]] : [];")
    out.append("")
    out.append(f"wall_height = {WALL_HEIGHT_FT:.6f};")
    out.append("")

    for label, outer, cavity in section_data:
        out.append(f"// Wall section {label.replace('_', ' to ')}")
        out.append(_scad_polygon(f"wall_{label}_outer", outer))
        out.append("")
        out.append(_scad_polygon(f"wall_{label}_cavity", cavity))
        out.append("")

    out.append("// --- Assembly ---")
    out.append("union() {")
    for label, _, _ in section_data:
        out.append("  linear_extrude(height = wall_height)")
        out.append("    difference() {")
        out.append(f"      wall_{label}_outer();")
        out.append(f"      wall_{label}_cavity();")
        out.append("    }")
    out.append("}")
    out.append("")

    with open(_OUT, "w") as f:
        f.write("\n".join(out))
    print(f"wrote {_OUT}")


if __name__ == "__main__":
    generate()
