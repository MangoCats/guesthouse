"""Generate flat_roof.scad — T-path shell centerline approach.

Each wall section has a single T-path tracing the center of the shell
material through both shells. OpenSCAD offsets this path by ±shell_t/2
to produce the shell boundary polygon.
All coordinates in feet.
"""
import math
import os

from floorplan.gen_floorplan import build_floorplan_data
from floorplan.constants import (WALL_OUTER, SHELL_THICKNESS,
                                  OPENING_INSIDE_RADIUS)
from shared.wall_shells import compute_inset_path, enumerate_wall_sections, lerp
from shared.types import ArcSeg

_DIR = os.path.dirname(os.path.abspath(__file__))
_OUT = os.path.join(_DIR, "flat_roof.scad")

WALL_HEIGHT_IN = 80.0  # inches
WALL_HEIGHT_FT = WALL_HEIGHT_IN / 12.0


# ── element types ────────────────────────────────────────────
# ("line", x1, y1, x2, y2)             — line segment
# ("arc", cx, cy, r, a1_deg, a2_deg)   — circular arc


def _seg_to_elem(seg, pts):
    """Convert a LineSeg or ArcSeg to a T-path element."""
    if isinstance(seg, ArcSeg):
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
    return ("line", pts[seg.start][0], pts[seg.start][1],
                    pts[seg.end][0], pts[seg.end][1])


def _f8f9_elems(pts, inset, R_turn):
    """F8-F9 corner override: straight south, arc 180→270, straight east."""
    F8, C8 = pts["F8"], pts["C8"]
    R_a8 = C8[0] - F8[0]
    se, sn = F8[0] - inset, F8[1]
    ee, en_ = C8[0], F8[1] - R_a8 - inset
    d = R_a8 + inset - R_turn
    acx, acy = se + R_turn, sn - d
    return [
        ("line", se, sn, se, acy),
        ("arc", acx, acy, R_turn, 180, 270),
        ("line", acx, acy - R_turn, ee, en_),
    ]


# ── boundary tracing ─────────────────────────────────────────

def _trace_elems(pts, segs, start_seg_idx, start_t, end_seg_idx, end_t,
                 R_out, seg_overrides=None):
    """Trace T-path elements between two opening boundaries."""
    n = len(segs)
    elems = []

    def slen(seg):
        a, b = pts[seg.start], pts[seg.end]
        return math.sqrt((b[0] - a[0]) ** 2 + (b[1] - a[1]) ** 2)

    if start_seg_idx == end_seg_idx:
        seg = segs[start_seg_idx]
        a, b = pts[seg.start], pts[seg.end]
        dt = R_out / slen(seg)
        p1 = lerp(a, b, start_t + dt)
        p2 = lerp(a, b, end_t - dt)
        elems.append(("line", p1[0], p1[1], p2[0], p2[1]))
        return elems

    # Start segment: partial line
    seg = segs[start_seg_idx]
    a, b = pts[seg.start], pts[seg.end]
    dt = R_out / slen(seg)
    p_start = lerp(a, b, start_t + dt)
    elems.append(("line", p_start[0], p_start[1], b[0], b[1]))

    # Intermediate full segments
    idx = (start_seg_idx + 1) % n
    while idx != end_seg_idx:
        seg = segs[idx]
        if seg_overrides and idx in seg_overrides:
            elems.extend(seg_overrides[idx])
        else:
            elems.append(_seg_to_elem(seg, pts))
        idx = (idx + 1) % n

    # End segment: partial line
    seg = segs[end_seg_idx]
    a, b = pts[seg.start], pts[seg.end]
    dt = R_out / slen(seg)
    p_end = lerp(a, b, end_t - dt)
    elems.append(("line", a[0], a[1], p_end[0], p_end[1]))

    return elems


def _rev_elem(e):
    """Reverse a single element."""
    if e[0] == "line":
        return ("line", e[3], e[4], e[1], e[2])
    _, cx, cy, r, a1, a2 = e
    return ("arc", cx, cy, r, a2, a1)


def _rev_elems(elems):
    """Reverse element list and each element's direction."""
    return [_rev_elem(e) for e in reversed(elems)]


# ── U-turn centerline ────────────────────────────────────────

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


def _uturn_center_elems(pts, outline_segs, inner_segs, seg_idx, t_param,
                         side, shell_t, R_in):
    """Compute centerline U-turn elements: F-arc, crossover, W-arc.

    Returns list of 3 elements (F-side arc, crossover line, W-side arc).
    These go from F-side to W-side. Reverse for W→F direction.
    """
    R_out = R_in + shell_t
    R_center = R_in + shell_t / 2.0

    seg = outline_segs[seg_idx]
    F_pt = lerp(pts[seg.start], pts[seg.end], t_param)
    iseg = inner_segs[seg_idx]
    W_pt = lerp(pts[iseg.start], pts[iseg.end], t_param)

    dx = pts[seg.end][0] - pts[seg.start][0]
    dy = pts[seg.end][1] - pts[seg.start][1]
    ln = math.sqrt(dx * dx + dy * dy)
    t_hat = (dx / ln, dy / ln)
    n_ext = (-t_hat[1], t_hat[0])

    if side == "start":
        open_dir = t_hat
    else:
        open_dir = (-t_hat[0], -t_hat[1])
    wall_dir = (-open_dir[0], -open_dir[1])
    n_int = (-n_ext[0], -n_ext[1])

    oc = (F_pt[0] - R_out * n_ext[0] + R_out * wall_dir[0],
          F_pt[1] - R_out * n_ext[1] + R_out * wall_dir[1])
    ic = (W_pt[0] + R_out * n_ext[0] + R_out * wall_dir[0],
          W_pt[1] + R_out * n_ext[1] + R_out * wall_dir[1])

    f_arc, _, f_ep = _qarc_elem(oc[0], oc[1], R_center, n_ext, open_dir)
    w_arc, w_sp, _ = _qarc_elem(ic[0], ic[1], R_center, open_dir, n_int)

    crossover = ("line", f_ep[0], f_ep[1], w_sp[0], w_sp[1])

    return [f_arc, crossover, w_arc]


# ── section assembly ──────────────────────────────────────────

def _build_tpath(pts, outline_segs, inner_segs, tf_segs, tw_segs,
                 start_op, end_op, shell_t, R_in, tw_ov):
    """Build T-path element list for one wall section."""
    R_out = R_in + shell_t

    # F-side run: trace TF segments between openings (building CW)
    f_elems = _trace_elems(pts, tf_segs, start_op.seg_idx, start_op.t_end,
                           end_op.seg_idx, end_op.t_start, R_out)

    # End U-turn (at end_op, side="start": opening ahead in CW direction)
    end_uturn = _uturn_center_elems(pts, outline_segs, inner_segs,
                                     end_op.seg_idx, end_op.t_start,
                                     "start", shell_t, R_in)

    # W-side run: trace TW segments, then reverse for T-path CW direction
    w_elems = _trace_elems(pts, tw_segs, start_op.seg_idx, start_op.t_end,
                           end_op.seg_idx, end_op.t_start, R_out,
                           seg_overrides=tw_ov)
    w_rev = _rev_elems(w_elems)

    # Start U-turn (at start_op, side="end": opening behind)
    # Computed F→W then reversed to get W→F direction
    start_uturn_base = _uturn_center_elems(pts, outline_segs, inner_segs,
                                            start_op.seg_idx, start_op.t_end,
                                            "end", shell_t, R_in)
    start_uturn = _rev_elems(start_uturn_base)

    return f_elems + end_uturn + w_rev + start_uturn


# ── OpenSCAD output ───────────────────────────────────────────

def _fmt_ft_in(ft):
    """Format feet as ft' inches\" with 4 decimal places on inches."""
    total_in = ft * 12
    whole_ft = int(total_in // 12)
    remaining_in = total_in - whole_ft * 12
    return f"{whole_ft}' {remaining_in:.4f}\""


def _scad_seg(elem):
    """Format a T-path element as SCAD array literal."""
    if elem[0] == "line":
        _, x1, y1, x2, y2 = elem
        return f"[0, {x1:.8f}, {y1:.8f}, {x2:.8f}, {y2:.8f}]"
    _, cx, cy, r, a1, a2 = elem
    return f"[1, {cx:.8f}, {cy:.8f}, {r:.8f}, {a1:.6f}, {a2:.6f}]"


def _seg_comment(elem):
    """Generate an inline comment for a T-path element."""
    if elem[0] == "line":
        _, x1, y1, x2, y2 = elem
        dE, dN = x2 - x1, y2 - y1
        length = math.sqrt(dE * dE + dN * dN)
        bearing = math.degrees(math.atan2(dE, dN)) % 360
        return f"// {_fmt_ft_in(length)} @ {bearing:.4f}deg"
    _, cx, cy, r, a1, a2 = elem
    sweep = a2 - a1
    direction = "CCW" if sweep > 0 else "CW"
    return f"// R {_fmt_ft_in(r)}, {direction} {abs(sweep):.4f}deg"


def generate():
    fp = build_floorplan_data()
    pts = fp.pts
    outline_segs = fp.outline_segs
    inner_segs = fp.inner_segs
    openings = fp.openings
    radii = fp.radii

    shell_t = SHELL_THICKNESS
    shell_half = shell_t / 2.0
    R_in = OPENING_INSIDE_RADIUS

    # Compute TF-series (shell centerline, shell_half from F-face)
    tf_pts, tf_segs = compute_inset_path(outline_segs, pts, radii,
                                          shell_half, "TF")
    pts.update(tf_pts)

    # Compute TW-series (shell centerline, shell_half from W-face)
    tw_pts, tw_segs = compute_inset_path(outline_segs, pts, radii,
                                          WALL_OUTER - shell_half, "TW")
    pts.update(tw_pts)

    # F8-F9 override for TW-side only (straight-arc-straight at 7" inset)
    tw_ov = {7: _f8f9_elems(pts, WALL_OUTER - shell_half,
                             OPENING_INSIDE_RADIUS + shell_half)}

    sections = enumerate_wall_sections(openings, outline_segs)
    sections = sections[-1:] + sections[:-1]

    section_data = []
    for start_op, end_op in sections:
        tpath = _build_tpath(
            pts, outline_segs, inner_segs, tf_segs, tw_segs,
            start_op, end_op, shell_t, R_in, tw_ov)
        label = f"{start_op.name}_{end_op.name}"
        section_data.append((label, tpath))

    out = []
    out.append("// flat_roof.scad - T-path shell centerline extrusion")
    out.append(f"// Wall height: {WALL_HEIGHT_IN:.0f}\" ({WALL_HEIGHT_FT:.4f} ft)")
    out.append("// Construction: 2\" outer shell / 4\" air gap / 2\" inner shell")
    out.append(f"// {len(section_data)} wall sections between openings")
    out.append("// Units: feet")
    out.append("// Generated by scad/gen_flat_roof.py")
    out.append("")
    out.append("// --- Helper functions ---")
    out.append("// Segment types: [0, x1, y1, x2, y2] = line")
    out.append("//                 [1, cx, cy, r, a1, a2] = arc (degrees)")
    out.append("")
    out.append("function arc_n(sweep) = max(4, round(abs(sweep) / 3));")
    out.append("")
    out.append("function line_pts(seg, d) =")
    out.append("  let(dx = seg[3]-seg[1], dy = seg[4]-seg[2],")
    out.append("      len = sqrt(dx*dx + dy*dy),")
    out.append("      nx = -dy/len, ny = dx/len)")
    out.append("  [[seg[1]+d*nx, seg[2]+d*ny], [seg[3]+d*nx, seg[4]+d*ny]];")
    out.append("")
    out.append("function arc_off_pts(seg, d) =")
    out.append("  let(cx = seg[1], cy = seg[2], r = seg[3],")
    out.append("      a1 = seg[4], a2 = seg[5],")
    out.append("      r_off = a2 > a1 ? r - d : r + d,")
    out.append("      n = arc_n(a2 - a1))")
    out.append("  [for (i=[0:n]) [cx + r_off*cos(a1 + i*(a2-a1)/n),")
    out.append("                   cy + r_off*sin(a1 + i*(a2-a1)/n)]];")
    out.append("")
    out.append("function seg_pts(seg, d) =")
    out.append("  seg[0] == 0 ? line_pts(seg, d) : arc_off_pts(seg, d);")
    out.append("")
    out.append("function tail(v) = len(v) > 1")
    out.append("  ? [for (i=[1:len(v)-1]) v[i]] : [];")
    out.append("")
    out.append("function shell_pts(path, d) =")
    out.append("  [for (i=[0:len(path)-1]) each")
    out.append("    (i == 0 ? seg_pts(path[i], d) : tail(seg_pts(path[i], d)))];")
    out.append("")
    out.append("module wall_shell(path, d) {")
    out.append("  ext = shell_pts(path, -d);")
    out.append("  gap = shell_pts(path, d);")
    out.append("  n_ext = len(ext);")
    out.append("  polygon(")
    out.append("    points = concat(ext, gap),")
    out.append("    paths = [")
    out.append("      [for (i=[0:n_ext-1]) i],")
    out.append("      [for (i=[0:len(gap)-1]) n_ext + i]")
    out.append("    ]")
    out.append("  );")
    out.append("}")
    out.append("")
    out.append(f"half_t = {shell_half:.6f};")
    out.append(f"wall_height = {WALL_HEIGHT_FT:.6f};")
    out.append("")

    # Pre-compute data parts to find alignment width
    all_data_parts = []
    for _, tpath in section_data:
        parts = []
        for i, elem in enumerate(tpath):
            comma = "," if i < len(tpath) - 1 else ""
            parts.append(f"  {_scad_seg(elem)}{comma}")
        all_data_parts.append(parts)
    max_data_width = max(len(p) for parts in all_data_parts for p in parts)

    for (label, tpath), parts in zip(section_data, all_data_parts):
        out.append(f"// T-path: wall section {label.replace('_', ' to ')}")
        out.append(f"t_{label} = [")
        for elem, data_part in zip(tpath, parts):
            pad = max(1, max_data_width + 1 - len(data_part))
            out.append(f"{data_part}{' ' * pad}{_seg_comment(elem)}")
        out.append("];")
        out.append("")

    out.append("// --- Assembly ---")
    out.append("union() {")
    for label, _ in section_data:
        out.append(f"  linear_extrude(height = wall_height)")
        out.append(f"    wall_shell(t_{label}, half_t);")
    out.append("}")
    out.append("")

    with open(_OUT, "w") as f:
        f.write("\n".join(out))
    print(f"wrote {_OUT}")


if __name__ == "__main__":
    generate()
