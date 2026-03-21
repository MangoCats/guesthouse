"""Generate flat_roof.scad — T-path shell centerline approach.

Each wall section has a single T-path tracing the center of the shell
material through both shells. OpenSCAD offsets this path by ±shell_t/2
to produce the shell boundary polygon.
All coordinates in feet.
"""
import math
import os
import sys

# Ensure project root is on sys.path for package imports
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

from floorplan.constants import (WALL_OUTER, SHELL_THICKNESS,
                                  OPENING_INSIDE_RADIUS)
from floorplan.roof import compute_roof_geometry, roof_polyline, roof_segments
from shared.wall_shells import compute_inset_path, enumerate_wall_sections, lerp
from shared.types import ArcSeg
from scad._common import (seg_to_elem as _seg_to_elem, f8f9_elems as _f8f9_elems,
                          trace_elems as _trace_elems, rev_elem as _rev_elem,
                          rev_elems as _rev_elems, qarc_elem as _qarc_elem)

_DIR = os.path.dirname(os.path.abspath(__file__))
_OUT = os.path.join(_DIR, "flat_roof.scad")

LOWER_HEIGHT_IN = 20.0  # 1'8" lower walls (doors O3, O6 only)
LOWER_HEIGHT_FT = LOWER_HEIGHT_IN / 12.0

WALL_HEIGHT_IN = 80.0   # 6'8" middle wall top (all openings)
WALL_HEIGHT_FT = WALL_HEIGHT_IN / 12.0
MIDDLE_HEIGHT_FT = (WALL_HEIGHT_IN - LOWER_HEIGHT_IN) / 12.0

UPPER_BASE_IN = 80.0    # bottom of upper wall band
UPPER_TOP_IN = 112.0    # 9'4" top of upper wall band
UPPER_BASE_FT = UPPER_BASE_IN / 12.0
UPPER_HEIGHT_FT = (UPPER_TOP_IN - UPPER_BASE_IN) / 12.0

ROOF_THICK_IN = 18.0   # minimum roof slab thickness (inches), at south edge
ROOF_THICK_FT = ROOF_THICK_IN / 12.0
ROOF_SLOPE_IN_PER_FT = 0.25   # 1/4" rise per foot of northing
ROOF_SLOPE = ROOF_SLOPE_IN_PER_FT / 12.0  # ft/ft

SEAM_SPACING_IN = 16.0    # standing seam on-center spacing (inches)
SEAM_SPACING_FT = SEAM_SPACING_IN / 12.0
SEAM_WIDTH_IN = 1.0       # seam rib width (inches)
SEAM_WIDTH_FT = SEAM_WIDTH_IN / 12.0
SEAM_HEIGHT_IN = 1.5      # seam rib height (inches)
SEAM_HEIGHT_FT = SEAM_HEIGHT_IN / 12.0



# ── U-turn centerline ────────────────────────────────────────


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
                 start_op, end_op, shell_t, R_in, tw_ov, full_wrap=False):
    """Build T-path element list for one wall section."""
    R_out = R_in + shell_t

    # F-side run: trace TF segments between openings (building CW)
    f_elems = _trace_elems(pts, tf_segs, start_op.seg_idx, start_op.t_end,
                           end_op.seg_idx, end_op.t_start, R_out,
                           full_wrap=full_wrap)

    # End U-turn (at end_op, side="start": opening ahead in CW direction)
    end_uturn = _uturn_center_elems(pts, outline_segs, inner_segs,
                                     end_op.seg_idx, end_op.t_start,
                                     "start", shell_t, R_in)

    # W-side run: trace TW segments, then reverse for T-path CW direction
    w_elems = _trace_elems(pts, tw_segs, start_op.seg_idx, start_op.t_end,
                           end_op.seg_idx, end_op.t_start, R_out,
                           seg_overrides=tw_ov, full_wrap=full_wrap)
    w_rev = _rev_elems(w_elems)

    # Start U-turn (at start_op, side="end": opening behind)
    # Computed F→W then reversed to get W→F direction
    start_uturn_base = _uturn_center_elems(pts, outline_segs, inner_segs,
                                            start_op.seg_idx, start_op.t_end,
                                            "end", shell_t, R_in)
    start_uturn = _rev_elems(start_uturn_base)

    return f_elems + end_uturn + w_rev + start_uturn


# ── OpenSCAD output ───────────────────────────────────────────

def _fmt_ft_in(ft, in_width=8):
    """Format feet as ft' inches\" with 4 decimal places on inches."""
    total_in = ft * 12
    whole_ft = int(total_in // 12)
    remaining_in = total_in - whole_ft * 12
    return f"{whole_ft:2d}' {remaining_in:{in_width}.4f}\""


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
        return f"// {_fmt_ft_in(length)} @ {bearing:8.4f}deg"
    _, cx, cy, r, a1, a2 = elem
    sweep = a2 - a1
    direction = "CCW" if sweep > 0 else "CW "
    return f"// {direction} {abs(sweep):8.4f}deg R {_fmt_ft_in(r, 7)}"


def generate(gd=None):
    if gd is None:
        from floorplan.gen_floorplan import build_floorplan_data
        fp = build_floorplan_data()
        pts = fp.pts
        outline_segs = fp.outline_segs
        inner_segs = fp.inner_segs
        openings = fp.openings
        radii = fp.radii
    else:
        pts = dict(gd.pts)
        outline_segs = gd.outline_segs
        inner_segs = gd.inner_segs
        openings = gd.openings
        radii = gd.radii

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
    _f8f9_idx = next((i for i, s in enumerate(outline_segs)
                      if s.start == "F8" and s.end == "F9"), None)
    tw_ov = {_f8f9_idx: _f8f9_elems(pts, WALL_OUTER - shell_half,
                             OPENING_INSIDE_RADIUS + shell_half)} if _f8f9_idx is not None else {}

    sections = enumerate_wall_sections(openings, outline_segs)
    sections = sections[-1:] + sections[:-1]

    section_data = []
    for start_op, end_op in sections:
        tpath = _build_tpath(
            pts, outline_segs, inner_segs, tf_segs, tw_segs,
            start_op, end_op, shell_t, R_in, tw_ov)
        label = f"{start_op.name}_{end_op.name}"
        section_data.append((label, tpath))

    # Door-only T-path: lower walls (0-20") with only O3, O6
    door_openings = [op for op in openings if op.name in ("O3", "O6")]
    lower_sections = enumerate_wall_sections(door_openings, outline_segs)
    lower_sections = lower_sections[-1:] + lower_sections[:-1]

    lower_section_data = []
    for start_op, end_op in lower_sections:
        tpath = _build_tpath(
            pts, outline_segs, inner_segs, tf_segs, tw_segs,
            start_op, end_op, shell_t, R_in, tw_ov)
        label = f"lower_{start_op.name}_{end_op.name}"
        lower_section_data.append((label, tpath))

    # Window openings: open in middle wall only (not doors O3/O6, not O4)
    window_names = {"O1", "O2", "O5", "O7", "O8", "O8a", "O9", "O10", "O11"}
    window_panels = []
    panel_half = 0.5 / 12.0  # 1" thick panel, half-thickness in feet
    for op in openings:
        if op.name not in window_names:
            continue
        seg = outline_segs[op.seg_idx]
        iseg = inner_segs[op.seg_idx]
        F_A, F_B = pts[seg.start], pts[seg.end]
        W_A, W_B = pts[iseg.start], pts[iseg.end]
        M_A = ((F_A[0] + W_A[0]) / 2, (F_A[1] + W_A[1]) / 2)
        M_B = ((F_B[0] + W_B[0]) / 2, (F_B[1] + W_B[1]) / 2)
        M_s = lerp(M_A, M_B, op.t_start)
        M_e = lerp(M_A, M_B, op.t_end)
        dx, dy = F_B[0] - F_A[0], F_B[1] - F_A[1]
        ln = math.sqrt(dx * dx + dy * dy)
        nx, ny = -dy / ln, dx / ln  # exterior normal
        window_panels.append((op.name, [
            (M_s[0] + nx * panel_half, M_s[1] + ny * panel_half),
            (M_e[0] + nx * panel_half, M_e[1] + ny * panel_half),
            (M_e[0] - nx * panel_half, M_e[1] - ny * panel_half),
            (M_s[0] - nx * panel_half, M_s[1] - ny * panel_half),
        ]))

    # Full-wall T-path: only O4 retained (upper wall band)
    o4_openings = [op for op in openings if op.name == "O4"]
    full_sections = enumerate_wall_sections(o4_openings, outline_segs)
    full_start, full_end = full_sections[0]
    full_tpath = _build_tpath(
        pts, outline_segs, inner_segs, tf_segs, tw_segs,
        full_start, full_end, shell_t, R_in, tw_ov, full_wrap=True)

    # Roof outline polygon and wedge parameters
    roof_geo = compute_roof_geometry(pts, radii)
    roof_pts = roof_polyline(roof_geo)
    roof_y_south = min(y for _, y in roof_pts)
    roof_y_north = max(y for _, y in roof_pts)
    roof_x_west = min(x for x, _ in roof_pts)
    roof_x_east = max(x for x, _ in roof_pts)
    # Flat bottom at wall_height; top slopes up toward north
    # At south edge: thickness = ROOF_THICK_FT (18")
    # At north edge (R3-R4): thickness = ROOF_THICK_FT + slope * delta_y
    roof_z_base = ROOF_THICK_FT - ROOF_SLOPE * roof_y_south
    max_roof_thick = ROOF_THICK_FT + ROOF_SLOPE * (roof_y_north - roof_y_south)
    max_roof_thick_in = max_roof_thick * 12.0

    out = []
    out.append("// flat_roof.scad - T-path shell centerline extrusion")
    out.append(f"// Lower walls:  0 to {LOWER_HEIGHT_IN:.0f}\" (1'8\") — doors O3, O6 only")
    out.append(f"// Middle walls: {LOWER_HEIGHT_IN:.0f}\" to {WALL_HEIGHT_IN:.0f}\" "
               f"(5'0\") — all openings")
    out.append(f"// Upper wall:   {UPPER_BASE_IN:.0f}\" to {UPPER_TOP_IN:.0f}\" "
               f"(9'4\") — O4 only")
    out.append(f"// Roof: {ROOF_THICK_IN:.0f}\"-{max_roof_thick_in:.1f}\" wedge slab "
               f"(1/4\"/ft slope N, {ROOF_THICK_IN:.0f}\" min at south)")
    out.append("// Construction: 2\" outer shell / 4\" air gap / 2\" inner shell")
    out.append(f"// {len(lower_section_data)} lower + {len(section_data)} middle "
               f"+ 1 upper wall sections")
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
    out.append(f"lower_height = {LOWER_HEIGHT_FT:.6f};")
    out.append(f"middle_height = {MIDDLE_HEIGHT_FT:.6f};")
    out.append(f"upper_base = {UPPER_BASE_FT:.6f};")
    out.append(f"upper_height = {UPPER_HEIGHT_FT:.6f};")
    out.append(f"roof_thick = {ROOF_THICK_FT:.6f};")
    out.append(f"max_roof_thick = {max_roof_thick:.6f};")
    out.append(f"seam_spacing = {SEAM_SPACING_FT:.6f};  // {SEAM_SPACING_IN:.0f}\" on center")
    out.append(f"seam_w = {SEAM_WIDTH_FT:.6f};  // {SEAM_WIDTH_IN:.0f}\" wide")
    out.append(f"seam_h = {SEAM_HEIGHT_FT:.6f};  // {SEAM_HEIGHT_IN:.1f}\" tall")
    out.append(f"roof_slope = {ROOF_SLOPE:.8f};  // {ROOF_SLOPE_IN_PER_FT}\" per ft")
    out.append(f"roof_z_base = {roof_z_base:.8f};")
    out.append("")

    # Roof outline as vector/arc segments
    roof_elems = roof_segments(roof_geo)

    # Collect all T-path data (lower + middle + full-wall + roof) for alignment
    all_entries = (list(lower_section_data) + list(section_data)
                   + [("full_O4", full_tpath)])

    # Pre-compute data parts to find alignment width
    all_data_parts = []
    for _, tpath in all_entries:
        parts = []
        for i, elem in enumerate(tpath):
            comma = "," if i < len(tpath) - 1 else ""
            parts.append(f"  {_scad_seg(elem)}{comma}")
        all_data_parts.append(parts)
    roof_parts = []
    for i, elem in enumerate(roof_elems):
        comma = "," if i < len(roof_elems) - 1 else ""
        roof_parts.append(f"  {_scad_seg(elem)}{comma}")
    max_data_width = max(len(p) for parts in all_data_parts + [roof_parts]
                         for p in parts)

    n_lower = len(lower_section_data)
    n_middle = len(section_data)

    # Emit lower wall section T-paths (doors O3, O6 only)
    for (label, tpath), parts in zip(lower_section_data,
                                      all_data_parts[:n_lower]):
        out.append(f"// T-path: lower wall {label.replace('lower_', '').replace('_', ' to ')}")
        out.append(f"t_{label} = [")
        for elem, data_part in zip(tpath, parts):
            pad = max(1, max_data_width + 1 - len(data_part))
            out.append(f"{data_part}{' ' * pad}{_seg_comment(elem)}")
        out.append("];")
        out.append("")

    # Emit middle wall section T-paths (all openings)
    for (label, tpath), parts in zip(section_data,
                                      all_data_parts[n_lower:n_lower + n_middle]):
        out.append(f"// T-path: wall section {label.replace('_', ' to ')}")
        out.append(f"t_{label} = [")
        for elem, data_part in zip(tpath, parts):
            pad = max(1, max_data_width + 1 - len(data_part))
            out.append(f"{data_part}{' ' * pad}{_seg_comment(elem)}")
        out.append("];")
        out.append("")

    # Emit full-wall T-path (upper band, O4 only)
    full_parts = all_data_parts[-1]
    out.append("// T-path: full wall (O4 only)")
    out.append("t_full_O4 = [")
    for elem, data_part in zip(full_tpath, full_parts):
        pad = max(1, max_data_width + 1 - len(data_part))
        out.append(f"{data_part}{' ' * pad}{_seg_comment(elem)}")
    out.append("];")
    out.append("")

    # Emit roof outline as vector/arc segments
    out.append("// Roof outline (R-series, 6\" overhang from F-face)")
    out.append("roof_outline = [")
    for elem, data_part in zip(roof_elems, roof_parts):
        pad = max(1, max_data_width + 1 - len(data_part))
        out.append(f"{data_part}{' ' * pad}{_seg_comment(elem)}")
    out.append("];")
    out.append("")

    out.append("// --- Assembly ---")
    out.append("wall_cream = [0.88, 0.82, 0.60];  // warm cream-yellow (match main house)")
    out.append("roof_teal = [0.10, 0.35, 0.33];  // dark teal-green metal (match main house)")
    out.append("color(wall_cream) union() {")
    out.append(f"  // Lower walls (0 to {LOWER_HEIGHT_IN:.0f}\", doors O3 and O6 only)")
    for label, _ in lower_section_data:
        out.append(f"  linear_extrude(height = lower_height)")
        out.append(f"    wall_shell(t_{label}, half_t);")
    out.append(f"  // Middle walls ({LOWER_HEIGHT_IN:.0f}\" to {WALL_HEIGHT_IN:.0f}\","
               f" all openings)")
    for label, _ in section_data:
        out.append(f"  translate([0, 0, lower_height])")
        out.append(f"    linear_extrude(height = middle_height)")
        out.append(f"      wall_shell(t_{label}, half_t);")
    out.append(f"  // Upper wall ({UPPER_BASE_IN:.0f}\" to {UPPER_TOP_IN:.0f}\", O4 only)")
    out.append("  translate([0, 0, upper_base])")
    out.append("    linear_extrude(height = upper_height)")
    out.append("      wall_shell(t_full_O4, half_t);")
    out.append("}")
    out.append(f"// Wedge roof slab ({ROOF_THICK_IN:.0f}\"-{max_roof_thick_in:.1f}\", "
               f"1/4\"/ft slope N)")
    out.append("color(roof_teal) {")
    out.append(f"  translate([0, 0, upper_base + upper_height])")
    out.append("    render() intersection() {")
    out.append("      linear_extrude(height = max_roof_thick + 0.1)")
    out.append("        polygon(points = shell_pts(roof_outline, 0));")
    out.append("      multmatrix([[1,0,0,0], [0,1,0,0],")
    out.append("                  [0, roof_slope, 1, roof_z_base], [0,0,0,1]])")
    out.append("        translate([-25, -20, -25])")
    out.append("          cube([50, 40, 25]);")
    out.append("    }")
    out.append(f"  // Standing seam ribs ({SEAM_SPACING_IN:.0f}\" o.c., "
               f"{SEAM_WIDTH_IN:.0f}\" wide, {SEAM_HEIGHT_IN:.1f}\" tall)")
    out.append(f"  translate([0, 0, upper_base + upper_height])")
    out.append(f"    for (x = [{roof_x_west + SEAM_SPACING_FT / 2:.6f} "
               f": seam_spacing : {roof_x_east:.6f}])")
    out.append("      intersection() {")
    out.append("        multmatrix([[1,0,0,0], [0,1,0,0],")
    out.append("                    [0, roof_slope, 1, roof_z_base], [0,0,0,1]])")
    out.append(f"          translate([x - seam_w/2, {roof_y_south:.6f}, 0])")
    out.append(f"            cube([seam_w, {roof_y_north - roof_y_south:.6f}, seam_h]);")
    out.append("        linear_extrude(height = max_roof_thick + seam_h + 0.01)")
    out.append("          polygon(points = shell_pts(roof_outline, 0));")
    out.append("      }")
    out.append("}")
    out.append("// Window panels (1\" opaque, middle wall openings only)")
    out.append("window_blue_grey = [0.80, 0.84, 0.90];")
    out.append("color(window_blue_grey) {")
    for name, poly in window_panels:
        pts_str = ", ".join(f"[{p[0]:.8f}, {p[1]:.8f}]" for p in poly)
        out.append(f"  // {name}")
        out.append(f"  translate([0, 0, lower_height])")
        out.append(f"    linear_extrude(height = middle_height)")
        out.append(f"      polygon(points = [{pts_str}]);")
    out.append("}")
    out.append("")

    with open(_OUT, "w") as f:
        f.write("\n".join(out))
    print(f"wrote {_OUT}")


if __name__ == "__main__":
    generate()
