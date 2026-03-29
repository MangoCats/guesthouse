"""Generate 2in12.scad — T-path shell centerline approach, 2:12 roof slope.

Each wall section has a single T-path tracing the center of the shell
material through both shells. OpenSCAD offsets this path by ±shell_t/2
to produce the shell boundary polygon.
All coordinates in feet.
"""
import os
import sys

# Ensure project root is on sys.path for package imports
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

from floorplan.constants import (WALL_OUTER, SHELL_THICKNESS,
                                  OPENING_INSIDE_RADIUS, OPENING_HEIGHT,
                                  ROOF_MIN_THICK, SHED_ROOF_SLOPE,
                                  SHED_ROOF_EAVE_ELEV, SEAM_SPACING,
                                  SEAM_WIDTH, SEAM_HEIGHT)
from floorplan.roof import compute_roof_geometry, roof_polyline, roof_segments
from shared.wall_shells import compute_inset_path
from scad._common import (seg_to_elem as _seg_to_elem, f8f9_elems as _f8f9_elems,
                          trace_elems as _trace_elems, rev_elems as _rev_elems,
                          uturn_center_elems as _uturn_center_elems,
                          build_tpath as _build_tpath,
                          fmt_ft_in as _fmt_ft_in, scad_seg as _scad_seg,
                          seg_comment as _seg_comment,
                          window_panel_poly as _window_panel_poly,
                          compute_wall_bands as _compute_wall_bands)

_DIR = os.path.dirname(os.path.abspath(__file__))
_OUT = os.path.join(_DIR, "2in12.scad")

# Fallback values (used when not in gd.constants; match floorplan/constants.py)
WALL_HEIGHT_FT = OPENING_HEIGHT  # for test introspection
ROOF_SLOPE = SHED_ROOF_SLOPE     # for test introspection



# ── OpenSCAD output ───────────────────────────────────────────


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

    _consts = gd.constants if gd is not None else {}
    shell_t = _consts.get('SHELL_THICKNESS', SHELL_THICKNESS)
    shell_half = shell_t / 2.0
    R_in = _consts.get('OPENING_INSIDE_RADIUS', OPENING_INSIDE_RADIUS)
    wall_outer = _consts.get('WALL_OUTER', WALL_OUTER)

    opening_h = _consts.get('OPENING_HEIGHT', OPENING_HEIGHT)
    roof_min_thick = _consts.get('ROOF_MIN_THICK', ROOF_MIN_THICK)
    roof_slope = _consts.get('SHED_ROOF_SLOPE', SHED_ROOF_SLOPE)
    roof_eave_elev = _consts.get('SHED_ROOF_EAVE_ELEV', SHED_ROOF_EAVE_ELEV)
    seam_spacing = _consts.get('SEAM_SPACING', SEAM_SPACING)
    seam_w = _consts.get('SEAM_WIDTH', SEAM_WIDTH)
    seam_h = _consts.get('SEAM_HEIGHT', SEAM_HEIGHT)

    # Compute TF-series (shell centerline, shell_half from F-face)
    tf_pts, tf_segs = compute_inset_path(outline_segs, pts, radii,
                                          shell_half, "TF")
    pts.update(tf_pts)

    # Compute TW-series (shell centerline, shell_half from W-face)
    tw_pts, tw_segs = compute_inset_path(outline_segs, pts, radii,
                                          wall_outer - shell_half, "TW")
    pts.update(tw_pts)

    # F8-F9 override for TW-side only (straight-arc-straight at 7" inset)
    _f8f9_idx = next((i for i, s in enumerate(outline_segs)
                      if s.start == "F8" and s.end == "F9"), None)
    tw_ov = {_f8f9_idx: _f8f9_elems(pts, wall_outer - shell_half,
                             R_in + shell_half)} if _f8f9_idx is not None else {}

    # Roof outline: T-path arc/line segments + sampled polygon for bounds
    # DB-driven when roof corners are configured; F-series fallback otherwise.
    _rc = getattr(gd, '_roof_corners_data', None) if gd is not None else None
    if _rc and _rc.get("corners"):
        from shared.roof_outline import db_roof_segments
        _corner_names = [c["center"] for c in _rc["corners"]]
        _corner_rad    = [c.get("radiused", False) for c in _rc["corners"]]
        _overhang      = float(_rc.get("overhang", 0.5))
        roof_elems = db_roof_segments(_corner_names, _corner_rad, pts, radii, _overhang)
        roof_pts = list(gd.roof_poly) if getattr(gd, 'roof_poly', None) else []
    elif gd is not None and getattr(gd, 'roof', None) and hasattr(gd.roof, 'r3_radius'):
        roof_elems = roof_segments(gd.roof)
        roof_pts = list(gd.roof_poly)
    else:
        roof_geo = compute_roof_geometry(pts, radii)
        roof_elems = roof_segments(roof_geo)
        roof_pts = roof_polyline(roof_geo)
    if not roof_pts:
        roof_pts = roof_polyline(compute_roof_geometry(pts, radii))
    min_roof_y = min(y for _, y in roof_pts)
    max_roof_y = max(y for _, y in roof_pts)
    min_roof_x = min(x for x, _ in roof_pts)
    max_roof_x = max(x for x, _ in roof_pts)
    # Reference elevation at south edge of roof outline (R19-R01 = low side of shed roof)
    ref_y = min_roof_y
    roof_z_offset = roof_eave_elev - roof_slope * ref_y
    max_roof_z = roof_slope * max_roof_y + roof_z_offset
    max_upper_h = max_roof_z - opening_h  # generous upper wall height
    # Compute intersection cube bounds from roof_elems (includes overhang).
    # For arcs use center ± radius as a conservative bound.
    _elem_xs = []
    _elem_ys = []
    for _e in roof_elems:
        if _e[0] == "line":
            _elem_xs += [_e[1], _e[3]]; _elem_ys += [_e[2], _e[4]]
        else:
            _elem_xs += [_e[1] - _e[3], _e[1] + _e[3]]
            _elem_ys += [_e[2] - _e[3], _e[2] + _e[3]]
    _margin = 2.0
    _cube_x_min = min(_elem_xs) - _margin
    _cube_x_max = max(_elem_xs) + _margin
    _cube_y_min = min(_elem_ys) - _margin
    _cube_y_max = max(_elem_ys) + _margin
    _cube_z_depth = 30.0  # generous depth in shear space

    # N-band walls: driven by per-opening bottom/top elevations
    bands = _compute_wall_bands(
        openings, outline_segs, max_roof_z,
        pts, tf_segs, tw_segs, inner_segs, shell_t, R_in, tw_ov)

    # Window panels: non-door openings with per-opening elevations
    panel_half = 0.5 / 12.0  # 1" thick panel, half-thickness in feet
    window_panels = [
        (op.name, _window_panel_poly(op, outline_segs, inner_segs, pts, panel_half),
         op.bottom_elev, op.top_elev - op.bottom_elev)
        for op in openings if op.opening_type != "door"
    ]

    out = []
    opening_in = opening_h * 12.0
    roof_min_in = roof_min_thick * 12.0
    out.append("// 2in12.scad - T-path shell centerline extrusion (2:12 slope)")
    out.append(f"// {len(bands)} elevation band(s): openings by per-opening bottom/top elevation")
    out.append(f"// Upper wall:   {opening_in:.0f}\" to sloped roof underside (full perimeter)")
    out.append(f"// Roof: {roof_min_in:.0f}\" slab, {roof_slope * 12:.0f}:12 slope N, "
               f"{roof_eave_elev:.1f}' at R19-R01")
    out.append("// Construction: 2\" outer shell / 4\" air gap / 2\" inner shell")
    out.append("// Units: feet")
    out.append("// Generated by scad/gen_2in12.py")
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
    out.append(f"roof_thick = {roof_min_thick:.6f};")
    out.append(f"seam_spacing = {seam_spacing:.6f};  // {seam_spacing * 12:.0f}\" on center")
    out.append(f"seam_w = {seam_w:.6f};  // {seam_w * 12:.0f}\" wide")
    out.append(f"seam_h = {seam_h:.6f};  // {seam_h * 12:.1f}\" tall")
    out.append(f"roof_slope = {roof_slope:.8f};  // {roof_slope * 12:.4f}\" per ft ({roof_slope * 12:.0f}:12)")
    out.append(f"roof_z_off = {roof_z_offset:.8f};")
    out.append("roof_shear = [[1,0,0,0], [0,1,0,0],")
    out.append("              [0, roof_slope, 1, roof_z_off], [0,0,0,1]];")
    out.append("")

    # Collect unique T-paths from bands: "full" label → emitted as t_full_upper
    seen_labels: dict[str, list] = {}
    for _z1, _z2, sec_data in bands:
        for label, tpath in sec_data:
            scad_label = "full_upper" if label == "full" else label
            if scad_label not in seen_labels:
                seen_labels[scad_label] = tpath
    all_entries = list(seen_labels.items())

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

    # Emit T-path data arrays
    for (scad_label, tpath), parts in zip(all_entries, all_data_parts):
        out.append(f"// T-path: {scad_label.replace('_', ' ')}")
        out.append(f"t_{scad_label} = [")
        for elem, data_part in zip(tpath, parts):
            pad = max(1, max_data_width + 1 - len(data_part))
            out.append(f"{data_part}{' ' * pad}{_seg_comment(elem)}")
        out.append("];")
        out.append("")

    # Emit roof outline as arc/line T-path segments
    out.append("// Roof outline (DB-driven, south edge = R19-R01 low side of shed roof)")
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
    for bi, (z1, z2, sec_data) in enumerate(bands):
        z1_in = z1 * 12.0
        z2_in = z2 * 12.0
        h = z2 - z1
        is_upper = all(lbl == "full" for lbl, _ in sec_data)
        active_names = [lbl for lbl, _ in sec_data if lbl != "full"]
        desc = "full perimeter" if not active_names else f"openings: {', '.join(active_names[:3])}"
        out.append(f"  // Band {bi}: {z1_in:.0f}\" to {z2_in:.0f}\" — {desc}")
        if is_upper:
            out.append("  render() intersection() {")
            for label, _ in sec_data:
                scad_label = "full_upper" if label == "full" else label
                if abs(z1) < 1e-9:
                    out.append(f"    linear_extrude(height = {h:.6f}, convexity = 10)")
                else:
                    out.append(f"    translate([0, 0, {z1:.6f}])")
                    out.append(f"      linear_extrude(height = {h:.6f}, convexity = 10)")
                out.append(f"        wall_shell(t_{scad_label}, half_t);")
            out.append("    multmatrix(roof_shear)")
            out.append(f"      translate([{_cube_x_min:.4f}, {_cube_y_min:.4f}, -{_cube_z_depth:.1f}])")
            out.append(f"        cube([{_cube_x_max - _cube_x_min:.4f}, {_cube_y_max - _cube_y_min:.4f}, {_cube_z_depth:.1f}]);")
            out.append("  }")
        else:
            for label, _ in sec_data:
                scad_label = "full_upper" if label == "full" else label
                if abs(z1) < 1e-9:
                    out.append(f"  linear_extrude(height = {h:.6f}, convexity = 10)")
                else:
                    out.append(f"  translate([0, 0, {z1:.6f}])")
                    out.append(f"    linear_extrude(height = {h:.6f}, convexity = 10)")
                out.append(f"      wall_shell(t_{scad_label}, half_t);")
    out.append("}")
    out.append(f"// Sloped roof slab ({roof_min_in:.0f}\", {roof_slope * 12:.0f}:12 slope N, "
               f"{roof_eave_elev:.1f}' at R19-R01)")
    out.append("color(roof_teal) {")
    out.append("  multmatrix(roof_shear)")
    out.append("    linear_extrude(height = roof_thick, convexity = 10)")
    out.append("      polygon(points = shell_pts(roof_outline, 0));")
    out.append(f"  // Standing seam ribs ({seam_spacing * 12:.0f}\" o.c., "
               f"{seam_w * 12:.0f}\" wide, {seam_h * 12:.1f}\" tall)")
    out.append("  multmatrix(roof_shear)")
    out.append(f"    for (x = [{min_roof_x + seam_spacing / 2:.6f} "
               f": seam_spacing : {max_roof_x:.6f}])")
    out.append("      intersection() {")
    out.append(f"        translate([x - seam_w/2, {min_roof_y:.6f}, roof_thick])")
    out.append(f"          cube([seam_w, {max_roof_y - min_roof_y:.6f}, seam_h]);")
    out.append("        linear_extrude(height = roof_thick + seam_h + 0.01, convexity = 10)")
    out.append("          polygon(points = shell_pts(roof_outline, 0));")
    out.append("      }")
    out.append("}")
    out.append("// Window panels (1\" opaque, per-opening elevation range)")
    out.append("window_blue_grey = [0.80, 0.84, 0.90];")
    out.append("color(window_blue_grey) {")
    for name, poly, bot, h in window_panels:
        pts_str = ", ".join(f"[{p[0]:.8f}, {p[1]:.8f}]" for p in poly)
        out.append(f"  // {name}")
        if abs(bot) < 1e-9:
            out.append(f"  linear_extrude(height = {h:.6f}, convexity = 10)")
        else:
            out.append(f"  translate([0, 0, {bot:.6f}])")
            out.append(f"    linear_extrude(height = {h:.6f}, convexity = 10)")
        out.append(f"      polygon(points = [{pts_str}]);")
    out.append("}")
    out.append("")

    with open(_OUT, "w") as f:
        f.write("\n".join(out))
    print(f"wrote {_OUT}")


if __name__ == "__main__":
    generate()
