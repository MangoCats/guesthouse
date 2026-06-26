"""Generate split2.scad — split single-pitch roof (MarkZ 'split2' style).

The roof is divided by a north-south seam (the centreline of the short N-S
exterior wall east of the office / half-bath) into two planar 6" sections:

* West of the seam: the existing 2:12 shed plane (z = roof_slope*y + roof_z_off).
* East of the seam: a single tilted plane (z = east_a + east_b*x + east_c*y)
  that coincides with the west plane along the seam (rising north at 2:12) and
  meets the closet's diagonal south wall face at a constant ("uniform") low
  elevation.

The upper wall is clipped to the lower envelope of the two planes (an
intersection with both shear half-spaces).  All roof sections are 6" thick.

Modelled on scad/gen_2in12.py; all coordinates in feet.
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

from floorplan.constants import (WALL_OUTER, SHELL_THICKNESS,
                                  OPENING_INSIDE_RADIUS, OPENING_HEIGHT,
                                  ROOF_MIN_THICK, SHED_ROOF_SLOPE,
                                  SHED_ROOF_EAVE_ELEV, SEAM_SPACING,
                                  SEAM_WIDTH, SEAM_HEIGHT)
from floorplan.roof import compute_roof_geometry, roof_polyline, roof_segments
from shared.wall_shells import compute_inset_path
from shared.split_roof import compute_east_plane, split_polygon_x
from scad._common import (f8f9_elems as _f8f9_elems,
                          scad_seg as _scad_seg, seg_comment as _seg_comment,
                          window_panel_poly as _window_panel_poly,
                          compute_wall_bands as _compute_wall_bands,
                          build_tw_overrides as _build_tw_overrides,
                          interior_wall_scad as _interior_wall_scad)

_DIR = os.path.dirname(os.path.abspath(__file__))
_OUT = os.path.join(_DIR, "split2.scad")

# This roof style makes every roof section 6" thick (overrides ROOF_MIN_THICK).
SPLIT2_ROOF_THICK = 6.0 / 12.0

WALL_HEIGHT_FT = OPENING_HEIGHT  # for test introspection
ROOF_SLOPE = SHED_ROOF_SLOPE     # for test introspection


def _roof_elems_and_poly(gd, pts, radii):
    """Roof outline arc/line elements + sampled polygon (DB-driven if present)."""
    _rc = getattr(gd, '_roof_corners_data', None) if gd is not None else None
    if _rc and _rc.get("corners"):
        from shared.roof_outline import db_roof_segments
        names = [c["center"] for c in _rc["corners"]]
        crad  = [c.get("radiused", False) for c in _rc["corners"]]
        csc   = [c.get("shortcut", False) for c in _rc["corners"]]
        oh    = float(_rc.get("overhang", 0.5))
        elems = db_roof_segments(names, crad, pts, radii, oh, corner_shortcut=csc)
        poly = list(gd.roof_poly) if getattr(gd, 'roof_poly', None) else []
    elif gd is not None and getattr(gd, 'roof', None) and hasattr(gd.roof, 'r3_radius'):
        elems = roof_segments(gd.roof)
        poly = list(gd.roof_poly)
    else:
        geo = compute_roof_geometry(pts, radii)
        elems = roof_segments(geo)
        poly = roof_polyline(geo)
    if not poly:
        poly = roof_polyline(compute_roof_geometry(pts, radii))
    return elems, poly


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
    roof_min_thick = SPLIT2_ROOF_THICK   # 6" — all sections in this style
    roof_slope = _consts.get('SHED_ROOF_SLOPE', SHED_ROOF_SLOPE)
    roof_eave_elev = _consts.get('SHED_ROOF_EAVE_ELEV', SHED_ROOF_EAVE_ELEV)
    seam_spacing = _consts.get('SEAM_SPACING', SEAM_SPACING)
    seam_w = _consts.get('SEAM_WIDTH', SEAM_WIDTH)
    seam_h = _consts.get('SEAM_HEIGHT', SEAM_HEIGHT)

    # Shell centerline paths (same as 2in12)
    tf_pts, tf_segs = compute_inset_path(outline_segs, pts, radii, shell_half, "TF")
    pts.update(tf_pts)
    tw_pts, tw_segs = compute_inset_path(outline_segs, pts, radii,
                                          wall_outer - shell_half, "TW")
    pts.update(tw_pts)

    _f8f9_idx = next((i for i, s in enumerate(outline_segs)
                      if s.start == "F8" and s.end == "F9"), None)
    tw_ov = {_f8f9_idx: _f8f9_elems(pts, wall_outer - shell_half,
                             R_in + shell_half)} if _f8f9_idx is not None else {}
    if gd is not None and gd.wall_overrides:
        _wov = _build_tw_overrides(gd.wall_overrides, inner_segs, pts)
        for _k, _v in _wov.items():
            tw_ov.setdefault(_k, _v)

    roof_elems, roof_pts = _roof_elems_and_poly(gd, pts, radii)
    min_roof_y = min(y for _, y in roof_pts)
    max_roof_y = max(y for _, y in roof_pts)
    min_roof_x = min(x for x, _ in roof_pts)
    max_roof_x = max(x for x, _ in roof_pts)

    # West plane (2:12) reference, exactly as the existing 2in12 roof.
    ref_y = min_roof_y
    roof_z_offset = roof_eave_elev - roof_slope * ref_y

    # East plane: solved from the seam + closet south wall (MarkZ point names).
    east = compute_east_plane(pts, slope=roof_slope, eave_elev=roof_eave_elev,
                              ref_y=ref_y, wall_outer=wall_outer)

    # Upper wall extrudes to the higher of the two planes' max, then is clipped.
    max_roof_z = roof_slope * max_roof_y + roof_z_offset
    # Cube bounds for the shear-clip half-spaces (cover the whole roof + margin).
    _elem_xs, _elem_ys = [], []
    for _e in roof_elems:
        if _e[0] == "line":
            _elem_xs += [_e[1], _e[3]]; _elem_ys += [_e[2], _e[4]]
        else:
            _elem_xs += [_e[1] - _e[3], _e[1] + _e[3]]
            _elem_ys += [_e[2] - _e[3], _e[2] + _e[3]]
    _margin = 2.0
    _cx0, _cx1 = min(_elem_xs) - _margin, max(_elem_xs) + _margin
    _cy0, _cy1 = min(_elem_ys) - _margin, max(_elem_ys) + _margin
    _cz = 30.0

    bands = _compute_wall_bands(
        openings, outline_segs, max_roof_z,
        pts, tf_segs, tw_segs, inner_segs, shell_t, R_in, tw_ov)

    panel_half = 0.5 / 12.0
    window_panels = [
        (op.name, _window_panel_poly(op, outline_segs, inner_segs, pts, panel_half),
         op.bottom_elev, op.top_elev - op.bottom_elev)
        for op in openings if op.opening_type != "door"
    ]

    # Split the sampled roof outline at the seam into west / east polygons.
    west_poly, east_poly = split_polygon_x(roof_pts, east.seam_x)

    out = []
    opening_in = opening_h * 12.0
    roof_min_in = roof_min_thick * 12.0
    out.append("// split2.scad - split single-pitch roof (west 2:12 + tilted east plane)")
    out.append(f"// Seam (N-S) at E={east.seam_x:.4f}; west of seam = 2:12, east = tilted plane")
    out.append(f"// East plane underside: z = {east.a:.5f} + ({east.b:.5f})*x + ({east.c:.5f})*y")
    out.append(f"// Level east low eave (closet S wall): {east.low_elev:.4f} ft")
    out.append(f"// Roof: {roof_min_in:.0f}\" slab, all sections")
    out.append("// Construction: 2\" outer shell / 4\" air gap / 2\" inner shell")
    out.append("// Units: feet — Generated by scad/gen_split2.py")
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
    out.append(f"seam_spacing = {seam_spacing:.6f};")
    out.append(f"seam_w = {seam_w:.6f};")
    out.append(f"seam_h = {seam_h:.6f};")
    out.append(f"roof_slope = {roof_slope:.8f};  // 2:12 N (west plane)")
    out.append(f"roof_z_off = {roof_z_offset:.8f};")
    out.append(f"east_a = {east.a:.8f};")
    out.append(f"east_b = {east.b:.8f};")
    out.append(f"east_c = {east.c:.8f};")
    out.append(f"seam_x = {east.seam_x:.8f};")
    out.append("// West plane shear: maps a level cube top to z = roof_slope*y + roof_z_off")
    out.append("west_shear = [[1,0,0,0], [0,1,0,0],")
    out.append("              [0, roof_slope, 1, roof_z_off], [0,0,0,1]];")
    out.append("// East plane shear: maps a level cube top to z = east_a + east_b*x + east_c*y")
    out.append("east_shear = [[1,0,0,0], [0,1,0,0],")
    out.append("              [east_b, east_c, 1, east_a], [0,0,0,1]];")
    out.append("")

    # Collect unique T-paths from bands
    seen_labels: dict[str, list] = {}
    for _z1, _z2, sec_data in bands:
        for label, tpath in sec_data:
            scad_label = "full_upper" if label == "full" else label
            if scad_label not in seen_labels:
                seen_labels[scad_label] = tpath
    all_entries = list(seen_labels.items())

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

    for (scad_label, tpath), parts in zip(all_entries, all_data_parts):
        out.append(f"// T-path: {scad_label.replace('_', ' ')}")
        out.append(f"t_{scad_label} = [")
        for elem, data_part in zip(tpath, parts):
            pad = max(1, max_data_width + 1 - len(data_part))
            out.append(f"{data_part}{' ' * pad}{_seg_comment(elem)}")
        out.append("];")
        out.append("")

    out.append("// Roof outline (full; split at seam_x for the two planes)")
    out.append("roof_outline = [")
    for elem, data_part in zip(roof_elems, roof_parts):
        pad = max(1, max_data_width + 1 - len(data_part))
        out.append(f"{data_part}{' ' * pad}{_seg_comment(elem)}")
    out.append("];")
    out.append("")

    def _poly_literal(poly):
        return "[" + ", ".join(f"[{p[0]:.6f}, {p[1]:.6f}]" for p in poly) + "]"

    out.append(f"roof_poly_west = {_poly_literal(west_poly)};")
    out.append(f"roof_poly_east = {_poly_literal(east_poly)};")
    out.append("")

    out.append("// --- Assembly ---")
    out.append("wall_cream = [0.88, 0.82, 0.60];")
    out.append("roof_teal = [0.10, 0.35, 0.33];")
    out.append("color(wall_cream) union() {")
    for bi, (z1, z2, sec_data) in enumerate(bands):
        z1_in = z1 * 12.0
        z2_in = z2 * 12.0
        h = z2 - z1
        is_upper = all(lbl == "full" for lbl, _ in sec_data)
        active = [lbl for lbl, _ in sec_data if lbl != "full"]
        desc = "full perimeter" if not active else f"openings: {', '.join(active[:3])}"
        out.append(f"  // Band {bi}: {z1_in:.0f}\" to {z2_in:.0f}\" — {desc}")
        if is_upper:
            # Upper wall clipped to the LOWER envelope of both planes:
            # intersection with both shear half-spaces.
            out.append("  render() intersection() {")
            for label, _ in sec_data:
                scad_label = "full_upper" if label == "full" else label
                if abs(z1) < 1e-9:
                    out.append(f"    linear_extrude(height = {h:.6f}, convexity = 10)")
                else:
                    out.append(f"    translate([0, 0, {z1:.6f}])")
                    out.append(f"      linear_extrude(height = {h:.6f}, convexity = 10)")
                out.append(f"        wall_shell(t_{scad_label}, half_t);")
            out.append("    multmatrix(west_shear)")
            out.append(f"      translate([{_cx0:.4f}, {_cy0:.4f}, -{_cz:.1f}])")
            out.append(f"        cube([{_cx1 - _cx0:.4f}, {_cy1 - _cy0:.4f}, {_cz:.1f}]);")
            out.append("    multmatrix(east_shear)")
            out.append(f"      translate([{_cx0:.4f}, {_cy0:.4f}, -{_cz:.1f}])")
            out.append(f"        cube([{_cx1 - _cx0:.4f}, {_cy1 - _cy0:.4f}, {_cz:.1f}]);")
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

    out.append(f"// Roof slabs ({roof_min_in:.0f}\" thick): west 2:12 plane, east tilted plane")
    out.append("color(roof_teal) {")
    out.append("  // West section (2:12), clipped to E <= seam_x")
    out.append("  multmatrix(west_shear)")
    out.append("    linear_extrude(height = roof_thick, convexity = 10)")
    out.append("      polygon(points = roof_poly_west);")
    out.append("  // East section (tilted plane), clipped to E >= seam_x")
    out.append("  multmatrix(east_shear)")
    out.append("    linear_extrude(height = roof_thick, convexity = 10)")
    out.append("      polygon(points = roof_poly_east);")
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

    # Interior walls: top follows the lower envelope of the two planes.
    iw_polys = getattr(gd, 'iw_polys', None) if gd is not None else None
    ro_polys = getattr(gd, 'ro_polys', None) if gd is not None else None
    if iw_polys:
        def _iw_top(name):
            poly = iw_polys[name]
            ax = sum(p[0] for p in poly) / len(poly)
            ay = sum(p[1] for p in poly) / len(poly)
            return min(roof_z_offset + roof_slope * ay, east.z(ax, ay))
        _interior_wall_scad(out, iw_polys, ro_polys, get_wall_top=_iw_top)

    with open(_OUT, "w", encoding="utf-8") as f:
        f.write("\n".join(out))
    print(f"wrote {_OUT}")


if __name__ == "__main__":
    generate()
