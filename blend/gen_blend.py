"""Generate Blender (.blend) models of the ADU, mirroring the .scad models.

This builds the same geometry the OpenSCAD / glTF pipeline produces — outer
wall shells, the sloped/wedge roof slab (with overhang), glazing panels and
interior walls — but as native Blender mesh objects with materials, saved to
.blend files.

Accuracy: the triangle geometry is reused verbatim from ``scad/gen_gltf.py``
(the same parser + mesh builders that drive ``opts/*.glb``).  The SCAD source
is generated fresh from the live editor DB (``app/adu.db``) for the currently
loaded model — using that model's own roof style — so the Blender model
matches it (and the SCAD/GLB pipeline) to the millimetre.

Coordinate system: SCAD/Blender are both Z-up (east X, north Y, up Z), so —
unlike the glTF export — no axis swap is applied.  Units are feet; the scene
is configured for imperial display (1 Blender unit = 1 foot).

A flat green ground plane covering the survey parcel (the property extent from
the site plans) is added beneath the building.  Its outline is the four parcel
corners from site/gen_site_plan.py, transformed from PDF/survey space back into
the building's FC-origin feet frame.  This is computed in the launching process
(which has PyMuPDF) and passed to headless Blender as a JSON argument, since
Blender's bundled Python lacks fitz.

Outputs: blend/<config_name>.blend (currently-loaded model) and
blend/<config_name>_2in12.blend (2:12 sloped-roof version) — e.g. MarkZ.blend
and MarkZ_2in12.blend.

Usage:
    python blend/gen_blend.py            # locates Blender and runs headless
    blender --background --python blend/gen_blend.py   # legacy: both roof types
"""

import os
import sys
import json
import math

_DIR = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_DIR)

# Legacy fallback (used only when Blender is invoked directly with no manifest):
# (display stem, source .scad file, roof_type for gen_gltf builders)
_MODELS = [
    ("flat_roof", "flat_roof.scad", "flat_roof"),
    ("2in12",     "2in12.scad",     "2in12"),
]

# Material palette — matches scad/gen_gltf.py (RGBA, roughness)
_MATERIALS = {
    "wall": ([0.88, 0.82, 0.60, 1.00], 0.85),  # warm cream-yellow shells
    "roof": ([0.10, 0.35, 0.33, 1.00], 0.55),  # dark teal-green metal slab
    "win":  ([0.80, 0.84, 0.90, 0.18], 0.05),  # blue-grey glazing (very transparent)
    "iw":   ([0.82, 0.72, 0.38, 1.00], 0.90),  # white-pine interior walls
    "grass": ([0.20, 0.50, 0.16, 1.00], 0.95),  # green ground plane (parcel)
    "fence_3rail": ([0.22, 0.12, 0.05, 1.00], 0.90),  # dark brown wood (3-rail)
    "fence_1rail": ([0.48, 0.30, 0.14, 1.00], 0.90),  # lighter brown wood (single-rail)
}

# Ground plane sits a hair below z=0 so it doesn't z-fight the wall base caps.
_GRASS_Z = -0.01  # feet

# ── Fence dimensions (feet) ──
_POST_W = 4.0 / 12.0       # 4" square posts
_RAIL_THICK = 2.0 / 12.0   # 2" rail thickness (perpendicular to fence plane)
_RAIL_H = 4.0 / 12.0       # 4" rail height (vertical)
_POST_SPACING = 10.0       # posts every 10'
_FENCE3_TOP = 4.0 + 8.0 / 12.0  # 3-rail fence top-of-top-rail = 4'8"
_FENCE1_TOP = 4.0               # single-rail fence top-of-top-rail = 4'0"


# ─────────────────────────  geometry (no bpy needed)  ────────────────────────

def _read_scad_file(scad_file):
    """Read SCAD text from the scad/ directory (legacy fallback path)."""
    p = os.path.join(_ROOT, "scad", scad_file)
    with open(p, encoding="utf-8", errors="replace") as f:
        return f.read()


def _build_meshes(scad_text, roof_type):
    """Parse SCAD text → ({key: (verts, indices)}, roof_spec, iw_walls, door_h).

    Walls and glazing reuse scad/gen_gltf.py's mesh builders verbatim (identical
    to the GLB geometry).  The roof and interior walls are returned as data so
    they can be built natively in Blender: the roof as a planar n-gon (gen_gltf's
    ear-clipping leaves a triangular gap on the repeated closing vertex), and the
    interior walls with their door openings cut out (the GLB build omits them).
    """
    if _ROOT not in sys.path:
        sys.path.insert(0, _ROOT)
    from scad import gen_gltf as g

    data = g._parse_scad_text(scad_text)

    wall_mesh = g.MeshBuilder()
    win_mesh  = g.MeshBuilder()

    s = data["scalars"]
    if roof_type == "2in12":
        half_t = s["half_t"]
        for name, z_base, height in data["walls"]:
            if name == "t_full_upper":
                continue  # clipped by roof slope — handled separately
            path = data["_arrays"].get(name)
            if path is None:
                continue
            outer = g.shell_pts(path, -half_t)
            inner = g.shell_pts(path, half_t)
            g.extrude_ring(outer, inner, wall_mesh,
                           z_bot=z_base, z_top=z_base + height)
        g._build_walls_2in12_upper(data, wall_mesh)
    else:
        g._build_walls(data, wall_mesh)

    g._build_windows(data, win_mesh)

    out = {}
    for key, mb in (("wall", wall_mesh), ("win", win_mesh)):
        verts, _norms, indices = mb.to_arrays()
        if verts is not None and len(verts):
            out[key] = (verts, indices)

    iw_walls, door_h = _parse_interior_walls(g, scad_text)
    return out, _roof_spec(g, data, roof_type), iw_walls, door_h


def _parse_interior_walls(g, scad_text):
    """Parse the interior-wall SCAD block into [{poly, z_top, openings}], door_h.

    Each wall is a ``linear_extrude(height=<float>)`` of a polygon; its door
    openings follow as ``linear_extrude(height = iw_door_h + ...)`` polygons
    (inside the same difference()).  Openings are grouped with the most recent
    wall, mirroring the SCAD emission order.
    """
    import ast
    import re

    pos = scad_text.find("// Interior walls")
    if pos < 0:
        return [], 82.5 / 12.0
    block = scad_text[pos:]

    m = re.search(r'iw_door_h\s*=\s*([\d.]+)', block)
    door_h = float(m.group(1)) if m else 82.5 / 12.0

    walls = []
    for em in re.finditer(r'linear_extrude\s*\(\s*height\s*=\s*([^,]+?)\s*,', block):
        height_expr = em.group(1)
        pm = re.search(r'polygon\s*\(\s*points\s*=\s*(\[)', block[em.end():])
        if not pm:
            continue
        start = em.end() + pm.start(1)
        arr = g._bracket_array(block, start)
        if arr is None:
            continue
        pts = ast.literal_eval(re.sub(r'//[^\n]*', '', arr))
        if "iw_door_h" in height_expr:
            if walls:
                walls[-1]["openings"].append(pts)
        else:
            walls.append({"poly": pts, "z_top": float(height_expr),
                          "openings": []})
    return walls, door_h


def _roof_spec(g, data, roof_type):
    """Roof outline polygon + bottom/top z-plane coefficients (z = a + b*x + c*y).

    Both roof caps are planar (the flat roof's wedge top and the 2:12 slab are
    tilted planes; the flat roof's underside is level), so the roof can be a
    clean planar n-gon prism.  Consecutive duplicate outline points — including
    the repeated closing vertex that breaks gen_gltf's ear-clipping — are removed.
    """
    poly = g.shell_pts(data["roof_outline"], 0)
    clean = []
    for p in poly:
        if not clean or math.hypot(p[0] - clean[-1][0],
                                   p[1] - clean[-1][1]) > 1e-9:
            clean.append((float(p[0]), float(p[1])))
    if len(clean) >= 2 and math.hypot(clean[0][0] - clean[-1][0],
                                      clean[0][1] - clean[-1][1]) < 1e-9:
        clean.pop()

    s = data["scalars"]
    if roof_type == "2in12":
        slope = s["roof_slope"]
        z_off = s["roof_z_off"]
        thick = s["roof_thick"]
        bot = (z_off, 0.0, slope)
        top = (z_off + thick, 0.0, slope)
    else:
        z_tr = s["upper_base"] + s["upper_height"]
        slope = s["roof_slope"]
        z_base_rel = s["roof_z_base"]
        bot = (z_tr, 0.0, 0.0)
        top = (z_tr + z_base_rel, 0.0, slope)
    return {"poly": clean, "bot": bot, "top": top}


# ───────────────────────────  Blender model build  ──────────────────────────

def _run_in_blender():
    import bpy  # noqa: F401  (only importable inside Blender)

    # The launching process passes a manifest JSON path after "--" describing
    # the currently-loaded model (config_name, roof_type, SCAD text, parcel).
    manifest = None
    argv = sys.argv
    if "--" in argv:
        extra = argv[argv.index("--") + 1:]
        if extra and os.path.isfile(extra[0]):
            with open(extra[0], encoding="utf-8") as f:
                manifest = json.load(f)

    if manifest and manifest.get("models"):
        parcel = manifest.get("parcel")
        footprint = manifest.get("footprint")
        for m in manifest["models"]:
            meshes, roof_spec, iw_walls, door_h = _build_meshes(
                m["scad_text"], m["roof_type"])
            _build_blend_file(m["name"], meshes, parcel, roof_spec,
                              iw_walls, door_h, footprint)
        return

    # Legacy fallback: build both roof types from the committed scad files
    # (no DB / manifest — e.g. when Blender is launched directly).
    for stem, scad_file, roof_type in _MODELS:
        meshes, roof_spec, iw_walls, door_h = _build_meshes(
            _read_scad_file(scad_file), roof_type)
        _build_blend_file(stem, meshes, None, roof_spec, iw_walls, door_h)


def _make_material(bpy, key):
    rgba, roughness = _MATERIALS[key]
    mat = bpy.data.materials.new(name=key)
    mat.use_nodes = True
    bsdf = mat.node_tree.nodes.get("Principled BSDF")
    if bsdf is not None:
        bsdf.inputs["Base Color"].default_value = rgba
        bsdf.inputs["Roughness"].default_value = roughness
        bsdf.inputs["Metallic"].default_value = 0.0
        # Alpha input name differs across versions; set defensively.
        if "Alpha" in bsdf.inputs:
            bsdf.inputs["Alpha"].default_value = rgba[3]
    if rgba[3] < 1.0:
        # Legacy EEVEE / older Blender.
        if hasattr(mat, "blend_method"):
            mat.blend_method = "BLEND"
        if hasattr(mat, "show_transparent_back"):
            mat.show_transparent_back = False
        # EEVEE Next (Blender 4.2+/5.x) uses a per-material render method.
        if hasattr(mat, "surface_render_method"):
            mat.surface_render_method = "BLENDED"
    # Viewport display colour (solid shading); alpha drives solid-mode opacity.
    mat.diffuse_color = rgba
    return mat


def _build_grass(bpy, coll, stem, parcel, footprint=None):
    """Add a flat green ground plane over the parcel, with the building footprint
    cut out so no lawn shows inside the exterior outline."""
    if not parcel or len(parcel) < 3:
        return
    import mathutils
    from mathutils.geometry import tessellate_polygon

    outer = [mathutils.Vector((float(x), float(y), 0.0)) for x, y in parcel]
    loops = [outer]
    if footprint and len(footprint) >= 3:
        loops.append([mathutils.Vector((float(x), float(y), 0.0))
                      for x, y in footprint])

    tris = tessellate_polygon(loops)
    verts = [(v.x, v.y, _GRASS_Z) for loop in loops for v in loop]
    name = f"{stem}_Grass"
    mesh = bpy.data.meshes.new(name)
    mesh.from_pydata(verts, [], tris)
    mesh.update()
    mesh.validate()
    mesh.materials.append(_make_material(bpy, "grass"))
    for poly in mesh.polygons:
        poly.use_smooth = False
    coll.objects.link(bpy.data.objects.new(name, mesh))


def _ccw(poly):
    """Return the polygon wound counter-clockwise (positive signed area).

    Prism faces are built assuming a CCW base so their normals point outward;
    a CW input would invert the cutter and the boolean would remove nothing.
    """
    n = len(poly)
    area2 = sum(poly[i][0] * poly[(i + 1) % n][1]
                - poly[(i + 1) % n][0] * poly[i][1] for i in range(n))
    return list(poly) if area2 >= 0 else list(reversed(poly))


def _prism_mesh(bpy, name, poly2d, z_bot, z_top):
    """Create a closed prism mesh from a 2D polygon extruded z_bot→z_top."""
    poly2d = _ccw(poly2d)
    n = len(poly2d)
    verts = [(float(x), float(y), z_bot) for x, y in poly2d]
    verts += [(float(x), float(y), z_top) for x, y in poly2d]
    faces = [tuple(range(n - 1, -1, -1)), tuple(range(n, 2 * n))]
    for i in range(n):
        j = (i + 1) % n
        faces.append((i, j, n + j, n + i))
    mesh = bpy.data.meshes.new(name)
    mesh.from_pydata(verts, [], faces)
    mesh.update()
    mesh.validate()
    return mesh


def _build_interior_walls_doors(bpy, coll, stem, iw_walls, door_h):
    """Build interior walls with their door openings cut out.

    Each wall is extruded floor→z_top; door openings are subtracted with a
    boolean (a cutter prism per opening, from just below the floor up to the
    door height) so the doorways read as real openings with a header above.
    """
    if not iw_walls:
        return
    mat = _make_material(bpy, "iw")
    for i, w in enumerate(iw_walls):
        poly = w["poly"]
        if len(poly) < 3:
            continue
        wmesh = _prism_mesh(bpy, f"{stem}_IW{i}", poly, 0.0, w["z_top"])
        wobj = bpy.data.objects.new(f"{stem}_InteriorWall{i}", wmesh)
        coll.objects.link(wobj)

        if w["openings"]:
            cv, cf = [], []
            for op in w["openings"]:
                op = _ccw(op)
                n = len(op)
                if n < 3:
                    continue
                b = len(cv)
                cv += [(float(x), float(y), -0.1) for x, y in op]
                cv += [(float(x), float(y), door_h) for x, y in op]
                cf.append(tuple(range(b + n - 1, b - 1, -1)))
                cf.append(tuple(range(b + n, b + 2 * n)))
                for k in range(n):
                    j = (k + 1) % n
                    cf.append((b + k, b + j, b + n + j, b + n + k))
            cmesh = bpy.data.meshes.new(f"{stem}_IWcut{i}")
            cmesh.from_pydata(cv, [], cf)
            cmesh.update()
            cmesh.validate()
            cobj = bpy.data.objects.new(f"{stem}_IWcut{i}", cmesh)
            coll.objects.link(cobj)

            mod = wobj.modifiers.new("door", "BOOLEAN")
            mod.operation = "DIFFERENCE"
            mod.object = cobj
            if hasattr(mod, "solver"):
                mod.solver = "EXACT"
            bpy.context.view_layer.objects.active = wobj
            bpy.ops.object.modifier_apply(modifier=mod.name)
            bpy.data.objects.remove(cobj, do_unlink=True)

        wobj.data.materials.append(mat)
        for p in wobj.data.polygons:
            p.use_smooth = False


def _box_into(verts, faces, p0, p1, half_w, z_bot, z_top):
    """Append an oriented box to verts/faces.

    The box's central axis runs from p0=(x,y) to p1=(x,y) in the XY plane,
    extruded ±half_w perpendicular to that axis and from z_bot to z_top.
    """
    x0, y0 = p0
    x1, y1 = p1
    dx, dy = x1 - x0, y1 - y0
    length = math.hypot(dx, dy)
    if length < 1e-9:
        ux, uy = 1.0, 0.0
    else:
        ux, uy = dx / length, dy / length
    px, py = -uy, ux  # perpendicular in XY plane

    corners_xy = [
        (x0 - half_w * px, y0 - half_w * py),  # start, -w
        (x1 - half_w * px, y1 - half_w * py),  # end,   -w
        (x1 + half_w * px, y1 + half_w * py),  # end,   +w
        (x0 + half_w * px, y0 + half_w * py),  # start, +w
    ]
    b = len(verts)
    for cx, cy in corners_xy:
        verts.append((cx, cy, z_bot))
    for cx, cy in corners_xy:
        verts.append((cx, cy, z_top))
    faces.extend([
        (b + 0, b + 1, b + 2, b + 3),  # bottom
        (b + 4, b + 5, b + 6, b + 7),  # top
        (b + 0, b + 1, b + 5, b + 4),  # sides
        (b + 1, b + 2, b + 6, b + 5),
        (b + 2, b + 3, b + 7, b + 6),
        (b + 3, b + 0, b + 4, b + 7),
    ])


def _build_fence(bpy, coll, stem, p_start, p_end, mat_key, n_rails, top_h, label):
    """Build a post-and-rail fence object along the segment p_start→p_end.

    n_rails horizontal rails are evenly spaced with the top rail's top at
    top_h; square posts run from the ground to top_h every _POST_SPACING feet
    (and at both endpoints).
    """
    x0, y0 = p_start[0], p_start[1]
    x1, y1 = p_end[0], p_end[1]
    dx, dy = x1 - x0, y1 - y0
    length = math.hypot(dx, dy)
    if length < 1e-9:
        return
    ux, uy = dx / length, dy / length

    verts, faces = [], []

    # Rails: tops evenly spaced; rail k (1..n) tops at top_h * k / n_rails.
    for k in range(1, n_rails + 1):
        rail_top = top_h * k / n_rails
        _box_into(verts, faces, (x0, y0), (x1, y1),
                  _RAIL_THICK / 2.0, rail_top - _RAIL_H, rail_top)

    # Posts every _POST_SPACING ft along the line, plus one at the far end.
    half = _POST_W / 2.0
    dists = []
    d = 0.0
    while d < length - 1e-6:
        dists.append(d)
        d += _POST_SPACING
    dists.append(length)
    for dd in dists:
        cx, cy = x0 + ux * dd, y0 + uy * dd
        _box_into(verts, faces,
                  (cx - ux * half, cy - uy * half),
                  (cx + ux * half, cy + uy * half),
                  half, 0.0, top_h)

    name = f"{stem}_{label}"
    mesh = bpy.data.meshes.new(name)
    mesh.from_pydata(verts, [], faces)
    mesh.update()
    mesh.validate()
    mesh.materials.append(_make_material(bpy, mat_key))
    for poly in mesh.polygons:
        poly.use_smooth = False
    coll.objects.link(bpy.data.objects.new(name, mesh))


def _build_roof_native(bpy, coll, stem, roof_spec):
    """Build the roof as a solid planar-capped prism (no triangulation gaps).

    Top and bottom are single n-gon faces (Blender triangulates planar n-gons
    cleanly); side quads connect them, so the roof reads solid from above and
    below across its whole footprint.
    """
    if not roof_spec:
        return
    poly = roof_spec["poly"]
    n = len(poly)
    if n < 3:
        return
    ba, bbx, bcy = roof_spec["bot"]
    ta, tbx, tcy = roof_spec["top"]
    verts = [(x, y, ba + bbx * x + bcy * y) for x, y in poly]
    verts += [(x, y, ta + tbx * x + tcy * y) for x, y in poly]
    faces = [tuple(range(n - 1, -1, -1)),   # bottom n-gon (downward normal)
             tuple(range(n, 2 * n))]        # top n-gon (upward normal)
    for i in range(n):
        j = (i + 1) % n
        faces.append((i, j, n + j, n + i))  # side quad
    name = f"{stem}_Roof"
    mesh = bpy.data.meshes.new(name)
    mesh.from_pydata(verts, [], faces)
    mesh.update()
    mesh.validate()
    mesh.materials.append(_make_material(bpy, "roof"))
    for p in mesh.polygons:
        p.use_smooth = False
    coll.objects.link(bpy.data.objects.new(name, mesh))


def _build_blend_file(stem, meshes, parcel=None, roof_spec=None,
                      iw_walls=None, door_h=None, footprint=None):
    import bpy

    # Start from a clean, empty scene.
    bpy.ops.wm.read_factory_settings(use_empty=True)

    # Don't leave .blend1 backup files behind when overwriting.
    bpy.context.preferences.filepaths.save_version = 0

    scene = bpy.context.scene
    scene.unit_settings.system = "IMPERIAL"
    scene.unit_settings.length_unit = "FEET"
    scene.unit_settings.scale_length = 0.3048  # 1 Blender unit == 1 foot

    coll = scene.collection

    _key_label = {"wall": "Walls", "win": "Windows"}

    for key in ("wall", "win"):
        if key not in meshes:
            continue
        verts, indices = meshes[key]
        faces = [tuple(int(i) for i in indices[t:t + 3])
                 for t in range(0, len(indices), 3)]
        vlist = [(float(v[0]), float(v[1]), float(v[2])) for v in verts]

        name = f"{stem}_{_key_label[key]}"
        mesh = bpy.data.meshes.new(name)
        mesh.from_pydata(vlist, [], faces)
        mesh.update()
        mesh.validate()

        mat = _make_material(bpy, key)
        mesh.materials.append(mat)

        obj = bpy.data.objects.new(name, mesh)
        coll.objects.link(obj)
        # Flat shading to match the faceted SCAD/GLB look.
        for poly in mesh.polygons:
            poly.use_smooth = False

    _build_roof_native(bpy, coll, stem, roof_spec)

    _build_interior_walls_doors(bpy, coll, stem, iw_walls, door_h)

    _build_grass(bpy, coll, stem, parcel, footprint)

    # Boundary fences (parcel order is [NW, NE, SE, SW]): the dark 3-rail fence
    # runs along the WEST boundary (NW–SW, the 275.08' line); the lighter
    # single-rail fence runs along the SOUTH boundary (SW–SE, the 216.73'
    # line).  They meet at the SW corner.
    if parcel and len(parcel) >= 4:
        nw, ne, se, sw = parcel[0], parcel[1], parcel[2], parcel[3]
        _build_fence(bpy, coll, stem, sw, nw, "fence_3rail",
                     3, _FENCE3_TOP, "Fence3Rail")
        _build_fence(bpy, coll, stem, se, sw, "fence_1rail",
                     1, _FENCE1_TOP, "Fence1Rail")

    out_path = os.path.join(_DIR, f"{stem}.blend")
    bpy.ops.wm.save_as_mainfile(filepath=out_path)
    print(f"  wrote {out_path}  ({len(coll.objects)} object(s))")


# ─────────────────────  currently-loaded model (launcher side)  ──────────────

def _parcel_from_gd(gd):
    """Survey parcel outline as [[E, N], ...] in FC-origin feet, from gd.

    Inverts site/gen_site_plan.py's building→PDF placement transform for the
    four parcel corners.  Returns None on failure so the model still builds
    without a grass plane.
    """
    try:
        import importlib.util
        # Load site/gen_site_plan.py by path (the dir name "site" clashes with
        # the stdlib site module, so a normal import is unsafe).
        sp_path = os.path.join(_ROOT, "site", "gen_site_plan.py")
        spec = importlib.util.spec_from_file_location("adu_gen_site_plan", sp_path)
        gsp = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(gsp)

        sp = gsp.build_site_plan_data(gd)
        f15 = sp.pts["_site_se_pt"]
        fpx, fpy = sp.f15_pdf
        scale = sp.SCALE
        rot = math.radians(sp.rotation_deg)
        cos_r, sin_r = math.cos(rot), math.sin(rot)

        def pdf_to_building(px, py):
            re = f15[0] + (px - fpx) / scale
            rn = f15[1] - (py - fpy) / scale
            de, dn = re - f15[0], rn - f15[1]
            return [f15[0] + de * cos_r + dn * sin_r,
                    f15[1] - de * sin_r + dn * cos_r]

        # Corners in true compass order: NW, NE, SE, SW.  (The survey plat is
        # drawn N-left / E-top, not N-up — see site/gen_site_plan.py.  The
        # building/FC frame uses true Easting/Northing, so these map directly
        # to the model's X=E / Y=N ground plane.)
        return [pdf_to_building(*c) for c in
                (gsp.CORNER_NW, gsp.CORNER_NE, gsp.CORNER_SE, gsp.CORNER_SW)]
    except Exception as exc:
        print(f"gen_blend: parcel computation failed ({exc}); no grass plane")
        return None


def _scad_text_from_gd(ggl, gd, roof_type):
    """Generate SCAD text for roof_type from gd, restoring the on-disk file.

    gen_gltf._generate_variant_scad writes scad/<roof_type>.scad as a side
    effect; we save and restore the committed file so a standalone run doesn't
    leave the working tree dirty.
    """
    scad_path = os.path.join(_ROOT, "scad", f"{roof_type}.scad")
    original = None
    if os.path.exists(scad_path):
        with open(scad_path, encoding="utf-8", errors="replace") as f:
            original = f.read()
    try:
        return ggl._generate_variant_scad(gd, roof_type)
    finally:
        if original is not None:
            with open(scad_path, "w", encoding="utf-8") as f:
                f.write(original)


def _prepare_from_db():
    """Build a manifest for the currently-loaded model from app/adu.db.

    Produces two .blend files from the same DB geometry: the currently-loaded
    model with its own roof (named <config_name>) and a 2:12 sloped-roof
    version (named <config_name>_2in12).

    Returns {models: [{name, roof_type, scad_text}, ...], parcel} or None on
    failure (caller then falls back to the committed-scad legacy build).
    """
    try:
        import sqlite3
        if _ROOT not in sys.path:
            sys.path.insert(0, _ROOT)
        from app.engine import build_generator_data_from_db
        from scad import gen_gltf as ggl

        db_path = os.path.join(_ROOT, "app", "adu.db")
        gd = build_generator_data_from_db(db_path)

        # Model identity + roof style from the live config table.
        con = sqlite3.connect(db_path)
        cfg = dict(con.execute("SELECT key, value FROM config").fetchall())
        con.close()
        config_name = cfg.get("config_name") or "model"
        roof_style = (cfg.get("roof_style") or "flat").lower()
        roof_type = "flat_roof" if roof_style.startswith("flat") else "2in12"

        # Currently-loaded model (its configured roof).
        models = [{"name": config_name, "roof_type": roof_type,
                   "scad_text": _scad_text_from_gd(ggl, gd, roof_type)}]
        # 2:12 sloped-roof version as its own file (skip if already 2:12).
        if roof_type != "2in12":
            models.append({"name": f"{config_name}_2in12", "roof_type": "2in12",
                           "scad_text": _scad_text_from_gd(ggl, gd, "2in12")})

        parcel = _parcel_from_gd(gd)
        # Building exterior outline (FC feet) — used to cut a hole in the grass
        # so the lawn doesn't show under the house.
        footprint = [[float(x), float(y)] for x, y in gd.outline_poly]
        return {"models": models, "parcel": parcel, "footprint": footprint}
    except Exception as exc:
        print(f"gen_blend: could not build from DB ({exc}); "
              "falling back to committed scad files")
        return None


# ─────────────────────────  headless re-launch shim  ─────────────────────────

def _find_blender():
    """Locate a blender executable; honour $BLENDER_EXE override."""
    env = os.environ.get("BLENDER_EXE")
    if env and os.path.isfile(env):
        return env
    import glob
    import shutil
    found = shutil.which("blender")
    if found:
        return found
    pats = [
        r"C:\Program Files\Blender Foundation\*\blender.exe",
        r"C:\Program Files (x86)\Blender Foundation\*\blender.exe",
    ]
    cands = []
    for pat in pats:
        cands.extend(glob.glob(pat))
    if cands:
        # Prefer the highest version directory.
        return sorted(cands)[-1]
    return None


def generate(gd=None):
    """Entry point.  Re-launches under Blender if bpy is unavailable.

    ``gd`` is accepted for gen_all-style compatibility but unused: the model is
    built from the live editor DB (app/adu.db) for the currently-loaded config.
    """
    try:
        import bpy  # noqa: F401
        _run_in_blender()
        return
    except ImportError:
        pass

    import subprocess
    import tempfile
    blender = _find_blender()
    if not blender:
        print("Blender not found (set $BLENDER_EXE); skipping .blend generation")
        return

    cmd = [blender, "--background", "--python", os.path.abspath(__file__)]

    # Build the currently-loaded model (SCAD + parcel) here, where the full app
    # stack and PyMuPDF are available, and hand it to Blender as a manifest.
    manifest = _prepare_from_db()
    manifest_path = None
    if manifest:
        fd, manifest_path = tempfile.mkstemp(prefix="adu_blend_", suffix=".json")
        with os.fdopen(fd, "w", encoding="utf-8") as f:
            json.dump(manifest, f)
        cmd += ["--", manifest_path]
        _names = ", ".join(f"{m['name']} ({m['roof_type']})"
                           for m in manifest["models"])
        print(f"gen_blend: building from app/adu.db: {_names}")
    else:
        print("gen_blend: building legacy flat_roof + 2in12 from committed scad")

    print(f"gen_blend: launching {blender} headless ...")
    try:
        subprocess.check_call(cmd, cwd=_ROOT)
    finally:
        if manifest_path and os.path.exists(manifest_path):
            os.remove(manifest_path)


if __name__ == "__main__":
    generate()
