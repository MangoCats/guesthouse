"""Generate Blender (.blend) models of the ADU, mirroring the .scad models.

This builds the same geometry the OpenSCAD / glTF pipeline produces — outer
wall shells, the sloped/wedge roof slab (with overhang), glazing panels and
interior walls — but as native Blender mesh objects with materials, saved to
.blend files.

Accuracy: the triangle geometry is reused verbatim from ``scad/gen_gltf.py``
(the same parser + mesh builders that drive ``opts/*.glb``), which in turn is
fed by the DB-generated ``scad/flat_roof.scad`` / ``scad/2in12.scad``.  No
geometry is re-derived here, so the Blender model matches the SCAD/GLB models
to the millimetre.

Coordinate system: SCAD/Blender are both Z-up (east X, north Y, up Z), so —
unlike the glTF export — no axis swap is applied.  Units are feet; the scene
is configured for imperial display (1 Blender unit = 1 foot).

A flat green ground plane covering the survey parcel (the property extent from
the site plans) is added beneath the building.  Its outline is the four parcel
corners from site/gen_site_plan.py, transformed from PDF/survey space back into
the building's FC-origin feet frame.  This is computed in the launching process
(which has PyMuPDF) and passed to headless Blender as a JSON argument, since
Blender's bundled Python lacks fitz.

Outputs: blend/flat_roof.blend, blend/2in12.blend

Usage:
    python blend/gen_blend.py            # locates Blender and runs headless
    blender --background --python blend/gen_blend.py
"""

import os
import sys
import json
import math

_DIR = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_DIR)

# (display stem, source .scad file, roof_type for gen_gltf builders)
_MODELS = [
    ("flat_roof", "flat_roof.scad", "flat_roof"),
    ("2in12",     "2in12.scad",     "2in12"),
]

# Material palette — matches scad/gen_gltf.py (RGBA, roughness)
_MATERIALS = {
    "wall": ([0.88, 0.82, 0.60, 1.00], 0.85),  # warm cream-yellow shells
    "roof": ([0.10, 0.35, 0.33, 1.00], 0.55),  # dark teal-green metal slab
    "win":  ([0.80, 0.84, 0.90, 0.90], 0.10),  # blue-grey glazing
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

def _build_meshes(scad_file, roof_type):
    """Parse a .scad file and return {key: (verts, indices)} per material.

    Reuses scad/gen_gltf.py's parser and mesh builders verbatim so the Blender
    geometry is identical to the GLB geometry.  ``verts`` is an (N,3) float
    array (SCAD Z-up coords, feet); ``indices`` is a flat triangle-index array.
    """
    if _ROOT not in sys.path:
        sys.path.insert(0, _ROOT)
    from scad import gen_gltf as g

    data = g._load_scad(scad_file)

    wall_mesh = g.MeshBuilder()
    roof_mesh = g.MeshBuilder()
    win_mesh  = g.MeshBuilder()
    iw_mesh   = g.MeshBuilder()

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
        g._build_roof_2in12(data, roof_mesh)
    else:
        g._build_walls(data, wall_mesh)
        g._build_roof_flat(data, roof_mesh)

    g._build_windows(data, win_mesh)
    g._build_interior_walls(data, iw_mesh)

    out = {}
    for key, mb in (("wall", wall_mesh), ("roof", roof_mesh),
                    ("win", win_mesh), ("iw", iw_mesh)):
        verts, _norms, indices = mb.to_arrays()
        if verts is not None and len(verts):
            out[key] = (verts, indices)
    return out


# ───────────────────────────  Blender model build  ──────────────────────────

def _run_in_blender():
    import bpy  # noqa: F401  (only importable inside Blender)

    # Parcel polygon (building/FC feet coords) is passed after "--" as a JSON
    # file path by the launching process; absent when run directly in Blender.
    parcel = None
    argv = sys.argv
    if "--" in argv:
        extra = argv[argv.index("--") + 1:]
        if extra and os.path.isfile(extra[0]):
            with open(extra[0]) as f:
                parcel = json.load(f)

    for stem, scad_file, roof_type in _MODELS:
        meshes = _build_meshes(scad_file, roof_type)
        _build_blend_file(stem, meshes, parcel)


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
        mat.blend_method = "BLEND"
    # Viewport display colour (solid shading).
    mat.diffuse_color = rgba
    return mat


def _build_grass(bpy, coll, stem, parcel):
    """Add a flat green ground-plane ngon covering the parcel outline."""
    if not parcel or len(parcel) < 3:
        return
    verts = [(float(x), float(y), _GRASS_Z) for x, y in parcel]
    face = [tuple(range(len(verts)))]
    name = f"{stem}_Grass"
    mesh = bpy.data.meshes.new(name)
    mesh.from_pydata(verts, [], face)
    mesh.update()
    mesh.validate()
    mesh.materials.append(_make_material(bpy, "grass"))
    for poly in mesh.polygons:
        poly.use_smooth = False
    coll.objects.link(bpy.data.objects.new(name, mesh))


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


def _build_blend_file(stem, meshes, parcel=None):
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

    _key_label = {"wall": "Walls", "roof": "Roof",
                  "win": "Windows", "iw": "InteriorWalls"}

    for key in ("wall", "roof", "win", "iw"):
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

    _build_grass(bpy, coll, stem, parcel)

    # Boundary fences: 3-rail along the SE→SW edge; single-rail along the
    # NE→SE edge (the empty side opposite the 3-rail). They meet at the SE
    # corner.  parcel order is [NW, SW, SE, NE].
    if parcel and len(parcel) >= 4:
        nw, sw, se, ne = parcel[0], parcel[1], parcel[2], parcel[3]
        _build_fence(bpy, coll, stem, se, sw, "fence_3rail",
                     3, _FENCE3_TOP, "Fence3Rail")
        _build_fence(bpy, coll, stem, ne, se, "fence_1rail",
                     1, _FENCE1_TOP, "Fence1Rail")

    out_path = os.path.join(_DIR, f"{stem}.blend")
    bpy.ops.wm.save_as_mainfile(filepath=out_path)
    print(f"  wrote {out_path}  ({len(coll.objects)} object(s))")


# ───────────────────────────  parcel (launcher side)  ───────────────────────

def _compute_parcel_poly():
    """Return the survey parcel outline as [[E, N], ...] in FC-origin feet.

    Inverts site/gen_site_plan.py's building→PDF placement transform for the
    four parcel corners.  Runs in the launching interpreter (needs PyMuPDF);
    returns None on any failure so the model still builds without grass.
    """
    try:
        import math
        import importlib.util
        if _ROOT not in sys.path:
            sys.path.insert(0, _ROOT)
        from app.engine import build_generator_data_from_db

        # Load site/gen_site_plan.py by path (the dir name "site" clashes with
        # the stdlib site module, so a normal import is unsafe).
        sp_path = os.path.join(_ROOT, "site", "gen_site_plan.py")
        spec = importlib.util.spec_from_file_location("adu_gen_site_plan", sp_path)
        gsp = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(gsp)

        gd = build_generator_data_from_db(os.path.join(_ROOT, "app", "adu.db"))
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

        # Order corners CCW around the parcel (NW→SW→SE→NE).
        return [pdf_to_building(*c) for c in
                (gsp.CORNER_NW, gsp.CORNER_SW, gsp.CORNER_SE, gsp.CORNER_NE)]
    except Exception as exc:
        print(f"gen_blend: parcel computation failed ({exc}); no grass plane")
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

    ``gd`` is accepted for gen_all-style compatibility but unused: geometry is
    read from the already-generated scad/*.scad files (which gen_all writes
    from the DB before this runs).
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

    # Compute the parcel here (PyMuPDF available) and hand it to Blender.
    parcel = _compute_parcel_poly()
    parcel_path = None
    if parcel:
        fd, parcel_path = tempfile.mkstemp(prefix="adu_parcel_", suffix=".json")
        with os.fdopen(fd, "w") as f:
            json.dump(parcel, f)
        cmd += ["--", parcel_path]

    print(f"gen_blend: launching {blender} headless ...")
    try:
        subprocess.check_call(cmd, cwd=_ROOT)
    finally:
        if parcel_path and os.path.exists(parcel_path):
            os.remove(parcel_path)


if __name__ == "__main__":
    generate()
