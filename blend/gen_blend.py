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

Outputs: blend/flat_roof.blend, blend/2in12.blend

Usage:
    python blend/gen_blend.py            # locates Blender and runs headless
    blender --background --python blend/gen_blend.py
"""

import os
import sys

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
}


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

    for stem, scad_file, roof_type in _MODELS:
        meshes = _build_meshes(scad_file, roof_type)
        _build_blend_file(stem, meshes)


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


def _build_blend_file(stem, meshes):
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

    out_path = os.path.join(_DIR, f"{stem}.blend")
    bpy.ops.wm.save_as_mainfile(filepath=out_path)
    tris = sum(len(meshes[k][1]) // 3 for k in meshes)
    print(f"  wrote {out_path}  ({tris:,} tris, "
          f"{len(meshes)} object(s))")


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
    blender = _find_blender()
    if not blender:
        print("Blender not found (set $BLENDER_EXE); skipping .blend generation")
        return
    print(f"gen_blend: launching {blender} headless ...")
    subprocess.check_call(
        [blender, "--background", "--python", os.path.abspath(__file__)],
        cwd=_ROOT,
    )


if __name__ == "__main__":
    generate()
