"""
Generate glTF binary (.glb) models from flat_roof.scad and 2in12.scad.

Each model is exported with four materials matching the OpenSCAD preview colours:
  wall_cream  [0.88, 0.82, 0.60] — outer/inner wall shells
  roof_teal   [0.10, 0.35, 0.33] — roof slab
  window_bg   [0.80, 0.84, 0.90] — glazing panels
  wood_tan    [0.82, 0.72, 0.38] — interior walls

Outputs: opts/flat_roof.glb, opts/2in12.glb, opts/models.js (base64 bundle)
"""

import ast, json, math, re, struct, os, sys
import numpy as np

_DIR  = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_DIR)
_OPTS = os.path.join(_ROOT, "opts")

WALL_CREAM = [0.88, 0.82, 0.60, 1.0]
ROOF_TEAL  = [0.10, 0.35, 0.33, 1.0]
WINDOW_BG  = [0.80, 0.84, 0.90, 0.90]
WOOD_TAN   = [0.82, 0.72, 0.38, 1.0]


# ─────────────────────────────  SCAD parsing  ────────────────────────────────

def _parse_scalars(text):
    """Return dict of simple numeric assignments: name = float;"""
    result = {}
    for m in re.finditer(r'^(\w+)\s*=\s*([-\d.]+)\s*;', text, re.MULTILINE):
        result[m.group(1)] = float(m.group(2))
    return result


def _extract_array(text, name):
    """Extract named array variable from SCAD text using bracket matching."""
    pat = re.compile(rf'\b{re.escape(name)}\s*=\s*(\[)')
    m = pat.search(text)
    if not m:
        return None
    start = m.start(1)
    depth = 0
    for i, c in enumerate(text[start:]):
        if c == '[':
            depth += 1
        elif c == ']':
            depth -= 1
            if depth == 0:
                end = start + i + 1
                break
    chunk = text[start:end]
    chunk = re.sub(r'//[^\n]*', '', chunk)   # strip line comments
    return ast.literal_eval(chunk)


def _bracket_array(text, start):
    """Return the string of the first [...] starting at 'start' using depth matching."""
    depth = 0
    for i, c in enumerate(text[start:]):
        if c == '[':
            depth += 1
        elif c == ']':
            depth -= 1
            if depth == 0:
                return text[start: start + i + 1]
    return None


def _parse_poly_extrusions(block, skip_door_cutouts=False):
    """
    Return list of (pts, z_base, height) for every
      [translate([0,0,Z])] linear_extrude(height=H) polygon(points=[...])
    in block.  If skip_door_cutouts, skip entries where the preceding context
    contains 'translate([0, 0, -0.001])'.
    """
    result = []
    # Find every polygon(points = [...]) call
    for pm in re.finditer(r'polygon\s*\(\s*points\s*=\s*(\[)', block):
        arr_str = _bracket_array(block, pm.start(1))
        if arr_str is None:
            continue
        pts = ast.literal_eval(re.sub(r'//[^\n]*', '', arr_str))

        # Look back up to 300 chars for height and z
        ctx = block[max(0, pm.start() - 300): pm.start()]

        if skip_door_cutouts and '-0.001' in ctx:
            continue

        leh = re.findall(r'linear_extrude\s*\(\s*height\s*=\s*([\d.]+)', ctx)
        if not leh:
            continue
        h = float(leh[-1])

        tz = re.findall(r'translate\s*\(\s*\[0\s*,\s*0\s*,\s*([\d.]+)\]', ctx)
        z  = float(tz[-1]) if tz else 0.0

        result.append((pts, z, h))
    return result


def _parse_assembly(text):
    """
    Parse the assembly block after '// --- Assembly ---'.

    Returns a dict with keys:
      'walls'   : list of (t_path_name, z_base, height)  — all in wall_cream
      'windows' : list of (poly_pts, z_base, height)     — in window_blue_grey
      'iw'      : list of (poly_pts, z_base, height)     — interior walls (wood_tan)
      'iw_door_h': float
      'roof_outline': list of SCAD segment arrays (for shell_pts)
      'scalars' : dict
    """
    s = _parse_scalars(text)

    asm_pos = text.find('// --- Assembly ---')
    asm = text[asm_pos:]

    # ── wall shells ──────────────────────────────────────────────────────────
    walls = []
    for m in re.finditer(r'wall_shell\((\w+),\s*half_t\)', asm):
        pos = m.start()
        ctx = asm[max(0, pos - 400):pos]
        tz  = re.findall(r'translate\(\[0\s*,\s*0\s*,\s*([\d.]+)\]', ctx)
        leh = re.findall(r'linear_extrude\(height\s*=\s*([\d.]+)', ctx)
        z   = float(tz[-1])  if tz  else 0.0
        h   = float(leh[-1]) if leh else 0.0
        walls.append((m.group(1), z, h))

    # ── window polygons (inside window_blue_grey block) ─────────────────────
    win_start = asm.find('window_blue_grey')
    win_end   = asm.find('// Interior walls')
    win_block = asm[win_start:win_end] if win_start >= 0 and win_end >= 0 else ''

    windows = []
    windows.extend(_parse_poly_extrusions(win_block))

    # ── interior walls (inside color([R,G,B]) union block) ───────────────────
    iw_start = asm.find('// Interior walls')
    iw_block = asm[iw_start:] if iw_start >= 0 else ''

    # Grab all polygon extrusions, skipping door-cutout ones (preceded by -0.001)
    iw = []
    for pts, z_base, height in _parse_poly_extrusions(iw_block, skip_door_cutouts=True):
        iw.append((pts, 0.0, height))

    roof_outline = _extract_array(text, 'roof_outline')

    return {
        'walls':        walls,
        'windows':      windows,
        'iw':           iw,
        'iw_door_h':    s.get('iw_door_h', 6.875),
        'roof_outline': roof_outline,
        'scalars':      s,
    }


# ──────────────────────────────  Shell geometry  ─────────────────────────────

_ARC_STEP_DEG = 10   # degrees per arc segment — coarser than SCAD (3°) but imperceptible on screen


def _arc_n(sweep_deg):
    return max(3, int(round(abs(sweep_deg) / _ARC_STEP_DEG)))


def _line_pts(seg, d):
    x1, y1, x2, y2 = seg[1], seg[2], seg[3], seg[4]
    dx, dy = x2 - x1, y2 - y1
    length = math.hypot(dx, dy)
    if length < 1e-12:
        return [[x1, y1], [x2, y2]]
    nx, ny = -dy / length, dx / length
    return [[x1 + d * nx, y1 + d * ny], [x2 + d * nx, y2 + d * ny]]


def _arc_pts(seg, d):
    cx, cy, r, a1, a2 = seg[1], seg[2], seg[3], seg[4], seg[5]
    r_off = r - d if a2 > a1 else r + d
    n = _arc_n(a2 - a1)
    return [
        [cx + r_off * math.cos(math.radians(a1 + i * (a2 - a1) / n)),
         cy + r_off * math.sin(math.radians(a1 + i * (a2 - a1) / n))]
        for i in range(n + 1)
    ]


def shell_pts(path, d):
    """Python equivalent of OpenSCAD shell_pts(path, d)."""
    pts = []
    for i, seg in enumerate(path):
        sp = _line_pts(seg, d) if seg[0] == 0 else _arc_pts(seg, d)
        pts.extend(sp if i == 0 else sp[1:])
    return pts


# ───────────────────────────────  Mesh builder  ──────────────────────────────

class MeshBuilder:
    """Accumulates flat-shaded triangles (each triangle owns its own 3 vertices)."""

    def __init__(self):
        self._verts = []
        self._norms = []
        self._idx   = []
        self._n     = 0

    def add_tri(self, p0, p1, p2):
        v0 = np.asarray(p0, dtype=np.float32)
        v1 = np.asarray(p1, dtype=np.float32)
        v2 = np.asarray(p2, dtype=np.float32)
        n  = np.cross(v1 - v0, v2 - v0)
        ln = float(np.linalg.norm(n))
        if ln < 1e-14:
            return
        n = (n / ln).astype(np.float32)
        base = self._n
        self._verts.extend([v0, v1, v2])
        self._norms.extend([n, n, n])
        self._idx.extend([base, base + 1, base + 2])
        self._n += 3

    def add_quad(self, p0, p1, p2, p3):
        self.add_tri(p0, p1, p2)
        self.add_tri(p0, p2, p3)

    def empty(self):
        return self._n == 0

    def to_arrays(self):
        if self.empty():
            return None, None, None
        verts   = np.array(self._verts, dtype=np.float32)
        norms   = np.array(self._norms, dtype=np.float32)
        indices = np.array(self._idx, dtype=np.uint32)
        return verts, norms, indices


# ─────────────────────────────  3-D primitives  ──────────────────────────────

def _stitch_caps(outer_3d, inner_3d, mesh):
    """
    Triangulate the annular (ring) top/bottom face between two 3-D contours.

    outer_3d and inner_3d must have the same length (guaranteed for shell_pts
    since arc_n() depends only on sweep angle, not radius).
    """
    n = len(outer_3d)
    for i in range(n):
        o0 = outer_3d[i]
        o1 = outer_3d[(i + 1) % n]
        i0 = inner_3d[i]
        i1 = inner_3d[(i + 1) % n]
        mesh.add_tri(o0, o1, i0)
        mesh.add_tri(o1, i1, i0)


def extrude_ring(outer_2d, inner_2d, mesh,
                 z_bot=0.0, z_top=None,
                 outer_top_fn=None, inner_top_fn=None):
    """
    Extrude an annular polygon into a solid band.

    z_top (uniform) or per-vertex top height functions may be supplied.
    """
    if outer_top_fn is None:
        ztop = float(z_top)
        outer_top_fn = lambda x, y: ztop
    if inner_top_fn is None:
        inner_top_fn = outer_top_fn

    ob = [[x, y, z_bot]                    for x, y in outer_2d]
    ib = [[x, y, z_bot]                    for x, y in inner_2d]
    ot = [[x, y, outer_top_fn(x, y)]       for x, y in outer_2d]
    it = [[x, y, inner_top_fn(x, y)]       for x, y in inner_2d]

    n = len(outer_2d)
    m = len(inner_2d)

    # Outer side walls
    for i in range(n):
        mesh.add_quad(ob[i], ob[(i + 1) % n], ot[(i + 1) % n], ot[i])

    # Inner side walls (reversed winding for inward-facing normals)
    for i in range(m):
        mesh.add_quad(ib[i], it[i], it[(i + 1) % m], ib[(i + 1) % m])

    # Bottom annular cap
    _stitch_caps(ob, ib, mesh)

    # Top annular cap (winding reversed relative to bottom)
    _stitch_caps(it, ot, mesh)


def extrude_polygon(pts_2d, mesh, z_bot=0.0, z_top=None,
                    top_fn=None, flip_top=False):
    """
    Extrude a simple (convex) polygon into a prism.
    Fan-triangulation from vertex 0 is used for the caps.
    top_fn(x, y) overrides z_top for per-vertex heights.
    """
    if top_fn is None:
        ztop = float(z_top)
        top_fn = lambda x, y: ztop

    bot = [[x, y, z_bot]        for x, y in pts_2d]
    top = [[x, y, top_fn(x, y)] for x, y in pts_2d]
    n   = len(pts_2d)

    # Side walls
    for i in range(n):
        mesh.add_quad(bot[i], bot[(i + 1) % n], top[(i + 1) % n], top[i])

    # Bottom cap (fan from vertex 0, winding: CCW when viewed from below)
    for i in range(1, n - 1):
        mesh.add_tri(bot[0], bot[i + 1], bot[i])

    # Top cap
    tris_top = [(top[0], top[i], top[i + 1]) for i in range(1, n - 1)]
    if flip_top:
        for a, b, c in tris_top:
            mesh.add_tri(a, c, b)
    else:
        for a, b, c in tris_top:
            mesh.add_tri(a, b, c)


# ──────────────────────────────  Model builders  ─────────────────────────────

def _build_walls(data, wall_mesh):
    """Add all outer wall shell extrusions to wall_mesh."""
    half_t = data['scalars']['half_t']
    arrays = data['_arrays']
    for name, z_base, height in data['walls']:
        path = arrays.get(name)
        if path is None:
            continue
        outer = shell_pts(path, -half_t)
        inner = shell_pts(path,  half_t)
        extrude_ring(outer, inner, wall_mesh,
                     z_bot=z_base, z_top=z_base + height)


def _build_walls_2in12_upper(data, wall_mesh):
    """
    Band 2 of the 2in12 model: upper wall clipped by the roof shear plane.
    z_top per point = min(6.666667 + 5.555556, roof_slope * y + roof_z_off)
    """
    s      = data['scalars']
    half_t = s['half_t']
    slope  = s['roof_slope']
    z_off  = s['roof_z_off']
    z_base = 6.666667
    z_max  = z_base + 5.555556

    def top_fn(x, y):
        return min(z_max, slope * y + z_off)

    arrays = data['_arrays']
    path   = arrays.get('t_full_upper')
    if path is None:
        return
    outer = shell_pts(path, -half_t)
    inner = shell_pts(path,  half_t)
    extrude_ring(outer, inner, wall_mesh,
                 z_bot=z_base,
                 outer_top_fn=top_fn, inner_top_fn=top_fn)


def _build_roof_flat(data, roof_mesh):
    """Flat (wedge) roof slab for flat_roof model."""
    s          = data['scalars']
    slope      = s['roof_slope']
    z_base_rel = s['roof_z_base']      # thickness at y = 0 minus slope correction
    z_translate = s['upper_base'] + s['upper_height']

    pts = shell_pts(data['roof_outline'], 0)

    def top_fn(x, y):
        return z_translate + z_base_rel + slope * y

    extrude_polygon(pts, roof_mesh,
                    z_bot=z_translate,
                    top_fn=top_fn)


def _build_roof_2in12(data, roof_mesh):
    """2:12 sloped roof slab."""
    s     = data['scalars']
    slope = s['roof_slope']
    z_off = s['roof_z_off']
    thick = s['roof_thick']

    pts = shell_pts(data['roof_outline'], 0)

    # Bottom face of slab at z = slope*y + z_off  (shear plane)
    # Top  face            at z = slope*y + z_off + thick
    def bot_fn(x, y):
        return slope * y + z_off

    def top_fn(x, y):
        return slope * y + z_off + thick

    bot3 = [[x, y, bot_fn(x, y)] for x, y in pts]
    top3 = [[x, y, top_fn(x, y)] for x, y in pts]
    n    = len(pts)

    # Side walls
    for i in range(n):
        roof_mesh.add_quad(bot3[i], bot3[(i + 1) % n],
                            top3[(i + 1) % n], top3[i])
    # Bottom cap
    for i in range(1, n - 1):
        roof_mesh.add_tri(bot3[0], bot3[i + 1], bot3[i])
    # Top cap
    for i in range(1, n - 1):
        roof_mesh.add_tri(top3[0], top3[i], top3[i + 1])


def _build_windows(data, win_mesh):
    for pts, z_base, height in data['windows']:
        extrude_polygon(pts, win_mesh, z_bot=z_base, z_top=z_base + height)


def _build_interior_walls(data, iw_mesh):
    for pts, z_base, height in data['iw']:
        extrude_polygon(pts, iw_mesh, z_bot=z_base, z_top=z_base + height)


# ──────────────────────────────  glTF / GLB  ─────────────────────────────────

_ROUGHNESS = {
    'wall': 0.85,
    'roof': 0.60,
    'win':  0.10,
    'iw':   0.90,
}


def _build_glb(primitives):
    """
    primitives: list of (verts_Nx3, norms_Nx3, indices_M, rgba_list, roughness)
    Returns bytes of a binary glTF (.glb) file.

    Coordinate system: SCAD X/Y/Z (east/north/up, Z-up).
    model-viewer is configured with camera-orbit to look from above-south-west.
    """
    bin_buf       = bytearray()
    buffer_views  = []
    accessors     = []
    materials     = []
    mesh_prims    = []

    ARRAY_BUF = 34962
    ELEM_BUF  = 34963

    def _add_bv(data_bytes, target):
        off = len(bin_buf)
        bin_buf.extend(data_bytes)
        while len(bin_buf) % 4:
            bin_buf.append(0)
        buffer_views.append({
            "buffer": 0,
            "byteOffset": off,
            "byteLength": len(data_bytes),
            "target": target,
        })
        return len(buffer_views) - 1

    for verts, norms, indices, rgba, roughness in primitives:
        if verts is None or len(verts) == 0:
            continue

        mat_idx = len(materials)
        alpha   = rgba[3]
        mat = {
            "pbrMetallicRoughness": {
                "baseColorFactor": list(rgba),
                "metallicFactor":  0.0,
                "roughnessFactor": roughness,
            },
            "doubleSided": True,
        }
        if alpha < 1.0:
            mat["alphaMode"] = "BLEND"
        materials.append(mat)

        # Convert SCAD Z-up (east, north, height) → glTF Y-up (east, height, -north)
        # so model-viewer shows the building upright by default.
        v = verts.copy().astype(np.float32)
        n = norms.copy().astype(np.float32)
        v[:, 1], v[:, 2] = verts[:, 2], -verts[:, 1]
        n[:, 1], n[:, 2] = norms[:, 2], -norms[:, 1]

        # POSITION
        vb   = v.tobytes()
        bv_v = _add_bv(vb, ARRAY_BUF)
        acc_v = len(accessors)
        accessors.append({
            "bufferView":    bv_v,
            "componentType": 5126,   # FLOAT
            "count":         len(v),
            "type":          "VEC3",
            "min":           v.min(axis=0).tolist(),
            "max":           v.max(axis=0).tolist(),
        })

        # NORMAL
        nb   = n.tobytes()
        bv_n = _add_bv(nb, ARRAY_BUF)
        acc_n = len(accessors)
        accessors.append({
            "bufferView":    bv_n,
            "componentType": 5126,
            "count":         len(n),
            "type":          "VEC3",
        })

        # INDICES
        ib   = indices.astype(np.uint32).tobytes()
        bv_i = _add_bv(ib, ELEM_BUF)
        acc_i = len(accessors)
        accessors.append({
            "bufferView":    bv_i,
            "componentType": 5125,   # UNSIGNED_INT
            "count":         len(indices),
            "type":          "SCALAR",
        })

        mesh_prims.append({
            "attributes": {"POSITION": acc_v, "NORMAL": acc_n},
            "indices":    acc_i,
            "material":   mat_idx,
            "mode":       4,   # TRIANGLES
        })

    gltf = {
        "asset":       {"version": "2.0", "generator": "ADU gen_gltf.py"},
        "scene":       0,
        "scenes":      [{"nodes": [0]}],
        "nodes":       [{"mesh": 0}],
        "meshes":      [{"primitives": mesh_prims}],
        "materials":   materials,
        "accessors":   accessors,
        "bufferViews": buffer_views,
        "buffers":     [{"byteLength": len(bin_buf)}],
    }

    json_bytes = json.dumps(gltf, separators=(',', ':')).encode('utf-8')
    while len(json_bytes) % 4:
        json_bytes += b' '
    while len(bin_buf) % 4:
        bin_buf.append(0)

    MAGIC   = 0x46546C67
    VERSION = 2
    total   = 12 + 8 + len(json_bytes) + 8 + len(bin_buf)

    out = bytearray()
    out += struct.pack('<III', MAGIC, VERSION, total)
    out += struct.pack('<II', len(json_bytes), 0x4E4F534A)   # "JSON"
    out += json_bytes
    out += struct.pack('<II', len(bin_buf), 0x004E4942)       # "BIN\0"
    out += bin_buf
    return bytes(out)


# ───────────────────────────────  Top-level  ─────────────────────────────────

def _load_scad(filename):
    """Read SCAD file and parse all geometry data."""
    path = os.path.join(_DIR, filename)
    with open(path, encoding='utf-8', errors='replace') as f:
        text = f.read()

    data = _parse_assembly(text)

    # Pre-load all t_* arrays and roof_outline
    arrays = {}
    for name, _, _ in data['walls']:
        if name not in arrays:
            arrays[name] = _extract_array(text, name)
    data['_arrays']      = arrays
    data['roof_outline'] = _extract_array(text, 'roof_outline')
    return data


def _write_model(scad_file, out_path, roof_builder):
    """Build and write a .glb for one model variant."""
    data = _load_scad(scad_file)
    s    = data['scalars']

    wall_mesh = MeshBuilder()
    roof_mesh = MeshBuilder()
    win_mesh  = MeshBuilder()
    iw_mesh   = MeshBuilder()

    # Outer walls (flat roof uses standard upper band; 2in12 overrides it)
    if scad_file == '2in12.scad':
        # Bands 0 and 1 are the same; band 2 (upper) is replaced by the
        # intersection-with-roof version
        for name, z_base, height in data['walls']:
            if name == 't_full_upper':
                continue          # handled by _build_walls_2in12_upper
            path = data['_arrays'].get(name)
            if path is None:
                continue
            half_t = s['half_t']
            outer  = shell_pts(path, -half_t)
            inner  = shell_pts(path,  half_t)
            extrude_ring(outer, inner, wall_mesh,
                         z_bot=z_base, z_top=z_base + height)
        _build_walls_2in12_upper(data, wall_mesh)
    else:
        _build_walls(data, wall_mesh)

    roof_builder(data, roof_mesh)
    _build_windows(data, win_mesh)
    _build_interior_walls(data, iw_mesh)

    primitives = []
    for mb, rgba, roughness in [
        (wall_mesh, WALL_CREAM, 0.85),
        (roof_mesh, ROOF_TEAL,  0.55),
        (win_mesh,  WINDOW_BG,  0.10),
        (iw_mesh,   WOOD_TAN,   0.90),
    ]:
        v, n, i = mb.to_arrays()
        if v is not None:
            primitives.append((v, n, i, rgba, roughness))

    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    data_bytes = _build_glb(primitives)
    with open(out_path, 'wb') as f:
        f.write(data_bytes)

    tri_count = sum(len(i) // 3 for _, _, i, _, _ in primitives)
    size_kb   = len(data_bytes) / 1024
    print(f"  wrote {out_path}  ({tri_count:,} tris, {size_kb:.0f} KB)")
    return data_bytes


def _write_models_js(flat_bytes, in12_bytes, out_path):
    """
    Write opts/models.js with both GLB files base64-encoded as JS globals.

    This lets index.html create Blob URLs at runtime — no fetch() required —
    so the page works when opened directly from the file system (file://).
    """
    import base64
    flat_b64 = base64.b64encode(flat_bytes).decode('ascii')
    in12_b64 = base64.b64encode(in12_bytes).decode('ascii')
    js = (
        "// Auto-generated by scad/gen_gltf.py - do not edit\n"
        "// GLB models embedded as base64 so the page works on file:// without a server.\n"
        f"window.GLB_FLAT_ROOF='{flat_b64}';\n"
        f"window.GLB_2IN12='{in12_b64}';\n"
    )
    with open(out_path, 'w', encoding='ascii') as f:
        f.write(js)
    size_kb = len(js) / 1024
    print(f"  wrote {out_path}  ({size_kb:.0f} KB)")


def generate(gd=None):
    """Generate both glTF models and the embedded JS bundle.  gd is unused."""
    os.makedirs(_OPTS, exist_ok=True)
    print("gen_gltf: generating glTF models...")
    flat_bytes = _write_model('flat_roof.scad',
                               os.path.join(_OPTS, 'flat_roof.glb'),
                               _build_roof_flat)
    in12_bytes = _write_model('2in12.scad',
                               os.path.join(_OPTS, '2in12.glb'),
                               _build_roof_2in12)
    _write_models_js(flat_bytes, in12_bytes,
                     os.path.join(_OPTS, 'models.js'))
    print("gen_gltf: done.")


if __name__ == '__main__':
    generate()
