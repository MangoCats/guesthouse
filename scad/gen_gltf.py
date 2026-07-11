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
    in block.  If skip_door_cutouts, skip entries where the 120 chars
    immediately preceding the polygon contain '-0.001' (door cutout marker).
    Uses a narrow window for the cutout check so the marker from a preceding
    difference() block doesn't bleed into the next wall's context.
    """
    result = []
    # Find every polygon(points = [...]) call
    for pm in re.finditer(r'polygon\s*\(\s*points\s*=\s*(\[)', block):
        arr_str = _bracket_array(block, pm.start(1))
        if arr_str is None:
            continue
        pts = ast.literal_eval(re.sub(r'//[^\n]*', '', arr_str))

        # Narrow context (120 chars) for door-cutout detection only
        near_ctx = block[max(0, pm.start() - 120): pm.start()]
        if skip_door_cutouts and '-0.001' in near_ctx:
            continue

        # Wider context (300 chars) for height and z-translate
        ctx = block[max(0, pm.start() - 300): pm.start()]

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


def _triangulate_polygon(pts_2d):
    """
    Ear-clipping triangulation for a simple polygon (convex or concave).
    Returns a list of (i, j, k) index triples whose winding matches the
    polygon orientation (CCW positive / CW negative area).
    """
    def _cross2(o, a, b):
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])

    def _pt_in_tri(p, a, b, c):
        d1 = _cross2(a, b, p)
        d2 = _cross2(b, c, p)
        d3 = _cross2(c, a, p)
        has_neg = (d1 < 0) or (d2 < 0) or (d3 < 0)
        has_pos = (d1 > 0) or (d2 > 0) or (d3 > 0)
        return not (has_neg and has_pos)

    pts = list(pts_2d)
    n = len(pts)
    if n < 3:
        return []
    if n == 3:
        return [(0, 1, 2)]

    # Signed area: positive = CCW, negative = CW
    area2 = sum(_cross2(pts[0], pts[i], pts[i + 1]) for i in range(1, n - 1))

    indices = list(range(n))
    triangles = []
    limit = n * n + 10  # safety iteration cap

    for _ in range(limit):
        m = len(indices)
        if m < 3:
            break
        if m == 3:
            triangles.append((indices[0], indices[1], indices[2]))
            break

        ear_found = False
        for i in range(m):
            pi = indices[(i - 1) % m]
            ci = indices[i]
            ni = indices[(i + 1) % m]
            a, b, c = pts[pi], pts[ci], pts[ni]

            # Ear tip must be convex w.r.t. polygon orientation
            cross = _cross2(a, b, c)
            if (area2 > 0 and cross <= 0) or (area2 < 0 and cross >= 0):
                continue

            # No other vertex may lie strictly inside the ear triangle
            is_ear = True
            for j in range(m):
                ji = indices[j]
                if ji in (pi, ci, ni):
                    continue
                if _pt_in_tri(pts[ji], a, b, c):
                    is_ear = False
                    break

            if is_ear:
                triangles.append((pi, ci, ni))
                indices.pop(i)
                ear_found = True
                break

        if not ear_found:
            break  # degenerate polygon — stop to avoid infinite loop

    return triangles


def extrude_polygon(pts_2d, mesh, z_bot=0.0, z_top=None,
                    top_fn=None, flip_top=False):
    """
    Extrude a simple polygon (convex or concave) into a prism.
    Ear-clipping triangulation is used for the caps.
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

    tris = _triangulate_polygon(pts_2d)

    # Bottom cap (reversed winding for downward-facing normal)
    for a, b, c in tris:
        mesh.add_tri(bot[a], bot[c], bot[b])

    # Top cap
    if flip_top:
        for a, b, c in tris:
            mesh.add_tri(top[a], top[c], top[b])
    else:
        for a, b, c in tris:
            mesh.add_tri(top[a], top[b], top[c])


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
    Band 5 of the 2in12 model: upper wall clipped by the roof shear plane.
    z_top per point = min(z_base + height, roof_slope * y + roof_z_off)
    where z_base and height come from the parsed 't_full_upper' wall entry.
    """
    s      = data['scalars']
    half_t = s['half_t']
    slope  = s['roof_slope']
    z_off  = s['roof_z_off']

    upper_entry = next((w for w in data['walls'] if w[0] == 't_full_upper'), None)
    if upper_entry is None:
        return
    _, z_base, height = upper_entry
    z_max = z_base + height

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

    tris = _triangulate_polygon(pts)
    # Bottom cap (reversed winding for downward-facing normal)
    for a, b, c in tris:
        roof_mesh.add_tri(bot3[a], bot3[c], bot3[b])
    # Top cap
    for a, b, c in tris:
        roof_mesh.add_tri(top3[a], top3[b], top3[c])


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

# Variant definitions: (display_name, config_db_stem, roof_type)
# roof_type must be 'flat_roof' or '2in12'
_VARIANTS = [
    ("Base",  "Base",  "2in12"),
    ("Mark6", "Mark6", "2in12"),
    ("Mark8", "Mark8", "flat_roof"),
    ("Mark9", "Mark9", "flat_roof"),
    ("MarkX", "MarkX", "flat_roof"),
]

_ROOF_BUILDERS = {
    "flat_roof": _build_roof_flat,
    "2in12":     _build_roof_2in12,
    "1in12":     _build_roof_2in12,  # same slab math; slope comes from scalars
}


def _parse_scad_text(text):
    """Parse SCAD text into a geometry data dict."""
    data = _parse_assembly(text)
    arrays = {}
    for name, _, _ in data['walls']:
        if name not in arrays:
            arrays[name] = _extract_array(text, name)
    data['_arrays']      = arrays
    data['roof_outline'] = _extract_array(text, 'roof_outline')
    return data


def _load_scad(filename):
    """Read a SCAD file from the scad/ directory and parse it."""
    path = os.path.join(_DIR, filename)
    with open(path, encoding='utf-8', errors='replace') as f:
        return _parse_scad_text(f.read())


def _build_model_bytes(data, roof_type):
    """Build GLB bytes from a parsed geometry data dict."""
    s            = data['scalars']
    roof_builder = _ROOF_BUILDERS[roof_type]

    wall_mesh = MeshBuilder()
    roof_mesh = MeshBuilder()
    win_mesh  = MeshBuilder()
    iw_mesh   = MeshBuilder()

    if roof_type in ('2in12', '1in12'):
        # Upper wall (band 2) is clipped by the roof slope
        for name, z_base, height in data['walls']:
            if name == 't_full_upper':
                continue          # handled separately below
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

    return _build_glb(primitives)


def _generate_variant_scad(var_gd, roof_type):
    """
    Call the appropriate SCAD generator with var_gd, read back the written
    SCAD text, and return it.  The generator writes to the normal .scad path.
    """
    if roof_type == 'flat_roof':
        import scad.gen_flat_roof as _gen
        scad_file = 'flat_roof.scad'
    elif roof_type == 'split2':
        import scad.gen_split2 as _gen
        scad_file = 'split2.scad'
    elif roof_type == '1in12':
        import scad.gen_1in12 as _gen
        scad_file = '1in12.scad'
    else:
        import scad.gen_2in12 as _gen
        scad_file = '2in12.scad'
    _gen.generate(var_gd)
    with open(os.path.join(_DIR, scad_file), encoding='utf-8', errors='replace') as f:
        return f.read()


def _write_models_js(variant_glbs, out_path):
    """
    Write opts/models.js with each variant's GLB base64-encoded as a JS global.

    This lets index.html create Blob URLs at runtime — no fetch() required —
    so the page works when opened directly from the file system (file://).
    """
    import base64
    lines = [
        "// Auto-generated by scad/gen_gltf.py - do not edit",
        "// GLB models embedded as base64 so the page works on file:// without a server.",
    ]
    for name, glb_bytes in variant_glbs.items():
        b64 = base64.b64encode(glb_bytes).decode('ascii')
        lines.append(f"window.GLB_{name.upper()}='{b64}';")
    js = "\n".join(lines) + "\n"
    with open(out_path, 'w', encoding='ascii') as f:
        f.write(js)
    print(f"  wrote {out_path}  ({len(js)/1024:.0f} KB)")


def generate(gd=None):
    """
    Generate one .glb per variant (Base/Mark6/Mark8/Mark9/MarkX) by running
    the appropriate SCAD generator against each variant's config DB, then
    bundle all models into opts/models.js as base64 globals.
    """
    sys.path.insert(0, _ROOT)
    os.makedirs(_OPTS, exist_ok=True)
    print("gen_gltf: generating per-variant glTF models...")

    configs_dir = os.path.join(_ROOT, "app", "configs")

    # Read current SCAD files so we can restore them after per-variant generation
    scad_originals = {}
    for scad_name in ("flat_roof.scad", "2in12.scad"):
        p = os.path.join(_DIR, scad_name)
        if os.path.exists(p):
            with open(p, encoding='utf-8', errors='replace') as f:
                scad_originals[scad_name] = f.read()

    try:
        from app.engine import build_generator_data_from_db
    except ImportError:
        build_generator_data_from_db = None

    variant_glbs = {}

    for variant_name, db_stem, roof_type in _VARIANTS:
        db_path = os.path.join(configs_dir, f"{db_stem}.db")
        scad_name = f"{roof_type}.scad"

        if build_generator_data_from_db and os.path.exists(db_path):
            try:
                var_gd    = build_generator_data_from_db(db_path)
                scad_text = _generate_variant_scad(var_gd, roof_type)
            except Exception as exc:
                print(f"  {variant_name}: SCAD generation failed ({exc}), using active-DB SCAD")
                scad_text = scad_originals.get(scad_name, "")
        elif scad_originals.get(scad_name):
            # Fallback: use the pre-generated SCAD from the active DB
            print(f"  {variant_name}: no DB at {db_path}, using active-DB SCAD")
            scad_text = scad_originals[scad_name]
        else:
            print(f"  {variant_name}: skipped (no DB and no SCAD)")
            continue

        data       = _parse_scad_text(scad_text)
        glb_bytes  = _build_model_bytes(data, roof_type)
        out_path   = os.path.join(_OPTS, f"{variant_name}.glb")
        with open(out_path, 'wb') as f:
            f.write(glb_bytes)
        tris = _count_tris(glb_bytes)
        print(f"  {variant_name} ({roof_type}): {tris:,} tris, "
              f"{len(glb_bytes)/1024:.0f} KB -> {os.path.basename(out_path)}")
        variant_glbs[variant_name] = glb_bytes

    # Restore SCAD files to their pre-loop state
    for scad_name, text in scad_originals.items():
        with open(os.path.join(_DIR, scad_name), 'w', encoding='utf-8') as f:
            f.write(text)

    _write_models_js(variant_glbs, os.path.join(_OPTS, 'models.js'))
    print("gen_gltf: done.")


def _count_tris(glb_bytes):
    """Quick triangle count from a GLB file's accessor metadata."""
    import struct as _s, json as _j
    jlen = _s.unpack_from('<I', glb_bytes, 12)[0]
    gltf = _j.loads(glb_bytes[20: 20 + jlen])
    total = 0
    for prim in gltf['meshes'][0]['primitives']:
        total += gltf['accessors'][prim['indices']]['count'] // 3
    return total


if __name__ == '__main__':
    generate()
