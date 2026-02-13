"""Generate outer wall detail SVG showing double-shell concrete construction.

Outer walls are double-shell 3D-printed concrete: two 2" shells separated
by a 4" air gap (8" total = WALL_OUTER). At openings, shells connect via
90-degree corner turns with configurable radii.

Outputs walls/walls.svg at 1:72 scale.
"""
import os, sys, math, subprocess, datetime
from typing import NamedTuple

# Ensure project root is on sys.path for package imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from shared.types import Point, LineSeg, ArcSeg, Segment
from shared.geometry import (
    segment_polyline, path_polygon, poly_area, arc_poly,
    compute_inner_walls, fmt_dist, left_norm,
)
from shared.svg import make_svg_transform, W, H
from floorplan.gen_floorplan import build_floorplan_data
from floorplan.layout import compute_interior_layout
from floorplan.constants import (
    WALL_OUTER, IW1_OFFSET_N, WALL_6IN, IW2_OFFSET_E,
    IW5_OFFSET_N,
)
from walls.constants import SHELL_THICKNESS, AIR_GAP, OPENING_INSIDE_RADIUS


# ============================================================
# Data types
# ============================================================

class WallOpening(NamedTuple):
    """Opening on an outline segment, parameterized along the segment."""
    name: str
    seg_idx: int    # index in outline_segs (0-based)
    t_start: float  # parametric position [0, 1] along the segment
    t_end: float    # parametric position [0, 1] along the segment


# ============================================================
# Shell path computation
# ============================================================

def _compute_inset_path(outline_segs, pts, radii, inset, prefix):
    """Compute a shell boundary path at given inset distance.

    Returns (new_pts_dict, new_segs) with point names using the given prefix
    (e.g., "S0".."S21" for prefix="S").
    """
    tmp_pts = dict(pts)
    tmp_segs = compute_inner_walls(outline_segs, tmp_pts, inset, radii)

    result_pts = {}
    for i in range(22):
        result_pts[f"{prefix}{i}"] = tmp_pts[f"W{i}"]

    result_segs = []
    for seg in tmp_segs:
        if isinstance(seg, LineSeg):
            s = prefix + seg.start[1:]
            e = prefix + seg.end[1:]
            result_segs.append(LineSeg(s, e))
        else:
            s = prefix + seg.start[1:]
            e = prefix + seg.end[1:]
            result_segs.append(ArcSeg(s, e, seg.center, seg.radius,
                                       seg.direction, seg.n_pts))
    return result_pts, result_segs


# ============================================================
# Opening computation
# ============================================================

def _seg_param(pts, seg, point):
    """Compute parametric position t of a point along a LineSeg."""
    A = pts[seg.start]
    B = pts[seg.end]
    dx = B[0] - A[0]
    dy = B[1] - A[1]
    if abs(dx) > abs(dy):
        return (point[0] - A[0]) / dx
    else:
        return (point[1] - A[1]) / dy


def _compute_openings(pts, outline_segs, layout):
    """Compute opening boundaries as WallOpening tuples.

    Each opening maps to a specific outline segment index and parametric range.
    """
    openings = []

    # Segment index lookup: outline_segs[i] has start = F{i}, end = F{(i+1)%22}
    # Build a map from (start, end) pair to index
    seg_map = {}
    for i, seg in enumerate(outline_segs):
        seg_map[(seg.start, seg.end)] = i

    # Openings numbered CW around the building outline (F1→F2→...→F21→F0)
    iw1_s = pts["W0"][1] + IW1_OFFSET_N
    iw1_n = iw1_s + WALL_6IN
    iw5_n = iw1_s - IW5_OFFSET_N
    iw5_s = iw5_n - WALL_6IN

    # --- O2: F1-F2, vertical, upper (near F2) — computed first, O1 depends on it ---
    idx = seg_map[("F1", "F2")]
    o2_n = pts["F2"][1] - 8.0 / 12.0
    o2_s = pts["F2"][1] - 33.0 / 12.0
    t1 = _seg_param(pts, outline_segs[idx], (pts["F2"][0], o2_s))
    t2 = _seg_param(pts, outline_segs[idx], (pts["F2"][0], o2_n))
    openings.append(WallOpening("O2", idx, min(t1, t2), max(t1, t2)))

    # --- O1: F1-F2, vertical, lower (south of IW1) ---
    o1_gap = o2_s - iw1_n
    o1_n = iw1_s - o1_gap
    o1_s = o1_n - 25.0 / 12.0
    t1 = _seg_param(pts, outline_segs[idx], (pts["F2"][0], o1_s))
    t2 = _seg_param(pts, outline_segs[idx], (pts["F2"][0], o1_n))
    openings.append(WallOpening("O1", idx, min(t1, t2), max(t1, t2)))

    # --- O3: F4-F5, vertical ---
    idx = seg_map[("F4", "F5")]
    o3_cn = (pts["F4"][1] + pts["F5"][1]) / 2
    o3_half = 16.0 / 12.0
    t1 = _seg_param(pts, outline_segs[idx], (pts["F4"][0], o3_cn - o3_half))
    t2 = _seg_param(pts, outline_segs[idx], (pts["F4"][0], o3_cn + o3_half))
    openings.append(WallOpening("O3", idx, min(t1, t2), max(t1, t2)))

    # --- O4: F6-F7, horizontal ---
    idx = seg_map[("F6", "F7")]
    iw2_w = layout.iw2_w
    o4_e = iw2_w - 18.0 / 12.0
    o4_w = o4_e - 9.0 / 12.0
    t1 = _seg_param(pts, outline_segs[idx], (o4_w, pts["F6"][1]))
    t2 = _seg_param(pts, outline_segs[idx], (o4_e, pts["F6"][1]))
    openings.append(WallOpening("O4", idx, min(t1, t2), max(t1, t2)))

    # --- O5: F9-F10, horizontal, big 6' opening ---
    idx = seg_map[("F9", "F10")]
    o5_w = pts["F9"][0] + 8.0 / 12.0
    o5_e = o5_w + 6.0
    t1 = _seg_param(pts, outline_segs[idx], (o5_w, pts["F9"][1]))
    t2 = _seg_param(pts, outline_segs[idx], (o5_e, pts["F9"][1]))
    openings.append(WallOpening("O5", idx, min(t1, t2), max(t1, t2)))

    # --- O6: F9-F10, horizontal, small opening ---
    o6_e = pts["F10"][0] - 4.0 / 12.0
    o6_w = pts["F10"][0] - 48.0 / 12.0
    t1 = _seg_param(pts, outline_segs[idx], (o6_w, pts["F9"][1]))
    t2 = _seg_param(pts, outline_segs[idx], (o6_e, pts["F9"][1]))
    openings.append(WallOpening("O6", idx, min(t1, t2), max(t1, t2)))

    # --- O7: F12-F13, diagonal ---
    idx = seg_map[("F12", "F13")]
    dE = pts["F13"][0] - pts["F12"][0]
    dN = pts["F13"][1] - pts["F12"][1]
    seg_len = math.sqrt(dE**2 + dN**2)
    o7_half = 36.0 / 12.0
    t1 = 0.5 - o7_half / seg_len
    t2 = 0.5 + o7_half / seg_len
    openings.append(WallOpening("O7", idx, t1, t2))

    # --- O8: F14-F15, vertical ---
    idx = seg_map[("F14", "F15")]
    o8_cn = (iw5_s + pts["F15"][1]) / 2
    o8_half = 12.5 / 12.0
    t1 = _seg_param(pts, outline_segs[idx], (pts["F15"][0], o8_cn - o8_half))
    t2 = _seg_param(pts, outline_segs[idx], (pts["F15"][0], o8_cn + o8_half))
    openings.append(WallOpening("O8", idx, min(t1, t2), max(t1, t2)))

    # --- O9: F18-F19, horizontal ---
    idx = seg_map[("F18", "F19")]
    o9_cn = (layout.bed_e + layout.iw4_w) / 2
    o9_half = 12.5 / 12.0
    t1 = _seg_param(pts, outline_segs[idx], (o9_cn - o9_half, pts["F18"][1]))
    t2 = _seg_param(pts, outline_segs[idx], (o9_cn + o9_half, pts["F18"][1]))
    openings.append(WallOpening("O9", idx, min(t1, t2), max(t1, t2)))

    # --- O10: F21-F0, horizontal (bed area) ---
    idx = seg_map[("F21", "F0")]
    o10_cn = (layout.bed_w + layout.iw3_e) / 2
    o10_half = 12.5 / 12.0
    t1 = _seg_param(pts, outline_segs[idx], (o10_cn - o10_half, pts["F0"][1]))
    t2 = _seg_param(pts, outline_segs[idx], (o10_cn + o10_half, pts["F0"][1]))
    openings.append(WallOpening("O10", idx, min(t1, t2), max(t1, t2)))

    # --- O11: F21-F0, horizontal (utility area) ---
    o11_cn = (layout.dryer_e + layout.ctr_w) / 2
    o11_half = 12.5 / 12.0
    t1 = _seg_param(pts, outline_segs[idx], (o11_cn - o11_half, pts["F0"][1]))
    t2 = _seg_param(pts, outline_segs[idx], (o11_cn + o11_half, pts["F0"][1]))
    openings.append(WallOpening("O11", idx, min(t1, t2), max(t1, t2)))

    return openings


def _openings_on_seg(openings, seg_idx):
    """Get openings on a given segment, sorted by t_start."""
    result = [o for o in openings if o.seg_idx == seg_idx]
    result.sort(key=lambda o: o.t_start)
    return result


def _solid_ranges(seg_openings):
    """Compute solid wall parametric ranges from sorted opening list.

    Returns list of (t_start, t_end) for solid wall sections.
    """
    ranges = []
    cursor = 0.0
    for o in seg_openings:
        if o.t_start > cursor + 1e-9:
            ranges.append((cursor, o.t_start))
        cursor = o.t_end
    if cursor < 1.0 - 1e-9:
        ranges.append((cursor, 1.0))
    return ranges


# ============================================================
# Polygon builders
# ============================================================

def _lerp(a, b, t):
    """Linear interpolation between two points."""
    return (a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1]))


def _arc_strip_poly(seg, pts, outer_prefix, inner_seg):
    """Build polygon for an arc shell strip (full segment).

    Returns list of (E, N) points forming the closed polygon.
    """
    outer_poly = segment_polyline(seg, pts)
    inner_poly = segment_polyline(inner_seg, pts)
    # Forward along outer, backward along inner
    return outer_poly + list(reversed(inner_poly))


def _line_strip_poly(pts, seg_start, seg_end, inner_start, inner_end):
    """Build polygon for a line shell strip (full segment or sub-range).

    4-point rectangle: start_outer, end_outer, end_inner, start_inner.
    """
    return [
        pts[seg_start], pts[seg_end],
        pts[inner_end], pts[inner_start],
    ]


def _partial_line_strip(pts, seg, inner_seg_start, inner_seg_end, t_start, t_end):
    """Build polygon for a partial line shell strip between t_start and t_end."""
    A_out = pts[seg.start]
    B_out = pts[seg.end]
    A_in = pts[inner_seg_start]
    B_in = pts[inner_seg_end]
    p1 = _lerp(A_out, B_out, t_start)
    p2 = _lerp(A_out, B_out, t_end)
    p3 = _lerp(A_in, B_in, t_end)
    p4 = _lerp(A_in, B_in, t_start)
    return [p1, p2, p3, p4]



def _uturn_arc_data(pts, outline_segs, inner_segs, seg_idx, t_param, side,
                    shell_t, R_in, wall_total, n_arc=12):
    """Compute U-turn arc point arrays at an opening boundary.

    Returns a dict with keys:
      'oc_F': F-face arc points (outer shell, outer face)
      'oc_S': S-face arc points (outer shell, inner face)
      'ic_W': W-face arc points (inner shell, inner face)
      'ic_G': G-face arc points (inner shell, outer face)

    Each arc goes from the shell face toward the cross-wall face:
      oc_F[0]/oc_S[0] = on F/S-face, R_out back from opening boundary
      oc_F[-1]/oc_S[-1] = on cross-wall face
      ic_W[0]/ic_G[0] = on cross-wall face
      ic_W[-1]/ic_G[-1] = on W/G-face, R_out back from opening boundary

    side: "start" means wall-to-opening transition (wall at t < t_param)
          "end" means opening-to-wall transition (wall at t > t_param)
    """
    R_out = R_in + shell_t
    seg = outline_segs[seg_idx]

    # Points at the boundary parameter on outer and inner faces
    F_A, F_B = pts[seg.start], pts[seg.end]
    F_pt = _lerp(F_A, F_B, t_param)

    inner_seg = inner_segs[seg_idx]
    W_A, W_B = pts[inner_seg.start], pts[inner_seg.end]
    W_pt = _lerp(W_A, W_B, t_param)

    # Tangent direction (along wall, CW traversal)
    dx, dy = F_B[0] - F_A[0], F_B[1] - F_A[1]
    t_len = math.sqrt(dx * dx + dy * dy)
    t_hat = (dx / t_len, dy / t_len)

    # Exterior normal (left of CW traversal direction)
    n_ext = (-t_hat[1], t_hat[0])

    # Direction toward the opening along the wall
    if side == "start":
        open_dir = t_hat
    else:
        open_dir = (-t_hat[0], -t_hat[1])
    wall_dir = (-open_dir[0], -open_dir[1])

    # --- Arc centers ---
    # Outer shell: R_out inward from F-face, R_out back from opening
    oc = (F_pt[0] - R_out * n_ext[0] + R_out * wall_dir[0],
          F_pt[1] - R_out * n_ext[1] + R_out * wall_dir[1])
    # Inner shell: R_out outward from W-face, R_out back from opening
    ic = (W_pt[0] + R_out * n_ext[0] + R_out * wall_dir[0],
          W_pt[1] + R_out * n_ext[1] + R_out * wall_dir[1])

    # Quarter-circle arc: center + R*(cos(a)*u0 + sin(a)*u1), a from 0 to pi/2
    def qarc(cx, cy, R, u0, u1):
        arc_pts = []
        for i in range(n_arc + 1):
            a = i * math.pi / (2 * n_arc)
            ca, sa = math.cos(a), math.sin(a)
            arc_pts.append((cx + R * (ca * u0[0] + sa * u1[0]),
                            cy + R * (ca * u0[1] + sa * u1[1])))
        return arc_pts

    # Outer shell arcs: from shell face (n_ext) to cross-wall (open_dir)
    oc_F = qarc(oc[0], oc[1], R_out, n_ext, open_dir)  # F-face arc
    oc_S = qarc(oc[0], oc[1], R_in, n_ext, open_dir)   # S-face arc

    # Inner shell arcs: from cross-wall (open_dir) to shell face (-n_ext)
    n_int = (-n_ext[0], -n_ext[1])
    ic_W = qarc(ic[0], ic[1], R_out, open_dir, n_int)   # W-face arc
    ic_G = qarc(ic[0], ic[1], R_in, open_dir, n_int)    # G-face arc

    return {'oc_F': oc_F, 'oc_S': oc_S, 'ic_W': ic_W, 'ic_G': ic_G}


def _uturn_polygon(pts, outline_segs, inner_segs, s_segs, g_segs,
                   seg_idx, t_param, side, shell_t, R_in, wall_total,
                   n_arc=12):
    """Build the U-turn polygon at an opening boundary.

    The turn connects the outer shell to the inner shell via two 90-degree
    arcs and a straight cross-wall section, curving toward building interior.

        F-face --,              R_out = R_in + shell_t
                 |
        S-face -,|              R_in
                ||
                || straight (wall_total - 2*(shell_t + R_in))
                ||
        G-face -'|              R_in
                 |
        W-face --'              R_out

    side: "start" means wall-to-opening transition (wall at t < t_param)
          "end" means opening-to-wall transition (wall at t > t_param)

    Returns a list of (E, N) points forming the U-turn polygon.
    """
    arcs = _uturn_arc_data(pts, outline_segs, inner_segs, seg_idx, t_param,
                           side, shell_t, R_in, wall_total, n_arc)

    # Assemble: outer profile forward, inner profile reversed
    poly = []
    poly.extend(arcs['oc_F'])           # F-face arc (shell -> cross-wall)
    # implicit straight: cross-wall outer face
    poly.extend(arcs['ic_W'])           # W-face arc (cross-wall -> shell)
    poly.extend(reversed(arcs['ic_G'])) # G-face arc reversed
    # implicit straight: cross-wall inner face
    poly.extend(reversed(arcs['oc_S'])) # S-face arc reversed
    return poly


# ============================================================
# Continuous outline builders
# ============================================================

def _trace_boundary_path(pts, segs, start_seg_idx, start_t, end_seg_idx,
                         end_t, R_out):
    """Trace a boundary path between two opening boundaries across segments.

    Traces CW from (start_seg_idx, start_t + delta_t) to
    (end_seg_idx, end_t - delta_t), spanning intermediate segments.

    segs: one of outline_segs/s_segs/g_segs/inner_segs (all 22 segments).
    start_t: parametric position of the starting opening's t_end.
    end_t: parametric position of the ending opening's t_start.
    R_out: trim distance in feet (R_in + shell_t) — converted to delta_t
           using each line segment's length.

    Returns list of (E, N) points along the boundary.
    """
    n_segs = len(segs)
    path = []

    if start_seg_idx == end_seg_idx:
        # Same segment — just two interpolated points
        seg = segs[start_seg_idx]
        A, B = pts[seg.start], pts[seg.end]
        seg_len = math.sqrt((B[0] - A[0])**2 + (B[1] - A[1])**2)
        dt = R_out / seg_len
        path.append(_lerp(A, B, start_t + dt))
        path.append(_lerp(A, B, end_t - dt))
        return path

    # Multi-segment: partial start + full intermediates + partial end

    # Start segment (from start_t + delta to segment end)
    seg = segs[start_seg_idx]
    A, B = pts[seg.start], pts[seg.end]
    seg_len = math.sqrt((B[0] - A[0])**2 + (B[1] - A[1])**2)
    dt = R_out / seg_len
    path.append(_lerp(A, B, start_t + dt))
    if isinstance(seg, ArcSeg):
        poly = segment_polyline(seg, pts)
        # Find closest point index to our start position and take rest
        # For arcs, the start_t + dt position is near the end of the arc
        # Since openings are only on LineSegs, this shouldn't happen
        path.append(B)
    else:
        path.append(B)

    # Intermediate full segments
    idx = (start_seg_idx + 1) % n_segs
    while idx != end_seg_idx:
        seg = segs[idx]
        if isinstance(seg, ArcSeg):
            poly = segment_polyline(seg, pts)
            path.extend(poly[1:])  # skip first (matches previous end)
        else:
            path.append(pts[seg.end])  # start matches previous end
        idx = (idx + 1) % n_segs

    # End segment (from segment start to end_t - delta)
    seg = segs[end_seg_idx]
    A, B = pts[seg.start], pts[seg.end]
    seg_len = math.sqrt((B[0] - A[0])**2 + (B[1] - A[1])**2)
    dt = R_out / seg_len
    path.append(_lerp(A, B, end_t - dt))

    return path


def _enumerate_wall_sections(openings, outline_segs):
    """Enumerate wall sections as (start_opening, end_opening) pairs.

    Each wall section is bounded by start_opening.t_end on one side and
    end_opening.t_start on the other, going CW around the building.
    Returns list of (start_op, end_op) tuples.
    """
    # Collect all opening boundaries in CW order
    ordered = []
    for seg_idx in range(len(outline_segs)):
        seg_ops = _openings_on_seg(openings, seg_idx)
        for op in seg_ops:
            ordered.append(op)

    # Wall sections go from each opening's t_end to the next opening's t_start
    sections = []
    n = len(ordered)
    for i in range(n):
        start_op = ordered[i]
        end_op = ordered[(i + 1) % n]
        sections.append((start_op, end_op))
    return sections


def _build_section_outlines(pts, outline_segs, inner_segs, s_segs, g_segs,
                            start_op, end_op, shell_t, R_in, wall_total,
                            n_arc=12):
    """Build outer and inner cavity outlines for one wall section.

    start_op: opening whose t_end starts this wall section
    end_op: opening whose t_start ends this wall section

    Returns (outer_path, cavity_path) as lists of (E, N) points.
    """
    R_out = R_in + shell_t

    # U-turn arc data at each end
    start_arcs = _uturn_arc_data(pts, outline_segs, inner_segs,
                                 start_op.seg_idx, start_op.t_end, "end",
                                 shell_t, R_in, wall_total, n_arc)
    end_arcs = _uturn_arc_data(pts, outline_segs, inner_segs,
                               end_op.seg_idx, end_op.t_start, "start",
                               shell_t, R_in, wall_total, n_arc)

    # Trace boundary paths between the two openings
    f_path = _trace_boundary_path(pts, outline_segs,
                                  start_op.seg_idx, start_op.t_end,
                                  end_op.seg_idx, end_op.t_start, R_out)
    s_path = _trace_boundary_path(pts, s_segs,
                                  start_op.seg_idx, start_op.t_end,
                                  end_op.seg_idx, end_op.t_start, R_out)
    g_path = _trace_boundary_path(pts, g_segs,
                                  start_op.seg_idx, start_op.t_end,
                                  end_op.seg_idx, end_op.t_start, R_out)
    w_path = _trace_boundary_path(pts, inner_segs,
                                  start_op.seg_idx, start_op.t_end,
                                  end_op.seg_idx, end_op.t_start, R_out)

    # --- Outer outline ---
    # F-face forward → end U-turn (F→W) → W-face backward → start U-turn (W→F)
    outer = list(f_path)
    outer.extend(end_arcs['oc_F'][1:])          # F-face arc at end
    outer.extend(end_arcs['ic_W'])              # cross-wall → W-face at end
    outer.extend(list(reversed(w_path)))        # W-face backward
    outer.extend(list(reversed(start_arcs['ic_W'])))  # W-face → cross-wall at start
    outer.extend(list(reversed(start_arcs['oc_F']))[1:])  # cross-wall → F-face at start

    # --- Cavity outline ---
    # S-face forward → end inner arcs → G-face backward → start inner arcs
    cavity = list(s_path)
    cavity.extend(end_arcs['oc_S'][1:])         # S-face arc at end
    cavity.append(end_arcs['ic_G'][0])          # cross-wall cavity face
    cavity.extend(end_arcs['ic_G'][1:])         # G-face arc at end
    cavity.extend(list(reversed(g_path)))       # G-face backward
    r_start_icG = list(reversed(start_arcs['ic_G']))
    cavity.extend(r_start_icG[1:])              # G-face arc at start (reversed)
    cavity.append(start_arcs['oc_S'][-1])       # cross-wall cavity face
    r_start_ocS = list(reversed(start_arcs['oc_S']))
    cavity.extend(r_start_ocS[1:])              # S-face arc at start (reversed)

    return outer, cavity


# ============================================================
# Data computation
# ============================================================

def build_wall_data():
    """Compute all geometry needed for the wall detail SVG."""
    fp_data = build_floorplan_data()
    pts = fp_data["pts"]
    to_svg = fp_data["to_svg"]
    outline_segs = fp_data["outline_segs"]
    inner_segs = fp_data["inner_segs"]
    radii = fp_data["radii"]
    inner_poly = fp_data["inner_poly"]

    # Compute S-series (2" inset = inner face of outer shell)
    s_pts, s_segs = _compute_inset_path(outline_segs, pts, radii,
                                         SHELL_THICKNESS, "S")
    pts.update(s_pts)

    # Compute G-series (6" inset = outer face of inner shell)
    g_pts, g_segs = _compute_inset_path(outline_segs, pts, radii,
                                         SHELL_THICKNESS + AIR_GAP, "G")
    pts.update(g_pts)

    # Compute interior layout (needed for opening positions)
    layout = compute_interior_layout(pts, inner_poly)

    # Compute openings
    openings = _compute_openings(pts, outline_segs, layout)

    # --- Page layout: 1:72 scale ---
    _f_svg = [to_svg(*pts[f"F{i}"]) for i in range(22)]
    _bldg_xmin = min(p[0] for p in _f_svg)
    _bldg_xmax = max(p[0] for p in _f_svg)
    _bldg_ymin = min(p[1] for p in _f_svg)
    _bldg_ymax = max(p[1] for p in _f_svg)
    _bldg_cx = (_bldg_xmin + _bldg_xmax) / 2

    _title_y = _bldg_ymin - 49

    _tb_w, _tb_h = 130, 80
    _tb_left = _bldg_xmax + 10
    _tb_right = _tb_left + _tb_w
    _tb_top = _title_y - 14 * 0.35
    _tb_bottom = _tb_top + _tb_h
    _tb_cx = (_tb_left + _tb_right) / 2

    _scale_ratio = 72
    _ft_per_inch = _scale_ratio / 12.0
    _svg_per_ft = to_svg(1, 0)[0] - to_svg(0, 0)[0]
    _fit_scale = 72.0 / (_ft_per_inch * _svg_per_ft)

    _margin_top = 36
    _vb_w = W / _fit_scale
    _vb_h = H / _fit_scale
    _cb_xmin = _bldg_xmin - 25
    _cb_xmax = _tb_right + 5
    _cb_cx = (_cb_xmin + _cb_xmax) / 2
    _cb_ymin = _title_y - 14 - 5
    _vb_x = _cb_cx - _vb_w / 2
    _vb_y = _cb_ymin - _margin_top / _fit_scale

    return {
        "pts": pts, "to_svg": to_svg,
        "outline_segs": outline_segs, "inner_segs": inner_segs,
        "s_segs": s_segs, "g_segs": g_segs,
        "radii": radii, "openings": openings,
        "outer_area": fp_data["outer_area"],
        "inner_area": fp_data["inner_area"],
        "vb_x": _vb_x, "vb_y": _vb_y, "vb_w": _vb_w, "vb_h": _vb_h,
        "title_x": _bldg_cx, "title_y": _title_y,
        "tb_left": _tb_left, "tb_right": _tb_right, "tb_top": _tb_top,
        "tb_bottom": _tb_bottom, "tb_w": _tb_w, "tb_h": _tb_h, "tb_cx": _tb_cx,
        "ft_per_inch": _ft_per_inch,
    }


# ============================================================
# SVG rendering
# ============================================================

def _svg_polygon(out, poly, to_svg, fill, stroke="#666", stroke_width="0.5"):
    """Render a polygon as an SVG element."""
    svg = " ".join(f"{to_svg(*p)[0]:.2f},{to_svg(*p)[1]:.2f}" for p in poly)
    out.append(f'<polygon points="{svg}" fill="{fill}" '
               f'stroke="{stroke}" stroke-width="{stroke_width}"/>')


def render_walls_svg(data):
    """Render the wall detail SVG. Returns SVG string."""
    pts = data["pts"]
    to_svg = data["to_svg"]
    outline_segs = data["outline_segs"]
    inner_segs = data["inner_segs"]
    s_segs = data["s_segs"]
    g_segs = data["g_segs"]
    openings = data["openings"]

    shell_t = SHELL_THICKNESS
    R_in = OPENING_INSIDE_RADIUS
    R_out = R_in + shell_t

    out = []
    out.append(f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}"'
               f' viewBox="{data["vb_x"]:.2f} {data["vb_y"]:.2f}'
               f' {data["vb_w"]:.2f} {data["vb_h"]:.2f}">')
    out.append(f'<rect x="{data["vb_x"]:.2f}" y="{data["vb_y"]:.2f}"'
               f' width="{data["vb_w"]:.2f}" height="{data["vb_h"]:.2f}" fill="white"/>')

    # Title
    out.append(f'<text x="{data["title_x"]:.1f}" y="{data["title_y"]:.1f}"'
               f' text-anchor="middle" font-family="Arial" font-size="14"'
               f' font-weight="bold">Outer Walls</text>')

    WALL_FILL = "rgba(180,180,180,0.5)"
    OPENING_FILL = "rgb(220,235,255)"

    # --- Draw wall sections ---
    for seg_idx in range(22):
        seg = outline_segs[seg_idx]
        inner_seg = inner_segs[seg_idx]
        s_seg = s_segs[seg_idx]
        g_seg = g_segs[seg_idx]

        seg_openings = _openings_on_seg(openings, seg_idx)

        if isinstance(seg, ArcSeg):
            # Arc segments have no openings — draw full strips
            # Outer shell: F-arc to S-arc
            outer_shell = _arc_strip_poly(seg, pts, "F", s_seg)
            _svg_polygon(out, outer_shell, to_svg, WALL_FILL, stroke="none")

            # Inner shell: G-arc to W-arc
            inner_shell = _arc_strip_poly(g_seg, pts, "G", inner_seg)
            _svg_polygon(out, inner_shell, to_svg, WALL_FILL, stroke="none")

        elif isinstance(seg, LineSeg):
            if not seg_openings:
                # No openings — draw full rectangle strips
                outer_strip = _line_strip_poly(pts, seg.start, seg.end,
                                               s_seg.start, s_seg.end)
                _svg_polygon(out, outer_strip, to_svg, WALL_FILL, stroke="none")

                inner_strip = _line_strip_poly(pts, g_seg.start, g_seg.end,
                                               inner_seg.start, inner_seg.end)
                _svg_polygon(out, inner_strip, to_svg, WALL_FILL, stroke="none")
            else:
                # Has openings — draw solid sections and U-turns
                solid_ranges = _solid_ranges(seg_openings)

                # Shrink ranges so shells end where U-turn arcs begin
                F_A, F_B = pts[seg.start], pts[seg.end]
                seg_len = math.sqrt((F_B[0]-F_A[0])**2 +
                                    (F_B[1]-F_A[1])**2)
                delta_t = R_out / seg_len
                adjusted = []
                for t_s, t_e in solid_ranges:
                    if t_s > 1e-9:   # borders an opening end
                        t_s += delta_t
                    if t_e < 1.0 - 1e-9:  # borders an opening start
                        t_e -= delta_t
                    if t_e > t_s + 1e-9:
                        adjusted.append((t_s, t_e))
                solid_ranges = adjusted

                for t_s, t_e in solid_ranges:
                    # Outer shell partial strip
                    outer_strip = _partial_line_strip(
                        pts, seg, s_seg.start, s_seg.end, t_s, t_e)
                    _svg_polygon(out, outer_strip, to_svg, WALL_FILL, stroke="none")

                    # Inner shell partial strip
                    inner_strip = _partial_line_strip_2(
                        pts, g_seg, inner_seg, t_s, t_e)
                    _svg_polygon(out, inner_strip, to_svg, WALL_FILL, stroke="none")

                # Draw U-turns at each opening boundary
                for op in seg_openings:
                    # U-turn at opening start (wall→opening transition)
                    uturn_start = _uturn_polygon(
                        pts, outline_segs, inner_segs, s_segs, g_segs,
                        seg_idx, op.t_start, "start", shell_t, R_in, WALL_OUTER)
                    _svg_polygon(out, uturn_start, to_svg, WALL_FILL,
                                 stroke="none")

                    # U-turn at opening end (opening→wall transition)
                    uturn_end = _uturn_polygon(
                        pts, outline_segs, inner_segs, s_segs, g_segs,
                        seg_idx, op.t_end, "end", shell_t, R_in, WALL_OUTER)
                    _svg_polygon(out, uturn_end, to_svg, WALL_FILL,
                                 stroke="none")

                # Draw opening regions
                for op in seg_openings:
                    F_A, F_B = pts[seg.start], pts[seg.end]
                    W_A = pts[inner_seg.start]
                    W_B = pts[inner_seg.end]
                    o_poly = [
                        _lerp(F_A, F_B, op.t_start),
                        _lerp(F_A, F_B, op.t_end),
                        _lerp(W_A, W_B, op.t_end),
                        _lerp(W_A, W_B, op.t_start),
                    ]
                    _svg_polygon(out, o_poly, to_svg, OPENING_FILL,
                                 stroke="#4682B4", stroke_width="0.5")

    # --- Continuous outlines per wall section ---
    sections = _enumerate_wall_sections(openings, outline_segs)
    for start_op, end_op in sections:
        outer_path, cavity_path = _build_section_outlines(
            pts, outline_segs, inner_segs, s_segs, g_segs,
            start_op, end_op, shell_t, R_in, WALL_OUTER)
        for path in [outer_path, cavity_path]:
            svg_pts = " ".join(
                f"{to_svg(*p)[0]:.2f},{to_svg(*p)[1]:.2f}" for p in path)
            out.append(f'<polygon points="{svg_pts}" fill="none" '
                       f'stroke="#999" stroke-width="0.3"/>')

    # --- Opening labels ---
    for op in openings:
        seg = outline_segs[op.seg_idx]
        inner_seg = inner_segs[op.seg_idx]
        t_mid = (op.t_start + op.t_end) / 2
        F_A, F_B = pts[seg.start], pts[seg.end]
        W_A, W_B = pts[inner_seg.start], pts[inner_seg.end]
        f_mid = _lerp(F_A, F_B, t_mid)
        w_mid = _lerp(W_A, W_B, t_mid)
        cx, cn = (f_mid[0] + w_mid[0]) / 2, (f_mid[1] + w_mid[1]) / 2
        sx, sy = to_svg(cx, cn)
        dE, dN = F_B[0] - F_A[0], F_B[1] - F_A[1]
        svg_angle = -math.degrees(math.atan2(dN, dE))
        if svg_angle > 90:
            svg_angle -= 180
        elif svg_angle < -90:
            svg_angle += 180
        rot = (f' transform="rotate({svg_angle:.1f},{sx:.1f},{sy:.1f})"'
               if abs(svg_angle) > 0.1 else "")
        out.append(f'<text x="{sx:.1f}" y="{sy:.1f}" text-anchor="middle"'
                   f' dominant-baseline="central" font-family="Arial"'
                   f' font-size="5" fill="#4682B4" font-weight="bold"'
                   f'{rot}>{op.name}</text>')

    # --- Title block ---
    out.append(f'<rect x="{data["tb_left"]:.1f}" y="{data["tb_top"]:.1f}"'
               f' width="{data["tb_w"]}" height="{data["tb_h"]}"'
               f' fill="white" stroke="#333" stroke-width="1"/>')
    out.append(f'<text x="{data["tb_cx"]:.1f}" y="{data["tb_top"]+14:.1f}"'
               f' text-anchor="middle" font-family="Arial" font-size="11"'
               f' font-weight="bold" fill="#333">'
               f'{data["outer_area"]:.2f} sq ft</text>')
    out.append(f'<text x="{data["tb_cx"]:.1f}" y="{data["tb_top"]+26:.1f}"'
               f' text-anchor="middle" font-family="Arial" font-size="8"'
               f' fill="#666">Exterior area</text>')

    _ratio = data["ft_per_inch"] * 12
    _scale_label = f'Scale 1:{_ratio:.1f} 1&#8243; = {fmt_dist(data["ft_per_inch"])}'
    out.append(f'<text x="{data["tb_cx"]:.1f}" y="{data["tb_top"]+40:.1f}"'
               f' text-anchor="middle" font-family="Arial" font-size="8"'
               f' fill="#666">{_scale_label}</text>')

    _now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    _git_desc = subprocess.check_output(
        ["git", "describe", "--always"], text=True).strip()
    out.append(f'<text x="{data["tb_cx"]:.1f}" y="{data["tb_top"]+54:.1f}"'
               f' text-anchor="middle" font-family="Arial" font-size="7.5"'
               f' fill="#999">Generated {_now}</text>')
    out.append(f'<text x="{data["tb_cx"]:.1f}" y="{data["tb_top"]+64:.1f}"'
               f' text-anchor="middle" font-family="Arial" font-size="7.5"'
               f' fill="#999">from {_git_desc}</text>')

    # Wall construction note
    out.append(f'<text x="{data["tb_cx"]:.1f}" y="{data["tb_top"]+76:.1f}"'
               f' text-anchor="middle" font-family="Arial" font-size="7"'
               f' fill="#999">2&#8243; shell / 4&#8243; gap / 2&#8243; shell</text>')

    out.append('</svg>')
    return "\n".join(out)


def _partial_line_strip_2(pts, g_seg, inner_seg, t_start, t_end):
    """Build inner shell strip for a partial line segment range."""
    G_A = pts[g_seg.start]
    G_B = pts[g_seg.end]
    W_A = pts[inner_seg.start]
    W_B = pts[inner_seg.end]
    p1 = _lerp(G_A, G_B, t_start)
    p2 = _lerp(G_A, G_B, t_end)
    p3 = _lerp(W_A, W_B, t_end)
    p4 = _lerp(W_A, W_B, t_start)
    return [p1, p2, p3, p4]


# ============================================================
# Main entry point
# ============================================================

if __name__ == "__main__":
    data = build_wall_data()
    svg_content = render_walls_svg(data)

    svg_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "walls.svg")
    with open(svg_path, "w") as f:
        f.write(svg_content)

    print(f"Wall detail written to {svg_path}")
    print(f"Shell: {SHELL_THICKNESS * 12:.0f}\" / "
          f"Gap: {AIR_GAP * 12:.0f}\" / "
          f"Shell: {SHELL_THICKNESS * 12:.0f}\"")
    print(f"Opening corner inside radius: {OPENING_INSIDE_RADIUS * 12:.1f}\"")
