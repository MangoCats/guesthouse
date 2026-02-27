"""Shell geometry functions for double-shell wall construction.

Computes inset paths, opening parametrics, shell strip polygons,
U-turn arcs, wall section enumeration, and section outlines.
Used by both floorplan/gen_floorplan.py and walls/gen_walls.py.
"""
import math

from shared.types import LineSeg, ArcSeg
from shared.geometry import compute_inner_walls, segment_polyline, GEOM_EPS

# Default arc point count for U-turn quarter-circle arcs
UTURN_ARC_PTS = 12

# ============================================================
# Shell path computation
# ============================================================

def compute_inset_path(outline_segs, pts, radii, inset, prefix
                       ) -> tuple[dict, list]:
    """Compute a shell boundary path at given inset distance.

    Returns (new_pts_dict, new_segs) with point names using the given prefix
    (e.g., "S0".."S20" for prefix="S").
    """
    tmp_pts = dict(pts)
    tmp_segs = compute_inner_walls(outline_segs, tmp_pts, inset, radii)

    result_pts = {}
    for seg in outline_segs:
        suffix = seg.end[1:]  # "0", "1", ..., "20a", "21"
        result_pts[f"{prefix}{suffix}"] = tmp_pts[f"W{suffix}"]

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


def openings_on_seg(openings, seg_idx) -> list:
    """Get openings on a given segment, sorted by t_start."""
    result = [o for o in openings if o.seg_idx == seg_idx]
    result.sort(key=lambda o: o.t_start)
    return result


def solid_ranges(seg_openings) -> list[tuple[float, float]]:
    """Compute solid wall parametric ranges from sorted opening list.

    Returns list of (t_start, t_end) for solid wall sections.
    """
    ranges = []
    cursor = 0.0
    for o in seg_openings:
        if o.t_start > cursor + GEOM_EPS:
            ranges.append((cursor, o.t_start))
        cursor = o.t_end
    if cursor < 1.0 - GEOM_EPS:
        ranges.append((cursor, 1.0))
    return ranges


# ============================================================
# Polygon builders
# ============================================================

def lerp(a, b, t) -> tuple[float, float]:
    """Linear interpolation between two points."""
    return (a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1]))


def arc_strip_poly(seg, pts, outer_prefix, inner_seg):
    """Build polygon for an arc shell strip (full segment).

    Returns list of (E, N) points forming the closed polygon.
    """
    outer_poly = segment_polyline(seg, pts)
    inner_poly = segment_polyline(inner_seg, pts)
    # Forward along outer, backward along inner
    return outer_poly + list(reversed(inner_poly))


def line_strip_poly(pts, seg_start, seg_end, inner_start, inner_end):
    """Build polygon for a line shell strip (full segment or sub-range).

    4-point rectangle: start_outer, end_outer, end_inner, start_inner.
    """
    return [
        pts[seg_start], pts[seg_end],
        pts[inner_end], pts[inner_start],
    ]


def partial_line_strip(pts, seg, inner_seg_start, inner_seg_end, t_start, t_end):
    """Build polygon for a partial line shell strip between t_start and t_end."""
    A_out = pts[seg.start]
    B_out = pts[seg.end]
    A_in = pts[inner_seg_start]
    B_in = pts[inner_seg_end]
    p1 = lerp(A_out, B_out, t_start)
    p2 = lerp(A_out, B_out, t_end)
    p3 = lerp(A_in, B_in, t_end)
    p4 = lerp(A_in, B_in, t_start)
    return [p1, p2, p3, p4]


def partial_line_strip_2(pts, g_seg, inner_seg, t_start, t_end):
    """Build inner shell strip for a partial line segment range."""
    G_A = pts[g_seg.start]
    G_B = pts[g_seg.end]
    W_A = pts[inner_seg.start]
    W_B = pts[inner_seg.end]
    p1 = lerp(G_A, G_B, t_start)
    p2 = lerp(G_A, G_B, t_end)
    p3 = lerp(W_A, W_B, t_end)
    p4 = lerp(W_A, W_B, t_start)
    return [p1, p2, p3, p4]


# ============================================================
# U-turn geometry
# ============================================================

def uturn_arc_data(pts, outline_segs, inner_segs, seg_idx, t_param, side,
                   shell_t, R_in, wall_total, n_arc=UTURN_ARC_PTS
                   ) -> dict[str, list[tuple[float, float]]]:
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
    F_pt = lerp(F_A, F_B, t_param)

    inner_seg = inner_segs[seg_idx]
    W_A, W_B = pts[inner_seg.start], pts[inner_seg.end]
    W_pt = lerp(W_A, W_B, t_param)

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


def uturn_polygon(pts, outline_segs, inner_segs, s_segs, g_segs,
                  seg_idx, t_param, side, shell_t, R_in, wall_total,
                  n_arc=UTURN_ARC_PTS) -> list[tuple[float, float]]:
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
    arcs = uturn_arc_data(pts, outline_segs, inner_segs, seg_idx, t_param,
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

def trace_boundary_path(pts, segs, start_seg_idx, start_t, end_seg_idx,
                        end_t, R_out, seg_overrides=None
                        ) -> list[tuple[float, float]]:
    """Trace a boundary path between two opening boundaries across segments.

    Traces CW from (start_seg_idx, start_t + delta_t) to
    (end_seg_idx, end_t - delta_t), spanning intermediate segments.

    segs: one of outline_segs/s_segs/g_segs/inner_segs (all 22 segments).
    start_t: parametric position of the starting opening's t_end.
    end_t: parametric position of the ending opening's t_start.
    R_out: trim distance in feet (R_in + shell_t) — converted to delta_t
           using each line segment's length.
    seg_overrides: optional dict mapping seg index to replacement polyline
                   (list of (E, N) points) for non-standard segment paths.

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
        path.append(lerp(A, B, start_t + dt))
        path.append(lerp(A, B, end_t - dt))
        return path

    # Multi-segment: partial start + full intermediates + partial end

    # Start segment (from start_t + delta to segment end)
    seg = segs[start_seg_idx]
    A, B = pts[seg.start], pts[seg.end]
    seg_len = math.sqrt((B[0] - A[0])**2 + (B[1] - A[1])**2)
    dt = R_out / seg_len
    path.append(lerp(A, B, start_t + dt))
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
        if seg_overrides and idx in seg_overrides:
            path.extend(seg_overrides[idx][1:])  # skip first (matches prev end)
        elif isinstance(seg, ArcSeg):
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
    path.append(lerp(A, B, end_t - dt))

    return path


def enumerate_wall_sections(openings, outline_segs) -> list[tuple]:
    """Enumerate wall sections as (start_opening, end_opening) pairs.

    Each wall section is bounded by start_opening.t_end on one side and
    end_opening.t_start on the other, going CW around the building.
    Returns list of (start_op, end_op) tuples.
    """
    # Collect all opening boundaries in CW order
    ordered = []
    for seg_idx in range(len(outline_segs)):
        seg_ops = openings_on_seg(openings, seg_idx)
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


def build_section_outlines(pts, outline_segs, inner_segs, s_segs, g_segs,
                           start_op, end_op, shell_t, R_in, wall_total,
                           n_arc=UTURN_ARC_PTS, g_seg_overrides=None,
                           w_seg_overrides=None
                           ) -> tuple[list[tuple[float, float]], list[tuple[float, float]]]:
    """Build outer and inner cavity outlines for one wall section.

    start_op: opening whose t_end starts this wall section
    end_op: opening whose t_start ends this wall section

    Returns (outer_path, cavity_path) as lists of (E, N) points.
    """
    R_out = R_in + shell_t

    # U-turn arc data at each end
    start_arcs = uturn_arc_data(pts, outline_segs, inner_segs,
                                start_op.seg_idx, start_op.t_end, "end",
                                shell_t, R_in, wall_total, n_arc)
    end_arcs = uturn_arc_data(pts, outline_segs, inner_segs,
                              end_op.seg_idx, end_op.t_start, "start",
                              shell_t, R_in, wall_total, n_arc)

    # Trace boundary paths between the two openings
    f_path = trace_boundary_path(pts, outline_segs,
                                 start_op.seg_idx, start_op.t_end,
                                 end_op.seg_idx, end_op.t_start, R_out)
    s_path = trace_boundary_path(pts, s_segs,
                                 start_op.seg_idx, start_op.t_end,
                                 end_op.seg_idx, end_op.t_start, R_out)
    g_path = trace_boundary_path(pts, g_segs,
                                 start_op.seg_idx, start_op.t_end,
                                 end_op.seg_idx, end_op.t_start, R_out,
                                 seg_overrides=g_seg_overrides)
    w_path = trace_boundary_path(pts, inner_segs,
                                 start_op.seg_idx, start_op.t_end,
                                 end_op.seg_idx, end_op.t_start, R_out,
                                 seg_overrides=w_seg_overrides)

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
