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

# Legacy default cap corner inside radius (10mm), used when no per-type radii
# are configured.
_DEFAULT_R_IN = 10.0 / 304.8


def opening_radii_config(constants) -> dict:
    """Build the per-type cap corner radii config from a constants dict.

    Returns {door_outer, door_inner, window_outer, window_inner} of cavity-side
    (inside) radii.  Each falls back to OPENING_INSIDE_RADIUS (legacy single
    value) then to the 10mm default.
    """
    base = (constants or {}).get("OPENING_INSIDE_RADIUS", _DEFAULT_R_IN)
    g = (constants or {}).get
    return {
        "door_outer":   g("OPENING_DOOR_OUTER_R", base),
        "door_inner":   g("OPENING_DOOR_INNER_R", base),
        "window_outer": g("OPENING_WINDOW_OUTER_R", base),
        "window_inner": g("OPENING_WINDOW_INNER_R", base),
    }


def opening_corner_radii(op, cfg) -> tuple[float, float]:
    """Return (R_in_outer_shell, R_in_inner_shell) for an opening.

    Doors use the door_* radii; every other opening type (window, casement,
    casement_r) uses the window_* radii.
    """
    if getattr(op, "opening_type", "window") == "door":
        return cfg["door_outer"], cfg["door_inner"]
    return cfg["window_outer"], cfg["window_inner"]

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


def solid_ranges_bordered(seg_openings) -> list[tuple]:
    """Like solid_ranges, but records the opening bordering each end.

    Returns list of (t_start, t_end, op_before, op_after) where op_before is the
    opening whose t_end == t_start (None at the segment start) and op_after is
    the opening whose t_start == t_end (None at the segment end).  Used to trim
    each solid shell strip by the correct per-opening cap radius at each end.
    """
    ranges = []
    cursor = 0.0
    op_before = None
    for o in seg_openings:
        if o.t_start > cursor + GEOM_EPS:
            ranges.append((cursor, o.t_start, op_before, o))
        cursor = o.t_end
        op_before = o
    if cursor < 1.0 - GEOM_EPS:
        ranges.append((cursor, 1.0, op_before, None))
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


def partial_line_strip(pts, outer_seg, inner_seg, t_start, t_end):
    """Build polygon for a partial line shell strip between t_start and t_end.

    outer_seg and inner_seg are LineSeg/ArcSeg objects with .start/.end names.
    """
    A_out = pts[outer_seg.start]
    B_out = pts[outer_seg.end]
    A_in = pts[inner_seg.start]
    B_in = pts[inner_seg.end]
    p1 = lerp(A_out, B_out, t_start)
    p2 = lerp(A_out, B_out, t_end)
    p3 = lerp(A_in, B_in, t_end)
    p4 = lerp(A_in, B_in, t_start)
    return [p1, p2, p3, p4]


# ============================================================
# U-turn geometry
# ============================================================

def uturn_arc_data(pts, outline_segs, inner_segs, seg_idx, t_param, side,
                   shell_t, R_in_oc, R_in_ic, wall_total, n_arc=UTURN_ARC_PTS
                   ) -> dict[str, list[tuple[float, float]]]:
    """Compute U-turn arc point arrays at an opening boundary.

    Returns a dict with keys:
      'oc_F': F-face arc points (outer shell, outer face)
      'oc_S': S-face arc points (outer shell, inner face)
      'ic_W': W-face arc points (inner shell, inner face)
      'ic_G': G-face arc points (inner shell, outer face)

    ``R_in_oc``/``R_in_ic`` are the cavity-side cap radii of the outer-shell and
    inner-shell corners respectively (independently adjustable per opening type);
    the corresponding shell-face radii are ``R_in_* + shell_t``.

    Each arc goes from the shell face toward the cross-wall face:
      oc_F[0]/oc_S[0] = on F/S-face, R_out_oc back from opening boundary
      oc_F[-1]/oc_S[-1] = on cross-wall face
      ic_W[0]/ic_G[0] = on cross-wall face
      ic_W[-1]/ic_G[-1] = on W/G-face, R_out_ic back from opening boundary

    side: "start" means wall-to-opening transition (wall at t < t_param)
          "end" means opening-to-wall transition (wall at t > t_param)
    """
    R_out_oc = R_in_oc + shell_t
    R_out_ic = R_in_ic + shell_t
    seg = outline_segs[seg_idx]

    # Points at the boundary parameter on outer and inner faces
    F_A, F_B = pts[seg.start], pts[seg.end]
    F_pt = lerp(F_A, F_B, t_param)

    inner_seg = inner_segs[seg_idx]
    W_A, W_B = pts[inner_seg.start], pts[inner_seg.end]
    W_pt = lerp(W_A, W_B, t_param)

    # Tangent direction (along wall, CW traversal)
    dx, dy = F_B[0] - F_A[0], F_B[1] - F_A[1]
    t_len = math.hypot(dx, dy)
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
    # Outer shell: R_out_oc inward from F-face, R_out_oc back from opening
    oc = (F_pt[0] - R_out_oc * n_ext[0] + R_out_oc * wall_dir[0],
          F_pt[1] - R_out_oc * n_ext[1] + R_out_oc * wall_dir[1])
    # Inner shell: R_out_ic outward from W-face, R_out_ic back from opening
    ic = (W_pt[0] + R_out_ic * n_ext[0] + R_out_ic * wall_dir[0],
          W_pt[1] + R_out_ic * n_ext[1] + R_out_ic * wall_dir[1])

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
    oc_F = qarc(oc[0], oc[1], R_out_oc, n_ext, open_dir)  # F-face arc
    oc_S = qarc(oc[0], oc[1], R_in_oc, n_ext, open_dir)   # S-face arc

    # Inner shell arcs: from cross-wall (open_dir) to shell face (-n_ext)
    n_int = (-n_ext[0], -n_ext[1])
    ic_W = qarc(ic[0], ic[1], R_out_ic, open_dir, n_int)   # W-face arc
    ic_G = qarc(ic[0], ic[1], R_in_ic, open_dir, n_int)    # G-face arc

    return {'oc_F': oc_F, 'oc_S': oc_S, 'ic_W': ic_W, 'ic_G': ic_G}


def uturn_polygon(pts, outline_segs, inner_segs, s_segs, g_segs,
                  seg_idx, t_param, side, shell_t, R_in_oc, R_in_ic, wall_total,
                  n_arc=UTURN_ARC_PTS) -> list[tuple[float, float]]:
    """Build the U-turn polygon at an opening boundary.

    The turn connects the outer shell to the inner shell via two 90-degree
    arcs and a straight cross-wall section, curving toward building interior.
    The outer-shell and inner-shell caps have independent cavity-side radii
    ``R_in_oc``/``R_in_ic`` (shell-face radii ``+ shell_t``).

        F-face --,              R_out_oc = R_in_oc + shell_t
                 |
        S-face -,|              R_in_oc
                ||
                || straight (wall_total - R_out_oc - R_out_ic)
                ||
        G-face -'|              R_in_ic
                 |
        W-face --'              R_out_ic = R_in_ic + shell_t

    side: "start" means wall-to-opening transition (wall at t < t_param)
          "end" means opening-to-wall transition (wall at t > t_param)

    Returns a list of (E, N) points forming the U-turn polygon.
    """
    arcs = uturn_arc_data(pts, outline_segs, inner_segs, seg_idx, t_param,
                          side, shell_t, R_in_oc, R_in_ic, wall_total, n_arc)

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
                        end_t, R_out_start, R_out_end=None, seg_overrides=None
                        ) -> list[tuple[float, float]]:
    """Trace a boundary path between two opening boundaries across segments.

    Traces CW from (start_seg_idx, start_t + delta_start) to
    (end_seg_idx, end_t - delta_end), spanning intermediate segments.

    segs: one of outline_segs/s_segs/g_segs/inner_segs (all 22 segments).
    start_t: parametric position of the starting opening's t_end.
    end_t: parametric position of the ending opening's t_start.
    R_out_start / R_out_end: trim distances in feet at the start / end ends
           (the bordering openings' cap outside radii for this shell face) —
           converted to delta_t using each line segment's length.  R_out_end
           defaults to R_out_start when omitted.
    seg_overrides: optional dict mapping seg index to replacement polyline
                   (list of (E, N) points) for non-standard segment paths.

    Returns list of (E, N) points along the boundary.
    """
    if R_out_end is None:
        R_out_end = R_out_start
    n_segs = len(segs)
    path = []

    if start_seg_idx == end_seg_idx:
        # Same segment — just two interpolated points
        seg = segs[start_seg_idx]
        A, B = pts[seg.start], pts[seg.end]
        seg_len = math.hypot(B[0] - A[0], B[1] - A[1])
        path.append(lerp(A, B, start_t + R_out_start / seg_len))
        path.append(lerp(A, B, end_t - R_out_end / seg_len))
        return path

    # Multi-segment: partial start + full intermediates + partial end

    # Start segment (from start_t + delta to segment end)
    seg = segs[start_seg_idx]
    A, B = pts[seg.start], pts[seg.end]
    seg_len = math.hypot(B[0] - A[0], B[1] - A[1])
    dt = R_out_start / seg_len
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
    seg_len = math.hypot(B[0] - A[0], B[1] - A[1])
    dt = R_out_end / seg_len
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
                           start_op, end_op, shell_t, radii_cfg, wall_total,
                           n_arc=UTURN_ARC_PTS, g_seg_overrides=None,
                           w_seg_overrides=None
                           ) -> tuple[list[tuple[float, float]], list[tuple[float, float]]]:
    """Build outer and inner cavity outlines for one wall section.

    start_op: opening whose t_end starts this wall section
    end_op: opening whose t_start ends this wall section
    radii_cfg: per-type cap radii config (see opening_radii_config); each U-turn
        uses the radii for its own opening's type, and the outer/inner shell
        boundary paths are trimmed by the matching outer-shell / inner-shell cap
        outside radius at each end.

    Returns (outer_path, cavity_path) as lists of (E, N) points.
    """
    s_oc, s_ic = opening_corner_radii(start_op, radii_cfg)
    e_oc, e_ic = opening_corner_radii(end_op, radii_cfg)
    # Outer-shell cap outside radii (F/S faces); inner-shell cap outside radii (G/W)
    s_Rout_oc, e_Rout_oc = s_oc + shell_t, e_oc + shell_t
    s_Rout_ic, e_Rout_ic = s_ic + shell_t, e_ic + shell_t

    # U-turn arc data at each end (each end uses its own opening's radii)
    start_arcs = uturn_arc_data(pts, outline_segs, inner_segs,
                                start_op.seg_idx, start_op.t_end, "end",
                                shell_t, s_oc, s_ic, wall_total, n_arc)
    end_arcs = uturn_arc_data(pts, outline_segs, inner_segs,
                              end_op.seg_idx, end_op.t_start, "start",
                              shell_t, e_oc, e_ic, wall_total, n_arc)

    # Trace boundary paths between the two openings.  Outer-shell faces (F, S)
    # are trimmed by the outer-shell cap radius; inner-shell faces (G, W) by the
    # inner-shell cap radius — independently at each end.
    f_path = trace_boundary_path(pts, outline_segs,
                                 start_op.seg_idx, start_op.t_end,
                                 end_op.seg_idx, end_op.t_start,
                                 s_Rout_oc, e_Rout_oc)
    s_path = trace_boundary_path(pts, s_segs,
                                 start_op.seg_idx, start_op.t_end,
                                 end_op.seg_idx, end_op.t_start,
                                 s_Rout_oc, e_Rout_oc)
    g_path = trace_boundary_path(pts, g_segs,
                                 start_op.seg_idx, start_op.t_end,
                                 end_op.seg_idx, end_op.t_start,
                                 s_Rout_ic, e_Rout_ic,
                                 seg_overrides=g_seg_overrides)
    w_path = trace_boundary_path(pts, inner_segs,
                                 start_op.seg_idx, start_op.t_end,
                                 end_op.seg_idx, end_op.t_start,
                                 s_Rout_ic, e_Rout_ic,
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
