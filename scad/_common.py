"""Shared helpers for SCAD generators."""
import math
from shared.types import ArcSeg, LineSeg
from shared.wall_shells import lerp


# ── element types ────────────────────────────────────────────
# ("line", x1, y1, x2, y2)             — line segment
# ("arc", cx, cy, r, a1_deg, a2_deg)   — circular arc


def seg_to_elem(seg, pts):
    """Convert a LineSeg or ArcSeg to a T-path element."""
    if isinstance(seg, ArcSeg):
        c = pts[seg.center]
        a1 = math.degrees(math.atan2(pts[seg.start][1] - c[1],
                                      pts[seg.start][0] - c[0]))
        a2 = math.degrees(math.atan2(pts[seg.end][1] - c[1],
                                      pts[seg.end][0] - c[0]))
        if seg.direction == "CW":
            sweep = (a1 - a2) % 360
            a2 = a1 - sweep
        else:
            sweep = (a2 - a1) % 360
            a2 = a1 + sweep
        return ("arc", c[0], c[1], seg.radius, a1, a2)
    return ("line", pts[seg.start][0], pts[seg.start][1],
                    pts[seg.end][0], pts[seg.end][1])


def f8f9_elems(pts, inset, R_turn):
    """F8-F9 corner override: straight south, arc 180→270, straight east."""
    F8, C8 = pts["F8"], pts["C8"]
    R_a8 = C8[0] - F8[0]
    se, sn = F8[0] - inset, F8[1]
    ee, en_ = C8[0], F8[1] - R_a8 - inset
    d = R_a8 + inset - R_turn
    acx, acy = se + R_turn, sn - d
    return [
        ("line", se, sn, se, acy),
        ("arc", acx, acy, R_turn, 180, 270),
        ("line", acx, acy - R_turn, ee, en_),
    ]


# ── boundary tracing ─────────────────────────────────────────

def trace_elems(pts, segs, start_seg_idx, start_t, end_seg_idx, end_t,
                R_out, seg_overrides=None, full_wrap=False):
    """Trace T-path elements between two opening boundaries."""
    n = len(segs)
    elems = []

    def slen(seg):
        a, b = pts[seg.start], pts[seg.end]
        return math.sqrt((b[0] - a[0]) ** 2 + (b[1] - a[1]) ** 2)

    if start_seg_idx == end_seg_idx and not full_wrap:
        seg = segs[start_seg_idx]
        a, b = pts[seg.start], pts[seg.end]
        dt = R_out / slen(seg)
        p1 = lerp(a, b, start_t + dt)
        p2 = lerp(a, b, end_t - dt)
        elems.append(("line", p1[0], p1[1], p2[0], p2[1]))
        return elems

    # Start segment: partial line
    seg = segs[start_seg_idx]
    a, b = pts[seg.start], pts[seg.end]
    dt = R_out / slen(seg)
    p_start = lerp(a, b, start_t + dt)
    elems.append(("line", p_start[0], p_start[1], b[0], b[1]))

    # Intermediate full segments
    idx = (start_seg_idx + 1) % n
    while idx != end_seg_idx:
        seg = segs[idx]
        if seg_overrides and idx in seg_overrides:
            elems.extend(seg_overrides[idx])
        else:
            elems.append(seg_to_elem(seg, pts))
        idx = (idx + 1) % n

    # End segment: partial line
    seg = segs[end_seg_idx]
    a, b = pts[seg.start], pts[seg.end]
    dt = R_out / slen(seg)
    p_end = lerp(a, b, end_t - dt)
    elems.append(("line", a[0], a[1], p_end[0], p_end[1]))

    return elems


def rev_elem(e):
    """Reverse a single element."""
    if e[0] == "line":
        return ("line", e[3], e[4], e[1], e[2])
    _, cx, cy, r, a1, a2 = e
    return ("arc", cx, cy, r, a2, a1)


def rev_elems(elems):
    """Reverse element list and each element's direction."""
    return [rev_elem(e) for e in reversed(elems)]


def qarc_elem(cx, cy, R, u0, u1):
    """Quarter-circle arc from direction u0 toward u1.

    Returns (element, start_pt, end_pt).
    """
    a1 = math.degrees(math.atan2(u0[1], u0[0]))
    a2 = math.degrees(math.atan2(u1[1], u1[0]))
    cross = u0[0] * u1[1] - u0[1] * u1[0]
    if cross > 0:          # CCW
        while a2 < a1:
            a2 += 360
        if a2 - a1 > 180:
            a2 -= 360
    elif cross < 0:        # CW
        while a2 > a1:
            a2 -= 360
        if a1 - a2 > 180:
            a2 += 360
    sp = (cx + R * u0[0], cy + R * u0[1])
    ep = (cx + R * u1[0], cy + R * u1[1])
    return ("arc", cx, cy, R, a1, a2), sp, ep


# ── U-turn and section assembly ───────────────────────────────

def uturn_center_elems(pts, outline_segs, inner_segs, seg_idx, t_param,
                       side, shell_t, R_in):
    """Compute centerline U-turn elements: F-arc, crossover, W-arc.

    Returns list of 3 elements (F-side arc, crossover line, W-side arc).
    These go from F-side to W-side. Reverse for W→F direction.
    """
    R_out = R_in + shell_t
    R_center = R_in + shell_t / 2.0

    seg = outline_segs[seg_idx]
    F_pt = lerp(pts[seg.start], pts[seg.end], t_param)
    iseg = inner_segs[seg_idx]
    W_pt = lerp(pts[iseg.start], pts[iseg.end], t_param)

    dx = pts[seg.end][0] - pts[seg.start][0]
    dy = pts[seg.end][1] - pts[seg.start][1]
    ln = math.sqrt(dx * dx + dy * dy)
    t_hat = (dx / ln, dy / ln)
    n_ext = (-t_hat[1], t_hat[0])

    if side == "start":
        open_dir = t_hat
    else:
        open_dir = (-t_hat[0], -t_hat[1])
    wall_dir = (-open_dir[0], -open_dir[1])
    n_int = (-n_ext[0], -n_ext[1])

    oc = (F_pt[0] - R_out * n_ext[0] + R_out * wall_dir[0],
          F_pt[1] - R_out * n_ext[1] + R_out * wall_dir[1])
    ic = (W_pt[0] + R_out * n_ext[0] + R_out * wall_dir[0],
          W_pt[1] + R_out * n_ext[1] + R_out * wall_dir[1])

    f_arc, _, f_ep = qarc_elem(oc[0], oc[1], R_center, n_ext, open_dir)
    w_arc, w_sp, _ = qarc_elem(ic[0], ic[1], R_center, open_dir, n_int)

    crossover = ("line", f_ep[0], f_ep[1], w_sp[0], w_sp[1])

    return [f_arc, crossover, w_arc]


def build_tpath(pts, outline_segs, inner_segs, tf_segs, tw_segs,
                start_op, end_op, shell_t, R_in, tw_ov, full_wrap=False):
    """Build T-path element list for one wall section."""
    R_out = R_in + shell_t

    # F-side run: trace TF segments between openings (building CW)
    f_elems = trace_elems(pts, tf_segs, start_op.seg_idx, start_op.t_end,
                          end_op.seg_idx, end_op.t_start, R_out,
                          full_wrap=full_wrap)

    # End U-turn (at end_op, side="start": opening ahead in CW direction)
    end_uturn = uturn_center_elems(pts, outline_segs, inner_segs,
                                   end_op.seg_idx, end_op.t_start,
                                   "start", shell_t, R_in)

    # W-side run: trace TW segments, then reverse for T-path CW direction
    w_elems = trace_elems(pts, tw_segs, start_op.seg_idx, start_op.t_end,
                          end_op.seg_idx, end_op.t_start, R_out,
                          seg_overrides=tw_ov, full_wrap=full_wrap)
    w_rev = rev_elems(w_elems)

    # Start U-turn (at start_op, side="end": opening behind)
    # Computed F→W then reversed to get W→F direction
    start_uturn_base = uturn_center_elems(pts, outline_segs, inner_segs,
                                          start_op.seg_idx, start_op.t_end,
                                          "end", shell_t, R_in)
    start_uturn = rev_elems(start_uturn_base)

    return f_elems + end_uturn + w_rev + start_uturn


# ── SCAD output formatting ────────────────────────────────────

def fmt_ft_in(ft, in_width=8):
    """Format feet as ft' inches\" with 4 decimal places on inches."""
    total_in = ft * 12
    whole_ft = int(total_in // 12)
    remaining_in = total_in - whole_ft * 12
    return f"{whole_ft:2d}' {remaining_in:{in_width}.4f}\""


def scad_seg(elem):
    """Format a T-path element as SCAD array literal."""
    if elem[0] == "line":
        _, x1, y1, x2, y2 = elem
        return f"[0, {x1:.8f}, {y1:.8f}, {x2:.8f}, {y2:.8f}]"
    _, cx, cy, r, a1, a2 = elem
    return f"[1, {cx:.8f}, {cy:.8f}, {r:.8f}, {a1:.6f}, {a2:.6f}]"


def seg_comment(elem):
    """Generate an inline comment for a T-path element."""
    if elem[0] == "line":
        _, x1, y1, x2, y2 = elem
        dE, dN = x2 - x1, y2 - y1
        length = math.sqrt(dE * dE + dN * dN)
        bearing = math.degrees(math.atan2(dE, dN)) % 360
        return f"// {fmt_ft_in(length)} @ {bearing:8.4f}deg"
    _, cx, cy, r, a1, a2 = elem
    sweep = a2 - a1
    direction = "CCW" if sweep > 0 else "CW "
    return f"// {direction} {abs(sweep):8.4f}deg R {fmt_ft_in(r, 7)}"


def window_panel_poly(op, outline_segs, inner_segs, pts, panel_half):
    """Compute a 1\"-panel polygon for a window/casement opening.

    Returns list of 4 (x, y) vertices at the shell centerline face.
    """
    seg = outline_segs[op.seg_idx]
    iseg = inner_segs[op.seg_idx]
    F_A, F_B = pts[seg.start], pts[seg.end]
    W_A, W_B = pts[iseg.start], pts[iseg.end]
    M_A = ((F_A[0] + W_A[0]) / 2, (F_A[1] + W_A[1]) / 2)
    M_B = ((F_B[0] + W_B[0]) / 2, (F_B[1] + W_B[1]) / 2)
    M_s = lerp(M_A, M_B, op.t_start)
    M_e = lerp(M_A, M_B, op.t_end)
    dx, dy = F_B[0] - F_A[0], F_B[1] - F_A[1]
    ln = math.sqrt(dx * dx + dy * dy)
    nx, ny = -dy / ln, dx / ln  # exterior normal
    return [
        (M_s[0] + nx * panel_half, M_s[1] + ny * panel_half),
        (M_e[0] + nx * panel_half, M_e[1] + ny * panel_half),
        (M_e[0] - nx * panel_half, M_e[1] - ny * panel_half),
        (M_s[0] - nx * panel_half, M_s[1] - ny * panel_half),
    ]


def _seg_start_bearing_deg(seg, pts):
    """Compute compass bearing (degrees, 0=N/90=E) at the start of an inner seg."""
    if isinstance(seg, ArcSeg):
        c = pts[seg.center]
        s = pts[seg.start]
        rx, ry = s[0] - c[0], s[1] - c[1]
        if seg.direction == "CW":
            dx, dy = ry, -rx  # CW tangent = right normal of radius
        else:
            dx, dy = -ry, rx  # CCW tangent = left normal
    else:
        s, e = pts[seg.start], pts[seg.end]
        dx, dy = e[0] - s[0], e[1] - s[1]
    return math.degrees(math.atan2(dx, dy)) % 360


def override_tw_elems(chain, w_start_pt, tw_start_pt, initial_bearing_deg):
    """Compute TW-series SCAD elements from a W-series override chain.

    Walks both W and TW paths simultaneously. For line sub-segments the
    same displacement is applied to both; for arc sub-segments the W arc
    center is used with R_TW = dist(tw_cur, center) so the TW arc is the
    parallel-offset curve of the W arc.

    Returns list of ("line",...) or ("arc",...) SCAD element tuples.
    """
    w_cur = w_start_pt
    tw_cur = tw_start_pt
    cur_bearing = initial_bearing_deg
    elems = []

    for sub in chain:
        seg_type = sub["seg_type"]
        if seg_type == "L":
            brg = sub["bearing"]
            dist = sub["distance"]
            brg_rad = math.radians(brg)
            dx = dist * math.sin(brg_rad)
            dy = dist * math.cos(brg_rad)
            w_next = (w_cur[0] + dx, w_cur[1] + dy)
            tw_next = (tw_cur[0] + dx, tw_cur[1] + dy)
            elems.append(("line", tw_cur[0], tw_cur[1], tw_next[0], tw_next[1]))
            w_cur = w_next
            tw_cur = tw_next
            cur_bearing = brg
        else:
            # Arc: CCW or CW
            radius = sub["radius"]
            sweep_deg = sub["sweep"]
            sweep_rad = math.radians(sweep_deg)

            # Travel direction vector from compass bearing
            dir_x = math.sin(math.radians(cur_bearing))
            dir_y = math.cos(math.radians(cur_bearing))

            # W arc center (left of travel for CCW, right for CW)
            if seg_type == "CCW":
                cx = w_cur[0] - dir_y * radius
                cy = w_cur[1] + dir_x * radius
            else:
                cx = w_cur[0] + dir_y * radius
                cy = w_cur[1] - dir_x * radius

            # TW arc: same center, radius = dist(tw_cur, center)
            r_tw = math.hypot(tw_cur[0] - cx, tw_cur[1] - cy)

            # SCAD arc angles from TW entry point
            a1 = math.degrees(math.atan2(tw_cur[1] - cy, tw_cur[0] - cx))
            if seg_type == "CCW":
                a2 = a1 + sweep_deg
            else:
                a2 = a1 - sweep_deg
            elems.append(("arc", cx, cy, r_tw, a1, a2))

            # Advance W along its arc
            w_entry = math.atan2(w_cur[1] - cy, w_cur[0] - cx)
            w_exit = w_entry + sweep_rad if seg_type == "CCW" else w_entry - sweep_rad
            w_cur = (cx + radius * math.cos(w_exit),
                     cy + radius * math.sin(w_exit))

            # Advance TW along its arc
            tw_entry = math.atan2(tw_cur[1] - cy, tw_cur[0] - cx)
            tw_exit = tw_entry + sweep_rad if seg_type == "CCW" else tw_entry - sweep_rad
            tw_cur = (cx + r_tw * math.cos(tw_exit),
                      cy + r_tw * math.sin(tw_exit))

            # Update compass bearing
            if seg_type == "CCW":
                cur_bearing = (cur_bearing - sweep_deg) % 360
            else:
                cur_bearing = (cur_bearing + sweep_deg) % 360

    return elems


def build_tw_overrides(overrides, inner_segs, pts):
    """Build tw_ov dict from W inner-wall override chains.

    For each override, computes TW SCAD elements that replace the default
    TW path for that segment index.  The TW start point is looked up as
    "TW" + (W-point numeric suffix), which must already be present in pts
    (added by compute_inset_path with prefix="TW").

    Returns dict {seg_idx: [elements]} suitable for merging into tw_ov.
    """
    result = {}
    n = len(inner_segs)
    for seg_idx, chain in overrides.items():
        if seg_idx >= n or not chain:
            continue
        seg = inner_segs[seg_idx]
        if seg.start not in pts:
            continue

        w_start = pts[seg.start]

        # TW point: same suffix as W point (W8 → TW8, W04 → TW04)
        tw_name = "TW" + seg.start[1:]
        if tw_name not in pts:
            continue
        tw_start = pts[tw_name]

        # Initial bearing: from first sub-seg if L, else compute from seg
        if chain[0]["seg_type"] == "L":
            initial_brg = chain[0]["bearing"]
        else:
            initial_brg = _seg_start_bearing_deg(seg, pts)

        elems = override_tw_elems(chain, w_start, tw_start, initial_brg)
        if elems:
            result[seg_idx] = elems

    return result


def compute_wall_bands(openings, outline_segs, upper_top,
                       pts, tf_segs, tw_segs, inner_segs,
                       shell_t, R_in, tw_ov):
    """Compute per-elevation wall bands driven by per-opening bottom/top elevations.

    Returns list of (z_start, z_end, section_data) where section_data is a list
    of (label, tpath) pairs ready for linear_extrude.  Bands with no active
    openings use the full-perimeter T-path (label="full", same tpath each time).

    The full-perimeter T-path is built once and reused for all solid bands.

    openings: list of WallOpening (must have .bottom_elev and .top_elev)
    upper_top: top of the wall structure (no openings above this)
    """
    from shared.wall_shells import enumerate_wall_sections

    # Build full-perimeter T-path (reused for solid bands and upper band)
    if openings:
        # Use first opening as seam anchor for full wrap
        all_secs = enumerate_wall_sections(openings, outline_segs)
        seam_op = (all_secs[-1:] + all_secs[:-1])[0][0]
    else:
        seam_op = None

    def _full_tpath():
        if seam_op is None:
            return []
        return build_tpath(pts, outline_segs, inner_segs, tf_segs, tw_segs,
                           seam_op, seam_op, shell_t, R_in, tw_ov, full_wrap=True)

    _cached_full = None

    def get_full_tpath():
        nonlocal _cached_full
        if _cached_full is None:
            _cached_full = _full_tpath()
        return _cached_full

    # Collect unique elevation breakpoints from openings
    breaks = sorted(set(
        [0.0] +
        [op.bottom_elev for op in openings] +
        [op.top_elev for op in openings]
    ))
    # Clamp breaks below upper_top
    breaks = [b for b in breaks if b < upper_top - 1e-9]

    bands = []
    for i, z1 in enumerate(breaks):
        z2 = breaks[i + 1] if i + 1 < len(breaks) else upper_top
        if z2 - z1 < 1e-9:
            continue
        # Active openings: those whose opening span fully covers this band
        active = [op for op in openings
                  if op.bottom_elev <= z1 + 1e-9 and op.top_elev >= z2 - 1e-9]
        if not active:
            bands.append((z1, z2, [("full", get_full_tpath())]))
            continue
        secs = enumerate_wall_sections(active, outline_segs)
        secs = secs[-1:] + secs[:-1]
        sec_data = []
        for start_op, end_op in secs:
            # full_wrap=True when start and end are the same opening (single
            # active opening in this band): the wall section must traverse the
            # entire perimeter instead of tracing a stub within one segment.
            wrap = (start_op is end_op)
            tp = build_tpath(pts, outline_segs, inner_segs, tf_segs, tw_segs,
                             start_op, end_op, shell_t, R_in, tw_ov,
                             full_wrap=wrap)
            sec_data.append((f"b{i}_{start_op.name}_{end_op.name}", tp))
        bands.append((z1, z2, sec_data))

    return bands


# ── Interior wall SCAD output ─────────────────────────────────

# White pine: pale, warm cream (Eastern White Pine)
IW_PINE_COLOR = "[0.91, 0.85, 0.65]"
IW_DOOR_H = 82.5 / 12.0  # 82.5" door height in feet


def _pts_str(poly):
    """Format polygon point list as SCAD array."""
    return ", ".join(f"[{p[0]:.8f}, {p[1]:.8f}]" for p in poly)


def interior_wall_scad(out, iw_polys, ro_polys, get_wall_top):
    """Append interior wall SCAD blocks to output list.

    Parameters:
        out          — list to append SCAD lines to
        iw_polys     — dict {wall_name: [[x,y], ...]} (4-vertex polygons)
        ro_polys     — list of rough-opening dicts with keys:
                        name, wall_name, poly
        get_wall_top — callable(wall_name) → top elevation in feet

    Each interior wall is extruded from z=0 to get_wall_top(name).
    Rough openings are subtracted from 0 to IW_DOOR_H (82.5").
    """
    if not iw_polys:
        return

    # Group openings by wall
    openings_by_wall = {}
    for ro in (ro_polys or []):
        openings_by_wall.setdefault(ro["wall_name"], []).append(ro)

    out.append("// Interior walls (white pine, openings at 82.5\")")
    out.append(f"iw_door_h = {IW_DOOR_H:.8f};  // 82.5\" door height")
    out.append(f"color({IW_PINE_COLOR}) union() {{")

    for name in sorted(iw_polys):
        poly = iw_polys[name]
        z_top = get_wall_top(name)
        wall_ros = openings_by_wall.get(name, [])

        out.append(f"  // {name}")
        if wall_ros:
            out.append(f"  difference() {{")
            out.append(f"    linear_extrude(height = {z_top:.6f}, convexity = 10)")
            out.append(f"      polygon(points = [{_pts_str(poly)}]);")
            for ro in wall_ros:
                out.append(f"    // {ro['name']}")
                out.append(f"    translate([0, 0, -0.001])")
                out.append(f"    linear_extrude(height = iw_door_h + 0.001, convexity = 10)")
                out.append(f"      polygon(points = [{_pts_str(ro['poly'])}]);")
            out.append(f"  }}")
        else:
            out.append(f"  linear_extrude(height = {z_top:.6f}, convexity = 10)")
            out.append(f"    polygon(points = [{_pts_str(poly)}]);")

    out.append("}")
    out.append("")
