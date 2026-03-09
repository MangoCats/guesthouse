"""Shared helpers for SCAD generators."""
import math
from shared.types import ArcSeg
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
