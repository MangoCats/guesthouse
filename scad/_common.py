"""Shared helpers for SCAD generators."""
import math
from shared.types import ArcSeg


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
