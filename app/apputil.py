"""Shared utilities for the app layer.

Small helpers used by both engine.py and variants.py to avoid
duplicating JSON serialisation logic within app/.
"""

from shared.types import LineSeg

# Arc discretisation constants — single source of truth for polygon
# approximations of arcs/circles within the app layer.
ARC_N_SEMICIRCLE = 32   # segments for semicircular arcs (bath sink bulge)
ARC_N_CIRCLE = 24       # segments for full circles (water heater, ET)


def point_to_list(pt):
    """Convert (e, n) tuple to [e, n] list for JSON."""
    return [round(pt[0], 6), round(pt[1], 6)]


def bbox_from_poly(poly):
    """Bounding box dict from polygon list of (e, n) tuples."""
    es = [p[0] for p in poly]
    ns = [p[1] for p in poly]
    return {"w": min(es), "s": min(ns), "e": max(es), "n": max(ns)}


def seg_to_dict(seg):
    """Convert LineSeg/ArcSeg to JSON-serialisable dict."""
    if isinstance(seg, LineSeg):
        return {"type": "line", "start": seg.start, "end": seg.end}
    return {
        "type": "arc", "start": seg.start, "end": seg.end,
        "center": seg.center, "radius": seg.radius,
        "direction": seg.direction, "n_pts": seg.n_pts,
    }
