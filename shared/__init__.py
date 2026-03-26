"""Shared types, geometry, survey computation, and SVG utilities."""

from .types import Point, LineSeg, ArcSeg, Segment
from .geometry import (
    GeometryError, GEOM_EPS, ARC_POLY_PTS, require_pts,
    left_norm, off_pt, line_isect, arc_poly,
    circle_circle_isect, line_circle_isect_min_t_gt, line_circle_isect_min_abs_t,
    poly_area, segment_polyline, path_polygon, arc_sweep_deg,
    brg_dist, fmt_brg, fmt_dist,
    horiz_isects, vert_isects,
    compute_inner_walls,
)
from .survey import compute_traverse, compute_three_arc, InsetResult, compute_inset, compute_pt1
from .svg import make_svg_transform, W, H
