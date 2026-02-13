"""Tests for shared/geometry.py pure functions."""
import math
import pytest
from shared.geometry import (
    GeometryError,
    left_norm, off_pt, line_isect, arc_poly,
    circle_circle_isect, line_circle_isect_min_t_gt, line_circle_isect_min_abs_t,
    poly_area, segment_polyline, path_polygon, arc_sweep_deg,
    brg_dist, fmt_brg, fmt_dist,
    horiz_isects, vert_isects,
)
from shared.types import LineSeg, ArcSeg


# --- left_norm ---

def test_left_norm_horizontal():
    n = left_norm((0, 0), (1, 0))
    assert abs(n[0] - 0.0) < 1e-12
    assert abs(n[1] - 1.0) < 1e-12


def test_left_norm_vertical():
    n = left_norm((0, 0), (0, 1))
    assert abs(n[0] - (-1.0)) < 1e-12
    assert abs(n[1] - 0.0) < 1e-12


# --- off_pt ---

def test_off_pt():
    p = off_pt((3, 4), (0, 1), 2.0)
    assert abs(p[0] - 3.0) < 1e-12
    assert abs(p[1] - 6.0) < 1e-12


# --- line_isect ---

def test_line_isect_perpendicular():
    # Horizontal line y=1 and vertical line x=2
    p = line_isect((0, 1), (1, 0), (2, 0), (0, 1))
    assert abs(p[0] - 2.0) < 1e-10
    assert abs(p[1] - 1.0) < 1e-10


def test_line_isect_parallel_raises():
    with pytest.raises(GeometryError, match="Parallel"):
        line_isect((0, 0), (1, 0), (0, 1), (1, 0))


# --- circle_circle_isect ---

def test_circle_circle_isect_unit_circles():
    # Two unit circles at (0,0) and (1,0), near point at (0.5, 1)
    p = circle_circle_isect((0, 0), 1.0, (1, 0), 1.0, (0.5, 1))
    assert abs(p[0] - 0.5) < 1e-10
    assert abs(p[1] - math.sqrt(3) / 2) < 1e-10


def test_circle_circle_isect_too_far():
    with pytest.raises(GeometryError, match="too far"):
        circle_circle_isect((0, 0), 1.0, (10, 0), 1.0, (5, 0))


def test_circle_circle_isect_contained():
    with pytest.raises(GeometryError, match="contained"):
        circle_circle_isect((0, 0), 5.0, (0.1, 0), 1.0, (0, 0))


# --- line_circle_isect_min_t_gt ---

def test_line_circle_isect_min_t_gt_basic():
    # Horizontal line y=0, circle at (3,0) r=1, t_min=0 → hit at x=2
    p = line_circle_isect_min_t_gt((0, 0), (1, 0), (3, 0), 1.0, 0)
    assert abs(p[0] - 2.0) < 1e-10
    assert abs(p[1] - 0.0) < 1e-10


def test_line_circle_isect_miss():
    with pytest.raises(GeometryError, match="misses circle"):
        line_circle_isect_min_t_gt((0, 0), (1, 0), (0, 10), 1.0, 0)


def test_line_circle_isect_no_candidate():
    # Line hits circle at t=-2 and t=-4, but t_min=0
    with pytest.raises(GeometryError, match="No intersection with t >"):
        line_circle_isect_min_t_gt((0, 0), (1, 0), (-3, 0), 1.0, 0)


# --- line_circle_isect_min_abs_t ---

def test_line_circle_isect_min_abs_t_basic():
    # Horizontal line y=0, circle at (5,0) r=1 → closest hit at x=4
    p = line_circle_isect_min_abs_t((0, 0), (1, 0), (5, 0), 1.0)
    assert abs(p[0] - 4.0) < 1e-10


# --- poly_area ---

def test_poly_area_unit_square():
    sq = [(0, 0), (1, 0), (1, 1), (0, 1)]
    assert abs(poly_area(sq) - 1.0) < 1e-12


def test_poly_area_triangle():
    tri = [(0, 0), (4, 0), (0, 3)]
    assert abs(poly_area(tri) - 6.0) < 1e-12


# --- brg_dist ---

def test_brg_dist_north():
    b, d = brg_dist((0, 0), (0, 1))
    assert abs(b - 0.0) < 1e-10
    assert abs(d - 1.0) < 1e-10


def test_brg_dist_east():
    b, d = brg_dist((0, 0), (1, 0))
    assert abs(b - 90.0) < 1e-10
    assert abs(d - 1.0) < 1e-10


def test_brg_dist_southwest():
    b, d = brg_dist((0, 0), (-1, -1))
    assert abs(b - 225.0) < 1e-10
    assert abs(d - math.sqrt(2)) < 1e-10


# --- fmt_brg / fmt_dist ---

def test_fmt_brg_zero():
    assert fmt_brg(0.0) == "0° 00' 00.0\""


def test_fmt_dist_one_foot():
    assert fmt_dist(1.0) == "1' 0\""


def test_fmt_dist_fractional():
    result = fmt_dist(2.5)  # 2'6"
    assert "2'" in result
    assert "6" in result


# --- horiz_isects / vert_isects ---

_SQUARE = [(0, 0), (10, 0), (10, 10), (0, 10)]


def test_horiz_isects_midline():
    xs = sorted(horiz_isects(_SQUARE, 5))
    assert len(xs) == 2
    assert abs(xs[0] - 0.0) < 1e-10
    assert abs(xs[1] - 10.0) < 1e-10


def test_vert_isects_midline():
    ns = sorted(vert_isects(_SQUARE, 5))
    assert len(ns) == 2
    assert abs(ns[0] - 0.0) < 1e-10
    assert abs(ns[1] - 10.0) < 1e-10


# --- segment_polyline / path_polygon / arc_sweep_deg ---

def test_segment_polyline_line():
    pts = {"A": (0, 0), "B": (1, 1)}
    poly = segment_polyline(LineSeg("A", "B"), pts)
    assert len(poly) == 2
    assert poly[0] == (0, 0)
    assert poly[1] == (1, 1)


def test_segment_polyline_arc():
    pts = {"A": (1, 0), "B": (0, 1), "C": (0, 0)}
    seg = ArcSeg("A", "B", "C", 1.0, "CCW", 20)
    poly = segment_polyline(seg, pts)
    assert len(poly) == 21
    # First point should be near A, last near B
    assert abs(poly[0][0] - 1.0) < 1e-10
    assert abs(poly[-1][1] - 1.0) < 1e-10


def test_arc_sweep_deg_quarter():
    pts = {"A": (1, 0), "B": (0, 1), "C": (0, 0)}
    seg = ArcSeg("A", "B", "C", 1.0, "CCW", 20)
    assert abs(arc_sweep_deg(seg, pts) - 90.0) < 1e-8
