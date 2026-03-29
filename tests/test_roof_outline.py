"""Tests for shared/roof_outline.py — DbRoofResult and related functions."""
import math
import pytest
from shared.roof_outline import (
    _ext_tangent_outside,
    _sample_arc_cw,
    compute_db_roof_outline,
    db_roof_segments,
    DbRoofResult,
)


# ---------------------------------------------------------------------------
# Shared test geometry: square building, 4 corners, all same radius
# ---------------------------------------------------------------------------

_SQ_PTS = {"C1": (5, 5), "C2": (5, -5), "C3": (-5, -5), "C4": (-5, 5)}
_SQ_RADII = {"R_a1": 1.0, "R_a2": 1.0, "R_a3": 1.0, "R_a4": 1.0}
_SQ_NAMES = ["C1", "C2", "C3", "C4"]
_OVERHANG = 0.5


# ---------------------------------------------------------------------------
# _ext_tangent_outside
# ---------------------------------------------------------------------------

class TestExtTangentOutside:
    def test_equal_radii_gives_parallel_tangent(self):
        """Equal-radius circles on x-axis: tangent points at top, normal=(0,1)."""
        ca, ra = (0.0, 0.0), 5.0
        cb, rb = (10.0, 0.0), 5.0
        pa, pb, nx, ny = _ext_tangent_outside(ca, ra, cb, rb)
        assert pa == pytest.approx((0.0, 5.0), abs=1e-9)
        assert pb == pytest.approx((10.0, 5.0), abs=1e-9)
        assert nx == pytest.approx(0.0, abs=1e-9)
        assert ny == pytest.approx(1.0, abs=1e-9)

    def test_tangent_points_lie_on_their_circles(self):
        """pa distance from ca == ra; pb distance from cb == rb."""
        ca, ra = (0.0, 0.0), 3.0
        cb, rb = (0.0, 12.0), 6.0
        pa, pb, nx, ny = _ext_tangent_outside(ca, ra, cb, rb)
        assert math.hypot(pa[0] - ca[0], pa[1] - ca[1]) == pytest.approx(ra, abs=1e-9)
        assert math.hypot(pb[0] - cb[0], pb[1] - cb[1]) == pytest.approx(rb, abs=1e-9)

    def test_normal_is_unit_vector(self):
        """Normal vector has magnitude 1."""
        _, _, nx, ny = _ext_tangent_outside((0.0, 0.0), 4.0, (15.0, 0.0), 7.0)
        assert math.hypot(nx, ny) == pytest.approx(1.0, abs=1e-9)

    def test_tangent_direction_perpendicular_to_normal(self):
        """Vector pa→pb is perpendicular to the outward normal."""
        ca, ra = (0.0, 0.0), 5.0
        cb, rb = (20.0, 10.0), 8.0
        pa, pb, nx, ny = _ext_tangent_outside(ca, ra, cb, rb)
        dx, dy = pb[0] - pa[0], pb[1] - pa[1]
        d = math.hypot(dx, dy)
        dot = (dx / d) * nx + (dy / d) * ny
        assert dot == pytest.approx(0.0, abs=1e-9)

    def test_coincident_centers_fallback(self):
        """Degenerate case (coincident centers) returns default normal (0, 1)."""
        pa, pb, nx, ny = _ext_tangent_outside((3.0, 3.0), 5.0, (3.0, 3.0), 7.0)
        assert nx == pytest.approx(0.0, abs=1e-9)
        assert ny == pytest.approx(1.0, abs=1e-9)

    def test_cos_t_clamped_when_radii_exceed_distance(self):
        """cos_t = (ra - rb) / d is clamped to [-1, 1]; no ValueError raised."""
        # ra - rb > d would give |cos_t| > 1 without clamping
        ca, ra = (0.0, 0.0), 10.0
        cb, rb = (1.0, 0.0), 1.0   # d=1, ra-rb=9 → cos_t=9 before clamping
        pa, pb, nx, ny = _ext_tangent_outside(ca, ra, cb, rb)
        assert math.hypot(nx, ny) == pytest.approx(1.0, abs=1e-9)


# ---------------------------------------------------------------------------
# _sample_arc_cw
# ---------------------------------------------------------------------------

class TestSampleArcCw:
    def test_all_points_lie_on_circle(self):
        """Every returned point is at distance `radius` from center."""
        center = (0.0, 0.0)
        r = 5.0
        pts = _sample_arc_cw(center, r, (0.0, 5.0), (5.0, 0.0), n=10)
        for pt in pts:
            assert math.hypot(pt[0] - center[0], pt[1] - center[1]) == pytest.approx(r, abs=1e-9)

    def test_endpoints_excluded(self):
        """Start and end points are not in the returned list."""
        start = (0.0, 5.0)
        end = (5.0, 0.0)
        pts = _sample_arc_cw((0.0, 0.0), 5.0, start, end, n=6)
        for pt in pts:
            assert pt != pytest.approx(start, abs=1e-9)
            assert pt != pytest.approx(end, abs=1e-9)

    def test_count_is_n_minus_one(self):
        """Returns exactly n-1 points (range(1, n))."""
        for n in (4, 10, 30):
            pts = _sample_arc_cw((0.0, 0.0), 5.0, (0.0, 5.0), (5.0, 0.0), n=n)
            assert len(pts) == n - 1

    def test_cw_ordering(self):
        """Arc from north (0,r) to east (r,0) with n=4: x values increase CW."""
        center = (0.0, 0.0)
        r = 5.0
        pts = _sample_arc_cw(center, r, (0.0, r), (r, 0.0), n=4)
        assert len(pts) == 3
        # Going CW from north to east: x increases (0 → r)
        xs = [p[0] for p in pts]
        assert xs[0] < xs[1] < xs[2]


# ---------------------------------------------------------------------------
# compute_db_roof_outline
# ---------------------------------------------------------------------------

class TestComputeDbRoofOutline:
    def test_raises_for_fewer_than_3_corners(self):
        with pytest.raises(ValueError):
            compute_db_roof_outline(
                ["C1", "C2"], [False, False],
                {"C1": (0, 0), "C2": (5, 0)},
                {"R_a1": 1.0, "R_a2": 1.0},
                overhang=0.5,
            )

    def test_returns_dbroofresult(self):
        result = compute_db_roof_outline(
            _SQ_NAMES, [False] * 4, _SQ_PTS, _SQ_RADII, _OVERHANG,
        )
        assert isinstance(result, DbRoofResult)

    def test_sharp_corners_one_vertex_per_corner(self):
        """4 sharp corners → polygon has exactly 4 vertices."""
        result = compute_db_roof_outline(
            _SQ_NAMES, [False] * 4, _SQ_PTS, _SQ_RADII, _OVERHANG,
        )
        assert len(result.poly) == 4

    def test_corner_pts_keys_match_corner_names(self):
        result = compute_db_roof_outline(
            _SQ_NAMES, [False] * 4, _SQ_PTS, _SQ_RADII, _OVERHANG,
        )
        assert set(result.corner_pts.keys()) == set(_SQ_NAMES)

    def test_corner_radii_equals_outline_plus_overhang(self):
        result = compute_db_roof_outline(
            _SQ_NAMES, [False] * 4, _SQ_PTS, _SQ_RADII, _OVERHANG,
        )
        for name in _SQ_NAMES:
            assert result.corner_radii[name] == pytest.approx(1.0 + _OVERHANG, abs=1e-9)

    def test_area_positive_and_larger_than_building(self):
        """Roof outline encloses more area than the 10×10 building interior."""
        result = compute_db_roof_outline(
            _SQ_NAMES, [False] * 4, _SQ_PTS, _SQ_RADII, _OVERHANG,
        )
        assert result.area > 0
        assert result.area > 100.0  # 10×10 building

    def test_r_series_names_map_to_c_series_centers(self):
        """R01 → looks up pts['C01'] via _cname."""
        pts_r = {"C01": (5, 5), "C02": (5, -5), "C03": (-5, -5)}
        radii_r = {"R_a01": 1.0, "R_a02": 1.0, "R_a03": 1.0}
        result = compute_db_roof_outline(
            ["R01", "R02", "R03"], [False, False, False],
            pts_r, radii_r, overhang=0.5,
        )
        assert isinstance(result, DbRoofResult)
        assert len(result.poly) == 3

    def test_radiused_corners_add_arc_samples(self):
        """Radiused corners insert arc points → poly longer than N."""
        result = compute_db_roof_outline(
            _SQ_NAMES, [True] * 4, _SQ_PTS, _SQ_RADII, _OVERHANG,
        )
        assert len(result.poly) > 4

    def test_radiused_corner_pts_are_float_tuples(self):
        """Label positions for radiused corners are (float, float) tuples."""
        result = compute_db_roof_outline(
            _SQ_NAMES, [True] * 4, _SQ_PTS, _SQ_RADII, _OVERHANG,
        )
        for name in _SQ_NAMES:
            pos = result.corner_pts[name]
            assert len(pos) == 2
            assert isinstance(pos[0], float)
            assert isinstance(pos[1], float)

    def test_sharp_vertices_are_outside_corner_circles(self):
        """Sharp corner vertices must be farther from corner center than r_roof."""
        result = compute_db_roof_outline(
            _SQ_NAMES, [False] * 4, _SQ_PTS, _SQ_RADII, _OVERHANG,
        )
        r_roof = 1.0 + _OVERHANG
        for i, name in enumerate(_SQ_NAMES):
            center = _SQ_PTS[name]
            vertex = result.poly[i]
            dist = math.hypot(vertex[0] - center[0], vertex[1] - center[1])
            assert dist >= r_roof - 1e-6


# ---------------------------------------------------------------------------
# db_roof_segments
# ---------------------------------------------------------------------------

class TestDbRoofSegments:
    def test_raises_for_fewer_than_3_corners(self):
        with pytest.raises(ValueError):
            db_roof_segments(
                ["C1", "C2"], [False, False],
                {"C1": (0, 0), "C2": (5, 0)},
                {"R_a1": 1.0, "R_a2": 1.0},
                overhang=0.5,
            )

    def test_returns_list(self):
        segs = db_roof_segments(
            _SQ_NAMES, [False] * 4, _SQ_PTS, _SQ_RADII, _OVERHANG,
        )
        assert isinstance(segs, list)

    def test_all_elements_are_line_or_arc(self):
        segs = db_roof_segments(
            _SQ_NAMES, [False] * 4, _SQ_PTS, _SQ_RADII, _OVERHANG,
        )
        for seg in segs:
            assert seg[0] in ("line", "arc")

    def test_all_sharp_corners_give_n_line_segments(self):
        """N sharp corners → exactly N line segments, no arcs."""
        n = 4
        segs = db_roof_segments(
            _SQ_NAMES, [False] * n, _SQ_PTS, _SQ_RADII, _OVERHANG,
        )
        assert len(segs) == n
        assert all(s[0] == "line" for s in segs)

    def test_line_segment_continuity(self):
        """End of each line segment is the start of the next (closed polygon)."""
        segs = db_roof_segments(
            _SQ_NAMES, [False] * 4, _SQ_PTS, _SQ_RADII, _OVERHANG,
        )
        n = len(segs)
        for i in range(n):
            end_x, end_y = segs[i][3], segs[i][4]
            start_x, start_y = segs[(i + 1) % n][1], segs[(i + 1) % n][2]
            assert end_x == pytest.approx(start_x, abs=1e-6)
            assert end_y == pytest.approx(start_y, abs=1e-6)
