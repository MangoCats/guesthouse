"""Tests for walls/gen_walls.py SVG generation."""
import math
import pytest
from shared.types import LineSeg, ArcSeg
from walls.constants import SHELL_THICKNESS, AIR_GAP, OPENING_INSIDE_RADIUS
from walls.gen_walls import (
    build_wall_data,
    render_walls_svg,
    _compute_inset_path,
    _solid_ranges,
    _openings_on_seg,
    _lerp,
    _line_strip_poly,
    _partial_line_strip,
)
from floorplan.openings import WallOpening


# ============================================================
# Helper unit tests
# ============================================================

class TestLerp:
    def test_at_zero(self):
        assert _lerp((0, 0), (10, 20), 0.0) == (0, 0)

    def test_at_one(self):
        assert _lerp((0, 0), (10, 20), 1.0) == (10, 20)

    def test_midpoint(self):
        result = _lerp((0, 0), (10, 20), 0.5)
        assert result == pytest.approx((5, 10))


class TestSolidRanges:
    def test_no_openings(self):
        assert _solid_ranges([]) == [(0.0, 1.0)]

    def test_one_opening_in_middle(self):
        op = WallOpening("O1", 0, 0.3, 0.7)
        ranges = _solid_ranges([op])
        assert len(ranges) == 2
        assert ranges[0] == pytest.approx((0.0, 0.3))
        assert ranges[1] == pytest.approx((0.7, 1.0))

    def test_opening_at_start(self):
        op = WallOpening("O1", 0, 0.0, 0.5)
        ranges = _solid_ranges([op])
        assert len(ranges) == 1
        assert ranges[0] == pytest.approx((0.5, 1.0))

    def test_two_openings(self):
        ops = [WallOpening("O1", 0, 0.2, 0.4), WallOpening("O2", 0, 0.6, 0.8)]
        ranges = _solid_ranges(ops)
        assert len(ranges) == 3
        assert ranges[0] == pytest.approx((0.0, 0.2))
        assert ranges[1] == pytest.approx((0.4, 0.6))
        assert ranges[2] == pytest.approx((0.8, 1.0))


class TestLineStripPoly:
    def test_returns_four_points(self):
        pts = {"A": (0.0, 0.0), "B": (10.0, 0.0), "C": (0.0, 1.0), "D": (10.0, 1.0)}
        result = _line_strip_poly(pts, "A", "B", "C", "D")
        assert len(result) == 4

    def test_correct_winding(self):
        pts = {"A": (0.0, 0.0), "B": (10.0, 0.0), "C": (0.0, 1.0), "D": (10.0, 1.0)}
        result = _line_strip_poly(pts, "A", "B", "C", "D")
        # Should be: A, B, D, C (outer start, outer end, inner end, inner start)
        assert result[0] == (0.0, 0.0)
        assert result[1] == (10.0, 0.0)
        assert result[2] == (10.0, 1.0)
        assert result[3] == (0.0, 1.0)


class TestPartialLineStrip:
    def test_full_range(self):
        pts = {"A": (0.0, 0.0), "B": (10.0, 0.0)}
        seg = LineSeg("A", "B")
        result = _partial_line_strip(pts, seg, "A", "B", 0.0, 1.0)
        assert len(result) == 4
        assert result[0] == pytest.approx((0.0, 0.0))
        assert result[1] == pytest.approx((10.0, 0.0))

    def test_sub_range(self):
        pts = {"A": (0.0, 0.0), "B": (10.0, 0.0), "C": (0.0, 1.0), "D": (10.0, 1.0)}
        seg = LineSeg("A", "B")
        result = _partial_line_strip(pts, seg, "C", "D", 0.2, 0.8)
        assert len(result) == 4
        assert result[0] == pytest.approx((2.0, 0.0))
        assert result[1] == pytest.approx((8.0, 0.0))
        assert result[2] == pytest.approx((8.0, 1.0))
        assert result[3] == pytest.approx((2.0, 1.0))


# ============================================================
# Integration tests
# ============================================================

@pytest.fixture(scope="module")
def wall_data():
    return build_wall_data()


@pytest.fixture(scope="module")
def rendered(wall_data):
    return render_walls_svg(wall_data)


class TestBuildWallData:
    def test_returns_namedtuple_with_expected_fields(self, wall_data):
        expected = {
            "pts", "to_svg", "outline_segs", "inner_segs",
            "s_segs", "g_segs", "radii", "openings",
            "vb_x", "vb_y", "vb_w", "vb_h",
        }
        assert expected.issubset(set(wall_data._fields))

    def test_s_series_points_exist(self, wall_data):
        pts = wall_data.pts
        for i in range(22):
            assert f"S{i}" in pts, f"Missing S{i}"

    def test_g_series_points_exist(self, wall_data):
        pts = wall_data.pts
        for i in range(22):
            assert f"G{i}" in pts, f"Missing G{i}"

    def test_shell_distances(self, wall_data):
        """Shell boundary distances from F-series should match expected insets."""
        pts = wall_data.pts
        for i in range(22):
            f_pt = pts[f"F{i}"]
            s_pt = pts[f"S{i}"]
            g_pt = pts[f"G{i}"]
            w_pt = pts[f"W{i}"]

            fs_dist = math.sqrt((f_pt[0] - s_pt[0])**2 + (f_pt[1] - s_pt[1])**2)
            fg_dist = math.sqrt((f_pt[0] - g_pt[0])**2 + (f_pt[1] - g_pt[1])**2)
            fw_dist = math.sqrt((f_pt[0] - w_pt[0])**2 + (f_pt[1] - w_pt[1])**2)

            assert fs_dist == pytest.approx(SHELL_THICKNESS, abs=0.01), \
                f"F{i}-S{i} distance {fs_dist} != {SHELL_THICKNESS}"
            assert fg_dist == pytest.approx(SHELL_THICKNESS + AIR_GAP, abs=0.01), \
                f"F{i}-G{i} distance {fg_dist} != {SHELL_THICKNESS + AIR_GAP}"
            assert fw_dist == pytest.approx(SHELL_THICKNESS * 2 + AIR_GAP, abs=0.01), \
                f"F{i}-W{i} distance {fw_dist} != {SHELL_THICKNESS * 2 + AIR_GAP}"

    def test_22_outline_segments(self, wall_data):
        assert len(wall_data.outline_segs) == 22

    def test_22_s_segments(self, wall_data):
        assert len(wall_data.s_segs) == 22

    def test_22_g_segments(self, wall_data):
        assert len(wall_data.g_segs) == 22

    def test_11_openings(self, wall_data):
        assert len(wall_data.openings) == 11

    def test_opening_names(self, wall_data):
        names = {o.name for o in wall_data.openings}
        expected = {f"O{i}" for i in range(1, 12)}
        assert names == expected

    def test_openings_on_line_segs(self, wall_data):
        """All openings should be on LineSeg segments."""
        outline_segs = wall_data.outline_segs
        for o in wall_data.openings:
            seg = outline_segs[o.seg_idx]
            assert isinstance(seg, LineSeg), \
                f"Opening {o.name} on non-line segment {type(seg)}"

    def test_opening_params_valid(self, wall_data):
        """Opening parametric ranges should be within [0, 1]."""
        for o in wall_data.openings:
            assert 0.0 <= o.t_start < o.t_end <= 1.0, \
                f"Opening {o.name}: t_start={o.t_start}, t_end={o.t_end}"


class TestRenderWallsSvg:
    def test_svg_envelope(self, rendered):
        assert rendered.strip().startswith("<svg")
        assert rendered.strip().endswith("</svg>")

    def test_contains_all_opening_labels(self, rendered):
        for i in range(1, 12):
            assert f">O{i}<" in rendered, f"Missing opening label O{i}"

    def test_wall_polygon_count(self, rendered):
        import re
        wall_fills = re.findall(r'fill="rgba\(180,180,180,0\.5\)"', rendered)
        assert len(wall_fills) == 88

    def test_opening_polygon_count(self, rendered):
        import re
        opening_fills = re.findall(r'fill="rgb\(220,235,255\)"', rendered)
        assert len(opening_fills) == 11

    def test_title_present(self, rendered):
        assert "Outer Walls" in rendered

    def test_wall_construction_note(self, rendered):
        assert "shell" in rendered and "gap" in rendered

    def test_scale_info_present(self, rendered):
        assert "Scale" in rendered
