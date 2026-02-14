"""Tests for floorplan/gen_floorplan.py SVG generation."""
import pytest
from shared.types import LineSeg, ArcSeg
from floorplan.gen_floorplan import (
    build_floorplan_data,
    render_floorplan_svg,
    dim_line_h, dim_line_v, wall_poly, stroke_segs,
)


# --- Mock transform for helper unit tests ---
def _mock_to_svg(e, n):
    return (e * 10, -n * 10)


# ============================================================
# Helper unit tests
# ============================================================

class TestDimLineH:
    def test_produces_line_and_text(self):
        out = []
        dim_line_h(out, 0, 5, 10, "10'", _mock_to_svg)
        joined = "\n".join(out)
        assert "<line" in joined
        assert "<text" in joined
        assert "10'" in joined

    def test_tick_marks(self):
        out = []
        dim_line_h(out, 0, 5, 10, "10'", _mock_to_svg)
        # main line + 2 tick marks + text = 4 elements
        assert len(out) == 4


class TestDimLineV:
    def test_produces_line_and_rotated_text(self):
        out = []
        dim_line_v(out, 5, 0, 10, "10'", _mock_to_svg)
        joined = "\n".join(out)
        assert "<line" in joined
        assert "rotate(-90" in joined
        assert "10'" in joined

    def test_tick_marks(self):
        out = []
        dim_line_v(out, 5, 0, 10, "10'", _mock_to_svg)
        # main line + 2 tick marks + text = 4 elements
        assert len(out) == 4


class TestWallPoly:
    def test_produces_polygon(self):
        out = []
        wall_poly(out, [(0, 0), (1, 0), (1, 1), (0, 1)], _mock_to_svg)
        assert len(out) == 2  # clipPath defs + polygon
        assert "<polygon" in out[1]
        assert "fill=" in out[1]

    def test_stroke_option(self):
        out_stroke = []
        wall_poly(out_stroke, [(0, 0), (1, 0), (1, 1)], _mock_to_svg, stroke=True)
        assert 'stroke="#666"' in out_stroke[1]
        assert 'clip-path=' in out_stroke[1]

        out_no_stroke = []
        wall_poly(out_no_stroke, [(0, 0), (1, 0), (1, 1)], _mock_to_svg, stroke=False)
        assert 'stroke="none"' in out_no_stroke[0]



class TestStrokeSegs:
    def test_line_seg(self):
        pts = {"A": (0.0, 0.0), "B": (5.0, 5.0)}
        segs = [LineSeg("A", "B")]
        out = []
        stroke_segs(out, segs, "red", 2, pts, _mock_to_svg)
        assert len(out) == 1
        assert "<line" in out[0]
        assert 'stroke="red"' in out[0]

    def test_arc_seg(self):
        pts = {"A": (1.0, 0.0), "B": (0.0, 1.0), "C": (0.0, 0.0)}
        segs = [ArcSeg("A", "B", "C", 1.0, "CCW", 60)]
        out = []
        stroke_segs(out, segs, "blue", 1.5, pts, _mock_to_svg)
        assert len(out) == 1
        assert "<polyline" in out[0]
        assert 'stroke="blue"' in out[0]


# ============================================================
# Integration tests
# ============================================================

@pytest.fixture(scope="module")
def floorplan_data():
    return build_floorplan_data()


@pytest.fixture(scope="module")
def rendered(floorplan_data):
    svg_content, inner_area, outer_area, vert_names = render_floorplan_svg(floorplan_data)
    return svg_content, inner_area, outer_area, vert_names


class TestBuildFloorplanData:
    def test_returns_dict_with_expected_keys(self, floorplan_data):
        expected = {"pts", "to_svg", "inner_area", "outer_area",
                    "outline_segs", "inner_segs", "outer_poly", "inner_poly",
                    "radii", "wall_t", "vb_x", "vb_y", "vb_w", "vb_h"}
        assert expected.issubset(floorplan_data.keys())

    def test_areas_positive(self, floorplan_data):
        assert floorplan_data["inner_area"] > 0
        assert floorplan_data["outer_area"] > 0

    def test_outer_area_greater_than_inner(self, floorplan_data):
        assert floorplan_data["outer_area"] > floorplan_data["inner_area"]

    def test_pts_contain_f_and_w_series(self, floorplan_data):
        pts = floorplan_data["pts"]
        for i in range(22):
            assert f"F{i}" in pts, f"Missing F{i}"
            assert f"W{i}" in pts, f"Missing W{i}"


class TestRenderFloorplanSvg:
    def test_svg_envelope(self, rendered):
        svg, _, _, _ = rendered
        assert svg.strip().startswith("<svg")
        assert svg.strip().endswith("</svg>")

    def test_vert_names_returned(self, rendered):
        _, _, _, vert_names = rendered
        assert len(vert_names) == 22, f"Expected 22 vertex names, got {len(vert_names)}"


    def test_contains_appliance_labels(self, rendered):
        svg, _, _, _ = rendered
        for label in ["DRYER", "WASHER", "COUNTER", "WH"]:
            assert label in svg, f"Missing appliance label {label}"

    def test_contains_bed_label(self, rendered):
        svg, _, _, _ = rendered
        assert "KING BED" in svg

    def test_contains_openings(self, rendered):
        svg, _, _, _ = rendered
        assert svg.count('fill="rgb(220,235,255)"') == 11, "Expected 11 opening polygons"

    def test_rendered_inner_area_less_than_data(self, rendered, floorplan_data):
        """Rendered inner area subtracts interior walls, so it's less than the raw polygon area."""
        _, inner_area, outer_area, _ = rendered
        assert inner_area < floorplan_data["inner_area"]
        assert inner_area > 0

    def test_outer_area_matches_data(self, rendered, floorplan_data):
        _, _, outer_area, _ = rendered
        assert outer_area == pytest.approx(floorplan_data["outer_area"])
