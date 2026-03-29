"""Tests for plumbing/gen_plumbing.py — boundary computation and SVG output."""
import math
import pytest
from conftest import _import_from

_plumbing = _import_from("plumbing", "gen_plumbing")
_compute_boundary_corners = _plumbing._compute_boundary_corners


# ---------------------------------------------------------------------------
# _compute_boundary_corners
# ---------------------------------------------------------------------------

class TestComputeBoundaryCorners:
    """Tests use pts from pts_with_outline which contains F-series and survey points."""

    @pytest.fixture(scope="class")
    def boundary(self, pts_with_outline):
        pts, _ = pts_with_outline
        return _compute_boundary_corners(pts)

    def test_returns_all_required_keys(self, boundary):
        assert set(boundary.keys()) == {"sw", "south_start", "west_start",
                                        "south_label", "west_label"}

    def test_south_label_value(self, boundary):
        assert boundary["south_label"] == "216.73'"

    def test_west_label_value(self, boundary):
        assert boundary["west_label"] == "275.08'"

    def test_sw_is_two_floats(self, boundary):
        sw = boundary["sw"]
        assert len(sw) == 2
        assert isinstance(sw[0], float)
        assert isinstance(sw[1], float)

    def test_south_start_is_two_floats(self, boundary):
        ss = boundary["south_start"]
        assert len(ss) == 2
        assert isinstance(ss[0], float)
        assert isinstance(ss[1], float)

    def test_west_start_is_two_floats(self, boundary):
        ws = boundary["west_start"]
        assert len(ws) == 2
        assert isinstance(ws[0], float)
        assert isinstance(ws[1], float)

    def test_south_start_is_east_of_sw(self, boundary):
        """South boundary runs NE from SW corner; clipped start should be east of sw."""
        assert boundary["south_start"][0] >= boundary["sw"][0] - 1e-6

    def test_west_start_is_north_of_sw(self, boundary):
        """West boundary runs north from SW corner; clipped start should be above sw."""
        assert boundary["west_start"][1] >= boundary["sw"][1] - 1e-6

    def test_sw_is_near_building_at_plausible_distance(self, boundary):
        """SW property corner should be within 200 ft of building origin in both axes."""
        sw = boundary["sw"]
        assert abs(sw[0]) < 200
        assert abs(sw[1]) < 200

    def test_south_start_not_coincident_with_sw(self, boundary):
        """south_start must differ from sw (boundary was clipped to building extent)."""
        sw = boundary["sw"]
        ss = boundary["south_start"]
        dist = math.hypot(ss[0] - sw[0], ss[1] - sw[1])
        assert dist > 0.1  # at least 0.1 ft apart


# ---------------------------------------------------------------------------
# SVG generation smoke test
# ---------------------------------------------------------------------------

class TestPlumbingSvgGeneration:
    """Smoke-test the full plumbing SVG pipeline using seed geometry."""

    @pytest.fixture(scope="class")
    def plumbing_svg(self, pts_with_outline):
        from floorplan.gen_floorplan import build_floorplan_data, render_floorplan_svg
        pts, _ = pts_with_outline
        data = build_floorplan_data()
        boundary = _compute_boundary_corners(data.pts)
        return render_floorplan_svg(
            data, room_title="Plumbing Plan",
            db=True, plumbing=True, boundary=boundary,
        )

    def test_svg_is_string(self, plumbing_svg):
        assert isinstance(plumbing_svg, str)

    def test_svg_starts_and_ends_correctly(self, plumbing_svg):
        stripped = plumbing_svg.strip()
        assert stripped.startswith("<svg")
        assert stripped.endswith("</svg>")

    def test_svg_contains_plumbing_title(self, plumbing_svg):
        assert "Plumbing Plan" in plumbing_svg
