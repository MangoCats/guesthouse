"""Tests for roof/gen_roof.py SVG generation."""
import re
import pytest
from floorplan.gen_floorplan import build_floorplan_data
from floorplan.roof import compute_roof_geometry
from conftest import _import_from

_mod = _import_from("roof", "gen_roof")
_svg_pts = _mod._svg_pts
render_roof_svg = _mod.render_roof_svg


@pytest.fixture(scope="module")
def fp_data():
    return build_floorplan_data()


@pytest.fixture(scope="module")
def roof_geo(fp_data):
    return compute_roof_geometry(fp_data.pts, fp_data.radii)


@pytest.fixture(scope="module")
def rendered(fp_data, roof_geo):
    return render_roof_svg(fp_data, roof_geo)


# ── unit tests ────────────────────────────────────────────────

class TestSvgPts:
    def test_single_point(self):
        to_svg = lambda e, n: (e * 10, -n * 10)
        result = _svg_pts([(1.0, 2.0)], to_svg)
        assert result == "10.00,-20.00"

    def test_multiple_points(self):
        to_svg = lambda e, n: (e, n)
        result = _svg_pts([(0.0, 0.0), (1.0, 2.0), (3.0, 4.0)], to_svg)
        parts = result.split(" ")
        assert len(parts) == 3


# ── integration tests ─────────────────────────────────────────

class TestRenderRoofSvg:
    def test_svg_envelope(self, rendered):
        assert rendered.startswith("<svg")
        assert rendered.rstrip().endswith("</svg>")

    def test_contains_title(self, rendered):
        assert "Roof Plan" in rendered

    def test_building_outline_polygon(self, rendered):
        assert 'fill="rgba(220,220,220,0.3)"' in rendered

    def test_roof_outline_polygon(self, rendered):
        assert 'stroke-dasharray="3,2"' in rendered

    def test_r_series_labels(self, rendered):
        for name in ["R1", "R5", "R6", "R7"]:
            assert f">{name}<" in rendered

    def test_r3_r4_fillet_labels(self, rendered):
        assert ">R3<" in rendered
        assert ">R4<" in rendered

    def test_title_block_present(self, rendered):
        assert "Area:" in rendered

    def test_area_positive(self, roof_geo):
        assert roof_geo.area > 0
