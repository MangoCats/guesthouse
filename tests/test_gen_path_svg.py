"""Tests for survey/gen_path_svg.py SVG generation."""
import sys, os
import pytest

# survey/ is not a package — add it to sys.path for import
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "survey"))

from gen_path_svg import compute_all, render_layer, build_outline_cfg


@pytest.fixture(scope="module")
def all_data():
    return compute_all()


class TestComputeAll:
    def test_returns_expected_keys(self, all_data):
        expected = {"pts", "pts_rot", "to_svg",
                    "outer_segs", "inset_segs", "outline_segs",
                    "outer_area", "inset_area", "outline_area",
                    "radii", "outer_poly", "inner_poly", "inner_segs", "layout"}
        assert expected.issubset(all_data.keys())

    def test_areas_positive(self, all_data):
        assert all_data["outer_area"] > 0
        assert all_data["inset_area"] > 0
        assert all_data["outline_area"] > 0

    def test_area_ordering(self, all_data):
        # Outer traverse > inset > outline (building envelope)
        assert all_data["outer_area"] > all_data["inset_area"]
        assert all_data["inset_area"] > all_data["outline_area"]


class TestBuildOutlineCfg:
    def test_returns_layer_config(self, all_data):
        cfg = build_outline_cfg(
            all_data["outline_segs"], all_data["pts"], all_data["radii"])
        assert cfg.opacity == 1.0
        assert isinstance(cfg.arc_styles, dict)

    def test_has_13_arc_styles(self, all_data):
        cfg = build_outline_cfg(
            all_data["outline_segs"], all_data["pts"], all_data["radii"])
        assert len(cfg.arc_styles) == 13


class TestRenderLayer:
    def test_produces_svg_content(self, all_data):
        cfg = build_outline_cfg(
            all_data["outline_segs"], all_data["pts"], all_data["radii"])
        lines = []
        render_layer(lines, all_data["outline_segs"], all_data["pts"],
                     cfg, all_data["to_svg"])
        joined = "\n".join(lines)
        assert "<polygon" in joined
        assert len(lines) > 0
