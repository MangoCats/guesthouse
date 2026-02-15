"""Tests for survey/gen_path_svg.py SVG generation."""
import sys, os
import pytest

# survey/ is not a package — add it to sys.path for import
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "survey"))

from gen_path_svg import (
    compute_all, render_layer, build_outline_cfg, render_floorplan,
    outer_cfg, inset_cfg,
)


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


class TestRenderLayerBranches:
    """Cover opacity wrapper and traverse overlay in render_layer."""

    def test_opacity_wraps_in_group(self, all_data):
        lines = []
        render_layer(lines, all_data["outer_segs"], all_data["pts_rot"],
                     outer_cfg, all_data["to_svg"])
        joined = "\n".join(lines)
        assert '<g opacity="' in joined
        assert '</g>' in joined

    def test_traverse_overlay_present(self, all_data):
        lines = []
        render_layer(lines, all_data["outer_segs"], all_data["pts_rot"],
                     outer_cfg, all_data["to_svg"])
        joined = "\n".join(lines)
        assert 'stroke-dasharray="4,4"' in joined

    def test_inset_layer_has_opacity(self, all_data):
        lines = []
        render_layer(lines, all_data["inset_segs"], all_data["pts_rot"],
                     inset_cfg, all_data["to_svg"])
        joined = "\n".join(lines)
        assert '<g opacity="' in joined


class TestRenderFloorplan:
    """Cover the render_floorplan function."""

    def test_produces_svg_content(self, all_data):
        lines = []
        render_floorplan(
            lines, all_data["to_svg"], all_data["pts"],
            all_data["outer_poly"], all_data["inner_poly"],
            all_data["inner_segs"], all_data["layout"])
        assert len(lines) > 0
        joined = "\n".join(lines)
        assert '<g opacity="0.5">' in joined
        assert '</g>' in joined

    def test_wall_band_path(self, all_data):
        lines = []
        render_floorplan(
            lines, all_data["to_svg"], all_data["pts"],
            all_data["outer_poly"], all_data["inner_poly"],
            all_data["inner_segs"], all_data["layout"])
        joined = "\n".join(lines)
        assert 'fill-rule="evenodd"' in joined

    def test_appliance_labels(self, all_data):
        lines = []
        render_floorplan(
            lines, all_data["to_svg"], all_data["pts"],
            all_data["outer_poly"], all_data["inner_poly"],
            all_data["inner_segs"], all_data["layout"])
        joined = "\n".join(lines)
        for label in ["DRYER", "WASHER", "COUNTER", "KING BED"]:
            assert label in joined, f"Missing label {label}"

    def test_room_labels(self, all_data):
        lines = []
        render_floorplan(
            lines, all_data["to_svg"], all_data["pts"],
            all_data["outer_poly"], all_data["inner_poly"],
            all_data["inner_segs"], all_data["layout"])
        joined = "\n".join(lines)
        assert "CLOSET" in joined
        assert "BEDROOM" in joined

    def test_iw_polygons(self, all_data):
        import re
        lines = []
        render_floorplan(
            lines, all_data["to_svg"], all_data["pts"],
            all_data["outer_poly"], all_data["inner_poly"],
            all_data["inner_segs"], all_data["layout"])
        joined = "\n".join(lines)
        iw_fills = re.findall(r'fill="rgba\(160,160,160,0\.5\)"', joined)
        assert len(iw_fills) >= 6
