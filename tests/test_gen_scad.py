"""Tests for scad/gen_flat_roof.py OpenSCAD generation (T-path approach)."""
import io
import pytest
from unittest.mock import patch
from conftest import _import_from

_mod = _import_from("scad", "gen_flat_roof")
_scad_seg = _mod._scad_seg
_rev_elem = _mod._rev_elem
generate = _mod.generate
WALL_HEIGHT_FT = _mod.WALL_HEIGHT_FT


# ── unit tests ────────────────────────────────────────────────

class TestScadSeg:
    def test_line_format(self):
        elem = ("line", 1.0, 2.0, 3.0, 4.0)
        result = _scad_seg(elem)
        assert result.startswith("[0,")
        assert "1.000000" in result
        assert "4.000000" in result

    def test_arc_format(self):
        elem = ("arc", 5.0, 6.0, 1.5, 0.0, 90.0)
        result = _scad_seg(elem)
        assert result.startswith("[1,")
        assert "5.000000" in result
        assert "1.500000" in result
        assert "90.000" in result


class TestRevElem:
    def test_reverse_line(self):
        elem = ("line", 1.0, 2.0, 3.0, 4.0)
        rev = _rev_elem(elem)
        assert rev == ("line", 3.0, 4.0, 1.0, 2.0)

    def test_reverse_arc(self):
        elem = ("arc", 5.0, 6.0, 1.5, 0.0, 90.0)
        rev = _rev_elem(elem)
        assert rev == ("arc", 5.0, 6.0, 1.5, 90.0, 0.0)


class TestWallHeight:
    def test_conversion(self):
        assert WALL_HEIGHT_FT == pytest.approx(80.0 / 12.0)


# ── integration test ──────────────────────────────────────────

class TestGenerate:
    def test_output_contains_tpath_data(self):
        buf = io.StringIO()
        with patch("builtins.open", return_value=buf):
            with patch.object(buf, "close"):  # prevent StringIO close
                generate()
        content = buf.getvalue()
        # 11 wall sections, each with a T-path variable
        assert "t_O11_O1 = [" in content
        assert "t_O1_O2 = [" in content
        assert "t_O10_O11 = [" in content

    def test_output_contains_assembly(self):
        buf = io.StringIO()
        with patch("builtins.open", return_value=buf):
            with patch.object(buf, "close"):
                generate()
        content = buf.getvalue()
        assert "linear_extrude" in content
        assert "union()" in content
        assert "wall_shell(t_O11_O1, half_t);" in content

    def test_output_contains_helpers(self):
        buf = io.StringIO()
        with patch("builtins.open", return_value=buf):
            with patch.object(buf, "close"):
                generate()
        content = buf.getvalue()
        assert "function shell_pts(" in content
        assert "function seg_pts(" in content
        assert "function line_pts(" in content
        assert "function arc_off_pts(" in content
        assert "function tail(" in content
        assert "module wall_shell(" in content

    def test_output_has_segment_types(self):
        buf = io.StringIO()
        with patch("builtins.open", return_value=buf):
            with patch.object(buf, "close"):
                generate()
        content = buf.getvalue()
        # Lines start with [0, and arcs with [1,
        assert "[0, " in content
        assert "[1, " in content

    def test_no_difference_in_assembly(self):
        buf = io.StringIO()
        with patch("builtins.open", return_value=buf):
            with patch.object(buf, "close"):
                generate()
        content = buf.getvalue()
        # T-path approach uses polygon with paths, no difference()
        assert "difference()" not in content

    def test_output_two_tier_walls(self):
        buf = io.StringIO()
        with patch("builtins.open", return_value=buf):
            with patch.object(buf, "close"):
                generate()
        content = buf.getvalue()
        assert "t_full_O4 = [" in content
        assert "upper_base" in content
        assert "upper_height" in content
        assert "wall_shell(t_full_O4, half_t);" in content

    def test_output_wedge_roof(self):
        buf = io.StringIO()
        with patch("builtins.open", return_value=buf):
            with patch.object(buf, "close"):
                generate()
        content = buf.getvalue()
        assert "translate([0, 0, upper_base + upper_height])" in content
        assert "polygon(points = shell_pts(roof_outline, 0))" in content
        assert "max_roof_thick" in content
        assert "roof_slope" in content
        assert "roof_z_base" in content
        assert "render() intersection()" in content
