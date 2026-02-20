"""Tests for scad/gen_flat_roof.py OpenSCAD generation."""
import io
import pytest
from unittest.mock import patch, mock_open
from conftest import _import_from

_mod = _import_from("scad", "gen_flat_roof")
_scad_polygon = _mod._scad_polygon
generate = _mod.generate
WALL_HEIGHT_FT = _mod.WALL_HEIGHT_FT


# ── unit tests ────────────────────────────────────────────────

class TestScadPolygon:
    def test_module_declaration(self):
        elems = [("pts", [(0.0, 0.0), (1.0, 0.0)])]
        result = _scad_polygon("test_wall", elems)
        assert "module test_wall()" in result

    def test_polygon_concat(self):
        elems = [("pts", [(0.0, 0.0), (1.0, 0.0)])]
        result = _scad_polygon("w", elems)
        assert "polygon(points = concat(" in result

    def test_arc_element(self):
        elems = [("pts", [(0.0, 0.0), (1.0, 0.0)]),
                 ("arc", 1.0, 1.0, 1.0, 0.0, 90.0)]
        result = _scad_polygon("a", elems)
        assert "arc_pts(" in result
        assert "tail(" in result

    def test_coordinate_precision(self):
        elems = [("pts", [(1.5, 2.5)])]
        result = _scad_polygon("p", elems)
        assert "[1.500000, 2.500000]" in result


class TestWallHeight:
    def test_conversion(self):
        assert WALL_HEIGHT_FT == pytest.approx(80.0 / 12.0)


# ── integration test ──────────────────────────────────────────

class TestGenerate:
    def test_output_contains_wall_sections(self):
        buf = io.StringIO()
        with patch("builtins.open", return_value=buf):
            with patch.object(buf, "close"):  # prevent StringIO close
                generate()
        content = buf.getvalue()
        # 11 wall sections between openings
        assert "module wall_O11_O1_outer()" in content
        assert "module wall_O11_O1_cavity()" in content
        assert "module wall_O1_O2_outer()" in content
        assert "module wall_O10_O11_outer()" in content

    def test_output_contains_assembly(self):
        buf = io.StringIO()
        with patch("builtins.open", return_value=buf):
            with patch.object(buf, "close"):
                generate()
        content = buf.getvalue()
        assert "linear_extrude" in content
        assert "union()" in content
        assert "difference()" in content
        assert "wall_O11_O1_outer();" in content
        assert "wall_O11_O1_cavity();" in content

    def test_output_contains_helpers(self):
        buf = io.StringIO()
        with patch("builtins.open", return_value=buf):
            with patch.object(buf, "close"):
                generate()
        content = buf.getvalue()
        assert "function arc_pts(" in content
        assert "function tail(" in content

    def test_output_uses_arcs(self):
        buf = io.StringIO()
        with patch("builtins.open", return_value=buf):
            with patch.object(buf, "close"):
                generate()
        content = buf.getvalue()
        assert "arc_pts(" in content
        assert "tail(" in content
        assert "concat(" in content
