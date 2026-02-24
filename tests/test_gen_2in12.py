"""Tests for scad/gen_2in12.py OpenSCAD generation (2:12 slope)."""
import re
import pytest
from conftest import _import_from, generate_scad_content

_mod = _import_from("scad", "gen_2in12")
generate = _mod.generate
WALL_HEIGHT_FT = _mod.WALL_HEIGHT_FT


# ── unit tests ────────────────────────────────────────────────

class TestConstants:
    def test_wall_height(self):
        assert WALL_HEIGHT_FT == pytest.approx(80.0 / 12.0)

    def test_slope(self):
        assert _mod.ROOF_SLOPE == pytest.approx(2.0 / 12.0)


# ── integration tests ────────────────────────────────────────

@pytest.fixture(scope="module")
def scad_content():
    return generate_scad_content(generate)


class TestGenerate2in12:
    def test_output_contains_tpath_data(self, scad_content):
        assert "t_O11_O1 = [" in scad_content
        assert "t_O1_O2 = [" in scad_content
        assert "t_O10_O11 = [" in scad_content

    def test_two_lower_sections(self, scad_content):
        matches = re.findall(r'^t_lower_O\d+_O\d+ = \[', scad_content, re.M)
        assert len(matches) == 2

    def test_eleven_middle_sections(self, scad_content):
        matches = re.findall(r'^t_O\d+_O\d+ = \[', scad_content, re.M)
        assert len(matches) == 11

    def test_output_contains_assembly(self, scad_content):
        assert "linear_extrude" in scad_content
        assert "union()" in scad_content
        assert "wall_shell(t_O11_O1, half_t);" in scad_content

    def test_output_contains_helpers(self, scad_content):
        assert "function shell_pts(" in scad_content
        assert "function seg_pts(" in scad_content
        assert "module wall_shell(" in scad_content

    def test_output_has_segment_types(self, scad_content):
        assert "[0, " in scad_content
        assert "[1, " in scad_content

    def test_three_tier_walls(self, scad_content):
        assert "lower_height" in scad_content
        assert "middle_height" in scad_content
        assert "t_full_O4 = [" in scad_content
        assert "upper_base" in scad_content
        assert "max_upper_h" in scad_content
        assert "wall_shell(t_full_O4, half_t);" in scad_content

    def test_2in12_specific_features(self, scad_content):
        assert "roof_shear" in scad_content
        assert "multmatrix(roof_shear)" in scad_content
        assert "roof_z_off" in scad_content

    def test_sloped_roof_slab(self, scad_content):
        assert "polygon(points = shell_pts(roof_outline, 0))" in scad_content
        assert "roof_thick" in scad_content

    def test_upper_wall_clipped_to_slope(self, scad_content):
        """Upper wall uses intersection to clip to sloped roof underside."""
        assert "render() intersection()" in scad_content

    def test_roof_outline_is_segments(self, scad_content):
        m = re.search(r'roof_outline = \[\n(.*?)\];',
                      scad_content, re.DOTALL)
        assert m, "roof_outline not found"
        body = m.group(1)
        data_lines = [l.strip() for l in body.split('\n')
                      if l.strip() and not l.strip().startswith('//')]
        assert len(data_lines) == 8
        for line in data_lines:
            assert line.startswith('[0,') or line.startswith('[1,'), \
                f"unexpected format: {line[:40]}"

    def test_no_difference_in_assembly(self, scad_content):
        assert "difference()" not in scad_content
