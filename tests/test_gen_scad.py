"""Tests for scad/gen_flat_roof.py OpenSCAD generation (T-path approach)."""
import math
import re
import pytest
from conftest import _import_from, generate_scad_content

from floorplan.gen_floorplan import build_floorplan_data
from floorplan.roof import compute_roof_geometry, roof_segments, roof_polyline
from shared.types import LineSeg, ArcSeg

_mod = _import_from("scad", "gen_flat_roof")
_scad_seg = _mod._scad_seg
_rev_elem = _mod._rev_elem
_seg_to_elem = _mod._seg_to_elem
_fmt_ft_in = _mod._fmt_ft_in
_seg_comment = _mod._seg_comment
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


class TestFmtFtIn:
    def test_exact_feet(self):
        result = _fmt_ft_in(2.0)
        assert result == " 2'   0.0000\""

    def test_half_foot(self):
        result = _fmt_ft_in(0.5)
        assert result == " 0'   6.0000\""

    def test_custom_width(self):
        result = _fmt_ft_in(0.25, in_width=7)
        assert result == " 0'  3.0000\""


class TestSegComment:
    def test_line_bearing_north(self):
        elem = ("line", 0.0, 0.0, 0.0, 1.0)  # due north
        comment = _seg_comment(elem)
        assert "0.0000deg" in comment
        assert "//" in comment

    def test_line_bearing_east(self):
        elem = ("line", 0.0, 0.0, 1.0, 0.0)  # due east
        comment = _seg_comment(elem)
        assert "90.0000deg" in comment

    def test_arc_cw_label(self):
        elem = ("arc", 0.0, 0.0, 1.0, 90.0, 0.0)  # a2 < a1 → CW
        comment = _seg_comment(elem)
        assert "CW" in comment

    def test_arc_ccw_label(self):
        elem = ("arc", 0.0, 0.0, 1.0, 0.0, 90.0)  # a2 > a1 → CCW
        comment = _seg_comment(elem)
        assert "CCW" in comment


class TestSegToElem:
    def test_line_seg(self):
        pts = {"A": (1.0, 2.0), "B": (3.0, 4.0)}
        seg = LineSeg("A", "B")
        result = _seg_to_elem(seg, pts)
        assert result == ("line", 1.0, 2.0, 3.0, 4.0)

    def test_arc_seg_cw(self):
        pts = {"A": (1.0, 0.0), "B": (0.0, -1.0), "C": (0.0, 0.0)}
        seg = ArcSeg("A", "B", "C", 1.0, "CW", 10)
        result = _seg_to_elem(seg, pts)
        assert result[0] == "arc"
        assert result[1] == pytest.approx(0.0)  # cx
        assert result[2] == pytest.approx(0.0)  # cy
        assert result[3] == pytest.approx(1.0)  # radius
        assert result[4] > result[5]             # CW: a1 > a2

    def test_arc_seg_ccw(self):
        pts = {"A": (1.0, 0.0), "B": (0.0, 1.0), "C": (0.0, 0.0)}
        seg = ArcSeg("A", "B", "C", 1.0, "CCW", 10)
        result = _seg_to_elem(seg, pts)
        assert result[0] == "arc"
        assert result[5] > result[4]             # CCW: a2 > a1


class TestWallHeight:
    def test_conversion(self):
        assert WALL_HEIGHT_FT == pytest.approx(80.0 / 12.0)


# ── roof_segments tests ──────────────────────────────────────

@pytest.fixture(scope="module")
def roof_geo():
    fp = build_floorplan_data()
    return compute_roof_geometry(fp.pts, fp.radii)


def _elem_start(elem):
    """Start point of an element."""
    if elem[0] == "line":
        return (elem[1], elem[2])
    _, cx, cy, r, a1, _a2 = elem
    return (cx + r * math.cos(math.radians(a1)),
            cy + r * math.sin(math.radians(a1)))


def _elem_end(elem):
    """End point of an element."""
    if elem[0] == "line":
        return (elem[3], elem[4])
    _, cx, cy, r, _a1, a2 = elem
    return (cx + r * math.cos(math.radians(a2)),
            cy + r * math.sin(math.radians(a2)))


class TestRoofSegments:
    def test_segment_count(self, roof_geo):
        segs = roof_segments(roof_geo)
        assert len(segs) == 8

    def test_segment_types(self, roof_geo):
        segs = roof_segments(roof_geo)
        lines = [s for s in segs if s[0] == "line"]
        arcs = [s for s in segs if s[0] == "arc"]
        assert len(lines) == 6
        assert len(arcs) == 2

    def test_arcs_are_cw(self, roof_geo):
        segs = roof_segments(roof_geo)
        for seg in segs:
            if seg[0] == "arc":
                a1, a2 = seg[4], seg[5]
                assert a2 < a1, f"arc a2={a2} should be < a1={a1} for CW"

    def test_arc_radii_match(self, roof_geo):
        segs = roof_segments(roof_geo)
        arcs = [s for s in segs if s[0] == "arc"]
        assert arcs[0][3] == pytest.approx(roof_geo.r3_radius)
        assert arcs[1][3] == pytest.approx(roof_geo.r4_radius)

    def test_arc_centers_match(self, roof_geo):
        segs = roof_segments(roof_geo)
        arcs = [s for s in segs if s[0] == "arc"]
        assert arcs[0][1] == pytest.approx(roof_geo.r3_center[0])
        assert arcs[0][2] == pytest.approx(roof_geo.r3_center[1])
        assert arcs[1][1] == pytest.approx(roof_geo.r4_center[0])
        assert arcs[1][2] == pytest.approx(roof_geo.r4_center[1])

    def test_chain_continuity(self, roof_geo):
        """End of each segment matches start of next (closed loop)."""
        segs = roof_segments(roof_geo)
        for i in range(len(segs)):
            ep = _elem_end(segs[i])
            sp = _elem_start(segs[(i + 1) % len(segs)])
            assert ep[0] == pytest.approx(sp[0], abs=1e-9), \
                f"seg {i}→{i+1} E gap: {ep[0]} vs {sp[0]}"
            assert ep[1] == pytest.approx(sp[1], abs=1e-9), \
                f"seg {i}→{i+1} N gap: {ep[1]} vs {sp[1]}"

    def test_agrees_with_polyline(self, roof_geo):
        """Sampled segment points should match roof_polyline output."""
        segs = roof_segments(roof_geo)
        poly = roof_polyline(roof_geo, n_arc=30)

        # Sample points from segments at the same resolution
        seg_pts = []
        for seg in segs:
            if seg[0] == "line":
                seg_pts.append((seg[1], seg[2]))
            else:
                _, cx, cy, r, a1, a2 = seg
                n = 30
                for i in range(n):
                    a = math.radians(a1 + (a2 - a1) * i / n)
                    seg_pts.append((cx + r * math.cos(a),
                                    cy + r * math.sin(a)))
        # Close: add last segment's endpoint
        last = segs[-1]
        seg_pts.append((last[3], last[4]))  # line endpoint

        # The polyline has R1..R7 (no closing R1); seg_pts starts at R1
        # and includes closing R1 at end. Compare the shared prefix.
        for i, (sp, pp) in enumerate(zip(seg_pts, poly)):
            assert sp[0] == pytest.approx(pp[0], abs=1e-6), \
                f"point {i} E: {sp[0]} vs {pp[0]}"
            assert sp[1] == pytest.approx(pp[1], abs=1e-6), \
                f"point {i} N: {sp[1]} vs {pp[1]}"


# ── T-path structural tests ─────────────────────────────────

class TestTpathStructure:
    @pytest.fixture(scope="class")
    def scad_content(self):
        return generate_scad_content(generate)

    def test_two_lower_sections(self, scad_content):
        matches = re.findall(r'^t_lower_O\d+_O\d+ = \[', scad_content, re.M)
        assert len(matches) == 2

    def test_twelve_middle_sections(self, scad_content):
        matches = re.findall(r'^t_O\w+_O\w+ = \[', scad_content, re.M)
        assert len(matches) == 12

    def test_roof_outline_is_segments(self, scad_content):
        """roof_outline uses [0,...] / [1,...] segment format, not [x,y]."""
        m = re.search(r'roof_outline = \[\n(.*?)\];',
                      scad_content, re.DOTALL)
        assert m, "roof_outline not found"
        body = m.group(1)
        # Every data line should start with [0, or [1,
        data_lines = [l.strip() for l in body.split('\n')
                      if l.strip() and not l.strip().startswith('//')]
        for line in data_lines:
            assert line.startswith('[0,') or line.startswith('[1,'), \
                f"unexpected format: {line[:40]}"

    def test_roof_outline_segment_count(self, scad_content):
        m = re.search(r'roof_outline = \[\n(.*?)\];',
                      scad_content, re.DOTALL)
        body = m.group(1)
        data_lines = [l.strip() for l in body.split('\n')
                      if l.strip() and not l.strip().startswith('//')]
        assert len(data_lines) == 8


# ── integration test ──────────────────────────────────────────

class TestGenerate:
    @pytest.fixture(scope="class")
    def scad_content(self):
        return generate_scad_content(generate)

    def test_output_contains_tpath_data(self, scad_content):
        # 12 wall sections, each with a T-path variable
        assert "t_O11_O1 = [" in scad_content
        assert "t_O1_O2 = [" in scad_content
        assert "t_O10_O11 = [" in scad_content

    def test_output_contains_assembly(self, scad_content):
        assert "linear_extrude" in scad_content
        assert "union()" in scad_content
        assert "wall_shell(t_O11_O1, half_t);" in scad_content

    def test_output_contains_helpers(self, scad_content):
        assert "function shell_pts(" in scad_content
        assert "function seg_pts(" in scad_content
        assert "function line_pts(" in scad_content
        assert "function arc_off_pts(" in scad_content
        assert "function tail(" in scad_content
        assert "module wall_shell(" in scad_content

    def test_output_has_segment_types(self, scad_content):
        # Lines start with [0, and arcs with [1,
        assert "[0, " in scad_content
        assert "[1, " in scad_content

    def test_no_difference_in_assembly(self, scad_content):
        # T-path approach uses polygon with paths, no difference()
        assert "difference()" not in scad_content

    def test_output_three_tier_walls(self, scad_content):
        assert "lower_height" in scad_content
        assert "middle_height" in scad_content
        assert "t_full_O4 = [" in scad_content
        assert "upper_base" in scad_content
        assert "upper_height" in scad_content
        assert "wall_shell(t_full_O4, half_t);" in scad_content

    def test_output_wedge_roof(self, scad_content):
        assert "translate([0, 0, upper_base + upper_height])" in scad_content
        assert "polygon(points = shell_pts(roof_outline, 0))" in scad_content
        assert "max_roof_thick" in scad_content
        assert "roof_slope" in scad_content
        assert "roof_z_base" in scad_content
        assert "render() intersection()" in scad_content
