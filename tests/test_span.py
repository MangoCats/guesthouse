"""Tests for span/gen_span.py, span/gen_span_min.py, span/gen_span_minmax.py."""
import math
import pytest
from conftest import _import_from

# Import generator modules (for _generate_svg / _compute_spans)
_span = _import_from("span", "gen_span")
_smin = _import_from("span", "gen_span_min")
_smm = _import_from("span", "gen_span_minmax")

# Import shared helpers from _common
from span._common import (
    rot_pt, rot_poly, seg_vert_isect,
    extract_iw_centerlines, max_span_at_angle,
    compute_rotation_data,
)


# ===================================================================
# Pure helper unit tests — via span._common
# ===================================================================

class TestRotPt:
    def test_identity(self):
        p = (3.0, 4.0)
        result = rot_pt(p, 0.0, 0.0, 1.0, 0.0)
        assert result[0] == pytest.approx(3.0)
        assert result[1] == pytest.approx(4.0)

    def test_90_degrees(self):
        # (1, 0) rotated 90° CCW around origin → (0, 1)
        result = rot_pt((1.0, 0.0), 0.0, 0.0, 0.0, 1.0)
        assert result[0] == pytest.approx(0.0, abs=1e-12)
        assert result[1] == pytest.approx(1.0, abs=1e-12)

    def test_around_offset_center(self):
        # (5, 3) rotated 180° around (3, 3) → (1, 3)
        result = rot_pt((5.0, 3.0), 3.0, 3.0, -1.0, 0.0)
        assert result[0] == pytest.approx(1.0, abs=1e-12)
        assert result[1] == pytest.approx(3.0, abs=1e-12)


class TestRotPoly:
    def test_rotates_all_points(self):
        poly = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]
        result = rot_poly(poly, 0.0, 0.0, 1.0, 0.0)  # identity
        assert len(result) == 4
        assert result[0][0] == pytest.approx(0.0)

    def test_empty_polygon(self):
        assert rot_poly([], 0.0, 0.0, 1.0, 0.0) == []


class TestSegVertIsect:
    def test_horizontal_crossing(self):
        # Segment (0,0)-(10,0), vertical at e=5 → n=0
        n = seg_vert_isect((0.0, 0.0), (10.0, 0.0), 5.0)
        assert n == pytest.approx(0.0)

    def test_diagonal_crossing(self):
        # Segment (0,0)-(10,10), vertical at e=5 → n=5
        n = seg_vert_isect((0.0, 0.0), (10.0, 10.0), 5.0)
        assert n == pytest.approx(5.0)

    def test_vertical_segment_returns_none(self):
        # Vertical segment (5,0)-(5,10) has de=0
        assert seg_vert_isect((5.0, 0.0), (5.0, 10.0), 5.0) is None

    def test_outside_segment_returns_none(self):
        # Segment (0,0)-(10,0), e=20 is outside
        assert seg_vert_isect((0.0, 0.0), (10.0, 0.0), 20.0) is None


# ===================================================================
# Integration tests using real geometry
# ===================================================================

class TestExtractIwCenterlines:
    def test_returns_three(self, span_geometry):
        _, _, _, inner_poly, layout, _ = span_geometry
        cls = extract_iw_centerlines(layout)
        assert len(cls) == 3

    def test_each_is_two_points(self, span_geometry):
        _, _, _, inner_poly, layout, _ = span_geometry
        cls = extract_iw_centerlines(layout)
        for cl in cls:
            assert len(cl) == 2
            assert len(cl[0]) == 2  # (E, N)


class TestComputeSpans:
    @pytest.fixture(scope="class")
    def spans_result(self, span_geometry):
        _, _, _, inner_poly, layout, _ = span_geometry
        return _span._compute_spans(inner_poly, layout)

    def test_four_equal_length_lists(self, spans_result):
        eastings, spans, south, north = spans_result
        assert len(eastings) == len(spans) == len(south) == len(north)

    def test_many_eastings(self, spans_result):
        assert len(spans_result[0]) > 100

    def test_max_span_in_range(self, spans_result):
        _, spans, _, _ = spans_result
        ms = max(spans)
        assert 20 < ms < 35  # building is ~27' N-S

    def test_no_negative_spans(self, spans_result):
        _, spans, south, north = spans_result
        assert all(s >= 0 for s in spans)
        assert all(s >= 0 for s in south)
        assert all(s >= 0 for s in north)

    def test_south_plus_north_leq_total(self, spans_result):
        _, spans, south, north = spans_result
        for s, ss, ns in zip(spans, south, north):
            if s > 0:
                assert ss <= s + 1e-9
                assert ns <= s + 1e-9


class TestMaxSpanAtAngle:
    def test_zero_rotation_positive(self, span_geometry):
        _, _, _, inner_poly, layout, _ = span_geometry
        iw_cls = extract_iw_centerlines(layout)
        cx = sum(p[0] for p in inner_poly) / len(inner_poly)
        cy = sum(p[1] for p in inner_poly) / len(inner_poly)
        ms = max_span_at_angle(inner_poly, iw_cls, 0.0, cx, cy)
        assert ms > 0

    def test_different_at_90(self, span_geometry):
        _, _, _, inner_poly, layout, _ = span_geometry
        iw_cls = extract_iw_centerlines(layout)
        cx = sum(p[0] for p in inner_poly) / len(inner_poly)
        cy = sum(p[1] for p in inner_poly) / len(inner_poly)
        ms0 = max_span_at_angle(inner_poly, iw_cls, 0.0, cx, cy)
        ms90 = max_span_at_angle(inner_poly, iw_cls, 90.0, cx, cy)
        assert ms0 != pytest.approx(ms90, abs=0.5)


class TestComputeRotationData:
    @pytest.fixture(scope="class")
    def rot_data(self, span_geometry):
        _, _, outer_poly, inner_poly, layout, roof_poly = span_geometry
        iw_cls = extract_iw_centerlines(layout)
        cx = sum(p[0] for p in inner_poly) / len(inner_poly)
        cy = sum(p[1] for p in inner_poly) / len(inner_poly)
        return compute_rotation_data(
            45.0, outer_poly, inner_poly, iw_cls, cx, cy, roof_poly)

    def test_expected_keys(self, rot_data):
        expected = {"angle", "eastings", "spans", "s_spans", "n_spans",
                    "r_outer", "r_inner", "r_cls", "max_span", "max_e",
                    "max_roof_span", "e_min", "e_max", "n_min", "n_max",
                    "roof_spans"}
        assert expected.issubset(rot_data.keys())

    def test_max_span_positive(self, rot_data):
        assert rot_data["max_span"] > 0

    def test_rotated_polys_same_length(self, rot_data, span_geometry):
        _, _, outer_poly, inner_poly, _, _ = span_geometry
        assert len(rot_data["r_outer"]) == len(outer_poly)
        assert len(rot_data["r_inner"]) == len(inner_poly)


# ===================================================================
# SVG generation smoke tests
# ===================================================================

class TestSpanSvg:
    @pytest.fixture(scope="class")
    def svg(self, span_geometry):
        pts, _, outer_poly, inner_poly, layout, roof_poly = span_geometry
        return _span._generate_svg(pts, outer_poly, inner_poly, layout, roof_poly)

    def test_svg_envelope(self, svg):
        assert svg.startswith("<svg")
        assert svg.rstrip().endswith("</svg>")

    def test_contains_title(self, svg):
        assert "Interior Span" in svg


class TestSpanMinSvg:
    def test_svg_envelope(self, span_geometry):
        pts, _, outer_poly, inner_poly, layout, roof_poly = span_geometry
        svg = _smin._generate_svg(pts, outer_poly, inner_poly, layout, roof_poly)
        assert svg.startswith("<svg")
        assert svg.rstrip().endswith("</svg>")
        assert "Minimum-Span" in svg


class TestSpanMinmaxSvg:
    def test_svg_envelope(self, span_geometry):
        pts, _, outer_poly, inner_poly, layout, roof_poly = span_geometry
        svg = _smm._generate_svg(pts, outer_poly, inner_poly, layout, roof_poly)
        assert svg.startswith("<svg")
        assert svg.rstrip().endswith("</svg>")
        assert "Rotation" in svg
