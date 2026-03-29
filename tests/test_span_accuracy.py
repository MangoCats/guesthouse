"""Span accuracy golden-value tests — numeric correctness of span computations."""
import pytest
from span._common import (
    max_span_at_angle, compute_rotation_data, find_min_span_angle,
    extract_iw_centerlines,
)


def _centroid(poly):
    """Compute centroid of a polygon."""
    n = len(poly)
    return sum(p[0] for p in poly) / n, sum(p[1] for p in poly) / n


# ---------------------------------------------------------------------------
# Golden value tests
# ---------------------------------------------------------------------------

class TestSpanGoldenValues:
    """Verify span magnitudes fall within expected ranges for the seed building."""

    @pytest.fixture(scope="class")
    def geom(self, span_geometry):
        _pts, _outline_segs, outline_poly, inner_poly, layout, roof_poly = span_geometry
        iw_cls = extract_iw_centerlines(layout)
        cx, cy = _centroid(inner_poly)
        return inner_poly, iw_cls, outline_poly, roof_poly, cx, cy

    def test_span_at_0_degrees_in_expected_range(self, geom):
        """N-S span at 0° rotation should be within the known building width."""
        inner_poly, iw_cls, _op, _rp, cx, cy = geom
        ms = max_span_at_angle(inner_poly, iw_cls, 0.0, cx, cy)
        assert 25.0 <= ms <= 35.0, f"Expected 25–35 ft N-S span at 0°, got {ms:.2f}"

    def test_span_at_90_degrees_in_expected_range(self, geom):
        """E-W span at 90° rotation should be in the expected range."""
        inner_poly, iw_cls, _op, _rp, cx, cy = geom
        ms = max_span_at_angle(inner_poly, iw_cls, 90.0, cx, cy)
        assert 25.0 <= ms <= 65.0, f"Expected 25–65 ft E-W span at 90°, got {ms:.2f}"

    def test_find_min_span_angle_returns_valid_range(self, geom):
        """Minimum span angle should be within [−10°, 90°] with span in [20, 40]."""
        inner_poly, iw_cls, _op, _rp, cx, cy = geom
        angle, span = find_min_span_angle(inner_poly, iw_cls, cx, cy)
        assert -10.0 <= angle <= 90.0, f"Min-span angle {angle:.2f}° out of expected range"
        assert 20.0 <= span <= 40.0, f"Min span {span:.2f} ft out of expected range"

    def test_min_span_leq_span_at_zero(self, geom):
        """Minimum span must be ≤ span at 0° (min-finding worked)."""
        inner_poly, iw_cls, _op, _rp, cx, cy = geom
        _angle, min_span = find_min_span_angle(inner_poly, iw_cls, cx, cy)
        span_at_0 = max_span_at_angle(inner_poly, iw_cls, 0.0, cx, cy)
        assert min_span <= span_at_0 + 1e-6

    def test_rotation_data_max_span_matches_direct_computation(self, geom):
        """compute_rotation_data max_span agrees with max_span_at_angle (within 1")."""
        inner_poly, iw_cls, outline_poly, _rp, cx, cy = geom
        data = compute_rotation_data(0.0, outline_poly, inner_poly, iw_cls, cx, cy)
        direct = max_span_at_angle(inner_poly, iw_cls, 0.0, cx, cy)
        assert abs(data["max_span"] - direct) <= 1.0 / 12.0  # within 1 inch

    def test_span_at_0_and_180_degrees_equal(self, geom):
        """Span at 0° and 180° should be equal (building mirrored same polygon)."""
        inner_poly, iw_cls, _op, _rp, cx, cy = geom
        ms0 = max_span_at_angle(inner_poly, iw_cls, 0.0, cx, cy)
        ms180 = max_span_at_angle(inner_poly, iw_cls, 180.0, cx, cy)
        assert ms0 == pytest.approx(ms180, abs=1.0 / 12.0)  # within 1 inch


# ---------------------------------------------------------------------------
# Span decomposition consistency
# ---------------------------------------------------------------------------

class TestSpanDecomposition:
    """Verify that span sub-components are internally consistent."""

    @pytest.fixture(scope="class")
    def rot_data(self, span_geometry):
        _pts, _outline_segs, outline_poly, inner_poly, layout, roof_poly = span_geometry
        iw_cls = extract_iw_centerlines(layout)
        cx, cy = _centroid(inner_poly)
        return compute_rotation_data(0.0, outline_poly, inner_poly, iw_cls, cx, cy)

    def test_parallel_list_lengths(self, rot_data):
        n = len(rot_data["eastings"])
        assert n == len(rot_data["spans"])
        assert n == len(rot_data["s_spans"])
        assert n == len(rot_data["n_spans"])

    def test_spans_non_negative(self, rot_data):
        assert all(s >= 0 for s in rot_data["spans"])

    def test_sub_spans_non_negative(self, rot_data):
        assert all(s >= 0 for s in rot_data["s_spans"])
        assert all(s >= 0 for s in rot_data["n_spans"])

    def test_sub_spans_sum_leq_total(self, rot_data):
        """South + north sub-spans must not exceed total span at each easting."""
        for total, s, n in zip(rot_data["spans"], rot_data["s_spans"], rot_data["n_spans"]):
            assert s + n <= total + 1e-6

    def test_max_e_within_bounds(self, rot_data):
        """Easting of max span must lie within the sampled easting range."""
        assert rot_data["e_min"] <= rot_data["max_e"] <= rot_data["e_max"]

    def test_expected_keys_present(self, rot_data):
        required = {
            "angle", "eastings", "spans", "s_spans", "n_spans",
            "max_span", "max_e", "e_min", "e_max",
        }
        assert required.issubset(rot_data.keys())
