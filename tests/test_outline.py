"""Tests for floorplan/geometry.py — outline geometry."""
import math
import pytest
from shared.types import LineSeg, ArcSeg
from shared.geometry import path_polygon, poly_area
from floorplan.geometry import OutlineGeometry


# Known-good F-series coordinates (regression snapshot, F1 = origin)
_EXPECTED_F = {
    "F1":  ( 0.0000000000,  0.0000000000),
    "F2":  (-0.7413070616,  0.8282364456),
    "F3":  (-0.7413070616, 12.9115697789),
    "F4":  (-0.6902252671, 13.8338712513),
    "F5":  ( 0.3562975575, 23.2525766731),
    "F6":  ( 2.6753596051, 25.3282364456),
    "F7":  ( 7.9253596051, 25.3282364456),
    "F8":  (10.2586929384, 22.9949031123),
    "F9":  (10.4253596051, 22.8282364456),
    "F10": (25.5920262717, 22.8282364456),
    "F11": (26.5962627758, 23.5988142175),
    "F11a": (28.8500897038, 25.3282364456),
    "F11b": (29.8500897038, 25.3282364456),
    "F12": (32.1039166318, 23.5988142175),
    "F13": (35.1732501348, 12.1439056391),
    "F14": (35.2586929384, 11.4949031123),
    "F15": (35.2586929384,  2.8282364456),
    "F16": (34.0220453027,  0.6862999095),
    "F17": (29.6919182838, -1.8137000905),
    "F18": (26.4895818403, -2.6717635544),
    "F19": (25.0920262717, -2.6717635544),
    "F20": (23.0061162653, -2.5562351406),
}


class TestOutlineGeometry:
    def test_returns_outline_geometry(self, outline_geo):
        assert isinstance(outline_geo, OutlineGeometry)

    def test_22_points(self, outline_geo):
        for i in [j for j in range(21) if j != 0]:
            assert f"F{i}" in outline_geo.fp_pts
        assert "F11a" in outline_geo.fp_pts
        assert "F11b" in outline_geo.fp_pts

    def test_22_segments(self, outline_geo):
        assert len(outline_geo.outline_segs) == 22

    def test_11_radii(self, outline_geo):
        assert len(outline_geo.radii) == 11
        for key, val in outline_geo.radii.items():
            assert key.startswith("R_a")
            assert val > 0, f"{key} = {val}"

    def test_segment_connectivity(self, outline_geo):
        """Each segment's end name matches the next segment's start name."""
        segs = outline_geo.outline_segs
        for i in range(len(segs)):
            j = (i + 1) % len(segs)
            assert segs[i].end == segs[j].start, (
                f"Segment {i} end={segs[i].end} != segment {j} start={segs[j].start}"
            )

    def test_segment_types(self, outline_geo):
        """Segments alternate between lines and arcs per the design."""
        for seg in outline_geo.outline_segs:
            assert isinstance(seg, (LineSeg, ArcSeg))

    def test_outline_area(self, outline_geo):
        poly = path_polygon(outline_geo.outline_segs, outline_geo.fp_pts)
        area = poly_area(poly)
        assert abs(area - 882.07) < 0.1

    @pytest.mark.parametrize("name,expected", list(_EXPECTED_F.items()))
    def test_f_series_regression(self, outline_geo, name, expected):
        """F-series coordinates match known-good values."""
        pt = outline_geo.fp_pts[name]
        assert abs(pt[0] - expected[0]) < 1e-6, f"{name} E: {pt[0]} != {expected[0]}"
        assert abs(pt[1] - expected[1]) < 1e-6, f"{name} N: {pt[1]} != {expected[1]}"
