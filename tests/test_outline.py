"""Tests for floorplan/geometry.py — outline geometry."""
import math
import pytest
from shared.types import LineSeg, ArcSeg
from shared.geometry import path_polygon, poly_area
from floorplan.geometry import OutlineGeometry


# Known-good F-series coordinates (regression snapshot)
_EXPECTED_F = {
    "F1":  ( 1.2413070616,  2.1717635544),
    "F2":  ( 0.5000000000,  3.0000000000),
    "F3":  ( 0.5000000000, 15.0833333333),
    "F4":  ( 0.6223253065, 17.2919638548),
    "F5":  ( 1.5259226813, 25.4243402275),
    "F6":  ( 3.8449847288, 27.5000000000),
    "F7":  ( 9.1666666667, 27.5000000000),
    "F8":  (11.5000000000, 25.1666666667),
    "F9":  (11.6666666667, 25.0000000000),
    "F10": (26.8333333333, 25.0000000000),
    "F11": (27.8375698374, 25.7705777719),
    "F11a": (30.0913967654, 27.5000000000),
    "F11b": (31.0913967654, 27.5000000000),
    "F12": (33.3452236934, 25.7705777719),
    "F13": (36.4145571964, 14.3156691935),
    "F14": (36.5000000000, 13.6666666667),
    "F15": (36.5000000000,  5.0000000000),
    "F16": (35.2633523643,  2.8580634639),
    "F17": (30.9332253454,  0.3580634639),
    "F18": (27.7308889019, -0.5000000000),
    "F19": (26.3333333333, -0.5000000000),
    "F20": (24.2474233269, -0.3844715862),
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
        assert abs(area - 882.88) < 0.1

    @pytest.mark.parametrize("name,expected", list(_EXPECTED_F.items()))
    def test_f_series_regression(self, outline_geo, name, expected):
        """F-series coordinates match known-good values."""
        pt = outline_geo.fp_pts[name]
        assert abs(pt[0] - expected[0]) < 1e-6, f"{name} E: {pt[0]} != {expected[0]}"
        assert abs(pt[1] - expected[1]) < 1e-6, f"{name} N: {pt[1]} != {expected[1]}"
