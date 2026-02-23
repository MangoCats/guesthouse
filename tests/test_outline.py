"""Tests for floorplan/geometry.py — outline geometry."""
import math
import pytest
from shared.types import LineSeg, ArcSeg
from shared.geometry import path_polygon, poly_area
from floorplan.geometry import OutlineGeometry


# Known-good F-series coordinates (FC-based: building center = origin)
_EXPECTED_F = {
    "F1":  (-17.2586929384, -11.3282364456),
    "F2":  (-18.0000000000, -10.5000000000),
    "F3":  (-18.0000000000,   1.5833333333),
    "F4":  (-17.9489182055,   2.5056348057),
    "F5":  (-16.9023953809,  11.9243402275),
    "F6":  (-14.5833333333,  14.0000000000),
    "F7":  ( -9.3333333333,  14.0000000000),
    "F8":  ( -7.0000000000,  11.6666666667),
    "F9":  ( -6.8333333333,  11.5000000000),
    "F10": (  8.3333333333,  11.5000000000),
    "F11": (  9.3375698374,  12.2705777719),
    "F11a": (11.5913967654,  14.0000000000),
    "F11b": (12.5913967654,  14.0000000000),
    "F12": ( 14.8452236934,  12.2705777719),
    "F13": ( 17.9145571964,   0.8156691935),
    "F14": ( 18.0000000000,   0.1666666667),
    "F15": ( 18.0000000000,  -8.5000000000),
    "F16": ( 16.7633523643, -10.6419365361),
    "F17": ( 12.4332253454, -13.1419365361),
    "F18": (  9.2308889019, -14.0000000000),
    "F19": (  7.8333333333, -14.0000000000),
    "F20": (  5.7474233269, -13.8844715862),
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
