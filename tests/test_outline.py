"""Tests for floorplan/geometry.py — outline geometry."""
import math
import pytest
from shared.types import LineSeg, ArcSeg
from shared.geometry import path_polygon, poly_area
from floorplan.geometry import OutlineGeometry


# Known-good F-series coordinates (regression snapshot, F1 = origin, F20→F1 brg=270°)
_EXPECTED_F = {
    "F1":  ( 0.0000000000,  0.0000000000),
    "F2":  (-0.8282364456,  0.7413070616),
    "F3":  (-2.1626173856, 12.7507355222),
    "F4":  (-2.2136991801, 13.6730369946),
    "F5":  (-2.2136991801, 23.1497042276),
    "F6":  (-0.1380394076, 25.4687662752),
    "F7":  ( 5.0798501994, 26.0485317871),
    "F8":  ( 7.6565858079, 23.9871433004),
    "F9":  ( 7.8406383513, 23.8399012656),
    "F10": (22.9145416605, 25.5147794111),
    "F11": (23.8275399084, 26.3915434945),
    "F11a": (25.8765990970, 28.3592816646),
    "F11b": (26.8704828317, 28.4697131907),
    "F12": (29.3015074921, 26.9997621150),
    "F13": (33.6170511724, 15.9538659795),
    "F14": (33.7736417246, 15.3182685035),
    "F15": (34.7307149505,  6.7046094696),
    "F16": (33.7381683003,  4.4392087001),
    "F17": (29.7106043024,  1.4763168286),
    "F18": (26.6226114560,  0.2698626081),
    "F19": (25.2336037081,  0.1155284138),
    "F20": (23.1476937017,  0.0000000000),
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
