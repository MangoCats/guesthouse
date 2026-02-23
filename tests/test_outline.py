"""Tests for floorplan/geometry.py — outline geometry."""
import math
import pytest
from shared.types import LineSeg, ArcSeg
from shared.geometry import path_polygon, poly_area
from floorplan.geometry import OutlineGeometry


# Known-good F-series coordinates (FC-based: building center = origin,
# rotated by COORD_ROTATION = arctan(14/85))
_EXPECTED_F = {
    "F1":  (-17.7822585106, -14.1308018679),
    "F2":  (-18.6155918439, -13.2974685346),
    "F5":  (-18.6155918439,   9.0189023597),
    "F6":  (-16.5399320714,  11.3379644073),
    "F7":  (-11.3220424644,  11.9177299192),
    "F8":  ( -8.7453068560,   9.8563414324),
    "F9":  ( -8.5612543125,   9.7090993977),
    "F10": (  6.5126489967,  11.3839775431),
    "F11": (  7.4256472445,  12.2607416266),
    "F11a": ( 9.4747064332,  14.2284797967),
    "F11b": (10.4685901679,  14.3389113227),
    "F12": ( 12.8996148283,  12.8689602470),
    "F13": ( 17.2151585085,   1.8230641115),
    "F14": ( 17.3717490607,   1.1874666355),
    "F15": ( 18.3288222867,  -7.4261923983),
    "F16": ( 17.3362756365,  -9.6915931678),
    "F17": ( 13.3087116385, -12.6544850393),
    "F18": ( 10.2207187921, -13.8609392599),
    "F19": (  8.8317110442, -14.0152734541),
    "F20": (  6.7458010378, -14.1308018679),
}


class TestOutlineGeometry:
    def test_returns_outline_geometry(self, outline_geo):
        assert isinstance(outline_geo, OutlineGeometry)

    def test_20_points(self, outline_geo):
        for i in [j for j in range(21) if j not in (0, 3, 4)]:
            assert f"F{i}" in outline_geo.fp_pts
        assert "F11a" in outline_geo.fp_pts
        assert "F11b" in outline_geo.fp_pts

    def test_20_segments(self, outline_geo):
        assert len(outline_geo.outline_segs) == 20

    def test_10_radii(self, outline_geo):
        assert len(outline_geo.radii) == 10
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
        assert abs(area - 891.73) < 0.1

    @pytest.mark.parametrize("name,expected", list(_EXPECTED_F.items()))
    def test_f_series_regression(self, outline_geo, name, expected):
        """F-series coordinates match known-good values."""
        pt = outline_geo.fp_pts[name]
        assert abs(pt[0] - expected[0]) < 1e-6, f"{name} E: {pt[0]} != {expected[0]}"
        assert abs(pt[1] - expected[1]) < 1e-6, f"{name} N: {pt[1]} != {expected[1]}"
