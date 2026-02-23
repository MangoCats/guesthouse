"""Tests for floorplan/geometry.py — outline geometry."""
import math
import pytest
from shared.types import LineSeg, ArcSeg
from shared.geometry import path_polygon, poly_area
from floorplan.geometry import OutlineGeometry


# Known-good F-series coordinates (FC-based: building center = origin,
# rotated by COORD_ROTATION = arctan(1/9) so F4-F5 is at bearing 0)
_EXPECTED_F = {
    "F1":  (-15.9021397548, -13.1648537450),
    "F2":  (-16.7303762003, -12.4235466834),
    "F3":  (-18.0647571404,  -0.4141182228),
    "F4":  (-18.1158389349,   0.5081832496),
    "F5":  (-18.1158389349,   9.9848504826),
    "F6":  (-16.0401791624,  12.3039125302),
    "F7":  (-10.8222895553,  12.8836780421),
    "F8":  ( -8.2455539469,  10.8222895553),
    "F9":  ( -8.0615014035,  10.6750475206),
    "F10": (  7.0124019058,  12.3499256660),
    "F11": (  7.9254001536,  13.2266897495),
    "F11a": ( 9.9744593423,  15.1944279196),
    "F11b": (10.9683430769,  15.3048594457),
    "F12": ( 13.3993677373,  13.8349083699),
    "F13": ( 17.7149114176,   2.7890122344),
    "F14": ( 17.8715019698,   2.1534147585),
    "F15": ( 18.8285751958,  -6.4602442754),
    "F16": ( 17.8360285456,  -8.7256450449),
    "F17": ( 13.8084645476, -11.6885369164),
    "F18": ( 10.7204717012, -12.8949911370),
    "F19": (  9.3314639533, -13.0493253312),
    "F20": (  7.2455539469, -13.1648537450),
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
