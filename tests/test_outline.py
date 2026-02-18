"""Tests for floorplan/geometry.py — outline geometry."""
import math
import pytest
from shared.types import LineSeg, ArcSeg
from shared.geometry import path_polygon, poly_area
from floorplan.geometry import OutlineGeometry


# Known-good F-series coordinates (regression snapshot)
_EXPECTED_F = {
    "F0":  ( 1.2407407407,  2.2134933417),
    "F1":  ( 0.5000000000,  3.0416666667),
    "F2":  ( 0.5000000000, 14.6666666667),
    "F3":  ( 0.5761162337, 15.5923323554),
    "F5":  ( 2.0313419786, 24.3811564601),
    "F6":  ( 4.3333333333, 26.3333333333),
    "F7":  ( 9.1666666667, 26.3333333333),
    "F8":  (11.5000000000, 24.0000000000),
    "F9":  (11.6666666667, 23.8333333333),
    "F10": (28.1666666667, 23.8333333333),
    "F11": (29.1709031708, 24.6039111052),
    "F12": (33.6785570268, 24.6039111052),
    "F13": (36.3082644586, 14.7897093608),
    "F14": (36.5000000000, 13.3333333333),
    "F15": (36.5000000000,  5.0000000000),
    "F16": (35.2633523643,  2.8580634639),
    "F17": (30.9332253454,  0.3580634639),
    "F18": (27.5000000000, -0.5000000000),
    "F19": (26.5000000000, -0.5000000000),
    "F20": (24.5280501561, -0.3901070016),
}


class TestOutlineGeometry:
    def test_returns_outline_geometry(self, outline_geo):
        assert isinstance(outline_geo, OutlineGeometry)

    def test_20_points(self, outline_geo):
        for i in [j for j in range(21) if j != 4]:
            assert f"F{i}" in outline_geo.fp_pts

    def test_20_segments(self, outline_geo):
        assert len(outline_geo.outline_segs) == 20

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
        assert abs(area - 837.77) < 0.1

    @pytest.mark.parametrize("name,expected", list(_EXPECTED_F.items()))
    def test_f_series_regression(self, outline_geo, name, expected):
        """F-series coordinates match known-good values."""
        pt = outline_geo.fp_pts[name]
        assert abs(pt[0] - expected[0]) < 1e-6, f"{name} E: {pt[0]} != {expected[0]}"
        assert abs(pt[1] - expected[1]) < 1e-6, f"{name} N: {pt[1]} != {expected[1]}"
