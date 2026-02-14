"""Tests for floorplan/geometry.py — outline geometry."""
import math
import pytest
from shared.types import LineSeg, ArcSeg
from shared.geometry import path_polygon, poly_area
from floorplan.geometry import OutlineGeometry


# Known-good F-series coordinates (regression snapshot)
_EXPECTED_F = {
    "F0":  ( 1.3333333333,  0.5000000000),
    "F1":  ( 0.5000000000,  1.3333333333),
    "F2":  ( 0.5000000000, 18.0000000000),
    "F3":  ( 1.0821917808, 19.5525114155),
    "F4":  ( 1.5000000000, 20.6666666667),
    "F5":  ( 1.5000000000, 24.0000000000),
    "F6":  ( 3.8333333333, 26.3333333333),
    "F7":  ( 9.1666666667, 26.3333333333),
    "F8":  (11.5000000000, 24.0000000000),
    "F9":  (11.6666666667, 23.8333333333),
    "F10": (29.1433271296, 23.8333333333),
    "F11": (29.3099937962, 24.0000000000),
    "F12": (33.8971540576, 24.6039111052),
    "F13": (36.0971159410, 16.3935415818),
    "F14": (36.5000000000, 13.3333333333),
    "F15": (36.5000000000,  5.0000000000),
    "F16": (35.2633523643,  2.8580634639),
    "F17": (30.9332253454,  0.3580634639),
    "F18": (27.7308889019, -0.5000000000),
    "F19": (21.5000000000, -0.5000000000),
    "F20": (19.9474885845,  0.0821917808),
    "F21": (18.8333333333,  0.5000000000),
}


class TestOutlineGeometry:
    def test_returns_outline_geometry(self, outline_geo):
        assert isinstance(outline_geo, OutlineGeometry)

    def test_22_points(self, outline_geo):
        for i in range(22):
            assert f"F{i}" in outline_geo.fp_pts

    def test_22_segments(self, outline_geo):
        assert len(outline_geo.outline_segs) == 22

    def test_13_radii(self, outline_geo):
        assert len(outline_geo.radii) == 13
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
        assert abs(area - 858.11) < 0.1

    @pytest.mark.parametrize("name,expected", list(_EXPECTED_F.items()))
    def test_f_series_regression(self, outline_geo, name, expected):
        """F-series coordinates match known-good values."""
        pt = outline_geo.fp_pts[name]
        assert abs(pt[0] - expected[0]) < 1e-6, f"{name} E: {pt[0]} != {expected[0]}"
        assert abs(pt[1] - expected[1]) < 1e-6, f"{name} N: {pt[1]} != {expected[1]}"
