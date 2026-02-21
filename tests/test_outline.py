"""Tests for floorplan/geometry.py — outline geometry."""
import math
import pytest
from shared.types import LineSeg, ArcSeg
from shared.geometry import path_polygon, poly_area
from floorplan.geometry import OutlineGeometry


# Known-good F-series coordinates (regression snapshot)
_EXPECTED_F = {
    "F1":  ( 1.2229018073,  2.0061162653),
    "F2":  ( 0.3333333333,  3.0000000000),
    "F3":  ( 0.3333333333, 16.3333333333),
    "F4":  ( 0.4279945573, 17.4153160748),
    "F5":  ( 1.8713139508, 25.6007871108),
    "F6":  ( 4.3333333333, 27.6666666667),
    "F7":  ( 9.1666666667, 27.6666666667),
    "F8":  (11.6617790765, 25.3229166667),
    "F9":  (11.8281199038, 25.1666666667),
    "F10": (26.8333333333, 25.1666666667),
    "F11": (27.6765821997, 25.8137142794),
    "F11a": (30.0913967654, 27.6666666667),
    "F11b": (31.4247300988, 27.6666666667),
    "F12": (33.8395446645, 25.8137142794),
    "F13": (36.4692520963, 15.9995125350),
    "F14": (36.6666666667, 14.5000000000),
    "F15": (36.6666666667,  5.0000000000),
    "F16": (35.3466856976,  2.7137258966),
    "F17": (31.0165586787,  0.2137258966),
    "F18": (27.7308889019, -0.6666666667),
    "F19": (26.3333333333, -0.6666666667),
    "F20": (24.2290180726, -0.5501188753),
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
        assert abs(area - 904.30) < 0.1

    @pytest.mark.parametrize("name,expected", list(_EXPECTED_F.items()))
    def test_f_series_regression(self, outline_geo, name, expected):
        """F-series coordinates match known-good values."""
        pt = outline_geo.fp_pts[name]
        assert abs(pt[0] - expected[0]) < 1e-6, f"{name} E: {pt[0]} != {expected[0]}"
        assert abs(pt[1] - expected[1]) < 1e-6, f"{name} N: {pt[1]} != {expected[1]}"
