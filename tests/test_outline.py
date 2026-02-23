"""Tests for floorplan/geometry.py — outline geometry."""
import math
import pytest
from shared.types import LineSeg, ArcSeg
from shared.geometry import path_polygon, poly_area
from floorplan.geometry import OutlineGeometry


# Known-good F-series coordinates (FC-based: building center = origin,
# F2.E = -18', F1.N = -13', F2 bearing 0 due north)
_EXPECTED_F = {
    "F1":  (-17.1666666667, -13.0000000000),
    "F2":  (-18.0000000000, -12.1666666667),
    "F5":  (-18.0000000000,  10.1497042276),
    "F6":  (-15.9243402275,  12.4687662752),
    "F7":  (-10.7064506205,  13.0485317871),
    "F8":  ( -8.1297150121,  10.9871433004),
    "F9":  ( -7.9456624686,  10.8399012656),
    "F10": (  7.1282408406,  12.5147794111),
    "F11": (  8.0412390885,  13.3915434945),
    "F11a": (10.0902982771,  15.3592816646),
    "F11b": (11.0841820118,  15.4697131907),
    "F12": ( 13.5152066722,  13.9997621150),
    "F13": ( 17.8307503524,   2.9538659795),
    "F14": ( 17.9873409046,   2.3182685035),
    "F15": ( 18.9444141306,  -6.2953905304),
    "F16": ( 17.9518674804,  -8.5607912999),
    "F17": ( 13.9243034824, -11.5236831714),
    "F18": ( 10.8363106361, -12.7301373919),
    "F19": (  9.4473028882, -12.8844715862),
    "F20": (  7.3613928818, -13.0000000000),
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

    def test_segment_index_mapping(self, outline_geo):
        """Segment start→end names match expected indices.

        Renderers use hardcoded indices (e.g. F8→F9 at [5]) for special
        corner handling. This test catches index drift if segments are
        added, removed, or reordered.
        """
        expected = [
            (0,  "F1",   "F2"),
            (1,  "F2",   "F5"),
            (2,  "F5",   "F6"),
            (3,  "F6",   "F7"),
            (4,  "F7",   "F8"),
            (5,  "F8",   "F9"),
            (6,  "F9",   "F10"),
            (7,  "F10",  "F11"),
            (8,  "F11",  "F11a"),
            (9,  "F11a", "F11b"),
            (10, "F11b", "F12"),
            (11, "F12",  "F13"),
            (12, "F13",  "F14"),
            (13, "F14",  "F15"),
            (14, "F15",  "F16"),
            (15, "F16",  "F17"),
            (16, "F17",  "F18"),
            (17, "F18",  "F19"),
            (18, "F19",  "F20"),
            (19, "F20",  "F1"),
        ]
        segs = outline_geo.outline_segs
        for idx, start, end in expected:
            assert segs[idx].start == start and segs[idx].end == end, (
                f"outline_segs[{idx}]: expected {start}→{end}, "
                f"got {segs[idx].start}→{segs[idx].end}"
            )

    @pytest.mark.parametrize("name,expected", list(_EXPECTED_F.items()))
    def test_f_series_regression(self, outline_geo, name, expected):
        """F-series coordinates match known-good values."""
        pt = outline_geo.fp_pts[name]
        assert abs(pt[0] - expected[0]) < 1e-6, f"{name} E: {pt[0]} != {expected[0]}"
        assert abs(pt[1] - expected[1]) < 1e-6, f"{name} N: {pt[1]} != {expected[1]}"
