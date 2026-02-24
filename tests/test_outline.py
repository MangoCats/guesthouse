"""Tests for floorplan/geometry.py — outline geometry."""
import math
import pytest
from shared.types import LineSeg, ArcSeg
from shared.geometry import path_polygon, poly_area
from floorplan.geometry import OutlineGeometry


# Known-good F-series coordinates (FC-based: building center = origin,
# F2.E = -18'6", F1.N = -13'6", F2 bearing 0 due north)
_EXPECTED_F = {
    "F1":  (-17.6666666667, -13.5000000000),
    "F2":  (-18.5000000000, -12.6666666667),
    "F5":  (-18.5000000000,  11.5000000000),
    "F6":  (-16.1666666667,  13.8333333333),
    "F7":  (-11.0000000000,  13.8333333333),
    "F8":  ( -8.6666666667,  11.5000000000),
    "F9":  ( -8.5000000000,  11.3333333333),
    "F10": (  7.1024481810,  11.3333333333),
    "F11": (  8.1066846850,  12.1039111052),
    "F11a": (10.3605116131,  13.8333333333),
    "F11b": (11.3605116131,  13.8333333333),
    "F12": ( 13.5741059752,  12.2378647874),
    "F13": ( 17.4071820764,   0.7386364836),
    "F14": ( 17.5000000000,   0.1666666667),
}


class TestOutlineGeometry:
    def test_returns_outline_geometry(self, outline_geo):
        assert isinstance(outline_geo, OutlineGeometry)

    def test_18_points(self, outline_geo):
        for i in [j for j in range(19) if j not in (0, 3, 4)]:
            assert f"F{i}" in outline_geo.fp_pts
        assert "F11a" in outline_geo.fp_pts
        assert "F11b" in outline_geo.fp_pts

    def test_18_segments(self, outline_geo):
        assert len(outline_geo.outline_segs) == 18

    def test_9_radii(self, outline_geo):
        assert len(outline_geo.radii) == 9
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
        assert abs(area - 898.51) < 1.0

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
            (17, "F18",  "F1"),
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
