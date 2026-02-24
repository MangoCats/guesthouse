"""Tests for floorplan/layout.py — interior layout."""
import math
import pytest
from floorplan.layout import InteriorLayout
from floorplan.constants import WALL_6IN, WALL_3IN, WALL_4IN, IW6_THICKNESS, IW2O_THICKNESS


def _wall_thickness(poly):
    """Perpendicular distance between the narrower pair of opposite sides.

    A wall quadrilateral has two parallel long sides (separated by the wall
    thickness) and two ends that may not be perpendicular to the long sides.
    This function returns the thickness — the smaller of the two opposite-side
    perpendicular distances.
    """
    def _perp_dist(p, a, b):
        dx, dy = b[0] - a[0], b[1] - a[1]
        length = math.sqrt(dx * dx + dy * dy)
        if length < 1e-12:
            return 0.0
        return abs(dx * (p[1] - a[1]) - dy * (p[0] - a[0])) / length
    # Perpendicular distance between sides 0→1 and 3→2
    d_01_32 = _perp_dist(poly[3], poly[0], poly[1])
    # Perpendicular distance between sides 1→2 and 0→3
    d_12_03 = _perp_dist(poly[0], poly[1], poly[2])
    return min(d_01_32, d_12_03)


class TestInteriorLayout:
    def test_returns_named_tuple(self, layout):
        assert isinstance(layout, InteriorLayout)

    def test_iw1_polygon(self, layout):
        assert len(layout.iw1.poly) == 4
        assert layout.iw1.n > layout.iw1.s

    def test_iw2_bounds(self, layout):
        assert layout.iw2.e > layout.iw2.w
        assert layout.iw2.n > layout.iw2.s

    def test_iw2o_bounds(self, layout):
        assert layout.iw2o.e > layout.iw2o.w
        assert layout.iw2o.n > layout.iw2o.s

    def test_iw2s_bounds(self, layout):
        assert layout.iw2s.e > layout.iw2s.w
        assert layout.iw2s.n > layout.iw2s.s

    def test_dryer_bounds(self, layout):
        assert layout.dryer.e > layout.dryer.w
        assert layout.dryer.n > layout.dryer.s

    def test_washer_above_dryer(self, layout):
        # Compare polygon centroids (rotation-safe, not BBox edges)
        wc_n = sum(p[1] for p in layout.washer.poly) / len(layout.washer.poly)
        dc_n = sum(p[1] for p in layout.dryer.poly) / len(layout.dryer.poly)
        assert wc_n > dc_n  # washer centroid north of dryer centroid

    def test_counter_bounds(self, layout):
        assert layout.ctr.e > layout.ctr.w
        assert layout.ctr.n > layout.ctr.s
        assert layout.ctr_nw_r == 0

    def test_wall_ordering_east(self, layout):
        """IW11 west of IW4."""
        assert layout.iw11.e < layout.iw4.w

    def test_bed_poly(self, layout):
        assert len(layout.bed.poly) == 4
        # Bed polygon centroid west of IW11
        bed_cx = sum(p[0] for p in layout.bed.poly) / 4
        assert bed_cx < layout.iw11.w

    def test_iw5_bounds(self, layout):
        assert layout.iw5.e > layout.iw5.w
        assert layout.iw5.n > layout.iw5.s

    def test_iw6_polygon(self, layout):
        assert len(layout.iw6.poly) == 4
        assert layout.iw6.n > layout.iw6.s

    def test_wall_thicknesses(self, layout):
        assert abs(WALL_6IN - 6.0 / 12.0) < 1e-12    # 6" = 0.5'
        assert abs(WALL_3IN - 3.0 / 12.0) < 1e-12   # 3" = 0.25'
        assert abs(WALL_4IN - 4.0 / 12.0) < 1e-12   # 4" ≈ 0.333'


# Expected thickness constant for each interior wall
_IW_EXPECTED_THICKNESS = {
    "iw1": WALL_6IN,
    "iw2": WALL_6IN,
    "iw2o": IW2O_THICKNESS,
    "iw2s": WALL_6IN,
    "iw3": WALL_4IN,
    "iw4": WALL_4IN,
    "iw5": WALL_3IN,
    "iw6": IW6_THICKNESS,
    "iw7": WALL_4IN,
    "iw8": WALL_6IN,
    "iw9": WALL_4IN,
    "iw11": WALL_4IN,
    "iw12": WALL_4IN,
    "iw16": WALL_4IN,
}


class TestWallThicknesses:
    """Verify that every interior wall has the correct thickness.

    Thickness is the perpendicular distance between the two parallel long
    sides of the wall quadrilateral.
    """

    @pytest.mark.parametrize("iw_name,expected", list(_IW_EXPECTED_THICKNESS.items()),
                             ids=list(_IW_EXPECTED_THICKNESS.keys()))
    def test_iw_thickness(self, layout, iw_name, expected):
        wall = getattr(layout, iw_name)
        thickness = _wall_thickness(wall.poly)
        assert thickness == pytest.approx(expected, abs=1e-10), \
            f"{iw_name.upper()} thickness {thickness * 12:.4f}\" != {expected * 12:.4f}\""
