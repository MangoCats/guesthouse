"""Tests for floorplan/layout.py — interior layout."""
import pytest
from floorplan.layout import InteriorLayout
from shared.geometry import poly_bbox


class TestInteriorLayout:
    def test_returns_named_tuple(self, layout):
        assert isinstance(layout, InteriorLayout)

    def test_iw1_polygon(self, layout):
        assert len(layout.iw1_poly) == 4
        assert layout.iw1_poly[3][1] > layout.iw1_poly[0][1]

    def test_iw2_bounds(self, layout):
        assert layout.iw2_poly[1][0] > layout.iw2_poly[0][0]
        assert layout.iw2_poly[3][1] > layout.iw2_poly[0][1]

    def test_dryer_bounds(self, layout):
        d = poly_bbox(layout.dryer_poly)
        assert d.e > d.w
        assert d.n > d.s

    def test_washer_above_dryer(self, layout):
        d = poly_bbox(layout.dryer_poly)
        w = poly_bbox(layout.washer_poly)
        assert w.s > d.n  # washer north of dryer (gap)

    def test_counter_bounds(self, layout):
        c = poly_bbox(layout.ctr_poly)
        assert c.e > c.w
        assert c.n > c.s
        assert layout.ctr_nw_r == 0

    def test_wall_ordering_east(self, layout):
        """IW11 west of IW4."""
        iw11_bb = poly_bbox(layout.iw11_poly)
        assert iw11_bb.e < layout.iw4_poly[0][0]

    def test_bed_poly(self, layout):
        assert len(layout.bed_poly) == 4
        # Bed polygon centroid west of IW11
        bed_cx = sum(p[0] for p in layout.bed_poly) / 4
        iw11_bb = poly_bbox(layout.iw11_poly)
        assert bed_cx < iw11_bb.w

    def test_iw5_bounds(self, layout):
        assert layout.iw5_poly[1][0] > layout.iw5_poly[0][0]
        assert layout.iw5_poly[3][1] > layout.iw5_poly[0][1]

    def test_iw6_polygon(self, layout):
        assert len(layout.iw6_poly) == 4
        assert layout.iw6_poly[3][1] > layout.iw6_poly[0][1]

    def test_wall_thicknesses(self, layout):
        assert abs(layout.iwt - 6.0 / 12.0) < 1e-12    # 6" = 0.5'
        assert abs(layout.iwt3 - 3.0 / 12.0) < 1e-12   # 3" = 0.25'
        assert abs(layout.iwt4 - 4.0 / 12.0) < 1e-12   # 4" ≈ 0.333'
