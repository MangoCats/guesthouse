"""Tests for floorplan/layout.py — interior layout."""
import pytest
from floorplan.layout import InteriorLayout


class TestInteriorLayout:
    def test_returns_named_tuple(self, layout):
        assert isinstance(layout, InteriorLayout)

    def test_iw1_polygon(self, layout):
        assert len(layout.iw1) == 4
        assert layout.iw1_n > layout.iw1_s

    def test_iw2_bounds(self, layout):
        assert layout.iw2.e > layout.iw2.w
        assert layout.iw2.n > layout.iw2.s

    def test_dryer_bounds(self, layout):
        assert layout.dryer.e > layout.dryer.w
        assert layout.dryer.n > layout.dryer.s

    def test_washer_above_dryer(self, layout):
        assert layout.washer.s > layout.dryer.n  # washer north of dryer (gap)

    def test_counter_bounds(self, layout):
        assert layout.ctr.e > layout.ctr.w
        assert layout.ctr.n > layout.ctr.s
        assert layout.ctr_nw_r == 0

    def test_wall_ordering_east(self, layout):
        """IW11 west of IW4."""
        assert layout.iw11.e < layout.iw4_w

    def test_bed_poly(self, layout):
        assert len(layout.bed_poly) == 4
        # Bed polygon centroid west of IW11
        bed_cx = sum(p[0] for p in layout.bed_poly) / 4
        assert bed_cx < layout.iw11.w

    def test_iw5_bounds(self, layout):
        assert layout.iw5.e > layout.iw5.w
        assert layout.iw5.n > layout.iw5.s

    def test_iw6_polygon(self, layout):
        assert len(layout.iw6_poly) == 4
        assert layout.iw6_n > layout.iw6_s

    def test_wall_thicknesses(self, layout):
        assert abs(layout.iwt - 6.0 / 12.0) < 1e-12    # 6" = 0.5'
        assert abs(layout.iwt3 - 3.0 / 12.0) < 1e-12   # 3" = 0.25'
        assert abs(layout.iwt4 - 4.0 / 12.0) < 1e-12   # 4" ≈ 0.333'
