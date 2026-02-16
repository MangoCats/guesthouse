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
        assert layout.ctr_nw_r > 0

    def test_wall_ordering_east(self, layout):
        """IW9 west of IW4 west of IW8."""
        assert layout.iw9.w < layout.iw4_w
        assert layout.iw4_e < layout.iw8_w

    def test_iw3_above_iw7(self, layout):
        """IW3 east face aligns with IW2 east, starts at IW7 L north face."""
        assert abs(layout.iw3.e - layout.iw2.e) < 1e-12  # east face = IW2 east face
        assert abs(layout.iw3.s - layout.iw7[4][1]) < 1e-12  # south = IW7 L north

    def test_iw9_bounds(self, layout):
        """IW9 at old IW3 position, south of IW7 L north face."""
        assert layout.iw9.e > layout.iw9.w
        assert layout.iw9.n > layout.iw9.s
        assert abs(layout.iw9.n - layout.iw3.s) < 1e-12  # IW9 top = IW3 bottom

    def test_iw7_polygon(self, layout):
        assert len(layout.iw7) == 6  # L-shape

    def test_iw8_polygon(self, layout):
        assert len(layout.iw8) == 6  # L-shape

    def test_bed_bounds(self, layout):
        assert layout.bed.e > layout.bed.w
        assert layout.bed.n > layout.bed.s
        # Bed center between IW9 and IW4
        assert layout.iw9.e < layout.bed_cx < layout.iw4_w

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
