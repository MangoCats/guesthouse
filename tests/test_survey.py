"""Tests for shared/survey.py — traverse, three-arc, inset."""
import math
import pytest
from shared.survey import compute_traverse, compute_three_arc, InsetResult, compute_inset


class TestComputeTraverse:
    def test_p3_is_origin(self, pts_base):
        assert pts_base["P3"] == (0, 0)

    def test_known_stations(self, pts_base):
        assert abs(pts_base["POB"][0] - 19.1177) < 1e-4
        assert abs(pts_base["POB"][1] - 32.9174) < 1e-3
        assert abs(pts_base["P2"][1] - 29.0) < 1e-4
        assert abs(pts_base["P4"][1] - 0.0) < 1e-4

    def test_p3_trav_tuple(self, traverse):
        _, p3_trav = traverse
        assert len(p3_trav) == 2
        assert isinstance(p3_trav[0], float)

    def test_five_stations(self, pts_base):
        for key in ["P3", "POB", "P2", "P4", "P5"]:
            assert key in pts_base


class TestComputeThreeArc:
    def test_positive_radii(self, arc_info):
        assert arc_info["R1"] > 0
        assert arc_info["R2"] > 0
        assert arc_info["R3"] > 0

    def test_known_radii(self, arc_info):
        assert abs(arc_info["R1"] - 10.0) < 1e-6
        assert abs(arc_info["R2"] - 12.5) < 1e-6
        assert abs(arc_info["R3"] - 11.0) < 1e-6

    def test_arc_centers_added(self, pts_base, arc_info):
        # arc_info fixture triggers compute_three_arc which mutates pts_base
        for key in ["T1", "TC1", "T2", "TC2", "PA", "T3", "TC3", "PX"]:
            assert key in pts_base, f"Missing {key} in pts"


class TestComputeInset:
    def test_returns_inset_result(self, inset_result):
        assert isinstance(inset_result, InsetResult)

    def test_inset_radii_offset(self, arc_info, inset_result):
        # CW arcs on the outer path: inset (interior side) adds delta
        assert abs(inset_result.R1i - (arc_info["R1"] + 0.5)) < 1e-10
        assert abs(inset_result.R2i - (arc_info["R2"] + 0.5)) < 1e-10
        assert abs(inset_result.R3i - (arc_info["R3"] + 0.5)) < 1e-10

    def test_pts_update_keys(self, inset_result):
        expected = {"PiOB", "Pi2", "Pi3", "Pi4", "Pi5", "Ti1", "Ti2", "Ti3", "PiX", "Ai2"}
        assert expected <= set(inset_result.pts_update.keys())

    def test_does_not_mutate_pts(self, pts_base, arc_info):
        """compute_inset should not mutate the input pts dict."""
        original_keys = set(pts_base.keys())
        # Note: pts_base was already mutated by arc_info fixture, but
        # compute_inset should not add anything new
        _ = compute_inset(
            pts_base, arc_info["R1"], arc_info["R2"], arc_info["R3"],
            arc_info["nE"], arc_info["nN"],
        )
        assert set(pts_base.keys()) == original_keys
