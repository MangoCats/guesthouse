"""Golden-gate identity tests: GeneratorData vs build_floorplan_data().

Verifies that GeneratorData (DB-driven) produces geometry identical
(within tolerance) to the hardcoded procedural build_floorplan_data().
"""
import math
import os
import sys
import pytest

_PROJECT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _PROJECT not in sys.path:
    sys.path.insert(0, _PROJECT)

from tests.test_zapp_conftest import fresh_db
from app.database import get_constants_dict
from app.gen_provider import build_generator_data
from floorplan.gen_floorplan import build_floorplan_data

TOL = 1e-8  # 1e-8 ft ≈ 0.0001 mil


def _pt_dist2(a, b):
    """Squared distance between two points."""
    return (a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2


@pytest.fixture
def reference():
    """Build reference FloorplanData from hardcoded modules."""
    return build_floorplan_data()


@pytest.fixture
def gen_data(fresh_db):
    """Build GeneratorData from a freshly-seeded database."""
    constants = get_constants_dict(fresh_db)
    return build_generator_data(constants)


# ── Points ──────────────────────────────────────────────────────────────

class TestPointIdentity:
    def test_f_series_points(self, reference, gen_data):
        """All F-series points match within tolerance."""
        for name in ("F1", "F2", "F5", "F6", "F7", "F8", "F9", "F10",
                      "F11", "F11a", "F11b", "F12", "F13", "F14", "F15",
                      "F16", "F17", "F18"):
            ref = reference.pts[name]
            gen = gen_data.pts[name]
            d2 = _pt_dist2(ref, gen)
            assert d2 < TOL, f"{name}: d²={d2}"

    def test_w_series_points(self, reference, gen_data):
        """All W-series (inner wall) points match."""
        for name, pt in reference.pts.items():
            if name.startswith("W"):
                gen = gen_data.pts.get(name)
                assert gen is not None, f"Missing {name} in GeneratorData"
                d2 = _pt_dist2(pt, gen)
                assert d2 < TOL, f"{name}: d²={d2}"

    def test_c_series_points(self, reference, gen_data):
        """All C-series (arc center) points match."""
        for name, pt in reference.pts.items():
            if name.startswith("C"):
                gen = gen_data.pts.get(name)
                assert gen is not None, f"Missing {name}"
                d2 = _pt_dist2(pt, gen)
                assert d2 < TOL, f"{name}: d²={d2}"

    def test_p_series_points(self, reference, gen_data):
        """Survey traverse points match."""
        for name in ("POB", "P2", "P3", "P4", "P5"):
            ref = reference.pts[name]
            gen = gen_data.pts[name]
            d2 = _pt_dist2(ref, gen)
            assert d2 < TOL, f"{name}: d²={d2}"


# ── Segments ────────────────────────────────────────────────────────────

class TestSegmentIdentity:
    def test_outline_seg_count(self, reference, gen_data):
        assert len(gen_data.outline_segs) == len(reference.outline_segs)

    def test_outline_seg_endpoints(self, reference, gen_data):
        for i, (ref, gen) in enumerate(zip(reference.outline_segs,
                                            gen_data.outline_segs)):
            assert ref.start == gen.start, f"seg {i} start: {ref.start} vs {gen.start}"
            assert ref.end == gen.end, f"seg {i} end: {ref.end} vs {gen.end}"

    def test_inner_seg_count(self, reference, gen_data):
        assert len(gen_data.inner_segs) == len(reference.inner_segs)

    def test_inner_seg_endpoints(self, reference, gen_data):
        for i, (ref, gen) in enumerate(zip(reference.inner_segs,
                                            gen_data.inner_segs)):
            assert ref.start == gen.start, f"seg {i}: {ref.start} vs {gen.start}"
            assert ref.end == gen.end, f"seg {i}: {ref.end} vs {gen.end}"


# ── Polygons ────────────────────────────────────────────────────────────

class TestPolygonIdentity:
    def test_outline_poly_count(self, reference, gen_data):
        assert len(gen_data.outline_poly) == len(reference.outer_poly)

    def test_outline_poly_points(self, reference, gen_data):
        for i, (ref, gen) in enumerate(zip(reference.outer_poly,
                                            gen_data.outline_poly)):
            d2 = _pt_dist2(ref, gen)
            assert d2 < TOL, f"outline_poly[{i}]: d²={d2}"

    def test_inner_poly_count(self, reference, gen_data):
        assert len(gen_data.inner_poly) == len(reference.inner_poly)

    def test_inner_poly_points(self, reference, gen_data):
        for i, (ref, gen) in enumerate(zip(reference.inner_poly,
                                            gen_data.inner_poly)):
            d2 = _pt_dist2(ref, gen)
            assert d2 < TOL, f"inner_poly[{i}]: d²={d2}"

    def test_areas(self, reference, gen_data):
        assert abs(gen_data.outer_area - reference.outer_area) < TOL
        assert abs(gen_data.inner_area - reference.inner_area) < TOL


# ── Radii ───────────────────────────────────────────────────────────────

class TestRadiiIdentity:
    def test_radii_keys(self, reference, gen_data):
        assert set(reference.radii.keys()) == set(gen_data.radii.keys())

    def test_radii_values(self, reference, gen_data):
        for name, val in reference.radii.items():
            assert abs(gen_data.radii[name] - val) < TOL, \
                f"{name}: {gen_data.radii[name]} vs {val}"


# ── Shell geometry (Phase 15-B) ─────────────────────────────────────────

class TestShellIdentity:
    def test_s_series_points(self, reference, gen_data):
        """S-series shell boundary points match."""
        for name, pt in reference.pts.items():
            if name.startswith("S"):
                gen = gen_data.pts.get(name)
                assert gen is not None, f"Missing {name}"
                d2 = _pt_dist2(pt, gen)
                assert d2 < TOL, f"{name}: d²={d2}"

    def test_g_series_points(self, reference, gen_data):
        """G-series shell boundary points match."""
        for name, pt in reference.pts.items():
            if name.startswith("G"):
                gen = gen_data.pts.get(name)
                assert gen is not None, f"Missing {name}"
                d2 = _pt_dist2(pt, gen)
                assert d2 < TOL, f"{name}: d²={d2}"

    def test_s_seg_count(self, reference, gen_data):
        assert len(gen_data.s_segs) == len(reference.s_segs)

    def test_g_seg_count(self, reference, gen_data):
        assert len(gen_data.g_segs) == len(reference.g_segs)

    def test_s_seg_endpoints(self, reference, gen_data):
        for i, (ref, gen) in enumerate(zip(reference.s_segs,
                                            gen_data.s_segs)):
            assert ref.start == gen.start, f"s_seg {i}: {ref.start} vs {gen.start}"
            assert ref.end == gen.end, f"s_seg {i}: {ref.end} vs {gen.end}"

    def test_g_seg_endpoints(self, reference, gen_data):
        for i, (ref, gen) in enumerate(zip(reference.g_segs,
                                            gen_data.g_segs)):
            assert ref.start == gen.start, f"g_seg {i}: {ref.start} vs {gen.start}"
            assert ref.end == gen.end, f"g_seg {i}: {ref.end} vs {gen.end}"

    def test_w_f8f9_poly(self, reference, gen_data):
        assert len(gen_data.w_f8f9_poly) == len(reference.w_f8f9_poly)
        for i, (ref, gen) in enumerate(zip(reference.w_f8f9_poly,
                                            gen_data.w_f8f9_poly)):
            d2 = _pt_dist2(ref, gen)
            assert d2 < TOL, f"w_f8f9[{i}]: d²={d2}"

    def test_g_f8f9_poly(self, reference, gen_data):
        assert len(gen_data.g_f8f9_poly) == len(reference.g_f8f9_poly)
        for i, (ref, gen) in enumerate(zip(reference.g_f8f9_poly,
                                            gen_data.g_f8f9_poly)):
            d2 = _pt_dist2(ref, gen)
            assert d2 < TOL, f"g_f8f9[{i}]: d²={d2}"

    def test_openings_count(self, reference, gen_data):
        assert len(gen_data.openings) == len(reference.openings)

    def test_openings_params(self, reference, gen_data):
        for i, (ref, gen) in enumerate(zip(reference.openings,
                                            gen_data.openings)):
            assert ref.name == gen.name, f"opening {i}: {ref.name} vs {gen.name}"
            assert ref.seg_idx == gen.seg_idx
            assert abs(ref.t_start - gen.t_start) < TOL
            assert abs(ref.t_end - gen.t_end) < TOL

    def test_wall_sections_count(self, reference, gen_data):
        ref_sections = len(gen_data.wall_sections)
        assert ref_sections == len(gen_data.openings)  # one section per opening


# ── Roof geometry (Phase 15-C) ──────────────────────────────────────────

class TestRoofIdentity:
    def test_roof_points(self, gen_data):
        """R-series points exist and are reasonable."""
        for name in ("R1", "R3s", "R3e", "R4s", "R4e", "R5", "R6", "R7"):
            assert name in gen_data.pts, f"Missing {name}"

    def test_roof_area_positive(self, gen_data):
        assert gen_data.roof.area > 0

    def test_roof_poly_nonempty(self, gen_data):
        assert len(gen_data.roof_poly) > 10

    def test_roof_vs_hardcoded(self, reference, gen_data):
        """Roof geometry from GeneratorData matches hardcoded computation."""
        from floorplan.roof import compute_roof_geometry, roof_polyline
        ref_roof = compute_roof_geometry(reference.pts, reference.radii)
        for name, pt in ref_roof.pts.items():
            gen = gen_data.roof.pts[name]
            d2 = _pt_dist2(pt, gen)
            assert d2 < TOL, f"roof {name}: d²={d2}"
        assert abs(gen_data.roof.area - ref_roof.area) < TOL


# ── Layout identity ─────────────────────────────────────────────────────

class TestLayoutIdentity:
    def test_layout_exists(self, gen_data):
        assert gen_data.layout is not None

    def test_wall_t(self, reference, gen_data):
        assert abs(gen_data.wall_t - reference.wall_t) < TOL
