"""Golden-gate identity tests: GeneratorData vs build_floorplan_data().

Verifies that GeneratorData (DB-driven) produces geometry identical
(within tolerance) to the hardcoded procedural build_floorplan_data().
"""
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


# ── Inner wall overrides (Phase 15½-B) ─────────────────────────────────

from app.gen_provider import walk_override_chain
from app.database import get_inner_wall_overrides


class TestWalkOverrideChain:
    def test_straight_line(self):
        """Walk a single line segment: bearing 0° (due north) for 1 ft."""
        chain = [{"seg_type": "L", "bearing": 0.0, "distance": 1.0}]
        poly = walk_override_chain(chain, (0.0, 0.0), 0.0)
        assert len(poly) == 2
        assert abs(poly[0][0]) < TOL
        assert abs(poly[0][1]) < TOL
        assert abs(poly[1][0]) < TOL
        assert abs(poly[1][1] - 1.0) < TOL

    def test_straight_east(self):
        """Walk bearing 90° (due east) for 2 ft."""
        chain = [{"seg_type": "L", "bearing": 90.0, "distance": 2.0}]
        poly = walk_override_chain(chain, (0.0, 0.0), 90.0)
        assert abs(poly[1][0] - 2.0) < TOL
        assert abs(poly[1][1]) < TOL

    def test_line_arc_line(self):
        """Walk a line-arc-line chain and verify endpoint continuity."""
        chain = [
            {"seg_type": "L", "bearing": 180.0, "distance": 1.0},
            {"seg_type": "CCW", "radius": 0.5, "sweep": 90.0, "n_pts": 20},
            {"seg_type": "L", "bearing": 90.0, "distance": 1.0},
        ]
        poly = walk_override_chain(chain, (0.0, 0.0), 180.0)
        # Should have 1 (start) + 1 (line end) + 20 (arc pts) + 1 (line end) = 23
        assert len(poly) == 23
        # Final point: went south 1, turned CCW 90° (now heading east), went east 1
        # After going south 1: (0, -1). CCW 90° turn with r=0.5 around center
        # that's LEFT of south = (0.5, -1): ends at (0.5, -1.5).
        # Then east 1: (1.5, -1.5)
        assert abs(poly[-1][0] - 1.5) < TOL
        assert abs(poly[-1][1] + 1.5) < TOL

    def test_cw_arc(self):
        """Walk a CW arc 90° and verify endpoint."""
        chain = [
            {"seg_type": "CW", "radius": 1.0, "sweep": 90.0, "n_pts": 20},
        ]
        # Heading north, CW turn 90° → end heading east
        # Center is RIGHT of north = (1, 0). Start at (0, 0).
        # After 90° CW: (1, 1)
        poly = walk_override_chain(chain, (0.0, 0.0), 0.0)
        assert len(poly) == 21  # 1 start + 20 arc
        assert abs(poly[-1][0] - 1.0) < TOL
        assert abs(poly[-1][1] - 1.0) < TOL


class TestOverrideIdentity:
    """Verify DB-seeded overrides produce geometry correctly (DB-only, no hardcoded fallback)."""

    @pytest.fixture
    def gen_data_no_overrides(self, fresh_db):
        """Build GeneratorData without any overrides."""
        constants = get_constants_dict(fresh_db)
        return build_generator_data(constants, overrides={})

    @pytest.fixture
    def gen_data_with_overrides(self, fresh_db):
        """Build GeneratorData using DB-seeded overrides."""
        constants = get_constants_dict(fresh_db)
        overrides = get_inner_wall_overrides(fresh_db)
        return build_generator_data(constants, overrides=overrides)

    def test_overrides_loaded(self, fresh_db):
        """DB seed creates seg_index=5 override."""
        overrides = get_inner_wall_overrides(fresh_db)
        assert 5 in overrides
        assert len(overrides[5]) == 3  # L, CCW, L

    def test_inner_poly_override_splices_arc(
            self, gen_data_no_overrides, gen_data_with_overrides):
        """DB override adds arc points to inner_poly; no-override baseline is shorter."""
        n_with = len(gen_data_with_overrides.inner_poly)
        n_without = len(gen_data_no_overrides.inner_poly)
        assert n_with > n_without, (
            f"override should add splice points: {n_with} <= {n_without}")

    def test_w_f8f9_poly_from_override(
            self, gen_data_no_overrides, gen_data_with_overrides):
        """w_f8f9_poly populated only when DB override exists for F8→W9."""
        assert len(gen_data_with_overrides.w_f8f9_poly) > 0
        assert len(gen_data_no_overrides.w_f8f9_poly) == 0

    def test_inner_area_with_override(self, gen_data_no_overrides, gen_data_with_overrides):
        """Override reduces inner area slightly (arc splice cuts the corner)."""
        diff = gen_data_no_overrides.inner_area - gen_data_with_overrides.inner_area
        assert 0.0 < diff < 0.5, f"area reduction {diff:.4f} sq ft out of expected range"


# ── compute_default_override ───────────────────────────────────────────

from app.gen_provider import compute_native_geometry, compute_default_override


class TestComputeDefaultOverride:
    @pytest.fixture
    def geom(self):
        from app.database import get_constants_dict as _gcd
        constants = _gcd()
        pts, outline_segs, inner_segs, radii = compute_native_geometry(constants)
        return pts, inner_segs, constants

    def test_line_segment_default(self, geom):
        """Default for a line segment is a single L chain."""
        pts, inner_segs, constants = geom
        # seg 1 is W2->W5 (LineSeg)
        chain = compute_default_override(1, inner_segs, pts, constants)
        assert len(chain) == 1
        assert chain[0]["seg_type"] == "L"
        assert chain[0]["distance"] > 0
        assert chain[0]["bearing"] is not None

    def test_arc_segment_default(self, geom):
        """Default for an arc segment is a single arc chain."""
        pts, inner_segs, constants = geom
        # seg 0 is W1->W2 (ArcSeg CW)
        chain = compute_default_override(0, inner_segs, pts, constants)
        assert len(chain) == 1
        assert chain[0]["seg_type"] in ("CW", "CCW")
        assert chain[0]["radius"] > 0
        assert chain[0]["sweep"] > 0

    def test_w8w9_default(self, geom):
        """W8-W9 default is L-CCW-L chain."""
        pts, inner_segs, constants = geom
        chain = compute_default_override(5, inner_segs, pts, constants)
        assert len(chain) == 3
        assert chain[0]["seg_type"] == "L"
        assert chain[1]["seg_type"] == "CCW"
        assert chain[2]["seg_type"] == "L"

    def test_invalid_index(self, geom):
        pts, inner_segs, constants = geom
        assert compute_default_override(99, inner_segs, pts, constants) == []


# ── Multi-segment spans ───────────────────────────────────────────────

from app.gen_provider import (
    compute_default_span_override, validate_override_endpoint,
    _seg_start_bearing, apply_overrides_to_poly, _walk_chain_exit,
)
from shared.geometry import path_polygon


class TestMultiSegmentSpan:
    @pytest.fixture
    def geom(self):
        from app.database import get_constants_dict as _gcd
        constants = _gcd()
        pts, outline_segs, inner_segs, radii = compute_native_geometry(constants)
        return pts, inner_segs, constants

    def test_span_default_concatenates(self, geom):
        """Default span override concatenates individual defaults."""
        pts, inner_segs, constants = geom
        # Span segs 1-3 (three consecutive segments)
        span_chain = compute_default_span_override(1, 3, inner_segs, pts, constants)
        individual = []
        for i in range(1, 4):
            individual.extend(compute_default_override(i, inner_segs, pts, constants))
        assert len(span_chain) == len(individual)
        for s, ind in zip(span_chain, individual):
            assert s["seg_type"] == ind["seg_type"]

    def test_span_default_invalid(self, geom):
        """Invalid span range returns empty list."""
        pts, inner_segs, constants = geom
        assert compute_default_span_override(5, 3, inner_segs, pts, constants) == []
        assert compute_default_span_override(-1, 2, inner_segs, pts, constants) == []

    def test_apply_span_override(self, geom):
        """Apply a multi-segment span override to inner_poly."""
        pts, inner_segs, constants = geom
        # Get default span chain for segments 1-3
        span_chain = compute_default_span_override(1, 3, inner_segs, pts, constants)
        # Tag with span_end
        for sub in span_chain:
            sub["span_end"] = 3
        inner_poly = path_polygon(inner_segs, pts)
        orig_len = len(inner_poly)
        apply_overrides_to_poly(inner_poly, inner_segs, pts, {1: span_chain})
        # Polygon should still be valid (may differ in length due to arc discretization)
        assert len(inner_poly) > 3

    def test_walk_chain_exit_line(self):
        """Exit bearing for a line segment equals the line's bearing."""
        exit_brg = _walk_chain_exit(
            [{"seg_type": "L", "bearing": 45.0, "distance": 1.0}], 0.0)
        assert abs(exit_brg - 45.0) < 0.01

    def test_walk_chain_exit_arc(self):
        """Exit bearing for a CW arc rotates by sweep."""
        exit_brg = _walk_chain_exit(
            [{"seg_type": "CW", "radius": 1.0, "sweep": 90.0}], 0.0)
        assert abs(exit_brg - 90.0) < 0.01

    def test_validate_endpoint_ok(self, geom):
        """Default chain passes endpoint validation."""
        pts, inner_segs, constants = geom
        seg = inner_segs[1]
        chain = compute_default_override(1, inner_segs, pts, constants)
        start_pt = pts[seg.start]
        end_pt = pts[seg.end]
        start_brg = _seg_start_bearing(seg, pts)
        next_brg = _seg_start_bearing(inner_segs[2], pts)
        result = validate_override_endpoint(
            chain, start_pt, start_brg, end_pt, next_brg)
        assert result["ok"]
        assert result["warnings"] == []

    def test_validate_endpoint_bad(self, geom):
        """Off-target chain produces warnings."""
        pts, inner_segs, constants = geom
        seg = inner_segs[1]
        chain = [{"seg_type": "L", "bearing": 0.0, "distance": 100.0}]
        start_pt = pts[seg.start]
        end_pt = pts[seg.end]
        start_brg = _seg_start_bearing(seg, pts)
        next_brg = _seg_start_bearing(inner_segs[2], pts)
        result = validate_override_endpoint(
            chain, start_pt, start_brg, end_pt, next_brg)
        assert not result["ok"]
        assert len(result["warnings"]) > 0


class TestOverrideOverlap:
    def test_no_overlap(self, fresh_db):
        """Non-overlapping override doesn't conflict."""
        from app.database import check_override_overlap
        conflicts = check_override_overlap(10, None, fresh_db)
        assert conflicts == []

    def test_overlap_single(self, fresh_db):
        """Span overlapping existing seg 5 override conflicts."""
        from app.database import check_override_overlap
        conflicts = check_override_overlap(3, 6, fresh_db)
        assert 5 in conflicts

    def test_self_upsert_no_conflict(self, fresh_db):
        """Upserting at same seg_index doesn't conflict with itself."""
        from app.database import check_override_overlap
        conflicts = check_override_overlap(5, None, fresh_db)
        assert conflicts == []

    def test_span_end_in_db(self, fresh_db):
        """Upserted span_end is stored and retrievable."""
        from app.database import upsert_inner_wall_override, get_inner_wall_override
        chain = [{"seg_type": "L", "bearing": 0.0, "distance": 1.0}]
        upsert_inner_wall_override(10, chain, span_end=12, db_path=fresh_db)
        stored = get_inner_wall_override(10, fresh_db)
        assert len(stored) == 1
        assert stored[0]["span_end"] == 12
