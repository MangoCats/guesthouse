"""Tests for app.engine — ENG-1 through ENG-10.

Verifies geometry computation, constant propagation, derived constants,
SVG file reading, and SVG generation.
"""
import os
import pytest
from app.database import (
    init_db, get_constants_dict, update_constant, reset_constants,
)
from app.engine import (
    compute_geometry, patch_constants, get_svg_content, generate_svg,
)

_PROJECT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

from tests.test_zapp_conftest import fresh_db, geometry


# ── ENG-1  Geometry Computation ──────────────────────────────────────

class TestENG1GeometryComputation:
    REQUIRED_KEYS = {
        "points", "outline_segments", "inner_segments", "outline_poly",
        "inner_poly", "interior_walls", "outer_openings", "rough_openings",
        "appliances", "furniture", "bbox",
    }

    def test_all_keys_present(self, geometry):
        assert self.REQUIRED_KEYS <= set(geometry.keys())

    def test_all_values_nonempty(self, geometry):
        for key in self.REQUIRED_KEYS:
            val = geometry[key]
            assert val, f"{key} is empty"


# ── ENG-2  Point Count ───────────────────────────────────────────────

class TestENG2PointCount:
    def test_at_least_50_points(self, geometry):
        assert len(geometry["points"]) >= 50

    def test_f_series_points(self, geometry):
        pts = geometry["points"]
        f_points = [k for k in pts if k.startswith("F")]
        assert len(f_points) >= 15, "Should have at least 15 F-series points"

    def test_w_series_points(self, geometry):
        pts = geometry["points"]
        w_points = [k for k in pts if k.startswith("W")]
        assert len(w_points) >= 15, "Should have at least 15 W-series points"

    def test_c_series_points(self, geometry):
        pts = geometry["points"]
        c_points = [k for k in pts if k.startswith("C")]
        assert len(c_points) >= 5, "Should have at least 5 C-series arc centers"


# ── ENG-3  Outline Segment Count ─────────────────────────────────────

class TestENG3SegmentCount:
    def test_outline_segment_count(self, geometry):
        assert len(geometry["outline_segments"]) == 18

    def test_inner_segment_count(self, geometry):
        assert len(geometry["inner_segments"]) == 18

    def test_segment_types_valid(self, geometry):
        for seg in geometry["outline_segments"]:
            assert seg["type"] in {"line", "arc"}


# ── ENG-4  Interior Wall Count ───────────────────────────────────────

class TestENG4InteriorWalls:
    EXPECTED = {"IW1", "IW2", "IW2O", "IW2S", "IW3", "IW4", "IW5",
                "IW6", "IW7", "IW8", "IW9", "IW11", "IW12"}

    def test_wall_names(self, geometry):
        assert set(geometry["interior_walls"].keys()) == self.EXPECTED

    def test_wall_has_poly_and_bbox(self, geometry):
        for name, wall in geometry["interior_walls"].items():
            assert "poly" in wall, f"{name} missing poly"
            assert "bbox" in wall, f"{name} missing bbox"
            assert len(wall["poly"]) == 4, f"{name} poly should have 4 corners"


# ── ENG-5  Opening Counts ────────────────────────────────────────────

class TestENG5OpeningCounts:
    def test_outer_opening_count(self, geometry):
        assert len(geometry["outer_openings"]) == 12

    def test_rough_opening_count(self, geometry):
        assert len(geometry["rough_openings"]) == 7

    def test_openings_have_names(self, geometry):
        for op in geometry["outer_openings"]:
            assert "name" in op
        for ro in geometry["rough_openings"]:
            assert "name" in ro


# ── ENG-6  Appliance and Furniture Counts ────────────────────────────

class TestENG6AppliancesFurniture:
    def test_appliance_keys(self, geometry):
        assert set(geometry["appliances"].keys()) == {"dryer", "washer", "counter"}

    def test_furniture_keys(self, geometry):
        assert set(geometry["furniture"].keys()) == {"bed", "dresser", "shelves"}

    def test_counter_has_clip(self, geometry):
        ctr = geometry["appliances"]["counter"]
        assert "clip" in ctr


# ── ENG-7  Constant Propagation ──────────────────────────────────────

class TestENG7ConstantPropagation:
    def test_bed_width_change_propagates(self, fresh_db):
        # Default geometry
        constants = get_constants_dict(fresh_db)
        geom1 = compute_geometry(constants)
        bed1 = geom1["furniture"]["bed"]
        w1 = bed1["bbox"]["e"] - bed1["bbox"]["w"]

        # Change BED_WIDTH to 80/12 (80 inches)
        update_constant("BED_WIDTH", 80.0 / 12.0, fresh_db)
        constants2 = get_constants_dict(fresh_db)
        geom2 = compute_geometry(constants2)
        bed2 = geom2["furniture"]["bed"]
        w2 = bed2["bbox"]["e"] - bed2["bbox"]["w"]

        assert abs(w2 - w1) > 0.01, "Bed width should change"
        # 80 inches is 6.6667 feet
        assert abs(w2 - 80.0 / 12.0) < 0.1 or abs((bed2["bbox"]["n"] - bed2["bbox"]["s"]) - 80.0 / 12.0) < 0.1


# ── ENG-8  Derived Constants ─────────────────────────────────────────

class TestENG8DerivedConstants:
    def test_derived_constants_recomputed(self, fresh_db):
        import floorplan.constants as mod
        constants = get_constants_dict(fresh_db)
        constants["WALL_OUTER"] = 10.0 / 12.0
        patch_constants(constants)
        expected_air_gap = 10.0 / 12.0 - 2 * mod.SHELL_THICKNESS
        assert abs(mod.AIR_GAP - expected_air_gap) < 1e-9

    def test_wall_extra_derived(self, fresh_db):
        import floorplan.constants as mod
        constants = get_constants_dict(fresh_db)
        constants["WALL_OUTER"] = 10.0 / 12.0
        patch_constants(constants)
        expected = 10.0 / 12.0 - 8.0 / 12.0
        assert abs(mod.WALL_EXTRA - expected) < 1e-9


# ── ENG-9  SVG File Reading ──────────────────────────────────────────

class TestENG9SVGReading:
    def test_existing_svg(self):
        svg_path = "floorplan/floorplan.svg"
        full = os.path.join(_PROJECT, svg_path)
        if not os.path.exists(full):
            pytest.skip("floorplan.svg not generated yet")
        content = get_svg_content(svg_path)
        assert content is not None
        assert content.strip().startswith("<")

    def test_nonexistent_svg(self):
        content = get_svg_content("nonexistent/does_not_exist.svg")
        assert content is None


# ── ENG-10  SVG Generation ───────────────────────────────────────────

class TestENG10SVGGeneration:
    def test_generate_floorplan(self):
        ok = generate_svg("floorplan", "floorplan/gen_floorplan.py")
        assert ok is True
        full = os.path.join(_PROJECT, "floorplan", "floorplan.svg")
        assert os.path.exists(full)

    def test_generate_nonexistent_script(self):
        ok = generate_svg("fake", "nonexistent/script.py")
        assert ok is False


# ── ENG-14  Variant Exclusion Filtering ───────────────────────────

class TestENG14VariantExclusions:
    def test_bare_excludes_iw6(self, fresh_db):
        """Bare variant omits IW6 from interior walls."""
        constants = get_constants_dict(fresh_db)
        geom = compute_geometry(constants, "bare")
        assert "IW6" not in geom["interior_walls"]

    def test_bare_excludes_ro5(self, fresh_db):
        """Bare variant omits RO5 from rough openings."""
        constants = get_constants_dict(fresh_db)
        geom = compute_geometry(constants, "bare")
        ro_names = {ro["name"] for ro in geom["rough_openings"]}
        assert "RO5" not in ro_names

    def test_standard_includes_both(self, fresh_db):
        """Standard variant includes IW6 and RO5."""
        constants = get_constants_dict(fresh_db)
        geom = compute_geometry(constants, "standard")
        assert "IW6" in geom["interior_walls"]
        ro_names = {ro["name"] for ro in geom["rough_openings"]}
        assert "RO5" in ro_names


# ── ENG-15  Room Label Computation ─────────────────────────────────

class TestENG15RoomLabels:
    EXPECTED_ROOMS = {
        "BEDROOM", "UTIL_N", "UTIL_S", "KITCHEN", "LIVING",
        "BATH", "OFFICE", "E CLOSET", "W CLOSET", "STORAGE", "WH",
    }

    def test_standard_has_11_labels(self, fresh_db):
        """Standard variant produces 11 room labels."""
        constants = get_constants_dict(fresh_db)
        geom = compute_geometry(constants, "standard")
        assert len(geom["room_labels"]) == 11

    def test_all_rooms_named(self, fresh_db):
        """All 11 expected room names are present."""
        constants = get_constants_dict(fresh_db)
        geom = compute_geometry(constants, "standard")
        names = {lbl["name"] for lbl in geom["room_labels"]}
        assert names == self.EXPECTED_ROOMS

    def test_labels_have_required_keys(self, fresh_db):
        """Each label has name, pos, and centroid keys."""
        constants = get_constants_dict(fresh_db)
        geom = compute_geometry(constants, "standard")
        for lbl in geom["room_labels"]:
            assert "name" in lbl, f"label missing 'name'"
            assert "pos" in lbl, f"{lbl.get('name')} missing 'pos'"
            assert "centroid" in lbl, f"{lbl.get('name')} missing 'centroid'"
            assert len(lbl["pos"]) == 2, f"{lbl['name']} pos not [E,N]"
            assert len(lbl["centroid"]) == 2, f"{lbl['name']} centroid not [E,N]"

    def test_label_positions_within_bbox(self, fresh_db):
        """All label positions fall within the building outline bbox."""
        constants = get_constants_dict(fresh_db)
        geom = compute_geometry(constants, "standard")
        ob = geom["bbox"]
        margin = 2.0
        for lbl in geom["room_labels"]:
            e, n = lbl["pos"]
            assert ob["w"] - margin < e < ob["e"] + margin, \
                f"{lbl['name']} pos E={e} outside bbox"
            assert ob["s"] - margin < n < ob["n"] + margin, \
                f"{lbl['name']} pos N={n} outside bbox"

    def test_sf_has_area_and_poly(self, fresh_db):
        """SF variant labels have area (positive number) and poly (list)."""
        constants = get_constants_dict(fresh_db)
        geom = compute_geometry(constants, "sf")
        for lbl in geom["room_labels"]:
            assert "area" in lbl, f"{lbl['name']} missing 'area'"
            assert lbl["area"] > 0, f"{lbl['name']} area not positive"
            assert "poly" in lbl, f"{lbl['name']} missing 'poly'"
            assert len(lbl["poly"]) >= 3, f"{lbl['name']} poly too short"

    def test_sf_has_sf_lines(self, fresh_db):
        """SF variant includes sf_lines with 3 partition lines."""
        constants = get_constants_dict(fresh_db)
        geom = compute_geometry(constants, "sf")
        assert "sf_lines" in geom
        assert len(geom["sf_lines"]) == 3
        for line in geom["sf_lines"]:
            assert "start" in line and "end" in line
            assert len(line["start"]) == 2 and len(line["end"]) == 2
