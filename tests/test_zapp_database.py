"""Tests for app.database — DB-1 through DB-8, DB-12, DB-13, shapes.

Verifies schema initialisation, seeding, and CRUD operations using an
isolated temporary database per test.
"""
import json
import os
import sqlite3
import pytest
from app.database import (
    init_db, get_db, get_all_constants, get_constants_dict,
    update_constant, update_constants_batch, get_categories,
    get_outline_chain, get_views, reset_constants,
    get_shapes, get_shape, get_variant_exclusions,
    get_room_label_offsets, set_room_label_offset,
    get_survey_legs, get_survey_config,
    update_survey_leg, update_survey_config, reset_survey,
    export_project, import_project,
)

# Re-use the fresh_db fixture from test_zapp_conftest.py
from tests.test_zapp_conftest import fresh_db


# ── DB-1  Schema Initialisation ─────────────────────────────────────

class TestDB1SchemaInit:
    def test_tables_created(self, fresh_db):
        """DB file contains all six tables."""
        with get_db(fresh_db) as conn:
            rows = conn.execute(
                "SELECT name FROM sqlite_master WHERE type='table' ORDER BY name"
            ).fetchall()
            names = {r["name"] for r in rows}
        assert {"constants", "outline_chain", "views",
                "shapes", "variant_exclusions", "room_label_offsets"} <= names

    def test_db_file_exists(self, fresh_db):
        assert os.path.exists(fresh_db)

    def test_wal_mode(self, fresh_db):
        with get_db(fresh_db) as conn:
            mode = conn.execute("PRAGMA journal_mode").fetchone()[0]
        assert mode == "wal"


# ── DB-2  Constants Seeding ──────────────────────────────────────────

class TestDB2ConstantsSeeding:
    def test_constant_count(self, fresh_db):
        """Seeded constants count >= 140."""
        constants = get_all_constants(fresh_db)
        assert len(constants) >= 140

    def test_constants_match_module(self, fresh_db):
        """Every constant in floorplan.constants is present in the database."""
        import floorplan.constants as mod
        module_names = {
            a for a in dir(mod)
            if a[0].isupper() and not a.startswith("_")
            and isinstance(getattr(mod, a), (int, float))
        }
        db_names = {c["name"] for c in get_all_constants(fresh_db)}
        missing = module_names - db_names
        assert missing == set(), f"Missing constants: {missing}"

    def test_constant_has_all_fields(self, fresh_db):
        """Each constant row has the required fields."""
        constants = get_all_constants(fresh_db)
        required = {"name", "value", "expr", "unit", "category", "description"}
        for c in constants:
            assert required <= set(c.keys()), f"Missing fields in {c['name']}"


# ── DB-3  Constant Categories ────────────────────────────────────────

class TestDB3ConstantCategories:
    EXPECTED_CATEGORIES = {
        "wall", "interior_wall", "opening", "appliance",
        "furniture", "fixture", "geometry", "construction", "misc",
    }

    def test_categories_match(self, fresh_db):
        cats = set(get_categories(fresh_db))
        assert cats == self.EXPECTED_CATEGORIES

    def test_each_constant_has_valid_category(self, fresh_db):
        for c in get_all_constants(fresh_db):
            assert c["category"] in self.EXPECTED_CATEGORIES, (
                f"{c['name']} has unexpected category {c['category']}"
            )


# ── DB-4  Outline Chain Seeding ──────────────────────────────────────

class TestDB4OutlineChain:
    def test_chain_count(self, fresh_db):
        chain = get_outline_chain(fresh_db)
        assert len(chain) == 18

    def test_seg_types_valid(self, fresh_db):
        chain = get_outline_chain(fresh_db)
        for seg in chain:
            assert seg["seg_type"] in {"L", "CW", "CCW"}

    def test_end_names_match_geometry(self, fresh_db):
        from floorplan.geometry import CHAIN_POINT_NAMES
        chain = get_outline_chain(fresh_db)
        db_names = [seg["end_name"] for seg in chain]
        assert db_names == list(CHAIN_POINT_NAMES)

    def test_seq_is_contiguous(self, fresh_db):
        chain = get_outline_chain(fresh_db)
        seqs = [seg["seq"] for seg in chain]
        assert seqs == list(range(18))


# ── DB-5  Views Seeding ──────────────────────────────────────────────

class TestDB5Views:
    def test_view_count(self, fresh_db):
        views = get_views(fresh_db)
        assert len(views) >= 11

    def test_view_has_required_fields(self, fresh_db):
        views = get_views(fresh_db)
        for v in views:
            assert "name" in v
            assert "label" in v
            assert "script" in v
            assert "svg_path" in v
            assert "category" in v

    def test_view_scripts_exist(self, fresh_db):
        project = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        views = get_views(fresh_db)
        for v in views:
            script_path = os.path.join(project, v["script"])
            assert os.path.exists(script_path), f"Missing script: {v['script']}"


# ── DB-6  Constant Update ────────────────────────────────────────────

class TestDB6ConstantUpdate:
    def test_update_returns_true(self, fresh_db):
        ok = update_constant("BED_WIDTH", 7.0, fresh_db)
        assert ok is True

    def test_updated_value_persists(self, fresh_db):
        update_constant("BED_WIDTH", 7.0, fresh_db)
        d = get_constants_dict(fresh_db)
        assert d["BED_WIDTH"] == 7.0

    def test_update_nonexistent_returns_false(self, fresh_db):
        ok = update_constant("DOES_NOT_EXIST", 1.0, fresh_db)
        assert ok is False


# ── DB-7  Batch Update ───────────────────────────────────────────────

class TestDB7BatchUpdate:
    def test_batch_updates_count(self, fresh_db):
        n = update_constants_batch(
            {"BED_WIDTH": 7.0, "BED_LENGTH": 8.0}, fresh_db
        )
        assert n == 2

    def test_batch_values_persist(self, fresh_db):
        update_constants_batch(
            {"BED_WIDTH": 7.0, "BED_LENGTH": 8.0}, fresh_db
        )
        d = get_constants_dict(fresh_db)
        assert d["BED_WIDTH"] == 7.0
        assert d["BED_LENGTH"] == 8.0


# ── DB-8  Reset Constants ────────────────────────────────────────────

class TestDB8ResetConstants:
    def test_reset_restores_original(self, fresh_db):
        import floorplan.constants as mod
        original = float(mod.BED_WIDTH)
        update_constant("BED_WIDTH", 999.0, fresh_db)
        reset_constants(fresh_db)
        d = get_constants_dict(fresh_db)
        assert abs(d["BED_WIDTH"] - original) < 1e-9

    def test_reset_preserves_count(self, fresh_db):
        before = len(get_all_constants(fresh_db))
        reset_constants(fresh_db)
        after = len(get_all_constants(fresh_db))
        assert after == before


# ── DB-12  Variant Exclusions Seeding ──────────────────────────────

class TestDB12VariantExclusions:
    def test_exclusion_count(self, fresh_db):
        """Variant exclusions table has 8 rows (walls, openings, dimensions)."""
        with get_db(fresh_db) as conn:
            count = conn.execute(
                "SELECT count(*) FROM variant_exclusions"
            ).fetchone()[0]
        assert count == 8

    def test_bare_exclusions(self, fresh_db):
        """Bare variant excludes IW6, RO5, dim12a, dim12b."""
        excl = get_variant_exclusions("bare", fresh_db)
        assert excl == {
            "wall": {"IW6"}, "rough_opening": {"RO5"},
            "dimension": {"dim12a", "dim12b"},
        }

    def test_standard_no_exclusions(self, fresh_db):
        """Standard variant has no exclusions."""
        excl = get_variant_exclusions("standard", fresh_db)
        assert excl == {}


# ── DB-13  Room Label Offsets ──────────────────────────────────────

class TestDB13RoomLabelOffsets:
    def test_initial_empty(self, fresh_db):
        """Room label offsets are empty by default."""
        offsets = get_room_label_offsets(fresh_db)
        assert offsets == {}

    def test_set_and_get(self, fresh_db):
        """Setting an offset makes it retrievable."""
        set_room_label_offset("BEDROOM", 0.5, -0.3, fresh_db)
        offsets = get_room_label_offsets(fresh_db)
        assert "BEDROOM" in offsets
        assert abs(offsets["BEDROOM"][0] - 0.5) < 1e-9
        assert abs(offsets["BEDROOM"][1] - (-0.3)) < 1e-9

    def test_upsert_overwrites(self, fresh_db):
        """Setting the same room twice overwrites the first value."""
        set_room_label_offset("OFFICE", 1.0, 2.0, fresh_db)
        set_room_label_offset("OFFICE", 3.0, 4.0, fresh_db)
        offsets = get_room_label_offsets(fresh_db)
        assert abs(offsets["OFFICE"][0] - 3.0) < 1e-9
        assert abs(offsets["OFFICE"][1] - 4.0) < 1e-9


# ── Shapes Table ───────────────────────────────────────────────────

class TestShapesTable:
    def test_shape_count(self, fresh_db):
        """Three shapes seeded: toilet, bath_sink, dining_table."""
        shapes = get_shapes(fresh_db)
        assert len(shapes) == 3
        names = {s["name"] for s in shapes}
        assert names == {"toilet", "bath_sink", "dining_table"}

    def test_get_shape_by_name(self, fresh_db):
        """get_shape returns the toilet shape with required fields."""
        shape = get_shape("toilet", fresh_db)
        assert shape is not None
        assert "poly_json" in shape
        assert "scale" in shape
        assert "origin" in shape
        assert shape["origin"] == "center"

    def test_get_nonexistent_shape(self, fresh_db):
        """get_shape returns None for unknown names."""
        assert get_shape("nonexistent", fresh_db) is None

    def test_shape_poly_valid_json(self, fresh_db):
        """Each shape's poly_json parses to a list of [x,y] pairs."""
        shapes = get_shapes(fresh_db)
        for s in shapes:
            pts = json.loads(s["poly_json"])
            assert isinstance(pts, list), f"{s['name']} poly is not a list"
            assert len(pts) >= 3, f"{s['name']} poly has < 3 points"
            for pt in pts:
                assert len(pt) == 2, f"{s['name']} point {pt} is not [x,y]"


# ── Phase 14-A  F2 Origin in DB ───────────────────────────────────────

class TestF2OriginConstants:
    def test_f2_constants_seeded(self, fresh_db):
        """F2_EASTING and F2_NORTHING are seeded in the constants table."""
        cd = get_constants_dict(fresh_db)
        assert "F2_EASTING" in cd, "F2_EASTING not seeded"
        assert "F2_NORTHING" in cd, "F2_NORTHING not seeded"

    def test_f2_default_values(self, fresh_db):
        """F2 constants have correct default values."""
        cd = get_constants_dict(fresh_db)
        assert cd["F2_EASTING"] == -18.5
        assert cd["F2_NORTHING"] == -13.5

    def test_f2_category_is_geometry(self, fresh_db):
        """F2 constants are in the geometry category."""
        constants = get_all_constants(fresh_db)
        f2_consts = {c["name"]: c for c in constants
                     if c["name"].startswith("F2_")}
        for name in ("F2_EASTING", "F2_NORTHING"):
            assert f2_consts[name]["category"] == "geometry"

    def test_f2_geometry_unchanged_with_defaults(self, fresh_db):
        """Geometry is identical whether F2 comes from DB or hardcoded."""
        from app.engine import compute_geometry
        from app.database import get_outline_chain
        cd = get_constants_dict(fresh_db)
        chain_rows = get_outline_chain(fresh_db)
        geo = compute_geometry(cd, chain_rows=chain_rows, db_path=fresh_db)
        # F2 position should match the default values
        f2 = geo["points"]["F2"]
        assert abs(f2[0] - (-18.5)) < 1e-6
        # F2_N = -13.5 + R_a1; check it's reasonable (positive, near 0)
        assert abs(f2[1]) < 50  # sanity check

    def test_f2_geometry_changes_with_edited_values(self, fresh_db):
        """Editing F2 constants changes geometry output."""
        from app.engine import compute_geometry
        from app.database import get_outline_chain
        cd = get_constants_dict(fresh_db)
        chain_rows = get_outline_chain(fresh_db)

        # Baseline
        geo_base = compute_geometry(cd, chain_rows=chain_rows, db_path=fresh_db)
        f2_base = geo_base["points"]["F2"]

        # Shift F2 easting by 1 ft
        cd_modified = dict(cd)
        cd_modified["F2_EASTING"] = -17.5
        geo_mod = compute_geometry(cd_modified, chain_rows=chain_rows, db_path=fresh_db)
        f2_mod = geo_mod["points"]["F2"]

        assert abs(f2_mod[0] - f2_base[0] - 1.0) < 1e-6


# ── Phase 14-B  Survey Data in DB ──────────────────────────────────────

class TestSurveyData:
    def test_survey_legs_seeded(self, fresh_db):
        """5 survey legs are seeded in order."""
        legs = get_survey_legs(fresh_db)
        assert len(legs) == 5
        assert [l["seq"] for l in legs] == [1, 2, 3, 4, 5]

    def test_survey_leg_labels(self, fresh_db):
        """Leg labels match traverse endpoints."""
        legs = get_survey_legs(fresh_db)
        labels = [l["label"] for l in legs]
        assert labels == ["P2", "P3", "P4", "P5", "POB"]

    def test_survey_config_seeded(self, fresh_db):
        """Survey config contains required keys."""
        config = get_survey_config(fresh_db)
        expected = {"FC_IN_P3_E", "FC_IN_P3_N", "COORD_ROTATION",
                    "P3_EASTING_OVERRIDE", "P2_P3_NORTHING_OFFSET"}
        assert expected <= set(config.keys())

    def test_survey_config_values(self, fresh_db):
        """Survey config has correct default values."""
        config = get_survey_config(fresh_db)
        assert abs(config["FC_IN_P3_E"] - 18.5141152720) < 1e-10
        assert abs(config["FC_IN_P3_N"] - 13.3968094375) < 1e-10
        assert abs(config["COORD_ROTATION"] - 0.0015153784) < 1e-10

    def test_update_survey_leg(self, fresh_db):
        """Updating a leg persists the change."""
        update_survey_leg(1, {"bearing_deg": 260}, fresh_db)
        legs = get_survey_legs(fresh_db)
        assert legs[0]["bearing_deg"] == 260

    def test_update_survey_config(self, fresh_db):
        """Updating survey config persists the change."""
        update_survey_config("FC_IN_P3_E", 20.0, fresh_db)
        config = get_survey_config(fresh_db)
        assert config["FC_IN_P3_E"] == 20.0

    def test_reset_survey(self, fresh_db):
        """Reset restores original leg and config values."""
        update_survey_leg(1, {"bearing_deg": 999}, fresh_db)
        update_survey_config("FC_IN_P3_E", 0.0, fresh_db)
        reset_survey(fresh_db)
        legs = get_survey_legs(fresh_db)
        assert legs[0]["bearing_deg"] == 257
        config = get_survey_config(fresh_db)
        assert abs(config["FC_IN_P3_E"] - 18.5141152720) < 1e-10

    def test_traverse_from_db_matches_hardcoded(self, fresh_db):
        """DB-driven traverse produces same results as hardcoded."""
        from shared.survey import compute_traverse
        from app.engine import _compute_traverse_from_db

        hardcoded = compute_traverse()
        from_db = _compute_traverse_from_db(fresh_db)

        for label in ("POB", "P2", "P3", "P4", "P5"):
            hc = hardcoded[label]
            db = from_db[label]
            assert abs(hc[0] - db[0]) < 1e-6, f"{label} easting mismatch"
            assert abs(hc[1] - db[1]) < 1e-6, f"{label} northing mismatch"


# ── Phase 14-D  Export / Import ────────────────────────────────────────

class TestExportImport:
    def test_export_has_required_keys(self, fresh_db):
        """Export contains all required top-level keys."""
        data = export_project(fresh_db)
        required = {"version", "exported_at", "constants", "outline_chain",
                    "elements", "element_formulas", "formula_deps", "doors",
                    "variants", "variant_exclusions", "survey_legs",
                    "survey_config", "plumbing_elements"}
        assert required <= set(data.keys())
        assert data["version"] == 1

    def test_export_constants_nonempty(self, fresh_db):
        """Export contains seeded constants."""
        data = export_project(fresh_db)
        assert len(data["constants"]) >= 140

    def test_export_outline_chain_nonempty(self, fresh_db):
        """Export contains outline chain rows."""
        data = export_project(fresh_db)
        assert len(data["outline_chain"]) >= 18

    def test_round_trip_constants(self, fresh_db):
        """Export→import round-trip preserves constants."""
        data = export_project(fresh_db)
        import_project(data, fresh_db)
        data2 = export_project(fresh_db)
        assert len(data2["constants"]) == len(data["constants"])
        orig = {c["name"]: c["value"] for c in data["constants"]}
        imported = {c["name"]: c["value"] for c in data2["constants"]}
        for name, val in orig.items():
            assert abs(imported[name] - val) < 1e-10, f"{name} value changed"

    def test_round_trip_geometry(self, fresh_db):
        """Export→import round-trip produces identical geometry."""
        from app.engine import compute_geometry

        cd = get_constants_dict(fresh_db)
        chain_rows = get_outline_chain(fresh_db)
        geo1 = compute_geometry(cd, chain_rows=chain_rows, db_path=fresh_db)

        data = export_project(fresh_db)
        import_project(data, fresh_db)

        cd2 = get_constants_dict(fresh_db)
        chain_rows2 = get_outline_chain(fresh_db)
        geo2 = compute_geometry(cd2, chain_rows=chain_rows2, db_path=fresh_db)

        # Compare all point positions
        for name, pt in geo1["points"].items():
            pt2 = geo2["points"].get(name)
            assert pt2 is not None, f"Point {name} missing after import"
            d2 = (pt[0] - pt2[0]) ** 2 + (pt[1] - pt2[1]) ** 2
            assert d2 < 1e-10, f"Point {name} moved: d²={d2}"

    def test_import_invalid_version(self, fresh_db):
        """Import rejects unsupported version."""
        data = export_project(fresh_db)
        data["version"] = 99
        with pytest.raises(ValueError, match="Unsupported"):
            import_project(data, fresh_db)

    def test_import_missing_keys(self, fresh_db):
        """Import rejects data with missing required keys."""
        with pytest.raises(ValueError, match="Missing"):
            import_project({"version": 1}, fresh_db)

    def test_import_cycle_detection(self, fresh_db):
        """Import rejects formula dependency cycles."""
        data = export_project(fresh_db)
        data["formula_deps"].append({
            "element_name": "CYCLE_A", "param_name": "position",
            "dep_type": "element", "dep_name": "CYCLE_B",
        })
        data["formula_deps"].append({
            "element_name": "CYCLE_B", "param_name": "position",
            "dep_type": "element", "dep_name": "CYCLE_A",
        })
        with pytest.raises(ValueError, match="cycle"):
            import_project(data, fresh_db)

    def test_import_no_state_change_on_error(self, fresh_db):
        """Failed import leaves DB unchanged."""
        before = export_project(fresh_db)
        bad_data = {"version": 99}
        try:
            import_project(bad_data, fresh_db)
        except ValueError:
            pass
        after = export_project(fresh_db)
        assert len(after["constants"]) == len(before["constants"])
