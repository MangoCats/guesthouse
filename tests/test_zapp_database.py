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
        """Variant exclusions table has exactly 4 rows."""
        with get_db(fresh_db) as conn:
            count = conn.execute(
                "SELECT count(*) FROM variant_exclusions"
            ).fetchone()[0]
        assert count == 4

    def test_bare_exclusions(self, fresh_db):
        """Bare variant excludes IW6 (wall) and RO5 (rough_opening)."""
        excl = get_variant_exclusions("bare", fresh_db)
        assert excl == {"wall": {"IW6"}, "rough_opening": {"RO5"}}

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
