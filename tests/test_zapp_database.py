"""Tests for app.database — DB-1 through DB-8.

Verifies schema initialisation, seeding, and CRUD operations using an
isolated temporary database per test.
"""
import os
import sqlite3
import pytest
from app.database import (
    init_db, get_db, get_all_constants, get_constants_dict,
    update_constant, update_constants_batch, get_categories,
    get_outline_chain, get_views, reset_constants,
)

# Re-use the fresh_db fixture from test_zapp_conftest.py
from tests.test_zapp_conftest import fresh_db


# ── DB-1  Schema Initialisation ─────────────────────────────────────

class TestDB1SchemaInit:
    def test_tables_created(self, fresh_db):
        """DB file contains the three core tables."""
        with get_db(fresh_db) as conn:
            rows = conn.execute(
                "SELECT name FROM sqlite_master WHERE type='table' ORDER BY name"
            ).fetchall()
            names = {r["name"] for r in rows}
        assert {"constants", "outline_chain", "views"} <= names

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
