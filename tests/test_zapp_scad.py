"""Tests for Phase 10b SCAD integration and database robustness."""
import os
import sqlite3
import sys
import pytest

_PROJECT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _PROJECT not in sys.path:
    sys.path.insert(0, _PROJECT)

from tests.test_zapp_conftest import fresh_db, app_client
from app.database import (
    get_config, set_config, get_all_config, get_views,
    validate_db, reset_db, init_db,
)


# =========================================================================
# TestConfig
# =========================================================================

class TestConfig:
    """Config table CRUD operations."""

    def test_default_roof_style(self, fresh_db):
        assert get_config("roof_style", fresh_db) == "flat"

    def test_set_config(self, fresh_db):
        set_config("roof_style", "2in12", fresh_db)
        assert get_config("roof_style", fresh_db) == "2in12"

    def test_get_unknown_key(self, fresh_db):
        assert get_config("nonexistent_key", fresh_db) is None

    def test_get_all_config(self, fresh_db):
        cfg = get_all_config(fresh_db)
        assert "roof_style" in cfg
        assert cfg["roof_style"] == "flat"

    def test_config_api_get(self, app_client):
        resp = app_client.get("/api/config/roof_style")
        assert resp.status_code == 200
        data = resp.get_json()
        assert data["key"] == "roof_style"
        assert data["value"] == "flat"

    def test_config_api_get_unknown(self, app_client):
        resp = app_client.get("/api/config/no_such_key")
        assert resp.status_code == 404

    def test_config_api_set(self, app_client):
        resp = app_client.put("/api/config/roof_style",
                              json={"value": "2in12"})
        assert resp.status_code == 200
        data = resp.get_json()
        assert data["ok"] is True
        # Verify it persisted
        resp2 = app_client.get("/api/config/roof_style")
        assert resp2.get_json()["value"] == "2in12"

    def test_config_api_set_missing_value(self, app_client):
        resp = app_client.put("/api/config/roof_style", json={})
        assert resp.status_code == 400


# =========================================================================
# TestScadViews
# =========================================================================

class TestScadViews:
    """SCAD views registered in database."""

    def test_scad_views_registered(self, fresh_db):
        views = get_views(fresh_db)
        names = {v["name"] for v in views}
        assert "3d_flat" in names
        assert "3d_2in12" in names
        assert "3views" in names

    def test_scad_views_category(self, fresh_db):
        views = get_views(fresh_db)
        scad = [v for v in views if v["category"] == "3d"]
        assert len(scad) == 3

    def test_views_api_includes_scad(self, app_client):
        resp = app_client.get("/api/views")
        data = resp.get_json()
        names = {v["name"] for v in data}
        assert "3d_flat" in names


# =========================================================================
# TestGenerateScad (SCAD-1)
# =========================================================================

class TestGenerateScad:
    """POST /api/generate-3d endpoint."""

    def test_generate_3d_endpoint(self, app_client):
        resp = app_client.post("/api/generate-3d")
        assert resp.status_code == 200
        data = resp.get_json()
        assert data["ok"] is True
        assert data["roof_style"] == "flat"
        assert "flat_roof" in data["output"]

    def test_generate_3d_2in12(self, app_client):
        # Set roof style to 2in12
        app_client.put("/api/config/roof_style", json={"value": "2in12"})
        resp = app_client.post("/api/generate-3d")
        data = resp.get_json()
        assert data["ok"] is True
        assert data["roof_style"] == "2in12"
        assert "2in12" in data["output"]

    def test_generate_3d_creates_file(self, app_client):
        app_client.post("/api/generate-3d")
        scad_path = os.path.join(_PROJECT, "scad", "flat_roof.scad")
        assert os.path.exists(scad_path)


# =========================================================================
# TestGenerateViews (SCAD-3)
# =========================================================================

class TestGenerateViews:
    """POST /api/generate-views endpoint."""

    def test_generate_views_endpoint(self, app_client):
        resp = app_client.post("/api/generate-views")
        data = resp.get_json()
        # May fail if OpenSCAD not installed — still check response format
        assert "ok" in data
        assert "roof_style" in data or "error" in data


# =========================================================================
# TestValidateDb
# =========================================================================

class TestValidateDb:
    """Database validation checks."""

    def test_healthy_db(self, fresh_db):
        ok, issues = validate_db(fresh_db)
        assert ok is True
        assert issues == []

    def test_missing_db(self, tmp_path):
        ok, issues = validate_db(str(tmp_path / "nonexistent.db"))
        assert ok is False
        assert any("does not exist" in i for i in issues)

    def test_empty_db(self, tmp_path):
        """Empty DB file has no tables."""
        empty_path = str(tmp_path / "empty.db")
        open(empty_path, "w").close()
        ok, issues = validate_db(empty_path)
        assert ok is False
        assert any("missing" in i for i in issues)

    def test_db_missing_table(self, fresh_db):
        """Drop a required table and validate."""
        conn = sqlite3.connect(fresh_db)
        conn.execute("DROP TABLE outline_chain")
        conn.commit()
        conn.close()
        ok, issues = validate_db(fresh_db)
        assert ok is False
        assert any("outline_chain" in i for i in issues)

    def test_db_empty_constants(self, fresh_db):
        """Empty constants table detected."""
        conn = sqlite3.connect(fresh_db)
        conn.execute("DELETE FROM constants")
        conn.commit()
        conn.close()
        ok, issues = validate_db(fresh_db)
        assert ok is False
        assert any("constants" in i for i in issues)

    def test_db_empty_outline(self, fresh_db):
        """Empty outline_chain table detected."""
        conn = sqlite3.connect(fresh_db)
        conn.execute("DELETE FROM outline_chain")
        conn.commit()
        conn.close()
        ok, issues = validate_db(fresh_db)
        assert ok is False
        assert any("outline chain" in i for i in issues)


# =========================================================================
# TestResetDb
# =========================================================================

class TestResetDb:
    """Database reset and recovery."""

    def test_reset_recovers_from_empty(self, tmp_path):
        """Reset creates a fresh healthy DB from an empty file."""
        db_path = str(tmp_path / "corrupt.db")
        open(db_path, "w").close()
        ok, issues = reset_db(db_path)
        assert ok is True
        assert issues == []

    def test_reset_recovers_from_missing_table(self, fresh_db):
        """Reset recovers from a dropped table."""
        conn = sqlite3.connect(fresh_db)
        conn.execute("DROP TABLE outline_chain")
        conn.commit()
        conn.close()
        ok, issues = reset_db(fresh_db)
        assert ok is True
        assert issues == []

    def test_reset_preserves_seed_data(self, tmp_path):
        """After reset, seed data is present."""
        db_path = str(tmp_path / "reset_test.db")
        init_db(db_path)
        reset_db(db_path)
        assert get_config("roof_style", db_path) == "flat"
        views = get_views(db_path)
        assert len(views) > 0

    def test_reset_api(self, app_client):
        """POST /api/reset-database returns ok."""
        resp = app_client.post("/api/reset-database")
        assert resp.status_code == 200
        data = resp.get_json()
        assert data["ok"] is True
        assert data["issues"] == []

    def test_geometry_after_reset(self, app_client):
        """Geometry loads after a database reset."""
        app_client.post("/api/reset-database")
        resp = app_client.get("/api/geometry")
        assert resp.status_code == 200
        data = resp.get_json()
        assert "outline" in data or "points" in data


# =========================================================================
# TestGeometryDbError
# =========================================================================

class TestGeometryDbError:
    """Geometry endpoint provides DB diagnostics on failure."""

    def test_geometry_db_issue_flag(self, tmp_path):
        """When DB is corrupt, geometry response includes db_issue."""
        db_path = str(tmp_path / "bad.db")
        init_db(db_path)
        # Corrupt: empty outline_chain
        conn = sqlite3.connect(db_path)
        conn.execute("DELETE FROM outline_chain")
        conn.commit()
        conn.close()
        from app.server import create_app
        app = create_app(db_path=db_path)
        app.config["TESTING"] = True
        with app.test_client() as client:
            resp = client.get("/api/geometry")
            assert resp.status_code == 500
            data = resp.get_json()
            assert "error" in data
