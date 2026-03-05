"""Tests for Phase 10b SCAD integration (SCAD-1, SCAD-2, SCAD-3)."""
import os
import sys
import pytest

_PROJECT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _PROJECT not in sys.path:
    sys.path.insert(0, _PROJECT)

from tests.test_zapp_conftest import fresh_db, app_client
from app.database import get_config, set_config, get_all_config, get_views


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
