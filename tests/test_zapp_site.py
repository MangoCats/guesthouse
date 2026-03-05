"""Tests for Phase 10c Site Plan integration (SITE-1 through SITE-4)."""
import os
import sys
import pytest

_PROJECT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _PROJECT not in sys.path:
    sys.path.insert(0, _PROJECT)

from tests.test_zapp_conftest import fresh_db, app_client
from app.database import get_config, set_config, get_views


# =========================================================================
# TestSiteConfig (SITE-1)
# =========================================================================

class TestSiteConfig:
    """Setback config values for site plan placement."""

    def test_default_setback_216(self, fresh_db):
        assert get_config("setback_216", fresh_db) == "11.0"

    def test_default_setback_275(self, fresh_db):
        assert get_config("setback_275", fresh_db) == "25.5"

    def test_set_setback(self, fresh_db):
        set_config("setback_216", "12.0", fresh_db)
        assert get_config("setback_216", fresh_db) == "12.0"

    def test_config_api_setback(self, app_client):
        resp = app_client.get("/api/config/setback_216")
        assert resp.status_code == 200
        data = resp.get_json()
        assert data["value"] == "11.0"

    def test_config_api_set_setback(self, app_client):
        resp = app_client.put("/api/config/setback_216",
                              json={"value": "12.5"})
        assert resp.status_code == 200
        resp2 = app_client.get("/api/config/setback_216")
        assert resp2.get_json()["value"] == "12.5"


# =========================================================================
# TestSiteViews
# =========================================================================

class TestSiteViews:
    """Site plan views registered and accessible."""

    def test_site_views_registered(self, fresh_db):
        views = get_views(fresh_db)
        names = {v["name"] for v in views}
        assert "site_plan_df" in names
        assert "site_plan_fs" in names

    def test_site_views_category(self, fresh_db):
        views = get_views(fresh_db)
        site_views = [v for v in views if v["category"] == "site"]
        assert len(site_views) == 2

    def test_pdf_file_endpoint(self, app_client):
        """PDF endpoint returns correct content type when file exists."""
        resp = app_client.get("/api/svg/site_plan_df/file")
        # File may not exist in test env, but endpoint should exist
        assert resp.status_code in (200, 404)
        if resp.status_code == 200:
            assert "pdf" in resp.content_type


# =========================================================================
# TestSurveyPoints (SITE-4)
# =========================================================================

class TestSurveyPoints:
    """Survey points API endpoint."""

    def test_survey_points_endpoint(self, app_client):
        resp = app_client.get("/api/survey-points")
        assert resp.status_code == 200
        data = resp.get_json()
        assert "points" in data
        assert "distances" in data
        assert "arc_radii" in data

    def test_survey_point_names(self, app_client):
        resp = app_client.get("/api/survey-points")
        data = resp.get_json()
        names = set(data["points"].keys())
        assert {"POB", "P2", "P3", "P4", "P5"}.issubset(names)

    def test_survey_distances(self, app_client):
        resp = app_client.get("/api/survey-points")
        data = resp.get_json()
        assert len(data["distances"]) == 5
        for d in data["distances"]:
            assert "from" in d and "to" in d and "distance" in d
            assert d["distance"] > 0

    def test_survey_arc_radii(self, app_client):
        resp = app_client.get("/api/survey-points")
        data = resp.get_json()
        radii = data["arc_radii"]
        assert radii["R1"] == 10.0
        assert radii["R2"] == 12.5
        assert radii["R3"] == 11.0

    def test_survey_points_coordinates(self, app_client):
        """P-series points have reasonable coordinates."""
        resp = app_client.get("/api/survey-points")
        data = resp.get_json()
        # P3 should be near the origin (offset by FC_IN_P3)
        p3 = data["points"]["P3"]
        assert len(p3) == 2
        # All points should have coordinates in reasonable range
        for name, pt in data["points"].items():
            assert -50 < pt[0] < 50, f"{name} easting out of range"
            assert -50 < pt[1] < 50, f"{name} northing out of range"


# =========================================================================
# TestSiteElements (SITE-2)
# =========================================================================

class TestSiteElements:
    """Site element CRUD for drainfield and annotations."""

    def test_create_drainfield(self, app_client):
        resp = app_client.post("/api/elements", json={
            "type": "site_element",
            "name": "test_drainfield",
            "properties": {
                "subtype": "drainfield",
                "width": 25.0, "height": 10.0,
                "label": "NEW DRAINFIELD",
                "x": 0, "y": 0, "rotation": 0,
            },
        })
        assert resp.status_code == 201
        data = resp.get_json()
        assert data["type"] == "site_element"
        assert data["name"] == "test_drainfield"

    def test_delete_drainfield(self, app_client):
        resp = app_client.post("/api/elements", json={
            "type": "site_element",
            "name": "test_df_del",
            "properties": {"subtype": "drainfield"},
        })
        eid = resp.get_json()["id"]
        resp2 = app_client.delete(f"/api/elements/{eid}")
        assert resp2.status_code == 200

    def test_create_annotation(self, app_client):
        resp = app_client.post("/api/elements", json={
            "type": "site_annotation",
            "name": "test_note",
            "properties": {
                "text": "PROPOSED ADU",
                "style": "text",
                "x": 0, "y": 0,
            },
        })
        assert resp.status_code == 201
        data = resp.get_json()
        assert data["type"] == "site_annotation"

    def test_update_annotation(self, app_client):
        resp = app_client.post("/api/elements", json={
            "type": "site_annotation",
            "name": "test_note_upd",
            "properties": {"text": "OLD TEXT"},
        })
        eid = resp.get_json()["id"]
        resp2 = app_client.put(f"/api/elements/{eid}", json={
            "properties": {"text": "NEW TEXT"},
        })
        assert resp2.status_code == 200


# =========================================================================
# TestGenerateSitePlan
# =========================================================================

class TestGenerateSitePlan:
    """POST /api/generate-site-plan endpoint."""

    def test_generate_endpoint(self, app_client):
        resp = app_client.post("/api/generate-site-plan")
        assert resp.status_code == 200
        data = resp.get_json()
        assert data["ok"] is True
        assert data["setback_216"] == "11.0"
        assert data["setback_275"] == "25.5"
