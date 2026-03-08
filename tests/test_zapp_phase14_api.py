"""Tests for Phase 14 API endpoints and undo — API-35 through API-43.

Covers survey data API (legs, config, reset), project export/import,
and undo/redo for all Phase 14 action types.
"""
import json
import pytest

from tests.test_zapp_conftest import fresh_db, app_client


# ── API-35  GET /api/survey/legs ───────────────────────────────────────

class TestAPI35GetSurveyLegs:
    def test_status_200(self, app_client):
        resp = app_client.get("/api/survey/legs")
        assert resp.status_code == 200

    def test_returns_5_legs(self, app_client):
        data = app_client.get("/api/survey/legs").get_json()
        assert len(data) == 5

    def test_leg_labels(self, app_client):
        data = app_client.get("/api/survey/legs").get_json()
        labels = [l["label"] for l in data]
        assert labels == ["P2", "P3", "P4", "P5", "POB"]

    def test_leg_has_required_keys(self, app_client):
        data = app_client.get("/api/survey/legs").get_json()
        required = {"seq", "bearing_deg", "bearing_min", "bearing_sec",
                    "distance_ft", "distance_inch", "label"}
        assert required <= set(data[0].keys())


# ── API-36  PUT /api/survey/legs/<seq> ─────────────────────────────────

class TestAPI36UpdateSurveyLeg:
    def test_update_bearing(self, app_client):
        resp = app_client.put("/api/survey/legs/1",
                              json={"bearing_deg": 260})
        assert resp.status_code == 200
        assert resp.get_json()["ok"] is True
        legs = app_client.get("/api/survey/legs").get_json()
        assert legs[0]["bearing_deg"] == 260

    def test_undo_restores_leg(self, app_client):
        # Modify
        app_client.put("/api/survey/legs/1", json={"bearing_deg": 260})
        # Undo
        resp = app_client.post("/api/undo")
        assert resp.status_code == 200
        # Verify restored
        legs = app_client.get("/api/survey/legs").get_json()
        assert legs[0]["bearing_deg"] == 257

    def test_redo_reapplies_leg(self, app_client):
        app_client.put("/api/survey/legs/1", json={"bearing_deg": 260})
        app_client.post("/api/undo")
        resp = app_client.post("/api/redo")
        assert resp.status_code == 200
        legs = app_client.get("/api/survey/legs").get_json()
        assert legs[0]["bearing_deg"] == 260


# ── API-37  GET /api/survey/config ─────────────────────────────────────

class TestAPI37GetSurveyConfig:
    def test_status_200(self, app_client):
        resp = app_client.get("/api/survey/config")
        assert resp.status_code == 200

    def test_has_required_keys(self, app_client):
        data = app_client.get("/api/survey/config").get_json()
        expected = {"FC_IN_P3_E", "FC_IN_P3_N", "COORD_ROTATION",
                    "P3_EASTING_OVERRIDE", "P2_P3_NORTHING_OFFSET"}
        assert expected <= set(data.keys())

    def test_config_values(self, app_client):
        data = app_client.get("/api/survey/config").get_json()
        assert abs(data["FC_IN_P3_E"] - 18.5141152720) < 1e-10


# ── API-38  PUT /api/survey/config/<key> ───────────────────────────────

class TestAPI38UpdateSurveyConfig:
    def test_update_config(self, app_client):
        resp = app_client.put("/api/survey/config/FC_IN_P3_E",
                              json={"value": 20.0})
        assert resp.status_code == 200
        config = app_client.get("/api/survey/config").get_json()
        assert config["FC_IN_P3_E"] == 20.0

    def test_undo_restores_config(self, app_client):
        app_client.put("/api/survey/config/FC_IN_P3_E",
                       json={"value": 20.0})
        app_client.post("/api/undo")
        config = app_client.get("/api/survey/config").get_json()
        assert abs(config["FC_IN_P3_E"] - 18.5141152720) < 1e-10

    def test_redo_reapplies_config(self, app_client):
        app_client.put("/api/survey/config/FC_IN_P3_E",
                       json={"value": 20.0})
        app_client.post("/api/undo")
        app_client.post("/api/redo")
        config = app_client.get("/api/survey/config").get_json()
        assert config["FC_IN_P3_E"] == 20.0


# ── API-39  POST /api/survey/reset ─────────────────────────────────────

class TestAPI39ResetSurvey:
    def test_reset_restores_defaults(self, app_client):
        app_client.put("/api/survey/legs/1", json={"bearing_deg": 999})
        app_client.put("/api/survey/config/FC_IN_P3_E",
                       json={"value": 0.0})
        resp = app_client.post("/api/survey/reset")
        assert resp.status_code == 200
        legs = app_client.get("/api/survey/legs").get_json()
        assert legs[0]["bearing_deg"] == 257
        config = app_client.get("/api/survey/config").get_json()
        assert abs(config["FC_IN_P3_E"] - 18.5141152720) < 1e-10

    def test_undo_reset(self, app_client):
        # Modify, then reset, then undo reset
        app_client.put("/api/survey/legs/1", json={"bearing_deg": 999})
        app_client.post("/api/survey/reset")
        app_client.post("/api/undo")
        legs = app_client.get("/api/survey/legs").get_json()
        assert legs[0]["bearing_deg"] == 999


# ── API-40  GET /api/project/export ────────────────────────────────────

class TestAPI40ProjectExport:
    def test_status_200(self, app_client):
        resp = app_client.get("/api/project/export")
        assert resp.status_code == 200

    def test_has_required_keys(self, app_client):
        data = app_client.get("/api/project/export").get_json()
        required = {"version", "exported_at", "constants", "outline_chain",
                    "elements", "element_formulas", "formula_deps", "doors",
                    "variants", "variant_exclusions", "survey_legs",
                    "survey_config", "plumbing_elements"}
        assert required <= set(data.keys())

    def test_version_is_1(self, app_client):
        data = app_client.get("/api/project/export").get_json()
        assert data["version"] == 1

    def test_constants_nonempty(self, app_client):
        data = app_client.get("/api/project/export").get_json()
        assert len(data["constants"]) >= 140


# ── API-41  POST /api/project/import ───────────────────────────────────

class TestAPI41ProjectImport:
    def test_round_trip(self, app_client):
        data = app_client.get("/api/project/export").get_json()
        resp = app_client.post("/api/project/import", json=data)
        assert resp.status_code == 200
        assert resp.get_json()["ok"] is True

    def test_bad_version(self, app_client):
        data = app_client.get("/api/project/export").get_json()
        data["version"] = 99
        resp = app_client.post("/api/project/import", json=data)
        assert resp.status_code == 400
        assert "Unsupported" in resp.get_json()["error"]

    def test_missing_keys(self, app_client):
        resp = app_client.post("/api/project/import",
                               json={"version": 1})
        assert resp.status_code == 400
        assert "Missing" in resp.get_json()["error"]

    def test_cycle_detection(self, app_client):
        data = app_client.get("/api/project/export").get_json()
        data["formula_deps"].extend([
            {"element_name": "CYC_A", "param_name": "position",
             "dep_type": "element", "dep_name": "CYC_B"},
            {"element_name": "CYC_B", "param_name": "position",
             "dep_type": "element", "dep_name": "CYC_A"},
        ])
        resp = app_client.post("/api/project/import", json=data)
        assert resp.status_code == 400
        assert "cycle" in resp.get_json()["error"]

    def test_undo_import(self, app_client):
        # Get baseline constant value
        data = app_client.get("/api/project/export").get_json()
        orig_count = len(data["constants"])
        # Modify a constant then export
        app_client.put("/api/constants/WALL_OUTER",
                       json={"value": 1.0})
        mod_data = app_client.get("/api/project/export").get_json()
        # Reset to get clean state, clear undo
        app_client.post("/api/constants/reset")
        # Import modified state
        app_client.post("/api/project/import", json=mod_data)
        # Verify import applied
        consts = app_client.get("/api/constants").get_json()
        wall_outer = next(c for c in consts if c["name"] == "WALL_OUTER")
        assert wall_outer["value"] == 1.0
        # Undo import
        app_client.post("/api/undo")
        consts = app_client.get("/api/constants").get_json()
        wall_outer = next(c for c in consts if c["name"] == "WALL_OUTER")
        assert wall_outer["value"] != 1.0

    def test_geometry_unchanged_after_round_trip(self, app_client):
        """Export→import produces identical geometry."""
        geo1 = app_client.get("/api/geometry").get_json()
        data = app_client.get("/api/project/export").get_json()
        app_client.post("/api/project/import", json=data)
        geo2 = app_client.get("/api/geometry").get_json()
        for name in geo1["points"]:
            p1 = geo1["points"][name]
            p2 = geo2["points"][name]
            d2 = (p1[0] - p2[0]) ** 2 + (p1[1] - p2[1]) ** 2
            assert d2 < 1e-10, f"Point {name} moved: d²={d2}"


# ── Full reset includes survey ─────────────────────────────────────────

class TestFullResetIncludesSurvey:
    def test_reset_restores_survey(self, app_client):
        """Full reset also resets survey data."""
        app_client.put("/api/survey/legs/1", json={"bearing_deg": 999})
        app_client.post("/api/constants/reset")
        legs = app_client.get("/api/survey/legs").get_json()
        assert legs[0]["bearing_deg"] == 257

    def test_undo_reset_restores_survey_edits(self, app_client):
        """Undoing full reset restores edited survey data."""
        app_client.put("/api/survey/legs/1", json={"bearing_deg": 999})
        app_client.post("/api/constants/reset")
        # Undo the full reset
        app_client.post("/api/undo")
        legs = app_client.get("/api/survey/legs").get_json()
        assert legs[0]["bearing_deg"] == 999
