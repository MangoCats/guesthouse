"""Tests for formula REST API endpoints (Phase 12b).

Covers: GET/PUT/DELETE /api/formulas, lock/unlock, deps, dependents.
"""
import json
import pytest

from tests.test_zapp_conftest import fresh_db, app_client


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _create_wall(client, name="TEST_IW"):
    """Create a wall element via the elements API."""
    return client.post(
        "/api/elements",
        data=json.dumps({"type": "wall", "name": name, "properties": {}}),
        content_type="application/json",
    )


def _put_formula(client, elem_name, param_name="position", formula=None):
    """PUT a formula via the API."""
    if formula is None:
        formula = {
            "type": "wall_rect",
            "anchor": [0, 0],
            "along": "east",
            "thickness_dir": "north",
            "thickness": 0.5,
            "end_mode": "fixed",
            "length": 5.0,
        }
    return client.put(
        f"/api/formulas/{elem_name}/{param_name}",
        data=json.dumps({"formula": formula}),
        content_type="application/json",
    )


# ---------------------------------------------------------------------------
# GET /api/formulas
# ---------------------------------------------------------------------------

class TestGetAllFormulas:
    def test_has_seeded_formulas(self, app_client):
        resp = app_client.get("/api/formulas")
        assert resp.status_code == 200
        data = resp.get_json()
        # 13 IW + 5 items + 12 outer openings + 7 rough openings = 37
        names = {f["element_name"] for f in data}
        assert "IW1" in names
        assert "DRYER" in names
        assert "O3" in names
        assert "RO1" in names
        assert len(names) >= 37

    def test_returns_formulas_after_upsert(self, app_client):
        _create_wall(app_client, "FW1")
        _put_formula(app_client, "FW1")
        data = app_client.get("/api/formulas").get_json()
        assert len(data) >= 1
        assert any(f["element_name"] == "FW1" for f in data)


# ---------------------------------------------------------------------------
# GET /api/formulas/<element_name>
# ---------------------------------------------------------------------------

class TestGetElementFormulas:
    def test_empty_for_no_formula(self, app_client):
        _create_wall(app_client, "FW2")
        resp = app_client.get("/api/formulas/FW2")
        assert resp.status_code == 200
        assert resp.get_json() == []

    def test_returns_formula(self, app_client):
        _create_wall(app_client, "FW3")
        _put_formula(app_client, "FW3")
        data = app_client.get("/api/formulas/FW3").get_json()
        assert len(data) == 1
        assert data[0]["element_name"] == "FW3"
        assert data[0]["param_name"] == "position"


# ---------------------------------------------------------------------------
# PUT /api/formulas/<element_name>/<param_name>
# ---------------------------------------------------------------------------

class TestPutFormula:
    def test_create_formula(self, app_client):
        _create_wall(app_client, "FW4")
        resp = _put_formula(app_client, "FW4")
        assert resp.status_code == 200
        data = resp.get_json()
        assert data["element_name"] == "FW4"
        assert data["param_name"] == "position"

    def test_update_formula(self, app_client):
        _create_wall(app_client, "FW5")
        _put_formula(app_client, "FW5")
        updated = {
            "type": "wall_rect",
            "anchor": [1, 1],
            "along": "north",
            "thickness_dir": "east",
            "thickness": 0.33,
            "end_mode": "fixed",
            "length": 8.0,
        }
        resp = _put_formula(app_client, "FW5", formula=updated)
        assert resp.status_code == 200
        stored = json.loads(resp.get_json()["formula_json"])
        assert stored["length"] == 8.0

    def test_missing_formula_returns_400(self, app_client):
        _create_wall(app_client, "FW6")
        resp = app_client.put(
            "/api/formulas/FW6/position",
            data=json.dumps({}),
            content_type="application/json",
        )
        assert resp.status_code == 400

    def test_creates_deps(self, app_client):
        _create_wall(app_client, "FW7")
        formula = {
            "type": "item_rect",
            "anchor": {"element": "OTHER", "corner": "NE"},
            "along": "east",
            "across": "north",
            "width": {"const": "BED_WIDTH"},
            "depth": 2.0,
            "anchor_corner": "sw",
        }
        _put_formula(app_client, "FW7", formula=formula)
        deps = app_client.get("/api/formulas/FW7/deps").get_json()
        dep_set = {(d["dep_type"], d["dep_name"]) for d in deps}
        assert ("element", "OTHER") in dep_set
        assert ("constant", "BED_WIDTH") in dep_set


# ---------------------------------------------------------------------------
# DELETE /api/formulas/<element_name>/<param_name>
# ---------------------------------------------------------------------------

class TestDeleteFormula:
    def test_delete_existing(self, app_client):
        _create_wall(app_client, "FW8")
        _put_formula(app_client, "FW8")
        resp = app_client.delete("/api/formulas/FW8/position")
        assert resp.status_code == 200
        assert resp.get_json()["ok"] is True
        # Verify gone
        data = app_client.get("/api/formulas/FW8").get_json()
        assert len(data) == 0

    def test_delete_nonexistent_returns_404(self, app_client):
        resp = app_client.delete("/api/formulas/NOSUCH/position")
        assert resp.status_code == 404


# ---------------------------------------------------------------------------
# POST /api/formulas/<element_name>/<param_name>/lock
# ---------------------------------------------------------------------------

class TestLockFormula:
    def test_lock(self, app_client):
        _create_wall(app_client, "FW9")
        _put_formula(app_client, "FW9")
        locked_val = {"poly": [[0, 0], [5, 0], [5, 0.5], [0, 0.5]]}
        resp = app_client.post(
            "/api/formulas/FW9/position/lock",
            data=json.dumps({"locked": True, "locked_value": locked_val}),
            content_type="application/json",
        )
        assert resp.status_code == 200
        # Verify locked state
        data = app_client.get("/api/formulas/FW9").get_json()
        assert data[0]["locked"] == 1

    def test_unlock(self, app_client):
        _create_wall(app_client, "FW10")
        _put_formula(app_client, "FW10")
        # Lock then unlock
        app_client.post(
            "/api/formulas/FW10/position/lock",
            data=json.dumps({"locked": True}),
            content_type="application/json",
        )
        resp = app_client.post(
            "/api/formulas/FW10/position/lock",
            data=json.dumps({"locked": False}),
            content_type="application/json",
        )
        assert resp.status_code == 200
        data = app_client.get("/api/formulas/FW10").get_json()
        assert data[0]["locked"] == 0

    def test_lock_nonexistent_returns_404(self, app_client):
        resp = app_client.post(
            "/api/formulas/NOSUCH/position/lock",
            data=json.dumps({"locked": True}),
            content_type="application/json",
        )
        assert resp.status_code == 404


# ---------------------------------------------------------------------------
# GET /api/formulas/<element_name>/deps
# ---------------------------------------------------------------------------

class TestFormulaDeps:
    def test_no_deps_for_literal_formula(self, app_client):
        _create_wall(app_client, "FW11")
        _put_formula(app_client, "FW11")
        deps = app_client.get("/api/formulas/FW11/deps").get_json()
        # wall_rect with literal anchor [0,0] has no element/constant deps
        dep_types = {d["dep_type"] for d in deps}
        assert "element" not in dep_types


# ---------------------------------------------------------------------------
# GET /api/formulas/<element_name>/dependents
# ---------------------------------------------------------------------------

class TestFormulaDependents:
    def test_dependents(self, app_client):
        _create_wall(app_client, "BASE_W")
        _create_wall(app_client, "DEP_W")
        _put_formula(app_client, "BASE_W")
        # DEP_W depends on BASE_W
        dep_formula = {
            "type": "item_rect",
            "anchor": {"element": "BASE_W", "corner": "NE"},
            "along": "east",
            "across": "north",
            "width": 2.0,
            "depth": 1.0,
            "anchor_corner": "sw",
        }
        _put_formula(app_client, "DEP_W", formula=dep_formula)
        resp = app_client.get("/api/formulas/BASE_W/dependents")
        assert resp.status_code == 200
        names = {d["element_name"] for d in resp.get_json()}
        assert "DEP_W" in names
