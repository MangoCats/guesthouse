"""Tests for elements and openings (Phase 3).

Covers: DB-9, API-20–22, API-24–26, DT-9–11, UNDO-4 (element types).
"""
import json

import pytest

from app.database import (
    get_db, init_db, get_all_elements, get_element, get_element_by_name,
    create_element, update_element, delete_element,
)
from app.elements import (
    get_elements_for_variant, get_controlling_constant, get_hosted_openings,
    IW_CONSTANT_MAP,
)
from app.undo import UndoManager

# Re-use fixtures from test_zapp_conftest.py
from tests.test_zapp_conftest import fresh_db, app_client  # noqa: F401


# ── DB-9: elements table ──────────────────────────────────────────────

class TestDB9Elements:
    """DB-9: The database SHALL maintain an elements table."""

    def test_elements_table_exists(self, fresh_db):
        with get_db(fresh_db) as conn:
            names = {r[0] for r in conn.execute(
                "SELECT name FROM sqlite_master WHERE type='table'"
            ).fetchall()}
        assert "elements" in names

    def test_13_walls_seeded(self, fresh_db):
        with get_db(fresh_db) as conn:
            count = conn.execute(
                "SELECT count(*) FROM elements WHERE type = 'wall'"
            ).fetchone()[0]
        assert count == 13

    def test_only_walls_and_labels_seeded(self, fresh_db):
        with get_db(fresh_db) as conn:
            count = conn.execute(
                "SELECT count(*) FROM elements WHERE type NOT IN ('wall', 'label')"
            ).fetchone()[0]
        assert count == 0

    def test_11_room_labels_seeded(self, fresh_db):
        with get_db(fresh_db) as conn:
            count = conn.execute(
                "SELECT count(*) FROM elements WHERE type = 'label'"
            ).fetchone()[0]
        assert count == 11

    def test_element_properties_json_valid(self, fresh_db):
        elements = get_all_elements(fresh_db)
        for e in elements:
            props = json.loads(e["properties"]) if isinstance(e["properties"], str) else e["properties"]
            assert isinstance(props, dict)
            if e["type"] == "wall":
                assert "thickness_constant" in props
                assert "orientation" in props


# ── Element CRUD ──────────────────────────────────────────────────────

class TestElementCRUD:
    """Unit tests for element CRUD operations."""

    def test_create_element(self, fresh_db):
        rec = create_element("furniture", "HAMPER", {"width": 31.5}, None, fresh_db)
        assert rec["id"] is not None
        assert rec["type"] == "furniture"
        assert rec["name"] == "HAMPER"

    def test_get_element_by_id(self, fresh_db):
        rec = create_element("furniture", "TEST_ITEM", {}, None, fresh_db)
        fetched = get_element(rec["id"], fresh_db)
        assert fetched is not None
        assert fetched["name"] == "TEST_ITEM"

    def test_update_element(self, fresh_db):
        rec = create_element("furniture", "UPD_ITEM", {"x": 1}, None, fresh_db)
        updated = update_element(rec["id"], {"properties": {"x": 2}}, fresh_db)
        props = json.loads(updated["properties"]) if isinstance(updated["properties"], str) else updated["properties"]
        assert props["x"] == 2

    def test_delete_element(self, fresh_db):
        rec = create_element("furniture", "DEL_ITEM", {}, None, fresh_db)
        deleted = delete_element(rec["id"], fresh_db)
        assert rec["id"] in deleted
        assert get_element(rec["id"], fresh_db) is None

    def test_delete_wall_cascades_opening(self, fresh_db):
        wall = create_element("wall", "TEST_WALL", {}, None, fresh_db)
        opening = create_element(
            "opening", "TEST_OPENING",
            {"host_wall": "TEST_WALL", "width": 36},
            None, fresh_db,
        )
        deleted = delete_element(wall["id"], fresh_db)
        assert wall["id"] in deleted
        assert opening["id"] in deleted
        assert get_element(opening["id"], fresh_db) is None


# ── Element business logic ────────────────────────────────────────────

class TestElementBusinessLogic:

    def test_get_elements_for_variant_all(self, fresh_db):
        # 13 IW walls + 11 room labels have variant=NULL, visible to all variants
        elems = get_elements_for_variant("standard", fresh_db)
        assert len(elems) == 24  # 13 walls + 11 labels

    def test_get_elements_for_variant_specific(self, fresh_db):
        create_element("furniture", "VARIANT_ITEM", {}, "standard", fresh_db)
        elems = get_elements_for_variant("standard", fresh_db)
        assert any(e["name"] == "VARIANT_ITEM" for e in elems)
        elems2 = get_elements_for_variant("minik", fresh_db)
        assert not any(e["name"] == "VARIANT_ITEM" for e in elems2)

    def test_controlling_constant(self):
        assert get_controlling_constant("IW1") == "IW1_OFFSET_FROM_W9"
        assert get_controlling_constant("IW2O") is None

    def test_hosted_openings(self):
        assert "RO1" in get_hosted_openings("IW1")
        assert get_hosted_openings("IW8") == []
        assert len(get_hosted_openings("IW9")) == 2


# ── API-20: POST /api/elements ────────────────────────────────────────

class TestAPI20CreateElement:
    """API-20: POST /api/elements SHALL create a new element."""

    def test_create_furniture_201(self, app_client):
        resp = app_client.post("/api/elements", json={
            "type": "furniture", "name": "HAMPER",
            "properties": {"width": 31.5, "depth": 19},
        })
        assert resp.status_code == 201
        data = resp.get_json()
        assert data["id"] is not None
        assert data["name"] == "HAMPER"

    def test_create_missing_type_400(self, app_client):
        resp = app_client.post("/api/elements", json={"name": "X"})
        assert resp.status_code == 400


# ── API-21: PUT /api/elements/<id> ────────────────────────────────────

class TestAPI21UpdateElement:
    """API-21: PUT /api/elements/<id> SHALL update element properties."""

    def test_update_properties_200(self, app_client):
        resp = app_client.post("/api/elements", json={
            "type": "furniture", "name": "UPD_TEST",
            "properties": {"x": 1},
        })
        eid = resp.get_json()["id"]
        resp = app_client.put(f"/api/elements/{eid}", json={
            "properties": {"x": 5.0, "y": -3.0},
        })
        assert resp.status_code == 200
        props = json.loads(resp.get_json()["properties"])
        assert props["x"] == 5.0

    def test_update_nonexistent_404(self, app_client):
        resp = app_client.put("/api/elements/99999", json={"properties": {}})
        assert resp.status_code == 404


# ── API-22: DELETE /api/elements/<id> ─────────────────────────────────

class TestAPI22DeleteElement:
    """API-22: DELETE /api/elements/<id> SHALL remove element and dependents."""

    def test_delete_element_200(self, app_client):
        resp = app_client.post("/api/elements", json={
            "type": "furniture", "name": "DEL_TEST",
        })
        eid = resp.get_json()["id"]
        resp = app_client.delete(f"/api/elements/{eid}")
        assert resp.status_code == 200
        assert eid in resp.get_json()["deleted"]

    def test_delete_cascade_removes_opening(self, app_client):
        # Create a wall
        resp = app_client.post("/api/elements", json={
            "type": "wall", "name": "CWALL",
        })
        wall_id = resp.get_json()["id"]
        # Create an opening hosted by that wall
        resp = app_client.post("/api/elements", json={
            "type": "opening", "name": "COPENING",
            "properties": {"host_wall": "CWALL", "width": 36},
        })
        op_id = resp.get_json()["id"]
        # Delete the wall — should cascade
        resp = app_client.delete(f"/api/elements/{wall_id}")
        deleted = resp.get_json()["deleted"]
        assert wall_id in deleted
        assert op_id in deleted


# ── API-24: POST /api/openings ────────────────────────────────────────

class TestAPI24CreateOpening:
    """API-24: POST /api/openings SHALL create a new opening."""

    def test_create_opening_201(self, app_client):
        resp = app_client.post("/api/openings", json={
            "name": "O8a", "segment": "F18-F1", "width": 19, "offset": 48,
        })
        assert resp.status_code == 201
        data = resp.get_json()
        assert data["type"] == "opening"
        assert data["name"] == "O8a"

    def test_create_opening_missing_segment_400(self, app_client):
        resp = app_client.post("/api/openings", json={"name": "X"})
        assert resp.status_code == 400


# ── API-25: PUT /api/openings/<name> ──────────────────────────────────

class TestAPI25UpdateOpening:
    """API-25: PUT /api/openings/<name> SHALL update opening properties."""

    def test_update_opening_width(self, app_client):
        app_client.post("/api/openings", json={
            "name": "TEST_O", "segment": "F1-F2", "width": 20,
        })
        resp = app_client.put("/api/openings/TEST_O", json={"width": 25})
        assert resp.status_code == 200
        props = json.loads(resp.get_json()["properties"])
        assert props["width"] == 25


# ── API-26: DELETE /api/openings/<name> ───────────────────────────────

class TestAPI26DeleteOpening:
    """API-26: DELETE /api/openings/<name> SHALL remove opening and door."""

    def test_delete_opening_200(self, app_client):
        app_client.post("/api/openings", json={
            "name": "DEL_O", "segment": "F1-F2", "width": 20,
        })
        resp = app_client.delete("/api/openings/DEL_O")
        assert resp.status_code == 200

    def test_delete_opening_removes_door(self, app_client):
        app_client.post("/api/openings", json={
            "name": "DOOR_O", "segment": "F1-F2", "width": 20,
        })
        app_client.post("/api/doors", json={
            "opening": "DOOR_O", "hinge_side": "east",
            "swing_direction": "south", "width": 36,
        })
        resp = app_client.delete("/api/openings/DOOR_O")
        assert resp.status_code == 200
        # Door should be gone
        door_resp = app_client.get("/api/doors")
        doors = door_resp.get_json()
        assert not any(d["opening_name"] == "DOOR_O" for d in doors)


# ── Undo integration (UNDO-4) ────────────────────────────────────────

class TestElementUndo:
    """UNDO-4: Undo works across element mutation types."""

    def test_undo_element_create(self, app_client):
        resp = app_client.post("/api/elements", json={
            "type": "furniture", "name": "UNDO_CREATE",
        })
        assert resp.status_code == 201
        eid = resp.get_json()["id"]
        # Undo should remove it
        resp = app_client.post("/api/undo")
        assert resp.status_code == 200
        # Verify removed
        resp = app_client.get("/api/elements")
        names = [e["name"] for e in resp.get_json()]
        assert "UNDO_CREATE" not in names

    def test_undo_element_delete(self, app_client):
        resp = app_client.post("/api/elements", json={
            "type": "furniture", "name": "UNDO_DEL",
        })
        eid = resp.get_json()["id"]
        app_client.delete(f"/api/elements/{eid}")
        # Undo should restore it
        resp = app_client.post("/api/undo")
        assert resp.status_code == 200
        resp = app_client.get("/api/elements")
        names = [e["name"] for e in resp.get_json()]
        assert "UNDO_DEL" in names
