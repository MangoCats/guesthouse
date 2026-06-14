"""Tests for doors (Phase 3).

Covers: DB-10, API-27–29, SEL-7, RT-5, UNDO-4 (door types).
"""
import json

import pytest

from app.database import (
    get_db, init_db, get_all_doors, get_door, create_door, update_door, delete_door,
    get_element_by_name, get_constants_dict,
)
from floorplan.constants import LOWER_WALL_HEIGHT
from app.doors import validate_door
from app.undo import UndoManager

# Re-use fixtures from test_zapp_conftest.py
from tests.test_zapp_conftest import fresh_db, app_client  # noqa: F401


# ── DB-10: doors table ────────────────────────────────────────────────

class TestDB10Doors:
    """DB-10: The database SHALL maintain a doors table."""

    def test_doors_table_exists(self, fresh_db):
        with get_db(fresh_db) as conn:
            names = {r[0] for r in conn.execute(
                "SELECT name FROM sqlite_master WHERE type='table'"
            ).fetchall()}
        assert "doors" in names

    def test_9_doors_seeded(self, fresh_db):
        doors = get_all_doors(fresh_db)
        assert len(doors) == 9
        names = {d["opening_name"] for d in doors}
        assert names == {"O3", "O6", "RO1", "RO2", "RO3", "RO4", "RO5", "RO6", "RO7"}

    def test_door_defaults_valid(self, fresh_db):
        doors = get_all_doors(fresh_db)
        for d in doors:
            assert d["hinge_side"] in ("east", "west", "north", "south")
            assert d["swing_direction"] in ("east", "west", "north", "south")
            assert d["door_type"] in ("single", "double", "hanging_slider")
            assert d["width"] > 0
        # RO6 and RO7 should be double doors
        ro6 = get_door("RO6", fresh_db)
        assert ro6["door_type"] == "double"
        ro7 = get_door("RO7", fresh_db)
        assert ro7["door_type"] == "double"


# ── Door CRUD ─────────────────────────────────────────────────────────

class TestDoorCRUD:
    """Unit tests for door CRUD operations."""

    def test_create_door(self, fresh_db):
        rec = create_door("TEST_RO", 42, "east", "south", "single", fresh_db)
        assert rec["id"] is not None
        assert rec["opening_name"] == "TEST_RO"
        assert rec["width"] == 42

    def test_get_door_by_opening(self, fresh_db):
        d = get_door("RO1", fresh_db)
        assert d is not None
        assert d["opening_name"] == "RO1"
        assert d["width"] == 36.0  # RO1_DOOR_WIDTH = 36/12 * 12 = 36"

    def test_update_door(self, fresh_db):
        updated = update_door("RO1", {"hinge_side": "west"}, fresh_db)
        assert updated["hinge_side"] == "west"

    def test_delete_door(self, fresh_db):
        result = delete_door("RO1", fresh_db)
        assert result is True
        assert get_door("RO1", fresh_db) is None


# ── Door / opening_type consistency ───────────────────────────────────
# An outer-wall opening must never be both a door (door record present) and
# a window/casement (its opening_type renders a void).  Door create/delete
# keeps opening_type in sync.

def _otype(name, db):
    return json.loads(get_element_by_name(name, db)["properties"]).get("opening_type")


class TestDoorOpeningTypeSync:

    def test_create_door_flips_window_to_door(self, fresh_db):
        assert _otype("O1", fresh_db) == "window"
        create_door("O1", 30, "east", "south", "single", fresh_db)
        assert _otype("O1", fresh_db) == "door"

    def test_delete_door_reverts_to_window(self, fresh_db):
        create_door("O1", 30, "east", "south", "single", fresh_db)
        delete_door("O1", fresh_db)
        assert _otype("O1", fresh_db) == "window"

    def test_interior_opening_untouched(self, fresh_db):
        # Rough openings carry no opening_type (always passages) and must not
        # gain one when a door is added/removed.
        before = json.loads(get_element_by_name("RO1", fresh_db)["properties"])
        assert "opening_type" not in before
        delete_door("RO1", fresh_db)
        create_door("RO1", 36, "east", "south", "single", fresh_db)
        after = json.loads(get_element_by_name("RO1", fresh_db)["properties"])
        assert "opening_type" not in after

    def test_create_door_on_missing_element_is_noop(self, fresh_db):
        # No element by this name — must not raise.
        rec = create_door("NO_SUCH_OPENING", 30, "east", "south", "single", fresh_db)
        assert rec["opening_name"] == "NO_SUCH_OPENING"

    def test_bottom_elev_follows_door_state(self, fresh_db):
        # Reopen to run the migration branch that seeds elevation constants.
        init_db(fresh_db)
        create_door("O1", 30, "east", "south", "single", fresh_db)
        assert get_constants_dict(fresh_db)["O1_BOTTOM_ELEV"] == 0.0
        delete_door("O1", fresh_db)
        assert get_constants_dict(fresh_db)["O1_BOTTOM_ELEV"] == LOWER_WALL_HEIGHT

    def test_seed_has_no_door_window_conflict(self, fresh_db):
        door_names = {d["opening_name"] for d in get_all_doors(fresh_db)}
        for name in door_names:
            el = get_element_by_name(name, fresh_db)
            if not el:
                continue
            props = json.loads(el["properties"])
            if "opening_type" in props:
                assert props["opening_type"] == "door", (
                    f"{name} has a door record but opening_type="
                    f"{props['opening_type']!r}")


# ── Door validation ───────────────────────────────────────────────────

class TestDoorValidation:

    def test_valid_door(self):
        assert validate_door("east", "south", "single") is None

    def test_invalid_hinge(self):
        assert validate_door("top", "south") is not None

    def test_invalid_type(self):
        assert validate_door("east", "south", "sliding") is not None


# ── API-27: POST /api/doors ───────────────────────────────────────────

class TestAPI27CreateDoor:
    """API-27: POST /api/doors SHALL create a door on an opening."""

    def test_create_door_201(self, app_client):
        resp = app_client.post("/api/doors", json={
            "opening": "TEST_NEW_RO",
            "hinge_side": "east", "swing_direction": "south",
            "width": 42, "door_type": "single",
        })
        assert resp.status_code == 201
        data = resp.get_json()
        assert data["opening_name"] == "TEST_NEW_RO"
        assert data["width"] == 42

    def test_create_duplicate_400(self, app_client):
        # RO1 already has a seeded door
        resp = app_client.post("/api/doors", json={
            "opening": "RO1",
            "hinge_side": "east", "swing_direction": "south",
            "width": 36,
        })
        assert resp.status_code == 400


# ── API-28: PUT /api/doors/<opening_name> ─────────────────────────────

class TestAPI28UpdateDoor:
    """API-28: PUT /api/doors/<opening_name> SHALL update door properties."""

    def test_update_hinge_side(self, app_client):
        resp = app_client.put("/api/doors/RO1", json={"hinge_side": "west"})
        assert resp.status_code == 200
        assert resp.get_json()["hinge_side"] == "west"

    def test_update_nonexistent_404(self, app_client):
        resp = app_client.put("/api/doors/NONEXISTENT", json={"width": 30})
        assert resp.status_code == 404


# ── API-29: DELETE /api/doors/<opening_name> ──────────────────────────

class TestAPI29DeleteDoor:
    """API-29: DELETE /api/doors/<opening_name> SHALL remove the door."""

    def test_delete_door_200(self, app_client):
        resp = app_client.delete("/api/doors/RO1")
        assert resp.status_code == 200
        assert resp.get_json()["ok"] is True
        # Verify gone
        doors_resp = app_client.get("/api/doors")
        names = [d["opening_name"] for d in doors_resp.get_json()]
        assert "RO1" not in names

    def test_delete_nonexistent_404(self, app_client):
        resp = app_client.delete("/api/doors/NONEXISTENT")
        assert resp.status_code == 404


# ── Door Undo (UNDO-4) ───────────────────────────────────────────────

class TestDoorUndo:
    """UNDO-4: Undo works for door mutations."""

    def test_undo_door_create(self, app_client):
        resp = app_client.post("/api/doors", json={
            "opening": "UNDO_DOOR",
            "hinge_side": "east", "swing_direction": "south",
            "width": 36,
        })
        assert resp.status_code == 201
        # Undo
        resp = app_client.post("/api/undo")
        assert resp.status_code == 200
        # Verify removed
        resp = app_client.get("/api/doors")
        names = [d["opening_name"] for d in resp.get_json()]
        assert "UNDO_DOOR" not in names

    def test_undo_door_update(self, app_client):
        # RO1 starts with hinge_side="east"
        resp = app_client.get("/api/doors")
        ro1 = next(d for d in resp.get_json() if d["opening_name"] == "RO1")
        original_hinge = ro1["hinge_side"]
        # Update
        app_client.put("/api/doors/RO1", json={"hinge_side": "west"})
        # Verify changed
        resp = app_client.get("/api/doors")
        ro1 = next(d for d in resp.get_json() if d["opening_name"] == "RO1")
        assert ro1["hinge_side"] == "west"
        # Undo
        resp = app_client.post("/api/undo")
        assert resp.status_code == 200
        # Verify restored
        resp = app_client.get("/api/doors")
        ro1 = next(d for d in resp.get_json() if d["opening_name"] == "RO1")
        assert ro1["hinge_side"] == original_hinge


# ── RT-5: element_changed event ──────────────────────────────────────

class TestRT5ElementChanged:
    """RT-5: Mutations broadcast element_changed event."""

    def test_element_changed_broadcast(self, app_client):
        # This is a structural test — we verify the endpoint works
        # and doesn't error.  Full SSE testing requires a WebSocket client.
        resp = app_client.post("/api/elements", json={
            "type": "furniture", "name": "SSE_TEST",
        })
        assert resp.status_code == 201
        # The broadcast happens but we can't observe it via test client.
        # Verify the element was created (side effect of the mutation).
        resp = app_client.get("/api/elements")
        names = [e["name"] for e in resp.get_json()]
        assert "SSE_TEST" in names
