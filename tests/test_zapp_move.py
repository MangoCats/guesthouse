"""Tests for Phase 4 — Move Tool (API-23, element_move undo, IW axis mapping).

Covers: API-23, TL-5 axis mapping, element_move undo action.
"""
import json

import pytest

from app.database import (
    get_db, get_constant_value, get_element, get_element_by_name,
    create_element, update_element,
)
from app.elements import (
    IW_CONSTANT_MAP, IW_MOVE_AXIS, compute_constant_delta,
)

from tests.test_zapp_conftest import fresh_db, app_client  # noqa: F401


# ── IW Axis Mapping Unit Tests ──────────────────────────────────────

class TestMoveAxisMapping:
    """Verify IW_MOVE_AXIS and compute_constant_delta."""

    def test_iw9_constant_fixed(self):
        """IW9 should map to IW3_OFFSET_IW9 (not the old IW9_OFFSET_O10)."""
        assert IW_CONSTANT_MAP["IW9"] == "IW3_OFFSET_IW9"

    def test_all_moveable_walls_have_axis(self):
        """Every IW with a non-None constant must have a non-None axis entry."""
        for name, const in IW_CONSTANT_MAP.items():
            if const is not None:
                axis_info = IW_MOVE_AXIS.get(name)
                assert axis_info is not None, f"{name} has constant but no axis"
                axis, sign = axis_info
                assert axis in ("x", "y"), f"{name} axis must be x or y"
                assert sign in (+1, -1), f"{name} sign must be +1 or -1"

    def test_nonmovable_walls_have_none_axis(self):
        """IW2O and IW8 should have None axis (not directly movable)."""
        assert IW_MOVE_AXIS["IW2O"] is None
        assert IW_MOVE_AXIS["IW8"] is None

    def test_compute_constant_delta_iw1_dy(self):
        """IW1 moves on y-axis with sign -1: dy=-0.5 → delta=+0.5."""
        result = compute_constant_delta("IW1", 0, -0.5)
        assert result is not None
        const_name, delta = result
        assert const_name == "IW1_OFFSET_FROM_W9"
        assert abs(delta - 0.5) < 1e-9

    def test_compute_constant_delta_iw2_dx(self):
        """IW2 moves on x-axis with sign +1: dx=1.0 → delta=+1.0."""
        result = compute_constant_delta("IW2", 1.0, 0)
        assert result is not None
        const_name, delta = result
        assert const_name == "IW2_DIST_W2W5"
        assert abs(delta - 1.0) < 1e-9

    def test_compute_constant_delta_nonmovable(self):
        """IW2O and IW8 should return None (not movable)."""
        assert compute_constant_delta("IW2O", 1.0, 1.0) is None
        assert compute_constant_delta("IW8", 1.0, 1.0) is None

    def test_compute_constant_delta_zero_move(self):
        """Zero movement on the relevant axis should return None."""
        # IW1 moves on y-axis; dx-only should give zero on y → None
        assert compute_constant_delta("IW1", 1.0, 0) is None

    def test_compute_constant_delta_unknown_wall(self):
        """Unknown wall name should return None."""
        assert compute_constant_delta("FAKE_WALL", 1.0, 1.0) is None


# ── API-23: POST /api/elements/<id>/move ─────────────────────────────

class TestAPI23Move:
    """API-23: POST /api/elements/<id>/move SHALL reposition an element."""

    def test_move_iw_wall_updates_constant(self, app_client):
        """Moving IW1 south by 0.5 ft should increase IW1_OFFSET_FROM_W9."""
        # Get IW1 element id
        resp = app_client.get("/api/elements")
        elements = resp.get_json()
        iw1 = next(e for e in elements if e["name"] == "IW1")

        resp = app_client.post(f"/api/elements/{iw1['id']}/move", json={
            "dx": 0, "dy": -0.5,
        })
        assert resp.status_code == 200
        data = resp.get_json()
        assert data["ok"] is True
        assert data["constant"] == "IW1_OFFSET_FROM_W9"
        # IW1: y axis, sign -1 → delta = -0.5 * -1 = +0.5
        assert abs(data["new_value"] - data["old_value"] - 0.5) < 1e-9

    def test_move_iw_wall_correct_axis_dx(self, app_client):
        """Moving IW2 east by 1.0 ft should increase IW2_DIST_W2W5."""
        resp = app_client.get("/api/elements")
        iw2 = next(e for e in resp.get_json() if e["name"] == "IW2")

        resp = app_client.post(f"/api/elements/{iw2['id']}/move", json={
            "dx": 1.0, "dy": 0,
        })
        assert resp.status_code == 200
        data = resp.get_json()
        assert data["constant"] == "IW2_DIST_W2W5"
        # IW2: x axis, sign +1 → delta = 1.0 * +1 = +1.0
        assert abs(data["new_value"] - data["old_value"] - 1.0) < 1e-9

    def test_move_nonmovable_wall_400(self, app_client):
        """Moving IW2O or IW8 should return 400 (not movable)."""
        resp = app_client.get("/api/elements")
        elements = resp.get_json()
        iw2o = next(e for e in elements if e["name"] == "IW2O")
        iw8 = next(e for e in elements if e["name"] == "IW8")

        resp = app_client.post(f"/api/elements/{iw2o['id']}/move", json={
            "dx": 1.0, "dy": 0,
        })
        assert resp.status_code == 400

        resp = app_client.post(f"/api/elements/{iw8['id']}/move", json={
            "dx": 0, "dy": 1.0,
        })
        assert resp.status_code == 400

    def test_move_custom_element_updates_offset(self, app_client):
        """Moving a custom element should update offset_x/offset_y."""
        # Create a furniture override element
        resp = app_client.post("/api/elements", json={
            "type": "furniture", "name": "MOVE_BED",
            "properties": {"offset_x": 0, "offset_y": 0, "source": "override"},
        })
        assert resp.status_code == 201
        eid = resp.get_json()["id"]

        resp = app_client.post(f"/api/elements/{eid}/move", json={
            "dx": 2.0, "dy": -1.5,
        })
        assert resp.status_code == 200
        data = resp.get_json()
        assert data["ok"] is True
        assert abs(data["offset_x"] - 2.0) < 1e-9
        assert abs(data["offset_y"] - (-1.5)) < 1e-9

    def test_move_nonexistent_404(self, app_client):
        """Moving a non-existent element should return 404."""
        resp = app_client.post("/api/elements/99999/move", json={
            "dx": 1.0, "dy": 0,
        })
        assert resp.status_code == 404

    def test_move_missing_payload_400(self, app_client):
        """Moving with dx=0 and dy=0 should return 400."""
        resp = app_client.get("/api/elements")
        iw1 = next(e for e in resp.get_json() if e["name"] == "IW1")

        resp = app_client.post(f"/api/elements/{iw1['id']}/move", json={
            "dx": 0, "dy": 0,
        })
        assert resp.status_code == 400

    def test_move_returns_can_undo(self, app_client):
        """Move response should include can_undo=True."""
        resp = app_client.get("/api/elements")
        iw1 = next(e for e in resp.get_json() if e["name"] == "IW1")

        resp = app_client.post(f"/api/elements/{iw1['id']}/move", json={
            "dx": 0, "dy": -1.0,
        })
        assert resp.status_code == 200
        data = resp.get_json()
        assert data["can_undo"] is True

    def test_move_cumulative_offset(self, app_client):
        """Multiple moves should accumulate offset values."""
        resp = app_client.post("/api/elements", json={
            "type": "furniture", "name": "CUM_TEST",
            "properties": {"offset_x": 0, "offset_y": 0},
        })
        eid = resp.get_json()["id"]

        app_client.post(f"/api/elements/{eid}/move", json={"dx": 1.0, "dy": 0})
        resp = app_client.post(f"/api/elements/{eid}/move", json={"dx": 0.5, "dy": -0.5})
        data = resp.get_json()
        assert abs(data["offset_x"] - 1.5) < 1e-9
        assert abs(data["offset_y"] - (-0.5)) < 1e-9


# ── Move Undo Tests ──────────────────────────────────────────────────

class TestMoveUndo:
    """Verify undo/redo of element_move actions."""

    def test_undo_iw_move_restores_constant(self, app_client):
        """Undo after moving IW1 should restore the original constant value."""
        # Get original constant value
        resp = app_client.get("/api/constants")
        orig = {c["name"]: c["value"] for c in resp.get_json()}
        orig_iw1 = orig["IW1_OFFSET_FROM_W9"]

        # Move IW1
        resp = app_client.get("/api/elements")
        iw1 = next(e for e in resp.get_json() if e["name"] == "IW1")
        app_client.post(f"/api/elements/{iw1['id']}/move", json={
            "dx": 0, "dy": -1.0,
        })

        # Verify constant changed
        resp = app_client.get("/api/constants")
        after = {c["name"]: c["value"] for c in resp.get_json()}
        assert after["IW1_OFFSET_FROM_W9"] != orig_iw1

        # Undo
        resp = app_client.post("/api/undo")
        assert resp.status_code == 200

        # Verify constant restored
        resp = app_client.get("/api/constants")
        restored = {c["name"]: c["value"] for c in resp.get_json()}
        assert abs(restored["IW1_OFFSET_FROM_W9"] - orig_iw1) < 1e-9

    def test_undo_custom_element_move(self, app_client):
        """Undo after moving a custom element should restore offset."""
        resp = app_client.post("/api/elements", json={
            "type": "furniture", "name": "UNDO_BED",
            "properties": {"offset_x": 0, "offset_y": 0},
        })
        eid = resp.get_json()["id"]

        # Move it
        resp = app_client.post(f"/api/elements/{eid}/move", json={
            "dx": 3.0, "dy": -2.0,
        })
        assert resp.status_code == 200
        assert abs(resp.get_json()["offset_x"] - 3.0) < 1e-9

        # Undo
        resp = app_client.post("/api/undo")
        assert resp.status_code == 200

        # Verify restored via elements list
        resp = app_client.get("/api/elements")
        el = next(e for e in resp.get_json() if e["name"] == "UNDO_BED")
        props = json.loads(el["properties"])
        assert abs(props["offset_x"]) < 1e-9
        assert abs(props["offset_y"]) < 1e-9


# ── Furniture Override Tests ──────────────────────────────────────────

class TestFurnitureOverride:
    """Verify furniture override creation and movement."""

    def test_create_furniture_override(self, app_client):
        """Creating a furniture override element should work."""
        resp = app_client.post("/api/elements", json={
            "type": "furniture", "name": "override_bed",
            "variant": "standard",
            "properties": {"offset_x": 0, "offset_y": 0, "source": "override"},
        })
        assert resp.status_code == 201
        data = resp.get_json()
        assert data["name"] == "override_bed"
        props = json.loads(data["properties"])
        assert props["source"] == "override"

    def test_move_furniture_override_updates_offset(self, app_client):
        """Moving a furniture override should update its offset."""
        resp = app_client.post("/api/elements", json={
            "type": "furniture", "name": "override_desk",
            "properties": {"offset_x": 0, "offset_y": 0, "source": "override"},
        })
        eid = resp.get_json()["id"]

        resp = app_client.post(f"/api/elements/{eid}/move", json={
            "dx": 1.5, "dy": -0.75,
        })
        assert resp.status_code == 200
        data = resp.get_json()
        assert abs(data["offset_x"] - 1.5) < 1e-9
        assert abs(data["offset_y"] - (-0.75)) < 1e-9
