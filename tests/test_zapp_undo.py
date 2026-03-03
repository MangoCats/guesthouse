"""Tests for the undo/redo system (Phase 2).

Covers: DB-11, API-30–31, UNDO-1–4.
"""
import json
import sqlite3

import pytest

from app.database import (
    get_db, init_db, get_constant_value, get_constants_dict,
    update_constant, update_constants_batch, reset_constants,
)
from app.undo import UndoManager

# Re-use fixtures from test_zapp_conftest.py
from tests.test_zapp_conftest import fresh_db, app_client  # noqa: F401


# ── DB-11: undo_history table ────────────────────────────────────────

class TestDB11UndoHistory:
    """DB-11: The database SHALL maintain an undo_history table."""

    def test_table_exists(self, fresh_db):
        with get_db(fresh_db) as conn:
            names = {r[0] for r in conn.execute(
                "SELECT name FROM sqlite_master WHERE type='table'"
            ).fetchall()}
        assert "undo_history" in names

    def test_table_empty_initially(self, fresh_db):
        with get_db(fresh_db) as conn:
            count = conn.execute("SELECT count(*) FROM undo_history").fetchone()[0]
        assert count == 0

    def test_entry_after_constant_update(self, fresh_db):
        mgr = UndoManager(fresh_db)
        mgr.record("constant_update", {"X": 1.0}, {"X": 2.0}, "test")
        with get_db(fresh_db) as conn:
            count = conn.execute("SELECT count(*) FROM undo_history").fetchone()[0]
        assert count == 1
        with get_db(fresh_db) as conn:
            row = conn.execute("SELECT * FROM undo_history").fetchone()
        assert row["action_type"] == "constant_update"
        assert json.loads(row["before_state"]) == {"X": 1.0}
        assert json.loads(row["after_state"]) == {"X": 2.0}


# ── UndoManager unit tests ──────────────────────────────────────────

class TestUndoManager:
    """Unit tests for UndoManager (UNDO-3, UNDO-4)."""

    def test_record_and_can_undo(self, fresh_db):
        mgr = UndoManager(fresh_db)
        assert not mgr.can_undo
        mgr.record("constant_update", {"X": 1.0}, {"X": 2.0})
        assert mgr.can_undo

    def test_undo_restores_value(self, fresh_db):
        old_val = get_constant_value("WALL_OUTER", fresh_db)
        update_constant("WALL_OUTER", old_val + 1.0, fresh_db)
        mgr = UndoManager(fresh_db)
        mgr.record(
            "constant_update",
            {"WALL_OUTER": old_val},
            {"WALL_OUTER": old_val + 1.0},
        )
        mgr.undo()
        assert get_constant_value("WALL_OUTER", fresh_db) == pytest.approx(old_val)

    def test_redo_reapplies_value(self, fresh_db):
        old_val = get_constant_value("WALL_OUTER", fresh_db)
        new_val = old_val + 1.0
        update_constant("WALL_OUTER", new_val, fresh_db)
        mgr = UndoManager(fresh_db)
        mgr.record("constant_update", {"WALL_OUTER": old_val}, {"WALL_OUTER": new_val})
        mgr.undo()
        assert get_constant_value("WALL_OUTER", fresh_db) == pytest.approx(old_val)
        mgr.redo()
        assert get_constant_value("WALL_OUTER", fresh_db) == pytest.approx(new_val)

    def test_undo_then_new_action_discards_redo(self, fresh_db):
        mgr = UndoManager(fresh_db)
        mgr.record("constant_update", {"X": 1.0}, {"X": 2.0})
        mgr.record("constant_update", {"X": 2.0}, {"X": 3.0})
        mgr.undo()  # position back to entry 0
        assert mgr.can_redo
        # New action discards the redo entry
        mgr.record("constant_update", {"X": 2.0}, {"X": 4.0})
        assert not mgr.can_redo

    def test_max_depth_50(self, fresh_db):
        mgr = UndoManager(fresh_db, max_depth=50)
        for i in range(60):
            mgr.record("constant_update", {"X": float(i)}, {"X": float(i + 1)})
        # Stack should be capped at 50
        count = 0
        while mgr.can_undo:
            mgr.undo()
            count += 1
        assert count == 50

    def test_undo_nothing_returns_none(self, fresh_db):
        mgr = UndoManager(fresh_db)
        assert mgr.undo() is None

    def test_redo_nothing_returns_none(self, fresh_db):
        mgr = UndoManager(fresh_db)
        assert mgr.redo() is None

    def test_batch_undo(self, fresh_db):
        vals = get_constants_dict(fresh_db)
        names = list(vals.keys())[:3]
        before = {n: vals[n] for n in names}
        after = {n: vals[n] + 100.0 for n in names}
        update_constants_batch(after, fresh_db)
        mgr = UndoManager(fresh_db)
        mgr.record("constant_batch", before, after)
        mgr.undo()
        for n in names:
            assert get_constant_value(n, fresh_db) == pytest.approx(before[n])

    def test_reset_undo(self, fresh_db):
        old_val = get_constant_value("WALL_OUTER", fresh_db)
        new_val = old_val + 5.0
        update_constant("WALL_OUTER", new_val, fresh_db)
        before = get_constants_dict(fresh_db)
        reset_constants(fresh_db)
        after = get_constants_dict(fresh_db)
        mgr = UndoManager(fresh_db)
        mgr.record("constant_reset", before, after)
        # After reset, WALL_OUTER should be back to original
        assert get_constant_value("WALL_OUTER", fresh_db) == pytest.approx(old_val)
        # Undo the reset — should restore the modified value
        mgr.undo()
        assert get_constant_value("WALL_OUTER", fresh_db) == pytest.approx(new_val)

    def test_persist_and_reload(self, fresh_db):
        mgr1 = UndoManager(fresh_db)
        mgr1.record("constant_update", {"X": 1.0}, {"X": 2.0}, "first")
        mgr1.record("constant_update", {"X": 2.0}, {"X": 3.0}, "second")
        # Create a new manager — should load from DB
        mgr2 = UndoManager(fresh_db)
        assert mgr2.can_undo
        entry = mgr2.undo()
        assert entry["description"] == "second"


# ── API tests (API-30, API-31) ───────────────────────────────────────

class TestAPI30Undo:
    """API-30: POST /api/undo SHALL revert the most recent edit."""

    def test_undo_constant_change(self, app_client):
        # Get original value
        resp = app_client.get("/api/constants")
        constants = resp.get_json()
        wt = next(c for c in constants if c["name"] == "WALL_OUTER")
        original = wt["value"]
        # Change it
        app_client.put(
            "/api/constants/WALL_OUTER",
            json={"value": original + 1.0},
        )
        # Verify changed
        resp = app_client.get("/api/constants")
        wt = next(c for c in resp.get_json() if c["name"] == "WALL_OUTER")
        assert wt["value"] == pytest.approx(original + 1.0)
        # Undo
        resp = app_client.post("/api/undo")
        assert resp.status_code == 200
        data = resp.get_json()
        assert data["ok"] is True
        assert data["action"] == "constant_update"
        # Verify restored
        resp = app_client.get("/api/constants")
        wt = next(c for c in resp.get_json() if c["name"] == "WALL_OUTER")
        assert wt["value"] == pytest.approx(original)

    def test_undo_returns_can_undo_can_redo(self, app_client):
        """Undo response includes can_undo/can_redo status."""
        resp = app_client.get("/api/constants")
        wt = next(c for c in resp.get_json() if c["name"] == "WALL_OUTER")
        original = wt["value"]
        # Two edits
        app_client.put("/api/constants/WALL_OUTER", json={"value": original + 1.0})
        app_client.put("/api/constants/WALL_OUTER", json={"value": original + 2.0})
        # Undo first — should still have more to undo, and can redo
        resp = app_client.post("/api/undo")
        data = resp.get_json()
        assert data["can_undo"] is True
        assert data["can_redo"] is True
        # Undo second — nothing left to undo, can still redo
        resp = app_client.post("/api/undo")
        data = resp.get_json()
        assert data["can_undo"] is False
        assert data["can_redo"] is True

    def test_undo_nothing_400(self, app_client):
        resp = app_client.post("/api/undo")
        assert resp.status_code == 400
        assert "nothing to undo" in resp.get_json()["error"]


class TestAPI31Redo:
    """API-31: POST /api/redo SHALL re-apply the most recently undone edit."""

    def test_redo_after_undo(self, app_client):
        resp = app_client.get("/api/constants")
        wt = next(c for c in resp.get_json() if c["name"] == "WALL_OUTER")
        original = wt["value"]
        new_val = original + 2.0
        app_client.put("/api/constants/WALL_OUTER", json={"value": new_val})
        app_client.post("/api/undo")
        # Redo
        resp = app_client.post("/api/redo")
        assert resp.status_code == 200
        data = resp.get_json()
        assert data["ok"] is True
        # Verify re-applied
        resp = app_client.get("/api/constants")
        wt = next(c for c in resp.get_json() if c["name"] == "WALL_OUTER")
        assert wt["value"] == pytest.approx(new_val)

    def test_redo_returns_can_undo_can_redo(self, app_client):
        """Redo response includes can_undo/can_redo status."""
        resp = app_client.get("/api/constants")
        wt = next(c for c in resp.get_json() if c["name"] == "WALL_OUTER")
        original = wt["value"]
        app_client.put("/api/constants/WALL_OUTER", json={"value": original + 1.0})
        app_client.post("/api/undo")
        # Redo — should be able to undo again, no more redo
        resp = app_client.post("/api/redo")
        data = resp.get_json()
        assert data["can_undo"] is True
        assert data["can_redo"] is False

    def test_redo_nothing_400(self, app_client):
        resp = app_client.post("/api/redo")
        assert resp.status_code == 400
        assert "nothing to redo" in resp.get_json()["error"]


# ── Integration tests (UNDO-1–4) ────────────────────────────────────

class TestUndoIntegration:
    """Integration tests for the full undo/redo cycle."""

    def test_50_level_depth(self, app_client):
        """UNDO-3: At least 50 levels of undo."""
        resp = app_client.get("/api/constants")
        wt = next(c for c in resp.get_json() if c["name"] == "WALL_OUTER")
        base = wt["value"]
        # Perform 55 edits
        for i in range(55):
            app_client.put(
                "/api/constants/WALL_OUTER",
                json={"value": base + i + 1},
            )
        # Undo 50 times — all should succeed
        for i in range(50):
            resp = app_client.post("/api/undo")
            assert resp.status_code == 200, f"Undo {i+1} failed"
        # 51st undo should still succeed (we have 55 entries, max_depth=50 trims oldest 5)
        # After 50 undos from 50-entry stack, nothing left
        resp = app_client.post("/api/undo")
        assert resp.status_code == 400

    def test_cross_type_undo(self, app_client):
        """UNDO-4: Undo works across different mutation types."""
        resp = app_client.get("/api/constants")
        constants = {c["name"]: c["value"] for c in resp.get_json()}
        wt_orig = constants["WALL_OUTER"]

        # 1. Single constant update
        app_client.put("/api/constants/WALL_OUTER", json={"value": wt_orig + 1.0})
        # 2. Batch update
        app_client.put("/api/constants/batch", json={
            "updates": {"WALL_OUTER": wt_orig + 2.0}
        })
        # 3. Reset
        app_client.post("/api/constants/reset")

        # Undo reset
        resp = app_client.post("/api/undo")
        assert resp.status_code == 200
        assert resp.get_json()["action"] == "constant_reset"
        resp = app_client.get("/api/constants")
        wt = next(c for c in resp.get_json() if c["name"] == "WALL_OUTER")
        assert wt["value"] == pytest.approx(wt_orig + 2.0)

        # Undo batch
        resp = app_client.post("/api/undo")
        assert resp.status_code == 200
        assert resp.get_json()["action"] == "constant_batch"
        resp = app_client.get("/api/constants")
        wt = next(c for c in resp.get_json() if c["name"] == "WALL_OUTER")
        assert wt["value"] == pytest.approx(wt_orig + 1.0)

        # Undo single
        resp = app_client.post("/api/undo")
        assert resp.status_code == 200
        assert resp.get_json()["action"] == "constant_update"
        resp = app_client.get("/api/constants")
        wt = next(c for c in resp.get_json() if c["name"] == "WALL_OUTER")
        assert wt["value"] == pytest.approx(wt_orig)

    def test_undo_reset_restores_all_values(self, app_client):
        """Undo of a reset restores all pre-reset constant values."""
        # Modify two constants
        resp = app_client.get("/api/constants")
        constants = {c["name"]: c["value"] for c in resp.get_json()}
        names = sorted(constants.keys())[:2]
        for n in names:
            app_client.put(f"/api/constants/{n}", json={"value": constants[n] + 99.0})
        # Reset
        app_client.post("/api/constants/reset")
        # Verify reset happened
        resp = app_client.get("/api/constants")
        for c in resp.get_json():
            if c["name"] in names:
                assert c["value"] == pytest.approx(constants[c["name"]])
        # Undo the reset (skip the 2 individual updates to get to reset)
        # Actually reset is the most recent action
        resp = app_client.post("/api/undo")
        assert resp.status_code == 200
        # Verify the modified values are restored
        resp = app_client.get("/api/constants")
        restored = {c["name"]: c["value"] for c in resp.get_json()}
        for n in names:
            assert restored[n] == pytest.approx(constants[n] + 99.0)
