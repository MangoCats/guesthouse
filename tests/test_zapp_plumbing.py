"""Tests for Phase 10d Plumbing Layout integration (PLUMB-1 through PLUMB-8)."""
import os
import sys
import pytest

_PROJECT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _PROJECT not in sys.path:
    sys.path.insert(0, _PROJECT)

from tests.test_zapp_conftest import fresh_db, app_client
from app.database import get_views
from tests.test_zapp_conftest import geometry
from app.database import get_constants_dict
from app.plumbing import (
    get_plumbing_elements, get_plumbing_element,
    create_plumbing_element, update_plumbing_element,
    delete_plumbing_element, FIXTURE_DEFS, PLUMBING_TYPES,
    compute_reference_plumbing, seed_reference_plumbing,
)


# =========================================================================
# TestPlumbingSchema (PLUMB-4)
# =========================================================================

class TestPlumbingSchema:
    """Plumbing elements table exists with correct structure."""

    def test_table_exists(self, fresh_db):
        from app.database import get_db
        with get_db(fresh_db) as conn:
            rows = conn.execute(
                "SELECT name FROM sqlite_master WHERE type='table' AND name='plumbing_elements'"
            ).fetchall()
        assert len(rows) == 1

    def test_columns(self, fresh_db):
        from app.database import get_db
        with get_db(fresh_db) as conn:
            info = conn.execute("PRAGMA table_info(plumbing_elements)").fetchall()
        col_names = {r["name"] for r in info}
        assert {"id", "type", "name", "path", "properties", "fixture"}.issubset(col_names)

    def test_valid_types(self):
        assert PLUMBING_TYPES == {"supply_pipe", "drain_pipe", "fitting", "fixture_connection"}


# =========================================================================
# TestPlumbingSeed (PLUMB-4)
# =========================================================================

class TestPlumbingSeed:
    """Default fixture connection records seeded correctly."""

    def test_seed_count(self, fresh_db):
        elems = get_plumbing_elements(fresh_db)
        assert len(elems) == 11

    def test_all_fixture_connections(self, fresh_db):
        elems = get_plumbing_elements(fresh_db)
        assert all(e["type"] == "fixture_connection" for e in elems)

    def test_fixture_names(self, fresh_db):
        elems = get_plumbing_elements(fresh_db)
        names = {e["name"] for e in elems}
        expected = {name for name, *_ in FIXTURE_DEFS}
        assert names == expected

    def test_cold_hot_drain_flags(self, fresh_db):
        elems = get_plumbing_elements(fresh_db)
        washer = next(e for e in elems if e["name"] == "Washer")
        assert washer["properties"]["cold"] is True
        assert washer["properties"]["hot"] is True
        assert washer["properties"]["drain"] is True

        fridge = next(e for e in elems if e["name"] == "Fridge")
        assert fridge["properties"]["cold"] is True
        assert fridge["properties"]["hot"] is False
        assert fridge["properties"]["drain"] is False


# =========================================================================
# TestPlumbingCRUD (PLUMB-4)
# =========================================================================

class TestPlumbingCRUD:
    """Direct CRUD operations on plumbing elements."""

    def test_create_supply_pipe(self, fresh_db):
        elem = create_plumbing_element(
            "supply_pipe", "test_pipe",
            path=[[0, 0], [5, 0], [5, 3]],
            properties={"hot_cold": "cold", "pipe_size": "0.75"},
            db_path=fresh_db,
        )
        assert elem["id"] > 0
        assert elem["type"] == "supply_pipe"
        assert elem["path"] == [[0, 0], [5, 0], [5, 3]]

    def test_create_drain_pipe(self, fresh_db):
        elem = create_plumbing_element(
            "drain_pipe", "test_drain",
            path=[[0, 0], [0, -5]],
            properties={"slope": "0.25 in/ft"},
            db_path=fresh_db,
        )
        assert elem["type"] == "drain_pipe"

    def test_create_fitting(self, fresh_db):
        elem = create_plumbing_element(
            "fitting", "test_tee",
            path=[[3, 2]],
            properties={"fitting_type": "tee", "rotation": 90},
            db_path=fresh_db,
        )
        assert elem["type"] == "fitting"

    def test_invalid_type_raises(self, fresh_db):
        with pytest.raises(ValueError, match="Invalid plumbing type"):
            create_plumbing_element("invalid", "bad", db_path=fresh_db)

    def test_read_by_id(self, fresh_db):
        elem = create_plumbing_element(
            "supply_pipe", "read_test", path=[[1, 2]], db_path=fresh_db)
        fetched = get_plumbing_element(elem["id"], fresh_db)
        assert fetched is not None
        assert fetched["name"] == "read_test"

    def test_update(self, fresh_db):
        elem = create_plumbing_element(
            "supply_pipe", "upd_test",
            properties={"hot_cold": "cold"}, db_path=fresh_db)
        updated = update_plumbing_element(
            elem["id"], {"properties": {"hot_cold": "hot"}}, fresh_db)
        assert updated["properties"]["hot_cold"] == "hot"

    def test_delete(self, fresh_db):
        elem = create_plumbing_element(
            "supply_pipe", "del_test", db_path=fresh_db)
        result = delete_plumbing_element(elem["id"], fresh_db)
        assert result == elem["id"]
        assert get_plumbing_element(elem["id"], fresh_db) is None

    def test_delete_nonexistent(self, fresh_db):
        result = delete_plumbing_element(99999, fresh_db)
        assert result is None


# =========================================================================
# TestPlumbingAPI (PLUMB-5)
# =========================================================================

class TestPlumbingAPI:
    """REST API endpoints for plumbing elements."""

    def test_get_list(self, app_client):
        resp = app_client.get("/api/plumbing")
        assert resp.status_code == 200
        data = resp.get_json()
        assert len(data) == 17  # 11 fixtures + 6 reference pipes

    def test_create(self, app_client):
        resp = app_client.post("/api/plumbing", json={
            "type": "supply_pipe",
            "name": "api_pipe",
            "path": [[0, 0], [3, 0]],
            "properties": {"hot_cold": "cold"},
        })
        assert resp.status_code == 201
        data = resp.get_json()
        assert data["type"] == "supply_pipe"
        assert data["name"] == "api_pipe"

    def test_update(self, app_client):
        resp = app_client.post("/api/plumbing", json={
            "type": "supply_pipe", "name": "api_upd",
        })
        eid = resp.get_json()["id"]
        resp2 = app_client.put(f"/api/plumbing/{eid}", json={
            "properties": {"hot_cold": "hot"},
        })
        assert resp2.status_code == 200
        assert resp2.get_json()["properties"]["hot_cold"] == "hot"

    def test_update_not_found(self, app_client):
        resp = app_client.put("/api/plumbing/99999", json={"name": "x"})
        assert resp.status_code == 404

    def test_delete(self, app_client):
        resp = app_client.post("/api/plumbing", json={
            "type": "drain_pipe", "name": "api_del",
        })
        eid = resp.get_json()["id"]
        resp2 = app_client.delete(f"/api/plumbing/{eid}")
        assert resp2.status_code == 200
        assert resp2.get_json()["deleted"] == eid

    def test_delete_not_found(self, app_client):
        resp = app_client.delete("/api/plumbing/99999")
        assert resp.status_code == 404


# =========================================================================
# TestPlumbingUndo (PLUMB-5)
# =========================================================================

class TestPlumbingUndo:
    """Undo/redo for plumbing mutations."""

    def test_create_undo(self, app_client):
        resp = app_client.post("/api/plumbing", json={
            "type": "supply_pipe", "name": "undo_pipe",
        })
        eid = resp.get_json()["id"]
        # Undo should remove it
        app_client.post("/api/undo")
        resp2 = app_client.get("/api/plumbing")
        names = {e["name"] for e in resp2.get_json()}
        assert "undo_pipe" not in names

    def test_delete_undo(self, app_client):
        resp = app_client.post("/api/plumbing", json={
            "type": "drain_pipe", "name": "undo_del",
        })
        eid = resp.get_json()["id"]
        app_client.delete(f"/api/plumbing/{eid}")
        # Undo should restore it
        app_client.post("/api/undo")
        resp2 = app_client.get("/api/plumbing")
        names = {e["name"] for e in resp2.get_json()}
        assert "undo_del" in names

    def test_update_undo(self, app_client):
        resp = app_client.post("/api/plumbing", json={
            "type": "fitting", "name": "undo_upd",
            "properties": {"fitting_type": "tee"},
        })
        eid = resp.get_json()["id"]
        app_client.put(f"/api/plumbing/{eid}", json={
            "properties": {"fitting_type": "elbow"},
        })
        # Undo should restore original
        app_client.post("/api/undo")
        resp2 = app_client.get("/api/plumbing")
        elem = next(e for e in resp2.get_json() if e["name"] == "undo_upd")
        assert elem["properties"]["fitting_type"] == "tee"


# =========================================================================
# TestPlumbingView (PLUMB-1)
# =========================================================================

class TestPlumbingView:
    """Plumbing is a layout variant with a view for SVG regeneration."""

    def test_plumbing_variant_registered(self, fresh_db):
        from app.database import get_variants
        variants = get_variants(fresh_db)
        names = {v["name"] for v in variants}
        assert "plumbing" in names

    def test_plumbing_view_for_regeneration(self, fresh_db):
        views = get_views(fresh_db)
        names = {v["name"] for v in views}
        assert "plumbing" in names

    def test_plumbing_edit_view_removed(self, fresh_db):
        views = get_views(fresh_db)
        names = {v["name"] for v in views}
        assert "plumbing_edit" not in names


# =========================================================================
# TestPlumbingReferenceSeed
# =========================================================================

class TestPlumbingReferenceSeed:
    """Reference plumbing configuration seeded into database."""

    def test_seed_creates_pipes(self, fresh_db, geometry):
        constants = get_constants_dict(fresh_db)
        wall_t = constants.get("WALL_OUTER", 10.0 / 12.0)
        seed_reference_plumbing(geometry, wall_t, fresh_db)
        elems = get_plumbing_elements(fresh_db)
        pipes = [e for e in elems if e["type"] in ("supply_pipe", "drain_pipe")]
        assert len(pipes) == 6

    def test_seed_pipe_types(self, fresh_db, geometry):
        constants = get_constants_dict(fresh_db)
        wall_t = constants.get("WALL_OUTER", 10.0 / 12.0)
        seed_reference_plumbing(geometry, wall_t, fresh_db)
        elems = get_plumbing_elements(fresh_db)
        supply = [e for e in elems if e["type"] == "supply_pipe"]
        drain = [e for e in elems if e["type"] == "drain_pipe"]
        assert len(supply) == 5
        assert len(drain) == 1

    def test_seed_pipe_paths_nonempty(self, fresh_db, geometry):
        constants = get_constants_dict(fresh_db)
        wall_t = constants.get("WALL_OUTER", 10.0 / 12.0)
        seed_reference_plumbing(geometry, wall_t, fresh_db)
        elems = get_plumbing_elements(fresh_db)
        pipes = [e for e in elems if e["type"] in ("supply_pipe", "drain_pipe")]
        for p in pipes:
            assert len(p["path"]) >= 2, f"{p['name']} has < 2 waypoints"

    def test_seed_fixture_positions(self, fresh_db, geometry):
        constants = get_constants_dict(fresh_db)
        wall_t = constants.get("WALL_OUTER", 10.0 / 12.0)
        seed_reference_plumbing(geometry, wall_t, fresh_db)
        elems = get_plumbing_elements(fresh_db)
        fixtures = [e for e in elems if e["type"] == "fixture_connection"]
        assert len(fixtures) == 11
        for f in fixtures:
            assert len(f["path"]) == 1, f"{f['name']} has no position"
            assert len(f["path"][0]) == 2, f"{f['name']} position not [E, N]"

    def test_seed_idempotent(self, fresh_db, geometry):
        constants = get_constants_dict(fresh_db)
        wall_t = constants.get("WALL_OUTER", 10.0 / 12.0)
        seed_reference_plumbing(geometry, wall_t, fresh_db)
        seed_reference_plumbing(geometry, wall_t, fresh_db)  # 2nd call
        elems = get_plumbing_elements(fresh_db)
        pipes = [e for e in elems if e["type"] in ("supply_pipe", "drain_pipe")]
        assert len(pipes) == 6  # no duplicates

    def test_compute_reference_returns_valid(self, geometry, fresh_db):
        constants = get_constants_dict(fresh_db)
        wall_t = constants.get("WALL_OUTER", 10.0 / 12.0)
        ref = compute_reference_plumbing(geometry, wall_t)
        assert len(ref["pipes"]) == 6
        assert len(ref["fixture_positions"]) == 11
        for name in ("Washer", "Toilet1", "Toilet2", "Util Sink",
                      "Bath Sink", "Fridge", "Shower", "Kitchen Sink",
                      "Dishwasher", "Ice Maker", "Water Heater"):
            assert name in ref["fixture_positions"]
