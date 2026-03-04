"""Tests for Phase 8: Labels, Dimensions, and Annotations.

Covers: TL-11–14, LABEL-1–4, DIS-7 (9 reqs).
"""
import json

import pytest

from app.database import (
    get_db, init_db, get_all_elements, get_element,
    create_element, update_element, delete_element,
)
from app.labels import (
    seed_room_labels, next_dimension_name, next_label_name,
    ROOM_LABEL_NAMES,
)
from app.engine import compute_geometry
from app.database import get_constants_dict

from tests.test_zapp_conftest import fresh_db, app_client  # noqa: F401


# ── Room label seeding ───────────────────────────────────────────────

class TestRoomLabelSeeding:
    """Room label elements are seeded on fresh DB."""

    def test_11_room_labels_exist(self, fresh_db):
        elements = get_all_elements(fresh_db)
        labels = [e for e in elements if e["type"] == "label"]
        assert len(labels) == 11

    def test_room_label_names(self, fresh_db):
        elements = get_all_elements(fresh_db)
        names = {e["name"] for e in elements if e["type"] == "label"}
        assert names == set(ROOM_LABEL_NAMES)

    def test_room_label_properties(self, fresh_db):
        elements = get_all_elements(fresh_db)
        for e in elements:
            if e["type"] != "label":
                continue
            props = json.loads(e["properties"]) if isinstance(e["properties"], str) else e["properties"]
            assert props["source"] == "room"
            assert "offset_e" in props
            assert "offset_n" in props
            assert "font_size" in props
            assert props["text"] == e["name"]

    def test_seed_idempotent(self, fresh_db):
        """Calling seed_room_labels again doesn't duplicate."""
        with get_db(fresh_db) as conn:
            seed_room_labels(conn)
        elements = get_all_elements(fresh_db)
        labels = [e for e in elements if e["type"] == "label"]
        assert len(labels) == 11

    def test_room_label_variant_null(self, fresh_db):
        """Room labels have variant=NULL (visible to all)."""
        elements = get_all_elements(fresh_db)
        for e in elements:
            if e["type"] == "label":
                assert e["variant"] is None


# ── Name generators ──────────────────────────────────────────────────

class TestNameGenerators:

    def test_next_dimension_name_empty(self):
        assert next_dimension_name([]) == "UD1"

    def test_next_dimension_name_increments(self):
        elems = [{"name": "UD1"}, {"name": "UD3"}, {"name": "BEDROOM"}]
        assert next_dimension_name(elems) == "UD4"

    def test_next_label_name_empty(self):
        assert next_label_name([]) == "UL1"

    def test_next_label_name_increments(self):
        elems = [{"name": "UL1"}, {"name": "UL2"}, {"name": "KITCHEN"}]
        assert next_label_name(elems) == "UL3"


# ── Dimension CRUD ───────────────────────────────────────────────────

class TestDimensionCRUD:
    """Dimension element CRUD operations."""

    def test_create_dimension(self, fresh_db):
        rec = create_element("dimension", "UD1", {
            "source": "user", "start": [0, 0], "end": [5, 0],
            "offset": 0, "label_rotation": "parallel",
        }, None, fresh_db)
        assert rec["type"] == "dimension"
        assert rec["name"] == "UD1"

    def test_update_dimension_offset(self, fresh_db):
        rec = create_element("dimension", "UD1", {
            "source": "user", "start": [0, 0], "end": [5, 0],
            "offset": 0, "label_rotation": "parallel",
        }, None, fresh_db)
        updated = update_element(rec["id"], {"properties": {
            "source": "user", "start": [0, 0], "end": [5, 0],
            "offset": 1.5, "label_rotation": "parallel",
        }}, fresh_db)
        props = json.loads(updated["properties"]) if isinstance(updated["properties"], str) else updated["properties"]
        assert props["offset"] == 1.5

    def test_delete_dimension(self, fresh_db):
        rec = create_element("dimension", "UD1", {
            "source": "user", "start": [0, 0], "end": [5, 0],
        }, None, fresh_db)
        deleted = delete_element(rec["id"], fresh_db)
        assert rec["id"] in deleted
        assert get_element(rec["id"], fresh_db) is None


# ── Label CRUD ───────────────────────────────────────────────────────

class TestLabelCRUD:
    """User label element CRUD operations."""

    def test_create_user_label(self, fresh_db):
        rec = create_element("label", "UL1", {
            "source": "user", "text": "MY LABEL",
            "position": [5.0, -3.0], "rotation": 0, "font_size": 0.25,
        }, None, fresh_db)
        assert rec["type"] == "label"
        assert rec["name"] == "UL1"

    def test_update_label_text(self, fresh_db):
        rec = create_element("label", "UL1", {
            "source": "user", "text": "ORIGINAL",
            "position": [0, 0], "font_size": 0.25,
        }, None, fresh_db)
        updated = update_element(rec["id"], {"properties": {
            "source": "user", "text": "EDITED",
            "position": [0, 0], "font_size": 0.25,
        }}, fresh_db)
        props = json.loads(updated["properties"]) if isinstance(updated["properties"], str) else updated["properties"]
        assert props["text"] == "EDITED"

    def test_update_label_font_size(self, fresh_db):
        rec = create_element("label", "UL1", {
            "source": "user", "text": "TEST",
            "position": [0, 0], "font_size": 0.25,
        }, None, fresh_db)
        updated = update_element(rec["id"], {"properties": {
            "source": "user", "text": "TEST",
            "position": [0, 0], "font_size": 0.5,
        }}, fresh_db)
        props = json.loads(updated["properties"]) if isinstance(updated["properties"], str) else updated["properties"]
        assert props["font_size"] == 0.5

    def test_delete_user_label(self, fresh_db):
        rec = create_element("label", "UL1", {
            "source": "user", "text": "DEL", "position": [0, 0],
        }, None, fresh_db)
        deleted = delete_element(rec["id"], fresh_db)
        assert rec["id"] in deleted


# ── API tests ────────────────────────────────────────────────────────

class TestDimensionAPI:
    """API CRUD for dimensions."""

    def test_create_dimension_201(self, app_client):
        resp = app_client.post("/api/elements", json={
            "type": "dimension", "name": "UD1",
            "properties": {
                "source": "user", "start": [0, 0], "end": [5, 0],
                "offset": 0, "label_rotation": "parallel",
            },
        })
        assert resp.status_code == 201
        assert resp.get_json()["name"] == "UD1"

    def test_update_dimension_200(self, app_client):
        resp = app_client.post("/api/elements", json={
            "type": "dimension", "name": "UD1",
            "properties": {"start": [0, 0], "end": [5, 0], "offset": 0},
        })
        eid = resp.get_json()["id"]
        resp = app_client.put(f"/api/elements/{eid}", json={
            "properties": {"start": [0, 0], "end": [5, 0], "offset": 2.0},
        })
        assert resp.status_code == 200

    def test_delete_dimension_200(self, app_client):
        resp = app_client.post("/api/elements", json={
            "type": "dimension", "name": "UD1",
            "properties": {"start": [0, 0], "end": [5, 0]},
        })
        eid = resp.get_json()["id"]
        resp = app_client.delete(f"/api/elements/{eid}")
        assert resp.status_code == 200

    def test_undo_dimension_create(self, app_client):
        resp = app_client.post("/api/elements", json={
            "type": "dimension", "name": "UD_UNDO",
            "properties": {"start": [0, 0], "end": [3, 0]},
        })
        assert resp.status_code == 201
        resp = app_client.post("/api/undo")
        assert resp.status_code == 200
        resp = app_client.get("/api/elements")
        names = [e["name"] for e in resp.get_json()]
        assert "UD_UNDO" not in names


class TestLabelAPI:
    """API CRUD for labels."""

    def test_create_user_label_201(self, app_client):
        resp = app_client.post("/api/elements", json={
            "type": "label", "name": "UL1",
            "properties": {
                "source": "user", "text": "HELLO",
                "position": [1, 2], "font_size": 0.3,
            },
        })
        assert resp.status_code == 201
        assert resp.get_json()["name"] == "UL1"

    def test_update_label_text_200(self, app_client):
        resp = app_client.post("/api/elements", json={
            "type": "label", "name": "UL1",
            "properties": {"source": "user", "text": "OLD", "position": [0, 0]},
        })
        eid = resp.get_json()["id"]
        resp = app_client.put(f"/api/elements/{eid}", json={
            "properties": {"source": "user", "text": "NEW", "position": [0, 0]},
        })
        assert resp.status_code == 200
        props = json.loads(resp.get_json()["properties"])
        assert props["text"] == "NEW"

    def test_undo_label_create(self, app_client):
        resp = app_client.post("/api/elements", json={
            "type": "label", "name": "UL_UNDO",
            "properties": {"source": "user", "text": "T", "position": [0, 0]},
        })
        assert resp.status_code == 201
        resp = app_client.post("/api/undo")
        assert resp.status_code == 200
        resp = app_client.get("/api/elements")
        names = [e["name"] for e in resp.get_json()]
        assert "UL_UNDO" not in names


# ── Geometry response ────────────────────────────────────────────────

class TestGeometryOutput:
    """Geometry response includes user_dimensions and label_elements."""

    def test_geometry_has_user_dimensions(self, fresh_db):
        constants = get_constants_dict(fresh_db)
        geom = compute_geometry(constants)
        assert "user_dimensions" in geom
        assert isinstance(geom["user_dimensions"], list)

    def test_geometry_has_label_elements(self, fresh_db):
        constants = get_constants_dict(fresh_db)
        geom = compute_geometry(constants)
        assert "label_elements" in geom
        assert len(geom["label_elements"]) == 11  # room labels

    def test_label_elements_have_centroid(self, fresh_db):
        constants = get_constants_dict(fresh_db)
        geom = compute_geometry(constants)
        for le in geom["label_elements"]:
            if le["properties"]["source"] == "room":
                assert "centroid" in le
                assert "pos" in le

    def test_geometry_user_dims_initially_empty(self, fresh_db):
        """No user dimensions exist in a fresh DB."""
        constants = get_constants_dict(fresh_db)
        geom = compute_geometry(constants)
        assert geom["user_dimensions"] == []
