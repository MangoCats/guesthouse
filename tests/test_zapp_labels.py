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
from app.engine import compute_geometry, _resolve_anchor, _face_midpoint
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


# ── Anchor resolution ───────────────────────────────────────────────

class TestFaceMidpoint:
    """Unit tests for _face_midpoint helper."""

    def test_south_face(self):
        poly = [[0, 0], [4, 0], [4, 2], [0, 2]]
        assert _face_midpoint(poly, "south") == [2.0, 0.0]

    def test_north_face(self):
        poly = [[0, 0], [4, 0], [4, 2], [0, 2]]
        assert _face_midpoint(poly, "north") == [2.0, 2.0]

    def test_east_face(self):
        poly = [[0, 0], [4, 0], [4, 2], [0, 2]]
        assert _face_midpoint(poly, "east") == [4.0, 1.0]

    def test_west_face(self):
        poly = [[0, 0], [4, 0], [4, 2], [0, 2]]
        assert _face_midpoint(poly, "west") == [0.0, 1.0]

    def test_invalid_face(self):
        poly = [[0, 0], [4, 0], [4, 2], [0, 2]]
        assert _face_midpoint(poly, "top") is None

    def test_short_poly(self):
        poly = [[0, 0], [4, 0]]
        assert _face_midpoint(poly, "south") is None


class TestResolveAnchor:
    """Unit tests for _resolve_anchor using mock geometry dicts."""

    MOCK_GEOM = {
        "points": {
            "W9": [5.0, 3.0],
            "F9": [6.0, 4.0],
        },
        "interior_walls": {
            "IW1": {"poly": [[0, 0], [10, 0], [10, 0.667], [0, 0.667]]},
            "IW2S": {"poly": [[1, 1], [1.5, 1], [1.5, 3], [1, 3]]},
        },
        "outer_openings": [
            {"name": "O6", "poly": [[2, 0], [5, 0], [5, 0.5], [2, 0.5]]},
        ],
        "rough_openings": [
            {"name": "RO1", "poly": [[3, 1], [4, 1], [4, 1.667], [3, 1.667]]},
        ],
    }

    def test_resolve_point(self):
        anchor = {"type": "point", "target": "W9"}
        assert _resolve_anchor(anchor, self.MOCK_GEOM) == [5.0, 3.0]

    def test_resolve_point_missing(self):
        anchor = {"type": "point", "target": "ZZZZ"}
        assert _resolve_anchor(anchor, self.MOCK_GEOM) is None

    def test_resolve_wall_face_north(self):
        anchor = {"type": "wall_face", "target": "IW1", "face": "north"}
        result = _resolve_anchor(anchor, self.MOCK_GEOM)
        assert result is not None
        assert abs(result[0] - 5.0) < 0.001  # midpoint of poly[2]→poly[3]
        assert abs(result[1] - 0.667) < 0.001

    def test_resolve_wall_face_south(self):
        anchor = {"type": "wall_face", "target": "IW1", "face": "south"}
        result = _resolve_anchor(anchor, self.MOCK_GEOM)
        assert result is not None
        assert abs(result[0] - 5.0) < 0.001
        assert abs(result[1] - 0.0) < 0.001

    def test_resolve_wall_face_missing_wall(self):
        anchor = {"type": "wall_face", "target": "IW99", "face": "north"}
        assert _resolve_anchor(anchor, self.MOCK_GEOM) is None

    def test_resolve_opening_face_outer(self):
        anchor = {"type": "opening_face", "target": "O6", "face": "south"}
        result = _resolve_anchor(anchor, self.MOCK_GEOM)
        assert result is not None
        assert abs(result[0] - 3.5) < 0.001
        assert abs(result[1] - 0.0) < 0.001

    def test_resolve_opening_face_rough(self):
        anchor = {"type": "opening_face", "target": "RO1", "face": "east"}
        result = _resolve_anchor(anchor, self.MOCK_GEOM)
        assert result is not None
        assert abs(result[0] - 4.0) < 0.001
        assert abs(result[1] - 1.3335) < 0.001

    def test_resolve_opening_face_missing(self):
        anchor = {"type": "opening_face", "target": "O99", "face": "south"}
        assert _resolve_anchor(anchor, self.MOCK_GEOM) is None

    def test_resolve_none_anchor(self):
        assert _resolve_anchor(None, self.MOCK_GEOM) is None

    def test_resolve_empty_anchor(self):
        assert _resolve_anchor({}, self.MOCK_GEOM) is None

    def test_resolve_unknown_type(self):
        anchor = {"type": "foobar", "target": "W9"}
        assert _resolve_anchor(anchor, self.MOCK_GEOM) is None


class TestDimensionAnchorCRUD:
    """Dimension elements with anchors via API."""

    def test_create_dimension_with_anchors(self, app_client):
        resp = app_client.post("/api/elements", json={
            "type": "dimension", "name": "UD1",
            "properties": {
                "source": "user", "start": [0, 0], "end": [5, 0],
                "offset": 0, "label_rotation": "parallel",
                "start_anchor": {"type": "point", "target": "W9"},
                "end_anchor": {"type": "wall_face", "target": "IW1", "face": "north"},
            },
        })
        assert resp.status_code == 201
        props = json.loads(resp.get_json()["properties"])
        assert props["start_anchor"]["type"] == "point"
        assert props["end_anchor"]["face"] == "north"

    def test_detach_anchor(self, app_client):
        resp = app_client.post("/api/elements", json={
            "type": "dimension", "name": "UD1",
            "properties": {
                "source": "user", "start": [0, 0], "end": [5, 0],
                "start_anchor": {"type": "point", "target": "W9"},
                "end_anchor": {"type": "point", "target": "F9"},
            },
        })
        eid = resp.get_json()["id"]
        # Detach start anchor
        resp = app_client.put(f"/api/elements/{eid}", json={
            "properties": {
                "source": "user", "start": [0, 0], "end": [5, 0],
                "end_anchor": {"type": "point", "target": "F9"},
            },
        })
        assert resp.status_code == 200
        props = json.loads(resp.get_json()["properties"])
        assert "start_anchor" not in props
        assert props["end_anchor"]["target"] == "F9"

    def test_dimension_without_anchors_unchanged(self, app_client):
        """Dimensions without anchors work as before."""
        resp = app_client.post("/api/elements", json={
            "type": "dimension", "name": "UD1",
            "properties": {
                "source": "user", "start": [1, 2], "end": [3, 4],
                "offset": 0.5,
            },
        })
        assert resp.status_code == 201
        props = json.loads(resp.get_json()["properties"])
        assert "start_anchor" not in props
        assert props["start"] == [1, 2]
