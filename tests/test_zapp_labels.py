"""Tests for Phase 8: Labels, Dimensions, Annotations, and Unified Dimensions.

Covers: TL-11–14, LABEL-1–4, DIS-7 (9 reqs), plus unified dimension tests.
"""
import json
import math

import pytest

from app.database import (
    get_db, init_db, get_all_elements, get_element,
    create_element, update_element, delete_element,
    reset_elements,
)
from app.labels import (
    seed_room_labels, seed_builtin_dimensions,
    next_dimension_name, next_label_name,
    ROOM_LABEL_NAMES, BUILTIN_DIMENSIONS,
)
from app.engine import (
    compute_geometry, _resolve_anchor, _face_midpoint,
    _resolve_point_spec, _resolve_dir_spec,
)
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
        geom = compute_geometry(constants, db_path=fresh_db)
        assert "user_dimensions" in geom
        assert isinstance(geom["user_dimensions"], list)

    def test_geometry_has_label_elements(self, fresh_db):
        constants = get_constants_dict(fresh_db)
        geom = compute_geometry(constants, db_path=fresh_db)
        assert "label_elements" in geom
        assert len(geom["label_elements"]) == 11  # room labels

    def test_label_elements_have_centroid(self, fresh_db):
        constants = get_constants_dict(fresh_db)
        geom = compute_geometry(constants, db_path=fresh_db)
        for le in geom["label_elements"]:
            if le["properties"]["source"] == "room":
                assert "centroid" in le
                assert "pos" in le

    def test_geometry_has_builtin_dims(self, fresh_db):
        """Fresh DB has seeded dimensions in user_dimensions (no source distinction)."""
        constants = get_constants_dict(fresh_db)
        geom = compute_geometry(constants, db_path=fresh_db)
        seeded = [d for d in geom["user_dimensions"]
                  if d["name"].startswith("dim")]
        # Standard variant: 20 normal dims (no dim12bare, dim20, dim21)
        assert len(seeded) == 20


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


# ── Builtin dimension seeding ──────────────────────────────────────

class TestBuiltinDimensionSeeding:
    """Built-in dimensions are seeded as elements on fresh DB."""

    def test_builtin_dims_count(self, fresh_db):
        elements = get_all_elements(fresh_db)
        dims = [e for e in elements if e["type"] == "dimension"]
        assert len(dims) == len(BUILTIN_DIMENSIONS)

    def test_builtin_dim_names(self, fresh_db):
        elements = get_all_elements(fresh_db)
        dim_names = {e["name"] for e in elements if e["type"] == "dimension"}
        expected = {d["name"] for d in BUILTIN_DIMENSIONS}
        assert dim_names == expected

    def test_builtin_dim_has_anchors(self, fresh_db):
        elements = get_all_elements(fresh_db)
        for e in elements:
            if e["type"] != "dimension":
                continue
            props = json.loads(e["properties"]) if isinstance(e["properties"], str) else e["properties"]
            assert "start_anchor" in props
            assert "end_anchor" in props

    def test_seed_idempotent(self, fresh_db):
        with get_db(fresh_db) as conn:
            seed_builtin_dimensions(conn)
        elements = get_all_elements(fresh_db)
        dims = [e for e in elements if e["type"] == "dimension"]
        assert len(dims) == len(BUILTIN_DIMENSIONS)

    def test_dim12a_variant_null(self, fresh_db):
        elements = get_all_elements(fresh_db)
        dim12a = next(e for e in elements if e["name"] == "dim12a")
        assert dim12a["variant"] is None

    def test_dim12bare_variant_bare(self, fresh_db):
        elements = get_all_elements(fresh_db)
        dim12bare = next(e for e in elements if e["name"] == "dim12bare")
        assert dim12bare["variant"] == "bare"

    def test_dim20_dim21_variant_bare(self, fresh_db):
        elements = get_all_elements(fresh_db)
        for name in ("dim20", "dim21"):
            dim = next(e for e in elements if e["name"] == name)
            assert dim["variant"] == "bare", f"{name} should be bare-only"

    def test_reset_restores_builtin_dims(self, fresh_db):
        """Deleting a builtin dim and resetting restores it."""
        elements = get_all_elements(fresh_db)
        dim01 = next(e for e in elements if e["name"] == "dim01")
        delete_element(dim01["id"], fresh_db)
        # Verify it's gone
        assert get_element(dim01["id"], fresh_db) is None
        # Reset restores it
        reset_elements(fresh_db)
        elements = get_all_elements(fresh_db)
        dim_names = {e["name"] for e in elements if e["type"] == "dimension"}
        assert "dim01" in dim_names


# ── Point spec resolution ──────────────────────────────────────────

MOCK_GEOM = {
    "points": {
        "W9": [5.0, 3.0], "F9": [6.0, 4.0], "F11": [8.0, 4.0],
        "W2": [-2.0, 0.0], "F5": [-2.0, 5.0],
        "W14": [3.0, -5.0], "W15": [4.0, -5.0],
        "W18": [-1.0, -3.0], "W1": [-1.0, 3.0],
        "C7": [4.0, 2.0],
    },
    "interior_walls": {
        "IW1": {"poly": [[0, 0], [10, 0], [10, 0.667], [0, 0.667]]},
        "IW2": {"poly": [[-2.0, 1.0], [-1.5, 1.0], [-1.5, 4.0], [-2.0, 4.0]]},
        "IW9": {"poly": [[1, 1], [1.5, 1], [1.5, 3], [1, 3]]},
        "IW11": {"poly": [[3, 1], [3.5, 1], [3.5, 3], [3, 3]]},
    },
    "outer_openings": [
        {"name": "O1", "poly": [[2, 0], [5, 0], [5, 0.5], [2, 0.5]]},
        {"name": "O6", "poly": [[2, 0], [5, 0], [5, 0.5], [2, 0.5]]},
    ],
    "rough_openings": [
        {"name": "RO1", "poly": [[3, 1], [4, 1], [4, 1.667], [3, 1.667]]},
    ],
    "radii": {"R_a7": 3.0, "R_a11": 4.0},
}


class TestResolvePointSpec:
    """Unit tests for _resolve_point_spec."""

    def test_string_point(self):
        assert _resolve_point_spec("W9", MOCK_GEOM) == [5.0, 3.0]

    def test_string_point_missing(self):
        assert _resolve_point_spec("ZZZZ", MOCK_GEOM) is None

    def test_face_mid(self):
        result = _resolve_point_spec({"face_mid": "IW1", "face": "north"}, MOCK_GEOM)
        assert result is not None
        assert abs(result[0] - 5.0) < 0.001
        assert abs(result[1] - 0.667) < 0.001

    def test_opening_face_mid(self):
        result = _resolve_point_spec(
            {"opening_face_mid": "O1", "face": "south"}, MOCK_GEOM)
        assert result is not None
        assert abs(result[0] - 3.5) < 0.001
        assert abs(result[1] - 0.0) < 0.001

    def test_opening_centroid(self):
        result = _resolve_point_spec({"opening_centroid": "O1"}, MOCK_GEOM)
        assert result is not None
        assert abs(result[0] - 3.5) < 0.001
        assert abs(result[1] - 0.25) < 0.001

    def test_midpoint(self):
        result = _resolve_point_spec({"midpoint": ["F9", "F11"]}, MOCK_GEOM)
        assert result is not None
        assert abs(result[0] - 7.0) < 0.001
        assert abs(result[1] - 4.0) < 0.001

    def test_offset(self):
        result = _resolve_point_spec(
            {"offset": "W9", "dir": "east", "dist": 2.0}, MOCK_GEOM)
        assert result is not None
        assert abs(result[0] - 7.0) < 0.001
        assert abs(result[1] - 3.0) < 0.001

    def test_offset_north(self):
        result = _resolve_point_spec(
            {"offset": "W9", "dir": "north", "dist": -1.0}, MOCK_GEOM)
        assert result is not None
        assert abs(result[0] - 5.0) < 0.001
        assert abs(result[1] - 2.0) < 0.001

    def test_arc_point_east(self):
        # Circle at C7=(4,2), R=3, reference at N=2 (same height as center)
        # dn = 0, de = sqrt(9-0) = 3, east side → point = (4+3, 2) = (7, 2)
        result = _resolve_point_spec({
            "arc_point": {
                "center": "C7", "radius_key": "R_a7",
                "reference": "C7",  # same point gives dn=0
                "side": "east",
            }
        }, MOCK_GEOM)
        assert result is not None
        assert abs(result[0] - 7.0) < 0.001
        assert abs(result[1] - 2.0) < 0.001

    def test_arc_point_west(self):
        result = _resolve_point_spec({
            "arc_point": {
                "center": "C7", "radius_key": "R_a7",
                "reference": "C7",
                "side": "west",
            }
        }, MOCK_GEOM)
        assert result is not None
        assert abs(result[0] - 1.0) < 0.001  # 4 - 3

    def test_none_spec(self):
        assert _resolve_point_spec(None, MOCK_GEOM) is None


# ── Direction spec resolution ──────────────────────────────────────

class TestResolveDirSpec:
    """Unit tests for _resolve_dir_spec."""

    def test_east(self):
        assert _resolve_dir_spec("east", MOCK_GEOM) == [1.0, 0.0]

    def test_north(self):
        assert _resolve_dir_spec("north", MOCK_GEOM) == [0.0, 1.0]

    def test_face_along(self):
        # IW1 south face: poly[0]=[0,0] → poly[1]=[10,0] → along = [1, 0]
        result = _resolve_dir_spec(
            {"face_along": "IW1", "face": "south"}, MOCK_GEOM)
        assert result is not None
        assert abs(result[0] - 1.0) < 0.001
        assert abs(result[1] - 0.0) < 0.001

    def test_face_perp(self):
        # IW1 south face: along = [1, 0], perp (right-hand) = [0, -1]
        result = _resolve_dir_spec(
            {"face_perp": "IW1", "face": "south"}, MOCK_GEOM)
        assert result is not None
        assert abs(result[0] - 0.0) < 0.001
        assert abs(result[1] - (-1.0)) < 0.001

    def test_segment(self):
        # W14=[3,-5] → W15=[4,-5] → along = [1, 0]
        result = _resolve_dir_spec(
            {"segment": ["W14", "W15"]}, MOCK_GEOM)
        assert result is not None
        assert abs(result[0] - 1.0) < 0.001
        assert abs(result[1] - 0.0) < 0.001

    def test_segment_perp(self):
        # W14→W15 along [1, 0], perp (right-hand) = [0, -1]
        result = _resolve_dir_spec(
            {"segment_perp": ["W14", "W15"]}, MOCK_GEOM)
        assert result is not None
        assert abs(result[0] - 0.0) < 0.001
        assert abs(result[1] - (-1.0)) < 0.001

    def test_none_spec(self):
        assert _resolve_dir_spec(None, MOCK_GEOM) is None

    def test_unknown_string(self):
        assert _resolve_dir_spec("diagonal", MOCK_GEOM) is None


# ── Line intersection anchor ──────────────────────────────────────

class TestLineIntersectionAnchor:
    """Tests for the line_intersection anchor type."""

    def test_perpendicular_lines(self):
        """Intersection of horizontal and vertical lines."""
        anchor = {
            "type": "line_intersection",
            "line1_point": "W9",      # [5, 3]
            "line1_dir": "east",       # horizontal at N=3
            "line2_point": "F9",       # [6, 4]
            "line2_dir": "north",      # vertical at E=6
        }
        result = _resolve_anchor(anchor, MOCK_GEOM)
        assert result is not None
        assert abs(result[0] - 6.0) < 0.001
        assert abs(result[1] - 3.0) < 0.001

    def test_face_based_intersection(self):
        """Intersection using face midpoints and face directions."""
        anchor = {
            "type": "line_intersection",
            "line1_point": {"face_mid": "IW1", "face": "north"},
            "line1_dir": {"face_along": "IW1", "face": "north"},
            "line2_point": "F9",
            "line2_dir": "north",
        }
        result = _resolve_anchor(anchor, MOCK_GEOM)
        assert result is not None
        # IW1 north face is at N=0.667, F9 at E=6.0
        assert abs(result[0] - 6.0) < 0.001
        assert abs(result[1] - 0.667) < 0.001

    def test_parallel_lines_returns_none(self):
        """Parallel lines should return None (no intersection)."""
        anchor = {
            "type": "line_intersection",
            "line1_point": "W9", "line1_dir": "east",
            "line2_point": "F9", "line2_dir": "east",
        }
        assert _resolve_anchor(anchor, MOCK_GEOM) is None

    def test_missing_point_returns_none(self):
        anchor = {
            "type": "line_intersection",
            "line1_point": "MISSING", "line1_dir": "east",
            "line2_point": "W9", "line2_dir": "north",
        }
        assert _resolve_anchor(anchor, MOCK_GEOM) is None


class TestComputedAnchor:
    """Tests for the computed anchor type."""

    def test_computed_midpoint(self):
        anchor = {
            "type": "computed",
            "spec": {"midpoint": ["W9", "F9"]},
        }
        result = _resolve_anchor(anchor, MOCK_GEOM)
        assert result is not None
        assert abs(result[0] - 5.5) < 0.001
        assert abs(result[1] - 3.5) < 0.001


# ── Variant filtering ──────────────────────────────────────────────

class TestDimensionVariantFiltering:
    """Variant filtering for builtin dimensions."""

    def test_standard_excludes_bare_only(self, fresh_db):
        constants = get_constants_dict(fresh_db)
        geom = compute_geometry(constants, variant="standard", db_path=fresh_db)
        dim_names = {d["name"] for d in geom["user_dimensions"]}
        assert "dim12a" in dim_names
        assert "dim12b" in dim_names
        assert "dim12bare" not in dim_names
        assert "dim20" not in dim_names
        assert "dim21" not in dim_names

    def test_bare_excludes_dim12a_dim12b(self, fresh_db):
        constants = get_constants_dict(fresh_db)
        geom = compute_geometry(constants, variant="bare", db_path=fresh_db)
        dim_names = {d["name"] for d in geom["user_dimensions"]}
        assert "dim12a" not in dim_names
        assert "dim12b" not in dim_names
        assert "dim12bare" in dim_names
        assert "dim20" in dim_names
        assert "dim21" in dim_names

    def test_bare_has_correct_count(self, fresh_db):
        """Bare variant: 20 normal dims - 2 (dim12a/b) + 3 (bare-only) = 21."""
        constants = get_constants_dict(fresh_db)
        geom = compute_geometry(constants, variant="bare", db_path=fresh_db)
        seeded = [d for d in geom["user_dimensions"]
                  if d["name"].startswith("dim")]
        assert len(seeded) == 21


class TestPropertiesVariantsFiltering:
    """properties.variants list overrides the variant column."""

    def test_variants_list_restricts_to_selected(self, fresh_db):
        create_element("dimension", "VF_TEST", {
            "source": "user", "start": [0, 0], "end": [5, 5],
            "variants": ["standard", "minik"],
        }, None, fresh_db)
        constants = get_constants_dict(fresh_db)
        std = compute_geometry(constants, "standard", db_path=fresh_db)
        mk = compute_geometry(constants, "minik", db_path=fresh_db)
        bare = compute_geometry(constants, "bare", db_path=fresh_db)
        std_names = {d["name"] for d in std["user_dimensions"]}
        mk_names = {d["name"] for d in mk["user_dimensions"]}
        bare_names = {d["name"] for d in bare["user_dimensions"]}
        assert "VF_TEST" in std_names
        assert "VF_TEST" in mk_names
        assert "VF_TEST" not in bare_names

    def test_empty_variants_list_hides_from_all(self, fresh_db):
        create_element("dimension", "VF_EMPTY", {
            "source": "user", "start": [0, 0], "end": [5, 5],
            "variants": [],
        }, None, fresh_db)
        constants = get_constants_dict(fresh_db)
        std = compute_geometry(constants, "standard", db_path=fresh_db)
        mk = compute_geometry(constants, "minik", db_path=fresh_db)
        assert "VF_EMPTY" not in {d["name"] for d in std["user_dimensions"]}
        assert "VF_EMPTY" not in {d["name"] for d in mk["user_dimensions"]}

    def test_variants_list_overrides_column(self, fresh_db):
        """If both properties.variants and variant column are set, list wins."""
        create_element("dimension", "VF_OVERRIDE", {
            "source": "user", "start": [0, 0], "end": [5, 5],
            "variants": ["bare"],
        }, "standard", fresh_db)  # column says standard, list says bare
        constants = get_constants_dict(fresh_db)
        std = compute_geometry(constants, "standard", db_path=fresh_db)
        bare = compute_geometry(constants, "bare", db_path=fresh_db)
        assert "VF_OVERRIDE" not in {d["name"] for d in std["user_dimensions"]}
        assert "VF_OVERRIDE" in {d["name"] for d in bare["user_dimensions"]}


# ── Geometry response with resolved builtin dims ───────────────────

class TestBuiltinDimGeometry:
    """Builtin dimensions resolve to valid coordinates."""

    def test_all_builtin_dims_resolve(self, fresh_db):
        """Every seeded dim has resolved start and end coordinates."""
        constants = get_constants_dict(fresh_db)
        geom = compute_geometry(constants, db_path=fresh_db)
        builtin = [d for d in geom["user_dimensions"]
                   if d["name"].startswith("dim")]
        for dim in builtin:
            p = dim["properties"]
            assert "start" in p, f"{dim['name']} missing start"
            assert "end" in p, f"{dim['name']} missing end"
            assert len(p["start"]) == 2, f"{dim['name']} start not [E,N]"
            assert len(p["end"]) == 2, f"{dim['name']} end not [E,N]"

    def test_dim01_has_positive_distance(self, fresh_db):
        constants = get_constants_dict(fresh_db)
        geom = compute_geometry(constants, db_path=fresh_db)
        dim01 = next(d for d in geom["user_dimensions"] if d["name"] == "dim01")
        p = dim01["properties"]
        dist = math.sqrt((p["end"][0] - p["start"][0])**2 +
                         (p["end"][1] - p["start"][1])**2)
        assert dist > 0.5  # should be several feet

    def test_radii_in_geometry(self, fresh_db):
        """Geometry response includes radii dict."""
        constants = get_constants_dict(fresh_db)
        geom = compute_geometry(constants, db_path=fresh_db)
        assert "radii" in geom
        assert "R_a7" in geom["radii"]

    def test_user_dim_still_works(self, fresh_db):
        """User-created dimensions still appear alongside builtins."""
        create_element("dimension", "UD1", {
            "source": "user", "start": [0, 0], "end": [5, 0],
            "offset": 0,
        }, None, fresh_db)
        constants = get_constants_dict(fresh_db)
        geom = compute_geometry(constants, db_path=fresh_db)
        names = {d["name"] for d in geom["user_dimensions"]}
        assert "UD1" in names
        assert "dim01" in names  # builtin still there
