"""Tests for Phase 7 tools: opening width editing, delete, draw, rotate, etc.

Covers: DT-7, SEL-10, TL-22, TL-23, SEL-4, TL-15–27.
"""
import json

import pytest

from app.database import get_db, init_db, get_constant_value
from tests.test_zapp_conftest import fresh_db, app_client  # noqa: F401


# ── DT-7 + SEL-10: Opening width editing ────────────────────────────

class TestDT7OpeningWidth:
    """DT-7: Opening width SHALL be editable via the controlling constant."""

    def test_outer_opening_width_constant_exists(self, fresh_db):
        """O1_WIDTH through O6_WIDTH exist in the database."""
        for name in ["O1_WIDTH", "O2_WIDTH", "O3_WIDTH", "O5_WIDTH", "O6_WIDTH"]:
            v = get_constant_value(name, fresh_db)
            assert v is not None, f"{name} should exist"
            assert v > 0

    def test_outer_opening_half_width_constant_exists(self, fresh_db):
        """O4_HALF_WIDTH through O11_HALF_WIDTH exist in the database."""
        for name in ["O4_HALF_WIDTH", "O7_HALF_WIDTH", "O8_HALF_WIDTH",
                      "O8A_HALF_WIDTH", "O9_HALF_WIDTH", "O10_HALF_WIDTH",
                      "O11_HALF_WIDTH"]:
            v = get_constant_value(name, fresh_db)
            assert v is not None, f"{name} should exist"
            assert v > 0

    def test_rough_opening_width_constants_exist(self, fresh_db):
        """IW-based RO width constants exist in the database."""
        for name in ["IW1_RO_WIDTH", "IW2_RO_WIDTH", "RO3_WIDTH",
                      "IW4_RO_WIDTH", "IW9_RO_WIDTH", "IW11_RO_WIDTH",
                      "IW6_RO_WIDTH"]:
            v = get_constant_value(name, fresh_db)
            assert v is not None, f"{name} should exist"
            assert v > 0

    def test_update_outer_opening_width_via_api(self, app_client):
        """PUT O1_WIDTH changes value and triggers geometry recompute."""
        new_val = 24.0 / 12.0  # 24 inches in feet
        resp = app_client.put(
            "/api/constants/O1_WIDTH",
            data=json.dumps({"value": new_val}),
            content_type="application/json",
        )
        assert resp.status_code == 200
        data = resp.get_json()
        assert data["ok"] is True
        assert abs(data["value"] - new_val) < 1e-9

    def test_update_half_width_via_api(self, app_client):
        """PUT O7_HALF_WIDTH changes value."""
        new_val = 40.0 / 12.0
        resp = app_client.put(
            "/api/constants/O7_HALF_WIDTH",
            data=json.dumps({"value": new_val}),
            content_type="application/json",
        )
        assert resp.status_code == 200
        data = resp.get_json()
        assert abs(data["value"] - new_val) < 1e-9

    def test_update_rough_opening_width_via_api(self, app_client):
        """PUT IW1_RO_WIDTH changes value."""
        new_val = 42.0 / 12.0
        resp = app_client.put(
            "/api/constants/IW1_RO_WIDTH",
            data=json.dumps({"value": new_val}),
            content_type="application/json",
        )
        assert resp.status_code == 200
        data = resp.get_json()
        assert abs(data["value"] - new_val) < 1e-9

    def test_width_change_affects_geometry(self, app_client):
        """Changing O3_WIDTH should change O3's polygon in geometry output."""
        # Get original geometry
        g1 = app_client.get("/api/geometry").get_json()
        o3_orig = [o for o in g1["outer_openings"] if o["name"] == "O3"][0]
        orig_poly = o3_orig["poly"]

        # Change O3 width
        new_val = 40.0 / 12.0
        app_client.put(
            "/api/constants/O3_WIDTH",
            data=json.dumps({"value": new_val}),
            content_type="application/json",
        )

        # Get updated geometry
        g2 = app_client.get("/api/geometry").get_json()
        o3_new = [o for o in g2["outer_openings"] if o["name"] == "O3"][0]
        new_poly = o3_new["poly"]

        # Polygon should differ
        assert orig_poly != new_poly

    def test_undo_restores_width(self, app_client):
        """Undo after width change restores original value."""
        # Get original
        orig = app_client.get("/api/constants").get_json()
        orig_val = [c for c in orig if c["name"] == "O1_WIDTH"][0]["value"]

        # Change
        app_client.put(
            "/api/constants/O1_WIDTH",
            data=json.dumps({"value": 3.0}),
            content_type="application/json",
        )

        # Undo
        resp = app_client.post("/api/undo")
        assert resp.status_code == 200

        # Verify restored
        after = app_client.get("/api/constants").get_json()
        after_val = [c for c in after if c["name"] == "O1_WIDTH"][0]["value"]
        assert abs(after_val - orig_val) < 1e-9


# ── TL-22 + TL-23: Delete element ───────────────────────────────────

class TestTL22Delete:
    """TL-22: Delete key SHALL remove selected elements."""

    def test_delete_custom_element(self, app_client):
        """Create a custom element and delete it via API."""
        # Create
        resp = app_client.post(
            "/api/elements",
            data=json.dumps({"type": "furniture", "name": "test_chair",
                             "properties": {"source": "placed"}}),
            content_type="application/json",
        )
        assert resp.status_code == 201
        eid = resp.get_json()["id"]

        # Delete
        resp = app_client.delete(f"/api/elements/{eid}")
        assert resp.status_code == 200
        assert resp.get_json()["ok"] is True

        # Verify gone
        resp = app_client.get("/api/elements")
        names = [e["name"] for e in resp.get_json()]
        assert "test_chair" not in names

    def test_delete_nonexistent_returns_404(self, app_client):
        resp = app_client.delete("/api/elements/99999")
        assert resp.status_code == 404

    def test_delete_and_undo_restores(self, app_client):
        """Undo after delete restores the element."""
        # Create
        resp = app_client.post(
            "/api/elements",
            data=json.dumps({"type": "furniture", "name": "undo_test",
                             "properties": {"source": "placed"}}),
            content_type="application/json",
        )
        eid = resp.get_json()["id"]

        # Delete
        app_client.delete(f"/api/elements/{eid}")

        # Undo
        resp = app_client.post("/api/undo")
        assert resp.status_code == 200

        # Verify restored
        resp = app_client.get("/api/elements")
        names = [e["name"] for e in resp.get_json()]
        assert "undo_test" in names


class TestTL23CascadeDelete:
    """TL-23: Deleting a wall SHALL cascade-delete hosted openings."""

    def test_cascade_deletes_hosted_opening(self, app_client):
        """Deleting an IW wall also removes its hosted opening element."""
        # Create a custom opening hosted by IW1
        resp = app_client.post(
            "/api/elements",
            data=json.dumps({
                "type": "opening", "name": "test_opening",
                "properties": {"host_wall": "IW1", "source": "custom"},
            }),
            content_type="application/json",
        )
        assert resp.status_code == 201

        # Get IW1's element ID
        resp = app_client.get("/api/elements")
        elements = resp.get_json()
        iw1 = [e for e in elements if e["name"] == "IW1"][0]

        # Delete IW1 — should cascade
        resp = app_client.delete(f"/api/elements/{iw1['id']}")
        assert resp.status_code == 200
        deleted_ids = resp.get_json()["deleted"]
        assert len(deleted_ids) >= 2  # IW1 + test_opening

        # Verify opening is gone
        resp = app_client.get("/api/elements")
        names = [e["name"] for e in resp.get_json()]
        assert "test_opening" not in names


# ── TL-15–17: Draw Wall tool ────────────────────────────────────────

class TestTL15DrawWall:
    """TL-15: Two-click wall drawing creates a custom wall element."""

    def test_create_drawn_wall(self, app_client):
        """POST a drawn wall with start/end/thickness/poly."""
        start = [0.0, 0.0]
        end = [5.0, 0.0]
        thickness = 4.0 / 12.0
        # Rectangle from start/end/thickness
        ht = thickness / 2
        poly = [
            [0.0, ht], [5.0, ht], [5.0, -ht], [0.0, -ht]
        ]
        resp = app_client.post(
            "/api/elements",
            data=json.dumps({
                "type": "wall", "name": "CW1",
                "properties": {
                    "source": "drawn",
                    "start": start, "end": end,
                    "thickness": thickness,
                    "poly": poly,
                },
            }),
            content_type="application/json",
        )
        assert resp.status_code == 201
        data = resp.get_json()
        assert data["name"] == "CW1"
        assert data["type"] == "wall"

    def test_drawn_wall_properties(self, app_client):
        """Created drawn wall has correct properties."""
        start = [1.0, 2.0]
        end = [4.0, 6.0]
        resp = app_client.post(
            "/api/elements",
            data=json.dumps({
                "type": "wall", "name": "CW_test",
                "properties": {
                    "source": "drawn",
                    "start": start, "end": end,
                    "thickness": 0.5,
                    "poly": [[0, 0], [1, 0], [1, 1], [0, 1]],
                },
            }),
            content_type="application/json",
        )
        assert resp.status_code == 201
        eid = resp.get_json()["id"]

        # Fetch and verify properties
        resp = app_client.get("/api/elements")
        el = [e for e in resp.get_json() if e["id"] == eid][0]
        props = json.loads(el["properties"]) if isinstance(el["properties"], str) else el["properties"]
        assert props["source"] == "drawn"
        assert props["start"] == start
        assert props["end"] == end
        assert abs(props["thickness"] - 0.5) < 1e-9

    def test_update_drawn_wall_thickness(self, app_client):
        """TL-16: PUT updates thickness and poly."""
        resp = app_client.post(
            "/api/elements",
            data=json.dumps({
                "type": "wall", "name": "CW_thick",
                "properties": {
                    "source": "drawn",
                    "start": [0, 0], "end": [3, 0],
                    "thickness": 4.0 / 12.0,
                    "poly": [[0, 0.167], [3, 0.167], [3, -0.167], [0, -0.167]],
                },
            }),
            content_type="application/json",
        )
        eid = resp.get_json()["id"]

        # Update thickness
        new_thick = 6.0 / 12.0
        new_poly = [[0, 0.25], [3, 0.25], [3, -0.25], [0, -0.25]]
        resp = app_client.put(
            f"/api/elements/{eid}",
            data=json.dumps({
                "properties": {
                    "source": "drawn",
                    "start": [0, 0], "end": [3, 0],
                    "thickness": new_thick,
                    "poly": new_poly,
                },
            }),
            content_type="application/json",
        )
        assert resp.status_code == 200
        props = resp.get_json()["properties"]
        if isinstance(props, str):
            props = json.loads(props)
        assert abs(props["thickness"] - new_thick) < 1e-9

    def test_drawn_wall_undo(self, app_client):
        """Undo after creating drawn wall removes it."""
        resp = app_client.post(
            "/api/elements",
            data=json.dumps({
                "type": "wall", "name": "CW_undo",
                "properties": {
                    "source": "drawn", "start": [0, 0], "end": [1, 0],
                    "thickness": 0.333, "poly": [[0, 0], [1, 0], [1, 0.333], [0, 0.333]],
                },
            }),
            content_type="application/json",
        )
        assert resp.status_code == 201

        # Undo
        resp = app_client.post("/api/undo")
        assert resp.status_code == 200

        # Verify removed
        resp = app_client.get("/api/elements")
        names = [e["name"] for e in resp.get_json()]
        assert "CW_undo" not in names


# ── TL-18–21: Add Element / Opening tools ────────────────────────────

class TestTL18PlaceFurniture:
    """TL-18: Placing furniture from catalog creates a DB element."""

    def test_create_placed_furniture(self, app_client):
        resp = app_client.post(
            "/api/elements",
            data=json.dumps({
                "type": "furniture", "name": "custom_bed",
                "properties": {
                    "source": "placed",
                    "center": [0.0, 0.0],
                    "width": 76.0 / 12.0, "depth": 80.0 / 12.0,
                    "rotation": 0,
                    "poly": [[-3.17, 3.33], [3.17, 3.33], [3.17, -3.33], [-3.17, -3.33]],
                    "shape": "rect",
                },
            }),
            content_type="application/json",
        )
        assert resp.status_code == 201
        assert resp.get_json()["type"] == "furniture"

    def test_create_placed_appliance(self, app_client):
        resp = app_client.post(
            "/api/elements",
            data=json.dumps({
                "type": "appliance", "name": "custom_stove",
                "properties": {
                    "source": "placed",
                    "center": [5.0, 5.0],
                    "width": 30.0 / 12.0, "depth": 25.0 / 12.0,
                    "rotation": 0,
                    "poly": [[3.75, 6.04], [6.25, 6.04], [6.25, 3.96], [3.75, 3.96]],
                    "shape": "rect",
                },
            }),
            content_type="application/json",
        )
        assert resp.status_code == 201
        assert resp.get_json()["type"] == "appliance"

    def test_create_placed_fixture(self, app_client):
        resp = app_client.post(
            "/api/elements",
            data=json.dumps({
                "type": "fixture", "name": "custom_toilet",
                "properties": {
                    "source": "placed",
                    "center": [-5.0, -5.0],
                    "width": 15.0 / 12.0, "depth": 28.0 / 12.0,
                    "rotation": 0,
                    "poly": [[-5.625, -3.833], [-4.375, -3.833],
                             [-4.375, -6.167], [-5.625, -6.167]],
                    "shape": "rect",
                },
            }),
            content_type="application/json",
        )
        assert resp.status_code == 201
        assert resp.get_json()["type"] == "fixture"

    def test_placed_element_undo(self, app_client):
        resp = app_client.post(
            "/api/elements",
            data=json.dumps({
                "type": "furniture", "name": "undo_place",
                "properties": {"source": "placed", "center": [0, 0],
                               "poly": [[0, 0], [1, 0], [1, 1], [0, 1]]},
            }),
            content_type="application/json",
        )
        assert resp.status_code == 201

        resp = app_client.post("/api/undo")
        assert resp.status_code == 200

        resp = app_client.get("/api/elements")
        names = [e["name"] for e in resp.get_json()]
        assert "undo_place" not in names


# ── TL-24: Rotate element ────────────────────────────────────────────

class TestTL24Rotate:
    """TL-24: R key SHALL open rotation dialog for placed/drawn elements."""

    def test_rotate_placed_element(self, app_client):
        """PUT with rotation updates properties."""
        # Create placed element
        resp = app_client.post(
            "/api/elements",
            data=json.dumps({
                "type": "furniture", "name": "rot_test",
                "properties": {
                    "source": "placed", "center": [0, 0],
                    "width": 2, "depth": 1, "rotation": 0,
                    "poly": [[-1, 0.5], [1, 0.5], [1, -0.5], [-1, -0.5]],
                },
            }),
            content_type="application/json",
        )
        eid = resp.get_json()["id"]

        # Update rotation to 90
        resp = app_client.put(
            f"/api/elements/{eid}",
            data=json.dumps({
                "properties": {
                    "source": "placed", "center": [0, 0],
                    "width": 2, "depth": 1, "rotation": 90,
                    "poly": [[-0.5, -1], [-0.5, 1], [0.5, 1], [0.5, -1]],
                },
            }),
            content_type="application/json",
        )
        assert resp.status_code == 200
        props = resp.get_json()["properties"]
        if isinstance(props, str):
            props = json.loads(props)
        assert props["rotation"] == 90

    def test_rotate_undo(self, app_client):
        """Undo after rotation restores original angle."""
        resp = app_client.post(
            "/api/elements",
            data=json.dumps({
                "type": "furniture", "name": "rot_undo",
                "properties": {"source": "placed", "rotation": 0,
                               "poly": [[0, 0], [1, 0], [1, 1], [0, 1]]},
            }),
            content_type="application/json",
        )
        eid = resp.get_json()["id"]

        # Rotate
        app_client.put(
            f"/api/elements/{eid}",
            data=json.dumps({
                "properties": {"source": "placed", "rotation": 45,
                               "poly": [[0, 0], [1, 0], [1, 1], [0, 1]]},
            }),
            content_type="application/json",
        )

        # Undo
        resp = app_client.post("/api/undo")
        assert resp.status_code == 200

        # Check restored
        resp = app_client.get("/api/elements")
        el = [e for e in resp.get_json() if e["name"] == "rot_undo"][0]
        props = json.loads(el["properties"]) if isinstance(el["properties"], str) else el["properties"]
        assert props["rotation"] == 0


# ── TL-26–28: Shape editor ───────────────────────────────────────────

class TestTL26ShapeEditor:
    """TL-26: Shape editor CRUD via API."""

    def test_get_shapes_returns_seeded(self, app_client):
        """GET /api/shapes returns seeded shapes."""
        resp = app_client.get("/api/shapes")
        assert resp.status_code == 200
        shapes = resp.get_json()
        assert isinstance(shapes, list)
        assert len(shapes) > 0  # seeded shapes exist

    def test_create_shape(self, app_client):
        """POST /api/shapes creates a new shape."""
        resp = app_client.post(
            "/api/shapes",
            data=json.dumps({
                "name": "test_triangle",
                "poly_json": [[0, 1], [1, -1], [-1, -1]],
                "description": "Unit triangle",
            }),
            content_type="application/json",
        )
        assert resp.status_code == 201
        data = resp.get_json()
        assert data["name"] == "test_triangle"

    def test_update_shape(self, app_client):
        """PUT /api/shapes/<name> updates poly."""
        # Create first
        app_client.post(
            "/api/shapes",
            data=json.dumps({
                "name": "test_update",
                "poly_json": [[0, 0], [1, 0], [0, 1]],
            }),
            content_type="application/json",
        )

        # Update
        resp = app_client.put(
            "/api/shapes/test_update",
            data=json.dumps({
                "poly_json": [[0, 0], [2, 0], [0, 2]],
                "description": "Updated",
            }),
            content_type="application/json",
        )
        assert resp.status_code == 200
        data = resp.get_json()
        assert data["description"] == "Updated"

    def test_update_nonexistent_returns_404(self, app_client):
        resp = app_client.put(
            "/api/shapes/nonexistent_xyz",
            data=json.dumps({"description": "nope"}),
            content_type="application/json",
        )
        assert resp.status_code == 404

    def test_shape_assignment_via_element(self, app_client):
        """TL-26: Assigning a shape to an element via properties update."""
        # Create shape
        app_client.post(
            "/api/shapes",
            data=json.dumps({
                "name": "diamond",
                "poly_json": [[0, 1], [1, 0], [0, -1], [-1, 0]],
            }),
            content_type="application/json",
        )

        # Create placed element
        resp = app_client.post(
            "/api/elements",
            data=json.dumps({
                "type": "furniture", "name": "shape_test",
                "properties": {
                    "source": "placed", "center": [0, 0],
                    "width": 2, "depth": 2, "rotation": 0,
                    "poly": [[-1, 1], [1, 1], [1, -1], [-1, -1]],
                    "shape": "rect",
                },
            }),
            content_type="application/json",
        )
        eid = resp.get_json()["id"]

        # Update with diamond shape
        resp = app_client.put(
            f"/api/elements/{eid}",
            data=json.dumps({
                "properties": {
                    "source": "placed", "center": [0, 0],
                    "width": 2, "depth": 2, "rotation": 0,
                    "poly": [[0, 1], [1, 0], [0, -1], [-1, 0]],
                    "shape": "diamond",
                },
            }),
            content_type="application/json",
        )
        assert resp.status_code == 200
        props = resp.get_json()["properties"]
        if isinstance(props, str):
            props = json.loads(props)
        assert props["shape"] == "diamond"


# =========================================================================
# DIS-7: User Dimensions Toggle
# =========================================================================

class TestDIS7UserDimsToggle:
    """DIS-7: Separate User Dims toggle."""

    def test_user_dims_checkbox_in_template(self, app_client):
        """index.html contains the show-user-dims checkbox."""
        resp = app_client.get("/")
        assert resp.status_code == 200
        html = resp.data.decode()
        assert 'id="show-user-dims"' in html

    def test_geometry_has_builtin_and_user_dims(self, app_client):
        """Geometry response contains both builtin and user-created dims."""
        # Create a user dimension
        app_client.post("/api/elements", json={
            "type": "dimension",
            "name": "test_dim1",
            "properties": {
                "start": [-10, -5], "end": [-5, -5], "offset": 0,
                "start_anchor": None, "end_anchor": None,
            },
        })
        resp = app_client.get("/api/geometry")
        data = resp.get_json()
        dims = data.get("user_dimensions", [])
        sources = set()
        for d in dims:
            p = d.get("properties", {})
            if isinstance(p, str):
                p = json.loads(p)
            sources.add(p.get("source", "user"))
        # Should have both builtin (seeded) and user-created
        assert "builtin" in sources
        assert "user" in sources or None in sources


# =========================================================================
# TL-17: Wall Endpoint Drag Handles
# =========================================================================

class TestTL17EndpointDrag:
    """TL-17: Updating drawn wall endpoints via API."""

    def _create_drawn_wall(self, client):
        """Helper: create a drawn wall and return its record."""
        start = [-10, -5]
        end = [-5, -5]
        thickness = 4.0 / 12
        dx = end[0] - start[0]
        dy = end[1] - start[1]
        length = (dx**2 + dy**2) ** 0.5
        px = -dy / length * (thickness / 2)
        py = dx / length * (thickness / 2)
        poly = [
            [start[0] + px, start[1] + py],
            [end[0] + px, end[1] + py],
            [end[0] - px, end[1] - py],
            [start[0] - px, start[1] - py],
        ]
        resp = client.post("/api/elements", json={
            "type": "wall", "name": "CW_test",
            "properties": {
                "source": "drawn", "start": start, "end": end,
                "thickness": thickness, "poly": poly,
            },
        })
        assert resp.status_code in (200, 201)
        return resp.get_json()

    def test_update_drawn_wall_endpoint(self, app_client):
        """PUT new start/end updates coordinates."""
        elem = self._create_drawn_wall(app_client)
        props = json.loads(elem["properties"]) if isinstance(elem["properties"], str) else elem["properties"]
        props["end"] = [-3, -5]  # extend the wall
        resp = app_client.put(f"/api/elements/{elem['id']}", json={
            "properties": props,
        })
        assert resp.status_code == 200
        up = resp.get_json()
        up_props = json.loads(up["properties"]) if isinstance(up["properties"], str) else up["properties"]
        assert abs(up_props["end"][0] - (-3)) < 1e-9
        assert abs(up_props["end"][1] - (-5)) < 1e-9

    def test_endpoint_update_preserves_thickness(self, app_client):
        """Updating endpoint does not change thickness."""
        elem = self._create_drawn_wall(app_client)
        props = json.loads(elem["properties"]) if isinstance(elem["properties"], str) else elem["properties"]
        orig_thickness = props["thickness"]
        props["start"] = [-12, -5]
        resp = app_client.put(f"/api/elements/{elem['id']}", json={
            "properties": props,
        })
        up = resp.get_json()
        up_props = json.loads(up["properties"]) if isinstance(up["properties"], str) else up["properties"]
        assert abs(up_props["thickness"] - orig_thickness) < 1e-9

    def test_endpoint_update_undo(self, app_client):
        """Undo after endpoint update restores original coordinates."""
        elem = self._create_drawn_wall(app_client)
        props = json.loads(elem["properties"]) if isinstance(elem["properties"], str) else elem["properties"]
        orig_end = props["end"][:]
        props["end"] = [-2, -3]
        app_client.put(f"/api/elements/{elem['id']}", json={
            "properties": props,
        })
        # Undo
        resp = app_client.post("/api/undo")
        assert resp.status_code == 200
        # Verify restored
        resp = app_client.get("/api/elements")
        elems = resp.get_json()
        restored = next(e for e in elems if e["name"] == "CW_test")
        rp = json.loads(restored["properties"]) if isinstance(restored["properties"], str) else restored["properties"]
        assert abs(rp["end"][0] - orig_end[0]) < 1e-9
        assert abs(rp["end"][1] - orig_end[1]) < 1e-9


# =========================================================================
# TL-21: Add Opening Tool
# =========================================================================

class TestTL21AddOpening:
    """TL-21: Creating openings on walls via API."""

    def _create_drawn_wall(self, client, name="CW_host"):
        """Create a drawn wall to host openings."""
        resp = client.post("/api/elements", json={
            "type": "wall", "name": name,
            "properties": {
                "source": "drawn",
                "start": [-10, 0], "end": [-5, 0],
                "thickness": 4.0 / 12,
                "poly": [[-10, 1/6], [-5, 1/6], [-5, -1/6], [-10, -1/6]],
            },
        })
        assert resp.status_code in (200, 201)
        return resp.get_json()

    def test_create_opening_on_wall(self, app_client):
        """POST opening with host_wall creates opening element."""
        self._create_drawn_wall(app_client)
        resp = app_client.post("/api/openings", json={
            "name": "UO1",
            "segment": "CW_host",
            "width": 32.0 / 12,
            "offset": 0,
            "host_wall": "CW_host",
        })
        assert resp.status_code == 201
        data = resp.get_json()
        assert data["type"] == "opening"
        assert data["name"] == "UO1"

    def test_opening_host_wall_reference(self, app_client):
        """Created opening has correct host_wall in properties."""
        self._create_drawn_wall(app_client)
        resp = app_client.post("/api/openings", json={
            "name": "UO2",
            "segment": "CW_host",
            "width": 24.0 / 12,
            "host_wall": "CW_host",
        })
        data = resp.get_json()
        props = json.loads(data["properties"]) if isinstance(data["properties"], str) else data["properties"]
        assert props["host_wall"] == "CW_host"

    def test_delete_host_wall_cascades_opening(self, app_client):
        """Deleting wall cascades to hosted openings."""
        wall = self._create_drawn_wall(app_client)
        app_client.post("/api/openings", json={
            "name": "UO3",
            "segment": "CW_host",
            "width": 24.0 / 12,
            "host_wall": "CW_host",
        })
        # Verify opening exists
        elems = app_client.get("/api/elements").get_json()
        assert any(e["name"] == "UO3" for e in elems)
        # Delete the wall
        resp = app_client.delete(f"/api/elements/{wall['id']}")
        assert resp.status_code == 200
        # Opening should be gone too
        elems = app_client.get("/api/elements").get_json()
        assert not any(e["name"] == "UO3" for e in elems)


# =========================================================================
# SEL-15: Constant Dependency Highlighting
# =========================================================================

class TestSEL15ConstantDeps:
    """SEL-15: Reverse constant → element dependency mapping."""

    def test_reverse_iw_map_completeness(self):
        """Every non-None IW_CONSTANT_MAP entry has a reverse mapping."""
        from app.elements import IW_CONSTANT_MAP, CONSTANT_TO_IW
        for iw, cname in IW_CONSTANT_MAP.items():
            if cname is None:
                continue
            assert cname in CONSTANT_TO_IW, f"{cname} missing from CONSTANT_TO_IW"
            assert iw in CONSTANT_TO_IW[cname], f"{iw} missing from CONSTANT_TO_IW[{cname}]"

    def test_specific_constant_mapping(self):
        """IW1_OFFSET_FROM_W9 maps back to IW1."""
        from app.elements import CONSTANT_TO_IW
        assert "IW1_OFFSET_FROM_W9" in CONSTANT_TO_IW
        assert "IW1" in CONSTANT_TO_IW["IW1_OFFSET_FROM_W9"]

    def test_hosted_openings_second_order(self):
        """IW9 hosts RO3 and RO7 (second-order deps)."""
        from app.elements import IW_HOSTED_OPENINGS
        hosted = IW_HOSTED_OPENINGS.get("IW9", [])
        assert "RO3" in hosted
        assert "RO7" in hosted
