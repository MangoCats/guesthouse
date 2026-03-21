"""Tests for Phase 7 tools: opening width editing, delete, draw, rotate, etc.

Covers: DT-7, SEL-10, TL-22, TL-23, SEL-4, TL-15–27.
"""
import json
import re

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


# ── TL-15–17: Add Interior Wall ──────────────────────────────────────

class TestTL15AddIWWall:
    """TL-15: Add > Wall wizard creates IWN walls via POST /api/interior-walls."""

    def _fixed_formula(self, anchor="W1"):
        return {
            "type": "wall_rect",
            "anchor": anchor,
            "along": {"segment": ["W1", "W18"]},
            "thickness_dir": {"perp": {"segment": ["W1", "W18"]}},
            "thickness": {"const": "WALL_4IN"},
            "end_mode": "fixed",
            "length": 8.0,
        }

    def test_create_iw_fixed_formula(self, app_client):
        """POST /api/interior-walls with fixed formula returns 201 IW13."""
        resp = app_client.post("/api/interior-walls", json={"formula": self._fixed_formula()})
        assert resp.status_code == 201
        data = resp.get_json()
        assert data["name"] == "IW13"
        assert data["type"] == "wall"

    def test_create_iw_intersect_formula(self, app_client):
        """Intersect-mode formula is stored correctly."""
        formula = {
            "type": "wall_rect",
            "anchor": "W2",
            "along": {"segment": ["W2", "W5"]},
            "thickness_dir": {"perp": {"segment": ["W2", "W5"]}},
            "thickness": {"const": "WALL_4IN"},
            "end_mode": "intersect",
            "end_target": "inner_poly",
            "select": "nearest",
        }
        resp = app_client.post("/api/interior-walls", json={"formula": formula})
        assert resp.status_code == 201
        name = resp.get_json()["name"]
        rows = app_client.get(f"/api/formulas/{name}").get_json()
        pos_row = next(r for r in rows if r["param_name"] == "position")
        stored = json.loads(pos_row["formula_json"]) if isinstance(
            pos_row["formula_json"], str) else pos_row["formula_json"]
        assert stored["end_mode"] == "intersect"

    def test_iw_naming_sequence(self, app_client):
        """Creating 3 walls produces IW13, IW14, IW15."""
        names = []
        for _ in range(3):
            resp = app_client.post("/api/interior-walls", json={"formula": self._fixed_formula()})
            assert resp.status_code == 201
            names.append(resp.get_json()["name"])
        assert names == ["IW13", "IW14", "IW15"]

    def test_iw_formula_deps(self, app_client):
        """Named anchor W1 creates a formula_deps row for W1."""
        resp = app_client.post("/api/interior-walls", json={"formula": self._fixed_formula("W1")})
        name = resp.get_json()["name"]
        deps = app_client.get(f"/api/formulas/{name}/deps").get_json()
        dep_names = [d["dep_name"] for d in deps]
        assert "W1" in dep_names


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
# TL-17: Interior Wall Move
# =========================================================================

class TestTL17IWWallMove:
    """TL-17: Move tool works on user-created IW walls with literal anchor."""

    def _create_literal_wall(self, client):
        """Create IW13 with a literal [x, y] anchor."""
        formula = {
            "type": "wall_rect",
            "anchor": [-5.0, 0.0],
            "along": [1.0, 0.0],
            "thickness_dir": [0.0, 1.0],
            "thickness": 4.0 / 12.0,
            "end_mode": "fixed",
            "length": 5.0,
        }
        resp = client.post("/api/interior-walls", json={"formula": formula})
        assert resp.status_code == 201
        return resp.get_json()

    def test_move_iw_literal_anchor(self, app_client):
        """Move IW wall with literal anchor translates anchor by dx/dy."""
        rec = self._create_literal_wall(app_client)
        eid, name = rec["id"], rec["name"]
        resp = app_client.post(f"/api/elements/{eid}/move", json={"dx": 1.0, "dy": 0.5})
        assert resp.status_code == 200
        assert resp.get_json()["ok"] is True
        rows = app_client.get(f"/api/formulas/{name}").get_json()
        pos = next(r for r in rows if r["param_name"] == "position")
        f = json.loads(pos["formula_json"]) if isinstance(pos["formula_json"], str) else pos["formula_json"]
        assert abs(f["anchor"][0] - (-4.0)) < 1e-9
        assert abs(f["anchor"][1] - 0.5) < 1e-9

    def test_move_iw_literal_anchor_undo(self, app_client):
        """Undo after moving IW wall with literal anchor restores anchor."""
        rec = self._create_literal_wall(app_client)
        eid, name = rec["id"], rec["name"]
        app_client.post(f"/api/elements/{eid}/move", json={"dx": 2.0, "dy": 0.0})
        resp = app_client.post("/api/undo")
        assert resp.status_code == 200
        rows = app_client.get(f"/api/formulas/{name}").get_json()
        pos = next(r for r in rows if r["param_name"] == "position")
        f = json.loads(pos["formula_json"]) if isinstance(pos["formula_json"], str) else pos["formula_json"]
        assert abs(f["anchor"][0] - (-5.0)) < 1e-9
        assert abs(f["anchor"][1] - 0.0) < 1e-9

    def test_move_iw_named_anchor_rejected(self, app_client):
        """IW wall with named-point anchor → move → 400."""
        formula = {
            "type": "wall_rect",
            "anchor": "W1",
            "along": {"segment": ["W1", "W18"]},
            "thickness_dir": {"perp": {"segment": ["W1", "W18"]}},
            "thickness": {"const": "WALL_4IN"},
            "end_mode": "fixed",
            "length": 8.0,
        }
        resp = app_client.post("/api/interior-walls", json={"formula": formula})
        eid = resp.get_json()["id"]
        resp = app_client.post(f"/api/elements/{eid}/move", json={"dx": 1.0, "dy": 0.0})
        assert resp.status_code == 400


# =========================================================================
# TL-21: Add Opening Tool
# =========================================================================

class TestTL21AddOpening:
    """TL-21: Creating openings on walls via API."""

    def _create_user_iw_wall(self, client):
        """Create a user IW wall to host openings via POST /api/interior-walls."""
        formula = {
            "type": "wall_rect",
            "anchor": [-10.0, 0.0],
            "along": [1.0, 0.0],
            "thickness_dir": [0.0, 1.0],
            "thickness": 4.0 / 12,
            "end_mode": "fixed",
            "length": 5.0,
        }
        resp = client.post("/api/interior-walls", json={"formula": formula})
        assert resp.status_code == 201
        return resp.get_json()

    def _create_opening(self, client, wall_name, width=32.0/12, gap=0.5):
        """Create a user RO opening on a wall via POST /api/openings."""
        return client.post("/api/openings", json={
            "wall": wall_name,
            "gap": gap,
            "width": width,
        })

    def test_create_opening_on_wall(self, app_client):
        """POST opening creates an RO-named formula-based opening element."""
        wall = self._create_user_iw_wall(app_client)
        resp = self._create_opening(app_client, wall["name"])
        assert resp.status_code == 201
        data = resp.get_json()
        assert data["type"] == "opening"
        assert re.match(r"^RO\d+$", data["name"]), f"expected RO{{N}}, got {data['name']}"

    def test_opening_auto_names_sequentially(self, app_client):
        """Two openings on the same wall get consecutive RO numbers."""
        wall = self._create_user_iw_wall(app_client)
        r1 = self._create_opening(app_client, wall["name"], gap=0.5).get_json()
        r2 = self._create_opening(app_client, wall["name"], gap=2.0).get_json()
        n1 = int(r1["name"][2:])
        n2 = int(r2["name"][2:])
        assert n2 == n1 + 1

    def test_opening_host_wall_reference(self, app_client):
        """Created opening has host_wall in properties for cascade-delete."""
        wall = self._create_user_iw_wall(app_client)
        host = wall["name"]
        data = self._create_opening(app_client, host).get_json()
        props = json.loads(data["properties"]) if isinstance(data["properties"], str) else data["properties"]
        assert props["host_wall"] == host
        assert props["wall_name"] == host

    def test_opening_has_formula(self, app_client):
        """Created opening stores a wall_opening formula with correct element refs."""
        wall = self._create_user_iw_wall(app_client)
        host = wall["name"]
        data = self._create_opening(app_client, host, width=2.667, gap=1.0).get_json()
        name = data["name"]
        resp = app_client.get(f"/api/formulas/{name}")
        assert resp.status_code == 200
        formulas = resp.get_json()
        pos = next((f for f in formulas if f["param_name"] == "position"), None)
        assert pos is not None
        fj = pos["formula_json"]
        if isinstance(fj, str):
            import json as _json; fj = _json.loads(fj)
        assert fj["type"] == "wall_opening"
        assert fj["outer_start"]["element"] == host
        assert abs(fj["gap"] - 1.0) < 1e-9
        assert abs(fj["width"] - 2.667) < 1e-6

    def test_opening_appears_in_geometry(self, app_client):
        """Created RO opening appears in /api/geometry rough_openings list."""
        wall = self._create_user_iw_wall(app_client)
        host = wall["name"]
        data = self._create_opening(app_client, host, width=2.0, gap=1.0).get_json()
        name = data["name"]
        geom = app_client.get("/api/geometry").get_json()
        ro_names = [r["name"] for r in geom.get("rough_openings", [])]
        assert name in ro_names, f"{name} not in rough_openings: {ro_names}"

    def test_delete_host_wall_cascades_opening(self, app_client):
        """Deleting wall cascades to hosted openings."""
        wall = self._create_user_iw_wall(app_client)
        host = wall["name"]
        opening = self._create_opening(app_client, host).get_json()
        ro_name = opening["name"]
        # Verify opening exists
        elems = app_client.get("/api/elements").get_json()
        assert any(e["name"] == ro_name for e in elems)
        # Delete the wall
        resp = app_client.delete(f"/api/elements/{wall['id']}")
        assert resp.status_code == 200
        # Opening should be gone too
        elems = app_client.get("/api/elements").get_json()
        assert not any(e["name"] == ro_name for e in elems)


# =========================================================================
# SEL-15: Constant Dependency Highlighting
# =========================================================================

class TestSEL15ConstantDeps:
    """SEL-15: Reverse constant → element dependency mapping."""

    def test_reverse_iw_map_completeness(self):
        """Every IW formula with a position_constant has a derivable reverse mapping."""
        from app.evaluator import get_iw_formulas
        formulas = get_iw_formulas()
        # Build reverse map from formulas
        constant_to_iw = {}
        for iw, formula in formulas.items():
            cname = formula.get("position_constant")
            if cname:
                constant_to_iw.setdefault(cname, []).append(iw)
        # Every position_constant must appear in the reverse map
        for iw, formula in formulas.items():
            cname = formula.get("position_constant")
            if cname is None:
                continue
            assert cname in constant_to_iw, f"{cname} missing from reverse map"
            assert iw in constant_to_iw[cname], f"{iw} missing from reverse map[{cname}]"

    def test_specific_constant_mapping(self):
        """IW1_OFFSET_FROM_W9 maps to IW1 in formula metadata."""
        from app.evaluator import get_iw_formulas
        formulas = get_iw_formulas()
        assert formulas["IW1"].get("position_constant") == "IW1_OFFSET_FROM_W9"

    def test_hosted_openings_second_order(self):
        """IW9 hosts RO3 and RO7 (second-order deps)."""
        from app.elements import IW_HOSTED_OPENINGS
        hosted = IW_HOSTED_OPENINGS.get("IW9", [])
        assert "RO3" in hosted
        assert "RO7" in hosted
