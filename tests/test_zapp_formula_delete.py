"""Tests for formula-aware element deletion with dependency re-basing (TL-25).

Covers: DELETE /api/formulas/<name>/element, inlining of element references
in dependent formulas, undo/redo of formula deletion.
"""
import json
import pytest

from tests.test_zapp_conftest import fresh_db, app_client


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _create_wall(client, name):
    return client.post(
        "/api/elements",
        data=json.dumps({"type": "wall", "name": name, "properties": {}}),
        content_type="application/json",
    )


def _put_formula(client, elem_name, formula, param_name="position"):
    return client.put(
        f"/api/formulas/{elem_name}/{param_name}",
        data=json.dumps({"formula": formula}),
        content_type="application/json",
    )


# ---------------------------------------------------------------------------
# DELETE /api/formulas/<element_name>/element
# ---------------------------------------------------------------------------

class TestFormulaDelete:
    def test_delete_element_removes_formula(self, app_client):
        """Deleting an element removes its formula."""
        _create_wall(app_client, "DEL1")
        formula = {
            "type": "wall_rect",
            "anchor": [0, 0], "along": "east",
            "thickness_dir": "north", "thickness": 0.5,
            "end_mode": "fixed", "length": 5.0,
        }
        _put_formula(app_client, "DEL1", formula)

        resp = app_client.delete("/api/formulas/DEL1/element")
        assert resp.status_code == 200
        assert resp.get_json()["ok"] is True

        # Formula is gone
        data = app_client.get("/api/formulas/DEL1").get_json()
        assert len(data) == 0

    def test_delete_nonexistent_returns_404(self, app_client):
        resp = app_client.delete("/api/formulas/NOSUCH/element")
        assert resp.status_code == 404

    def test_delete_rebases_dependents(self, app_client):
        """Deleting an element inlines its position in dependent formulas."""
        _create_wall(app_client, "BASE")
        _create_wall(app_client, "DEP")

        base_formula = {
            "type": "wall_rect",
            "anchor": [1, 2], "along": "east",
            "thickness_dir": "north", "thickness": 0.5,
            "end_mode": "fixed", "length": 5.0,
        }
        _put_formula(app_client, "BASE", base_formula)

        dep_formula = {
            "type": "item_rect",
            "anchor": {"element": "BASE", "corner": "NE"},
            "along": "east", "across": "north",
            "width": 2.0, "depth": 1.0,
            "anchor_corner": "sw",
        }
        _put_formula(app_client, "DEP", dep_formula)

        # Verify DEP depends on BASE
        deps = app_client.get("/api/formulas/DEP/deps").get_json()
        dep_names = {d["dep_name"] for d in deps if d["dep_type"] == "element"}
        assert "BASE" in dep_names

        # Delete BASE
        resp = app_client.delete("/api/formulas/BASE/element")
        assert resp.status_code == 200
        result = resp.get_json()
        assert "DEP" in result["rebased"]

        # DEP's formula should no longer reference BASE
        dep_formulas = app_client.get("/api/formulas/DEP").get_json()
        assert len(dep_formulas) == 1
        updated = json.loads(dep_formulas[0]["formula_json"])
        # The anchor should now be a literal [E, N] instead of {"element": "BASE", ...}
        anchor = updated["anchor"]
        assert isinstance(anchor, list), f"Expected literal point, got {anchor}"
        assert len(anchor) == 2

        # DEP should no longer depend on BASE
        deps_after = app_client.get("/api/formulas/DEP/deps").get_json()
        dep_names_after = {d["dep_name"] for d in deps_after if d["dep_type"] == "element"}
        assert "BASE" not in dep_names_after

    def test_delete_undo_restores_formula(self, app_client):
        """Undoing a formula delete restores the original formula."""
        _create_wall(app_client, "UND1")
        formula = {
            "type": "wall_rect",
            "anchor": [3, 4], "along": "east",
            "thickness_dir": "north", "thickness": 0.5,
            "end_mode": "fixed", "length": 5.0,
        }
        _put_formula(app_client, "UND1", formula)

        # Delete
        app_client.delete("/api/formulas/UND1/element")
        assert len(app_client.get("/api/formulas/UND1").get_json()) == 0

        # Undo
        resp = app_client.post("/api/undo")
        assert resp.status_code == 200

        # Formula is back
        data = app_client.get("/api/formulas/UND1").get_json()
        assert len(data) == 1
        restored = json.loads(data[0]["formula_json"])
        assert restored["anchor"] == [3, 4]

    def test_delete_undo_restores_dependent_formulas(self, app_client):
        """Undoing a formula delete restores dependent formulas to pre-inlined state."""
        _create_wall(app_client, "UB")
        _create_wall(app_client, "UD")

        base_formula = {
            "type": "wall_rect",
            "anchor": [0, 0], "along": "east",
            "thickness_dir": "north", "thickness": 0.5,
            "end_mode": "fixed", "length": 5.0,
        }
        _put_formula(app_client, "UB", base_formula)

        dep_formula = {
            "type": "item_rect",
            "anchor": {"element": "UB", "corner": "SE"},
            "along": "east", "across": "north",
            "width": 2.0, "depth": 1.0,
            "anchor_corner": "sw",
        }
        _put_formula(app_client, "UD", dep_formula)

        # Delete UB (inlines into UD)
        app_client.delete("/api/formulas/UB/element")
        dep_after = json.loads(
            app_client.get("/api/formulas/UD").get_json()[0]["formula_json"]
        )
        assert isinstance(dep_after["anchor"], list)  # inlined

        # Undo
        app_client.post("/api/undo")

        # UD's formula should reference UB again
        dep_restored = json.loads(
            app_client.get("/api/formulas/UD").get_json()[0]["formula_json"]
        )
        assert dep_restored["anchor"] == {"element": "UB", "corner": "SE"}

    def test_delete_redo_re_inlines_dependents(self, app_client):
        """Redo after undo re-inlines dependent formulas and re-deletes."""
        _create_wall(app_client, "RB")
        _create_wall(app_client, "RD")

        base_formula = {
            "type": "wall_rect",
            "anchor": [0, 0], "along": "east",
            "thickness_dir": "north", "thickness": 0.5,
            "end_mode": "fixed", "length": 5.0,
        }
        _put_formula(app_client, "RB", base_formula)

        dep_formula = {
            "type": "item_rect",
            "anchor": {"element": "RB", "corner": "SE"},
            "along": "east", "across": "north",
            "width": 2.0, "depth": 1.0,
            "anchor_corner": "sw",
        }
        _put_formula(app_client, "RD", dep_formula)

        # Delete RB (inlines into RD)
        app_client.delete("/api/formulas/RB/element")

        # Undo
        app_client.post("/api/undo")
        # RD should reference RB again
        dep_restored = json.loads(
            app_client.get("/api/formulas/RD").get_json()[0]["formula_json"]
        )
        assert dep_restored["anchor"] == {"element": "RB", "corner": "SE"}

        # Redo
        resp = app_client.post("/api/redo")
        assert resp.status_code == 200

        # RB formula should be gone again
        assert len(app_client.get("/api/formulas/RB").get_json()) == 0

        # RD's formula should be inlined again (literal point, not reference)
        dep_redone = json.loads(
            app_client.get("/api/formulas/RD").get_json()[0]["formula_json"]
        )
        assert isinstance(dep_redone["anchor"], list), \
            f"Expected inlined point after redo, got {dep_redone['anchor']}"

    def test_seeded_element_delete(self, app_client):
        """Deleting a seeded element like IW1 works and returns rebased list."""
        # IW1 should have dependents
        resp = app_client.delete("/api/formulas/IW1/element")
        assert resp.status_code == 200
        data = resp.get_json()
        assert data["ok"] is True
        # IW1 formula should be gone
        assert len(app_client.get("/api/formulas/IW1").get_json()) == 0

    def test_delete_variant_specific_element(self, app_client):
        """Deleting a variant-specific element (e.g. loveseat) works."""
        # loveseat is seeded with variant="standard"
        resp = app_client.delete(
            "/api/formulas/loveseat/element?variant=standard"
        )
        assert resp.status_code == 200
        data = resp.get_json()
        assert data["ok"] is True
        # Formula should be gone
        formulas = app_client.get(
            "/api/formulas/loveseat?variant=standard"
        ).get_json()
        assert len(formulas) == 0


# ---------------------------------------------------------------------------
# Evaluator: inline_element_refs unit tests
# ---------------------------------------------------------------------------

class TestInlineElementRefs:
    def test_inline_element_corner(self):
        from app.evaluator import inline_element_refs
        elem_data = {"poly": [[1, 2], [3, 2], [3, 4], [1, 4]]}
        formula = {"anchor": {"element": "X", "corner": "SW"}}
        result = inline_element_refs(formula, "X", elem_data)
        assert result["anchor"] == [1, 2]

    def test_inline_element_centroid(self):
        from app.evaluator import inline_element_refs
        elem_data = {"poly": [[0, 0], [4, 0], [4, 4], [0, 4]]}
        formula = {"center": {"element_centroid": "X"}}
        result = inline_element_refs(formula, "X", elem_data)
        assert result["center"] == [2, 2]

    def test_inline_face_mid(self):
        from app.evaluator import inline_element_refs
        elem_data = {"poly": [[0, 0], [4, 0], [4, 2], [0, 2]]}
        formula = {"ref_point": {"face_mid": "X", "face": "south"}}
        result = inline_element_refs(formula, "X", elem_data)
        assert result["ref_point"] == [2, 0]

    def test_inline_face_along(self):
        from app.evaluator import inline_element_refs
        elem_data = {"poly": [[0, 0], [4, 0], [4, 2], [0, 2]]}
        formula = {"dir": {"face_along": "X", "face": "south"}}
        result = inline_element_refs(formula, "X", elem_data)
        assert result["dir"] == [1, 0]

    def test_no_match_unchanged(self):
        from app.evaluator import inline_element_refs
        elem_data = {"poly": [[0, 0], [1, 0], [1, 1], [0, 1]]}
        formula = {"anchor": {"element": "OTHER", "corner": "SW"}}
        result = inline_element_refs(formula, "X", elem_data)
        assert result["anchor"] == {"element": "OTHER", "corner": "SW"}

    def test_nested_inline(self):
        from app.evaluator import inline_element_refs
        elem_data = {"poly": [[0, 0], [4, 0], [4, 2], [0, 2]]}
        formula = {
            "type": "item_rect",
            "anchor": {
                "offset": {"element": "X", "corner": "NE"},
                "dir": "east",
                "dist": 1.0,
            },
        }
        result = inline_element_refs(formula, "X", elem_data)
        assert result["anchor"]["offset"] == [4, 2]
