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


def _put_formula(client, elem_name, formula, param_name="position",
                  variant=None):
    body = {"formula": formula}
    if variant is not None:
        body["variant"] = variant
    return client.put(
        f"/api/formulas/{elem_name}/{param_name}",
        data=json.dumps(body),
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

    def test_delete_element_without_formula(self, app_client):
        """Deleting an element that has a DB record but no formula succeeds."""
        # Create element with no formula (e.g. legacy seeded item)
        resp = app_client.post(
            "/api/elements",
            data=json.dumps({
                "type": "furniture", "name": "NOFMLA",
                "properties": {"width": 3.0, "depth": 2.0},
            }),
            content_type="application/json",
        )
        assert resp.status_code in (200, 201)

        # Verify element exists
        elems = app_client.get("/api/elements").get_json()
        assert any(e["name"] == "NOFMLA" for e in elems)

        # Delete should succeed (not return 404 "no formula found")
        resp = app_client.delete("/api/formulas/NOFMLA/element")
        assert resp.status_code == 200
        assert resp.get_json()["ok"] is True

        # Element is gone
        elems = app_client.get("/api/elements").get_json()
        assert not any(e["name"] == "NOFMLA" for e in elems)

        # Undo restores the element
        resp = app_client.post("/api/undo")
        assert resp.status_code == 200
        elems = app_client.get("/api/elements").get_json()
        assert any(e["name"] == "NOFMLA" for e in elems)

    def test_delete_preserves_other_variant_formulas(self, app_client):
        """Deleting an element in one variant preserves other variants' formulas.

        Regression: delete_element cascade-wiped all variant formulas.
        Deleting fridge in 'standard' should keep the minik formula intact.
        """
        # Create element with a formula for two variants
        app_client.post(
            "/api/elements",
            data=json.dumps({
                "type": "appliance", "name": "gadget",
                "properties": {"width": 2.0, "depth": 1.5,
                               "variants": ["standard", "minik"]},
            }),
            content_type="application/json",
        )
        formula_std = {
            "type": "item_rect",
            "anchor": [1, 2], "along": [1, 0], "across": [0, 1],
            "width": 2.0, "depth": 1.5, "anchor_corner": "center",
        }
        formula_mk = {
            "type": "item_rect",
            "anchor": [3, 4], "along": [1, 0], "across": [0, 1],
            "width": 1.5, "depth": 1.0, "anchor_corner": "center",
        }
        _put_formula(app_client, "gadget", formula_std, variant="standard")
        _put_formula(app_client, "gadget", formula_mk, variant="minik")

        # Delete in standard variant
        resp = app_client.delete(
            "/api/formulas/gadget/element?variant=standard")
        assert resp.status_code == 200

        # Element record should still exist (minik formula remains)
        elems = app_client.get("/api/elements").get_json()
        assert any(e["name"] == "gadget" for e in elems), \
            "element record should survive when other variant formulas remain"

        # Minik formula should still be there
        # (GET without variant returns only variant=NULL rows,
        #  so query with variant=minik to include it)
        formulas = app_client.get(
            "/api/formulas/gadget?variant=minik").get_json()
        assert len(formulas) > 0, "minik formula should survive"

    def test_delete_rebases_dependents(self, app_client):
        """Deleting an element structurally rebases dependent formulas."""
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

        # DEP's formula should no longer reference BASE — anchor is now a
        # structural expression (offset spec), not a literal point
        dep_formulas = app_client.get("/api/formulas/DEP").get_json()
        assert len(dep_formulas) == 1
        updated = json.loads(dep_formulas[0]["formula_json"])
        anchor = updated["anchor"]
        # wall_rect NE = offset(offset(anchor, length, along), thickness, thick_dir)
        assert isinstance(anchor, dict), f"Expected rebased spec, got {anchor}"
        assert "offset" in anchor

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

        # Delete UB (rebases UD)
        app_client.delete("/api/formulas/UB/element")
        dep_after = json.loads(
            app_client.get("/api/formulas/UD").get_json()[0]["formula_json"]
        )
        assert dep_after["anchor"] != {"element": "UB", "corner": "SE"}  # rebased

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

        # RD's formula should be rebased again (not referencing RB)
        dep_redone = json.loads(
            app_client.get("/api/formulas/RD").get_json()[0]["formula_json"]
        )
        assert dep_redone["anchor"] != {"element": "RB", "corner": "SE"}, \
            f"Expected rebased formula after redo, got {dep_redone['anchor']}"

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

    def test_delete_rebases_unrebaseable_to_four_corner(self, app_client):
        """Deleting dining_table converts chairs to four_corner literals."""
        # Verify chairs reference dining_table before deletion
        c1 = app_client.get("/api/formulas/dining_chair_1").get_json()
        c2 = app_client.get("/api/formulas/dining_chair_2").get_json()
        assert len(c1) > 0, "dining_chair_1 formula should exist"
        c1_f = json.loads(c1[0]["formula_json"])
        assert c1_f.get("table") == "dining_table"

        resp = app_client.delete(
            "/api/formulas/dining_table/element?variant=standard")
        assert resp.status_code == 200
        data = resp.get_json()
        assert data["ok"] is True
        assert "dining_chair_1" in data.get("rebased", [])
        assert "dining_chair_2" in data.get("rebased", [])

        # Chair formulas should still exist but be four_corner literals
        c1_after = app_client.get("/api/formulas/dining_chair_1").get_json()
        c2_after = app_client.get("/api/formulas/dining_chair_2").get_json()
        assert len(c1_after) > 0, "dining_chair_1 should still exist"
        assert len(c2_after) > 0, "dining_chair_2 should still exist"
        c1_new = json.loads(c1_after[0]["formula_json"])
        assert c1_new["type"] == "four_corner"
        assert isinstance(c1_new["sw"], list), "corners should be literals"
        assert len(c1_new["sw"]) == 2
        # Should not reference dining_table anymore
        assert "dining_table" not in json.dumps(c1_new)

    def test_four_corner_fallback_undo_restores_original(self, app_client):
        """Undoing a four_corner fallback restores the original formula."""
        app_client.delete(
            "/api/formulas/dining_table/element?variant=standard")
        app_client.post("/api/undo")

        # Table and chairs should all be restored with original formulas
        tbl = app_client.get(
            "/api/formulas/dining_table?variant=standard").get_json()
        assert len(tbl) > 0, "dining_table formula should be restored"
        c1 = app_client.get(
            "/api/formulas/dining_chair_1?variant=standard").get_json()
        assert len(c1) > 0, "dining_chair_1 should be restored"
        c1_f = json.loads(c1[0]["formula_json"])
        assert c1_f.get("table") == "dining_table", \
            "chair should reference dining_table again after undo"

    def test_four_corner_fallback_redo(self, app_client):
        """Redo after undo re-applies four_corner conversion."""
        app_client.delete(
            "/api/formulas/dining_table/element?variant=standard")
        app_client.post("/api/undo")
        app_client.post("/api/redo")

        # Table should be gone
        tbl = app_client.get(
            "/api/formulas/dining_table?variant=standard").get_json()
        assert len(tbl) == 0
        # Chairs should be four_corner
        c1 = app_client.get("/api/formulas/dining_chair_1").get_json()
        assert len(c1) > 0, "dining_chair_1 should still exist"
        c1_f = json.loads(c1[0]["formula_json"])
        assert c1_f["type"] == "four_corner"


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


# ---------------------------------------------------------------------------
# Evaluator: rebase_element_refs unit tests
# ---------------------------------------------------------------------------

class TestRebaseElementRefs:
    """Structural rebasing replaces element references with the deleted
    element's own formula sub-expressions instead of literal coordinates."""

    def test_rebase_four_corner_element_corner(self):
        """Corner ref → deleted element's corner sub-expression."""
        from app.evaluator import rebase_element_refs
        elem_formula = {
            "type": "four_corner",
            "sw": {"offset": "W2", "dist": 3.0, "dir": "east"},
            "se": {"offset": "W2", "dist": 5.0, "dir": "east"},
            "ne": [10, 10],
            "nw": [0, 10],
        }
        dep_formula = {"anchor": {"element": "X", "corner": "SW"}}
        result = rebase_element_refs(dep_formula, "X", elem_formula)
        # Should get the sw spec, not a literal coordinate
        assert result["anchor"] == {"offset": "W2", "dist": 3.0, "dir": "east"}

    def test_rebase_four_corner_named_corner(self):
        """Corner string names (SW/SE/NE/NW) all resolve correctly."""
        from app.evaluator import rebase_element_refs
        elem_formula = {
            "type": "four_corner",
            "sw": [0, 0], "se": [4, 0], "ne": [4, 2], "nw": [0, 2],
        }
        for corner_name, expected in [("SW", [0, 0]), ("SE", [4, 0]),
                                       ("NE", [4, 2]), ("NW", [0, 2])]:
            dep = {"pt": {"element": "X", "corner": corner_name}}
            result = rebase_element_refs(dep, "X", elem_formula)
            assert result["pt"] == expected

    def test_rebase_four_corner_integer_corner(self):
        """Corner integer indices (0-3) resolve correctly."""
        from app.evaluator import rebase_element_refs
        elem_formula = {
            "type": "four_corner",
            "sw": "A", "se": "B", "ne": "C", "nw": "D",
        }
        dep = {"pt": {"element": "X", "corner": 2}}
        result = rebase_element_refs(dep, "X", elem_formula)
        assert result["pt"] == "C"  # corner 2 = NE

    def test_rebase_face_mid(self):
        """face_mid → midpoint of two corner specs."""
        from app.evaluator import rebase_element_refs
        elem_formula = {
            "type": "four_corner",
            "sw": "P1", "se": "P2", "ne": "P3", "nw": "P4",
        }
        dep = {"ref": {"face_mid": "X", "face": "south"}}
        result = rebase_element_refs(dep, "X", elem_formula)
        assert result["ref"] == {"midpoint": ["P1", "P2"]}

    def test_rebase_face_along(self):
        """face_along → segment direction from corner specs."""
        from app.evaluator import rebase_element_refs
        elem_formula = {
            "type": "four_corner",
            "sw": "P1", "se": "P2", "ne": "P3", "nw": "P4",
        }
        dep = {"dir": {"face_along": "X", "face": "east"}}
        result = rebase_element_refs(dep, "X", elem_formula)
        assert result["dir"] == {"segment": ["P2", "P3"]}

    def test_rebase_face_perp(self):
        """face_perp → segment_perp from corner specs."""
        from app.evaluator import rebase_element_refs
        elem_formula = {
            "type": "four_corner",
            "sw": "P1", "se": "P2", "ne": "P3", "nw": "P4",
        }
        dep = {"dir": {"face_perp": "X", "face": "north"}}
        result = rebase_element_refs(dep, "X", elem_formula)
        assert result["dir"] == {"segment_perp": ["P3", "P4"]}

    def test_rebase_element_centroid(self):
        """element_centroid → nested midpoints of all 4 corners."""
        from app.evaluator import rebase_element_refs
        elem_formula = {
            "type": "four_corner",
            "sw": "A", "se": "B", "ne": "C", "nw": "D",
        }
        dep = {"center": {"element_centroid": "X"}}
        result = rebase_element_refs(dep, "X", elem_formula)
        assert result["center"] == {
            "midpoint": [{"midpoint": ["A", "B"]}, {"midpoint": ["C", "D"]}]
        }

    def test_rebase_wall_rect_corners(self):
        """wall_rect formula derives structural corner specs."""
        from app.evaluator import rebase_element_refs
        elem_formula = {
            "type": "wall_rect",
            "anchor": [1, 2], "along": "east",
            "thickness_dir": "north", "thickness": 0.5,
            "end_mode": "fixed", "length": 5.0,
        }
        # NE = offset(offset(anchor, length, along), thickness, thick_dir)
        dep = {"pt": {"element": "X", "corner": "NE"}}
        result = rebase_element_refs(dep, "X", elem_formula)
        ne = result["pt"]
        assert isinstance(ne, dict) and "offset" in ne

    def test_rebase_item_rect_corners(self):
        """item_rect formula derives structural corner specs."""
        from app.evaluator import rebase_element_refs
        elem_formula = {
            "type": "item_rect",
            "anchor": {"element": "PARENT", "corner": "SW"},
            "along": "east", "across": "north",
            "width": 2.0, "depth": 1.0,
            "anchor_corner": "sw",
        }
        # SW = anchor (since anchor_corner is "sw")
        dep = {"pt": {"element": "X", "corner": "SW"}}
        result = rebase_element_refs(dep, "X", elem_formula)
        # Should be the anchor spec from elem_formula, referencing PARENT
        assert result["pt"] == {"element": "PARENT", "corner": "SW"}

    def test_rebase_no_match_unchanged(self):
        """References to other elements are left unchanged."""
        from app.evaluator import rebase_element_refs
        elem_formula = {
            "type": "four_corner",
            "sw": "A", "se": "B", "ne": "C", "nw": "D",
        }
        dep = {"pt": {"element": "OTHER", "corner": "SW"}}
        result = rebase_element_refs(dep, "X", elem_formula)
        assert result["pt"] == {"element": "OTHER", "corner": "SW"}

    def test_rebase_deep_copy(self):
        """Rebased specs are deep copies — no shared mutation."""
        from app.evaluator import rebase_element_refs
        sw_spec = {"offset": "W2", "dist": 3.0, "dir": "east"}
        elem_formula = {
            "type": "four_corner",
            "sw": sw_spec, "se": "B", "ne": "C", "nw": "D",
        }
        dep = {"a": {"element": "X", "corner": "SW"},
               "b": {"element": "X", "corner": "SW"}}
        result = rebase_element_refs(dep, "X", elem_formula)
        # Modify one — the other must not change
        result["a"]["dist"] = 999
        assert result["b"]["dist"] == 3.0

    def test_rebase_falls_back_to_literal_for_unsupported_type(self):
        """Unsupported formula types fall back to literal inlining."""
        from app.evaluator import rebase_element_refs
        elem_formula = {"type": "item_circle", "center": [5, 5], "radius": 1.0}
        elem_data = {"poly": [[4, 5], [5, 4], [6, 5], [5, 6]]}  # 4 approx pts
        dep = {"pt": {"element": "X", "corner": "SW"}}
        result = rebase_element_refs(dep, "X", elem_formula, elem_data)
        # Should fall back to literal inlining from elem_data
        assert result["pt"] == [4, 5]

    def test_rebase_preserves_geometry(self):
        """Rebasing produces identical geometry to the original."""
        from app.evaluator import FormulaEvaluator, rebase_element_refs
        import json

        base_formula = {
            "type": "four_corner",
            "sw": [1, 2], "se": [5, 2], "ne": [5, 4], "nw": [1, 4],
        }
        dep_formula = {
            "type": "item_rect",
            "anchor": {"element": "BASE", "corner": "NE"},
            "along": "east", "across": "north",
            "width": 2.0, "depth": 1.0,
            "anchor_corner": "sw",
        }

        # Evaluate with both elements present
        ev1 = FormulaEvaluator({}, {}, [], {})
        ev1.add_formula("BASE", "position", base_formula)
        ev1.add_formula("DEP", "position", dep_formula)
        ev1.topo_sort()
        ev1.evaluate_all()
        poly_before = ev1.elements["DEP"]["poly"]

        # Rebase DEP and evaluate without BASE
        rebased = rebase_element_refs(dep_formula, "BASE", base_formula)
        ev2 = FormulaEvaluator({}, {}, [], {})
        ev2.add_formula("DEP", "position", rebased)
        ev2.topo_sort()
        ev2.evaluate_all()
        poly_after = ev2.elements["DEP"]["poly"]

        for i, (pb, pa) in enumerate(zip(poly_before, poly_after)):
            assert abs(pb[0] - pa[0]) < 1e-9, \
                f"Corner {i} E mismatch: {pb[0]} vs {pa[0]}"
            assert abs(pb[1] - pa[1]) < 1e-9, \
                f"Corner {i} N mismatch: {pb[1]} vs {pa[1]}"
