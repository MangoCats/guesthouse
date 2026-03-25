"""Tests for Phase 5 — Outline Chain Editing.

Covers: ENG-11, API-16–19, outline undo actions, solver cross-validation.
"""
import math

import pytest

from app.database import (
    get_outline_chain, get_outline_chain_row, update_outline_segment,
    get_constants_dict,
)
from app.outline_solver import (
    ChainEntry, SolverResult, WalkResult,
    chain_offset, solve_closure,
    db_rows_to_chain, walk_chain,
    solve_for_constraint, flex_specs_from_chain_rows,
)

from tests.test_zapp_conftest import fresh_db, app_client  # noqa: F401


# ── Solver Unit Tests ─────────────────────────────────────────────────

class TestChainOffset:
    """Verify chain_offset for basic segments."""

    def test_empty_chain(self):
        """Empty chain returns (0, 0, start_brg)."""
        dE, dN, brg = chain_offset([], start_brg=0.5)
        assert dE == 0.0
        assert dN == 0.0
        assert brg == 0.5

    def test_line_north(self):
        """Single northward line: bearing 0, distance 10."""
        seg = ChainEntry("L", 10.0, None, None, None, 60, "P1")
        dE, dN, brg = chain_offset([seg], start_brg=0.0)
        assert abs(dE) < 1e-12
        assert abs(dN - 10.0) < 1e-12
        assert abs(brg) < 1e-12

    def test_line_east(self):
        """Single eastward line: bearing π/2, distance 5."""
        seg = ChainEntry("L", 5.0, None, None, None, 60, "P1")
        dE, dN, brg = chain_offset([seg], start_brg=math.pi / 2)
        assert abs(dE - 5.0) < 1e-12
        assert abs(dN) < 1e-12

    def test_cw_arc_90(self):
        """90° CW arc of radius 1 from bearing 0."""
        seg = ChainEntry("CW", None, 1.0, math.pi / 2, "C", 20, "P1")
        dE, dN, brg = chain_offset([seg], start_brg=0.0)
        # CW 90° from bearing 0: should end up 1 unit east, 1 unit north
        assert abs(dE - 1.0) < 1e-12
        assert abs(dN - 1.0) < 1e-12
        assert abs(brg - math.pi / 2) < 1e-12

    def test_ccw_arc_90(self):
        """90° CCW arc of radius 1 from bearing 0."""
        seg = ChainEntry("CCW", None, 1.0, math.pi / 2, "C", 20, "P1")
        dE, dN, brg = chain_offset([seg], start_brg=0.0)
        # CCW 90° from bearing 0: should end up 1 unit west, 1 unit north
        assert abs(dE - (-1.0)) < 1e-12
        assert abs(dN - 1.0) < 1e-12
        assert abs(brg - (-math.pi / 2)) < 1e-12


class TestSolverCrossValidation:
    """Verify app solver matches floorplan/geometry.py for default chain."""

    def test_solve_closure_matches_floorplan(self, fresh_db):
        """Default DB chain must produce bit-identical closure distances."""
        from floorplan.geometry import OUTLINE_CHAIN, CHAIN_POINT_NAMES
        import floorplan.constants as fc

        # Get default chain from DB
        chain_rows = get_outline_chain(fresh_db)
        chain = db_rows_to_chain(chain_rows)
        R_a1 = fc.CORNER_SW_R

        result = solve_closure(chain, R_a1)

        assert result.valid, f"Default chain should be valid, error={result.closure_error}"

        # Compare with floorplan's solved values
        # OUTLINE_CHAIN[0] = ("L", d_F2_F5)
        fp_d_F2_F5 = OUTLINE_CHAIN[0][1]
        # OUTLINE_CHAIN[16] = ("L", d_F18_F1)
        fp_d_F18_F1 = OUTLINE_CHAIN[16][1]

        d_F2_F5 = result.solved_values[0][1]
        d_F18_F1 = result.solved_values[len(chain) - 2][1]
        sweep_closure = result.solved_values[len(chain) - 1][1]
        assert abs(d_F2_F5 - fp_d_F2_F5) < 1e-9, \
            f"d_F2_F5: app={d_F2_F5} vs fp={fp_d_F2_F5}"
        assert abs(d_F18_F1 - fp_d_F18_F1) < 1e-9, \
            f"d_F18_F1: app={d_F18_F1} vs fp={fp_d_F18_F1}"
        # Default closure arc sweep should be π/2 (90°)
        assert abs(sweep_closure - math.pi / 2) < 1e-9, \
            f"sweep_closure: app={sweep_closure} vs expected={math.pi / 2}"

    def test_sweep_closure_adjusts_for_modified_sweep(self, fresh_db):
        """Changing an arc sweep changes the closure arc sweep to compensate."""
        chain_rows = get_outline_chain(fresh_db)
        chain = db_rows_to_chain(chain_rows)
        import floorplan.constants as fc
        R_a1 = fc.CORNER_SW_R

        # Default closure sweep
        result0 = solve_closure(chain, R_a1)
        assert result0.valid

        # Reduce F5→F6 arc (seq 1) sweep by 10°
        delta = math.radians(10)
        chain = list(chain)
        old_sweep = chain[1].sweep
        chain[1] = chain[1]._replace(sweep=old_sweep - delta)

        result1 = solve_closure(chain, R_a1)
        assert result1.valid, "Modified sweep should still close"
        # Closure sweep should increase by the same 10°
        n = len(chain)
        sweep0 = result0.solved_values[n - 1][1]
        sweep1 = result1.solved_values[n - 1][1]
        assert abs(sweep1 - (sweep0 + delta)) < 1e-9, \
            f"sweep_closure should increase by {delta:.4f}: " \
            f"got {sweep1:.6f}, expected {sweep0 + delta:.6f}"

    def test_d_F2_F5_positive(self, fresh_db):
        """Solved d_F2_F5 must be positive."""
        chain = db_rows_to_chain(get_outline_chain(fresh_db))
        import floorplan.constants as fc
        result = solve_closure(chain, fc.CORNER_SW_R)
        assert result.solved_values[0][1] > 0

    def test_d_F18_F1_positive(self, fresh_db):
        """Solved d_F18_F1 must be positive."""
        chain = db_rows_to_chain(get_outline_chain(fresh_db))
        import floorplan.constants as fc
        result = solve_closure(chain, fc.CORNER_SW_R)
        assert result.solved_values[len(chain) - 2][1] > 0

    def test_walk_chain_matches_floorplan_points(self, fresh_db):
        """All F-series points from app solver match floorplan/geometry.py."""
        from floorplan.geometry import walk_outline_chain, OUTLINE_CHAIN
        import floorplan.constants as fc

        # Get floorplan reference points
        fp_pts = walk_outline_chain()

        # Get app solver points
        chain_rows = get_outline_chain(fresh_db)
        chain = db_rows_to_chain(chain_rows)
        R_a1 = fc.CORNER_SW_R
        sr = solve_closure(chain, R_a1)
        chain = list(chain)
        for seq, (param, value) in sr.solved_values.items():
            chain[seq] = chain[seq]._replace(**{param: value})

        start_E = -18.5
        start_N = -13.5 + R_a1
        wr = walk_chain(chain, start_E, start_N)

        # Compare all F-series points
        for name in ["F1", "F2", "F5", "F6", "F7", "F8", "F9", "F10",
                      "F11", "F11a", "F11b", "F12", "F13", "F14", "F15",
                      "F16", "F17", "F18"]:
            fp = fp_pts[name]
            ap = wr.points[name]
            assert abs(fp[0] - ap[0]) < 1e-9, \
                f"{name} E: fp={fp[0]} vs app={ap[0]}"
            assert abs(fp[1] - ap[1]) < 1e-9, \
                f"{name} N: fp={fp[1]} vs app={ap[1]}"

    def test_solve_closure_modified_radius(self, fresh_db):
        """Changing one arc radius still produces a valid closure."""
        chain_rows = get_outline_chain(fresh_db)
        chain = db_rows_to_chain(chain_rows)

        # Modify the F13->F14 arc radius (seq 11) slightly
        chain = list(chain)
        old_r = chain[11].radius
        chain[11] = chain[11]._replace(radius=old_r + 0.1)

        import floorplan.constants as fc
        result = solve_closure(chain, fc.CORNER_SW_R)
        assert result.valid, "Modified radius should still close"
        assert result.solved_values[0][1] > 0
        assert result.solved_values[len(chain) - 2][1] > 0

    def test_solve_closure_impossible_params(self):
        """Impossible parameters produce valid=False."""
        # Create a tiny chain that can't possibly close
        chain = [
            ChainEntry("L", 1.0, None, None, None, 60, "F5"),
            ChainEntry("L", 1000.0, None, None, None, 60, "F18"),
            ChainEntry("L", 1.0, None, None, None, 60, "F1"),
            ChainEntry("CW", None, 0.001, math.pi / 2, "C1", 20, "F2"),
        ]
        result = solve_closure(chain, 0.001)
        # Either invalid or the distances are negative
        if result.valid:
            assert result.solved_values[0][1] < 0 or result.solved_values[2][1] < 0


# ── Engine Integration Tests ──────────────────────────────────────────

class TestEngineIntegration:
    """Verify engine uses DB chain correctly (ENG-11)."""

    def test_engine_with_chain_rows(self, fresh_db):
        """compute_geometry with chain_rows produces valid geometry."""
        from app.engine import compute_geometry
        constants = get_constants_dict(fresh_db)
        chain_rows = get_outline_chain(fresh_db)
        geom = compute_geometry(constants, "standard", chain_rows, db_path=fresh_db)

        assert "points" in geom
        assert "outline_segments" in geom
        assert "F1" in geom["points"]
        assert "F2" in geom["points"]
        assert len(geom["outline_segments"]) == 18

    def test_engine_default_chain_matches(self, fresh_db):
        """Default DB chain produces identical geometry to module-scope."""
        from app.engine import compute_geometry
        constants = get_constants_dict(fresh_db)

        # Without chain_rows (old path)
        geom_old = compute_geometry(constants, "standard", None, db_path=fresh_db)
        # With chain_rows (new path)
        chain_rows = get_outline_chain(fresh_db)
        geom_new = compute_geometry(constants, "standard", chain_rows, db_path=fresh_db)

        # Compare F-series points
        for name in ["F1", "F2", "F5", "F6", "F7", "F8", "F9", "F10",
                      "F11", "F11a", "F11b", "F12", "F13", "F14", "F15",
                      "F16", "F17", "F18"]:
            old_pt = geom_old["points"][name]
            new_pt = geom_new["points"][name]
            assert abs(old_pt[0] - new_pt[0]) < 1e-9, \
                f"{name} E mismatch: {old_pt[0]} vs {new_pt[0]}"
            assert abs(old_pt[1] - new_pt[1]) < 1e-9, \
                f"{name} N mismatch: {old_pt[1]} vs {new_pt[1]}"

    def test_engine_modified_chain_differs(self, fresh_db):
        """Modified chain produces different geometry."""
        from app.engine import compute_geometry
        constants = get_constants_dict(fresh_db)
        chain_rows = get_outline_chain(fresh_db)

        geom_default = compute_geometry(constants, "standard", chain_rows, db_path=fresh_db)

        # Modify F9->F10 distance (seq 5, a line — seq 4 is the F8->F9 arc)
        chain_rows_mod = get_outline_chain(fresh_db)
        for row in chain_rows_mod:
            if row["seq"] == 5:
                row["distance"] = (row["distance"] or 0) + 1.0  # add 1 foot
                break

        geom_modified = compute_geometry(constants, "standard", chain_rows_mod, db_path=fresh_db)

        # F10 should be at a different position
        f10_default = geom_default["points"]["F10"]
        f10_modified = geom_modified["points"]["F10"]
        dist = math.sqrt(
            (f10_default[0] - f10_modified[0])**2 +
            (f10_default[1] - f10_modified[1])**2
        )
        assert dist > 0.5, f"F10 should have moved significantly, dist={dist}"


# ── API Tests ─────────────────────────────────────────────────────────

class TestAPI16UpdateOutline:
    """PUT /api/outline/<seq> tests."""

    def test_update_arc_radius(self, app_client):
        """Change an arc radius, verify closure valid."""
        # Seq 11 is F13->F14 arc
        resp = app_client.put("/api/outline/11",
                              json={"dist_or_radius": 2.0})
        assert resp.status_code == 200
        data = resp.get_json()
        assert data["ok"] is True
        assert data["closure_valid"] is True
        # Pivot-aware solver assigns flex to the section containing the edit;
        # verify at least one solved value exists with a positive distance.
        solved = data["solved_values"]
        assert solved, "Expected at least one solved flex value"
        dist_values = [v["value"] for v in solved.values()
                       if v["param"] == "distance"]
        assert any(d > 0 for d in dist_values), \
            "At least one solved distance flex should be positive"

    def test_update_line_distance(self, app_client):
        """Change a line distance, verify closure valid."""
        # Seq 5 is F9->F10 line (seq 4 is the F8->F9 arc)
        resp = app_client.put("/api/outline/5",
                              json={"dist_or_radius": 16.0})
        assert resp.status_code == 200
        data = resp.get_json()
        assert data["ok"] is True
        assert data["closure_valid"] is True

    def test_update_arc_sweep(self, app_client):
        """Change an arc sweep angle, verify closure is valid and a sweep flex adjusts."""
        # Seq 1 is F5->F6 arc (CW)
        resp = app_client.put("/api/outline/1",
                              json={"sweep": math.pi / 2 + 0.01})
        assert resp.status_code == 200
        data = resp.get_json()
        assert data["ok"] is True
        assert data["closure_valid"] is True
        # Pivot-aware solver assigns sweep flex to the section containing seq 1;
        # verify at least one sweep-flex value was solved.
        solved = data["solved_values"]
        assert solved, "Expected at least one solved flex value"
        sweep_values = [v for v in solved.values() if v["param"] == "sweep"]
        assert sweep_values, "Expected at least one solved sweep flex"

    def test_reject_solved_distance_seq0(self, app_client):
        """Cannot directly edit solved distance at seq 0."""
        resp = app_client.put("/api/outline/0",
                              json={"dist_or_radius": 5.0})
        assert resp.status_code == 400
        assert "solved" in resp.get_json()["error"].lower()

    def test_reject_solved_distance_seq16(self, app_client):
        """Cannot directly edit solved distance at seq N-2."""
        # Get chain length
        chain_resp = app_client.get("/api/outline")
        chain = chain_resp.get_json()
        solved_seq = len(chain) - 2  # second to last
        resp = app_client.put(f"/api/outline/{solved_seq}",
                              json={"dist_or_radius": 5.0})
        assert resp.status_code == 400

    def test_reject_closure_arc_sweep(self, app_client):
        """Cannot directly edit the closure arc sweep (last segment)."""
        chain_resp = app_client.get("/api/outline")
        chain = chain_resp.get_json()
        closure_seq = len(chain) - 1
        resp = app_client.put(f"/api/outline/{closure_seq}",
                              json={"sweep": math.pi / 3})
        assert resp.status_code == 400
        assert "solved" in resp.get_json()["error"].lower()

    def test_sweep_change_updates_closure_arc(self, app_client):
        """Modifying one arc sweep adjusts the section sweep-flex arc in DB."""
        # Get initial chain
        chain0 = app_client.get("/api/outline").get_json()

        # Change F5->F6 arc (seq 1) sweep by -0.01 rad
        resp = app_client.put("/api/outline/1",
                              json={"sweep": chain0[1]["sweep"] - 0.01})
        assert resp.status_code == 200
        data = resp.get_json()
        assert data["ok"] is True
        assert data["closure_valid"] is True

        # Pivot-aware solver assigns sweep flex to Section A (seq 9 by default).
        # Seq 9 sweep should have increased by ~0.01 to compensate.
        chain1 = app_client.get("/api/outline").get_json()
        sweep9_before = chain0[9]["sweep"]
        sweep9_after = chain1[9]["sweep"]
        assert abs(sweep9_after - (sweep9_before + 0.01)) < 1e-9, \
            f"Section A sweep flex (seq 9): {sweep9_after} vs expected {sweep9_before + 0.01}"

    def test_not_found(self, app_client):
        """Non-existent seq returns 404."""
        resp = app_client.put("/api/outline/999",
                              json={"dist_or_radius": 5.0})
        assert resp.status_code == 404

    def test_geometry_reflects_change(self, app_client):
        """After outline edit, geometry endpoint returns changed data."""
        # Get original F10 position
        geom_before = app_client.get("/api/geometry").get_json()
        f10_before = geom_before["points"]["F10"]

        # Change F9->F10 distance (seq 5)
        resp = app_client.put("/api/outline/5",
                              json={"dist_or_radius": 16.0})
        assert resp.status_code == 200

        # Get new F10 position
        geom_after = app_client.get("/api/geometry").get_json()
        f10_after = geom_after["points"]["F10"]

        # Should be different
        dist = math.sqrt(
            (f10_before[0] - f10_after[0])**2 +
            (f10_before[1] - f10_after[1])**2
        )
        assert dist > 0.1


class TestAPI17ValidateOutline:
    """POST /api/outline/validate tests."""

    def test_validate_current_chain(self, app_client):
        """Validate current (default) chain returns valid=true."""
        resp = app_client.post("/api/outline/validate",
                               json={"changes": {}})
        assert resp.status_code == 200
        data = resp.get_json()
        assert data["valid"] is True
        assert data["closure_error"] < 1e-9

    def test_validate_with_change(self, app_client):
        """Validate with a reasonable change returns valid."""
        resp = app_client.post("/api/outline/validate",
                               json={"changes": {"5": {"dist_or_radius": 16.0}}})
        assert resp.status_code == 200
        data = resp.get_json()
        assert data["valid"] is True


class TestAPI18AddPoint:
    """POST /api/outline/add-point tests."""

    def test_add_point_line_split(self, app_client):
        """Splitting a line segment adds two rows (L + CW arc + L)."""
        chain_before = app_client.get("/api/outline").get_json()
        n_before = len(chain_before)

        resp = app_client.post("/api/outline/add-point",
                               json={
                                   "after_seq": 5,
                                   "end_name": "F9a",
                               })
        assert resp.status_code == 200
        data = resp.get_json()
        assert data["ok"] is True
        assert data["closure_valid"] is True
        assert len(data["chain"]) == n_before + 2
        assert data["added"] == ["F9a", "F9ab"]
        # Middle segment should be a CW arc with 4 ft radius
        mid = next(r for r in data["chain"] if r["end_name"] == "F9ab")
        assert mid["seg_type"] == "CW"
        assert abs(mid["radius"] - 4.0) < 1e-6

    def test_add_point_arc_split(self, app_client):
        """Splitting an arc segment adds two rows (arc + L + arc)."""
        chain_before = app_client.get("/api/outline").get_json()
        n_before = len(chain_before)

        # Seq 1 is F5->F6 arc
        resp = app_client.post("/api/outline/add-point",
                               json={
                                   "after_seq": 1,
                                   "end_name": "F5a",
                               })
        assert resp.status_code == 200
        data = resp.get_json()
        assert data["ok"] is True
        assert len(data["chain"]) == n_before + 2
        assert data["added"] == ["F5a", "F5ab"]
        # Middle segment should be a 1" line
        mid = next(r for r in data["chain"] if r["end_name"] == "F5ab")
        assert mid["seg_type"] == "L"
        assert abs(mid["distance"] - 1.0 / 12.0) < 1e-6

    def test_add_point_duplicate_name_rejected(self, app_client):
        """Duplicate point name is rejected."""
        resp = app_client.post("/api/outline/add-point",
                               json={
                                   "after_seq": 5,
                                   "end_name": "F10",  # already exists
                               })
        assert resp.status_code == 400
        assert "already exists" in resp.get_json()["error"]

    def test_add_point_missing_params(self, app_client):
        """Missing required params returns 400."""
        resp = app_client.post("/api/outline/add-point",
                               json={"after_seq": 5})
        assert resp.status_code == 400


class TestAPI19DeletePoint:
    """DELETE /api/outline/<seq> tests."""

    def test_delete_point(self, app_client):
        """Delete a mid-chain line segment."""
        # First add a point so we have something safe to remove
        app_client.post("/api/outline/add-point",
                        json={"after_seq": 5, "end_name": "F9a"})

        chain_before = app_client.get("/api/outline").get_json()
        n_before = len(chain_before)

        # The new point is at seq 6 (F9->F10 line split at seq 5)
        resp = app_client.delete("/api/outline/6")
        assert resp.status_code == 200
        data = resp.get_json()
        assert data["ok"] is True
        assert data["closure_valid"] is True

    def test_delete_closure_arc_rejected(self, app_client):
        """Cannot delete the closure arc (last segment)."""
        chain = app_client.get("/api/outline").get_json()
        last_seq = len(chain) - 1
        resp = app_client.delete(f"/api/outline/{last_seq}")
        assert resp.status_code == 400
        assert "closure" in resp.get_json()["error"].lower()

    def test_delete_invalid_seq(self, app_client):
        """Invalid seq returns 400."""
        resp = app_client.delete("/api/outline/999")
        assert resp.status_code == 400


# ── Undo Tests ────────────────────────────────────────────────────────

class TestOutlineUndo:
    """Verify undo/redo for outline edits."""

    def test_undo_outline_update(self, app_client):
        """Update then undo restores original chain."""
        chain_before = app_client.get("/api/outline").get_json()

        # Edit seq 5 distance (F9->F10 line)
        app_client.put("/api/outline/5",
                       json={"dist_or_radius": 16.0})

        # Undo
        resp = app_client.post("/api/undo")
        assert resp.status_code == 200
        data = resp.get_json()
        assert data["ok"] is True

        # Chain should be restored
        chain_after = app_client.get("/api/outline").get_json()
        for a, b in zip(chain_before, chain_after):
            if a["distance"] is not None:
                assert abs((a["distance"] or 0) - (b["distance"] or 0)) < 1e-9
            if a["radius"] is not None:
                assert abs((a["radius"] or 0) - (b["radius"] or 0)) < 1e-9

    def test_undo_outline_add_point(self, app_client):
        """Add point then undo restores original chain."""
        chain_before = app_client.get("/api/outline").get_json()
        n_before = len(chain_before)

        app_client.post("/api/outline/add-point",
                        json={"after_seq": 5, "end_name": "F9a"})

        # Verify chain grew by 2 (two new points added)
        chain_mid = app_client.get("/api/outline").get_json()
        assert len(chain_mid) == n_before + 2

        # Undo
        resp = app_client.post("/api/undo")
        assert resp.status_code == 200

        chain_after = app_client.get("/api/outline").get_json()
        assert len(chain_after) == n_before

    def test_redo_outline_update(self, app_client):
        """Undo then redo re-applies edit."""
        # Edit seq 5 (F9->F10 line)
        app_client.put("/api/outline/5",
                       json={"dist_or_radius": 16.0})

        # Get edited chain
        chain_edited = app_client.get("/api/outline").get_json()
        dist_edited = next(r for r in chain_edited if r["seq"] == 5)["distance"]

        # Undo
        app_client.post("/api/undo")

        # Redo
        resp = app_client.post("/api/redo")
        assert resp.status_code == 200

        chain_redo = app_client.get("/api/outline").get_json()
        dist_redo = next(r for r in chain_redo if r["seq"] == 5)["distance"]
        assert abs(dist_edited - dist_redo) < 1e-9

    def test_undo_outline_remove_point(self, app_client):
        """Remove then undo restores chain."""
        # Add a point first (split F9->F10 line at seq 5, new row at seq 6)
        app_client.post("/api/outline/add-point",
                        json={"after_seq": 5, "end_name": "F9a"})

        chain_with_extra = app_client.get("/api/outline").get_json()
        n_with_extra = len(chain_with_extra)

        # Delete it
        app_client.delete("/api/outline/6")
        assert len(app_client.get("/api/outline").get_json()) < n_with_extra

        # Undo the delete
        app_client.post("/api/undo")
        chain_restored = app_client.get("/api/outline").get_json()
        assert len(chain_restored) == n_with_extra
