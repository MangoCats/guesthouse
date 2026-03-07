"""Tests for module-state isolation: compute_geometry does NOT mutate floorplan modules.

Phase 12g removed patch_constants + importlib.reload from compute_geometry.
These tests verify that floorplan module state is never modified by the app.
"""
import pytest
from app.database import init_db, get_constants_dict, update_constant
from app.engine import compute_geometry
from tests.test_zapp_conftest import (
    _FLOORPLAN_MODULES, _snapshot_modules, _restore_modules,
)


class TestModuleIsolation:
    """Verify that compute_geometry does not mutate floorplan modules."""

    def test_snapshot_captures_bed_width_in_constants(self):
        import floorplan.constants as c
        snap = _snapshot_modules()
        assert "BED_WIDTH" in snap["floorplan.constants"]
        assert snap["floorplan.constants"]["BED_WIDTH"] == c.BED_WIDTH

    def test_snapshot_captures_bed_width_in_layout(self):
        import floorplan.layout as l
        snap = _snapshot_modules()
        assert "BED_WIDTH" in snap["floorplan.layout"]
        assert snap["floorplan.layout"]["BED_WIDTH"] == l.BED_WIDTH

    def test_snapshot_captures_closure_distances(self):
        import floorplan.geometry as g
        snap = _snapshot_modules()
        assert "_d_F2_F5" in snap["floorplan.geometry"]
        assert snap["floorplan.geometry"]["_d_F2_F5"] == g._d_F2_F5
        assert "_d_F18_F1" in snap["floorplan.geometry"]

    def test_constants_and_layout_agree(self):
        """floorplan.constants.BED_WIDTH == floorplan.layout.BED_WIDTH."""
        import floorplan.constants as c
        import floorplan.layout as l
        assert c.BED_WIDTH == l.BED_WIDTH

    def test_compute_geometry_does_not_patch_layout(self, tmp_path):
        """compute_geometry must NOT modify floorplan.layout module attrs."""
        import floorplan.layout as l
        orig = l.BED_WIDTH

        db = str(tmp_path / "test.db")
        init_db(db)
        update_constant("BED_WIDTH", orig + 2.0, db)
        compute_geometry(get_constants_dict(db), db_path=db)

        # layout must be unchanged — formulas handle constant propagation
        assert abs(l.BED_WIDTH - orig) < 1e-9

    def test_modules_unchanged_after_compute(self, tmp_path):
        """Every numeric attr in every module unchanged after compute_geometry."""
        snap_before = _snapshot_modules()

        db = str(tmp_path / "test.db")
        init_db(db)
        update_constant("BED_WIDTH", 5.0, db)
        compute_geometry(get_constants_dict(db), db_path=db)

        snap_after = _snapshot_modules()

        for mod_name in _FLOORPLAN_MODULES:
            for attr, orig_val in snap_before[mod_name].items():
                current_val = snap_after[mod_name].get(attr)
                assert current_val is not None, (
                    f"{mod_name}.{attr} missing after compute_geometry"
                )
                assert abs(current_val - orig_val) < 1e-12, (
                    f"{mod_name}.{attr}: expected={orig_val}, got={current_val}"
                )
