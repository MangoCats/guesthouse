"""Tests for module-state save/restore used by app test fixtures.

Verifies that after compute_geometry patches module-level constants,
the save/restore mechanism fully restores the state of all four
floorplan modules so other tests are unaffected.
"""
import pytest
from app.database import init_db, get_constants_dict, update_constant
from app.engine import compute_geometry
from tests.test_zapp_conftest import (
    _FLOORPLAN_MODULES, _snapshot_modules, _restore_modules,
)


class TestModuleIsolation:
    """Verify that save/restore round-trips perfectly."""

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

    def test_compute_geometry_patches_layout(self, tmp_path):
        """compute_geometry reloads layout, changing its BED_WIDTH."""
        import floorplan.layout as l
        snap = _snapshot_modules()
        orig = l.BED_WIDTH

        db = str(tmp_path / "test.db")
        init_db(db)
        update_constant("BED_WIDTH", orig + 2.0, db)
        compute_geometry(get_constants_dict(db), db_path=db)

        # layout now has the patched value
        assert abs(l.BED_WIDTH - (orig + 2.0)) < 1e-9

        # Restore ALL modules
        _restore_modules(snap)
        assert abs(l.BED_WIDTH - orig) < 1e-9

    def test_full_round_trip(self, tmp_path):
        """Every numeric attr in every module restored after compute_geometry."""
        snap_before = _snapshot_modules()

        db = str(tmp_path / "test.db")
        init_db(db)
        update_constant("BED_WIDTH", 5.0, db)
        compute_geometry(get_constants_dict(db), db_path=db)

        _restore_modules(snap_before)
        snap_after = _snapshot_modules()

        for mod_name in _FLOORPLAN_MODULES:
            for attr, orig_val in snap_before[mod_name].items():
                restored_val = snap_after[mod_name].get(attr)
                assert restored_val is not None, (
                    f"{mod_name}.{attr} missing after restore"
                )
                assert abs(restored_val - orig_val) < 1e-12, (
                    f"{mod_name}.{attr}: orig={orig_val}, restored={restored_val}"
                )

    def test_class_identity_stable_after_restore(self, tmp_path):
        """Restore does not reload modules — class objects stay the same."""
        import floorplan.layout as l
        snap = _snapshot_modules()

        db = str(tmp_path / "test.db")
        init_db(db)
        compute_geometry(get_constants_dict(db), db_path=db)
        reloaded_class_id = id(l.InteriorLayout)

        _restore_modules(snap)
        assert id(l.InteriorLayout) == reloaded_class_id
