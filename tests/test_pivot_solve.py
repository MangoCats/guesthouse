"""Regression tests for anchor-relative pivot solve (Phase altF).

Verifies that:
- anchor/pivot positions are stored and retrieved correctly
- default anchor/pivot is seeded at DB init
- Section A points remain fixed when a Section B segment is edited
- Section B points move correctly after a Section B edit
- The global rotated walk and sub-chain B walk use consistent bearings

Chain structure (fresh DB, 18 segs):
  seq0..seq17 ending at F5,F6,F7,F8,F9,F10,F11,F11a,F11b,F12,F13,F14,F15,F16,F17,F18,F1,F2

Test setup uses anchor=F7 (seq2), pivot=F14 (seq11):
  Section A = seqs [3..11] → points F8..F14
  Section B = seqs [12..17,0,1,2] → points F15..F18,F1,F2,F5,F6,F7
"""
import math
import pytest

from app.database import (
    get_outline_chain, get_outline_anchor_pivot, get_outline_anchor_pos,
    set_outline_anchor_pivot, clear_outline_pivot,
    get_constant_value,
)
from app.outline_solver import (
    db_rows_to_chain, walk_chain, point_name_to_seq, section_seqs,
)

from tests.test_zapp_conftest import fresh_db, app_client  # noqa: F401

# Anchor/pivot point names valid in the standard 18-seg fresh DB chain
_ANCHOR = "F7"    # end of seq 2
_PIVOT  = "F14"   # end of seq 11
# A Section B arc to edit: seq 17 is CW arc ending at F2
_EDIT_SEQ = 17
_EDIT_DELTA_RAD = math.radians(5.0)

# Section A points that must not move
_SECTION_A_POINTS = ["F8", "F9", "F10", "F11", "F11a", "F11b",
                      "F12", "F13", "F14"]
# Section B points (should move after the edit)
_SECTION_B_POINT = "F2"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _rotated_walk(chain_rows, db_path):
    """Walk chain from stored anchor position (rotated), return WalkResult."""
    chain = db_rows_to_chain(chain_rows)
    n = len(chain)
    pos = get_outline_anchor_pos(db_path)
    anchor_name, _ = get_outline_anchor_pivot(db_path)
    anchor_pt_seq = point_name_to_seq(chain, anchor_name)
    a_start = (anchor_pt_seq + 1) % n
    E, N, brg = pos
    rotated = [chain[(a_start + i) % n] for i in range(n)]
    return walk_chain(rotated, E, N, brg)


# ---------------------------------------------------------------------------
# Test 1: anchor position stored and retrieved
# ---------------------------------------------------------------------------

class TestAnchorPosStorage:
    def test_stored_and_retrieved(self, fresh_db):
        """set_outline_anchor_pivot stores E/N/brg; get_outline_anchor_pos retrieves."""
        set_outline_anchor_pivot("F9", "F12", 3.14, -7.5, 0.25, fresh_db)
        anchor, pivot = get_outline_anchor_pivot(fresh_db)
        pos = get_outline_anchor_pos(fresh_db)

        assert anchor == "F9"
        assert pivot == "F12"
        assert pos is not None
        E, N, brg = pos
        assert abs(E - 3.14) < 1e-12
        assert abs(N - (-7.5)) < 1e-12
        assert abs(brg - 0.25) < 1e-12

    def test_returns_none_when_not_set(self, fresh_db):
        """get_outline_anchor_pos returns None if coords keys are missing."""
        from app.database import get_db
        with get_db(fresh_db) as conn:
            conn.execute("DELETE FROM config WHERE key IN "
                         "('outline_anchor_E', 'outline_anchor_N', 'outline_anchor_brg')")
        assert get_outline_anchor_pos(fresh_db) is None


# ---------------------------------------------------------------------------
# Test 2: defaults seeded at init
# ---------------------------------------------------------------------------

class TestDefaultAnchorSeed:
    def test_defaults_seeded_at_init(self, fresh_db):
        """Fresh DB has anchor, pivot, and anchor coords set."""
        anchor, pivot = get_outline_anchor_pivot(fresh_db)
        pos = get_outline_anchor_pos(fresh_db)

        assert anchor is not None, "anchor should be seeded"
        assert pivot is not None, "pivot should be seeded"
        assert pos is not None, "anchor coords should be seeded"

        E, N, brg = pos
        R_a1 = float(get_constant_value("CORNER_SW_R", fresh_db) or 1.0)
        expected_E = float(get_constant_value("F2_EASTING", fresh_db) or -18.5)
        expected_N = float(get_constant_value("F2_NORTHING", fresh_db) or -13.5) + R_a1

        assert abs(E - expected_E) < 1e-6
        assert abs(N - expected_N) < 1e-6
        assert abs(brg) < 1e-9  # bearing north = 0

    def test_default_anchor_is_chain_wrap_point(self, fresh_db):
        """Default anchor = chain[-1].end_name (chain wrap point)."""
        chain_rows = get_outline_chain(fresh_db)
        anchor, _ = get_outline_anchor_pivot(fresh_db)
        assert anchor == chain_rows[-1]["end_name"]

    def test_default_pivot_is_midpoint(self, fresh_db):
        """Default pivot = chain[n//2].end_name."""
        chain_rows = get_outline_chain(fresh_db)
        n = len(chain_rows)
        _, pivot = get_outline_anchor_pivot(fresh_db)
        assert pivot == chain_rows[n // 2]["end_name"]


# ---------------------------------------------------------------------------
# Test 3: clear_outline_pivot re-seeds defaults
# ---------------------------------------------------------------------------

class TestClearPivot:
    def test_clear_reseeds_defaults(self, fresh_db):
        """clear_outline_pivot restores default anchor/pivot/coords."""
        set_outline_anchor_pivot("F9", "F12", 3.14, -7.5, 0.25, fresh_db)
        clear_outline_pivot(fresh_db)

        anchor, pivot = get_outline_anchor_pivot(fresh_db)
        pos = get_outline_anchor_pos(fresh_db)

        assert anchor is not None
        assert pivot is not None
        assert pos is not None

        chain_rows = get_outline_chain(fresh_db)
        assert anchor == chain_rows[-1]["end_name"]

    def test_clear_anchor_pos_not_none(self, fresh_db):
        """After clear, anchor_pos is not None (defaults restored)."""
        set_outline_anchor_pivot("F9", "F12", 1.0, 2.0, 0.5, fresh_db)
        clear_outline_pivot(fresh_db)
        assert get_outline_anchor_pos(fresh_db) is not None


# ---------------------------------------------------------------------------
# Test 4+5: Section A fixed, Section B moves after Section B edit
# ---------------------------------------------------------------------------

class TestPivotSolveConsistency:

    def _set_anchor_pivot(self, client):
        r = client.put("/api/outline/pivot",
                       json={"anchor": _ANCHOR, "pivot": _PIVOT},
                       content_type="application/json")
        assert r.status_code == 200, f"PUT pivot failed: {r.data}"

    def _edit_sec_b_seg(self, client, fresh_db, delta_rad=_EDIT_DELTA_RAD):
        chain_rows = get_outline_chain(fresh_db)
        seg = next(r for r in chain_rows if r["seq"] == _EDIT_SEQ)
        new_sweep = seg["sweep"] + delta_rad
        r = client.put(f"/api/outline/{_EDIT_SEQ}",
                       json={"sweep": new_sweep},
                       content_type="application/json")
        assert r.status_code == 200, f"PUT seg{_EDIT_SEQ} failed: {r.data}"

    def test_section_a_points_fixed_after_section_b_edit(self, app_client, fresh_db):
        """Section A points (F8..F14) must not move when a Section B arc is edited."""
        self._set_anchor_pivot(app_client)

        # Record Section A point positions before the edit
        wr_before = _rotated_walk(get_outline_chain(fresh_db), fresh_db)
        pts_before = {name: wr_before.points[name]
                      for name in _SECTION_A_POINTS
                      if name in wr_before.points}
        assert pts_before, "Section A points should exist in walk result"

        # Edit a Section B arc
        self._edit_sec_b_seg(app_client, fresh_db)

        # Get positions after the edit
        wr_after = _rotated_walk(get_outline_chain(fresh_db), fresh_db)

        for name, (E_before, N_before) in pts_before.items():
            pt_after = wr_after.points.get(name)
            assert pt_after is not None, f"{name} not found after edit"
            assert abs(pt_after[0] - E_before) < 1e-6, \
                f"{name} E shifted by {pt_after[0] - E_before:.9f} after Section B edit"
            assert abs(pt_after[1] - N_before) < 1e-6, \
                f"{name} N shifted by {pt_after[1] - N_before:.9f} after Section B edit"

    def test_anchor_point_fixed_after_section_b_edit(self, app_client, fresh_db):
        """Anchor point (F7) must not move when a Section B arc is edited."""
        self._set_anchor_pivot(app_client)

        wr_before = _rotated_walk(get_outline_chain(fresh_db), fresh_db)
        f7_before = wr_before.points.get(_ANCHOR)
        assert f7_before is not None

        self._edit_sec_b_seg(app_client, fresh_db)

        wr_after = _rotated_walk(get_outline_chain(fresh_db), fresh_db)
        f7_after = wr_after.points.get(_ANCHOR)
        assert f7_after is not None

        assert abs(f7_after[0] - f7_before[0]) < 1e-6, \
            f"{_ANCHOR} E drifted by {f7_after[0] - f7_before[0]:.9f}"
        assert abs(f7_after[1] - f7_before[1]) < 1e-6, \
            f"{_ANCHOR} N drifted by {f7_after[1] - f7_before[1]:.9f}"

    def test_section_b_point_moves_after_edit(self, app_client, fresh_db):
        """A Section B point (F2) should change position after a Section B arc edit."""
        self._set_anchor_pivot(app_client)

        wr_before = _rotated_walk(get_outline_chain(fresh_db), fresh_db)
        f2_before = wr_before.points.get(_SECTION_B_POINT)
        assert f2_before is not None

        # Use a larger delta to ensure a detectable position change
        self._edit_sec_b_seg(app_client, fresh_db, delta_rad=math.radians(10.0))

        wr_after = _rotated_walk(get_outline_chain(fresh_db), fresh_db)
        f2_after = wr_after.points.get(_SECTION_B_POINT)
        assert f2_after is not None

        moved = (abs(f2_after[0] - f2_before[0]) > 1e-4 or
                 abs(f2_after[1] - f2_before[1]) > 1e-4)
        assert moved, f"{_SECTION_B_POINT} (Section B) should move after a Section B arc edit"

    def test_bearing_consistency_after_edit(self, app_client, fresh_db):
        """After a Section B edit, rotated-walk and sub-chain B walk agree at each point."""
        self._set_anchor_pivot(app_client)
        self._edit_sec_b_seg(app_client, fresh_db)

        chain_rows = get_outline_chain(fresh_db)
        chain = db_rows_to_chain(chain_rows)
        n = len(chain)
        pos = get_outline_anchor_pos(fresh_db)
        anchor_name, pivot_name = get_outline_anchor_pivot(fresh_db)
        anchor_pt_seq = point_name_to_seq(chain, anchor_name)
        pivot_pt_seq = point_name_to_seq(chain, pivot_name)
        a_start = (anchor_pt_seq + 1) % n
        p_start = (pivot_pt_seq + 1) % n

        E, N, brg = pos
        rotated = [chain[(a_start + i) % n] for i in range(n)]
        wr_global = walk_chain(rotated, E, N, brg)

        # Compute bearing at pivot in the rotated walk
        pivot_rotated_idx = (p_start - a_start) % n
        cur_brg = brg
        for i, seg in enumerate(rotated):
            if i == pivot_rotated_idx:
                break
            if seg.seg_type == "CW":
                cur_brg += seg.sweep
            elif seg.seg_type == "CCW":
                cur_brg -= seg.sweep
        brg_at_pivot = cur_brg

        # Walk sub-chain B from pivot
        b_seqs = []
        i = p_start
        while i != a_start:
            b_seqs.append(i)
            i = (i + 1) % n
        subchain = [chain[s] for s in b_seqs]
        pivot_pos = wr_global.points.get(pivot_name)
        assert pivot_pos is not None

        sub_wr = walk_chain(subchain, pivot_pos[0], pivot_pos[1], brg_at_pivot)

        # All Section B end-points should match between global and sub-chain walks
        for seg in subchain:
            name = seg.end_name
            g_pos = wr_global.points.get(name)
            s_pos = sub_wr.points.get(name)
            if g_pos is not None and s_pos is not None:
                assert abs(g_pos[0] - s_pos[0]) < 1e-6, \
                    f"{name} E: global={g_pos[0]:.9f}, sub={s_pos[0]:.9f}"
                assert abs(g_pos[1] - s_pos[1]) < 1e-6, \
                    f"{name} N: global={g_pos[1]:.9f}, sub={s_pos[1]:.9f}"
