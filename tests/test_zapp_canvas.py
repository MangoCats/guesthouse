"""Tests for Phase 6 — Enhanced Canvas Rendering (door arcs, clearance zones)."""
import math
import sys
import os

_PROJECT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _PROJECT not in sys.path:
    sys.path.insert(0, _PROJECT)

from tests.test_zapp_conftest import fresh_db, app_client, geometry  # noqa: F401
from app.database import (
    get_constants_dict, get_outline_chain, get_all_doors,
    create_door, delete_door,
)
from app.engine import compute_geometry


def _geom_with_doors(fresh_db):
    """Helper: compute geometry with door data."""
    constants = get_constants_dict(fresh_db)
    chain_rows = get_outline_chain(fresh_db)
    doors = get_all_doors(fresh_db)
    return compute_geometry(constants, "standard", chain_rows, doors_data=doors)


class TestDoorArcs:
    """Door arc geometry in engine output."""

    def test_door_arcs_in_geometry(self, fresh_db):
        g = _geom_with_doors(fresh_db)
        assert "door_arcs" in g
        assert len(g["door_arcs"]) == 9  # O3, O6, RO1-RO7

    def test_door_arc_ro1_single(self, fresh_db):
        g = _geom_with_doors(fresh_db)
        ro1 = next(a for a in g["door_arcs"] if a["opening_name"] == "RO1")
        assert ro1["door_type"] == "single"
        assert len(ro1["leaves"]) == 1
        leaf = ro1["leaves"][0]
        assert len(leaf["hinge"]) == 2
        assert len(leaf["tip"]) == 2
        assert len(leaf["arc_pts"]) == 21

    def test_door_arc_ro7_double(self, fresh_db):
        g = _geom_with_doors(fresh_db)
        ro7 = next(a for a in g["door_arcs"] if a["opening_name"] == "RO7")
        assert ro7["door_type"] == "double"
        assert len(ro7["leaves"]) == 2
        # Both leaves should have 21 arc points
        for leaf in ro7["leaves"]:
            assert len(leaf["arc_pts"]) == 21

    def test_door_arc_ro6_double(self, fresh_db):
        g = _geom_with_doors(fresh_db)
        ro6 = next(a for a in g["door_arcs"] if a["opening_name"] == "RO6")
        assert ro6["door_type"] == "double"
        assert len(ro6["leaves"]) == 2

    def test_door_arc_o3(self, fresh_db):
        g = _geom_with_doors(fresh_db)
        o3 = next(a for a in g["door_arcs"] if a["opening_name"] == "O3")
        assert o3["door_type"] == "single"
        assert len(o3["leaves"]) == 1

    def test_door_arc_o6(self, fresh_db):
        g = _geom_with_doors(fresh_db)
        o6 = next(a for a in g["door_arcs"] if a["opening_name"] == "O6")
        assert o6["door_type"] == "single"
        assert len(o6["leaves"]) == 1

    def test_door_arc_no_door(self, fresh_db):
        """Openings without doors should not have arc entries."""
        g = _geom_with_doors(fresh_db)
        names = {a["opening_name"] for a in g["door_arcs"]}
        # O1, O2, O4, O5, O7-O11 have no doors
        for name in ("O1", "O2", "O4", "O5", "O7", "O8", "O9", "O10", "O11"):
            assert name not in names

    def test_door_arc_points_count(self, fresh_db):
        """All door arc polylines should have 21 points (0..20)."""
        g = _geom_with_doors(fresh_db)
        for arc in g["door_arcs"]:
            for leaf in arc["leaves"]:
                assert len(leaf["arc_pts"]) == 21

    def test_door_arc_radius_matches_width(self, fresh_db):
        """Arc radius (distance from hinge to first arc point) should match door width."""
        g = _geom_with_doors(fresh_db)
        ro1 = next(a for a in g["door_arcs"] if a["opening_name"] == "RO1")
        leaf = ro1["leaves"][0]
        hx, hy = leaf["hinge"]
        p0 = leaf["arc_pts"][0]
        dist = math.sqrt((p0[0] - hx)**2 + (p0[1] - hy)**2)
        # RO1 door width is 36" = 3.0 ft
        assert abs(dist - 3.0) < 0.01


class TestClearanceZones:
    """Clearance zone geometry in engine output."""

    def test_clearance_zones_in_geometry(self, fresh_db):
        g = _geom_with_doors(fresh_db)
        assert "clearance_zones" in g
        assert len(g["clearance_zones"]) >= 1

    def test_dresser_clearance(self, fresh_db):
        g = _geom_with_doors(fresh_db)
        dresser = next(
            (z for z in g["clearance_zones"] if z["name"] == "dresser_clearance"),
            None)
        assert dresser is not None
        assert len(dresser["poly"]) == 4
        assert dresser["style"] == "dashed"

    def test_clearance_zone_bare(self, fresh_db):
        """Bare variant should have no clearance zones."""
        constants = get_constants_dict(fresh_db)
        chain_rows = get_outline_chain(fresh_db)
        doors = get_all_doors(fresh_db)
        g = compute_geometry(constants, "bare", chain_rows, doors_data=doors)
        assert g["clearance_zones"] == []


class TestDoorArcAPI:
    """Door mutations update arc geometry."""

    def test_door_arc_after_door_create(self, fresh_db):
        """Creating a door on an opening should add an arc entry."""
        # O1 doesn't have a door by default
        g1 = _geom_with_doors(fresh_db)
        o1_before = [a for a in g1["door_arcs"] if a["opening_name"] == "O1"]
        assert len(o1_before) == 0

        # Add a door to O1
        create_door("O1", 30, "north", "west", "single", fresh_db)
        g2 = _geom_with_doors(fresh_db)
        o1_after = [a for a in g2["door_arcs"] if a["opening_name"] == "O1"]
        assert len(o1_after) == 1

    def test_door_arc_after_door_delete(self, fresh_db):
        """Deleting a door should remove its arc entry."""
        g1 = _geom_with_doors(fresh_db)
        ro1_before = [a for a in g1["door_arcs"] if a["opening_name"] == "RO1"]
        assert len(ro1_before) == 1

        delete_door("RO1", fresh_db)
        g2 = _geom_with_doors(fresh_db)
        ro1_after = [a for a in g2["door_arcs"] if a["opening_name"] == "RO1"]
        assert len(ro1_after) == 0

    def test_door_invalidation_via_api(self, app_client):
        """PUT door hinge triggers geometry recomputation with new arcs."""
        # Use RO1: along=(1,0), so east/west flip moves hinge between ends
        r1 = app_client.get("/api/geometry?variant=standard")
        g1 = r1.get_json()
        ro1_1 = next(a for a in g1["door_arcs"] if a["opening_name"] == "RO1")
        hinge1 = ro1_1["leaves"][0]["hinge"]

        # Flip hinge from east to west
        r2 = app_client.put("/api/doors/RO1",
                            json={"hinge_side": "west"})
        assert r2.status_code == 200

        # Get new geometry — hinge should have moved
        r3 = app_client.get("/api/geometry?variant=standard")
        g3 = r3.get_json()
        ro1_3 = next(a for a in g3["door_arcs"] if a["opening_name"] == "RO1")
        hinge3 = ro1_3["leaves"][0]["hinge"]

        # Hinge position should differ after flip
        dist = math.sqrt((hinge3[0] - hinge1[0])**2 + (hinge3[1] - hinge1[1])**2)
        assert dist > 0.1  # should move significantly
