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
    return compute_geometry(constants, "standard", chain_rows,
                            doors_data=doors, db_path=fresh_db)


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
        g = compute_geometry(constants, "bare", chain_rows,
                             doors_data=doors, db_path=fresh_db)
        assert g["clearance_zones"] == []


class TestApplianceClearanceZones:
    """Clearance zones for stove, D/W, hamper from variant_items metadata."""

    def test_stove_clearance(self, fresh_db):
        g = _geom_with_doors(fresh_db)
        stove = next(
            (z for z in g["clearance_zones"] if z["name"] == "stove_clearance"),
            None)
        assert stove is not None
        assert len(stove["poly"]) == 4
        assert stove["style"] == "dashed"

    def test_dishwasher_clearance(self, fresh_db):
        g = _geom_with_doors(fresh_db)
        dw = next(
            (z for z in g["clearance_zones"] if z["name"] == "dishwasher_clearance"),
            None)
        assert dw is not None
        assert len(dw["poly"]) == 4

    def test_hamper_clearance(self, fresh_db):
        g = _geom_with_doors(fresh_db)
        hm = next(
            (z for z in g["clearance_zones"] if z["name"] == "hamper_clearance"),
            None)
        assert hm is not None
        assert len(hm["poly"]) == 4

    def test_clearance_count_standard(self, fresh_db):
        """Standard variant: dresser + stove + D/W + hamper = 4 clearance zones."""
        g = _geom_with_doors(fresh_db)
        assert len(g["clearance_zones"]) == 4

    def test_minik_no_stove_dw_clearance(self, fresh_db):
        """Minik variant has no stove or D/W, so fewer clearance zones."""
        constants = get_constants_dict(fresh_db)
        chain_rows = get_outline_chain(fresh_db)
        doors = get_all_doors(fresh_db)
        g = compute_geometry(constants, "minik", chain_rows,
                             doors_data=doors, db_path=fresh_db)
        names = {z["name"] for z in g["clearance_zones"]}
        assert "stove_clearance" not in names
        assert "dishwasher_clearance" not in names
        # Still has dresser and hamper
        assert "dresser_clearance" in names
        assert "hamper_clearance" in names

    def test_clearance_extends_outward(self, fresh_db):
        """Clearance polygon should extend away from item centroid."""
        g = _geom_with_doors(fresh_db)
        stove_cz = next(z for z in g["clearance_zones"]
                        if z["name"] == "stove_clearance")
        stove_item = g["variant_items"]["stove"]
        # Clearance poly[2] and poly[3] should be farther from stove centroid
        # than poly[0] and poly[1]
        cx = sum(p[0] for p in stove_item["poly"]) / len(stove_item["poly"])
        cy = sum(p[1] for p in stove_item["poly"]) / len(stove_item["poly"])
        face_mid = ((stove_cz["poly"][0][0] + stove_cz["poly"][1][0]) / 2,
                    (stove_cz["poly"][0][1] + stove_cz["poly"][1][1]) / 2)
        ext_mid = ((stove_cz["poly"][2][0] + stove_cz["poly"][3][0]) / 2,
                   (stove_cz["poly"][2][1] + stove_cz["poly"][3][1]) / 2)
        d_face = math.sqrt((face_mid[0] - cx)**2 + (face_mid[1] - cy)**2)
        d_ext = math.sqrt((ext_mid[0] - cx)**2 + (ext_mid[1] - cy)**2)
        assert d_ext > d_face


class TestApplianceDoors:
    """Appliance door arcs (fridge, washer, dryer, microwave)."""

    def test_appliance_doors_in_geometry(self, fresh_db):
        g = _geom_with_doors(fresh_db)
        assert "appliance_doors" in g
        names = {d["item_name"] for d in g["appliance_doors"]}
        assert "dryer" in names
        assert "washer" in names
        assert "fridge" in names
        assert "microwave" in names

    def test_appliance_door_count_standard(self, fresh_db):
        """Standard variant: dryer, washer, fridge, microwave = 4."""
        g = _geom_with_doors(fresh_db)
        assert len(g["appliance_doors"]) == 4

    def test_appliance_door_arc_points(self, fresh_db):
        """All appliance door arcs should have 21 points."""
        g = _geom_with_doors(fresh_db)
        for ad in g["appliance_doors"]:
            assert len(ad["arc_pts"]) == 21, f"{ad['item_name']} has {len(ad['arc_pts'])} points"

    def test_appliance_door_radius(self, fresh_db):
        """Arc radius should match door width (distance from hinge to arc[0])."""
        g = _geom_with_doors(fresh_db)
        for ad in g["appliance_doors"]:
            hx, hy = ad["hinge"]
            p0 = ad["arc_pts"][0]
            dist = math.sqrt((p0[0] - hx)**2 + (p0[1] - hy)**2)
            # Width stored in variant_items door metadata
            item = g["variant_items"][ad["item_name"]]
            expected_w = item["door"]["width"]
            assert abs(dist - expected_w) < 0.01, (
                f"{ad['item_name']}: arc radius {dist:.4f} != width {expected_w:.4f}")

    def test_appliance_door_tip_on_arc(self, fresh_db):
        """Tip should coincide with first arc point (open position)."""
        g = _geom_with_doors(fresh_db)
        for ad in g["appliance_doors"]:
            tip = ad["tip"]
            arc0 = ad["arc_pts"][0]
            dist = math.sqrt((tip[0] - arc0[0])**2 + (tip[1] - arc0[1])**2)
            assert dist < 0.01, f"{ad['item_name']}: tip-arc[0] dist = {dist:.6f}"

    def test_appliance_door_hinge_on_poly(self, fresh_db):
        """Hinge should be at the specified polygon vertex."""
        g = _geom_with_doors(fresh_db)
        for ad in g["appliance_doors"]:
            item = g["variant_items"][ad["item_name"]]
            idx = item["door"]["hinge_idx"]
            poly_pt = item["poly"][idx]
            dist = math.sqrt((ad["hinge"][0] - poly_pt[0])**2 +
                             (ad["hinge"][1] - poly_pt[1])**2)
            assert dist < 0.001, f"{ad['item_name']}: hinge not at poly[{idx}]"

    def test_minik_fridge_door(self, fresh_db):
        """Minik fridge door should still be present with different hinge."""
        constants = get_constants_dict(fresh_db)
        chain_rows = get_outline_chain(fresh_db)
        doors = get_all_doors(fresh_db)
        g = compute_geometry(constants, "minik", chain_rows,
                             doors_data=doors, db_path=fresh_db)
        fridge_doors = [d for d in g["appliance_doors"] if d["item_name"] == "fridge"]
        assert len(fridge_doors) == 1
        # Minik fridge hinge is at idx 1 (SE corner)
        item = g["variant_items"]["fridge"]
        assert item["door"]["hinge_idx"] == 1


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
