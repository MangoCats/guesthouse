"""Tests for variant furniture/appliance computation and API."""
import math
import os
import sys
import pytest

_PROJECT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _PROJECT not in sys.path:
    sys.path.insert(0, _PROJECT)

from tests.test_zapp_conftest import fresh_db, app_client, geometry
from app.variants import compute_variant_items, VARIANTS
from app.database import get_constants_dict
from app.engine import compute_geometry


# =========================================================================
# Fixtures
# =========================================================================

@pytest.fixture
def variant_items_standard(fresh_db):
    """Standard variant items from a fresh database."""
    constants = get_constants_dict(fresh_db)
    geom = compute_geometry(constants, "standard")
    return geom["variant_items"]


@pytest.fixture
def all_variant_items(fresh_db):
    """Dict mapping variant name to its items from a fresh database."""
    constants = get_constants_dict(fresh_db)
    result = {}
    for v in VARIANTS:
        geom = compute_geometry(constants, v)
        result[v] = geom["variant_items"]
    return result


# =========================================================================
# TestVariantRegistry
# =========================================================================

class TestVariantRegistry:
    def test_all_five_variants_defined(self):
        assert set(VARIANTS.keys()) == {"standard", "minik", "daybed", "bare", "sf"}

    def test_variant_labels_non_empty(self):
        for name, info in VARIANTS.items():
            assert "label" in info
            assert len(info["label"]) > 0, f"{name} has empty label"

    def test_variant_flags_are_dicts(self):
        for name, info in VARIANTS.items():
            assert isinstance(info["flags"], dict), f"{name} flags not a dict"


# =========================================================================
# TestComputeVariantItems
# =========================================================================

class TestComputeVariantItems:
    def test_standard_has_many_items(self, variant_items_standard):
        """Standard variant should have ~25+ items."""
        assert len(variant_items_standard) >= 20, \
            f"Expected >=20 items, got {len(variant_items_standard)}: {list(variant_items_standard.keys())}"

    def test_standard_has_expected_items(self, variant_items_standard):
        expected = {"dryer", "washer", "hamper", "counter", "water_heater",
                    "toilet_s", "toilet_n", "util_sink", "bath_sink",
                    "stove", "dishwasher", "kitchen_sink", "fridge",
                    "ice_maker", "work_counter", "microwave", "north_counter",
                    "coffee_maker", "dining_table", "dining_chair_1", "dining_chair_2",
                    "bed", "dresser", "shelves", "loveseat", "et", "loveseat2",
                    "chair", "ottoman", "desk", "desk_chair"}
        actual = set(variant_items_standard.keys())
        missing = expected - actual
        assert not missing, f"Missing items: {missing}"

    def test_minik_has_minik_items(self, all_variant_items):
        mk = all_variant_items["minik"]
        assert "sofa" in mk, "minik should have sofa"
        assert "rocker" in mk, "minik should have rocker"
        assert "loveseat" not in mk, "minik should not have loveseat"
        assert "loveseat2" not in mk, "minik should not have loveseat2"
        assert "stove" not in mk, "minik should not have stove"
        assert "dishwasher" not in mk, "minik should not have dishwasher"
        assert "cooktop" in mk, "minik should have cooktop"
        assert "toaster" in mk, "minik should have toaster"

    def test_daybed_has_daybed_items(self, all_variant_items):
        db = all_variant_items["daybed"]
        assert "daybed" in db, "daybed variant should have daybed"
        assert "shelves2" in db, "daybed variant should have shelves2"
        assert "et_east" in db, "daybed variant should have et_east"
        assert "et_west" in db, "daybed variant should have et_west"
        assert "rocker" in db, "daybed variant should have rocker"
        assert "loveseat" not in db, "daybed variant should not have loveseat"
        assert "sofa" not in db, "daybed variant should not have sofa"

    def test_bare_has_no_items(self, all_variant_items):
        assert len(all_variant_items["bare"]) == 0

    def test_sf_has_no_items(self, all_variant_items):
        assert len(all_variant_items["sf"]) == 0


# =========================================================================
# TestItemGeometry
# =========================================================================

class TestItemGeometry:
    def test_all_items_have_poly(self, variant_items_standard):
        for name, item in variant_items_standard.items():
            assert "poly" in item, f"{name} missing poly"
            assert len(item["poly"]) >= 3, f"{name} poly has < 3 points"

    def test_all_items_have_bbox(self, variant_items_standard):
        for name, item in variant_items_standard.items():
            assert "bbox" in item, f"{name} missing bbox"
            b = item["bbox"]
            assert b["e"] >= b["w"], f"{name} bbox east < west"
            assert b["n"] >= b["s"], f"{name} bbox north < south"

    def test_all_items_have_type(self, variant_items_standard):
        valid_types = {"appliance", "furniture", "fixture"}
        for name, item in variant_items_standard.items():
            assert item["type"] in valid_types, f"{name} has invalid type {item['type']}"

    def test_water_heater_is_circle(self, variant_items_standard):
        wh = variant_items_standard.get("water_heater")
        assert wh is not None, "water_heater missing"
        assert wh["shape"] == "circle"
        assert "center" in wh
        assert "radius" in wh
        assert wh["radius"] > 0

    def test_bed_position_matches_layout(self, fresh_db):
        """Bed poly from variant items should match layout.bed.poly."""
        constants = get_constants_dict(fresh_db)
        geom = compute_geometry(constants, "standard")
        vi_bed = geom["variant_items"]["bed"]
        layout_bed = geom["furniture"]["bed"]
        # Compare first point of each
        vi_p0 = vi_bed["poly"][0]
        ly_p0 = layout_bed["poly"][0]
        assert abs(vi_p0[0] - ly_p0[0]) < 1e-4, "bed E mismatch"
        assert abs(vi_p0[1] - ly_p0[1]) < 1e-4, "bed N mismatch"

    def test_desk_along_w16_w17(self, variant_items_standard):
        desk = variant_items_standard.get("desk")
        assert desk is not None, "desk missing"
        # Desk should have a reasonable size (60" x 30" = 5' x 2.5')
        b = desk["bbox"]
        w = b["e"] - b["w"]
        h = b["n"] - b["s"]
        assert max(w, h) > 3.0, "desk should be at least 3 ft in one dimension"

    def test_chair_has_four_corners(self, variant_items_standard):
        ch = variant_items_standard.get("chair")
        assert ch is not None, "chair missing"
        assert len(ch["poly"]) == 4, "chair should have 4 corners"

    def test_et_is_circle_in_standard(self, variant_items_standard):
        et = variant_items_standard.get("et")
        assert et is not None, "et missing in standard"
        assert et["shape"] == "circle"

    def test_all_items_within_outline(self, fresh_db):
        """Every item bbox should be within outline bbox (with margin)."""
        constants = get_constants_dict(fresh_db)
        geom = compute_geometry(constants, "standard")
        ob = geom["bbox"]
        margin = 2.0  # 2 ft margin for tolerance
        for name, item in geom["variant_items"].items():
            b = item["bbox"]
            assert b["w"] > ob["w"] - margin, f"{name} west {b['w']} outside outline {ob['w']}"
            assert b["e"] < ob["e"] + margin, f"{name} east {b['e']} outside outline {ob['e']}"
            assert b["s"] > ob["s"] - margin, f"{name} south {b['s']} outside outline {ob['s']}"
            assert b["n"] < ob["n"] + margin, f"{name} north {b['n']} outside outline {ob['n']}"


# =========================================================================
# TestVariantAPI
# =========================================================================

class TestVariantAPI:
    def test_variants_endpoint(self, app_client):
        resp = app_client.get("/api/variants")
        assert resp.status_code == 200
        data = resp.get_json()
        assert len(data) == 5
        names = {v["name"] for v in data}
        assert names == {"standard", "minik", "daybed", "bare", "sf"}

    def test_geometry_default_variant_is_standard(self, app_client):
        resp = app_client.get("/api/geometry")
        assert resp.status_code == 200
        data = resp.get_json()
        assert data["variant"] == "standard"
        assert "variant_items" in data
        assert len(data["variant_items"]) >= 20

    def test_geometry_with_variant_param(self, app_client):
        resp = app_client.get("/api/geometry?variant=bare")
        assert resp.status_code == 200
        data = resp.get_json()
        assert data["variant"] == "bare"
        assert len(data["variant_items"]) == 0

    def test_geometry_minik_variant(self, app_client):
        resp = app_client.get("/api/geometry?variant=minik")
        assert resp.status_code == 200
        data = resp.get_json()
        assert data["variant"] == "minik"
        assert "sofa" in data["variant_items"]
        assert "loveseat" not in data["variant_items"]

    def test_unknown_variant_returns_standard(self, app_client):
        resp = app_client.get("/api/geometry?variant=nonexistent")
        assert resp.status_code == 200
        data = resp.get_json()
        assert data["variant"] == "standard"

    def test_variant_items_in_response_structure(self, app_client):
        resp = app_client.get("/api/geometry")
        data = resp.get_json()
        assert "available_variants" in data
        assert set(data["available_variants"]) == {"standard", "minik", "daybed", "bare", "sf"}
        # Each item should have required keys
        for name, item in data["variant_items"].items():
            assert "type" in item
            assert "poly" in item
            assert "bbox" in item
            assert "label" in item
            assert "shape" in item


# =========================================================================
# TestDimensionData
# =========================================================================

class TestDimensionData:
    def test_geometry_has_dimensions(self, fresh_db):
        """Geometry result should include dimension line data via user_dimensions."""
        constants = get_constants_dict(fresh_db)
        geom = compute_geometry(constants, "standard", db_path=fresh_db)
        assert "user_dimensions" in geom
        builtin = [d for d in geom["user_dimensions"]
                   if d["properties"].get("source") == "builtin"]
        assert len(builtin) >= 18

    def test_dimension_has_endpoints(self, fresh_db):
        """Each builtin dimension should have start, end endpoints."""
        constants = get_constants_dict(fresh_db)
        geom = compute_geometry(constants, "standard", db_path=fresh_db)
        import math
        for dim in geom["user_dimensions"]:
            p = dim["properties"]
            if p.get("source") != "builtin":
                continue
            name = dim["name"]
            assert "start" in p and "end" in p, f"{name} missing endpoints"
            assert len(p["start"]) == 2 and len(p["end"]) == 2, \
                f"{name} endpoints not [E,N] pairs"
            dist = math.sqrt((p["end"][0] - p["start"][0])**2 +
                             (p["end"][1] - p["start"][1])**2)
            assert dist > 0, f"{name} has non-positive dist"

    def test_bare_variant_has_extra_dims(self, fresh_db):
        """Bare variant should include dim20 and dim21."""
        constants = get_constants_dict(fresh_db)
        geom = compute_geometry(constants, "bare", db_path=fresh_db)
        dim_names = {d["name"] for d in geom["user_dimensions"]}
        assert "dim20" in dim_names
        assert "dim21" in dim_names

    def test_dimension_api_response(self, app_client):
        """API should return dimensions in geometry response."""
        resp = app_client.get("/api/geometry")
        data = resp.get_json()
        assert "user_dimensions" in data
        builtin = [d for d in data["user_dimensions"]
                   if d["properties"].get("source") == "builtin"]
        assert len(builtin) >= 18
        # Spot-check one dimension
        import math
        dim01 = next((d for d in builtin if d["name"] == "dim01"), None)
        assert dim01 is not None
        p = dim01["properties"]
        dist = math.sqrt((p["end"][0] - p["start"][0])**2 +
                         (p["end"][1] - p["start"][1])**2)
        assert dist > 0
