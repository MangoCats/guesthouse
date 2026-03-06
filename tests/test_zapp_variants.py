"""Tests for variant furniture/appliance computation and API."""
import json
import math
import os
import sys
import pytest

_PROJECT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _PROJECT not in sys.path:
    sys.path.insert(0, _PROJECT)

from tests.test_zapp_conftest import fresh_db, app_client, geometry
from app.variants import (compute_variant_items, VARIANTS,
                          _resolve_product_url, _PRODUCT_URLS_BASE,
                          _PRODUCT_URLS_VARIANT, get_variant_flags)
from app.database import (get_constants_dict, get_variants, get_variant,
                          get_variant_by_id, update_variant)
from app.engine import compute_geometry


# =========================================================================
# Fixtures
# =========================================================================

@pytest.fixture
def variant_items_standard(fresh_db):
    """Standard variant items from a fresh database."""
    constants = get_constants_dict(fresh_db)
    geom = compute_geometry(constants, "standard", db_path=fresh_db)
    return geom["variant_items"]


@pytest.fixture
def all_variant_items(fresh_db):
    """Dict mapping variant name to its items from a fresh database."""
    constants = get_constants_dict(fresh_db)
    result = {}
    for v in VARIANTS:
        geom = compute_geometry(constants, v, db_path=fresh_db)
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
        geom = compute_geometry(constants, "standard", db_path=fresh_db)
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
        geom = compute_geometry(constants, "standard", db_path=fresh_db)
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
# TestProductUrls — _resolve_product_url and variant item URLs
# =========================================================================

class TestProductUrls:
    """Product URL resolution from _PRODUCT_URLS_BASE and _PRODUCT_URLS_VARIANT."""

    def test_base_url_returns_for_any_variant(self):
        """Items in _PRODUCT_URLS_BASE get the same URL regardless of flags."""
        url = _resolve_product_url("hamper", False, False)
        assert url and "homedepot.com" in url
        assert _resolve_product_url("hamper", True, False) == url
        assert _resolve_product_url("hamper", False, True) == url

    def test_variant_url_minik(self):
        """Fridge in minik gets the IKEA URL, not the Lowe's default."""
        url = _resolve_product_url("fridge", True, False)
        assert url and "ikea.com" in url

    def test_variant_url_default(self):
        """Fridge in standard (not minik, not db) gets the default Lowe's URL."""
        url = _resolve_product_url("fridge", False, False)
        assert url and "lowes.com" in url

    def test_variant_url_db_flag(self):
        """Rocker with db flag gets the db-specific URL."""
        url = _resolve_product_url("rocker", False, True)
        assert url and "ikea.com" in url

    def test_minik_only_items_absent_in_standard(self):
        """Dryer/washer have URLs only for minik; standard returns None."""
        assert _resolve_product_url("dryer", False, False) is None
        assert _resolve_product_url("dryer", True, False) is not None

    def test_unknown_item_returns_none(self):
        assert _resolve_product_url("nonexistent_xyz", False, False) is None

    def test_standard_items_have_urls(self, variant_items_standard):
        """Items with base URLs should have product_url in variant output."""
        for name in ("hamper", "kitchen_sink", "shelves", "dining_table"):
            item = variant_items_standard.get(name)
            assert item is not None, f"{name} missing"
            assert "product_url" in item, f"{name} has no product_url"
            assert item["product_url"].startswith("http"), f"{name} URL invalid"

    def test_minik_fridge_url_differs_from_standard(self, all_variant_items):
        """Fridge URL should differ between standard and minik."""
        std_url = all_variant_items["standard"]["fridge"]["product_url"]
        mk_url = all_variant_items["minik"]["fridge"]["product_url"]
        assert std_url != mk_url

    def test_minik_has_dryer_washer_urls(self, all_variant_items):
        mk = all_variant_items["minik"]
        assert "product_url" in mk["dryer"]
        assert "product_url" in mk["washer"]

    def test_standard_dryer_washer_no_url(self, variant_items_standard):
        """Dryer/washer in standard have no product URL (minik-only)."""
        assert "product_url" not in variant_items_standard.get("dryer", {})
        assert "product_url" not in variant_items_standard.get("washer", {})

    def test_api_returns_product_urls(self, app_client):
        """Geometry API response includes product_url on items."""
        resp = app_client.get("/api/geometry?variant=standard")
        data = resp.get_json()
        shelves = data["variant_items"].get("shelves")
        assert shelves is not None
        assert "product_url" in shelves


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


# =========================================================================
# TestVariantsTable — Phase 11a: database-driven variant registry
# =========================================================================

class TestVariantsTable:
    def test_variants_table_exists(self, fresh_db):
        """Schema should include variants table."""
        import sqlite3
        conn = sqlite3.connect(fresh_db)
        conn.row_factory = sqlite3.Row
        rows = conn.execute(
            "SELECT name FROM sqlite_master WHERE type='table' AND name='variants'"
        ).fetchall()
        conn.close()
        assert len(rows) == 1

    def test_builtin_variants_seeded(self, fresh_db):
        """5 built-in variants should be seeded."""
        variants = get_variants(fresh_db)
        assert len(variants) == 5
        names = {v["name"] for v in variants}
        assert names == {"standard", "minik", "daybed", "bare", "sf"}

    def test_builtin_labels(self, fresh_db):
        """Built-in variants have correct labels."""
        variants = get_variants(fresh_db)
        label_map = {v["name"]: v["label"] for v in variants}
        assert label_map["standard"] == "Standard"
        assert label_map["minik"] == "Small Kitchen"
        assert label_map["daybed"] == "Daybed"
        assert label_map["bare"] == "Room Dimensions"
        assert label_map["sf"] == "Square Footage"

    def test_builtin_flags(self, fresh_db):
        """Built-in variants have correct flags."""
        v = get_variant("minik", fresh_db)
        assert v["flags"] == {"minik": True}
        v = get_variant("daybed", fresh_db)
        assert v["flags"] == {"db": True}
        v = get_variant("bare", fresh_db)
        assert v["flags"] == {"bare": True}
        v = get_variant("standard", fresh_db)
        assert v["flags"] == {}

    def test_builtin_flag_set(self, fresh_db):
        """All built-ins have is_builtin=True."""
        for v in get_variants(fresh_db):
            assert v["is_builtin"] is True, f"{v['name']} not builtin"

    def test_get_variant_by_name(self, fresh_db):
        """get_variant returns correct record."""
        v = get_variant("sf", fresh_db)
        assert v is not None
        assert v["name"] == "sf"
        assert v["label"] == "Square Footage"
        assert v["flags"] == {"sf": True}

    def test_get_variant_by_id(self, fresh_db):
        """get_variant_by_id returns correct record."""
        v = get_variant("standard", fresh_db)
        v2 = get_variant_by_id(v["id"], fresh_db)
        assert v2["name"] == "standard"

    def test_get_variant_unknown_returns_none(self, fresh_db):
        """get_variant for unknown name returns None."""
        assert get_variant("nonexistent", fresh_db) is None

    def test_update_layer_config(self, fresh_db):
        """update_variant updates layer_config JSON."""
        v = get_variant("standard", fresh_db)
        cfg = {"points": False, "labels": True, "dims": True}
        updated = update_variant(v["id"], {"layer_config": cfg}, fresh_db)
        assert updated["layer_config"] == cfg
        # Re-read to confirm persistence
        v2 = get_variant("standard", fresh_db)
        assert v2["layer_config"] == cfg

    def test_update_preserves_other_fields(self, fresh_db):
        """Partial update doesn't clobber name, flags, etc."""
        v = get_variant("minik", fresh_db)
        update_variant(v["id"], {"layer_config": {"grid": True}}, fresh_db)
        v2 = get_variant("minik", fresh_db)
        assert v2["name"] == "minik"
        assert v2["label"] == "Small Kitchen"
        assert v2["flags"] == {"minik": True}
        assert v2["layer_config"] == {"grid": True}

    def test_default_layer_config_empty(self, fresh_db):
        """New DB has empty layer_config for all variants."""
        for v in get_variants(fresh_db):
            assert v["layer_config"] == {}, f"{v['name']} has non-empty layer_config"


# =========================================================================
# TestVariantsAPI_DB — Phase 11a: API tests for DB-driven variants
# =========================================================================

class TestVariantsAPI_DB:
    def test_get_variants_from_db(self, app_client):
        """GET /api/variants returns 5 variants with all fields."""
        resp = app_client.get("/api/variants")
        assert resp.status_code == 200
        data = resp.get_json()
        assert len(data) == 5
        names = {v["name"] for v in data}
        assert names == {"standard", "minik", "daybed", "bare", "sf"}
        # Each variant should have DB fields
        for v in data:
            assert "id" in v
            assert "label" in v
            assert "flags" in v
            assert "layer_config" in v
            assert "is_builtin" in v

    def test_put_layer_config(self, app_client):
        """PUT /api/variants/<id> updates layer_config."""
        # Get standard variant id
        resp = app_client.get("/api/variants")
        variants = resp.get_json()
        std = next(v for v in variants if v["name"] == "standard")
        # Update layer config
        cfg = {"points": False, "furniture": True, "grid": True}
        resp = app_client.put(f"/api/variants/{std['id']}",
                              data=json.dumps({"layer_config": cfg}),
                              content_type="application/json")
        assert resp.status_code == 200
        updated = resp.get_json()
        assert updated["layer_config"] == cfg

    def test_put_layer_config_undo(self, app_client):
        """Undo restores previous layer_config."""
        resp = app_client.get("/api/variants")
        std = next(v for v in resp.get_json() if v["name"] == "standard")
        # Set initial config
        cfg1 = {"dims": True}
        app_client.put(f"/api/variants/{std['id']}",
                       data=json.dumps({"layer_config": cfg1}),
                       content_type="application/json")
        # Set new config
        cfg2 = {"dims": False, "grid": True}
        app_client.put(f"/api/variants/{std['id']}",
                       data=json.dumps({"layer_config": cfg2}),
                       content_type="application/json")
        # Verify cfg2 is current
        resp = app_client.get("/api/variants")
        std2 = next(v for v in resp.get_json() if v["name"] == "standard")
        assert std2["layer_config"] == cfg2
        # Undo
        resp = app_client.post("/api/undo")
        assert resp.status_code == 200
        # Verify cfg1 is restored
        resp = app_client.get("/api/variants")
        std3 = next(v for v in resp.get_json() if v["name"] == "standard")
        assert std3["layer_config"] == cfg1

    def test_put_nonexistent_variant_404(self, app_client):
        """PUT with bad id returns 404."""
        resp = app_client.put("/api/variants/9999",
                              data=json.dumps({"layer_config": {}}),
                              content_type="application/json")
        assert resp.status_code == 404

    def test_geometry_accepts_db_variants(self, app_client):
        """Geometry endpoint validates variant via DB."""
        resp = app_client.get("/api/geometry?variant=minik")
        assert resp.status_code == 200
        data = resp.get_json()
        assert data["variant"] == "minik"

    def test_geometry_unknown_variant_fallback(self, app_client):
        """Unknown variant falls back to standard via DB lookup."""
        resp = app_client.get("/api/geometry?variant=nonexistent")
        assert resp.status_code == 200
        data = resp.get_json()
        assert data["variant"] == "standard"


# =========================================================================
# TestVariantFlagsLookup — Phase 11a: get_variant_flags
# =========================================================================

class TestVariantFlagsLookup:
    def test_get_variant_flags_builtin(self, fresh_db):
        """get_variant_flags returns correct flags from DB."""
        flags = get_variant_flags("minik", fresh_db)
        assert flags == {"minik": True}
        flags = get_variant_flags("daybed", fresh_db)
        assert flags == {"db": True}

    def test_get_variant_flags_standard_empty(self, fresh_db):
        """Standard variant has empty flags."""
        flags = get_variant_flags("standard", fresh_db)
        assert flags == {}

    def test_get_variant_flags_unknown_fallback(self, fresh_db):
        """Unknown variant name falls back to standard flags."""
        flags = get_variant_flags("nonexistent", fresh_db)
        assert flags == {}

    def test_compute_variant_items_uses_db_flags(self, fresh_db):
        """compute_variant_items works with DB-driven flags."""
        constants = get_constants_dict(fresh_db)
        geom = compute_geometry(constants, "minik", db_path=fresh_db)
        assert "sofa" in geom["variant_items"]
        assert "loveseat" not in geom["variant_items"]
