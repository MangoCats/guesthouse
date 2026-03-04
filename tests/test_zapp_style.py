"""Tests for Phase 9: Styling and Product Links.

Covers: STYLE-1–4, LINK-1–2, SEL-12, CV-12 (8 reqs).
"""
import json

import pytest

from app.style import (
    TYPE_DEFAULTS, STYLE_KEYS, VALID_STROKE_STYLES,
    validate_color, validate_opacity, validate_stroke_style,
    validate_stroke_width, validate_style_props,
    get_defaults, resolve_style,
)
from app.database import (
    get_all_elements, create_element, update_element, get_element, init_db,
)

from tests.test_zapp_conftest import fresh_db, app_client  # noqa: F401


# ── Style defaults ──────────────────────────────────────────────────

class TestStyleDefaults:

    def test_all_types_have_defaults(self):
        for t in ("appliance", "furniture", "fixture", "wall", "opening",
                   "dimension", "label"):
            assert t in TYPE_DEFAULTS

    def test_default_keys(self):
        for t, d in TYPE_DEFAULTS.items():
            for key in STYLE_KEYS:
                assert key in d, f"{t} missing {key}"

    def test_appliance_defaults_match_css(self):
        d = TYPE_DEFAULTS["appliance"]
        assert d["fill_color"] == "#2a3a4a"
        assert d["stroke_color"] == "#4682B4"
        assert d["stroke_width"] == 0.02

    def test_furniture_defaults_match_css(self):
        d = TYPE_DEFAULTS["furniture"]
        assert d["fill_color"] == "#2a3a2a"
        assert d["stroke_color"] == "#5a8a5a"

    def test_fixture_defaults_match_css(self):
        d = TYPE_DEFAULTS["fixture"]
        assert d["fill_color"] == "#3a2a3a"
        assert d["stroke_color"] == "#8a5a8a"

    def test_wall_defaults_match_css(self):
        d = TYPE_DEFAULTS["wall"]
        assert d["fill_color"] == "#334"
        assert d["stroke_color"] == "#556"

    def test_opening_defaults_match_css(self):
        d = TYPE_DEFAULTS["opening"]
        assert d["fill_color"] == "#2a4a6a"
        assert d["stroke_color"] == "#4488cc"

    def test_get_defaults_returns_copy(self):
        d1 = get_defaults("appliance")
        d2 = get_defaults("appliance")
        assert d1 == d2
        d1["fill_color"] = "#000"
        assert d2["fill_color"] != "#000"

    def test_get_defaults_unknown_type(self):
        d = get_defaults("unknown_type")
        assert d["opacity"] == 100
        assert d["fill_color"] is None


# ── Validation ──────────────────────────────────────────────────────

class TestStyleValidation:

    def test_valid_hex_3(self):
        assert validate_color("#abc") is True

    def test_valid_hex_6(self):
        assert validate_color("#aabbcc") is True

    def test_none_color_valid(self):
        assert validate_color(None) is True

    def test_invalid_color_name(self):
        assert validate_color("red") is False

    def test_invalid_color_no_hash(self):
        assert validate_color("aabbcc") is False

    def test_invalid_color_bad_chars(self):
        assert validate_color("#gg0000") is False

    def test_invalid_color_int(self):
        assert validate_color(42) is False

    def test_opacity_zero(self):
        assert validate_opacity(0) is True

    def test_opacity_hundred(self):
        assert validate_opacity(100) is True

    def test_opacity_fifty(self):
        assert validate_opacity(50) is True

    def test_opacity_float(self):
        assert validate_opacity(33.5) is True

    def test_opacity_negative(self):
        assert validate_opacity(-1) is False

    def test_opacity_over(self):
        assert validate_opacity(101) is False

    def test_opacity_string(self):
        assert validate_opacity("50") is False

    def test_stroke_style_solid(self):
        assert validate_stroke_style("solid") is True

    def test_stroke_style_dashed(self):
        assert validate_stroke_style("dashed") is True

    def test_stroke_style_dotted(self):
        assert validate_stroke_style("dotted") is True

    def test_stroke_style_invalid(self):
        assert validate_stroke_style("dash") is False

    def test_stroke_style_none(self):
        assert validate_stroke_style(None) is False

    def test_stroke_width_zero(self):
        assert validate_stroke_width(0) is True

    def test_stroke_width_positive(self):
        assert validate_stroke_width(0.02) is True

    def test_stroke_width_negative(self):
        assert validate_stroke_width(-0.01) is False

    def test_stroke_width_string(self):
        assert validate_stroke_width("0.02") is False

    def test_validate_style_props_valid(self):
        ok, err = validate_style_props({
            "fill_color": "#ff6600",
            "stroke_color": "#cc4400",
            "stroke_width": 0.03,
            "stroke_style": "dashed",
            "opacity": 80,
        })
        assert ok is True
        assert err is None

    def test_validate_style_props_empty(self):
        ok, err = validate_style_props({})
        assert ok is True

    def test_validate_style_props_bad_color(self):
        ok, err = validate_style_props({"fill_color": "red"})
        assert ok is False
        assert "fill_color" in err

    def test_validate_style_props_bad_opacity(self):
        ok, err = validate_style_props({"opacity": 200})
        assert ok is False
        assert "opacity" in err

    def test_validate_style_props_not_dict(self):
        ok, err = validate_style_props("not a dict")
        assert ok is False

    def test_validate_view_overrides_valid(self):
        ok, err = validate_style_props({
            "view_overrides": {"bare": {"opacity": 20}},
        })
        assert ok is True

    def test_validate_view_overrides_invalid_opacity(self):
        ok, err = validate_style_props({
            "view_overrides": {"bare": {"opacity": 200}},
        })
        assert ok is False
        assert "bare" in err
        assert "opacity" in err

    def test_validate_view_overrides_invalid_color(self):
        ok, err = validate_style_props({
            "view_overrides": {"sf": {"fill_color": "red"}},
        })
        assert ok is False
        assert "sf" in err

    def test_validate_view_overrides_not_dict(self):
        ok, err = validate_style_props({"view_overrides": "not a dict"})
        assert ok is False
        assert "view_overrides" in err

    def test_validate_view_overrides_value_not_dict(self):
        ok, err = validate_style_props({"view_overrides": {"bare": "bad"}})
        assert ok is False
        assert "bare" in err

    def test_validate_view_overrides_empty(self):
        ok, err = validate_style_props({"view_overrides": {}})
        assert ok is True

    def test_validate_product_url_valid(self):
        ok, err = validate_style_props({"product_url": "https://example.com"})
        assert ok is True

    def test_validate_product_url_http(self):
        ok, err = validate_style_props({"product_url": "http://example.com"})
        assert ok is True

    def test_validate_product_url_invalid_scheme(self):
        ok, err = validate_style_props({"product_url": "ftp://example.com"})
        assert ok is False
        assert "product_url" in err

    def test_validate_product_url_no_scheme(self):
        ok, err = validate_style_props({"product_url": "example.com"})
        assert ok is False

    def test_validate_product_url_null(self):
        ok, err = validate_style_props({"product_url": None})
        assert ok is True

    def test_validate_product_url_empty_string(self):
        ok, err = validate_style_props({"product_url": ""})
        assert ok is True


# ── Resolution ──────────────────────────────────────────────────────

class TestStyleResolution:

    def test_resolve_defaults_only(self):
        result = resolve_style("appliance", None)
        assert result == TYPE_DEFAULTS["appliance"]

    def test_resolve_defaults_empty_props(self):
        result = resolve_style("furniture", {})
        assert result == TYPE_DEFAULTS["furniture"]

    def test_resolve_base_override(self):
        result = resolve_style("appliance", {"fill_color": "#ff0000"})
        assert result["fill_color"] == "#ff0000"
        assert result["stroke_color"] == "#4682B4"  # unchanged default

    def test_resolve_partial_override(self):
        result = resolve_style("furniture", {"opacity": 50})
        assert result["opacity"] == 50
        assert result["fill_color"] == "#2a3a2a"

    def test_resolve_view_override(self):
        props = {
            "opacity": 80,
            "view_overrides": {
                "bare": {"opacity": 20},
            },
        }
        result = resolve_style("appliance", props, variant="bare")
        assert result["opacity"] == 20  # view override wins
        assert result["fill_color"] == "#2a3a4a"  # default

    def test_resolve_view_override_partial(self):
        props = {
            "fill_color": "#ff0000",
            "view_overrides": {
                "sf": {"opacity": 50},
            },
        }
        result = resolve_style("furniture", props, variant="sf")
        assert result["fill_color"] == "#ff0000"  # base override preserved
        assert result["opacity"] == 50  # view override applied

    def test_resolve_no_matching_view(self):
        props = {
            "opacity": 80,
            "view_overrides": {
                "bare": {"opacity": 20},
            },
        }
        result = resolve_style("appliance", props, variant="standard")
        assert result["opacity"] == 80  # base override, not view

    def test_resolve_no_variant_param(self):
        props = {
            "opacity": 80,
            "view_overrides": {
                "bare": {"opacity": 20},
            },
        }
        result = resolve_style("appliance", props, variant=None)
        assert result["opacity"] == 80

    def test_resolve_unknown_type(self):
        result = resolve_style("unknown", {"fill_color": "#abc"})
        assert result["fill_color"] == "#abc"
        assert result["opacity"] == 100


# ── Product URL in DB ───────────────────────────────────────────────

class TestProductUrl:

    def test_url_stored_in_properties(self, fresh_db):
        elem = create_element("furniture", "test_item",
                              {"source": "placed", "product_url": "https://example.com"},
                              db_path=fresh_db)
        props = json.loads(elem["properties"]) if isinstance(elem["properties"], str) else elem["properties"]
        assert props["product_url"] == "https://example.com"

    def test_url_persists_across_reads(self, fresh_db):
        elem = create_element("furniture", "test_item",
                              {"source": "placed", "product_url": "https://example.com"},
                              db_path=fresh_db)
        read_back = get_element(elem["id"], fresh_db)
        props = json.loads(read_back["properties"]) if isinstance(read_back["properties"], str) else read_back["properties"]
        assert props["product_url"] == "https://example.com"

    def test_url_cleared_with_null(self, fresh_db):
        elem = create_element("furniture", "test_item",
                              {"source": "placed", "product_url": "https://example.com"},
                              db_path=fresh_db)
        update_element(elem["id"], {"properties": {"source": "placed", "product_url": None}},
                       db_path=fresh_db)
        read_back = get_element(elem["id"], fresh_db)
        props = json.loads(read_back["properties"]) if isinstance(read_back["properties"], str) else read_back["properties"]
        assert props.get("product_url") is None


# ── Style via API ───────────────────────────────────────────────────

class TestStyleAPI:

    def test_put_style_properties(self, app_client):
        # Create element
        resp = app_client.post("/api/elements", json={
            "type": "furniture", "name": "styled_item",
            "properties": {"source": "placed", "poly": [[0,0],[1,0],[1,1],[0,1]]},
        })
        assert resp.status_code == 201
        elem_id = resp.get_json()["id"]
        # Update with style
        resp = app_client.put(f"/api/elements/{elem_id}", json={
            "properties": {
                "source": "placed", "poly": [[0,0],[1,0],[1,1],[0,1]],
                "fill_color": "#ff0000", "opacity": 75,
            },
        })
        assert resp.status_code == 200
        props = resp.get_json()["properties"]
        if isinstance(props, str):
            props = json.loads(props)
        assert props["fill_color"] == "#ff0000"
        assert props["opacity"] == 75

    def test_put_view_overrides(self, app_client):
        resp = app_client.post("/api/elements", json={
            "type": "appliance", "name": "ov_item",
            "properties": {"source": "placed", "poly": [[0,0],[1,0],[1,1],[0,1]]},
        })
        elem_id = resp.get_json()["id"]
        resp = app_client.put(f"/api/elements/{elem_id}", json={
            "properties": {
                "source": "placed", "poly": [[0,0],[1,0],[1,1],[0,1]],
                "view_overrides": {"bare": {"opacity": 20}},
            },
        })
        assert resp.status_code == 200
        props = resp.get_json()["properties"]
        if isinstance(props, str):
            props = json.loads(props)
        assert props["view_overrides"]["bare"]["opacity"] == 20

    def test_put_product_url(self, app_client):
        resp = app_client.post("/api/elements", json={
            "type": "fixture", "name": "url_item",
            "properties": {"source": "placed", "poly": [[0,0],[1,0],[1,1],[0,1]]},
        })
        elem_id = resp.get_json()["id"]
        resp = app_client.put(f"/api/elements/{elem_id}", json={
            "properties": {
                "source": "placed", "poly": [[0,0],[1,0],[1,1],[0,1]],
                "product_url": "https://www.ikea.com/product/123",
            },
        })
        assert resp.status_code == 200
        props = resp.get_json()["properties"]
        if isinstance(props, str):
            props = json.loads(props)
        assert props["product_url"] == "https://www.ikea.com/product/123"

    def test_style_in_geometry_response(self, app_client):
        """Style properties pass through to geometry for dimension elements."""
        # Create a dimension with style properties
        resp = app_client.post("/api/elements", json={
            "type": "dimension", "name": "styled_dim",
            "properties": {
                "source": "user",
                "start": [0, 0], "end": [5, 0], "offset": 0.5,
                "stroke_color": "#ff6600", "opacity": 60,
            },
        })
        assert resp.status_code == 201
        # Fetch geometry
        resp = app_client.get("/api/geometry?variant=standard")
        assert resp.status_code == 200
        g = resp.get_json()
        dims = g.get("user_dimensions", [])
        styled = [d for d in dims if d["name"] == "styled_dim"]
        assert len(styled) == 1
        props = styled[0]["properties"]
        assert props["stroke_color"] == "#ff6600"
        assert props["opacity"] == 60


# ── Per-view resolution end-to-end ──────────────────────────────────

class TestPerViewStyleResolution:

    def test_resolve_standard_no_override(self):
        props = {"fill_color": "#ff0000"}
        r = resolve_style("furniture", props, variant="standard")
        assert r["fill_color"] == "#ff0000"
        assert r["opacity"] == 100

    def test_resolve_bare_with_override(self):
        props = {
            "fill_color": "#ff0000",
            "view_overrides": {"bare": {"opacity": 20}},
        }
        r = resolve_style("furniture", props, variant="bare")
        assert r["opacity"] == 20
        assert r["fill_color"] == "#ff0000"

    def test_override_does_not_affect_base(self):
        props = {
            "fill_color": "#ff0000",
            "opacity": 80,
            "view_overrides": {"bare": {"opacity": 20}},
        }
        r_std = resolve_style("furniture", props, variant="standard")
        r_bare = resolve_style("furniture", props, variant="bare")
        assert r_std["opacity"] == 80
        assert r_bare["opacity"] == 20

    def test_multiple_view_overrides(self):
        props = {
            "view_overrides": {
                "bare": {"opacity": 20},
                "sf": {"fill_color": "#cccccc", "opacity": 50},
            },
        }
        r_bare = resolve_style("furniture", props, variant="bare")
        r_sf = resolve_style("furniture", props, variant="sf")
        assert r_bare["opacity"] == 20
        assert r_bare["fill_color"] == "#2a3a2a"  # default
        assert r_sf["opacity"] == 50
        assert r_sf["fill_color"] == "#cccccc"
