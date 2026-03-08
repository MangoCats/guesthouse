"""Tests for app.server REST API — API-1 through API-15, API-32, API-33.

Verifies all REST endpoints using the Flask test client with an isolated
database.
"""
import json
import pytest

from tests.test_zapp_conftest import fresh_db, app_client


# ── API-1  GET /api/constants ────────────────────────────────────────

class TestAPI1GetConstants:
    def test_status_200(self, app_client):
        resp = app_client.get("/api/constants")
        assert resp.status_code == 200

    def test_returns_json_array(self, app_client):
        data = app_client.get("/api/constants").get_json()
        assert isinstance(data, list)
        assert len(data) >= 140

    def test_constant_has_required_keys(self, app_client):
        data = app_client.get("/api/constants").get_json()
        required = {"name", "value", "expr", "unit", "category", "description"}
        assert required <= set(data[0].keys())


# ── API-2  GET /api/constants?category=wall ──────────────────────────

class TestAPI2ConstantsByCategory:
    def test_filters_by_category(self, app_client):
        data = app_client.get("/api/constants?category=wall").get_json()
        assert len(data) > 0
        for c in data:
            assert c["category"] == "wall"

    def test_empty_category_returns_all(self, app_client):
        all_data = app_client.get("/api/constants").get_json()
        no_filter = app_client.get("/api/constants?category=").get_json()
        # Empty string doesn't match any category, returns empty
        # or returns all — depends on implementation
        # The code does `if cat:` so empty string is falsy, returns all
        assert len(no_filter) == len(all_data)


# ── API-3  GET /api/constants/categories ─────────────────────────────

class TestAPI3Categories:
    def test_status_200(self, app_client):
        resp = app_client.get("/api/constants/categories")
        assert resp.status_code == 200

    def test_returns_array(self, app_client):
        data = app_client.get("/api/constants/categories").get_json()
        assert isinstance(data, list)

    def test_contains_expected_categories(self, app_client):
        data = app_client.get("/api/constants/categories").get_json()
        cats = set(data)
        assert "wall" in cats
        assert "opening" in cats
        assert "furniture" in cats


# ── API-4  PUT /api/constants/<name> ─────────────────────────────────

class TestAPI4UpdateConstant:
    def test_update_returns_200(self, app_client):
        resp = app_client.put(
            "/api/constants/BED_WIDTH",
            data=json.dumps({"value": 7.0}),
            content_type="application/json",
        )
        assert resp.status_code == 200

    def test_update_response_fields(self, app_client):
        resp = app_client.put(
            "/api/constants/BED_WIDTH",
            data=json.dumps({"value": 7.0}),
            content_type="application/json",
        )
        data = resp.get_json()
        assert data["ok"] is True
        assert data["name"] == "BED_WIDTH"
        assert data["value"] == 7.0


# ── API-5  PUT /api/constants/<name> -- validation ───────────────────

class TestAPI5UpdateValidation:
    def test_missing_value_returns_400(self, app_client):
        resp = app_client.put(
            "/api/constants/BED_WIDTH",
            data=json.dumps({}),
            content_type="application/json",
        )
        assert resp.status_code == 400

    def test_invalid_value_returns_400(self, app_client):
        resp = app_client.put(
            "/api/constants/BED_WIDTH",
            data=json.dumps({"value": "abc"}),
            content_type="application/json",
        )
        assert resp.status_code == 400

    def test_nonexistent_constant_returns_404(self, app_client):
        resp = app_client.put(
            "/api/constants/NONEXISTENT_CONSTANT_XYZ",
            data=json.dumps({"value": 1.0}),
            content_type="application/json",
        )
        assert resp.status_code == 404


# ── API-6  PUT /api/constants/batch ──────────────────────────────────

class TestAPI6BatchUpdate:
    def test_batch_update(self, app_client):
        resp = app_client.put(
            "/api/constants/batch",
            data=json.dumps({"updates": {"BED_WIDTH": 7.0, "BED_LENGTH": 8.0}}),
            content_type="application/json",
        )
        assert resp.status_code == 200
        data = resp.get_json()
        assert data["ok"] is True
        assert data["changed"] == 2

    def test_batch_invalid_values_returns_400(self, app_client):
        resp = app_client.put(
            "/api/constants/batch",
            data=json.dumps({"updates": {"BED_WIDTH": "bad"}}),
            content_type="application/json",
        )
        assert resp.status_code == 400


# ── API-7  POST /api/constants/reset ─────────────────────────────────

class TestAPI7ResetConstants:
    def test_reset_returns_ok(self, app_client):
        resp = app_client.post("/api/constants/reset")
        assert resp.status_code == 200
        assert resp.get_json()["ok"] is True

    def test_reset_restores_values(self, app_client):
        # Get original
        orig = app_client.get("/api/constants").get_json()
        orig_bed = next(c for c in orig if c["name"] == "BED_WIDTH")

        # Modify
        app_client.put(
            "/api/constants/BED_WIDTH",
            data=json.dumps({"value": 999.0}),
            content_type="application/json",
        )

        # Reset
        app_client.post("/api/constants/reset")

        # Verify restored
        after = app_client.get("/api/constants").get_json()
        after_bed = next(c for c in after if c["name"] == "BED_WIDTH")
        assert abs(after_bed["value"] - orig_bed["value"]) < 1e-9


# ── API-8  GET /api/geometry ─────────────────────────────────────────

class TestAPI8GetGeometry:
    def test_status_200(self, app_client):
        resp = app_client.get("/api/geometry")
        assert resp.status_code == 200

    def test_geometry_has_required_keys(self, app_client):
        data = app_client.get("/api/geometry").get_json()
        for key in ("points", "outline_poly", "interior_walls", "bbox"):
            assert key in data, f"Missing key: {key}"


# ── API-10  GET /api/outline ─────────────────────────────────────────

class TestAPI10GetOutline:
    def test_status_200(self, app_client):
        resp = app_client.get("/api/outline")
        assert resp.status_code == 200

    def test_returns_18_segments(self, app_client):
        data = app_client.get("/api/outline").get_json()
        assert len(data) == 18

    def test_segment_has_required_keys(self, app_client):
        data = app_client.get("/api/outline").get_json()
        for seg in data:
            assert "seq" in seg
            assert "seg_type" in seg
            assert "end_name" in seg


# ── API-11  GET /api/views ───────────────────────────────────────────

class TestAPI11GetViews:
    def test_status_200(self, app_client):
        resp = app_client.get("/api/views")
        assert resp.status_code == 200

    def test_returns_at_least_11_views(self, app_client):
        data = app_client.get("/api/views").get_json()
        assert len(data) >= 11

    def test_view_has_required_keys(self, app_client):
        data = app_client.get("/api/views").get_json()
        for v in data:
            assert "name" in v
            assert "label" in v
            assert "script" in v
            assert "svg_path" in v
            assert "category" in v


# ── API-12  GET /api/svg/<view_name> ─────────────────────────────────

class TestAPI12GetSVG:
    def test_unknown_view_returns_404(self, app_client):
        resp = app_client.get("/api/svg/nonexistent")
        assert resp.status_code == 404

    def test_existing_view_returns_svg(self, app_client):
        import os
        svg_path = os.path.join(
            os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
            "floorplan", "floorplan.svg",
        )
        if not os.path.exists(svg_path):
            pytest.skip("floorplan.svg not yet generated")
        resp = app_client.get("/api/svg/floorplan")
        assert resp.status_code == 200
        assert resp.content_type.startswith("image/svg+xml")


# ── API-13  GET /api/svg/<view_name>/file ────────────────────────────

class TestAPI13GetSVGFile:
    def test_unknown_view_returns_404(self, app_client):
        resp = app_client.get("/api/svg/nonexistent/file")
        assert resp.status_code == 404

    def test_existing_view_serves_file(self, app_client):
        import os
        svg_path = os.path.join(
            os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
            "floorplan", "floorplan.svg",
        )
        if not os.path.exists(svg_path):
            pytest.skip("floorplan.svg not yet generated")
        resp = app_client.get("/api/svg/floorplan/file")
        assert resp.status_code == 200


# ── API-14  POST /api/regenerate ─────────────────────────────────────

class TestAPI14Regenerate:
    def test_regenerate_single_view(self, app_client):
        resp = app_client.post(
            "/api/regenerate",
            data=json.dumps({"view": "floorplan"}),
            content_type="application/json",
        )
        assert resp.status_code == 200
        data = resp.get_json()
        assert data["ok"] is True
        assert data["view"] == "floorplan"

    def test_regenerate_unknown_view_returns_404(self, app_client):
        resp = app_client.post(
            "/api/regenerate",
            data=json.dumps({"view": "nonexistent"}),
            content_type="application/json",
        )
        assert resp.status_code == 404

    def test_regenerate_all(self, app_client):
        resp = app_client.post("/api/regenerate")
        assert resp.status_code == 200
        data = resp.get_json()
        assert data["ok"] is True
        assert "results" in data


# ── API-15  Geometry Cache Invalidation ──────────────────────────────

class TestAPI15CacheInvalidation:
    def test_constant_change_invalidates_cache(self, app_client):
        # Prime cache
        geom1 = app_client.get("/api/geometry").get_json()
        bed1 = geom1["furniture"]["bed"]

        # Read current BED_WIDTH and change to something very different
        constants = app_client.get("/api/constants").get_json()
        cur_width = next(c["value"] for c in constants if c["name"] == "BED_WIDTH")
        new_width = cur_width + 2.0  # add 2 feet — clearly different

        app_client.put(
            "/api/constants/BED_WIDTH",
            data=json.dumps({"value": new_width}),
            content_type="application/json",
        )

        # Geometry should reflect the change
        geom2 = app_client.get("/api/geometry").get_json()
        bed2 = geom2["furniture"]["bed"]
        # At least one dimension should differ
        w1 = bed1["bbox"]["e"] - bed1["bbox"]["w"]
        w2 = bed2["bbox"]["e"] - bed2["bbox"]["w"]
        n1 = bed1["bbox"]["n"] - bed1["bbox"]["s"]
        n2 = bed2["bbox"]["n"] - bed2["bbox"]["s"]
        assert abs(w2 - w1) > 0.01 or abs(n2 - n1) > 0.01


# ── API-32  GET /api/events (SSE) ────────────────────────────────────

class TestAPI32SSE:
    def test_event_stream_content_type(self, app_client):
        resp = app_client.get("/api/events")
        assert "text/event-stream" in resp.content_type


class TestAPIVersion:
    def test_version_returns_git_and_started(self, app_client):
        resp = app_client.get("/api/version")
        assert resp.status_code == 200
        data = resp.get_json()
        assert "git" in data
        assert "started" in data
        assert len(data["started"]) == 19  # YYYY-MM-DD HH:MM:SS


# ── API-33  GET / (root) ────────────────────────────────────────────

class TestAPI33Root:
    def test_root_returns_200(self, app_client):
        resp = app_client.get("/")
        assert resp.status_code == 200

    def test_root_contains_adu_editor(self, app_client):
        resp = app_client.get("/")
        assert b"ADU Editor" in resp.data


# ── API-9  GET /api/geometry -- error handling ─────────────────────

class TestAPI9GeometryError:
    def test_geometry_error_returns_500(self, app_client, monkeypatch):
        """Geometry computation failure returns HTTP 500 with error message."""
        def _raise(*a, **kw):
            raise RuntimeError("test error")

        # Patch where the name is looked up (server.py's module globals)
        import app.server
        monkeypatch.setattr(app.server, "compute_geometry", _raise)
        resp = app_client.get("/api/geometry")
        assert resp.status_code == 500
        data = resp.get_json()
        assert "error" in data


# ── Variant SVG Serving ────────────────────────────────────────────

class TestVariantSVGServing:
    def test_floorplan_standard_serves_base_svg(self, app_client):
        """Standard variant serves the base floorplan.svg."""
        import os
        svg_path = os.path.join(
            os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
            "floorplan", "floorplan.svg",
        )
        if not os.path.exists(svg_path):
            pytest.skip("floorplan.svg not yet generated")
        resp = app_client.get("/api/svg/floorplan?variant=standard")
        assert resp.status_code == 200
        assert resp.content_type.startswith("image/svg+xml")

    def test_floorplan_minik_attempts_variant_svg(self, app_client):
        """Minik variant attempts to serve floorplan_minik.svg."""
        import os
        svg_path = os.path.join(
            os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
            "floorplan", "floorplan_minik.svg",
        )
        if not os.path.exists(svg_path):
            # If the variant SVG doesn't exist, the API returns 404 with a
            # message about regeneration — this is correct behaviour
            resp = app_client.get("/api/svg/floorplan?variant=minik")
            assert resp.status_code == 404
        else:
            resp = app_client.get("/api/svg/floorplan?variant=minik")
            assert resp.status_code == 200
            assert resp.content_type.startswith("image/svg+xml")


# ── Inner Wall Overrides API (Phase 15½-C) ─────────────────────────────

class TestInnerWallOverridesAPI:
    def test_get_all(self, app_client):
        resp = app_client.get("/api/inner-wall-overrides")
        assert resp.status_code == 200
        data = resp.get_json()
        assert "5" in data  # seeded override for seg 5

    def test_get_single(self, app_client):
        resp = app_client.get("/api/inner-wall-overrides/5")
        assert resp.status_code == 200
        chain = resp.get_json()
        assert len(chain) == 3
        assert chain[0]["seg_type"] == "L"
        assert chain[1]["seg_type"] == "CCW"
        assert chain[2]["seg_type"] == "L"

    def test_get_missing(self, app_client):
        resp = app_client.get("/api/inner-wall-overrides/99")
        assert resp.status_code == 404

    def test_upsert(self, app_client):
        chain = [
            {"seg_type": "L", "bearing": 180.0, "distance": 0.5},
            {"seg_type": "CCW", "radius": 0.2, "sweep": 90.0, "n_pts": 10},
            {"seg_type": "L", "bearing": 90.0, "distance": 0.5},
        ]
        resp = app_client.put("/api/inner-wall-overrides/5",
                              json={"chain": chain})
        assert resp.status_code == 200
        data = resp.get_json()
        assert data["ok"]
        assert len(data["chain"]) == 3

    def test_upsert_validation(self, app_client):
        # Missing chain
        resp = app_client.put("/api/inner-wall-overrides/5", json={})
        assert resp.status_code == 400
        # Invalid seg_type
        resp = app_client.put("/api/inner-wall-overrides/5",
                              json={"chain": [{"seg_type": "X"}]})
        assert resp.status_code == 400
        # Line missing bearing
        resp = app_client.put("/api/inner-wall-overrides/5",
                              json={"chain": [{"seg_type": "L", "distance": 1.0}]})
        assert resp.status_code == 400

    def test_delete(self, app_client):
        resp = app_client.delete("/api/inner-wall-overrides/5")
        assert resp.status_code == 200
        # Now it's gone
        resp = app_client.get("/api/inner-wall-overrides/5")
        assert resp.status_code == 404

    def test_delete_missing(self, app_client):
        resp = app_client.delete("/api/inner-wall-overrides/99")
        assert resp.status_code == 404
