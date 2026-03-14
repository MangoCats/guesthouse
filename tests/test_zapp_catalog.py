"""Tests for catalog_items table and /api/catalog endpoints.

Covers: catalog seeding, CRUD operations, catalog delete endpoint,
element deletion with catalog removal.
"""
import json
import pytest

from tests.test_zapp_conftest import fresh_db, app_client


class TestCatalogSeeding:
    """catalog_items table is populated during init_db."""

    def test_catalog_items_seeded(self, app_client):
        """Seed database has catalog items."""
        resp = app_client.get("/api/catalog")
        assert resp.status_code == 200
        data = resp.get_json()
        # Should have all three categories
        assert "furniture" in data
        assert "appliance" in data
        assert "fixture" in data
        total = len(data["furniture"]) + len(data["appliance"]) + len(data["fixture"])
        assert total > 0, "catalog should contain seeded items"

    def test_catalog_has_known_items(self, app_client):
        """Known seed items appear in the catalog."""
        resp = app_client.get("/api/catalog")
        data = resp.get_json()
        all_keys = set()
        for cat_list in data.values():
            for item in cat_list:
                all_keys.add(item["key"])
        # Check a few known items from _VARIANT_ITEMS
        for key in ("dryer", "washer", "bed", "toilet_s", "bath_sink"):
            assert key in all_keys, f"expected {key} in catalog"

    def test_catalog_item_has_dimensions(self, app_client):
        """Seeded catalog items have width/depth or radius."""
        resp = app_client.get("/api/catalog")
        data = resp.get_json()
        all_items = []
        for cat_list in data.values():
            all_items.extend(cat_list)
        # At least some items should have non-zero dimensions
        has_dims = [i for i in all_items
                    if (i.get("width") and i["width"] > 0)
                    or (i.get("radius") and i["radius"] > 0)]
        assert len(has_dims) > 0, "some catalog items should have dimensions"

    def test_catalog_item_fields(self, app_client):
        """Each catalog item has required fields."""
        resp = app_client.get("/api/catalog")
        data = resp.get_json()
        for cat_type, items in data.items():
            for item in items:
                assert "key" in item
                assert "label" in item
                assert "type" in item
                assert "shape" in item


class TestCatalogCRUD:
    """CRUD operations on catalog_items via database functions."""

    def test_get_all_catalog_items(self, fresh_db):
        from app.database import get_all_catalog_items
        items = get_all_catalog_items(fresh_db)
        assert len(items) > 0

    def test_get_catalog_item(self, fresh_db):
        from app.database import get_catalog_item
        item = get_catalog_item("dryer", fresh_db)
        assert item is not None
        assert item["key"] == "dryer"
        assert item["item_type"] == "appliance"

    def test_get_nonexistent_catalog_item(self, fresh_db):
        from app.database import get_catalog_item
        item = get_catalog_item("nonexistent_item_xyz", fresh_db)
        assert item is None

    def test_create_catalog_item(self, fresh_db):
        from app.database import create_catalog_item, get_catalog_item
        data = {
            "key": "test_widget",
            "item_type": "furniture",
            "label": "Test Widget",
            "width": 2.5,
            "depth": 1.5,
            "shape": "rect",
            "variants": "[]",
        }
        create_catalog_item(data, fresh_db)
        item = get_catalog_item("test_widget", fresh_db)
        assert item is not None
        assert item["label"] == "Test Widget"
        assert item["width"] == 2.5

    def test_delete_catalog_item(self, fresh_db):
        from app.database import delete_catalog_item, get_catalog_item
        # Verify it exists first
        assert get_catalog_item("dryer", fresh_db) is not None
        delete_catalog_item("dryer", fresh_db)
        assert get_catalog_item("dryer", fresh_db) is None

    def test_delete_nonexistent_catalog_item(self, fresh_db):
        """Deleting a non-existent catalog item doesn't error."""
        from app.database import delete_catalog_item
        # Should not raise
        delete_catalog_item("nonexistent_xyz", fresh_db)


class TestCatalogEndpoints:
    """HTTP endpoints for catalog operations."""

    def test_get_catalog(self, app_client):
        resp = app_client.get("/api/catalog")
        assert resp.status_code == 200
        data = resp.get_json()
        assert isinstance(data, dict)
        assert all(isinstance(v, list) for v in data.values())

    def test_delete_catalog_item_endpoint(self, app_client):
        """DELETE /api/catalog/<key> removes the item."""
        # Verify dryer exists in catalog
        resp = app_client.get("/api/catalog")
        data = resp.get_json()
        appliance_keys = {i["key"] for i in data.get("appliance", [])}
        assert "dryer" in appliance_keys

        # Delete it
        resp = app_client.delete("/api/catalog/dryer")
        assert resp.status_code == 200
        assert resp.get_json()["ok"] is True

        # Verify it's gone
        resp = app_client.get("/api/catalog")
        data = resp.get_json()
        appliance_keys = {i["key"] for i in data.get("appliance", [])}
        assert "dryer" not in appliance_keys

    def test_delete_nonexistent_catalog_returns_404(self, app_client):
        resp = app_client.delete("/api/catalog/nonexistent_xyz")
        assert resp.status_code == 404

    def test_catalog_delete_does_not_affect_elements(self, app_client):
        """Removing a catalog item doesn't delete placed element instances."""
        # Verify dryer element exists
        elems = app_client.get("/api/elements").get_json()
        has_dryer = any(e["name"] == "dryer" for e in elems)
        assert has_dryer, "dryer element should exist in seed DB"

        # Delete from catalog
        resp = app_client.delete("/api/catalog/dryer")
        assert resp.status_code == 200

        # Element should still exist
        elems = app_client.get("/api/elements").get_json()
        assert any(e["name"] == "dryer" for e in elems), \
            "dryer element should remain after catalog deletion"
