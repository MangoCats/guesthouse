"""Tests for Phase 10a Analysis features (ANALYSIS-1, ANALYSIS-2, ANALYSIS-3)."""
import os
import sys
import pytest

_PROJECT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _PROJECT not in sys.path:
    sys.path.insert(0, _PROJECT)

from tests.test_zapp_conftest import fresh_db, app_client, geometry
from app.database import get_constants_dict
from app.engine import compute_geometry


# =========================================================================
# TestRoomAreas (ANALYSIS-3)
# =========================================================================

class TestRoomAreas:
    """Room areas should be present for all variants, not just SF."""

    def test_standard_has_areas(self, fresh_db):
        constants = get_constants_dict(fresh_db)
        geom = compute_geometry(constants, "standard", db_path=fresh_db)
        for lbl in geom["room_labels"]:
            assert "area" in lbl, f"{lbl['name']} missing area"
            assert isinstance(lbl["area"], float), f"{lbl['name']} area not float"
            assert lbl["area"] > 0, f"{lbl['name']} area <= 0"

    def test_bare_has_areas(self, fresh_db):
        constants = get_constants_dict(fresh_db)
        geom = compute_geometry(constants, "bare", db_path=fresh_db)
        for lbl in geom["room_labels"]:
            assert "area" in lbl, f"{lbl['name']} missing area in bare"

    def test_sf_has_areas_and_polys(self, fresh_db):
        """SF variant should have both area values AND highlight polygons."""
        constants = get_constants_dict(fresh_db)
        geom = compute_geometry(constants, "sf", db_path=fresh_db)
        for lbl in geom["room_labels"]:
            assert "area" in lbl, f"{lbl['name']} missing area in sf"
            assert "poly" in lbl, f"{lbl['name']} missing poly in sf"

    def test_standard_no_polys(self, fresh_db):
        """Non-SF variants should NOT have highlight polygons."""
        constants = get_constants_dict(fresh_db)
        geom = compute_geometry(constants, "standard", db_path=fresh_db)
        for lbl in geom["room_labels"]:
            assert "poly" not in lbl, f"{lbl['name']} has poly in standard"

    def test_api_returns_areas(self, app_client):
        resp = app_client.get("/api/geometry?variant=standard")
        data = resp.get_json()
        for lbl in data["room_labels"]:
            assert "area" in lbl, f"{lbl['name']} missing area via API"

    def test_area_sum_reasonable(self, fresh_db):
        """Total area should be between 500 and 1200 SF."""
        constants = get_constants_dict(fresh_db)
        geom = compute_geometry(constants, "standard", db_path=fresh_db)
        total = sum(lbl["area"] for lbl in geom["room_labels"])
        assert 500 < total < 1200, f"Total area {total} out of expected range"


# =========================================================================
# TestSpanData (ANALYSIS-1)
# =========================================================================

class TestSpanData:
    """Span data API returns valid arrays."""

    def test_span_data_endpoint(self, app_client):
        resp = app_client.get("/api/span-data")
        assert resp.status_code == 200
        data = resp.get_json()
        assert "eastings" in data
        assert "spans" in data
        assert "south_spans" in data
        assert "north_spans" in data

    def test_span_data_arrays_same_length(self, app_client):
        data = app_client.get("/api/span-data").get_json()
        n = len(data["eastings"])
        assert n > 0
        assert len(data["spans"]) == n
        assert len(data["south_spans"]) == n
        assert len(data["north_spans"]) == n

    def test_span_data_values_positive(self, app_client):
        data = app_client.get("/api/span-data").get_json()
        assert max(data["spans"]) > 10, "Max span should be > 10 ft"

    def test_span_data_eastings_monotonic(self, app_client):
        data = app_client.get("/api/span-data").get_json()
        for i in range(1, len(data["eastings"])):
            assert data["eastings"][i] >= data["eastings"][i - 1], \
                "Eastings not monotonically increasing"


# =========================================================================
# TestSpanRotation (ANALYSIS-2)
# =========================================================================

class TestSpanRotation:
    """Span rotation analysis returns min/max with angles."""

    def test_span_rotation_endpoint(self, app_client):
        resp = app_client.get("/api/span-rotation")
        assert resp.status_code == 200
        data = resp.get_json()
        assert "min_angle" in data
        assert "min_span" in data
        assert "max_angle" in data
        assert "max_span" in data
        assert "data" in data

    def test_min_less_than_max(self, app_client):
        data = app_client.get("/api/span-rotation").get_json()
        assert data["min_span"] <= data["max_span"]

    def test_data_array_has_entries(self, app_client):
        data = app_client.get("/api/span-rotation").get_json()
        assert len(data["data"]) >= 30, "Should have 35 angle samples"
        for entry in data["data"]:
            assert len(entry) == 2, "Each entry should be [angle, span]"
            assert 0 < entry[0] <= 175
            assert entry[1] > 0

    def test_min_span_positive(self, app_client):
        data = app_client.get("/api/span-rotation").get_json()
        assert data["min_span"] > 0
        assert data["max_span"] > 0
