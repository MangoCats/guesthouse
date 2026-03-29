"""Tests for app/engine.py helper functions:
_face_midpoint, _face_vertices, _find_wall_poly, _find_opening_poly,
_resolve_point_spec.
"""
import pytest
from app.engine import (
    _face_midpoint,
    _face_vertices,
    _find_wall_poly,
    _find_opening_poly,
    _resolve_point_spec,
)


# ---------------------------------------------------------------------------
# Shared test geometry
# ---------------------------------------------------------------------------

# Axis-aligned rectangle: south=0→1, east=1→2, north=2→3, west=3→0
_RECT = [[0, 0], [10, 0], [10, 10], [0, 10]]

# Minimal geometry dict for resolve_point_spec tests
_GEOM = {
    "points": {"F1": [1.0, 2.0]},
    "interior_walls": {
        "IW1": {"poly": [[0, 0], [10, 0], [10, 5], [0, 5]]},
    },
    "outer_openings": [
        {"name": "O1", "poly": [[4, 0], [6, 0], [6, 3], [4, 3]]},
    ],
    "rough_openings": [
        {"name": "RO1", "poly": [[1, 0], [4, 0], [4, 2], [1, 2]]},
    ],
}


# ---------------------------------------------------------------------------
# _face_midpoint
# ---------------------------------------------------------------------------

class TestFaceMidpoint:
    def test_south_face(self):
        result = _face_midpoint(_RECT, "south")
        assert result == pytest.approx([5.0, 0.0])

    def test_east_face(self):
        result = _face_midpoint(_RECT, "east")
        assert result == pytest.approx([10.0, 5.0])

    def test_north_face(self):
        result = _face_midpoint(_RECT, "north")
        assert result == pytest.approx([5.0, 10.0])

    def test_west_face(self):
        result = _face_midpoint(_RECT, "west")
        assert result == pytest.approx([0.0, 5.0])

    def test_invalid_face_name_returns_none(self):
        assert _face_midpoint(_RECT, "top") is None

    def test_poly_too_short_returns_none(self):
        assert _face_midpoint([[0, 0], [1, 0], [1, 1]], "south") is None


# ---------------------------------------------------------------------------
# _face_vertices
# ---------------------------------------------------------------------------

class TestFaceVertices:
    def test_south_vertices(self):
        a, b = _face_vertices(_RECT, "south")
        assert a == [0, 0]
        assert b == [10, 0]

    def test_east_vertices(self):
        a, b = _face_vertices(_RECT, "east")
        assert a == [10, 0]
        assert b == [10, 10]

    def test_invalid_face_returns_none_pair(self):
        a, b = _face_vertices(_RECT, "bad")
        assert a is None
        assert b is None


# ---------------------------------------------------------------------------
# _find_wall_poly
# ---------------------------------------------------------------------------

class TestFindWallPoly:
    def test_found_returns_poly(self):
        geom = {"interior_walls": {"IW1": {"poly": [[0, 0], [1, 0], [1, 1], [0, 1]]}}}
        poly = _find_wall_poly(geom, "IW1")
        assert poly == [[0, 0], [1, 0], [1, 1], [0, 1]]

    def test_not_found_returns_none(self):
        geom = {"interior_walls": {"IW1": {"poly": [[0, 0], [1, 0], [1, 1], [0, 1]]}}}
        assert _find_wall_poly(geom, "IW99") is None

    def test_empty_geom_returns_none(self):
        assert _find_wall_poly({}, "IW1") is None


# ---------------------------------------------------------------------------
# _find_opening_poly
# ---------------------------------------------------------------------------

class TestFindOpeningPoly:
    def test_found_in_outer_openings(self):
        poly = _find_opening_poly(_GEOM, "O1")
        assert poly == [[4, 0], [6, 0], [6, 3], [4, 3]]

    def test_found_in_rough_openings(self):
        poly = _find_opening_poly(_GEOM, "RO1")
        assert poly == [[1, 0], [4, 0], [4, 2], [1, 2]]

    def test_not_found_returns_none(self):
        assert _find_opening_poly(_GEOM, "O99") is None

    def test_entry_without_poly_key_skipped(self):
        geom = {
            "outer_openings": [{"name": "O1"}],  # no "poly" key
            "rough_openings": [],
        }
        assert _find_opening_poly(geom, "O1") is None


# ---------------------------------------------------------------------------
# _resolve_point_spec
# ---------------------------------------------------------------------------

class TestResolvePointSpec:
    def test_none_returns_none(self):
        assert _resolve_point_spec(None, _GEOM) is None

    def test_string_returns_named_point(self):
        result = _resolve_point_spec("F1", _GEOM)
        assert result == pytest.approx([1.0, 2.0])

    def test_missing_string_returns_none(self):
        assert _resolve_point_spec("MISSING", _GEOM) is None

    def test_face_mid_spec(self):
        """face_mid returns midpoint of named wall face."""
        result = _resolve_point_spec({"face_mid": "IW1", "face": "south"}, _GEOM)
        # IW1 south face: poly[0]=[0,0] → poly[1]=[10,0], midpoint=[5,0]
        assert result == pytest.approx([5.0, 0.0])

    def test_midpoint_spec(self):
        """midpoint averages two resolved point specs."""
        spec = {"midpoint": ["F1", {"face_mid": "IW1", "face": "south"}]}
        result = _resolve_point_spec(spec, _GEOM)
        # F1=[1,2], south mid=[5,0] → avg=[3,1]
        assert result == pytest.approx([3.0, 1.0])

    def test_opening_centroid_spec(self):
        """opening_centroid returns centroid of opening polygon."""
        spec = {"opening_centroid": "O1"}
        result = _resolve_point_spec(spec, _GEOM)
        # O1 poly: [[4,0],[6,0],[6,3],[4,3]] → centroid=(5, 1.5)
        assert result == pytest.approx([5.0, 1.5])

    def test_offset_spec(self):
        """offset moves base point by dist * direction."""
        spec = {"offset": "F1", "dir": "north", "dist": 3.0}
        result = _resolve_point_spec(spec, _GEOM)
        # F1=[1,2], north=[0,1], dist=3 → [1, 5]
        assert result == pytest.approx([1.0, 5.0])

    def test_face_mid_missing_wall_returns_none(self):
        spec = {"face_mid": "MISSING_WALL", "face": "south"}
        assert _resolve_point_spec(spec, _GEOM) is None

    def test_opening_face_mid_spec(self):
        """opening_face_mid returns midpoint of opening face."""
        spec = {"opening_face_mid": "O1", "face": "south"}
        result = _resolve_point_spec(spec, _GEOM)
        # O1 south: poly[0]=[4,0]→poly[1]=[6,0], midpoint=[5,0]
        assert result == pytest.approx([5.0, 0.0])
