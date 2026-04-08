"""Geometry output tests for FormulaEvaluator.

All assertions are on ev.elements[name] output keys (poly, bbox, etc.),
never on internal _eval_* method names.  Tests survive refactoring of
the dispatch layer or method names.

Coverage: all 12 element types.
  DB-free (no fixture): wall_rect, item_rect, item_circle, four_corner,
    wall_opening, toilet_shape, bath_sink_shape, dining_triangle,
    dining_chair, ellipse_rect, area_poly
  Requires DB (fresh_db): shape_transform
"""
import math
import pytest
from tests.test_zapp_conftest import fresh_db  # noqa: F401

from app.evaluator import FormulaEvaluator


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_ev(**overrides):
    """Create a FormulaEvaluator with minimal base context."""
    constants = {
        "WALL_T": 8 / 12,
        "CORNER_SW_R": 10 / 12,
        "SHELL_THICKNESS": 3.5 / 12,
        "AIR_GAP": 1 / 12,
        "OPENING_INSIDE_RADIUS": 1 / 12,
        **overrides.get("constants", {}),
    }
    pts = overrides.get("points", {
        "W1": [0, 0], "W2": [10, 0], "W3": [10, 10], "W4": [0, 10],
        "W5": [5, 0], "W9": [0, 5], "W10": [10, 5],
    })
    inner_poly = overrides.get("inner_poly", [
        [0, 0], [10, 0], [10, 10], [0, 10],
    ])
    radii = overrides.get("radii", {})
    return FormulaEvaluator(constants, pts, inner_poly, radii)


def _bbox_area(bbox):
    return (bbox["e"] - bbox["w"]) * (bbox["n"] - bbox["s"])


# ---------------------------------------------------------------------------
# wall_rect
# ---------------------------------------------------------------------------

class TestWallRectGeometry:
    def _formula(self, **kw):
        f = {
            "type": "wall_rect",
            "anchor": [0, 0],
            "along": "east",
            "thickness_dir": "north",
            "thickness": 8 / 12,
            "length": 10.0,
        }
        f.update(kw)
        return f

    def test_poly_has_four_vertices(self):
        ev = _make_ev()
        result = ev._evaluate_formula("wr", self._formula())
        assert result is not None
        assert len(result["poly"]) == 4

    def test_bbox_keys_present(self):
        ev = _make_ev()
        result = ev._evaluate_formula("wr", self._formula())
        assert {"w", "s", "e", "n"} <= set(result["bbox"])

    def test_length_matches_bbox_width(self):
        ev = _make_ev()
        result = ev._evaluate_formula("wr", self._formula(length=7.5))
        bb = result["bbox"]
        assert bb["e"] - bb["w"] == pytest.approx(7.5, abs=1e-9)

    def test_thickness_matches_bbox_height(self):
        ev = _make_ev()
        t = 6 / 12
        result = ev._evaluate_formula("wr", self._formula(thickness=t))
        bb = result["bbox"]
        assert bb["n"] - bb["s"] == pytest.approx(t, abs=1e-9)

    def test_anchor_is_sw_corner(self):
        ev = _make_ev()
        result = ev._evaluate_formula("wr", self._formula(anchor=[2.0, 3.0]))
        sw = result["poly"][0]
        assert sw == pytest.approx([2.0, 3.0], abs=1e-9)


# ---------------------------------------------------------------------------
# item_rect
# ---------------------------------------------------------------------------

class TestItemRectGeometry:
    def _formula(self, **kw):
        f = {
            "type": "item_rect",
            "anchor": [0, 0],
            "along": "east",
            "across": "north",
            "width": 2.0,
            "depth": 1.5,
        }
        f.update(kw)
        return f

    def test_poly_has_four_vertices(self):
        ev = _make_ev()
        result = ev._evaluate_formula("ir", self._formula())
        assert result is not None
        assert len(result["poly"]) == 4

    def test_bbox_matches_width(self):
        ev = _make_ev()
        result = ev._evaluate_formula("ir", self._formula(width=3.0))
        bb = result["bbox"]
        assert bb["e"] - bb["w"] == pytest.approx(3.0, abs=1e-9)

    def test_bbox_matches_depth(self):
        ev = _make_ev()
        result = ev._evaluate_formula("ir", self._formula(depth=2.5))
        bb = result["bbox"]
        assert bb["n"] - bb["s"] == pytest.approx(2.5, abs=1e-9)

    def test_center_key_at_midpoint(self):
        ev = _make_ev()
        result = ev._evaluate_formula("ir", self._formula(width=2.0, depth=1.0))
        assert "center" in result
        assert result["center"] == pytest.approx([1.0, 0.5], abs=1e-9)

    def test_width_and_depth_keys_in_result(self):
        ev = _make_ev()
        result = ev._evaluate_formula("ir", self._formula(width=3.0, depth=2.0))
        assert result["width"] == pytest.approx(3.0, abs=1e-9)
        assert result["depth"] == pytest.approx(2.0, abs=1e-9)

    def test_ne_corner_anchor(self):
        ev = _make_ev()
        # NE anchor at [4, 3], width=2, depth=1 → SW at [2, 2]
        result = ev._evaluate_formula("ir", self._formula(
            anchor=[4.0, 3.0], width=2.0, depth=1.0, anchor_corner="ne",
        ))
        bb = result["bbox"]
        assert bb["w"] == pytest.approx(2.0, abs=1e-9)
        assert bb["s"] == pytest.approx(2.0, abs=1e-9)


# ---------------------------------------------------------------------------
# item_circle
# ---------------------------------------------------------------------------

class TestItemCircleGeometry:
    def test_has_many_vertices(self):
        ev = _make_ev()
        result = ev._evaluate_formula("ic", {
            "type": "item_circle", "center": [5, 5], "radius": 1.0,
        })
        assert result is not None
        assert len(result["poly"]) >= 8

    def test_center_key_correct(self):
        ev = _make_ev()
        result = ev._evaluate_formula("ic", {
            "type": "item_circle", "center": [5, 5], "radius": 1.0,
        })
        assert "center" in result
        assert result["center"] == pytest.approx([5.0, 5.0], abs=1e-9)

    def test_radius_key_correct(self):
        ev = _make_ev()
        result = ev._evaluate_formula("ic", {
            "type": "item_circle", "center": [0, 0], "radius": 1.5,
        })
        assert "radius" in result
        assert result["radius"] == pytest.approx(1.5, abs=1e-9)

    def test_all_vertices_on_circle(self):
        ev = _make_ev()
        cx, cy, r = 3.0, 4.0, 2.0
        result = ev._evaluate_formula("ic", {
            "type": "item_circle", "center": [cx, cy], "radius": r,
        })
        for p in result["poly"]:
            d = math.sqrt((p[0] - cx) ** 2 + (p[1] - cy) ** 2)
            assert d == pytest.approx(r, abs=1e-9)

    def test_bbox_approx_diameter(self):
        ev = _make_ev()
        r = 1.0
        result = ev._evaluate_formula("ic", {
            "type": "item_circle", "center": [0, 0], "radius": r,
        })
        bb = result["bbox"]
        assert bb["e"] - bb["w"] == pytest.approx(2 * r, abs=1e-3)
        assert bb["n"] - bb["s"] == pytest.approx(2 * r, abs=1e-3)


# ---------------------------------------------------------------------------
# four_corner
# ---------------------------------------------------------------------------

class TestFourCornerGeometry:
    def _formula(self):
        return {
            "type": "four_corner",
            "sw": [0, 0], "se": [4, 0], "ne": [4, 3], "nw": [0, 3],
        }

    def test_poly_has_four_vertices(self):
        ev = _make_ev()
        result = ev._evaluate_formula("fc", self._formula())
        assert result is not None
        assert len(result["poly"]) == 4

    def test_poly_matches_input_corners(self):
        ev = _make_ev()
        sw, se, ne, nw = [0, 0], [4, 0], [4, 3], [0, 3]
        result = ev._evaluate_formula("fc", {
            "type": "four_corner", "sw": sw, "se": se, "ne": ne, "nw": nw,
        })
        poly = result["poly"]
        assert poly[0] == pytest.approx(sw, abs=1e-9)
        assert poly[1] == pytest.approx(se, abs=1e-9)
        assert poly[2] == pytest.approx(ne, abs=1e-9)
        assert poly[3] == pytest.approx(nw, abs=1e-9)

    def test_center_is_sw_ne_midpoint(self):
        ev = _make_ev()
        result = ev._evaluate_formula("fc", self._formula())
        assert "center" in result
        assert result["center"] == pytest.approx([2.0, 1.5], abs=1e-9)

    def test_bbox_covers_all_corners(self):
        ev = _make_ev()
        result = ev._evaluate_formula("fc", {
            "type": "four_corner",
            "sw": [1, 2], "se": [5, 2], "ne": [5, 6], "nw": [1, 6],
        })
        bb = result["bbox"]
        assert bb["w"] == pytest.approx(1.0, abs=1e-9)
        assert bb["s"] == pytest.approx(2.0, abs=1e-9)
        assert bb["e"] == pytest.approx(5.0, abs=1e-9)
        assert bb["n"] == pytest.approx(6.0, abs=1e-9)


# ---------------------------------------------------------------------------
# wall_opening
# ---------------------------------------------------------------------------

class TestWallOpeningGeometry:
    def _base_formula(self, **kw):
        f = {
            "type": "wall_opening",
            "outer_start": [0, 0],
            "outer_end": [10, 0],
            "inner_start": [0, 1],
            "inner_end": [10, 1],
            "width": 3.0,
            "centered": True,
            "center_t": 0.5,
        }
        f.update(kw)
        return f

    def test_poly_has_four_vertices(self):
        ev = _make_ev()
        result = ev._evaluate_formula("wo", self._base_formula())
        assert result is not None
        assert len(result["poly"]) == 4

    def test_bbox_has_positive_area(self):
        ev = _make_ev()
        result = ev._evaluate_formula("wo", self._base_formula())
        assert _bbox_area(result["bbox"]) > 0

    def test_centered_opening_width(self):
        ev = _make_ev()
        # Centered 3-ft opening on 10-ft wall → outer edge width = 3
        result = ev._evaluate_formula("wo", self._base_formula(width=3.0))
        outer_a, outer_b = result["poly"][0], result["poly"][1]
        w = math.dist(outer_a, outer_b)
        assert w == pytest.approx(3.0, abs=1e-9)

    def test_gap_from_start_positions_opening(self):
        ev = _make_ev()
        # gap=1 from start, width=2 → opening E 1.0..3.0
        result = ev._evaluate_formula("wo", {
            "type": "wall_opening",
            "outer_start": [0, 0], "outer_end": [10, 0],
            "inner_start": [0, 1], "inner_end": [10, 1],
            "width": 2.0, "from_end": False, "gap": 1.0,
        })
        outer_start = result["poly"][0]
        assert outer_start[0] == pytest.approx(1.0, abs=1e-9)


# ---------------------------------------------------------------------------
# toilet_shape
# ---------------------------------------------------------------------------

class TestToiletShapeGeometry:
    def _formula(self, **kw):
        f = {
            "type": "toilet_shape",
            "center": [0, 0],
            "facing_dir": "north",
            "width_dir": "east",
        }
        f.update(kw)
        return f

    def test_poly_has_many_vertices(self):
        ev = _make_ev()
        result = ev._evaluate_formula("toilet", self._formula())
        assert result is not None
        assert len(result["poly"]) > 4

    def test_center_key_correct(self):
        ev = _make_ev()
        result = ev._evaluate_formula("toilet", self._formula(center=[3, 4]))
        assert "center" in result
        assert result["center"] == pytest.approx([3.0, 4.0], abs=1e-9)

    def test_bbox_has_positive_area(self):
        ev = _make_ev()
        result = ev._evaluate_formula("toilet", self._formula())
        assert _bbox_area(result["bbox"]) > 0

    def test_facing_north_extends_northward(self):
        ev = _make_ev()
        result = ev._evaluate_formula("toilet", self._formula(center=[0, 0]))
        # Bowl projects outward (north of anchor)
        assert result["bbox"]["n"] > 0

    def test_poly_translates_with_center(self):
        ev = _make_ev()
        r0 = ev._evaluate_formula("t0", self._formula(center=[0, 0]))
        r1 = ev._evaluate_formula("t1", self._formula(center=[5, 3]))
        # Every vertex should shift by [5, 3]
        for p0, p1 in zip(r0["poly"], r1["poly"]):
            assert p1[0] == pytest.approx(p0[0] + 5, abs=1e-9)
            assert p1[1] == pytest.approx(p0[1] + 3, abs=1e-9)


# ---------------------------------------------------------------------------
# bath_sink_shape
# ---------------------------------------------------------------------------

class TestBathSinkShapeGeometry:
    def _formula(self, **kw):
        f = {
            "type": "bath_sink_shape",
            "anchor": [0, 0],
            "along": "east",
            "outward": "north",
            "length": 2.0,
            "depth": 1.0,
        }
        f.update(kw)
        return f

    def test_poly_has_more_than_four_vertices(self):
        ev = _make_ev()
        result = ev._evaluate_formula("sink", self._formula())
        assert result is not None
        assert len(result["poly"]) > 4

    def test_bbox_has_positive_area(self):
        ev = _make_ev()
        result = ev._evaluate_formula("sink", self._formula())
        assert _bbox_area(result["bbox"]) > 0

    def test_center_is_anchor(self):
        ev = _make_ev()
        result = ev._evaluate_formula("sink", self._formula(anchor=[1.5, 2.5]))
        assert "center" in result
        assert result["center"] == pytest.approx([1.5, 2.5], abs=1e-9)

    def test_bbox_width_matches_length(self):
        ev = _make_ev()
        result = ev._evaluate_formula("sink", self._formula(length=2.0))
        bb = result["bbox"]
        # along=east → e-w = length
        assert bb["e"] - bb["w"] == pytest.approx(2.0, abs=1e-9)

    def test_bbox_height_matches_depth(self):
        ev = _make_ev()
        result = ev._evaluate_formula("sink", self._formula(depth=1.5))
        bb = result["bbox"]
        # outward=north → n-s = depth
        assert bb["n"] - bb["s"] == pytest.approx(1.5, abs=1e-9)


# ---------------------------------------------------------------------------
# dining_triangle
# ---------------------------------------------------------------------------

class TestDiningTriangleGeometry:
    def _formula(self, **kw):
        f = {
            "type": "dining_triangle",
            "base_center": [0, 0],
            "toward_apex": "north",
            "along_base": "east",
            "base_width": 3.0,
            "height": 2.6,
            "apex_radius": 0.2,
            "fillet_radius": 0.1,
        }
        f.update(kw)
        return f

    def test_poly_has_vertices(self):
        ev = _make_ev()
        result = ev._evaluate_formula("tri", self._formula())
        assert result is not None
        assert len(result["poly"]) > 0

    def test_ne_side_key_is_two_points(self):
        ev = _make_ev()
        result = ev._evaluate_formula("tri", self._formula())
        assert "ne_side" in result
        assert len(result["ne_side"]) == 2

    def test_nw_side_key_is_two_points(self):
        ev = _make_ev()
        result = ev._evaluate_formula("tri", self._formula())
        assert "nw_side" in result
        assert len(result["nw_side"]) == 2

    def test_base_center_key_correct(self):
        ev = _make_ev()
        result = ev._evaluate_formula("tri", self._formula())
        assert "base_center" in result
        assert result["base_center"] == pytest.approx([0.0, 0.0], abs=1e-9)

    def test_bbox_has_positive_area(self):
        ev = _make_ev()
        result = ev._evaluate_formula("tri", self._formula())
        assert _bbox_area(result["bbox"]) > 0

    def test_apex_north_of_base(self):
        ev = _make_ev()
        result = ev._evaluate_formula("tri", self._formula())
        # Apex points north → bbox extends north of base_center
        assert result["bbox"]["n"] > 0

    def test_bbox_symmetric_east_west(self):
        ev = _make_ev()
        result = ev._evaluate_formula("tri", self._formula(base_center=[0, 0]))
        bb = result["bbox"]
        # Symmetric triangle → |east| ≈ |west|
        assert abs(bb["e"]) == pytest.approx(abs(bb["w"]), abs=1e-6)


# ---------------------------------------------------------------------------
# dining_chair (depends on dining_triangle via topo_sort)
# ---------------------------------------------------------------------------

class TestDiningChairGeometry:
    def _setup_ev(self, side="ne_side"):
        ev = _make_ev()
        ev.add_formula("tri_table", "position", {
            "type": "dining_triangle",
            "base_center": [0, 0],
            "toward_apex": "north",
            "along_base": "east",
            "base_width": 3.0,
            "height": 2.6,
            "apex_radius": 0.2,
            "fillet_radius": 0.1,
        })
        ev.add_formula("chair", "position", {
            "type": "dining_chair",
            "table": "tri_table",
            "side": side,
            "chair_width": 16 / 12,
            "chair_depth": 14 / 12,
            "gap": 2 / 12,
        })
        ev.topo_sort()
        ev.evaluate_all()
        return ev

    def test_poly_has_four_vertices(self):
        ev = self._setup_ev()
        assert "chair" in ev.elements
        assert len(ev.elements["chair"]["poly"]) == 4

    def test_chair_evaluated_after_triangle(self):
        ev = self._setup_ev()
        # Both elements should be present
        assert "tri_table" in ev.elements
        assert "chair" in ev.elements

    def test_chair_has_center_key(self):
        ev = self._setup_ev()
        assert "center" in ev.elements["chair"]

    def test_chair_ne_side_position(self):
        ev = self._setup_ev(side="ne_side")
        chair_c = ev.elements["chair"]["center"]
        # ne_side chair should be east of and/or north of base center
        assert chair_c[0] > 0 or chair_c[1] > 0

    def test_missing_table_returns_none(self):
        ev = _make_ev()
        result = ev._evaluate_formula("chair", {
            "type": "dining_chair",
            "table": "nonexistent_table",
            "side": "ne_side",
            "chair_width": 16 / 12,
            "chair_depth": 14 / 12,
        })
        assert result is None


# ---------------------------------------------------------------------------
# ellipse_rect
# ---------------------------------------------------------------------------

class TestEllipseRectGeometry:
    def _formula(self, **kw):
        f = {
            "type": "ellipse_rect",
            "anchor": [0, 0],
            "along": "east",
            "outward": "north",
            "rx": 1.5,
            "ry": 0.75,
        }
        f.update(kw)
        return f

    def test_poly_has_four_vertices(self):
        ev = _make_ev()
        result = ev._evaluate_formula("er", self._formula())
        assert result is not None
        assert len(result["poly"]) == 4

    def test_bbox_width_is_2rx(self):
        ev = _make_ev()
        rx = 1.5
        result = ev._evaluate_formula("er", self._formula(rx=rx))
        bb = result["bbox"]
        assert bb["e"] - bb["w"] == pytest.approx(2 * rx, abs=1e-9)

    def test_bbox_height_is_2ry(self):
        ev = _make_ev()
        ry = 0.75
        result = ev._evaluate_formula("er", self._formula(ry=ry))
        bb = result["bbox"]
        assert bb["n"] - bb["s"] == pytest.approx(2 * ry, abs=1e-9)

    def test_center_is_anchor(self):
        ev = _make_ev()
        result = ev._evaluate_formula("er", self._formula(anchor=[3.0, 4.0]))
        assert "center" in result
        assert result["center"] == pytest.approx([3.0, 4.0], abs=1e-9)

    def test_anchor_is_bbox_center(self):
        ev = _make_ev()
        result = ev._evaluate_formula("er", self._formula(anchor=[2.0, 5.0], rx=1.0, ry=0.5))
        bb = result["bbox"]
        cx = (bb["w"] + bb["e"]) / 2
        cy = (bb["s"] + bb["n"]) / 2
        assert cx == pytest.approx(2.0, abs=1e-9)
        assert cy == pytest.approx(5.0, abs=1e-9)


# ---------------------------------------------------------------------------
# area_poly
# ---------------------------------------------------------------------------

class TestAreaPolyGeometry:
    def test_poly_matches_vertices(self):
        ev = _make_ev()
        verts = [[0, 0], [4, 0], [4, 3], [0, 3]]
        result = ev._evaluate_formula("area", {
            "type": "area_poly", "vertices": verts,
        })
        assert result is not None
        assert len(result["poly"]) == 4
        for i, v in enumerate(verts):
            assert result["poly"][i] == pytest.approx(v, abs=1e-9)

    def test_area_rectangle(self):
        ev = _make_ev()
        # 4×3 rectangle → area = 12 sq ft
        result = ev._evaluate_formula("area", {
            "type": "area_poly",
            "vertices": [[0, 0], [4, 0], [4, 3], [0, 3]],
        })
        assert "area" in result
        assert result["area"] == pytest.approx(12.0, abs=1e-9)

    def test_area_right_triangle(self):
        ev = _make_ev()
        # Base=3, height=4 → area = 6
        result = ev._evaluate_formula("area", {
            "type": "area_poly",
            "vertices": [[0, 0], [3, 0], [0, 4]],
        })
        assert result["area"] == pytest.approx(6.0, abs=1e-9)

    def test_too_few_vertices_returns_none(self):
        ev = _make_ev()
        result = ev._evaluate_formula("area", {
            "type": "area_poly", "vertices": [[0, 0], [1, 0]],
        })
        assert result is None

    def test_bbox_covers_polygon(self):
        ev = _make_ev()
        result = ev._evaluate_formula("area", {
            "type": "area_poly",
            "vertices": [[1, 2], [5, 2], [5, 7], [1, 7]],
        })
        bb = result["bbox"]
        assert bb["w"] == pytest.approx(1.0, abs=1e-9)
        assert bb["s"] == pytest.approx(2.0, abs=1e-9)
        assert bb["e"] == pytest.approx(5.0, abs=1e-9)
        assert bb["n"] == pytest.approx(7.0, abs=1e-9)

    def test_named_point_vertices(self):
        ev = _make_ev()
        # Use named points from the evaluator's pts dict
        result = ev._evaluate_formula("area", {
            "type": "area_poly",
            "vertices": ["W1", "W2", "W3", "W4"],
        })
        assert result is not None
        # W1=[0,0], W2=[10,0], W3=[10,10], W4=[0,10] → area=100
        assert result["area"] == pytest.approx(100.0, abs=1e-9)


# ---------------------------------------------------------------------------
# shape_transform (requires DB — shape polygon data in shapes table)
# ---------------------------------------------------------------------------

class TestShapeTransformGeometry:
    def test_poly_has_vertices(self, fresh_db):
        ev = _make_ev()
        ev._db_path = fresh_db
        result = ev._evaluate_formula("st", {
            "type": "shape_transform",
            "shape_name": "toilet",
            "center": [5, 5],
            "rotation_deg": 0,
        })
        assert result is not None
        assert len(result["poly"]) >= 3

    def test_bbox_is_non_degenerate(self, fresh_db):
        ev = _make_ev()
        ev._db_path = fresh_db
        result = ev._evaluate_formula("st", {
            "type": "shape_transform",
            "shape_name": "toilet",
            "center": [5, 5],
            "rotation_deg": 0,
        })
        bb = result["bbox"]
        assert bb["e"] > bb["w"]
        assert bb["n"] > bb["s"]

    def test_center_matches_formula(self, fresh_db):
        ev = _make_ev()
        ev._db_path = fresh_db
        result = ev._evaluate_formula("st", {
            "type": "shape_transform",
            "shape_name": "toilet",
            "center": [3.0, 7.0],
            "rotation_deg": 0,
        })
        assert result["center"] == pytest.approx([3.0, 7.0], abs=1e-9)

    def test_rotation_swaps_bbox_dimensions(self, fresh_db):
        ev = _make_ev()
        ev._db_path = fresh_db
        r0 = ev._evaluate_formula("st0", {
            "type": "shape_transform",
            "shape_name": "toilet",
            "center": [5, 5],
            "rotation_deg": 0,
        })
        r90 = ev._evaluate_formula("st90", {
            "type": "shape_transform",
            "shape_name": "toilet",
            "center": [5, 5],
            "rotation_deg": 90,
        })
        w0 = r0["bbox"]["e"] - r0["bbox"]["w"]
        h0 = r0["bbox"]["n"] - r0["bbox"]["s"]
        w90 = r90["bbox"]["e"] - r90["bbox"]["w"]
        h90 = r90["bbox"]["n"] - r90["bbox"]["s"]
        # 90° rotation swaps width and height
        assert w90 == pytest.approx(h0, abs=1e-6)
        assert h90 == pytest.approx(w0, abs=1e-6)

    def test_unknown_shape_returns_none(self, fresh_db):
        ev = _make_ev()
        ev._db_path = fresh_db
        result = ev._evaluate_formula("st", {
            "type": "shape_transform",
            "shape_name": "nonexistent_shape_xyz",
            "center": [0, 0],
        })
        assert result is None

    def test_bath_sink_shape_also_seeded(self, fresh_db):
        ev = _make_ev()
        ev._db_path = fresh_db
        result = ev._evaluate_formula("st", {
            "type": "shape_transform",
            "shape_name": "bath_sink",
            "center": [0, 0],
            "rotation_deg": 0,
        })
        assert result is not None
        assert len(result["poly"]) >= 3
