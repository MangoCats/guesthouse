"""Tests for the parametric formula evaluator (Phase 12a).

Covers: FormulaEvaluator point/dir/length resolution, wall_rect and
item_rect evaluation, topological sort, cycle detection, dependency
extraction, and database formula CRUD.
"""
import json
import math
import pytest
from tests.test_zapp_conftest import fresh_db  # noqa: F401

from app.evaluator import (
    FormulaEvaluator,
    CycleError,
    extract_deps,
    _line_poly_intersections,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_evaluator(**kwargs):
    """Create an evaluator with minimal base context."""
    constants = kwargs.get("constants", {"WALL_6IN": 0.5, "WALL_4IN": 1/3})
    points = kwargs.get("points", {
        "W1": [0, 0], "W2": [10, 0], "W3": [10, 10], "W4": [0, 10],
        "W5": [10, 0], "W9": [0, 5], "W10": [10, 5],
    })
    inner_poly = kwargs.get("inner_poly", [
        [0, 0], [10, 0], [10, 10], [0, 10],
    ])
    radii = kwargs.get("radii", {})
    return FormulaEvaluator(constants, points, inner_poly, radii)


# ---------------------------------------------------------------------------
# Point resolution tests
# ---------------------------------------------------------------------------

class TestResolvePoint:
    """Test point spec resolution."""

    def test_named_point(self):
        ev = _make_evaluator()
        assert ev.resolve_point("W1") == [0, 0]
        assert ev.resolve_point("W2") == [10, 0]

    def test_unknown_point_returns_none(self):
        ev = _make_evaluator()
        assert ev.resolve_point("NONEXISTENT") is None

    def test_literal_list(self):
        ev = _make_evaluator()
        assert ev.resolve_point([3.5, 7.2]) == [3.5, 7.2]

    def test_none_returns_none(self):
        ev = _make_evaluator()
        assert ev.resolve_point(None) is None

    def test_midpoint(self):
        ev = _make_evaluator()
        result = ev.resolve_point({"midpoint": ["W1", "W2"]})
        assert result == pytest.approx([5, 0])

    def test_offset_literal_dist(self):
        ev = _make_evaluator()
        result = ev.resolve_point({"offset": "W1", "dir": "east", "dist": 3.0})
        assert result == pytest.approx([3, 0])

    def test_offset_const_dist(self):
        ev = _make_evaluator()
        result = ev.resolve_point(
            {"offset": "W1", "dir": "north", "dist": {"const": "WALL_6IN"}}
        )
        assert result == pytest.approx([0, 0.5])

    def test_line_intersection(self):
        ev = _make_evaluator()
        result = ev.resolve_point({
            "type": "line_intersection",
            "line1_point": "W1", "line1_dir": "east",
            "line2_point": [5, -1], "line2_dir": "north",
        })
        assert result == pytest.approx([5, 0])

    def test_element_corner_int(self):
        ev = _make_evaluator()
        ev.elements["IW1"] = {"poly": [[1, 2], [3, 2], [3, 4], [1, 4]]}
        assert ev.resolve_point({"element": "IW1", "corner": 0}) == [1, 2]
        assert ev.resolve_point({"element": "IW1", "corner": 2}) == [3, 4]

    def test_element_corner_named(self):
        ev = _make_evaluator()
        ev.elements["IW1"] = {"poly": [[1, 2], [3, 2], [3, 4], [1, 4]]}
        assert ev.resolve_point({"element": "IW1", "corner": "SW"}) == [1, 2]
        assert ev.resolve_point({"element": "IW1", "corner": "NE"}) == [3, 4]

    def test_element_corner_missing_element(self):
        ev = _make_evaluator()
        assert ev.resolve_point({"element": "MISSING", "corner": 0}) is None

    def test_face_mid(self):
        ev = _make_evaluator()
        ev.elements["IW1"] = {"poly": [[0, 0], [4, 0], [4, 2], [0, 2]]}
        result = ev.resolve_point({"face_mid": "IW1", "face": "south"})
        assert result == pytest.approx([2, 0])
        result = ev.resolve_point({"face_mid": "IW1", "face": "north"})
        assert result == pytest.approx([2, 2])

    def test_line_poly_isect_nearest(self):
        ev = _make_evaluator()
        result = ev.resolve_point({
            "line_poly_isect": {
                "origin": [5, 5], "dir": "east", "poly": "inner_poly",
            }
        })
        assert result is not None
        assert result[0] == pytest.approx(10)  # east wall of square
        assert result[1] == pytest.approx(5)

    def test_arc_point(self):
        ev = _make_evaluator(
            points={"C1": [0, 0], "REF": [0, 5]},
            radii={"R1": 5.0},
        )
        result = ev.resolve_point({
            "arc_point": {
                "center": "C1", "radius_key": "R1",
                "reference": "REF", "side": "east",
            }
        })
        assert result is not None
        # Point on circle at northing=5 → (0 + sqrt(25-25), 5) = (0, 5)
        # Actually ref northing is 5, center at 0,0, R=5
        # dn = ref[1] - center[1] = 5
        # disc = 25 - 25 = 0, de = 0
        assert result == pytest.approx([0, 5], abs=1e-9)


# ---------------------------------------------------------------------------
# Direction resolution tests
# ---------------------------------------------------------------------------

class TestResolveDir:
    """Test direction spec resolution."""

    def test_cardinal_directions(self):
        ev = _make_evaluator()
        assert ev.resolve_dir("east") == [1, 0]
        assert ev.resolve_dir("north") == [0, 1]
        assert ev.resolve_dir("west") == [-1, 0]
        assert ev.resolve_dir("south") == [0, -1]

    def test_none_returns_none(self):
        ev = _make_evaluator()
        assert ev.resolve_dir(None) is None

    def test_segment(self):
        ev = _make_evaluator()
        result = ev.resolve_dir({"segment": ["W1", "W2"]})
        assert result == pytest.approx([1, 0])  # W1(0,0) → W2(10,0) = east

    def test_segment_perp(self):
        ev = _make_evaluator()
        result = ev.resolve_dir({"segment_perp": ["W1", "W2"]})
        # Perpendicular of east = [0, -1] (right-hand rule: dy/len, -dx/len)
        assert result == pytest.approx([0, -1])

    def test_face_along(self):
        ev = _make_evaluator()
        ev.elements["IW1"] = {"poly": [[0, 0], [4, 0], [4, 2], [0, 2]]}
        result = ev.resolve_dir({"face_along": "IW1", "face": "south"})
        assert result == pytest.approx([1, 0])  # south face goes east

    def test_face_perp(self):
        ev = _make_evaluator()
        ev.elements["IW1"] = {"poly": [[0, 0], [4, 0], [4, 2], [0, 2]]}
        result = ev.resolve_dir({"face_perp": "IW1", "face": "south"})
        # perp of east-going face: [0, -1]
        assert result == pytest.approx([0, -1])

    def test_literal_dir(self):
        ev = _make_evaluator()
        assert ev.resolve_dir([0.707, 0.707]) == pytest.approx([0.707, 0.707])

    def test_neg_dir(self):
        ev = _make_evaluator()
        assert ev.resolve_dir({"neg": "east"}) == [-1.0, 0.0]
        assert ev.resolve_dir({"neg": "north"}) == [0.0, -1.0]

    def test_perp_dir(self):
        ev = _make_evaluator()
        # perp of east [1,0] → [0, -1] (CW rotation)
        assert ev.resolve_dir({"perp": "east"}) == pytest.approx([0, -1])
        # perp of north [0,1] → [1, 0]
        assert ev.resolve_dir({"perp": "north"}) == pytest.approx([1, 0])

    def test_neg_segment(self):
        ev = _make_evaluator()
        # neg of W1→W2 (east) = west
        result = ev.resolve_dir({"neg": {"segment": ["W1", "W2"]}})
        assert result == pytest.approx([-1, 0])


# ---------------------------------------------------------------------------
# Length resolution tests
# ---------------------------------------------------------------------------

class TestResolveLength:
    """Test length expression resolution."""

    def test_literal(self):
        ev = _make_evaluator()
        assert ev.resolve_length(3.5) == 3.5
        assert ev.resolve_length(0) == 0.0

    def test_const_ref(self):
        ev = _make_evaluator()
        assert ev.resolve_length({"const": "WALL_6IN"}) == 0.5

    def test_unknown_const(self):
        ev = _make_evaluator()
        assert ev.resolve_length({"const": "NONEXISTENT"}) is None

    def test_none(self):
        ev = _make_evaluator()
        assert ev.resolve_length(None) is None

    def test_neg(self):
        ev = _make_evaluator()
        assert ev.resolve_length({"neg": 3.0}) == -3.0
        assert ev.resolve_length({"neg": {"const": "WALL_6IN"}}) == -0.5

    def test_add(self):
        ev = _make_evaluator()
        assert ev.resolve_length({"add": [1.0, 2.0, 3.0]}) == 6.0
        assert ev.resolve_length({"add": [{"const": "WALL_6IN"}, 1.0]}) == 1.5

    def test_sub(self):
        ev = _make_evaluator()
        assert ev.resolve_length({"sub": [5.0, 2.0]}) == 3.0
        assert ev.resolve_length({"sub": [{"const": "WALL_6IN"}, 0.25]}) == 0.25

    def test_mul(self):
        ev = _make_evaluator()
        assert ev.resolve_length({"mul": [3.0, 4.0]}) == 12.0
        assert ev.resolve_length({"mul": [{"const": "WALL_6IN"}, 2.0]}) == 1.0

    def test_dist(self):
        ev = _make_evaluator()
        assert ev.resolve_length({"dist": ["W1", "W2"]}) == pytest.approx(10.0)
        assert ev.resolve_length({"dist": [[0, 0], [3, 4]]}) == pytest.approx(5.0)

    def test_nested_arithmetic(self):
        ev = _make_evaluator()
        # IW7_OFFSET + WALL_4IN pattern
        result = ev.resolve_length({"add": [6.0, {"const": "WALL_4IN"}]})
        assert result == pytest.approx(6.0 + 1/3)


# ---------------------------------------------------------------------------
# Wall rect evaluation tests
# ---------------------------------------------------------------------------

class TestWallRect:
    """Test wall_rect formula evaluation."""

    def test_fixed_length_wall(self):
        ev = _make_evaluator()
        formula = {
            "type": "wall_rect",
            "anchor": [1, 0],
            "along": "east",
            "thickness_dir": "north",
            "thickness": 0.5,
            "end_mode": "fixed",
            "length": 4.0,
        }
        result = ev._evaluate_formula("test_wall", formula)
        assert result is not None
        poly = result["poly"]
        assert len(poly) == 4
        # SW = anchor = [1, 0]
        assert poly[0] == pytest.approx([1, 0])
        # SE = [5, 0]
        assert poly[1] == pytest.approx([5, 0])
        # NE = [5, 0.5]
        assert poly[2] == pytest.approx([5, 0.5])
        # NW = [1, 0.5]
        assert poly[3] == pytest.approx([1, 0.5])

    def test_const_thickness(self):
        ev = _make_evaluator()
        formula = {
            "type": "wall_rect",
            "anchor": [0, 0],
            "along": "east",
            "thickness_dir": "north",
            "thickness": {"const": "WALL_6IN"},
            "end_mode": "fixed",
            "length": 3.0,
        }
        result = ev._evaluate_formula("test_wall", formula)
        assert result is not None
        assert result["poly"][2] == pytest.approx([3, 0.5])

    def test_intersect_mode(self):
        ev = _make_evaluator()  # inner_poly = [0,0], [10,0], [10,10], [0,10]
        formula = {
            "type": "wall_rect",
            "anchor": [5, 3],
            "along": "east",
            "thickness_dir": "north",
            "thickness": 0.5,
            "end_mode": "intersect",
            "end_target": "inner_poly",
        }
        result = ev._evaluate_formula("test_wall", formula)
        assert result is not None
        # Wall should extend to inner poly edges
        poly = result["poly"]
        # Near end should be at x=0, far end at x=10
        assert poly[0][0] == pytest.approx(0)  # near (west wall)
        assert poly[1][0] == pytest.approx(10)  # far (east wall)


# ---------------------------------------------------------------------------
# Item rect evaluation tests
# ---------------------------------------------------------------------------

class TestItemRect:
    """Test item_rect formula evaluation."""

    def test_sw_anchor(self):
        ev = _make_evaluator()
        formula = {
            "type": "item_rect",
            "anchor": [2, 3],
            "along": "east",
            "across": "north",
            "width": 4.0,
            "depth": 2.0,
            "anchor_corner": "sw",
        }
        result = ev._evaluate_formula("bed", formula)
        assert result is not None
        poly = result["poly"]
        assert poly[0] == pytest.approx([2, 3])   # SW
        assert poly[1] == pytest.approx([6, 3])   # SE
        assert poly[2] == pytest.approx([6, 5])   # NE
        assert poly[3] == pytest.approx([2, 5])   # NW

    def test_ne_anchor(self):
        ev = _make_evaluator()
        formula = {
            "type": "item_rect",
            "anchor": [6, 5],
            "along": "east",
            "across": "north",
            "width": 4.0,
            "depth": 2.0,
            "anchor_corner": "ne",
        }
        result = ev._evaluate_formula("item", formula)
        assert result is not None
        poly = result["poly"]
        # NE is anchor, so SW = anchor - width*along - depth*across = [2, 3]
        assert poly[0] == pytest.approx([2, 3])

    def test_const_dimensions(self):
        ev = _make_evaluator(constants={"BED_W": 5.0, "BED_D": 6.5})
        formula = {
            "type": "item_rect",
            "anchor": [0, 0],
            "along": "east",
            "across": "north",
            "width": {"const": "BED_W"},
            "depth": {"const": "BED_D"},
            "anchor_corner": "sw",
        }
        result = ev._evaluate_formula("bed", formula)
        assert result is not None
        assert result["poly"][1] == pytest.approx([5, 0])
        assert result["poly"][2] == pytest.approx([5, 6.5])


# ---------------------------------------------------------------------------
# Item circle evaluation tests
# ---------------------------------------------------------------------------

class TestItemCircle:
    """Test item_circle formula evaluation."""

    def test_basic_circle(self):
        ev = _make_evaluator()
        formula = {
            "type": "item_circle",
            "center": [5, 5],
            "radius": 1.0,
        }
        result = ev._evaluate_formula("wh", formula)
        assert result is not None
        assert result["center"] == [5, 5]
        assert result["radius"] == 1.0
        assert len(result["poly"]) == 24
        # First point should be at angle=0 → (6, 5)
        assert result["poly"][0] == pytest.approx([6, 5])

    def test_const_radius(self):
        ev = _make_evaluator(constants={"WH_RADIUS": 0.75})
        formula = {
            "type": "item_circle",
            "center": [0, 0],
            "radius": {"const": "WH_RADIUS"},
        }
        result = ev._evaluate_formula("wh", formula)
        assert result is not None
        assert result["radius"] == 0.75


# ---------------------------------------------------------------------------
# Topological sort tests
# ---------------------------------------------------------------------------

class TestTopoSort:
    """Test topological sorting and cycle detection."""

    def test_independent_formulas(self):
        ev = _make_evaluator()
        ev.add_formula("A", "position", {"type": "item_rect",
                       "anchor": [0, 0], "along": "east", "across": "north",
                       "width": 1, "depth": 1, "anchor_corner": "sw"})
        ev.add_formula("B", "position", {"type": "item_rect",
                       "anchor": [5, 0], "along": "east", "across": "north",
                       "width": 1, "depth": 1, "anchor_corner": "sw"})
        ev.topo_sort()
        assert len(ev.eval_order) == 2

    def test_dependency_order(self):
        ev = _make_evaluator()
        # A depends on nothing
        ev.add_formula("A", "position", {"type": "item_rect",
                       "anchor": [0, 0], "along": "east", "across": "north",
                       "width": 1, "depth": 1, "anchor_corner": "sw"})
        # B depends on A
        ev.add_formula("B", "position", {"type": "item_rect",
                       "anchor": {"element": "A", "corner": "NE"},
                       "along": "east", "across": "north",
                       "width": 1, "depth": 1, "anchor_corner": "sw"})
        ev.topo_sort()
        order = [k[0] for k in ev.eval_order]
        assert order.index("A") < order.index("B")

    def test_chain_dependency(self):
        ev = _make_evaluator()
        ev.add_formula("A", "position", {"type": "item_rect",
                       "anchor": [0, 0], "along": "east", "across": "north",
                       "width": 1, "depth": 1, "anchor_corner": "sw"})
        ev.add_formula("B", "position", {"type": "item_rect",
                       "anchor": {"element": "A", "corner": "NE"},
                       "along": "east", "across": "north",
                       "width": 1, "depth": 1, "anchor_corner": "sw"})
        ev.add_formula("C", "position", {"type": "item_rect",
                       "anchor": {"element": "B", "corner": "NE"},
                       "along": "east", "across": "north",
                       "width": 1, "depth": 1, "anchor_corner": "sw"})
        ev.topo_sort()
        order = [k[0] for k in ev.eval_order]
        assert order.index("A") < order.index("B") < order.index("C")

    def test_cycle_detection(self):
        ev = _make_evaluator()
        # A depends on B, B depends on A
        ev.add_formula("A", "position", {"type": "item_rect",
                       "anchor": {"element": "B", "corner": 0},
                       "along": "east", "across": "north",
                       "width": 1, "depth": 1, "anchor_corner": "sw"})
        ev.add_formula("B", "position", {"type": "item_rect",
                       "anchor": {"element": "A", "corner": 0},
                       "along": "east", "across": "north",
                       "width": 1, "depth": 1, "anchor_corner": "sw"})
        with pytest.raises(CycleError):
            ev.topo_sort()


# ---------------------------------------------------------------------------
# Full evaluation tests
# ---------------------------------------------------------------------------

class TestEvaluateAll:
    """Test end-to-end formula evaluation."""

    def test_chained_items(self):
        ev = _make_evaluator()
        ev.add_formula("A", "position", {"type": "item_rect",
                       "anchor": [0, 0], "along": "east", "across": "north",
                       "width": 2, "depth": 1, "anchor_corner": "sw"})
        ev.add_formula("B", "position", {"type": "item_rect",
                       "anchor": {"element": "A", "corner": "SE"},
                       "along": "east", "across": "north",
                       "width": 3, "depth": 1, "anchor_corner": "sw"})
        ev.evaluate_all()
        assert "A" in ev.elements
        assert "B" in ev.elements
        # B's SW should be A's SE = [2, 0]
        assert ev.elements["B"]["poly"][0] == pytest.approx([2, 0])
        # B's SE = [5, 0]
        assert ev.elements["B"]["poly"][1] == pytest.approx([5, 0])

    def test_locked_formula(self):
        ev = _make_evaluator()
        locked_value = {"poly": [[10, 10], [12, 10], [12, 12], [10, 12]],
                        "bbox": {"w": 10, "s": 10, "e": 12, "n": 12}}
        ev.add_formula("A", "position", {"type": "item_rect",
                       "anchor": [0, 0], "along": "east", "across": "north",
                       "width": 2, "depth": 2, "anchor_corner": "sw"},
                       locked=True, locked_value=locked_value)
        ev.evaluate_all()
        # Should use locked value, not evaluate the formula
        assert ev.elements["A"]["poly"][0] == [10, 10]

    def test_wall_with_element_reference(self):
        ev = _make_evaluator()
        # First: create a reference wall
        ev.add_formula("IW1", "position", {"type": "wall_rect",
                       "anchor": [0, 3], "along": "east",
                       "thickness_dir": "north", "thickness": 0.5,
                       "end_mode": "fixed", "length": 10})
        # Second: item anchored to IW1's NE corner
        ev.add_formula("bed", "position", {"type": "item_rect",
                       "anchor": {"element": "IW1", "corner": "NE"},
                       "along": "east", "across": "north",
                       "width": 5, "depth": 6.5, "anchor_corner": "sw"})
        ev.evaluate_all()
        # IW1 NE = [10, 3.5]
        assert ev.elements["bed"]["poly"][0] == pytest.approx([10, 3.5])


# ---------------------------------------------------------------------------
# Dependency extraction tests
# ---------------------------------------------------------------------------

class TestExtractDeps:
    """Test static dependency extraction from formula JSON."""

    def test_element_ref(self):
        deps = extract_deps({"element": "IW1", "corner": "SW"})
        assert ("element", "IW1") in deps

    def test_face_mid_ref(self):
        deps = extract_deps({"face_mid": "IW2", "face": "north"})
        assert ("element", "IW2") in deps

    def test_const_ref(self):
        deps = extract_deps({"const": "WALL_6IN"})
        assert ("constant", "WALL_6IN") in deps

    def test_point_ref(self):
        deps = extract_deps("W9")
        assert ("point", "W9") in deps

    def test_nested_deps(self):
        formula = {
            "type": "wall_rect",
            "anchor": {"offset": "W9",
                       "dir": {"segment_perp": ["W9", "W10"]},
                       "dist": {"const": "IW1_OFFSET"}},
            "along": {"segment": ["W9", "W10"]},
            "thickness_dir": "north",
            "thickness": {"const": "WALL_6IN"},
            "end_mode": "fixed",
            "length": 5.0,
        }
        deps = extract_deps(formula)
        assert ("point", "W9") in deps
        assert ("point", "W10") in deps
        assert ("constant", "IW1_OFFSET") in deps
        assert ("constant", "WALL_6IN") in deps

    def test_element_deps_in_formula(self):
        formula = {
            "type": "item_rect",
            "anchor": {"element": "IW1", "corner": "NE"},
            "along": {"face_along": "IW2", "face": "south"},
            "across": {"face_perp": "IW2", "face": "south"},
            "width": {"const": "BED_W"},
            "depth": {"const": "BED_D"},
            "anchor_corner": "sw",
        }
        deps = extract_deps(formula)
        assert ("element", "IW1") in deps
        assert ("element", "IW2") in deps
        assert ("constant", "BED_W") in deps
        assert ("constant", "BED_D") in deps

    def test_empty_formula(self):
        deps = extract_deps({})
        assert len(deps) == 0

    def test_literal_has_no_deps(self):
        deps = extract_deps(3.5)
        assert len(deps) == 0


# ---------------------------------------------------------------------------
# Dependency query tests
# ---------------------------------------------------------------------------

class TestDependencyQueries:
    """Test get_dependents and get_dependencies."""

    def test_get_dependents(self):
        ev = _make_evaluator()
        ev.add_formula("A", "position", {"type": "item_rect",
                       "anchor": [0, 0], "along": "east", "across": "north",
                       "width": 1, "depth": 1, "anchor_corner": "sw"})
        ev.add_formula("B", "position", {"type": "item_rect",
                       "anchor": {"element": "A", "corner": "NE"},
                       "along": "east", "across": "north",
                       "width": 1, "depth": 1, "anchor_corner": "sw"})
        ev.add_formula("C", "position", {"type": "item_rect",
                       "anchor": {"element": "A", "corner": "SE"},
                       "along": "east", "across": "north",
                       "width": 1, "depth": 1, "anchor_corner": "sw"})
        dependents = ev.get_dependents("A")
        assert "B" in dependents
        assert "C" in dependents
        assert "A" not in dependents

    def test_get_dependencies(self):
        ev = _make_evaluator()
        ev.add_formula("A", "position", {"type": "item_rect",
                       "anchor": [0, 0], "along": "east", "across": "north",
                       "width": 1, "depth": 1, "anchor_corner": "sw"})
        ev.add_formula("B", "position", {"type": "item_rect",
                       "anchor": {"element": "A", "corner": "NE"},
                       "along": "east", "across": "north",
                       "width": 1, "depth": 1, "anchor_corner": "sw"})
        ev.add_formula("C", "position", {"type": "item_rect",
                       "anchor": {"element": "B", "corner": "NE"},
                       "along": "east", "across": "north",
                       "width": 1, "depth": 1, "anchor_corner": "sw"})
        deps = ev.get_dependencies("C")
        assert "B" in deps
        assert "A" in deps  # transitive

    def test_transitive_dependents(self):
        ev = _make_evaluator()
        ev.add_formula("A", "position", {"type": "item_rect",
                       "anchor": [0, 0], "along": "east", "across": "north",
                       "width": 1, "depth": 1, "anchor_corner": "sw"})
        ev.add_formula("B", "position", {"type": "item_rect",
                       "anchor": {"element": "A", "corner": 1},
                       "along": "east", "across": "north",
                       "width": 1, "depth": 1, "anchor_corner": "sw"})
        ev.add_formula("C", "position", {"type": "item_rect",
                       "anchor": {"element": "B", "corner": 1},
                       "along": "east", "across": "north",
                       "width": 1, "depth": 1, "anchor_corner": "sw"})
        dependents = ev.get_dependents("A")
        assert "B" in dependents
        assert "C" in dependents  # transitive


# ---------------------------------------------------------------------------
# Line-polygon intersection helper tests
# ---------------------------------------------------------------------------

class TestLinePolyIntersections:
    """Test the line-polygon intersection utility."""

    def test_horizontal_through_square(self):
        poly = [[0, 0], [10, 0], [10, 10], [0, 10]]
        result = _line_poly_intersections([5, 5], [1, 0], poly)
        assert len(result) == 1
        assert result[0] == pytest.approx([10, 5])

    def test_both_directions(self):
        poly = [[0, 0], [10, 0], [10, 10], [0, 10]]
        # East hits east wall
        east = _line_poly_intersections([5, 5], [1, 0], poly)
        assert len(east) == 1
        assert east[0][0] == pytest.approx(10)
        # West hits west wall
        west = _line_poly_intersections([5, 5], [-1, 0], poly)
        assert len(west) == 1
        assert west[0][0] == pytest.approx(0)

    def test_no_intersection(self):
        poly = [[0, 0], [10, 0], [10, 10], [0, 10]]
        # Shooting away from polygon
        result = _line_poly_intersections([15, 5], [1, 0], poly)
        assert len(result) == 0


# ---------------------------------------------------------------------------
# Database formula CRUD tests
# ---------------------------------------------------------------------------

class TestFormulaCRUD:
    """Test formula database operations."""

    def test_upsert_and_get_formula(self, fresh_db):
        from app.database import upsert_formula, get_element_formulas
        from app.database import get_db
        # First create an element to reference
        with get_db(fresh_db) as conn:
            conn.execute(
                "INSERT INTO elements (type, name, properties) VALUES (?, ?, ?)",
                ("wall", "TEST_WALL", "{}"),
            )
        formula = {"type": "wall_rect", "anchor": [0, 0], "along": "east",
                   "thickness_dir": "north", "thickness": 0.5,
                   "end_mode": "fixed", "length": 5.0}
        result = upsert_formula("TEST_WALL", "position", formula, db_path=fresh_db)
        assert result is not None
        assert result["element_name"] == "TEST_WALL"
        assert result["param_name"] == "position"
        parsed = json.loads(result["formula_json"])
        assert parsed["type"] == "wall_rect"

        # Get it back
        formulas = get_element_formulas("TEST_WALL", db_path=fresh_db)
        assert len(formulas) == 1
        assert formulas[0]["param_name"] == "position"

    def test_upsert_updates_existing(self, fresh_db):
        from app.database import upsert_formula, get_element_formulas, get_db
        with get_db(fresh_db) as conn:
            conn.execute(
                "INSERT INTO elements (type, name, properties) VALUES (?, ?, ?)",
                ("wall", "TW2", "{}"),
            )
        formula1 = {"type": "wall_rect", "anchor": [0, 0], "along": "east",
                    "thickness_dir": "north", "thickness": 0.5,
                    "end_mode": "fixed", "length": 5.0}
        upsert_formula("TW2", "position", formula1, db_path=fresh_db)
        formula2 = {"type": "wall_rect", "anchor": [1, 1], "along": "north",
                    "thickness_dir": "east", "thickness": 0.5,
                    "end_mode": "fixed", "length": 8.0}
        upsert_formula("TW2", "position", formula2, db_path=fresh_db)
        formulas = get_element_formulas("TW2", db_path=fresh_db)
        assert len(formulas) == 1
        parsed = json.loads(formulas[0]["formula_json"])
        assert parsed["length"] == 8.0

    def test_delete_formula(self, fresh_db):
        from app.database import upsert_formula, delete_formula, get_element_formulas, get_db
        with get_db(fresh_db) as conn:
            conn.execute(
                "INSERT INTO elements (type, name, properties) VALUES (?, ?, ?)",
                ("wall", "TW3", "{}"),
            )
        formula = {"type": "item_rect", "anchor": [0, 0], "along": "east",
                   "across": "north", "width": 1, "depth": 1,
                   "anchor_corner": "sw"}
        upsert_formula("TW3", "position", formula, db_path=fresh_db)
        assert delete_formula("TW3", "position", db_path=fresh_db)
        assert len(get_element_formulas("TW3", db_path=fresh_db)) == 0

    def test_formula_lock(self, fresh_db):
        from app.database import upsert_formula, set_formula_lock, get_element_formulas, get_db
        with get_db(fresh_db) as conn:
            conn.execute(
                "INSERT INTO elements (type, name, properties) VALUES (?, ?, ?)",
                ("wall", "TW4", "{}"),
            )
        formula = {"type": "item_rect", "anchor": [0, 0], "along": "east",
                   "across": "north", "width": 1, "depth": 1,
                   "anchor_corner": "sw"}
        upsert_formula("TW4", "position", formula, db_path=fresh_db)
        locked_val = {"poly": [[0, 0], [1, 0], [1, 1], [0, 1]]}
        assert set_formula_lock("TW4", "position", True, locked_val, db_path=fresh_db)
        formulas = get_element_formulas("TW4", db_path=fresh_db)
        assert formulas[0]["locked"] == 1
        assert json.loads(formulas[0]["locked_value"]) == locked_val

    def test_get_all_formulas(self, fresh_db):
        from app.database import upsert_formula, get_all_formulas, get_db
        with get_db(fresh_db) as conn:
            conn.execute(
                "INSERT INTO elements (type, name, properties) VALUES (?, ?, ?)",
                ("wall", "TW5", "{}"),
            )
            conn.execute(
                "INSERT INTO elements (type, name, properties) VALUES (?, ?, ?)",
                ("furniture", "TW6", "{}"),
            )
        f1 = {"type": "wall_rect", "anchor": [0, 0], "along": "east",
              "thickness_dir": "north", "thickness": 0.5,
              "end_mode": "fixed", "length": 5.0}
        f2 = {"type": "item_rect", "anchor": [0, 0], "along": "east",
              "across": "north", "width": 1, "depth": 1,
              "anchor_corner": "sw"}
        upsert_formula("TW5", "position", f1, db_path=fresh_db)
        upsert_formula("TW6", "position", f2, db_path=fresh_db)
        all_formulas = get_all_formulas(db_path=fresh_db)
        assert len(all_formulas) >= 2

    def test_rebuild_formula_deps(self, fresh_db):
        from app.database import rebuild_formula_deps, get_formula_deps
        deps = [("element", "IW1"), ("constant", "WALL_6IN"), ("point", "W9")]
        rebuild_formula_deps("TEST", "position", deps, db_path=fresh_db)
        stored = get_formula_deps("TEST", "position", db_path=fresh_db)
        dep_set = {(d["dep_type"], d["dep_name"]) for d in stored}
        assert ("element", "IW1") in dep_set
        assert ("constant", "WALL_6IN") in dep_set
        assert ("point", "W9") in dep_set

    def test_get_dependents(self, fresh_db):
        from app.database import rebuild_formula_deps, get_dependents as db_get_dependents
        rebuild_formula_deps("bed", "position",
                             [("element", "IW1")], db_path=fresh_db)
        rebuild_formula_deps("dresser", "position",
                             [("element", "IW1")], db_path=fresh_db)
        result = db_get_dependents("IW1", dep_type="element", db_path=fresh_db)
        names = {r["element_name"] for r in result}
        assert "bed" in names
        assert "dresser" in names


# ---------------------------------------------------------------------------
# Evaluator loading from DB test
# ---------------------------------------------------------------------------

class TestEvaluatorDB:
    """Test evaluator loading formulas from database."""

    def test_load_and_evaluate(self, fresh_db):
        from app.database import upsert_formula, get_db
        with get_db(fresh_db) as conn:
            conn.execute(
                "INSERT INTO elements (type, name, properties) VALUES (?, ?, ?)",
                ("furniture", "TEST_BED", "{}"),
            )
        formula = {"type": "item_rect", "anchor": [0, 0], "along": "east",
                   "across": "north", "width": 5, "depth": 6.5,
                   "anchor_corner": "sw"}
        upsert_formula("TEST_BED", "position", formula, db_path=fresh_db)

        ev = FormulaEvaluator(
            constants={},
            base_points={},
            inner_poly=[],
            radii={},
        )
        ev.load_formulas_from_db(db_path=fresh_db)
        ev.evaluate_all()
        assert "TEST_BED" in ev.elements
        assert ev.elements["TEST_BED"]["poly"][0] == pytest.approx([0, 0])
        assert ev.elements["TEST_BED"]["poly"][1] == pytest.approx([5, 0])


# ---------------------------------------------------------------------------
# Variant item constants in DB tests
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# IW formula verification tests (Phase 12c)
# ---------------------------------------------------------------------------

class TestIWFormulas:
    """Verify formula-computed IW walls match procedural layout output."""

    @pytest.fixture
    def geom_and_evaluator(self, fresh_db):
        """Compute procedural geometry and build a formula evaluator."""
        from app.database import get_constants_dict
        from app.engine import compute_geometry
        from app.evaluator import FormulaEvaluator, get_iw_formulas

        constants = get_constants_dict(fresh_db)
        geom = compute_geometry(constants, db_path=fresh_db)

        base_points = geom.get("points", {})
        inner = [(p[0], p[1]) for p in (geom.get("inner_poly") or [])]
        radii = geom.get("radii", {})

        ev = FormulaEvaluator(constants, base_points, inner, radii)
        formulas = get_iw_formulas()
        # Add formulas respecting dependency order
        for name, formula in formulas.items():
            ev.add_formula(name, "position", formula)
        ev.topo_sort()
        ev.evaluate_all()

        return geom, ev

    def _assert_poly_match(self, geom, ev, iw_name, atol=0.01):
        """Assert formula polygon matches procedural polygon within tolerance."""
        proc = geom["interior_walls"][iw_name]
        proc_poly = proc["poly"]
        form = ev.elements.get(iw_name)
        assert form is not None, f"Formula did not produce result for {iw_name}"
        form_poly = form["poly"]
        assert len(form_poly) == len(proc_poly), \
            f"{iw_name}: poly length mismatch {len(form_poly)} vs {len(proc_poly)}"
        for i, (fp, pp) in enumerate(zip(form_poly, proc_poly)):
            assert fp[0] == pytest.approx(pp[0], abs=atol), \
                f"{iw_name} corner {i} E: formula={fp[0]:.4f} proc={pp[0]:.4f}"
            assert fp[1] == pytest.approx(pp[1], abs=atol), \
                f"{iw_name} corner {i} N: formula={fp[1]:.4f} proc={pp[1]:.4f}"

    def test_iw1(self, geom_and_evaluator):
        geom, ev = geom_and_evaluator
        self._assert_poly_match(geom, ev, "IW1")

    def test_iw2(self, geom_and_evaluator):
        geom, ev = geom_and_evaluator
        self._assert_poly_match(geom, ev, "IW2")

    def test_iw2o(self, geom_and_evaluator):
        geom, ev = geom_and_evaluator
        self._assert_poly_match(geom, ev, "IW2O")

    def test_iw2s(self, geom_and_evaluator):
        geom, ev = geom_and_evaluator
        self._assert_poly_match(geom, ev, "IW2S")

    def test_iw3(self, geom_and_evaluator):
        geom, ev = geom_and_evaluator
        self._assert_poly_match(geom, ev, "IW3")

    def test_iw4(self, geom_and_evaluator):
        geom, ev = geom_and_evaluator
        self._assert_poly_match(geom, ev, "IW4")

    def test_iw5(self, geom_and_evaluator):
        geom, ev = geom_and_evaluator
        self._assert_poly_match(geom, ev, "IW5")

    def test_iw6(self, geom_and_evaluator):
        geom, ev = geom_and_evaluator
        self._assert_poly_match(geom, ev, "IW6")

    def test_iw7(self, geom_and_evaluator):
        geom, ev = geom_and_evaluator
        self._assert_poly_match(geom, ev, "IW7")

    def test_iw8(self, geom_and_evaluator):
        geom, ev = geom_and_evaluator
        self._assert_poly_match(geom, ev, "IW8")

    def test_iw9(self, geom_and_evaluator):
        geom, ev = geom_and_evaluator
        self._assert_poly_match(geom, ev, "IW9")

    def test_iw11(self, geom_and_evaluator):
        geom, ev = geom_and_evaluator
        self._assert_poly_match(geom, ev, "IW11")

    def test_iw12(self, geom_and_evaluator):
        geom, ev = geom_and_evaluator
        self._assert_poly_match(geom, ev, "IW12")

    def test_all_13_walls_evaluated(self, geom_and_evaluator):
        _, ev = geom_and_evaluator
        expected = {"IW1", "IW2", "IW2O", "IW2S", "IW3", "IW4",
                    "IW5", "IW6", "IW7", "IW8", "IW9", "IW11", "IW12"}
        assert expected <= set(ev.elements.keys())


class TestLayoutItemFormulas:
    """Verify formula-computed layout items match procedural output."""

    @pytest.fixture
    def geom_and_evaluator(self, fresh_db):
        """Compute procedural geometry and build evaluator with IW + layout formulas."""
        from app.database import get_constants_dict
        from app.engine import compute_geometry
        from app.evaluator import FormulaEvaluator, get_iw_formulas, get_layout_item_formulas

        constants = get_constants_dict(fresh_db)
        geom = compute_geometry(constants, db_path=fresh_db)

        base_points = geom.get("points", {})
        inner = [(p[0], p[1]) for p in (geom.get("inner_poly") or [])]
        radii = geom.get("radii", {})

        ev = FormulaEvaluator(constants, base_points, inner, radii)
        # Add IW formulas (layout items depend on IW walls)
        for name, formula in get_iw_formulas().items():
            ev.add_formula(name, "position", formula)
        # Add layout item formulas
        for name, formula in get_layout_item_formulas().items():
            ev.add_formula(name, "position", formula)
        ev.topo_sort()
        ev.evaluate_all()

        return geom, ev

    def _assert_item_match(self, geom, ev, item_name, section, atol=0.01):
        """Assert formula polygon matches procedural polygon."""
        proc = geom[section][item_name.lower()]
        proc_poly = proc["poly"]
        form = ev.elements.get(item_name)
        assert form is not None, f"Formula did not produce result for {item_name}"
        form_poly = form["poly"]
        assert len(form_poly) == len(proc_poly), \
            f"{item_name}: poly length mismatch {len(form_poly)} vs {len(proc_poly)}"
        for i, (fp, pp) in enumerate(zip(form_poly, proc_poly)):
            assert fp[0] == pytest.approx(pp[0], abs=atol), \
                f"{item_name} corner {i} E: formula={fp[0]:.4f} proc={pp[0]:.4f}"
            assert fp[1] == pytest.approx(pp[1], abs=atol), \
                f"{item_name} corner {i} N: formula={fp[1]:.4f} proc={pp[1]:.4f}"

    def test_dryer(self, geom_and_evaluator):
        geom, ev = geom_and_evaluator
        self._assert_item_match(geom, ev, "dryer", "appliances")

    def test_washer(self, geom_and_evaluator):
        geom, ev = geom_and_evaluator
        self._assert_item_match(geom, ev, "washer", "appliances")

    def test_counter(self, geom_and_evaluator):
        geom, ev = geom_and_evaluator
        self._assert_item_match(geom, ev, "counter", "appliances")

    def test_dresser(self, geom_and_evaluator):
        geom, ev = geom_and_evaluator
        self._assert_item_match(geom, ev, "dresser", "furniture")

    def test_shelves(self, geom_and_evaluator):
        geom, ev = geom_and_evaluator
        self._assert_item_match(geom, ev, "shelves", "furniture")


class TestOuterOpeningFormulas:
    """Verify formula-computed outer openings match procedural output."""

    @pytest.fixture
    def geom_and_evaluator(self, fresh_db):
        from app.database import get_constants_dict
        from app.engine import compute_geometry
        from app.evaluator import (FormulaEvaluator, get_iw_formulas,
                                   get_layout_item_formulas,
                                   get_outer_opening_formulas)

        constants = get_constants_dict(fresh_db)
        geom = compute_geometry(constants, db_path=fresh_db)

        base_points = geom.get("points", {})
        inner = [(p[0], p[1]) for p in (geom.get("inner_poly") or [])]
        radii = geom.get("radii", {})

        ev = FormulaEvaluator(constants, base_points, inner, radii)
        for name, formula in get_iw_formulas().items():
            ev.add_formula(name, "position", formula)
        for name, formula in get_layout_item_formulas().items():
            ev.add_formula(name, "position", formula)
        for name, formula in get_outer_opening_formulas().items():
            ev.add_formula(name, "position", formula)
        ev.topo_sort()
        ev.evaluate_all()

        return geom, ev

    def _assert_opening_match(self, geom, ev, name, atol=0.01):
        proc_list = geom.get("outer_openings", [])
        proc = next((o for o in proc_list if o["name"] == name), None)
        assert proc is not None, f"No procedural opening {name}"
        form = ev.elements.get(name)
        assert form is not None, f"Formula did not produce result for {name}"
        proc_poly = proc["poly"]
        form_poly = form["poly"]
        assert len(form_poly) == len(proc_poly), \
            f"{name}: poly length mismatch"
        for i, (fp, pp) in enumerate(zip(form_poly, proc_poly)):
            assert fp[0] == pytest.approx(pp[0], abs=atol), \
                f"{name} corner {i} E: formula={fp[0]:.4f} proc={pp[0]:.4f}"
            assert fp[1] == pytest.approx(pp[1], abs=atol), \
                f"{name} corner {i} N: formula={fp[1]:.4f} proc={pp[1]:.4f}"

    def test_o1(self, geom_and_evaluator):
        self._assert_opening_match(*geom_and_evaluator, "O1")

    def test_o2(self, geom_and_evaluator):
        self._assert_opening_match(*geom_and_evaluator, "O2")

    def test_o3(self, geom_and_evaluator):
        self._assert_opening_match(*geom_and_evaluator, "O3")

    def test_o4(self, geom_and_evaluator):
        self._assert_opening_match(*geom_and_evaluator, "O4")

    def test_o5(self, geom_and_evaluator):
        self._assert_opening_match(*geom_and_evaluator, "O5")

    def test_o6(self, geom_and_evaluator):
        self._assert_opening_match(*geom_and_evaluator, "O6")

    def test_o7(self, geom_and_evaluator):
        self._assert_opening_match(*geom_and_evaluator, "O7")

    def test_o8(self, geom_and_evaluator):
        self._assert_opening_match(*geom_and_evaluator, "O8")

    def test_o8a(self, geom_and_evaluator):
        self._assert_opening_match(*geom_and_evaluator, "O8a")

    def test_o9(self, geom_and_evaluator):
        self._assert_opening_match(*geom_and_evaluator, "O9")

    def test_o10(self, geom_and_evaluator):
        self._assert_opening_match(*geom_and_evaluator, "O10")

    def test_o11(self, geom_and_evaluator):
        self._assert_opening_match(*geom_and_evaluator, "O11")


class TestRoughOpeningFormulas:
    """Verify formula-computed rough openings match procedural output."""

    @pytest.fixture
    def geom_and_evaluator(self, fresh_db):
        from app.database import get_constants_dict
        from app.engine import compute_geometry
        from app.evaluator import (FormulaEvaluator, get_iw_formulas,
                                   get_layout_item_formulas,
                                   get_rough_opening_formulas)

        constants = get_constants_dict(fresh_db)
        geom = compute_geometry(constants, db_path=fresh_db)

        base_points = geom.get("points", {})
        inner = [(p[0], p[1]) for p in (geom.get("inner_poly") or [])]
        radii = geom.get("radii", {})

        ev = FormulaEvaluator(constants, base_points, inner, radii)
        for name, formula in get_iw_formulas().items():
            ev.add_formula(name, "position", formula)
        for name, formula in get_layout_item_formulas().items():
            ev.add_formula(name, "position", formula)
        for name, formula in get_rough_opening_formulas().items():
            ev.add_formula(name, "position", formula)
        ev.topo_sort()
        ev.evaluate_all()

        return geom, ev

    def _assert_ro_match(self, geom, ev, name, atol=0.01):
        proc_list = geom.get("rough_openings", [])
        proc = next((o for o in proc_list if o["name"] == name), None)
        assert proc is not None, f"No procedural rough opening {name}"
        form = ev.elements.get(name)
        assert form is not None, f"Formula did not produce result for {name}"
        proc_poly = proc["poly"]
        form_poly = form["poly"]
        assert len(form_poly) == len(proc_poly), \
            f"{name}: poly length mismatch"
        for i, (fp, pp) in enumerate(zip(form_poly, proc_poly)):
            assert fp[0] == pytest.approx(pp[0], abs=atol), \
                f"{name} corner {i} E: formula={fp[0]:.4f} proc={pp[0]:.4f}"
            assert fp[1] == pytest.approx(pp[1], abs=atol), \
                f"{name} corner {i} N: formula={fp[1]:.4f} proc={pp[1]:.4f}"

    def test_ro1(self, geom_and_evaluator):
        self._assert_ro_match(*geom_and_evaluator, "RO1")

    def test_ro2(self, geom_and_evaluator):
        self._assert_ro_match(*geom_and_evaluator, "RO2")

    def test_ro3(self, geom_and_evaluator):
        self._assert_ro_match(*geom_and_evaluator, "RO3")

    def test_ro4(self, geom_and_evaluator):
        self._assert_ro_match(*geom_and_evaluator, "RO4")

    def test_ro5(self, geom_and_evaluator):
        self._assert_ro_match(*geom_and_evaluator, "RO5")

    def test_ro6(self, geom_and_evaluator):
        self._assert_ro_match(*geom_and_evaluator, "RO6")

    def test_ro7(self, geom_and_evaluator):
        self._assert_ro_match(*geom_and_evaluator, "RO7")


class TestVariantItemConstants:
    """Test that the 24 variant item constants are seeded in the DB."""

    def test_all_24_constants_exist(self, fresh_db):
        from app.database import get_constants_dict
        constants = get_constants_dict(fresh_db)
        expected = [
            "HAMPER_W", "HAMPER_D", "MINIK_APPL_W", "MINIK_APPL_D",
            "MICROWAVE_W", "MICROWAVE_D", "COFFEE_W", "COFFEE_D",
            "COOKTOP_W", "COOKTOP_D", "TOASTER_W", "TOASTER_D",
            "DINING_CHAIR_W", "DINING_CHAIR_D", "DINING_TBL_BASE",
            "DINING_TBL_H", "DAYBED_W", "DAYBED_D", "WORK_CTR_W",
            "WORK_CTR_D", "STD_FRIDGE_W", "STD_FRIDGE_D",
            "SOFA_FULL_W", "SOFA_FULL_D",
        ]
        for name in expected:
            assert name in constants, f"Missing constant: {name}"

    def test_hamper_w_value(self, fresh_db):
        from app.database import get_constants_dict
        constants = get_constants_dict(fresh_db)
        assert constants["HAMPER_W"] == pytest.approx(31.5 / 12.0)

    def test_constants_match_variants_py(self, fresh_db):
        """All 24 DB constants match the hardcoded variants.py values."""
        from app.database import get_constants_dict
        from app.variants import (
            HAMPER_W, HAMPER_D, MINIK_APPL_W, MINIK_APPL_D,
            MICROWAVE_W, MICROWAVE_D, COFFEE_W, COFFEE_D,
            COOKTOP_W, COOKTOP_D, TOASTER_W, TOASTER_D,
            DINING_CHAIR_W, DINING_CHAIR_D, DINING_TBL_BASE,
            DINING_TBL_H, DAYBED_W, DAYBED_D, WORK_CTR_W,
            WORK_CTR_D, STD_FRIDGE_W, STD_FRIDGE_D,
            SOFA_FULL_W, SOFA_FULL_D,
        )
        constants = get_constants_dict(fresh_db)
        expected = {
            "HAMPER_W": HAMPER_W, "HAMPER_D": HAMPER_D,
            "MINIK_APPL_W": MINIK_APPL_W, "MINIK_APPL_D": MINIK_APPL_D,
            "MICROWAVE_W": MICROWAVE_W, "MICROWAVE_D": MICROWAVE_D,
            "COFFEE_W": COFFEE_W, "COFFEE_D": COFFEE_D,
            "COOKTOP_W": COOKTOP_W, "COOKTOP_D": COOKTOP_D,
            "TOASTER_W": TOASTER_W, "TOASTER_D": TOASTER_D,
            "DINING_CHAIR_W": DINING_CHAIR_W, "DINING_CHAIR_D": DINING_CHAIR_D,
            "DINING_TBL_BASE": DINING_TBL_BASE, "DINING_TBL_H": DINING_TBL_H,
            "DAYBED_W": DAYBED_W, "DAYBED_D": DAYBED_D,
            "WORK_CTR_W": WORK_CTR_W, "WORK_CTR_D": WORK_CTR_D,
            "STD_FRIDGE_W": STD_FRIDGE_W, "STD_FRIDGE_D": STD_FRIDGE_D,
            "SOFA_FULL_W": SOFA_FULL_W, "SOFA_FULL_D": SOFA_FULL_D,
        }
        for name, expected_val in expected.items():
            assert constants[name] == pytest.approx(expected_val), \
                f"{name}: DB={constants[name]} vs py={expected_val}"

    def test_constant_categories(self, fresh_db):
        from app.database import get_all_constants
        all_consts = get_all_constants(fresh_db)
        by_name = {c["name"]: c for c in all_consts}
        assert by_name["HAMPER_W"]["category"] == "furniture"
        assert by_name["MICROWAVE_W"]["category"] == "appliance"
        assert by_name["STD_FRIDGE_W"]["category"] == "appliance"
        assert by_name["SOFA_FULL_W"]["category"] == "furniture"


# ---------------------------------------------------------------------------
# DB-seeded IW formula tests
# ---------------------------------------------------------------------------

class TestSeededIWFormulas:
    """Verify IW formulas are seeded into element_formulas + formula_deps."""

    def test_all_formulas_seeded(self, fresh_db):
        from app.database import get_all_formulas
        formulas = get_all_formulas(db_path=fresh_db)
        names = {f["element_name"] for f in formulas}
        # 13 IW walls + 5 layout items + 12 outer openings + 7 rough openings = 37
        expected_iw = {"IW1", "IW2", "IW2S", "IW2O", "IW3", "IW4", "IW5",
                       "IW6", "IW7", "IW8", "IW9", "IW11", "IW12"}
        expected_items = {"dryer", "washer", "counter", "dresser", "shelves"}
        expected_openings = {"O1", "O2", "O3", "O4", "O5", "O6", "O7",
                             "O8", "O8a", "O9", "O10", "O11"}
        expected_ro = {"RO1", "RO2", "RO3", "RO4", "RO5", "RO6", "RO7"}
        assert expected_iw.issubset(names)
        assert expected_items.issubset(names)
        assert expected_openings.issubset(names)
        assert expected_ro.issubset(names)
        assert len(names) >= 37

    def test_iw1_has_deps(self, fresh_db):
        from app.database import get_formula_deps
        deps = get_formula_deps("IW1", "position", db_path=fresh_db)
        dep_set = {(d["dep_type"], d["dep_name"]) for d in deps}
        # IW1 depends on W-series points and WALL_6IN constant
        assert ("constant", "WALL_6IN") in dep_set
        assert ("point", "W9") in dep_set or ("point", "W10") in dep_set

    def test_iw4_depends_on_iw11_and_iw12(self, fresh_db):
        from app.database import get_formula_deps
        deps = get_formula_deps("IW4", "position", db_path=fresh_db)
        dep_set = {(d["dep_type"], d["dep_name"]) for d in deps}
        assert ("element", "IW11") in dep_set
        assert ("element", "IW12") in dep_set

    def test_iw12_depends_on_iw11_not_iw4(self, fresh_db):
        from app.database import get_formula_deps
        deps = get_formula_deps("IW12", "position", db_path=fresh_db)
        dep_set = {(d["dep_type"], d["dep_name"]) for d in deps}
        assert ("element", "IW11") in dep_set
        # IW12 inlines IW4.SW computation, so no direct dep on IW4
        assert ("element", "IW4") not in dep_set

    def test_seeding_is_idempotent(self, fresh_db):
        from app.database import get_all_formulas, get_db, seed_iw_formulas
        before = len(get_all_formulas(db_path=fresh_db))
        with get_db(fresh_db) as conn:
            seed_iw_formulas(conn)
        after = len(get_all_formulas(db_path=fresh_db))
        assert before == after


# ---------------------------------------------------------------------------
# Variant item formula tests (Phase 12e)
# ---------------------------------------------------------------------------

class TestVariantItemFormulas:
    """Verify formula-computed variant items match procedural output."""

    TOLERANCE = 0.02  # 0.02 ft = 0.24" — slightly looser for complex items

    @pytest.fixture(params=["standard", "minik", "daybed"])
    def variant_geom(self, request, fresh_db):
        """Compute procedural geometry and build evaluator for each variant."""
        variant = request.param
        from app.database import get_constants_dict
        from app.engine import compute_geometry
        from app.evaluator import (
            FormulaEvaluator, get_iw_formulas, get_layout_item_formulas,
            get_variant_item_formulas,
        )

        constants = get_constants_dict(fresh_db)
        geom = compute_geometry(constants, variant=variant, db_path=fresh_db)

        base_points = geom.get("points", {})
        inner = [(p[0], p[1]) for p in (geom.get("inner_poly") or [])]
        radii = geom.get("radii", {})

        ev = FormulaEvaluator(constants, base_points, inner, radii)
        # Add IW + layout formulas (variant items depend on them)
        for name, formula in get_iw_formulas().items():
            ev.add_formula(name, "position", formula)
        for name, formula in get_layout_item_formulas().items():
            ev.add_formula(name, "position", formula)
        # Add variant item formulas matching this variant
        for name, var, formula in get_variant_item_formulas():
            if var is None or var == variant:
                ev.add_formula(name, "position", formula)
        ev.topo_sort()
        ev.evaluate_all()

        return variant, geom, ev

    def _assert_poly_match(self, proc_poly, form_poly, name, atol=None):
        if atol is None:
            atol = self.TOLERANCE
        assert len(form_poly) == len(proc_poly), \
            f"{name}: poly length mismatch {len(form_poly)} vs {len(proc_poly)}"
        for i, (fp, pp) in enumerate(zip(form_poly, proc_poly)):
            fp_e = fp[0] if isinstance(fp, (list, tuple)) else fp
            fp_n = fp[1] if isinstance(fp, (list, tuple)) else fp
            pp_e = pp[0] if isinstance(pp, (list, tuple)) else pp
            pp_n = pp[1] if isinstance(pp, (list, tuple)) else pp
            assert fp_e == pytest.approx(pp_e, abs=atol), \
                f"{name} vertex {i} E: formula={fp_e:.6f} proc={pp_e:.6f}"
            assert fp_n == pytest.approx(pp_n, abs=atol), \
                f"{name} vertex {i} N: formula={fp_n:.6f} proc={pp_n:.6f}"

    def _check_item(self, variant_geom, item_name):
        variant, geom, ev = variant_geom
        vi = geom.get("variant_items", {})
        if item_name not in vi:
            pytest.skip(f"{item_name} not in {variant} variant")
        proc_poly = vi[item_name]["poly"]
        form = ev.elements.get(item_name)
        assert form is not None, \
            f"Formula did not produce result for {item_name} ({variant})"
        self._assert_poly_match(proc_poly, form["poly"],
                                f"{item_name}({variant})")

    # All-variant items
    def test_hamper(self, variant_geom):
        self._check_item(variant_geom, "hamper")

    def test_water_heater(self, variant_geom):
        variant, geom, ev = variant_geom
        vi = geom.get("variant_items", {})
        if "water_heater" not in vi:
            pytest.skip("no water_heater")
        proc = vi["water_heater"]
        form = ev.elements.get("water_heater")
        assert form is not None
        # Check center and radius match
        assert form["center"][0] == pytest.approx(proc["center"][0], abs=0.01)
        assert form["center"][1] == pytest.approx(proc["center"][1], abs=0.01)
        assert form["radius"] == pytest.approx(proc["radius"], abs=1e-6)

    def test_toilet_s(self, variant_geom):
        self._check_item(variant_geom, "toilet_s")

    def test_toilet_n(self, variant_geom):
        self._check_item(variant_geom, "toilet_n")

    def test_util_sink(self, variant_geom):
        self._check_item(variant_geom, "util_sink")

    def test_bath_sink(self, variant_geom):
        self._check_item(variant_geom, "bath_sink")

    def test_kitchen_sink(self, variant_geom):
        self._check_item(variant_geom, "kitchen_sink")

    def test_ice_maker(self, variant_geom):
        self._check_item(variant_geom, "ice_maker")

    def test_coffee_maker(self, variant_geom):
        self._check_item(variant_geom, "coffee_maker")

    def test_dining_table(self, variant_geom):
        self._check_item(variant_geom, "dining_table")

    def test_dining_chair_1(self, variant_geom):
        self._check_item(variant_geom, "dining_chair_1")

    def test_dining_chair_2(self, variant_geom):
        self._check_item(variant_geom, "dining_chair_2")

    def test_chair(self, variant_geom):
        self._check_item(variant_geom, "chair")

    def test_ottoman(self, variant_geom):
        self._check_item(variant_geom, "ottoman")

    def test_desk(self, variant_geom):
        self._check_item(variant_geom, "desk")

    def test_desk_chair(self, variant_geom):
        self._check_item(variant_geom, "desk_chair")

    # Standard/daybed-only items
    def test_stove(self, variant_geom):
        self._check_item(variant_geom, "stove")

    def test_dishwasher(self, variant_geom):
        self._check_item(variant_geom, "dishwasher")

    def test_fridge(self, variant_geom):
        self._check_item(variant_geom, "fridge")

    def test_work_counter(self, variant_geom):
        self._check_item(variant_geom, "work_counter")

    def test_microwave(self, variant_geom):
        self._check_item(variant_geom, "microwave")

    def test_north_counter(self, variant_geom):
        self._check_item(variant_geom, "north_counter")

    # Minik-only items
    def test_kitchen_counter(self, variant_geom):
        self._check_item(variant_geom, "kitchen_counter")

    def test_cooktop(self, variant_geom):
        self._check_item(variant_geom, "cooktop")

    def test_toaster(self, variant_geom):
        self._check_item(variant_geom, "toaster")

    def test_sofa(self, variant_geom):
        self._check_item(variant_geom, "sofa")

    # Standard-only items
    def test_loveseat(self, variant_geom):
        self._check_item(variant_geom, "loveseat")

    def test_et(self, variant_geom):
        self._check_item(variant_geom, "et")

    def test_loveseat2(self, variant_geom):
        self._check_item(variant_geom, "loveseat2")

    # Daybed-only items
    def test_shelves2(self, variant_geom):
        self._check_item(variant_geom, "shelves2")

    def test_et_east(self, variant_geom):
        self._check_item(variant_geom, "et_east")

    def test_daybed(self, variant_geom):
        self._check_item(variant_geom, "daybed")

    def test_et_west(self, variant_geom):
        self._check_item(variant_geom, "et_west")

    # Rocker (minik and daybed variants)
    def test_rocker(self, variant_geom):
        self._check_item(variant_geom, "rocker")


class TestVariantItemFormulaCounts:
    """Verify expected formula counts."""

    def test_formula_definitions_count(self):
        from app.evaluator import get_variant_item_formulas
        formulas = get_variant_item_formulas()
        # Should have formulas for ~30+ items (some with variant-specific records)
        assert len(formulas) >= 40

    def test_formula_definitions_structure(self):
        from app.evaluator import get_variant_item_formulas
        for name, variant, formula in get_variant_item_formulas():
            assert isinstance(name, str) and len(name) > 0
            assert variant is None or isinstance(variant, str)
            assert isinstance(formula, dict)
            assert "type" in formula

    def test_seeded_formulas_include_variant_items(self, fresh_db):
        from app.database import get_all_formulas
        formulas = get_all_formulas(variant="standard", db_path=fresh_db)
        names = {f["element_name"] for f in formulas}
        # Check a few expected variant items are seeded
        assert "hamper" in names
        assert "kitchen_sink" in names
        assert "stove" in names
        assert "fridge" in names

    def test_minik_formulas_loaded(self, fresh_db):
        from app.database import get_all_formulas
        formulas = get_all_formulas(variant="minik", db_path=fresh_db)
        names = {f["element_name"] for f in formulas}
        assert "kitchen_counter" in names
        assert "hamper" in names
        # stove should NOT be in minik
        has_stove = any(f["element_name"] == "stove" and f["variant"] == "minik"
                        for f in formulas)
        assert not has_stove

    def test_daybed_formulas_loaded(self, fresh_db):
        from app.database import get_all_formulas
        formulas = get_all_formulas(variant="daybed", db_path=fresh_db)
        names = {f["element_name"] for f in formulas}
        assert "shelves2" in names
        assert "daybed" in names
        assert "et_east" in names


class TestNewFormulaTypes:
    """Test new formula type evaluators (Phase 12e)."""

    def test_toilet_shape_basic(self):
        ev = _make_evaluator()
        formula = {
            "type": "toilet_shape",
            "center": [5, 5],
            "facing_dir": "north",
            "width_dir": "east",
        }
        result = ev._eval_toilet_shape(formula)
        assert result is not None
        assert "poly" in result
        assert len(result["poly"]) == len([p for p in result["poly"]])
        assert len(result["poly"]) == 34

    def test_bath_sink_shape_basic(self):
        ev = _make_evaluator()
        formula = {
            "type": "bath_sink_shape",
            "anchor": [5, 5],
            "along": "east",
            "outward": "north",
            "length": 2.0,
            "depth": 1.5,
        }
        result = ev._eval_bath_sink_shape(formula)
        assert result is not None
        assert "poly" in result
        # Should have rect corners + semicircular arc points
        assert len(result["poly"]) > 4

    def test_dining_triangle_basic(self):
        ev = _make_evaluator()
        formula = {
            "type": "dining_triangle",
            "base_center": [5, 2],
            "toward_apex": "south",
            "along_base": "east",
            "base_width": 2.5,
            "height": 3.0,
            "apex_radius": 1.0,
            "fillet_radius": 0.5,
        }
        result = ev._eval_dining_triangle("test_table", formula)
        assert result is not None
        assert "poly" in result
        assert len(result["poly"]) > 4  # Triangle with arcs has many points
        assert "ne_side" in result
        assert "nw_side" in result
        assert "base_center" in result

    def test_dining_chair_basic(self):
        ev = _make_evaluator()
        # First create a dining table
        table_formula = {
            "type": "dining_triangle",
            "base_center": [5, 2],
            "toward_apex": "south",
            "along_base": "east",
            "base_width": 2.5,
            "height": 3.0,
            "apex_radius": 1.0,
            "fillet_radius": 0.5,
        }
        ev.elements["dining_table"] = ev._eval_dining_triangle("dining_table",
                                                                 table_formula)
        chair_formula = {
            "type": "dining_chair",
            "table": "dining_table",
            "side": "ne_side",
            "chair_width": 1.5,
            "chair_depth": 1.75,
            "gap": 2.0 / 12.0,
        }
        result = ev._eval_dining_chair(chair_formula)
        assert result is not None
        assert len(result["poly"]) == 4

    def test_ellipse_rect(self):
        ev = _make_evaluator()
        formula = {
            "type": "ellipse_rect",
            "anchor": [5, 5],
            "along": "east",
            "outward": "north",
            "rx": 1.0,
            "ry": 0.75,
        }
        result = ev._eval_ellipse_rect(formula)
        assert result is not None
        assert len(result["poly"]) == 4
        # Check bbox dimensions
        assert result["bbox"]["e"] - result["bbox"]["w"] == pytest.approx(2.0)
        assert result["bbox"]["n"] - result["bbox"]["s"] == pytest.approx(1.5)

    def test_rotated_dir(self):
        ev = _make_evaluator()
        # Rotate east (1,0) by -45° should give (cos(-45), sin(-45)) ≈ (0.707, -0.707)
        d = ev.resolve_dir({"rotated": "east", "angle": -math.pi / 4})
        assert d[0] == pytest.approx(math.sqrt(2) / 2, abs=1e-10)
        assert d[1] == pytest.approx(-math.sqrt(2) / 2, abs=1e-10)

    def test_element_centroid(self):
        ev = _make_evaluator()
        ev.elements["test_elem"] = {"poly": [[0, 0], [2, 0], [2, 2], [0, 2]]}
        result = ev.resolve_point({"element_centroid": "test_elem"})
        assert result == pytest.approx([1, 1])

    def test_ray_circle_isect(self):
        ev = _make_evaluator()
        # Ray from (0,0) going east, target at (5,3), distance = 5
        # |t*(1,0) - (5,3)|² = 25
        # (t-5)² + 9 = 25 → (t-5)² = 16 → t = 5±4
        # farthest: t=9, nearest: t=1
        result = ev.resolve_point({"ray_circle_isect": {
            "origin": [0, 0],
            "dir": "east",
            "target": [5, 3],
            "distance": 5.0,
            "select": "farthest",
        }})
        assert result[0] == pytest.approx(9.0, abs=1e-10)
        assert result[1] == pytest.approx(0.0, abs=1e-10)

        result2 = ev.resolve_point({"ray_circle_isect": {
            "origin": [0, 0],
            "dir": "east",
            "target": [5, 3],
            "distance": 5.0,
            "select": "nearest",
        }})
        assert result2[0] == pytest.approx(1.0, abs=1e-10)

    def test_radius_key_length(self):
        ev = _make_evaluator(radii={"R_a7": 15.5})
        val = ev.resolve_length({"radius_key": "R_a7"})
        assert val == pytest.approx(15.5)
