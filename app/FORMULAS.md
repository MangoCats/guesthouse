# ADU Editor — Formula Language Reference

Element geometry in the ADU editor is defined by **JSON formulas** stored in
the `element_formulas` table and evaluated by
`app/evaluator.py:FormulaEvaluator` in topological dependency order.  This
document is the reference for the formula language: the top-level formula
types and the three spec grammars (point, direction, length) they are built
from.

Source of truth for behaviour: the point/direction/length spec grammar is
implemented once in **`app/specs.py`** (shared by the FormulaEvaluator and by
`app/engine.py`'s dimension-anchor resolution via a `SpecContext` adapter);
the formula types (`_eval_*` methods) live in `app/evaluator.py`.  If this
document and the code disagree, the code wins — update this file.

## Evaluation model

- Each formula row is keyed `(element_name, param_name, variant)`.
  `param_name` is usually `position` or `poly`; `variant` NULL applies to
  all variants.
- Dependencies (element → element/point/constant) are recorded in
  `formula_deps` and drive a topological sort in `evaluate_all()`.
  Cycles raise `CycleError`.
- A **locked** formula (`locked=1`) returns its frozen `locked_value`
  instead of re-evaluating; upstream changes stop propagating through it.
- Formula results are dicts with at least `poly` (vertex list, building
  coordinates in feet) and `bbox` (`{w, s, e, n}`); some types add
  `center`, `radius`, `width`, `depth`, `rotation`, `area`.
- Rectangular polygons follow the `[SW, SE, NE, NW]` corner convention; face
  names map `south=(0,1)`, `east=(1,2)`, `north=(2,3)`, `west=(3,0)`.
- Available named points during evaluation: F-series, W-series, S/G-series,
  arc centers, FC, plus all previously-evaluated element polygons.

---

## Formula types

### wall_rect — interior wall rectangle

```json
{"type": "wall_rect",
 "anchor": <point>,          // start corner (SW)
 "along": <dir>,             // wall run direction
 "thickness_dir": <dir>,     // outer→inner face direction
 "thickness": <length>,
 "end_mode": "fixed" | "intersect" | "to_line",
 // end_mode-specific:
 "length": <length>,                             // fixed
 "end_target": "inner_poly", "select": "nearest", // intersect (far end on inner wall polygon)
 "end_line_point": <point>, "end_line_dir": <dir> // to_line
}
```

Real example (user wall `IW39`):

```json
{"type": "wall_rect",
 "anchor": {"offset": "W01", "dir": {"segment": ["W01", "W02"]}, "dist": 25.667},
 "along": {"perp": {"segment": ["W01", "W02"]}},
 "thickness_dir": {"perp": {"perp": {"segment": ["W01", "W02"]}}},
 "thickness": 0.5, "end_mode": "fixed", "length": 2.833}
```

### item_rect — furniture/appliance rectangle

```json
{"type": "item_rect",
 "anchor": <point>,
 "anchor_corner": "sw" | "se" | "nw" | "ne" | "center",   // default "sw"
 "along": <dir>, "across": <dir>,
 "width": <length>, "depth": <length>,
 "rotation_deg": <length>}   // optional, rotates around center
```

Result includes `center`, `width`, `depth`, `rotation` (world angle =
along-vector angle + `rotation_deg`).

### item_circle — circular item (water heater, end table)

```json
{"type": "item_circle", "center": <point>, "radius": <length>, "n_pts": 24}
```

### four_corner — four independent point specs

```json
{"type": "four_corner", "sw": <point>, "se": <point>, "ne": <point>, "nw": <point>}
```

Also the "frozen" form produced when an item is moved to literal
coordinates (each corner a bare `[E, N]` pair).

### wall_opening — opening spanning outer→inner wall face

```json
{"type": "wall_opening",
 "outer_start": <point>, "outer_end": <point>,   // outer face segment
 "inner_start": <point>, "inner_end": <point>,   // inner face segment
 "width": <length>,
 // one of four positioning modes:
 "from_end": true|false, "gap": <length>,        // 1. gap from segment end/start
 "ref_point": <point>, "ref_offset": <length>,   // 2. projection of a reference point
 "centered": true, "center_t": 0.5,              // 3. centered at parametric t
 "center_refs": [<point>, <point>],              // 4. centered between two projections
 "poly_order": ...}                               // optional vertex ordering
```

Parametric t-values are computed on the outer segment and applied to both
faces.  Real example (`RO11` in wall `IW15`):

```json
{"type": "wall_opening",
 "outer_start": {"element": "IW15", "corner": 0},
 "outer_end":   {"element": "IW15", "corner": 1},
 "inner_start": {"element": "IW15", "corner": 3},
 "inner_end":   {"element": "IW15", "corner": 2},
 "from_end": false, "gap": 1.000134, "width": 5.166667}
```

### Shape formulas

| Type | Purpose | Key fields |
|------|---------|-----------|
| `toilet_shape` | Fixed toilet SVG polygon transformed to building coords | anchor/direction fields |
| `bath_sink_shape` | D-shaped sink: rectangular base + semicircular outward bulge | wall anchor, length, depth |
| `dining_triangle` | Triangle table with apex arc + fillet corners; stores side endpoints for chairs | corner/arc fields |
| `dining_chair` | Chair at midpoint of a dining-table side, facing outward | table reference, side |
| `ellipse_rect` | Bounding rect of an ellipse on a wall face | `anchor`, `along`, `outward`, `rx`, `ry` |
| `shape_transform` | Named polygon from the `shapes` table placed by center + rotation (+ optional width/depth scaling) | `shape_name`, `center`, `rotation_deg`, `width`, `depth` |
| `area_poly` | Arbitrary polygon (room/SF areas) with vertex list and optional arc bulges | `vertices` (point specs), `arc_adjustments` (`from_idx`, `to_idx`, `center`, `sign`) |

`area_poly` real example (water-heater room `WH`):

```json
{"type": "area_poly",
 "vertices": [{"element": "IW2S", "corner": 2}, "W7", "W8", "W9",
              {"type": "line_intersection",
               "line1_point": {"element": "IW2S", "corner": 2}, "line1_dir": "north",
               "line2_point": "W9", "line2_dir": "east"}],
 "arc_adjustments": [{"from_idx": 1, "to_idx": 2, "center": "C7", "sign": 1},
                     {"from_idx": 2, "to_idx": 3, "center": "C8", "sign": -1}]}
```

---

## Point specs (`app/specs.py:resolve_point`)

The three spec grammars below are implemented once in `app/specs.py` and
used by both formulas and dimension anchors.  Which names resolve depends on
the caller's context: in formulas, every polygon lookup is against evaluated
elements; in dimension anchors, `face_*` specs see interior walls and
`opening_*` specs see outer/rough openings.

| Spec | Meaning |
|------|---------|
| `[E, N]` | Literal coordinates (feet) |
| `"F9"`, `"W18"`, … | Named geometry point |
| `{"ref": <point>, "dAlong": d, "dPerp": d}` | Offset in the local frame of a named chain point (tangent / 90° CCW of tangent) |
| `{"element": name, "corner": 0-3 or "SW"/"SE"/"NE"/"NW"}` | Element polygon corner |
| `{"element_centroid": name}` | Vertex-average centroid of an element polygon |
| `{"face_mid": name, "face": f}` | Midpoint of a named face of an element |
| `{"midpoint": [<point>, <point>]}` | Midpoint of two recursively-resolved specs |
| `{"offset": <point>, "dir": <dir>, "dist": <length>}` | Base + distance along direction |
| `{"type": "line_intersection", "line1_point", "line1_dir", "line2_point", "line2_dir"}` | Intersection of two infinite lines |
| `{"line_poly_isect": {"origin", "dir", "poly": "inner_poly", "select": "nearest"\|"farthest"}}` | Ray–polygon intersection |
| `{"arc_point": {"center", "radius_key", "reference", "side": "east"\|"west"}}` | Point on an arc circle at the reference's northing |
| `{"ray_circle_isect": {"origin", "dir", "target", "distance", "select"}}` | Point on ray at given distance from target |

(Additionally, `engine.py`'s dimension-anchor resolution supports
`{"opening_face_mid"}` and `{"opening_centroid"}` specs for opening
polygons.)

## Direction specs (`resolve_dir`)

| Spec | Meaning |
|------|---------|
| `[dE, dN]` | Literal vector |
| `"east" / "west" / "north" / "south"` | Cardinal unit vectors |
| `{"ref_tangent": name, "angle_offset_deg": a}` | Chain-point tangent, optionally rotated |
| `{"face_along": name, "face": f}` | Unit vector along an element face |
| `{"face_perp": name, "face": f}` | Right-hand perpendicular of a face |
| `{"segment": [<point>, <point>]}` | Unit vector A→B |
| `{"segment_perp": [<point>, <point>]}` | Right-hand perpendicular of A→B |
| `{"neg": <dir>}` | Negation |
| `{"rotated": <dir>, "angle": radians}` | Rotation by angle (radians, CCW) |
| `{"perp": <dir>}` | 90° CW rotation (`[dy, -dx]`) |

## Length specs (`resolve_length`)

| Spec | Meaning |
|------|---------|
| number | Literal feet |
| `{"const": name}` | Constant lookup (DB `constants` table) |
| `{"radius_key": name}` | Arc radius lookup (e.g. `"R_a5"`) |
| `{"neg": s}` / `{"add": [..]}` / `{"sub": [a, b]}` / `{"mul": [..]}` / `{"div": [a, b]}` | Arithmetic |
| `{"proj": {"target", "anchor", "dir"}}` | Signed projection of (target − anchor) onto dir |
| `{"dist": [<point>, <point>]}` | Euclidean distance |

---

## Editing formulas

- REST: formula endpoints in `app/server.py` (get/set per element+param,
  lock/unlock); `GET /api/deps/graph` returns the full dependency DAG.
- UI: properties panel → formula section (type, deps, lock/unlock with
  padlock indicators); selecting an element highlights upstream (blue) and
  downstream (orange) dependencies on the canvas.
- Moving a furniture item with the Move tool typically rewrites its formula
  to a literal `four_corner` (locked position); walls move by adjusting
  their controlling constant or anchor.
