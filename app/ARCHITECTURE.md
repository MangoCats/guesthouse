# ADU Editor — Architecture

## Overview

The ADU Editor is a single-page Flask application backed by SQLite.  It
provides an interactive canvas and data tables for viewing and editing
the parametric building design.  All element geometry (walls, openings,
furniture, appliances) is computed by the FormulaEvaluator from
database-stored JSON formulas (Phase 12h: procedural baselines eliminated).

### DB Authority

**The SQLite database (`app/adu.db`) is the sole source of truth for all
building geometry and dimensions.**

The `floorplan/` Python modules (`constants.py`, `geometry.py`, `layout.py`,
`openings.py`) are **seed/reference sources only**.  They are used in exactly
one context: populating a fresh database via File → Reset Database.  They are
never consulted during normal operation, SVG generation, or the `gen_all.py`
regeneration workflow.

Implications:
- `gen_all.py` builds `GeneratorData` from the DB and uses in-process
  generation for all supported scripts.  Scripts without in-process handlers
  (survey, SCAD) fall back to subprocess — acceptable because those derive
  geometry from survey traversal / 3D model parameters, not DB formulas.
- `compute_geometry()` must always receive `chain_rows` from the DB so that
  F-series/W-series points reflect the live outline, not hardcoded constants.
- `GeneratorData.iw_polys` provides DB-driven interior wall polygons for span
  generators; `gd.layout` (from `compute_interior_layout()`) is a legacy
  fallback retained only for the non-DB standalone path.

```
Browser (HTML/CSS/JS)
    │  REST API (/api/*)
    │  SSE (/api/events)
    ▼
app/server.py          Flask routes, SSE broadcast, geometry cache
    │
    ├─ app/database.py     SQLite schema, seeding, CRUD
    │      └─ floorplan/constants.py    (seed source)
    │      └─ floorplan/geometry.py     (outline chain source)
    │
    ├─ app/elements.py     Element business logic, IW→constant mapping
    │      └─ app/database.py     (element CRUD)
    │
    ├─ app/doors.py        Door validation (hinge/swing/type)
    │      └─ app/database.py     (door CRUD)
    │
    ├─ app/undo.py         Undo/redo manager (50-level command stack)
    │      └─ app/database.py     (state application via CRUD functions)
    │
    ├─ app/engine.py       Geometry computation + in-process generator dispatch
    │      ├─ app/gen_provider.py       (shared native geometry computation)
    │      ├─ floorplan/geometry.py     (outline, F-series)
    │      ├─ floorplan/gen_floorplan.py (dimension endpoints)
    │      ├─ shared/geometry.py        (inner walls, polygons)
    │      ├─ shared/survey.py          (traverse, inset)
    │      ├─ app/database.py           (element metadata, formulas)
    │      └─ generators (11 scripts)   (in-process render via Phase 19)
    │
    ├─ app/gen_provider.py Generator data provider (Phase 15)
    │      ├─ shared/geometry.py        (inner walls, F8-F9 polyline)
    │      ├─ shared/wall_shells.py     (S/G-series, wall sections)
    │      ├─ floorplan/roof.py         (R-series roof geometry)
    │      └─ floorplan/openings.py     (wall openings for sections)
    │
    └─ app/static/          HTML template, JS, CSS (no build step)
```

All dependencies flow downward.  No circular imports exist.  `floorplan/`
and `shared/` never import from `app/`.

---

## Modules

### app/database.py — Persistence

Seventeen SQLite tables:

| Table | Rows | Purpose |
|-------|------|---------|
| `constants` | 143+ | Named numeric constants (value, unit, category) |
| `outline_chain` | 18 | Closed outline segments (line/arc definitions) |
| `views` | 11 | Registered SVG generators and output paths |
| `shapes` | ~15 | Complex item shapes (polygon coordinate lists) |
| `variants` | 5+ | Variant definitions (name, label, flags, layer_config, is_builtin) |
| `variant_exclusions` | 4 | Per-variant element hiding (wall/opening exclusions) |
| `room_label_offsets` | 0 | User-adjusted room label positions (offset from centroid) |
| `undo_history` | 0–50 | Serialised undo/redo entries (action type, before/after state) |
| `elements` | ~59+ | All design elements seeded (walls, openings, variant items) + user-added |
| `doors` | 9+ | Door configurations per opening (width, hinge side, swing direction, type) |
| `element_formulas` | ~100+ | JSON position/poly formulas per element (locked, variant-scoped) |
| `formula_deps` | ~200+ | Formula dependency edges (element→element, element→point, element→constant) |
| `config` | 2+ | Application config key-value pairs (default_variant, etc.) |
| `plumbing_elements` | 17+ | Supply/drain pipes, fittings, fixture connections |
| `survey_legs` | 5 | Traverse legs (bearing deg/min/sec, distance ft/inch, label) |
| `survey_config` | 5 | Survey constants (FC_IN_P3, COORD_ROTATION, corrections) |
| `inner_wall_overrides` | 3+ | Parametric override chains for W-series inner wall segments (seg_index, span_end, sub_seq, seg_type L/CW/CCW, bearing, distance, radius, sweep, n_pts) |

**Seeding** — On first run, `init_db()` creates the schema and populates
tables from source:

- **constants**: parsed from `floorplan/constants.py` source lines via
  regex, extracting name, expression, unit, and comment.  Category
  assigned by name prefix (`WALL_*` → wall, `IW_*` → interior_wall,
  `O*_WIDTH` → opening, etc.)
- **outline_chain**: copied from `OUTLINE_CHAIN` and `CHAIN_POINT_NAMES`
  in `floorplan/geometry.py`
- **views**: hardcoded list of 11 generator registrations
- **shapes**: hardcoded complex item shapes (polygon coordinates)
- **variant_exclusions**: per-variant wall/opening exclusion rules
  (bare and sf variants exclude IW6 and RO5)
- **room_label_offsets**: empty by default; populated when the user
  moves a room label from its computed centroid

**CRUD functions**:

| Function | Purpose |
|----------|---------|
| `get_constants_dict()` | `{name: value}` dict for engine |
| `get_constant_value(name)` | Single constant value (for undo capture) |
| `get_all_constants()` | Full rows for UI table |
| `update_constant()` | Single constant update |
| `update_constants_batch()` | Multi-constant transaction |
| `reset_constants()` | Clear and re-seed from source |
| `reset_elements()` | Clear elements/doors, re-seed IW walls + doors |
| `restore_elements(elems, doors)` | Full replace from snapshot (undo) |
| `get_categories()` | Distinct category list |
| `get_outline_chain()` | All segment rows (18 default) |
| `get_outline_chain_row(seq)` | Single segment by seq |
| `update_outline_segment(seq, updates)` | Update segment fields |
| `insert_outline_segment(seq, row_data)` | Insert segment, shift seqs |
| `delete_outline_segment(seq)` | Delete segment, renumber |
| `restore_outline_chain(snapshot)` | Full chain replace (undo/rollback) |
| `get_views()` | Enabled view definitions |
| `get_shapes()` | All shape rows |
| `get_shape(name)` | Single shape by name |
| `get_variant_exclusions(variant)` | `{element_type: {names}}` for a variant |
| `get_room_label_offsets()` | `{name: (offset_e, offset_n)}` dict |
| `set_room_label_offset(name, e, n)` | Upsert a room label offset |
| `get_all_elements()` | All element rows |
| `get_element(id)` | Single element by ID |
| `get_element_by_name(name)` | Single element by name |
| `create_element(type, name, props, variant)` | Insert element, return row |
| `create_element_raw(record)` | Insert element from full dict (undo re-insert) |
| `update_element(id, fields)` | Update element fields |
| `delete_element(id)` | Delete element + cascade openings hosted on it |
| `get_all_doors()` | All door rows |
| `get_door(opening_name)` | Single door by opening name |
| `create_door(opening, width, hinge, swing, type)` | Insert door, return row |
| `create_door_raw(record)` | Insert door from full dict (undo re-insert) |
| `update_door(opening_name, fields)` | Update door fields |
| `delete_door(opening_name)` | Delete door by opening name |
| `get_variants()` | All variant rows (5 built-in + user-created) |
| `get_variant(name)` | Single variant by name |
| `get_variant_by_id(id)` | Single variant by ID |
| `update_variant(id, updates)` | Update variant fields (layer_config, label) |
| `create_variant(name, label, source, flags)` | Create user variant, return row |
| `create_variant_raw(record)` | Insert variant from full dict (undo re-insert) |
| `delete_variant(id)` | Delete variant by ID |
| `clone_variant_exclusions(source, target)` | Copy exclusion rows for new variant |
| `delete_variant_exclusions(variant)` | Remove all exclusion rows for a variant |
| `clone_variant_elements(source, target)` | Add target to element visibility lists |
| `unclone_variant_elements(target)` | Remove target from element visibility lists |
| `get_inner_wall_overrides()` | `{seg_index: [chain dicts]}` for all overrides |
| `get_inner_wall_override(seg_index)` | Chain dicts for single segment |
| `upsert_inner_wall_override(seg_index, chain, span_end)` | Create/update override chain |
| `delete_inner_wall_override(seg_index)` | Remove override for segment |
| `reset_inner_wall_overrides()` | Re-seed defaults (W8-W9) |
| `check_override_overlap(seg_index, span_end)` | Validate no overlapping spans |
| `snapshot_inner_wall_overrides()` | Capture full state for undo |
| `restore_inner_wall_overrides(snapshot)` | Restore from undo snapshot |

Connection management uses a context manager with WAL journaling and
foreign keys enabled.

### app/elements.py — Element Business Logic

Maps interior walls to their controlling constants and hosted openings.

| Data / Function | Purpose |
|-----------------|---------|
| `IW_CONSTANT_MAP` | Dict mapping IW name → controlling constant name (e.g., `"IW1"` → `"IW1_OFFSET_FROM_W9"`) |
| `IW_MOVE_AXIS` | Dict mapping IW name → `(axis, sign)` or `None`; axis is `"x"`/`"y"`, sign is `+1`/`-1` |
| `IW_HOSTED_OPENINGS` | Dict mapping IW name → list of hosted RO names (e.g., `"IW9"` → `["RO3", "RO7"]`) |
| `compute_constant_delta(iw_name, dx, dy)` | Translate world-coordinate move to `(constant_name, delta)` or `None` |
| `get_elements_for_variant(variant)` | Return elements visible to a variant (variant=NULL or matching) |
| `get_controlling_constant(iw_name)` | Return the constant that controls an IW's position |
| `get_hosted_openings(iw_name)` | Return list of RO names hosted by an IW |

### app/doors.py — Door Validation

Validates door parameters against allowed values.

| Data / Function | Purpose |
|-----------------|---------|
| `VALID_SIDES` | Set of allowed hinge/swing values: `{east, west, north, south}` |
| `VALID_TYPES` | Set of allowed door types: `{single, double}` |
| `validate_door(hinge, swing, type)` | Return error string or `None` |

### app/gen_provider.py — Generator Data Provider (Phase 15 + 15½)

Provides `GeneratorData`, a typed container of native Python geometry objects
for SVG generators.  Replaces direct imports from `floorplan/` modules.
Also implements the inner wall override engine (Phase 15½).

| Function / Class | Purpose |
|------------------|---------|
| `compute_native_geometry(constants, chain_rows, db_path)` | Shared steps 1-3: survey traverse, outline, inner walls.  Returns `(pts, outline_segs, inner_segs, radii)` as native objects.  Used by both `compute_geometry()` and `GeneratorData`. |
| `build_generator_data(constants, chain_rows, db_path)` | Main entry point — builds `GeneratorData` from DB state.  When `db_path` is provided, calls `compute_geometry()` once with `chain_rows` to populate both `openings` (DB-driven outer openings) and `iw_polys` (DB-driven interior wall polygons for span generators). |
| `GeneratorData` | Contains `pts`, `outline_segs`, `inner_segs`, `radii`, `outline_poly`, `inner_poly`, `s_segs`, `g_segs`, `w_f8f9_poly`, `g_f8f9_poly`, `openings`, `wall_sections`, `roof`, `roof_poly`, `layout`, `iw_polys`, `constants`, `wall_t`, `outer_area`, `inner_area`.  **`iw_polys`** is a `{name: poly}` dict of formula-evaluated interior wall polygons (DB-driven); span generators prefer this over `layout.iwN.poly`.  **`layout`** is from `compute_interior_layout()` (hardcoded seed values) and is retained for the standalone subprocess path only. |
| `walk_override_chain(chain, start_pt, start_bearing)` | Parametric chain walk: line (bearing+distance) and arc (radius+sweep CW/CCW) segments. Returns polyline. |
| `compute_default_override(seg_index, inner_segs, pts, constants)` | Computes default parametric chain for a single inner wall segment. |
| `compute_default_span_override(seg_index, span_end, ...)` | Computes default chain for a multi-segment span. |
| `apply_overrides_to_poly(inner_poly, inner_segs, pts, overrides)` | Splices override polylines into `inner_poly` in-place, replacing default offset geometry. Processes in descending index order. |

### app/engine.py — Computation

**`compute_geometry(constants_dict, variant="standard", chain_rows=None, doors_data=None)`**
orchestrates the full pipeline:

1. **Survey traverse** → F-series alignment (via `compute_native_geometry()`)
2. **Outline geometry** → 20 F-series points, 18 segments, arc radii
3. **Inner walls** → 18 W-series inset segments, closed polygon
3b. **Inner wall overrides** → `apply_overrides_to_poly()` splices
   DB-stored override chains (from `inner_wall_overrides` table) into
   `inner_poly`, replacing default offset geometry with parametric
   line/arc chains.  Supports single-segment and multi-segment spans.
4. **Formula evaluation** → `FormulaEvaluator` evaluates all DB-stored
   formulas in topological order, producing element polygons
5. **Build elements from formulas** → `_build_elements_from_formulas()`
   queries DB element metadata (label, type, shape, door config,
   clearance config, product URLs, variant lists) and combines with
   FormulaEvaluator output to build interior walls, openings, and
   variant items dicts
6. **Variant exclusions** → filter interior walls and rough openings
   per variant (e.g., bare/sf exclude IW6 and RO5)
7. **Door arcs** → structural door arcs (from `doors_data`), appliance
   door arcs (from DB element door configs), clearance zones (from DB
   element clearance configs)
8. **Room labels** → area-weighted centroids for 11 rooms (BEDROOM,
   UTIL_N, UTIL_S, KITCHEN, LIVING, BATH, OFFICE, E CLOSET,
   W CLOSET, STORAGE, WH), with DB-stored offsets applied; for SF
   variant, includes area values and highlight polygons
9. **Dimensions** → all dimension elements (builtin + user) from the
   `elements` table are resolved via `_resolve_anchor()` to produce
   18–22 dimension line endpoint pairs with distances

Returns a JSON-serialisable dict with all computed geometry.  Also
provides `generate_svg_db()` (in-process generator dispatch with
DB-built `GeneratorData`, subprocess fallback), `generate_svg()` (legacy
subprocess path), and `get_svg_content()` (reads SVG files from disk).

**DB-driven regeneration (Phase 19 + 19½)** — `build_generator_data_from_db(db_path)`
constructs a `GeneratorData` from the current database state.
`_run_generator_inprocess(script_path, gd)` dispatches to 11 generator
render functions (floorplan variants including plumbing, roof, walls, span,
site plan, SCAD, survey).  The plumbing SVG is a floorplan variant
(`render_floorplan_svg_db()` with `variant="plumbing"`), not a standalone
generator.  `generate_svg_db(view_name, script_path, gd)` tries in-process
first, falls back to subprocess for scripts without handlers (gen_views,
gen_line_drawings, gen_3views).  All regeneration API endpoints in
`server.py` use this path so generated SVGs reflect DB state.

**Door arc computation** — `_compute_door_arcs(opening_polys, doors_data)`
resolves abstract hinge/swing directions (east/west/north/south) to
concrete positions by projecting onto opening polygon geometry.  For each
door, computes hinge position, open-tip position, and a 21-point 90° arc
polyline.  Double doors produce two leaves hinged at opposite ends.

**Clearance zone computation** — `_compute_clearance_zones(variant_items)`
reads clearance metadata from DB element properties (face vertex indices +
distance).  Computes outward-extending rectangles using centroid-based
direction check.  No hardcoded clearance data remains in code.

**Appliance door computation** — `_compute_appliance_doors(variant_items)`
reads intrinsic door configs from DB element properties (hinge_idx +
target_idx) and computes arcs via `_swing_arc()`.  Propagates `stacked`
flag for SVG paint-order correctness.

**Formula-only geometry (Phase 12h)** — All element geometry is computed
by `FormulaEvaluator` from database-stored JSON formulas.  The engine no
longer calls `compute_interior_layout()`, `compute_outer_openings()`,
`compute_rough_openings()`, or `compute_variant_items()`.  Element
metadata (label, type, shape, door/clearance configs, product URLs,
variant lists) is stored in the `elements` table and queried by
`_build_elements_from_formulas()`.  The procedural `floorplan/` modules
provide outline geometry only.  `patch_constants()` is no longer called
from span analysis (Phase 14-C); span functions use `compute_geometry()`
result directly.  Survey traverse is computed from DB tables when
`db_path` is provided (Phase 14-B).

**Derived constants** — `_derive_constant()` computes `WALL_EXTRA`,
`CORNER_SW_R`, and other derived values from the constants dict without
modifying any module state.

**Radii exposure** — `result["radii"]` is included in the geometry
response, exposing all arc radii (R_a1, R_a5, ..., R_a19).  This is
needed by the `arc_point` point spec to resolve dimension anchors that
reference points on arc segments.

**Anchor resolution system** — `_resolve_anchor(anchor, geometry_result)`
resolves a dimension anchor dict to `[E, N]` coordinates.  Supports
five anchor types:

| Anchor type | Resolution |
|-------------|------------|
| `point` | Named geometry point from `result["points"]` |
| `wall_face` | Midpoint of a named face (south/east/north/west) on an interior wall polygon |
| `opening_face` | Midpoint of a named face on an outer or rough opening polygon |
| `line_intersection` | Intersection of two infinite lines, each defined by a point spec + direction spec |
| `computed` | Arbitrary point resolved via `_resolve_point_spec` |

**Point spec resolution** — `_resolve_point_spec(spec, geom)` resolves
a point specification to `[E, N]` coordinates.  Point specs can be:

| Spec form | Resolution |
|-----------|------------|
| `str` | Named point from `geom["points"]` (e.g. `"F9"`, `"W18"`) |
| `{face_mid: T, face: F}` | Midpoint of wall face |
| `{opening_face_mid: T, face: F}` | Midpoint of opening face |
| `{opening_centroid: T}` | Vertex-average centroid of opening polygon |
| `{midpoint: [A, B]}` | Midpoint of two recursively resolved point specs |
| `{offset: base, dir: D, dist: N}` | Base point + distance along direction |
| `{arc_point: {center, radius_key, reference, side}}` | Point on arc circle at reference northing |

**Direction spec resolution** — `_resolve_dir_spec(spec, geom)` resolves
a direction specification to a `[dE, dN]` unit vector:

| Spec form | Resolution |
|-----------|------------|
| `"east"` | `[1, 0]` |
| `"north"` | `[0, 1]` |
| `{face_along: T, face: F}` | Unit vector along wall face direction |
| `{face_perp: T, face: F}` | Right-hand perpendicular of wall face |
| `{segment: [A, B]}` | Unit vector from point A to point B |
| `{segment_perp: [A, B]}` | Right-hand perpendicular of A-to-B direction |

**Unified variant filtering** — When resolving dimension and label
elements, `properties["variants"]` (a list of variant names) takes
precedence over the DB `variant` column.  If `variants` is present in
properties, the element appears only in those listed variants.  If
absent, the DB `variant` column is used (NULL = all variants, specific
value = that variant only).  This allows builtin dimensions like
`dim12bare` to target specific variants via the DB column, while
user-created elements can use the properties list for multi-variant
visibility.

**App solver bypass (Phase 5)** — When `chain_rows` is provided (always
in the web editor), `compute_geometry()` uses `app/outline_solver.py`
instead of `compute_outline_geometry()`.  The app solver walks the DB
chain, solves closure, and produces F-series points and outline segments
directly.  This makes the DB chain authoritative and allows user edits
to propagate through the full geometry pipeline.

### app/outline_solver.py — Closure Solver (Phase 5)

Pure-math reimplementation of the outline closure solver from
`floorplan/geometry.py`, with zero imports from `floorplan/`.

| Function | Purpose |
|----------|---------|
| `chain_offset(chain, start_brg)` | Walk chain entries, return (dE, dN, exit_brg) |
| `solve_closure(chain, R_a1)` | Solve d_F2_F5, d_F18_F1, and sweep_closure for chain closure |
| `db_rows_to_chain(rows)` | Convert DB row dicts to ChainEntry NamedTuples |
| `walk_chain(chain, F2_E, F2_N)` | Full point generation → WalkResult(points, radii) |
| `check_closure(chain, flex_specs)` | Non-mutating closure check → {valid, closure_error, ...} |
| `solve_for_constraint(...)` | Secant method for target distance constraints |

Cross-validation tests verify bit-identical results with `floorplan/geometry.py`
for default chain parameters.

### app/variants.py — Variant Registry (149 lines)

Contains the `VARIANTS` dict (legacy reference for variant flags),
`get_variant_flags()` (reads from DB with dict fallback), dimension
constants used by formula evaluation, and product URL lookup (legacy
reference).  All procedural positioning math was removed in Phase 12h
(835 lines deleted).

**Variant registry** — six built-in variants stored in the `variants` DB
table (Phase 11 + Phase 20).  User-created variants clone from a source
variant and inherit its flags, exclusions, element visibility, and formulas.

| Variant | Label | Items | Built-in |
|---------|-------|-------|----------|
| standard | Standard | ~31 (full set) | Yes |
| minik | Small Kitchen | ~22 (cooktop, toaster, sofa, no stove/dishwasher) | Yes |
| daybed | Daybed | ~24 (daybed, shelves2, no loveseat/sofa) | Yes |
| bare | Room Dimensions | 0 (walls only, IW6/RO5 excluded) | Yes |
| sf | Square Footage | 0 (walls only, IW6/RO5 excluded; adds SF labels + highlight polygons) | Yes |
| plumbing | Plumbing | standard items (plumbing canvas overlay) | Yes |
| (user) | (custom) | Matches source variant | No |

Each item is a dict with `type` (appliance/furniture/fixture), `poly`
(coordinate list), `bbox`, `label`, `shape` (rect/circle), and for
circles: `center` and `radius`.  Items may also carry optional metadata:
`door` (hinge_idx, target_idx) for appliance door arcs,
`clearance` (face vertex indices, distance) for clearance zones, and
`stacked` (boolean) for items rendered above their parent counter.
All metadata is stored in the `elements` table properties JSON and
queried by `_build_elements_from_formulas()` — no metadata remains
hardcoded in `variants.py`.

### app/server.py — HTTP & SSE

**Flask app factory** `create_app(db_path)` registers all routes and
initialises the database.

**Geometry cache** — dict keyed by variant name, each entry holding
computed geometry and a dirty flag.  Protected by a threading lock.
`_invalidate()` marks all entries dirty and broadcasts an SSE event.

**SSE** — `GET /api/events` returns a `text/event-stream` response.
Each connected client gets a `queue.Queue`; `_broadcast()` pushes
messages to all queues.  Events: `constants_changed`, `geometry_changed`,
`svg_updated`, `element_changed`, `outline_changed`, `undo_status`, `connected`.  Keepalive every 30 seconds.

**Floorplan variant mapping** — When the Floorplan view is requested
with a `?variant=` parameter, the server maps variant names to
SVG file suffixes (standard → `floorplan.svg`, minik →
`floorplan_minik.svg`, daybed → `floorplan_db.svg`, bare →
`floorplan_bare.svg`, sf → `floorplan_sf.svg`).

**API endpoints** (70 total):

| Method | Path | Purpose |
|--------|------|---------|
| GET | `/` | Serve index.html |
| GET | `/api/constants` | List constants (optional `?category=`) |
| GET | `/api/constants/categories` | Category list |
| PUT | `/api/constants/<name>` | Update one constant |
| PUT | `/api/constants/batch` | Batch update |
| POST | `/api/constants/reset` | Reset to source defaults |
| POST | `/api/undo` | Undo last mutation |
| POST | `/api/redo` | Redo last undone mutation |
| GET | `/api/geometry` | Computed geometry (`?variant=`) |
| GET | `/api/variants` | Variant list from DB (name, label, flags, layer_config) |
| PUT | `/api/variants/<id>` | Update variant (layer_config, label) |
| POST | `/api/variants` | Create user variant (clone from source) |
| DELETE | `/api/variants/<id>` | Delete user variant (built-in protected) |
| GET | `/api/shapes` | List all shapes |
| POST | `/api/shapes` | Create shape |
| PUT | `/api/shapes/<name>` | Update shape |
| GET | `/api/outline` | Outline chain (18 segments) |
| GET | `/api/views` | Registered views |
| GET | `/api/svg/<name>` | SVG content |
| GET | `/api/svg/<name>/file` | File download |
| POST | `/api/regenerate` | Run generator scripts |
| GET | `/api/events` | SSE stream |
| GET | `/api/elements` | List elements (optional `?variant=`) |
| POST | `/api/elements` | Create element (API-20) |
| PUT | `/api/elements/<id>` | Update element (API-21) |
| DELETE | `/api/elements/<id>` | Delete element + cascade (API-22) |
| POST | `/api/elements/<id>/move` | Move element: constant-based (seeded IW), anchor-translation (user IW), or offset (API-23) |
| POST | `/api/interior-walls` | Create user IW wall with `wall_rect` formula; auto-names IW13, IW14… |
| GET | `/api/version` | Server git describe + start time |
| POST | `/api/openings` | Create opening (API-24) |
| PUT | `/api/openings/<name>` | Update opening (API-25) |
| DELETE | `/api/openings/<name>` | Delete opening + door (API-26) |
| GET | `/api/doors` | List doors |
| POST | `/api/doors` | Create door (API-27) |
| PUT | `/api/doors/<opening_name>` | Update door (API-28) |
| DELETE | `/api/doors/<opening_name>` | Delete door (API-29) |
| PUT | `/api/outline/<seq>` | Update outline segment (API-16) |
| POST | `/api/outline/validate` | Dry-run chain validation (API-17) |
| POST | `/api/outline/add-point` | Insert F-point by splitting (API-18) |
| DELETE | `/api/outline/<seq>` | Remove point from chain (API-19) |
| GET | `/api/span-data` | Span analysis data (ANALYSIS-1) |
| GET | `/api/span-rotation` | Span rotation min/max (ANALYSIS-2) |
| GET | `/api/config/<key>` | Get config value (SCAD-2, SITE-1) |
| PUT | `/api/config/<key>` | Set config value (SCAD-2, SITE-1) |
| POST | `/api/reset-database` | Reset DB to defaults |
| GET | `/api/survey-points` | P-series survey points (SITE-4) |
| GET | `/api/survey/legs` | Survey traverse legs |
| PUT | `/api/survey/legs/<seq>` | Update a survey leg |
| GET | `/api/survey/config` | Survey configuration |
| PUT | `/api/survey/config/<key>` | Update survey config value |
| POST | `/api/survey/reset` | Reset survey to defaults |
| GET | `/api/project/export` | Export full project as JSON |
| POST | `/api/project/import` | Import project from JSON |
| POST | `/api/generate-site-plan` | Regenerate site plan PDFs (SITE-1) |
| POST | `/api/generate-3d` | Generate SCAD 3D model (SCAD-1) |
| POST | `/api/generate-views` | Generate multi-view PDF (SCAD-3) |
| GET | `/api/plumbing` | List all plumbing elements (PLUMB-5) |
| POST | `/api/plumbing` | Create plumbing element (PLUMB-5) |
| PUT | `/api/plumbing/<id>` | Update plumbing element (PLUMB-5) |
| DELETE | `/api/plumbing/<id>` | Delete plumbing element (PLUMB-5) |
| GET | `/api/deps/graph` | Full formula dependency DAG (nodes + edges) (FORM-20) |

### app/plumbing.py — Plumbing CRUD Module

Manages the `plumbing_elements` table: supply pipes, drain pipes, fittings,
and fixture connections.  Includes reference plumbing computation and seeding.

- `PLUMBING_TYPES`: `supply_pipe`, `drain_pipe`, `fitting`, `fixture_connection`
- `FIXTURE_DEFS`: 11 seeded fixture connections with cold/hot/drain flags
- `get_plumbing_elements()`, `get_plumbing_element(id)` — read
- `create_plumbing_element(type, name, path, properties, fixture)` — create
- `update_plumbing_element(id, updates)`, `delete_plumbing_element(id)` — mutate
- `create_plumbing_raw(record)` — raw re-insert for undo
- `seed_plumbing(conn)` — seed 11 default fixture connections (idempotent)
- `compute_reference_plumbing(geom, wall_t)` — compute reference pipe paths
  and fixture positions from geometry (replicates gen_floorplan.py routing)
- `seed_reference_plumbing(geom, wall_t, db_path)` — seed reference pipes
  into DB if none exist (called on startup and database reset)

### app/undo.py — Undo/Redo Manager

Command-pattern undo system with 50-level depth (UNDO-3).

**`UndoManager(db_path, max_depth=50)`** — maintains an in-memory stack
of undo entries with a position pointer.  Entries are persisted to the
`undo_history` table.

| Method | Purpose |
|--------|---------|
| `record(action_type, before_state, after_state, desc)` | Append entry, trim redo, enforce depth |
| `undo()` | Apply `before_state`, return entry or `None` |
| `redo()` | Apply `after_state`, return entry or `None` |
| `can_undo` / `can_redo` | Boolean properties |

**State dispatch (`_apply`):** Dispatches by `action_type`:

| Action type | Undo behaviour |
|-------------|---------------|
| `constant_update`, `constant_batch`, `constant_reset` | Apply `{name: value}` dict via `update_constants_batch()` |
| `element_create` | Delete the created element by ID |
| `element_delete` | Re-insert full element record(s) via `create_element_raw()` |
| `element_update` | Restore old field values via `update_element()` |
| `door_create` | Delete the created door by opening name |
| `door_delete` | Re-insert full door record via `create_door_raw()` |
| `door_update` | Restore old field values via `update_door()` |
| `opening_create` | Delete the created opening (element) by ID |
| `opening_delete` | Re-insert full opening record(s) via `create_element_raw()` |
| `opening_update` | Restore old field values via `update_element()` |
| `element_move` (constant) | Restore constant value via `update_constant()` |
| `element_move` (position) | Restore element properties via `update_element()` |
| `outline_update` | Restore full chain snapshot via `restore_outline_chain()` |
| `outline_add_point` | Restore full chain snapshot (N rows → N-1) |
| `outline_remove_point` | Restore full chain snapshot (N rows → N+1) |
| `full_reset` | Restore both constants dict and outline chain snapshot |

**Lifecycle:**
- On startup: loads stack from DB, sets position to end
- On record: truncates redo entries, appends, trims oldest if > 50, persists
- On undo/redo: applies state, adjusts position (no DB write needed — the
  database state changes happen through the apply function)

### app/labels.py — Label & Dimension Helpers (Phase 8)

Room label seeding, builtin dimension seeding, and auto-name generators
for user labels/dimensions.

| Data / Function | Purpose |
|-----------------|---------|
| `ROOM_LABEL_NAMES` | List of 11 room label names (BEDROOM, UTIL_N, UTIL_S, KITCHEN, LIVING, BATH, OFFICE, E CLOSET, W CLOSET, STORAGE, WH) |
| `BUILTIN_DIMENSIONS` | List of 22 anchor-based dimension specs.  Each entry is a dict with `name`, `variant` (None or specific), `start_anchor`, and `end_anchor`.  Anchors use the five anchor types supported by `_resolve_anchor()` in `engine.py` (point, wall_face, opening_face, line_intersection, computed).  Shared sub-expressions (e.g. `_F9F11_MID_OFFSET`, `_W18W1_DIR`) are factored out as module-level constants. |
| `seed_room_labels(conn)` | Create 11 room label elements (type `'label'`, source `'room'`) if not present; migrate offsets from legacy `room_label_offsets` table |
| `seed_builtin_dimensions(conn)` | Create 22 builtin dimension elements (type `'dimension'`, source `'builtin'`) if not present.  Each dimension stores its start/end anchors in properties JSON for resolution during `compute_geometry()` |
| `_next_auto_name(elements, prefix, pattern)` | Shared helper: scan element names matching `pattern`, return next sequential name with `prefix` (e.g. `UD4`, `UL3`) |
| `next_dimension_name(elements)` | Return next auto-name for a user dimension (`UD1`, `UD2`, ...) via `_next_auto_name` |
| `next_label_name(elements)` | Return next auto-name for a user label (`UL1`, `UL2`, ...) via `_next_auto_name` |

### app/style.py — Per-Element Style Management (Phase 9)

Defaults, validation, and per-view resolution for element visual properties
(fill colour, stroke colour/width/style, opacity).  Style properties are
stored in the element `properties` JSON column — no schema changes required.

| Data / Function | Purpose |
|-----------------|---------|
| `TYPE_DEFAULTS` | Dict mapping element type → default style dict (fill_color, stroke_color, stroke_width, stroke_style, opacity) matching the CSS class defaults in `app.css` |
| `STYLE_KEYS` | Tuple of the five style property names |
| `VALID_STROKE_STYLES` | Tuple: `"solid"`, `"dashed"`, `"dotted"` |
| `validate_color(value)` | Return `True` if value is `None` or valid CSS hex colour (`#rgb` / `#rrggbb`) |
| `validate_opacity(value)` | Return `True` if value is a number in [0, 100] |
| `validate_stroke_style(value)` | Return `True` if value is a recognised stroke style |
| `validate_stroke_width(value)` | Return `True` if value is a non-negative number |
| `validate_style_props(props)` | Validate all style keys present in a properties dict; returns `(ok, error_msg)` |
| `get_defaults(element_type)` | Return a copy of the default style dict for the given type |
| `resolve_style(element_type, props, variant=None)` | Three-layer merge: type defaults → base element props → `view_overrides[variant]`.  Returns a fully-populated style dict |

**Per-view overrides** are stored in `properties.view_overrides` as a dict
keyed by variant name, each value a partial style dict.  The frontend resolves
styles client-side using the same merge order.

### app/parse_input.py — Unit-Aware Dimension Parser

Parses user-entered dimension strings and returns a value in feet.
Supports feet, inches, centimetres, millimetres, and metres, with
flexible whitespace and foot-mark/inch-mark symbols.

| Data / Function | Purpose |
|-----------------|---------|
| `_CONVERSIONS` | Dict mapping unit names to feet conversion factors (`ft`, `in`, `cm`, `mm`, `m` and long forms) |
| `_TOKEN_RE` | Regex pattern matching number + optional unit (foot-mark `'`, inch-mark `"`, or word unit) |
| `parse_dimension(text)` | Parse a dimension string, return value in feet or `None`.  Rules: bare number = feet; number + unit = converted; two bare numbers = feet + inches; fractions (`1/3`) evaluated as division; up to 3 tokens summed |

### app/apputil.py — Shared Serialisation Helpers

Utility functions and constants used by `engine.py`, `variants.py`, and
`database.py` for converting geometry objects to JSON-serialisable dicts
and standardising polygon approximations:

| Constant / Function | Purpose |
|---------------------|---------|
| `ARC_N_SEMICIRCLE` | Arc discretisation: 32 segments for semicircular arcs (bath sink bulge) |
| `ARC_N_CIRCLE` | Arc discretisation: 24 segments for full circles (water heater, ET) |
| `point_to_list(pt)` | `(E, N)` tuple → `[E, N]` list |
| `bbox_from_poly(poly)` | Polygon → `{w, s, e, n}` bounding box |
| `seg_to_dict(seg)` | `LineSeg`/`ArcSeg` → JSON-serialisable dict |

### app/static/js/dialogs.js — Dialog Framework

Generic modal dialog system.

| Function | Purpose |
|----------|---------|
| `Dialog.show(opts)` | Show dialog with title, fields, onSubmit/onCancel. Supports `customContent` (DOM node) and `presetButtons` (target + values array) |
| `Dialog.close()` | Remove overlay |
| `parseOffsetString(str)` | Parse "6in east" → `{dx, dy}` in feet |

### app/static/js/tools.js — Move Tools

Client-side move tool and shared utilities.

| Object/Function | Purpose |
|-----------------|---------|
| `IW_MOVE_AXIS` | Client-side mirror of `elements.py` axis mapping |
| `IW_HOSTED_OPENINGS` | Client-side mirror for cascade delete warnings |
| `MoveTool` | State machine: active, startWorld, ghost, targets, origTransforms |
| `moveToolMouseDown/Move/Up` | Move tool mouse handlers |
| `commitMove(targets, dx, dy)` | POST move for each target; auto-create override for furniture |
| `showOffsetDialog()` | Show offset dialog (Enter key trigger) |
| `findElementRecord(type, name)` | Look up DB record from App.state.elements |

### app/static/js/catalog.js — Catalog & Placement Tool

Hardcoded catalog data and click-to-place element creation.

| Object/Function | Purpose |
|-----------------|---------|
| `CATALOG` | Furniture (8), appliance (6), fixture (3) item arrays |
| `PlaceTool` | Placement state: active, itemTemplate, itemType |
| `showCatalog(type)` | Open catalog grid dialog for item type |
| `startPlacement(item, type)` | Enter placement mode |
| `placeToolMouseDown(e)` | Canvas click creates element at position |
| `placeElement(wx, wy)` | POST placed element with `rectPoly()` polygon |

### app/static/js/shape-editor.js — Shape Editor

Interactive shape editing dialog with SVG canvas and API integration.

| Object/Function | Purpose |
|-----------------|---------|
| `ShapeEditor.shapes` | Cached shapes from API |
| `loadShapes()` | GET /api/shapes into cache |
| `showShapeEditor(name, cb)` | Open modal with draggable vertex handles |
| `parseSvgPolygon(text)` | Extract points from SVG polygon/path markup |
| `addShapePicker(tbody, rec, props)` | Shape dropdown in properties panel |

### app/static/js/app.js — Client

Single-file client application (~2100 lines).  No build step or
framework.

**API helper** — `apiFetch(url, opts)` wraps `fetch()` and throws on
non-ok responses, extracting error messages from JSON bodies.  All data
loading and mutation calls use this wrapper so HTTP errors surface as
readable messages rather than silent JSON parse failures.

**State** — `App.state` object holds: current geometry, constants array,
views list, active view/tool, zoom/pan, display toggles, variant
selection, sort/filter state, selections array (multi-select).

**Initialisation** (on DOMContentLoaded):
1. `cacheElements()` — populate `App.els` with ~30 DOM references
2. `setupEventListeners()` — attach handlers to buttons, inputs, canvas
3. `connectSSE()` — open EventSource to `/api/events`
4. `loadViews()` → `renderViewTabs()` — build tab bar
5. `loadConstants()` → `renderConstantsTable()` — populate right panel
6. `loadGeometry()` → `renderCanvas()` — draw interactive view

**Canvas rendering** — `renderCanvas()` calls layer functions in order:
outline → inner walls → interior walls → openings → doors →
furniture → clearance zones → room labels → points → dimensions.
Each function creates SVG elements (`<polygon>`, `<circle>`, `<line>`,
`<text>`) in the corresponding `<g>` layer.  Display toggles control
which layers are populated.  Stacked variant items (MICRO, coffee
maker, etc.) are sorted to render after their parent counters so they
appear on top in the SVG paint order.  Stacked appliance door arcs
(e.g., microwave) are rendered inside `renderFurniture()` after all
item polygons, rather than in `renderDoors()`, so they appear above
their parent counter in SVG paint order.  Move offsets are computed
once per render cycle via `itemOverrides()` and passed to
`renderDoors()`, `renderFurniture()`, and `renderClearanceZones()`.

**Room labels** — `renderRoomLabels(g)` renders room name text at
centroid positions (with DB offsets applied).  For the SF variant,
it also renders dashed partition lines and clickable area polygons
that highlight on click.

**Units formatting** — `fmtFtIn(ft)` converts feet to `X' Y.YY"`
format with trailing zeroes removed, matching `shared/geometry.py`'s
`fmt_dist()`.  Used in all value displays: coordinate readout, property
panel, constants table, measurement tool, dimension labels, and data
tables.

### app/templates/index.html — Layout

Single-page five-region layout:

```
┌──────────────────────────────────────────────┐
│  Menu bar (32px)                             │
├────────┬─────────────────────────┬───────────┤
│  Tool  │  View tabs + variant    │  Panel    │
│  palet │  selector               │  tabs     │
│  te    ├─────────────────────────┤           │
│ (140px)│  Canvas / SVG viewer    │  Props /  │
│        │                         │  Consts / │
│        │                         │  Outline /│
│        ├─────────────────────────┤  Openings │
│        │  Status bar (24px)      │  (320px)  │
└────────┴─────────────────────────┴───────────┘
```

The SVG canvas defines 12 layer groups: outline, inner, walls, openings,
doors, furniture, clearance, rooms, points, labels, dims, measure.

### app/static/css/app.css — Styling

Dark theme (Catppuccin palette) via CSS custom properties.  18 colour
variables, 5 layout dimension variables.  Major style sections: menu bar,
tool palette, viewport, tabs, tables, right panel, SVG canvas classes
(outline/wall/opening/furniture/point/dimension/measure/room-label
fills and strokes), scrollbar, toast notifications, room label
highlights and SF partition lines.  Responsive breakpoint at
1000px collapses tool labels and narrows the right panel.

---

## Sources of Truth

### Current State (Phases 0–1)

| Data | Source of truth | Editable via app? |
|------|----------------|-------------------|
| Named constants | `constants` table (seeded from `floorplan/constants.py`) | Yes — inline editing |
| Outline chain | `outline_chain` table (seeded from `floorplan/geometry.py`) | Editable: PUT, add-point, delete; closure re-solved automatically |
| Interior walls | Computed by `floorplan/layout.py` from constants | Indirectly (edit constants) |
| Openings | Computed by `floorplan/openings.py` from constants | Indirectly (edit constants) |
| Variant items | Computed by `app/variants.py` from layout + constants | Indirectly (edit constants) |
| Dimension lines | Computed by `floorplan/gen_floorplan.py` from layout | Indirectly (edit constants) |
| SVG views | Generated by scripts listed in `views` table | Regenerate via menu |
| Variant exclusions | `variant_exclusions` table (seeded with bare/sf rules) | Read-only |
| Room label offsets | `room_label_offsets` table | Yes — move operation |
| Item shapes | `shapes` table (seeded from hardcoded data) | Read-only |
| Plumbing elements | `plumbing_elements` table (seeded with reference pipes/fixtures) | Yes — full CRUD via API |
| Survey legs | `survey_legs` table (5 legs seeded from `shared/survey.py`) | Yes — PUT endpoint |
| Survey config | `survey_config` table (FC_IN_P3, COORD_ROTATION, corrections) | Yes — PUT endpoint |

The constants table is the single editable root for building geometry.
Every other geometric value is deterministically derived from it through
the computation pipeline.  Plumbing elements are independently editable.

### Evolution Through Phases

The sources of truth evolve as phases are completed:

| Phase | Change to data model |
|-------|---------------------|
| 3 | `elements` and `doors` tables added.  Interior walls seeded as DB entities.  User-added custom elements stored with absolute positions.  Engine-computed items (furniture/appliances) overlaid with DB-stored custom items on the canvas. |
| 5 | `outline_chain` becomes editable.  DB chain is authoritative — the engine uses DB-stored chain parameters, not `floorplan/geometry.py`'s hardcoded chain. |
| 8 | Room labels stored as `elements` (type `'label'`).  `room_label_offsets` table subsumed — offsets and rotation stored per-element.  Auto-computed centroids remain the default position; DB stores offset + rotation from centroid. |
| 10 | `plumbing_elements` table added.  Pipes, fittings, and fixture connections stored as DB entities with full CRUD API.  Reference plumbing (6 pipes, 11 fixture positions) seeded on startup from geometry computation.  Re-seeded on database reset. |
| 12a | `element_formulas` and `formula_deps` tables added.  `FormulaEvaluator` class with topo sort, cycle detection, `wall_rect`/`item_rect`/`item_circle` evaluators.  24 variant item constants seeded into `constants` table. |
| 12b | 7 formula REST API endpoints.  Properties panel formula section (type, deps, lock/unlock). |
| 12c | Evaluator extensions: `four_corner` type, `proj`/`dist`/`neg`/`add`/`sub`/`mul` length specs, `neg`/`perp` dir specs.  All 13 IW wall formulas written, verified to 1e-9 ft vs procedural, seeded into DB.  Hybrid engine active (formula results override procedural). |
| 12d | `wall_opening` formula type (4 positioning modes: gap/ref_point/centered/center_refs; 4 poly_order options).  5 layout item formulas, 12 outer opening formulas, 7 rough opening formulas — total 37 formulas seeded.  Formula overrides preserve extra fields (counter clip). |
| 12e | 5 new formula types (`toilet_shape`, `bath_sink_shape`, `dining_triangle`, `dining_chair`, `ellipse_rect`).  New specs: `element_centroid`, `ray_circle_isect` point specs; `rotated` dir spec; `radius_key` length spec.  ~50 variant item formulas across 3 variants (standard, minik, daybed).  Hybrid engine overrides variant items preserving metadata (door, clearance, product_url). |
| 12f | Lock/unlock UI: padlock icon on locked canvas elements, lock/unlock button in properties panel with emoji indicators.  Formula dependency highlighting: blue (upstream) and orange (downstream) on element select.  Lock/unlock undo/redo (`formula_lock` action type).  `formula_locked` SSE event.  `GET /api/deps/graph` endpoint (full DAG as nodes + edges).  `locked_elements` list in geometry response.  `get_all_formula_deps()` database function. |
| 12g | Cutover: removed `patch_constants()`/`importlib.reload()` from `compute_geometry()`, added outer/rough opening formula overrides, BED formula (`div` length spec), lifted NF-4. |
| 12h | Eliminated procedural baselines: seeded ~59 element metadata rows (label, type, shape, door/clearance configs, product URLs, variant lists) into `elements` table.  Rewrote `compute_geometry()` as formula-only — removed all procedural calls.  Removed 835 lines of dead code from `variants.py` (984→149 lines).  Product URLs moved to DB.  Variant formula cloning on custom variant creation. |
| 14 | Fully DB-driven outline & export/import.  F2 origin (`F2_EASTING`, `F2_NORTHING`) as DB constants.  `survey_legs` and `survey_config` tables with 5 API endpoints.  Span analysis refactored to use `compute_geometry()` result — `patch_constants()` no longer called from span functions.  Full project export/import with outline closure and DAG cycle validation. |

### Current Architecture (Phase 14 Complete)

The database is the **sole** authoritative source for all design data.
Every element — interior walls, openings, furniture, appliances,
fixtures, labels, dimensions — is a database entity with parametric
position formulas and metadata.  `compute_geometry()` no longer imports
from `floorplan/layout.py` or `floorplan/openings.py`.  "Reset to
Defaults" regenerates the entire database from seed sources.  See
CHARTER.md § Design Principles for the full parametric dependency model.

---

## Computation Flow

```
constants table + element_formulas table + elements table
    │
    ▼
compute_outline_geometry()     → F-series points, segments, radii
    │
    ▼
compute_inner_walls()          → W-series inset path
    │
    ▼
apply_overrides_to_poly()      → splice DB-stored inner wall overrides
    │                             (line/arc chains replacing offset segments)
    ▼
FormulaEvaluator               → evaluate all formulas in topo order
    │                             (walls, openings, variant items)
    ▼
_build_elements_from_formulas()  → query DB element metadata (label,
    │                               type, shape, door/clearance configs,
    │                               product URLs, variant lists) and
    │                               build interior_walls, openings,
    │                               variant_items dicts
    ▼
get_variant_exclusions()       → filter IW/RO per variant
    │
    ├──► _compute_door_arcs()           → structural door arcs (9 openings)
    │    _compute_appliance_doors()     → appliance door arcs (from DB configs)
    │    _compute_clearance_zones()     → clearance zone polygons (from DB configs)
    │
    ├──► _compute_room_labels()         → 11 room centroids + areas
    │
    ├──► resolve dimensions (builtin + user) from elements table
    │         _resolve_anchor()          → resolve start/end anchors
    │         _resolve_point_spec()      → resolve point specs (string, face_mid, etc.)
    │         _resolve_dir_spec()        → resolve direction specs (east, north, etc.)
    │
    └──► resolve label elements          → room labels + user labels
              │
              ▼
         JSON response → browser → canvas rendering
```

Changing any constant triggers a full recomputation.  The geometry cache
avoids redundant computation when the same variant is requested without
intervening constant changes.

---

## Real-Time Update Cycle

1. User edits a constant in the browser
2. `PUT /api/constants/<name>` updates the database
3. Server marks all geometry caches dirty
4. Server broadcasts `geometry_changed` via SSE
5. Client receives event, calls `loadGeometry()`
6. `GET /api/geometry?variant=<current>` triggers recomputation
7. Server returns new geometry JSON
8. Client re-renders canvas with updated positions

The same cycle applies to batch updates and constant resets.  SVG
regeneration follows a parallel path: `POST /api/regenerate` runs
generator scripts as subprocesses and broadcasts `svg_updated`.

---

## Test Infrastructure

Tests live in `tests/test_zapp_*.py` (the `zapp` prefix sorts them after
the pre-existing geometry/layout tests alphabetically, ensuring app
tests run last).  Shared fixtures in `tests/test_zapp_conftest.py`
(imported explicitly by app test files, not auto-discovered by pytest)
provide:

- **`fresh_db`** — isolated SQLite database in `tmp_path`, with
  module-level numeric snapshot/restore across `floorplan.constants`,
  `.geometry`, `.layout`, `.openings` to prevent test pollution
- **`app_client`** — Flask test client backed by the fresh database
- **`geometry`** — pre-computed geometry fixture

The snapshot/restore mechanism was originally necessary because
`compute_geometry()` mutated module-level state via `patch_constants()`
+ `importlib.reload()`.  Since Phase 14-C, no engine function calls
`patch_constants()` — it remains only for SVG generator endpoints.
The snapshot mechanism is retained for test isolation.

---

## NF-4 Constraint: Lifted (Phase 12g)

NF-4 has been lifted.  The FormulaEvaluator is the sole source of element
geometry — `compute_geometry()` no longer patches `floorplan.constants` or
reloads modules.  The database is the authoritative source for constants,
element formulas, and element metadata.

**No remaining duplication (Phase 12h).**  `app/variants.py` has been
reduced from 984 to 149 lines.  All procedural positioning math has been
removed.  Element metadata (labels, door configs, clearance configs,
product URLs, shapes) is stored in the `elements` table and queried by
`_build_elements_from_formulas()`.  `compute_geometry()` no longer imports
from `floorplan/layout.py` or `floorplan/openings.py`.

**Module patching.**  `patch_constants()` is retained for SVG generator
functions that call into `floorplan/` modules, but is no longer used by
the span analysis functions (Phase 14-C refactored them to use
`compute_geometry()` result directly).

---

## Roadmap

The charter describes a full parametric editor; the current
implementation is a **parametric viewer with constant editing, undo,
element/door CRUD, move tool, outline chain editing, enhanced canvas
rendering, element creation/deletion/rotation tools, and label/dimension
annotations** (Phases 0–8 complete) — 197 of 226 requirements are
implemented across 629 app tests (1215 total).  Phase 1 established
automated test coverage for all implemented server-side requirements.
Phase 2 added undo/redo infrastructure.  Phase 3 added elements and
doors as first-class database objects with full CRUD APIs.  Phase 4
added the move tool with drag-to-reposition, ghost preview, axis
constraints, grid snap, multi-select group move, and offset dialog.
Phase 5 added the outline closure solver, editable chain parameters,
and canvas F-point dragging.  Phase 6 added door arc rendering
(structural + appliance), clearance zones, and display toggles.
Phase 7 added draw wall tool, catalog placement, delete with cascade,
rotate, multi-select, opening width editing, and shape editor.
Phase 8 added user dimension tool, label tool, room label migration
to elements table, context menu for dimension rotation, inline label
editing, and font size control.
Phase 9 added per-element styling (fill/stroke/opacity), product URL
links, and view override controls.
Phase 10 added domain views: site plan (10a), 3D SCAD (10b), analysis
(10c), and plumbing layout (10d) with interactive pipe/fixture editing,
reference plumbing seeding, and full CRUD for plumbing elements.
Phase 11 added variant selection UI and polish.
Phase 12 (a–h) implemented the formula-driven architecture: evaluator,
schema, constant consolidation (12a), formula REST API (12b), IW wall
formulas (12c), layout item formulas (12d), opening formulas (12e),
variant item formulas (12f), cutover (12g) — removing
`patch_constants`/`importlib.reload` from `compute_geometry()` and
lifting NF-4, and procedural baseline elimination (12h) — seeding all
element metadata into DB, removing all procedural calls from the engine,
and reducing `variants.py` from 984 to 149 lines.

The development arc followed three stages:

1. **Phases 0–1:** Parametric viewer.  Constants-only editing,
   read-only chain, engine-computed geometry.
2. **Phases 2–11:** Progressive element CRUD, canvas tools, and domain
   views.  The database gradually became authoritative (chain editing
   in Phase 5, element storage in Phase 3, room labels in Phase 8).
3. **Phase 12 (complete):** Parametric dependency system replaced
   procedural computation.  All element geometry is computed by the
   FormulaEvaluator from database-stored JSON formulas.  All element
   metadata stored in DB.  NF-4 lifted.  No procedural baselines remain.

Phase 15½ generalised the hardcoded W8-W9 inner wall corner treatment
into a database-driven system.  `inner_wall_overrides` table stores
parametric line/arc chains keyed by inner wall segment index.
`walk_override_chain()` evaluates chains; `apply_overrides_to_poly()`
splices results into `inner_poly`.  Supports single-segment and
multi-segment spans with endpoint validation and overlap detection.
Editor UI provides a dialog with compute-default, click-to-define mode,
and live canvas preview.

Phase 19½ migrated the plumbing SVG from a standalone generator to a
floorplan variant.  `render_floorplan_svg_db()` handles
`variant="plumbing"` with plumbing pipes, ghosted furniture/openings,
and a supplies table.  The plumbing view seed points to
`floorplan/gen_floorplan.py`.

Phase 20 made plumbing a layout variant in the variant selector (like
standard/minik/daybed/bare/sf).  The standalone Plumbing Edit and
Plumbing tabs were removed.  Interactive view renders the plumbing
canvas when variant is plumbing; Floorplan view loads the plumbing SVG.

See `app/ROADMAP.md` for the complete development plan
covering all remaining requirements, phase dependencies, new file
inventory, test growth estimates, anticipated challenges, and cutover
criteria.
