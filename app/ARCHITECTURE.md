# ADU Editor — Architecture

## Overview

The ADU Editor is a single-page Flask application backed by SQLite.  It
provides an interactive canvas and data tables for viewing and editing
the parametric building design defined in the project's `floorplan/`,
`shared/`, and related packages.  The app imports from but never modifies
those packages (NF-4 constraint).

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
    ├─ app/engine.py       Geometry computation orchestrator
    │      ├─ floorplan/geometry.py     (outline, F-series)
    │      ├─ floorplan/layout.py       (interior walls, rooms)
    │      ├─ floorplan/openings.py     (O1–O11, RO1–RO7)
    │      ├─ floorplan/gen_floorplan.py (dimension endpoints)
    │      ├─ shared/geometry.py        (inner walls, polygons)
    │      ├─ shared/survey.py          (traverse, inset)
    │      └─ app/variants.py           (variant furniture)
    │             └─ shared/geometry.py  (seg_vecs, offset_pt, line_isect)
    │
    └─ app/static/          HTML template, JS, CSS (no build step)
```

All dependencies flow downward.  No circular imports exist.  `floorplan/`
and `shared/` never import from `app/`.

---

## Modules

### app/database.py — Persistence

Nine SQLite tables:

| Table | Rows | Purpose |
|-------|------|---------|
| `constants` | 143 | Named numeric constants (value, unit, category) |
| `outline_chain` | 18 | Closed outline segments (line/arc definitions) |
| `views` | 11 | Registered SVG generators and output paths |
| `shapes` | ~15 | Complex item shapes (polygon coordinate lists) |
| `variant_exclusions` | 4 | Per-variant element hiding (wall/opening exclusions) |
| `room_label_offsets` | 0 | User-adjusted room label positions (offset from centroid) |
| `undo_history` | 0–50 | Serialised undo/redo entries (action type, before/after state) |
| `elements` | 13+ | Interior walls (seeded) + user-added elements (type, name, properties JSON, variant) |
| `doors` | 7+ | Door configurations per opening (width, hinge side, swing direction, type) |

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
| `get_categories()` | Distinct category list |
| `get_outline_chain()` | 18 segment rows |
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

### app/engine.py — Computation

**`compute_geometry(constants_dict, variant="standard")`** orchestrates
the full pipeline:

1. **Patch** `floorplan.constants` with database values, then reload
   `floorplan.geometry`, `.layout`, `.openings` so module-scope code
   re-executes with patched values
2. **Survey traverse** → F-series alignment
3. **Outline geometry** → 20 F-series points, 18 segments, arc radii
4. **Inner walls** → 18 W-series inset segments, closed polygon
5. **Interior layout** → 13 interior walls (IW1–IW9, IW11–IW12, IW2O, IW2S; no IW10),
   rooms, appliance/furniture placements
6. **Variant exclusions** → filter interior walls and rough openings
   per variant (e.g., bare/sf exclude IW6 and RO5)
7. **Openings** → 12 outer (O1–O11, O8A), 7 rough (RO1–RO7)
8. **Variant items** → 20–31 furniture/appliance items per variant
9. **Room labels** → area-weighted centroids for 11 rooms (BEDROOM,
   UTIL_N, UTIL_S, KITCHEN, LIVING, BATH, OFFICE, E CLOSET,
   W CLOSET, STORAGE, WH), with DB-stored offsets applied; for SF
   variant, includes area values and highlight polygons
10. **Dimensions** → 18–22 dimension line endpoint pairs with distances

Returns a JSON-serialisable dict with all computed geometry.  Also
provides `generate_svg()` (runs generator scripts via subprocess) and
`get_svg_content()` (reads SVG files from disk).

**Module reloading** — The engine uses `importlib.reload()` on four
floorplan modules after patching constants.  This is necessary because
those modules use `from floorplan.constants import ...` at module scope;
a simple `setattr` on the constants module would leave the other modules
with stale values.  The reload sequence is ordered to respect import
dependencies: constants → geometry → layout → openings.

**Derived constants** — Five constants are recomputed from others after
patching: `WALL_EXTRA`, `AIR_GAP`, `DOOR_FLAT_FACE`,
`F8F9_INNER_TURN_R`, `CORNER_SW_R`.

### app/variants.py — Variant Furniture

Replicates positioning math from `gen_floorplan.py`'s `_render_appliances()`,
`_render_kitchen()`, and `_render_furniture()` functions.  Uses the same
`seg_vecs()` / `offset_pt()` / `line_isect()` utilities from
`shared/geometry.py`.

**Variant registry** — five variants with boolean flags:

| Variant | Label | Items |
|---------|-------|-------|
| standard | Standard | ~31 (full set) |
| minik | Small Kitchen | ~22 (cooktop, toaster, sofa, no stove/dishwasher) |
| daybed | Daybed | ~24 (daybed, shelves2, no loveseat/sofa) |
| bare | Room Dimensions | 0 (walls only, IW6/RO5 excluded) |
| sf | Square Footage | 0 (walls only, IW6/RO5 excluded; adds SF labels + highlight polygons) |

Each item is a dict with `type` (appliance/furniture/fixture), `poly`
(coordinate list), `bbox`, `label`, `shape` (rect/circle), and for
circles: `center` and `radius`.

### app/server.py — HTTP & SSE

**Flask app factory** `create_app(db_path)` registers all routes and
initialises the database.

**Geometry cache** — dict keyed by variant name, each entry holding
computed geometry and a dirty flag.  Protected by a threading lock.
`_invalidate()` marks all entries dirty and broadcasts an SSE event.

**SSE** — `GET /api/events` returns a `text/event-stream` response.
Each connected client gets a `queue.Queue`; `_broadcast()` pushes
messages to all queues.  Events: `constants_changed`, `geometry_changed`,
`svg_updated`, `element_changed`, `undo_status`, `connected`.  Keepalive every 30 seconds.

**Floorplan variant mapping** — When the Floorplan view is requested
with a `?variant=` parameter, the server maps variant names to
SVG file suffixes (standard → `floorplan.svg`, minik →
`floorplan_minik.svg`, daybed → `floorplan_db.svg`, bare →
`floorplan_bare.svg`, sf → `floorplan_sf.svg`).

**API endpoints** (29 total):

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
| GET | `/api/variants` | Variant names and labels |
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
| POST | `/api/elements/<id>/move` | Move element: constant-based (IW) or offset (API-23) |
| GET | `/api/version` | Server git describe + start time |
| POST | `/api/openings` | Create opening (API-24) |
| PUT | `/api/openings/<name>` | Update opening (API-25) |
| DELETE | `/api/openings/<name>` | Delete opening + door (API-26) |
| GET | `/api/doors` | List doors |
| POST | `/api/doors` | Create door (API-27) |
| PUT | `/api/doors/<opening_name>` | Update door (API-28) |
| DELETE | `/api/doors/<opening_name>` | Delete door (API-29) |

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

**Lifecycle:**
- On startup: loads stack from DB, sets position to end
- On record: truncates redo entries, appends, trims oldest if > 50, persists
- On undo/redo: applies state, adjusts position (no DB write needed — the
  database state changes happen through the apply function)

### app/apputil.py — Shared Serialisation Helpers

Utility functions used by both `engine.py` and `variants.py` for
converting geometry objects to JSON-serialisable dicts:

| Function | Purpose |
|----------|---------|
| `point_to_list(pt)` | `(E, N)` tuple → `[E, N]` list |
| `bbox_from_poly(poly)` | Polygon → `{w, s, e, n}` bounding box |
| `seg_to_dict(seg)` | `LineSeg`/`ArcSeg` → JSON-serialisable dict |

### app/static/js/dialogs.js — Dialog Framework

Generic modal dialog system.

| Function | Purpose |
|----------|---------|
| `Dialog.show(opts)` | Show dialog with title, fields, onSubmit/onCancel |
| `Dialog.close()` | Remove overlay |
| `parseOffsetString(str)` | Parse "6in east" → `{dx, dy}` in feet |

### app/static/js/tools.js — Move Tool

Client-side move tool with drag, ghost preview, axis constraints, and snap.

| Object/Function | Purpose |
|-----------------|---------|
| `IW_MOVE_AXIS` | Client-side mirror of `elements.py` axis mapping |
| `MoveTool` | State machine: active, startWorld, ghost, targets, origTransforms |
| `moveToolMouseDown(e)` | Start drag, create ghost clones |
| `moveToolMouseMove(e)` | Update ghost positions (shift-constrain, grid snap) |
| `moveToolMouseUp(e)` | Remove ghosts, commit move via API |
| `commitMove(targets, dx, dy)` | POST move for each target; auto-create override for furniture |
| `showOffsetDialog()` | Show offset dialog (Enter key trigger) |
| `findElementRecord(type, name)` | Look up DB record from App.state.elements |

### app/static/js/app.js — Client

Single-file client application (~1600 lines).  No build step or
framework.

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
outline → inner walls → interior walls → openings → furniture →
room labels → points → dimensions.  Each function creates SVG
elements (`<polygon>`, `<circle>`, `<line>`, `<text>`) in the
corresponding `<g>` layer.  Display toggles control which layers
are populated.  Stacked variant items (MICRO, coffee maker, etc.)
are sorted to render after their parent counters so they appear
on top in the SVG paint order.

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

The SVG canvas defines 10 layer groups: outline, inner, walls, openings,
furniture, rooms, points, labels, dims, measure.

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
| Outline chain | `outline_chain` table (seeded from `floorplan/geometry.py`) | Read-only (editable in Phase 5) |
| Interior walls | Computed by `floorplan/layout.py` from constants | Indirectly (edit constants) |
| Openings | Computed by `floorplan/openings.py` from constants | Indirectly (edit constants) |
| Variant items | Computed by `app/variants.py` from layout + constants | Indirectly (edit constants) |
| Dimension lines | Computed by `floorplan/gen_floorplan.py` from layout | Indirectly (edit constants) |
| SVG views | Generated by scripts listed in `views` table | Regenerate via menu |
| Variant exclusions | `variant_exclusions` table (seeded with bare/sf rules) | Read-only |
| Room label offsets | `room_label_offsets` table | Yes — move operation |
| Item shapes | `shapes` table (seeded from hardcoded data) | Read-only |

The constants table is the single editable root.  Every other geometric
value is deterministically derived from it through the computation
pipeline.

### Evolution Through Phases

The sources of truth evolve as phases are completed:

| Phase | Change to data model |
|-------|---------------------|
| 3 | `elements` and `doors` tables added.  Interior walls seeded as DB entities.  User-added custom elements stored with absolute positions.  Engine-computed items (furniture/appliances) overlaid with DB-stored custom items on the canvas. |
| 5 | `outline_chain` becomes editable.  DB chain is authoritative — the engine uses DB-stored chain parameters, not `floorplan/geometry.py`'s hardcoded chain. |
| 8 | Room labels stored as `elements` (type `'label'`).  `room_label_offsets` table subsumed — offsets and rotation stored per-element.  Auto-computed centroids remain the default position; DB stores offset + rotation from centroid. |
| 12 | **Cutover.**  All positioning becomes formula-driven.  Constants become DB-stored values (no longer Python module attributes).  Element positions defined by parametric formulas referencing other elements and/or constants.  Existing scripts retained only as seed sources for "Reset to Defaults." |

### Target Architecture (Phase 12+)

After cutover, the database is the **sole** authoritative source for all
design data.  Every element — exterior walls, interior walls, openings,
furniture, appliances, fixtures, labels, dimensions — is a database
entity with parametric position formulas.  "Reset to Defaults" regenerates
the entire database from the existing generator scripts' output,
reproducing the reference design.  See CHARTER.md § Design Principles
for the full parametric dependency model.

---

## Computation Flow

```
constants table
    │
    ▼
patch_constants()          Monkey-patch floorplan.constants module
    │
    ▼
importlib.reload()         Reload geometry, layout, openings modules
    │
    ├──► compute_outline_geometry()     → F-series points, segments, radii
    │
    ├──► compute_inner_walls()          → W-series inset path
    │
    ├──► compute_interior_layout()      → 13 interior walls, rooms
    │
    ├──► get_variant_exclusions()       → filter IW/RO per variant
    │
    ├──► compute_outer_openings()       → 12 outer openings
    │    compute_rough_openings()       → 7 rough openings
    │
    ├──► compute_variant_items()        → 20–31 furniture/appliance items
    │
    ├──► _compute_room_labels()         → 11 room centroids + areas
    │
    └──► compute_dimension_endpoints()  → 18–22 dimension line pairs
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

The snapshot/restore mechanism is necessary because `compute_geometry()`
mutates module-level state via `patch_constants()` + `importlib.reload()`.
Without it, one test's constant changes would leak into subsequent tests.

---

## NF-4 Constraint: No Modification of Existing Packages

The app imports from but never modifies `shared/`, `floorplan/`, `walls/`,
`span/`, `survey/`, `roof/`, `site/`, `scad/`, or `plumbing/`.  This
constraint applies during Phases 0–12.  It is lifted at cutover, when
the database becomes the sole authoritative source and code duplication
is consolidated (see CHARTER.md § Transition from Principles 4 → 5).

**Consequence: intentional duplication.**  `app/variants.py` replicates
~700 lines of positioning math from `floorplan/gen_floorplan.py`'s
`_render_appliances()`, `_render_kitchen()`, and `_render_furniture()`
functions.  It also carries 24 hardcoded item-dimension constants
(hamper, microwave, dining table, etc.) that duplicate values scattered
through the generator.  This duplication is deliberate: the existing
scripts are the reference implementation, and the app must reproduce
their output without modifying them.  At cutover, these constants move
into the database as the single source, and the shared positioning math
is extracted into functions callable by both the SVG renderer and the
app engine.

**Consequence: module reloading.**  The engine uses
`importlib.reload()` and `patch_constants()` to inject database values
into the floorplan modules at runtime (see Computation Flow above).
Five derived constants (`WALL_EXTRA`, `AIR_GAP`, `DOOR_FLAT_FACE`,
`F8F9_INNER_TURN_R`, `CORNER_SW_R`) are recomputed in `engine.py`
after patching because they depend on other constants and their
derivation formulas cannot be imported without modifying source.  This
will be unnecessary once the app owns the constants directly (Phase 12).

---

## Roadmap

The charter describes a full parametric editor; the current
implementation is a **parametric viewer with constant editing, undo, and
element/door CRUD** (Phases 0–3 complete) — 145 of 226 requirements are
implemented across 253 app tests (839 total).  Phase 1 established
automated test coverage for all implemented server-side requirements.
Phase 2 added undo/redo infrastructure.  Phase 3 added elements and doors
as first-class database objects with full CRUD APIs.  Next phase: Phase 4
(move tool).

The development arc follows three stages:

1. **Phases 0–1 (complete):** Parametric viewer.  Constants-only
   editing, read-only chain, engine-computed geometry.
2. **Phases 2–11:** Progressive element CRUD, canvas tools, and domain
   views.  The database gradually becomes authoritative (chain editing
   in Phase 5, element storage in Phase 3, room labels in Phase 8).
   Existing scripts remain the reference for default-value output.
3. **Phase 12 (cutover):** Parametric dependency system replaces
   procedural computation.  All constants and positions become
   database-driven formulas.  NF-4 is lifted.  Code duplication
   consolidated.

See `app/ROADMAP.md` for the complete 13-phase development plan
covering all remaining requirements, phase dependencies, new file
inventory, test growth estimates, anticipated challenges, and cutover
criteria.
