# ADU Editor — Development Roadmap

This roadmap traces the path from the current **parametric viewer with constant
editing** to the charter's goal: a **full parametric editor** where every
element can be created, moved, edited, deleted, and persisted from the browser
— all without AI assistance.

NF-4 has been lifted (Phase 12g cutover complete).  The FormulaEvaluator
is the sole source of element geometry.  Phase 12h eliminated all procedural
element baselines — `compute_geometry()` is now formula-only.

Phase 14 is complete: F2 origin, survey traverse, and span analysis are
fully DB-driven; `patch_constants()` is no longer called from span analysis;
full project export/import with validation is available.

---

## Current State (Phase 14 complete)

**275 of 275 requirements implemented.**  ~970 app tests, ~580 pre-existing
tests (~1550 total).  All implemented requirements have automated test coverage.

| Capability | Status |
|------------|--------|
| Constants: view, edit, persist, reset | Done |
| Outline chain: view, edit, add/remove, closure solver | Done |
| Geometry engine: full recomputation | Done |
| Interactive canvas: outline, walls, openings, furniture, points, dims | Done |
| Five layout variants | Done |
| 15 view tabs (SVG + PDF + PNG + plumbing canvas) with zoom/pan | Done |
| Properties panel with related constants | Done |
| Constants table: sort, filter, inline edit, category colours | Done |
| Openings table: outer + rough | Done |
| REST API: 50 endpoints, SSE | Done |
| Real-time update cycle | Done |
| Feet-inches display (NF-6) | Done |
| Unit-aware dimension input parser (CT-7a–j, CT-8) | Done |
| Room labels with area-weighted centroids (CV-8 partial, ENG-12 partial) | Done |
| SF variant: room areas, partition lines, clickable highlights | Done |
| Realistic item shapes: toilet, bath sink, dining table | Done |
| Item labels on canvas | Done |
| Stacked item rendering (microwave, cooktop, toaster, coffee maker) | Done |
| Database-driven variant exclusions (bare/sf hide IW6, RO5) | Done |
| Layout selector on Floorplan SVG tab (variant-specific SVGs) | Done |
| DB tables: shapes, variant_exclusions, room_label_offsets, undo_history | Done |
| Undo/redo: 50-level, Ctrl+Z / Ctrl+Shift+Z, cross-type | Done |
| Elements table: 13 IW seeds, CRUD, variant filtering | Done |
| Doors table: 9 seeds (7 RO + O3, O6), CRUD, validation | Done |
| Element/door/opening undo support (9 action types) | Done |
| Elements tab: interior walls table, furniture/appliance table | Done |
| Door properties in selection panel (hinge, swing, type) | Done |
| Furniture/appliance properties in selection panel | Done |
| element_changed SSE event on mutations | Done |
| Move tool: drag IW walls (constant-based) | Done |
| Move tool: drag furniture/appliance (offset override) | Done |
| Move tool: ghost preview, shift-constrain, grid snap | Done |
| Move tool: offset dialog (Enter key, parsed input) | Done |
| Move tool: multi-select (Ctrl+Click) + group move | Done |
| Move tool: undo/redo for all move types | Done |
| Dialog framework (dialogs.js) | Done |
| Outline solver (app/outline_solver.py) — closure, walk, validate | Done |
| Outline APIs: PUT, validate, add-point, DELETE (API-16–19) | Done |
| Outline undo: 3 new action types (12 total) | Done |
| Outline table: editable cells, closure indicator, +/- toolbar (DT-2–4) | Done |
| Canvas: F-point drag (OE-1), arc handle (OE-2) | Done |
| Engine: chain_rows parameter, app solver bypass (ENG-11) | Done |
| Door arc rendering: structural (9 openings) + appliance (4 items) | Done |
| Double door rendering (RO6, RO7 two-leaf arcs) | Done |
| Clearance zones: dresser, stove, D/W, hamper (4 per standard variant) | Done |
| Display toggles: Doors (DIS-10), Clearance (DIS-9) | Done |
| Appliance door arcs: fridge, washer, dryer, microwave | Done |
| Stacked appliance doors render above counters (SVG paint order) | Done |
| Door arcs + clearance zones follow move offsets | Done |
| Delete key: remove selected elements with confirmation | Done |
| Cascading delete: walls → hosted openings/doors | Done |
| Multi-select: Shift+Click and rubber-band drag-select | Done |
| Opening width editing: inline in table + properties panel | Done |
| Draw Wall tool: two-click placement, thickness editing | Done |
| Add Element tools: furniture/appliance/fixture from catalog | Done |
| Catalog dialog with grid layout, click-to-place mode | Done |
| Rotate tool: R key, rotation dialog with presets | Done |
| Shape editor: vertex editing, SVG import, shape assignment | Done |
| Shape API: GET/POST/PUT /api/shapes endpoints | Done |
| Dimension line tool: two-click placement with ft-in labels | Done |
| Room label tool: click-to-place, inline text editing, font size | Done |
| Anchored dimensions: snap-to-geometry (point, wall_face, opening_face) | Done |
| Unified builtin dimensions: 22 seeded anchored elements (line_intersection) | Done |
| Dimension/label selectability: click to select, properties panel | Done |
| Per-dimension style: solid/dashed visual control | Done |
| Multi-layout visibility: per-element variant checkbox list | Done |
| Cache invalidation on element/opening/door CRUD | Done |
| Delete button in properties panel | Done |
| Per-element fill colour, stroke colour/width/style (STYLE-1, STYLE-2) | Done |
| Per-element opacity slider (STYLE-3) | Done |
| Per-view style overrides (STYLE-4) | Done |
| Product URL field + Open button (LINK-1, SEL-12) | Done |
| SVG link wrapping for elements with URLs (LINK-2) | Done |
| Canvas link icon overlay (CV-12) | Done |
| Style module: defaults, validation, resolution (app/style.py) | Done |
| Analysis: room areas toggle, span tooltip, rotation min/max | Done |
| 3D/SCAD: config table, roof style, generate-3d/views endpoints | Done |
| Site plan: PDF view tabs, setback config, survey points, elements | Done |
| DB error handling: validation, reset, error banner UI | Done |
| Plumbing canvas: building ghost overlay, pipe/fitting/fixture rendering | Done |
| Plumbing drawing tools: cold supply, hot supply, drain (polyline) | Done |
| Plumbing fitting placement: tee, elbow, valve with rotation | Done |
| Plumbing fixtures table panel with cold/hot/drain indicators | Done |
| Plumbing DB: plumbing_elements table, CRUD API, undo/redo | Done |
| DB-driven variants: `variants` table, 5 built-in seeds, CRUD API | Done |
| Per-variant layer config: toggle persistence per variant | Done |
| User-defined variants: create (clone from source), delete, undo/redo | Done |
| Variant element/exclusion cloning on create, cleanup on delete | Done |
| Dynamic variant dropdown (populated from DB, + / × buttons) | Done |
| Furniture property editing: editable Width/Depth for placed items | Done |
| Separate User Dims toggle (DIS-7) | Done |
| Wall endpoint drag handles for drawn walls (TL-17) | Done |
| Add Opening placement tool for walls (TL-21) | Done |
| Constant dependency highlighting: white first-order, pink second-order (SEL-15) | Done |
| Formula evaluator: topo sort, cycle detection, wall_rect/item_rect/item_circle | Done |
| Formula API: 7 REST endpoints, properties panel formula section | Done |
| IW wall formulas: 13 walls, verified to 1e-9 ft vs procedural | Done |
| Layout + opening formulas: 37 formulas (5 layout, 12 outer, 7 rough openings) | Done |
| Variant item formulas: ~50 formulas across 3 variants, 5 new formula types | Done |
| Lock/unlock UI: padlock icon on canvas, lock button in properties panel | Done |
| Formula dependency highlighting: blue upstream, orange downstream on select | Done |
| Lock/unlock undo/redo and formula_locked SSE event | Done |
| Dependency graph endpoint: GET /api/deps/graph (full DAG) | Done |
| DB-seeded element metadata: label, type, shape, door/clearance configs | Done |
| Formula-only engine: no procedural calls, all elements from formulas + DB | Done |
| Product URLs in DB element properties (editable via API) | Done |
| Variant formula cloning on custom variant creation | Done |
| F2 origin from DB constants (F2_EASTING, F2_NORTHING) | Done |
| Survey data in DB (survey_legs, survey_config tables) | Done |
| DB-driven traverse: `_compute_traverse_from_db()` replaces hardcoded | Done |
| Span analysis without `patch_constants()` | Done |
| Survey API: 5 endpoints (GET/PUT legs, GET/PUT config, POST reset) | Done |
| Full project export (`GET /api/project/export`) | Done |
| Full project import with validation (`POST /api/project/import`) | Done |

**What's missing:** Phases 15–18 (DB-driven SVG generators), and electrical
layout (aspirational).  The database is the sole authoritative source for
all element data, survey traverse, and outline origin.  No procedural element
baselines remain in `compute_geometry()`.  Span analysis uses geometry
result directly without module patching.  SVG generators still read from
hardcoded Python modules; edits made in the interactive editor are not
reflected in generated SVGs (Phases 15–18 scope).

---

## Phase 1 — Foundation Test Coverage (Complete)

**Goal:** Achieve 100% automated test coverage for all 121 implemented
requirements before any new development.  This is the safety net.

**Completed.** 22 new tests added (774 total, up from 752).

**Coverage audit results:**
- **Already covered:** DB-2–8, ENG-1–3, ENG-6–10, ENG-13, API-1–8, API-10–15,
  API-32–34, CT-7a–j, CT-8, module isolation, variant items, dimensions
- **Extended:** DB-1 (all 6 tables), test_zapp_database.py (DB-12, DB-13,
  shapes), test_zapp_engine.py (ENG-14, ENG-15), test_zapp_api.py (API-9,
  variant SVG paths)
- **Frontend-only (manual acceptance):** ~55 requirements (GEN, UI, CV, DIS,
  SEL, TL, CT-1–6, CT-9–10, DT-1/5/6/8, RT, APP, NF-1–2/5–6) verified by
  visual inspection, not automatable server-side
- **Verified by infrastructure:** NF-3 (all 774 tests pass), NF-4 (git diff)

**Deviation from estimate:** 22 new tests vs. estimated 30. The difference is
because DB-1 was modified in place (not a new test), and several requirements
initially flagged as gaps (ENG-13, API-34) were already covered in
test_zapp_variants.py.

**Dependencies:** None.

---

## Phase 2 — Undo/Redo System (Complete)

**Goal:** Before adding mutation capabilities, establish the undo/redo
infrastructure so every future mutation is automatically reversible.

**Requirements:** DB-11, API-30–31, UNDO-1–4 (7 reqs)

**Completed.** 20 new tests added (794 total, up from 774).

**Implementation:**
- `undo_history` table added to `app/database.py` (DB-11)
- `app/undo.py` — `UndoManager` class with command-pattern record/undo/redo,
  50-level depth, DB persistence, state dispatch via `update_constants_batch`
- `POST /api/undo` and `POST /api/redo` endpoints in `app/server.py` (API-30–31)
- All three constant mutation endpoints (`PUT /api/constants/<name>`,
  `PUT /api/constants/batch`, `POST /api/constants/reset`) wrapped with
  before/after state capture and undo recording
- `Ctrl+Z` / `Ctrl+Shift+Z` keyboard shortcuts in `app.js` (UNDO-1–2)
- `constants_changed` SSE event added for real-time constant table refresh
- Cross-type undo verified (constant_update, constant_batch, constant_reset)

**New files:** `app/undo.py`, `tests/test_zapp_undo.py`
**Modified:** `app/database.py`, `app/server.py`, `app/static/js/app.js`

**Dependencies:** None.

---

## Phase 3 — Elements & Doors Database (Complete)

**Goal:** Extend the schema to represent interior elements and doors as
first-class database objects, with full CRUD APIs.

**Requirements:** DB-9–10, API-20–22, API-24–29, DT-9–11, SEL-7–8, RT-5 (17 reqs)
*(API-23 — element move — is assigned to Phase 4.  UI-7 — Elements tab —
is counted in Phase 0 as partial; this phase completes it.)*

**Completed.** 45 new tests added (839 total, up from 794).

**Implementation:**
- `elements` and `doors` tables added to `app/database.py` (DB-9, DB-10)
- 13 IW records seeded; 7 RO door records seeded with widths from constants
- `app/elements.py` — IW→constant mapping, variant-aware queries, hosted
  openings lookup
- `app/doors.py` — door validation (hinge side, swing direction, door type)
- `app/undo.py` — extended `_apply()` with 9 new action types:
  element_create/delete/update, door_create/delete/update,
  opening_create/delete/update
- 11 new API endpoints in `app/server.py`: elements CRUD (4), openings
  CRUD (3), doors CRUD (4), all with undo recording + SSE broadcast
- Elements tab in `index.html` with interior walls table (DT-9) and
  furniture/appliance table (DT-11)
- `app.js`: `loadElements()`, `updateElementsTable()`, door property
  editing (SEL-7), furniture/appliance properties (SEL-8),
  `element_changed` SSE listener (RT-5)

**New files:** `app/elements.py`, `app/doors.py`, `tests/test_zapp_elements.py`,
`tests/test_zapp_doors.py`
**Modified:** `app/database.py`, `app/server.py`, `app/undo.py`,
`app/templates/index.html`, `app/static/js/app.js`

**Deviation from estimate:** 45 new tests vs. estimated ~38.  Additional
tests cover element business logic (variant filtering, constant mapping,
hosted openings) and door validation.

**Key design — dual-source element model (transitional):**

During Phases 3–11, two categories of elements coexist:

1. **Engine-computed elements** — the 13 interior walls, furniture,
   appliances, and fixtures are computed by the engine from constants and
   layout logic per-variant.  Moving a base element (e.g., IW1) means
   changing its controlling constant via a hand-maintained mapping table
   in `app/elements.py` (e.g., `"IW1"` → `IW1_OFFSET_FROM_W9`).

2. **User-added (custom) elements** — stored in the `elements` table
   with absolute positions in `properties` JSON.  These are overlaid on
   the engine output during canvas rendering.

**Variant membership:** Custom elements have a nullable `variant` column.
`NULL` means the element appears in all variants.  A non-null value
(e.g., `"standard"`) restricts it to that variant only.

**Canvas merging:** The Elements panel (DT-11) and canvas merge both
sources — engine-computed items are shown as read-only rows (not directly
editable, but respond to constant changes), while custom items are shown
as fully editable rows.  Custom items display a visual indicator (e.g.,
badge or icon) to distinguish them from engine-computed items.

**Conflict avoidance:** If a custom element shares a name with an
engine-computed item (e.g., user adds a custom "stove"), the custom item
takes precedence and the engine-computed item is hidden for that variant.

**Migration to Phase 12:** This dual-source model is transitional.  At
Phase 12 (cutover), all engine-computed items are migrated to
database-stored elements with parametric formulas, eliminating the
distinction.  See CHARTER.md § Design Principle 5.

**Dependencies:** Phase 2 (undo recording wraps all mutations).

---

## Phase 4 — Move Tool (Complete)

**Goal:** Enable drag-to-reposition for walls, furniture, and openings.

**Requirements:** TL-5–10, API-23 (7 reqs)

**Completed.** 20 new tests added (862 total, up from 842).

**Implementation:**
- `app/elements.py` — fixed `IW_CONSTANT_MAP["IW9"]` (was `IW9_OFFSET_O10`,
  now `IW3_OFFSET_IW9`), added `IW_MOVE_AXIS` dict mapping each IW wall to
  its move axis and sign, added `compute_constant_delta(iw_name, dx, dy)`
- `app/undo.py` — added `element_move` action type to `_apply()` dispatch,
  handles both `move_type="constant"` and `move_type="position"` states
- `app/server.py` — `POST /api/elements/<id>/move` endpoint: IW walls
  translate dx/dy to constant changes; custom elements update offset
  properties; anchor-based format supported
- `app/static/js/tools.js` — `MoveTool` state machine with drag, ghost
  preview (TL-5), axis constraint for IW walls, shift-constrain (TL-7),
  grid snap to 1" (TL-8), group move from multi-selection (TL-9), offset
  dialog via Enter key (TL-6)
- `app/static/js/dialogs.js` — generic modal dialog framework with
  `Dialog.show()`, `Dialog.close()`, `parseOffsetString()` for text input
  (e.g. "6in east", "2ft north")
- `app/static/js/app.js` — multi-select (Ctrl+Click toggles selections,
  SEL-4 subset for TL-9), mouse handler hooks for move tool, Escape
  cancels drag, furniture override rendering merge (applies offset_x/
  offset_y from DB override elements to engine-computed variant items)
- `app/static/css/app.css` — `.move-ghost`, `.multi-selected`,
  `.dialog-overlay/.dialog` styles (Catppuccin palette)
- `app/templates/index.html` — script tags for dialogs.js and tools.js

**Bug fix:** `IW_CONSTANT_MAP["IW9"]` pointed to `"IW9_OFFSET_O10"` which
is defined in `floorplan/constants.py` but never used in `layout.py`. The
actual controlling constant is `IW3_OFFSET_IW9` (layout.py:160).

**New files:** `app/static/js/tools.js`, `app/static/js/dialogs.js`,
`tests/test_zapp_move.py`
**Modified:** `app/elements.py`, `app/undo.py`, `app/server.py`,
`app/static/js/app.js`, `app/static/css/app.css`, `app/templates/index.html`

**Dependencies:** Phase 3 (element CRUD, constant-to-element mapping).

---

## Phase 5 — Outline Chain Editing (Complete)

**Goal:** Make the building outline editable — change arc radii, segment
lengths, add/remove points — with automatic closure re-solving.

**Requirements:** ENG-11, API-16–19, OE-1–3, DT-2–4 (11 reqs)

**Actual:** 37 new tests, 899 total (310 app + 589 pre-existing).  All 11
requirements implemented.  OE-3 (constraint dialog) deferred to future phase
as the solver infrastructure is in place but the UI dialog is not yet needed.

**Work:**
- Create `app/outline_solver.py` — parallel closure solver reimplemented from
  `floorplan/geometry.py`'s `_chain_offset` logic (pure math, no import from
  `floorplan/`).  This is the most architecturally significant new module.
  - `chain_offset()` — walk chain entries from origin
  - `solve_closure()` — solve for `d_F2_F5`, `d_F18_F1`, and `sweep_closure`
    (3 variables: two line distances plus the closure arc's sweep angle for
    full positional + angular closure)
  - `validate_chain()` — check proposed modifications
  - `solve_for_constraint()` — secant method solver for target distances
- Add mutation APIs: `PUT /api/outline/<seq>` (API-16),
  `POST /api/outline/validate` (API-17), `POST /api/outline/add-point` (API-18),
  `DELETE /api/outline/<seq>` (API-19)
- Editable outline table cells (DT-2), closure indicator (DT-3), add/remove
  buttons (DT-4)
- F-point dragging on canvas (OE-1), arc radius handle (OE-2), constraint-based
  editing dialog (OE-3)
- Engine integration: solver writes modified chain to DB, injects solved
  distances as constants before `compute_geometry()` (ENG-11)
- All chain parameter types become editable: distances, bearings, arc radii,
  sweep angles — both graphically (canvas handles) and textually (table cells)

**New files:** `app/outline_solver.py`, `app/static/js/outline-editor.js`,
`tests/test_zapp_outline.py`

**Chain authority:** After Phase 5, the `outline_chain` database table is
authoritative for all chain parameters.  The engine reads chain data from the
DB and injects it (including solved closure distances) as patched constants
before reloading floorplan modules.  The floorplan module's hardcoded chain
definition is bypassed at runtime but retained as the seed source for
"Reset to Defaults" (which resets both constants and chain).

**Two solvers:** `floorplan/geometry.py`'s module-scope solver and
`app/outline_solver.py` must produce bit-identical results for **default**
chain parameters.  Cross-validation tests enforce this.  When the user
modifies chain parameters beyond defaults, the app solver is authoritative
and d² regression tests are expected to fail for the modified values (they
still pass for default values after reset).

**Dependencies:** Phase 2 (undo), Phase 3 (element CRUD for cascading effects).

---

## Phase 6 — Enhanced Canvas Rendering (Complete)

**Goal:** Add door arcs, clearance zones, and display toggles to the
interactive canvas.

**Requirements:** CV-7, CV-11, DIS-9–10, SEL-11, DOOR-1–4 (9 reqs)

**Completed.** 28 new tests added (927 total, up from 899).  All 9
requirements implemented.

**Implementation:**
- `app/engine.py` — `_compute_door_arcs()` computes structural door arcs
  from DB door data (hinge/swing resolution via wall unit vectors);
  `_compute_clearance_zones()` computes clearance polygons from variant item
  metadata (face vertices + distance); `_compute_appliance_doors()` computes
  appliance door arcs from variant item door metadata (hinge vertex, open/closed
  direction vectors).  All three integrated into `compute_geometry()` with
  `doors_data` parameter.
- `app/database.py` — O3 and O6 added to `_DOOR_SEED` (9 doors total).
  Door seeding changed from `INSERT OR REPLACE` to `INSERT OR IGNORE` to
  preserve user edits across app restarts.
- `app/variants.py` — Door metadata (hinge_idx, width, open_dir, closed_dir)
  added to dryer, washer, fridge, microwave.  Clearance metadata (face,
  distance) added to stove, dishwasher, hamper.  Stacked flag propagated
  for SVG paint-order correctness.
- `app/server.py` — `_get_geometry()` loads doors and passes `doors_data`
  to `compute_geometry()`
- `app/static/js/app.js` — `renderDoors()` renders structural door arcs
  (hinge circle + line + dashed arc polyline) and non-stacked appliance
  door arcs; `renderClearanceZones()` renders dashed clearance polygons;
  stacked appliance doors rendered in `renderFurniture()` above counter
  items.  `renderApplDoor()` shared helper eliminates rendering duplication.
  `itemOverrides()` computed once per render cycle and passed to all
  offset-aware render functions.
- `app/templates/index.html` — Doors and Clearance toggle checkboxes,
  `#layer-doors` and `#layer-clearance` SVG groups
- `app/static/css/app.css` — `.door-line`, `.door-arc`, `.door-hinge`,
  `.appl-door-line`, `.appl-door-arc`, `.clearance-zone` styles

**New files:** `tests/test_zapp_canvas.py`
**Modified:** `app/engine.py`, `app/database.py`, `app/variants.py`,
`app/server.py`, `app/static/js/app.js`, `app/static/css/app.css`,
`app/templates/index.html`

**Key design — server-side door arc computation:**
Door arc geometry (hinge position, open-tip position, arc polyline with
21 points) is computed server-side in `engine.py` and returned as part of
the geometry JSON.  The client just renders the pre-computed points.  This
follows the existing pattern: the server computes all geometry, the client
renders it.

**Key design — stacked appliance door paint order:**
Appliance doors flagged as `stacked: true` (e.g., microwave) are rendered
in `renderFurniture()` after all item polygons, rather than in
`renderDoors()`.  This ensures the door arc appears above the counter
polygon it sits on in SVG paint order.

**Dependencies:** Phase 3 (doors data), Phase 5 (outline chain editing).

---

## Phase 7 — Draw, Add, Delete, Rotate, and Shape Editor Tools (Complete)

**Goal:** Full element creation, deletion, and shape customisation from the
canvas.

**Requirements:** TL-15–27, SEL-4, SEL-10, DT-7 (16 reqs)
**Actual:** 16 requirements implemented, 27 new tests (954 total)

**Implemented:**
- **Opening width editing** (DT-7, SEL-10): `findWidthConstant()` maps opening
  names to their controlling constants (handles `_WIDTH`, `_HALF_WIDTH`, and
  IW-based RO patterns). Editable in both openings table and properties panel.
- **Delete tool** (TL-22, TL-23): Delete/Backspace key triggers
  `deleteSelectedElements()` with confirmation. Server-side cascade deletes
  hosted openings and doors when a wall is deleted. Client-side
  `IW_HOSTED_OPENINGS` map provides cascade warnings.
- **Multi-select** (SEL-4): Shift+Click added alongside existing Ctrl+Click.
  Rubber-band drag-select: mousedown on empty space starts selection rectangle,
  mouseup computes world-coordinate bbox and selects intersecting elements via
  SVG `getBBox()`.
- **Draw Wall tool** (TL-15–17): `DrawWallTool` state machine in tools.js.
  Two-click placement creates custom wall element (`source: "drawn"`) with
  `start/end/thickness/poly` properties. Preview line during draw. Editable
  thickness in properties panel recomputes wall polygon. `W` keyboard shortcut.
- **Add Element tools** (TL-18–20): `catalog.js` with hardcoded catalog data
  (8 furniture, 6 appliance, 3 fixture items). `showCatalog()` opens grid
  dialog, click enters placement mode, canvas click creates element via API.
  Placed elements rendered in `renderFurniture()` from `App.state.elements`.
- **Rotate tool** (TL-24): `R` key opens rotation dialog with angle input and
  preset buttons (0/90/180/270). `rotatedRectPoly()` recomputes polygon from
  center/width/depth/angle. Works for placed and drawn elements.
- **Shape editor** (TL-25–27): `shape-editor.js` with interactive SVG vertex
  editing, drag handles, vertex add/remove. Shape CRUD API (`GET/POST/PUT
  /api/shapes`). Shape assignment dropdown in properties panel for placed
  elements. SVG import: `parseSvgPolygon()` parses `<polygon>` and simple
  `<path>` elements.

**New files:** `app/static/js/catalog.js` (150 lines), `app/static/js/shape-editor.js`
(280 lines), `tests/test_zapp_tools.py` (27 tests)
**Modified:** `app/database.py` (shape CRUD: create/update/delete),
`app/server.py` (3 shape API routes), `app/static/js/app.js` (+400 lines:
keyboard shortcuts, delete, rotate, rubber-band, custom element rendering,
findWidthConstant, bboxFromPoly), `app/static/js/tools.js` (DrawWallTool,
IW_HOSTED_OPENINGS), `app/static/js/dialogs.js` (customContent, presetButtons),
`app/static/css/app.css` (rubber-band, draw-preview, catalog, preset styles),
`app/templates/index.html` (Wall button, Add section, script tags)

**Post-phase fixes:**
- Fixed race condition: `loadElements()` now awaits before first
  `loadGeometry()` call so item overrides are applied on initial render.
  SSE `element_changed` handler now re-renders the canvas.
- Fixed Reset to Defaults: now clears override elements (move offsets),
  placed items, and drawn walls via `reset_elements()`.  Elements and doors
  snapshots included in undo state for full_reset.  Client reloads elements
  and geometry after reset (no page reload required).

**Deferred to future phases:**
- TL-17 endpoint drag handles (infrastructure exists but interactive handles
  not yet rendered on drawn wall selection)
- TL-21 Add Opening tool (click-on-wall placement with segment detection)

---

## Phase 8 — Labels, Dimensions, and Annotations ✓

**Goal:** User-placeable room labels, custom dimension lines, and annotation
editing.  Unifies the room label system under the `elements` table.

**Requirements:** TL-11–14, LABEL-1–4, DIS-7, DIM-1–5, VAR-1–3, SEL-13–14 (19 reqs) — all implemented

**Completed.** 90 new tests (1044 total, up from 954).  27 base tests in
`tests/test_zapp_labels.py`, 63 additional tests across existing test files
for post-phase enhancements.

Post-phase fixes:
- Fixed variant exclusion seeding on DB upgrade path (pre-existing issue from Phase 7).

Post-phase enhancements:
- Anchored dimensions: snap-to-geometry with point, wall_face, opening_face anchor types
- Unified dimensions: all 22 builtin dimensions converted to seeded anchored elements with line_intersection anchor type
- Dimension selectability: click to select, properties panel shows source/anchors/style
- Per-dimension style: solid/dashed visual style control
- Multi-layout visibility: properties.variants checkbox list for per-element layout control
- Cache invalidation: element/opening/door CRUD endpoints now invalidate geometry cache
- Delete button in properties panel

**Work:**
- Create `app/labels.py` — label and dimension line management using the
  `elements` table (type `'label'` or `'dimension'`)
- **Dimension Line tool** (TL-11): click two points to create dimension with
  extension lines, arrowheads, and ft-in label
- Dimension line repositioning (TL-12), label rotation (TL-13), deletion (TL-14)
- **Room Label tool** (LABEL-1): click to place, type label text
- Label dragging (LABEL-2), inline text editing (LABEL-3), font size (LABEL-4)
- Dimensions toggle (DIS-7): controls user-created persistent dimension
  annotations (distinct from the existing computed dimensions toggle in CV-9,
  which controls engine-computed dimension lines)

**Room label model:** Phase 8 unifies room labels under the `elements`
table, superseding the `room_label_offsets` table from Phase 0.  Each room
label is stored as an element with type `'label'`.  The model:

1. **Default position:** the area-weighted centroid of the room's boundary
   polygon (auto-computed from room geometry, not stored).
2. **Stored state:** offset from centroid `(de, dn)` and rotation in degrees.
   Defaults are `(0, 0)` offset and `0°` rotation.
3. **Seeding / reset:** "Reset to Defaults" back-computes label positions
   from the existing SVG generation scripts' output and stores the
   corresponding offsets and rotations.  At initial startup (no DB), the
   same seed process populates default label elements.
4. **Migration:** the `room_label_offsets` table is deprecated.  Existing
   offsets are migrated to label elements on schema upgrade.  The 11
   auto-computed room labels become editable label elements; users can
   also add additional custom labels.

**New files:** `app/labels.py`, `tests/test_zapp_labels.py`

**Dependencies:** Phase 3 (elements table), Phase 7 (draw tool patterns).

**Deferred dimension editing features:**
- TL-17D endpoint drag handles for dimensions: render draggable handles at
  start/end points; on drag-release, snap to nearest geometry point and
  re-anchor the endpoint (updating `start_anchor`/`end_anchor` properties)
- Editable Start/End coordinate fields in properties panel: text inputs with
  ft-in parsing that update absolute coordinates (detaching any existing
  anchor on that endpoint)
- Reattach anchor: right-click endpoint handle → "Attach to nearest" snaps
  the endpoint to geometry and sets an anchor (complement to existing
  "Detach" context menu items)

---

## Phase 9 — Styling and Product Links (Complete)

**Goal:** Per-element visual customisation and product URL attachments.

**Requirements:** STYLE-1–4, LINK-1–2, SEL-12, CV-12 (8 reqs) — all implemented

**Implemented:**
- `app/style.py` — TYPE_DEFAULTS (7 element types), validation (colour, opacity,
  stroke), resolve_style() with 3-layer merge (defaults → base → view overrides)
- Properties panel: colour picker (STYLE-1), stroke colour/style (STYLE-2),
  opacity range slider (STYLE-3), per-view override section (STYLE-4)
- Product URL text field with "Open" button (SEL-12, LINK-1)
- SVG `<a xlink:href>` wrapping on interactive canvas (LINK-2)
- Canvas link icon overlay at item corner (CV-12)
- Inline style application in renderFurniture(), renderUserDimensions(),
  renderUserLabels() — overrides CSS class defaults
- Reset buttons to revert individual style properties to defaults
- Override auto-creation for engine-computed items (same pattern as move tool)

**New files:** `app/style.py`, `tests/test_zapp_style.py` (69 tests)

**Dependencies:** Phase 3 (elements), Phase 6 (canvas rendering).

### Post-phase enhancements
- None yet.  Style system designed for extensibility — new element types can
  add entries to TYPE_DEFAULTS and get full style editing automatically.

---

## Phase 10 — Domain Views: Site Plan, 3D, Analysis, Plumbing

**Goal:** Interactive editing for specialised views, wrapping existing
generators via the regeneration API.

**Requirements:** SITE-1–4, SCAD-1–3, ANALYSIS-1–3, PLUMB-1–8 (18 reqs)

**Work:**

### Site Plan (SITE-1–4) ✅ Phase 10c — Complete
- Setback distance inputs stored in config table, editable in properties panel (SITE-1) ✅
- Drainfield element tool — custom elements with type `'site_element'` via File menu (SITE-2) ✅
- Text annotation tool — `'site_annotation'` elements via File menu (SITE-3) ✅
- P-series survey points with distances in properties panel (SITE-4) ✅
- PDF view tabs for site_plan_df, site_plan_fs, and 3views
- PNG view tabs for 3D renders (flat roof, 2:12 roof)
- `compute_survey_points()` in engine.py, `/api/survey-points` endpoint
- `POST /api/generate-site-plan` endpoint with setback config response
- DB error handling: validation, `/api/reset-database`, error banner UI
- 30 new tests in `test_zapp_site.py` (18) and `test_zapp_scad.py` (+12)

### 3D / SCAD (SCAD-1–3) ✅ Phase 10b — Complete
- `POST /api/generate-3d` endpoint wrapping SCAD generator (SCAD-1) ✅
- Roof style selection stored in `config` table (SCAD-2) ✅
- Multi-view PDF generation via `POST /api/generate-views` (SCAD-3) ✅
- Config table for app-level settings, File menu UI with dropdown
- 15 new tests in `test_zapp_scad.py`

### Analysis (ANALYSIS-1–3) ✅ Phase 10a — Complete
- Span hover tooltip on SVG view (ANALYSIS-1) ✅
- Span rotation min/max in properties panel (ANALYSIS-2) ✅
- Room area display for all variants with toggle (ANALYSIS-3) ✅
- 14 tests in `test_zapp_analysis.py`

### Plumbing Layout (PLUMB-1–11) ✅ Phase 10d — Complete
- Plumbing interactive canvas with building ghost overlay, zoom/pan/selection (PLUMB-1) ✅
- Supply line drawing tool with hot (red) / cold (blue) polylines (PLUMB-2) ✅
- Fixtures/supplies table in Plumbing panel tab (PLUMB-3) ✅
- `plumbing_elements` table: pipes, fittings, fixture connections (PLUMB-4) ✅
- CRUD API: GET/POST/PUT/DELETE `/api/plumbing`, SSE broadcast, undo/redo (PLUMB-5) ✅
- Drain line drawing tool with brown polylines and slope annotations (PLUMB-6) ✅
- Fitting placement: tee, elbow, valve with rotation (PLUMB-7/8) ✅
- Fixture property editing: select fixture to edit cold/hot/drain flags (PLUMB-9) ✅
- Add/remove fixture connections via Place Fixture tool (PLUMB-10) ✅
- Fixture/pipe properties panel with name, flags, position, delete (PLUMB-11) ✅
- 11 seeded fixture connections with cold/hot/drain flags
- Plumbing tool palette: Cold, Hot, Drain, Fitting, Fixture (visible only in plumbing_edit view)
- `isCanvasView()` helper replaces hardcoded `=== "interactive"` checks
- Reference plumbing configuration seeded into database (6 pipes, 11 fixture positions)
- Database reset re-seeds reference plumbing automatically
- Move tool works for plumbing elements in plumbing_edit view
- Buried pipes rendered as dashed lines (`buried` property)
- Plumbing undo/redo updates canvas display
- Tool legend colors driven from single `PLUMBING_COLORS` constant
- Display toggles work correctly in plumbing_edit view
- Appliance door arcs respect Furniture toggle
- 32 tests in `test_zapp_plumbing.py`

**Further development (beyond Phase 10):** The plumbing layout will expand to
include service location indicators (hot, cold, and drain at each fixture),
full drain line routing from fixtures through the building to the exterior,
and site-level plumbing elements: septic tank placement, drainfield layout,
and the well-to-building supply run.  These extensions build on the Phase 10
foundation and will be specified as additional PLUMB requirements when the
base layout is complete.

**New files:** `app/plumbing.py`, `tests/test_zapp_plumbing.py`

**Dependencies:** Phase 3 (elements table), Phase 7 (draw/add tools).

---

## Phase 11 — View Variants and Polish ✅ Complete

**Goal:** User-defined layout variants with per-view layer configuration.

**Requirements:** UI-5–6, SEL-8a (3 reqs — all implemented)

**Implemented:**
- `variants` table in database (id, name, label, source_variant, flags,
  layer_config JSON, is_builtin flag).  5 built-in variants seeded.
- Phase 11a — Database-driven variants + per-variant layer config:
  - Variant definitions moved from hardcoded `VARIANTS` dict to DB
  - GET/PUT /api/variants endpoints for reading and updating layer_config
  - Dynamic variant dropdown (populated from DB, replaces hardcoded HTML)
  - Layer toggle state persisted per variant (saveCurrentLayerConfig /
    restoreLayerConfig cycle on variant switch)
  - `get_variant_flags()` reads from DB with dict fallback
  - `variant_update` undo action type
- Phase 11b — User-defined variant creation/deletion:
  - POST /api/variants creates named variant cloned from source
  - DELETE /api/variants/<id> removes user variant (built-in protected)
  - `clone_variant_exclusions()` / `clone_variant_elements()` copy source
    data to new variant on creation
  - `unclone_variant_elements()` / `delete_variant_exclusions()` clean up
    on deletion
  - "+" / "×" buttons in variant selector UI
  - `variant_create` / `variant_delete` undo action types (bidirectional)
  - SVG suffix fallback: user variants use source variant's floorplan SVG
- SEL-8a — Furniture property editing:
  - Width/Depth editable for placed elements in Properties panel
  - `handleElementPropEdit()` updates properties via PUT API with poly
    recomputation
  - Display uses stored dimensions (not bbox) for correct rotated values

**Tests:** 44 new variant tests (79 total in test_zapp_variants.py).
1516 tests pass.

---

## Phase 12 — Parametric Dependencies and Cutover (Charter Principle 5)

**Goal:** The charter's most ambitious feature and the project's end state —
every element's position stored as an editable formula referencing other
elements, with dependency chains displayed both as a formula table and
graphically highlighted in the drawing.  This phase also performs cutover:
NF-4 is lifted, code duplication is consolidated, and the database becomes
the sole authoritative source for all design data.

This is the largest architectural addition and should be designed as a
standalone specification before implementation begins.  Given its scope, it
will likely be implemented as multiple sub-phases (12a, 12b, 12c...) once
the design spec is finalised.

**Requirements:** Design specification for Charter Principle 5

SEL-15 (Constant Dependency Highlighting) was the first concrete requirement
for this phase and is now implemented: when a constant is focused in the
Properties panel, all geometry elements whose position depends on that constant
are highlighted on the canvas (first-order in white, downstream in pink).

**Architectural requirements:**

1. **Formula storage** — Each element's positioning parameters stored as
   expressions referencing other elements and constants (e.g., "IW3 spans
   from IW2.east_face to IW1.north_face")
2. **Dependency graph** — Directed acyclic graph of element references,
   persisted in the database
3. **Topological evaluator** — Sorts elements by dependency order, evaluates
   formulas to produce concrete positions
4. **Lock/unlock** — Any parametric dependency can be "locked" to freeze it
   at its current computed value, converting it to a fixed position so
   upstream changes no longer propagate through it
5. **Dependency chain UI** — When an element is selected, its dependency
   chain is displayed as a formula table and graphically highlighted on the
   canvas
6. **Design-first workflow** — The system supports building a design from
   any starting point (e.g., place a bed, then define walls relative to it;
   or define the exterior shell first and position everything inward)

**Cutover tasks:**

1. All engine-computed elements (interior walls, furniture, appliances,
   fixtures) are migrated to database-stored elements with parametric
   formulas, eliminating the Phase 3–11 dual-source model
2. Constants move from Python module attributes to database-stored values;
   the `constants` table evolves from a flat value store to a formula-aware
   parameter store
3. The 24 duplicated dimension constants in `app/variants.py` and the ~700
   lines of replicated positioning math are consolidated — the database
   formulas become the single source
4. All hardcoded metadata moves to the database: product URLs (currently
   in `variants.py`), style defaults (currently in `style.py` and
   `app.js`), door/clearance metadata, item dimensions, and variant
   configurations.  No design content remains in code
5. NF-4 is lifted; existing generator scripts are retained as seed sources
   for "Reset to Defaults" (which regenerates the entire database from
   their output — positions, metadata, URLs, styles, and all)
6. All five built-in variants and user-defined variants render correctly
   from database-driven formulas
7. d² regression tests pass for default database values (matching original
   script output)

**Impact:** Replaces the current procedural Python computation with a
data-driven evaluation model.  The positioning logic and all design
metadata currently baked into `floorplan/layout.py`, `floorplan/
gen_floorplan.py`, and `app/variants.py` are encoded as database-stored
records and formulas.  After cutover, the existing scripts are no longer
the authoritative source for any design content — they are seed-only.

**Dependencies:** All prior phases (the formula system builds on the element
CRUD, canvas rendering, and property editing infrastructure).

### Phase 12 Sub-Phase Status

| Sub-Phase | Status | Summary |
|-----------|--------|---------|
| **12a** | Complete | FormulaEvaluator, schema (`element_formulas`, `formula_deps`), 24 variant constants in DB, topo sort, cycle detection, dependency extraction, wall_rect/item_rect/item_circle evaluators |
| **12b** | Complete | 7 formula REST API endpoints (GET/PUT/DELETE formulas, lock/unlock, deps, dependents), properties panel formula section |
| **12c** | Complete | Evaluator extensions (`four_corner`, `proj`/`dist`/`neg`/`add`/`sub`/`mul` length specs, `neg`/`perp` dir specs), all 13 IW wall formulas written and verified to 1e-9 ft vs procedural, seeded into DB, hybrid engine active |
| **12d** | Complete | `wall_opening` formula type (4 positioning modes, 4 poly_order options), 5 layout item formulas (DRYER, WASHER, COUNTER, DRESSER, SHELVES), 12 outer opening formulas (O1-O11, O8a), 7 rough opening formulas (RO1-RO7), all 37 formulas seeded in DB |
| **12e** | Complete | 5 new formula types (`toilet_shape`, `bath_sink_shape`, `dining_triangle`, `dining_chair`, `ellipse_rect`), new specs (`element_centroid`, `ray_circle_isect`, `rotated`, `radius_key`), ~50 variant item formulas across 3 variants, hybrid engine overrides variant items preserving metadata |
| **12f** | Complete | Lock/unlock UI (padlock icon on canvas, lock/unlock button in properties panel), formula dependency highlighting (blue upstream, orange downstream on element select), lock/unlock undo/redo, `formula_locked` SSE event, `GET /api/deps/graph` endpoint (full DAG with nodes and edges), `locked_elements` in geometry response |
| **12g** | Complete | Cutover: removed `patch_constants()`/`importlib.reload()` from `compute_geometry()`, added outer/rough opening formula overrides, BED formula (with `div` length spec), lifted NF-4, updated isolation tests |
| **12h** | Complete | Eliminated procedural element baselines: seeded element metadata into DB (~59 elements), rewrote `compute_geometry()` as formula-only (removed `compute_interior_layout`/`compute_outer_openings`/`compute_rough_openings`/`compute_variant_items` calls), moved product URLs to DB element properties, removed 835 lines of dead code from `variants.py` (984→149 lines), clone variant formulas on custom variant creation |

---

## Phase 14 — Fully DB-Driven Outline & Export/Import

**Goal:** Make the entire building definition — outline origin, survey
traverse, and span analysis — fully database-driven and editable, then
provide a single-file export/import for the complete project state.

Currently the outline chain segments are DB-driven, but several upstream
inputs remain hardcoded: the F2 origin position, survey traverse legs, and
span analysis relies on `patch_constants()` + module reloading.  Phase 14
eliminates these last procedural holdouts.

**Note:** Phase 14 is expected to be implemented before Phase 13.

### Sub-phases

| Sub-phase | Status | Summary |
|-----------|--------|---------|
| **14-A** | Complete | F2 origin in DB: `F2_EASTING` / `F2_NORTHING` constants (category "geometry"), read via `constants_dict.get()` in `compute_geometry()` |
| **14-B** | Complete | Survey data in DB: `survey_legs` and `survey_config` tables seeded from `shared/survey.py`, `_compute_traverse_from_db()` reads from DB, 5 API endpoints |
| **14-C** | Complete | Span analysis without module patching: `compute_span_data()` and `compute_span_rotation()` use `compute_geometry()` result directly, no `patch_constants()` calls |
| **14-D** | Complete | Full export/import: `GET /api/project/export` returns JSON (13 data sections), `POST /api/project/import` validates closure + DAG cycles, undo support |

### Phase 14-A: F2 Origin in DB

Add two constants (`F2_EASTING` = −18.5, `F2_NORTHING` = −13.5, category
"outline", unit "ft") to the `constants` table.  Replace hardcoded values
in `compute_geometry()` with lookups from `constants_dict`.  The outline
solver `walk_chain()` already accepts F2_E / F2_N as parameters.

### Phase 14-B: Survey Data in DB

New tables:

```sql
CREATE TABLE IF NOT EXISTS survey_legs (
    seq       INTEGER PRIMARY KEY,
    bearing   REAL NOT NULL,      -- degrees
    distance  REAL NOT NULL,      -- feet
    label     TEXT                -- POB, P2, P3, P4, P5
);

CREATE TABLE IF NOT EXISTS survey_config (
    key   TEXT PRIMARY KEY,
    value TEXT NOT NULL            -- JSON value
);
-- Keys: FC_IN_P3_E, FC_IN_P3_N, COORD_ROTATION
```

Seed from `shared/survey.py` hardcoded values.  `compute_traverse()` reads
from DB instead of hardcoded arrays.  `FC_IN_P3` and `COORD_ROTATION` from
`survey_config`.

### Phase 14-C: Span Analysis Without Module Patching

Remove last `patch_constants()` usage.  `compute_span_analysis()` and
`compute_3d_model()` extract measurements from `compute_geometry()` result
directly instead of patching floorplan modules and reloading.

### Phase 14-D: Full Export/Import

**Export** (`GET /api/project/export`): returns JSON containing all DB
tables — constants, outline_chain, survey_legs, survey_config, elements,
element_formulas, formula_deps, doors, variants, variant_exclusions,
plumbing_elements.

**Import** (`POST /api/project/import`): validates structure, outline
closure, formula DAG (no cycles), then replaces all data in a transaction
with undo support (before-state = full export of current).

**Reset** (`POST /api/reset`): extend to also reset survey data and
re-seed all metadata.

### Files affected

- `app/database.py` — new tables, seeding, bulk export/import ops
- `app/engine.py` — read F2 from constants, survey from DB, remove patching
- `app/server.py` — export/import endpoints
- `shared/survey.py` — remains as seed data source (read-only)

---

## Phase 15 — Generator Data Provider (Planned)

**Goal:** Build a data-access layer that provides SVG generators with all the
geometry they need from the database and engine, replacing direct imports from
hardcoded Python modules (`floorplan/layout.py`, `floorplan/openings.py`,
`floorplan/constants.py`).  This is the foundation for Phases 16–18.

**Status:** Planned

**Motivation:** Today every generator imports `build_floorplan_data()` and/or
individual `floorplan/` modules, which read hardcoded Python source.  Edits
made in the interactive editor (database) are invisible to generators.  This
phase creates a single entry point that generators call instead, backed by the
DB-stored constants, formulas, and outline chain.

**Scope:**

| Sub-phase | Summary |
|-----------|---------|
| **15-A** | `app/gen_provider.py` — data provider module exposing a `GeneratorData` object that wraps `compute_geometry()` output plus shell geometry, roof geometry, and survey data |
| **15-B** | Extend `compute_geometry()` (or add companion functions) to expose shell boundary paths (S/G-series), wall section enumeration, and U-turn arc polygons currently computed only inside `shared/wall_shells.py` |
| **15-C** | Extend to expose roof geometry (R-series corners, roof polyline, roof area) currently computed only inside `floorplan/roof.py` |
| **15-D** | Golden-gate identity tests: for a freshly-seeded (unedited) database, verify that `GeneratorData` returns values identical (within tolerance) to the hardcoded module outputs |

### Design: GeneratorData

```python
class GeneratorData:
    """All geometry a generator needs, sourced from the database."""
    # Outline
    points: dict           # F/W/C/P/TC/R-series named points
    outline_segs: list     # LineSeg/ArcSeg chain
    outline_poly: list     # [E, N] vertices
    inner_segs: list       # inner wall boundary segments
    inner_poly: list       # inner wall boundary polygon
    radii: dict            # R_a-series arc radii
    constants: dict        # all constants from DB

    # Interior
    interior_walls: dict   # IW name → poly/bbox
    outer_openings: list   # O1-O11 with polys
    rough_openings: list   # RO1-RO7 with polys
    variant_items: dict    # furniture/appliance/fixture geometry
    door_arcs: list        # structural door arcs
    appliance_doors: list  # appliance door arcs
    clearance_zones: list  # clearance zone polys
    room_labels: list      # room name/pos/area

    # Construction detail (Phase 15-B)
    shell_paths: dict      # S/G-series inset paths per shell
    wall_sections: list    # enumerated wall sections with openings
    uturn_polys: list      # U-turn arc polygon fill regions

    # Roof (Phase 15-C)
    roof_points: dict      # R-series corner points
    roof_poly: list        # roof outline polygon
    roof_area: float       # roof area in sq ft

    # Survey (requires Phase 14-B)
    survey_points: dict    # P-series traverse points
    survey_legs: list      # bearing/distance legs
```

### Files affected

- **New:** `app/gen_provider.py` — wraps engine + shell/roof/survey computation
- **New:** `tests/test_gen_provider.py` — identity tests vs. hardcoded modules
- **Modified:** `app/engine.py` — expose shell/roof helper functions
- **Modified:** `shared/wall_shells.py` — factor out reusable computation (if needed)

**Dependencies:** Phase 14 (DB-driven outline and survey data).

---

## Phase 16 — Floorplan & Roof Generator Migration (Planned)

**Goal:** Migrate `gen_floorplan.py` and `gen_roof.py` to read all data from
the database via `GeneratorData`, so that edits made in the interactive editor
are reflected in generated SVGs.  The initial (unedited seed) database must
produce SVG output identical to the current hardcoded generators.

**Status:** Planned

**Scope:**

| Sub-phase | Summary |
|-----------|---------|
| **16-A** | Refactor `gen_floorplan.py`: replace `build_floorplan_data()` and direct `floorplan/` imports with `GeneratorData` lookups.  Preserve all rendering logic (SVG paths, labels, dimension lines, colours) |
| **16-B** | Refactor `gen_roof.py`: replace `build_floorplan_data()` + `compute_roof_geometry()` imports with `GeneratorData.roof_*` fields |
| **16-C** | SVG identity tests: binary-diff (or pixel-diff) of seed-DB-generated SVGs against current hardcoded SVGs, tolerating only non-semantic whitespace differences |
| **16-D** | App integration: add "Regenerate SVG" button/endpoint that runs the migrated generator against the live database, so the Floorplan SVG tab reflects the current editor state |

### Key challenge

`gen_floorplan.py` is the master generator (~1000+ lines) with extensive
rendering logic.  The migration must change only the data source, not the
rendering.  A phased approach within 16-A is recommended: first replace the
data acquisition calls, then run identity tests, then clean up.

### Files affected

- **Modified:** `floorplan/gen_floorplan.py` — data source switch
- **Modified:** `roof/gen_roof.py` — data source switch
- **Modified:** `app/server.py` — regenerate endpoint (16-D)
- **New:** `tests/test_gen_identity.py` — SVG identity comparison tests

**Dependencies:** Phase 15 (GeneratorData provider).

---

## Phase 17 — Wall Detail Generator Migration (Planned)

**Goal:** Migrate `gen_walls.py` to read from the database via
`GeneratorData`, including the shell boundary paths and wall section
enumeration added in Phase 15-B.

**Status:** Planned

**Scope:**

| Sub-phase | Summary |
|-----------|---------|
| **17-A** | Refactor `gen_walls.py`: replace all `floorplan/` and `walls/constants.py` imports with `GeneratorData` lookups.  Shell geometry (S/G-series paths, U-turn arcs) comes from `GeneratorData.shell_paths` / `wall_sections` / `uturn_polys` |
| **17-B** | SVG identity tests for `walls.svg` and `all_walls.svg` |
| **17-C** | App integration: wall detail SVG regeneration from live DB state |

### Key challenge

The wall detail generator is the most construction-specific view, relying on
shell inset paths, U-turn arc polygons, and per-section opening cutouts.  All
of this geometry must be exposed through Phase 15-B before this phase can
begin.  The rendering logic itself (SVG path construction, colour scheme,
annotations) stays unchanged.

### Files affected

- **Modified:** `walls/gen_walls.py` — data source switch
- **Modified:** `walls/constants.py` — may become a thin re-export or removed
- **Extended:** `tests/test_gen_identity.py` — wall SVG identity tests

**Dependencies:** Phase 15-B (shell geometry in GeneratorData), Phase 16
(pattern established).

---

## Phase 18 — Remaining Generator Migration (Planned)

**Goal:** Migrate all remaining SVG/PDF generators to read from the database,
completing the transition to fully editable views.

**Status:** Planned

**Scope:**

| Sub-phase | Generator | Output files | Notes |
|-----------|-----------|-------------|-------|
| **18-A** | `span/gen_span.py`, `span/gen_span_minmax.py`, `span/gen_span_min.py` | `span/*.svg` | Span analysis reads interior wall positions and outline geometry; both now available from `GeneratorData` |
| **18-B** | `site/gen_site_plan.py` | `site/site_plan_*.pdf` | Requires Phase 14-B survey data in DB.  Site-specific calibration constants (PDF coordinates, property line distances) move to DB constants or a `site_config` table |
| **18-C** | `plumbing/gen_plumbing.py` | `plumbing/plumbing.svg` | Unify plumbing into the layout model: plumbing becomes a layout variant (selectable alongside standard/minik/daybed/bare/sf) rather than a separate top-level view tab.  The Plumbing SVG is generated from DB `plumbing_elements` data overlaid on the building ghost (same data the Plumbing Edit canvas already renders).  The standalone `_render_plumbing_path()` in gen_floorplan.py (~235 lines of hardcoded pipe routing) is deleted — the DB-seeded reference pipes replace it.  The current separate "Plumbing" SVG tab and "Plumbing Edit" canvas tab merge into a single Floorplan view with plumbing as the selected layout |
| **18-D** | `scad/gen_flat_roof.py`, `scad/gen_2in12.py`, `scad/gen_views.py`, `scad/gen_line_drawings.py` | `scad/*.scad`, `scad/*.svg` | 3D model generators; outline + roof + shell geometry from `GeneratorData`.  Wall heights and roof slope constants move to DB |
| **18-E** | `survey/gen_path_svg.py` | `survey/path_area.svg` | Survey path overlay; requires Phase 14-B survey data |
| **18-F** | SVG identity tests for all migrated generators |
| **18-G** | Procedural module sunset: mark `floorplan/layout.py`, `floorplan/openings.py`, and rendering helpers in `floorplan/gen_floorplan.py` as deprecated seed-only sources.  `gen_all.py` switches to DB-backed generation |

### End state

After Phase 18-G:
- `python gen_all.py` reads from the database (seeded on first run if empty)
- All SVG/PDF outputs reflect the current editor state
- Edits in the interactive editor are visible in every generated view
- The hardcoded `floorplan/` modules remain as seed data sources only
- A freshly-seeded database produces output identical to the pre-migration
  hardcoded generators

### Files affected

- **Modified:** All generator scripts in `span/`, `site/`, `plumbing/`, `scad/`, `survey/`
- **Modified:** `gen_all.py` — DB-backed generation pipeline
- **Extended:** `tests/test_gen_identity.py` — full identity suite
- **Deprecated:** `floorplan/layout.py`, `floorplan/openings.py` (seed-only)

**Dependencies:** Phase 14 (survey data in DB), Phase 15 (GeneratorData),
Phase 16 (migration pattern established).

---

## Phase 13 — Electrical Layout (Aspirational)

**Goal:** An interactive electrical layout, following the same pattern as the
plumbing layout (Phase 10): a full interactive canvas with database-stored
elements, CRUD APIs, and domain-specific drawing tools.

Unlike plumbing, no reference `electrical.svg` generator exists in the
codebase — the electrical layout will be built from scratch within the editor,
without a pre-existing generator to replicate or wrap.

**Scope (to be specified):**
- `electrical_elements` database table for circuits, wiring runs, devices
  (outlets, switches, lights, panels), and connections
- CRUD API endpoints for electrical elements
- Interactive canvas with building outline overlay showing electrical elements
- Circuit drawing tool: draw wiring runs between devices, assign to circuits
  with amperage and breaker size
- Device placement tools: outlets (standard, GFCI, 240V), switches (single,
  3-way, dimmer), light fixtures (ceiling, wall, recessed), smoke/CO detectors
- Panel schedule: breaker panel layout with circuit assignments, load
  calculations, and wire gauge recommendations
- Service entry: main panel location, meter base, service lateral routing
- Conduit and wire routing with bend radius constraints
- NEC compliance annotations (e.g., outlet spacing, GFCI requirements in
  wet locations, AFCI requirements in living spaces)

**New files (anticipated):** `app/electrical.py`,
`app/static/js/electrical-editor.js`, `tests/test_zapp_electrical.py`

**Dependencies:** Phase 10 (plumbing layout establishes the pattern for
domain-specific interactive layouts), Phase 3 (elements table), Phase 7
(draw/add tools).

---

## Phase Dependency Graph

```
Phase 1 (Test Coverage) ─────────────────────────────┐
                                                       │
Phase 2 (Undo/Redo) ─────────┐                        │
                              │                        │
Phase 3 (Elements & Doors) ◄──┘                        │
    │         │                                        │
    │         ├── Phase 4 (Move Tool)                  │
    │         │                                        │
    │         ├── Phase 5 (Outline Editing)             │
    │         │         :                              │
    │         └── Phase 6 (Canvas Rendering) ◄── 5?    │
    │                   │                              │
    │                   │                              │
    ├── Phase 7 (Draw/Add/Delete) ◄── 3, 4, 6         │
    │                                                  │
    ├── Phase 8 (Labels & Dimensions) ◄── 3, 7        │
    │                                                  │
    ├── Phase 9 (Styling & Links) ◄── 3, 6            │
    │                                                  │
    ├── Phase 10 (Domain Views) ◄── 3, 7              │
    │         │                                        │
    │         └── Phase 13 (Electrical) ◄── 10, 18     │
    │                                                  │
    ├── Phase 11 (View Variants & Polish) ◄── all     │
    │                                                  │
    └── Phase 12 (Parametric Dependencies) ◄── all ───┘
              │
              └── Phase 14 (DB-Driven Outline & Export) ◄── 12
                        │
                        └── Phase 15 (Generator Data Provider) ◄── 14
                                  │
                                  ├── Phase 16 (Floorplan & Roof Gen) ◄── 15
                                  │         │
                                  │         ├── Phase 17 (Wall Detail Gen) ◄── 15-B, 16
                                  │         │
                                  │         └── Phase 18 (Remaining Gens) ◄── 14, 15, 16
                                  │
                                  └── Phase 13 (Electrical) ◄── 10, 18
```

(`:` and `?` denote optional dependencies — Phase 6 can start without Phase 5
but room boundaries adapt to outline changes once Phase 5 is complete.)

Phases 1 and 2 can proceed in parallel.  Phases 4, 5, and 6 can proceed in
parallel after Phase 3.  Phase 14 follows Phase 12.  Phases 15–18 form a
sequential chain migrating SVG generators to DB-backed data.  Phase 13
(Electrical) follows Phase 18 (all generators DB-backed).

---

## Phase Completion Protocol

Each phase is considered complete only after:

1. All phase requirements pass automated tests
2. All 1516+ tests continue to pass (`python -m pytest tests/ -x -q`)
3. All SVGs regenerate successfully (`python gen_all.py`)
4. User acknowledgement that all phase goals are met with no known outstanding
   issues

**Before proceeding to the next phase**, perform these maintenance steps:

1. **Update this ROADMAP.md** — move the completed phase into a "Completed"
   section, record actual test count, note any deviations from the plan, and
   update the Current State summary at the top
2. **Update ARCHITECTURE.md** — reflect any new modules, tables, API endpoints,
   or architectural patterns introduced by the phase
3. **Update REQUIREMENTS.md** — mark implemented requirements (remove **(NEW)**
   tag, update acceptance results if applicable)
4. **Verify undo coverage (UNDO-4)** — if the phase introduced new mutation
   types (element add/move/delete, door changes, outline edits, label changes,
   styling, etc.), verify that each mutation type is covered by undo/redo
   tests.  This is required by UNDO-4 and applies to every phase from Phase 3
   onward.
5. **Commit the documentation update** as a separate commit after the phase
   implementation commit(s)

This ensures the documentation always reflects the true state of the project
and that fresh context windows have accurate information.

**Database policy during development:** Schema changes (e.g., new tables in
Phase 3) do not require migration logic.  During development, databases are
wiped and recreated from defaults when the schema changes.  Database migration
will become a concern only after initial development is complete and users have
persistent data worth preserving.

---

## Requirement-to-Phase Mapping

| Phase | Requirement IDs | Count |
|-------|----------------|-------|
| 0 (done) | DB-1–8, DB-12–13, ENG-1–10, ENG-12–15, API-1–15, API-32–34, GEN-1–4, UI-1–4, UI-7 (partial), UI-8, CV-1–6, CV-8 (partial), CV-9, CV-10, CV-10a, CV-13–17, DIS-1–6, SEL-1–3, SEL-5–6, SEL-9, TL-1–4, CT-1–10, CT-7a–j, DT-1, DT-5–6, DT-8, RT-1–4, APP-1–4, NF-1–6 | 121 |
| 1 (done) | (test coverage for above — 22 new tests) | 0 new reqs |
| 2 | DB-11, API-30–31, UNDO-1–4 | 7 |
| 3 | DB-9–10, API-20–22, API-24–29, DT-9–11, SEL-7–8, RT-5 | 17 |
| 4 | TL-5–10, API-23 | 7 |
| 5 | ENG-11, API-16–19, OE-1–3, DT-2–4 | 11 |
| 6 | CV-7, CV-11, DIS-9–10, SEL-11, DOOR-1–4 | 9 |
| 7 | TL-15–27, SEL-4, SEL-10, DT-7 | 16 |
| 8 (done) | TL-11–14, LABEL-1–4, DIS-7, DIM-1–5, VAR-1–3, SEL-13–14 | 19 |
| 9 | STYLE-1–4, LINK-1–2, SEL-12, CV-12 | 8 |
| 10 | SITE-1–4, SCAD-1–3, ANALYSIS-1–3, PLUMB-1–11 | 21 |
| 11 | UI-5–6, SEL-8a | 3 |
| 12 | SEL-15, Charter Principle 5 | 1 + design spec |
| 14 | (DB-driven outline, survey, export/import — to be specified) | TBD |
| 15 | (generator data provider — to be specified) | TBD |
| 16 | (floorplan & roof generator migration — to be specified) | TBD |
| 17 | (wall detail generator migration — to be specified) | TBD |
| 18 | (remaining generator migration — to be specified) | TBD |
| 13 | (aspirational — to be specified) | TBD |
| **Total** | | **240 + Phase 13–18 TBD** |

---

## New Files by Phase

Files already created during Phase 0 work: `app/apputil.py`,
`app/parse_input.py`, `tests/test_zapp_parse_input.py`.

| Phase | New Python files | New JS files | New test files |
|-------|-----------------|-------------|---------------|
| 1 (done) | — | — | — (extended existing test files) |
| 2 | undo.py | — | test_zapp_undo.py |
| 3 | elements.py, doors.py | — | test_zapp_elements.py, test_zapp_doors.py |
| 4 | — | tools.js, dialogs.js | test_zapp_move.py |
| 5 | outline_solver.py | outline-editor.js | test_zapp_outline.py |
| 6 (done) | — | — | test_zapp_canvas.py |
| 7 | — | selection.js | test_zapp_tools.py |
| 8 | labels.py | — | test_zapp_labels.py |
| 9 | style.py | — | test_zapp_style.py |
| 10 | plumbing.py | site-plan.js, plumbing-editor.js | test_zapp_plumbing.py, test_zapp_views.py |
| 11 | — | — | test_zapp_variants_ext.py |
| 12 | dependencies.py | dependency-graph.js | test_zapp_deps.py |
| 14 | survey_db.py | — | test_zapp_survey.py |
| 15 | gen_provider.py | — | test_gen_provider.py |
| 16 | — | — | test_gen_identity.py |
| 17 | — | — | (extends test_gen_identity.py) |
| 18 | — | — | (extends test_gen_identity.py) |
| 13 | electrical.py | electrical-editor.js | test_zapp_electrical.py |

---

## Estimated Test Growth

| Phase | New tests | Cumulative |
|-------|-----------|-----------|
| 0 | — | 752 |
| 1 (done) | 22 | 774 |
| 2 (done) | 20 | 794 |
| 3 (done) | 45 | 839 |
| 4 (done) | 20 | 862 |
| 5 (done) | 37 | 899 |
| 6 (done) | 28 | 927 |
| 7 (done) | 27 | 954 |
| 8 (done) | 90 | 1044 |
| 9 (done) | 17 | 1061 |
| 10a (done) | 14 | 1075 |
| 10b (done) | 15 | 1090 |
| 10c (done) | 30 | 1120 |
| 10 (actual) | 63 (across 10a-c + DB error handling) | 1183 |
| 11 (actual) | 76 | 1259 |
| 12 (actual) | 257 | 1516 |
| 14 | ~30 | ~1546 |
| 15 | ~25 | ~1571 |
| 16 | ~20 | ~1591 |
| 17 | ~10 | ~1601 |
| 18 | ~30 | ~1631 |
| 13 | ~20 | ~1651 |

---

## Anticipated Challenges

1. **Closure Solver Reimplementation (Phase 5)** — The `_chain_offset` function
   in `floorplan/geometry.py` runs at module import.  The parallel solver in
   `app/outline_solver.py` must produce bit-identical results.  Cross-validation
   tests between the two are critical.

2. **Constant-to-Element Mapping (Phase 3–4)** — Mapping a drag on IW1 back to
   `IW1_OFFSET_FROM_W9` requires a hand-maintained table in `app/elements.py`.
   This table must stay in sync with `floorplan/layout.py`'s computation logic.

3. **Module Reloading + Outline Mutation (Phase 5)** — After outline
   modification, the engine must inject the solver's computed `d_F2_F5` and
   `d_F18_F1` as patched constants before reloading modules.  The module-scope
   solver in `floorplan/geometry.py` will otherwise run with stale values.

4. **Frontend Module Splitting (Phase 4+)** — `app.js` is ~1440 lines (grew
   from 1100 with room labels, SVG zoom/pan, item labels, and SF rendering).
   Adding tools, canvas rendering, and selection will push past 3000 without
   splitting.  New JS files (`tools.js`, `canvas.js`, `selection.js`, etc.)
   share the global `App` namespace via additional `<script>` tags.

5. **Dual-Source Element Model (Phase 3–11)** — Engine-computed items and
   DB-stored custom items coexist on the canvas.  The merging logic must
   handle variant membership (custom elements can be all-variant or
   single-variant), name conflicts (custom items shadow engine-computed
   items), and clear UI indicators so users know which items are editable
   vs. computed.  This transitional model is eliminated at Phase 12
   cutover when all items become database-stored elements.

6. **Site Plan / Plumbing Drainfield Overlap (Phase 10)** — The drainfield
   appears in both the site plan (SITE-2, as a placement rectangle in the
   `elements` table) and the plumbing system (future PLUMB extensions, as
   a piped element in `plumbing_elements`).  These are two views of the
   same physical object.  The design should store the drainfield once (in
   `elements`) and render it in both views, so moving it in one view
   updates the other.

7. **Parametric Dependencies (Phase 12)** — Fundamentally replaces procedural
   computation with data-driven evaluation.  Requires its own design spec
   before implementation.  Given its scope (formula storage, DAG evaluation,
   element migration, cutover consolidation), it will likely be implemented
   as multiple sub-phases.  This is the largest planned change by far.

8. **SVG Identity Verification (Phases 15–18)** — The seed-database SVG
   output must be identical to the current hardcoded-module SVG output.
   Floating-point precision, SVG path discretization (arc-to-polyline
   step counts), and string formatting differences can cause false
   negatives.  Identity tests need a tolerance-aware comparison strategy
   (e.g., numeric path coordinate tolerance of ±0.001 ft, ignoring
   whitespace-only differences).

9. **Generator Rendering/Data Separation (Phase 15–16)** — Generators
   currently interleave data acquisition and SVG rendering in a single
   script.  The migration must cleanly separate "get data" (now from DB)
   from "render SVG" (unchanged).  Some generators (notably
   `gen_floorplan.py`) compute intermediate derived values inline during
   rendering; these must be identified and either moved into the engine
   or preserved as local rendering logic.

---

## Cutover Status (Phase 12h Complete)

The cutover criteria have been fully met:

1. ✅ All 258 requirements pass automated tests (~1500 tests total)
2. ✅ All element positions are database-stored JSON formulas evaluated
   by the FormulaEvaluator in topological dependency order
3. ✅ All 24 dimension constants consolidated in the database; derived
   constants (`WALL_EXTRA`, `CORNER_SW_R`) computed from dict values
4. ✅ `app/variants.py` reduced to 149 lines (registry + flags only);
   all procedural positioning math removed (835 lines deleted)
5. ✅ d² regression tests (281 golden-gate checks) pass unchanged
6. ✅ All five layout variants render from formula-driven geometry
7. ✅ NF-4 lifted; `compute_geometry()` no longer patches modules
8. ✅ All element metadata (labels, types, shapes, door/clearance configs,
   product URLs, variant lists) stored in `elements` table — no
   hardcoded metadata remains in code
9. ✅ `compute_geometry()` no longer imports from `floorplan/layout.py`
   or `floorplan/openings.py` — formula-only engine
10. ✅ Variant formula cloning on custom variant creation
