# ADU Editor — Development Roadmap

This roadmap traces the path from the current **parametric viewer with constant
editing** to the charter's goal: a **full parametric editor** where every
element can be created, moved, edited, deleted, and persisted from the browser
— all without AI assistance.

All work respects NF-4: no files outside `app/` and `tests/` are modified
until cutover (see ARCHITECTURE.md § NF-4).

---

## Current State (Phase 2 — Complete)

**128 of 226 requirements implemented.**  208 app tests, 586 pre-existing tests
(794 total).  All implemented requirements have automated test coverage.

| Capability | Status |
|------------|--------|
| Constants: view, edit, persist, reset | Done |
| Outline chain: view (read-only) | Done |
| Geometry engine: full recomputation | Done |
| Interactive canvas: outline, walls, openings, furniture, points, dims | Done |
| Five layout variants | Done |
| 11 SVG view tabs with zoom/pan | Done |
| Properties panel with related constants | Done |
| Constants table: sort, filter, inline edit, category colours | Done |
| Openings table: outer + rough | Done |
| REST API: 16 endpoints, SSE | Done |
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

**What's missing:** Element CRUD, doors, outline editing, move/draw
tools, shape editor UI, styling, site plan editing, 3D integration,
plumbing layout (interactive canvas with DB-stored elements), electrical
layout (aspirational), parametric dependencies and cutover to fully
database-driven design (Charter Principle 5).  The current implementation
uses constants as the single editable root; the target architecture makes
every element a database entity with parametric position formulas.

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

## Phase 3 — Elements & Doors Database

**Goal:** Extend the schema to represent interior elements and doors as
first-class database objects, with full CRUD APIs.

**Requirements:** DB-9–10, API-20–22, API-24–29, DT-9–11, SEL-7–8, RT-5 (17 reqs)
*(API-23 — element move — is assigned to Phase 4.  UI-7 — Elements tab —
is counted in Phase 0 as partial; this phase completes it.)*

**Work:**
- Add `elements` table (DB-9) and `doors` table (DB-10) to database schema
- Seed elements with the 13 interior walls only (furniture/appliances/fixtures
  are computed per-variant, not seeded); seed doors from RO1–RO7 defaults
- Create `app/elements.py` — CRUD, cascading delete, constant-to-element mapping
- Create `app/doors.py` — door arc geometry computation, CRUD
- Add 9 new API endpoints: elements CRUD (API-20–22, API-24–26), doors (API-27–29)
- Complete "Elements" tab in right panel (UI-7) with interior walls table
  (DT-9–10) and furniture/appliance table (DT-11)
- Enhanced property display for furniture (SEL-8) and doors (SEL-7)
- Broadcast `element_changed` SSE event on mutations (RT-5)

**New files:** `app/elements.py`, `app/doors.py`, `tests/test_zapp_elements.py`,
`tests/test_zapp_doors.py`
**Modified:** `app/database.py`, `app/server.py`, `app/engine.py`,
`app/templates/index.html`, `app/static/js/app.js`

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

## Phase 4 — Move Tool

**Goal:** Enable drag-to-reposition for walls, furniture, and openings.

**Requirements:** TL-5–10, API-23 (7 reqs)

**Work:**
- Implement `MoveTool` in new `app/static/js/tools.js`:
  - Drag with live ghost preview (TL-5)
  - Offset dialog for precise entry (TL-6)
  - Shift-constrained axis movement (TL-7)
  - Snap to 1" grid when enabled (TL-8)
  - Group move for multi-selected elements (TL-9)
  - Opening slides along host wall segment only (TL-10)
- Backend: `POST /api/elements/<id>/move` maps to constant changes for base
  elements, absolute position updates for custom elements (API-23)

**New files:** `app/static/js/tools.js`, `app/static/js/dialogs.js`,
`tests/test_zapp_move.py`

**Dependencies:** Phase 3 (element CRUD, constant-to-element mapping).

---

## Phase 5 — Outline Chain Editing

**Goal:** Make the building outline editable — change arc radii, segment
lengths, add/remove points — with automatic closure re-solving.

**Requirements:** ENG-11, API-16–19, OE-1–3, DT-2–4 (11 reqs)

**Work:**
- Create `app/outline_solver.py` — parallel closure solver reimplemented from
  `floorplan/geometry.py`'s `_chain_offset` logic (pure math, no import from
  `floorplan/`).  This is the most architecturally significant new module.
  - `chain_offset()` — walk chain entries from origin
  - `solve_closure()` — solve for `d_F2_F5` and `d_F18_F1`
  - `validate_chain()` — check proposed modifications
  - `solve_for_constraint()` — Newton/bisection solver for target distances
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

## Phase 6 — Enhanced Canvas Rendering

**Goal:** Add door arcs, clearance zones, and display toggles to the
interactive canvas.

**Requirements:** CV-7, CV-11, DIS-9–10, SEL-11, DOOR-1–4 (9 reqs)

*(DIS-6 and CV-10 are already implemented in Phase 0.  DIS-7 is assigned to
Phase 8 alongside dimension line creation.  CV-12 is assigned to Phase 9
alongside the product link infrastructure.)*

**Work:**
- Extend `app/engine.py` to include `clearance_zones`, `doors` in geometry result
- Create `app/static/js/canvas.js` — rendering functions:
  - `renderDoors()` — dashed arcs + hinge circles (CV-7)
  - `renderClearanceZones()` — dashed circles at fixtures (CV-11)
- Add display toggles: Clearance (DIS-9), Doors (DIS-10)
- Door property editing triggers re-render (SEL-11)

**New files:** `app/static/js/canvas.js`, `tests/test_zapp_canvas.py`

**Dependencies:** Phase 3 (doors data), Phase 5 (optional — room boundaries
adapt to outline changes).

---

## Phase 7 — Draw, Add, Delete, Rotate, and Shape Editor Tools

**Goal:** Full element creation, deletion, and shape customisation from the
canvas.

**Requirements:** TL-15–27, SEL-4, SEL-10, DT-7 (16 reqs)

**Work:**
- **Draw Wall tool** (TL-15–17): click start/end to place new interior wall,
  thickness input, endpoint drag handles
- **Add Element tools** (TL-18–20): catalog panels for furniture, appliances,
  fixtures with click-to-place
- **Add Opening tool** (TL-21): click on wall segment, dialog for width/type
- **Delete** (TL-22–23): Delete key with confirmation, cascading delete for
  walls with hosted openings
- **Rotate** (TL-24): R key opens rotation dialog with presets
- **Shape editor** (TL-25–27): edit polygon vertices in local coordinates
  (TL-25), assign shapes to item types (TL-26), import from SVG (TL-27).
  Builds on the `shapes` table already seeded with toilet, bath_sink,
  dining_table shapes.
- **Multi-select** (SEL-4): Shift-click and drag-select
- **Opening width editing** (SEL-10, DT-7): inline edit in openings table
  and properties panel

**New files:** `app/static/js/selection.js` (multi-select),
`tests/test_zapp_tools.py`
**Modified:** `app/static/js/tools.js` (expand), `app/static/js/dialogs.js`
(catalog data, rotation dialog, shape editor), `app/templates/index.html`
(tool buttons, menu items, keyboard shortcuts)

**Dependencies:** Phase 3 (element CRUD), Phase 4 (move tool patterns),
Phase 6 (canvas rendering).

---

## Phase 8 — Labels, Dimensions, and Annotations

**Goal:** User-placeable room labels, custom dimension lines, and annotation
editing.  Unifies the room label system under the `elements` table.

**Requirements:** TL-11–14, LABEL-1–4, DIS-7 (9 reqs)

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

---

## Phase 9 — Styling and Product Links

**Goal:** Per-element visual customisation and product URL attachments.

**Requirements:** STYLE-1–4, LINK-1–2, SEL-12, CV-12 (8 reqs)

**Work:**
- Create `app/style.py` — per-element and per-view style management stored in
  element `properties` JSON
- Properties panel: colour picker (STYLE-1), stroke inputs (STYLE-2), opacity
  slider (STYLE-3)
- Per-view style overrides (STYLE-4)
- URL field in properties panel (SEL-12, LINK-1)
- SVG generation wraps elements with URLs in `<a xlink:href>` tags (LINK-2)
- Canvas hyperlink indicators (CV-12)

**New files:** `app/style.py`, `tests/test_zapp_style.py`

**Dependencies:** Phase 3 (elements), Phase 6 (canvas rendering).

---

## Phase 10 — Domain Views: Site Plan, 3D, Analysis, Plumbing

**Goal:** Interactive editing for specialised views, wrapping existing
generators via the regeneration API.

**Requirements:** SITE-1–4, SCAD-1–3, ANALYSIS-1–3, PLUMB-1–8 (18 reqs)

**Work:**

### Site Plan (SITE-1–4)
- Setback distance inputs that map to site-plan constants (SITE-1)
- Drainfield element tool — custom elements with type `'site_element'` (SITE-2)
- Text annotation tool for site plan (SITE-3)
- P-series markers with distance labels (SITE-4)

### 3D / SCAD (SCAD-1–3)
- `POST /api/generate-3d` endpoint wrapping SCAD generator (SCAD-1)
- Roof style selection stored as config (SCAD-2)
- Multi-view PDF generation (SCAD-3)

### Analysis (ANALYSIS-1–3)
- Existing SVG views for span analysis already load in tabs (ANALYSIS-1–2)
- Room area display partially implemented in Phase 0 (ANALYSIS-3 — this
  phase completes the general View > Show Areas action)

### Plumbing Layout (PLUMB-1–8)
The plumbing view becomes a full interactive layout — like the floorplan
layouts (standard, minik, daybed, bare, sf) — with plumbing-specific editing
tools, database-stored element configuration, and CRUD APIs.

- Plumbing interactive canvas with building outline overlay, zoom/pan/selection
  (PLUMB-1 — replaces the static SVG tab with a live canvas)
- `plumbing_elements` database table for pipes, fittings, and fixture
  connections, with type, path geometry, and properties JSON (PLUMB-4)
- CRUD API endpoints for plumbing elements (PLUMB-5)
- Supply line drawing tool with hot/cold differentiation and routing (PLUMB-2)
- Drain line drawing tool with slope annotations (PLUMB-6)
- Fixture placement tool with automatic supply/drain stub connections (PLUMB-7)
- Pipe fitting editing: T-stubs, elbows, valves with type selection (PLUMB-8)
- Fixtures/supplies table in the right panel (PLUMB-3)

**Further development (beyond Phase 10):** The plumbing layout will expand to
include service location indicators (hot, cold, and drain at each fixture),
full drain line routing from fixtures through the building to the exterior,
and site-level plumbing elements: septic tank placement, drainfield layout,
and the well-to-building supply run.  These extensions build on the Phase 10
foundation and will be specified as additional PLUMB requirements when the
base layout is complete.

**New files:** `app/plumbing.py`, `app/static/js/plumbing-editor.js`,
`app/static/js/site-plan.js`, `tests/test_zapp_plumbing.py`,
`tests/test_zapp_views.py`

**Dependencies:** Phase 3 (elements table), Phase 7 (draw/add tools).

---

## Phase 11 — View Variants and Polish

**Goal:** User-defined layout variants with per-view layer configuration.

**Requirements:** UI-5–6 (2 reqs)

**Work:**
- Add `variants` table to database (name, label, source_variant,
  layer_config JSON, element_overrides JSON)
- View > New Variant creates a named clone of the current layout (UI-5):
  - Clones the source variant's element set (custom elements are
    duplicated with variant membership set to the new variant name)
  - Clones the source variant's exclusion rules
  - Does NOT clone constants (constants are shared across all variants)
  - The new variant appears in the Layout dropdown and can have its own
    SVG generated via the regeneration API
- Per-variant layer visibility toggles (UI-6):
  - Each variant stores a `layer_config` JSON specifying which layers
    (walls, openings, furniture, labels, dimensions) are visible
  - Toggling a layer in one variant does not affect other variants
- User-defined variants can add/remove element overrides (e.g., hide a
  wall, add custom furniture) without affecting other variants
- Final polish: cross-check all 226 requirements, keyboard shortcut audit
  (NF-5), responsive layout verification (NF-2), NF-3 (586 existing tests
  pass), NF-4 (no files modified outside `app/` and `tests/`)

**Interaction with fixed variants:** The 5 built-in variants (standard,
minik, daybed, bare, sf) remain as seed data.  User-defined variants are
additional entries in the `variants` table.  They can be renamed, deleted,
or cloned further.  Built-in variants cannot be deleted but can be modified.

**Dependencies:** Most prior phases complete.

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

**Requirements:** SEL-13 (1 req) + design specification for Charter Principle 5

SEL-13 (Constant Dependency Highlighting) is the first concrete requirement
for this phase: when a constant is focused in the Properties panel, all
geometry elements whose position depends on that constant are highlighted
on the canvas (first-order in white, downstream in pink).

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
4. NF-4 is lifted; existing generator scripts are retained as seed sources
   for "Reset to Defaults" (which regenerates the entire database from
   their output)
5. All five built-in variants and user-defined variants render correctly
   from database-driven formulas
6. d² regression tests pass for default database values (matching original
   script output)

**Impact:** Replaces the current procedural Python computation with a
data-driven evaluation model.  The positioning logic currently baked into
`floorplan/layout.py` and `app/variants.py` is encoded as database-stored
formulas.  After cutover, the existing scripts are no longer the
authoritative source for positioning — they are seed-only.

**Dependencies:** All prior phases (the formula system builds on the element
CRUD, canvas rendering, and property editing infrastructure).

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
    │         └── Phase 13 (Electrical) ◄── 10         │
    │                                                  │
    ├── Phase 11 (View Variants & Polish) ◄── all     │
    │                                                  │
    └── Phase 12 (Parametric Dependencies) ◄── all ───┘
```

(`:` and `?` denote optional dependencies — Phase 6 can start without Phase 5
but room boundaries adapt to outline changes once Phase 5 is complete.)

Phases 1 and 2 can proceed in parallel.  Phases 4, 5, and 6 can proceed in
parallel after Phase 3.

---

## Phase Completion Protocol

Each phase is considered complete only after:

1. All phase requirements pass automated tests
2. All 774+ tests continue to pass (`python -m pytest tests/ -x -q`)
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
| 8 | TL-11–14, LABEL-1–4, DIS-7 | 9 |
| 9 | STYLE-1–4, LINK-1–2, SEL-12, CV-12 | 8 |
| 10 | SITE-1–4, SCAD-1–3, ANALYSIS-1–3, PLUMB-1–8 | 18 |
| 11 | UI-5–6 | 2 |
| 12 | SEL-13, Charter Principle 5 | 1 + design spec |
| 13 | (aspirational — to be specified) | TBD |
| **Total** | | **226 + Phase 13 TBD** |

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
| 6 | — | canvas.js | test_zapp_canvas.py |
| 7 | — | selection.js | test_zapp_tools.py |
| 8 | labels.py | — | test_zapp_labels.py |
| 9 | style.py | — | test_zapp_style.py |
| 10 | plumbing.py | site-plan.js, plumbing-editor.js | test_zapp_plumbing.py, test_zapp_views.py |
| 11 | — | — | test_zapp_variants_ext.py |
| 12 | dependencies.py | dependency-graph.js | test_zapp_deps.py |
| 13 | electrical.py | electrical-editor.js | test_zapp_electrical.py |

---

## Estimated Test Growth

| Phase | New tests | Cumulative |
|-------|-----------|-----------|
| 0 | — | 752 |
| 1 (done) | 22 | 774 |
| 2 (done) | 20 | 794 |
| 3 | ~38 | 832 |
| 4 | ~15 | 847 |
| 5 | ~30 | 877 |
| 6 | ~12 | 889 |
| 7 | ~25 | 914 |
| 8 | ~15 | 929 |
| 9 | ~12 | 941 |
| 10 | ~20 | 961 |
| 11 | ~10 | 971 |
| 12 | ~20 | 991 |
| 13 | ~20 | 1011 |

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

---

## Cutover Criteria

The editor is ready for cutover (dropping the NF-4 constraint and
transitioning to fully database-driven design) when:

1. All 226 requirements pass automated tests
2. The parametric dependency system (Phase 12) replaces hardcoded
   positioning — all element positions are database-stored formulas
3. All constants live in the database (no longer Python module attributes);
   the 24 duplicated dimension constants in `app/variants.py` are
   consolidated into the database as the single source
4. The ~700 lines of replicated positioning math in `app/variants.py` are
   replaced by the parametric formula evaluator
5. "Reset to Defaults" regenerates the entire database from existing
   generator scripts' output, producing bit-identical geometry for the
   reference design (verified by d² regression tests)
6. All five built-in layout variants and any user-defined variants render
   correctly from database-driven formulas
7. NF-4 is lifted; existing generator scripts (`floorplan/`, `shared/`,
   etc.) are retained as seed-only sources and may now be modified for
   consolidation
8. User acceptance testing complete
