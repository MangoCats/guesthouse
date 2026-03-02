# ADU Editor — Development Roadmap

This roadmap traces the path from the current **parametric viewer with constant
editing** to the charter's goal: a **full parametric editor** where every
element can be created, moved, edited, deleted, and persisted from the browser
— all without AI assistance.

All work respects NF-4: no files outside `app/` and `tests/` are modified
until cutover (see ARCHITECTURE.md § NF-4).

---

## Current State (Phase 0 — Complete)

**80 of 189 requirements implemented.**  117 app tests, 586 pre-existing tests
(703 total).

| Capability | Status |
|------------|--------|
| Constants: view, edit, persist, reset | Done |
| Outline chain: view (read-only) | Done |
| Geometry engine: full recomputation | Done |
| Interactive canvas: outline, walls, openings, furniture, points, dims | Done |
| Five layout variants | Done |
| 11 SVG view tabs | Done |
| Properties panel with related constants | Done |
| Constants table: sort, filter, inline edit, category colours | Done |
| Openings table: outer + rough | Done |
| REST API: 14 endpoints, SSE | Done |
| Real-time update cycle | Done |
| Feet-inches display (NF-6) | Done |

**What's missing:** Element CRUD, doors, undo/redo, outline editing, move/draw
tools, room areas, labels, styling, site plan editing, 3D integration,
plumbing editing, parametric dependencies (Charter Principle 5).

---

## Phase 1 — Foundation Test Coverage

**Goal:** Achieve 100% automated test coverage for all 80 existing requirements
before any new development.  This is the safety net.

**Requirements:** DB-1–8, ENG-1–10, ENG-13, API-1–15, API-32–34, GEN-1–4,
UI-1–4, UI-7–8, CV-1–6, CV-9, CV-13–17, DIS-1–5, SEL-1–3, SEL-5–6, SEL-9,
TL-1–4, CT-1–10, DT-1, DT-5–6, DT-8, RT-1–4, APP-1–4, NF-1–6 (80 reqs)

**Work:**
- Create `tests/test_zapp_database.py` — schema, seeding, CRUD (DB-1–8)
- Create `tests/test_zapp_engine_full.py` — geometry keys, counts, propagation,
  derived constants (ENG-1–10, ENG-13)
- Create `tests/test_zapp_api_full.py` — all 14 endpoints (API-1–15, API-32–34)
- Extend existing test files as needed for UI/canvas requirements (manual +
  programmatic acceptance where possible)

**Deliverable:** Every existing requirement has at least one automated test.
~50 new tests.

**Dependencies:** None.

---

## Phase 2 — Undo/Redo System

**Goal:** Before adding mutation capabilities, establish the undo/redo
infrastructure so every future mutation is automatically reversible.

**Requirements:** DB-11, API-30–31, UNDO-1–4 (6 reqs)

**Work:**
- Add `undo_history` table to `app/database.py` (DB-11)
- Create `app/undo.py` — `UndoManager` class with command-pattern record/undo/redo
- Add `POST /api/undo` and `POST /api/redo` endpoints (API-30–31)
- Wrap existing constant mutation endpoints with undo recording
- Add `Ctrl+Z` / `Ctrl+Shift+Z` handlers to `app.js` (UNDO-1–2)
- 50-level depth (UNDO-3), cross-type undo (UNDO-4)

**New files:** `app/undo.py`, `tests/test_zapp_undo.py`
**Modified:** `app/database.py`, `app/server.py`, `app/static/js/app.js`

**Dependencies:** None (can run in parallel with Phase 1).

---

## Phase 3 — Elements & Doors Database

**Goal:** Extend the schema to represent interior elements and doors as
first-class database objects, with full CRUD APIs.

**Requirements:** DB-9–10, API-20–29, DT-9–11, SEL-7–8, RT-5, UI-7 (21 reqs)

**Work:**
- Add `elements` table (DB-9) and `doors` table (DB-10) to database schema
- Seed elements from `InteriorLayout` fields (13 walls, 6 furniture/appliances,
  fixtures) and doors from RO1–RO7 defaults
- Create `app/elements.py` — CRUD, cascading delete, constant-to-element mapping
- Create `app/doors.py` — door arc geometry computation, CRUD
- Add 10 new API endpoints: elements CRUD (API-20–26), doors (API-27–29)
- Add "Elements" tab to right panel (UI-7 completion) with interior walls table
  (DT-9–10) and furniture/appliance table (DT-11)
- Enhanced property display for furniture (SEL-8) and doors (SEL-7)
- Broadcast `element_changed` SSE event on mutations (RT-5)

**New files:** `app/elements.py`, `app/doors.py`, `tests/test_zapp_elements.py`,
`tests/test_zapp_doors.py`
**Modified:** `app/database.py`, `app/server.py`, `app/engine.py`,
`app/templates/index.html`, `app/static/js/app.js`

**Key design:** Base elements (the 13 walls, furniture from `floorplan/layout.py`)
are positioned by constants.  Moving them = changing the controlling constants.
Custom (user-added) elements store absolute positions in the `properties` JSON and
are overlaid on the engine output.  A hand-maintained mapping table in
`app/elements.py` maps element names to their controlling constants (e.g.,
`"IW1"` → `IW1_OFFSET_FROM_W9`).

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

**Requirements:** ENG-11, API-16–19, OE-1–3, DT-2–4 (10 reqs)

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

**New files:** `app/outline_solver.py`, `app/static/js/outline-editor.js`,
`tests/test_zapp_outline.py`

**Challenge:** The module-scope closure solver in `floorplan/geometry.py` runs
on import.  After outline modification, we must inject our solver's computed
`d_F2_F5` and `d_F18_F1` as patched constants before the engine reloads modules.
Extensive cross-validation tests between the two solvers are critical.

**Dependencies:** Phase 2 (undo), Phase 3 (element CRUD for cascading effects).

---

## Phase 6 — Enhanced Canvas Rendering

**Goal:** Add door arcs, room labels, area labels, clearance zones, and
hyperlink indicators to the interactive canvas.

**Requirements:** CV-7–8, CV-10–12, DIS-6–10, ENG-12, SEL-11 (13 reqs)

**Work:**
- Create `app/room_areas.py` — room boundary polygon construction + shoelace
  area computation + centroid calculation (ENG-12)
- Extend `app/engine.py` to include `room_areas`, `room_centroids`,
  `clearance_zones`, `doors` in geometry result
- Create `app/static/js/canvas.js` — rendering functions:
  - `renderDoors()` — dashed arcs + hinge circles (CV-7)
  - `renderRoomLabels()` — room names at centroids (CV-8)
  - `renderAreaLabels()` — "125 sf" at centroids (CV-10)
  - `renderClearanceZones()` — dashed circles at fixtures (CV-11)
  - `renderHyperlinkIndicators()` — link icons (CV-12)
- Add SVG layers and display toggles: Room Labels (DIS-6), Dimensions (DIS-7),
  Areas (DIS-8), Clearance (DIS-9), Doors (DIS-10)
- Door property editing triggers re-render (SEL-11)

**New files:** `app/room_areas.py`, `app/static/js/canvas.js`,
`tests/test_zapp_areas.py`

**Dependencies:** Phase 3 (doors data), Phase 5 (optional — room boundaries
adapt to outline changes).

---

## Phase 7 — Draw, Add, Delete, and Rotate Tools

**Goal:** Full element creation and deletion from the canvas.

**Requirements:** TL-15–24, SEL-4, SEL-10, DT-7 (14 reqs)

**Work:**
- **Draw Wall tool** (TL-15–17): click start/end to place new interior wall,
  thickness input, endpoint drag handles
- **Add Element tools** (TL-18–20): catalog panels for furniture, appliances,
  fixtures with click-to-place
- **Add Opening tool** (TL-21): click on wall segment, dialog for width/type
- **Delete** (TL-22–23): Delete key with confirmation, cascading delete for
  walls with hosted openings
- **Rotate** (TL-24): R key opens rotation dialog with presets
- **Multi-select** (SEL-4): Shift-click and drag-select
- **Opening width editing** (SEL-10, DT-7): inline edit in openings table
  and properties panel

**New files:** `app/static/js/selection.js` (multi-select),
`tests/test_zapp_tools.py`
**Modified:** `app/static/js/tools.js` (expand), `app/static/js/dialogs.js`
(catalog data, rotation dialog), `app/templates/index.html` (tool buttons,
menu items, keyboard shortcuts)

**Dependencies:** Phase 3 (element CRUD), Phase 4 (move tool patterns),
Phase 6 (canvas rendering).

---

## Phase 8 — Labels, Dimensions, and Annotations

**Goal:** User-placeable room labels, custom dimension lines, and annotation
editing.

**Requirements:** TL-11–14, LABEL-1–4, DIS-7 (8 reqs)

**Work:**
- Create `app/labels.py` — label and dimension line management using the
  `elements` table (type `'label'` or `'dimension'`)
- **Dimension Line tool** (TL-11): click two points to create dimension with
  extension lines, arrowheads, and ft-in label
- Dimension line repositioning (TL-12), label rotation (TL-13), deletion (TL-14)
- **Room Label tool** (LABEL-1): click to place, type label text
- Label dragging (LABEL-2), inline text editing (LABEL-3), font size (LABEL-4)
- Persistent dimensions toggle (DIS-7)

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

**Requirements:** SITE-1–4, SCAD-1–3, ANALYSIS-1–3, PLUMB-1–3 (13 reqs)

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
- Room area display implemented in Phase 6 (ANALYSIS-3)

### Plumbing (PLUMB-1–3)
- Plumbing SVG view already loads in tabs (PLUMB-1)
- Supply line drawing tool with hot/cold differentiation (PLUMB-2)
- Fixtures/supplies table (PLUMB-3)

**New files:** `app/static/js/site-plan.js`, `tests/test_zapp_views.py`

**Dependencies:** Phase 3 (elements table), Phase 7 (draw/add tools).

---

## Phase 11 — View Variants and Polish

**Goal:** User-defined layout variants with per-view layer configuration.

**Requirements:** UI-5–6 (2 reqs)

**Work:**
- Add `variants` table to database (name, label, source, layer_config JSON)
- View > New Variant creates a named clone of the current layout (UI-5)
- Per-variant layer visibility toggles (UI-6)
- Final polish: cross-check all 189 requirements, keyboard shortcut audit (NF-5),
  responsive layout verification (NF-2), NF-3 (586 existing tests pass),
  NF-4 (no files modified outside `app/` and `tests/`)

**Dependencies:** Most prior phases complete.

---

## Phase 12 — Parametric Dependencies (Charter Principle 5)

**Goal:** The charter's most ambitious feature — every element's position
stored as an editable formula referencing other elements, with dependency
chains displayed both as a formula table and graphically highlighted in the
drawing.

This is the largest architectural addition and should be designed as a
standalone specification before implementation begins.

**Architectural requirements:**

1. **Formula storage** — Each element's positioning parameters stored as
   expressions referencing other elements and constants (e.g., "IW3 spans
   from IW2.east_face to IW1.north_face")
2. **Dependency graph** — Directed acyclic graph of element references,
   persisted in the database
3. **Topological evaluator** — Sorts elements by dependency order, evaluates
   formulas to produce concrete positions
4. **Dependency chain UI** — When an element is selected, its dependency chain
   is displayed as a formula table and graphically highlighted on the canvas

**Impact:** Replaces the current procedural Python computation with a
data-driven evaluation model.  This is a fundamental shift — the positioning
logic currently baked into `floorplan/layout.py` and `app/variants.py` would
be encoded as database-stored formulas.

**Dependencies:** All prior phases (the formula system builds on the element
CRUD, canvas rendering, and property editing infrastructure).

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
    │         │                                        │
    │         └── Phase 6 (Canvas Rendering)            │
    │                   │                              │
    │                   │                              │
    ├── Phase 7 (Draw/Add/Delete) ◄── 3, 4, 6         │
    │                                                  │
    ├── Phase 8 (Labels & Dimensions) ◄── 3, 7        │
    │                                                  │
    ├── Phase 9 (Styling & Links) ◄── 3, 6            │
    │                                                  │
    ├── Phase 10 (Domain Views) ◄── 3, 7              │
    │                                                  │
    ├── Phase 11 (View Variants & Polish) ◄── all     │
    │                                                  │
    └── Phase 12 (Parametric Dependencies) ◄── all ───┘
```

Phases 1 and 2 can proceed in parallel.  Phases 4, 5, and 6 can proceed in
parallel after Phase 3.

---

## Phase Completion Protocol

Each phase is considered complete only after:

1. All phase requirements pass automated tests
2. All 703+ tests continue to pass (`python -m pytest tests/ -x -q`)
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
4. **Commit the documentation update** as a separate commit after the phase
   implementation commit(s)

This ensures the documentation always reflects the true state of the project
and that fresh context windows have accurate information.

---

## Requirement-to-Phase Mapping

| Phase | Requirement IDs | Count |
|-------|----------------|-------|
| 0 (done) | DB-1–8, ENG-1–10, ENG-13, API-1–15, API-32–34, GEN-1–4, UI-1–4, UI-7–8, CV-1–6, CV-9, CV-13–17, DIS-1–5, SEL-1–3, SEL-5–6, SEL-9, TL-1–4, CT-1–10, DT-1, DT-5–6, DT-8, RT-1–4, APP-1–4, NF-1–6 | 80 |
| 1 | (test coverage for above) | 0 new |
| 2 | DB-11, API-30–31, UNDO-1–4 | 6 |
| 3 | DB-9–10, API-20–29, DT-9–11, SEL-7–8, RT-5, UI-7 | 21 |
| 4 | TL-5–10, API-23 | 7 |
| 5 | ENG-11, API-16–19, OE-1–3, DT-2–4 | 10 |
| 6 | CV-7–8, CV-10–12, DIS-6–10, ENG-12, SEL-11 | 13 |
| 7 | TL-15–24, SEL-4, SEL-10, DT-7 | 14 |
| 8 | TL-11–14, LABEL-1–4, DIS-7 | 8 |
| 9 | STYLE-1–4, LINK-1–2, SEL-12, CV-12 | 8 |
| 10 | SITE-1–4, SCAD-1–3, ANALYSIS-1–3, PLUMB-1–3 | 13 |
| 11 | UI-5–6 | 2 |
| 12 | Charter Principle 5 | (design spec) |
| **Total** | | **189 + Principle 5** |

---

## New Files by Phase

| Phase | New Python files | New JS files | New test files |
|-------|-----------------|-------------|---------------|
| 1 | — | — | test_zapp_database.py, test_zapp_engine_full.py, test_zapp_api_full.py |
| 2 | undo.py | — | test_zapp_undo.py |
| 3 | elements.py, doors.py | — | test_zapp_elements.py, test_zapp_doors.py |
| 4 | — | tools.js, dialogs.js | test_zapp_move.py |
| 5 | outline_solver.py | outline-editor.js | test_zapp_outline.py |
| 6 | room_areas.py | canvas.js | test_zapp_areas.py |
| 7 | — | selection.js | test_zapp_tools.py |
| 8 | labels.py | — | test_zapp_labels.py |
| 9 | style.py | — | test_zapp_style.py |
| 10 | — | site-plan.js | test_zapp_views.py |
| 11 | — | — | test_zapp_variants_ext.py |
| 12 | dependencies.py | dependency-graph.js | test_zapp_deps.py |

---

## Estimated Test Growth

| Phase | New tests | Cumulative |
|-------|-----------|-----------|
| 0 (current) | — | 703 |
| 1 | ~50 | 753 |
| 2 | ~20 | 773 |
| 3 | ~30 | 803 |
| 4 | ~15 | 818 |
| 5 | ~30 | 848 |
| 6 | ~15 | 863 |
| 7 | ~25 | 888 |
| 8 | ~15 | 903 |
| 9 | ~12 | 915 |
| 10 | ~12 | 927 |
| 11 | ~10 | 937 |
| 12 | ~20 | 957 |

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

4. **Frontend Module Splitting (Phase 4+)** — `app.js` is 1100 lines.  Adding
   tools, canvas rendering, and selection will push past 3000 without splitting.
   New JS files (`tools.js`, `canvas.js`, `selection.js`, etc.) share the global
   `App` namespace via additional `<script>` tags.

5. **Custom Element Overlay (Phase 3)** — The engine computes base geometry from
   constants.  Custom elements exist only in the DB and are overlaid on the
   result.  This dual-source model (constants-driven base + DB-driven overlay)
   must avoid duplication or conflict.

6. **Parametric Dependencies (Phase 12)** — Fundamentally replaces procedural
   computation with data-driven evaluation.  Requires its own design spec before
   implementation.  This is the largest planned change by far.

---

## Cutover Criteria

The editor is ready for cutover (dropping the NF-4 constraint and consolidating
code across `app/` and existing packages) when:

1. All 189 requirements pass automated tests
2. All five layout variants render identically to reference SVGs
3. Parametric dependency system (Phase 12) replaces hardcoded positioning
4. All 24 duplicated dimension constants in `app/variants.py` have been moved
   to `floorplan/constants.py` (or database) as single source
5. The ~700 lines of positioning math in `app/variants.py` have been extracted
   into a shared function callable by both SVG renderer and app engine
6. All generator scripts produce identical output from database-driven constants
7. User acceptance testing complete
