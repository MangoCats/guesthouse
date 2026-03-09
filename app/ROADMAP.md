# ADU Editor — Development Roadmap

This roadmap traces the path from the current **parametric viewer with constant
editing** to the charter's goal: a **full parametric editor** where every
element can be created, moved, edited, deleted, and persisted from the browser
— all without AI assistance.

NF-4 has been lifted (Phase 12g cutover complete).  The FormulaEvaluator
is the sole source of element geometry.  Phase 12h eliminated all procedural
element baselines — `compute_geometry()` is now formula-only.

Phases 15½ and 16–19 are complete: inner wall segment overrides are
DB-driven and editable from the UI; all SVG/PDF generators accept
`GeneratorData` as an optional parameter, and the app's regeneration
endpoints pass a DB-built `GeneratorData` in-process.  Standalone
execution (via `gen_all.py`) constructs `GeneratorData` from the hardcoded
procedural modules.

---

## Current State (Phase 19 + 15½ complete)

**285 of 285 requirements implemented.**  ~1100 app tests, ~590 pre-existing
tests (~1695 total).  All implemented requirements have automated test coverage.

| Capability | Status |
|------------|--------|
| Constants: view, edit, persist, reset | Done |
| Outline chain: view, edit, add/remove, closure solver | Done |
| Geometry engine: full recomputation (formula-only) | Done |
| Interactive canvas: outline, walls, openings, furniture, points, dims | Done |
| Five layout variants + user-defined variants | Done |
| 15 view tabs (SVG + PDF + PNG + plumbing canvas) with zoom/pan | Done |
| Properties panel with related constants | Done |
| Constants table: sort, filter, inline edit, category colours | Done |
| REST API: 70 endpoints, SSE | Done |
| Undo/redo: 50-level, Ctrl+Z / Ctrl+Shift+Z, cross-type | Done |
| Elements table: 13 IW seeds, ~59 total elements, CRUD | Done |
| Doors table: 9 seeds (7 RO + O3, O6), CRUD, validation | Done |
| Move tool: drag walls/furniture, ghost preview, multi-select, grid snap | Done |
| Outline editing: F-point drag, arc handle, closure solver | Done |
| Door arcs, clearance zones, display toggles | Done |
| Draw Wall, Add Element, Delete, Rotate, Shape Editor tools | Done |
| Dimension lines, room labels, anchored dimensions | Done |
| Per-element styling (fill, stroke, opacity) and product URLs | Done |
| Site plan, 3D/SCAD, span analysis, plumbing canvas | Done |
| Formula evaluator: topo sort, cycle detection, ~50 formula types | Done |
| Lock/unlock, dependency highlighting, dependency graph API | Done |
| DB-seeded element metadata, product URLs, variant formulas | Done |
| F2 origin, survey data, span analysis — all DB-driven | Done |
| Full project export/import | Done |
| Generator data provider (`GeneratorData`) | Done |
| All generators accept `GeneratorData` (Phases 16–18) | Done |
| DB-driven regeneration: in-process dispatch (Phase 19) | Done |
| Inner wall segment overrides: DB, engine, API, editor UI (Phase 15½) | Done |

**What's missing:** Electrical layout (Phase 13, aspirational).

---

## Completed Phases (1–19)

### Phase 1 — Foundation Test Coverage
22 new tests (774 total).  100% coverage audit for 121 requirements.

### Phase 2 — Undo/Redo System
DB-11, API-30–31, UNDO-1–4 (7 reqs).  `app/undo.py`, 50-level depth,
command-pattern, 20 new tests (794 total).

### Phase 3 — Elements & Doors Database
DB-9–10, API-20–29, DT-9–11, SEL-7–8, RT-5 (17 reqs).  `elements` and
`doors` tables, 13 IW seeds, 9 door seeds, CRUD APIs, undo support.
45 new tests (839 total).

### Phase 4 — Move Tool
TL-5–10, API-23 (7 reqs).  Drag walls/furniture, ghost preview, shift-
constrain, grid snap, multi-select group move, offset dialog.
20 new tests (862 total).

### Phase 5 — Outline Chain Editing
ENG-11, API-16–19, OE-1–3, DT-2–4 (11 reqs).  `app/outline_solver.py`
parallel closure solver, chain mutation APIs, F-point drag, arc handle.
37 new tests (899 total).

### Phase 6 — Enhanced Canvas Rendering
CV-7, CV-11, DIS-9–10, SEL-11, DOOR-1–4 (9 reqs).  Door arcs (structural
+ appliance), clearance zones, display toggles.  28 new tests (927 total).

### Phase 7 — Draw, Add, Delete, Rotate, Shape Editor
TL-15–27, SEL-4, SEL-10, DT-7 (16 reqs).  Draw Wall tool, catalog dialog,
delete with cascade, rubber-band select, rotate, shape editor with SVG import.
27 new tests (954 total).

### Phase 8 — Labels, Dimensions, Annotations
TL-11–14, LABEL-1–4, DIS-7, DIM-1–5, VAR-1–3, SEL-13–14 (19 reqs).
Dimension lines, room labels, anchored dimensions, multi-layout visibility.
90 new tests (1044 total).

### Phase 9 — Styling and Product Links
STYLE-1–4, LINK-1–2, SEL-12, CV-12 (8 reqs).  `app/style.py` with 3-layer
merge, colour/stroke/opacity editing, product URLs, SVG link wrapping.
69 new tests.

### Phase 10 — Domain Views
SITE-1–4, SCAD-1–3, ANALYSIS-1–3, PLUMB-1–11 (21 reqs across 4 sub-phases).
Site plan PDF, 3D/SCAD generation, span analysis, plumbing interactive canvas
with drawing tools.  ~95 new tests.

### Phase 11 — View Variants and Polish
UI-5–6, SEL-8a (3 reqs).  DB-driven variants table, user-defined variant
creation/deletion with cloning, per-variant layer config persistence.
44 new tests.

### Phase 12 — Parametric Dependencies and Cutover
8 sub-phases (12a–12h).  FormulaEvaluator with topo sort, cycle detection,
~50 formula types.  All ~59 elements migrated to DB-stored formulas.  NF-4
lifted.  `variants.py` reduced from 984→149 lines.  257 new tests (1516 total).

### Phase 14 — DB-Driven Outline & Export/Import
4 sub-phases (14-A–D).  F2 origin in DB, survey data in DB, span analysis
without module patching, full project export/import with validation.

### Phase 15 — Generator Data Provider
`app/gen_provider.py` exposing `GeneratorData` with outline, interior, shell,
roof, and survey geometry.  Golden-gate identity tests vs hardcoded modules.

### Phase 16 — Floorplan & Roof Generator Migration
`gen_floorplan.py` and `gen_roof.py` accept `GeneratorData`.
`build_floorplan_data(gd=None)` bridges to `FloorplanData`.

### Phase 17 — Wall Detail Generator Migration
`gen_walls.py` accepts `GeneratorData` via `build_wall_data(gd=None)`.

### Phase 18 — Remaining Generator Migration
`span/_common.py`, `site/gen_site_plan.py`, `scad/gen_flat_roof.py`,
`scad/gen_2in12.py` all accept `GeneratorData`.

### Phase 19 — DB-Driven Regeneration
In-process generator dispatch replaces subprocess execution.
`build_generator_data_from_db()` builds `GeneratorData` from DB state;
`_run_generator_inprocess()` dispatches to 11 generator render functions;
`generate_svg_db()` wraps with subprocess fallback.  All regeneration
endpoints (`/api/regenerate`, `/api/generate-site-plan`, `/api/generate-3d`)
pass DB-built `GeneratorData`.  17 new tests.

### Phase 15½ — Inner Wall Segment Overrides
4 sub-phases (15½-A–D).  `inner_wall_overrides` DB table with
`walk_override_chain()` parametric engine.  Generalises the hardcoded
W8-W9 corner treatment into a data-driven system supporting single-segment
and multi-segment span overrides.  Editor UI with live canvas preview,
compute-default, and click-to-define mode.  Endpoint validation (position
and bearing tolerance), overlap detection.  ~28 new tests (1695 total).

---

## Phase 13 — Electrical Layout (Aspirational)

**Goal:** An interactive electrical layout, following the plumbing layout
pattern (Phase 10): interactive canvas with database-stored elements, CRUD
APIs, and domain-specific drawing tools.

**Status:** Aspirational — no reference generator exists; built from scratch.

**Scope (to be specified):**
- `electrical_elements` database table for circuits, wiring runs, devices
- CRUD API endpoints, undo/redo
- Interactive canvas with building outline overlay
- Circuit drawing tool, device placement tools (outlets, switches, lights)
- Panel schedule with load calculations
- NEC compliance annotations

**Dependencies:** Phase 10 (plumbing pattern), Phase 18 (all generators
DB-backed).

---

## Phase Dependency Graph

```
Phases 1–19 + 15½ (all complete)
    │
    └── Phase 13 (Electrical) ◄── 10, 18
```

---

## Phase Completion Protocol

Each phase is considered complete only after:

1. All phase requirements pass automated tests
2. All 1695+ tests continue to pass (`python -m pytest tests/ -x -q`)
3. All SVGs regenerate successfully (`python gen_all.py`)
4. User acknowledgement that all phase goals are met

**Before proceeding to the next phase**, update this ROADMAP.md,
`app/ARCHITECTURE.md`, and `app/REQUIREMENTS.md` to reflect the new state.

**Database policy:** Schema changes do not require migration logic during
development.  Databases are wiped and recreated from defaults when the
schema changes.
