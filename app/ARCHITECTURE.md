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

Three SQLite tables:

| Table | Rows | Purpose |
|-------|------|---------|
| `constants` | 143 | Named numeric constants (value, unit, category) |
| `outline_chain` | 18 | Closed outline segments (line/arc definitions) |
| `views` | 11 | Registered SVG generators and output paths |

**Seeding** — On first run, `init_db()` creates the schema and populates
all three tables from source:

- **constants**: parsed from `floorplan/constants.py` source lines via
  regex, extracting name, expression, unit, and comment.  Category
  assigned by name prefix (`WALL_*` → wall, `IW_*` → interior_wall,
  `O*_WIDTH` → opening, etc.)
- **outline_chain**: copied from `OUTLINE_CHAIN` and `CHAIN_POINT_NAMES`
  in `floorplan/geometry.py`
- **views**: hardcoded list of 11 generator registrations

**CRUD functions**:

| Function | Purpose |
|----------|---------|
| `get_constants_dict()` | `{name: value}` dict for engine |
| `get_all_constants()` | Full rows for UI table |
| `update_constant()` | Single constant update |
| `update_constants_batch()` | Multi-constant transaction |
| `reset_constants()` | Clear and re-seed from source |
| `get_categories()` | Distinct category list |
| `get_outline_chain()` | 18 segment rows |
| `get_views()` | Enabled view definitions |

Connection management uses a context manager with WAL journaling and
foreign keys enabled.

### app/engine.py — Computation

**`compute_geometry(constants_dict, variant="standard")`** orchestrates
the full pipeline:

1. **Patch** `floorplan.constants` with database values, then reload
   `floorplan.geometry`, `.layout`, `.openings` so module-scope code
   re-executes with patched values
2. **Survey traverse** → F-series alignment
3. **Outline geometry** → 20 F-series points, 18 segments, arc radii
4. **Inner walls** → 18 W-series inset segments, closed polygon
5. **Interior layout** → 13 interior walls (IW1–IW12, IW2O, IW2S),
   rooms, appliance/furniture placements
6. **Openings** → 12 outer (O1–O11, O8A), 7 rough (RO1–RO7)
7. **Variant items** → 20–31 furniture/appliance items per variant
8. **Dimensions** → 18–22 dimension line endpoint pairs with distances

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
| bare | Room Dimensions | 0 (walls only) |
| sf | Square Footage | 0 (walls only) |

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
messages to all queues.  Events: `geometry_changed`, `svg_updated`,
`connected`.  Keepalive every 30 seconds.

**API endpoints** (14 total):

| Method | Path | Purpose |
|--------|------|---------|
| GET | `/` | Serve index.html |
| GET | `/api/constants` | List constants (optional `?category=`) |
| GET | `/api/constants/categories` | Category list |
| PUT | `/api/constants/<name>` | Update one constant |
| PUT | `/api/constants/batch` | Batch update |
| POST | `/api/constants/reset` | Reset to source defaults |
| GET | `/api/geometry` | Computed geometry (`?variant=`) |
| GET | `/api/variants` | Variant names and labels |
| GET | `/api/outline` | Outline chain (18 segments) |
| GET | `/api/views` | Registered views |
| GET | `/api/svg/<name>` | SVG content |
| GET | `/api/svg/<name>/file` | File download |
| POST | `/api/regenerate` | Run generator scripts |
| GET | `/api/events` | SSE stream |

### app/static/js/app.js — Client

Single-file client application (~1100 lines).  No build step or
framework.

**State** — `App.state` object holds: current geometry, constants array,
views list, active view/tool, zoom/pan, display toggles, variant
selection, sort/filter state.

**Initialisation** (on DOMContentLoaded):
1. `cacheElements()` — populate `App.els` with ~30 DOM references
2. `setupEventListeners()` — attach handlers to buttons, inputs, canvas
3. `connectSSE()` — open EventSource to `/api/events`
4. `loadViews()` → `renderViewTabs()` — build tab bar
5. `loadConstants()` → `renderConstantsTable()` — populate right panel
6. `loadGeometry()` → `renderCanvas()` — draw interactive view

**Canvas rendering** — `renderCanvas()` calls layer functions in order:
outline → inner walls → interior walls → openings → furniture →
points → dimensions.  Each function creates SVG elements (`<polygon>`,
`<circle>`, `<line>`, `<text>`) in the corresponding `<g>` layer.
Display toggles control which layers are populated.

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

The SVG canvas defines 9 layer groups: outline, inner, walls, openings,
furniture, points, labels, dims, measure.

### app/static/css/app.css — Styling

Dark theme (Catppuccin palette) via CSS custom properties.  18 colour
variables, 5 layout dimension variables.  Major style sections: menu bar,
tool palette, viewport, tabs, tables, right panel, SVG canvas classes
(outline/wall/opening/furniture/point/dimension/measure fills and
strokes), scrollbar, toast notifications.  Responsive breakpoint at
1000px collapses tool labels and narrows the right panel.

---

## Sources of Truth

| Data | Source of truth | Editable via app? |
|------|----------------|-------------------|
| Named constants | `constants` table (seeded from `floorplan/constants.py`) | Yes — inline editing |
| Outline chain | `outline_chain` table (seeded from `floorplan/geometry.py`) | Read-only |
| Interior layout | Computed by `floorplan/layout.py` from constants | Indirectly (edit constants) |
| Openings | Computed by `floorplan/openings.py` from constants | Indirectly (edit constants) |
| Variant items | Computed by `app/variants.py` from layout + constants | Indirectly (edit constants) |
| Dimension lines | Computed by `floorplan/gen_floorplan.py` from layout | Indirectly (edit constants) |
| SVG views | Generated by scripts listed in `views` table | Regenerate via menu |

The constants table is the single editable root.  Every other geometric
value is deterministically derived from it through the computation
pipeline.

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
    ├──► compute_outer_openings()       → 12 outer openings
    │    compute_rough_openings()       → 7 rough openings
    │
    ├──► compute_variant_items()        → 20–31 furniture/appliance items
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

Tests live in `tests/test_zapp_*.py` (prefixed to sort after existing
tests).  Shared fixtures in `tests/test_zapp_conftest.py` provide:

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
constraint applies until the editor has achieved 100% functional
completeness and has been approved for cutover to database-only data
sources.

**Consequence: intentional duplication.**  `app/variants.py` replicates
~700 lines of positioning math from `floorplan/gen_floorplan.py`'s
`_render_appliances()`, `_render_kitchen()`, and `_render_furniture()`
functions.  It also carries 24 hardcoded item-dimension constants
(hamper, microwave, dining table, etc.) that duplicate values scattered
through the generator.  This duplication is deliberate: the existing
scripts are the reference implementation, and the app must reproduce
their output without modifying them.  At cutover, these constants will
be consolidated into a single source (likely `floorplan/constants.py`)
and the shared positioning math extracted into a function callable by
both the SVG renderer and the app engine.

**Consequence: module reloading.**  The engine uses
`importlib.reload()` and `patch_constants()` to inject database values
into the floorplan modules at runtime (see Computation Flow above).
Five derived constants (`WALL_EXTRA`, `AIR_GAP`, `DOOR_FLAT_FACE`,
`F8F9_INNER_TURN_R`, `CORNER_SW_R`) are recomputed in `engine.py`
after patching because they depend on other constants and their
derivation formulas cannot be imported without modifying source.  This
will be unnecessary once the app owns the constants directly.

---

## Roadmap: Current State vs Charter Goals

The charter describes a full parametric editor; the current
implementation is a **parametric viewer with constant editing**.  This
section documents what is implemented, what is planned, and the
architectural additions each planned capability requires.

### Implemented

| Capability | Charter section | Status |
|------------|----------------|--------|
| Named constants: view, edit, persist, reset | Data Layer | Done |
| Outline chain: view (read-only) | Data Layer | Done |
| Geometry engine: full recomputation from constants | Data Layer | Done |
| Interactive canvas: outline, walls, openings, furniture, points, dims | Interactive Canvas | Done |
| Five layout variants switchable from UI | Variant Selector | Done |
| SVG view tabs: 11 generated views | SVG View Tabs | Done |
| Properties panel with related constants | Right Panel | Done |
| Constants table: sort, filter, inline edit, category colours | Right Panel | Done |
| Openings table: outer + rough with computed widths | Right Panel | Done |
| REST API: constants CRUD, geometry, variants, SVG, SSE | REST API | Done |
| Real-time update cycle: edit → recompute → broadcast → re-render | Real-Time | Done |
| Feet-inches display (NF-6) across all value surfaces | Units Display | Done |

### Not Yet Implemented

| Capability | Req IDs | Architectural impact |
|------------|---------|---------------------|
| **`elements` table** — interior walls, furniture, appliances, fixtures as database rows with JSON properties | DB-9 | New table, seed logic, CRUD API (API-20–26) |
| **`doors` table** — door configs per opening (hinge side, swing, type) | DB-10 | New table, seed logic, API (API-27–29), canvas rendering (CV-7) |
| **Undo/redo** — action history with before/after state | DB-11 | New table, API (API-30–31), UI (UNDO-1–4) |
| **Outline chain editing** — add/remove/reorder segments | API-16–19, ENG-11 | Mutation API, engine re-solve, d² test updates |
| **Element CRUD** — create, move, delete interior elements via UI | OE-1–3, TL-5–24 | Canvas tools, hit testing, drag/drop, constraint solver |
| **Parametric dependency tracking** (Charter Principle 5) | — | Major addition: dependency graph, formula storage per element, topological evaluation, dependency chain visualisation |
| **Room area computation and labels** | ENG-12, LABEL-1–4 | Polygon area calculation, label positioning |
| **Styling** — user-customisable colours, line weights | STYLE-1–4 | Preferences table or CSS variable API |
| **Site plan integration** — interactive site plan editing | SITE-1–4 | Extend canvas to site-plan mode |
| **3D/SCAD integration** — trigger SCAD generation from UI | SCAD-1–3 | Subprocess management, 3D viewer |
| **Analysis views** — span analysis from UI | ANALYSIS-1–3 | Parameterise span scripts |
| **Plumbing plan** — interactive plumbing editing | PLUMB-1–3 | Element layer for plumbing fixtures |

### Parametric Dependencies (Charter Principle 5)

The charter envisions that every element's position is stored as an
editable formula referencing other elements and constants, with
dependency chains displayed both as a formula table and as graphical
highlights.  The current architecture has no support for this:

- Positioning is computed procedurally in Python code, not stored as
  data.
- The three-table schema (constants, outline_chain, views) has no place
  to encode "IW3 spans from IW2.east_face to IW1.north_face".
- No dependency graph or topological sort exists.

Implementing this requires:

1. An `elements` table where each row stores positioning parameters as
   structured data (reference element, face, offset formula).
2. A dependency resolver that topologically sorts elements and evaluates
   formulas in order.
3. UI for displaying and editing dependency chains.
4. Canvas overlay for highlighting the dependency graph of the selected
   element.

This is the largest planned architectural change and should be designed
before incremental implementation begins.
