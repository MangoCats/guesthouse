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
