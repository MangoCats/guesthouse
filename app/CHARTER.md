# ADU Editor — Application Charter

## Purpose

This application provides an interactive web-based editor for the ADU
(Accessory Dwelling Unit) building design that has been developed across
869 commits of AI-assisted iterative design work in the project's other
folders (`shared/`, `floorplan/`, `walls/`, `span/`, `survey/`, `roof/`,
`site/`, `scad/`, `plumbing/`).

The editor's goal is to enable a user to **deterministically create,
edit, view, and persist** the same building designs that were previously
developed only through interactive chat with an AI agent.  At full
maturity (see ROADMAP.md Phase 12), **no content is hardcoded** — every
aspect of the design lives in the database: walls, openings, furniture
and appliance symbols and placements, dimension lines, product URLs,
style defaults, item dimensions, variant membership, and the parametric
dependencies among them.  All of this can be created, edited, stored,
and retrieved from the database without requiring further AI interaction.

The database is initially populated with the default content defined by
the existing generator scripts (`floorplan/`, `shared/`, etc.) and can
be restored to that state at any time via "Reset to Defaults."

## What the Project Built via Chat

The ~870 commits (as of initial charter, Mar 2026) span approximately
three weeks of intensive AI-assisted design (Feb 10 – Mar 1, 2026) and
cover every aspect of a curved-wall residential building, from outline
geometry through construction documents and site planning.  Many commits
touch multiple areas; the keyword-based counts below reflect this
overlap.

### Outline Geometry (~334 commits)

The building outline is a closed 20-point chain of line segments and
tangent arcs (F-series points F1–F20, with arcs at convex and concave
corners).  Work included defining the chain walk from F2, solving closure
and tangency constraints, adjusting arc radii and sweep angles, renumber-
ing points, and maintaining d² regression tests that lock the geometry
against unintended drift.

### Walls & Interior Layout (~267 commits)

Thirteen interior walls (IW1–IW9, IW11–IW12 plus IW2O, IW2S; no IW10)
divide the interior
into rooms: bedroom, office, bathroom, utility, kitchen, living, and
storage areas.  Work included wall-relative positioning of all interior
elements, parameterised wall thickness, double-shell construction
modelling with cavity and air gap, shell inset paths, U-turn arcs at
openings, and per-section wall enumeration for the construction detail
drawings.

### Furniture & Appliances (~215 commits)

Full furnishing across five layout variants (Standard, Small Kitchen,
Daybed, Room Dimensions, Square Footage).  Items include bed, dresser,
shelves, desk and chair, loveseat/sofa/daybed, dining table and chairs,
rocker, ottoman, stove, refrigerator, dishwasher, kitchen and bath
sinks, toilets, washer, dryer, water heater, counters, microwave,
coffee maker, cooktop, toaster, ice maker, and hamper.  All placed via
wall-relative vector math using the same `seg_vecs()` / `offset_pt()` /
`line_isect()` utilities from `shared/geometry.py`.

### Openings & Doors (~173 commits)

Twelve outer openings (O1–O11 plus O8A) and seven rough openings
(RO1–RO7) with door swings, hinge positions, jamb blocks, and casement
illustrations.  Opening widths are parameterised constants stored in
`floorplan/constants.py`.

### Dimension Lines & Annotations (~138 commits)

Twenty-three dimension lines (dim01–dim22, some with bare-only
variants) computed wall-relative via `compute_dimension_endpoints()`.
Room area labels, room name labels, font sizing, and label positioning
for construction document readability.

### SVG Generation & Rendering (~88 commits)

Multi-view SVG pipeline producing 11 registered views: Floorplan, Wall
Detail, All Walls, Span Analysis, Min Span, Span vs Rotation, Survey
Path, Roof, Plumbing, and two Site Plan variants.  Rendering utilities
in `shared/svg.py` with consistent transforms, page sizing, and a DRY
helper architecture.

### Span Analysis (~60 commits)

Interior N–S span measurement graph, span-vs-rotation analysis across
5–175 degrees, and minimum-span rotation finder for structural
evaluation.

### Site Plan & Survey (~51 commits)

Property boundary placement, survey traverse alignment (P-series to
F-series), setback and constraint annotations, drainfield
visualisation, parcel corner markers, and site plan PDF publishing.

### SCAD / 3D Model & Roof (~50 commits)

OpenSCAD wall extrusion with shell cavities, three-tier wall sections
(below doors, mid-height, full-height), flat and 2:12-slope roof
variants, and rendered PNG/PDF output.

### Plumbing (~27 commits)

Cold and hot water supply lines, fixture connections, T-stubs, waste
lines, well-to-building routing, and a water services table.

### Code Refactoring (~237 commits)

Module decomposition, type definitions (NamedTuple, BBox, LineSeg,
ArcSeg), DRY extraction of helper functions, rotation-invariant
placement conversion, wall-relative positioning refactoring, and package
reorganisation into the current `shared/` → `floorplan/` → `walls/` etc.
dependency hierarchy.

### Testing (~53 commits)

586 pre-existing tests covering geometry, layout, openings, walls, SVG
generation, site plan, and d² regression.  The app layer adds 515 more
tests covering database operations, engine computation, API endpoints,
variant items, dimensions, labels, styling, and undo/redo (1113 total).

## What the App Implements (as of Phase 0)

The editor makes the following subset of the above functionality
available through a browser-based GUI backed by an SQLite database.
See `ROADMAP.md` for the 13-phase plan to reach full charter capability.

### Data Layer
- **143 named constants** seeded from `floorplan/constants.py`, editable
  and persistable, with category assignment and parametric propagation
- **18 outline chain segments** defining the building boundary
- **Geometry engine** that recomputes all points, walls, openings,
  variant items, and dimension lines from the stored constants

### Interactive Canvas
- Building outline, 18 inner wall segments, 13 interior walls
- 19 openings (12 outer + 7 rough)
- Full variant furniture/appliance set (20–31 items per variant)
- 60+ labelled points (F-series, W-series, C-series)
- 18+ dimension lines with tick marks and feet-inches labels
- Zoom, pan, fit-to-window, and coordinate tracking
- Display toggles for points, labels, dimensions, grid, openings,
  furniture

### Variant Selector
- Five layout variants switchable via dropdown: Standard, Small Kitchen,
  Daybed, Room Dimensions, Square Footage
- Each variant shows the correct furniture and appliance set

### SVG View Tabs
- Embedded viewing of all 11 generated SVGs (floorplan, walls, span,
  survey, roof, plumbing, site plan)
- Regenerate individual or all views from the menu

### Right Panel
- **Properties tab**: selected element details in NF-6 compliant
  feet-inches format, with related constants
- **Constants tab**: sortable, filterable, inline-editable table of all
  143 constants with category colour coding
- **Outline tab**: 18-row segment chain table
- **Openings tab**: outer and rough openings with computed widths

### REST API
- Constants CRUD with batch update and reset
- Geometry computation with variant parameter
- Variants endpoint listing all five layouts
- SVG file serving and regeneration
- SSE event stream for real-time geometry change notifications
- Cache invalidation across all variants on constant changes

### Units Display
- All dimensions shown in feet and inches with inches to two decimal
  places and trailing zeroes removed (e.g. `5' 3.5"` not `5' 3.50"`)
- Consistent across coordinate display, property panel, constants table,
  measurement tool, dimension labels, and data tables

## Design Principles

1. **No modification of existing code (NF-4)** — the app imports from
   but never modifies `shared/`, `floorplan/`, `walls/`, or other
   existing packages.  This constraint holds during Phases 0–12;
   see Principle 4 for why and Principle 5 for when it ends.

2. **Parametric from database** — every geometric decision traces back
   to named values in the database; changing a value recomputes the
   entire design.  During Phases 0–11 the editable root is the
   `constants` table (seeded from `floorplan/constants.py`).  After
   Phase 12 the editable roots are per-element parametric formulas
   stored in the database, with constants as one possible input.

3. **Deterministic without AI** — once the database state is set, the
   geometry is fully determined by the computation pipeline; no AI
   agent is needed to reproduce or edit the design.

4. **Same math, different interface (transitional)** — during Phases
   0–11, while NF-4 is in effect, the app replicates the same
   positioning logic from the SVG generators (wall-relative vectors,
   line intersections, arc tangency) rather than introducing new
   geometric approaches.  The existing scripts remain the reference
   implementation; the app must reproduce their output for default
   values.  This principle is retired at cutover (see Principle 5).

5. **Database-driven parametric dependencies (target architecture)** —
   the charter's end state.  At full maturity (Phase 12 and beyond),
   **every** element — exterior walls, interior walls, openings,
   furniture, appliances, fixtures, dimension lines, labels — is a
   database-stored entity whose position and size are defined by
   editable formulas referencing other entities and/or constants.
   **All content** is database-stored: positions, dimensions, product
   URLs, visual styles, variant membership, door/clearance metadata,
   and inter-element relationships.  Nothing is hardcoded in Python or
   JavaScript — the code provides computation and rendering, while the
   database provides all design content.

   The parametric model supports:
   - **Fixed positions** — an element can be locked to absolute
     coordinates, independent of other elements
   - **Relative positions** — an element's position defined as an
     offset, bearing, or parametric expression referencing one or more
     other elements (e.g., "6 inches south of IW2 east face",
     "midpoint between W9 and W10")
   - **Mixed** — some parameters fixed, others relative (e.g., a wall
     with a fixed bearing but endpoints anchored to other walls)
   - **Lock/unlock** — any parametric dependency can be "locked" to
     freeze it at its current computed value, converting it to a fixed
     position so that upstream changes no longer propagate through it

   The dependency graph is a DAG persisted in the database.  A
   topological evaluator resolves formulas in dependency order to
   produce concrete positions.  When an element is selected, its
   dependency chain is displayed both as a formula table and
   graphically highlighted on the canvas.

   **Design-first workflow:** the system supports building a design
   from any starting point — place a bed, then define walls relative
   to it, then anchor openings to walls.  Or define the exterior shell
   first and position everything inward.  Dependencies can later be
   edited, reversed, or locked as the design matures.

   After Phase 12, the existing generator scripts (`floorplan/`,
   `shared/`, etc.) are no longer the authoritative source for any
   design content.  They remain only as seed sources: "Reset to
   Defaults" regenerates the entire database from their output —
   constants, element positions, product URLs, style defaults, item
   dimensions, variant configurations — to reproduce the reference
   design.  NF-4 is lifted, and code duplication between `app/` and
   the existing packages is consolidated.

## Transition from Principles 4 → 5

During development (Phases 0–11), Principle 4 governs: the app
replicates existing script logic and must produce bit-identical results
for default values.  Principle 5 describes the target architecture that
replaces Principle 4 at cutover.

The transition happened in two stages:

1. **Phases 0–11 (NF-4 in effect):** the database stored constants and
   element overrides; the computation pipeline patched existing modules
   and reloaded them.  Existing scripts were the reference implementation.
   d² regression tests enforced bit-identical output for default values.

2. **Phase 12 (cutover — complete):** the parametric dependency system
   replaced procedural computation.  Constants are database-stored values
   (not Python module attributes).  The FormulaEvaluator computes all
   element geometry from database-stored JSON formulas.  NF-4 is lifted.
   The existing scripts are retained as seed sources and for SVG output.