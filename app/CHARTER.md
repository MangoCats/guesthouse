# ADU Editor — Application Charter

## Purpose

This application provides an interactive web-based editor for the ADU
(Accessory Dwelling Unit) building design that has been developed across
869 commits of AI-assisted iterative design work in the project's other
folders (`shared/`, `floorplan/`, `walls/`, `span/`, `survey/`, `roof/`,
`site/`, `scad/`, `plumbing/`).

The editor's goal is to enable a user to **deterministically create,
edit, view, and persist** the same building designs that were previously
developed only through interactive chat with an AI agent.  All walls,
openings, furniture and appliance symbols and placements, dimension
lines, and tables of base data — as well as the parametric dependencies
among them — can be created, edited, stored, and retrieved from a
database without requiring further AI interaction.

## What the Project Built via Chat

The 869 commits span approximately three weeks of intensive AI-assisted
design (Feb 10 – Mar 1, 2026) and cover every aspect of a curved-wall
residential building, from outline geometry through construction
documents and site planning.  Many commits touch multiple areas; the
keyword-based counts below reflect this overlap.

### Outline Geometry (~334 commits)

The building outline is a closed 20-point chain of line segments and
tangent arcs (F-series points F1–F20, with arcs at convex and concave
corners).  Work included defining the chain walk from F2, solving closure
and tangency constraints, adjusting arc radii and sweep angles, renumber-
ing points, and maintaining d² regression tests that lock the geometry
against unintended drift.

### Walls & Interior Layout (~267 commits)

Thirteen interior walls (IW1–IW12 plus IW2O, IW2S) divide the interior
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
generation, site plan, and d² regression.  The app layer adds 117 more
tests covering database operations, engine computation, API endpoints,
variant items, and dimension data (703 total).

## What the App Implements

The editor makes the following subset of the above functionality
available through a browser-based GUI backed by an SQLite database:

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

1. **No modification of existing code** — the app imports from but never
   modifies `shared/`, `floorplan/`, `walls/`, or other existing
   packages (NF-4)
2. **Parametric from constants** — every geometric decision traces back
   to named constants in the database; changing a constant recomputes
   the entire design
3. **Deterministic without AI** — once the constants and chain are set,
   the geometry is fully determined by the computation pipeline; no AI
   agent is needed to reproduce or edit the design
4. **Same math, different interface** — the app replicates the same
   positioning logic from the SVG generators (wall-relative vectors,
   line intersections, arc tangency) rather than introducing new
   geometric approaches
