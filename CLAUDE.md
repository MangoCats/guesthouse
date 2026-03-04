# ADU Project

Curved-wall building outline geometry and floorplan SVG generation.

## Coordinate System
- FC (building center) = origin, Easting (E) / Northing (N), units in feet
- P3 = (-18.5, -13.5), defined relative to FC
- F-series outline: chain walk from F2 bearing north is the single source of truth
- Dimensions typically specified in inches, converted via `/ 12.0`

## Project Structure

```
shared/              — Common types, geometry, survey computation, SVG utilities
  types.py           — Point, BBox, LineSeg, ArcSeg, Segment
  geometry.py        — Pure geometry functions, path ops, polygon utilities, compute_inner_walls
  wall_shells.py     — Shell inset paths, U-turn arcs/polygons, wall section enumeration
  survey.py          — compute_traverse, compute_three_arc, compute_inset
  svg.py             — make_svg_transform, W/H page constants

floorplan/           — Building design: single source of truth for geometry and layout
  constants.py       — Named physical dimension constants (wall thicknesses, room sizes,
                       shell thickness, air gap, opening corner radius, etc.)
  geometry.py        — compute_outline_geometry → F-series points, segments, radii
  layout.py          — compute_interior_layout → InteriorLayout (rooms, appliances, furniture)
  openings.py        — compute_outer_openings, compute_rough_openings (single source for O1-O11, RO1-RO7)
  gen_floorplan.py   — Detailed floorplan SVG renderer. Outputs floorplan/floorplan.svg

walls/               — Outer wall construction detail drawing
  constants.py       — Re-exports SHELL_THICKNESS, AIR_GAP, OPENING_INSIDE_RADIUS from floorplan.constants
  gen_walls.py       — Double-shell wall detail SVG renderer. Outputs walls/walls.svg, walls/all_walls.svg

span/                — N-S span measurement analysis
  gen_span.py        — Interior span graph + plan view. Outputs span/span.svg
  gen_span_minmax.py — Span-vs-rotation analysis (5-175°). Outputs span/span_minmax.svg
  gen_span_min.py    — Minimum-span rotation finder. Outputs span/span_min.svg

survey/              — Survey scripts and data (not a Python package)
  gen_path_svg.py    — Outline + inset path SVG with labels. Outputs survey/path_area.svg
  compute_path.py    — Diagnostic/computation script
  rough_survey.txt   — Raw field measurements
  distances.md       — Theoretical distances
  adjust_pentagon.py — Least-squares survey adjustment

roof/                — Roof outline generation
  gen_roof.py        — Roof outline SVG with R-series corners. Outputs roof/roof.svg

site/                — Site plan generation
  gen_site_plan.py   — Building on survey parcel PDF. Outputs site/site_plan_df.pdf, site_plan_fs.pdf

scad/                — OpenSCAD 3D model and line drawing generation
  gen_flat_roof.py   — Flat roof SCAD model
  gen_2in12.py       — 2:12 slope roof SCAD model
  gen_views.py       — Multi-view rendered output
  gen_line_drawings.py — Line drawing SVGs from SCAD

plumbing/            — Plumbing plan generation
  gen_plumbing.py    — Water supply/drain plan SVG. Outputs plumbing/plumbing.svg
```

## Dependency Graph

```
survey/gen_path_svg.py ──→ floorplan/ ──→ shared/
                       └──→ shared/
floorplan/gen_floorplan.py ──→ floorplan/ ──→ shared/
                           └──→ shared/
walls/gen_walls.py ──→ walls/ ──→ floorplan/ ──→ shared/
                   └──→ shared/
span/gen_span*.py ──→ floorplan/ ──→ shared/
                  └──→ shared/
roof/gen_roof.py ──→ floorplan/ ──→ shared/
site/gen_site_plan.py ──→ floorplan/ ──→ shared/
scad/gen_*.py ──→ floorplan/ ──→ shared/
              └──→ shared/
plumbing/gen_plumbing.py ──→ floorplan/ ──→ shared/
survey/compute_path.py ──→ shared/
```

No circular dependencies. floorplan/ never imports from survey/, walls/, span/, roof/, site/, scad/, or plumbing/.

## Traversal Conventions
- **Survey traverse** (POB→P2→P3→P4→P5→POB): **CCW** as viewed from above
- **Construction outline** (F1→F2→...→F20→F1): **CW** as viewed from above (opposite of survey). Interior is on the **right** side of the traversal direction
- `left_norm(p1, p2)` returns the left perpendicular of the direction p1→p2. For CW traversal, left = exterior. Code uses `_wt = -wall_t` to negate and offset inward
- Individual arc `direction` ("CW"/"CCW" in `ArcSeg`) refers to each arc's own sweep direction, not the outline traversal. CW arcs (convex corners) get inner radius `R - wall_t`; CCW arcs (concave corners) get `R + wall_t`

## Key Patterns
- Outline points: F-series (`F1`, `F2`, `F5`..`F20`, `F11a`, `F11b` — no F3/F4), primary naming; U-series derived as aliases in survey/gen_path_svg.py
- Inner wall points: W-series (`W1`, `W2`, `W5`..`W20`, `W11a`, `W11b`), 8" inset from outline, matching F-series numbering
- Shell boundary points: S-series / G-series, matching F-series numbering. Computed via `shared/wall_shells.py:compute_inset_path` with custom inset + prefix rename
- Arc centers: C-series by lower point number (`C1`, `C5`, `C7`, `C8`, `C10`, `C11`, `C11a`, `C13`, `C15`, `C17`, `C19`); radii: R_a-series (`R_a1`, `R_a5`, ..., `R_a19`). `R_a1` = `CORNER_SW_R` (design constant); F20→F1 distance derived from closure
- Traverse arc centers: `TC1`, `TC2`, `TC3` (outer/inset path)
- `outline_segs`: list of 20 `LineSeg`/`ArcSeg` defining the closed outline path (CW traversal: F1→F2→F5→...→F11→F11a→F11b→F12→...→F20→F1)
- All radii in `OutlineGeometry.radii` dict; passed to `compute_inner_walls`
- Arc tangency: `|center1 - center2| = R1 + R2` (external)
- Physical constants defined once in `floorplan/constants.py` — no magic numbers in geometry/layout code

## F-Series Geometry Principle

The F-series outline is defined purely by:
1. **Start point**: F2 position (offset from FC) and bearing
2. **Segment definitions**: each segment is either a line (bearing, length) or an arc (radius, sweep, CW/CCW direction)
3. **Tangency**: automatic — each segment starts at the position and bearing where the previous segment ended
4. **Closure**: d_F2_F5 and d_F18_F1 are solved so the chain closes back to F2

**That's it.** No other runtime constraints exist in the geometry code. When a change request specifies additional constraints (e.g., "keep F9-F20 fixed", "maintain clear span"), those constraints are used to DERIVE new segment parameter values. Once the values are verified to satisfy the change request plus tangency and closure, the derived values are hardcoded as segment definitions and the change-request constraints are discarded. The solver in `geometry.py` must never accumulate constraints from past change requests.

## Workflow
- After each successful (as determined by passing of all tests) completed request, always: `git commit -a -m "<summary>"` then `python gen_all.py` to regenerate all SVGs.  summary shall be 25 words or less, and shall not include "Co-Authored-By: Claude".
- Outline geometry lives in `floorplan/geometry.py`; dimension constants in `floorplan/constants.py`
- Interior layout (rooms, furniture) lives in `floorplan/layout.py`
- Opening positions (O1-O11, RO1-RO7) live in `floorplan/openings.py` — single source of truth
- Wall construction constants (shell thickness, air gap, opening radius) defined in `floorplan/constants.py`, re-exported by `walls/constants.py`
- Shell geometry utilities (inset paths, U-turn polygons, wall sections) live in `shared/wall_shells.py`
- Pure geometry utilities (intersections, polygon ops) live in `shared/geometry.py`
- After geometry changes, update d² regression tests: `python tests/update_d2.py` (see HOWTO §14)

## App (Interactive Editor)

```
app/                 — Flask web editor: interactive canvas + constants editing
  server.py          — Flask routes, SSE broadcast, geometry cache
  engine.py          — Geometry computation orchestrator (patches floorplan.constants)
  database.py        — SQLite schema, seeding from floorplan/ sources, CRUD
  variants.py        — Variant furniture/appliance positioning (replicates gen_floorplan.py math)
  apputil.py         — Shared JSON serialisation helpers (point_to_list, bbox_from_poly, seg_to_dict)
  templates/         — index.html (single-page layout)
  static/js/app.js   — Client application (~1440 lines, no build step)
  static/css/app.css — Dark theme (Catppuccin palette)
  ARCHITECTURE.md    — Detailed module docs, computation flow, roadmap
  REQUIREMENTS.md    — 236 testable requirements (213 implemented, 23 planned)
  CHARTER.md         — Purpose, history, design principles
```

### NF-4: No Modification of Existing Packages
The app SHALL NOT modify files in `shared/`, `floorplan/`, `walls/`, `span/`,
`survey/`, `roof/`, `site/`, `scad/`, or `plumbing/`.  This constraint holds
until the editor achieves 100% functional completeness and is approved for
cutover.  Until then the existing generator scripts are the reference
implementations — the app must reproduce their output from database-driven
constants without changing them.

**Duplication between `app/` and existing packages is intentional.**
`app/variants.py` replicates positioning math from `floorplan/gen_floorplan.py`
and carries 24 hardcoded item-dimension constants that duplicate values in the
generator.  Do not "fix" this by modifying files outside `app/` — it will be
consolidated at cutover.  See `app/ARCHITECTURE.md` § NF-4 and `app/ROADMAP.md`
for the 12-phase development plan.

### Phase Completion Protocol
After each roadmap phase is complete (user acknowledges all goals met, no
outstanding issues), update `app/ROADMAP.md`, `app/ARCHITECTURE.md`, and
`app/REQUIREMENTS.md` to reflect the new state before proceeding to the next
phase.  See `app/ROADMAP.md` § Phase Completion Protocol for details.

### App Dependency Graph
```
app/server.py ──→ app/database.py ──→ floorplan/constants.py (seed source)
              └──→ app/engine.py  ──→ floorplan/ (geometry, layout, openings)
                                  └──→ shared/ (geometry, survey, types)
                                  └──→ app/variants.py ──→ shared/geometry.py
```
`floorplan/` and `shared/` never import from `app/`.

## HOWTO Reference
See `HOWTO.md` for step-by-step instructions on common tasks (adding dimension lines, walls, openings, appliances, identifying wall faces). Consult it before researching the codebase from scratch. If you complete a complex task not covered there, add a section to it.
