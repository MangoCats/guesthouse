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
  roof_outline.py    — Roof outline computation (R-series corners from roof_corners config)
  split_roof.py      — Split-roof (split1/split2) geometry helpers
  survey.py          — compute_traverse, compute_three_arc, compute_inset
  svg.py             — make_svg_transform, W/H page constants, git_describe

floorplan/           — Reference geometry and seed source (NOT the live source of truth — see DB Authority)
  constants.py       — Seed values for DB constants table (wall thicknesses, room sizes,
                       shell thickness, air gap, opening corner radius, etc.)
  geometry.py        — compute_outline_geometry → F-series points, segments, radii (seed for outline_chain)
  layout.py          — compute_interior_layout → InteriorLayout (seed/reference only; app uses DB formulas)
  openings.py        — compute_outer_openings, compute_rough_openings (seed/reference only; app uses DB formulas)
  roof.py            — R-series roof geometry (used by app/gen_provider.py)
  gen_floorplan.py   — Floorplan SVG rendering engine. The DB path (app/db_render.py)
                       renders all variants: floorplan.svg, _minik, _db, _bare, _sf, _plumbing

walls/               — Outer wall construction detail drawing
  gen_walls.py       — Double-shell wall detail SVG renderer. Outputs walls/walls.svg, walls/all_walls.svg

span/                — N-S span measurement analysis
  _common.py         — Shared span computation helpers
  gen_span.py        — Interior span graph + plan view. Outputs span/span.svg
  gen_span_minmax.py — Span-vs-rotation analysis (5-175°). Outputs span/span_minmax.svg
  gen_span_min.py    — Minimum-span rotation finder. Outputs span/span_min.svg

survey/              — Survey scripts and data (not a Python package)
  gen_path_svg.py    — Outline + inset path SVG with labels. Outputs survey/path_area.svg
  gen_path_svg_wo.py — Variant: outline-only path drawing. Outputs survey/path_area_wo.svg
  gen_path_svg_ks.py — Variant: known-scale path drawing. Outputs survey/path_area_ks.svg
  gen_points_pdf.py  — Survey points PDF. Outputs survey/points.pdf
  compute_path.py    — Diagnostic/computation script
  rough_survey.txt   — Raw field measurements
  distances.md, updated_distances.md — Theoretical distances
  adjust_pentagon.py — Least-squares survey adjustment

roof/                — Roof outline generation
  gen_roof.py        — Roof outline SVG with R-series corners. Outputs roof/roof.svg

site/                — Site plan generation
  gen_site_plan.py   — Building on survey parcel PDF. Outputs site/site_plan_df.pdf, site_plan_fs.pdf

scad/                — OpenSCAD 3D model and line drawing generation
  _common.py         — Shared SCAD generation helpers
  gen_flat_roof.py   — Flat roof SCAD model
  gen_2in12.py       — 2:12 slope roof SCAD model
  gen_1in12.py       — 1:12 slope roof SCAD model
  gen_split1.py      — Split roof (1:12) SCAD model
  gen_split2.py      — Split roof (2:12) SCAD model
  gen_gltf.py        — GLB (glTF) 3D model export; SCAD parser + mesh builders
                       reused by blend/gen_blend.py
  gen_views.py       — Multi-view rendered output
  gen_line_drawings.py — Line drawing SVGs from SCAD

blend/               — Blender (.blend) 3D model generation from the live DB.
                       See blend/README.md for the full pipeline (roof variants,
                       manifest protocol, ground plane, fences, camera animation).

plumbing/            — Plumbing plan generation
  gen_plumbing.py    — Standalone water supply/drain plan SVG. Outputs plumbing/plumbing.svg.
                       The app also renders plumbing as a floorplan variant
                       (floorplan/floorplan_plumbing.svg) — see PLUMBING_PLAN.md
```

### Additional directories and root scripts

```
app/                 — Interactive editor (see § App below and app/ARCHITECTURE.md)
tests/               — Pytest suite (geometry/layout regression + app tests)
elevations/          — gen_elevations.py: elevation drawings
sketch/              — gen_sketch_nnw.py: perspective sketch from reference photo
opt_rot/             — C17 sweep-angle analysis (see opt_rot/ABOUT_OPT_ROT.md)
ani/                 — Floorplan-history animation tooling (see ani/README.md)
options/, opts/      — Published option-comparison pages and GLB model exports (gen_options.py)
main/                — Reference site photos
X/                   — ARCHIVAL: frozen snapshot of generated outputs; not live code
gen_all.py           — Regenerate all outputs from the DB (see GENERATORS.md)
gen_3views.py        — 3-view PDF layout. Outputs 3views.pdf
gen_options.py       — Builds options/ and opts/ comparison pages
gen_overlay.py       — Overlay comparison drawing
run_app.py           — Launch the interactive editor (Flask, port 5000)
recover_markx*.py    — One-off DB recovery scripts (historical)
publish.sh           — Publishing helper
```

See `GENERATORS.md` for the authoritative generator → output → data-source table.

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
blend/gen_blend.py ──→ scad/gen_gltf.py ──→ floorplan/ ──→ shared/
gen_all.py ──→ app/engine.py (DB-driven in-process generation)
```

No circular dependencies. floorplan/ never imports from survey/, walls/, span/, roof/, site/, scad/, or plumbing/.

## Traversal Conventions
- **Survey traverse** (POB→P2→P3→P4→P5→POB): **CCW** as viewed from above
- **Construction outline** (F1→F2→F5→...→F18→F1): **CW** as viewed from above (opposite of survey). Interior is on the **right** side of the traversal direction
- `left_norm(p1, p2)` returns the left perpendicular of the direction p1→p2. For CW traversal, left = exterior. Code uses `_wt = -wall_t` to negate and offset inward
- Individual arc `direction` ("CW"/"CCW" in `ArcSeg`) refers to each arc's own sweep direction, not the outline traversal. CW arcs (convex corners) get inner radius `R - wall_t`; CCW arcs (concave corners) get `R + wall_t`

## Key Patterns
- Outline points: F-series (`F1`, `F2`, `F5`..`F18`, `F11a`, `F11b` — no F3/F4, no F19/F20), primary naming; U-series derived as aliases in survey/gen_path_svg.py. The live DB chain may add or rename points (e.g. `W01`-style names in some saved configs) — the `outline_chain` table is authoritative
- Inner wall points: W-series (`W1`, `W2`, `W5`..`W18`, `W11a`, `W11b`), inset from outline by `WALL_OUTER` (8" seed default; live value in DB constants), matching F-series numbering
- Shell boundary points: S-series / G-series, matching F-series numbering. Computed via `shared/wall_shells.py:compute_inset_path` with custom inset + prefix rename
- Arc centers: C-series by lower point number (`C1`, `C5`, `C7`, `C8`, `C10`, `C11a`, `C11`, `C13`, `C15`, `C17`); radii: R_a-series (`R_a1`, `R_a5`, `R_a7`, `R_a8`, `R_a10`, `R_a11`, `R_a13`, `R_a15`, `R_a17`). `R_a1` = `CORNER_SW_R` (design constant); `d_F2_F5` and `d_F18_F1` derived from closure
- Traverse arc centers: `TC1`, `TC2`, `TC3` (outer/inset path)
- `outline_segs`: `LineSeg`/`ArcSeg` list defining the closed outline path (18 segments in the seed chain; the live DB chain may differ — CW traversal: F2→F5→F6→...→F11→F11a→F11b→F12→...→F18→F1→F2)
- All radii in `OutlineGeometry.radii` dict; passed to `compute_inner_walls`
- Arc tangency: `|center1 - center2| = R1 + R2` (external)
- Physical constants defined in `floorplan/constants.py` (seed/reference) and in the DB `constants` table (live/authoritative). **When docs or code comments cite a physical number, note whether it is the seed default or the live DB value** — e.g. `WALL_OUTER` is 8" in the seed but may be larger in the loaded config

## F-Series Geometry Principle

The F-series outline is defined purely by:
1. **Start point**: F2 position (offset from FC) and bearing
2. **Segment definitions**: each segment is either a line (bearing, length) or an arc (radius, sweep, CW/CCW direction)
3. **Tangency**: automatic — each segment starts at the position and bearing where the previous segment ended
4. **Closure**: d_F2_F5 and d_F18_F1 are solved so the chain closes back to F2

**That's it.** No other runtime constraints exist in the geometry code. When a change request specifies additional constraints (e.g., "keep F9-F18 fixed", "maintain clear span"), those constraints are used to DERIVE new segment parameter values. Once the values are verified to satisfy the change request plus tangency and closure, the derived values are hardcoded as segment definitions and the change-request constraints are discarded. The solver in `geometry.py` must never accumulate constraints from past change requests.

## DB Authority

**The database is the sole source of truth for all building geometry.**
`floorplan/constants.py`, `floorplan/geometry.py`, `floorplan/layout.py`, and
`floorplan/openings.py` are **seed/reference sources only** — used exclusively
to populate a fresh DB via File → Reset Database.  They are never consulted
during normal app operation or SVG generation.

The only legitimate use of these hardcoded modules outside of seeding is:
- `floorplan/geometry.py`: outline geometry for standalone reference scripts
  (`gen_all.py` now uses the DB for all scripts with in-process handlers;
  only survey/SCAD scripts without handlers still use the procedural path)
- `floorplan/gen_floorplan.py` and other generators: still used as the rendering
  engine, but always fed DB-sourced data via `GeneratorData` / `compute_geometry()`

## Workflow
- After each successful (as determined by passing of all tests) completed request, always: `git commit -a -m "<summary>"` then `python gen_all.py` to regenerate all SVGs.  summary shall be 25 words or less, and shall not include "Co-Authored-By: Claude".
- `gen_all.py` reads the DB (`app/adu.db`) and uses in-process generation for all supported scripts; survey and SCAD scripts without in-process handlers fall back to subprocess (uses `floorplan/constants.py` as reference — acceptable since those scripts derive geometry from survey traversal / 3D model parameters, not DB formulas).
- Outline geometry lives in `floorplan/geometry.py` (seed/reference); live outline in the DB `outline_chain` table
- Interior layout (rooms, furniture) lives in `floorplan/layout.py` (seed/reference only; app uses DB formulas)
- Opening positions (O1-O11 + O8A, RO1-RO7) live in `floorplan/openings.py` (seed/reference only; app uses DB formulas)
- Wall construction constants (shell thickness, air gap, opening radius) defined in `floorplan/constants.py` (seed) and the DB `constants` table (live)
- Shell geometry utilities (inset paths, U-turn polygons, wall sections) live in `shared/wall_shells.py`
- Pure geometry utilities (intersections, polygon ops) live in `shared/geometry.py`
- After geometry changes, update d² regression tests: `python tests/update_d2.py` (see HOWTO §14)
- **Doc update protocol (root docs):** when a change renames geometry points, adds/removes a generator script, or adds/changes API endpoints or DB tables, update the affected sections of `CLAUDE.md`, `HOWTO.md`, `COORDINATE_SYSTEMS.md`, `GENERATORS.md`, and `app/DATABASE.md` in the same commit. Avoid hardcoding volatile counts (endpoints, tests, line counts) in prose — reference the source instead.

## App (Interactive Editor)

```
app/                 — Flask web editor: interactive canvas + full parametric editing
  server.py          — Flask routes (REST + SSE), geometry cache, config save/load
  engine.py          — Geometry computation orchestrator (formula-only since Phase 12g;
                       does NOT patch floorplan.constants)
  evaluator.py       — FormulaEvaluator: topo-sorted evaluation of DB-stored JSON
                       formulas (see app/FORMULAS.md)
  database.py        — SQLite schema + CRUD (see app/DATABASE.md)
  database_seeds.py  — Seed data (variants, views, shapes, exclusions)
  db_render.py       — DB-path floorplan SVG renderer (all variants incl. plumbing)
  gen_provider.py    — GeneratorData provider + inner wall override engine
  outline_solver.py  — DB-chain closure solver (no floorplan/ imports)
  elements.py, doors.py, labels.py, style.py, plumbing.py, undo.py,
  variants.py, parse_input.py, apputil.py — supporting modules (see ARCHITECTURE.md)
  configs/           — Named saved configurations (File → Save As), e.g. Mark10.db
  templates/         — index.html (single-page layout)
  static/js/         — Client (app.js + tools, catalog, dialogs, shape-editor,
                       outline-editor; no build step)
  static/css/app.css — Dark theme (Catppuccin palette)
  ARCHITECTURE.md    — Detailed module docs, computation flow, roadmap
  REQUIREMENTS.md    — Testable requirements (see its summary table for counts)
  DATABASE.md        — Database schema reference
  FORMULAS.md        — Formula language reference
  CHARTER.md         — Purpose, history, design principles (historical snapshot)
```

### NF-4: Lifted (Phase 12g Cutover Complete)
NF-4 has been lifted.  The app's FormulaEvaluator is the sole source of element
geometry — `compute_geometry()` no longer patches `floorplan.constants` or
reloads modules.  The `floorplan/` modules are seed/reference sources only.
Modifications to `shared/` and `floorplan/` are permitted when needed.

The generator scripts (`gen_floorplan.py`, `gen_walls.py`, etc.) remain the
rendering engines but are always invoked with DB-sourced data via
`GeneratorData` / `compute_geometry()`.  They never drive geometry themselves.

### Phase Completion Protocol
After each roadmap phase is complete (user acknowledges all goals met, no
outstanding issues), update `app/ROADMAP.md`, `app/ARCHITECTURE.md`, and
`app/REQUIREMENTS.md` to reflect the new state before proceeding to the next
phase.  See `app/ROADMAP.md` § Phase Completion Protocol for details.

### App Dependency Graph
```
app/server.py ──→ app/database.py ──→ floorplan/constants.py (seed source)
              └──→ app/engine.py  ──→ floorplan/ (geometry only)
                                  └──→ shared/ (geometry, survey, types)
                                  └──→ app/variants.py
```
`floorplan/` and `shared/` never import from `app/`.

## HOWTO Reference
See `HOWTO.md` for step-by-step instructions on common tasks (adding dimension lines, walls, openings, appliances, identifying wall faces). Consult it before researching the codebase from scratch. If you complete a complex task not covered there, add a section to it.
