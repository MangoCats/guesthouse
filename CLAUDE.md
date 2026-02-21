# Hut2 Project

Curved-wall building outline geometry and floorplan SVG generation.

## Coordinate System
- F1 = origin, E / N axes, units in feet; rotated ~6.34° from survey E/N so bearing F20→F1 = 270°
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
  openings.py        — compute_outer_openings, compute_rough_openings (single source for O1-O11, RO1-RO5)
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
survey/compute_path.py ──→ shared/
```

No circular dependencies. floorplan/ never imports from survey/, walls/, or span/.

## Traversal Conventions
- **Survey traverse** (POB→P2→P3→P4→P5→POB): **CCW** as viewed from above
- **Construction outline** (F1→F2→...→F20→F1): **CW** as viewed from above (opposite of survey). Interior is on the **right** side of the traversal direction
- `left_norm(p1, p2)` returns the left perpendicular of the direction p1→p2. For CW traversal, left = exterior. Code uses `_wt = -wall_t` to negate and offset inward
- Individual arc `direction` ("CW"/"CCW" in `ArcSeg`) refers to each arc's own sweep direction, not the outline traversal. CW arcs (convex corners) get inner radius `R - wall_t`; CCW arcs (concave corners) get `R + wall_t`

## Key Patterns
- Outline points: F-series (`F1`..`F20`, `F11a`, `F11b`), primary naming; U-series derived as aliases in survey/gen_path_svg.py
- Inner wall points: W-series (`W1`..`W20`, `W11a`, `W11b`), 8" inset from outline, matching F-series numbering
- Shell boundary points: S-series (`S1`..`S20`, `S11a`, `S11b`) = 2" inset (inner face of outer shell); G-series (`G1`..`G20`, `G11a`, `G11b`) = 6" inset (outer face of inner shell). Computed via `shared/wall_shells.py:compute_inset_path` with custom inset + prefix rename
- Arc centers: C-series by lower point number (`C1`, `C3`, `C5`, `C7`, `C8`, `C10`, `C11`, `C11a`, `C13`, `C15`, `C17`, `C19`); radii: R_a-series (`R_a1`, `R_a3`, ..., `R_a19`)
- Traverse arc centers: `TC1`, `TC2`, `TC3` (outer/inset path)
- `outline_segs`: list of 22 `LineSeg`/`ArcSeg` defining the closed outline path (CW traversal: F1→F2→...→F11→F11a→F11b→F12→...→F20→F1)
- All radii in `OutlineGeometry.radii` dict; passed to `compute_inner_walls`
- Arc tangency: `|center1 - center2| = R1 + R2` (external)
- Physical constants defined once in `floorplan/constants.py` — no magic numbers in geometry/layout code

## Workflow
- After each successful (as determined by passing of all tests) completed request, always: `git commit -a -m "<summary>"` then `python gen_all.py` to regenerate all SVGs.  summary shall be 25 words or less, and shall not include "Co-Authored-By: Claude".
- Outline geometry lives in `floorplan/geometry.py`; dimension constants in `floorplan/constants.py`
- Interior layout (rooms, furniture) lives in `floorplan/layout.py`
- Opening positions (O1-O11, RO1-RO5) live in `floorplan/openings.py` — single source of truth
- Wall construction constants (shell thickness, air gap, opening radius) defined in `floorplan/constants.py`, re-exported by `walls/constants.py`
- Shell geometry utilities (inset paths, U-turn polygons, wall sections) live in `shared/wall_shells.py`
- Pure geometry utilities (intersections, polygon ops) live in `shared/geometry.py`

## HOWTO Reference
See `HOWTO.md` for step-by-step instructions on common tasks (adding dimension lines, walls, openings, appliances, identifying wall faces). Consult it before researching the codebase from scratch. If you complete a complex task not covered there, add a section to it.
