# Hut2 Project

Curved-wall building outline geometry and floorplan SVG generation.

## Coordinate System
- P3 = origin, Easting (E) / Northing (N), units in feet
- Dimensions typically specified in inches, converted via `/ 12.0`

## Project Structure

```
shared/              — Common types, geometry, survey computation, SVG utilities
  types.py           — Point, LineSeg, ArcSeg, Segment
  geometry.py        — Pure geometry functions, path ops, polygon utilities, compute_inner_walls
  survey.py          — compute_traverse, compute_three_arc, compute_inset
  svg.py             — make_svg_transform, W/H page constants

floorplan/           — Building design: single source of truth for geometry and layout
  constants.py       — Named physical dimension constants (wall thicknesses, room sizes, etc.)
  geometry.py        — compute_outline_geometry → F-series points, segments, radii
  layout.py          — compute_interior_layout (rooms, appliances, furniture)
  gen_floorplan.py   — Detailed floorplan SVG renderer. Outputs floorplan/floorplan.svg

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
survey/compute_path.py ──→ shared/
```

No circular dependencies. floorplan/ never imports from survey/.

## Key Patterns
- Outline points: F-series (`F0`..`F21`), primary naming; U-series derived as aliases in survey/gen_path_svg.py
- Inner wall points: W-series (`W0`..`W21`), 8" inset from outline, matching F-series numbering
- Arc centers: C-series by lower point number (`C0`, `C2`, `C3`, `C5`, `C7`, `C8`, `C10`, `C11`, `C13`, `C15`, `C17`, `C19`, `C20`); radii: R_a-series (`R_a0`, `R_a2`, ..., `R_a20`)
- Traverse arc centers: `TC1`, `TC2`, `TC3` (outer/inset path)
- `outline_segs`: list of `LineSeg`/`ArcSeg` defining the closed outline path (CCW traversal: F0→F1→...→F21→F0)
- All radii in `OutlineGeometry.radii` dict; passed to `compute_inner_walls`
- Arc tangency: `|center1 - center2| = R1 + R2` (external); inner wall radii: CW arcs get `-wall_t`, CCW get `+wall_t` (CCW traversal convention)
- `_inner_point` in shared/geometry.py uses `_wt = -wall_t` to place inner walls on the interior (left) side of CCW traversal
- Physical constants defined once in `floorplan/constants.py` — no magic numbers in geometry/layout code

## Workflow
- After geometry changes, verify by running: `python survey/gen_path_svg.py` then `python floorplan/gen_floorplan.py`
- Outline geometry lives in `floorplan/geometry.py`; dimension constants in `floorplan/constants.py`
- Interior layout (rooms, furniture) lives in `floorplan/layout.py`
- Pure geometry utilities (intersections, polygon ops) live in `shared/geometry.py`
