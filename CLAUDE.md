# Hut2 Project

Curved-wall building outline geometry and floorplan SVG generation.

## Coordinate System
- P3 = origin, Easting (E) / Northing (N), units in feet
- Dimensions typically specified in inches, converted via `/ 12.0`

## Files
- `gen_path_svg.py` — Outline geometry (`pts` dict, `outline_segs`), radii, SVG rendering. Outputs `path_area.svg`
- `survey.py` — Shared: `LineSeg`/`ArcSeg`, `compute_inner_walls`, `compute_interior_layout`, polygon utilities
- `floorplan/gen_floorplan.py` — Detailed floorplan SVG. Imports from both above. Outputs `floorplan/floorplan.svg`
- `compute_path.py` — Diagnostics/computation script

## Key Patterns
- Outline points: U-series (`U0`..`U21`), sequential; display labels F-series (`F0`..`F21`)
- Inner wall points: W-series (`W0`..`W21`), 8" inset from outline, matching U-series numbering
- Arc centers: C-series (`Cf`, `Cf2`, `Cf4`, `Cc1`, `Cc2`, `Ct1`..`Ct4`, `Ct6a`, `Cw1`..`Cw3`) — unchanged
- `outline_segs`: list of `LineSeg`/`ArcSeg` defining the closed outline path (CCW traversal: U0→U1→...→U21→U0)
- Radii exported from `gen_path_svg` and passed via `_radii` dict to `compute_inner_walls`
- Arc tangency: `|center1 - center2| = R1 + R2` (external); inner wall radii: CW arcs get `-wall_t`, CCW get `+wall_t` (CCW traversal convention)
- `_inner_point` in survey.py uses `_wt = -wall_t` to place inner walls on the interior (left) side of CCW traversal

## Workflow
- Outline geometry changes typically touch all 3 files: gen_path_svg (geometry + SVG config), survey (inner_segs), gen_floorplan (imports + _radii)
- After changes, verify by running: `python gen_path_svg.py` then `python floorplan/gen_floorplan.py`
- When adding new radii: export from gen_path_svg, add to _radii dicts in both gen_path_svg and gen_floorplan, add to survey docstring
