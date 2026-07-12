# Generator Inventory

Authoritative table of generator scripts, their outputs, and their data
source.  Maintained alongside `gen_all.py` — when adding a generator, update
this file (see HOWTO §15).

**Data sources:**

- **DB (in-process)** — run by `gen_all.py` via
  `app/engine.py:_run_generator_inprocess()` with `GeneratorData` built from
  the live database (`app/adu.db`).  Reflects the loaded configuration.
- **Seed (subprocess)** — run as `python <script>`; derives geometry from the
  seed `floorplan/` modules, survey traversal, or 3D model parameters rather
  than DB formulas.

## Run by `python gen_all.py`

### In-process (DB-driven)

| Script | Outputs |
|--------|---------|
| `floorplan/gen_floorplan.py` (rendered by `app/db_render.py`) | `floorplan/floorplan.svg`, `floorplan_minik.svg`, `floorplan_db.svg`, `floorplan_bare.svg`, `floorplan_sf.svg`, `floorplan_plumbing.svg` |
| `walls/gen_walls.py` | `walls/walls.svg`, `walls/all_walls.svg` |
| `span/gen_span.py` | `span/span.svg` |
| `span/gen_span_minmax.py` | `span/span_minmax.svg` |
| `span/gen_span_min.py` | `span/span_min.svg` |
| `roof/gen_roof.py` | `roof/roof.svg` |
| `site/gen_site_plan.py` | `site/site_plan_df.pdf`, `site/site_plan_fs.pdf` (+ PNGs) |
| `scad/gen_flat_roof.py` | `scad/flat_roof.scad` (+ rendered PNGs via gen_views) |
| `scad/gen_2in12.py` | `scad/2in12.scad` (+ rendered PNGs via gen_views) |
| `survey/gen_path_svg.py` | `survey/path_area.svg` |

### Subprocess (seed-driven)

| Script | Outputs |
|--------|---------|
| `survey/gen_path_svg_wo.py` | `survey/path_area_wo.svg` |
| `survey/gen_path_svg_ks.py` | `survey/path_area_ks.svg` |
| `survey/gen_points_pdf.py` | `survey/points.pdf` |
| `scad/gen_views.py` | `scad/*_patio.png`, `*_corner.png`, `*_bumpout.png` (per roof style) |
| `scad/gen_line_drawings.py` | `scad/*_ld.png` line drawings |
| `gen_3views.py` | `3views.pdf` (+ PNG) |
| `plumbing/gen_plumbing.py` | `plumbing/plumbing.svg` (standalone; the DB path renders `floorplan/floorplan_plumbing.svg` instead) |
| `scad/gen_gltf.py` | `opts/*.glb` glTF models |
| `blend/gen_blend.py` | `blend/<config>*.blend` Blender models (see `blend/README.md`; requires Blender) |

After the scripts run, `gen_all.py` renders the site-plan PDFs and
`3views.pdf` to 1920px PNGs (requires PyMuPDF).

## Registered app views (`views` table, 14 rows)

Floorplan, Wall Detail, All Walls, Span Analysis, Min Span, Span vs Rotation
(disabled by default), Survey Path, Roof, Plumbing (floorplan variant), Site
Plan DF/FS, 3D Flat Roof, 3D 2:12 Roof, 3-View Layout.  The app's
Regenerate menu drives the same in-process/subprocess dispatch as
`gen_all.py`.

## Standalone-only generators (not in gen_all.py)

| Script | Outputs | Notes |
|--------|---------|-------|
| `scad/gen_1in12.py`, `scad/gen_split1.py`, `scad/gen_split2.py` | `scad/1in12.scad`, `split1.scad`, `split2.scad` | Alternate roof styles; also generated on demand by gen_gltf/gen_blend |
| `elevations/gen_elevations.py` | elevation drawings | |
| `sketch/gen_sketch_nnw.py` | perspective sketch | |
| `opt_rot/gen_opt_rot.py` | `opt_rot/path_area_*.svg` | C17 sweep analysis (see `opt_rot/ABOUT_OPT_ROT.md`) |
| `gen_options.py` | `options/index.html`, `opts/index.html` + GLBs | Option comparison pages |
| `gen_overlay.py` | overlay comparison drawing | |
| `ani/make_video.py` | floorplan history MP4 | See `ani/README.md` |
