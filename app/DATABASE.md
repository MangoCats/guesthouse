# ADU Editor — Database Schema Reference

The SQLite database `app/adu.db` is the **sole source of truth** for all
building geometry and design content (CLAUDE.md § DB Authority).  This file
documents the 19 tables and the JSON property shapes stored inside them.

Regenerate the authoritative column list at any time:

```
python -c "import sqlite3; [print(r[0], '\n') for r in sqlite3.connect('app/adu.db').execute(\"SELECT sql FROM sqlite_master WHERE type='table'\")]"
```

Schema DDL lives in `app/database.py:init_db()`; seed data in
`app/database.py` + `app/database_seeds.py`.  **Database policy:** schema
changes do not require migrations — the DB is wiped and reseeded (File →
Reset Database) when the schema changes.  Named snapshots are stored as
whole-file copies in `app/configs/<name>.db`.

## Geometry roots

| Table | Key columns | Purpose |
|-------|-------------|---------|
| `constants` | `name` PK, `value` REAL, `expr`, `unit` ('ft'), `category`, `description` | Named numeric constants — the editable roots (wall thicknesses, offsets, opening widths). Seeded by regex-parsing `floorplan/constants.py`. |
| `outline_chain` | `seq` PK, `seg_type` ('L'/'CW'/'CCW'), `distance`, `radius`, `sweep_name`, `sweep` (radians), `center_name`, `n_pts`, `end_name`, `flex` ('distance'/'radius'/'sweep'/NULL), `bearing_flex` | The exterior outline chain, walked and closure-solved by `app/outline_solver.py`. `flex` marks the solver-adjusted parameters. |
| `survey_legs` | `seq` PK, `bearing_deg/min/sec`, `distance_ft`, `distance_inch`, `label` | The 5-leg survey traverse (POB→P2→P3→P4→P5→POB). |
| `survey_config` | `key` PK, `value` (JSON) | FC_IN_P3, COORD_ROTATION, per-leg corrections. |

## Elements and formulas

| Table | Key columns | Purpose |
|-------|-------------|---------|
| `elements` | `id` PK, `type`, `name` UNIQUE, `properties` (JSON), `variant` (NULL = all) | Every design element: walls, openings, furniture, appliances, fixtures, clearances, labels, dimensions. |
| `element_formulas` | `element_name`, `param_name`, `formula_json`, `locked`, `locked_value`, `variant`; UNIQUE(element, param, variant) | Parametric position/poly formulas — see `app/FORMULAS.md`. |
| `formula_deps` | `element_name`, `param_name`, `dep_type` ('element'/'point'/'constant'), `dep_name` | Dependency edges driving topological evaluation and the `/api/deps/graph` DAG. |
| `doors` | `id` PK, `opening_name` UNIQUE, `width` (inches), `hinge_side`, `swing_direction`, `door_type` ('single'/'double') | Door configuration per opening. |
| `inner_wall_overrides` | (`seg_index`, `sub_seq`) PK, `span_end`, `seg_type` ('L'/'CW'/'CCW'), `bearing` (deg), `distance`, `radius`, `sweep` (deg), `n_pts` | Parametric line/arc chains spliced into the W-series inner wall polygon (Phase 15½). |
| `shapes` | `name` PK, `poly_json`, `scale`, `origin` ('center'/'corner'), `width_key`, `depth_key`, `description` | Custom item outline polygons (toilet, sinks, dining table) referenced by `shape_transform` formulas. |
| `catalog_items` | `key` PK, `item_type`, `label`, `width`, `depth`, `radius`, `shape`, `door` (JSON), `clearance` (JSON), `product_url` (JSON), `variants` (JSON array), `stacked`, `clip_to_inner` | Placement-catalog templates for the Add Element tool. |

### `elements.type` values

`wall`, `opening`, `furniture`, `appliance`, `fixture`, `clearance`,
`label`, `dimension`.

### `elements.properties` JSON keys (by element kind)

All optional unless noted; unknown keys are preserved on update.

- **Common/style:** `fill_color`, `stroke_color`, `stroke_width`,
  `stroke_style` ('solid'/'dashed'/'dotted'), `opacity` (0–100),
  `view_overrides` (`{variant: {partial style}}`) — see `app/style.py`.
- **Items (furniture/appliance/fixture):** `label`, `item_type`, `shape`
  ('rect'/'circle'/shape name), `door` (`{hinge_idx, target_idx}` or
  variant-keyed dict of the same), `clearance` (face vertex indices +
  distance), `product_url` (string or `{variant|"default": url}`),
  `variants` (list of variant names — takes precedence over the `variant`
  column), `stacked` (render above parent counter), `clip_to_inner`.
- **Dimensions:** `start_anchor`, `end_anchor` (anchor specs — see
  `app/ARCHITECTURE.md` § Anchor resolution), `label_prefix`, rotation.
- **Labels:** text, offset from computed centroid, rotation, font size;
  room labels carry `source: 'room'`.

## Application state

| Table | Key columns | Purpose |
|-------|-------------|---------|
| `variants` | `id` PK, `name` UNIQUE, `label`, `source_variant`, `flags` (JSON), `layer_config` (JSON), `is_builtin` | 6 built-in variants (standard, minik, daybed, bare, sf, plumbing) + user-created clones. |
| `variant_exclusions` | (`variant`, `element_type`, `element_name`) PK | Per-variant hiding (e.g. bare/sf exclude IW6, RO5). |
| `views` | `name` PK, `label`, `script`, `svg_path`, `category`, `enabled` | Registered generator views (14 seeded) — see `GENERATORS.md`. |
| `config` | `key` PK, `value` | App key-value config (`default_variant`, `config_name`, roof style, …). |
| `undo_history` | `id` PK, `timestamp`, `action_type`, `before_state` (JSON), `after_state` (JSON), `description` | Persisted undo/redo stack (50 deep). |
| `plumbing_elements` | `id` PK, `type` ('supply_pipe'/'drain_pipe'/'fitting'/'fixture_connection'), `name` UNIQUE, `path` (JSON waypoints), `properties` (JSON), `fixture`, `sort_order` | Plumbing layout entities. |
| `roof_corners` | `seq` PK, `center_name`, `radiused`, `shortcut` | Roof perimeter corner treatment config (R-series). |
| `room_label_offsets` | `room_name` PK, `offset_e`, `offset_n` | Legacy — room label offsets now live on label elements; retained for migration. |

## Related files

- `app/adu.db` — the live database (binary; committed to git as the current
  working state)
- `app/configs/*.db` — named saved configurations (File → Save As / Load)
- `adu.db` (repo root) and `app/*.db` besides `adu.db` — historical/scratch
  snapshots, not read by the app
