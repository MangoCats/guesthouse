# ADU Editor — Requirements

Testable requirements for the interactive building editor application.
Each requirement has a unique ID, a description, and acceptance criteria
that can be verified by automated or manual testing.

Requirements are organized into a modular tree that groups functionality
both by **implementation layer** (data, API, UI) and by **user experience**
(what the user can do). New requirements identified from analysis of the
full commit history are marked **(NEW)**.

---

## 1  Data Layer

### 1.1  Database Schema & Seeding

#### DB-1  Schema Initialisation
The application SHALL create an SQLite database with tables `constants`,
`outline_chain`, `views`, `shapes`, `variant_exclusions`,
`room_label_offsets`, `undo_history`, `elements`, and `doors` when
launched for the first time.

**Acceptance:** Start with no `app/adu.db` file. Run `python run_app.py
--no-browser`. Verify the file is created and contains all six tables.

#### DB-2  Constants Seeding
On first initialisation the database SHALL contain every uppercase numeric
constant defined in `floorplan/constants.py`.

**Acceptance:** Compare the set of constant names in the database against
`dir(floorplan.constants)` filtered to uppercase float/int attributes.
Count must be >= 140.

#### DB-3  Constant Categories
Each constant SHALL be assigned to exactly one category from the set
{wall, interior_wall, opening, appliance, furniture, fixture, geometry,
construction, misc} based on its name prefix.

**Acceptance:** Query `SELECT DISTINCT category FROM constants` and verify
the result set equals the set above.

#### DB-4  Outline Chain Seeding
On first initialisation the database SHALL contain exactly 18 outline
chain rows matching the `OUTLINE_CHAIN` list in `floorplan/geometry.py`.

**Acceptance:** `SELECT count(*) FROM outline_chain` returns 18.
Each row's `seg_type` is one of {'L', 'CW', 'CCW'}. Each `end_name`
matches the corresponding entry in `CHAIN_POINT_NAMES`.

#### DB-5  Views Seeding
On first initialisation the database SHALL register at least 11 views
covering the floorplan, walls, span, survey, roof, plumbing, and site
plan generators.

**Acceptance:** `SELECT count(*) FROM views WHERE enabled = 1` >= 11.
Each row has a non-empty `script` path pointing to an existing Python file.

#### DB-6  Constant Update
`update_constant(name, value)` SHALL change the stored value and return
`True`. A subsequent `get_constants_dict()` SHALL reflect the new value.

**Acceptance:** Update `BED_WIDTH` to `7.0`. Read it back. Verify it
equals `7.0`. Reset constants afterwards.

#### DB-7  Batch Update
`update_constants_batch(updates)` SHALL update all named constants in a
single transaction and return the number of rows changed.

**Acceptance:** Batch-update `BED_WIDTH` and `BED_LENGTH` to new values.
Verify both changed. Return value equals 2.

#### DB-8  Reset Constants
`reset_constants()` SHALL delete all rows from the `constants` table and
re-seed from `floorplan/constants.py`, restoring original values.

**Acceptance:** Modify `BED_WIDTH`. Call `reset_constants()`. Read
`BED_WIDTH`. Verify it equals the value in the Python source module.

#### DB-9  Elements Table
The `elements` table SHALL store interior walls and user-added elements
with columns: `id`, `type` (wall / furniture / appliance / fixture /
clearance / label / dimension), `name`, `properties` (JSON), `variant`
(nullable — for floorplan variant membership).

**Variant membership:** A `NULL` variant means the element appears in all
variants.  A non-null value (e.g., `"standard"`) restricts the element to
that variant only.

Seeding populates rows for the 13 interior walls only (IW1–IW9, IW11–IW12,
IW2O, IW2S; no IW10).  Furniture, appliances, and fixtures are **not**
seeded — they are computed by the engine per-variant from constants and
layout logic (see ENG-13).  User-created elements (custom furniture,
additional walls, etc.) are stored here with absolute positions in
`properties` JSON.

**Dual-source model (transitional, Phases 3–11):** The canvas merges
engine-computed items (from constants/layout logic) with DB-stored custom
elements.  If a custom element shares a name with an engine-computed item,
the custom element takes precedence and the engine-computed item is hidden
for that variant.  This dual-source model is eliminated at Phase 12
cutover, when all items become database-stored elements with parametric
formulas.

**Acceptance:** Table exists with the specified columns. After seeding,
`SELECT count(*) FROM elements WHERE type = 'wall'` returns 13.
`SELECT count(*) FROM elements WHERE type != 'wall'` returns 0.

#### DB-10  Doors Table
The `doors` table SHALL store door configurations with columns: `id`,
`opening_name`, `width`, `hinge_side`, `swing_direction`, `door_type`
(single / double).

**Acceptance:** Table exists. Seeding populates rows for all openings that
have doors (RO1-RO7, appliance doors).

#### DB-11  Undo History Table
The database SHALL maintain an `undo_history` table storing serialised
snapshots of changed state, enabling undo and redo operations.

**Acceptance:** After an edit, `SELECT count(*) FROM undo_history` >= 1.
Each row contains a `timestamp`, `action_type`, and `before_state` /
`after_state` JSON.

#### DB-12  Variant Exclusions Seeding
On first initialisation the `variant_exclusions` table SHALL contain
exclusion rules for the `bare` and `sf` variants, hiding interior wall
`IW6` and rough opening `RO5` from those layouts.

**Acceptance:** `SELECT count(*) FROM variant_exclusions` returns 4.
`get_variant_exclusions("bare")` returns `{"wall": {"IW6"},
"rough_opening": {"RO5"}}`. The `standard` variant returns empty sets.

#### DB-13  Room Label Offsets
The `room_label_offsets` table SHALL store per-room (E, N) offsets from
computed centroids. Default is (0, 0) for all rooms (no rows).
`set_room_label_offset(name, e, n)` SHALL upsert a row.
`get_room_label_offsets()` SHALL return a dict of `{name: (e, n)}`.

**Acceptance:** Initially `get_room_label_offsets()` returns `{}`.
After `set_room_label_offset("BEDROOM", 0.5, -0.3)`,
`get_room_label_offsets()` returns `{"BEDROOM": (0.5, -0.3)}`.

### 1.2  Geometry Engine

#### ENG-1  Geometry Computation
`compute_geometry(constants_dict)` SHALL return a dict containing at least
the keys: `points`, `outline_segments`, `inner_segments`, `outline_poly`,
`inner_poly`, `interior_walls`, `outer_openings`, `rough_openings`,
`appliances`, `furniture`, `bbox`.

**Acceptance:** Call `compute_geometry` with the default constants dict.
Verify all listed keys exist and have non-empty values.

#### ENG-2  Point Count
The computed geometry SHALL include at least 50 named points covering
F-series, W-series, C-series, and survey points.

**Acceptance:** `len(result["points"]) >= 50`.

#### ENG-3  Outline Segment Count
The computed geometry SHALL include exactly 18 outline segments and 18
inner wall segments.

**Acceptance:** `len(result["outline_segments"]) == 18` and
`len(result["inner_segments"]) == 18`.

#### ENG-4  Interior Wall Count
The computed geometry SHALL include 13 interior walls for the standard
variant: IW1, IW2, IW2O, IW2S, IW3, IW4, IW5, IW6, IW7, IW8, IW9,
IW11, IW12. Variants with exclusion rules (bare, sf) SHALL have fewer
walls per their `variant_exclusions` entries (see ENG-14).

**Acceptance:** For standard: `set(result["interior_walls"].keys())`
equals the full 13-wall set. For bare/sf: `"IW6"` is absent.

#### ENG-5  Opening Counts
The computed geometry SHALL include exactly 12 outer openings and 7 rough
openings.

**Acceptance:** `len(result["outer_openings"]) == 12` and
`len(result["rough_openings"]) == 7`.

#### ENG-6  Appliance and Furniture Counts
The `appliances` dict SHALL contain keys `dryer`, `washer`, `counter`.
The `furniture` dict SHALL contain keys `bed`, `dresser`, `shelves`.

**Acceptance:** Verify both key sets.

#### ENG-7  Constant Propagation
Changing a constant in the database SHALL produce different computed
geometry when `compute_geometry` is called again.

**Acceptance:** Set `BED_WIDTH` to `80/12`. Recompute geometry. Verify the
bed bounding box width is 80 inches (approx. 6.6667 ft). Reset and verify
original width is 76 inches.

#### ENG-8  Derived Constants
When `WALL_OUTER` is patched, the engine SHALL recompute the derived
constants `WALL_EXTRA`, `AIR_GAP`, `DOOR_FLAT_FACE`, `F8F9_INNER_TURN_R`,
and `CORNER_SW_R` before computing geometry.

**Acceptance:** Patch `WALL_OUTER` to 10/12. Verify `AIR_GAP` equals
`WALL_OUTER - 2 * SHELL_THICKNESS`.

#### ENG-9  SVG File Reading
`get_svg_content(svg_path)` SHALL return the file contents as a string
when the file exists, and `None` when it does not.

**Acceptance:** Read `floorplan/floorplan.svg` (which exists after
`gen_all.py`). Verify result starts with `<` (XML). Read
`nonexistent.svg`. Verify result is `None`.

#### ENG-10  SVG Generation
`generate_svg(view_name, script_path)` SHALL run the generator script and
return `True` on success, `False` on failure or timeout.

**Acceptance:** Call with `("floorplan", "floorplan/gen_floorplan.py")`.
Verify it returns `True` and the SVG file is updated.

#### ENG-11  Outline Chain Mutation
The engine SHALL accept modified outline chain parameters (radius, sweep,
bearing, length) from the `outline_chain` database table and re-solve for
closure/tangency before recomputing geometry.  After Phase 5, the DB chain
is authoritative — the engine reads chain data from the DB rather than
from `floorplan/geometry.py`'s hardcoded chain.  All chain parameter types
(distances, bearings, arc radii, sweep angles) are editable.  "Reset to
Defaults" resets the chain (and all constants) to values seeded from the
existing scripts.  The app's outline solver and `floorplan/geometry.py`'s
solver must produce bit-identical results for default chain parameters;
for user-modified parameters, the app solver is authoritative.

**Acceptance:** Change F13-F14 arc radius from 28 to 30 inches. Engine
re-solves closure distances. Computed geometry reflects the new radius.
All tangency constraints hold.

#### ENG-12  Room Area Computation
The engine SHALL compute room areas (in square feet) for each enclosed
region defined by interior walls and the building outline. UTIL SHALL
be split into UTIL_N and UTIL_S by the E-W partition line.

**Acceptance:** `result["room_labels"]` (SF variant) contains entries for
BEDROOM, OFFICE, BATH, UTIL_N, UTIL_S, KITCHEN, LIVING, E CLOSET,
W CLOSET, STORAGE, WH with positive numeric `area` values.

#### ENG-13  Variant Item Computation
`compute_geometry(constants_dict, variant)` SHALL accept a `variant`
parameter (one of `standard`, `minik`, `daybed`, `bare`, `sf`) and return
a `variant_items` dict containing all furniture, appliances, and fixtures
for that variant. The standard variant SHALL produce at least 20 items
including `dryer`, `washer`, `hamper`, `counter`, `water_heater`,
`toilet_s`, `toilet_n`, `util_sink`, `bath_sink`, `stove`, `dishwasher`,
`kitchen_sink`, `fridge`, `ice_maker`, `work_counter`, `microwave`,
`north_counter`, `coffee_maker`, `dining_table`, `dining_chair_1`,
`dining_chair_2`, `bed`, `dresser`, `shelves`, `loveseat`, `et`,
`loveseat2`, `chair`, `ottoman`, `desk`, `desk_chair`. The `bare` and
`sf` variants SHALL produce zero items. The `minik` variant SHALL include
`sofa`, `rocker`, `cooktop`, `toaster` and exclude `loveseat`, `stove`,
`dishwasher`. The `daybed` variant SHALL include `daybed`, `shelves2`,
`et_east`, `et_west`, `rocker` and exclude `loveseat`, `sofa`.

Each item SHALL have keys: `type` (one of `appliance`, `furniture`,
`fixture`), `poly` (list of [E,N] coordinate pairs with >= 3 points),
`bbox` ({`w`, `s`, `e`, `n`}), `label` (display name), `shape` (`rect`
or `circle`). Circle items SHALL additionally have `center` ([E,N]) and
`radius` (positive float).

The response SHALL also include `variant` (the active variant name) and
`available_variants` (list of all variant names).

**Acceptance:** Call `compute_geometry` with default constants and each
variant name. Verify standard has >= 20 items with all required keys.
Verify bare and sf have 0 items. Verify minik contains `sofa` but not
`loveseat`. Verify daybed contains `daybed` but not `sofa`. Verify every
item bbox is within the building outline bbox (with 2 ft margin).

#### ENG-14  Variant Exclusion Filtering
When computing geometry for a variant with exclusion rules, the engine
SHALL omit the excluded interior walls from `interior_walls` and excluded
rough openings from `rough_openings`.

**Acceptance:** `compute_geometry(constants, "bare")` returns
`interior_walls` without key `"IW6"` and `rough_openings` without key
`"RO5"`. `compute_geometry(constants, "standard")` includes both.

#### ENG-15  Room Label Computation
The engine SHALL compute room labels for 11 rooms: BEDROOM, UTIL_N,
UTIL_S, KITCHEN, LIVING, BATH, OFFICE, E CLOSET, W CLOSET, STORAGE, WH.
Each label SHALL have `name`, `pos` (area-weighted centroid plus DB
offset), and `centroid` (raw centroid). For the SF variant, each label
SHALL additionally include `area` (numeric sf) and `poly` (room boundary
polygon for highlight rendering). The response SHALL include
`room_labels` (list of label dicts) and optionally `sf_lines` (dashed
partition line endpoints for SF layout).

**Acceptance:** `compute_geometry(constants, "standard")["room_labels"]`
has 11 entries. Each has `name`, `pos`, `centroid` keys. All `pos`
values fall within the building outline bbox.
`compute_geometry(constants, "sf")["room_labels"]` entries additionally
have `area` and `poly` keys. `sf_lines` is present and non-empty.

### 1.3  File Generation

#### GEN-1  Regenerate Single View
File > Regenerate Current View SHALL run the generator script for the
active SVG view and reload its content in the viewport.

**Acceptance:** Switch to the Floorplan tab. Click Regenerate Current
View. The SVG content refreshes. A toast confirms completion.

#### GEN-2  Regenerate All
File > Regenerate All SHALL run all registered generator scripts and
broadcast SVG update events.

**Acceptance:** Click Regenerate All. All generator scripts run. Toasts
appear for each updated view.

#### GEN-3  Reset to Defaults
File > Reset to Defaults SHALL prompt for confirmation, then reset all
constants and the outline chain, and reload the Constants and Outline tables.

**Acceptance:** Click Reset to Defaults. A confirmation dialog appears.
Confirm. Constants table refreshes with original values. Outline table
refreshes with seed values from floorplan/geometry.py.

#### GEN-4  Interactive View Not Regenerable
Clicking Regenerate Current View while on the Interactive tab SHALL show
a toast explaining that the interactive view cannot be regenerated.

**Acceptance:** Select Interactive tab. Click Regenerate Current View.
Toast reads "Cannot regenerate interactive view".

### 1.4  Anchored Dimensions

Dimension elements support anchor-based endpoint resolution, where each
endpoint is defined by an anchor spec that resolves to absolute
coordinates during geometry computation. Builtin dimensions are seeded
as database elements and re-anchored on every geometry recomputation.

#### DIM-1  Dimension Anchor Properties
Dimension elements SHALL support `start_anchor` and `end_anchor`
properties that resolve to absolute coordinates during geometry
computation. Each anchor is a JSON object specifying the anchor type
and its parameters. The resolved coordinates are stored as `start` and
`end` ([E, N]) in the dimension's properties after computation.

**Acceptance:** Query `SELECT properties FROM elements WHERE type =
'dimension'`. Each row's JSON contains `start_anchor` and `end_anchor`
objects. After geometry computation, each also contains `start` and
`end` coordinate arrays that match the resolved anchor positions.
**Tested:** Implemented in `app/engine.py` dimension anchor resolver.

#### DIM-2  Anchor Type Repertoire
The anchor resolver SHALL support the following anchor types:
- `point`: resolves to a named geometry point (F-series, W-series, etc.)
- `wall_face`: resolves to a point on a specified face (N/S/E/W) of a
  named interior wall
- `opening_face`: resolves to a point on a specified face of a named
  opening
- `line_intersection`: resolves to the intersection of two infinite
  lines (see DIM-4)
- `computed`: resolves to a previously computed coordinate stored in the
  anchor spec

**Acceptance:** Create dimension elements using each anchor type. After
geometry computation, verify that the resolved `start` and `end`
coordinates match the expected positions for each anchor type.
**Tested:** Implemented in `app/engine.py` `resolve_anchor()`.

#### DIM-3  Builtin Dimension Seeding
Builtin dimensions SHALL be seeded as elements with `type='dimension'`
and `source='builtin'` during database initialisation and reset. The
seeded dimensions correspond to the standard measurement annotations
produced by the existing floorplan generator scripts. On "Reset to
Defaults", all builtin dimension elements are deleted and re-seeded
from the canonical definition list.

**Acceptance:** After a fresh database init, `SELECT count(*) FROM
elements WHERE type = 'dimension' AND source = 'builtin'` returns
>= 18. After modifying a builtin dimension and clicking "Reset to
Defaults", the dimension reverts to its seeded anchor configuration.
**Tested:** Implemented in `app/database.py` seeding and
`app/engine.py` dimension computation.

#### DIM-4  Line Intersection Anchor
The `line_intersection` anchor type SHALL compute the intersection of
two infinite lines. Each line is defined by a point spec (a named
geometry point or [E, N] coordinate) and a direction spec (a bearing
angle in degrees, a named geometry point defining the direction from
the line's point, or an axis keyword like `"E-W"` or `"N-S"`). The
resolver returns the intersection coordinate, or skips the dimension
if the lines are parallel.

**Acceptance:** Create a dimension with a `line_intersection` anchor
where line 1 passes through W5 bearing east and line 2 passes through
W9 bearing north. The resolved coordinate is the intersection of those
two infinite lines. Verify the coordinate is correct by manual
calculation.
**Tested:** Implemented in `app/engine.py` `resolve_anchor()`.

#### DIM-5  Dimension Visual Style
Dimension elements SHALL support a `dim_style` property controlling
per-dimension visual appearance. Supported styles are `"solid"` (the
default, rendering a solid line with tick marks) and `"dashed"`
(rendering a dashed line). The style is stored in the element's
`properties` JSON and applied during canvas rendering.

**Acceptance:** Create two dimensions: one with `dim_style: "solid"` and
one with `dim_style: "dashed"`. On the canvas, the solid dimension
renders with a continuous line; the dashed dimension renders with a
dashed stroke pattern. Verify via SVG inspection that the dashed
dimension's `<line>` element has a `stroke-dasharray` attribute.
**Tested:** Implemented in `app/static/js/app.js` dimension rendering.

### 1.5  Variant Filtering

Elements support fine-grained variant membership control via a
`properties.variants` array, complementing the top-level `variant`
column in the `elements` table.

#### VAR-1  Element Variant Array
Elements SHALL support a `properties.variants` array specifying which
layout variants the element appears in. The array contains variant name
strings (e.g., `["standard", "minik"]`). When `properties.variants` is
absent or null, the element's visibility is governed by the top-level
`variant` column (DB-9 rules: `NULL` = all variants, non-null = that
variant only).

**Acceptance:** Create an element with `properties.variants =
["standard", "daybed"]`. When viewing the `standard` variant, the
element is visible. When viewing the `minik` variant, the element is
hidden. When viewing the `daybed` variant, the element is visible.
**Tested:** Implemented in `app/engine.py` variant filtering logic.

#### VAR-2  Properties Variants Precedence
When `properties.variants` is set (non-null, non-empty array), it
SHALL take precedence over the top-level `variant` column for
determining element visibility. The `variant` column value is ignored
for filtering purposes when `properties.variants` is present.

**Acceptance:** Create an element with `variant = "standard"` and
`properties.variants = ["minik", "daybed"]`. When viewing the
`standard` variant, the element is hidden (despite the `variant` column
matching). When viewing the `minik` variant, the element is visible.
**Tested:** Implemented in `app/engine.py` variant filtering logic.

#### VAR-3  Layout Checkbox UI
The Properties panel for elements SHALL display layout checkboxes
allowing users to select any subset of layout variants in which the
element should appear. The checkboxes SHALL reflect the current
`properties.variants` array (or all-checked if the array is absent/null
and `variant` is NULL). Toggling a checkbox SHALL update the
`properties.variants` array via `PUT /api/elements/<id>` and trigger
geometry recomputation for affected variants.

**Acceptance:** Select a furniture element. The Properties panel shows
checkboxes for Standard, Small Kitchen, Daybed, Room Dimensions, and
Square Footage. Uncheck "Daybed". The element disappears when switching
to the Daybed variant. Re-check "Daybed". The element reappears.
**Tested:** Implemented in `app/static/js/app.js` properties panel.

---

## 2  REST API

### 2.1  Constants API

#### API-1  GET /api/constants
SHALL return HTTP 200 with a JSON array of constant objects. Each object
SHALL have keys: `name`, `value`, `expr`, `unit`, `category`,
`description`.

**Acceptance:** Response is a JSON array with length >= 140. First element
has all six keys.

#### API-2  GET /api/constants?category=wall
SHALL return only constants whose `category` equals the query parameter.

**Acceptance:** Every object in the response has `category == "wall"`.

#### API-3  GET /api/constants/categories
SHALL return HTTP 200 with a JSON array of distinct category strings.

**Acceptance:** Response is a JSON array containing at least the strings
`"wall"`, `"opening"`, `"furniture"`.

#### API-4  PUT /api/constants/<name>
SHALL accept a JSON body `{"value": <number>}`, update the constant, and
return `{"ok": true, "name": ..., "value": ...}`.

**Acceptance:** PUT `BED_WIDTH` with `{"value": 7.0}`. Response status is
200, `ok` is `true`, `value` is `7.0`.

#### API-5  PUT /api/constants/<name> -- validation
SHALL return HTTP 400 if `value` is missing or not a valid number.
SHALL return HTTP 404 if the constant name does not exist.

**Acceptance:** PUT with `{}` body returns 400. PUT with
`{"value": "abc"}` returns 400. PUT to `/api/constants/NONEXISTENT`
returns 404.

#### API-6  PUT /api/constants/batch
SHALL accept `{"updates": {"A": 1.0, "B": 2.0}}` and return
`{"ok": true, "changed": N}` where N is the count of rows updated.

**Acceptance:** Batch-update two known constants. Verify `changed == 2`.

#### API-7  POST /api/constants/reset
SHALL reset all constants to their original source values and return
`{"ok": true}`.

**Acceptance:** Modify a constant. POST reset. GET the constant. Value
matches the original.

### 2.2  Geometry & Views API

#### API-8  GET /api/geometry
SHALL return HTTP 200 with a JSON object containing the complete computed
geometry (see ENG-1 for required keys). SHALL accept an optional
`?variant=` query parameter (one of `standard`, `minik`, `daybed`,
`bare`, `sf`; default `standard`). An unrecognised variant name SHALL
fall back to `standard`. The response SHALL include `variant` (the active
variant name), `variant_items` (dict of items for that variant), and
`available_variants` (list of all variant names).

**Acceptance:** GET `/api/geometry` returns JSON with keys `points`,
`outline_poly`, `interior_walls`, `bbox`, `variant`, `variant_items`,
`available_variants`. `variant` equals `"standard"` when no param given.
GET `/api/geometry?variant=bare` returns `variant == "bare"` with 0
variant_items. GET `/api/geometry?variant=nonexistent` returns
`variant == "standard"`.

#### API-9  GET /api/geometry -- error handling
SHALL return HTTP 500 with `{"error": "..."}` if geometry computation
fails.

**Acceptance:** (Edge case -- verified by code inspection of the try/except
block.)

#### API-10  GET /api/outline
SHALL return HTTP 200 with a JSON array of 18 outline chain segment
objects.

**Acceptance:** Response length is 18. Each object has keys `seq`,
`seg_type`, `end_name`.

#### API-11  GET /api/views
SHALL return HTTP 200 with a JSON array of enabled view objects. Each has
keys `name`, `label`, `script`, `svg_path`, `category`.

**Acceptance:** Response length >= 11. Each element has all five keys.

#### API-12  GET /api/svg/<view_name>
SHALL return the SVG file content with MIME type `image/svg+xml` when the
file exists. SHALL return HTTP 404 when the SVG has not been generated or
the view name is unknown.

**Acceptance:** GET `/api/svg/floorplan` returns 200 with content starting
with `<`. GET `/api/svg/nonexistent` returns 404.

#### API-13  GET /api/svg/<view_name>/file
SHALL serve the file as an attachment. SVG files use MIME `image/svg+xml`;
PDF files use MIME `application/pdf`.

**Acceptance:** GET `/api/svg/floorplan/file` returns 200. Content-Type
header contains `svg`.

#### API-14  POST /api/regenerate
With no body: SHALL regenerate all enabled views and return
`{"ok": true, "results": {...}}`.
With body `{"view": "floorplan"}`: SHALL regenerate only that view.

**Acceptance:** POST with `{"view": "floorplan"}` returns 200 with
`ok == true` and `view == "floorplan"`. POST with no body returns 200 with
`results` dict.

#### API-15  Geometry Cache Invalidation
After a successful PUT to `/api/constants/<name>`, the geometry cache
SHALL be marked dirty for all variants so the next GET `/api/geometry`
recomputes.

**Acceptance:** GET geometry (cache populated). PUT a constant. GET
geometry again. Verify the response reflects the changed constant.

#### API-34  GET /api/variants
SHALL return HTTP 200 with a JSON array of variant objects. Each object
SHALL have keys `name` (string identifier) and `label` (human-readable
display name). The array SHALL contain exactly 5 entries: `standard`,
`minik`, `daybed`, `bare`, `sf`.

**Acceptance:** GET `/api/variants` returns 200 with a JSON array of
length 5. The set of `name` values equals
`{"standard", "minik", "daybed", "bare", "sf"}`. Each entry has a
non-empty `label`.

### 2.3  Outline Chain Mutation API

After Phase 5, the `outline_chain` database table is the authoritative
source for all chain parameters.  The engine reads chain data from the DB
and injects it as patched constants before reloading floorplan modules.
All chain parameter types (distances, bearings, arc radii, sweep angles)
are editable.  "Reset to Defaults" resets the chain to the values seeded
from `floorplan/geometry.py`.

#### API-16  PUT /api/outline/<seq>
SHALL accept a JSON body with optional keys `dist_or_radius`, `sweep`,
`bearing` and update the specified outline chain segment in the database.
SHALL trigger closure re-solve and geometry recomputation.

**Acceptance:** PUT `/api/outline/5` with `{"dist_or_radius": 30.0}`.
Response contains the updated segment and `closure_valid: true`. GET
`/api/geometry` reflects the changed arc radius.

#### API-17  POST /api/outline/validate
SHALL accept a proposed set of outline chain parameter changes and return
whether closure/tangency constraints are satisfied without committing.

**Acceptance:** POST with a valid parameter set returns
`{"valid": true, "closure_error": 0.0}`. POST with an impossible set
returns `{"valid": false, "closure_error": <nonzero>}`.

#### API-18  POST /api/outline/add-point
SHALL insert a new F-point into the outline chain at the specified
position (splitting an existing segment).

**Acceptance:** POST to add F11a between F11 and F12. Chain length
increases by 1. Closure is re-solved.

#### API-19  DELETE /api/outline/<seq>
SHALL remove a point from the outline chain and re-solve for closure.

**Acceptance:** DELETE `/api/outline/3` (removing F4). Chain length
decreases by 1. Closure is re-solved.

### 2.4  Element CRUD API

#### API-20  POST /api/elements
SHALL create a new element (interior wall, furniture, appliance, fixture)
and return the created element with its assigned ID.

**Acceptance:** POST `{"type": "furniture", "name": "HAMPER",
"properties": {"width": 18, "depth": 14}}`. Response status 201 with
assigned `id`.

#### API-21  PUT /api/elements/<id>
SHALL update the properties of an existing element (position, size,
rotation, name).

**Acceptance:** PUT with `{"properties": {"x": 5.0, "y": -3.0}}`.
Response 200 with updated element.

#### API-22  DELETE /api/elements/<id>
SHALL remove an element and any dependent elements (e.g., removing a wall
removes its openings). SHALL return the list of removed element IDs.

**Acceptance:** DELETE an interior wall that hosts a rough opening. Both
the wall and its opening are removed. Response includes both IDs.

#### API-23  POST /api/elements/<id>/move
SHALL accept `{"dx": <feet>, "dy": <feet>}` or
`{"anchor": "<element>", "face": "<side>", "offset": <inches>}` and
reposition the element.

**Acceptance:** POST move IW1 by `{"dx": 0.5, "dy": 0}` (6 inches east).
Geometry recomputes with IW1 at the new position.

#### API-24  POST /api/openings
SHALL create a new opening on a specified wall segment with position and
width parameters.

**Acceptance:** POST `{"name": "O8a", "segment": "F18-F1", "width": 19,
"offset": 48}`. Response 201. Opening appears in computed geometry.

#### API-25  PUT /api/openings/<name>
SHALL update opening width, position, or type.

**Acceptance:** PUT O8 with `{"width": 25}`. Geometry recomputes with
the new opening width.

#### API-26  DELETE /api/openings/<name>
SHALL remove an opening and its associated door (if any).

**Acceptance:** DELETE O8a. Opening no longer appears in computed geometry.

### 2.5  Door API

#### API-27  POST /api/doors
SHALL create a door on an opening with hinge side, swing direction, width,
and type (single/double).

**Acceptance:** POST `{"opening": "RO1", "hinge_side": "east",
"swing_direction": "south", "width": 42, "door_type": "single"}`.
Response 201. Door appears in computed geometry.

#### API-28  PUT /api/doors/<opening_name>
SHALL update door hinge side, swing direction, or type.

**Acceptance:** PUT RO1 door with `{"hinge_side": "west"}`. Door arc
flips in computed geometry.

#### API-29  DELETE /api/doors/<opening_name>
SHALL remove the door from an opening.

**Acceptance:** DELETE RO1 door. Opening remains but door arc is removed.

### 2.6  Undo/Redo API

#### API-30  POST /api/undo
SHALL revert the most recent edit operation and return the restored state.

**Acceptance:** Edit a constant. POST undo. Constant returns to its
previous value. Geometry recomputes.

#### API-31  POST /api/redo
SHALL re-apply the most recently undone operation.

**Acceptance:** Undo an edit. POST redo. The edit is re-applied. Geometry
recomputes.

### 2.7  Real-Time Events

#### API-32  GET /api/events (SSE)
SHALL return a streaming response with MIME type `text/event-stream`.
The first event SHALL be `event: connected`.

**Acceptance:** Connect to the endpoint. Read at least the first event.
Verify it is `event: connected\ndata: {}\n\n`.

#### API-33  GET / (root)
The root route SHALL return HTTP 200 with HTML containing the string
"ADU Editor".

**Acceptance:** `GET /` returns status 200. Body contains `"ADU Editor"`.

---

## 3  User Interface -- Layout & Navigation

### 3.1  Page Structure

#### UI-1  Five-Region Layout
The page SHALL display a menu bar (top), tool palette (left), viewport
(centre), property/data panel (right), and status bar (bottom).

**Acceptance:** Load the page. All five regions are visible and occupy
their expected positions.

#### UI-2  Menu Bar
The menu bar SHALL have menus: File, Edit, View, Tools. Each opens a
dropdown on hover.

**Acceptance:** Hover over each menu name. A dropdown with at least one
action button appears.

### 3.2  View Management

#### UI-3  View Tabs
The viewport area SHALL display a tab bar. The first tab SHALL be
"Interactive". Additional tabs SHALL correspond to registered SVG views.

**Acceptance:** Load the page. The Interactive tab is visible and active.
At least one additional tab (e.g. "Floorplan") is visible.

#### UI-4  Tab Switching
Clicking a view tab SHALL switch the viewport content. Clicking
"Interactive" shows the SVG canvas. Clicking a generated view tab loads
and displays the corresponding SVG file.

**Acceptance:** Click "Interactive" -- canvas is visible, SVG container is
hidden. Click "Floorplan" -- canvas is hidden, SVG container shows
floorplan SVG or a "not generated" message.

#### UI-5  Create View Variant **(NEW)**
View > New Variant SHALL create a named copy of the current floorplan
layout as a separate variant (e.g., "daybed", "minik").

**Acceptance:** Create a variant named "test". A new tab appears. The
variant initially matches the source layout. Edits to the variant do not
affect the original.

#### UI-6  Per-View Layer Configuration **(NEW)**
Each view variant SHALL have independently configurable layer visibility
(walls, openings, furniture, labels, dimensions).

**Acceptance:** In variant "bare", disable furniture and labels. The
generated SVG for "bare" omits those layers. The main floorplan still
shows them.

#### UI-8  Variant Selector
The view tab bar SHALL contain a "Layout" dropdown selector visible when
the Interactive or Floorplan tab is active. The selector SHALL list all
5 floorplan variants (Standard, Small Kitchen, Daybed, Room Dimensions,
Square Footage). On the Interactive tab, selecting a variant SHALL
reload the geometry with that variant's furniture/appliance set. On the
Floorplan tab, selecting a variant SHALL load the corresponding
variant-specific SVG file (standard → `floorplan.svg`, minik →
`floorplan_minik.svg`, daybed → `floorplan_db.svg`, bare →
`floorplan_bare.svg`, sf → `floorplan_sf.svg`). The selector SHALL be
hidden when any other SVG view tab is active.

**Acceptance:** Load the page with the Interactive tab active. The
"Layout" dropdown is visible with 5 options. Select "Small Kitchen" --
the canvas re-renders with minik furniture (sofa, rocker, cooktop
visible; loveseat, stove absent). Switch to the Floorplan tab -- the
dropdown remains visible. Select "Daybed" -- the floorplan SVG switches
to `floorplan_db.svg`. Switch to a non-floorplan SVG view tab -- the
dropdown is hidden.

### 3.3  Right Panel

#### UI-7  Right Panel Tabs
The right panel SHALL have tabs: Properties, Constants, Outline, Openings,
Elements. Clicking a tab shows only that panel's content.

**Acceptance:** Click each tab. The corresponding panel content is visible;
all others are hidden.

---

## 4  Interactive Canvas

### 4.1  Rendering

#### CV-1  Outline Rendering
The interactive canvas SHALL render the building outline as a filled
polygon with a stroke border.

**Acceptance:** Load the page. A polygon shape is visible in the canvas
area. The polygon has >= 40 vertices (arcs are polygonised).

#### CV-2  Interior Wall Rendering
The canvas SHALL render all 13 interior walls as filled polygons.

**Acceptance:** Count elements in `#layer-walls`. There are 13 polygon
elements.

#### CV-3  Opening Rendering
When the "Openings" toggle is checked, the canvas SHALL render all 12
outer openings and 7 rough openings.

**Acceptance:** Count elements in `#layer-openings`. There are 19 polygon
elements.

#### CV-4  Furniture Rendering
When the "Furniture" toggle is checked, the canvas SHALL render all
variant items (appliances, furniture, and fixtures) for the currently
selected variant. The standard variant SHALL render at least 20 items.
Items with `shape == "rect"` SHALL be rendered as SVG `<polygon>`
elements. Items with `shape == "circle"` SHALL be rendered as SVG
`<circle>` elements. Each item type SHALL have a distinct CSS class:
`item-appliance`, `item-furniture`, `item-fixture`. Stacked items
(e.g., MICRO on counter, coffee maker on counter) SHALL be rendered
after their parent items so they appear on top in SVG paint order.

**Acceptance:** Select the standard variant. Count elements in
`#layer-furniture`. There are >= 20 elements (polygons and circles).
Select the bare variant. The layer contains 0 elements. Select minik.
The layer contains elements for `sofa` and `cooktop` but not `loveseat`
or `stove`.

#### CV-5  Point Markers
When the "Points" toggle is checked, the canvas SHALL render circle
markers for F-series, W-series, and C-series points.

**Acceptance:** `#layer-points` contains circle elements. F-series circles
use the accent colour; W-series circles use green.

#### CV-6  Point Labels
When the "Labels" toggle is checked, the canvas SHALL render text labels
adjacent to each point marker showing the point name.

**Acceptance:** `#layer-labels` contains text elements with content such
as `"F1"`, `"W2"`, `"C5"`.

#### CV-7  Door Arc Rendering
When a door is configured on an opening, the canvas SHALL render the door
swing arc showing hinge position and sweep direction.

**Acceptance:** An opening with a configured door displays a dashed arc
from closed to open position. Hinge point is marked with a small circle.
**Tested:** `test_zapp_canvas.py::TestDoorArcs` (8 tests).

#### CV-8  Room Labels
When the "Room Names" toggle is checked, the canvas SHALL render room
name labels for all 11 rooms (BEDROOM, UTIL_N, UTIL_S, KITCHEN, LIVING,
BATH, OFFICE, E CLOSET, W CLOSET, STORAGE, WH) positioned at their
area-weighted centroids (plus any stored DB offset and rotation).

**Phases 0–7:** Labels are auto-computed from room polygons.  The
`room_label_offsets` table stores `(de, dn)` offsets from centroids.

**Phase 8 transition:** Labels migrate to the `elements` table (type
`'label'`) with offset, rotation, font size, and text stored per-element.
The `room_label_offsets` table is deprecated.  The 11 auto-computed room
labels become editable label elements.  Users can also add additional
custom labels beyond the 11 defaults (LABEL-1 through LABEL-4).
"Reset to Defaults" back-computes label positions from the existing SVG
generation scripts' output and stores the corresponding offsets and
rotations.

**Acceptance:** `#layer-rooms` contains text elements. Each label is
positioned within its room boundary polygon. All 11 rooms have labels.

#### CV-9  Dimension Lines
When the "Dimensions" toggle is checked, the canvas SHALL render all
dimension lines computed by `compute_dimension_endpoints()` (at least 18
for standard variant). Each dimension SHALL be rendered as a line
segment in `#layer-dims` with perpendicular tick marks at each endpoint
and a rotated text label at the midpoint showing the distance in
feet-inches format (NF-6 compliant). Label rotation SHALL be normalised
to [-90, 90) degrees for readability.

The geometry API response SHALL include a `dimensions` dict keyed by
dimension name, each containing `A` and `B` endpoint coordinates and
`dist` (distance in feet).

**Acceptance:** Check the "Dims" toggle. `#layer-dims` contains >= 54
child elements (18 dimensions x 3 elements each: 1 line + 2 ticks,
plus 18 text labels). Each text label shows a value in `X' Y.YY"`
format. GET `/api/geometry` response includes `dimensions` with >= 18
entries, each having `A`, `B`, and `dist` keys.

#### CV-10  Area Labels
When the SF layout variant is selected, the canvas SHALL render room area
labels (in square feet) below each room name label. Clicking an area
label SHALL toggle highlighting of the corresponding room polygon.

**Acceptance:** Select the SF variant. `#layer-rooms` contains area text
elements showing values like "125 sf". Click an area label — the room
polygon highlights with a semi-transparent fill and stroke. Click again
— the highlight is removed.

#### CV-10a  SF Partition Lines
When the SF layout variant is selected, the canvas SHALL render dashed
partition lines showing the boundaries used for area calculation.

**Acceptance:** Select the SF variant. `#layer-rooms` contains `<line>`
elements with class `sf-partition` and dashed stroke style.

#### CV-11  Clearance Zones
The canvas SHALL render clearance zones as dashed polygons extending
from item faces when the "Clearance" toggle is checked.

**Acceptance:** Clearance zones appear at dresser, stove, dishwasher,
and hamper locations as dashed rectangles extending from the item face.
**Tested:** `test_zapp_canvas.py::TestClearanceZones` (3 tests),
`TestApplianceClearanceZones` (7 tests).

#### CV-12  Hyperlink Indicators
Elements with attached product URLs SHALL display a small link icon
overlay on the canvas.

**Acceptance:** A furniture item with a URL shows a link icon at its
corner. Clicking the icon opens the URL in a new tab.

**Tested:** `test_zapp_style.py::TestProductUrl` (URL storage),
`test_zapp_style.py::TestStyleAPI::test_put_product_url` (API round-trip).

### 4.2  Canvas Navigation

#### CV-13  Fit to Window
On initial load or when the user presses F or selects View > Fit to
Window, the building SHALL be centred and scaled to fill the viewport with
a margin.

**Acceptance:** Load the page. The outline is fully visible within the
viewport. No scrolling is needed.

#### CV-14  Mouse Wheel Zoom
Scrolling the mouse wheel SHALL zoom the canvas toward or away from the
cursor position.

**Acceptance:** Place cursor over the building. Scroll up -- the building
gets larger. Scroll down -- it gets smaller. The area under the cursor
remains roughly stationary.

#### CV-15  Pan with Middle Mouse Button
Pressing and dragging with the middle mouse button SHALL pan the canvas.

**Acceptance:** Middle-click and drag. The building moves with the cursor.

#### CV-16  Pan Tool
When the Pan tool is active, left-click drag SHALL pan the canvas. The
cursor SHALL be `grab` (at rest) or `grabbing` (while dragging).

**Acceptance:** Select Pan tool (click button or press H). Left-drag
pans the canvas. Cursor changes appropriately.

#### CV-17  Coordinate Display
As the mouse moves over the viewport, the status bar SHALL display the
current cursor position in world coordinates (feet and inches).

**Acceptance:** Move mouse. `#coord-display` updates continuously showing
E and N values.

### 4.3  Display Toggles

#### DIS-1  Points Toggle
Unchecking the "Points" checkbox SHALL hide all point markers and
re-render the canvas without them.

**Acceptance:** Uncheck Points. `#layer-points` is empty. Check Points.
Markers reappear.

#### DIS-2  Labels Toggle
Unchecking the "Labels" checkbox SHALL hide all point labels.

**Acceptance:** Uncheck Labels. `#layer-labels` is empty.

#### DIS-3  Grid Toggle
Checking the "Grid" checkbox SHALL show a 1-foot grid overlay.

**Acceptance:** Check Grid. The `#grid-rect` element's `visibility`
attribute changes to `visible`.

#### DIS-4  Openings Toggle
Unchecking the "Openings" checkbox SHALL hide all opening polygons.

**Acceptance:** Uncheck Openings. `#layer-openings` is empty.

#### DIS-5  Furniture Toggle
Unchecking the "Furniture" checkbox SHALL hide all appliance and furniture
polygons.

**Acceptance:** Uncheck Furniture. `#layer-furniture` is empty.

#### DIS-6  Room Names Toggle
Checking the "Room Names" checkbox SHALL show/hide room name labels,
area labels (SF variant), and SF partition lines. The checkbox is
checked by default.

**Acceptance:** Toggle the "Room Names" checkbox. `#layer-rooms` is
populated when checked, empty when unchecked. On the SF variant, area
labels and dashed partition lines also appear/disappear with this toggle.

#### DIS-7  User Dimensions Toggle **(NEW)**
Checking the "User Dims" checkbox SHALL show/hide user-created persistent
dimension annotations (placed via TL-11).  This is distinct from the
existing "Dims" toggle (CV-9), which controls engine-computed dimension
lines.  Both toggles operate independently.

**Acceptance:** Create a user dimension via TL-11.  Toggle "User Dims"
off — user-created dimension disappears; engine-computed dimensions
(controlled by CV-9 toggle) remain visible.  Toggle "User Dims" on —
user dimension reappears.

#### DIS-9  Clearance Toggle
Checking the "Clearance" checkbox SHALL show/hide clearance zones.

**Acceptance:** Toggle checkbox. Clearance zone polygons appear/disappear.
**Tested:** Client-side toggle (manual verification).

#### DIS-10  Doors Toggle
Checking the "Doors" checkbox SHALL show/hide door swing arcs.

**Acceptance:** Toggle checkbox. Door arc elements appear/disappear.
**Tested:** Client-side toggle (manual verification).

#### DIS-11  Open Links Toggle
The "Open Links" checkbox SHALL be checked by default. When unchecked,
canvas elements SHALL NOT be wrapped in `<a xlink:href>` tags and link
icon overlays (CV-12) SHALL be hidden, suppressing link-opening on click
in the interactive view. The Properties panel "Open" button is unaffected.

**Acceptance:** Uncheck "Open Links". Click an item with a product URL —
no new tab opens; link icon is not visible. Check "Open Links" — link
icons reappear and clicking opens the URL.
**Tested:** Client-side toggle (manual verification).

#### DIS-12  Areas Toggle
The "Areas" checkbox SHALL be unchecked by default. When checked, room
area labels (in square feet) SHALL appear within each room boundary on
the interactive canvas. On the SF variant, areas are always shown
regardless of this toggle.

**Acceptance:** Check "Areas" on Standard variant. Room areas appear.
Uncheck — areas disappear. Switch to SF — areas always shown.
**Tested:** `test_zapp_analysis.py::TestRoomAreas` (engine emits areas
for all variants).

---

## 5  Selection & Properties

### 5.1  Selection

#### SEL-1  Click-to-Select
Clicking on a wall, opening, appliance, or furniture polygon SHALL select
that element. The element SHALL be visually highlighted and the Properties
panel SHALL populate with its details.

**Acceptance:** Click on an interior wall polygon. It receives the CSS
class `selected-highlight`. The Properties panel title shows the wall
name.

#### SEL-2  Point Selection
Clicking on a point marker SHALL select the point and show its Easting
and Northing in feet and inches.

**Acceptance:** Click on a point circle. Properties panel shows E, N,
and their inch equivalents.

#### SEL-3  Clear Selection
Clicking on the viewport background or pressing Escape SHALL deselect any
selected element and clear the Properties panel.

**Acceptance:** Select an element. Click background. The
`selected-highlight` class is removed. Properties panel shows the empty
state.

#### SEL-4  Multi-Select
Shift-clicking or drag-selecting SHALL allow multiple elements to be
selected simultaneously for group operations (move, delete).

**Acceptance:** Shift-click or Ctrl+Click adds to multi-selection. Rubber-band
drag-select: click empty space + drag creates selection rectangle; mouseup
selects intersecting elements. All selected elements highlighted yellow.
**Tested:** Client-side (rubber-band logic + shift-click handler in app.js).

### 5.2  Property Display

#### SEL-5  Wall Properties
When a wall is selected, the Properties panel SHALL display its width
(inches), height (inches), and bounding box coordinates (feet). It SHALL
also list related constants with editable input fields.

**Acceptance:** Select IW1. Panel shows Width, Height, West, South, East,
North values. Constants starting with "IW1" appear as editable rows.

#### SEL-6  Opening Properties
When an opening is selected, the Properties panel SHALL display its name,
segment reference, width, and actual computed width.

**Acceptance:** Select O3. Panel shows Name, Segment (F2-F5), and Width
values.

#### SEL-7  Door Properties
When an opening with a door is selected, the Properties panel SHALL
display door width, hinge side, swing direction, and door type, all
editable.

**Acceptance:** Select RO1. Panel shows Door section with Hinge Side
dropdown (east/west/north/south), Swing Direction dropdown, and Type
dropdown (single/double).

#### SEL-8  Furniture Properties
When a furniture or appliance item is selected, the Properties panel
SHALL display its width, depth, rotation angle, position (wall-relative
offset), and product URL (if any).

**Acceptance:** Select BED. Panel shows Width, Depth (inches), Rotation
(degrees), Wall Reference, Offset, and URL field.

### 5.3  Property Editing

#### SEL-9  Inline Constant Editing
Editing a related-constant value in the Properties panel and pressing
Enter SHALL send a PUT request to update the constant and trigger geometry
recomputation.

**Acceptance:** Select a wall. Edit a related constant input. A PUT
request is sent. A success toast appears. The geometry is reloaded.

#### SEL-10  Opening Width Editing
The opening width field in the Properties panel SHALL be editable.
Changing it SHALL update the opening width constant and trigger geometry
recomputation.

**Acceptance:** Select opening in properties panel. `findWidthConstant()` maps
opening name to controlling constant (O1_WIDTH, O7_HALF_WIDTH, IW1_RO_WIDTH,
etc). Editable input calls `handleConstantEdit()` on change.
**Tested:** test_zapp_tools.py::TestDT7OpeningWidth (8 tests).

#### SEL-11  Door Property Editing
Changing a door property (hinge side, swing direction) in the Properties
panel SHALL update the door configuration and re-render the door arc.

**Acceptance:** Select RO3. Change hinge side from "south" to "north" via
dropdown. Door arc flips to the opposite side.
**Tested:** `test_zapp_canvas.py::TestDoorArcAPI::test_door_invalidation_via_api`.

#### SEL-12  Product URL Field
When a furniture, appliance, or fixture item is selected, the Properties
panel SHALL display a "Link" text field. The field SHALL show the URL
currently associated with the item, or be blank if no URL is associated.
The field SHALL be editable: the user can type or paste a URL and press
Enter to save it. When the field contains a valid URL (starts with
`http://` or `https://`), an "Open" button SHALL appear beside the field.
Clicking "Open" SHALL open the URL in a new browser tab.

**Acceptance:** Select FRIDGE. Link field is blank. Paste
`https://example.com/fridge`. Press Enter. URL is saved. "Open" button
appears. Click "Open". A new browser tab navigates to that URL. Select
FRIDGE again after reload — the URL persists. Clear the field and press
Enter — URL is removed, "Open" button disappears.

**Tested:** `test_zapp_style.py::TestStyleAPI::test_put_product_url`,
`test_zapp_style.py::TestProductUrl`.

#### SEL-13  Delete Button
When an element has a database record (i.e., is stored in the `elements`
table, not purely engine-computed), the Properties panel SHALL include a
Delete button. Clicking the button SHALL prompt for confirmation and then
remove the element via `DELETE /api/elements/<id>`.

**Acceptance:** Select a user-placed furniture item. The Properties panel
shows a Delete button. Click Delete. Confirmation dialog appears. Confirm.
The element is removed from the canvas and the database. Select an
engine-computed item (e.g., standard-variant `bed` with no DB record).
No Delete button appears.
**Tested:** Implemented in `app.js` properties panel rendering.

#### SEL-14  Dimension and Label Selectability
Dimension elements and label elements SHALL be selectable via click on
the interactive canvas. Clicking a dimension line or label SHALL select
it, apply the `selected-highlight` visual indicator, and populate the
Properties panel with its properties (name, type, anchors, style, text).

**Acceptance:** Click on a dimension line in `#layer-dims`. The dimension
receives `selected-highlight`. Properties panel shows dimension name,
start/end anchors, distance, and dim_style. Click on a room label in
`#layer-rooms`. The label receives `selected-highlight`. Properties
panel shows label name, text, position offset, rotation, and font size.
**Tested:** Implemented in `app.js` click handler and properties panel.

#### SEL-15  Constant Dependency Highlighting **(NEW)**
When a constant is focused in the Properties panel constant list (e.g.,
selecting `COUNTER_GAP` while viewing the counter appliance properties),
all geometry elements whose position or size depends on that constant
SHALL be highlighted on the canvas. The first-order dependent element
(the element directly controlled by that constant) SHALL be highlighted
in white. Second-order and downstream dependencies (elements whose
position is indirectly affected because they are anchored to the
first-order element) SHALL be highlighted in pink.

**Acceptance:** Select the counter appliance. In the Properties panel,
click the `COUNTER_GAP` constant row. The counter polygon highlights
white on the canvas. Adjacent elements whose positions shift when
`COUNTER_GAP` changes (e.g., the work counter, microwave) highlight
pink. Clicking away or selecting a different constant clears the
highlights.

---

## 6  Tools

### 6.1  Tool Palette

#### TL-1  Tool Buttons
The left palette SHALL display buttons for Select, Pan, Measure, Move,
Dimension, Draw Wall, and Add Opening tools. Exactly one tool SHALL be
active at any time.

**Acceptance:** Page loads with Select active. Click Pan -- only Pan is
highlighted.

#### TL-2  Tool Keyboard Shortcuts
The keyboard shortcuts V (Select), H (Pan), M (Measure), G (Move), D
(Dimension), W (Draw Wall), and F (Fit) SHALL activate the corresponding
tool or action when no text input is focused.

**Acceptance:** Press V -- Select tool is active. Press H -- Pan tool is
active. Press F -- viewport fits to window.

### 6.2  Select Tool

Covered by SEL-1 through SEL-4. When the Select tool is active,
clicking selects elements and the Properties panel populates.

### 6.3  Pan Tool

Covered by CV-16. When the Pan tool is active, left-click drag pans the
canvas.

### 6.4  Measure Tool

#### TL-3  Measure Tool
When the Measure tool is active, clicking and dragging SHALL draw a
dashed red line between two points. A label SHALL show the distance in
feet and inches. The status bar SHALL display the distance on mouse-up.

**Acceptance:** Select Measure tool. Click-drag between two points on the
building. A red dashed line and distance label appear. The `#measure-info`
element shows the distance.

#### TL-4  Escape Clears Measure
Pressing Escape while a measurement is displayed SHALL clear the
measurement line and info.

**Acceptance:** Measure between two points. Press Escape. Line and info
text are gone.

### 6.5  Move Tool

#### TL-5  Move by Drag
When the Move tool is active and an element is selected, clicking and
dragging SHALL reposition the element with a live preview. Releasing the
mouse SHALL commit the new position.

**Acceptance:** Select Move tool. Click on the BED furniture item. Drag
it 6 inches east. Release. Geometry recomputes with the bed at the new
position. Canvas re-renders.

#### TL-6  Move by Offset Dialog
When the Move tool is active, pressing Enter or double-clicking SHALL
open an offset dialog where the user can type a precise offset (e.g.,
"6in east", "2ft north").

**Acceptance:** Select Move tool. Select IW1. Press Enter. Dialog appears
with dx/dy input fields. Enter "6" in the East field. Click OK. IW1
moves 6 inches east.

#### TL-7  Constrained Movement
While dragging with the Move tool, holding Shift SHALL constrain movement
to the horizontal or vertical axis.

**Acceptance:** Select Move tool. Shift-drag an element. Movement is
restricted to a single axis (whichever had the larger initial delta).

#### TL-8  Snap to Grid
While dragging with the Move tool, the element SHALL snap to 1-inch grid
positions when the Grid toggle is enabled.

**Acceptance:** Enable Grid. Move tool. Drag an element. It snaps to the
nearest inch boundary.

#### TL-9  Group Move
When multiple elements are selected (via SEL-4), the Move tool SHALL move
all selected elements by the same offset, preserving their relative
positions.

**Acceptance:** Multi-select IW3, IW7, IW9. Move tool. Drag 2 inches
west. All three walls move 2 inches west. Their relative spacing is
unchanged.

#### TL-10  Move Opening Along Wall
When an opening is selected with the Move tool, dragging SHALL slide the
opening along its host wall segment (not in arbitrary directions).

**Acceptance:** Select O5. Move tool. Drag. O5 slides along the F8-F9
wall segment. Its perpendicular distance from the wall does not change.

### 6.6  Dimension Line Tool

#### TL-11  Place Dimension Line
When the Dimension tool is active, clicking two reference points or
element faces SHALL create a persistent dimension line with extension
lines, arrowheads, and a distance label.

**Acceptance:** Select Dimension tool. Click on IW2 east face. Click on
W14. A dimension line appears showing the distance in feet-inches format.

#### TL-12  Move Dimension Line
With the Select tool, clicking and dragging a dimension line SHALL
reposition it (offset perpendicular to the measurement axis).

**Acceptance:** Select a dimension line. Drag it 2 inches east. The
dimension line moves; its measured value does not change.

#### TL-13  Rotate Dimension Label
Right-clicking a dimension line SHALL offer a context menu with rotation
options (horizontal, vertical, parallel to measurement, perpendicular).

**Acceptance:** Right-click a dimension line. Select "Rotate 90 CW". The
label rotates.

#### TL-14  Delete Dimension Line
Selecting a dimension line and pressing Delete SHALL remove it.

**Acceptance:** Select a dimension line. Press Delete. The dimension line
is removed from the canvas.

### 6.7  Draw Wall Tool

#### TL-15  Draw Interior Wall
When the Draw Wall tool is active, clicking a start point and end point
SHALL create a new interior wall with a default thickness (4 inches).

**Acceptance:** `W` key activates draw-wall tool. First click sets start,
second click creates wall via `createDrawnWall()`. Element stored with
`source: "drawn"`, `start/end/thickness/poly` properties. Auto-named CW1, CW2...
**Tested:** test_zapp_tools.py::TestTL15DrawWall (4 tests).

#### TL-16  Wall Thickness Input
After placing a wall, the Properties panel SHALL show a thickness input
that can be changed before or after placement.

**Acceptance:** Select drawn wall. Properties panel shows editable Thickness
input. Changing value calls `wallPoly()` to recompute polygon, PUTs updated
properties to API.
**Tested:** test_zapp_tools.py::TestTL15DrawWall::test_update_drawn_wall_thickness.

#### TL-17  Wall Endpoint Editing
Selecting an existing wall with the Select tool and dragging an endpoint
handle SHALL extend or shorten the wall.

**Acceptance:** Infrastructure exists (start/end stored in properties).
Interactive drag handles deferred to future phase.

### 6.8  Add Element Tools

#### TL-18  Add Furniture from Catalog
Add > Furniture SHALL open a catalog dialog listing standard items with
dimensions. Clicking an item enters placement mode where the next
canvas click positions the item.

**Acceptance:** Click "Furniture" in Add palette section. Catalog grid shows
8 items (bed, dresser, shelves, loveseat, desk, chair, sofa, rocker). Click
item, then click canvas. Element created via API with `source: "placed"`.
**Tested:** test_zapp_tools.py::TestTL18PlaceFurniture (4 tests).

#### TL-19  Add Appliance from Catalog
Add > Appliance SHALL open a catalog listing standard appliances.

**Acceptance:** Catalog shows 6 appliances (washer, dryer, stove, dishwasher,
ice maker, kitchen sink). Click to place creates appliance element.
**Tested:** test_zapp_tools.py::TestTL18PlaceFurniture::test_create_placed_appliance.

#### TL-20  Add Fixture from Catalog
Add > Fixture SHALL list fixtures (toilet, bathtub, vanity).

**Acceptance:** Catalog shows 3 fixtures. Click to place creates fixture element.
**Tested:** test_zapp_tools.py::TestTL18PlaceFurniture::test_create_placed_fixture.

#### TL-21  Add Opening
Add Opening SHALL enter opening placement mode on a wall segment.

**Acceptance:** Deferred to future phase. Opening creation via API already
works (Phase 3 element CRUD); wall-click placement needs segment detection.

### 6.9  Delete Tool

#### TL-22  Delete Selected Element
Pressing Delete or Backspace with an element selected SHALL remove the
element after a confirmation prompt.

**Acceptance:** Select element. Press Delete/Backspace. `deleteSelectedElements()`
shows `confirm()` listing targets. On confirm, calls `DELETE /api/elements/<id>`
for each. Supports single and multi-selection.
**Tested:** test_zapp_tools.py::TestTL22Delete (3 tests).

#### TL-23  Cascading Delete
Deleting a wall SHALL also delete any openings or doors hosted on that
wall.

**Acceptance:** Delete IW1. Server-side cascade removes hosted openings and
doors. Client-side `IW_HOSTED_OPENINGS` map provides cascade warning in
confirmation dialog.
**Tested:** test_zapp_tools.py::TestTL23CascadeDelete (1 test).

### 6.10  Rotate Tool

#### TL-24  Rotate Element
When a placed/drawn element is selected, pressing R SHALL open a
rotation dialog with angle input and preset buttons (0/90/180/270).

**Acceptance:** Select placed element. Press R. `showRotationDialog()` opens
Dialog with angle input + presets. Submit updates element properties with
new rotation and recomputed polygon via `rotatedRectPoly()`.
**Tested:** test_zapp_tools.py::TestTL24Rotate (2 tests).

### 6.11  Shape Editor

#### TL-25  Shape Editor Dialog
The application SHALL provide a shape editor dialog for creating and
modifying item shapes stored in the `shapes` database table.
The editor SHALL support adding/moving polygon vertices and previewing
the resulting shape in real time.

**Acceptance:** `showShapeEditor()` opens modal with inline SVG canvas
showing polygon with draggable vertex handles. Vertex coordinate list
updates in real time. Add/remove vertex buttons. Shape saved via
`POST/PUT /api/shapes`.
**Tested:** test_zapp_tools.py::TestTL25ShapeEditor (5 tests).

#### TL-26  Shape Assignment
The Properties panel for placed elements SHALL include a Shape dropdown
listing all shapes from the database. Selecting a shape transforms its
polygon to the element's position/rotation.

**Acceptance:** `addShapePicker()` renders `<select>` with rect + DB shapes.
Change triggers PUT with transformed polygon.
**Tested:** test_zapp_tools.py::TestTL25ShapeEditor::test_shape_assignment_via_element.

#### TL-27  Shape Import from SVG
The shape editor SHALL support importing a polygon outline from an SVG
`<polygon>` or `<path>` element via paste.

**Acceptance:** "Import SVG" button opens textarea dialog. `parseSvgPolygon()`
extracts points from `<polygon points="...">` or `<path d="M...L...">`.
Parsed vertices populate the shape editor canvas.

---

## 7  Constants Table

#### CT-1  Table Population
The Constants panel SHALL display a table of all constants with columns:
Name, Value, Unit, Category.

**Acceptance:** Switch to the Constants tab. The table contains >= 140
rows.

#### CT-2  Value Formatting
Constants with unit `ft` SHALL display their value in inches with a `"`
suffix. Other units SHALL display 6 decimal places.

**Acceptance:** A constant like `BED_WIDTH = 76/12 ft` displays as `76"`.

#### CT-3  Category Filtering
Selecting a category from the dropdown SHALL show only constants in that
category.

**Acceptance:** Select "wall". Only constants with category "wall" are
visible.

#### CT-4  Search Filtering
Typing in the search box SHALL filter constants by name or description
(case-insensitive).

**Acceptance:** Type "bed". Only constants whose name or description
contains "bed" are shown.

#### CT-5  Column Sorting
Clicking a sortable column header SHALL sort the table by that column.
Clicking again SHALL reverse the sort direction.

**Acceptance:** Click "Name" header. Constants are sorted A-Z. Click
again. Sorted Z-A.

#### CT-6  Inline Editing
Each value cell SHALL contain an editable input. Changing the value and
pressing Enter SHALL send a PUT request to the server.

**Acceptance:** Change BED_WIDTH to `80"` and press Enter. A PUT is sent.
Toast confirms the update.

#### CT-7  Unit-Aware Value Parsing

The input parser SHALL accept numeric values with optional unit suffixes
and convert them to feet for storage. The parser is used by both the
constants table and any other dimension input field in the application.

##### CT-7a  Raw Number
A number with no unit suffix SHALL be interpreted as feet.

**Acceptance:** Enter `6.5`. The value sent to the server is `6.5`.

##### CT-7b  Feet Suffixes
A number followed by `'`, `ft`, or `feet` SHALL be interpreted as feet.
Whitespace between the number and suffix is optional.

**Acceptance:** Enter `6.5'`. Value sent is `6.5`. Enter `6.5 ft`.
Value sent is `6.5`. Enter `6.5feet`. Value sent is `6.5`.

##### CT-7c  Inch Suffixes
A number followed by `"`, `in`, or `inches` SHALL be interpreted as
inches and converted to feet by dividing by 12.

**Acceptance:** Enter `80"`. Value sent is `80 / 12 ≈ 6.6667`.
Enter `6.5 in`. Value sent is `6.5 / 12 ≈ 0.5417`.
Enter `78inches`. Value sent is `78 / 12 = 6.5`.

##### CT-7d  Centimetre Suffixes
A number followed by `cm` or `centimeters` SHALL be interpreted as
centimetres and converted to feet (÷ 30.48).

**Acceptance:** Enter `30.48cm`. Value sent is `1.0`.

##### CT-7e  Millimetre Suffixes
A number followed by `mm` or `millimeters` SHALL be interpreted as
millimetres and converted to feet (÷ 304.8).

**Acceptance:** Enter `304.8mm`. Value sent is `1.0`.

##### CT-7f  Metre Suffixes
A number followed by `m` or `meters` SHALL be interpreted as metres
and converted to feet (÷ 0.3048).

**Acceptance:** Enter `0.3048m`. Value sent is `1.0`.

##### CT-7g  Two Bare Numbers as Feet-Inches
When two numbers appear without recognised unit suffixes, the first
SHALL be interpreted as feet and the second as inches.

**Acceptance:** Enter `6 6`. Value sent is `6 + 6/12 = 6.5`.
Enter `5 0`. Value sent is `5.0`.

##### CT-7h  Multi-Token Summation
Up to three number-unit tokens may appear on a single input line.
The final value SHALL be the sum of all tokens converted to feet.
Tokens may appear in any order.

**Acceptance:** Enter `6' 6"`. Value sent is `6.5`.
Enter `6.7in 6' 2cm`. Value sent is
`6 + 6.7/12 + 2/30.48 ≈ 6.6242`.

##### CT-7i  Fourth Token Ignored
The fourth and subsequent number-unit tokens on a single input line
SHALL be silently ignored.

**Acceptance:** Enter `1' 2" 3cm 999ft`. Value sent is
`1 + 2/12 + 3/30.48 ≈ 1.2651` (the `999ft` is ignored).

##### CT-7j  Whitespace Flexibility
Whitespace between a number and its unit suffix SHALL be optional.
Whitespace SHALL separate consecutive tokens.

**Acceptance:** Enter `6'6"`. Value sent is `6.5`.
Enter `80 "`. Value sent is `80/12 ≈ 6.6667`.

#### CT-8  Value Parsing -- Fractions
Entering a value like `1/3` SHALL be parsed as `0.333...`.

**Acceptance:** Enter `1/3`. The value sent is approximately `0.3333`.

#### CT-9  Category Colour Coding
Each category cell SHALL be colour-coded (e.g. wall = blue, furniture =
green, opening = yellow).

**Acceptance:** Inspect category cells. Different categories have visually
distinct colours.

#### CT-10  Changed Value Indicator
When an input value differs from its original loaded value, the input
SHALL receive a visual indicator (yellow border).

**Acceptance:** Change a value without pressing Enter. The input border
turns yellow.

---

## 8  Data Tables

### 8.1  Outline Chain Table

#### DT-1  Outline Chain Display
The Outline panel SHALL display a table with columns: #, Type, Dist/R,
Sweep, End. It SHALL contain 18 rows.

**Acceptance:** Switch to the Outline tab. Table has 18 rows. Line
segments show distance; arc segments show radius and sweep.

#### DT-2  Outline Chain Editing
Clicking a cell in the Dist/R, Sweep, or Bearing columns SHALL make it
editable. Pressing Enter SHALL commit the change, trigger closure
re-solve, and recompute geometry.

**Acceptance:** Click the Dist/R cell for segment 5 (an arc). Change the
value from 28 to 30. Press Enter. Engine re-solves closure. Canvas
re-renders with the updated outline.

#### DT-3  Outline Closure Indicator
The Outline panel SHALL display a closure status indicator showing whether
the current chain parameters produce a valid closed outline. Red = open,
green = closed.

**Acceptance:** After a valid edit, indicator is green. After setting an
impossible radius, indicator turns red with an error distance.

#### DT-4  Add/Remove Chain Segment
The Outline panel SHALL have "+" and "-" buttons to insert a new F-point
(splitting a segment) or remove an existing F-point.

**Acceptance:** Click "+" on segment 11-12. A dialog prompts for the new
point name (e.g., F11a) and segment type. Confirm. A new row appears.
Chain re-solves.

### 8.2  Openings Tables

#### DT-5  Outer Openings Table
The Openings panel SHALL display outer openings (O1-O11, O8a) with
columns: Name, Segment, Width (inches).

**Acceptance:** Panel shows 12 rows with correct opening names.

#### DT-6  Rough Openings Table
The Openings panel SHALL display rough openings (RO1-RO7) with columns:
Name, Wall, Width (inches), Orientation.

**Acceptance:** Panel shows 7 rows. Each row has a wall name (e.g. IW1)
and orientation (H, V, or R).

#### DT-7  Openings Table Editing
Width cells in the openings tables SHALL be editable. Changing a width
SHALL update the opening constant and trigger geometry recomputation.

**Acceptance:** `updateOpeningsTable()` creates `<input>` elements for width
cells where `findWidthConstant()` maps the opening name to a constant.
Change triggers `handleConstantEdit()`. Both outer and rough openings tables
support inline width editing.
**Tested:** test_zapp_tools.py::TestDT7OpeningWidth (8 tests).

#### DT-8  Table Row Selection
Clicking a row in the openings tables SHALL select the corresponding
element on the canvas and show its properties.

**Acceptance:** Click a row in the outer openings table. The
corresponding polygon on the canvas receives the `selected-highlight`
class.

### 8.3  Interior Walls Table

#### DT-9  Interior Walls Table
The Elements panel SHALL display an interior walls table with columns:
Name, Thickness (inches), Length (ft-in), Orientation, Hosted Openings.

**Acceptance:** Table shows 13 rows matching ENG-4 wall names. Thickness
values are in inches. Length is formatted as ft'in".

#### DT-10  Interior Walls Table Editing
Clicking a wall row SHALL select it on the canvas. Thickness and length
cells SHALL be editable.

**Acceptance:** Click IW8 row. IW8 highlights on canvas. Change thickness
from 4 to 6. Press Enter. Geometry recomputes.

### 8.4  Furniture/Appliance Table

#### DT-11  Furniture/Appliance Table
The Elements panel SHALL display a furniture/appliance table with columns:
Name, Type, Width, Depth, Rotation, URL.

**Acceptance:** Table shows rows for all furniture and appliance items.
Items with product URLs show clickable links.

---

## 9  Element Operations

### 9.1  Outline Editing

#### OE-1  Drag F-Points
When the Select tool is active and the "Points" toggle is on, clicking
and dragging an F-series point SHALL reshape the building outline. The
engine SHALL re-solve the outline chain parameters to match the dragged
position while maintaining tangency and closure constraints.

**Acceptance:** Drag F10 2 inches west. Engine derives new segment
parameters. Outline re-renders. Tangency holds at all arc junctions.

#### OE-2  Arc Radius Handle
When an arc segment is selected, a radius adjustment handle SHALL appear.
Dragging the handle SHALL change the arc radius with live preview.

**Acceptance:** Select the F13-F14 arc. A handle appears on the arc
midpoint. Drag inward. Radius decreases. Release. Engine re-solves.

#### OE-3  Constraint-Based Outline Editing
Edit > Set Constraint SHALL allow the user to specify a target distance
between two named points (e.g., "F6-F7 = 5'3""). The engine SHALL solve
for chain parameters that satisfy the constraint while maintaining
closure.

**Acceptance:** Set constraint "W20-W1 to W9-W10 clear span = 23'6"".
Engine adjusts arc radii/lengths to satisfy. Result is verified by
measuring the span in the computed geometry.

### 9.2  Door Configuration

#### DOOR-1  Add Door to Opening
Right-clicking an opening SHALL offer "Add Door". A dialog SHALL prompt
for door width, hinge side, swing direction, and type (single/double).

**Acceptance:** Right-click RO1. Select "Add Door". Dialog appears. Set
hinge = east, swing = south, type = single. Confirm. Door arc renders on
the canvas.
**Tested:** `test_zapp_canvas.py::TestDoorArcAPI::test_door_arc_after_door_create`.

#### DOOR-2  Edit Door Hinge and Swing
Selecting a door arc and clicking a "Flip Hinge" button SHALL move the
hinge to the opposite side. "Flip Swing" SHALL reverse the swing
direction.

**Acceptance:** Select RO3 door. Click "Flip Hinge". Door hinge moves
from south to north. Door arc re-renders.
**Tested:** `test_zapp_canvas.py::TestDoorArcAPI::test_door_invalidation_via_api`.

#### DOOR-3  Appliance Door Arcs
Appliances (fridge, washer, dryer, microwave) SHALL render door arcs
with hinge corner and swing direction defined by variant item metadata.

**Acceptance:** Fridge, washer, dryer, and microwave display door arcs.
Stacked appliance doors (microwave) render above their parent counter.
**Tested:** `test_zapp_canvas.py::TestApplianceDoors` (7 tests).

#### DOOR-4  Double Door
Openings SHALL support double-door configuration with two leaves, each
with independent hinge side.

**Acceptance:** Configure RO7 as a double door. Two door arcs render,
one hinged on each side of the opening.
**Tested:** `test_zapp_canvas.py::TestDoorArcs::test_door_arc_ro7_double`,
`test_door_arc_ro6_double`.

### 9.3  Hyperlinks

#### LINK-1  Attach Product URL
Each furniture, appliance, fixture, and opening element SHALL have an
optional product URL property.

**Acceptance:** Select SHELVES. In Properties panel, enter an IKEA URL.
Save. The URL persists across page reloads.

**Tested:** `test_zapp_style.py::TestProductUrl`,
`test_zapp_style.py::TestStyleAPI::test_put_product_url`.

#### LINK-2  Clickable SVG Links
Elements with product URLs SHALL be wrapped in `<a xlink:href="...">` tags
in generated SVGs, making them clickable in browsers and PDF viewers.

**Acceptance:** Assign a URL to FRIDGE. Regenerate floorplan SVG. Open the
SVG in a browser. Clicking FRIDGE opens the product URL.

**Tested:** Interactive canvas wraps elements with `<a xlink:href>` tags.
`test_zapp_style.py::TestStyleAPI::test_put_product_url`.

### 9.4  Room Labels & Annotations

Phase 8 unifies room labels under the `elements` table (type `'label'`),
superseding the `room_label_offsets` table from Phase 0.  Each room label
is stored as an element whose default position is the area-weighted
centroid of its room boundary polygon.  The database stores an offset
`(de, dn)` from that centroid and a rotation angle (defaults: `(0, 0)`
offset, `0°` rotation).  "Reset to Defaults" back-computes label
positions from the existing SVG generation scripts' output.  The 11
auto-computed room labels become editable label elements; users can also
add additional custom labels.

#### LABEL-1  Add Room Label
Edit > Add Room Label SHALL allow the user to place a room name label by
clicking a position on the canvas and typing the label text (e.g.,
"BEDROOM", "OFFICE").  The label SHALL be stored as an element in the
`elements` table with type `'label'`.

**Acceptance:** Add Room Label. Click in the bedroom area. Type "BEDROOM".
Press Enter. Label appears on the canvas and persists in the database
(`SELECT * FROM elements WHERE type = 'label' AND name = 'BEDROOM'`
returns a row).

#### LABEL-2  Move Label
With the Select tool, labels SHALL be draggable to reposition them.
Dragging SHALL update the stored offset from the room centroid.

**Acceptance:** Select the "OFFICE" room label. Drag it 1 foot north.
Label position updates. The stored offset reflects the new position.
Generated SVGs reflect the new position.

#### LABEL-3  Edit Label Text
Double-clicking a label SHALL open an inline text editor to change the
label text.

**Acceptance:** Double-click "BATH" label. Text becomes editable. Change
to "BATHROOM". Press Enter. Label updates.

#### LABEL-4  Label Font Size
Each label SHALL have a configurable font size property, editable in the
Properties panel.

**Acceptance:** Select a room label. Properties panel shows Font Size
field. Change from 10pt to 14pt. Label re-renders at the new size.

---

## 10  Styling (Phase 9)

#### STYLE-1  Element Fill Colour
Each element type (wall, opening, furniture, etc.) SHALL have a
configurable fill colour, editable via the Properties panel or a colour
picker.

**Acceptance:** Select a wall. Open colour picker. Change fill from
default grey to light brown. Wall re-renders with the new fill colour.

**Tested:** `test_zapp_style.py::TestStyleDefaults` (defaults match CSS),
`test_zapp_style.py::TestStyleAPI::test_put_style_properties`.

#### STYLE-2  Element Stroke Properties
Each element SHALL have configurable stroke colour, width, and style
(solid, dashed, dotted).

**Acceptance:** Select a dimension line. Change stroke style to dashed.
Line re-renders as dashed.

**Tested:** `test_zapp_style.py::TestStyleValidation` (stroke_style,
stroke_width validation),
`test_zapp_style.py::TestStyleResolution::test_resolve_base_override`.

#### STYLE-3  Element Opacity
Each element SHALL have a configurable opacity (0-100%), editable in the
Properties panel.

**Acceptance:** Select interior dim lines. Set opacity to 50%. Lines
render at half opacity.

**Tested:** `test_zapp_style.py::TestStyleValidation` (opacity range),
`test_zapp_style.py::TestStyleResolution::test_resolve_partial_override`.

#### STYLE-4  Per-View Styling
View variants SHALL support per-view element styling (e.g., doors at 20%
opacity in the plumbing view).

**Acceptance:** In plumbing view variant, set door opacity to 20%. Doors
render faintly. Main floorplan view still shows doors at 100%.

**Tested:** `test_zapp_style.py::TestStyleResolution::test_resolve_view_override`,
`test_zapp_style.py::TestPerViewStyleResolution`,
`test_zapp_style.py::TestStyleAPI::test_put_view_overrides`.

---

## 11  Site Plan **(NEW)**

#### SITE-1  Structure Placement
The site plan view SHALL allow positioning the building outline on the
survey parcel by specifying setback distances from property boundaries.

**Acceptance:** Open site plan view. Enter north setback = 11.5 feet.
Building repositions. Setback distance labels update.

#### SITE-2  Drainfield Operations
The site plan SHALL support adding, positioning, and sizing drainfield
rectangles.  The drainfield SHALL be stored as a single element in the
`elements` table (type `'site_element'`), shared across both the site
plan and plumbing views — moving it in one view updates the other.

**Acceptance:** Tools > Add Drainfield. Click on the site plan. Enter
dimensions 25x10 feet. Drainfield rectangle appears. Drag to reposition.

#### SITE-3  Site Annotations
The site plan SHALL support text annotations, distance labels, and
directional arrows (e.g., "FRONT" arrow).

**Acceptance:** Add a text annotation "PROPOSED 950SF MAX ADU". Position
it on the site plan. Add a "FRONT" arrow pointing north.

#### SITE-4  Parcel Corner Markers
The site plan SHALL render survey point markers (P-series) with distance
labels between corners.

**Acceptance:** Site plan shows P2, P3, P4, P5, POB markers with
connecting boundary lines and distance labels.

---

## 12  3D Model **(NEW)**

#### SCAD-1  SCAD Generation ✅
Tools > Generate 3D Model SHALL run the OpenSCAD generator and produce a
3D model file.

**Acceptance:** Click Generate 3D Model. SCAD file is generated. Success
toast appears.

*Implemented: `POST /api/generate-3d` runs selected roof script via
`generate_svg()`. Tests: `test_zapp_scad.py::TestGenerateScad`.*

#### SCAD-2  Roof Style Selection ✅
The 3D model settings SHALL allow selecting roof style (flat, 2:12 slope)
and overhang distance.

**Acceptance:** Select "2:12 slope" in 3D settings. Regenerate. SCAD file
contains sloped roof geometry.

*Implemented: `config` table stores `roof_style`. File menu has roof style
dropdown. Overhang editable via `ROOF_OVERHANG` constant.
Tests: `test_zapp_scad.py::TestConfig`.*

#### SCAD-3  Multi-View Layout ✅
Tools > Generate Views SHALL produce a multi-view PDF (3views.pdf) with
floor plan and elevation views.

**Acceptance:** Click Generate Views. PDF is created with 4 panels
(floorplan + 3 elevations).

*Implemented: `POST /api/generate-views` runs SCAD generator, view renderer,
line drawings, and 3views composer. Tests: `test_zapp_scad.py::TestGenerateViews`.*

---

## 13  Analysis

#### ANALYSIS-1  Span Analysis View
Tools > Span Analysis SHALL display the N-S interior span graph showing
span distance vs. position.

**Acceptance:** Open Span Analysis. Graph shows span curve. Mouse hover
displays span value at each position.

**Tested:** `test_zapp_analysis.py::TestSpanData` (API returns valid
arrays, monotonic eastings, positive spans).

#### ANALYSIS-2  Span vs. Rotation
Tools > Span vs. Rotation SHALL display span measurements across rotation
angles (5-175 degrees).

**Acceptance:** Graph shows span vs. rotation with min/max markers.
Properties panel shows min/max span and rotation angles.

**Tested:** `test_zapp_analysis.py::TestSpanRotation` (API returns
min/max angles and spans, data array with 35 angle samples).

#### ANALYSIS-3  Room Area Display
The "Areas" display toggle (DIS-12) SHALL compute and display room areas
in square feet for each enclosed region. On the SF variant, areas are
always shown.

**Acceptance:** Room areas appear within room boundaries. Total area sums
correctly. Values update when walls are moved.

**Tested:** `test_zapp_analysis.py::TestRoomAreas` (areas present for
all variants, SF has polys, total area in range).

---

## 14  Plumbing Layout **(NEW)**

The plumbing view SHALL be a full interactive layout — like the floorplan
layouts (standard, minik, daybed, bare, sf) — with plumbing-specific editing
tools, database-stored element configuration, and CRUD APIs.  Plumbing
elements include pipes (supply and drain), fittings (T-stubs, elbows, valves),
and fixture connections.

#### PLUMB-1  Plumbing Interactive Canvas
The plumbing layout SHALL render the building outline with all plumbing
elements on an interactive canvas (not a static SVG), supporting the same
zoom, pan, and selection interactions as the floorplan canvas.

**Acceptance:** Switch to Plumbing tab. Canvas shows the building outline
with fixtures (toilets, sinks, washer), supply lines, and drain lines.
Pan and zoom work. Clicking a pipe or fixture selects it and shows
properties in the right panel.

#### PLUMB-2  Supply Line Routing **(NEW)**
The plumbing plan editor SHALL allow drawing supply line paths from the
water source to fixtures, with hot (red) and cold (blue) line
differentiation.

**Acceptance:** Select supply line tool. Draw a path from Well to sink.
Line renders in blue. Switch to hot line. Draw parallel path. Renders red.

#### PLUMB-3  Plumbing Fixtures Table **(NEW)**
The plumbing view SHALL display a fixtures/supplies table listing each
fixture with its supply and drain connections.

**Acceptance:** Table shows rows for each plumbing fixture with type,
supply pipe size, and drain pipe size.

#### PLUMB-4  Plumbing Elements Database **(NEW)**
Plumbing elements (pipes, fittings, fixture connections) SHALL be stored
in a `plumbing_elements` table with columns for element type, geometry
(path coordinates), properties (pipe size, material, hot/cold), and
fixture associations.

**Acceptance:** Create a supply line via the drawing tool. Query
`SELECT * FROM plumbing_elements WHERE type = 'supply_pipe'` — returns the
created pipe with path coordinates and properties.  Close and reopen the
app — the pipe persists.

#### PLUMB-5  Plumbing CRUD API **(NEW)**
The API SHALL provide CRUD endpoints for plumbing elements:
`GET /api/plumbing`, `POST /api/plumbing`, `PUT /api/plumbing/<id>`,
`DELETE /api/plumbing/<id>`.

**Acceptance:** `POST /api/plumbing` with pipe JSON returns 201. `GET`
returns the created element.  `PUT` updates properties.  `DELETE` removes
it.  Each mutation broadcasts an SSE event.

#### PLUMB-6  Drain Line Routing **(NEW)**
The plumbing editor SHALL allow drawing drain line paths from fixtures to
the waste outlet, with slope annotations.

**Acceptance:** Select drain line tool. Draw a path from toilet to waste
outlet.  Line renders in the drain colour.  Slope annotation displays
grade (e.g., "1/4 in/ft").

#### PLUMB-7  Fixture Placement Tool **(NEW)**
The plumbing editor SHALL provide a fixture placement tool for adding
plumbing fixtures (toilets, sinks, washer, water heater) with automatic
supply and drain stub connections.

**Acceptance:** Select fixture tool, choose "sink" from catalog.  Click
to place on canvas.  Fixture renders with supply and drain stubs.  Stubs
are selectable and connectable to supply/drain lines.

#### PLUMB-8  Pipe Fitting Editing **(NEW)**
The plumbing editor SHALL support placing and editing pipe fittings
(T-stubs, elbows, valves) at pipe junctions and along pipe runs.

**Acceptance:** Select fitting tool, choose "T-stub."  Click on a pipe
intersection.  Fitting renders at the junction.  Properties panel shows
fitting type and size.  Changing fitting type re-renders.

---

## 15  Undo/Redo

#### UNDO-1  Undo Last Action
Edit > Undo or Ctrl+Z SHALL revert the most recent edit operation
(constant change, element move, add, delete).

**Acceptance:** Move IW1 6 inches east. Press Ctrl+Z. IW1 returns to its
original position. Geometry recomputes.

#### UNDO-2  Redo Last Undo
Edit > Redo or Ctrl+Shift+Z SHALL re-apply the most recently undone
operation.

**Acceptance:** Undo the IW1 move. Press Ctrl+Shift+Z. IW1 moves back to
the edited position.

#### UNDO-3  Undo History Depth
The undo history SHALL support at least 50 levels of undo.

**Acceptance:** Perform 50 edits. Undo all 50. Each undo correctly
restores the previous state.

#### UNDO-4  Undo Across All Mutation Operations
Undo SHALL work for all mutation operations introduced in any phase:
constant edits, element moves, adds, deletes, door changes, outline edits,
label changes, styling changes, site plan edits, and variant operations.
Each phase that introduces new mutation types SHALL extend undo coverage
to include those mutations.  This is enforced by the Phase Completion
Protocol (see ROADMAP.md), which requires undo verification for all new
mutation types as a gate for phase completion from Phase 3 onward.

**Acceptance:** Add a door, move a wall, change a constant. Undo three
times. All three operations are correctly reverted. As new mutation types
are added in later phases, each SHALL be verified undoable.

---

## 16  Real-Time Updates

#### RT-1  SSE Connection
On page load the browser SHALL establish a Server-Sent Events connection
to `/api/events` and display a green "Connected" indicator.

**Acceptance:** Load the page. The connection status shows "Connected"
with a green background.

#### RT-2  Geometry Change Notification
When a constant is updated via the API, the server SHALL broadcast a
`geometry_changed` event. All connected browsers SHALL reload geometry
and re-render the canvas.

**Acceptance:** Open two browser tabs. In tab A, edit a constant. Tab B
automatically re-renders with the new geometry.

#### RT-3  SVG Update Notification
When a view is regenerated via the API, the server SHALL broadcast an
`svg_updated` event with the view name. If a browser is viewing that SVG,
it SHALL reload automatically.

**Acceptance:** Open the app viewing the Floorplan SVG tab. Trigger
regeneration via the File menu. A toast appears and the SVG reloads.

#### RT-4  Reconnection
If the SSE connection drops, the browser SHALL attempt to reconnect after
3 seconds and restore the "Connected" indicator on success.

**Acceptance:** Stop and restart the server. Within a few seconds the
browser reconnects and the status returns to "Connected".

#### RT-5  Element Change Notification
When an element is added, moved, deleted, or modified via the API, the
server SHALL broadcast an `element_changed` event so all connected
browsers update in real time.

**Acceptance:** Open two tabs. In tab A, add a furniture item. Tab B
shows the new item without manual refresh.

---

## 17  Application

### 17.1  Launcher

#### APP-1  Default Start
`python run_app.py` SHALL start the server on `127.0.0.1:5000` and open
the default browser.

**Acceptance:** Run the command. Server log shows "Running on
http://127.0.0.1:5000". Browser opens to that URL.

#### APP-2  Custom Port
`python run_app.py --port 8080` SHALL start the server on port 8080.

**Acceptance:** Run the command. Server log shows port 8080.

#### APP-3  No-Browser Flag
`python run_app.py --no-browser` SHALL start the server without opening
a browser.

**Acceptance:** Run the command. No browser window opens. Server is
accessible via HTTP.

#### APP-4  Database Auto-Creation
If `app/adu.db` does not exist, the server SHALL create and seed it
automatically on startup.

**Acceptance:** Delete `app/adu.db`. Start the server. Verify the file is
created with seeded data.

### 17.2  Non-Functional

#### NF-1  Dark Theme
The UI SHALL use a dark colour scheme with light text on dark backgrounds.

**Acceptance:** Visual inspection -- background is dark (#1e1e2e), text is
light.

#### NF-2  Responsive Layout
At viewport widths below 1000 px, the tool palette SHALL collapse to icon-
only mode and the right panel SHALL narrow.

**Acceptance:** Resize window below 1000 px. Tool button labels disappear.
Panel width shrinks.

#### NF-3  Existing Tests Unaffected
All pre-existing tests in `tests/` SHALL continue to pass after the app
is added.

**Acceptance:** `python -m pytest tests/ -x -q` reports all pre-existing
tests (i.e. non-`test_zapp_*` tests) passed, 0 failed.

#### NF-4  No Modification of Existing Code
The app SHALL not modify any files in `shared/`, `floorplan/`, `walls/`,
`span/`, `survey/`, `roof/`, `site/`, `scad/`, or `plumbing/`.

This constraint applies until the editor has achieved 100% functional
completeness and has been approved for cutover to database-only data
sources.  Until then, the existing generator scripts serve as reference
implementations to verify the app produces identical output.  Code
duplication between `app/` and the existing packages (e.g.,
`app/variants.py` replicating positioning math from
`floorplan/gen_floorplan.py`) is intentional during this period and will
be consolidated at cutover.

**Acceptance:** `git diff --name-only` shows changes only in `app/`,
`tests/`, `run_app.py`, `pyproject.toml`, and `.gitignore`.

#### NF-5  Keyboard Shortcut Safety
Keyboard shortcuts SHALL NOT fire when the user is typing in an `<input>`
or `<select>` element.

**Acceptance:** Focus a constant value input. Press V. The tool does not
change; the character is typed into the input instead.

#### NF-6  Default Display Units
Unless otherwise specified, all displayed dimensions SHALL use feet and
inches. Inch values SHALL be displayed to two decimal places of precision
with trailing zeroes removed (e.g., `5' 3.5"` not `5' 3.50"`; `12' 0"`
not `12' 0.00"`). All angles SHALL be displayed in degrees to four
decimal places of precision with trailing zeroes removed and a degree
symbol suffix (e.g., `45.5°` not `45.5000°`; `90°` not `90.0000°`).

**Acceptance:** Inspect dimension labels, property panel values, and
measurement tool output. All use feet-inches format. A value of exactly
3 inches displays as `3"`, not `3.00"`. A value of 7.50 inches displays
as `7.5"`. Inspect the Outline table Sweep column. Arc sweep values
display with degree symbol and trailing zeroes removed (e.g.,
`32.5921°`).

---

## Appendix A: Requirement Cross-Reference by User Operation

This table maps common user operations (identified from the full commit
history) to the requirements that enable each operation through the GUI.

| Operation | Commits | Requirements |
|-----------|---------|-------------|
| Move element (wall/furniture/opening) | ~146 | TL-5, TL-6, TL-7, TL-8, TL-9, TL-10, API-23 |
| Edit outline geometry (F-series chain) | ~79 | OE-1, OE-2, OE-3, DT-2, DT-3, DT-4, API-16..19 |
| Add/edit dimension lines & labels | ~78 | TL-11..14, CV-9, DIS-7, DIM-1..5, SEL-14, LABEL-1..4 |
| Add/resize/rotate furniture | ~67 | TL-18..20, TL-24, SEL-8, DT-11 |
| Add/edit/remove interior walls | ~52 | TL-15..17, TL-22, TL-23, DT-9, DT-10 |
| Add/edit openings | ~42 | TL-21, DT-7, SEL-10, API-24..26 |
| Configure doors (hinge/swing) | ~35 | DOOR-1..4, SEL-7, SEL-11, CV-7, API-27..29 |
| Edit site plan | ~36 | SITE-1..4 |
| 3D/SCAD model | ~31 | SCAD-1..3 |
| View variants | ~24 | UI-5, UI-6, UI-8, ENG-13, API-8, API-34, CV-4, VAR-1..3 |
| Plumbing layout | ~22 | PLUMB-1..8 |
| Span/area analysis | ~22 | ANALYSIS-1..3, DIS-12 |
| Resize elements | ~20 | SEL-10, DT-7, DT-10 |
| Element styling (colour/opacity) | ~18 | STYLE-1..4 |
| Product hyperlinks | ~14 | LINK-1, LINK-2, SEL-12, CV-12, DIS-11 |
| Delete elements | ~8 | TL-22, TL-23, SEL-13, API-22, API-26, API-29 |
| Undo/redo | implicit | UNDO-1..4, API-30, API-31 |

## Appendix B: Requirements Summary

"Implemented" = working with test coverage (Phases 0–9).
"Planned" = future phases, marked **(NEW)** on the requirement
line or inherited from a **(NEW)** section/subsection heading.

| Section | Implemented | Planned | Total |
|---------|-------------|---------|-------|
| 1 Data Layer | 40 | 0 | 40 |
| 2 REST API | 34 | 0 | 34 |
| 3 UI Layout | 6 | 2 | 8 |
| 4 Canvas | 26 | 1 | 27 |
| 5 Selection | 14 | 1 | 15 |
| 6 Tools | 27 | 0 | 27 |
| 7 Constants | 20 | 0 | 20 |
| 8 Data Tables | 11 | 0 | 11 |
| 9 Element Ops | 12 | 1 | 13 |
| 10 Styling | 4 | 0 | 4 |
| 11 Site Plan | 0 | 4 | 4 |
| 12 3D Model | 0 | 3 | 3 |
| 13 Analysis | 0 | 3 | 3 |
| 14 Plumbing | 0 | 8 | 8 |
| 15 Undo/Redo | 4 | 0 | 4 |
| 16 Real-Time | 5 | 0 | 5 |
| 17 Application | 10 | 0 | 10 |
| **Total** | **213** | **23** | **236** |

CT-7 (Unit-Aware Value Parsing) is counted as one requirement alongside
its 10 sub-requirements CT-7a through CT-7j, which are also counted
individually.  This gives 11 items from the CT-7 family: the parent
requirement plus its 10 lettered children.
