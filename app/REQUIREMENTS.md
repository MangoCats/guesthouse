# ADU Editor — Requirements

Testable requirements for the interactive building editor application.
Each requirement has a unique ID, a description, and acceptance criteria
that can be verified by automated or manual testing.

Requirements are organized into a modular tree that groups functionality
both by **implementation layer** (data, API, UI) and by **user experience**
(what the user can do). New requirements identified from analysis of the
full 864-commit history are marked **(NEW)**.

---

## 1  Data Layer

### 1.1  Database Schema & Seeding

#### DB-1  Schema Initialisation
The application SHALL create an SQLite database with tables `constants`,
`outline_chain`, `views`, `elements`, and `doors` when launched for the
first time.

**Acceptance:** Start with no `app/adu.db` file. Run `python run_app.py
--no-browser`. Verify the file is created and contains all five tables.

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

#### DB-9  Elements Table **(NEW)**
The `elements` table SHALL store interior walls, furniture, appliances,
fixtures, and clearance zones with columns: `id`, `type` (wall /
furniture / appliance / fixture / clearance), `name`, `properties` (JSON),
`variant` (nullable — for floorplan variant membership).

**Acceptance:** Table exists with the specified columns. Seeding populates
rows for all 13 interior walls, 6 furniture/appliance items, and fixtures.

#### DB-10  Doors Table **(NEW)**
The `doors` table SHALL store door configurations with columns: `id`,
`opening_name`, `width`, `hinge_side`, `swing_direction`, `door_type`
(single / double).

**Acceptance:** Table exists. Seeding populates rows for all openings that
have doors (RO1-RO7, appliance doors).

#### DB-11  Undo History Table **(NEW)**
The database SHALL maintain an `undo_history` table storing serialised
snapshots of changed state, enabling undo and redo operations.

**Acceptance:** After an edit, `SELECT count(*) FROM undo_history` >= 1.
Each row contains a `timestamp`, `action_type`, and `before_state` /
`after_state` JSON.

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
The computed geometry SHALL include exactly 13 interior walls: IW1, IW2,
IW2O, IW2S, IW3, IW4, IW5, IW6, IW7, IW8, IW9, IW11, IW12.

**Acceptance:** `set(result["interior_walls"].keys())` equals the expected
set.

#### ENG-5  Opening Counts
The computed geometry SHALL include exactly 12 outer openings and 7 rough
openings.

**Acceptance:** `len(result["outer_openings"]) == 12` and
`len(result["rough_openings"]) == 7`.

#### ENG-6  Appliance and Furniture Counts
The `appliances` dict SHALL contain keys `dryer`, `washer`, `counter`.
The `furniture` dict SHALL contain keys `bed`, `dresser`, `shelves`.

**Acceptance:** Verify both key sets.

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

#### ENG-11  Outline Chain Mutation **(NEW)**
The engine SHALL accept modified outline chain parameters (radius, sweep,
bearing, length) and re-solve for closure/tangency before recomputing
geometry.

**Acceptance:** Change F13-F14 arc radius from 28 to 30 inches. Engine
re-solves closure distances. Computed geometry reflects the new radius.
All tangency constraints hold.

#### ENG-12  Room Area Computation **(NEW)**
The engine SHALL compute room areas (in square feet) for each enclosed
region defined by interior walls and the building outline.

**Acceptance:** `result["room_areas"]` contains entries for BEDROOM,
OFFICE, BATH, UTIL, KITCHEN with positive numeric values.

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
constants and reload the Constants table.

**Acceptance:** Click Reset to Defaults. A confirmation dialog appears.
Confirm. Constants table refreshes with original values.

#### GEN-4  Interactive View Not Regenerable
Clicking Regenerate Current View while on the Interactive tab SHALL show
a toast explaining that the interactive view cannot be regenerated.

**Acceptance:** Select Interactive tab. Click Regenerate Current View.
Toast reads "Cannot regenerate interactive view".

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

### 2.3  Outline Chain Mutation API **(NEW)**

#### API-16  PUT /api/outline/<seq>
SHALL accept a JSON body with optional keys `dist_or_radius`, `sweep`,
`bearing` and update the specified outline chain segment. SHALL trigger
closure re-solve and geometry recomputation.

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

### 2.4  Element CRUD API **(NEW)**

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

### 2.5  Door API **(NEW)**

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

### 2.6  Undo/Redo API **(NEW)**

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
the Interactive tab is active. The selector SHALL list all 5 floorplan
variants (Standard, Small Kitchen, Daybed, Room Dimensions, Square
Footage). Selecting a variant SHALL reload the geometry with that
variant's furniture/appliance set. The selector SHALL be hidden when an
SVG view tab is active.

**Acceptance:** Load the page with the Interactive tab active. The
"Layout" dropdown is visible with 5 options. Select "Small Kitchen" --
the canvas re-renders with minik furniture (sofa, rocker, cooktop
visible; loveseat, stove absent). Switch to an SVG view tab -- the
dropdown is hidden. Switch back to Interactive -- the dropdown reappears
with the previously selected variant.

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
`item-appliance`, `item-furniture`, `item-fixture`.

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

#### CV-7  Door Arc Rendering **(NEW)**
When a door is configured on an opening, the canvas SHALL render the door
swing arc showing hinge position and sweep direction.

**Acceptance:** An opening with a configured door displays a dashed arc
from closed to open position. Hinge point is marked with a small circle.

#### CV-8  Room Labels **(NEW)**
When the "Room Labels" toggle is checked, the canvas SHALL render room
name labels (BEDROOM, OFFICE, BATH, UTIL, KITCHEN) centred within their
respective enclosed regions.

**Acceptance:** `#layer-room-labels` contains text elements. Each label
is positioned within its room boundary polygon.

#### CV-9  Dimension Lines **(NEW)**
When the "Dimensions" toggle is checked, the canvas SHALL render
persistent dimension annotations showing distances between reference
points/faces with extension lines, arrowheads, and value labels.

**Acceptance:** `#layer-dimensions` contains dimension line groups. Each
group has two extension lines, a measurement line, and a text label
showing the distance in feet-inches format.

#### CV-10  Area Labels **(NEW)**
When the "Areas" toggle is checked, the canvas SHALL render room area
labels (in square feet) within each enclosed region.

**Acceptance:** `#layer-areas` contains text elements showing values like
"125 sf" positioned within room boundaries.

#### CV-11  Clearance Zones **(NEW)**
The canvas SHALL render clearance circles (WW-series) as dashed circles
at fixture locations when the "Clearance" toggle is checked.

**Acceptance:** Clearance circles appear at toilet, sink, and bed
locations with the specified radii.

#### CV-12  Hyperlink Indicators **(NEW)**
Elements with attached product URLs SHALL display a small link icon
overlay on the canvas.

**Acceptance:** A furniture item with a URL shows a link icon at its
corner. Clicking the icon opens the URL in a new tab.

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

#### DIS-6  Room Labels Toggle **(NEW)**
Checking the "Room Labels" checkbox SHALL show/hide room name labels.

**Acceptance:** Toggle checkbox. Room label text elements appear/disappear.

#### DIS-7  Dimensions Toggle **(NEW)**
Checking the "Dimensions" checkbox SHALL show/hide persistent dimension
annotations.

**Acceptance:** Toggle checkbox. Dimension line groups appear/disappear.

#### DIS-8  Areas Toggle **(NEW)**
Checking the "Areas" checkbox SHALL show/hide room area (sf) labels.

**Acceptance:** Toggle checkbox. Area label text elements appear/disappear.

#### DIS-9  Clearance Toggle **(NEW)**
Checking the "Clearance" checkbox SHALL show/hide fixture clearance
circles.

**Acceptance:** Toggle checkbox. WW-series dashed circles appear/disappear.

#### DIS-10  Doors Toggle **(NEW)**
Checking the "Doors" checkbox SHALL show/hide door swing arcs.

**Acceptance:** Toggle checkbox. Door arc elements appear/disappear.

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

#### SEL-4  Multi-Select **(NEW)**
Shift-clicking or drag-selecting SHALL allow multiple elements to be
selected simultaneously for group operations (move, delete).

**Acceptance:** Shift-click three interior walls. All three are highlighted.
Properties panel shows "3 elements selected" with common properties.

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

#### SEL-7  Door Properties **(NEW)**
When an opening with a door is selected, the Properties panel SHALL
display door width, hinge side, swing direction, and door type, all
editable.

**Acceptance:** Select RO1. Panel shows Door section with Hinge Side
dropdown (east/west/north/south), Swing Direction dropdown, and Type
dropdown (single/double).

#### SEL-8  Furniture Properties **(NEW)**
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

#### SEL-10  Opening Width Editing **(NEW)**
The opening width field in the Properties panel SHALL be editable.
Changing it SHALL update the opening width constant and trigger geometry
recomputation.

**Acceptance:** Select O8. Change width from 19 to 25 inches. Press
Enter. Opening re-renders at the new width.

#### SEL-11  Door Property Editing **(NEW)**
Changing a door property (hinge side, swing direction) in the Properties
panel SHALL update the door configuration and re-render the door arc.

**Acceptance:** Select RO3. Change hinge side from "south" to "north" via
dropdown. Door arc flips to the opposite side.

#### SEL-12  Product URL Editing **(NEW)**
Each furniture, appliance, and fixture element SHALL have an editable URL
field in the Properties panel for linking to product pages.

**Acceptance:** Select SHELVES. Enter an IKEA URL in the URL field. Press
Enter. The URL is saved. The generated SVG wraps the element in an
`<a>` tag linking to that URL.

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

### 6.5  Move Tool **(NEW)**

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

#### TL-9  Group Move **(NEW)**
When multiple elements are selected (via SEL-4), the Move tool SHALL move
all selected elements by the same offset, preserving their relative
positions.

**Acceptance:** Multi-select IW3, IW7, IW9. Move tool. Drag 2 inches
west. All three walls move 2 inches west. Their relative spacing is
unchanged.

#### TL-10  Move Opening Along Wall **(NEW)**
When an opening is selected with the Move tool, dragging SHALL slide the
opening along its host wall segment (not in arbitrary directions).

**Acceptance:** Select O5. Move tool. Drag. O5 slides along the F8-F9
wall segment. Its perpendicular distance from the wall does not change.

### 6.6  Dimension Line Tool **(NEW)**

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

### 6.7  Draw Wall Tool **(NEW)**

#### TL-15  Draw Interior Wall
When the Draw Wall tool is active, clicking a start point and end point
SHALL create a new interior wall with a default thickness (4 inches).

**Acceptance:** Select Draw Wall tool. Click on W18-W19 midpoint. Click
on IW1 south face. A new interior wall appears connecting those points.
The wall has 4-inch thickness.

#### TL-16  Wall Thickness Input
After placing a wall, the Properties panel SHALL show a thickness input
that can be changed before or after placement.

**Acceptance:** Draw a new wall. Properties panel shows Thickness = 4
inches. Change to 8 inches. Wall re-renders at the new thickness.

#### TL-17  Wall Endpoint Editing
Selecting an existing wall with the Select tool and dragging an endpoint
handle SHALL extend or shorten the wall.

**Acceptance:** Select IW4. Drag its south endpoint 6 inches north. IW4
shortens. Geometry recomputes.

### 6.8  Add Element Tools **(NEW)**

#### TL-18  Add Furniture from Catalog
Tools > Add Furniture SHALL open a catalog panel listing standard items
(bed, dresser, shelves, sofa, table, chairs, hamper, etc.) with icons
and dimensions. Clicking an item enters placement mode where the next
canvas click positions the item.

**Acceptance:** Open Tools > Add Furniture. Select "HAMPER". Click on the
canvas at a position. HAMPER appears at that location. Properties panel
shows its dimensions.

#### TL-19  Add Appliance from Catalog
Tools > Add Appliance SHALL open a catalog listing standard appliances
(fridge, stove, washer, dryer, dishwasher, microwave, etc.).

**Acceptance:** Open Tools > Add Appliance. Select "MICRO". Click on the
canvas. MICRO appears at the clicked position.

#### TL-20  Add Fixture from Catalog
Tools > Add Fixture SHALL list fixtures (toilet, sink, bathtub).

**Acceptance:** Select TOILET. Click on canvas. Toilet appears.

#### TL-21  Add Opening
Tools > Add Opening SHALL enter opening placement mode. Clicking on a
wall segment SHALL place a new opening at that position. A dialog
prompts for width and type (window / rough opening / casement).

**Acceptance:** Select Add Opening. Click on the F8-F9 wall segment. A
dialog prompts for width (default 19 inches) and type. Confirm. New
opening appears on the wall. Engine recomputes geometry.

### 6.9  Delete Tool **(NEW)**

#### TL-22  Delete Selected Element
Pressing Delete or Backspace with an element selected SHALL remove the
element after a confirmation prompt.

**Acceptance:** Select HAMPER. Press Delete. Confirmation dialog appears.
Confirm. HAMPER is removed from the canvas and database.

#### TL-23  Cascading Delete
Deleting a wall SHALL also delete any openings or doors hosted on that
wall.

**Acceptance:** Delete IW9. RO3 (hosted on IW9) is also deleted. Both
disappear from the canvas.

### 6.10  Rotate Tool **(NEW)**

#### TL-24  Rotate Element
When an element is selected, Edit > Rotate or pressing R SHALL open a
rotation input. The user can type an angle or select from presets
(0, 90, 180, 270 degrees).

**Acceptance:** Select SHELVES. Press R. Rotation dialog appears with
preset buttons and angle input. Select 90. SHELVES rotates 90 degrees.

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

#### CT-7  Value Parsing -- Inches
Entering a value like `35"` SHALL be parsed as `35 / 12.0` feet.

**Acceptance:** Enter `80"` for a ft-unit constant. The value sent to the
server is `80 / 12 = 6.6667`.

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

#### DT-2  Outline Chain Editing **(NEW)**
Clicking a cell in the Dist/R, Sweep, or Bearing columns SHALL make it
editable. Pressing Enter SHALL commit the change, trigger closure
re-solve, and recompute geometry.

**Acceptance:** Click the Dist/R cell for segment 5 (an arc). Change the
value from 28 to 30. Press Enter. Engine re-solves closure. Canvas
re-renders with the updated outline.

#### DT-3  Outline Closure Indicator **(NEW)**
The Outline panel SHALL display a closure status indicator showing whether
the current chain parameters produce a valid closed outline. Red = open,
green = closed.

**Acceptance:** After a valid edit, indicator is green. After setting an
impossible radius, indicator turns red with an error distance.

#### DT-4  Add/Remove Chain Segment **(NEW)**
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

#### DT-7  Openings Table Editing **(NEW)**
Width cells in the openings tables SHALL be editable. Changing a width
SHALL update the opening constant and trigger geometry recomputation.

**Acceptance:** Click the Width cell for O8. Change from 19 to 25. Press
Enter. O8 re-renders at 25 inches wide.

#### DT-8  Table Row Selection
Clicking a row in the openings tables SHALL select the corresponding
element on the canvas and show its properties.

**Acceptance:** Click a row in the outer openings table. The
corresponding polygon on the canvas receives the `selected-highlight`
class.

### 8.3  Interior Walls Table **(NEW)**

#### DT-9  Interior Walls Table
The Elements panel SHALL display an interior walls table with columns:
Name, Thickness (inches), Length (ft-in), Orientation, Hosted Openings.

**Acceptance:** Table shows 13 rows matching ENG-4 wall names. Thickness
values are in inches. Length is formatted as ft'in".

#### DT-10  Interior Walls Table Editing **(NEW)**
Clicking a wall row SHALL select it on the canvas. Thickness and length
cells SHALL be editable.

**Acceptance:** Click IW8 row. IW8 highlights on canvas. Change thickness
from 4 to 6. Press Enter. Geometry recomputes.

### 8.4  Furniture/Appliance Table **(NEW)**

#### DT-11  Furniture/Appliance Table
The Elements panel SHALL display a furniture/appliance table with columns:
Name, Type, Width, Depth, Rotation, URL.

**Acceptance:** Table shows rows for all furniture and appliance items.
Items with product URLs show clickable links.

---

## 9  Element Operations

### 9.1  Outline Editing

#### OE-1  Drag F-Points **(NEW)**
When the Select tool is active and the "Points" toggle is on, clicking
and dragging an F-series point SHALL reshape the building outline. The
engine SHALL re-solve the outline chain parameters to match the dragged
position while maintaining tangency and closure constraints.

**Acceptance:** Drag F10 2 inches west. Engine derives new segment
parameters. Outline re-renders. Tangency holds at all arc junctions.

#### OE-2  Arc Radius Handle **(NEW)**
When an arc segment is selected, a radius adjustment handle SHALL appear.
Dragging the handle SHALL change the arc radius with live preview.

**Acceptance:** Select the F13-F14 arc. A handle appears on the arc
midpoint. Drag inward. Radius decreases. Release. Engine re-solves.

#### OE-3  Constraint-Based Outline Editing **(NEW)**
Edit > Set Constraint SHALL allow the user to specify a target distance
between two named points (e.g., "F6-F7 = 5'3""). The engine SHALL solve
for chain parameters that satisfy the constraint while maintaining
closure.

**Acceptance:** Set constraint "W20-W1 to W9-W10 clear span = 23'6"".
Engine adjusts arc radii/lengths to satisfy. Result is verified by
measuring the span in the computed geometry.

### 9.2  Door Configuration

#### DOOR-1  Add Door to Opening **(NEW)**
Right-clicking an opening SHALL offer "Add Door". A dialog SHALL prompt
for door width, hinge side, swing direction, and type (single/double).

**Acceptance:** Right-click RO1. Select "Add Door". Dialog appears. Set
hinge = east, swing = south, type = single. Confirm. Door arc renders on
the canvas.

#### DOOR-2  Edit Door Hinge and Swing **(NEW)**
Selecting a door arc and clicking a "Flip Hinge" button SHALL move the
hinge to the opposite side. "Flip Swing" SHALL reverse the swing
direction.

**Acceptance:** Select RO3 door. Click "Flip Hinge". Door hinge moves
from south to north. Door arc re-renders.

#### DOOR-3  Appliance Door Arcs **(NEW)**
Appliances (fridge, washer, dryer, microwave) SHALL support door arc
configuration with hinge corner and swing direction.

**Acceptance:** Select FRIDGE. Add door with hinge at NE corner, swing
NW. A door arc renders showing the fridge door swing.

#### DOOR-4  Double Door **(NEW)**
Openings SHALL support double-door configuration with two leaves, each
with independent hinge side.

**Acceptance:** Configure RO7 as a double door. Two door arcs render,
one hinged on each side of the opening.

### 9.3  Hyperlinks

#### LINK-1  Attach Product URL **(NEW)**
Each furniture, appliance, fixture, and opening element SHALL have an
optional product URL property.

**Acceptance:** Select SHELVES. In Properties panel, enter an IKEA URL.
Save. The URL persists across page reloads.

#### LINK-2  Clickable SVG Links **(NEW)**
Elements with product URLs SHALL be wrapped in `<a xlink:href="...">` tags
in generated SVGs, making them clickable in browsers and PDF viewers.

**Acceptance:** Assign a URL to FRIDGE. Regenerate floorplan SVG. Open the
SVG in a browser. Clicking FRIDGE opens the product URL.

### 9.4  Room Labels & Annotations

#### LABEL-1  Add Room Label **(NEW)**
Edit > Add Room Label SHALL allow the user to place a room name label by
clicking a position on the canvas and typing the label text (e.g.,
"BEDROOM", "OFFICE").

**Acceptance:** Add Room Label. Click in the bedroom area. Type "BEDROOM".
Press Enter. Label appears on the canvas and persists in the database.

#### LABEL-2  Move Label **(NEW)**
With the Select tool, labels SHALL be draggable to reposition them.

**Acceptance:** Select the "OFFICE" room label. Drag it 1 foot north.
Label position updates. Generated SVGs reflect the new position.

#### LABEL-3  Edit Label Text **(NEW)**
Double-clicking a label SHALL open an inline text editor to change the
label text.

**Acceptance:** Double-click "BATH" label. Text becomes editable. Change
to "BATHROOM". Press Enter. Label updates.

#### LABEL-4  Label Font Size **(NEW)**
Each label SHALL have a configurable font size property, editable in the
Properties panel.

**Acceptance:** Select a room label. Properties panel shows Font Size
field. Change from 10pt to 14pt. Label re-renders at the new size.

---

## 10  Styling **(NEW)**

#### STYLE-1  Element Fill Colour
Each element type (wall, opening, furniture, etc.) SHALL have a
configurable fill colour, editable via the Properties panel or a colour
picker.

**Acceptance:** Select a wall. Open colour picker. Change fill from
default grey to light brown. Wall re-renders with the new fill colour.

#### STYLE-2  Element Stroke Properties
Each element SHALL have configurable stroke colour, width, and style
(solid, dashed, dotted).

**Acceptance:** Select a dimension line. Change stroke style to dashed.
Line re-renders as dashed.

#### STYLE-3  Element Opacity
Each element SHALL have a configurable opacity (0-100%), editable in the
Properties panel.

**Acceptance:** Select interior dim lines. Set opacity to 50%. Lines
render at half opacity.

#### STYLE-4  Per-View Styling **(NEW)**
View variants SHALL support per-view element styling (e.g., doors at 20%
opacity in the plumbing view).

**Acceptance:** In plumbing view variant, set door opacity to 20%. Doors
render faintly. Main floorplan view still shows doors at 100%.

---

## 11  Site Plan **(NEW)**

#### SITE-1  Structure Placement
The site plan view SHALL allow positioning the building outline on the
survey parcel by specifying setback distances from property boundaries.

**Acceptance:** Open site plan view. Enter north setback = 11.5 feet.
Building repositions. Setback distance labels update.

#### SITE-2  Drainfield Operations
The site plan SHALL support adding, positioning, and sizing drainfield
rectangles.

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

#### SCAD-1  SCAD Generation
Tools > Generate 3D Model SHALL run the OpenSCAD generator and produce a
3D model file.

**Acceptance:** Click Generate 3D Model. SCAD file is generated. Success
toast appears.

#### SCAD-2  Roof Style Selection
The 3D model settings SHALL allow selecting roof style (flat, 2:12 slope)
and overhang distance.

**Acceptance:** Select "2:12 slope" in 3D settings. Regenerate. SCAD file
contains sloped roof geometry.

#### SCAD-3  Multi-View Layout
Tools > Generate Views SHALL produce a multi-view PDF (3views.pdf) with
floor plan and elevation views.

**Acceptance:** Click Generate Views. PDF is created with 4 panels
(floorplan + 3 elevations).

---

## 13  Analysis **(NEW)**

#### ANALYSIS-1  Span Analysis View
Tools > Span Analysis SHALL display the N-S interior span graph showing
span distance vs. position.

**Acceptance:** Open Span Analysis. Graph shows span curve. Mouse hover
displays span value at each position.

#### ANALYSIS-2  Span vs. Rotation
Tools > Span vs. Rotation SHALL display span measurements across rotation
angles (5-175 degrees).

**Acceptance:** Graph shows span vs. rotation with min/max markers.

#### ANALYSIS-3  Room Area Display
View > Show Areas SHALL compute and display room areas in square feet for
each enclosed region.

**Acceptance:** Room areas appear within room boundaries. Total area sums
correctly. Values update when walls are moved.

---

## 14  Plumbing Plan **(NEW)**

#### PLUMB-1  Plumbing View
The plumbing plan view SHALL display the building outline with fixture
locations, supply lines, and drain lines.

**Acceptance:** Switch to Plumbing tab. SVG shows fixtures (toilets, sinks,
washer) with blue supply lines and coloured drain lines.

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

---

## 15  Undo/Redo **(NEW)**

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

#### UNDO-4  Undo Across Element Types
Undo SHALL work for all element operations: constant edits, element moves,
adds, deletes, door changes, outline edits, label changes.

**Acceptance:** Add a door, move a wall, change a constant. Undo three
times. All three operations are correctly reverted.

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

#### RT-5  Element Change Notification **(NEW)**
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

**Acceptance:** `python -m pytest tests/ -x -q` reports 586 passed, 0
failed.

#### NF-4  No Modification of Existing Code
The app SHALL not modify any files in `shared/`, `floorplan/`, `walls/`,
`span/`, `survey/`, `roof/`, `site/`, `scad/`, or `plumbing/`.

**Acceptance:** `git diff --name-only` shows changes only in `app/`,
`run_app.py`, `pyproject.toml`, and `.gitignore`.

#### NF-5  Keyboard Shortcut Safety
Keyboard shortcuts SHALL NOT fire when the user is typing in an `<input>`
or `<select>` element.

**Acceptance:** Focus a constant value input. Press V. The tool does not
change; the character is typed into the input instead.

#### NF-6  Default Display Units
Unless otherwise specified, all displayed dimensions SHALL use feet and
inches. Inch values SHALL be displayed to two decimal places of precision
with trailing zeroes removed (e.g., `5' 3.5"` not `5' 3.50"`; `12' 0"`
not `12' 0.00"`).

**Acceptance:** Inspect dimension labels, property panel values, and
measurement tool output. All use feet-inches format. A value of exactly
3 inches displays as `3"`, not `3.00"`. A value of 7.50 inches displays
as `7.5"`.

---

## Appendix A: Requirement Cross-Reference by User Operation

This table maps common user operations (identified from the full 864-commit
history) to the requirements that enable each operation through the GUI.

| Operation | Commits | Requirements |
|-----------|---------|-------------|
| Move element (wall/furniture/opening) | ~146 | TL-5, TL-6, TL-7, TL-8, TL-9, TL-10, API-23 |
| Edit outline geometry (F-series chain) | ~79 | OE-1, OE-2, OE-3, DT-2, DT-3, DT-4, API-16..19 |
| Add/edit dimension lines & labels | ~78 | TL-11..14, CV-9, DIS-7, LABEL-1..4 |
| Add/resize/rotate furniture | ~67 | TL-18..20, TL-24, SEL-8, DT-11 |
| Add/edit/remove interior walls | ~52 | TL-15..17, TL-22, TL-23, DT-9, DT-10 |
| Add/edit openings | ~42 | TL-21, DT-7, SEL-10, API-24..26 |
| Configure doors (hinge/swing) | ~35 | DOOR-1..4, SEL-7, SEL-11, CV-7, API-27..29 |
| Edit site plan | ~36 | SITE-1..4 |
| 3D/SCAD model | ~31 | SCAD-1..3 |
| View variants | ~24 | UI-5, UI-6, UI-8, ENG-13, API-8, API-34, CV-4 |
| Plumbing plan | ~22 | PLUMB-1..3 |
| Span/area analysis | ~22 | ANALYSIS-1..3 |
| Resize elements | ~20 | SEL-10, DT-7, DT-10 |
| Element styling (colour/opacity) | ~18 | STYLE-1..4 |
| Product hyperlinks | ~14 | LINK-1, LINK-2, SEL-12, CV-12 |
| Delete elements | ~8 | TL-22, TL-23, API-22, API-26, API-29 |
| Undo/redo | implicit | UNDO-1..4, API-30, API-31 |

## Appendix B: Requirements Summary

| Section | Existing | New | Total |
|---------|----------|-----|-------|
| 1 Data Layer | 10 | 6 | 16 |
| 2 REST API | 17 | 17 | 34 |
| 3 UI Layout | 5 | 4 | 9 |
| 4 Canvas | 11 | 10 | 21 |
| 5 Selection | 6 | 8 | 14 |
| 6 Tools | 4 | 20 | 24 |
| 7 Constants | 10 | 0 | 10 |
| 8 Data Tables | 4 | 7 | 11 |
| 9 Element Ops | 0 | 14 | 14 |
| 10 Styling | 0 | 4 | 4 |
| 11 Site Plan | 0 | 4 | 4 |
| 12 3D Model | 0 | 3 | 3 |
| 13 Analysis | 0 | 3 | 3 |
| 14 Plumbing | 0 | 3 | 3 |
| 15 Undo/Redo | 0 | 4 | 4 |
| 16 Real-Time | 4 | 1 | 5 |
| 17 Application | 9 | 1 | 10 |
| **Total** | **80** | **109** | **189** |
