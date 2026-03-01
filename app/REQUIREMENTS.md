# ADU Editor — Requirements

Testable requirements for the interactive building editor application.
Each requirement has a unique ID, a category, a description, and acceptance
criteria that can be verified by automated or manual testing.

---

## 1  Database Layer

### DB-1  Schema Initialisation
The application SHALL create an SQLite database with tables `constants`,
`outline_chain`, and `views` when launched for the first time.

**Acceptance:** Start with no `app/adu.db` file. Run `python run_app.py
--no-browser`. Verify the file is created and contains all three tables.

### DB-2  Constants Seeding
On first initialisation the database SHALL contain every uppercase numeric
constant defined in `floorplan/constants.py`.

**Acceptance:** Compare the set of constant names in the database against
`dir(floorplan.constants)` filtered to uppercase float/int attributes.
Count must be ≥ 140.

### DB-3  Constant Categories
Each constant SHALL be assigned to exactly one category from the set
{wall, interior_wall, opening, appliance, furniture, fixture, geometry,
construction, misc} based on its name prefix.

**Acceptance:** Query `SELECT DISTINCT category FROM constants` and verify
the result set equals the set above.

### DB-4  Outline Chain Seeding
On first initialisation the database SHALL contain exactly 18 outline
chain rows matching the `OUTLINE_CHAIN` list in `floorplan/geometry.py`.

**Acceptance:** `SELECT count(*) FROM outline_chain` returns 18.
Each row's `seg_type` is one of {'L', 'CW', 'CCW'}. Each `end_name`
matches the corresponding entry in `CHAIN_POINT_NAMES`.

### DB-5  Views Seeding
On first initialisation the database SHALL register at least 11 views
covering the floorplan, walls, span, survey, roof, plumbing, and site
plan generators.

**Acceptance:** `SELECT count(*) FROM views WHERE enabled = 1` ≥ 11.
Each row has a non-empty `script` path pointing to an existing Python file.

### DB-6  Constant Update
`update_constant(name, value)` SHALL change the stored value and return
`True`. A subsequent `get_constants_dict()` SHALL reflect the new value.

**Acceptance:** Update `BED_WIDTH` to `7.0`. Read it back. Verify it
equals `7.0`. Reset constants afterwards.

### DB-7  Batch Update
`update_constants_batch(updates)` SHALL update all named constants in a
single transaction and return the number of rows changed.

**Acceptance:** Batch-update `BED_WIDTH` and `BED_LENGTH` to new values.
Verify both changed. Return value equals 2.

### DB-8  Reset Constants
`reset_constants()` SHALL delete all rows from the `constants` table and
re-seed from `floorplan/constants.py`, restoring original values.

**Acceptance:** Modify `BED_WIDTH`. Call `reset_constants()`. Read
`BED_WIDTH`. Verify it equals the value in the Python source module.

---

## 2  Geometry Engine

### ENG-1  Geometry Computation
`compute_geometry(constants_dict)` SHALL return a dict containing at least
the keys: `points`, `outline_segments`, `inner_segments`, `outline_poly`,
`inner_poly`, `interior_walls`, `outer_openings`, `rough_openings`,
`appliances`, `furniture`, `bbox`.

**Acceptance:** Call `compute_geometry` with the default constants dict.
Verify all listed keys exist and have non-empty values.

### ENG-2  Point Count
The computed geometry SHALL include at least 50 named points covering
F-series, W-series, C-series, and survey points.

**Acceptance:** `len(result["points"]) >= 50`.

### ENG-3  Outline Segment Count
The computed geometry SHALL include exactly 18 outline segments and 18
inner wall segments.

**Acceptance:** `len(result["outline_segments"]) == 18` and
`len(result["inner_segments"]) == 18`.

### ENG-4  Interior Wall Count
The computed geometry SHALL include exactly 13 interior walls: IW1, IW2,
IW2O, IW2S, IW3, IW4, IW5, IW6, IW7, IW8, IW9, IW11, IW12.

**Acceptance:** `set(result["interior_walls"].keys())` equals the expected
set.

### ENG-5  Opening Counts
The computed geometry SHALL include exactly 12 outer openings and 7 rough
openings.

**Acceptance:** `len(result["outer_openings"]) == 12` and
`len(result["rough_openings"]) == 7`.

### ENG-6  Appliance and Furniture Counts
The `appliances` dict SHALL contain keys `dryer`, `washer`, `counter`.
The `furniture` dict SHALL contain keys `bed`, `dresser`, `shelves`.

**Acceptance:** Verify both key sets.

### ENG-7  Constant Propagation
Changing a constant in the database SHALL produce different computed
geometry when `compute_geometry` is called again.

**Acceptance:** Set `BED_WIDTH` to `80/12`. Recompute geometry. Verify the
bed bounding box width is 80" (≈ 6.6667 ft). Reset and verify original
width is 76".

### ENG-8  Derived Constants
When `WALL_OUTER` is patched, the engine SHALL recompute the derived
constants `WALL_EXTRA`, `AIR_GAP`, `DOOR_FLAT_FACE`, `F8F9_INNER_TURN_R`,
and `CORNER_SW_R` before computing geometry.

**Acceptance:** Patch `WALL_OUTER` to 10"/12. Verify `AIR_GAP` equals
`WALL_OUTER − 2 × SHELL_THICKNESS`.

### ENG-9  SVG File Reading
`get_svg_content(svg_path)` SHALL return the file contents as a string
when the file exists, and `None` when it does not.

**Acceptance:** Read `floorplan/floorplan.svg` (which exists after
`gen_all.py`). Verify result starts with `<` (XML). Read
`nonexistent.svg`. Verify result is `None`.

### ENG-10  SVG Generation
`generate_svg(view_name, script_path)` SHALL run the generator script and
return `True` on success, `False` on failure or timeout.

**Acceptance:** Call with `("floorplan", "floorplan/gen_floorplan.py")`.
Verify it returns `True` and the SVG file is updated.

---

## 3  REST API

### API-1  GET /
The root route SHALL return HTTP 200 with HTML containing the string
"ADU Editor".

**Acceptance:** `GET /` returns status 200. Body contains `"ADU Editor"`.

### API-2  GET /api/constants
SHALL return HTTP 200 with a JSON array of constant objects. Each object
SHALL have keys: `name`, `value`, `expr`, `unit`, `category`,
`description`.

**Acceptance:** Response is a JSON array with length ≥ 140. First element
has all six keys.

### API-3  GET /api/constants?category=wall
SHALL return only constants whose `category` equals the query parameter.

**Acceptance:** Every object in the response has `category == "wall"`.

### API-4  GET /api/constants/categories
SHALL return HTTP 200 with a JSON array of distinct category strings.

**Acceptance:** Response is a JSON array containing at least the strings
`"wall"`, `"opening"`, `"furniture"`.

### API-5  PUT /api/constants/<name>
SHALL accept a JSON body `{"value": <number>}`, update the constant, and
return `{"ok": true, "name": ..., "value": ...}`.

**Acceptance:** PUT `BED_WIDTH` with `{"value": 7.0}`. Response status is
200, `ok` is `true`, `value` is `7.0`.

### API-6  PUT /api/constants/<name> — validation
SHALL return HTTP 400 if `value` is missing or not a valid number.
SHALL return HTTP 404 if the constant name does not exist.

**Acceptance:** PUT with `{}` body returns 400. PUT with
`{"value": "abc"}` returns 400. PUT to `/api/constants/NONEXISTENT`
returns 404.

### API-7  PUT /api/constants/batch
SHALL accept `{"updates": {"A": 1.0, "B": 2.0}}` and return
`{"ok": true, "changed": N}` where N is the count of rows updated.

**Acceptance:** Batch-update two known constants. Verify `changed == 2`.

### API-8  POST /api/constants/reset
SHALL reset all constants to their original source values and return
`{"ok": true}`.

**Acceptance:** Modify a constant. POST reset. GET the constant. Value
matches the original.

### API-9  GET /api/geometry
SHALL return HTTP 200 with a JSON object containing the complete computed
geometry (see ENG-1 for required keys).

**Acceptance:** Response is JSON with keys `points`, `outline_poly`,
`interior_walls`, `bbox`.

### API-10  GET /api/geometry — error handling
SHALL return HTTP 500 with `{"error": "..."}` if geometry computation
fails.

**Acceptance:** (Edge case — verified by code inspection of the try/except
block.)

### API-11  GET /api/outline
SHALL return HTTP 200 with a JSON array of 18 outline chain segment
objects.

**Acceptance:** Response length is 18. Each object has keys `seq`,
`seg_type`, `end_name`.

### API-12  GET /api/views
SHALL return HTTP 200 with a JSON array of enabled view objects. Each has
keys `name`, `label`, `script`, `svg_path`, `category`.

**Acceptance:** Response length ≥ 11. Each element has all five keys.

### API-13  GET /api/svg/<view_name>
SHALL return the SVG file content with MIME type `image/svg+xml` when the
file exists. SHALL return HTTP 404 when the SVG has not been generated or
the view name is unknown.

**Acceptance:** GET `/api/svg/floorplan` returns 200 with content starting
with `<`. GET `/api/svg/nonexistent` returns 404.

### API-14  GET /api/svg/<view_name>/file
SHALL serve the file as an attachment. SVG files use MIME `image/svg+xml`;
PDF files use MIME `application/pdf`.

**Acceptance:** GET `/api/svg/floorplan/file` returns 200. Content-Type
header contains `svg`.

### API-15  POST /api/regenerate
With no body: SHALL regenerate all enabled views and return
`{"ok": true, "results": {...}}`.
With body `{"view": "floorplan"}`: SHALL regenerate only that view.

**Acceptance:** POST with `{"view": "floorplan"}` returns 200 with
`ok == true` and `view == "floorplan"`. POST with no body returns 200 with
`results` dict.

### API-16  GET /api/events (SSE)
SHALL return a streaming response with MIME type `text/event-stream`.
The first event SHALL be `event: connected`.

**Acceptance:** Connect to the endpoint. Read at least the first event.
Verify it is `event: connected\ndata: {}\n\n`.

### API-17  Geometry Cache Invalidation
After a successful PUT to `/api/constants/<name>`, the geometry cache
SHALL be marked dirty so the next GET `/api/geometry` recomputes.

**Acceptance:** GET geometry (cache populated). PUT a constant. GET
geometry again. Verify the response reflects the changed constant.

---

## 4  Frontend — Layout & Navigation

### UI-1  Page Structure
The page SHALL display a menu bar (top), tool palette (left), viewport
(centre), property/data panel (right), and status bar (bottom).

**Acceptance:** Load the page. All five regions are visible and occupy
their expected positions.

### UI-2  Menu Bar
The menu bar SHALL have menus: File, Edit, View, Tools. Each opens a
dropdown on hover.

**Acceptance:** Hover over each menu name. A dropdown with at least one
action button appears.

### UI-3  View Tabs
The viewport area SHALL display a tab bar. The first tab SHALL be
"Interactive". Additional tabs SHALL correspond to registered SVG views.

**Acceptance:** Load the page. The Interactive tab is visible and active.
At least one additional tab (e.g. "Floorplan") is visible.

### UI-4  Tab Switching
Clicking a view tab SHALL switch the viewport content. Clicking
"Interactive" shows the SVG canvas. Clicking a generated view tab loads
and displays the corresponding SVG file.

**Acceptance:** Click "Interactive" — canvas is visible, SVG container is
hidden. Click "Floorplan" — canvas is hidden, SVG container shows
floorplan SVG or a "not generated" message.

### UI-5  Right Panel Tabs
The right panel SHALL have tabs: Properties, Constants, Outline, Openings.
Clicking a tab shows only that panel's content.

**Acceptance:** Click each tab. The corresponding panel content is visible;
all others are hidden.

---

## 5  Frontend — Interactive Canvas

### CV-1  Outline Rendering
The interactive canvas SHALL render the building outline as a filled
polygon with a stroke border.

**Acceptance:** Load the page. A polygon shape is visible in the canvas
area. The polygon has ≥ 40 vertices (arcs are polygonised).

### CV-2  Interior Wall Rendering
The canvas SHALL render all 13 interior walls as filled polygons.

**Acceptance:** Count elements in `#layer-walls`. There are 13 polygon
elements.

### CV-3  Opening Rendering
When the "Openings" toggle is checked, the canvas SHALL render all 12
outer openings and 7 rough openings.

**Acceptance:** Count elements in `#layer-openings`. There are 19 polygon
elements.

### CV-4  Furniture Rendering
When the "Furniture" toggle is checked, the canvas SHALL render appliances
(dryer, washer, counter) and furniture (bed, dresser, shelves).

**Acceptance:** Count elements in `#layer-furniture`. There are 6 polygon
elements.

### CV-5  Point Markers
When the "Points" toggle is checked, the canvas SHALL render circle
markers for F-series, W-series, and C-series points.

**Acceptance:** `#layer-points` contains circle elements. F-series circles
use the accent colour; W-series circles use green.

### CV-6  Point Labels
When the "Labels" toggle is checked, the canvas SHALL render text labels
adjacent to each point marker showing the point name.

**Acceptance:** `#layer-labels` contains text elements with content such
as `"F1"`, `"W2"`, `"C5"`.

### CV-7  Fit to Window
On initial load or when the user presses F or selects View > Fit to
Window, the building SHALL be centred and scaled to fill the viewport with
a margin.

**Acceptance:** Load the page. The outline is fully visible within the
viewport. No scrolling is needed.

### CV-8  Mouse Wheel Zoom
Scrolling the mouse wheel SHALL zoom the canvas toward or away from the
cursor position.

**Acceptance:** Place cursor over the building. Scroll up — the building
gets larger. Scroll down — it gets smaller. The area under the cursor
remains roughly stationary.

### CV-9  Pan with Middle Mouse Button
Pressing and dragging with the middle mouse button SHALL pan the canvas.

**Acceptance:** Middle-click and drag. The building moves with the cursor.

### CV-10  Pan Tool
When the Pan tool is active, left-click drag SHALL pan the canvas. The
cursor SHALL be `grab` (at rest) or `grabbing` (while dragging).

**Acceptance:** Select Pan tool (click button or press H). Left-drag
pans the canvas. Cursor changes appropriately.

### CV-11  Coordinate Display
As the mouse moves over the viewport, the status bar SHALL display the
current cursor position in world coordinates (feet and inches).

**Acceptance:** Move mouse. `#coord-display` updates continuously showing
E and N values.

---

## 6  Frontend — Selection & Properties

### SEL-1  Click-to-Select
Clicking on a wall, opening, appliance, or furniture polygon SHALL select
that element. The element SHALL be visually highlighted and the Properties
panel SHALL populate with its details.

**Acceptance:** Click on an interior wall polygon. It receives the CSS
class `selected-highlight`. The Properties panel title shows the wall
name.

### SEL-2  Point Selection
Clicking on a point marker SHALL select the point and show its Easting
and Northing in feet and inches.

**Acceptance:** Click on a point circle. Properties panel shows E, N,
and their inch equivalents.

### SEL-3  Clear Selection
Clicking on the viewport background or pressing Escape SHALL deselect any
selected element and clear the Properties panel.

**Acceptance:** Select an element. Click background. The
`selected-highlight` class is removed. Properties panel shows the empty
state.

### SEL-4  Wall Properties
When a wall is selected, the Properties panel SHALL display its width
(inches), height (inches), and bounding box coordinates (feet). It SHALL
also list related constants with editable input fields.

**Acceptance:** Select IW1. Panel shows Width, Height, West, South, East,
North values. Constants starting with "IW1" appear as editable rows.

### SEL-5  Opening Properties
When an opening is selected, the Properties panel SHALL display its name,
segment reference, width, and actual computed width.

**Acceptance:** Select O3. Panel shows Name, Segment (F2–F5), and Width
values.

### SEL-6  Inline Constant Editing
Editing a related-constant value in the Properties panel and pressing
Enter SHALL send a PUT request to update the constant and trigger geometry
recomputation.

**Acceptance:** Select a wall. Edit a related constant input. A PUT
request is sent. A success toast appears. The geometry is reloaded.

---

## 7  Frontend — Constants Table

### CT-1  Table Population
The Constants panel SHALL display a table of all constants with columns:
Name, Value, Unit, Category.

**Acceptance:** Switch to the Constants tab. The table contains ≥ 140
rows.

### CT-2  Value Formatting
Constants with unit `ft` SHALL display their value in inches with a `"`
suffix. Other units SHALL display 6 decimal places.

**Acceptance:** A constant like `BED_WIDTH = 76/12 ft` displays as `76"`.

### CT-3  Category Filtering
Selecting a category from the dropdown SHALL show only constants in that
category.

**Acceptance:** Select "wall". Only constants with category "wall" are
visible.

### CT-4  Search Filtering
Typing in the search box SHALL filter constants by name or description
(case-insensitive).

**Acceptance:** Type "bed". Only constants whose name or description
contains "bed" are shown.

### CT-5  Column Sorting
Clicking a sortable column header SHALL sort the table by that column.
Clicking again SHALL reverse the sort direction.

**Acceptance:** Click "Name" header. Constants are sorted A–Z. Click
again. Sorted Z–A.

### CT-6  Inline Editing
Each value cell SHALL contain an editable input. Changing the value and
pressing Enter SHALL send a PUT request to the server.

**Acceptance:** Change BED_WIDTH to `80"` and press Enter. A PUT is sent.
Toast confirms the update.

### CT-7  Value Parsing — Inches
Entering a value like `35"` SHALL be parsed as `35 / 12.0` feet.

**Acceptance:** Enter `80"` for a ft-unit constant. The value sent to the
server is `80 / 12 ≈ 6.6667`.

### CT-8  Value Parsing — Fractions
Entering a value like `1/3` SHALL be parsed as `0.333...`.

**Acceptance:** Enter `1/3`. The value sent is approximately `0.3333`.

### CT-9  Category Colour Coding
Each category cell SHALL be colour-coded (e.g. wall = blue, furniture =
green, opening = yellow).

**Acceptance:** Inspect category cells. Different categories have visually
distinct colours.

### CT-10  Changed Value Indicator
When an input value differs from its original loaded value, the input
SHALL receive a visual indicator (yellow border).

**Acceptance:** Change a value without pressing Enter. The input border
turns yellow.

---

## 8  Frontend — Data Tables

### DT-1  Outline Chain Table
The Outline panel SHALL display a table with columns: #, Type, Dist/R,
Sweep, End. It SHALL contain 18 rows.

**Acceptance:** Switch to the Outline tab. Table has 18 rows. Line
segments show distance; arc segments show radius and sweep.

### DT-2  Outer Openings Table
The Openings panel SHALL display outer openings (O1–O11, O8a) with
columns: Name, Segment, Width (inches).

**Acceptance:** Panel shows 12 rows with correct opening names.

### DT-3  Rough Openings Table
The Openings panel SHALL display rough openings (RO1–RO7) with columns:
Name, Wall, Width (inches), Orientation.

**Acceptance:** Panel shows 7 rows. Each row has a wall name (e.g. IW1)
and orientation (H, V, or R).

### DT-4  Table Row Selection
Clicking a row in the openings tables SHALL select the corresponding
element on the canvas and show its properties.

**Acceptance:** Click a row in the outer openings table. The
corresponding polygon on the canvas receives the `selected-highlight`
class.

---

## 9  Frontend — Tools

### TL-1  Tool Palette
The left palette SHALL display buttons for Select, Pan, Measure, and Move
tools. Exactly one tool SHALL be active at any time.

**Acceptance:** Page loads with Select active. Click Pan — only Pan is
highlighted.

### TL-2  Tool Keyboard Shortcuts
The keyboard shortcuts V (Select), H (Pan), M (Measure), G (Move), and
F (Fit) SHALL activate the corresponding tool or action when no text input
is focused.

**Acceptance:** Press V — Select tool is active. Press H — Pan tool is
active. Press F — viewport fits to window.

### TL-3  Measure Tool
When the Measure tool is active, clicking and dragging SHALL draw a
dashed red line between two points. A label SHALL show the distance in
feet and inches. The status bar SHALL display the distance on mouse-up.

**Acceptance:** Select Measure tool. Click-drag between two points on the
building. A red dashed line and distance label appear. The `#measure-info`
element shows the distance.

### TL-4  Escape Clears Measure
Pressing Escape while a measurement is displayed SHALL clear the
measurement line and info.

**Acceptance:** Measure between two points. Press Escape. Line and info
text are gone.

---

## 10  Frontend — Display Toggles

### DIS-1  Points Toggle
Unchecking the "Points" checkbox SHALL hide all point markers and
re-render the canvas without them.

**Acceptance:** Uncheck Points. `#layer-points` is empty. Check Points.
Markers reappear.

### DIS-2  Labels Toggle
Unchecking the "Labels" checkbox SHALL hide all point labels.

**Acceptance:** Uncheck Labels. `#layer-labels` is empty.

### DIS-3  Grid Toggle
Checking the "Grid" checkbox SHALL show a 1-foot grid overlay.

**Acceptance:** Check Grid. The `#grid-rect` element's `visibility`
attribute changes to `visible`.

### DIS-4  Openings Toggle
Unchecking the "Openings" checkbox SHALL hide all opening polygons.

**Acceptance:** Uncheck Openings. `#layer-openings` is empty.

### DIS-5  Furniture Toggle
Unchecking the "Furniture" checkbox SHALL hide all appliance and furniture
polygons.

**Acceptance:** Uncheck Furniture. `#layer-furniture` is empty.

---

## 11  Real-Time Updates

### RT-1  SSE Connection
On page load the browser SHALL establish a Server-Sent Events connection
to `/api/events` and display a green "Connected" indicator.

**Acceptance:** Load the page. The connection status shows "Connected"
with a green background.

### RT-2  Geometry Change Notification
When a constant is updated via the API, the server SHALL broadcast a
`geometry_changed` event. All connected browsers SHALL reload geometry
and re-render the canvas.

**Acceptance:** Open two browser tabs. In tab A, edit a constant. Tab B
automatically re-renders with the new geometry.

### RT-3  SVG Update Notification
When a view is regenerated via the API, the server SHALL broadcast an
`svg_updated` event with the view name. If a browser is viewing that SVG,
it SHALL reload automatically.

**Acceptance:** Open the app viewing the Floorplan SVG tab. Trigger
regeneration via the File menu. A toast appears and the SVG reloads.

### RT-4  Reconnection
If the SSE connection drops, the browser SHALL attempt to reconnect after
3 seconds and restore the "Connected" indicator on success.

**Acceptance:** Stop and restart the server. Within a few seconds the
browser reconnects and the status returns to "Connected".

---

## 12  Product File Generation

### GEN-1  Regenerate Single View
File > Regenerate Current View SHALL run the generator script for the
active SVG view and reload its content in the viewport.

**Acceptance:** Switch to the Floorplan tab. Click Regenerate Current
View. The SVG content refreshes. A toast confirms completion.

### GEN-2  Regenerate All
File > Regenerate All SHALL run all registered generator scripts and
broadcast SVG update events.

**Acceptance:** Click Regenerate All. All generator scripts run. Toasts
appear for each updated view.

### GEN-3  Reset to Defaults
File > Reset to Defaults SHALL prompt for confirmation, then reset all
constants and reload the Constants table.

**Acceptance:** Click Reset to Defaults. A confirmation dialog appears.
Confirm. Constants table refreshes with original values.

### GEN-4  Interactive View Not Regenerable
Clicking Regenerate Current View while on the Interactive tab SHALL show
a toast explaining that the interactive view cannot be regenerated.

**Acceptance:** Select Interactive tab. Click Regenerate Current View.
Toast reads "Cannot regenerate interactive view".

---

## 13  Application Launcher

### APP-1  Default Start
`python run_app.py` SHALL start the server on `127.0.0.1:5000` and open
the default browser.

**Acceptance:** Run the command. Server log shows "Running on
http://127.0.0.1:5000". Browser opens to that URL.

### APP-2  Custom Port
`python run_app.py --port 8080` SHALL start the server on port 8080.

**Acceptance:** Run the command. Server log shows port 8080.

### APP-3  No-Browser Flag
`python run_app.py --no-browser` SHALL start the server without opening
a browser.

**Acceptance:** Run the command. No browser window opens. Server is
accessible via HTTP.

### APP-4  Database Auto-Creation
If `app/adu.db` does not exist, the server SHALL create and seed it
automatically on startup.

**Acceptance:** Delete `app/adu.db`. Start the server. Verify the file is
created with seeded data.

---

## 14  Non-Functional

### NF-1  Dark Theme
The UI SHALL use a dark colour scheme with light text on dark backgrounds.

**Acceptance:** Visual inspection — background is dark (#1e1e2e), text is
light.

### NF-2  Responsive Layout
At viewport widths below 1000 px, the tool palette SHALL collapse to icon-
only mode and the right panel SHALL narrow.

**Acceptance:** Resize window below 1000 px. Tool button labels disappear.
Panel width shrinks.

### NF-3  Existing Tests Unaffected
All pre-existing tests in `tests/` SHALL continue to pass after the app
is added.

**Acceptance:** `python -m pytest tests/ -x -q` reports 586 passed, 0
failed.

### NF-4  No Modification of Existing Code
The app SHALL not modify any files in `shared/`, `floorplan/`, `walls/`,
`span/`, `survey/`, `roof/`, `site/`, `scad/`, or `plumbing/`.

**Acceptance:** `git diff --name-only` shows changes only in `app/`,
`run_app.py`, `pyproject.toml`, and `.gitignore`.

### NF-5  Keyboard Shortcut Safety
Keyboard shortcuts SHALL NOT fire when the user is typing in an `<input>`
or `<select>` element.

**Acceptance:** Focus a constant value input. Press V. The tool does not
change; the character is typed into the input instead.
