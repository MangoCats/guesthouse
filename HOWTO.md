# HOWTO — Common Tasks Reference

Step-by-step instructions for common but complex tasks in the Hut2 project. Consult this before researching the codebase from scratch.

## Table of Contents

1. [Adding a Dimension Line](#1-adding-a-dimension-line)
2. [Identifying Wall Faces and Coordinates](#2-identifying-wall-faces-and-coordinates)
3. [Adding an Interior Wall](#3-adding-an-interior-wall)
4. [Adding an Opening](#4-adding-an-opening)
5. [Adding an Appliance or Furniture Item](#5-adding-an-appliance-or-furniture-item)
6. [Wall Construction Detail Drawing](#6-wall-construction-detail-drawing)
7. [Verifying Changes](#7-verifying-changes)
8. [Adding Text Captions to site_plan.pdf](#8-adding-text-captions-to-site_planpdf)
9. [Contributing to This Document](#9-contributing-to-this-document)

---

## 1. Adding a Dimension Line

**File:** `floorplan/gen_floorplan.py` (dimension lines section, after wall definitions, before openings)

### Helper functions

Two helpers are defined at the top of `gen_floorplan.py`:

- `dim_line_h(out, e1, n, e2, label)` — Horizontal (W-E) dimension line at a fixed northing `n`, between eastings `e1` and `e2`. Tick marks are vertical. Label is centered above.
- `dim_line_v(out, e, n1, n2, label)` — Vertical (S-N) dimension line at a fixed easting `e`, between northings `n1` and `n2`. Tick marks are horizontal. Label is rotated 90° to the left.

### Coordinates

All dimension helper parameters use **survey coordinates** (Easting/Northing in feet), not SVG pixels. The `to_svg()` conversion is handled internally.

### Label formatting

Use `fmt_dist(distance_in_feet)` from `shared/geometry.py` to format distances as `X' Y"`. The distance argument should be positive.

### Step-by-step

1. **Identify the two endpoints** in survey coordinates:
   - Perimeter wall outer face: `pts["F<n>"][0]` (easting) or `pts["F<n>"][1]` (northing)
   - Perimeter wall inner face: `pts["W<n>"][0]` or `pts["W<n>"][1]`
   - Interior wall faces: use `layout.<field>` (e.g., `layout.iw3.w`, `layout.iw3.e`, `layout.iw4_w`, `layout.iw1_s`)

2. **Determine the offset position** (where the line is placed, perpendicular to the measurement):
   - For `dim_line_h`: choose a northing `n` (e.g., `layout.ctr.n + layout.iwt3 + 1.0` for "1' north of IW7 north face")
   - For `dim_line_v`: choose an easting `e` (e.g., `layout.iw3.e + 2.0` for "2' east of IW3 east face")

3. **Add the call** in the dimension lines section (between wall rendering and openings):
   ```python
   # Description of what the dimension measures
   dim_line_h(out, start_e, n, end_e, fmt_dist(end_e - start_e))
   ```

4. **Prefixed labels** — For labeled dimensions (e.g., closets, storage), pass a formatted string:
   ```python
   dim_line_h(out, e1, n, e2, f"CLOSET {fmt_dist(e2 - e1)}")
   ```

### Example (real code)

```python
# F2-F3 east face to IW3 west face, 1' north of IW7 north face
dim_n = layout.ctr.n + layout.iwt3 + 1.0
dim_line_h(out, pts["W3"][0], dim_n, layout.iw3.w, fmt_dist(layout.iw3.w - pts["W3"][0]), to_svg)
```

---

## 2. Identifying Wall Faces and Coordinates

Understanding which coordinate to use for "east face of X wall" is the most common source of confusion.

### Perimeter walls (8" thick, F-series outer / W-series inner)

The outline traverses CW (as viewed from above): F1 → F2 → ... → F20 → F1. The interior is on the **right** side.

| Wall side of building | Outer (exterior) face | Inner (interior) face |
|-|-|-|
| **West** (F2-F3, F4-F5) | `pts["F<n>"][0]` (smaller easting) | `pts["W<n>"][0]` (larger easting) |
| **East** (F14-F15) | `pts["F<n>"][0]` (larger easting) | `pts["W<n>"][0]` (smaller easting) |
| **North** (F6-F7) | `pts["F<n>"][1]` (larger northing) | `pts["W<n>"][1]` (smaller northing) |
| **South** (F18-F19, F20-F1) | `pts["F<n>"][1]` (smaller northing) | `pts["W<n>"][1]` (larger northing) |

**Key insight:** For walls on the west side (like F2-F3), the "east face" is the **inner** face at `pts["W<n>"]`, not the F-series point.

### Interior walls

Interior wall positions are computed in `floorplan/layout.py` and returned in the `InteriorLayout` NamedTuple. Access them via `layout.<field>`:

| Wall | Type | West face | East face | South face | North face |
|-|-|-|-|-|-|
| **IW1** (horizontal, 6") | polygon | — | — | `layout.iw1_s` | `layout.iw1_n` |
| **IW2** (vertical, 6") | BBox | `layout.iw2.w` | `layout.iw2.e` | `layout.iw2.s` | `layout.iw2.n` |
| **IW3** (perpendicular to W20-W1, 4") | BBox+poly | `layout.iw3.w` | `layout.iw3.e` | `layout.iw3.s` | `layout.iw3.n` |
| **IW4** (vertical, 4") | scalars | `layout.iw4_w` | `layout.iw4_e` | `layout.wall_south_n` | `layout.iw1_s` |
| **IW5** (horizontal, 3") | BBox | `layout.iw5.w` | `layout.iw5.e` | `layout.iw5.s` | `layout.iw5.n` |
| **IW6** (horizontal, 1") | polygon | — | — | `layout.iw6_s` | `layout.iw6_n` |
| **IW7** (parallel to W20-W1, 3") | BBox+poly | `layout.iw7.w` | `layout.iw7.e` | `layout.iw7.s` | `layout.iw7.n` |
| **IW8** (horizontal, 6") | BBox | `layout.iw8.w` | `layout.iw8.e` | `layout.iw8.s` | `layout.iw8.n` |
| **IW9** (perpendicular to W20-W1, 3") | BBox+poly | `layout.iw9.w` | `layout.iw9.e` | `layout.iw9.s` | `layout.iw9.n` |
| **IW11** (N-S, 4") | BBox+poly | `layout.iw11.w` | `layout.iw11.e` | `layout.iw11.s` | `layout.iw11.n` |
| **IW12** (perpendicular to IW11, 4") | BBox+poly | `layout.iw12.w` | `layout.iw12.e` | `layout.iw12.s` | `layout.iw12.n` |
| **IW14** (parallel to IW12, 3") | BBox+poly | `layout.iw14.w` | `layout.iw14.e` | `layout.iw14.s` | `layout.iw14.n` |
| **IW15** (N-S, 4") | BBox | `layout.iw15.w` | `layout.iw15.e` | `layout.iw15.s` | `layout.iw15.n` |
| **IW16** (N-S, 4") | polygon | — | — | — | — |

BBox-type walls use `.w`, `.s`, `.e`, `.n` accessors. "BBox+poly" walls have both a `BBox` field and a `_poly` field (`list[Point]`) for the actual polygon (which may differ from the BBox on curved walls). IW1 and IW6 are pure polygon walls. IW4 uses individual scalar fields.

### Room-relative references

| Reference | Variable | Notes |
|-|-|-|
| Counter east edge | `layout.ctr.e` | |
| Counter north edge | `layout.ctr.n` | |
| Counter south edge | `layout.ctr.s` | Same as `pts["W1"][1]` |
| Bedroom center E-W | `(layout.iw3.e + layout.iw4_w) / 2` | |
| Inner south wall | `pts["W1"][1]` | |
| Inner west wall | `pts["W2"][0]` | |

---

## 3. Adding an Interior Wall

Interior wall positions are computed in `floorplan/layout.py` and rendered in `floorplan/gen_floorplan.py` (`_render_walls` function).

### Step-by-step

1. **Add a constant** for any new offset/dimension in `floorplan/constants.py`.

2. **Add fields** to the `InteriorLayout` NamedTuple in `floorplan/layout.py`:
   ```python
   mywall: BBox   # for rectangular walls
   # or
   mywall: list[Point]  # for L-shaped or irregular walls
   ```

3. **Compute the position** in `compute_interior_layout()` in `layout.py`:
   ```python
   mywall_w = <west face easting>
   mywall_e = mywall_w + WALL_3IN  # use named constants
   mywall_s = <south face northing>
   mywall_n = <north face northing>
   ```

4. **Render** in `_render_walls()` in `gen_floorplan.py`:
   ```python
   wall_poly(out, [(mywall.w, mywall.s), (mywall.e, mywall.s),
                   (mywall.e, mywall.n), (mywall.w, mywall.n)], to_svg)
   ```

5. **Add the polygon to `compute_iw_area()`** in `gen_floorplan.py` so its area is subtracted from the interior area.

6. **Add tests** in `tests/test_layout.py`.

### Wall thickness constants

Accessed via layout fields or from `floorplan/constants.py`:
- `WALL_3IN` = 3" (0.25') — also `layout.iwt3`
- `WALL_4IN` = 4" (0.333') — also `layout.iwt4`
- `WALL_6IN` = 6" (0.5') — also `layout.iwt`

---

## 4. Adding an Opening

**File:** `floorplan/openings.py` — single source of truth for all opening positions.

Openings are rendered as light-blue rectangles (`rgb(220,235,255)`) with `#4682B4` stroke. They cut through the wall from outer (F-series) to inner (W-series) face.

### Types

- `OuterOpening(name, seg_start, seg_end, poly)` — perimeter wall opening (O1-O11)
- `RoughOpening(name, bbox, wall_name, orientation)` — interior wall rough opening (RO1-RO5)
- `WallOpening(name, seg_idx, t_start, t_end)` — parametric form used by `walls/gen_walls.py`

### Adding a perimeter opening

1. Add the opening definition in `compute_outer_openings()` in `floorplan/openings.py`.
2. The function returns `OuterOpening` objects with 4-vertex polygons in survey coords.
3. Both `gen_floorplan.py` and `gen_walls.py` consume openings from this single source — no need to update two places.

### Adding a rough opening (interior wall)

1. Add the opening in `compute_rough_openings()` in `floorplan/openings.py`.
2. Update rendering in `_render_walls()` in `gen_floorplan.py` and `_render_interior_walls()` in `walls/gen_walls.py`.

### Tests

Opening tests are in `tests/test_gen_floorplan.py` and `tests/test_gen_walls.py` (11 outer openings, 5 rough openings, segment index validity, parametric ranges).

---

## 5. Adding an Appliance or Furniture Item

**File:** `floorplan/gen_floorplan.py` (appliances section)

### Rectangle items (washer, dryer, bed)

1. Compute bounding box in survey coords: `item_w`, `item_e`, `item_s`, `item_n`.
2. Convert to SVG and render:
   ```python
   sx1, sy1 = to_svg(item_w, item_n)   # SVG top-left = survey NW corner
   sx2, sy2 = to_svg(item_e, item_s)   # SVG bottom-right = survey SE corner
   sw = sx2 - sx1; sh = sy2 - sy1
   out.append(f'<rect x="{sx1:.1f}" y="{sy1:.1f}" width="{sw:.1f}" height="{sh:.1f}"'
              f' fill="rgba(100,150,200,0.2)" stroke="#4682B4" stroke-width="0.8"/>')
   ```
3. Add a centered label.

**Note:** SVG y-axis is inverted from northing. `to_svg(e, n)` handles this, but the NW corner (max northing) maps to the SVG top-left (min y).

### Circular items (water heater)

Use `<circle>` with `cx`, `cy` from `to_svg()` and radius converted via the scale factor.

---

## 6. Wall Construction Detail Drawing

**Files:** `walls/gen_walls.py`, `shared/wall_shells.py`, `floorplan/constants.py`

The wall detail drawing (`walls/walls.svg`) shows the double-shell 3D-printed concrete outer wall construction at 1:72 scale.

### Wall construction model

- **Outer shell**: 2" thick concrete (F-series to S-series boundary)
- **Air gap**: 4" between shells
- **Inner shell**: 2" thick concrete (G-series to W-series boundary)
- **Total**: 8" (`WALL_OUTER`)

Four concentric boundary paths trace the building perimeter:

| Path | Point series | Inset from F | Description |
|-|-|-|-|
| Outer face of outer shell | F-series | 0" | Existing `outline_segs` |
| Inner face of outer shell | S-series | 2" | `compute_inset_path(..., SHELL_THICKNESS, "S")` |
| Outer face of inner shell | G-series | 6" | `compute_inset_path(..., SHELL_THICKNESS + AIR_GAP, "G")` |
| Inner face of inner shell | W-series | 8" | Existing `inner_segs` |

### Construction constants

Defined in `floorplan/constants.py` (re-exported by `walls/constants.py`):

- `SHELL_THICKNESS` = 2/12 ft (2")
- `AIR_GAP` = 4/12 ft (4")
- `OPENING_INSIDE_RADIUS` = 10/304.8 ft (10mm)

### Opening U-turn corners

At each opening boundary, the shells connect via 90-degree corner turns:

- **Inside radius** (`R_in`): `OPENING_INSIDE_RADIUS` (10mm)
- **Outside radius** (`R_out`): `R_in + SHELL_THICKNESS`
- The turned outside face is flush with the opening boundary
- `uturn_polygon()` in `shared/wall_shells.py` builds the U-turn as a single closed polygon using quarter-circle arcs

### Modifying wall constants

1. Edit values in `floorplan/constants.py`
2. Run `python walls/gen_walls.py` to regenerate
3. Run `python -m pytest tests/test_gen_walls.py` to verify

### Adding/modifying openings

Opening positions are defined in `floorplan/openings.py` (single source of truth). `walls/gen_walls.py` imports openings via `outer_to_wall_openings()` which converts them to parametric `WallOpening` objects. Each opening maps to a parametric range `(t_start, t_end)` along its outline segment. See Section 4 for details.

---

## 7. Verifying Changes

After any geometry or layout change, regenerate and inspect all SVGs:

```bash
python gen_all.py
```

This captures `git describe --always --dirty=-DEV` once into `.git_describe`, runs all six generators (survey, floorplan, walls, and three span scripts) using that cached value, then deletes the cache. This ensures all title blocks embed the same version string even though writing the first SVG dirties the working tree.

Individual scripts can also be run standalone — they fall back to a live `git describe` if the cache file is absent:

```bash
python survey/gen_path_svg.py
python floorplan/gen_floorplan.py
python walls/gen_walls.py
python span/gen_span.py
python span/gen_span_minmax.py
python span/gen_span_min.py
```

The floorplan script prints:
- Outer/inner/wall areas (sanity check for area subtraction)
- All F-series and W-series point coordinates (verify geometry)

The walls script prints:
- Shell and gap dimensions
- Opening corner radius

Open `floorplan/floorplan.svg`, `survey/path_area.svg`, `walls/walls.svg`, and `span/span.svg` to visually inspect.

---

## 8. Adding Text Captions to site_plan.pdf

**File:** `site/gen_site_plan.py`

The site plan overlays the building outline on a copy of `site/site_survey.pdf` using pymupdf (fitz). Adding rotated text captions has several pitfalls.

### Use `page.insert_text` with `morph`, not `TextWriter`

`TextWriter.write_text()` with `morph` produces incorrect positions on pages copied from other PDFs (the text is displaced hundreds of points from the intended location). Use `page.insert_text()` with its `morph` parameter instead — this works correctly on copied pages.

```python
page.insert_text(
    fitz.Point(x, y),
    text, fontname="helv", fontsize=8, color=(0, 0, 0),
    morph=(fitz.Point(pivot_x, pivot_y), fitz.Matrix(angle_degrees)))
```

### Morph rotation convention

`fitz.Matrix(angle)` rotates anti-clockwise in the visual (on-screen) sense. On the survey-copied page:

- Positive angle = visually anti-clockwise (text tilts upward to the right)
- Negative angle = visually clockwise (text tilts downward to the right)

To match the survey's perpendicular-distance captions (46.7', 39.5', etc.), compute the perpendicular direction angle and **negate** it:

```python
perp_deg = math.degrees(math.atan2(dy_pdf, dx_pdf))  # PDF coords (y-down)
morph=(pivot, fitz.Matrix(-perp_deg))                 # negate for correct visual rotation
```

### Centering text at a point

To center text horizontally at `(cx, cy)`:

```python
tw = fitz.get_text_length(text, fontname="helv", fontsize=fs)
page.insert_text(
    fitz.Point(cx - tw / 2.0, cy + fs / 3.0),  # offset left by half width, down by ~1/3 font size
    text, ...,
    morph=(fitz.Point(cx, cy), fitz.Matrix(angle)))
```

The morph pivot should be at the intended center `(cx, cy)`. The text origin is offset by half the text width; the morph rotation keeps the center fixed.

### Coordinate system

All positions in `gen_site_plan.py` use PDF coordinates (x-right, y-down from top-left). The survey scale is 2.4 PDF points per foot (72 pts/inch ÷ 30 ft/inch). Building coordinates are transformed via `building_to_pdf()`.

---

## 9. Contributing to This Document

**For future agents:** If you encounter a complex task that required significant codebase research and is not already covered here, please add a new section documenting the procedure. Follow the existing format:

1. Add an entry to the Table of Contents.
2. Write a section with:
   - Which file(s) to edit
   - The relevant helper functions or patterns
   - Step-by-step instructions
   - A concrete code example from the codebase
3. Keep instructions concise and focused on the "how", not the "why" — the architecture is documented in `CLAUDE.md`.
