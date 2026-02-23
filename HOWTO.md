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
9. [Getting Rotation Direction Right](#9-getting-rotation-direction-right)
10. [Changing Exterior Wall Thickness](#10-changing-exterior-wall-thickness)
11. [Contributing to This Document](#11-contributing-to-this-document)
12. [Rotation-Invariant Placement](#12-rotation-invariant-placement)

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
   - Interior wall faces: use `layout.<field>` (e.g., `layout.iw3.w`, `layout.iw3.e`, `layout.iw4.w`, `layout.iw1.s`)

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

### Perimeter walls (WALL_OUTER thick, F-series outer / W-series inner)

The outline traverses CW (as viewed from above): F1 → F2 → ... → F20 → F1. The interior is on the **right** side.

| Wall side of building | Outer (exterior) face | Inner (interior) face |
|-|-|-|
| **West** (F2-F3, F4-F5) | `pts["F<n>"][0]` (smaller easting) | `pts["W<n>"][0]` (larger easting) |
| **East** (F14-F15) | `pts["F<n>"][0]` (larger easting) | `pts["W<n>"][0]` (smaller easting) |
| **North** (F6-F7) | `pts["F<n>"][1]` (larger northing) | `pts["W<n>"][1]` (smaller northing) |
| **South** (F18-F19, F20-F1) | `pts["F<n>"][1]` (smaller northing) | `pts["W<n>"][1]` (larger northing) |

**Key insight:** For walls on the west side (like F2-F3), the "east face" is the **inner** face at `pts["W<n>"]`, not the F-series point.

### Interior walls

Interior wall positions are computed in `floorplan/layout.py` and returned in the `InteriorLayout` NamedTuple. Every wall, appliance, and furniture item is a `Wall` NamedTuple (defined in `shared/types.py`) with both polygon corners and axis-aligned bounding box:

```python
class Wall(NamedTuple):
    poly: list[Point]  # [SW, SE, NE, NW] corners
    w: float; s: float; e: float; n: float  # axis-aligned bounding box
```

Access via `layout.<field>.w`, `.s`, `.e`, `.n` for bounding box, or `layout.<field>.poly` for the corner polygon.

| Wall | Thickness | West face | East face | South face | North face | Polygon |
|-|-|-|-|-|-|-|
| **IW1** (horizontal) | 6" | `layout.iw1.w` | `layout.iw1.e` | `layout.iw1.s` | `layout.iw1.n` | `layout.iw1.poly` |
| **IW2** (vertical) | 6" | `layout.iw2.w` | `layout.iw2.e` | `layout.iw2.s` | `layout.iw2.n` | `layout.iw2.poly` |
| **IW3** (perp. to W20-W1) | 4" | `layout.iw3.w` | `layout.iw3.e` | `layout.iw3.s` | `layout.iw3.n` | `layout.iw3.poly` |
| **IW4** (vertical) | 4" | `layout.iw4.w` | `layout.iw4.e` | `layout.iw4.s` | `layout.iw4.n` | `layout.iw4.poly` |
| **IW5** (horizontal) | 3" | `layout.iw5.w` | `layout.iw5.e` | `layout.iw5.s` | `layout.iw5.n` | `layout.iw5.poly` |
| **IW6** (horizontal) | 1" | `layout.iw6.w` | `layout.iw6.e` | `layout.iw6.s` | `layout.iw6.n` | `layout.iw6.poly` |
| **IW7** (par. to W20-W1) | 3" | `layout.iw7.w` | `layout.iw7.e` | `layout.iw7.s` | `layout.iw7.n` | `layout.iw7.poly` |
| **IW8** (horizontal) | 6" | `layout.iw8.w` | `layout.iw8.e` | `layout.iw8.s` | `layout.iw8.n` | `layout.iw8.poly` |
| **IW9** (perp. to W20-W1) | 3" | `layout.iw9.w` | `layout.iw9.e` | `layout.iw9.s` | `layout.iw9.n` | `layout.iw9.poly` |
| **IW11** (N-S) | 4" | `layout.iw11.w` | `layout.iw11.e` | `layout.iw11.s` | `layout.iw11.n` | `layout.iw11.poly` |
| **IW12** (perp. to IW11) | 4" | `layout.iw12.w` | `layout.iw12.e` | `layout.iw12.s` | `layout.iw12.n` | `layout.iw12.poly` |
| **IW14** (par. to IW12) | 3" | `layout.iw14.w` | `layout.iw14.e` | `layout.iw14.s` | `layout.iw14.n` | `layout.iw14.poly` |
| **IW15** (N-S) | 4" | `layout.iw15.w` | `layout.iw15.e` | `layout.iw15.s` | `layout.iw15.n` | `layout.iw15.poly` |
| **IW16** (N-S) | 4" | `layout.iw16.w` | `layout.iw16.e` | `layout.iw16.s` | `layout.iw16.n` | `layout.iw16.poly` |

Appliances and furniture also use `Wall`: `layout.dryer`, `layout.washer`, `layout.ctr`, `layout.dresser`, `layout.bed`.

### Room-relative references

| Reference | Variable | Notes |
|-|-|-|
| Counter east edge | `layout.ctr.e` | |
| Counter north edge | `layout.ctr.n` | |
| Counter south edge | `layout.ctr.s` | Same as `pts["W1"][1]` |
| Bedroom center E-W | `(layout.iw3.e + layout.iw4.w) / 2` | |
| Inner south wall | `pts["W1"][1]` | |
| Inner west wall | `pts["W2"][0]` | |

---

## 3. Adding an Interior Wall

Interior wall positions are computed in `floorplan/layout.py` and rendered in `floorplan/gen_floorplan.py` (`_render_walls` function).

### Step-by-step

1. **Add a constant** for any new offset/dimension in `floorplan/constants.py`.

2. **Add a field** to the `InteriorLayout` NamedTuple in `floorplan/layout.py`:
   ```python
   mywall: Wall
   ```

3. **Compute the position** in `compute_interior_layout()` in `layout.py` using `_make_wall`:
   ```python
   mywall_sw = <southwest corner>
   mywall_se = <southeast corner>
   mywall_ne = <northeast corner>
   mywall_nw = <northwest corner>
   # _make_wall computes the bounding box automatically
   ```
   Then in the return statement: `mywall=_make_wall([mywall_sw, mywall_se, mywall_ne, mywall_nw])`

4. **Render** in `_render_walls()` in `gen_floorplan.py`:
   ```python
   wall_poly(out, layout.mywall.poly, to_svg)
   ```

5. **Add the polygon to `compute_iw_area()`** in `gen_floorplan.py` so its area is subtracted from the interior area.

6. **Add tests** in `tests/test_layout.py`.

### Wall thickness constants

Import from `floorplan/constants.py`:
- `WALL_3IN` = 3" (0.25')
- `WALL_4IN` = 4" (1/3')
- `WALL_6IN` = 6" (0.5')

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

**Important:** All positions must be defined relative to wall segments using direction vectors, never using raw coordinate indexing (`pts[...][0]`, `pts[...][1]`) or hardcoded angles. See [Section 12](#12-rotation-invariant-placement) for the full pattern.

### Rectangle items (washer, dryer, bed)

1. Compute the item's anchor point using `seg_vecs()` + `offset_pt()` relative to the nearest wall segment or interior wall polygon face.
2. Build the bounding box from the anchor using direction vectors:
   ```python
   item_nw = offset_pt(anchor, depth, wall_inward)
   item_se = offset_pt(anchor, -width, wall_along)
   ```
3. Convert to SVG and render:
   ```python
   sx1, sy1 = to_svg(item_nw[0], item_nw[1])
   sx2, sy2 = to_svg(item_se[0], item_se[1])
   sw = sx2 - sx1; sh = sy2 - sy1
   out.append(f'<rect x="{sx1:.1f}" y="{sy1:.1f}" width="{sw:.1f}" height="{sh:.1f}"'
              f' fill="rgba(100,150,200,0.2)" stroke="#4682B4" stroke-width="0.8"/>')
   ```
4. Add a centered label.

**Note:** SVG y-axis is inverted from northing. `to_svg(e, n)` handles this, but the NW corner (max northing) maps to the SVG top-left (min y).

### Circular items (water heater)

Use `<circle>` with `cx`, `cy` from `to_svg()` and radius converted via the scale factor.

### Rotated items (chair, loveseat, rocker)

Compute the SVG rotation angle from a wall direction vector using `_svg_angle()`, then apply an SVG `rotate()` transform. See [Section 12](#12-rotation-invariant-placement) for details.

---

## 6. Wall Construction Detail Drawing

**Files:** `walls/gen_walls.py`, `shared/wall_shells.py`, `floorplan/constants.py`

The wall detail drawing (`walls/walls.svg`) shows the double-shell 3D-printed concrete outer wall construction at 1:72 scale.

### Wall construction model

- **Outer shell**: `SHELL_THICKNESS` (2") concrete (F-series to S-series boundary)
- **Air gap**: `AIR_GAP` (computed from `WALL_OUTER - 2*SHELL_THICKNESS`)
- **Inner shell**: `SHELL_THICKNESS` (2") concrete (G-series to W-series boundary)
- **Total**: `WALL_OUTER` (adjustable 8"–12", see [Section 10](#10-changing-exterior-wall-thickness))

Four concentric boundary paths trace the building perimeter:

| Path | Point series | Inset from F | Description |
|-|-|-|-|
| Outer face of outer shell | F-series | 0 | Existing `outline_segs` |
| Inner face of outer shell | S-series | `SHELL_THICKNESS` | `compute_inset_path(..., SHELL_THICKNESS, "S")` |
| Outer face of inner shell | G-series | `SHELL_THICKNESS + AIR_GAP` | `compute_inset_path(..., SHELL_THICKNESS + AIR_GAP, "G")` |
| Inner face of inner shell | W-series | `WALL_OUTER` | Existing `inner_segs` |

### Construction constants

Defined in `floorplan/constants.py` (re-exported by `walls/constants.py`):

- `SHELL_THICKNESS` = 2/12 ft (2") — constant regardless of wall thickness
- `AIR_GAP` = `WALL_OUTER - 2*SHELL_THICKNESS` — computed, increases with wall thickness
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

## 9. Getting Rotation Direction Right

Rotation direction errors are common when adding door swings, casement windows, or any arc that should swing "inward" vs "outward" relative to the building. The issue is that the math is correct in isolation but the sign is backwards for the intended physical direction.

### The pattern

A typical swing illustration has:

- **Hinge point** — fixed corner
- **Closed direction** — unit vector from hinge toward the closed-end corner
- **Cross direction** — unit vector perpendicular to the wall, pointing inward or outward
- **Sweep angle** — how far to rotate (90° for doors, 45° for casement windows)

The rotation angle sign determines whether the free end swings inward or outward.

### How to determine the sign

1. Compute the **closed direction** (`c_dir`) and the **target direction** (the direction you want the swing to go — e.g., inward for doors, outward for casements).

2. Take the 2D cross product to find the rotation sense:
   ```python
   cross = c_dir[0] * target_dir[1] - c_dir[1] * target_dir[0]
   ```

3. **Positive cross** → the target is CCW from the closed direction → use `rot_sign = +1`.
   **Negative cross** → the target is CW → use `rot_sign = -1`.

4. Apply the rotation:
   ```python
   angle = rot_sign * sweep_angle
   open_dir = (c_dir[0] * cos(angle) - c_dir[1] * sin(angle),
               c_dir[0] * sin(angle) + c_dir[1] * cos(angle))
   ```

### Common mistake: confusing inward and outward

The "inward direction" is computed as the vector from the outer face midpoint to the inner face midpoint:

```python
outer_mid = ((poly[0][0] + poly[1][0]) / 2, (poly[0][1] + poly[1][1]) / 2)
inner_mid = ((poly[2][0] + poly[3][0]) / 2, (poly[2][1] + poly[3][1]) / 2)
inward_dir = normalize(inner_mid - outer_mid)
```

- To swing **inward** (into the building, e.g., interior doors): use `inward_dir` as the target direction.
- To swing **outward** (away from the building, e.g., casement windows): **negate the sign**. Either use `-inward_dir` as the target, or flip the `rot_sign`:
  ```python
  rot_sign = -1 if cross > 0 else 1   # outward = opposite of inward
  ```

### Quick sanity check

Before committing, verify the swing direction by reasoning about one concrete case. For example, on the east wall (F14-F15):

- Closed direction is **north** (along the wall)
- Inward direction is **west** (toward building center)
- Outward direction is **east** (away from building)
- A casement that should swing east: the free end should rotate CW from north → rot_sign should be **-1**

If the cross product of (north, east) = `(0)(0) - (1)(1) = -1` (negative), then the outward formula `rot_sign = -1 if cross > 0 else 1` gives `+1`. But wait — that rotates CCW from north, which goes west (inward), not east. So for outward you need the opposite: `rot_sign = -1 if cross > 0 else 1` applied to the **inward** cross gives the correct outward rotation. The key insight: compute the cross product against the **inward** direction, then negate.

### Existing examples

- **Interior doors** (RO3, RO4, RO5 in `_render_openings`): swing inward, 90° arcs, use `JAMB_COLOR`.
- **Casement windows** (O8, O9, O10 in `_render_openings`): swing outward, 45° arcs, use `OPENING_STROKE` color. The rotation sign is flipped relative to the inward cross product.

---

## 10. Changing Exterior Wall Thickness

**File:** `floorplan/constants.py` (line 10)

The exterior wall thickness is parameterized. To change it, edit the single constant `WALL_OUTER` (valid range: 8"–12"). All dependent geometry updates automatically.

### What to change

Edit one line in `floorplan/constants.py`:

```python
WALL_OUTER = 10.0 / 12.0          # 10" outer wall (adjustable: 8"-12")
```

Change `10.0` to the desired thickness in inches (e.g., `9.0` for 9").

### What happens automatically

The private constant `_WE = WALL_OUTER - 8.0 / 12.0` (wall extra beyond the 8" baseline) drives 11 derived constants:

| Constant | Formula | Why |
|-|-|-|
| `AIR_GAP` | `WALL_OUTER - 2*SHELL_THICKNESS` | Shell thickness is fixed; gap absorbs the change |
| `CORNER_NE_R` | `10"/12 + _WE` | Arc radii grow to keep arc centers (and W-series) fixed |
| `CORNER_NW_R` | `28"/12 + _WE` | Same |
| `UPPER_E_R` | `28"/12 + _WE` | Same |
| `ARC_180_R` | `28"/12 + _WE` | Same |
| `ARC_F3_R` | `base + _WE` | Same |
| `ARC_F13_R` | `base + _WE` | Same |
| `ARC_F19_R` | `base + _WE` | Same |
| `F6_HEIGHT` | `base + 2*_WE` | F1 moves south and F6 moves north, so height changes by 2x |
| `SOUTH_WALL_N` | `base - _WE` | South wall face moves south |
| `F15_OFFSET_E` | `base + _WE` | East wall moves east |

Additionally, anchor offsets (`_WE` applied to Pi2, Pi3, PiX, Pi5) in `gen_floorplan.py`, `gen_path_svg.py`, and `tests/conftest.py` shift the outline outward while keeping W-series (inner wall) coordinates constant.

### Key invariant

**W-series points do not move** when wall thickness changes. The inner wall geometry is fixed by the survey inset path. Changing `WALL_OUTER` only moves F-series points outward (or inward), preserving all interior dimensions.

### After changing

1. Run tests: `python -m pytest tests/ -x`
2. **Update regression test expectations** in `tests/test_outline.py` (`_EXPECTED_F` coordinates, area) and `tests/test_gen_walls.py` (wall polygon count) to match the new geometry.
3. Regenerate SVGs: `python gen_all.py`

### Wall polygon count

The wall polygon count (tested in `test_gen_walls.py::test_wall_polygon_count`) may change at different thicknesses. Wider walls increase the U-turn trim distance (`R_out / seg_len`), which can eliminate small solid wall ranges between openings. To investigate, count polygons per segment:

```python
from shared.wall_shells import openings_on_seg, solid_ranges
# For each LineSeg with openings, check if adjusted ranges survive:
# t_e_adj > t_s_adj + 1e-9 after delta_t = R_out / seg_len trim
```

---

## 12. Rotation-Invariant Placement

All positions in the floorplan are defined relative to wall segments using direction vectors. This ensures geometry remains correct if the building outline is rotated or if wall angles change.

### Core functions (`shared/geometry.py`)

```python
def seg_vecs(p1, p2):
    """Along-direction and CW-inward normal for segment p1->p2."""
    dx, dy = p2[0] - p1[0], p2[1] - p1[1]
    length = math.sqrt(dx*dx + dy*dy)
    along = (dx/length, dy/length)
    inward = (dy/length, -dx/length)    # CW-right perpendicular
    return along, inward

def offset_pt(origin, dist, direction):
    """Offset point by dist along direction vector."""
    return (origin[0] + dist*direction[0], origin[1] + dist*direction[1])
```

- `seg_vecs(p1, p2)` returns two unit vectors: `along` (direction from p1 to p2) and `inward` (CW-right perpendicular, which points into the building for the CW outline traversal).
- `offset_pt(origin, dist, direction)` moves a point by `dist` feet along a direction vector.
- `line_isect(p1, d1, p2, d2)` finds the intersection of two lines defined by point + direction.

### Positioning pattern

**Step 1: Get wall direction vectors**

```python
# From W-series perimeter wall points:
w9w10_al, w9w10_in = seg_vecs(pts["W9"], pts["W10"])   # north wall

# From interior wall polygon faces:
_iw4_w_al, _iw4_w_out = seg_vecs(layout.iw4.poly[3], layout.iw4.poly[0])  # IW4 west face
```

Interior wall polygons are ordered `[SW, SE, NE, NW]`. Face vectors:
- **West face**: `seg_vecs(poly[3], poly[0])` — NW to SW, outward = west
- **East face**: `seg_vecs(poly[1], poly[2])` — SE to NE, outward = east
- **North face**: `seg_vecs(poly[3], poly[2])` — NW to NE, CW = south; negate for outward = north
- **South face**: `seg_vecs(poly[0], poly[1])` — SW to SE, CW = north; negate for outward = south

**Step 2: Position using offsets**

```python
# Single-wall anchor: offset from a wall point
hamper_pos = offset_pt(pts["W2"], 2.0/12.0, w2w3_in)    # 2" inward from W2

# Two-wall anchor: offset from wall intersection corner
corner = line_isect(layout.iw4.poly[3], _iw4_w_al,
                    layout.iw1.poly[3], _iw1_n_al)
item_pos = offset_pt(offset_pt(corner, d_w, _iw4_w_out), d_n, _iw1_n_out)
```

**Step 3: Compute SVG rotation from wall direction**

```python
def _svg_angle(along):
    """SVG rotation angle (degrees, CW) for a direction vector."""
    return -math.degrees(math.atan2(along[1], along[0]))

desk_angle = _svg_angle(w16w17_al)   # desk aligns with W16-W17 wall
```

### Local anchor helpers

Each render function defines local helpers for its reference walls:

- `_render_kitchen` uses `_iwp(d_along, d_inward)` — offset from IW1/W9 intersection along the north wall.
- `_render_furniture` uses `_lwp(d_w, d_n)` — offset from IW4/IW1 corner along their outward directions; and `_nwp(d_along, d_inward)` — offset from W9 along the north wall.
- `_render_appliances` uses direct `offset_pt` calls from wall polygon faces.

### What NOT to do

- **Raw coordinate indexing**: `pts["W9"][0]`, `pts["W9"][1]` — fragile if walls are not axis-aligned.
- **Hardcoded angles**: `rotate(30, ...)` — derive from `_svg_angle(wall_along)` instead.
- **BBox face references for positioning**: `layout.iw5.w`, `layout.iw5.n` — use `seg_vecs` on the polygon face instead. This applies to ALL position definitions including dimension line endpoints.
- **Cardinal direction arithmetic**: `item_e = wall_e + gap` — use `offset_pt(wall_pt, gap, wall_outward)`.

### Existing examples

| Item | Reference wall | Pattern |
|------|---------------|---------|
| Kitchen appliances | W9-W10 (north wall) | `_iwp(d_along, d_inward)` from IW1/W9 corner |
| Sofa, loveseat | IW4 west + IW1 north | `_lwp(d_w, d_n)` from IW4/IW1 corner |
| Desk | W16-W17 | `_svg_angle(w16w17_al)` for rotation |
| Chair | W11-W12 chord | `_svg_angle(w12w13_al) - 45` for rotation |
| Water heater | IW2 east face | Line-circle intersection along `_iw2_e_al` |
| Washer/dryer | W2-W3 + W9-W10 | `offset_pt` along wall vectors |

---

## 11. Contributing to This Document

**For future agents:** If you encounter a complex task that required significant codebase research and is not already covered here, please add a new section documenting the procedure. Follow the existing format:

1. Add an entry to the Table of Contents.
2. Write a section with:
   - Which file(s) to edit
   - The relevant helper functions or patterns
   - Step-by-step instructions
   - A concrete code example from the codebase
3. Keep instructions concise and focused on the "how", not the "why" — the architecture is documented in `CLAUDE.md`.
