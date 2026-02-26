# opt_rot: C17 Sweep Angle Analysis

## What this does

`gen_opt_rot.py` generates SVG drawings of the F-series building outline at
different C17 arc sweep angles, overlaid on the survey outer path and 6" inset
path. Each SVG shows how the outline shape and alignment change as C17 varies.

## Background

The F-series outline has one shape degree of freedom when:
- C15 and C17 sweep angles are free (constrained only by C15 + C17 = 90°)
- The F14-F15 and F18-F1 line lengths are free
- Rigid-body placement (rotation + translation) is free

All other segment lengths, radii, and sweep angles are locked.

### Closure constraint

Bearing closure forces C15 + C17 = π/2, so the single shape parameter is C17.
Given a C17 sweep angle, the two free line lengths are determined analytically:

    d_F14_F15 = d_F2_F5 + dN14 + R_a1 - R - L * cos(π/2 - C17)
    d_F18_F1  = dE14 - R - L * sin(π/2 - C17) - R_a1

where dE14/dN14 are the offsets from walking the locked chain F5→F14, R is the
C15/C17 arc radius (1.8087'), L is the F16-F17 line length (5'), and R_a1 is
the F1-F2 corner radius.

### Alignment constraints

The building outline must satisfy 4 alignment constraints from the survey:

1. **T3 on F18-F1**: survey point T3 lies on the F18→F1 line
2. **F12-F13 tangent to TC1 arc**: the signed distance from TC1 to line
   F12→F13 equals the TC1 arc radius R1
3. **F16 on Pi5-PiX line**: outline point F16 lies on the Pi5→PiX line
4. **F17 on Pi5-PiX line**: outline point F17 lies on the Pi5→PiX line

With 4 unknowns (rotation, dx, dy, C17) and 4 constraints, the system is
fully constrained. The exact solution is C17 ≈ 36.0805° — the current design.

### What the sweep shows

At angles other than the exact solution, the script finds the best-fit
rigid-body placement (3 unknowns, 4 constraints → least squares). The RMS
residual measures how far each angle is from satisfying all constraints
simultaneously. The exact solution has RMS ≈ 0.

## Output files

Each analysis produces `opt_rot/path_area_{label}.svg` containing:
- Title bar with C17/C15 angles and RMS residual
- Subtitle with derived lengths, area, rotation, and translation
- Individual constraint residuals (green if RMS < 0.01', red otherwise)
- Outer survey path at 20% opacity
- Inset path at 20% opacity
- F-series outline at 100% with point labels
- Enclosed area label centered on the drawing
- North arrow and legend

### Current analyses

| Label | C17 angle | Source |
|-------|-----------|--------|
| `300`–`420` | 30.0°–42.0° in 1° steps | Integer degree sweep |
| `712` | atan(7/12) ≈ 30.2564° | Named analysis |
| `812` | atan(8/12) ≈ 33.6901° | Named analysis |

## How to add a new analysis

Edit `gen_opt_rot.py` and add one line to the `angles` list (around line 302):

```python
angles = [(x / 10.0, f"{x:03d}") for x in range(300, 421, 10)]
angles.append((math.degrees(math.atan(7.0 / 12.0)), "712"))
angles.append((math.degrees(math.atan(8.0 / 12.0)), "812"))
angles.append((YOUR_ANGLE_IN_DEGREES, "YOUR_LABEL"))        # <-- add here
```

**Parameters:**
- First element: the C17 sweep angle in degrees (float)
- Second element: a filename-safe label string. The output file will be
  `path_area_{label}.svg`

**For named atan angles**, also add an entry to the `_atan_labels` dict
(around line 360) so the SVG title shows the symbolic form:

```python
_atan_labels = {"712": "atan(7/12)", "812": "atan(8/12)", "NEW": "atan(a/b)"}
```

Then regenerate:

```
python opt_rot/gen_opt_rot.py
```

### Constraints on valid angles

- Both d_F14_F15 and d_F18_F1 must be positive. Angles that produce a
  non-positive length are skipped automatically with a "SKIPPED" message.
- The valid range is approximately 26°–50° for C17.

## Dependencies

- `scipy` (for `least_squares` and `fsolve`): `pip install -e '.[adjust]'`
- All project packages (`shared/`, `floorplan/`, `survey/`)

## Running

```
python opt_rot/gen_opt_rot.py
```

Prints a table of results to stdout and writes SVG files to `opt_rot/`.
