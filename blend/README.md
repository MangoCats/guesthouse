# blend/ — Blender Model Generation

`gen_blend.py` builds native Blender (`.blend`) models of the
currently-loaded editor model (`app/adu.db`).

## Pipeline

1. **Launcher (plain Python):** builds `GeneratorData` from the live DB,
   reads `config_name` and `roof_style` from the `config` table, and
   generates SCAD text fresh for each roof variant (restoring the on-disk
   `scad/*.scad` files afterward so the working tree stays clean).
2. **Manifest:** hands headless Blender a manifest of models + survey parcel
   + building footprint + walk-around camera path.  The launcher re-launches
   itself under Blender if invoked with plain `python` (`$BLENDER_EXE`
   overrides the auto-detected `blender.exe`).  Invoking Blender directly
   with no manifest falls back to legacy `flat_roof` + `2in12` builds from
   the committed SCAD files.
3. **Mesh build (inside Blender):** reuses `scad/gen_gltf.py`'s SCAD parser
   and mesh builders so geometry matches the `.glb` models exactly, except:
   - the roof is rebuilt natively as a planar n-gon prism (gen_gltf's
     ear-clipping leaves a triangular gap on the repeated closing vertex);
   - interior walls are built per-wall with their door openings boolean-cut
     to door height (the GLB build omits them);
   - window glazing gets a very transparent material (alpha 0.18, EEVEE
     Next BLENDED).

## Roof variants / outputs

One `.blend` per roof style: the configured `roof_style` under the plain
config name, plus each of the other styles suffixed — for a config named
`Mark10` with a flat roof:

| File | Roof |
|------|------|
| `Mark10.blend` | configured style (e.g. flat) |
| `Mark10_2in12.blend` | 2:12 shed |
| `Mark10_1in12.blend` | 1:12 shed |
| `Mark10_split2.blend` | split: west 2:12 + tilted east plane |
| `Mark10_split1.blend` | split: west 1:12 + tilted east plane (7'6" closet eave) |

## Scene extras

- **Ground plane:** green grass covering the survey parcel (site-plan
  property extent), computed by inverting `site/gen_site_plan.py`'s
  building→PDF transform (needs PyMuPDF), with the building exterior outline
  cut out so no lawn shows under the house.
- **Boundary fences:** post-and-rail — dark-brown 3-rail (4'8") along the
  west property boundary (275.08' line), lighter-brown single-rail (4'0")
  along the south boundary (216.73' line), meeting at the SW corner, 4"
  posts every 10'.  (Survey plat is N-left/E-top, not N-up — see
  `site/gen_site_plan.py`.)
- **Camera animation:** the 2:12 model gets a 1800-frame walk-around
  animation (30 fps): camera circles the building clockwise from the 42"
  door side along an outline-offset path, 5'6" above ground, kept pointed at
  the origin (5'6") by a Track-To constraint, vertical FOV 105°.
  (`MarkZ_2in12_walkaround.mp4` is a rendered example.)

## Running

```
python blend/gen_blend.py        # or via python gen_all.py
```

Requires Blender (auto-detected or `$BLENDER_EXE`) and PyMuPDF.
