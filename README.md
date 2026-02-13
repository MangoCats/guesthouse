# Hut2

Curved-wall building outline geometry and floorplan SVG generation.

## Setup

Requires Python 3.12+.

```
pip install -e .
```

This installs the `shared` and `floorplan` packages in editable mode so all scripts can import them.

To also install test dependencies:

```
pip install -e ".[test]"
```

The `survey/adjust_pentagon.py` script requires numpy and scipy:

```
pip install -e ".[adjust]"
```

## Usage

Regenerate all SVG files:

```
python gen_all.py
```

This captures `git describe` once before generating, so all title blocks show the same version even though writing the first SVG dirties the working tree.

Individual scripts can also be run standalone (they fall back to a live `git describe`):

```
python survey/gen_path_svg.py        # → survey/path_area.svg
python floorplan/gen_floorplan.py    # → floorplan/floorplan.svg
python walls/gen_walls.py            # → walls/walls.svg, walls/all_walls.svg
```

## Tests

```
pytest
```

## Project Structure

```
shared/              Common types, geometry, survey computation, SVG utilities
floorplan/           Building design: geometry, layout, constants, SVG renderer
walls/               Outer wall construction detail drawing
survey/              Survey scripts and data
tests/               Unit tests
```

See `CLAUDE.md` for detailed file descriptions, dependency graph, and naming conventions.
