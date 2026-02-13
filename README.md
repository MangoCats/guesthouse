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

Generate the survey path SVG (outline + inset path with labels):

```
python survey/gen_path_svg.py
```

Output: `survey/path_area.svg`

Generate the detailed floorplan SVG (walls, rooms, appliances, dimensions):

```
python floorplan/gen_floorplan.py
```

Output: `floorplan/floorplan.svg`

## Tests

```
pytest
```

## Project Structure

```
shared/              Common types, geometry, survey computation, SVG utilities
floorplan/           Building design: geometry, layout, constants, SVG renderer
survey/              Survey scripts and data
tests/               Unit tests
```

See `CLAUDE.md` for detailed file descriptions, dependency graph, and naming conventions.
