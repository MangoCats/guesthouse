"""SQLite database layer for the ADU building definition.

Stores all building constants, outline chain definitions, and entity metadata
so the building can be edited interactively and regenerated from the database.
"""
import ast
import os
import re
import sqlite3
from contextlib import contextmanager

_DIR = os.path.dirname(os.path.abspath(__file__))
_PROJECT = os.path.dirname(_DIR)
DB_PATH = os.path.join(_PROJECT, "app", "adu.db")

# ---------------------------------------------------------------------------
# Connection helper
# ---------------------------------------------------------------------------

@contextmanager
def get_db(db_path=None):
    """Yield an SQLite connection with WAL mode and row factory."""
    conn = sqlite3.connect(db_path or DB_PATH)
    conn.row_factory = sqlite3.Row
    conn.execute("PRAGMA journal_mode=WAL")
    conn.execute("PRAGMA foreign_keys=ON")
    try:
        yield conn
        conn.commit()
    except Exception:
        conn.rollback()
        raise
    finally:
        conn.close()


# ---------------------------------------------------------------------------
# Schema
# ---------------------------------------------------------------------------

_SCHEMA = """\
CREATE TABLE IF NOT EXISTS constants (
    name        TEXT PRIMARY KEY,
    value       REAL NOT NULL,
    expr        TEXT,               -- original Python expression
    unit        TEXT DEFAULT 'ft',
    category    TEXT DEFAULT 'misc',
    description TEXT DEFAULT ''
);

CREATE TABLE IF NOT EXISTS outline_chain (
    seq         INTEGER PRIMARY KEY,
    seg_type    TEXT NOT NULL,       -- 'L', 'CW', 'CCW'
    distance    REAL,               -- for lines
    radius      REAL,               -- for arcs
    sweep_name  TEXT,               -- symbolic sweep name (e.g. '_PI_2')
    sweep       REAL,               -- sweep in radians
    center_name TEXT,               -- arc center point name
    n_pts       INTEGER DEFAULT 60, -- arc discretisation
    end_name    TEXT NOT NULL        -- produced point name
);

CREATE TABLE IF NOT EXISTS views (
    name        TEXT PRIMARY KEY,
    label       TEXT NOT NULL,
    script      TEXT NOT NULL,       -- relative path to generator script
    svg_path    TEXT NOT NULL,       -- relative path to output SVG/PDF
    category    TEXT DEFAULT 'design',
    enabled     INTEGER DEFAULT 1
);
"""


def init_db(db_path=None):
    """Create tables and seed data if the database does not exist."""
    db_path = db_path or DB_PATH
    fresh = not os.path.exists(db_path)
    with get_db(db_path) as conn:
        conn.executescript(_SCHEMA)
        if fresh:
            _seed_constants(conn)
            _seed_outline_chain(conn)
            _seed_views(conn)


# ---------------------------------------------------------------------------
# Seed: constants from floorplan/constants.py
# ---------------------------------------------------------------------------

def _categorize(name: str) -> str:
    """Assign a category based on the constant name."""
    n = name.upper()
    if any(n.startswith(p) for p in ("WALL_", "SHELL_", "AIR_")):
        return "wall"
    if any(n.startswith(p) for p in ("IW", "RO")):
        return "interior_wall"
    if n.startswith("O") and (n[1:2].isdigit() or n.startswith("O8") or n.startswith("OPENING")):
        return "opening"
    if any(n.startswith(p) for p in ("APPLIANCE", "COUNTER", "KITCHEN", "DW_", "STOVE", "FRIDGE", "MINIK", "ICE_")):
        return "appliance"
    if any(n.startswith(p) for p in ("BED_", "DRESSER", "SHELVES", "LOVESEAT", "DESK", "CHAIR_",
                                      "OTTOMAN", "SOFA_", "ROCKER", "ET_")):
        return "furniture"
    if any(n.startswith(p) for p in ("WH_", "TOILET", "SINK", "BATH_SINK")):
        return "fixture"
    if any(n.startswith(p) for p in ("CORNER_", "F11AB", "F2_", "F10_", "ROOF_")):
        return "geometry"
    if any(n.startswith(p) for p in ("DOOR_", "F8F9", "JAMB", "STD_")):
        return "construction"
    return "misc"


def _parse_constants_source():
    """Parse floorplan/constants.py to extract constants with expressions and comments."""
    src_path = os.path.join(_PROJECT, "floorplan", "constants.py")
    with open(src_path, "r") as f:
        source = f.read()

    results = []
    for line in source.splitlines():
        line = line.strip()
        if not line or line.startswith("#") or line.startswith("assert") or line.startswith('"'):
            continue
        m = re.match(r'^([A-Z_][A-Z0-9_]*)\s*=\s*(.+?)(?:\s*#\s*(.*))?$', line)
        if m:
            name, expr, comment = m.group(1), m.group(2).strip(), (m.group(3) or "").strip()
            results.append((name, expr, comment))
    return results


def _seed_constants(conn):
    """Extract constants from the Python module and insert into the database."""
    import importlib
    import floorplan.constants as mod
    # Reload from disk so patched-in values from engine.patch_constants()
    # are discarded and the original source values are used.
    importlib.reload(mod)

    parsed = _parse_constants_source()
    name_to_info = {name: (expr, comment) for name, expr, comment in parsed}

    for attr in sorted(dir(mod)):
        if attr.startswith("_") or not attr[0].isupper():
            continue
        val = getattr(mod, attr)
        if not isinstance(val, (int, float)):
            continue
        expr, desc = name_to_info.get(attr, (str(val), ""))
        cat = _categorize(attr)
        unit = "cm" if "CM" in attr else "ft"
        conn.execute(
            "INSERT OR REPLACE INTO constants (name, value, expr, unit, category, description) "
            "VALUES (?, ?, ?, ?, ?, ?)",
            (attr, float(val), expr, unit, cat, desc),
        )


# ---------------------------------------------------------------------------
# Seed: outline chain
# ---------------------------------------------------------------------------

def _seed_outline_chain(conn):
    """Seed the outline chain from floorplan/geometry.py definitions."""
    from floorplan.geometry import OUTLINE_CHAIN, CHAIN_POINT_NAMES
    import math

    sweep_names = {
        math.pi / 2: "_PI_2",
        5 * math.pi / 12: "_5PI_12",
        math.pi / 12: "_PI_12",
    }

    for seq, (entry, end_name) in enumerate(zip(OUTLINE_CHAIN, CHAIN_POINT_NAMES)):
        if entry[0] == "L":
            conn.execute(
                "INSERT INTO outline_chain (seq, seg_type, distance, end_name) VALUES (?, ?, ?, ?)",
                (seq, "L", entry[1], end_name),
            )
        else:
            direction, radius, sweep = entry[0], entry[1], entry[2]
            center_name = entry[3]
            n_pts = entry[4]
            s_name = sweep_names.get(sweep, f"{sweep:.12f}")
            conn.execute(
                "INSERT INTO outline_chain "
                "(seq, seg_type, radius, sweep_name, sweep, center_name, n_pts, end_name) "
                "VALUES (?, ?, ?, ?, ?, ?, ?, ?)",
                (seq, direction, radius, s_name, sweep, center_name, n_pts, end_name),
            )


# ---------------------------------------------------------------------------
# Seed: available views
# ---------------------------------------------------------------------------

def _seed_views(conn):
    """Register all generator scripts and their outputs."""
    views = [
        ("floorplan", "Floorplan", "floorplan/gen_floorplan.py", "floorplan/floorplan.svg", "design"),
        ("walls", "Wall Detail", "walls/gen_walls.py", "walls/walls.svg", "construction"),
        ("all_walls", "All Walls", "walls/gen_walls.py", "walls/all_walls.svg", "construction"),
        ("span", "Span Analysis", "span/gen_span.py", "span/span.svg", "analysis"),
        ("span_min", "Min Span", "span/gen_span_min.py", "span/span_min.svg", "analysis"),
        ("span_minmax", "Span vs Rotation", "span/gen_span_minmax.py", "span/span_minmax.svg", "analysis"),
        ("path_area", "Survey Path", "survey/gen_path_svg.py", "survey/path_area.svg", "survey"),
        ("roof", "Roof", "roof/gen_roof.py", "roof/roof.svg", "design"),
        ("plumbing", "Plumbing", "plumbing/gen_plumbing.py", "plumbing/plumbing.svg", "design"),
        ("site_plan_df", "Site Plan (DF)", "site/gen_site_plan.py", "site/site_plan_df.pdf", "site"),
        ("site_plan_fs", "Site Plan (FS)", "site/gen_site_plan.py", "site/site_plan_fs.pdf", "site"),
    ]
    for name, label, script, svg_path, cat in views:
        conn.execute(
            "INSERT OR REPLACE INTO views (name, label, script, svg_path, category, enabled) "
            "VALUES (?, ?, ?, ?, ?, 1)",
            (name, label, script, svg_path, cat),
        )


# ---------------------------------------------------------------------------
# CRUD operations
# ---------------------------------------------------------------------------

def get_all_constants(db_path=None):
    """Return all constants as a list of dicts."""
    with get_db(db_path) as conn:
        rows = conn.execute(
            "SELECT name, value, expr, unit, category, description FROM constants ORDER BY category, name"
        ).fetchall()
        return [dict(r) for r in rows]


def get_constants_dict(db_path=None):
    """Return constants as a simple name→value dict for the geometry engine."""
    with get_db(db_path) as conn:
        rows = conn.execute("SELECT name, value FROM constants").fetchall()
        return {r["name"]: r["value"] for r in rows}


def update_constant(name, value, db_path=None):
    """Update a single constant value. Returns True if changed."""
    with get_db(db_path) as conn:
        cur = conn.execute(
            "UPDATE constants SET value = ?, expr = ? WHERE name = ?",
            (float(value), str(value), name),
        )
        return cur.rowcount > 0


def update_constants_batch(updates: dict, db_path=None):
    """Update multiple constants at once. Returns number of rows changed."""
    changed = 0
    with get_db(db_path) as conn:
        for name, value in updates.items():
            cur = conn.execute(
                "UPDATE constants SET value = ?, expr = ? WHERE name = ?",
                (float(value), str(value), name),
            )
            changed += cur.rowcount
    return changed


def get_categories(db_path=None):
    """Return distinct constant categories."""
    with get_db(db_path) as conn:
        rows = conn.execute("SELECT DISTINCT category FROM constants ORDER BY category").fetchall()
        return [r["category"] for r in rows]


def get_outline_chain(db_path=None):
    """Return outline chain as list of dicts."""
    with get_db(db_path) as conn:
        rows = conn.execute("SELECT * FROM outline_chain ORDER BY seq").fetchall()
        return [dict(r) for r in rows]


def get_views(db_path=None):
    """Return all registered views."""
    with get_db(db_path) as conn:
        rows = conn.execute("SELECT * FROM views WHERE enabled = 1 ORDER BY category, name").fetchall()
        return [dict(r) for r in rows]


def reset_constants(db_path=None):
    """Reset all constants to their original values from source code."""
    db_path = db_path or DB_PATH
    with get_db(db_path) as conn:
        conn.execute("DELETE FROM constants")
        _seed_constants(conn)
