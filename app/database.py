"""SQLite database layer for the ADU building definition.

Stores all building constants, outline chain definitions, and entity metadata
so the building can be edited interactively and regenerated from the database.
"""
import ast
import json
import os
import re
import sqlite3
from contextlib import contextmanager

from app.apputil import ARC_N_SEMICIRCLE

_DIR = os.path.dirname(os.path.abspath(__file__))
_PROJECT = os.path.dirname(_DIR)
DB_PATH = os.path.join(_PROJECT, "app", "adu.db")

# Inner wall override column list — shared across SELECT/INSERT queries
_IW_OV_COLS = ("seg_index, span_end, sub_seq, seg_type, bearing, "
               "distance, radius, sweep, n_pts")

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
    end_name    TEXT NOT NULL,       -- produced point name
    flex        TEXT DEFAULT NULL,   -- solved param: 'distance','radius','sweep', or NULL
    bearing_flex INTEGER DEFAULT 0   -- 1 = line bearing may rotate (opt-in); 0 = fixed (default)
);

CREATE TABLE IF NOT EXISTS views (
    name        TEXT PRIMARY KEY,
    label       TEXT NOT NULL,
    script      TEXT NOT NULL,       -- relative path to generator script
    svg_path    TEXT NOT NULL,       -- relative path to output SVG/PDF
    category    TEXT DEFAULT 'design',
    enabled     INTEGER DEFAULT 1
);

CREATE TABLE IF NOT EXISTS shapes (
    name        TEXT PRIMARY KEY,
    poly_json   TEXT NOT NULL,       -- JSON: [[dx, dy], ...] in local coords
    scale       REAL DEFAULT 1.0,    -- conversion factor from poly units to feet
    origin      TEXT DEFAULT 'center', -- 'center' or 'corner'
    width_key   TEXT,                -- constant name for width dimension, if any
    depth_key   TEXT,                -- constant name for depth dimension, if any
    description TEXT DEFAULT ''
);

CREATE TABLE IF NOT EXISTS variant_exclusions (
    variant      TEXT NOT NULL,       -- variant name (e.g. 'bare', 'sf')
    element_type TEXT NOT NULL,       -- 'wall', 'opening', 'rough_opening'
    element_name TEXT NOT NULL,       -- e.g. 'IW6', 'RO5', 'O3'
    PRIMARY KEY (variant, element_type, element_name)
);

CREATE TABLE IF NOT EXISTS room_label_offsets (
    room_name   TEXT PRIMARY KEY,     -- e.g. 'BEDROOM', 'UTIL_N'
    offset_e    REAL DEFAULT 0.0,     -- easting offset from centroid (feet)
    offset_n    REAL DEFAULT 0.0      -- northing offset from centroid (feet)
);

CREATE TABLE IF NOT EXISTS undo_history (
    id          INTEGER PRIMARY KEY AUTOINCREMENT,
    timestamp   REAL NOT NULL,
    action_type TEXT NOT NULL,        -- e.g. 'constant_update', 'constant_batch', 'constant_reset'
    before_state TEXT NOT NULL,       -- JSON serialised previous state
    after_state  TEXT NOT NULL,       -- JSON serialised new state
    description TEXT DEFAULT ''
);

CREATE TABLE IF NOT EXISTS elements (
    id         INTEGER PRIMARY KEY AUTOINCREMENT,
    type       TEXT NOT NULL,          -- wall, furniture, appliance, fixture, opening, clearance, label, dimension
    name       TEXT NOT NULL UNIQUE,
    properties TEXT DEFAULT '{}',      -- JSON: metadata, dimensions, position
    variant    TEXT                    -- NULL = all variants; non-null restricts to one
);

CREATE TABLE IF NOT EXISTS doors (
    id              INTEGER PRIMARY KEY AUTOINCREMENT,
    opening_name    TEXT NOT NULL UNIQUE,     -- e.g. 'RO1'
    width           REAL NOT NULL,            -- door width in inches
    hinge_side      TEXT NOT NULL,            -- 'east', 'west', 'north', 'south'
    swing_direction TEXT NOT NULL,            -- 'east', 'west', 'north', 'south'
    door_type       TEXT NOT NULL DEFAULT 'single'  -- 'single' or 'double'
);

CREATE TABLE IF NOT EXISTS config (
    key   TEXT PRIMARY KEY,
    value TEXT NOT NULL
);

CREATE TABLE IF NOT EXISTS plumbing_elements (
    id         INTEGER PRIMARY KEY AUTOINCREMENT,
    type       TEXT NOT NULL,
    name       TEXT NOT NULL UNIQUE,
    path       TEXT DEFAULT '[]',
    properties TEXT DEFAULT '{}',
    fixture    TEXT,
    sort_order INTEGER DEFAULT 0
);

CREATE TABLE IF NOT EXISTS variants (
    id              INTEGER PRIMARY KEY AUTOINCREMENT,
    name            TEXT NOT NULL UNIQUE,
    label           TEXT NOT NULL,
    source_variant  TEXT,
    flags           TEXT DEFAULT '{}',
    layer_config    TEXT DEFAULT '{}',
    is_builtin      INTEGER DEFAULT 0
);

CREATE TABLE IF NOT EXISTS survey_legs (
    seq            INTEGER PRIMARY KEY,
    bearing_deg    INTEGER NOT NULL,
    bearing_min    INTEGER NOT NULL,
    bearing_sec    INTEGER NOT NULL,
    distance_ft    INTEGER NOT NULL,
    distance_inch  REAL NOT NULL,
    label          TEXT               -- endpoint label: P2, P3, P4, P5, POB
);

CREATE TABLE IF NOT EXISTS survey_config (
    key   TEXT PRIMARY KEY,
    value TEXT NOT NULL               -- JSON-encoded value
);

CREATE TABLE IF NOT EXISTS inner_wall_overrides (
    seg_index   INTEGER NOT NULL,   -- span start inner_seg index (0-based)
    span_end    INTEGER,            -- span end inner_seg index (NULL = single segment)
    sub_seq     INTEGER NOT NULL,   -- ordering within the override chain
    seg_type    TEXT NOT NULL,       -- 'L' (line), 'CW', or 'CCW' (arc)
    bearing     REAL,               -- degrees, compass convention (line segs)
    distance    REAL,               -- feet (line segments)
    radius      REAL,               -- feet (arc segments)
    sweep       REAL,               -- degrees (arc segments)
    n_pts       INTEGER DEFAULT 20, -- arc discretization point count
    PRIMARY KEY (seg_index, sub_seq)
);

CREATE TABLE IF NOT EXISTS element_formulas (
    id           INTEGER PRIMARY KEY AUTOINCREMENT,
    element_name TEXT NOT NULL COLLATE NOCASE,  -- FK to elements.name
    param_name   TEXT NOT NULL,          -- 'position', 'poly', 'center', etc.
    formula_json TEXT NOT NULL,          -- JSON formula spec
    locked       INTEGER DEFAULT 0,     -- 1 = frozen at locked_value
    locked_value TEXT,                   -- JSON: frozen computed result
    variant      TEXT,                   -- NULL = all variants
    UNIQUE(element_name, param_name, variant)
);

CREATE TABLE IF NOT EXISTS formula_deps (
    element_name TEXT NOT NULL COLLATE NOCASE,  -- element with the formula
    param_name   TEXT NOT NULL,          -- which parameter formula
    dep_type     TEXT NOT NULL,          -- 'element', 'point', 'constant'
    dep_name     TEXT NOT NULL COLLATE NOCASE,  -- name of dependency target
    UNIQUE(element_name, param_name, dep_type, dep_name)
);

CREATE TABLE IF NOT EXISTS catalog_items (
    key         TEXT PRIMARY KEY,
    item_type   TEXT NOT NULL,          -- 'furniture', 'appliance', 'fixture'
    label       TEXT NOT NULL,
    width       REAL,                   -- default width in feet (NULL for circles)
    depth       REAL,                   -- default depth in feet (NULL for circles)
    radius      REAL,                   -- radius in feet (for circle shapes)
    shape       TEXT NOT NULL DEFAULT 'rect',
    door        TEXT,                   -- JSON: door config or variant-keyed dict
    clearance   TEXT,                   -- JSON: clearance config
    product_url TEXT,                   -- JSON: string or variant-keyed dict
    variants    TEXT NOT NULL DEFAULT '[]',  -- JSON array of variant names
    stacked     INTEGER DEFAULT 0,
    clip_to_inner INTEGER DEFAULT 0
);
"""


def _create_catalog_table(conn):
    """Create catalog_items table if it doesn't exist."""
    conn.execute("""
        CREATE TABLE IF NOT EXISTS catalog_items (
            key         TEXT PRIMARY KEY,
            item_type   TEXT NOT NULL,
            label       TEXT NOT NULL,
            width       REAL,
            depth       REAL,
            radius      REAL,
            shape       TEXT NOT NULL DEFAULT 'rect',
            door        TEXT,
            clearance   TEXT,
            product_url TEXT,
            variants    TEXT NOT NULL DEFAULT '[]',
            stacked     INTEGER DEFAULT 0,
            clip_to_inner INTEGER DEFAULT 0
        )
    """)


def init_db(db_path=None):
    """Create tables and seed data if the database does not exist."""
    db_path = db_path or DB_PATH
    fresh = not os.path.exists(db_path)
    with get_db(db_path) as conn:
        conn.executescript(_SCHEMA)
        if fresh:
            _seed_constants(conn)
            _seed_variant_item_constants(conn)
            _seed_outline_chain(conn)
            _seed_survey(conn)
            _seed_inner_wall_overrides(conn)
            _seed_views(conn)
            _seed_shapes(conn)
            _seed_variant_exclusions(conn)
            _seed_variants(conn)
            _seed_elements(conn)
            _seed_doors(conn)
            _seed_config(conn)
            from app.labels import seed_room_labels, seed_builtin_dimensions
            seed_room_labels(conn)
            seed_builtin_dimensions(conn)
            from app.plumbing import seed_plumbing
            seed_plumbing(conn)
            seed_iw_formulas(conn)
        else:
            # Ensure variant item constants exist (Phase 12a upgrade)
            _seed_variant_item_constants(conn)
            # Ensure all seed doors exist (handles additions like O3, O6)
            _seed_doors(conn)
            # Ensure variant exclusions exist (Phase 7 upgrade)
            _seed_variant_exclusions(conn)
            # Ensure variant definitions exist (Phase 11a upgrade)
            _seed_variants(conn)
            # Ensure config defaults exist (Phase 10b upgrade)
            _seed_config(conn)
            # Ensure room label elements exist (Phase 8 upgrade)
            from app.labels import seed_room_labels, seed_builtin_dimensions
            seed_room_labels(conn)
            # Ensure builtin dimension elements exist (unified dims upgrade)
            seed_builtin_dimensions(conn)
            # Ensure plumbing fixture records exist (Phase 10d upgrade)
            from app.plumbing import seed_plumbing
            seed_plumbing(conn)
            # Ensure IW wall formulas exist (Phase 12c upgrade)
            seed_iw_formulas(conn)
            # Ensure survey data exists (Phase 14-B upgrade)
            _seed_survey(conn)
            # Add span_end column if missing (Phase 15½ multi-segment upgrade)
            cols = {r[1] for r in conn.execute(
                "PRAGMA table_info(inner_wall_overrides)").fetchall()}
            if "span_end" not in cols:
                conn.execute("ALTER TABLE inner_wall_overrides "
                             "ADD COLUMN span_end INTEGER")
            # Ensure inner wall overrides exist (Phase 15½ upgrade)
            _seed_inner_wall_overrides(conn)
            # Add sort_order column if missing (plumbing supplies upgrade)
            pe_cols = {r[1] for r in conn.execute(
                "PRAGMA table_info(plumbing_elements)").fetchall()}
            if "sort_order" not in pe_cols:
                conn.execute("ALTER TABLE plumbing_elements "
                             "ADD COLUMN sort_order INTEGER DEFAULT 0")
            # Add flex column to outline_chain if missing
            oc_cols = {r[1] for r in conn.execute(
                "PRAGMA table_info(outline_chain)").fetchall()}
            if "flex" not in oc_cols:
                conn.execute("ALTER TABLE outline_chain "
                             "ADD COLUMN flex TEXT DEFAULT NULL")
                # Set defaults: first line, second-to-last line, last arc
                rows = conn.execute(
                    "SELECT seq FROM outline_chain ORDER BY seq"
                ).fetchall()
                if rows:
                    n = len(rows)
                    conn.execute(
                        "UPDATE outline_chain SET flex = 'distance' "
                        "WHERE seq = 0")
                    conn.execute(
                        "UPDATE outline_chain SET flex = 'distance' "
                        "WHERE seq = ?", (n - 2,))
                    conn.execute(
                        "UPDATE outline_chain SET flex = 'sweep' "
                        "WHERE seq = ?", (n - 1,))
            if "bearing_flex" not in oc_cols:
                conn.execute("ALTER TABLE outline_chain "
                             "ADD COLUMN bearing_flex INTEGER DEFAULT 0")
            # Migrate element_formulas/formula_deps to COLLATE NOCASE
            _migrate_nocase_formulas(conn)
            # Migrate existing placed items: create formulas if missing
            _migrate_placed_item_formulas(conn)
            # Add prop_constants to IW wall records if missing
            _migrate_iw_prop_constants(conn)
    # Seed anchor/pivot defaults after DB is committed (needs constants table)
    if fresh:
        _seed_default_anchor_pivot(db_path)
    else:
        # Migrate: seed anchor coords if missing (new in this phase)
        if get_outline_anchor_pos(db_path) is None:
            _seed_default_anchor_pivot(db_path)
    # Seed catalog items after DB is committed (needs geometry computation)
    _seed_catalog_items_post(db_path)


def _seed_catalog_items_post(db_path):
    """Seed catalog_items after DB is committed (separate connection)."""
    with get_db(db_path) as conn:
        _seed_catalog_items(conn, db_path)
        # Fix placed items with 0x0 dimensions from catalog
        _fix_zero_dim_placed_items(conn)


def _fix_zero_dim_placed_items(conn):
    """Fix placed items that have width=0/depth=0 by looking up catalog dims."""
    rows = conn.execute(
        "SELECT e.id, e.name, e.properties, ef.formula_json "
        "FROM elements e "
        "LEFT JOIN element_formulas ef "
        "  ON ef.element_name = e.name AND ef.param_name = 'position' "
        "WHERE e.type IN ('appliance', 'furniture', 'fixture')"
    ).fetchall()
    for r in rows:
        props = json.loads(r["properties"]) if isinstance(
            r["properties"], str) else (r["properties"] or {})
        cat_key = props.get("catalog_key")
        if not cat_key:
            continue
        pw = props.get("width", 0) or 0
        pd = props.get("depth", 0) or 0
        if pw > 0 and pd > 0:
            continue
        # Look up catalog dimensions
        cat = conn.execute(
            "SELECT width, depth FROM catalog_items WHERE key = ?",
            (cat_key,)).fetchone()
        if not cat or not cat["width"] or not cat["depth"]:
            continue
        cw, cd = cat["width"], cat["depth"]
        # Fix element properties
        props["width"] = cw
        props["depth"] = cd
        conn.execute("UPDATE elements SET properties = ? WHERE id = ?",
                     (json.dumps(props), r["id"]))
        # Fix formula
        if r["formula_json"]:
            fj = json.loads(r["formula_json"]) if isinstance(
                r["formula_json"], str) else r["formula_json"]
            if fj and fj.get("type") == "item_rect":
                fj["width"] = cw
                fj["depth"] = cd
                conn.execute(
                    "UPDATE element_formulas SET formula_json = ? "
                    "WHERE element_name = ? AND param_name = 'position'",
                    (json.dumps(fj), r["name"]))


def _migrate_nocase_formulas(conn):
    """Migrate element_formulas/formula_deps to COLLATE NOCASE columns.

    SQLite doesn't support ALTER COLUMN, so we recreate the tables with
    COLLATE NOCASE on element_name/dep_name.  This makes all lookups and
    UNIQUE constraints case-insensitive, fixing the case mismatch between
    formula names (e.g. BED) and UI element names (e.g. bed).
    """
    # Check if migration is needed by inspecting column collation
    info = conn.execute("PRAGMA table_info(element_formulas)").fetchall()
    # PRAGMA table_info doesn't expose COLLATE — check the CREATE TABLE SQL
    sql = conn.execute(
        "SELECT sql FROM sqlite_master WHERE name = 'element_formulas'"
    ).fetchone()
    if sql and "COLLATE NOCASE" in sql[0]:
        return  # Already migrated

    # --- element_formulas ---
    conn.execute("""
        CREATE TABLE IF NOT EXISTS element_formulas_new (
            id           INTEGER PRIMARY KEY AUTOINCREMENT,
            element_name TEXT NOT NULL COLLATE NOCASE,
            param_name   TEXT NOT NULL,
            formula_json TEXT NOT NULL,
            locked       INTEGER DEFAULT 0,
            locked_value TEXT,
            variant      TEXT,
            UNIQUE(element_name, param_name, variant)
        )
    """)
    conn.execute("""
        INSERT OR IGNORE INTO element_formulas_new
            (id, element_name, param_name, formula_json, locked, locked_value, variant)
        SELECT id, element_name, param_name, formula_json, locked, locked_value, variant
        FROM element_formulas
    """)
    conn.execute("DROP TABLE element_formulas")
    conn.execute("ALTER TABLE element_formulas_new RENAME TO element_formulas")

    # --- formula_deps ---
    conn.execute("""
        CREATE TABLE IF NOT EXISTS formula_deps_new (
            element_name TEXT NOT NULL COLLATE NOCASE,
            param_name   TEXT NOT NULL,
            dep_type     TEXT NOT NULL,
            dep_name     TEXT NOT NULL COLLATE NOCASE,
            UNIQUE(element_name, param_name, dep_type, dep_name)
        )
    """)
    conn.execute("""
        INSERT OR IGNORE INTO formula_deps_new
            (element_name, param_name, dep_type, dep_name)
        SELECT element_name, param_name, dep_type, dep_name
        FROM formula_deps
    """)
    conn.execute("DROP TABLE formula_deps")
    conn.execute("ALTER TABLE formula_deps_new RENAME TO formula_deps")

    # Remove orphaned uppercase layout_item element records (BED, DRYER, etc.)
    # that collide case-insensitively with lowercase variant items.
    _OLD_LAYOUT_ITEMS = ("BED", "COUNTER", "DRESSER", "DRYER", "SHELVES", "WASHER")
    for name in _OLD_LAYOUT_ITEMS:
        conn.execute(
            "DELETE FROM elements WHERE name = ? AND "
            "json_extract(properties, '$.layout_item') = 1",
            (name,),
        )


def _migrate_placed_item_formulas(conn):
    """Create formulas for placed items that don't have them yet."""
    # Shapes that have entries in the shapes table (custom polygons)
    shape_table_names = {r[0] for r in conn.execute(
        "SELECT name FROM shapes").fetchall()}

    rows = conn.execute(
        "SELECT name, properties FROM elements "
        "WHERE json_extract(properties, '$.source') = 'placed'"
    ).fetchall()
    for row in rows:
        name = row["name"]
        existing = conn.execute(
            "SELECT id FROM element_formulas WHERE element_name = ?",
            (name,),
        ).fetchone()
        if existing:
            continue
        props = json.loads(row["properties"]) if isinstance(
            row["properties"], str) else (row["properties"] or {})
        center = props.get("center", [0, 0])
        shape = props.get("shape", "rect")
        rotation = props.get("rotation", 0)
        # Resolve shape name: prefer catalog_key for shapes table lookup
        catalog_key = props.get("catalog_key")
        shape_key = shape
        if shape not in ("rect", "circle"):
            if catalog_key and catalog_key in shape_table_names:
                shape_key = catalog_key
        # Dedicated formula types for specific catalog items
        if catalog_key == "dining_table":
            import math as _math
            rad = rotation * _math.pi / 180
            cos_r, sin_r = _math.cos(rad), _math.sin(rad)
            formula = {
                "type": "dining_triangle",
                "base_center": center,
                "toward_apex": [sin_r, -cos_r],
                "along_base": [cos_r, sin_r],
                "base_width": {"const": "DINING_TBL_BASE"},
                "height": {"const": "DINING_TBL_H"},
                "apex_radius": 1.0,
                "fillet_radius": 0.5,
            }
        elif shape == "toilet":
            import math as _math
            rad = rotation * _math.pi / 180
            cos_r, sin_r = _math.cos(rad), _math.sin(rad)
            formula = {
                "type": "toilet_shape",
                "center": center,
                "facing_dir": [-sin_r, cos_r],
                "width_dir": [cos_r, sin_r],
            }
        elif shape == "bath_sink":
            import math as _math
            rad = rotation * _math.pi / 180
            cos_r, sin_r = _math.cos(rad), _math.sin(rad)
            formula = {
                "type": "bath_sink_shape",
                "anchor": center,
                "along": [cos_r, sin_r],
                "outward": [-sin_r, cos_r],
                "length": {"const": "BATH_SINK_LENGTH"},
                "depth": {"const": "BATH_SINK_DEPTH"},
            }
        elif shape == "circle":
            radius = props.get("radius") or props.get("width", 1) / 2
            formula = {"type": "item_circle", "center": center,
                       "radius": radius}
        elif shape_key in shape_table_names:
            # Use shape_transform for items with custom shape polygons
            formula = {
                "type": "shape_transform",
                "shape_name": shape_key,
                "center": center,
                "rotation_deg": rotation,
            }
            if props.get("width"):
                formula["width"] = props["width"]
            if props.get("depth"):
                formula["depth"] = props["depth"]
        else:
            w = props.get("width", 1)
            d = props.get("depth", 1)
            formula = {
                "type": "item_rect", "anchor": center,
                "along": [1, 0], "across": [0, 1],
                "width": w, "depth": d,
                "anchor_corner": "center",
                "rotation_deg": rotation,
            }
        conn.execute(
            "INSERT INTO element_formulas "
            "(element_name, param_name, formula_json, variant) "
            "VALUES (?, 'position', ?, NULL)",
            (name, json.dumps(formula)),
        )


def _migrate_iw_prop_constants(conn):
    """Sync prop_constants for all IW wall records to the canonical IW_PROP_CONSTANTS
    mapping (upgrade migration — runs on every startup for existing DBs)."""
    from app.elements import IW_PROP_CONSTANTS
    rows = conn.execute(
        "SELECT id, name, properties FROM elements "
        "WHERE type='wall' AND name GLOB 'IW*'"
    ).fetchall()
    for row in rows:
        name = row["name"]
        try:
            props = json.loads(row["properties"]) if row["properties"] else {}
        except Exception:
            props = {}
        expected = IW_PROP_CONSTANTS.get(name, {})
        if props.get("prop_constants") != expected:
            props["prop_constants"] = expected
            conn.execute(
                "UPDATE elements SET properties=? WHERE id=?",
                (json.dumps(props), row["id"]),
            )


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
# Seed: variant item dimension constants (from app/variants.py)
# ---------------------------------------------------------------------------

# These 24 constants are hardcoded in gen_floorplan.py and replicated in
# app/variants.py but not defined in floorplan/constants.py.  Phase 12a
# moves them into the database as first-class constants.
_VARIANT_ITEM_CONSTANTS = [
    # (name, value_inches, category, description)
    ("HAMPER_W",         31.5,   "furniture",  "Hamper width"),
    ("HAMPER_D",         19.0,   "furniture",  "Hamper depth"),
    ("MINIK_APPL_W",     32.0,   "appliance",  "Small kitchen appliance width"),
    ("MINIK_APPL_D",     27.0,   "appliance",  "Small kitchen appliance depth"),
    ("MICROWAVE_W",      19.5,   "appliance",  "Microwave width"),
    ("MICROWAVE_D",      16.625, "appliance",  "Microwave depth"),
    ("COFFEE_W",          7.2,   "appliance",  "Coffee maker width"),
    ("COFFEE_D",          9.2,   "appliance",  "Coffee maker depth"),
    ("COOKTOP_W",        13.4,   "appliance",  "Cooktop width"),
    ("COOKTOP_D",        20.5,   "appliance",  "Cooktop depth"),
    ("TOASTER_W",        13.7,   "appliance",  "Toaster width"),
    ("TOASTER_D",        12.5,   "appliance",  "Toaster depth"),
    ("DINING_CHAIR_W",   18.0,   "furniture",  "Dining chair width"),
    ("DINING_CHAIR_D",   21.0,   "furniture",  "Dining chair depth"),
    ("DINING_TBL_BASE",  31.5,   "furniture",  "Dining table base width"),
    ("DINING_TBL_H",     35.25,  "furniture",  "Dining table height (plan depth)"),
    ("DAYBED_W",         86.0,   "furniture",  "Daybed width"),
    ("DAYBED_D",         43.0,   "furniture",  "Daybed depth"),
    ("WORK_CTR_W",       60.0,   "furniture",  "Work counter width"),
    ("WORK_CTR_D",       18.0,   "furniture",  "Work counter depth"),
    ("STD_FRIDGE_W",     32.75,  "appliance",  "Standard fridge width"),
    ("STD_FRIDGE_D",     35.0,   "appliance",  "Standard fridge depth"),
    ("SOFA_FULL_W",      80.75,  "furniture",  "Full sofa width"),
    ("SOFA_FULL_D",      34.625, "furniture",  "Full sofa depth"),
]


def _seed_variant_item_constants(conn):
    """Seed dimension constants for variant items (Phase 12a)."""
    for name, value_in, category, description in _VARIANT_ITEM_CONSTANTS:
        value_ft = value_in / 12.0
        expr = f"{value_in} / 12.0"
        conn.execute(
            "INSERT OR IGNORE INTO constants "
            "(name, value, expr, unit, category, description) "
            "VALUES (?, ?, ?, ?, ?, ?)",
            (name, value_ft, expr, "ft", category, description),
        )

    # F2 origin position (Phase 14-A) — already in feet
    for name, value, description in [
        ("F2_EASTING",  -18.5, "F2 easting offset from FC"),
        ("F2_NORTHING", -13.5, "F2 northing offset from FC (before R_a1)"),
    ]:
        conn.execute(
            "INSERT OR IGNORE INTO constants "
            "(name, value, expr, unit, category, description) "
            "VALUES (?, ?, ?, ?, ?, ?)",
            (name, value, str(value), "ft", "geometry", description),
        )


# ---------------------------------------------------------------------------
# Seed: survey data (Phase 14-B)
# ---------------------------------------------------------------------------

# Raw traverse legs from shared/survey.py:_accumulate_legs()
_SURVEY_LEGS = [
    # (seq, deg, min, sec, ft, inch, label)
    (1, 257, 53, 45, 19, 1.0,  "P2"),
    (2, 180, 54, 31, 26, 11.0, "P3"),
    (3,  93, 36,  7, 31, 10.5, "P4"),
    (4,  56, 36, 31, 13,  2.5, "P5"),
    (5, 317, 11, 44, 34, 11.5, "POB"),
]

# Survey configuration constants
_SURVEY_CONFIG = {
    "FC_IN_P3_E": 18.5141152720,
    "FC_IN_P3_N": 13.3968094375,
    "COORD_ROTATION": 0.0015153784,
    "P3_EASTING_OVERRIDE": -19.1177,
    "P2_P3_NORTHING_OFFSET": 29.0,
}


def _seed_inner_wall_overrides(conn):
    """Seed the W8-W9 straight-arc-straight override for segment 5 (F8→F9).

    Computes parametric values from the default geometry: two straight
    segments flanking a 90° CCW arc at the concave F8-F9 corner.

    Only seeds if no overrides exist yet — avoids clobbering databases
    where the override lives at a different seg_index (e.g. Mark2).
    """
    count = conn.execute(
        "SELECT COUNT(*) FROM inner_wall_overrides").fetchone()[0]
    if count > 0:
        return
    import math
    from floorplan.geometry import compute_outline_geometry
    from floorplan.constants import (
        WALL_OUTER, SHELL_THICKNESS, OPENING_INSIDE_RADIUS,
    )
    g = compute_outline_geometry()
    pts = g.fp_pts
    F8, C7 = pts["F8"], pts["C7"]
    F9, F10 = pts["F9"], pts["F10"]

    # CW traversal direction at F8 (tangent at exit of C7 arc)
    r8x, r8y = F8[0] - C7[0], F8[1] - C7[1]
    r8_len = math.hypot(r8x, r8y)
    dir_f8 = (r8y / r8_len, -r8x / r8_len)

    # CW traversal direction at F9 (F9→F10 line)
    d9x, d9y = F10[0] - F9[0], F10[1] - F9[1]
    d9_len = math.hypot(d9x, d9y)
    dir_f9 = (d9x / d9_len, d9y / d9_len)

    brg_f8 = math.degrees(math.atan2(dir_f8[0], dir_f8[1]))
    brg_f9 = math.degrees(math.atan2(dir_f9[0], dir_f9[1]))

    R_turn = OPENING_INSIDE_RADIUS + SHELL_THICKNESS

    # Inset points (W8, W9)
    ins_f8 = (dir_f8[1], -dir_f8[0])  # right of CW = inward
    ins_f9 = (dir_f9[1], -dir_f9[0])
    W8 = (F8[0] + WALL_OUTER * ins_f8[0], F8[1] + WALL_OUTER * ins_f8[1])
    W9 = (F9[0] + WALL_OUTER * ins_f9[0], F9[1] + WALL_OUTER * ins_f9[1])

    # Tangent point offsets: arc center is R_turn to the LEFT of each direction
    left_f8 = (-dir_f8[1], dir_f8[0])
    left_f9 = (-dir_f9[1], dir_f9[0])

    # Arc center via intersection
    P1 = (W8[0] + R_turn * left_f8[0], W8[1] + R_turn * left_f8[1])
    P2 = (W9[0] + R_turn * left_f9[0], W9[1] + R_turn * left_f9[1])
    dp = (P2[0] - P1[0], P2[1] - P1[1])
    cross = dir_f8[0] * dir_f9[1] - dir_f8[1] * dir_f9[0]
    t = (dp[0] * dir_f9[1] - dp[1] * dir_f9[0]) / cross
    arc_tangent1 = (W8[0] + t * dir_f8[0], W8[1] + t * dir_f8[1])

    # Re-derive from end side
    t2_dp = (P1[0] - P2[0], P1[1] - P2[1])
    t2_cross = dir_f9[0] * dir_f8[1] - dir_f9[1] * dir_f8[0]
    t2 = (t2_dp[0] * dir_f8[1] - t2_dp[1] * dir_f8[0]) / t2_cross
    arc_tangent2 = (W9[0] + t2 * dir_f9[0], W9[1] + t2 * dir_f9[1])

    dist_start = math.hypot(arc_tangent1[0] - W8[0], arc_tangent1[1] - W8[1])
    dist_end = math.hypot(W9[0] - arc_tangent2[0], W9[1] - arc_tangent2[1])

    # Arc sweep (CCW)
    entry = math.atan2(-left_f8[1], -left_f8[0])
    exit_ = math.atan2(-left_f9[1], -left_f9[0])
    sweep_rad = (exit_ - entry) % (2 * math.pi)
    sweep_deg = math.degrees(sweep_rad)

    chain = [
        (0, "L", brg_f8, dist_start, None, None, 20),
        (1, "CCW", None, None, R_turn, sweep_deg, 20),
        (2, "L", brg_f9, dist_end, None, None, 20),
    ]
    for sub_seq, seg_type, bearing, distance, radius, sweep, n_pts in chain:
        conn.execute(
            f"INSERT OR IGNORE INTO inner_wall_overrides "
            f"({_IW_OV_COLS}) "
            "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)",
            (5, None, sub_seq, seg_type, bearing, distance, radius, sweep,
             n_pts),
        )


def _seed_survey(conn):
    """Seed survey traverse legs and config from hardcoded values."""
    for seq, deg, mn, sec, ft, inch, label in _SURVEY_LEGS:
        conn.execute(
            "INSERT OR IGNORE INTO survey_legs "
            "(seq, bearing_deg, bearing_min, bearing_sec, "
            "distance_ft, distance_inch, label) "
            "VALUES (?, ?, ?, ?, ?, ?, ?)",
            (seq, deg, mn, sec, ft, inch, label),
        )
    for key, value in _SURVEY_CONFIG.items():
        conn.execute(
            "INSERT OR IGNORE INTO survey_config (key, value) VALUES (?, ?)",
            (key, json.dumps(value)),
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

    # Set default flex segments: first line, second-to-last line, closure arc
    n = len(OUTLINE_CHAIN)
    conn.execute("UPDATE outline_chain SET flex = 'distance' WHERE seq = 0")
    conn.execute("UPDATE outline_chain SET flex = 'distance' WHERE seq = ?",
                 (n - 2,))
    conn.execute("UPDATE outline_chain SET flex = 'sweep' WHERE seq = ?",
                 (n - 1,))


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
        ("plumbing", "Plumbing", "floorplan/gen_floorplan.py", "floorplan/floorplan_plumbing.svg", "design"),
        ("site_plan_df", "Site Plan (DF)", "site/gen_site_plan.py", "site/site_plan_df.pdf", "site"),
        ("site_plan_fs", "Site Plan (FS)", "site/gen_site_plan.py", "site/site_plan_fs.pdf", "site"),
        ("3d_flat", "3D Flat Roof", "scad/gen_flat_roof.py", "scad/flat_roof_patio.png", "3d"),
        ("3d_2in12", "3D 2:12 Roof", "scad/gen_2in12.py", "scad/2in12_patio.png", "3d"),
        ("3views", "3-View Layout", "gen_3views.py", "3views.pdf", "3d"),
    ]
    for name, label, script, svg_path, cat in views:
        conn.execute(
            "INSERT OR REPLACE INTO views (name, label, script, svg_path, category, enabled) "
            "VALUES (?, ?, ?, ?, ?, 1)",
            (name, label, script, svg_path, cat),
        )


# ---------------------------------------------------------------------------
# Seed: config
# ---------------------------------------------------------------------------

def _seed_config(conn):
    """Seed default configuration values."""
    defaults = [
        ("roof_style", "flat"),
        ("setback_216", "11.0"),
        ("setback_275", "25.5"),
    ]
    for key, value in defaults:
        conn.execute(
            "INSERT OR IGNORE INTO config (key, value) VALUES (?, ?)",
            (key, value),
        )


# ---------------------------------------------------------------------------
# Seed: item shapes
# ---------------------------------------------------------------------------

def _seed_shapes(conn):
    """Register built-in special item shapes."""

    import math

    # Toilet plan-view polygon (from gen_floorplan.py _TOILET_SVG).
    # Local frame: dx = across width (centered on 0), dy = depth from wall.
    # Scale: SVG units → feet via 10.0 / 30.48.
    toilet_pts = [
        [-1.905, 0], [-1.905, 2.032], [-0.841, 2.032],
        [-1.078, 2.224], [-1.267, 2.455], [-1.408, 2.719],
        [-1.495, 3.005], [-1.524, 3.302],
        [-1.732, 5.461], [-1.699, 5.799], [-1.600, 6.124],
        [-1.440, 6.423], [-1.225, 6.686], [-0.962, 6.901],
        [-0.663, 7.061], [-0.338, 7.160], [0, 7.193],
        [0.338, 7.160], [0.663, 7.061], [0.962, 6.901],
        [1.225, 6.686], [1.440, 6.423], [1.600, 6.124],
        [1.699, 5.799], [1.732, 5.461],
        [1.524, 3.302], [1.495, 3.005], [1.408, 2.719],
        [1.267, 2.455], [1.078, 2.224], [0.847, 2.035],
        [0.841, 2.032], [1.905, 2.032], [1.905, 0],
    ]
    conn.execute(
        "INSERT OR REPLACE INTO shapes "
        "(name, poly_json, scale, origin, description) VALUES (?, ?, ?, ?, ?)",
        ("toilet", json.dumps(toilet_pts), 10.0 / 30.48, "center",
         "Realistic toilet plan view; center = wall face at centerline, "
         "dy = outward from wall"),
    )

    # Bath sink (Tripoli wall-mount, 33-7/8 x 18-3/4) with semicircular bulge.
    # Local frame: dx = along wall (centered on 0), dy = outward from wall.
    # Dimensions in feet (scale = 1.0).
    bs_length = 33.875 / 12.0
    bs_depth = 18.75 / 12.0
    half_len = bs_length / 2
    quarter_len = bs_length / 4
    rect_depth = bs_depth - quarter_len
    bs_pts = [[half_len, 0], [-half_len, 0], [-half_len, rect_depth]]
    for i in range(ARC_N_SEMICIRCLE + 1):
        t = math.pi - math.pi * i / ARC_N_SEMICIRCLE
        bs_pts.append([round(math.cos(t) * quarter_len, 6),
                       round(rect_depth + math.sin(t) * quarter_len, 6)])
    bs_pts.append([half_len, rect_depth])
    conn.execute(
        "INSERT OR REPLACE INTO shapes "
        "(name, poly_json, scale, origin, width_key, depth_key, description) "
        "VALUES (?, ?, ?, ?, ?, ?, ?)",
        ("bath_sink", json.dumps(bs_pts), 1.0, "center",
         "BATH_SINK_LENGTH", "BATH_SINK_DEPTH",
         "Tripoli wall-mount sink; rectangle with semicircular bulge "
         "from central 50% of far edge"),
    )

    # Dining table (Oscar triangle set): base + two tangent lines + apex arc
    # + corner fillets.  Stored as pre-computed polygon in local frame
    # (dx = along base centered on 0, dy = toward apex).
    # Dimensions in feet (scale = 1.0).
    tbl_base = 31.5 / 12.0
    tbl_h = 35.25 / 12.0
    apex_r = 12.0 / 12.0
    fillet_r = 6.0 / 12.0
    half_base = tbl_base / 2
    apex_y = tbl_h
    arc_cy = apex_y - apex_r
    # Tangent from NE corner to apex arc
    dx_r = half_base
    dy_r = -arc_cy
    dist_r = math.sqrt(dx_r ** 2 + dy_r ** 2)
    angle_cp = math.atan2(dy_r, dx_r)
    delta_a = math.acos(apex_r / dist_r)
    alpha_r = angle_cp - delta_a
    tr_x = apex_r * math.cos(alpha_r)
    tr_y = arc_cy + apex_r * math.sin(alpha_r)
    tl_x = -tr_x
    tl_y = tr_y
    # NE fillet
    d_base = (-1, 0)
    dtr_len = math.sqrt((tr_x - half_base) ** 2 + tr_y ** 2)
    d_tang = ((tr_x - half_base) / dtr_len, tr_y / dtr_len)
    cos_th = d_base[0] * d_tang[0] + d_base[1] * d_tang[1]
    half_ang = math.acos(max(-1, min(1, cos_th))) / 2
    fillet_dist = fillet_r / math.sin(half_ang)
    bis = (d_base[0] + d_tang[0], d_base[1] + d_tang[1])
    bis_len = math.sqrt(bis[0] ** 2 + bis[1] ** 2)
    bis = (bis[0] / bis_len, bis[1] / bis_len)
    fc_ne = (half_base + fillet_dist * bis[0], fillet_dist * bis[1])
    f_ne_base_x = fc_ne[0]
    f_ne_base_y = 0
    v_ne = (fc_ne[0] - half_base, fc_ne[1])
    t_proj = v_ne[0] * d_tang[0] + v_ne[1] * d_tang[1]
    f_ne_tang = (half_base + t_proj * d_tang[0], t_proj * d_tang[1])
    f_nw_base_x = -f_ne_base_x
    f_nw_tang = (-f_ne_tang[0], f_ne_tang[1])
    fc_nw = (-fc_ne[0], fc_ne[1])
    # Build polygon
    tbl_pts = [[round(-f_ne_base_x, 6), 0], [round(f_ne_base_x, 6), 0]]
    # NE fillet arc
    a0 = math.atan2(f_ne_base_y - fc_ne[1], f_ne_base_x - fc_ne[0])
    a1 = math.atan2(f_ne_tang[1] - fc_ne[1], f_ne_tang[0] - fc_ne[0])
    sw = (a1 - a0) % (2 * math.pi)
    if sw > math.pi:
        sw -= 2 * math.pi
    for i in range(1, 8):
        t = a0 + sw * i / 8
        tbl_pts.append([round(fc_ne[0] + fillet_r * math.cos(t), 6),
                        round(fc_ne[1] + fillet_r * math.sin(t), 6)])
    tbl_pts.append([round(f_ne_tang[0], 6), round(f_ne_tang[1], 6)])
    tbl_pts.append([round(tr_x, 6), round(tr_y, 6)])
    # Apex arc
    a0 = math.atan2(tr_y - arc_cy, tr_x)
    a1 = math.atan2(tl_y - arc_cy, tl_x)
    sw = (a1 - a0) % (2 * math.pi)
    if sw > math.pi:
        sw -= 2 * math.pi
    for i in range(1, 16):
        t = a0 + sw * i / 16
        tbl_pts.append([round(apex_r * math.cos(t), 6),
                        round(arc_cy + apex_r * math.sin(t), 6)])
    tbl_pts.append([round(tl_x, 6), round(tl_y, 6)])
    tbl_pts.append([round(f_nw_tang[0], 6), round(f_nw_tang[1], 6)])
    # NW fillet arc
    a0 = math.atan2(f_nw_tang[1] - fc_nw[1], f_nw_tang[0] - fc_nw[0])
    a1 = math.atan2(-fc_nw[1], f_nw_base_x - fc_nw[0])
    sw = (a1 - a0) % (2 * math.pi)
    if sw > math.pi:
        sw -= 2 * math.pi
    for i in range(1, 8):
        t = a0 + sw * i / 8
        tbl_pts.append([round(fc_nw[0] + fillet_r * math.cos(t), 6),
                        round(fc_nw[1] + fillet_r * math.sin(t), 6)])
    conn.execute(
        "INSERT OR REPLACE INTO shapes "
        "(name, poly_json, scale, origin, description) VALUES (?, ?, ?, ?, ?)",
        ("dining_table", json.dumps(tbl_pts), 1.0, "center",
         "Oscar triangle dining table; base centered at origin, "
         "apex arc at dy=height, 6\" corner fillets"),
    )


# ---------------------------------------------------------------------------
# Seed: variant exclusions
# ---------------------------------------------------------------------------

def _seed_variant_exclusions(conn):
    """Register elements hidden in specific layout variants."""
    exclusions = [
        # bare (Room Dimensions): IW6 and RO5 are omitted per gen_floorplan.py
        ("bare", "wall", "IW6"),
        ("bare", "rough_opening", "RO5"),
        # bare: dim12a/dim12b replaced by dim12bare
        ("bare", "dimension", "dim12a"),
        ("bare", "dimension", "dim12b"),
        # sf (Square Footage): same exclusions as bare
        ("sf", "wall", "IW6"),
        ("sf", "rough_opening", "RO5"),
        ("sf", "dimension", "dim12a"),
        ("sf", "dimension", "dim12b"),
    ]
    for variant, etype, ename in exclusions:
        conn.execute(
            "INSERT OR REPLACE INTO variant_exclusions "
            "(variant, element_type, element_name) VALUES (?, ?, ?)",
            (variant, etype, ename),
        )


# ---------------------------------------------------------------------------
# Seed: variants
# ---------------------------------------------------------------------------

_VARIANT_SEEDS = [
    ("standard", "Standard",       "{}",                 1),
    ("minik",    "Small Kitchen",   '{"minik": true}',   1),
    ("daybed",   "Daybed",          '{"db": true}',      1),
    ("bare",     "Room Dimensions", '{"bare": true}',    1),
    ("sf",       "Square Footage",  '{"sf": true}',      1),
    ("plumbing", "Plumbing",        '{"plumbing": true}', 1),
]


def _seed_variants(conn):
    """Seed built-in variant definitions."""
    for name, label, flags, is_builtin in _VARIANT_SEEDS:
        conn.execute(
            "INSERT OR IGNORE INTO variants (name, label, flags, is_builtin) "
            "VALUES (?, ?, ?, ?)",
            (name, label, flags, is_builtin),
        )


# ---------------------------------------------------------------------------
# Seed: elements (13 interior walls)
# ---------------------------------------------------------------------------

# IW name → (thickness constant, orientation)
_IW_SEED = [
    ("IW1",  "WALL_6IN",        "H"),
    ("IW2",  "WALL_6IN",        "V"),
    ("IW2O", "IW2O_THICKNESS",  "R"),
    ("IW2S", "WALL_6IN",        "V"),
    ("IW3",  "WALL_4IN",        "V"),
    ("IW4",  "WALL_4IN",        "V"),
    ("IW5",  "WALL_4IN",        "H"),
    ("IW6",  "IW6_THICKNESS",   "H"),
    ("IW7",  "WALL_4IN",        "H"),
    ("IW8",  "WALL_4IN",        "V"),
    ("IW9",  "WALL_4IN",        "V"),
    ("IW11", "WALL_4IN",        "V"),
    ("IW12", "WALL_4IN",        "H"),
]


# ---------------------------------------------------------------------------
# Variant items master list (module-level, used by both _seed_elements and
# _seed_catalog_items)
# ---------------------------------------------------------------------------
_VI_ALL = ["standard", "minik", "daybed"]
_VI_ALL_P = ["standard", "minik", "daybed", "plumbing"]
_VI_STD_DB = ["standard", "daybed"]
_VI_STD_DB_P = ["standard", "daybed", "plumbing"]

_VARIANT_ITEMS = [
    # --- Utility / laundry ---
    ("dryer", "appliance", {
        "label": "DRYER", "item_type": "appliance", "shape": "rect",
        "door": {"hinge_idx": 1, "target_idx": 2},
        "product_url": {"minik": "https://www.lowes.com/pd/Electrolux-8-cu-ft-Stackable-Steam-Cycle-Electric-Dryer-Titanium-ENERGY-STAR/5015416377"},
        "variants": _VI_ALL,
    }, None),
    ("washer", "appliance", {
        "label": "WASHER", "item_type": "appliance", "shape": "rect",
        "door": {"hinge_idx": 2, "target_idx": 1},
        "product_url": {"minik": "https://www.lowes.com/pd/Electrolux-Smartboost-Optic-Whites-and-Pure-Rinse-4-5-cu-ft-High-Efficiency-Stackable-Steam-Cycle-Front-Load-Washer-Titanium-ENERGY-STAR/5015416375"},
        "variants": _VI_ALL_P,
    }, None),
    ("hamper", "appliance", {
        "label": "HAMPER", "item_type": "appliance", "shape": "rect",
        "clearance": {"face": [3, 2], "distance": 19.0 / 12.0},
        "product_url": "https://www.homedepot.com/p/Casual-Home-Eco-Home-Laundry-Prep-Hamper-761-30/307595219",
        "variants": _VI_ALL,
    }, None),
    ("water_heater", "appliance", {
        "label": "WH", "item_type": "appliance", "shape": "circle",
        "variants": _VI_ALL_P,
    }, None),
    ("counter", "appliance", {
        "label": "COUNTER", "item_type": "appliance", "shape": "rect",
        "clip_to_inner": True,
        "variants": _VI_ALL,
    }, None),
    # --- Toilets & sinks ---
    ("toilet_s", "fixture", {
        "label": "TOILET", "item_type": "fixture", "shape": "toilet",
        "variants": _VI_ALL_P,
    }, None),
    ("toilet_n", "fixture", {
        "label": "TOILET", "item_type": "fixture", "shape": "toilet",
        "variants": _VI_ALL_P,
    }, None),
    ("util_sink", "fixture", {
        "label": "SINK", "item_type": "fixture", "shape": "rect",
        "product_url": "https://www.magnushomeproducts.com/products/24-petten-matte-gray-vitreous-china-console-sink-with-black-powdercoat-steel-stand-and-shelves",
        "variants": _VI_ALL_P,
    }, None),
    ("bath_sink", "fixture", {
        "label": "BATH SINK", "item_type": "fixture", "shape": "bath_sink",
        "product_url": "https://www.magnushomeproducts.com/products/tripoli-vitreous-china-wall-mount-bathroom-sink",
        "variants": _VI_ALL_P,
    }, None),
    ("kitchen_sink", "fixture", {
        "label": "SINK", "item_type": "fixture", "shape": "rect",
        "product_url": "https://www.webstaurantstore.com/advance-tabco-fs1181824l-45-fabricated-one-compartment-sink-with-24-left-drainboard-18-x-18-x-14-bowl/109FS1L241818.html",
        "variants": _VI_ALL_P,
    }, None),
    # --- Kitchen: standard/daybed ---
    ("stove", "appliance", {
        "label": "STOVE", "item_type": "appliance", "shape": "rect",
        "clearance": {"face": [0, 1], "distance": 24.0 / 12.0},
        "variants": _VI_STD_DB,
    }, None),
    ("dishwasher", "appliance", {
        "label": "D/W", "item_type": "appliance", "shape": "rect",
        "clearance": {"face": [0, 1], "distance": 31.0 / 12.0},
        "variants": _VI_STD_DB_P,
    }, None),
    ("north_counter", "appliance", {
        "label": "COUNTER", "item_type": "appliance", "shape": "rect",
        "product_url": "https://www.webstaurantstore.com/regency-spec-line-30-x-36-14-gauge-stainless-steel-commercial-work-table-with-4-backsplash-and-undershelf/600TSSB3036S.html",
        "variants": _VI_STD_DB,
    }, None),
    ("work_counter", "appliance", {
        "label": "COUNTER", "item_type": "appliance", "shape": "rect",
        "product_url": "https://www.webstaurantstore.com/table-s-s-18x60-s-s-under/600TS1860S.html",
        "variants": _VI_STD_DB,
    }, None),
    # --- Kitchen: minik ---
    ("kitchen_counter", "appliance", {
        "label": "COUNTER", "item_type": "appliance", "shape": "rect",
        "product_url": "https://www.webstaurantstore.com/regency-spec-line-30-x-72-14-gauge-stainless-steel-commercial-work-table-with-4-backsplash-and-undershelf/600TSSB3072S.html",
        "variants": ["minik"],
    }, None),
    ("cooktop", "appliance", {
        "label": "COOKTOP", "item_type": "appliance", "shape": "rect",
        "stacked": True,
        "product_url": "https://www.homedepot.com/p/Empava-Portable-13-4-in-Induction-Electric-Cooktop-in-Black-with-2-Elements-EMPV-ID12/313815692",
        "variants": ["minik"],
    }, None),
    ("toaster", "appliance", {
        "label": "TOASTER", "item_type": "appliance", "shape": "rect",
        "stacked": True,
        "product_url": "https://www.amazon.com/Roter-Mond-Stainless-Independent-Removable/dp/B0CGTQZTDZ?th=1",
        "variants": ["minik"],
    }, None),
    # --- Kitchen: all variants ---
    ("fridge", "appliance", {
        "label": "FRIDGE", "item_type": "appliance", "shape": "rect",
        "door": {
            "standard": {"hinge_idx": 3, "target_idx": 2},
            "minik": {"hinge_idx": 1, "target_idx": 0},
            "daybed": {"hinge_idx": 3, "target_idx": 2},
        },
        "product_url": {"minik": "https://www.ikea.com/us/en/p/bergsnaes-bottom-freezer-refrigerator-stainless-steel-color-60607883/", "default": "https://www.lowes.com/pd/LG-25-5-cu-ft-Bottom-Freezer-Refrigerator-with-Ice-Maker-Fingerprint-Resistant-Printproof-Stainless-Steel-ENERGY-STAR/1002543648"},
        "variants": _VI_ALL_P,
    }, None),
    ("ice_maker", "appliance", {
        "label": "ICE", "item_type": "appliance", "shape": "rect",
        "product_url": "https://www.homedepot.com/p/EUHOMY-17-3-in-100-lb-24H-Full-Ice-Sizes-Commercial-Ice-Maker-in-Black-33-lb-Storage-Bin-Ice-Full-Alert-and-Auto-Cleaning-CIM001-100BL-E/337185876",
        "variants": _VI_ALL_P,
    }, None),
    ("microwave", "appliance", {
        "label": "MICRO", "item_type": "appliance", "shape": "rect",
        "stacked": True,
        "door": {
            "standard": {"hinge_idx": 2, "target_idx": 3},
            "minik": {"hinge_idx": 0, "target_idx": 1},
            "daybed": {"hinge_idx": 2, "target_idx": 3},
        },
        "product_url": "https://www.ikea.com/us/en/p/gatebo-microwave-oven-with-air-fryer-function-ikea-500-black-70603506/",
        "variants": _VI_ALL,
    }, None),
    ("coffee_maker", "appliance", {
        "label": "C", "item_type": "appliance", "shape": "rect",
        "stacked": True,
        "product_url": "https://www.amazon.com/Holstein-Housewares-HH-0914701E-5-Cup-Coffee/dp/B08HSRCC4T/?th=1",
        "variants": _VI_ALL,
    }, None),
    # --- Dining ---
    ("dining_table", "furniture", {
        "label": "TABLE", "item_type": "furniture", "shape": "triangle",
        "product_url": "https://www.homedepot.com/pep/NEW-CLASSIC-HOME-FURNISHINGS-New-Classic-Furniture-Oscar-3-Piece-Wood-Top-Triangle-Dining-Set-Walnut-40-1651-D2C/327836175",
        "variants": _VI_ALL,
    }, None),
    ("dining_chair_1", "furniture", {
        "label": "CHAIR", "item_type": "furniture", "shape": "rect",
        "variants": _VI_ALL,
    }, None),
    ("dining_chair_2", "furniture", {
        "label": "CHAIR", "item_type": "furniture", "shape": "rect",
        "variants": _VI_ALL,
    }, None),
    # --- Bedroom ---
    ("bed", "furniture", {
        "label": "KING BED", "item_type": "furniture", "shape": "rect",
        "variants": _VI_ALL,
    }, None),
    ("dresser", "furniture", {
        "label": "DRESSER", "item_type": "furniture", "shape": "rect",
        "clearance": {"face": [0, 1], "distance": 15.0 / 12.0},
        "variants": _VI_ALL,
    }, None),
    ("shelves", "furniture", {
        "label": "SHELVES", "item_type": "furniture", "shape": "rect",
        "product_url": "https://www.ikea.com/us/en/p/kallax-shelving-unit-with-underframe-white-stained-oak-effect-black-s49442718/",
        "variants": _VI_ALL,
    }, None),
    # --- Living: chair + ottoman (all) ---
    ("chair", "furniture", {
        "label": "CHAIR", "item_type": "furniture", "shape": "rect",
        "product_url": "https://www.ikea.com/us/en/p/havberg-swivel-easy-chair-and-footstool-grann-bomstad-golden-brown-s59485321/",
        "variants": _VI_ALL,
    }, None),
    ("ottoman", "furniture", {
        "label": "OTTO", "item_type": "furniture", "shape": "rect",
        "product_url": "https://www.ikea.com/us/en/p/havberg-swivel-easy-chair-and-footstool-grann-bomstad-golden-brown-s59485321/",
        "variants": _VI_ALL,
    }, None),
    # --- Living: standard seating ---
    ("loveseat", "furniture", {
        "label": "LOVESEAT", "item_type": "furniture", "shape": "rect",
        "product_url": "https://www.ikea.com/us/en/p/saltsjoebaden-loveseat-tonerud-red-brown-s59579188/",
        "variants": ["standard"],
    }, "standard"),
    ("et", "furniture", {
        "label": "ET", "item_type": "furniture", "shape": "circle",
        "product_url": "https://www.ikea.com/us/en/p/listerby-side-table-oak-veneer-30515314/",
        "variants": ["standard"],
    }, "standard"),
    ("loveseat2", "furniture", {
        "label": "LOVESEAT", "item_type": "furniture", "shape": "rect",
        "product_url": "https://www.ikea.com/us/en/p/saltsjoebaden-loveseat-tonerud-red-brown-s59579188/",
        "variants": ["standard"],
    }, "standard"),
    # --- Living: minik seating ---
    ("sofa", "furniture", {
        "label": "SOFA", "item_type": "furniture", "shape": "rect",
        "product_url": "https://www.ikea.com/us/en/p/saltsjoebaden-3-seat-sofa-gunnared-light-green-s89599953/",
        "variants": ["minik"],
    }, "minik"),
    ("rocker", "furniture", {
        "label": "ROCKER", "item_type": "furniture", "shape": "rect",
        "product_url": {"minik": "https://www.ikea.com/us/en/p/poaeng-rocking-chair-brown-gunnared-beige-s39502048/", "db": "https://www.ikea.com/us/en/p/poaeng-rocking-chair-brown-gunnared-beige-s39502048/"},
        "variants": ["minik", "daybed"],
    }, None),
    # --- Living: daybed seating ---
    ("shelves2", "furniture", {
        "label": "SHELVES", "item_type": "furniture", "shape": "rect",
        "product_url": "https://www.ikea.com/us/en/p/kallax-shelving-unit-with-underframe-white-stained-oak-effect-black-s49442718/",
        "variants": ["daybed"],
    }, "daybed"),
    ("et_east", "furniture", {
        "label": "ET", "item_type": "furniture", "shape": "circle",
        "product_url": "https://www.ikea.com/us/en/p/listerby-side-table-oak-veneer-30515314/",
        "variants": ["daybed"],
    }, "daybed"),
    ("daybed", "furniture", {
        "label": "DAYBED", "item_type": "furniture", "shape": "rect",
        "variants": ["daybed"],
    }, "daybed"),
    ("et_west", "furniture", {
        "label": "ET", "item_type": "furniture", "shape": "circle",
        "product_url": "https://www.ikea.com/us/en/p/listerby-side-table-oak-veneer-30515314/",
        "variants": ["daybed"],
    }, "daybed"),
    # --- Office ---
    ("desk", "furniture", {
        "label": "DESK", "item_type": "furniture", "shape": "rect",
        "variants": _VI_ALL,
    }, None),
    ("desk_chair", "furniture", {
        "label": "CHAIR", "item_type": "furniture", "shape": "rect",
        "product_url": "https://www.amazon.com/BESTFAIR-Ergonomic-Office-Chair-Adjustable/dp/B0FDQDMP2D?th=1",
        "variants": _VI_ALL,
    }, None),
]


def _seed_elements(conn):
    """Seed the elements table with interior walls, openings, and variant items."""

    # --- Interior walls (13) ---
    from app.elements import IW_PROP_CONSTANTS
    for name, thickness_const, orientation in _IW_SEED:
        props = json.dumps({
            "thickness_constant": thickness_const,
            "orientation": orientation,
            "prop_constants": IW_PROP_CONSTANTS.get(name, {}),
        })
        conn.execute(
            "INSERT OR REPLACE INTO elements (type, name, properties, variant) "
            "VALUES (?, ?, ?, ?)",
            ("wall", name, props, None),
        )

    # --- Outer openings (O1-O11, O8a) ---
    _OUTER_OPENINGS = [
        ("O1",  {"seg_start": "F2",  "seg_end": "F5",  "opening_type": "window"}),
        ("O2",  {"seg_start": "F2",  "seg_end": "F5",  "opening_type": "window"}),
        ("O3",  {"seg_start": "F2",  "seg_end": "F5",  "opening_type": "door"}),
        ("O4",  {"seg_start": "F6",  "seg_end": "F7",  "opening_type": "window"}),
        ("O5",  {"seg_start": "F9",  "seg_end": "F10", "opening_type": "window"}),
        ("O6",  {"seg_start": "F9",  "seg_end": "F10", "opening_type": "door"}),
        ("O7",  {"seg_start": "F12", "seg_end": "F13", "opening_type": "window"}),
        ("O8",  {"seg_start": "F14", "seg_end": "F15", "opening_type": "casement"}),
        ("O8a", {"seg_start": "F18", "seg_end": "F1",  "opening_type": "window"}),
        ("O9",  {"seg_start": "F18", "seg_end": "F1",  "opening_type": "casement"}),
        ("O10", {"seg_start": "F18", "seg_end": "F1",  "opening_type": "casement"}),
        ("O11", {"seg_start": "F18", "seg_end": "F1",  "opening_type": "window"}),
    ]
    for name, props in _OUTER_OPENINGS:
        conn.execute(
            "INSERT OR REPLACE INTO elements (type, name, properties, variant) "
            "VALUES (?, ?, ?, ?)",
            ("opening", name, json.dumps(props), None),
        )

    # --- Rough openings (RO1-RO7) ---
    _ROUGH_OPENINGS = [
        ("RO1", {"wall_name": "IW1",  "orientation": "H"}),
        ("RO2", {"wall_name": "IW11", "orientation": "V"}),
        ("RO3", {"wall_name": "IW9",  "orientation": "V"}),
        ("RO4", {"wall_name": "IW2O", "orientation": "V"}),
        ("RO5", {"wall_name": "IW6",  "orientation": "H"}),
        ("RO6", {"wall_name": "IW11", "orientation": "V"}),
        ("RO7", {"wall_name": "IW9",  "orientation": "V"}),
    ]
    for name, props in _ROUGH_OPENINGS:
        conn.execute(
            "INSERT OR REPLACE INTO elements (type, name, properties, variant) "
            "VALUES (?, ?, ?, ?)",
            ("opening", name, json.dumps(props), None),
        )

    for name, elem_type, props, variant in _VARIANT_ITEMS:
        conn.execute(
            "INSERT OR REPLACE INTO elements (type, name, properties, variant) "
            "VALUES (?, ?, ?, ?)",
            (elem_type, name, json.dumps(props), variant),
        )

    # Layout item placeholders (BED, COUNTER, etc.) were removed — layout
    # formulas now use lowercase names matching _VARIANT_ITEMS directly.


# ---------------------------------------------------------------------------
# Seed: catalog items
# ---------------------------------------------------------------------------

def _seed_catalog_items(conn, db_path=None):
    """Seed catalog_items from _VARIANT_ITEMS with computed dimensions.

    Uses INSERT OR IGNORE so existing entries are not overwritten on upgrade.
    Dimensions come from geometry computation at seed time.
    """
    # Compute dimensions from geometry for all variants
    dims = {}  # key -> {width, depth} or {radius}
    constants = {}
    try:
        constants = {r[0]: r[1] for r in conn.execute(
            "SELECT name, value FROM constants").fetchall()}
    except Exception:
        pass
    try:
        from app.engine import compute_geometry
        for v in ("standard", "minik", "daybed"):
            from app.engine import compute_geometry
            geom = compute_geometry(constants, variant=v, db_path=db_path)
            vi = geom.get("variant_items", {})
            for name, item in vi.items():
                if name in dims:
                    continue
                b = item.get("bbox", {})
                shape = item.get("shape", "rect")
                if shape == "circle":
                    dims[name] = {"radius": round(item.get("radius", 0), 6)}
                else:
                    dims[name] = {
                        "width": round(b.get("e", 0) - b.get("w", 0), 6),
                        "depth": round(b.get("n", 0) - b.get("s", 0), 6),
                    }
    except Exception:
        pass  # Dimensions will be NULL if geometry computation fails

    # Resolve dimensions from element formulas' constant references
    # for any items not already in dims (geometry may miss catalog templates)
    def _eval_expr(expr):
        """Evaluate a simple constant expression (const, mul, add, sub)."""
        if isinstance(expr, (int, float)):
            return expr
        if isinstance(expr, dict):
            if "const" in expr:
                return constants.get(expr["const"])
            if "mul" in expr:
                result = 1.0
                for operand in expr["mul"]:
                    v = _eval_expr(operand)
                    if v is None:
                        return None
                    result *= v
                return result
            if "add" in expr:
                result = 0.0
                for operand in expr["add"]:
                    v = _eval_expr(operand)
                    if v is None:
                        return None
                    result += v
                return result
            if "sub" in expr:
                parts = expr["sub"]
                if len(parts) >= 2:
                    a = _eval_expr(parts[0])
                    b = _eval_expr(parts[1])
                    if a is not None and b is not None:
                        return a - b
        return None

    try:
        rows = conn.execute(
            "SELECT element_name, formula_json FROM element_formulas "
            "WHERE param_name = 'position'"
        ).fetchall()
        for r in rows:
            name = r[0].lower()
            if name in dims:
                continue
            fj = json.loads(r[1]) if isinstance(r[1], str) else r[1]
            if not fj:
                continue
            ftype = fj.get("type", "")

            # item_rect: width/depth as const refs or literals
            if ftype == "item_rect":
                w_ref = fj.get("width")
                d_ref = fj.get("depth")
                w_val = _eval_expr(w_ref)
                d_val = _eval_expr(d_ref)
                if w_val and d_val:
                    dims[name] = {"width": round(w_val, 6),
                                  "depth": round(d_val, 6)}

            # item_circle: radius as const ref or expression
            elif ftype == "item_circle":
                r_val = _eval_expr(fj.get("radius"))
                if r_val and r_val > 0:
                    dims[name] = {"radius": round(r_val, 6)}

            # dining_chair: chair_width/chair_depth
            elif ftype == "dining_chair":
                w_val = _eval_expr(fj.get("chair_width"))
                d_val = _eval_expr(fj.get("chair_depth"))
                if w_val and d_val:
                    dims[name] = {"width": round(w_val, 6),
                                  "depth": round(d_val, 6)}

            # dining_triangle: base_width/height
            elif ftype == "dining_triangle":
                w_val = _eval_expr(fj.get("base_width") or
                                   fj.get("base_side"))
                d_val = _eval_expr(fj.get("height"))
                if w_val and d_val:
                    dims[name] = {"width": round(w_val, 6),
                                  "depth": round(d_val, 6)}

            # toilet_shape: use TOILET_WIDTH / (TOILET_WIDTH + TOILET_TANK_DEPTH)
            elif ftype == "toilet_shape":
                tw = constants.get("TOILET_WIDTH")
                td = constants.get("TOILET_TANK_DEPTH")
                if tw and td:
                    dims[name] = {"width": round(tw, 6),
                                  "depth": round(tw + td, 6)}

            # ellipse_rect: rx/ry → width=2*rx, depth=2*ry
            elif ftype == "ellipse_rect":
                rx = _eval_expr(fj.get("rx"))
                ry = _eval_expr(fj.get("ry"))
                if rx and ry:
                    dims[name] = {"width": round(rx * 2, 6),
                                  "depth": round(ry * 2, 6)}

            # bath_sink_shape: BATH_SINK_LENGTH / BATH_SINK_DEPTH
            elif ftype == "bath_sink_shape":
                bl = constants.get("BATH_SINK_LENGTH")
                bd = _eval_expr(fj.get("depth")) or constants.get(
                    "BATH_SINK_DEPTH")
                if bl and bd:
                    dims[name] = {"width": round(bl, 6),
                                  "depth": round(bd, 6)}

            # four_corner: look up known dimension constants by item name
            elif ftype == "four_corner":
                _fc_dims = {
                    "counter": ("COUNTER_LENGTH", "COUNTER_DEPTH"),
                    "chair": ("CHAIR_WIDTH", "CHAIR_DEPTH"),
                    "ottoman": ("OTTOMAN_SIZE", "OTTOMAN_SIZE"),
                    "loveseat": ("LOVESEAT_LENGTH", "LOVESEAT_WIDTH"),
                    "sofa": ("SOFA_WIDTH", "SOFA_DEPTH"),
                    "rocker": ("ROCKER_WIDTH", "ROCKER_DEPTH"),
                    "daybed": ("DAYBED_W", "DAYBED_D"),
                }
                pair = _fc_dims.get(name)
                if pair:
                    w_val = constants.get(pair[0])
                    d_val = constants.get(pair[1])
                    if w_val and d_val:
                        dims[name] = {"width": round(w_val, 6),
                                      "depth": round(d_val, 6)}

    except Exception:
        pass

    for name, elem_type, props, _variant in _VARIANT_ITEMS:
        d = dims.get(name, {})
        product_url = props.get("product_url")
        door = props.get("door")
        clearance = props.get("clearance")
        variants = props.get("variants", [])
        conn.execute(
            "INSERT OR IGNORE INTO catalog_items "
            "(key, item_type, label, width, depth, radius, shape, "
            " door, clearance, product_url, variants, stacked, clip_to_inner) "
            "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
            (name, elem_type, props.get("label", name.upper()),
             d.get("width"), d.get("depth"), d.get("radius"),
             props.get("shape", "rect"),
             json.dumps(door) if door else None,
             json.dumps(clearance) if clearance else None,
             json.dumps(product_url) if product_url else None,
             json.dumps(variants),
             1 if props.get("stacked") else 0,
             1 if props.get("clip_to_inner") else 0),
        )

    # Backfill: update existing entries that have NULL dimensions
    for name in dims:
        d = dims[name]
        if "radius" in d:
            conn.execute(
                "UPDATE catalog_items SET radius = ? "
                "WHERE key = ? AND radius IS NULL",
                (d["radius"], name))
        elif "width" in d and "depth" in d:
            conn.execute(
                "UPDATE catalog_items SET width = ?, depth = ? "
                "WHERE key = ? AND (width IS NULL OR depth IS NULL)",
                (d["width"], d["depth"], name))


# ---------------------------------------------------------------------------
# Seed: doors (O3, O6, RO1–RO7 defaults)
# ---------------------------------------------------------------------------

# (opening_name, door_width_constant, hinge_side, swing_direction, door_type)
_DOOR_SEED = [
    ("O3",  "O3_DOOR_WIDTH",  "north", "east",  "single"),
    ("O6",  "O6_DOOR_WIDTH",  "east",  "south", "single"),
    ("RO1", "RO1_DOOR_WIDTH", "east",  "south", "single"),
    ("RO2", "RO2_DOOR_WIDTH", "north", "east",  "single"),
    ("RO3", "RO3_DOOR_WIDTH", "south", "west",  "single"),
    ("RO4", "RO4_DOOR_WIDTH", "south", "west",  "single"),
    ("RO5", "RO5_DOOR_WIDTH", "east",  "north", "single"),
    ("RO6", "RO6_DOOR_WIDTH", "west",  "west",  "double"),
    ("RO7", "RO7_DOOR_WIDTH", "east",  "east",  "double"),
]


def _seed_doors(conn):
    """Seed the doors table with default configurations for O3, O6, RO1–RO7."""
    import importlib
    import floorplan.constants as mod
    importlib.reload(mod)

    for opening_name, width_const, hinge, swing, dtype in _DOOR_SEED:
        width_ft = getattr(mod, width_const, 3.0)
        width_in = round(width_ft * 12.0, 2)
        conn.execute(
            "INSERT OR IGNORE INTO doors "
            "(opening_name, width, hinge_side, swing_direction, door_type) "
            "VALUES (?, ?, ?, ?, ?)",
            (opening_name, width_in, hinge, swing, dtype),
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


def get_constant_value(name, db_path=None):
    """Return the current value of a single constant, or None if not found."""
    with get_db(db_path) as conn:
        row = conn.execute(
            "SELECT value FROM constants WHERE name = ?", (name,)
        ).fetchone()
        return row["value"] if row else None


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


def get_outline_chain_row(seq, db_path=None):
    """Return a single outline chain row by seq, or None."""
    with get_db(db_path) as conn:
        row = conn.execute(
            "SELECT * FROM outline_chain WHERE seq = ?", (seq,)
        ).fetchone()
        return dict(row) if row else None


def update_outline_segment(seq, updates, db_path=None):
    """Update outline chain segment fields.  Returns updated row or None."""
    allowed = {"distance", "radius", "sweep", "sweep_name", "seg_type",
               "center_name", "n_pts", "end_name", "bearing_flex"}
    sets = []
    vals = []
    for k, v in updates.items():
        if k in allowed:
            sets.append(f"{k} = ?")
            vals.append(v)
    if not sets:
        return get_outline_chain_row(seq, db_path)
    vals.append(seq)
    with get_db(db_path) as conn:
        conn.execute(
            f"UPDATE outline_chain SET {', '.join(sets)} WHERE seq = ?", vals
        )
    return get_outline_chain_row(seq, db_path)


def insert_outline_segment(seq, row_data, db_path=None):
    """Insert a new outline chain segment at seq.

    All segments with seq >= target are shifted up by 1.
    Returns the inserted row.
    """
    with get_db(db_path) as conn:
        # Negate affected seqs to avoid UNIQUE collision, then flip back +1
        conn.execute(
            "UPDATE outline_chain SET seq = -(seq + 1) WHERE seq >= ?", (seq,)
        )
        conn.execute(
            "UPDATE outline_chain SET seq = -seq WHERE seq < 0"
        )
        conn.execute(
            "INSERT INTO outline_chain "
            "(seq, seg_type, distance, radius, sweep_name, sweep, "
            "center_name, n_pts, end_name, flex, bearing_flex) "
            "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
            (seq, row_data["seg_type"], row_data.get("distance"),
             row_data.get("radius"), row_data.get("sweep_name"),
             row_data.get("sweep"), row_data.get("center_name"),
             row_data.get("n_pts", 60), row_data["end_name"],
             row_data.get("flex"), row_data.get("bearing_flex", 0)),
        )
    return get_outline_chain_row(seq, db_path)


def delete_outline_segment(seq, db_path=None):
    """Delete an outline chain segment and renumber.

    Returns the deleted row dict (for undo) or None.
    """
    with get_db(db_path) as conn:
        row = conn.execute(
            "SELECT * FROM outline_chain WHERE seq = ?", (seq,)
        ).fetchone()
        if not row:
            return None
        deleted = dict(row)
        conn.execute("DELETE FROM outline_chain WHERE seq = ?", (seq,))
        # Negate affected seqs to avoid UNIQUE collision, then flip back -1
        conn.execute(
            "UPDATE outline_chain SET seq = -seq WHERE seq > ?", (seq,)
        )
        conn.execute(
            "UPDATE outline_chain SET seq = (-seq) - 1 WHERE seq < 0"
        )
    return deleted


def reset_outline_chain(db_path=None):
    """Reset outline chain to seed values from floorplan/geometry.py."""
    db_path = db_path or DB_PATH
    with get_db(db_path) as conn:
        conn.execute("DELETE FROM outline_chain")
        _seed_outline_chain(conn)


def restore_outline_chain(snapshot, db_path=None):
    """Restore outline_chain table from a full snapshot (for undo/rollback)."""
    with get_db(db_path) as conn:
        conn.execute("DELETE FROM outline_chain")
        for row in snapshot:
            conn.execute(
                "INSERT INTO outline_chain "
                "(seq, seg_type, distance, radius, sweep_name, sweep, "
                "center_name, n_pts, end_name, flex, bearing_flex) "
                "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
                (row["seq"], row["seg_type"], row.get("distance"),
                 row.get("radius"), row.get("sweep_name"),
                 row.get("sweep"), row.get("center_name"),
                 row.get("n_pts", 60), row["end_name"],
                 row.get("flex"), row.get("bearing_flex", 0)),
            )


def get_views(db_path=None):
    """Return all registered views."""
    with get_db(db_path) as conn:
        rows = conn.execute("SELECT * FROM views WHERE enabled = 1 ORDER BY category, name").fetchall()
        return [dict(r) for r in rows]


# ---------------------------------------------------------------------------
# CRUD: config
# ---------------------------------------------------------------------------

def get_config(key, db_path=None):
    """Return a config value by key, or None if not set."""
    with get_db(db_path) as conn:
        row = conn.execute("SELECT value FROM config WHERE key = ?", (key,)).fetchone()
        return row["value"] if row else None


def set_config(key, value, db_path=None):
    """Set a config value (insert or update)."""
    with get_db(db_path) as conn:
        conn.execute(
            "INSERT INTO config (key, value) VALUES (?, ?) "
            "ON CONFLICT(key) DO UPDATE SET value = ?",
            (key, str(value), str(value)),
        )


def get_all_config(db_path=None):
    """Return all config as a dict."""
    with get_db(db_path) as conn:
        rows = conn.execute("SELECT key, value FROM config ORDER BY key").fetchall()
        return {r["key"]: r["value"] for r in rows}


def get_outline_anchor_pivot(db_path=None):
    """Return (anchor_name, pivot_name) from config, or (None, None)."""
    anchor = get_config("outline_anchor", db_path)
    pivot = get_config("outline_pivot", db_path)
    return anchor, pivot


def get_outline_anchor_pos(db_path=None):
    """Return (E, N, brg) of stored anchor start position, or None."""
    E = get_config("outline_anchor_E", db_path)
    N = get_config("outline_anchor_N", db_path)
    brg = get_config("outline_anchor_brg", db_path)
    if E is None or N is None or brg is None:
        return None
    return (float(E), float(N), float(brg))


def set_outline_anchor_pivot(anchor, pivot, anchor_E, anchor_N, anchor_brg,
                             db_path=None, user_set=False):
    """Store anchor/pivot names and anchor absolute position in config.

    user_set=True marks this as an explicit user choice (shown in UI as active).
    Defaults seeded by _seed_default_anchor_pivot pass user_set=False so the
    UI correctly shows the 'Set anchor/pivot' flow rather than 'Clear'.
    """
    set_config("outline_anchor", anchor, db_path)
    set_config("outline_pivot", pivot, db_path)
    set_config("outline_anchor_E", str(anchor_E), db_path)
    set_config("outline_anchor_N", str(anchor_N), db_path)
    set_config("outline_anchor_brg", str(anchor_brg), db_path)
    set_config("outline_pivot_user_set", "1" if user_set else "0", db_path)


def get_outline_pivot_user_set(db_path=None):
    """Return True if the anchor/pivot was explicitly set by the user."""
    v = get_config("outline_pivot_user_set", db_path)
    return v == "1"


def _seed_default_anchor_pivot(db_path):
    """Seed anchor/pivot defaults from constants table.

    Anchor = end_name of last chain segment (chain wrap point).
    Pivot  = end_name of middle chain segment.
    Anchor coords = (F2_EASTING, F2_NORTHING + CORNER_SW_R, 0.0).
    """
    import floorplan.constants as fc
    chain_rows = get_outline_chain(db_path)
    if not chain_rows:
        return
    n = len(chain_rows)
    anchor_name = chain_rows[-1]["end_name"]
    pivot_name = chain_rows[n // 2]["end_name"]
    R_a1 = get_constant_value("CORNER_SW_R", db_path) or fc.CORNER_SW_R
    anchor_E = get_constant_value("F2_EASTING", db_path) or -18.5
    anchor_N = get_constant_value("F2_NORTHING", db_path) or -13.5
    # F2_NORTHING is stored before the R_a1 offset (same convention as before)
    anchor_N = float(anchor_N) + float(R_a1)
    set_outline_anchor_pivot(anchor_name, pivot_name,
                             float(anchor_E), anchor_N, 0.0, db_path)


def clear_outline_pivot(db_path=None):
    """Reset anchor/pivot to chain-default positions.

    Pivot is always active; 'clearing' means restoring defaults, not removing.
    Flex columns are also reset to auto-assign defaults.
    """
    db_path = db_path or DB_PATH
    with get_db(db_path) as conn:
        conn.execute("DELETE FROM config WHERE key IN "
                     "('outline_anchor', 'outline_pivot', "
                     "'outline_anchor_E', 'outline_anchor_N', "
                     "'outline_anchor_brg', 'outline_pivot_user_set')")
        conn.execute("UPDATE outline_chain SET flex = NULL")
    _seed_default_anchor_pivot(db_path)


def get_shapes(db_path=None):
    """Return all shape definitions as a list of dicts."""
    with get_db(db_path) as conn:
        rows = conn.execute("SELECT * FROM shapes ORDER BY name").fetchall()
        return [dict(r) for r in rows]


def get_shape(name, db_path=None):
    """Return a single shape definition by name, or None."""
    with get_db(db_path) as conn:
        row = conn.execute("SELECT * FROM shapes WHERE name = ?", (name,)).fetchone()
        return dict(row) if row else None


def create_shape(name, poly_json, scale=1.0, origin="center",
                 width_key=None, depth_key=None, description="", db_path=None):
    """Create a new shape. Returns the created record dict."""

    with get_db(db_path) as conn:
        conn.execute(
            "INSERT INTO shapes (name, poly_json, scale, origin, width_key, depth_key, description) "
            "VALUES (?, ?, ?, ?, ?, ?, ?)",
            (name, json.dumps(poly_json) if isinstance(poly_json, list) else poly_json,
             scale, origin, width_key, depth_key, description),
        )
    return get_shape(name, db_path)


def update_shape(name, **kwargs):
    """Update a shape by name. Accepted kwargs: poly_json, scale, origin,
    width_key, depth_key, description. Returns True if updated."""

    db_path = kwargs.pop("db_path", None)
    sets, vals = [], []
    for key in ("poly_json", "scale", "origin", "width_key", "depth_key", "description"):
        if key in kwargs:
            v = kwargs[key]
            if key == "poly_json" and isinstance(v, list):
                v = json.dumps(v)
            sets.append(f"{key} = ?")
            vals.append(v)
    if not sets:
        return False
    vals.append(name)
    with get_db(db_path) as conn:
        cur = conn.execute(f"UPDATE shapes SET {', '.join(sets)} WHERE name = ?", vals)
        return cur.rowcount > 0


def delete_shape(name, db_path=None):
    """Delete a shape by name. Returns True if deleted."""
    with get_db(db_path) as conn:
        cur = conn.execute("DELETE FROM shapes WHERE name = ?", (name,))
        return cur.rowcount > 0


def get_variant_exclusions(variant, db_path=None):
    """Return excluded element names for a variant, grouped by type.

    Returns dict like: {"wall": {"IW6"}, "rough_opening": {"RO5"}}
    """
    with get_db(db_path) as conn:
        try:
            rows = conn.execute(
                "SELECT element_type, element_name FROM variant_exclusions "
                "WHERE variant = ?", (variant,)
            ).fetchall()
        except sqlite3.OperationalError:
            return {}
    result = {}
    for r in rows:
        result.setdefault(r["element_type"], set()).add(r["element_name"])
    return result


# ---------------------------------------------------------------------------
# Variants CRUD
# ---------------------------------------------------------------------------

def _variant_row_to_dict(row):
    """Convert a variant DB row to a dict with parsed JSON fields."""
    return {
        "id": row["id"],
        "name": row["name"],
        "label": row["label"],
        "source_variant": row["source_variant"],
        "flags": json.loads(row["flags"]) if row["flags"] else {},
        "layer_config": json.loads(row["layer_config"]) if row["layer_config"] else {},
        "is_builtin": bool(row["is_builtin"]),
    }


def get_variants(db_path=None):
    """Return list of all variant dicts."""
    with get_db(db_path) as conn:
        rows = conn.execute(
            "SELECT id, name, label, source_variant, flags, layer_config, is_builtin "
            "FROM variants ORDER BY id"
        ).fetchall()
    return [_variant_row_to_dict(r) for r in rows]


def get_variant(name, db_path=None):
    """Return a single variant dict by name, or None."""
    with get_db(db_path) as conn:
        row = conn.execute(
            "SELECT id, name, label, source_variant, flags, layer_config, is_builtin "
            "FROM variants WHERE name = ?", (name,)
        ).fetchone()
    return _variant_row_to_dict(row) if row else None


def get_variant_by_id(variant_id, db_path=None):
    """Return a single variant dict by id, or None."""
    with get_db(db_path) as conn:
        row = conn.execute(
            "SELECT id, name, label, source_variant, flags, layer_config, is_builtin "
            "FROM variants WHERE id = ?", (variant_id,)
        ).fetchone()
    return _variant_row_to_dict(row) if row else None


def update_variant(variant_id, updates, db_path=None):
    """Update a variant record. Returns updated dict or None."""
    allowed = {"label", "flags", "layer_config"}
    with get_db(db_path) as conn:
        row = conn.execute(
            "SELECT id FROM variants WHERE id = ?", (variant_id,)
        ).fetchone()
        if not row:
            return None
        sets = []
        vals = []
        for k, v in updates.items():
            if k not in allowed:
                continue
            if k in ("flags", "layer_config") and not isinstance(v, str):
                v = json.dumps(v)
            sets.append(f"{k} = ?")
            vals.append(v)
        if sets:
            vals.append(variant_id)
            conn.execute(
                f"UPDATE variants SET {', '.join(sets)} WHERE id = ?",
                vals,
            )
        row = conn.execute(
            "SELECT id, name, label, source_variant, flags, layer_config, is_builtin "
            "FROM variants WHERE id = ?", (variant_id,)
        ).fetchone()
    return _variant_row_to_dict(row) if row else None


def create_variant(name, label, source_variant, flags, db_path=None):
    """Create a new user-defined variant. Returns the created variant dict."""
    if isinstance(flags, dict):
        flags = json.dumps(flags)
    with get_db(db_path) as conn:
        conn.execute(
            "INSERT INTO variants (name, label, source_variant, flags, is_builtin) "
            "VALUES (?, ?, ?, ?, 0)",
            (name, label, source_variant, flags),
        )
        row = conn.execute(
            "SELECT id, name, label, source_variant, flags, layer_config, is_builtin "
            "FROM variants WHERE name = ?", (name,)
        ).fetchone()
    return _variant_row_to_dict(row) if row else None


def create_variant_raw(record, db_path=None):
    """Re-insert a variant record (for undo). Restores id if possible."""
    flags = record.get("flags", {})
    if isinstance(flags, dict):
        flags = json.dumps(flags)
    layer_config = record.get("layer_config", {})
    if isinstance(layer_config, dict):
        layer_config = json.dumps(layer_config)
    with get_db(db_path) as conn:
        conn.execute(
            "INSERT INTO variants (id, name, label, source_variant, flags, layer_config, is_builtin) "
            "VALUES (?, ?, ?, ?, ?, ?, ?)",
            (record["id"], record["name"], record["label"],
             record.get("source_variant"), flags, layer_config,
             1 if record.get("is_builtin") else 0),
        )


def delete_variant(variant_id, db_path=None):
    """Delete a variant by id. Returns deleted id or None."""
    with get_db(db_path) as conn:
        row = conn.execute(
            "SELECT id FROM variants WHERE id = ?", (variant_id,)
        ).fetchone()
        if not row:
            return None
        conn.execute("DELETE FROM variants WHERE id = ?", (variant_id,))
    return variant_id


def clone_variant_exclusions(source, target, db_path=None):
    """Copy variant_exclusions from source variant to target variant."""
    with get_db(db_path) as conn:
        rows = conn.execute(
            "SELECT element_type, element_name FROM variant_exclusions "
            "WHERE variant = ?", (source,)
        ).fetchall()
        for r in rows:
            conn.execute(
                "INSERT OR IGNORE INTO variant_exclusions "
                "(variant, element_type, element_name) VALUES (?, ?, ?)",
                (target, r["element_type"], r["element_name"]),
            )


def delete_variant_exclusions(variant, db_path=None):
    """Remove all variant_exclusions for a variant."""
    with get_db(db_path) as conn:
        conn.execute(
            "DELETE FROM variant_exclusions WHERE variant = ?", (variant,)
        )


def set_variant_exclusion(variant, element_type, element_name,
                          excluded, db_path=None):
    """Add or remove a single variant exclusion.

    Returns True if state changed, False if already in desired state.
    """
    with get_db(db_path) as conn:
        exists = conn.execute(
            "SELECT 1 FROM variant_exclusions "
            "WHERE variant = ? AND element_type = ? AND element_name = ?",
            (variant, element_type, element_name),
        ).fetchone()
        if excluded and not exists:
            conn.execute(
                "INSERT INTO variant_exclusions "
                "(variant, element_type, element_name) VALUES (?, ?, ?)",
                (variant, element_type, element_name),
            )
            return True
        if not excluded and exists:
            conn.execute(
                "DELETE FROM variant_exclusions "
                "WHERE variant = ? AND element_type = ? AND element_name = ?",
                (variant, element_type, element_name),
            )
            return True
        return False


def _parse_props(raw):
    """Parse properties JSON, handling double-encoding."""
    if not raw:
        return {}
    parsed = json.loads(raw)
    if isinstance(parsed, str):
        parsed = json.loads(parsed)
    return parsed if isinstance(parsed, dict) else {}


def clone_variant_elements(source, target, db_path=None):
    """Add target variant to properties.variants lists for elements visible in source.

    For elements with properties.variants containing source, add target.
    For elements with variant column == source, convert to properties.variants list.
    For elements with variant IS NULL and no properties.variants, they're already
    visible in all variants — no change needed.
    """
    with get_db(db_path) as conn:
        rows = conn.execute(
            "SELECT id, properties, variant FROM elements"
        ).fetchall()
        for r in rows:
            props = _parse_props(r["properties"])
            variants_list = props.get("variants")
            if variants_list is not None:
                # Explicit list: add target if source is in the list
                if source in variants_list and target not in variants_list:
                    variants_list.append(target)
                    props["variants"] = variants_list
                    conn.execute(
                        "UPDATE elements SET properties = ? WHERE id = ?",
                        (json.dumps(props), r["id"]),
                    )
            elif r["variant"] == source:
                # Single variant column: convert to properties.variants list
                props["variants"] = [source, target]
                conn.execute(
                    "UPDATE elements SET properties = ?, variant = NULL WHERE id = ?",
                    (json.dumps(props), r["id"]),
                )
            # variant IS NULL and no properties.variants → visible everywhere, no change

        # Clone variant-specific formulas: copy formulas with variant=source
        # as variant=target so the cloned variant has the same element geometry.
        source_formulas = conn.execute(
            "SELECT element_name, param_name, formula_json, locked, locked_value "
            "FROM element_formulas WHERE variant = ?",
            (source,),
        ).fetchall()
        for f in source_formulas:
            # Check if target already has this formula
            existing = conn.execute(
                "SELECT id FROM element_formulas "
                "WHERE element_name = ? AND param_name = ? AND variant = ?",
                (f["element_name"], f["param_name"], target),
            ).fetchone()
            if not existing:
                conn.execute(
                    "INSERT INTO element_formulas "
                    "(element_name, param_name, formula_json, variant, locked, locked_value) "
                    "VALUES (?, ?, ?, ?, ?, ?)",
                    (f["element_name"], f["param_name"], f["formula_json"],
                     target, f["locked"], f["locked_value"]),
                )


def unclone_variant_elements(target, db_path=None):
    """Remove target variant from properties.variants lists on all elements."""
    with get_db(db_path) as conn:
        rows = conn.execute(
            "SELECT id, properties FROM elements"
        ).fetchall()
        for r in rows:
            props = _parse_props(r["properties"])
            variants_list = props.get("variants")
            if variants_list is not None and target in variants_list:
                variants_list.remove(target)
                props["variants"] = variants_list if variants_list else None
                conn.execute(
                    "UPDATE elements SET properties = ? WHERE id = ?",
                    (json.dumps(props), r["id"]),
                )


def get_room_label_offsets(db_path=None):
    """Return dict of room label offsets: {name: (offset_e, offset_n)}."""
    with get_db(db_path) as conn:
        try:
            rows = conn.execute(
                "SELECT room_name, offset_e, offset_n FROM room_label_offsets"
            ).fetchall()
        except sqlite3.OperationalError:
            return {}
    return {r["room_name"]: (r["offset_e"], r["offset_n"]) for r in rows}


def set_room_label_offset(room_name, offset_e, offset_n, db_path=None):
    """Set the label offset for a room (relative to centroid)."""
    with get_db(db_path) as conn:
        conn.execute(
            "INSERT INTO room_label_offsets (room_name, offset_e, offset_n) "
            "VALUES (?, ?, ?) "
            "ON CONFLICT(room_name) DO UPDATE SET offset_e=?, offset_n=?",
            (room_name, offset_e, offset_n, offset_e, offset_n),
        )


def validate_db(db_path=None):
    """Check database health. Returns (ok, issues) tuple.

    ok: True if DB is healthy
    issues: list of strings describing problems found
    """
    db_path = db_path or DB_PATH
    issues = []
    if not os.path.exists(db_path):
        return False, ["database file does not exist"]
    try:
        with get_db(db_path) as conn:
            # Check required tables exist
            rows = conn.execute(
                "SELECT name FROM sqlite_master WHERE type='table'"
            ).fetchall()
            existing = {r["name"] for r in rows}
            required = {"constants", "outline_chain", "views", "elements",
                        "doors", "config", "variants"}
            for tbl in sorted(required - existing):
                issues.append(f"{tbl} table missing")
            if issues:
                return False, issues
            # Check critical tables have data
            for tbl, label in [("constants", "constants"),
                               ("outline_chain", "outline chain"),
                               ("views", "views")]:
                count = conn.execute(f"SELECT COUNT(*) AS n FROM {tbl}").fetchone()["n"]
                if count == 0:
                    issues.append(f"{label} table is empty")
    except Exception as exc:
        issues.append(f"database error: {exc}")
    return (len(issues) == 0), issues


def reset_db(db_path=None):
    """Delete and recreate the database with fresh seed data.

    Returns (ok, issues) from validate_db on the new database.
    """
    db_path = db_path or DB_PATH
    if os.path.exists(db_path):
        os.remove(db_path)
    # Also remove WAL/SHM files if present
    for suffix in ("-wal", "-shm"):
        wal = db_path + suffix
        if os.path.exists(wal):
            os.remove(wal)
    init_db(db_path)
    return validate_db(db_path)


def reset_constants(db_path=None):
    """Reset all constants to their original values from source code."""
    db_path = db_path or DB_PATH
    with get_db(db_path) as conn:
        conn.execute("DELETE FROM constants")
        _seed_constants(conn)
        _seed_variant_item_constants(conn)


def reset_elements(db_path=None):
    """Reset elements and doors to seed state.

    Deletes all non-seeded elements (overrides, placed items, drawn walls,
    user labels, user dimensions) and re-seeds IW walls, doors, and room labels.
    """
    db_path = db_path or DB_PATH
    with get_db(db_path) as conn:
        conn.execute("DELETE FROM elements")
        conn.execute("DELETE FROM doors")
        _seed_elements(conn)
        _seed_doors(conn)
        from app.labels import seed_room_labels, seed_builtin_dimensions
        seed_room_labels(conn)
        seed_builtin_dimensions(conn)


def restore_elements(elements, doors, db_path=None):
    """Restore elements and doors tables from a snapshot (for undo/redo)."""

    db_path = db_path or DB_PATH
    with get_db(db_path) as conn:
        conn.execute("DELETE FROM elements")
        conn.execute("DELETE FROM doors")
        for e in elements:
            props = e.get("properties", "{}")
            if isinstance(props, dict):
                props = json.dumps(props)
            conn.execute(
                "INSERT INTO elements (id, type, name, properties, variant) "
                "VALUES (?, ?, ?, ?, ?)",
                (e["id"], e["type"], e["name"], props, e.get("variant")),
            )
        for d in doors:
            conn.execute(
                "INSERT INTO doors (opening_name, width, hinge_side, "
                "swing_direction, door_type) VALUES (?, ?, ?, ?, ?)",
                (d["opening_name"], d["width"], d["hinge_side"],
                 d["swing_direction"], d["door_type"]),
            )


# ---------------------------------------------------------------------------
# CRUD: elements
# ---------------------------------------------------------------------------

def get_all_elements(db_path=None):
    """Return all elements as a list of dicts."""
    with get_db(db_path) as conn:
        rows = conn.execute(
            "SELECT id, type, name, properties, variant FROM elements ORDER BY type, name"
        ).fetchall()
        return [dict(r) for r in rows]


def get_element(element_id, db_path=None):
    """Return a single element by ID, or None."""
    with get_db(db_path) as conn:
        row = conn.execute(
            "SELECT id, type, name, properties, variant FROM elements WHERE id = ?",
            (element_id,),
        ).fetchone()
        return dict(row) if row else None


def get_element_by_name(name, db_path=None):
    """Return a single element by name, or None."""
    with get_db(db_path) as conn:
        row = conn.execute(
            "SELECT id, type, name, properties, variant FROM elements WHERE name = ?",
            (name,),
        ).fetchone()
        return dict(row) if row else None


def create_element(type_, name, properties=None, variant=None, db_path=None):
    """Create a new element. Returns the created record dict with assigned id."""

    props = json.dumps(properties or {})
    with get_db(db_path) as conn:
        cur = conn.execute(
            "INSERT INTO elements (type, name, properties, variant) VALUES (?, ?, ?, ?)",
            (type_, name, props, variant),
        )
        return {
            "id": cur.lastrowid, "type": type_, "name": name,
            "properties": props, "variant": variant,
        }


def create_element_raw(record, db_path=None):
    """Re-insert an element from a full record dict (used by undo)."""
    with get_db(db_path) as conn:
        conn.execute(
            "INSERT INTO elements (id, type, name, properties, variant) VALUES (?, ?, ?, ?, ?)",
            (record["id"], record["type"], record["name"],
             record["properties"], record.get("variant")),
        )


def update_element(element_id, updates, db_path=None):
    """Update element fields.  Returns updated record or None if not found."""

    allowed = {"type", "name", "properties", "variant"}
    sets = []
    vals = []
    for k, v in updates.items():
        if k in allowed:
            if k == "properties" and isinstance(v, dict):
                v = json.dumps(v)
            sets.append(f"{k} = ?")
            vals.append(v)
    if not sets:
        return get_element(element_id, db_path)
    vals.append(element_id)
    with get_db(db_path) as conn:
        conn.execute(
            f"UPDATE elements SET {', '.join(sets)} WHERE id = ?", vals
        )
    return get_element(element_id, db_path)


def delete_element(element_id, db_path=None):
    """Delete an element by ID.  Returns list of deleted IDs (includes cascade)."""
    deleted = []
    with get_db(db_path) as conn:
        # Check for hosted openings if deleting a wall
        row = conn.execute("SELECT type, name FROM elements WHERE id = ?", (element_id,)).fetchone()
        if not row:
            return deleted
        if row["type"] == "wall":
            # Cascade: delete openings that reference this wall
            hosted = conn.execute(
                "SELECT id FROM elements WHERE type = 'opening' AND "
                "json_extract(properties, '$.host_wall') = ?",
                (row["name"],),
            ).fetchall()
            for h in hosted:
                # Also remove door for the opening
                opening_row = conn.execute(
                    "SELECT name FROM elements WHERE id = ?", (h["id"],)
                ).fetchone()
                if opening_row:
                    conn.execute(
                        "DELETE FROM doors WHERE opening_name = ?",
                        (opening_row["name"],),
                    )
                conn.execute("DELETE FROM elements WHERE id = ?", (h["id"],))
                deleted.append(h["id"])
        # Clean up door record if deleting an opening
        if row["type"] == "opening":
            conn.execute("DELETE FROM doors WHERE opening_name = ?",
                         (row["name"],))
        # Clean up formulas and deps for the deleted element
        conn.execute("DELETE FROM element_formulas WHERE element_name = ?",
                     (row["name"],))
        conn.execute("DELETE FROM formula_deps WHERE element_name = ?",
                     (row["name"],))
        conn.execute("DELETE FROM elements WHERE id = ?", (element_id,))
        deleted.append(element_id)
    return deleted


# ---------------------------------------------------------------------------
# CRUD: doors
# ---------------------------------------------------------------------------

def get_all_doors(db_path=None):
    """Return all doors as a list of dicts."""
    with get_db(db_path) as conn:
        rows = conn.execute(
            "SELECT id, opening_name, width, hinge_side, swing_direction, door_type "
            "FROM doors ORDER BY opening_name"
        ).fetchall()
        return [dict(r) for r in rows]


def get_door(opening_name, db_path=None):
    """Return a door by opening name, or None."""
    with get_db(db_path) as conn:
        row = conn.execute(
            "SELECT id, opening_name, width, hinge_side, swing_direction, door_type "
            "FROM doors WHERE opening_name = ?",
            (opening_name,),
        ).fetchone()
        return dict(row) if row else None


def create_door(opening_name, width, hinge_side, swing_direction, door_type="single", db_path=None):
    """Create a new door. Returns the created record dict."""
    with get_db(db_path) as conn:
        cur = conn.execute(
            "INSERT INTO doors (opening_name, width, hinge_side, swing_direction, door_type) "
            "VALUES (?, ?, ?, ?, ?)",
            (opening_name, float(width), hinge_side, swing_direction, door_type),
        )
        return {
            "id": cur.lastrowid, "opening_name": opening_name,
            "width": float(width), "hinge_side": hinge_side,
            "swing_direction": swing_direction, "door_type": door_type,
        }


def create_door_raw(record, db_path=None):
    """Re-insert a door from a full record dict (used by undo)."""
    with get_db(db_path) as conn:
        conn.execute(
            "INSERT INTO doors (id, opening_name, width, hinge_side, swing_direction, door_type) "
            "VALUES (?, ?, ?, ?, ?, ?)",
            (record["id"], record["opening_name"], record["width"],
             record["hinge_side"], record["swing_direction"], record["door_type"]),
        )


def update_door(opening_name, updates, db_path=None):
    """Update door fields.  Returns updated record or None if not found."""
    allowed = {"width", "hinge_side", "swing_direction", "door_type"}
    sets = []
    vals = []
    for k, v in updates.items():
        if k in allowed:
            if k == "width":
                v = float(v)
            sets.append(f"{k} = ?")
            vals.append(v)
    if not sets:
        return get_door(opening_name, db_path)
    vals.append(opening_name)
    with get_db(db_path) as conn:
        conn.execute(
            f"UPDATE doors SET {', '.join(sets)} WHERE opening_name = ?", vals
        )
    return get_door(opening_name, db_path)


def delete_door(opening_name, db_path=None):
    """Delete a door by opening name.  Returns True if deleted."""
    with get_db(db_path) as conn:
        cur = conn.execute(
            "DELETE FROM doors WHERE opening_name = ?", (opening_name,)
        )
        return cur.rowcount > 0


# ---------------------------------------------------------------------------
# CRUD: element_formulas
# ---------------------------------------------------------------------------

def get_element_formulas(element_name, variant=None, db_path=None):
    """Return all formulas for an element as list of dicts."""
    with get_db(db_path) as conn:
        if variant is None:
            rows = conn.execute(
                "SELECT * FROM element_formulas WHERE element_name = ? "
                "AND variant IS NULL ORDER BY param_name",
                (element_name,),
            ).fetchall()
        else:
            rows = conn.execute(
                "SELECT * FROM element_formulas WHERE element_name = ? "
                "AND (variant IS NULL OR variant = ?) ORDER BY param_name",
                (element_name, variant),
            ).fetchall()
        return [dict(r) for r in rows]


def get_all_formulas(variant=None, db_path=None):
    """Return all formulas, optionally filtered by variant."""
    with get_db(db_path) as conn:
        if variant is None:
            rows = conn.execute(
                "SELECT * FROM element_formulas ORDER BY element_name, param_name"
            ).fetchall()
        else:
            rows = conn.execute(
                "SELECT * FROM element_formulas "
                "WHERE variant IS NULL OR variant = ? "
                "ORDER BY element_name, param_name",
                (variant,),
            ).fetchall()
        return [dict(r) for r in rows]


def upsert_formula(element_name, param_name, formula_json, variant=None,
                   db_path=None):
    """Insert or update a formula.  Returns the formula row dict."""
    fj = json.dumps(formula_json) if not isinstance(formula_json, str) else formula_json
    with get_db(db_path) as conn:
        # SQLite treats NULLs as distinct in UNIQUE, so handle manually
        if variant is None:
            existing = conn.execute(
                "SELECT id FROM element_formulas "
                "WHERE element_name = ? AND param_name = ? AND variant IS NULL",
                (element_name, param_name),
            ).fetchone()
        else:
            existing = conn.execute(
                "SELECT id FROM element_formulas "
                "WHERE element_name = ? AND param_name = ? AND variant = ?",
                (element_name, param_name, variant),
            ).fetchone()
        if existing:
            conn.execute(
                "UPDATE element_formulas SET formula_json = ? WHERE id = ?",
                (fj, existing["id"]),
            )
        else:
            conn.execute(
                "INSERT INTO element_formulas "
                "(element_name, param_name, formula_json, variant) "
                "VALUES (?, ?, ?, ?)",
                (element_name, param_name, fj, variant),
            )
    formulas = get_element_formulas(element_name, variant, db_path)
    return next((f for f in formulas if f["param_name"] == param_name), None)


def delete_formula(element_name, param_name, variant=None, db_path=None):
    """Delete a formula.  Returns True if deleted."""
    with get_db(db_path) as conn:
        if variant is None:
            cur = conn.execute(
                "DELETE FROM element_formulas "
                "WHERE element_name = ? AND param_name = ? AND variant IS NULL",
                (element_name, param_name),
            )
        else:
            cur = conn.execute(
                "DELETE FROM element_formulas "
                "WHERE element_name = ? AND param_name = ? AND variant = ?",
                (element_name, param_name, variant),
            )
        return cur.rowcount > 0


def set_formula_lock(element_name, param_name, locked, locked_value=None,
                     variant=None, db_path=None):
    """Set lock state on a formula.  Returns True if updated."""
    lv = json.dumps(locked_value) if locked_value is not None else None
    with get_db(db_path) as conn:
        if variant is None:
            cur = conn.execute(
                "UPDATE element_formulas SET locked = ?, locked_value = ? "
                "WHERE element_name = ? AND param_name = ? AND variant IS NULL",
                (1 if locked else 0, lv, element_name, param_name),
            )
        else:
            cur = conn.execute(
                "UPDATE element_formulas SET locked = ?, locked_value = ? "
                "WHERE element_name = ? AND param_name = ? AND variant = ?",
                (1 if locked else 0, lv, element_name, param_name, variant),
            )
        return cur.rowcount > 0


# ---------------------------------------------------------------------------
# CRUD: formula_deps
# ---------------------------------------------------------------------------

def get_formula_deps(element_name, param_name=None, db_path=None):
    """Return dependencies for an element (optionally specific param)."""
    with get_db(db_path) as conn:
        if param_name:
            rows = conn.execute(
                "SELECT * FROM formula_deps "
                "WHERE element_name = ? AND param_name = ?",
                (element_name, param_name),
            ).fetchall()
        else:
            rows = conn.execute(
                "SELECT * FROM formula_deps WHERE element_name = ?",
                (element_name,),
            ).fetchall()
        return [dict(r) for r in rows]


def get_dependents(dep_name, dep_type=None, db_path=None):
    """Return elements that depend on dep_name."""
    with get_db(db_path) as conn:
        if dep_type:
            rows = conn.execute(
                "SELECT DISTINCT element_name, param_name FROM formula_deps "
                "WHERE dep_name = ? AND dep_type = ?",
                (dep_name, dep_type),
            ).fetchall()
        else:
            rows = conn.execute(
                "SELECT DISTINCT element_name, param_name FROM formula_deps "
                "WHERE dep_name = ?",
                (dep_name,),
            ).fetchall()
        return [dict(r) for r in rows]


def get_all_formula_deps(db_path=None):
    """Return all formula dependency edges."""
    with get_db(db_path) as conn:
        rows = conn.execute(
            "SELECT element_name, param_name, dep_type, dep_name "
            "FROM formula_deps ORDER BY element_name, param_name"
        ).fetchall()
        return [dict(r) for r in rows]


def rebuild_formula_deps(element_name, param_name, deps, db_path=None):
    """Replace all deps for a formula with a new set.

    deps: list of (dep_type, dep_name) tuples.
    """
    with get_db(db_path) as conn:
        conn.execute(
            "DELETE FROM formula_deps WHERE element_name = ? AND param_name = ?",
            (element_name, param_name),
        )
        for dep_type, dep_name in deps:
            conn.execute(
                "INSERT OR IGNORE INTO formula_deps "
                "(element_name, param_name, dep_type, dep_name) "
                "VALUES (?, ?, ?, ?)",
                (element_name, param_name, dep_type, dep_name),
            )


# ---------------------------------------------------------------------------
# Catalog items CRUD
# ---------------------------------------------------------------------------

def get_all_catalog_items(db_path=None):
    """Return all catalog items as list of dicts."""
    with get_db(db_path) as conn:
        try:
            rows = conn.execute(
                "SELECT * FROM catalog_items ORDER BY item_type, label, key"
            ).fetchall()
        except Exception:
            # Table missing — create and seed it
            _create_catalog_table(conn)
            _seed_catalog_items(conn, db_path)
            conn.commit()
            rows = conn.execute(
                "SELECT * FROM catalog_items ORDER BY item_type, label, key"
            ).fetchall()
        return [dict(r) for r in rows]


def get_catalog_item(key, db_path=None):
    """Return a single catalog item by key, or None."""
    with get_db(db_path) as conn:
        row = conn.execute(
            "SELECT * FROM catalog_items WHERE key = ?", (key,)
        ).fetchone()
        return dict(row) if row else None


def create_catalog_item(data, db_path=None):
    """Insert a new catalog item. Returns the created record."""
    with get_db(db_path) as conn:
        conn.execute(
            "INSERT INTO catalog_items "
            "(key, item_type, label, width, depth, radius, shape, "
            " door, clearance, product_url, variants, stacked, clip_to_inner) "
            "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
            (data["key"], data["item_type"], data["label"],
             data.get("width"), data.get("depth"), data.get("radius"),
             data.get("shape", "rect"),
             data.get("door"), data.get("clearance"),
             data.get("product_url"), data.get("variants", "[]"),
             data.get("stacked", 0), data.get("clip_to_inner", 0)),
        )
        return get_catalog_item(data["key"], db_path)


def update_catalog_item(key, data, db_path=None):
    """Update specific fields of a catalog item."""
    allowed = {"label", "width", "depth", "radius", "shape", "door",
               "clearance", "product_url", "variants", "stacked",
               "clip_to_inner", "item_type"}
    sets = []
    vals = []
    for k, v in data.items():
        if k in allowed:
            sets.append(f"{k} = ?")
            vals.append(v)
    if not sets:
        return
    vals.append(key)
    with get_db(db_path) as conn:
        conn.execute(
            f"UPDATE catalog_items SET {', '.join(sets)} WHERE key = ?",
            vals,
        )


def delete_catalog_item(key, db_path=None):
    """Delete a catalog item by key."""
    with get_db(db_path) as conn:
        conn.execute("DELETE FROM catalog_items WHERE key = ?", (key,))


def ensure_catalog_item(data, db_path=None):
    """Insert a catalog item if it doesn't already exist (INSERT OR IGNORE)."""
    with get_db(db_path) as conn:
        conn.execute(
            "INSERT OR IGNORE INTO catalog_items "
            "(key, item_type, label, width, depth, radius, shape, "
            " door, clearance, product_url, variants, stacked, clip_to_inner) "
            "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
            (data["key"], data["item_type"], data["label"],
             data.get("width"), data.get("depth"), data.get("radius"),
             data.get("shape", "rect"),
             data.get("door"), data.get("clearance"),
             data.get("product_url"), data.get("variants", "[]"),
             data.get("stacked", 0), data.get("clip_to_inner", 0)),
        )


def seed_iw_formulas(conn):
    """Seed interior wall + layout item + variant item formulas into element_formulas + formula_deps.

    Idempotent: skips elements that already have a position formula for the
    given variant.
    """
    from app.evaluator import (get_iw_formulas, get_layout_item_formulas,
                               get_outer_opening_formulas,
                               get_rough_opening_formulas,
                               get_variant_item_formulas, extract_deps)

    # Non-variant formulas (variant=NULL): IW walls, layout items, openings
    formulas = {**get_iw_formulas(), **get_layout_item_formulas(),
                **get_outer_opening_formulas(), **get_rough_opening_formulas()}
    for wall_name, formula in formulas.items():
        existing = conn.execute(
            "SELECT id FROM element_formulas "
            "WHERE element_name = ? AND param_name = 'position' AND variant IS NULL",
            (wall_name,),
        ).fetchone()
        if existing:
            continue
        fj = json.dumps(formula)
        conn.execute(
            "INSERT INTO element_formulas "
            "(element_name, param_name, formula_json, variant) "
            "VALUES (?, 'position', ?, NULL)",
            (wall_name, fj),
        )
        deps = extract_deps(formula)
        for dep_type, dep_name in deps:
            conn.execute(
                "INSERT OR IGNORE INTO formula_deps "
                "(element_name, param_name, dep_type, dep_name) "
                "VALUES (?, 'position', ?, ?)",
                (wall_name, dep_type, dep_name),
            )

    # Variant item formulas (Phase 12e)
    for elem_name, variant, formula in get_variant_item_formulas():
        if variant is None:
            existing = conn.execute(
                "SELECT id FROM element_formulas "
                "WHERE element_name = ? AND param_name = 'position' AND variant IS NULL",
                (elem_name,),
            ).fetchone()
        else:
            existing = conn.execute(
                "SELECT id FROM element_formulas "
                "WHERE element_name = ? AND param_name = 'position' AND variant = ?",
                (elem_name, variant),
            ).fetchone()
        if existing:
            continue
        fj = json.dumps(formula)
        conn.execute(
            "INSERT INTO element_formulas "
            "(element_name, param_name, formula_json, variant) "
            "VALUES (?, 'position', ?, ?)",
            (elem_name, fj, variant),
        )
        deps = extract_deps(formula)
        for dep_type, dep_name in deps:
            conn.execute(
                "INSERT OR IGNORE INTO formula_deps "
                "(element_name, param_name, dep_type, dep_name) "
                "VALUES (?, 'position', ?, ?)",
                (elem_name, dep_type, dep_name),
            )


# ---------------------------------------------------------------------------
# Survey CRUD (Phase 14-B)
# ---------------------------------------------------------------------------

def get_survey_legs(db_path=None):
    """Return all survey legs ordered by seq."""
    with get_db(db_path or DB_PATH) as conn:
        rows = conn.execute(
            "SELECT seq, bearing_deg, bearing_min, bearing_sec, "
            "distance_ft, distance_inch, label "
            "FROM survey_legs ORDER BY seq"
        ).fetchall()
        return [dict(r) for r in rows]


def get_survey_config(db_path=None):
    """Return survey config as a dict of key→value."""
    with get_db(db_path or DB_PATH) as conn:
        rows = conn.execute("SELECT key, value FROM survey_config").fetchall()
        return {r["key"]: json.loads(r["value"]) for r in rows}


def update_survey_leg(seq, data, db_path=None):
    """Update a survey leg. data is a dict with optional keys:
    bearing_deg, bearing_min, bearing_sec, distance_ft, distance_inch, label."""
    allowed = {"bearing_deg", "bearing_min", "bearing_sec",
               "distance_ft", "distance_inch", "label"}
    updates = {k: v for k, v in data.items() if k in allowed}
    if not updates:
        return
    cols = ", ".join(f"{k} = ?" for k in updates)
    vals = list(updates.values()) + [seq]
    with get_db(db_path or DB_PATH) as conn:
        conn.execute(f"UPDATE survey_legs SET {cols} WHERE seq = ?", vals)


def update_survey_config(key, value, db_path=None):
    """Update a survey config value."""
    with get_db(db_path or DB_PATH) as conn:
        conn.execute(
            "INSERT OR REPLACE INTO survey_config (key, value) VALUES (?, ?)",
            (key, json.dumps(value)),
        )


def reset_survey(db_path=None):
    """Reset survey data to defaults."""
    with get_db(db_path or DB_PATH) as conn:
        conn.execute("DELETE FROM survey_legs")
        conn.execute("DELETE FROM survey_config")
        _seed_survey(conn)


def restore_survey(legs, config, db_path=None):
    """Restore survey legs and config from snapshot (for undo/redo)."""
    with get_db(db_path or DB_PATH) as conn:
        conn.execute("DELETE FROM survey_legs")
        conn.execute("DELETE FROM survey_config")
        for leg in legs:
            conn.execute(
                "INSERT INTO survey_legs "
                "(seq, bearing_deg, bearing_min, bearing_sec, "
                "distance_ft, distance_inch, label) "
                "VALUES (?, ?, ?, ?, ?, ?, ?)",
                (leg["seq"], leg["bearing_deg"], leg["bearing_min"],
                 leg["bearing_sec"], leg["distance_ft"],
                 leg["distance_inch"], leg.get("label")),
            )
        for key, value in config.items():
            conn.execute(
                "INSERT INTO survey_config (key, value) VALUES (?, ?)",
                (key, json.dumps(value)),
            )


# ---------------------------------------------------------------------------
# Project Export / Import (Phase 14-D)
# ---------------------------------------------------------------------------

def export_project(db_path=None):
    """Export the entire project state as a JSON-serialisable dict."""
    import datetime
    from app.plumbing import get_plumbing_elements

    db_path = db_path or DB_PATH
    with get_db(db_path) as conn:
        constants = [dict(r) for r in conn.execute(
            "SELECT name, value, expr, unit, category, description "
            "FROM constants ORDER BY name"
        ).fetchall()]

        outline_chain = [dict(r) for r in conn.execute(
            "SELECT seq, seg_type, distance, radius, sweep_name, sweep, "
            "center_name, n_pts, end_name, flex FROM outline_chain ORDER BY seq"
        ).fetchall()]

        elements = [dict(r) for r in conn.execute(
            "SELECT id, type, name, properties, variant "
            "FROM elements ORDER BY id"
        ).fetchall()]

        element_formulas = [dict(r) for r in conn.execute(
            "SELECT element_name, param_name, formula_json, locked, "
            "locked_value, variant FROM element_formulas ORDER BY id"
        ).fetchall()]

        formula_deps = [dict(r) for r in conn.execute(
            "SELECT element_name, param_name, dep_type, dep_name "
            "FROM formula_deps ORDER BY element_name, param_name"
        ).fetchall()]

        doors = [dict(r) for r in conn.execute(
            "SELECT opening_name, width, hinge_side, swing_direction, door_type "
            "FROM doors ORDER BY opening_name"
        ).fetchall()]

        variants = [dict(r) for r in conn.execute(
            "SELECT name, label, source_variant, flags, layer_config, is_builtin "
            "FROM variants ORDER BY id"
        ).fetchall()]

        variant_exclusions = [dict(r) for r in conn.execute(
            "SELECT variant, element_type, element_name "
            "FROM variant_exclusions ORDER BY variant, element_type, element_name"
        ).fetchall()]

        survey_legs = [dict(r) for r in conn.execute(
            "SELECT seq, bearing_deg, bearing_min, bearing_sec, "
            "distance_ft, distance_inch, label "
            "FROM survey_legs ORDER BY seq"
        ).fetchall()]

        survey_config = {}
        for r in conn.execute("SELECT key, value FROM survey_config").fetchall():
            survey_config[r["key"]] = json.loads(r["value"])

        outline_anchor = None
        outline_pivot = None
        for r in conn.execute(
                "SELECT key, value FROM config "
                "WHERE key IN ('outline_anchor', 'outline_pivot')").fetchall():
            if r["key"] == "outline_anchor":
                outline_anchor = r["value"]
            elif r["key"] == "outline_pivot":
                outline_pivot = r["value"]

    plumbing = get_plumbing_elements(db_path)

    with get_db(db_path) as conn:
        inner_wall_overrides = [dict(r) for r in conn.execute(
            f"SELECT {_IW_OV_COLS} FROM inner_wall_overrides "
            "ORDER BY seg_index, sub_seq"
        ).fetchall()]

    return {
        "version": 1,
        "exported_at": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "constants": constants,
        "outline_chain": outline_chain,
        "elements": elements,
        "element_formulas": element_formulas,
        "formula_deps": formula_deps,
        "doors": doors,
        "variants": variants,
        "variant_exclusions": variant_exclusions,
        "survey_legs": survey_legs,
        "survey_config": survey_config,
        "inner_wall_overrides": inner_wall_overrides,
        "plumbing_elements": plumbing,
        "outline_anchor": outline_anchor,
        "outline_pivot": outline_pivot,
    }


def import_project(data, db_path=None):
    """Import a project from an exported dict, replacing all data.

    Validates structure before importing. Raises ValueError on invalid data.
    """
    db_path = db_path or DB_PATH

    # Validate required keys
    required = {"version", "constants", "outline_chain", "elements",
                "element_formulas", "formula_deps", "doors", "variants",
                "variant_exclusions", "survey_legs", "survey_config"}
    missing = required - set(data.keys())
    if missing:
        raise ValueError(f"Missing required keys: {missing}")

    if data["version"] != 1:
        raise ValueError(f"Unsupported export version: {data['version']}")

    # Validate outline closure
    if data["outline_chain"]:
        from app.outline_solver import db_rows_to_chain, solve_closure
        chain = db_rows_to_chain(data["outline_chain"])
        # Find CORNER_SW_R from imported constants
        cdict = {c["name"]: c["value"] for c in data["constants"]}
        r_a1 = cdict.get("CORNER_SW_R", 10.0 / 12.0 + (8.0 / 12.0 - 8.0 / 12.0))
        result = solve_closure(chain, r_a1)
        if not result.valid:
            raise ValueError(
                f"Imported outline chain does not close: "
                f"error={result.closure_error:.6f}")

    # Validate formula DAG (no cycles)
    if data["element_formulas"]:
        from app.evaluator import extract_deps
        # Build adjacency from formula_deps
        dep_graph = {}
        for d in data.get("formula_deps", []):
            if d["dep_type"] == "element":
                dep_graph.setdefault(d["element_name"], set()).add(d["dep_name"])
        # Topological sort check
        visited = set()
        in_stack = set()
        def _has_cycle(node):
            if node in in_stack:
                return True
            if node in visited:
                return False
            visited.add(node)
            in_stack.add(node)
            for dep in dep_graph.get(node, ()):
                if _has_cycle(dep):
                    return True
            in_stack.discard(node)
            return False
        for node in dep_graph:
            if _has_cycle(node):
                raise ValueError(f"Formula dependency cycle detected at {node}")

    with get_db(db_path) as conn:
        # Clear all tables
        for table in ("constants", "outline_chain", "elements", "doors",
                      "element_formulas", "formula_deps", "variants",
                      "variant_exclusions", "survey_legs", "survey_config",
                      "inner_wall_overrides", "plumbing_elements",
                      "room_label_offsets"):
            conn.execute(f"DELETE FROM {table}")

        # Import constants
        for c in data["constants"]:
            conn.execute(
                "INSERT INTO constants (name, value, expr, unit, category, description) "
                "VALUES (?, ?, ?, ?, ?, ?)",
                (c["name"], c["value"], c.get("expr"), c.get("unit", "ft"),
                 c.get("category", "misc"), c.get("description", "")),
            )

        # Import outline chain
        for row in data["outline_chain"]:
            conn.execute(
                "INSERT INTO outline_chain "
                "(seq, seg_type, distance, radius, sweep_name, sweep, "
                "center_name, n_pts, end_name, flex) "
                "VALUES (?,?,?,?,?,?,?,?,?,?)",
                (row["seq"], row["seg_type"], row.get("distance"),
                 row.get("radius"), row.get("sweep_name"), row.get("sweep"),
                 row.get("center_name"), row.get("n_pts", 60), row["end_name"],
                 row.get("flex")),
            )
        # If imported data has no flex designations, set defaults
        has_flex = any(r.get("flex") for r in data["outline_chain"])
        if not has_flex and data["outline_chain"]:
            n = len(data["outline_chain"])
            conn.execute(
                "UPDATE outline_chain SET flex = 'distance' WHERE seq = 0")
            conn.execute(
                "UPDATE outline_chain SET flex = 'distance' WHERE seq = ?",
                (n - 2,))
            conn.execute(
                "UPDATE outline_chain SET flex = 'sweep' WHERE seq = ?",
                (n - 1,))

        # Import elements
        for e in data["elements"]:
            props = e.get("properties", "{}")
            if isinstance(props, dict):
                props = json.dumps(props)
            conn.execute(
                "INSERT INTO elements (id, type, name, properties, variant) "
                "VALUES (?, ?, ?, ?, ?)",
                (e.get("id"), e["type"], e["name"], props, e.get("variant")),
            )

        # Import element formulas
        for f in data["element_formulas"]:
            fj = f["formula_json"]
            if isinstance(fj, dict):
                fj = json.dumps(fj)
            lv = f.get("locked_value")
            if isinstance(lv, (dict, list)):
                lv = json.dumps(lv)
            conn.execute(
                "INSERT INTO element_formulas "
                "(element_name, param_name, formula_json, locked, locked_value, variant) "
                "VALUES (?, ?, ?, ?, ?, ?)",
                (f["element_name"], f["param_name"], fj,
                 f.get("locked", 0), lv, f.get("variant")),
            )

        # Import formula deps
        for d in data["formula_deps"]:
            conn.execute(
                "INSERT OR IGNORE INTO formula_deps "
                "(element_name, param_name, dep_type, dep_name) "
                "VALUES (?, ?, ?, ?)",
                (d["element_name"], d["param_name"], d["dep_type"], d["dep_name"]),
            )

        # Import doors
        for d in data["doors"]:
            conn.execute(
                "INSERT INTO doors "
                "(opening_name, width, hinge_side, swing_direction, door_type) "
                "VALUES (?, ?, ?, ?, ?)",
                (d["opening_name"], d["width"], d["hinge_side"],
                 d["swing_direction"], d.get("door_type", "single")),
            )

        # Import variants
        for v in data["variants"]:
            flags = v.get("flags", "{}")
            if isinstance(flags, dict):
                flags = json.dumps(flags)
            lc = v.get("layer_config", "{}")
            if isinstance(lc, dict):
                lc = json.dumps(lc)
            conn.execute(
                "INSERT INTO variants "
                "(name, label, source_variant, flags, layer_config, is_builtin) "
                "VALUES (?, ?, ?, ?, ?, ?)",
                (v["name"], v["label"], v.get("source_variant"),
                 flags, lc, v.get("is_builtin", 0)),
            )

        # Import variant exclusions
        for ve in data["variant_exclusions"]:
            conn.execute(
                "INSERT OR IGNORE INTO variant_exclusions "
                "(variant, element_type, element_name) VALUES (?, ?, ?)",
                (ve["variant"], ve["element_type"], ve["element_name"]),
            )

        # Import survey legs
        for leg in data["survey_legs"]:
            conn.execute(
                "INSERT INTO survey_legs "
                "(seq, bearing_deg, bearing_min, bearing_sec, "
                "distance_ft, distance_inch, label) VALUES (?,?,?,?,?,?,?)",
                (leg["seq"], leg["bearing_deg"], leg["bearing_min"],
                 leg["bearing_sec"], leg["distance_ft"],
                 leg["distance_inch"], leg.get("label")),
            )

        # Import survey config
        for key, value in data.get("survey_config", {}).items():
            conn.execute(
                "INSERT INTO survey_config (key, value) VALUES (?, ?)",
                (key, json.dumps(value)),
            )

        # Import inner wall overrides
        for ov in data.get("inner_wall_overrides", []):
            conn.execute(
                "INSERT INTO inner_wall_overrides "
                f"({_IW_OV_COLS}) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)",
                (ov["seg_index"], ov.get("span_end"), ov["sub_seq"],
                 ov["seg_type"], ov.get("bearing"), ov.get("distance"),
                 ov.get("radius"), ov.get("sweep"),
                 ov.get("n_pts", 20)),
            )

        # Import plumbing elements
        for p in data.get("plumbing_elements", []):
            path = p.get("path", "[]")
            if isinstance(path, list):
                path = json.dumps(path)
            props = p.get("properties", "{}")
            if isinstance(props, dict):
                props = json.dumps(props)
            conn.execute(
                "INSERT INTO plumbing_elements "
                "(type, name, path, properties, fixture, sort_order) "
                "VALUES (?, ?, ?, ?, ?, ?)",
                (p["type"], p["name"], path, props, p.get("fixture"),
                 p.get("sort_order", 0)),
            )

        # Import outline anchor/pivot config
        conn.execute("DELETE FROM config WHERE key IN "
                     "('outline_anchor', 'outline_pivot')")
        if data.get("outline_anchor"):
            conn.execute(
                "INSERT INTO config (key, value) VALUES (?, ?)",
                ("outline_anchor", data["outline_anchor"]))
        if data.get("outline_pivot"):
            conn.execute(
                "INSERT INTO config (key, value) VALUES (?, ?)",
                ("outline_pivot", data["outline_pivot"]))


# ---------------------------------------------------------------------------
# Inner Wall Overrides CRUD (Phase 15½)
# ---------------------------------------------------------------------------

def get_inner_wall_overrides(db_path=None):
    """Return all inner wall overrides grouped by seg_index.

    Returns dict {seg_index: [sub-segment dicts ordered by sub_seq]}.
    Each sub-segment dict includes span_end (None for single-segment overrides).
    """
    with get_db(db_path or DB_PATH) as conn:
        rows = conn.execute(
            f"SELECT {_IW_OV_COLS} FROM inner_wall_overrides "
            "ORDER BY seg_index, sub_seq"
        ).fetchall()
    result = {}
    for r in rows:
        d = dict(r)
        si = d["seg_index"]
        result.setdefault(si, []).append(d)
    return result


def get_inner_wall_override(seg_index, db_path=None):
    """Return the override chain for a single segment, or empty list."""
    with get_db(db_path or DB_PATH) as conn:
        rows = conn.execute(
            f"SELECT {_IW_OV_COLS} FROM inner_wall_overrides "
            "WHERE seg_index = ? ORDER BY sub_seq",
            (seg_index,),
        ).fetchall()
        return [dict(r) for r in rows]


def upsert_inner_wall_override(seg_index, chain, span_end=None, db_path=None):
    """Set the override chain for a segment or span (replaces existing).

    When span_end is not None, the override spans inner_segs[seg_index]
    through inner_segs[span_end] inclusive.
    """
    with get_db(db_path or DB_PATH) as conn:
        conn.execute("DELETE FROM inner_wall_overrides WHERE seg_index = ?",
                     (seg_index,))
        for sub_seq, seg in enumerate(chain):
            conn.execute(
                "INSERT INTO inner_wall_overrides "
                f"({_IW_OV_COLS}) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)",
                (seg_index, span_end, sub_seq, seg["seg_type"],
                 seg.get("bearing"), seg.get("distance"),
                 seg.get("radius"), seg.get("sweep"),
                 seg.get("n_pts", 20)),
            )


def check_override_overlap(seg_index, span_end, db_path=None):
    """Check if a proposed override span overlaps any existing overrides.

    Returns list of overlapping seg_index values, or empty list if no overlap.
    Ignores the override at seg_index itself (for upsert replacement).
    """
    new_start = seg_index
    new_end = span_end if span_end is not None else seg_index
    with get_db(db_path or DB_PATH) as conn:
        rows = conn.execute(
            "SELECT DISTINCT seg_index, span_end FROM inner_wall_overrides "
            "WHERE seg_index != ?", (seg_index,)
        ).fetchall()
    conflicts = []
    for r in rows:
        ex_start = r["seg_index"]
        ex_end = r["span_end"] if r["span_end"] is not None else ex_start
        # Ranges overlap if new_start <= ex_end and ex_start <= new_end
        if new_start <= ex_end and ex_start <= new_end:
            conflicts.append(ex_start)
    return conflicts


def reset_inner_wall_overrides(db_path=None):
    """Re-seed inner wall overrides from defaults."""
    with get_db(db_path or DB_PATH) as conn:
        conn.execute("DELETE FROM inner_wall_overrides")
        _seed_inner_wall_overrides(conn)


def delete_inner_wall_override(seg_index, db_path=None):
    """Remove the override for a segment (reverts to default computation)."""
    with get_db(db_path or DB_PATH) as conn:
        conn.execute("DELETE FROM inner_wall_overrides WHERE seg_index = ?",
                     (seg_index,))


def snapshot_inner_wall_overrides(db_path=None):
    """Snapshot all overrides for undo state capture."""
    with get_db(db_path or DB_PATH) as conn:
        rows = conn.execute(
            f"SELECT {_IW_OV_COLS} FROM inner_wall_overrides "
            "ORDER BY seg_index, sub_seq"
        ).fetchall()
        return [dict(r) for r in rows]


def restore_inner_wall_overrides(snapshot, db_path=None):
    """Restore all overrides from a snapshot (for undo/redo)."""
    with get_db(db_path or DB_PATH) as conn:
        conn.execute("DELETE FROM inner_wall_overrides")
        for ov in snapshot:
            conn.execute(
                "INSERT INTO inner_wall_overrides "
                f"({_IW_OV_COLS}) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)",
                (ov["seg_index"], ov.get("span_end"), ov["sub_seq"],
                 ov["seg_type"], ov.get("bearing"), ov.get("distance"),
                 ov.get("radius"), ov.get("sweep"),
                 ov.get("n_pts", 20)),
            )
