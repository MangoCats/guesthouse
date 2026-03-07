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
    fixture    TEXT
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

CREATE TABLE IF NOT EXISTS element_formulas (
    id           INTEGER PRIMARY KEY AUTOINCREMENT,
    element_name TEXT NOT NULL,          -- FK to elements.name
    param_name   TEXT NOT NULL,          -- 'position', 'poly', 'center', etc.
    formula_json TEXT NOT NULL,          -- JSON formula spec
    locked       INTEGER DEFAULT 0,     -- 1 = frozen at locked_value
    locked_value TEXT,                   -- JSON: frozen computed result
    variant      TEXT,                   -- NULL = all variants
    UNIQUE(element_name, param_name, variant)
);

CREATE TABLE IF NOT EXISTS formula_deps (
    element_name TEXT NOT NULL,          -- element with the formula
    param_name   TEXT NOT NULL,          -- which parameter formula
    dep_type     TEXT NOT NULL,          -- 'element', 'point', 'constant'
    dep_name     TEXT NOT NULL,          -- name of dependency target
    UNIQUE(element_name, param_name, dep_type, dep_name)
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
            _seed_variant_item_constants(conn)
            _seed_outline_chain(conn)
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
        ("plumbing_edit", "Plumbing Edit", "", "plumbing_edit", "design"),
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


def _seed_elements(conn):
    """Seed the elements table with the 13 interior walls."""

    for name, thickness_const, orientation in _IW_SEED:
        props = json.dumps({
            "thickness_constant": thickness_const,
            "orientation": orientation,
        })
        conn.execute(
            "INSERT OR REPLACE INTO elements (type, name, properties, variant) "
            "VALUES (?, ?, ?, ?)",
            ("wall", name, props, None),
        )


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
               "center_name", "n_pts", "end_name"}
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
            "center_name, n_pts, end_name) "
            "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)",
            (seq, row_data["seg_type"], row_data.get("distance"),
             row_data.get("radius"), row_data.get("sweep_name"),
             row_data.get("sweep"), row_data.get("center_name"),
             row_data.get("n_pts", 60), row_data["end_name"]),
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
                "center_name, n_pts, end_name) "
                "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)",
                (row["seq"], row["seg_type"], row.get("distance"),
                 row.get("radius"), row.get("sweep_name"),
                 row.get("sweep"), row.get("center_name"),
                 row.get("n_pts", 60), row["end_name"]),
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


def add_variant_exclusion(variant, element_type, element_name, db_path=None):
    """Add a single variant exclusion."""
    with get_db(db_path) as conn:
        conn.execute(
            "INSERT OR IGNORE INTO variant_exclusions "
            "(variant, element_type, element_name) VALUES (?, ?, ?)",
            (variant, element_type, element_name),
        )


def remove_variant_exclusion(variant, element_type, element_name, db_path=None):
    """Remove a single variant exclusion."""
    with get_db(db_path) as conn:
        conn.execute(
            "DELETE FROM variant_exclusions "
            "WHERE variant = ? AND element_type = ? AND element_name = ?",
            (variant, element_type, element_name),
        )


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
