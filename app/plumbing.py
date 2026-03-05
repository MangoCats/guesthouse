"""Plumbing elements CRUD and seeding for the ADU Editor.

Manages supply pipes, drain pipes, fittings, and fixture connections
stored in the plumbing_elements table.
"""
import json

from app.database import get_db

# Valid plumbing element types
PLUMBING_TYPES = {"supply_pipe", "drain_pipe", "fitting", "fixture_connection"}

# Fixture definitions: (name, cold, hot, drain)
FIXTURE_DEFS = [
    ("Washer", True, True, True),
    ("Toilet1", True, False, True),
    ("Toilet2", True, False, True),
    ("Util Sink", True, True, True),
    ("Bath Sink", True, True, True),
    ("Fridge", True, False, False),
    ("Shower", True, True, True),
    ("Kitchen Sink", True, True, True),
    ("Dishwasher", False, True, True),
    ("Ice Maker", True, False, False),
]


def _row_to_dict(row):
    """Convert a sqlite3.Row to a plain dict with parsed JSON fields."""
    d = dict(row)
    for key in ("path", "properties"):
        if key in d and isinstance(d[key], str):
            d[key] = json.loads(d[key])
    return d


def get_plumbing_elements(db_path=None):
    """Return all plumbing elements as a list of dicts."""
    with get_db(db_path) as conn:
        rows = conn.execute(
            "SELECT * FROM plumbing_elements ORDER BY id"
        ).fetchall()
    return [_row_to_dict(r) for r in rows]


def get_plumbing_element(element_id, db_path=None):
    """Return a single plumbing element by ID, or None."""
    with get_db(db_path) as conn:
        row = conn.execute(
            "SELECT * FROM plumbing_elements WHERE id = ?", (element_id,)
        ).fetchone()
    return _row_to_dict(row) if row else None


def create_plumbing_element(type_, name, path=None, properties=None,
                            fixture=None, db_path=None):
    """Create a plumbing element. Returns the created record dict."""
    if type_ not in PLUMBING_TYPES:
        raise ValueError(f"Invalid plumbing type: {type_}")
    path_json = json.dumps(path or [])
    props_json = json.dumps(properties or {})
    with get_db(db_path) as conn:
        cur = conn.execute(
            "INSERT INTO plumbing_elements (type, name, path, properties, fixture) "
            "VALUES (?, ?, ?, ?, ?)",
            (type_, name, path_json, props_json, fixture),
        )
        return {
            "id": cur.lastrowid, "type": type_, "name": name,
            "path": path or [], "properties": properties or {},
            "fixture": fixture,
        }


def create_plumbing_raw(record, db_path=None):
    """Re-insert a plumbing element from a full record dict (used by undo)."""
    path_json = json.dumps(record.get("path", []))
    props_json = record.get("properties", "{}")
    if isinstance(props_json, dict):
        props_json = json.dumps(props_json)
    with get_db(db_path) as conn:
        conn.execute(
            "INSERT INTO plumbing_elements (id, type, name, path, properties, fixture) "
            "VALUES (?, ?, ?, ?, ?, ?)",
            (record["id"], record["type"], record["name"],
             path_json, props_json, record.get("fixture")),
        )


def update_plumbing_element(element_id, updates, db_path=None):
    """Update plumbing element fields. Returns updated record or None."""
    allowed = {"type", "name", "path", "properties", "fixture"}
    sets = []
    vals = []
    for k, v in updates.items():
        if k in allowed:
            if k in ("path", "properties") and isinstance(v, (dict, list)):
                v = json.dumps(v)
            sets.append(f"{k} = ?")
            vals.append(v)
    if not sets:
        return get_plumbing_element(element_id, db_path)
    vals.append(element_id)
    with get_db(db_path) as conn:
        conn.execute(
            f"UPDATE plumbing_elements SET {', '.join(sets)} WHERE id = ?", vals
        )
    return get_plumbing_element(element_id, db_path)


def delete_plumbing_element(element_id, db_path=None):
    """Delete a plumbing element by ID. Returns the deleted ID or None."""
    with get_db(db_path) as conn:
        row = conn.execute(
            "SELECT id FROM plumbing_elements WHERE id = ?", (element_id,)
        ).fetchone()
        if not row:
            return None
        conn.execute(
            "DELETE FROM plumbing_elements WHERE id = ?", (element_id,)
        )
    return element_id


def seed_plumbing(conn):
    """Seed default fixture connection records."""
    existing = conn.execute(
        "SELECT COUNT(*) as cnt FROM plumbing_elements"
    ).fetchone()["cnt"]
    if existing > 0:
        return
    for name, cold, hot, drain in FIXTURE_DEFS:
        props = json.dumps({"cold": cold, "hot": hot, "drain": drain})
        conn.execute(
            "INSERT INTO plumbing_elements (type, name, path, properties, fixture) "
            "VALUES (?, ?, '[]', ?, ?)",
            ("fixture_connection", name, props, name),
        )
