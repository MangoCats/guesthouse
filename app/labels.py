"""Room label seeding and label/dimension name generators.

Manages the migration of room labels from the legacy `room_label_offsets`
table into the `elements` table (type='label', source='room').
"""
import json
import re


# The 11 room names, matching _compute_room_labels in engine.py
ROOM_LABEL_NAMES = [
    "BEDROOM", "UTIL_N", "UTIL_S", "KITCHEN", "LIVING",
    "BATH", "OFFICE", "E CLOSET", "W CLOSET", "STORAGE", "WH",
]


def seed_room_labels(conn):
    """Create room label elements if they don't already exist.

    Migrates any offsets from the legacy room_label_offsets table.
    Called from init_db() and reset_elements().
    """
    # Check if room label elements already exist
    count = conn.execute(
        "SELECT COUNT(*) FROM elements WHERE type = 'label'"
    ).fetchone()[0]
    if count >= len(ROOM_LABEL_NAMES):
        return

    # Load legacy offsets if available
    offsets = {}
    try:
        rows = conn.execute(
            "SELECT room_name, offset_e, offset_n FROM room_label_offsets"
        ).fetchall()
        offsets = {r["room_name"]: (r["offset_e"], r["offset_n"]) for r in rows}
    except Exception:
        pass  # Table may not exist in test DBs

    for name in ROOM_LABEL_NAMES:
        # Skip if this specific label already exists
        existing = conn.execute(
            "SELECT id FROM elements WHERE type = 'label' AND name = ?",
            (name,),
        ).fetchone()
        if existing:
            continue
        de, dn = offsets.get(name, (0.0, 0.0))
        props = json.dumps({
            "source": "room",
            "room_name": name,
            "text": name,
            "offset_e": de,
            "offset_n": dn,
            "rotation": 0,
            "font_size": 0.25,
        })
        conn.execute(
            "INSERT INTO elements (type, name, properties, variant) "
            "VALUES (?, ?, ?, ?)",
            ("label", name, props, None),
        )


def next_dimension_name(elements):
    """Return the next auto-name for a user dimension (UD1, UD2, ...)."""
    max_n = 0
    for e in elements:
        m = re.match(r"^UD(\d+)$", e.get("name", ""))
        if m:
            n = int(m.group(1))
            if n > max_n:
                max_n = n
    return f"UD{max_n + 1}"


def next_label_name(elements):
    """Return the next auto-name for a user label (UL1, UL2, ...)."""
    max_n = 0
    for e in elements:
        m = re.match(r"^UL(\d+)$", e.get("name", ""))
        if m:
            n = int(m.group(1))
            if n > max_n:
                max_n = n
    return f"UL{max_n + 1}"
