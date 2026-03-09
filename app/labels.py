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


# ---------------------------------------------------------------------------
# Built-in dimensions: anchor specs for all 22 dimension lines
# ---------------------------------------------------------------------------
# Each dimension has start_anchor and end_anchor that resolve to [E, N]
# coordinates via _resolve_anchor() in engine.py.  The new "line_intersection"
# anchor type computes the intersection of two infinite lines, each defined
# by a point spec + direction spec.

# Shared sub-expressions used by multiple dimensions
_F9F11_MID_OFFSET = {"offset": {"midpoint": ["F9", "F11"]}, "dir": "east", "dist": 1.0}
_W18W1_DIR = {"segment": ["W18", "W1"]}
_W18W1_PERP = {"segment_perp": ["W18", "W1"]}
_STOR_REF = {"midpoint": [{"face_mid": "IW5", "face": "north"},
                           {"face_mid": "IW1", "face": "south"}]}
_DIM12_REF = {"offset": "F6", "dir": "east", "dist": 1.0}
_DIM14_ARC_REF = {"offset": "F7", "dir": "north", "dist": -4.0 / 12.0}

BUILTIN_DIMENSIONS = [
    # dim01: IW1-north → W9 (vertical probe east of F9-F11 midpoint)
    {
        "name": "dim01", "variant": None,
        "start_anchor": {
            "type": "line_intersection",
            "line1_point": {"face_mid": "IW1", "face": "north"},
            "line1_dir": {"face_along": "IW1", "face": "north"},
            "line2_point": _F9F11_MID_OFFSET, "line2_dir": "north",
        },
        "end_anchor": {
            "type": "line_intersection",
            "line1_point": "W9", "line1_dir": "east",
            "line2_point": _F9F11_MID_OFFSET, "line2_dir": "north",
        },
    },
    # dim02: IW2-east → W13-W12 (horizontal probe at F12-F13 midpoint)
    {
        "name": "dim02", "variant": None,
        "start_anchor": {
            "type": "line_intersection",
            "line1_point": {"face_mid": "IW2", "face": "east"},
            "line1_dir": {"face_along": "IW2", "face": "east"},
            "line2_point": {"midpoint": ["F12", "F13"]}, "line2_dir": "east",
        },
        "end_anchor": {
            "type": "line_intersection",
            "line1_point": "W13", "line1_dir": {"segment": ["W13", "W12"]},
            "line2_point": {"midpoint": ["F12", "F13"]}, "line2_dir": "east",
        },
    },
    # dim03: IW12-south-mid → W18-W1 (rotated, along IW11 east direction)
    {
        "name": "dim03", "variant": None, "label_prefix": "CLOSET ",
        "start_anchor": {"type": "wall_face", "target": "IW12", "face": "south"},
        "end_anchor": {
            "type": "line_intersection",
            "line1_point": {"face_mid": "IW12", "face": "south"},
            "line1_dir": {"face_along": "IW11", "face": "east"},
            "line2_point": "W18", "line2_dir": _W18W1_DIR,
        },
    },
    # dim04: IW7-south-mid → W18-W1 (rotated, along IW7 west direction)
    {
        "name": "dim04", "variant": None, "label_prefix": "CLOSET ",
        "start_anchor": {"type": "wall_face", "target": "IW7", "face": "south"},
        "end_anchor": {
            "type": "line_intersection",
            "line1_point": {"face_mid": "IW7", "face": "south"},
            "line1_dir": {"face_along": "IW7", "face": "west"},
            "line2_point": "W18", "line2_dir": _W18W1_DIR,
        },
    },
    # dim05: W2 → IW3-west (parallel to W18-W1 at IW3 west midpoint)
    {
        "name": "dim05", "variant": None,
        "start_anchor": {
            "type": "line_intersection",
            "line1_point": "W2", "line1_dir": "north",
            "line2_point": {"face_mid": "IW3", "face": "west"},
            "line2_dir": _W18W1_DIR,
        },
        "end_anchor": {"type": "wall_face", "target": "IW3", "face": "west"},
    },
    # dim06: IW4-east → O8-north (perp to W14-W15 at O8 center)
    {
        "name": "dim06", "variant": None,
        "start_anchor": {
            "type": "line_intersection",
            "line1_point": {"face_mid": "IW4", "face": "east"},
            "line1_dir": {"face_along": "IW4", "face": "east"},
            "line2_point": {"opening_face_mid": "O8", "face": "north"},
            "line2_dir": {"segment_perp": ["W14", "W15"]},
        },
        "end_anchor": {"type": "opening_face", "target": "O8", "face": "north"},
    },
    # dim07: IW11-east → W14-W15 (horizontal at storage midpoint)
    {
        "name": "dim07", "variant": None, "label_prefix": "STORAGE ",
        "start_anchor": {
            "type": "line_intersection",
            "line1_point": {"face_mid": "IW11", "face": "east"},
            "line1_dir": {"face_along": "IW11", "face": "east"},
            "line2_point": _STOR_REF, "line2_dir": "east",
        },
        "end_anchor": {
            "type": "line_intersection",
            "line1_point": "W14", "line1_dir": {"segment": ["W14", "W15"]},
            "line2_point": _STOR_REF, "line2_dir": "east",
        },
    },
    # dim08: O1-north-face → IW9-west (horizontal at O1 centroid northing)
    {
        "name": "dim08", "variant": None,
        "start_anchor": {
            "type": "line_intersection",
            "line1_point": {"opening_face_mid": "O1", "face": "north"},
            "line1_dir": "north",
            "line2_point": {"opening_centroid": "O1"}, "line2_dir": "east",
        },
        "end_anchor": {
            "type": "line_intersection",
            "line1_point": {"face_mid": "IW9", "face": "west"},
            "line1_dir": {"face_along": "IW9", "face": "west"},
            "line2_point": {"opening_centroid": "O1"}, "line2_dir": "east",
        },
    },
    # dim09: W2 → IW2-west (horizontal at O2 centroid)
    {
        "name": "dim09", "variant": None,
        "start_anchor": {
            "type": "line_intersection",
            "line1_point": "W2", "line1_dir": "north",
            "line2_point": {"opening_centroid": "O2"}, "line2_dir": "east",
        },
        "end_anchor": {
            "type": "line_intersection",
            "line1_point": {"face_mid": "IW2", "face": "west"},
            "line1_dir": {"face_along": "IW2", "face": "west"},
            "line2_point": {"opening_centroid": "O2"}, "line2_dir": "east",
        },
    },
    # dim10: W2 → IW2S-west (horizontal at F5 northing)
    {
        "name": "dim10", "variant": None,
        "start_anchor": {
            "type": "line_intersection",
            "line1_point": "W2", "line1_dir": "north",
            "line2_point": "F5", "line2_dir": "east",
        },
        "end_anchor": {
            "type": "line_intersection",
            "line1_point": {"face_mid": "IW2S", "face": "west"},
            "line1_dir": {"face_along": "IW2S", "face": "west"},
            "line2_point": "F5", "line2_dir": "east",
        },
    },
    # dim11: IW5-south → W18 (vertical at F18 easting)
    {
        "name": "dim11", "variant": None,
        "start_anchor": {
            "type": "line_intersection",
            "line1_point": {"face_mid": "IW5", "face": "south"},
            "line1_dir": {"face_along": "IW5", "face": "south"},
            "line2_point": "F18", "line2_dir": "north",
        },
        "end_anchor": {
            "type": "line_intersection",
            "line1_point": "W18", "line1_dir": "east",
            "line2_point": "F18", "line2_dir": "north",
        },
    },
    # dim12a: IW6-north → W6 (vertical, normal variants only)
    {
        "name": "dim12a", "variant": None,
        "start_anchor": {
            "type": "line_intersection",
            "line1_point": {"face_mid": "IW6", "face": "north"},
            "line1_dir": {"face_along": "IW6", "face": "north"},
            "line2_point": _DIM12_REF, "line2_dir": "north",
        },
        "end_anchor": {
            "type": "line_intersection",
            "line1_point": "W6", "line1_dir": "east",
            "line2_point": _DIM12_REF, "line2_dir": "north",
        },
    },
    # dim12b: IW8-north → IW6-south (vertical, normal variants only)
    {
        "name": "dim12b", "variant": None,
        "start_anchor": {
            "type": "line_intersection",
            "line1_point": {"face_mid": "IW8", "face": "north"},
            "line1_dir": {"face_along": "IW8", "face": "north"},
            "line2_point": _DIM12_REF, "line2_dir": "north",
        },
        "end_anchor": {
            "type": "line_intersection",
            "line1_point": {"face_mid": "IW6", "face": "south"},
            "line1_dir": {"face_along": "IW6", "face": "south"},
            "line2_point": _DIM12_REF, "line2_dir": "north",
        },
    },
    # dim12bare: IW8-north → W6 (vertical, bare variant only)
    {
        "name": "dim12bare", "variant": "bare",
        "start_anchor": {
            "type": "line_intersection",
            "line1_point": {"face_mid": "IW8", "face": "north"},
            "line1_dir": {"face_along": "IW8", "face": "north"},
            "line2_point": {"opening_centroid": "O4"}, "line2_dir": "north",
        },
        "end_anchor": {
            "type": "line_intersection",
            "line1_point": "W6", "line1_dir": "east",
            "line2_point": {"opening_centroid": "O4"}, "line2_dir": "north",
        },
    },
    # dim13: F18 → F6 external (vertical west of F2)
    {
        "name": "dim13", "variant": None,
        "start_anchor": {
            "type": "line_intersection",
            "line1_point": "F18", "line1_dir": "east",
            "line2_point": {"offset": "F2", "dir": "east", "dist": -2.7},
            "line2_dir": "north",
        },
        "end_anchor": {
            "type": "line_intersection",
            "line1_point": "F6", "line1_dir": "east",
            "line2_point": {"offset": "F2", "dir": "east", "dist": -2.7},
            "line2_dir": "north",
        },
    },
    # dim14: Arc width (horizontal, F7-side arc → F11a-side arc)
    {
        "name": "dim14", "variant": None,
        "start_anchor": {
            "type": "line_intersection",
            "line1_point": {"arc_point": {
                "center": "C7", "radius_key": "R_a7",
                "reference": _DIM14_ARC_REF, "side": "east",
            }},
            "line1_dir": "north",
            "line2_point": {"offset": "F6", "dir": "north", "dist": 1.0},
            "line2_dir": "east",
        },
        "end_anchor": {
            "type": "line_intersection",
            "line1_point": {"arc_point": {
                "center": "C11a", "radius_key": "R_a11",
                "reference": _DIM14_ARC_REF, "side": "west",
            }},
            "line1_dir": "north",
            "line2_point": {"offset": "F6", "dir": "north", "dist": 1.0},
            "line2_dir": "east",
        },
    },
    # dim15: F2 → F15 external (horizontal south of F18)
    {
        "name": "dim15", "variant": None,
        "start_anchor": {
            "type": "line_intersection",
            "line1_point": "F2", "line1_dir": "north",
            "line2_point": {"offset": "F18", "dir": "north", "dist": -3.0},
            "line2_dir": "east",
        },
        "end_anchor": {
            "type": "line_intersection",
            "line1_point": "F15", "line1_dir": "north",
            "line2_point": {"offset": "F18", "dir": "north", "dist": -3.0},
            "line2_dir": "east",
        },
    },
    # dim17: O10-north-face → IW1-south (rotated, perp to W18-W1)
    {
        "name": "dim17", "variant": None,
        "start_anchor": {"type": "opening_face", "target": "O10", "face": "north"},
        "end_anchor": {
            "type": "line_intersection",
            "line1_point": {"opening_face_mid": "O10", "face": "north"},
            "line1_dir": _W18W1_PERP,
            "line2_point": {"face_mid": "IW1", "face": "south"},
            "line2_dir": {"face_along": "IW1", "face": "south"},
        },
    },
    # dim18: IW9-east-mid → IW11-west (rotated, perp to IW9 east face)
    {
        "name": "dim18", "variant": None,
        "start_anchor": {"type": "wall_face", "target": "IW9", "face": "east"},
        "end_anchor": {
            "type": "line_intersection",
            "line1_point": {"face_mid": "IW9", "face": "east"},
            "line1_dir": {"face_perp": "IW9", "face": "east"},
            "line2_point": {"face_mid": "IW11", "face": "west"},
            "line2_dir": {"face_along": "IW11", "face": "west"},
        },
    },
    # dim19: O11-center → IW8-south (vertical at O11 centroid easting)
    {
        "name": "dim19", "variant": None,
        "start_anchor": {
            "type": "line_intersection",
            "line1_point": {"opening_centroid": "O11"}, "line1_dir": "north",
            "line2_point": {"opening_face_mid": "O11", "face": "north"},
            "line2_dir": "east",
        },
        "end_anchor": {
            "type": "line_intersection",
            "line1_point": {"face_mid": "IW8", "face": "south"},
            "line1_dir": {"face_along": "IW8", "face": "south"},
            "line2_point": {"opening_centroid": "O11"}, "line2_dir": "north",
        },
    },
    # dim22: IW12-north-mid → IW5-south (rotated, perp to W18-W1)
    {
        "name": "dim22", "variant": None,
        "start_anchor": {"type": "wall_face", "target": "IW12", "face": "north"},
        "end_anchor": {
            "type": "line_intersection",
            "line1_point": {"face_mid": "IW12", "face": "north"},
            "line1_dir": _W18W1_PERP,
            "line2_point": {"face_mid": "IW5", "face": "south"},
            "line2_dir": {"face_along": "IW5", "face": "south"},
        },
    },
    # dim20: IW2-east → W14 (horizontal at W14 northing, bare only)
    {
        "name": "dim20", "variant": "bare",
        "start_anchor": {
            "type": "line_intersection",
            "line1_point": {"face_mid": "IW2", "face": "east"},
            "line1_dir": {"face_along": "IW2", "face": "east"},
            "line2_point": "W14", "line2_dir": "east",
        },
        "end_anchor": {"type": "point", "target": "W14"},
    },
    # dim21: IW1-north → midpoint(W11a, W11b) (vertical, bare only)
    {
        "name": "dim21", "variant": "bare",
        "start_anchor": {
            "type": "line_intersection",
            "line1_point": {"face_mid": "IW1", "face": "north"},
            "line1_dir": {"face_along": "IW1", "face": "north"},
            "line2_point": {"midpoint": ["W11a", "W11b"]}, "line2_dir": "north",
        },
        "end_anchor": {
            "type": "computed",
            "spec": {"midpoint": ["W11a", "W11b"]},
        },
    },
]


def seed_builtin_dimensions(conn):
    """Create built-in dimension elements if they don't already exist.

    Called from init_db() and reset_elements().
    """
    count = conn.execute(
        "SELECT COUNT(*) FROM elements WHERE type = 'dimension' "
        "AND json_extract(properties, '$.source') = 'builtin'"
    ).fetchone()[0]
    if count >= len(BUILTIN_DIMENSIONS):
        return

    for dim in BUILTIN_DIMENSIONS:
        existing = conn.execute(
            "SELECT id FROM elements WHERE type = 'dimension' AND name = ?",
            (dim["name"],),
        ).fetchone()
        if existing:
            continue
        props_dict = {
            "source": "builtin",
            "start_anchor": dim["start_anchor"],
            "end_anchor": dim["end_anchor"],
            "offset": 0,
            "label_rotation": "parallel",
        }
        if dim.get("label_prefix"):
            props_dict["label_prefix"] = dim["label_prefix"]
        props = json.dumps(props_dict)
        conn.execute(
            "INSERT INTO elements (type, name, properties, variant) "
            "VALUES (?, ?, ?, ?)",
            ("dimension", dim["name"], props, dim["variant"]),
        )


def _next_auto_name(elements, prefix, pattern):
    """Return next auto-name for a given prefix (e.g. 'UD', 'UL')."""
    max_n = 0
    for e in elements:
        m = re.match(pattern, e.get("name", ""))
        if m:
            n = int(m.group(1))
            if n > max_n:
                max_n = n
    return f"{prefix}{max_n + 1}"


def next_dimension_name(elements):
    """Return the next auto-name for a user dimension (UD1, UD2, ...)."""
    return _next_auto_name(elements, "UD", r"^UD(\d+)$")


def next_label_name(elements):
    """Return the next auto-name for a user label (UL1, UL2, ...)."""
    return _next_auto_name(elements, "UL", r"^UL(\d+)$")
