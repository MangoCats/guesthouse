"""Element business logic for the ADU Editor.

Provides higher-level operations on the elements table, including
variant-aware queries, cascading deletes, and the IW→constant mapping
that connects interior wall positions to their controlling constants.
"""
import json

from app.database import (
    get_db, get_all_elements, get_element, get_element_by_name,
    create_element, update_element, delete_element,
)

# ---------------------------------------------------------------------------
# IW → controlling constant mapping (Phase 3–11 transitional)
# ---------------------------------------------------------------------------

IW_CONSTANT_MAP = {
    "IW1":  "IW1_OFFSET_FROM_W9",
    "IW2":  "IW2_DIST_W2W5",
    "IW2O": None,                   # oblique — connects IW2 to IW2S, not directly movable
    "IW2S": "IW2S_W2REF_OFFSET",
    "IW3":  "IW3_DIST_W2W5",
    "IW4":  "IW4_GAP_IW11",
    "IW5":  "IW5_S_OFFSET_FROM_IW1",
    "IW6":  "IW6_OFFSET_FROM_W6",
    "IW7":  "IW7_OFFSET_FROM_W18W1",
    "IW8":  None,                   # derived from IW1 + IW2 positions
    "IW9":  "IW3_OFFSET_IW9",
    "IW11": "IW9_IW11_GAP",
    "IW12": "IW12_S_OFFSET_W18W1",
}

# ---------------------------------------------------------------------------
# IW → move axis mapping (Phase 4)
# ---------------------------------------------------------------------------
# Each moveable IW wall has one degree of freedom controlled by its constant.
# axis: "x" or "y" — the world-coordinate axis the constant affects.
# sign: +1 or -1 — delta_constant = drag_delta_on_axis * sign.
#   e.g. IW1 sign=-1: dragging south (dy<0) increases the constant.

IW_MOVE_AXIS = {
    "IW1":  ("y", -1),   # IW1_OFFSET_FROM_W9: +const = further south
    "IW2":  ("x", +1),   # IW2_DIST_W2W5: +const = further east
    "IW2O": None,         # oblique — not directly movable
    "IW2S": ("x", +1),   # IW2S_W2REF_OFFSET: +const = further east
    "IW3":  ("x", +1),   # IW3_DIST_W2W5: +const = further east
    "IW4":  ("x", +1),   # IW4_GAP_IW11: +const = further east
    "IW5":  ("y", -1),   # IW5_S_OFFSET_FROM_IW1: +const = further south
    "IW6":  ("y", -1),   # IW6_OFFSET_FROM_W6: +const = further south
    "IW7":  ("y", +1),   # IW7_OFFSET_FROM_W18W1: +const = further north
    "IW8":  None,         # derived — not directly movable
    "IW9":  ("x", +1),   # IW3_OFFSET_IW9: +const = further east
    "IW11": ("x", +1),   # IW9_IW11_GAP: +const = further east
    "IW12": ("y", +1),   # IW12_S_OFFSET_W18W1: +const = further north
}


def compute_constant_delta(iw_name, dx, dy):
    """Translate a world-coordinate move (dx, dy) into a constant change.

    Returns (constant_name, delta) or None if the wall is not movable.
    """
    axis_info = IW_MOVE_AXIS.get(iw_name)
    if axis_info is None:
        return None
    const_name = IW_CONSTANT_MAP.get(iw_name)
    if const_name is None:
        return None
    axis, sign = axis_info
    drag = dx if axis == "x" else dy
    delta = drag * sign
    if delta == 0:
        return None
    return (const_name, delta)


# Hosted rough openings per interior wall
IW_HOSTED_OPENINGS = {
    "IW1":  ["RO1"],
    "IW2":  [],
    "IW2O": ["RO4"],
    "IW2S": [],
    "IW3":  [],
    "IW4":  [],
    "IW5":  [],
    "IW6":  ["RO5"],
    "IW7":  [],
    "IW8":  [],
    "IW9":  ["RO3", "RO7"],
    "IW11": ["RO2", "RO6"],
    "IW12": [],
}


def get_elements_for_variant(variant=None, db_path=None):
    """Return elements visible to a variant (variant=NULL or matching)."""
    with get_db(db_path) as conn:
        rows = conn.execute(
            "SELECT id, type, name, properties, variant FROM elements "
            "WHERE variant IS NULL OR variant = ? "
            "ORDER BY type, name",
            (variant,),
        ).fetchall()
        return [dict(r) for r in rows]


def get_controlling_constant(iw_name):
    """Return the constant name that controls an IW's position, or None."""
    return IW_CONSTANT_MAP.get(iw_name)


def get_hosted_openings(iw_name):
    """Return list of rough opening names hosted by an interior wall."""
    return IW_HOSTED_OPENINGS.get(iw_name, [])
