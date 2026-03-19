"""Element business logic for the ADU Editor.

Provides higher-level operations on the elements table, including
variant-aware queries, cascading deletes, and the IW→constant mapping
that connects interior wall positions to their controlling constants.
"""
from app.database import get_db

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


# Reverse of IW_CONSTANT_MAP: constant name → list of IW names
CONSTANT_TO_IW = {}
for _iw, _cname in IW_CONSTANT_MAP.items():
    if _cname:
        CONSTANT_TO_IW.setdefault(_cname, []).append(_iw)
del _iw, _cname

# Span-dimension constants: IW wall name → (property_label, constant_name).
# Only for walls where the span is controlled by a single constant.
# "width"  = the bbox E-W extent (horizontal span)
# "height" = the bbox N-S extent (vertical span)
# Note: some span constants are shared with another wall's "position" constant
# (e.g. IW3_OFFSET_IW9 positions IW9 and sets IW7's width; IW4_GAP_IW11
#  positions IW4 and sets IW12's width). That coupling is by design.
_IW_SPAN_CONST = {
    "IW2":  ("height", "IW2_LENGTH"),      # N-S span (V wall)
    "IW2S": ("height", "IW2S_LENGTH"),     # N-S span (V wall)
    "IW7":  ("width",  "IW3_OFFSET_IW9"), # E-W span = gap IW3.SE→IW9.SW
    "IW12": ("width",  "IW4_GAP_IW11"),   # E-W span = gap IW11.SE→IW4.SW
}

# Combined property map: IW wall name → {property_label: constant_name}.
# Stored in elements.properties["prop_constants"] so the UI is DB-driven.
# Seeding source only — the live value lives in the DB.
IW_PROP_CONSTANTS = {}
for _iw, _cname in IW_CONSTANT_MAP.items():
    if _cname:
        IW_PROP_CONSTANTS.setdefault(_iw, {})["position"] = _cname
for _iw, (_prop, _cname) in _IW_SPAN_CONST.items():
    IW_PROP_CONSTANTS.setdefault(_iw, {})[_prop] = _cname
del _iw, _cname, _prop


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
