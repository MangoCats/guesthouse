"""Element business logic for the ADU Editor.

Provides higher-level operations on the elements table, including
variant-aware queries, cascading deletes, and hosted-opening lookups.

Position/move metadata (formerly IW_CONSTANT_MAP / IW_MOVE_AXIS) now lives
inside each wall formula as position_constant / position_axis / position_sign
fields, read at runtime by app/server.py via get_iw_formulas().
"""
from app.database import get_db

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


def get_hosted_openings(iw_name):
    """Return list of rough opening names hosted by an interior wall."""
    return IW_HOSTED_OPENINGS.get(iw_name, [])
