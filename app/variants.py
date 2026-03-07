"""Variant registry and item dimension constants.

All element geometry is computed by the FormulaEvaluator from database-stored
JSON formulas.  This module provides:
  - VARIANTS dict and get_variant_flags() for variant metadata
  - Item dimension constants (used by evaluator formula seeding)
"""

# ---------------------------------------------------------------------------
# Variant registry
# ---------------------------------------------------------------------------
VARIANTS = {
    "standard": {"label": "Standard", "flags": {}},
    "minik":    {"label": "Small Kitchen", "flags": {"minik": True}},
    "daybed":   {"label": "Daybed", "flags": {"db": True}},
    "bare":     {"label": "Room Dimensions", "flags": {"bare": True}},
    "sf":       {"label": "Square Footage", "flags": {"sf": True}},
}


def get_variant_flags(variant, db_path=None):
    """Return flags dict for a variant, from DB with fallback to VARIANTS dict."""
    try:
        from app.database import get_variant
        v = get_variant(variant, db_path)
        if v and v.get("flags"):
            return v["flags"]
    except Exception:
        pass
    info = VARIANTS.get(variant, VARIANTS["standard"])
    return info["flags"]


# ---------------------------------------------------------------------------
# Item dimensions hardcoded in gen_floorplan.py (not in floorplan/constants.py)
# ---------------------------------------------------------------------------
HAMPER_W = 31.5 / 12.0
HAMPER_D = 19.0 / 12.0
MINIK_APPL_W = 32.0 / 12.0
MINIK_APPL_D = 27.0 / 12.0
MICROWAVE_W = 19.5 / 12.0
MICROWAVE_D = 16.625 / 12.0
COFFEE_W = 7.2 / 12.0
COFFEE_D = 9.2 / 12.0
COOKTOP_W = 13.4 / 12.0
COOKTOP_D = 20.5 / 12.0
TOASTER_W = 13.7 / 12.0
TOASTER_D = 12.5 / 12.0
DINING_CHAIR_W = 18.0 / 12.0
DINING_CHAIR_D = 21.0 / 12.0
DINING_TBL_BASE = 31.5 / 12.0
DINING_TBL_H = 35.25 / 12.0
DAYBED_W = 86.0 / 12.0
DAYBED_D = 43.0 / 12.0
WORK_CTR_W = 60.0 / 12.0
WORK_CTR_D = 18.0 / 12.0
STD_FRIDGE_W = 32.75 / 12.0
STD_FRIDGE_D = 35.0 / 12.0
SOFA_FULL_W = 80.75 / 12.0
SOFA_FULL_D = 34.625 / 12.0

