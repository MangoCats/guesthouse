"""Variant registry, item dimension constants, and product URL lookup.

All element geometry is computed by the FormulaEvaluator from database-stored
JSON formulas.  This module provides:
  - VARIANTS dict and get_variant_flags() for variant metadata
  - Item dimension constants (used by evaluator formula seeding)
  - Product URL lookup (_resolve_product_url)
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

# ---------------------------------------------------------------------------
# Product URLs — replicated from floorplan/gen_floorplan.py per NF-4
# ---------------------------------------------------------------------------
# Base URLs shared across all variants
_PRODUCT_URLS_BASE = {
    "hamper": "https://www.homedepot.com/p/Casual-Home-Eco-Home-Laundry-Prep-Hamper-761-30/307595219",
    "kitchen_sink": "https://www.webstaurantstore.com/advance-tabco-fs1181824l-45-fabricated-one-compartment-sink-with-24-left-drainboard-18-x-18-x-14-bowl/109FS1L241818.html",
    "ice_maker": "https://www.homedepot.com/p/EUHOMY-17-3-in-100-lb-24H-Full-Ice-Sizes-Commercial-Ice-Maker-in-Black-33-lb-Storage-Bin-Ice-Full-Alert-and-Auto-Cleaning-CIM001-100BL-E/337185876",
    "microwave": "https://www.ikea.com/us/en/p/gatebo-microwave-oven-with-air-fryer-function-ikea-500-black-70603506/",
    "coffee_maker": "https://www.amazon.com/Holstein-Housewares-HH-0914701E-5-Cup-Coffee/dp/B08HSRCC4T/?th=1",
    "dining_table": "https://www.homedepot.com/pep/NEW-CLASSIC-HOME-FURNISHINGS-New-Classic-Furniture-Oscar-3-Piece-Wood-Top-Triangle-Dining-Set-Walnut-40-1651-D2C/327836175",
    "shelves": "https://www.ikea.com/us/en/p/kallax-shelving-unit-with-underframe-white-stained-oak-effect-black-s49442718/",
    "chair": "https://www.ikea.com/us/en/p/havberg-swivel-easy-chair-and-footstool-grann-bomstad-golden-brown-s59485321/",
    "ottoman": "https://www.ikea.com/us/en/p/havberg-swivel-easy-chair-and-footstool-grann-bomstad-golden-brown-s59485321/",
    "desk_chair": "https://www.amazon.com/BESTFAIR-Ergonomic-Office-Chair-Adjustable/dp/B0FDQDMP2D?th=1",
    "bath_sink": "https://www.magnushomeproducts.com/products/tripoli-vitreous-china-wall-mount-bathroom-sink",
    "util_sink": "https://www.magnushomeproducts.com/products/24-petten-matte-gray-vitreous-china-console-sink-with-black-powdercoat-steel-stand-and-shelves",
}
# Variant-specific URLs (keyed by item name → {variant_flag → url})
_PRODUCT_URLS_VARIANT = {
    "dryer": {
        "minik": "https://www.lowes.com/pd/Electrolux-8-cu-ft-Stackable-Steam-Cycle-Electric-Dryer-Titanium-ENERGY-STAR/5015416377",
    },
    "washer": {
        "minik": "https://www.lowes.com/pd/Electrolux-Smartboost-Optic-Whites-and-Pure-Rinse-4-5-cu-ft-High-Efficiency-Stackable-Steam-Cycle-Front-Load-Washer-Titanium-ENERGY-STAR/5015416375",
    },
    "fridge": {
        "minik": "https://www.ikea.com/us/en/p/bergsnaes-bottom-freezer-refrigerator-stainless-steel-color-60607883/",
        "default": "https://www.lowes.com/pd/LG-25-5-cu-ft-Bottom-Freezer-Refrigerator-with-Ice-Maker-Fingerprint-Resistant-Printproof-Stainless-Steel-ENERGY-STAR/1002543648",
    },
    "work_counter": {
        "default": "https://www.webstaurantstore.com/table-s-s-18x60-s-s-under/600TS1860S.html",
    },
    "kitchen_counter": {
        "minik": "https://www.webstaurantstore.com/regency-spec-line-30-x-72-14-gauge-stainless-steel-commercial-work-table-with-4-backsplash-and-undershelf/600TSSB3072S.html",
    },
    "north_counter": {
        "default": "https://www.webstaurantstore.com/regency-spec-line-30-x-36-14-gauge-stainless-steel-commercial-work-table-with-4-backsplash-and-undershelf/600TSSB3036S.html",
    },
    "cooktop": {
        "minik": "https://www.homedepot.com/p/Empava-Portable-13-4-in-Induction-Electric-Cooktop-in-Black-with-2-Elements-EMPV-ID12/313815692",
    },
    "toaster": {
        "minik": "https://www.amazon.com/Roter-Mond-Stainless-Independent-Removable/dp/B0CGTQZTDZ?th=1",
    },
    "sofa": {
        "minik": "https://www.ikea.com/us/en/p/saltsjoebaden-3-seat-sofa-gunnared-light-green-s89599953/",
    },
    "rocker": {
        "minik": "https://www.ikea.com/us/en/p/poaeng-rocking-chair-brown-gunnared-beige-s39502048/",
        "db": "https://www.ikea.com/us/en/p/poaeng-rocking-chair-brown-gunnared-beige-s39502048/",
    },
    "shelves2": {
        "db": "https://www.ikea.com/us/en/p/kallax-shelving-unit-with-underframe-white-stained-oak-effect-black-s49442718/",
    },
    "loveseat": {
        "default": "https://www.ikea.com/us/en/p/saltsjoebaden-loveseat-tonerud-red-brown-s59579188/",
    },
    "loveseat2": {
        "default": "https://www.ikea.com/us/en/p/saltsjoebaden-loveseat-tonerud-red-brown-s59579188/",
    },
    "et": {
        "default": "https://www.ikea.com/us/en/p/listerby-side-table-oak-veneer-30515314/",
    },
    "et_east": {
        "db": "https://www.ikea.com/us/en/p/listerby-side-table-oak-veneer-30515314/",
    },
    "et_west": {
        "db": "https://www.ikea.com/us/en/p/listerby-side-table-oak-veneer-30515314/",
    },
}


def _resolve_product_url(name, minik, db):
    """Return the product URL for an item, or None."""
    url = _PRODUCT_URLS_BASE.get(name)
    if url:
        return url
    variant_map = _PRODUCT_URLS_VARIANT.get(name)
    if not variant_map:
        return None
    if minik and "minik" in variant_map:
        return variant_map["minik"]
    if db and "db" in variant_map:
        return variant_map["db"]
    return variant_map.get("default")

