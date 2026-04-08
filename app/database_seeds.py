"""Static seed data for the ADU database.

All module-level literal data used exclusively by database.py seed functions
lives here so that database.py can focus on schema and CRUD operations.
"""

# ---------------------------------------------------------------------------
# Variant item dimension constants
# ---------------------------------------------------------------------------

# (name, value_inches, category, description)
_VARIANT_ITEM_CONSTANTS = [
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
    ("NORDVIKEN_L",      82.625, "furniture",  "NORDVIKEN extendable table core length"),
    ("NORDVIKEN_W",      41.375, "furniture",  "NORDVIKEN extendable table width"),
    ("SKOGSTA_L",        47.25,  "furniture",  "SKOGSTA bench with storage length"),
    ("SKOGSTA_W",        18.125, "furniture",  "SKOGSTA bench with storage width"),
    ("STOCKHOLM_SOFA_L", 83.125, "furniture",  "STOCKHOLM sofa length"),
    ("STOCKHOLM_SOFA_D", 34.625, "furniture",  "STOCKHOLM sofa depth"),
    ("TONSTAD_W",        18.5,   "furniture",  "TONSTAD chair width"),
    ("TONSTAD_D",        22.0,   "furniture",  "TONSTAD chair depth"),
]


# ---------------------------------------------------------------------------
# Survey seed data
# ---------------------------------------------------------------------------

# Raw traverse legs from shared/survey.py:_accumulate_legs()
# (seq, deg, min, sec, ft, inch, label)
_SURVEY_LEGS = [
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


# ---------------------------------------------------------------------------
# Variant definitions
# ---------------------------------------------------------------------------

# (name, label, flags_json, is_builtin)
_VARIANT_SEEDS = [
    ("standard", "Standard",       "{}",                 1),
    ("minik",    "Small Kitchen",   '{"minik": true}',   1),
    ("daybed",   "Daybed",          '{"db": true}',      1),
    ("bare",     "Room Dimensions", '{"bare": true}',    1),
    ("sf",       "Square Footage",  '{"sf": true}',      1),
    ("plumbing", "Plumbing",        '{"plumbing": true}', 1),
]


# ---------------------------------------------------------------------------
# Interior wall seed definitions
# ---------------------------------------------------------------------------

# (name, thickness_constant, orientation)
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
# Variant items master list
# (used by _seed_elements and _seed_catalog_items in database.py)
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
    ("dryer_sm", "appliance", {
        "label": "DRYER", "item_type": "appliance", "shape": "rect",
        "width": 32.0 / 12.0, "depth": 27.0 / 12.0,
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
    ("washer_sm", "appliance", {
        "label": "WASHER", "item_type": "appliance", "shape": "rect",
        "width": 32.0 / 12.0, "depth": 27.0 / 12.0,
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
    ("nordviken", "furniture", {
        "label": "TABLE",
        "item_type": "furniture",
        "shape": "rect",
        "width": 82.625 / 12.0,
        "depth": 41.375 / 12.0,
        "product_url": "https://www.ikea.com/us/en/p/nordviken-extendable-table-antique-stain-00488543/#content",
        # Clearance on each short end represents 15 9/16" extension leaf
        "clearance": [
            {"face": [1, 2], "distance": 15.5625 / 12.0},
            {"face": [3, 0], "distance": 15.5625 / 12.0},
        ],
        "variants": _VI_ALL,
    }, None),
    ("tonstad_chair", "furniture", {
        "label": "CHAIR",
        "item_type": "furniture",
        "shape": "rect",
        "width": 18.5 / 12.0,
        "depth": 22.0 / 12.0,
        "product_url": "https://www.ikea.com/us/en/p/tonstad-chair-bomstad-golden-brown-brown-oak-effect-s29602171/#content",
        "variants": _VI_ALL,
    }, None),
    ("skogsta_bench", "furniture", {
        "label": "BENCH",
        "item_type": "furniture",
        "shape": "rect",
        "width": 47.25 / 12.0,
        "depth": 18.125 / 12.0,
        "product_url": "https://www.ikea.com/us/en/p/skogsta-bench-with-storage-acacia-20611882/",
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
    ("stockholm_sofa", "furniture", {
        "label": "SOFA",
        "item_type": "furniture",
        "shape": "rect",
        "width": 83.125 / 12.0,
        "depth": 34.625 / 12.0,
        "product_url": "https://www.ikea.com/us/en/p/stockholm-sofa-seglora-natural-20245049/",
        "variants": _VI_ALL,
    }, None),
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
    # --- Accessibility ---
    ("turning_circle", "furniture", {
        "label": "TURNING CIRCLE", "item_type": "furniture", "shape": "circle",
        "radius": 2.5,
        "variants": _VI_ALL,
    }, None),
    ("three_feet", "furniture", {
        "label": "THREE FEET", "item_type": "furniture", "shape": "circle",
        "radius": 1.5,
        "variants": _VI_ALL,
    }, None),
]


# ---------------------------------------------------------------------------
# Door defaults
# ---------------------------------------------------------------------------

# (opening_name, width_inches, hinge_side, swing_direction, door_type)
# Width values match the corresponding *_DOOR_WIDTH constants in floorplan/constants.py.
_DOOR_SEED = [
    ("O3",  30.0, "north", "east",  "single"),
    ("O6",  42.0, "east",  "south", "single"),
    ("RO1", 36.0, "east",  "south", "single"),
    ("RO2", 36.0, "north", "east",  "single"),
    ("RO3", 36.0, "south", "west",  "single"),
    ("RO4", 36.0, "south", "west",  "single"),
    ("RO5", 36.0, "east",  "north", "single"),
    ("RO6", 24.0, "west",  "west",  "double"),
    ("RO7", 24.0, "east",  "east",  "double"),
]
