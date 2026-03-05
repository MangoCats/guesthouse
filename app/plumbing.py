"""Plumbing elements CRUD and seeding for the ADU Editor.

Manages supply pipes, drain pipes, fittings, and fixture connections
stored in the plumbing_elements table.
"""
import json
import math

from app.database import get_db
from shared.geometry import seg_vecs, offset_pt, line_isect

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


# ---------------------------------------------------------------------------
# Reference plumbing computation — replicates _render_plumbing_path() from
# floorplan/gen_floorplan.py using geometry dict data.
# Constants hardcoded per NF-4 duplication principle.
# ---------------------------------------------------------------------------
_WH_RADIUS = 14.0 / 12.0
_KITCHEN_SINK_WIDTH = 45.0 / 12.0
_DW_WIDTH = 28.0 / 12.0
_STOVE_WIDTH = 30.0 / 12.0
_NORTH_CTR_LENGTH = 36.0 / 12.0
_KITCHEN_APPL_GAP = 3.0 / 12.0
_SVG_TO_FT = 10.0 / 30.48
_APPL_D = 27.0 / 12.0


def _poly_centroid(poly):
    """Centroid of a polygon given as list of [E, N]."""
    n = len(poly)
    if n == 0:
        return [0, 0]
    return [sum(p[0] for p in poly) / n, sum(p[1] for p in poly) / n]


def compute_reference_plumbing(geom, wall_t):
    """Compute reference plumbing pipe paths and fixture positions from geometry.

    Replicates the pipe routing logic from _render_plumbing_path() in
    floorplan/gen_floorplan.py, translating layout/pts references to the
    geometry dict returned by compute_geometry().

    Returns dict with keys:
      "pipes": list of {name, type, path, properties}
      "fixture_positions": dict of fixture_name → [E, N]
    """
    pts = geom["points"]
    iw = geom["interior_walls"]
    appl = geom["appliances"]
    vi = geom.get("variant_items", {})

    iw8 = iw["IW8"]["poly"]
    iw2 = iw["IW2"]["poly"]
    iw2o = iw["IW2O"]["poly"]
    iw2s = iw["IW2S"]["poly"]

    # --- Washer center (replicates _render_plumbing_path logic) ---
    w2w5_al, w2w5_in = seg_vecs(pts["W2"], pts["W5"])
    dryer_sw = appl["dryer"]["poly"][0]
    dryer_nw = offset_pt(dryer_sw, _APPL_D, w2w5_al)
    washer_sw = offset_pt(dryer_nw, 1.0 / 12.0, w2w5_al)
    washer_nw = offset_pt(washer_sw, _APPL_D, w2w5_al)
    washer_cn = (washer_sw[1] + washer_nw[1]) / 2

    # Reference points
    f2_e = pts["F2"][0]
    f1_n = pts["F1"][1]
    f14_e = pts["F14"][0]

    # ===== Green utility path (drain) =====
    wp1 = [f2_e, washer_cn]
    wp2 = [f2_e - 24.0 / 12.0, washer_cn]
    wp3 = [f2_e - 24.0 / 12.0, f1_n - 24.0 / 12.0]
    wp4 = [f14_e, f1_n - 24.0 / 12.0]

    # ===== Blue cold supply line =====
    gw = 1.5 / 12.0   # 1.5" offset from exterior wall
    gi = 1.0 / 12.0   # 1.0" offset from interior walls

    # Outside waypoints
    bp1 = [f14_e + 5.0, f1_n - 24.0 / 12.0]
    bp2 = [f2_e - 24.0 / 12.0, f1_n - 24.0 / 12.0]
    bp3 = [f2_e - 24.0 / 12.0, washer_cn]

    # Inside building
    bp4 = [pts["W2"][0] + gw, washer_cn]
    bp5 = [pts["W2"][0] + gw, iw8[3][1] - gi]
    bp6 = [iw2[0][0] + gi, iw8[3][1] - gi]

    # IW2 west face: poly[0]→poly[3]
    iw2_w_al, iw2_w_in = seg_vecs(iw2[0], iw2[3])
    iw2_anchor = offset_pt(iw2[0], gi, iw2_w_in)

    # IW2o oblique: poly[1]→poly[2]
    iw2o_al, iw2o_in = seg_vecs(iw2o[1], iw2o[2])
    iw2o_anchor = offset_pt(iw2o[1], gi, iw2o_in)
    bp7 = list(line_isect(iw2_anchor, iw2_w_al, iw2o_anchor, iw2o_al))

    # IW2s west face: poly[0]→poly[3]
    iw2s_w_al, iw2s_w_in = seg_vecs(iw2s[0], iw2s[3])
    iw2s_anchor = offset_pt(iw2s[0], gi, iw2s_w_in)
    bp8 = list(line_isect(iw2o_anchor, iw2o_al, iw2s_anchor, iw2s_w_al))

    # W9-W10 north wall
    w9w10_al, w9w10_in = seg_vecs(pts["W9"], pts["W10"])
    w9_inset = offset_pt(pts["W9"], gi, w9w10_in)
    bp9 = list(line_isect(iw2s_anchor, iw2s_w_al, w9_inset, w9w10_al))

    # Kitchen appliance chain distances
    iw2s_ne = iw2s[2]
    iw2_d = ((iw2s_ne[0] - pts["W9"][0]) * w9w10_al[0] +
             (iw2s_ne[1] - pts["W9"][1]) * w9w10_al[1])
    dw_d = (iw2_d + _NORTH_CTR_LENGTH + _KITCHEN_APPL_GAP +
            _STOVE_WIDTH + _KITCHEN_APPL_GAP + 2.0 / 12.0 +
            _KITCHEN_SINK_WIDTH + _KITCHEN_APPL_GAP)
    ice_d = dw_d + _DW_WIDTH + 2.0 / 12.0
    bp10 = list(offset_pt(w9_inset, ice_d, w9w10_al))

    blue_path = [bp1, bp2, bp3, bp4, bp5, bp6, bp7, bp8, bp9, bp10]

    # ===== Blue WH T-stub =====
    iw2s_e_al, iw2s_e_out = seg_vecs(iw2s[1], iw2s[2])
    wh_ref = offset_pt(iw2s[2], _WH_RADIUS, iw2s_e_out)
    wh_tan_r = (geom["radii"]["R_a7"] - wall_t) - _WH_RADIUS
    c7 = pts["C7"]
    wh_d = (wh_ref[0] - c7[0], wh_ref[1] - c7[1])
    wh_d_al = wh_d[0] * iw2s_e_al[0] + wh_d[1] * iw2s_e_al[1]
    wh_d2 = wh_d[0] ** 2 + wh_d[1] ** 2
    wh_t = -wh_d_al + math.sqrt(wh_tan_r ** 2 - wh_d2 + wh_d_al ** 2)
    wh_ctr = offset_pt(wh_ref, wh_t, iw2s_e_al)

    wh_bl = list(line_isect(w9_inset, w9w10_al, (wh_ctr[0], 0), (0, 1)))
    wh_bl_n = [wh_bl[0], wh_bl[1] + 1.0]

    # ===== Red hot supply west =====
    h = 1.5 / 12.0   # center-to-center offset from blue
    gwh = gw + h
    gih = gi + h

    rp1 = [pts["W2"][0] + gwh, washer_cn]
    rp2 = [pts["W2"][0] + gwh, iw8[3][1] - gih]
    rp3 = [iw2[0][0] + gih, iw8[3][1] - gih]

    iw2_anchor_r = offset_pt(iw2[0], gih, iw2_w_in)
    iw2o_anchor_r = offset_pt(iw2o[1], gih, iw2o_in)
    rp4 = list(line_isect(iw2_anchor_r, iw2_w_al, iw2o_anchor_r, iw2o_al))

    iw2s_anchor_r = offset_pt(iw2s[0], gih, iw2s_w_in)
    rp5 = list(line_isect(iw2o_anchor_r, iw2o_al, iw2s_anchor_r, iw2s_w_al))

    w9_inset_r = offset_pt(pts["W9"], gih, w9w10_in)
    rp6 = list(line_isect(iw2s_anchor_r, iw2s_w_al, w9_inset_r, w9w10_al))

    rp7 = list(offset_pt(w9_inset_r, dw_d + _DW_WIDTH / 2, w9w10_al))

    # WH break points
    wh_brk_w = list(line_isect(w9_inset_r, w9w10_al,
                               (wh_ctr[0] - h, 0), (0, 1)))
    wh_brk_e = list(line_isect(w9_inset_r, w9w10_al,
                               (wh_ctr[0] + h, 0), (0, 1)))
    wh_brk_wn = [wh_brk_w[0], wh_bl_n[1]]
    wh_brk_en = [wh_brk_e[0], wh_bl_n[1]]

    red_w_path = [rp1, rp2, rp3, rp4, rp5, rp6, wh_brk_w, wh_brk_wn]
    red_e_path = [wh_brk_en, wh_brk_e, rp7]

    # ===== 2nd red hot line (bath sink to kitchen sink) =====
    gihh = gi + 2 * h

    # North toilet east side
    iw8_al_p, _ = seg_vecs(iw8[0], iw8[1])
    iw8_n_ref = iw8[3]
    dryer_poly = appl["dryer"]["poly"]
    dryer_cx = sum(p[0] for p in dryer_poly) / 4
    dryer_cy = sum(p[1] for p in dryer_poly) / 4
    d_toilet_n_al = ((dryer_cx - iw8_n_ref[0]) * iw8_al_p[0] +
                     (dryer_cy - iw8_n_ref[1]) * iw8_al_p[1])
    toilet_n_ctr = offset_pt(iw8_n_ref, d_toilet_n_al - 4.0 / 12.0, iw8_al_p)
    toilet_n_east = offset_pt(toilet_n_ctr, 1.905 * _SVG_TO_FT, iw8_al_p)

    rr1 = [toilet_n_east[0], iw8[3][1] - gihh]
    rr2 = [iw2[0][0] + gihh, iw8[3][1] - gihh]

    iw2_anchor_rr = offset_pt(iw2[0], gihh, iw2_w_in)
    iw2o_anchor_rr = offset_pt(iw2o[1], gihh, iw2o_in)
    rr3 = list(line_isect(iw2_anchor_rr, iw2_w_al, iw2o_anchor_rr, iw2o_al))

    iw2s_anchor_rr = offset_pt(iw2s[0], gihh, iw2s_w_in)
    rr4 = list(line_isect(iw2o_anchor_rr, iw2o_al, iw2s_anchor_rr, iw2s_w_al))

    w9_inset_rr = offset_pt(pts["W9"], gihh, w9w10_in)
    rr5 = list(line_isect(iw2s_anchor_rr, iw2s_w_al, w9_inset_rr, w9w10_al))

    ks_ctr_d = dw_d - _KITCHEN_APPL_GAP - _KITCHEN_SINK_WIDTH / 2
    rr6 = list(offset_pt(w9_inset_rr, ks_ctr_d, w9w10_al))

    # T-connections to 1st red line
    rr0 = [toilet_n_east[0], iw8[3][1] - gih]
    rr7 = list(offset_pt(w9_inset_r, ks_ctr_d, w9w10_al))

    red2_path = [rr0, rr1, rr2, rr3, rr4, rr5, rr6, rr7]

    # ===== Build pipe records =====
    pipes = [
        {"name": "Utility Path", "type": "drain_pipe",
         "path": [wp1, wp2, wp3, wp4],
         "properties": {"slope": "0.25 in/ft"}},
        {"name": "Cold Supply Main", "type": "supply_pipe",
         "path": blue_path,
         "properties": {"hot_cold": "cold"}},
        {"name": "Cold WH Stub", "type": "supply_pipe",
         "path": [wh_bl, wh_bl_n],
         "properties": {"hot_cold": "cold"}},
        {"name": "Hot Supply West", "type": "supply_pipe",
         "path": red_w_path,
         "properties": {"hot_cold": "hot"}},
        {"name": "Hot Supply East", "type": "supply_pipe",
         "path": red_e_path,
         "properties": {"hot_cold": "hot"}},
        {"name": "Hot Supply 2nd", "type": "supply_pipe",
         "path": red2_path,
         "properties": {"hot_cold": "hot"}},
    ]

    # ===== Fixture positions from variant_items centroids =====
    fixture_map = {
        "Washer": "washer",
        "Toilet1": "toilet_s",
        "Toilet2": "toilet_n",
        "Util Sink": "util_sink",
        "Bath Sink": "bath_sink",
        "Fridge": "fridge",
        "Kitchen Sink": "kitchen_sink",
        "Dishwasher": "dishwasher",
        "Ice Maker": "ice_maker",
    }
    fixture_positions = {}
    for fixture_name, item_key in fixture_map.items():
        # Check variant_items first, then appliances
        item = vi.get(item_key) or appl.get(item_key)
        if item and "poly" in item:
            fixture_positions[fixture_name] = _poly_centroid(item["poly"])

    # Shower: approximate center of shower room (IW2S east face midpoint)
    iw2s_e_mid = [(iw2s[1][0] + iw2s[2][0]) / 2,
                  (iw2s[1][1] + iw2s[2][1]) / 2]
    # Offset east into the shower room by ~1.5 ft
    fixture_positions["Shower"] = list(offset_pt(iw2s_e_mid, 1.5, iw2s_e_out))

    return {"pipes": pipes, "fixture_positions": fixture_positions}


def seed_reference_plumbing(geom, wall_t, db_path=None):
    """Seed reference plumbing pipes and fixture positions into the database.

    Only runs if no supply_pipe or drain_pipe records exist yet.
    Updates fixture_connection positions from computed geometry.
    """
    existing = get_plumbing_elements(db_path)
    if any(e["type"] in ("supply_pipe", "drain_pipe") for e in existing):
        return  # Already seeded or user has drawn pipes

    ref = compute_reference_plumbing(geom, wall_t)

    # Create pipe records
    for pipe in ref["pipes"]:
        create_plumbing_element(
            pipe["type"], pipe["name"], pipe["path"],
            pipe["properties"], db_path=db_path,
        )

    # Update fixture_connection positions
    for e in existing:
        if e["type"] == "fixture_connection" and e["name"] in ref["fixture_positions"]:
            pos = ref["fixture_positions"][e["name"]]
            update_plumbing_element(e["id"], {"path": [pos]}, db_path)
