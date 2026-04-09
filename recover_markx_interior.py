"""
Recover MarkX interior walls (IW), outer openings (O), and rough openings (RO)
from X/walls/all_walls.svg and X/floorplan/floorplan_bare.svg.

Sources used:
  all_walls.svg poly[111-119]  — IW body 4-corner polygons
  all_walls.svg poly[120-125]  — RO cutout 4-corner polygons
  all_walls.svg poly[14,15,32,33,50,59,72,89,90] — glazing (O-series positions/widths)
  all_walls.svg text table     — IW names, lengths, RO widths
  walls.svg text labels        — O-series name→label-position mapping

Strategy:
  IW walls  → four_corner formula with literal [E,N] corners (sw/se/ne/nw)
  O openings → wall_opening formula: ref_point=[label E,N], width from glazing polygon
  RO openings → wall_opening formula: outer/inner from parent IW, ref_point=center,
                width from table
"""

import math, json, sqlite3, sys, os
import xml.etree.ElementTree as ET

_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _DIR)
import shared.svg as svg_mod

DST_DB           = os.path.join(_DIR, "app", "configs", "MarkX.db")
ALL_WALLS_SVG    = os.path.join(_DIR, "X", "walls", "all_walls.svg")
WALLS_SVG        = os.path.join(_DIR, "X", "walls", "walls.svg")
FLOORPLAN_SVG    = os.path.join(_DIR, "X", "floorplan", "floorplan_bare.svg")

_s, _px, _py = svg_mod._s, svg_mod._px, svg_mod._py

def svg_to_en(x, y):
    return (x - _px) / _s, (_py - y) / _s

def get_polygon_corners(svg_path, index):
    """Extract E,N corners from a <polygon> by index."""
    tree = ET.parse(svg_path)
    polys = list(tree.getroot().iter("{http://www.w3.org/2000/svg}polygon"))
    p = polys[index]
    tokens = p.get("points", "").split()
    corners = []
    for tok in tokens:
        xy = tok.split(",")
        if len(xy) == 2:
            corners.append(svg_to_en(float(xy[0]), float(xy[1])))
    return corners

def poly_centroid(corners):
    cx = sum(e for e,n in corners) / len(corners)
    cn = sum(n for e,n in corners) / len(corners)
    return cx, cn

def assign_sw_se_ne_nw(corners):
    """
    Assign SW/SE/NE/NW to a 4-corner polygon using the polygon vertex order.
    The SVG polygon vertices run around the perimeter; we preserve that order.
    This avoids duplicate-assignment bugs on thin diagonal rectangles where two
    corners share the same geographic quadrant.

    Convention: corners[0]→sw, corners[1]→se, corners[2]→ne, corners[3]→nw.
    The RO formula then references the outer face as corner-index 0→1 and the
    inner face as corner-index 3→2 (same direction, opposite side).
    """
    assert len(corners) == 4
    return {
        "SW": list(corners[0]),
        "SE": list(corners[1]),
        "NE": list(corners[2]),
        "NW": list(corners[3]),
    }

def proj_t(seg_start, seg_end, pt):
    """Parametric t of pt projected onto segment start→end."""
    dx, dy = seg_end[0]-seg_start[0], seg_end[1]-seg_start[1]
    L2 = dx*dx + dy*dy
    if L2 < 1e-12:
        return 0.0
    return ((pt[0]-seg_start[0])*dx + (pt[1]-seg_start[1])*dy) / L2

def seg_len(a, b):
    return math.hypot(b[0]-a[0], b[1]-a[1])

def opening_width_from_glazing(corners):
    """
    Glazing polygon is a 4-pt rectangle spanning the full wall thickness.
    The opening width is the diagonal across the opening (the long direction).
    For axis-aligned openings: max of width/height.
    For diagonal openings: diagonal of bbox.
    """
    es = [e for e,n in corners]
    ns = [n for e,n in corners]
    bw = (max(es)-min(es))*12  # inches
    bh = (max(ns)-min(ns))*12
    # The opening extent is max of the two dimensions (thin dimension = wall thickness)
    # But for diagonal openings both bw and bh are non-zero; use the full perimeter diagonal.
    # Compute all 4 edge lengths and take the longer pair average.
    edges = [seg_len(corners[i], corners[(i+1)%4]) for i in range(4)]
    edges.sort(reverse=True)
    # Two longer edges = the opening width direction
    return (edges[0] + edges[1]) / 2 * 12  # inches

# ─────────────────────────────────────────────────────────
# Parsed data from SVGs (extracted offline, stored as constants)
# ─────────────────────────────────────────────────────────

# IW body polygon corners from all_walls.svg [111-119]
# Keyed by IW name; corners in polygon order (not yet SW/SE/NE/NW sorted)
IW_POLY_INDEX = {
    "IW15": 111,
    "IW17": 112,
    "IW25": 113,
    "IW31": 114,
    "IW32": 115,
    "IW35": 116,
    "IW36": 117,
    "IW37": 118,
    "IW38": 119,
}

# IW metadata from all_walls.svg table
IW_META = {
    #  name   thk_in  orientation  len_in
    "IW15": {"thk": 6, "orientation": "R", "len_in": 156.0},   # 13'0"
    "IW17": {"thk": 6, "orientation": "R", "len_in": 156.0},   # 13'0"
    "IW25": {"thk": 6, "orientation": "R", "len_in": 173.0},   # 14'5"
    "IW31": {"thk": 6, "orientation": "H", "len_in":  34.0},   # 2'10"
    "IW32": {"thk": 6, "orientation": "R", "len_in":  64.0},   # 5'4"
    "IW35": {"thk": 6, "orientation": "H", "len_in":  89.0},   # 7'5"
    "IW36": {"thk": 6, "orientation": "R", "len_in":  74.5},   # 6'2.5"
    "IW37": {"thk": 6, "orientation": "R", "len_in":  98.0},   # 8'2"
    "IW38": {"thk": 6, "orientation": "R", "len_in":  40.0},   # 3'4"
}

# RO cutout polygon indices and assignments
RO_DATA = [
    # name   poly_idx  parent_iw  width_in
    ("RO11",  120,     "IW15",     62),
    ("RO12",  121,     "IW15",     62),
    ("RO8",   122,     "IW17",     38),
    ("RO17",  123,     "IW25",     36),
    ("RO18",  124,     "IW35",     38),
    ("RO19",  125,     "IW37",     36),
]

# Outer opening glazing polygon indices and label positions from walls.svg
# Format: (O_name, glazing_poly_index_in_all_walls, label_E, label_N, opening_type)
O_DATA = [
    ("O29",  14,   -17.671,   9.736,  "casement"),  # west wall upper
    ("O33",  15,   -17.671,  -8.347,  "casement"),  # west wall lower
    ("O22",  32,    -9.342,  11.940,  "window"),     # north wall left
    ("O25",  33,    -3.637,  11.940,  "window"),     # north wall right
    ("O26",  50,     7.300,  12.172,  "window"),     # NE area
    ("O21",  59,    16.582,   1.106,  "window"),     # SE wall
    ("O30",  72,     5.104, -13.983,  "window"),     # south area
    ("O32",  89,    -4.461,  -8.759,  "window"),     # SW area
    ("O31",  90,    -8.990, -11.787,  "window"),     # SW area
]

# ─────────────────────────────────────────────────────────
# Outer wall segment mapping for O-series
# Derived from MarkX outline_segs (line segments only, arcs excluded)
# Format: (outer_start_name, outer_end_name, inner_start_name, inner_end_name)
# ─────────────────────────────────────────────────────────

OUTER_SEGS = [
    # West wall
    ("F01", "F02", "W01", "W02"),
    # North wall segments
    ("F05", "F06", "W05", "W06"),
    ("F07", "F08", "W07", "W08"),
    ("F09", "F10", "W09", "W10"),
    ("F11", "F14", "W11", "W14"),
    # SE segments
    ("F15", "F16", "W15", "W16"),
    ("F17", "F18", "W17", "W18"),
    ("F19", "F20", "W19", "W20"),
    ("F21", "F22", "W21", "W22"),
    # SW segments
    ("F23", "F24", "W23", "W24"),
    ("F24b","F26",  "W24b","W26"),
]


def find_best_outer_seg(label_e, label_n, pts):
    """
    Project the opening label position onto each line segment.
    Return the segment where t ∈ [0,1] and perpendicular distance is minimal.
    """
    best = None
    best_score = 9e9
    for os_name, oe_name, is_name, ie_name in OUTER_SEGS:
        o_s = pts[os_name]
        o_e = pts[oe_name]
        t = proj_t(o_s, o_e, (label_e, label_n))
        if not (-0.5 <= t <= 1.5):
            continue
        # Perpendicular distance from label to segment
        dx, dy = o_e[0]-o_s[0], o_e[1]-o_s[1]
        L = math.hypot(dx, dy)
        if L < 1e-9:
            continue
        perp = abs((label_e-o_s[0])*(-dy) + (label_n-o_s[1])*dx) / L
        score = perp + max(0, t-1)*L + max(0, -t)*L
        if score < best_score:
            best_score = score
            best = (os_name, oe_name, is_name, ie_name, t)
    return best


def make_wall_opening_formula(o_s, o_e, i_s, i_e, pts_dict, glazing_corners):
    """
    Build a wall_opening formula using ref_point mode.
    The ref_point is the midpoint of the glazing polygon projected to center on the segment.
    """
    gx, gn = poly_centroid(glazing_corners)
    # Compute width from glazing polygon diagonal
    width_ft = opening_width_from_glazing(glazing_corners) / 12.0

    # Project glazing centroid onto outer segment to get t_center
    seg_s = pts_dict[o_s]
    seg_e = pts_dict[o_e]
    t_center = proj_t(seg_s, seg_e, (gx, gn))
    seg_length = seg_len(seg_s, seg_e)

    # Use gap from start (from_end=False)
    gap_ft = t_center * seg_length - width_ft / 2.0

    return {
        "type": "wall_opening",
        "outer_start": o_s,
        "outer_end":   o_e,
        "inner_start": i_s,
        "inner_end":   i_e,
        "from_end":    False,
        "gap":         round(gap_ft, 6),
        "width":       round(width_ft, 6),
    }


def make_ro_formula(parent_iw_name, iw_raw_corners, ro_corners, width_in):
    """
    Build a wall_opening formula for an RO using the parent IW polygon.
    Uses integer corner indices (0,1 = outer face; 3,2 = inner face, same direction)
    so the formula is unambiguous regardless of geographic orientation.
    Width comes from the all_walls table (authoritative).
    """
    ro_cx, ro_cn = poly_centroid(ro_corners)

    # Outer face: corners[0] → corners[1]
    outer_start_pt = iw_raw_corners[0]
    outer_end_pt   = iw_raw_corners[1]

    ro_width_ft = width_in / 12.0

    # Project RO centroid onto the outer face
    t_center = proj_t(outer_start_pt, outer_end_pt, (ro_cx, ro_cn))
    seg_length = seg_len(outer_start_pt, outer_end_pt)
    gap_ft = t_center * seg_length - ro_width_ft / 2.0

    return {
        "type":        "wall_opening",
        "outer_start": {"element": parent_iw_name, "corner": 0},
        "outer_end":   {"element": parent_iw_name, "corner": 1},
        "inner_start": {"element": parent_iw_name, "corner": 3},
        "inner_end":   {"element": parent_iw_name, "corner": 2},
        "from_end":    False,
        "gap":         round(gap_ft, 6),
        "width":       round(ro_width_ft, 6),
    }


def build_iw_four_corner(corners_named):
    """Build a four_corner formula with literal [E,N] corners."""
    return {
        "type": "four_corner",
        "sw": corners_named["SW"],
        "se": corners_named["SE"],
        "ne": corners_named["NE"],
        "nw": corners_named["NW"],
    }


def write_interior(db_path, pts_dict):
    conn = sqlite3.connect(db_path)

    # ── 1. Delete all existing IW/RO/O elements and their formulas ──
    for prefix in ("IW", "RO", "O"):
        conn.execute(
            "DELETE FROM element_formulas WHERE element_name LIKE ?", (prefix+"%",))
        conn.execute(
            "DELETE FROM formula_deps WHERE element_name LIKE ?", (prefix+"%",))
        conn.execute(
            "DELETE FROM elements WHERE name LIKE ?", (prefix+"%",))
    conn.execute("DELETE FROM doors")
    conn.execute("DELETE FROM inner_wall_overrides")
    print("Cleared existing IW/RO/O elements, doors, inner_wall_overrides.")

    # ── 2. Insert IW walls ──
    print("\nInserting IW walls:")
    iw_raw_corners_cache = {}  # polygon vertex order, for RO projection
    iw_corners_cache = {}      # named corners, for four_corner formula
    for iw_name, poly_idx in IW_POLY_INDEX.items():
        meta = IW_META[iw_name]
        raw_corners = get_polygon_corners(ALL_WALLS_SVG, poly_idx)
        corners = assign_sw_se_ne_nw(raw_corners)
        iw_raw_corners_cache[iw_name] = raw_corners
        iw_corners_cache[iw_name] = corners

        thk_const = "WALL_6IN" if meta["thk"] == 6 else "WALL_4IN"
        props = json.dumps({
            "thickness_constant": thk_const,
            "orientation": meta["orientation"],
            "prop_constants": {},
        })
        conn.execute(
            "INSERT OR REPLACE INTO elements (type, name, properties) VALUES (?,?,?)",
            ("wall", iw_name, props))

        formula = build_iw_four_corner(corners)
        conn.execute(
            "INSERT OR REPLACE INTO element_formulas (element_name, param_name, formula_json)"
            " VALUES (?,?,?)",
            (iw_name, "position", json.dumps(formula)))

        cx = (corners["SW"][0] + corners["NE"][0]) / 2
        cn = (corners["SW"][1] + corners["NE"][1]) / 2
        print(f"  {iw_name}: center=({cx:.3f},{cn:.3f})"
              f"  SW={corners['SW']}  NE={corners['NE']}")

    # ── 3. Insert RO openings ──
    print("\nInserting RO openings:")
    for ro_name, poly_idx, parent_iw, width_in in RO_DATA:
        ro_corners = get_polygon_corners(ALL_WALLS_SVG, poly_idx)
        parent_raw = iw_raw_corners_cache[parent_iw]
        formula = make_ro_formula(parent_iw, parent_raw, ro_corners, width_in)

        props = json.dumps({"wall_name": parent_iw, "orientation": "R"})
        conn.execute(
            "INSERT OR REPLACE INTO elements (type, name, properties) VALUES (?,?,?)",
            ("opening", ro_name, props))
        conn.execute(
            "INSERT OR REPLACE INTO element_formulas (element_name, param_name, formula_json)"
            " VALUES (?,?,?)",
            (ro_name, "position", json.dumps(formula)))

        ro_cx, ro_cn = poly_centroid(ro_corners)
        print(f"  {ro_name}: parent={parent_iw} width={width_in}\" center=({ro_cx:.3f},{ro_cn:.3f})")

    # ── 4. Insert O outer openings ──
    print("\nInserting O outer openings:")
    for o_name, glazing_idx, lbl_e, lbl_n, o_type in O_DATA:
        glazing_corners = get_polygon_corners(ALL_WALLS_SVG, glazing_idx)
        seg = find_best_outer_seg(lbl_e, lbl_n, pts_dict)
        if seg is None:
            print(f"  {o_name}: WARNING — could not find matching outer segment")
            continue
        os_name, oe_name, is_name, ie_name, t = seg
        formula = make_wall_opening_formula(
            os_name, oe_name, is_name, ie_name, pts_dict, glazing_corners)

        props = json.dumps({
            "seg_start": os_name,
            "seg_end":   oe_name,
            "opening_type": o_type,
        })
        conn.execute(
            "INSERT OR REPLACE INTO elements (type, name, properties) VALUES (?,?,?)",
            ("opening", o_name, props))
        conn.execute(
            "INSERT OR REPLACE INTO element_formulas (element_name, param_name, formula_json)"
            " VALUES (?,?,?)",
            (o_name, "position", json.dumps(formula)))

        width_in = opening_width_from_glazing(glazing_corners)
        gx, gn = poly_centroid(glazing_corners)
        print(f"  {o_name}: seg={os_name}->{oe_name}  t={t:.3f}  width={width_in:.1f}\"  "
              f"center=({gx:.3f},{gn:.3f})")

    conn.commit()
    conn.close()
    print(f"\nWrote interior geometry to {db_path}")


def main():
    # Load MarkX geometry to get F-series and W-series points
    from app.database import get_constants_dict, get_outline_chain
    from app.engine import compute_geometry

    db_path = DST_DB
    constants = get_constants_dict(db_path)
    chain_rows = get_outline_chain(db_path)
    print("Computing MarkX geometry for segment endpoints...")
    gd = compute_geometry(constants, chain_rows=chain_rows, db_path=db_path)
    pts = gd.get("points", {})

    # Verify key points exist
    missing = [n for n in ["F01","F02","F05","F06","W01","W02","W05","W06"] if n not in pts]
    if missing:
        print(f"WARNING: missing points: {missing}")

    write_interior(db_path, pts)
    print("\nDone. Reload MarkX in the app to verify interior geometry.")


if __name__ == "__main__":
    main()
