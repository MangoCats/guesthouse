"""
Recover MarkX.db furniture, appliance, and fixture placements from X/floorplan/floorplan.svg.

Strategy:
  1. Extract every furniture polygon/circle from the backup SVG (fill='rgba(100,150,200,0.2)')
  2. Convert SVG pixel coords to E/N feet using the same calibration as shared/svg.py
  3. Map each shape to its element name by label text proximity
  4. Build appropriate formula:
       - 4-pt rectangles  -> four_corner with literal [E,N] vertices (exact rotation)
       - Circles (ET, WH, clearances) -> item_circle with literal center + radius
       - Toilet (34-pt)   -> toilet_shape with literal center, facing_dir, width_dir
       - Bath sink (37-pt) -> bath_sink_shape with literal anchor, along, outward
  5. Also update element properties where MarkX used different labels ('DESK', 'OTTO')
     and variant lists where MarkX showed elements not visible in adu standard variant.
  6. Update element_formulas for all variant rows that match the standard variant shapes.
"""

import math, sqlite3, json, sys, os
import xml.etree.ElementTree as ET

_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _DIR)
import shared.svg as svg_mod

SVG_PATH = os.path.join(_DIR, "X", "floorplan", "floorplan.svg")
DST_DB   = os.path.join(_DIR, "app", "configs", "MarkX.db")

# ---------------------------------------------------------------------------
# Coordinate helpers
# ---------------------------------------------------------------------------

_s, _px, _py = svg_mod._s, svg_mod._px, svg_mod._py

def to_en(x, y):
    return (x - _px) / _s, (_py - y) / _s

def poly_centroid(pts_en):
    cx = sum(e for e,n in pts_en) / len(pts_en)
    cy = sum(n for e,n in pts_en) / len(pts_en)
    return cx, cy

def edge_len(a, b):
    return math.hypot(b[0]-a[0], b[1]-a[1])

# ---------------------------------------------------------------------------
# Parse SVG
# ---------------------------------------------------------------------------

def parse_svg():
    tree = ET.parse(SVG_PATH)
    root = tree.getroot()
    ns = '{http://www.w3.org/2000/svg}'

    polys   = list(root.iter(f'{ns}polygon'))
    texts   = list(root.iter(f'{ns}text'))
    circles = list(root.iter(f'{ns}circle'))

    # Furniture polygons: fill='rgba(100,150,200,0.2)' stroke='#4682B4'
    furn_polys = []
    for i, p in enumerate(polys):
        if p.get('fill','') == 'rgba(100,150,200,0.2)':
            pts_svg = [tuple(map(float, pt.split(',')))
                       for pt in p.get('points','').split()]
            pts_en = [to_en(x,y) for x,y in pts_svg]
            furn_polys.append({'svgi': i, 'pts': pts_en,
                               'npts': len(pts_en),
                               'center': poly_centroid(pts_en)})

    # Text labels with positions
    labels = []
    for t in texts:
        x, y = float(t.get('x', 0)), float(t.get('y', 0))
        e, n = to_en(x, y)
        txt = (t.text or '').strip()
        labels.append({'pos': (e, n), 'text': txt})

    # Circles (ET, WH, three_feet, turning_circle)
    circ_shapes = []
    for c in circles:
        cx, cy = float(c.get('cx', 0)), float(c.get('cy', 0))
        r_px   = float(c.get('r', 0))
        e, n   = to_en(cx, cy)
        r_ft   = r_px / _s
        if c.get('fill','') == 'rgba(100,150,200,0.2)':
            circ_shapes.append({'center': (e, n), 'radius': r_ft})

    return furn_polys, labels, circ_shapes

# ---------------------------------------------------------------------------
# Nearest-label matching
# ---------------------------------------------------------------------------

def nearest_label(center, labels, used, candidates=None):
    """Return text of nearest unused label to center, optionally restricted to candidates."""
    best_d, best_l = 1e9, None
    for l in labels:
        if id(l) in used:
            continue
        if candidates and l['text'] not in candidates:
            continue
        d = math.hypot(center[0]-l['pos'][0], center[1]-l['pos'][1])
        if d < best_d:
            best_d, best_l = d, l
    if best_l:
        used.add(id(best_l))
    return best_l['text'] if best_l else '?'

# ---------------------------------------------------------------------------
# Formula builders
# ---------------------------------------------------------------------------

def four_corner_formula(pts_en):
    """
    Build a four_corner formula from 4 extracted EN polygon vertices.
    Polygon vertex order from SVG (CW traversal in SVG, but N-up flip makes it CCW in EN):
    [0]=one corner, [1]=next CW in SVG=CCW in EN...
    We pass all 4 vertices as literals; the evaluator uses sw/se/ne/nw keys.
    Use the same vertex->corner assignment as recover_markx_interior.py:
    index 0->SW, 1->SE, 2->NE, 3->NW (direct polygon order).
    """
    corners = [list(map(lambda v: round(v, 6), pt)) for pt in pts_en]
    return {
        "type": "four_corner",
        "sw": corners[0],
        "se": corners[1],
        "ne": corners[2],
        "nw": corners[3],
    }

def item_circle_formula(center, radius):
    return {
        "type": "item_circle",
        "center": [round(center[0], 6), round(center[1], 6)],
        "radius": round(radius, 6),
    }

def toilet_formula(pts_en):
    """
    Extract toilet orientation from its 34-pt polygon.
    Canonical shape: back wall at y=0 with pts[0]=[-1.905,0] and pts[-1]=[+1.905,0].
    In the rendered polygon the back corners are pts[0] and pts[-1].
    The facing direction (tank->bowl) is the perpendicular pointing toward the bowl.
    """
    p0, pm1 = pts_en[0], pts_en[-1]
    back_dir = (pm1[0]-p0[0], pm1[1]-p0[1])
    back_len = math.hypot(*back_dir)
    back_dir_n = (back_dir[0]/back_len, back_dir[1]/back_len)
    back_mid = ((p0[0]+pm1[0])/2, (p0[1]+pm1[1])/2)

    # Two candidate perpendiculars; bowl is at canonical y ~ 2.224 * scale = 0.729 ft
    cx, cy = poly_centroid(pts_en)
    perp1 = (-back_dir_n[1],  back_dir_n[0])
    perp2 = ( back_dir_n[1], -back_dir_n[0])
    v = (cx-back_mid[0], cy-back_mid[1])
    facing = perp1 if (v[0]*perp1[0]+v[1]*perp1[1]) > 0 else perp2

    # The toilet_shape evaluator uses 'center' as the origin for the _TOILET_SVG
    # coordinate system: local (dx=0, dy=0) maps to center.  In _TOILET_SVG,
    # the back wall is at dy=0 and the bowl extends in the +dy direction.
    # So the formula 'center' is exactly the midpoint of the back (tank) wall edge.
    return {
        "type": "toilet_shape",
        "center": [round(back_mid[0], 6), round(back_mid[1], 6)],
        "facing_dir": [round(facing[0], 6), round(facing[1], 6)],
        "width_dir":  [round(back_dir_n[0], 6), round(back_dir_n[1], 6)],
    }

def bath_sink_formula(pts_en, sink_type):
    """
    Extract bath sink orientation from its 37-pt polygon.
    Canonical shape: back edge at y=0, pts[0]=[+width/2, 0], pts[1]=[-width/2, 0].
    Anchor = center of back edge.
    along = direction from pts[0] to pts[1] (along the wall).
    outward = perpendicular pointing toward bowl (toward centroid).
    """
    p0, p1 = pts_en[0], pts_en[1]
    back_dir = (p1[0]-p0[0], p1[1]-p0[1])
    back_len = math.hypot(*back_dir)
    along_n  = (back_dir[0]/back_len, back_dir[1]/back_len)
    anchor   = ((p0[0]+p1[0])/2, (p0[1]+p1[1])/2)

    cx, cy = poly_centroid(pts_en)
    perp1 = (-along_n[1],  along_n[0])
    perp2 = ( along_n[1], -along_n[0])
    v = (cx-anchor[0], cy-anchor[1])
    outward = perp1 if (v[0]*perp1[0]+v[1]*perp1[1]) > 0 else perp2

    pfx = "BATH_SINK_L" if sink_type == "bath_sink_l" else "BATH_SINK"
    return {
        "type": "bath_sink_shape",
        "anchor":  [round(anchor[0], 6), round(anchor[1], 6)],
        "along":   [round(along_n[0], 6), round(along_n[1], 6)],
        "outward": [round(outward[0], 6), round(outward[1], 6)],
        "length":  {"const": f"{pfx}_LENGTH"},
        "depth":   {"const": f"{pfx}_DEPTH"},
    }

# ---------------------------------------------------------------------------
# Assignment map: label text -> candidate element names (in priority order)
# When a label appears multiple times we assign in spatial order.
# ---------------------------------------------------------------------------

# Manual mapping built from spatial analysis of backup SVG.
# Each entry: (poly_svgi, element_name, formula_type_hint)
# Determined by matching label text proximity and spatial context.
POLY_TO_ELEM = {
    139: ("bath_sink",      "bath_sink_shape"),   # BATH SINK, 37-pt
    140: ("bath_sink_l",    "bath_sink_shape"),   # BATH SINK L, 37-pt
    141: ("bed",            "four_corner"),        # KING BED, ang=33.8
    142: ("chair",          "four_corner"),        # CHAIR near bedroom
    143: ("north_counter",  "four_corner"),        # COUNTER, 6x2, axis-aligned
    144: ("work_counter",   "four_corner"),        # DESK, 5x2.5 - MarkX label
    145: ("tonstad_chair",  "four_corner"),        # CHAIR near desk
    146: ("dresser",        "four_corner"),        # DRESSER, ang=33.8
    147: ("dryer",          "four_corner"),        # DRYER, axis-aligned
    148: ("fridge",         "four_corner"),        # FRIDGE, axis-aligned
    149: ("kitchen_sink",   "four_corner"),        # SINK (kitchen area)
    150: ("loveseat",       "four_corner"),        # LOVESEAT, ang=45
    151: ("dining_table",   "four_corner"),        # TABLE
    152: ("counter",        "four_corner"),        # COUNTER (stove area)
    153: ("nordviken",      "four_corner"),        # OTTO - square ottoman; MarkX label
    154: ("skogsta_bench",  "four_corner"),        # BENCH
    155: ("stockholm_sofa", "four_corner"),        # SOFA, ang=140
    156: ("stove",          "four_corner"),        # STOVE
    157: ("toilet_s",       "toilet_shape"),       # TOILET (south bath, ang~33.8)
    158: ("toilet_n",       "toilet_shape"),       # TOILET (north bath)
    159: ("dining_chair_1", "four_corner"),        # CHAIR at dining table
    160: ("dining_chair_2", "four_corner"),        # CHAIR at dining table
    161: ("loveseat2",      "four_corner"),        # CHAIR at dining table (MarkX repurposed)
    162: ("hamper",         "four_corner"),        # CHAIR at dining table (MarkX repurposed)
    163: ("washer",         "four_corner"),        # WASHER, ang=123.7
    164: ("microwave",      "four_corner"),        # MICRO
}

# Circles: matched by proximity to text labels
# (6.622,2.436) r=0.824 -> ET  (standard variant)
# (0.223,10.336) r=0.824 -> ET (second ET in MarkX, use et_east mapped to standard)
# (-15.912,13.407) r=1.167 -> WH
# (-4.452,-6.203) r=1.501 -> THREE FEET
# (2.539,-7.549) r=2.496 -> TURNING CIRCLE
CIRCLE_ELEM = [
    # (approx_center_E, approx_center_N, element_name, variants_override)
    (6.622,   2.436,   "et",              None),
    (0.223,   10.336,  "et_east",         ["standard", "daybed"]),  # MarkX showed in standard
    (-15.912, 13.407,  "water_heater",    None),
    (-4.452,  -6.203,  "three_feet",      None),
    (2.539,   -7.549,  "turning_circle",  None),
]

# Label overrides: elements whose 'label' property differs from adu.db in MarkX
# key = element_name, value = new label string
LABEL_OVERRIDES = {
    "work_counter": "DESK",    # MarkX showed DESK instead of COUNTER
    "nordviken":    "OTTO",    # MarkX showed OTTO (ottoman) instead of TABLE
    "loveseat2":    "CHAIR",   # repurposed as dining chair in MarkX
    "hamper":       "CHAIR",   # repurposed as dining chair in MarkX
}

# Variants override: elements whose variants list differs from adu.db in MarkX standard view
# key = element_name, value = new variants list
VARIANTS_OVERRIDES = {
    "et_east":        ["standard", "daybed"],
    "loveseat2":      ["standard", "minik", "daybed"],
    "hamper":         ["standard", "minik", "daybed"],
    "three_feet":     ["standard", "bare"],    # MarkX shows in standard view
    "turning_circle": ["standard", "bare"],    # MarkX shows in standard view
}

# Elements whose existing formula_json variant column must be NULL so the formula
# loads for all variants (including 'standard').  For elements whose existing rows
# are all variant-specific (e.g. variant='daybed'), we also insert a variant=NULL row.
FORCE_NULL_VARIANT = {"et_east", "three_feet", "turning_circle"}

# ---------------------------------------------------------------------------
# Build formulas from SVG
# ---------------------------------------------------------------------------

def build_formulas(furn_polys, circ_shapes):
    """Return dict: element_name -> formula dict."""
    # Index furn_polys by svgi for O(1) lookup
    by_svgi = {fp['svgi']: fp for fp in furn_polys}
    formulas = {}

    # --- Polygon elements ---
    for svgi, (elem_name, ftype) in POLY_TO_ELEM.items():
        fp = by_svgi.get(svgi)
        if fp is None:
            print(f"WARNING: poly[{svgi}] not found in SVG")
            continue
        pts = fp['pts']
        npts = fp['npts']

        if ftype == "four_corner":
            assert npts == 4, f"Expected 4-pt polygon for {elem_name}, got {npts}"
            formula = four_corner_formula(pts)
        elif ftype == "toilet_shape":
            assert npts == 34, f"Expected 34-pt polygon for {elem_name}, got {npts}"
            formula = toilet_formula(pts)
        elif ftype == "bath_sink_shape":
            assert npts == 37, f"Expected 37-pt polygon for {elem_name}, got {npts}"
            formula = bath_sink_formula(pts, elem_name)
        else:
            raise ValueError(f"Unknown formula type {ftype} for {elem_name}")

        formulas[elem_name] = formula
        print(f"  {elem_name:20s} <- poly[{svgi}] {ftype} npts={npts}")

    # --- Circle elements ---
    for approx_e, approx_n, elem_name, _ in CIRCLE_ELEM:
        # Find the matching circle by proximity to expected center
        best = min(circ_shapes,
                   key=lambda c: math.hypot(c['center'][0]-approx_e,
                                            c['center'][1]-approx_n))
        dist = math.hypot(best['center'][0]-approx_e, best['center'][1]-approx_n)
        if dist > 0.5:
            print(f"WARNING: {elem_name}: no circle within 0.5ft of ({approx_e},{approx_n}), "
                  f"nearest is {dist:.2f}ft away")
        formula = item_circle_formula(best['center'], best['radius'])
        formulas[elem_name] = formula
        print(f"  {elem_name:20s} <- circle center=({best['center'][0]:.3f},{best['center'][1]:.3f}) r={best['radius']:.3f}ft")

    return formulas

# ---------------------------------------------------------------------------
# Write to DB
# ---------------------------------------------------------------------------

def write_furniture(formulas):
    conn = sqlite3.connect(DST_DB)

    # Update element_formulas: for each element, update ALL variant rows to the
    # new literal formula (position is the same regardless of variant since we're
    # hardcoding the EN coordinates).
    # For elements that only appear in one variant (e.g. et=standard only),
    # we update just that row.
    updated = 0
    for elem_name, formula in formulas.items():
        fj = json.dumps(formula)
        rows = conn.execute(
            "SELECT id, variant FROM element_formulas WHERE element_name=? AND param_name='position'",
            (elem_name,)
        ).fetchall()
        if not rows:
            # Insert a new row with variant=NULL
            conn.execute(
                "INSERT INTO element_formulas (element_name, param_name, formula_json, variant) "
                "VALUES (?, 'position', ?, NULL)",
                (elem_name, fj)
            )
            print(f"  INSERT {elem_name} (no existing formula, variant=NULL)")
            updated += 1
        else:
            for row_id, variant in rows:
                conn.execute(
                    "UPDATE element_formulas SET formula_json=? WHERE id=?",
                    (fj, row_id)
                )
            updated += len(rows)
            # For elements that must be visible in 'standard' variant, ensure there
            # is a variant=NULL (universal) row in addition to any variant-specific rows.
            # get_all_formulas for variant='standard' loads only NULL or 'standard' rows.
            if elem_name in FORCE_NULL_VARIANT:
                null_row = conn.execute(
                    "SELECT id FROM element_formulas WHERE element_name=? AND param_name='position' AND variant IS NULL",
                    (elem_name,)
                ).fetchone()
                if not null_row:
                    conn.execute(
                        "INSERT INTO element_formulas (element_name, param_name, formula_json, variant) "
                        "VALUES (?, 'position', ?, NULL)",
                        (elem_name, fj)
                    )
                    print(f"  INSERT {elem_name} variant=NULL (for standard visibility)")
                    updated += 1

    print(f"Updated {updated} element_formula rows")

    # Apply label overrides
    for elem_name, new_label in LABEL_OVERRIDES.items():
        row = conn.execute(
            "SELECT properties FROM elements WHERE name=?", (elem_name,)
        ).fetchone()
        if row:
            props = json.loads(row[0])
            old_label = props.get('label', '')
            props['label'] = new_label
            conn.execute(
                "UPDATE elements SET properties=? WHERE name=?",
                (json.dumps(props), elem_name)
            )
            print(f"  label override: {elem_name}: {old_label!r} -> {new_label!r}")

    # Apply variants overrides
    for elem_name, new_variants in VARIANTS_OVERRIDES.items():
        row = conn.execute(
            "SELECT properties FROM elements WHERE name=?", (elem_name,)
        ).fetchone()
        if row:
            props = json.loads(row[0])
            old_v = props.get('variants', [])
            props['variants'] = new_variants
            conn.execute(
                "UPDATE elements SET properties=? WHERE name=?",
                (json.dumps(props), elem_name)
            )
            print(f"  variants override: {elem_name}: {old_v} -> {new_variants}")

    conn.commit()
    conn.close()
    print(f"Done. Updated {len(formulas)} elements in {DST_DB}")

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print(f"Parsing {SVG_PATH}...")
    furn_polys, labels, circ_shapes = parse_svg()
    print(f"  Found {len(furn_polys)} furniture polygons, {len(circ_shapes)} circles")

    print("\nBuilding formulas:")
    formulas = build_formulas(furn_polys, circ_shapes)

    print(f"\nWriting {len(formulas)} element formulas to {DST_DB}...")
    write_furniture(formulas)

if __name__ == "__main__":
    main()
