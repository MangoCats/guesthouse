"""
Recover MarkX.db outline_chain from X/survey/path_area.svg backup.

The survey SVG was generated in-process from the MarkX DB and contains:
  - 12 F-series line segments as <line> elements [14]-[25] (stroke=#333)
  - 12 F-series arc segments as <polyline> elements [8]-[19] (stroke=#333, sw=2.0, 13pts)
  - Bearing/distance text labels in chain order (F-series labels only, after the 7 survey traverse labels)
  - Arc R/sweep text labels in chain order

Strategy:
  1. Parse the SVG to verify the 24-segment interleaved chain (verified: perfect endpoint closure)
  2. Use the embedded text labels for exact numeric values (more accurate than pixel-fitted values)
  3. Assign F-series end_names by SVG rendering order (confirmed matches chain order)
  4. Determine CW/CCW direction from circle-fit cross-product on arc polylines
  5. Copy adu.db to MarkX.db, replace outline_chain, clear undo_history
"""

import math, shutil, sqlite3, os, sys
import xml.etree.ElementTree as ET
import numpy as np

_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _DIR)
import shared.svg as svg_mod

SVG_PATH  = os.path.join(_DIR, "X", "survey", "path_area.svg")
SRC_DB    = os.path.join(_DIR, "app", "adu.db")
DST_DB    = os.path.join(_DIR, "app", "configs", "MarkX.db")

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def svg_to_en(x, y):
    _s, _px, _py = svg_mod._s, svg_mod._px, svg_mod._py
    return (x - _px) / _s, (_py - y) / _s

def parse_dist_ft(txt):
    """Parse "25' 4\"" or "0' 2\"" → feet as float."""
    txt = txt.strip().replace('"', '').replace("'", "' ")
    parts = txt.split("'")
    feet = int(parts[0].strip())
    inches = float(parts[1].strip()) if len(parts) > 1 and parts[1].strip() else 0.0
    return feet + inches / 12.0

def parse_bearing(txt):
    """Parse "270.00°" or "37.98°" → float degrees."""
    return float(txt.replace('°', '').strip())

def parse_r_inches(txt):
    """Parse "R=16\"" → float inches."""
    return float(txt.replace('R=', '').replace('"', '').strip())

def parse_sweep_deg(txt):
    """Parse "90.0°" → float degrees."""
    return float(txt.replace('°', '').strip())

def fit_circle(pts):
    xs = np.array([p[0] for p in pts]); ys = np.array([p[1] for p in pts])
    A = np.column_stack([2*xs, 2*ys, np.ones(len(pts))])
    b = xs**2 + ys**2
    result, _, _, _ = np.linalg.lstsq(A, b, rcond=None)
    cx, cy = result[0], result[1]
    r = math.sqrt(result[2] + cx**2 + cy**2)
    return cx, cy, r

def arc_direction(pts_en):
    """Determine CW or CCW from midpoint cross-product (standard E,N coords, N-up)."""
    start, mid, end = pts_en[0], pts_en[len(pts_en)//2], pts_en[-1]
    v1 = (mid[0]-start[0], mid[1]-start[1])
    v2 = (end[0]-mid[0],   end[1]-mid[1])
    cross = v1[0]*v2[1] - v1[1]*v2[0]
    return 'CCW' if cross > 0 else 'CW'

def sweep_name_for(deg):
    """Return symbolic sweep_name string for common angles, else numeric."""
    if abs(deg - 90.0) < 0.1:   return '_PI_2'
    if abs(deg - 180.0) < 0.1:  return '_PI'
    return str(math.radians(deg))

# ---------------------------------------------------------------------------
# Parse SVG
# ---------------------------------------------------------------------------

def parse_svg():
    tree = ET.parse(SVG_PATH)
    root = tree.getroot()
    ns = '{http://www.w3.org/2000/svg}'

    lines_elem    = list(root.iter(f'{ns}line'))
    polylines_all = list(root.iter(f'{ns}polyline'))
    g_blocks      = list(root.iter(f'{ns}g'))
    texts         = list(root.iter(f'{ns}text'))

    # --- Line segments: elements [14]-[25] (stroke=#333, F-series outline) ---
    fseg_lines = []
    for ln in lines_elem[14:26]:
        x1, y1 = float(ln.get('x1')), float(ln.get('y1'))
        x2, y2 = float(ln.get('x2')), float(ln.get('y2'))
        fseg_lines.append({'start': svg_to_en(x1,y1), 'end': svg_to_en(x2,y2)})

    # --- Arc polylines: [8]-[19] (stroke=#333, sw=2.0, 13 pts) ---
    fseg_arcs = []
    for pl in polylines_all[8:20]:
        pts_svg = [tuple(map(float, p.split(','))) for p in pl.get('points','').split()]
        pts_en  = [svg_to_en(x,y) for x,y in pts_svg]
        cx, cy, r = fit_circle(pts_en)
        direction = arc_direction(pts_en)
        fseg_arcs.append({
            'center':    (cx, cy),
            'radius_fit': r,
            'start':     pts_en[0],
            'end':       pts_en[-1],
            'direction': direction,
        })

    # --- Bearing/dist text labels: the first 7 are survey traverse; F-series start at [7] ---
    import re
    brg_dist_labels = []  # list of (bearing_deg, dist_ft) in chain order
    for g in g_blocks:
        t = g.get('transform','')
        if 'translate' in t and 'rotate' in t:
            children = list(g)
            if len(children) == 2:
                t1 = (children[0].text or '').strip()
                t2 = (children[1].text or '').strip()
                # F-series labels have decimal degree format "270.00°"
                if re.match(r'\d+\.\d+°', t1):
                    brg_dist_labels.append((parse_bearing(t1), parse_dist_ft(t2)))

    assert len(brg_dist_labels) == 12, f"Expected 12 F-series line labels, got {len(brg_dist_labels)}"

    # --- Arc R/sweep text labels ---
    arc_labels = []  # list of (r_ft, sweep_deg) in chain order
    text_vals  = [(float(t.get('x',0)), float(t.get('y',0)), (t.text or '').strip()) for t in texts]
    i = 0
    while i < len(text_vals):
        _, _, txt = text_vals[i]
        if txt.startswith('R='):
            nxt_txt = text_vals[i+1][2] if i+1 < len(text_vals) else '?'
            r_in    = parse_r_inches(txt)
            sw_deg  = parse_sweep_deg(nxt_txt)
            arc_labels.append((r_in / 12.0, sw_deg))
            i += 2
        else:
            i += 1

    assert len(arc_labels) == 12, f"Expected 12 arc labels, got {len(arc_labels)}"

    return fseg_lines, fseg_arcs, brg_dist_labels, arc_labels

# ---------------------------------------------------------------------------
# Build chain rows
# ---------------------------------------------------------------------------

# End-names for each of the 24 segments in chain order:
# Chain: L→Ar→L→Ar→...  (L14, Arc8, L15, Arc9, ..., L25, Arc19)
# Verified by SVG position matching.
END_NAMES = [
    'F26',  # L14 end (270° westward, SW approach)
    'F01',  # Arc8 end  (CW 90°, SW corner)
    'F02',  # L15 end  (0° north, west wall)
    'F03',  # Arc9 end  (CW 180°, NW U-turn)
    'F04',  # L16 end  (180° south, tiny step)
    'F05',  # Arc10 end (CCW 90°)
    'F06',  # L17 end  (90° east, north shelf)
    'F07',  # Arc11 end (CCW 52°)
    'F08',  # L18 end  (37.98°)
    'F09',  # Arc12 end (CW 52°)
    'F10',  # L19 end  (90°, tiny step, essentially same as F09)
    'F11',  # Arc13 end (CW 37.9°)
    'F14',  # L20 end  (127.88°, SE wall approach)
    'F15',  # Arc14 end (CW 18.4°, large R)
    'F16',  # L21 end  (146.31°, SE wall)
    'F17',  # Arc15 end (CW 90°, SE corner)
    'F18',  # L22 end  (236.31°, south wall going west)
    'F19',  # Arc16 end (CW 71.6°)
    'F20',  # L23 end  (307.87°)
    'F21',  # Arc17 end (CW 18.4°, large R)
    'F22',  # L24 end  (326.31°)
    'F23',  # Arc18 end (CCW 90°)
    'F24',  # L25 end  (236.31°)
    'F24b', # Arc19 end (CW 33.7°) — chain starts from F24b
]

# Center names: C + start_point_name for each arc.
# Arc start = end_name of previous segment.
# prev end_names for arcs (odd indices): F26,F02,F03,F06,F07,F09,F14,F16,F18,F20,F22,F24
CENTER_NAMES = ['C26','C02','C04','C06','C08','C10','C14','C16','C18','C20','C22','C24']

# flex: west wall (L15, seq=2) adjusts distance for closure;
#       closing arc (Arc19, seq=23) adjusts sweep for closure.
FLEX = {2: 'distance', 23: 'sweep'}

def build_chain_rows(fseg_lines, fseg_arcs, brg_dist_labels, arc_labels):
    rows = []
    arc_idx = 0
    line_idx = 0
    name_idx = 0
    center_idx = 0

    for i in range(24):
        end_name = END_NAMES[i]
        flex     = FLEX.get(i, None)

        if i % 2 == 0:
            # Line segment
            brg_deg, dist_ft = brg_dist_labels[line_idx]
            line_idx += 1
            rows.append({
                'seq':          i,
                'seg_type':     'L',
                'distance':     dist_ft,
                'radius':       None,
                'sweep_name':   None,
                'sweep':        None,
                'center_name':  None,
                'n_pts':        60,
                'end_name':     end_name,
                'flex':         flex,
                'bearing_flex': 0,
                '_bearing_deg': brg_deg,  # for verification only
            })
        else:
            # Arc segment
            r_ft, sweep_deg = arc_labels[arc_idx]
            direction = fseg_arcs[arc_idx]['direction']
            arc_idx += 1
            center_name = CENTER_NAMES[center_idx]; center_idx += 1
            rows.append({
                'seq':          i,
                'seg_type':     direction,   # 'CW' or 'CCW'
                'distance':     None,
                'radius':       r_ft,
                'sweep_name':   sweep_name_for(sweep_deg),
                'sweep':        math.radians(sweep_deg),
                'center_name':  center_name,
                'n_pts':        20,
                'end_name':     end_name,
                'flex':         flex,
                'bearing_flex': 0,
                '_sweep_deg':   sweep_deg,  # for verification only
            })

    return rows

# ---------------------------------------------------------------------------
# Verify chain closure using the extracted endpoint geometry
# ---------------------------------------------------------------------------

def verify_closure(fseg_lines, fseg_arcs):
    chain = []
    for i in range(12):
        chain.append(fseg_lines[i])
        chain.append(fseg_arcs[i])
    max_gap = 0.0
    for i in range(24):
        end   = chain[i]['end']
        start = chain[(i+1) % 24]['start']
        gap   = math.hypot(end[0]-start[0], end[1]-start[1])
        max_gap = max(max_gap, gap)
    return max_gap

# ---------------------------------------------------------------------------
# Write to DB
# ---------------------------------------------------------------------------

def write_db(rows):
    print(f"Copying {SRC_DB} -> {DST_DB}")
    shutil.copy2(SRC_DB, DST_DB)

    conn = sqlite3.connect(DST_DB)
    conn.execute("DELETE FROM outline_chain")
    conn.execute("DELETE FROM undo_history")
    # Clear stale baseline anchor config (copied from adu.db, references 'F2'/'F12'
    # which don't exist in MarkX's chain — let solver use the closure fallback).
    # Set the correct anchor so the chain walk starts at F24b heading 270° (west).
    # F24b position is extracted from the SVG arc endpoint; brg in radians.
    # Without this the fallback uses F2_EASTING/F2_NORTHING (-18.5, ~-12.2) and
    # start_brg=0 (north) — wrong position and wrong orientation.
    F24b_E   = -11.408783
    F24b_N   = -13.391088
    F24b_brg = 3 * math.pi / 2  # 270° west, the bearing of the first line (F24b→F26)
    anchor_rows = [
        ('outline_anchor',          'F24b'),
        ('outline_pivot',           'F11'),   # arbitrary; used by interactive editor
        ('outline_anchor_E',        str(F24b_E)),
        ('outline_anchor_N',        str(F24b_N)),
        ('outline_anchor_brg',      str(F24b_brg)),
        ('outline_pivot_user_set',  '1'),
    ]
    for key, val in anchor_rows:
        conn.execute("DELETE FROM config WHERE key=?", (key,))
        conn.execute("INSERT INTO config (key, value) VALUES (?,?)", (key, val))

    for r in rows:
        conn.execute("""
            INSERT INTO outline_chain
                (seq, seg_type, distance, radius, sweep_name, sweep,
                 center_name, n_pts, end_name, flex, bearing_flex)
            VALUES (?,?,?,?,?,?,?,?,?,?,?)
        """, (
            r['seq'], r['seg_type'], r['distance'], r['radius'],
            r['sweep_name'], r['sweep'], r['center_name'], r['n_pts'],
            r['end_name'], r['flex'], r['bearing_flex'],
        ))

    conn.commit()
    conn.close()
    print(f"Wrote {len(rows)} rows to outline_chain in {DST_DB}")

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print("Parsing SVG…")
    fseg_lines, fseg_arcs, brg_dist_labels, arc_labels = parse_svg()

    max_gap = verify_closure(fseg_lines, fseg_arcs)
    print(f"Chain closure gap: {max_gap:.6f} ft (expect ~0)")

    rows = build_chain_rows(fseg_lines, fseg_arcs, brg_dist_labels, arc_labels)

    print("\n=== Proposed outline_chain (24 segments) ===")
    print(f"{'seq':3s} {'type':4s} {'dist_ft':9s} {'R_in':7s} {'sweep°':8s} {'dir':4s} {'center':6s} {'end':6s} {'flex'}")
    for r in rows:
        if r['seg_type'] == 'L':
            print(f"{r['seq']:3d}  L     {r['distance']:8.4f}  {'':7s} {'':8s}      {'':6s} {r['end_name']:6s}  {r['flex'] or ''}")
        else:
            sw = math.degrees(r['sweep'])
            print(f"{r['seq']:3d}  {r['seg_type']:4s}  {'':9s} {r['radius']*12:6.2f}\"  {sw:7.1f}°  {'':4s} {r['center_name']:6s} {r['end_name']:6s}  {r['flex'] or ''}")

    print()
    write_db(rows)
    print("Done. Load MarkX variant in the app to verify outline geometry.")

if __name__ == "__main__":
    main()
