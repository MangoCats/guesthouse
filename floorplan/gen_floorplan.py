"""Generate floorplan SVG with 8" wall inset from the outline path.

Runs gen_path_svg.py to obtain outline geometry, then computes an 8" inset.
Points C0-C15 correspond to outline points O0-O15.

Inset rules (CCW path, interior to left):
  - CW arcs  (center outside shape): increase radius by wall thickness
  - CCW arcs (center inside shape):  decrease radius by wall thickness
  - Line segments: offset inward along left normal
"""
import math, sys, os, io

# --- Load outline geometry by executing gen_path_svg.py ---
_parent = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_old_stdout = sys.stdout
sys.stdout = io.StringIO()
_ns = {}
exec(open(os.path.join(_parent, "gen_path_svg.py")).read(), _ns)
sys.stdout = _old_stdout

pts          = _ns["pts"]
outline_segs = _ns["outline_segs"]
LineSeg      = _ns["LineSeg"]
ArcSeg       = _ns["ArcSeg"]
left_norm    = _ns["left_norm"]
off_pt       = _ns["off_pt"]
to_svg       = _ns["to_svg"]
segment_polyline = _ns["segment_polyline"]
path_polygon     = _ns["path_polygon"]
poly_area        = _ns["poly_area"]
W, H         = _ns["W"], _ns["H"]

# Outline arc radii
R_fillet  = _ns["R_fillet"]     # Cf   CCW (inside)
R3i       = _ns["R3i"]          # C3   CW  (outside)
R_f_pox   = _ns["R_f_pox"]     # Cf3  CCW (inside)
R_f_po5   = _ns["R_f_po5"]     # Cf4  CCW (inside)
R1i       = _ns["R1i"]          # C1   CW  (outside)
R_turn3   = _ns["R_turn3"]     # Ct3  CCW (inside)
R_turn2   = _ns["R_turn2"]     # Ct2  CW  (outside)
R_turn1   = _ns["R_turn1"]     # Ct1  CCW (inside)
R_fillet2 = _ns["R_fillet2"]   # Cf2  CCW (inside)

# --- Wall thickness ---
wall_t = 8.0 / 12.0  # 8 inches in feet

# --- Compute inner wall points ---
# At each junction, tangency is preserved by the inset:
#   - Line offset and arc radius adjustment keep center-to-line distance = adjusted radius
#   - Arc-arc external tangency is preserved since +wall_t and -wall_t cancel

def _inner_point(seg_before, seg_after):
    """Inner wall point at junction of seg_before and seg_after."""
    b_line = isinstance(seg_before, LineSeg)
    a_line = isinstance(seg_after, LineSeg)

    if not b_line and not a_line:
        # Arc-Arc: new tangent point on the line between centers
        c1 = pts[seg_before.center]; c2 = pts[seg_after.center]
        r1 = (seg_before.radius + wall_t) if seg_before.direction == "CW" \
             else (seg_before.radius - wall_t)
        dx = c2[0] - c1[0]; dy = c2[1] - c1[1]
        d = math.sqrt(dx*dx + dy*dy)
        return (c1[0] + r1*dx/d, c1[1] + r1*dy/d)

    # Line-Arc or Arc-Line: perpendicular foot from arc center onto offset line
    line_seg = seg_before if b_line else seg_after
    arc_seg  = seg_after  if b_line else seg_before
    center = pts[arc_seg.center]
    S = pts[line_seg.start]; E = pts[line_seg.end]
    D = (E[0]-S[0], E[1]-S[1])
    LN = left_norm(S, E)
    P = off_pt(S, LN, wall_t)
    t = ((center[0]-P[0])*D[0] + (center[1]-P[1])*D[1]) / (D[0]**2 + D[1]**2)
    return (P[0]+t*D[0], P[1]+t*D[1])

for i in range(16):
    seg_b = outline_segs[i]
    seg_a = outline_segs[(i+1) % 16]
    n = int(seg_b.end[1:])   # "O7" -> 7
    pts[f"W{n}"] = _inner_point(seg_b, seg_a)

# --- Inner wall segments (same topology, adjusted radii) ---
inner_segs = [
    LineSeg("W2",  "W1"),
    ArcSeg("W1",  "W0",  "Cf",  R_fillet  - wall_t, "CCW", 20),
    LineSeg("W0",  "W15"),
    ArcSeg("W15", "W14", "C3",  R3i       + wall_t, "CW",  60),
    ArcSeg("W14", "W13", "Cf3", R_f_pox   - wall_t, "CCW", 20),
    LineSeg("W13", "W12"),
    ArcSeg("W12", "W11", "Cf4", R_f_po5   - wall_t, "CCW", 20),
    LineSeg("W11", "W10"),
    ArcSeg("W10", "W9",  "C1",  R1i       + wall_t, "CW",  60),
    LineSeg("W9",  "W8"),
    ArcSeg("W8",  "W7",  "Ct3", R_turn3   - wall_t, "CCW", 20),
    LineSeg("W7",  "W6"),
    ArcSeg("W6",  "W5",  "Ct2", R_turn2   + wall_t, "CW",  20),
    ArcSeg("W5",  "W4",  "Ct1", R_turn1   - wall_t, "CCW", 20),
    LineSeg("W4",  "W3"),
    ArcSeg("W3",  "W2",  "Cf2", R_fillet2 - wall_t, "CCW", 20),
]

outer_poly = path_polygon(outline_segs, pts)
inner_poly = path_polygon(inner_segs, pts)
outer_area = poly_area(outer_poly)
inner_area = poly_area(inner_poly)

# --- SVG ---
out = []
out.append(f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}">')
out.append(f'<rect width="{W}" height="{H}" fill="white"/>')
out.append('<defs>')
out.append('  <marker id="ah" markerWidth="8" markerHeight="6" refX="8" refY="3" orient="auto">'
           '<polygon points="0 0, 8 3, 0 6" fill="#333"/></marker>')
out.append('</defs>')
out.append(f'<text x="{W/2}" y="30" text-anchor="middle" font-family="Arial" font-size="14"'
           f' font-weight="bold">Floorplan &#8212; 8&#8243; Walls</text>')

# Wall fill: outer gray, inner white cutout
outer_svg = " ".join(f"{to_svg(*p)[0]:.1f},{to_svg(*p)[1]:.1f}" for p in outer_poly)
inner_rev = list(reversed(inner_poly))
inner_svg = " ".join(f"{to_svg(*p)[0]:.1f},{to_svg(*p)[1]:.1f}" for p in inner_rev)
out.append(f'<polygon points="{outer_svg}" fill="rgba(160,160,160,0.35)" stroke="none"/>')
out.append(f'<polygon points="{inner_svg}" fill="white" stroke="none"/>')

# Outer outline strokes
for seg in outline_segs:
    if isinstance(seg, LineSeg):
        sx1, sy1 = to_svg(*pts[seg.start]); sx2, sy2 = to_svg(*pts[seg.end])
        out.append(f'<line x1="{sx1:.1f}" y1="{sy1:.1f}" x2="{sx2:.1f}" y2="{sy2:.1f}"'
                   f' stroke="#333" stroke-width="1.5"/>')
    else:
        poly = segment_polyline(seg, pts)
        svg_p = " ".join(f"{to_svg(e,n)[0]:.1f},{to_svg(e,n)[1]:.1f}" for e,n in poly)
        out.append(f'<polyline points="{svg_p}" fill="none" stroke="#333"'
                   f' stroke-width="1.5" stroke-linecap="round"/>')

# Inner outline strokes
for seg in inner_segs:
    if isinstance(seg, LineSeg):
        sx1, sy1 = to_svg(*pts[seg.start]); sx2, sy2 = to_svg(*pts[seg.end])
        out.append(f'<line x1="{sx1:.1f}" y1="{sy1:.1f}" x2="{sx2:.1f}" y2="{sy2:.1f}"'
                   f' stroke="#666" stroke-width="1.0"/>')
    else:
        poly = segment_polyline(seg, pts)
        svg_p = " ".join(f"{to_svg(e,n)[0]:.1f},{to_svg(e,n)[1]:.1f}" for e,n in poly)
        out.append(f'<polyline points="{svg_p}" fill="none" stroke="#666"'
                   f' stroke-width="1.0" stroke-linecap="round"/>')

# C-point labels (outer vertices displayed as C0-C15)
vs_map = _ns["outline_cfg"].vertex_styles
for i in range(16):
    o_name = f"O{i}"
    sx, sy = to_svg(*pts[o_name])
    out.append(f'<circle cx="{sx:.1f}" cy="{sy:.1f}" r="2.5" fill="#333"/>')
    if o_name in vs_map:
        vs = vs_map[o_name]
        out.append(f'<text x="{sx+vs.dx:.1f}" y="{sy+vs.dy:.1f}" text-anchor="{vs.anchor}"'
                   f' font-family="Arial" font-size="9" font-weight="bold"'
                   f' fill="#333">C{i}</text>')

# Interior area label
cx_o = sum(pts[f"O{i}"][0] for i in range(16)) / 16
cy_o = sum(pts[f"O{i}"][1] for i in range(16)) / 16
sx, sy = to_svg(cx_o, cy_o)
out.append(f'<text x="{sx:.1f}" y="{sy:.1f}" text-anchor="middle" font-family="Arial"'
           f' font-size="11" fill="#333" font-weight="bold">{inner_area:.2f} sq ft</text>')
out.append(f'<text x="{sx:.1f}" y="{sy+13:.1f}" text-anchor="middle" font-family="Arial"'
           f' font-size="8" fill="#666">(Interior area)</text>')

# North arrow
out.append('<line x1="742" y1="560" x2="742" y2="524" stroke="#333" stroke-width="2"'
           ' marker-end="url(#ah)"/>')
out.append('<text x="742" y="518" text-anchor="middle" font-family="Arial"'
           ' font-size="13" font-weight="bold">N</text>')

# Footer
out.append(f'<text x="{W/2}" y="{H-2}" text-anchor="middle" font-family="Arial" font-size="7.5"'
           f' fill="#999">Wall thickness: 8&#8243; &#8226;'
           f' Outer: {outer_area:.2f} sq ft &#8226; Interior: {inner_area:.2f} sq ft</text>')
out.append('</svg>')

svg_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "floorplan.svg")
with open(svg_path, "w") as f:
    f.write("\n".join(out))

print(f"Floorplan written to {svg_path}")
print(f"Outer area:    {outer_area:.2f} sq ft")
print(f"Interior area: {inner_area:.2f} sq ft")
print(f"Wall area:     {outer_area - inner_area:.2f} sq ft")
print()
for i in range(16):
    o = pts[f"O{i}"]; w = pts[f"W{i}"]
    print(f"  C{i:<2d} ({o[0]:8.4f}, {o[1]:8.4f})  ->  inner ({w[0]:8.4f}, {w[1]:8.4f})")
