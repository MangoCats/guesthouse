"""DB-driven floorplan SVG renderer.

Renders floorplan SVGs entirely from compute_geometry() output (DB state).
No hardcoded element geometry — all interior walls, openings, variant items,
dimensions, and room labels come from the database via FormulaEvaluator.

Outer wall shell rendering (double-shell strips, U-turns, section outlines)
is delegated to the existing _render_walls() in gen_floorplan.py which
already uses GeneratorData (non-hardcoded) for shell geometry.

This module implements the ARCHITECTURAL MANDATE: no generated floorplan SVG
shall ever derive element geometry from hardcoded modules.
"""
import math
import datetime

from shared.svg import make_svg_transform, W, H
from shared.geometry import fmt_dist, seg_vecs, offset_pt, poly_area

# --- Style constants (match gen_floorplan.py) ---
APPL_FILL = 'rgba(100,150,200,0.2)'
APPL_STROKE = '#4682B4'
APPL_SW = '0.8'
WALL_FILL = 'rgba(160,160,160,0.35)'
WALL_STROKE = '#666'
WALL_SW = '1.0'
OPENING_FILL = 'rgb(220,235,255)'
OPENING_STROKE = '#4682B4'
JAMB_COLOR = 'darkred'
DIM_COLOR = '#999'

# Jamb width in feet (1" = 1/12')
JAMB_WIDTH = 1.0 / 12.0


# ---------------------------------------------------------------------------
# SVG helpers (reused from gen_floorplan.py patterns)
# ---------------------------------------------------------------------------

def _svg_angle(along):
    """SVG rotation angle (degrees, CW-positive) for a direction vector."""
    return -math.degrees(math.atan2(along[1], along[0]))


def _poly_svg(pts, to_svg):
    """Convert list of (E,N) points to SVG polygon points string."""
    return " ".join(f"{to_svg(p[0], p[1])[0]:.1f},{to_svg(p[0], p[1])[1]:.1f}"
                    for p in pts)


def _wall_poly(out, poly, to_svg, stroke=True):
    """Render wall polygon with standard gray fill and optional clipped stroke."""
    svg = _poly_svg(poly, to_svg)
    if stroke:
        cid = f"wc{len(out)}"
        out.append(f'<defs><clipPath id="{cid}"><polygon points="{svg}"/></clipPath></defs>')
        out.append(f'<polygon points="{svg}" fill="{WALL_FILL}"'
                   f' stroke="{WALL_STROKE}" stroke-width="1.6" clip-path="url(#{cid})"/>')
    else:
        out.append(f'<polygon points="{svg}" fill="{WALL_FILL}" stroke="none"/>')


def _wall_stroke_line(out, p_start, p_end, half_sw, to_svg):
    """Render a wall face stroke line inset by half_sw."""
    _al, _out = seg_vecs(p_start, p_end)
    _inw = (-_out[0], -_out[1])
    _p1 = offset_pt(p_start, half_sw, _inw)
    _p2 = offset_pt(p_end, half_sw, _inw)
    sx1, sy1 = to_svg(*_p1)
    sx2, sy2 = to_svg(*_p2)
    out.append(f'<line x1="{sx1:.1f}" y1="{sy1:.1f}" x2="{sx2:.1f}" y2="{sy2:.1f}"'
               f' stroke="{WALL_STROKE}" stroke-width="{WALL_SW}"/>')


def _jamb_poly(out, j1, j2, along_dir, to_svg):
    """Render a jamb polygon from edge points j1, j2 + JAMB_WIDTH along along_dir."""
    jn = (along_dir[0] * JAMB_WIDTH, along_dir[1] * JAMB_WIDTH)
    j_poly = [j1, j2, (j2[0] + jn[0], j2[1] + jn[1]),
              (j1[0] + jn[0], j1[1] + jn[1])]
    jp = _poly_svg(j_poly, to_svg)
    out.append(f'<polygon points="{jp}" fill="{JAMB_COLOR}" stroke="none"/>')


def _appl_poly(out, corners, to_svg, label=None, href=None,
               fill=APPL_FILL, stroke=APPL_STROKE, sw=APPL_SW,
               font_size="7", dash=False, text_rot=None):
    """Render appliance/furniture polygon with optional label and link."""
    if href:
        out.append(f'<a href="{href}" target="_blank">')
    pts_svg = _poly_svg(corners, to_svg)
    attrs = ' stroke-dasharray="4,3"' if dash else ''
    out.append(f'<polygon points="{pts_svg}" fill="{fill}" stroke="{stroke}" '
               f'stroke-width="{sw}"{attrs}/>')
    if label:
        cx = sum(p[0] for p in corners) / len(corners)
        cy = sum(p[1] for p in corners) / len(corners)
        scx, scy = to_svg(cx, cy)
        rot = f' transform="rotate({text_rot:.1f},{scx:.1f},{scy+3:.1f})"' if text_rot is not None else ''
        out.append(f'<text x="{scx:.1f}" y="{scy+3:.1f}" text-anchor="middle" '
                   f'font-family="Arial" font-size="{font_size}" fill="{stroke}"{rot}>{label}</text>')
    if href:
        out.append('</a>')


def _appl_circle(out, center, radius, to_svg, label=None, href=None,
                 fill=APPL_FILL, stroke=APPL_STROKE, sw=APPL_SW,
                 font_size="7"):
    """Render circular appliance/furniture element."""
    if href:
        out.append(f'<a href="{href}" target="_blank">')
    cx, cy = to_svg(center[0], center[1])
    # Scale radius from world to SVG
    sx0, _ = to_svg(0, 0)
    sx1, _ = to_svg(radius, 0)
    r_svg = abs(sx1 - sx0)
    out.append(f'<circle cx="{cx:.1f}" cy="{cy:.1f}" r="{r_svg:.1f}" '
               f'fill="{fill}" stroke="{stroke}" stroke-width="{sw}"/>')
    if label:
        out.append(f'<text x="{cx:.1f}" y="{cy + 3:.1f}" text-anchor="middle" '
                   f'font-family="Arial" font-size="{font_size}" fill="{stroke}">{label}</text>')
    if href:
        out.append('</a>')


# ---------------------------------------------------------------------------
# Interior wall rendering from DB data
# ---------------------------------------------------------------------------

def render_interior_walls_db(out, geom, to_svg):
    """Render interior walls from compute_geometry() output.

    Each IW is rendered as a filled polygon with clipped strokes on all faces.
    Rough openings split their host wall into segments with jamb blocks.
    """
    interior_walls = geom.get("interior_walls", {})
    rough_openings = geom.get("rough_openings", [])

    # Build lookup: wall_name -> list of rough openings
    ro_by_wall = {}
    for ro in rough_openings:
        wn = ro.get("wall_name", "")
        ro_by_wall.setdefault(wn, []).append(ro)

    # SVG scale for half stroke width
    sx0, _ = to_svg(0, 0)
    sx1, _ = to_svg(1, 0)
    svg_per_ft = abs(sx1 - sx0)
    half_sw = 0.5 / svg_per_ft if svg_per_ft > 0 else 0

    for iw_name, iw_data in sorted(interior_walls.items()):
        poly = iw_data["poly"]
        if len(poly) < 4:
            continue
        # Convert poly points to tuples
        poly = [(p[0], p[1]) for p in poly]

        wall_ros = ro_by_wall.get(iw_name, [])
        if not wall_ros:
            # Simple wall — no openings
            _wall_poly(out, poly, to_svg, stroke=False)
            # Stroke all 4 faces
            for i in range(len(poly)):
                _wall_stroke_line(out, poly[i], poly[(i + 1) % len(poly)], half_sw, to_svg)
        else:
            # Wall with rough openings — split and render segments
            _render_split_wall(out, iw_name, poly, wall_ros, half_sw, to_svg)


def _render_split_wall(out, iw_name, poly, rough_openings, half_sw, to_svg):
    """Render a wall split by rough openings.

    For each rough opening, determines which face it intersects and splits
    the wall polygon accordingly, rendering jamb blocks at the split points.
    """
    if len(poly) < 4:
        return

    # For each RO, find where it intersects this wall
    # RO poly has 4 vertices; its bbox tells us the cut location
    for ro in rough_openings:
        ro_poly = ro.get("poly", [])
        if len(ro_poly) < 4:
            continue
        ro_poly = [(p[0], p[1]) for p in ro_poly]

        # Determine orientation from RO properties
        orientation = ro.get("orientation", "H")

        # Find which wall edge the RO intersects by projecting RO center
        # onto wall edges
        ro_cx = sum(p[0] for p in ro_poly) / len(ro_poly)
        ro_cy = sum(p[1] for p in ro_poly) / len(ro_poly)

        # Simple approach: render the wall polygon minus the RO area
        # by splitting into sub-polygons on each side of the RO
        _render_wall_with_opening(out, poly, ro_poly, orientation, half_sw, to_svg)

    # If multiple ROs on same wall, they've each been rendered
    # For now, handle single RO case (most common)


def _render_wall_with_opening(out, wall_poly, ro_poly, orientation, half_sw, to_svg):
    """Render wall polygon split by a single rough opening.

    Determines the split axis from orientation (H=horizontal wall split
    vertically by RO, V=vertical wall split horizontally by RO).
    """
    wp = wall_poly
    rp = ro_poly

    # Get wall and RO bounding boxes
    w_es = [p[0] for p in wp]
    w_ns = [p[1] for p in wp]
    r_es = [p[0] for p in rp]
    r_ns = [p[1] for p in rp]

    if orientation == "H":
        # Horizontal wall (runs E-W), RO splits it vertically
        # Wall vertices: SW=0, SE=1, NE=2, NW=3 (typical)
        # Split into west piece and east piece
        # West piece: wall vertices west of RO
        # East piece: wall vertices east of RO

        # Find RO west and east edges on wall's south and north faces
        ro_w = min(r_es)
        ro_e = max(r_es)

        # Build west sub-wall: from wall west end to RO west edge
        # Build east sub-wall: from RO east edge to wall east end
        # This works for any wall vertex ordering if we know the split coords

        # Use RO poly vertices directly as the split boundary
        # RO poly: [SW, SE, NE, NW] typically
        # Sort RO vertices into west pair and east pair
        ro_sorted_e = sorted(rp, key=lambda p: p[0])
        ro_west_pair = sorted(ro_sorted_e[:2], key=lambda p: p[1])  # [south, north]
        ro_east_pair = sorted(ro_sorted_e[2:], key=lambda p: p[1])  # [south, north]

        # Sort wall vertices into west pair and east pair
        w_sorted_e = sorted(wp, key=lambda p: p[0])
        w_west_pair = sorted(w_sorted_e[:2], key=lambda p: p[1])  # [south, north]
        w_east_pair = sorted(w_sorted_e[2:], key=lambda p: p[1])  # [south, north]

        # West sub-wall: w_west_south, ro_west_south, ro_west_north, w_west_north
        west_poly = [w_west_pair[0], ro_west_pair[0], ro_west_pair[1], w_west_pair[1]]
        # East sub-wall: ro_east_south, w_east_south, w_east_north, ro_east_north
        east_poly = [ro_east_pair[0], w_east_pair[0], w_east_pair[1], ro_east_pair[1]]

        for sub in [west_poly, east_poly]:
            if _poly_has_area(sub):
                _wall_poly(out, sub, to_svg, stroke=False)
                for i in range(len(sub)):
                    _wall_stroke_line(out, sub[i], sub[(i + 1) % len(sub)], half_sw, to_svg)

        # Jamb blocks at RO edges
        # Along direction = wall long axis (E-W)
        wall_al = _unit_vec(w_west_pair[0], w_east_pair[0])
        neg_al = (-wall_al[0], -wall_al[1])
        _jamb_poly(out, ro_west_pair[1], ro_west_pair[0], wall_al, to_svg)
        _jamb_poly(out, ro_east_pair[0], ro_east_pair[1], neg_al, to_svg)

    else:
        # Vertical wall (runs N-S), RO splits it horizontally
        ro_s = min(r_ns)
        ro_n = max(r_ns)

        ro_sorted_n = sorted(rp, key=lambda p: p[1])
        ro_south_pair = sorted(ro_sorted_n[:2], key=lambda p: p[0])
        ro_north_pair = sorted(ro_sorted_n[2:], key=lambda p: p[0])

        w_sorted_n = sorted(wp, key=lambda p: p[1])
        w_south_pair = sorted(w_sorted_n[:2], key=lambda p: p[0])
        w_north_pair = sorted(w_sorted_n[2:], key=lambda p: p[0])

        south_poly = [w_south_pair[0], w_south_pair[1], ro_south_pair[1], ro_south_pair[0]]
        north_poly = [ro_north_pair[0], ro_north_pair[1], w_north_pair[1], w_north_pair[0]]

        for sub in [south_poly, north_poly]:
            if _poly_has_area(sub):
                _wall_poly(out, sub, to_svg, stroke=False)
                for i in range(len(sub)):
                    _wall_stroke_line(out, sub[i], sub[(i + 1) % len(sub)], half_sw, to_svg)

        wall_al = _unit_vec(w_south_pair[0], w_north_pair[0])
        neg_al = (-wall_al[0], -wall_al[1])
        _jamb_poly(out, ro_south_pair[0], ro_south_pair[1], wall_al, to_svg)
        _jamb_poly(out, ro_north_pair[1], ro_north_pair[0], neg_al, to_svg)


def _poly_has_area(poly):
    """Check if polygon has non-trivial area."""
    if len(poly) < 3:
        return False
    # Quick cross product check
    area = 0
    for i in range(len(poly)):
        j = (i + 1) % len(poly)
        area += poly[i][0] * poly[j][1] - poly[j][0] * poly[i][1]
    return abs(area) > 1e-10


def _unit_vec(p1, p2):
    """Unit direction vector from p1 to p2."""
    dx = p2[0] - p1[0]
    dy = p2[1] - p1[1]
    d = math.sqrt(dx * dx + dy * dy)
    if d < 1e-12:
        return (1, 0)
    return (dx / d, dy / d)


# ---------------------------------------------------------------------------
# Variant items rendering from DB data
# ---------------------------------------------------------------------------

def render_variant_items_db(out, geom, to_svg):
    """Render all variant items (furniture, appliances, fixtures) from DB data.

    Each item is rendered as a polygon or circle with label and optional
    product URL link.  Stacked items are rendered last (SVG paint order).
    """
    variant_items = geom.get("variant_items", {})

    # Sort: non-stacked first, stacked last
    entries = sorted(variant_items.items(),
                     key=lambda kv: (1 if kv[1].get("stacked") else 0, kv[0]))

    for name, item in entries:
        shape = item.get("shape", "rect")
        label = item.get("label", name.upper())
        href = item.get("product_url")

        if shape == "circle":
            center = item.get("center")
            radius = item.get("radius")
            if center and radius:
                _appl_circle(out, center, radius, to_svg,
                             label=label, href=href, font_size="7")
        else:
            poly = item.get("poly", [])
            if not poly:
                continue
            poly = [(p[0], p[1]) for p in poly]
            # Compute text rotation from polygon edge direction
            text_rot = None
            if len(poly) >= 2:
                # Use edge 0->1 direction for rotation
                dx = poly[1][0] - poly[0][0]
                dy = poly[1][1] - poly[0][1]
                if abs(dx) > 1e-6 or abs(dy) > 1e-6:
                    angle = _svg_angle((dx, dy))
                    # Only rotate if significantly non-axis-aligned
                    if abs(angle) > 5 and abs(abs(angle) - 90) > 5 and abs(abs(angle) - 180) > 5:
                        text_rot = angle

            _appl_poly(out, poly, to_svg, label=label, href=href,
                       font_size="7", text_rot=text_rot)


# ---------------------------------------------------------------------------
# Door arcs rendering from DB data
# ---------------------------------------------------------------------------

def render_door_arcs_db(out, geom, to_svg):
    """Render door swing arcs from compute_geometry() output.

    Handles both structural doors (door_arcs) and appliance doors
    (appliance_doors).
    """
    # Structural door arcs (O3, O6, RO1-RO7)
    for da in geom.get("door_arcs", []):
        for leaf in da.get("leaves", []):
            hinge = leaf["hinge"]
            tip = leaf["tip"]
            arc_pts = leaf["arc_pts"]

            hx, hy = to_svg(hinge[0], hinge[1])
            tx, ty = to_svg(tip[0], tip[1])

            # Hinge-to-tip line
            out.append(f'<line x1="{hx:.1f}" y1="{hy:.1f}" '
                       f'x2="{tx:.1f}" y2="{ty:.1f}" '
                       f'stroke="{JAMB_COLOR}" stroke-width="1.0"/>')
            # Swing arc
            if arc_pts:
                pts_str = " ".join(
                    f"{to_svg(p[0], p[1])[0]:.1f},{to_svg(p[0], p[1])[1]:.1f}"
                    for p in arc_pts)
                out.append(f'<polyline points="{pts_str}" fill="none" '
                           f'stroke="{JAMB_COLOR}" stroke-width="0.5"/>')

    # Appliance door arcs (fridge, washer, dryer, microwave)
    for ad in geom.get("appliance_doors", []):
        hinge = ad["hinge"]
        tip = ad["tip"]
        arc_pts = ad["arc_pts"]

        hx, hy = to_svg(hinge[0], hinge[1])
        tx, ty = to_svg(tip[0], tip[1])

        out.append(f'<line x1="{hx:.1f}" y1="{hy:.1f}" '
                   f'x2="{tx:.1f}" y2="{ty:.1f}" '
                   f'stroke="{APPL_STROKE}" stroke-width="0.8"/>')
        if arc_pts:
            pts_str = " ".join(
                f"{to_svg(p[0], p[1])[0]:.1f},{to_svg(p[0], p[1])[1]:.1f}"
                for p in arc_pts)
            out.append(f'<polyline points="{pts_str}" fill="none" '
                       f'stroke="{APPL_STROKE}" stroke-width="0.5"/>')


# ---------------------------------------------------------------------------
# Casement window rendering from DB data
# ---------------------------------------------------------------------------

def render_casement_windows_db(out, geom, to_svg, shell_thickness):
    """Render casement window swing arcs for openings with opening_type=casement.

    Casement windows swing outward at 45 degrees (vs doors at 90 degrees).
    """
    for oo in geom.get("outer_openings", []):
        if oo.get("opening_type") != "casement":
            continue
        poly = oo.get("poly", [])
        if len(poly) < 4:
            continue
        poly = [(p[0], p[1]) for p in poly]

        # Inward direction (outer face midpoint toward inner face midpoint)
        omid = ((poly[0][0] + poly[1][0]) / 2, (poly[0][1] + poly[1][1]) / 2)
        imid = ((poly[2][0] + poly[3][0]) / 2, (poly[2][1] + poly[3][1]) / 2)
        iE = imid[0] - omid[0]
        iN = imid[1] - omid[1]
        ilen = math.sqrt(iE ** 2 + iN ** 2)
        if ilen < 1e-12:
            continue
        idir = (iE / ilen, iN / ilen)

        # Opening width direction (along outer face)
        dE = poly[1][0] - poly[0][0]
        dN = poly[1][1] - poly[0][1]
        wlen = math.sqrt(dE ** 2 + dN ** 2)
        if wlen < 1e-12:
            continue

        # Hinge at one end, offset inward by shell thickness
        # Determine hinge end based on opening name convention
        name = oo.get("name", "")
        # O9 hinges at poly[1] side, O8 and O10 hinge at poly[0] side
        if name == "O9":
            hinge_idx, close_idx = 1, 0
        else:
            hinge_idx, close_idx = 0, 1

        hinge = (poly[hinge_idx][0] + shell_thickness * idir[0],
                 poly[hinge_idx][1] + shell_thickness * idir[1])
        close = (poly[close_idx][0] + shell_thickness * idir[0],
                 poly[close_idx][1] + shell_thickness * idir[1])

        cdE = close[0] - hinge[0]
        cdN = close[1] - hinge[1]
        clen = math.sqrt(cdE ** 2 + cdN ** 2)
        if clen < 1e-12:
            continue
        cdir = (cdE / clen, cdN / clen)

        # Swing outward (away from interior)
        cross = cdir[0] * idir[1] - cdir[1] * idir[0]
        rsign = -1 if cross > 0 else 1
        open_angle = rsign * math.pi / 4  # 45 degrees

        # Open tip position
        cos_a = math.cos(open_angle)
        sin_a = math.sin(open_angle)
        odir = (cdir[0] * cos_a - cdir[1] * sin_a,
                cdir[0] * sin_a + cdir[1] * cos_a)
        tip = (hinge[0] + clen * odir[0], hinge[1] + clen * odir[1])

        hx, hy = to_svg(*hinge)
        tx, ty = to_svg(*tip)

        # Product URL from opening properties
        href = oo.get("product_url", "")

        line_el = (f'<line x1="{hx:.1f}" y1="{hy:.1f}" '
                   f'x2="{tx:.1f}" y2="{ty:.1f}" '
                   f'stroke="{OPENING_STROKE}" stroke-width="1.0"/>')

        # Arc from open to closed
        n_arc = 10
        arc_pts_strs = []
        for i in range(n_arc + 1):
            a = open_angle * (1 - i / n_arc)
            ca = math.cos(a)
            sa = math.sin(a)
            d = (cdir[0] * ca - cdir[1] * sa,
                 cdir[0] * sa + cdir[1] * ca)
            pt = (hinge[0] + clen * d[0], hinge[1] + clen * d[1])
            sx, sy = to_svg(*pt)
            arc_pts_strs.append(f"{sx:.1f},{sy:.1f}")
        arc_el = (f'<polyline points="{" ".join(arc_pts_strs)}" fill="none" '
                  f'stroke="{OPENING_STROKE}" stroke-width="0.5"/>')

        if href:
            out.append(f'<a href="{href}">{line_el}{arc_el}</a>')
        else:
            out.append(f'{line_el}{arc_el}')


# ---------------------------------------------------------------------------
# Clearance zones rendering from DB data
# ---------------------------------------------------------------------------

def render_clearance_zones_db(out, geom, to_svg):
    """Render clearance zones from compute_geometry() output."""
    for cz in geom.get("clearance_zones", []):
        poly = cz.get("poly", [])
        if not poly:
            continue
        poly = [(p[0], p[1]) for p in poly]
        pts_svg = _poly_svg(poly, to_svg)
        out.append(f'<polygon points="{pts_svg}" fill="none" '
                   f'stroke="{APPL_STROKE}" stroke-width="0.5" '
                   f'stroke-dasharray="4,3"/>')


# ---------------------------------------------------------------------------
# Dimension lines rendering from DB data
# ---------------------------------------------------------------------------

def render_dimensions_db(out, geom, to_svg):
    """Render dimension lines from compute_geometry() output."""
    for dim in geom.get("user_dimensions", []):
        props = dim.get("properties", {})
        start = props.get("start")
        end = props.get("end")
        if not start or not end:
            continue

        p1 = (start[0], start[1])
        p2 = (end[0], end[1])

        # Compute distance
        dx = p2[0] - p1[0]
        dy = p2[1] - p1[1]
        dist = math.sqrt(dx * dx + dy * dy)
        if dist < 1e-6:
            continue

        # Format label
        label_prefix = props.get("label_prefix", "")
        label = label_prefix + fmt_dist(dist)

        # Optional label positioning point
        label_pt = None
        offset_data = props.get("offset")
        if offset_data and "offset" in offset_data:
            label_pt = tuple(offset_data["offset"])

        # Use the same rendering as gen_floorplan's _rotated_dim
        _rotated_dim(out, p1, p2, label, to_svg, label_pt=label_pt)


def _rotated_dim(out, p1, p2, label, to_svg, label_pt=None):
    """Rotated dimension line with tick marks and label."""
    sx1, sy1 = to_svg(*p1)
    sx2, sy2 = to_svg(*p2)
    sdx = sx2 - sx1
    sdy = sy2 - sy1
    slen = math.sqrt(sdx ** 2 + sdy ** 2)
    if slen < 1e-6:
        return
    px = -sdy / slen
    py = sdx / slen
    tk = 4

    out.append(f'<line x1="{sx1:.1f}" y1="{sy1:.1f}" x2="{sx2:.1f}" y2="{sy2:.1f}" '
               f'stroke="{DIM_COLOR}" stroke-width="0.8"/>')
    for sx, sy in [(sx1, sy1), (sx2, sy2)]:
        out.append(f'<line x1="{sx - tk * px:.1f}" y1="{sy - tk * py:.1f}" '
                   f'x2="{sx + tk * px:.1f}" y2="{sy + tk * py:.1f}" '
                   f'stroke="{DIM_COLOR}" stroke-width="0.8"/>')

    if label_pt is not None:
        lsx, lsy = to_svg(*label_pt)
        ux, uy = sdx / slen, sdy / slen
        t = (lsx - sx1) * ux + (lsy - sy1) * uy
        lmx = sx1 + t * ux
        lmy = sy1 + t * uy
    else:
        lmx = (sx1 + sx2) / 2
        lmy = (sy1 + sy2) / 2

    ang = math.degrees(math.atan2(sdy, sdx))
    if ang >= 90:
        ang -= 180
    elif ang < -90:
        ang += 180

    ang_rad = math.radians(ang)
    lx = lmx + 3 * math.sin(ang_rad)
    ly = lmy - 3 * math.cos(ang_rad)
    out.append(f'<text x="{lx:.1f}" y="{ly:.1f}" text-anchor="middle" font-family="Arial" '
               f'font-size="8" fill="{DIM_COLOR}" transform="rotate({ang:.1f},{lx:.1f},{ly:.1f})">'
               f'{label}</text>')


# ---------------------------------------------------------------------------
# Room labels rendering from DB data
# ---------------------------------------------------------------------------

def render_room_labels_db(out, geom, to_svg):
    """Render room labels from compute_geometry() output."""
    for rl in geom.get("room_labels", []):
        pos = rl.get("pos")
        if not pos:
            continue
        name = rl.get("name", "")

        # Look up font size from label_elements
        font_size = 8  # default
        for le in geom.get("label_elements", []):
            if le.get("name") == name:
                props = le.get("properties", {})
                fs = props.get("font_size")
                if fs and fs > 0:
                    # Convert from feet to SVG px
                    sx0, _ = to_svg(0, 0)
                    sx1, _ = to_svg(fs, 0)
                    font_size = abs(sx1 - sx0)
                break

        sx, sy = to_svg(pos[0], pos[1])
        # Display text from label element or room name
        display_text = name
        for le in geom.get("label_elements", []):
            if le.get("name") == name:
                props = le.get("properties", {})
                display_text = props.get("text", name)
                break

        out.append(f'<text x="{sx:.1f}" y="{sy:.1f}" text-anchor="middle" '
                   f'font-family="Arial" font-size="{font_size:.1f}" '
                   f'fill="#666">{display_text}</text>')


# ---------------------------------------------------------------------------
# Main DB-driven SVG renderer
# ---------------------------------------------------------------------------

def render_floorplan_svg_db(geom, data, room_title="Parent Suite",
                            shell_thickness=None):
    """Render complete floorplan SVG from DB-driven geometry.

    Args:
        geom: compute_geometry() output dict (all element geometry from DB)
        data: FloorplanData (page layout, outer wall shell geometry)
        room_title: title string for the SVG
        shell_thickness: shell thickness in feet (for casement window hinges)

    Returns:
        SVG string
    """
    from floorplan.gen_floorplan import (
        _render_walls, _render_title_block, _render_sf_extras,
        _render_plumbing_path, _render_supplies_table,
        compute_iw_area, git_describe,
    )
    from floorplan.constants import SHELL_THICKNESS

    if shell_thickness is None:
        shell_thickness = SHELL_THICKNESS

    pts = data.pts
    to_svg = data.to_svg
    layout = data.layout  # still needed for outer wall shell rendering

    vb_x, vb_y, vb_w, vb_h = data.vb_x, data.vb_y, data.vb_w, data.vb_h
    page_w, page_h = W, H

    out = []
    out.append(f'<svg xmlns="http://www.w3.org/2000/svg" width="{page_w:.0f}" height="{page_h:.0f}"'
               f' viewBox="{vb_x:.2f} {vb_y:.2f} {vb_w:.2f} {vb_h:.2f}">')
    out.append(f'<rect x="{vb_x:.2f}" y="{vb_y:.2f}" width="{vb_w:.2f}" height="{vb_h:.2f}" fill="white"/>')
    out.append('<defs>')
    out.append('  <marker id="ah" markerWidth="8" markerHeight="6" refX="8" refY="3" orient="auto">'
               '<polygon points="0 0, 8 3, 0 6" fill="#333"/></marker>')
    out.append('</defs>')
    out.append(f'<text x="{data.title_x:.1f}" y="{data.title_y:.1f}" text-anchor="middle" font-family="Arial" font-size="14"'
               f' font-weight="bold">{room_title}</text>')

    # Outer wall shells (uses GeneratorData — already non-hardcoded)
    # skip_interior_walls=True: interior walls are rendered by render_interior_walls_db
    variant = geom.get("variant", "standard")
    bare = variant in ("bare", "sf")
    _render_walls(out, data, layout, bare=bare, skip_interior_walls=True)

    # DB-driven interior walls (replaces hardcoded IW rendering in _render_walls)
    render_interior_walls_db(out, geom, to_svg)

    # DB-driven variant items (replaces _render_appliances, _render_kitchen, _render_furniture)
    if not bare:
        render_variant_items_db(out, geom, to_svg)

    # Door arcs (replaces hardcoded _render_openings door swing code)
    render_door_arcs_db(out, geom, to_svg)

    # Casement windows (replaces hardcoded casement rendering)
    render_casement_windows_db(out, geom, to_svg, shell_thickness)

    # Clearance zones
    if not bare:
        render_clearance_zones_db(out, geom, to_svg)

    # Dimensions
    out.append('<g opacity="0.5">')
    render_dimensions_db(out, geom, to_svg)
    out.append('</g>')

    # Room labels
    render_room_labels_db(out, geom, to_svg)

    # Title block
    # Compute inner area from DB geometry
    inner_poly = geom.get("inner_poly", [])
    inner_area = poly_area([(p[0], p[1]) for p in inner_poly]) if inner_poly else data.inner_area
    # Subtract IW area from DB interior walls
    iw_area = 0
    for iw_name, iw_data in geom.get("interior_walls", {}).items():
        poly = iw_data.get("poly", [])
        if len(poly) >= 3:
            iw_area += abs(poly_area([(p[0], p[1]) for p in poly]))
    net_area = inner_area - iw_area
    _render_title_block(out, data, net_area)

    out.append('</svg>')
    return "\n".join(out)
