"""Generate 4 architectural elevation PDFs for building department submission.

Each PDF is US Letter landscape (11" x 8.5") showing an orthographic
elevation of the 2:12 roofed structure from one cardinal direction.

Outputs:
  elevations/elevation_south.pdf
  elevations/elevation_north.pdf
  elevations/elevation_east.pdf
  elevations/elevation_west.pdf

Convention: "South Elevation" shows the south-facing wall as seen by a
viewer standing to the south looking north, per standard architectural
practice.
"""
import math
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

import fitz  # PyMuPDF

from shared.types import ArcSeg
from floorplan.geometry import compute_outline_geometry
from floorplan.constants import WALL_OUTER, ROOF_OVERHANG
from floorplan.roof import compute_roof_geometry, roof_polyline
from shared.geometry import compute_inner_walls, path_polygon

_DIR = os.path.dirname(os.path.abspath(__file__))

# ── Height constants (matching scad/gen_2in12.py) ──────────────────────
LOWER_HEIGHT_FT = 20.0 / 12.0       # 1'8"  — bottom of window openings
WALL_HEIGHT_FT = 80.0 / 12.0        # 6'8"  — top of middle walls / head of openings
ROOF_THICK_FT = 18.0 / 12.0         # 1'6"  — roof slab thickness
ROOF_SLOPE = 2.0 / 12.0             # 2:12 slope (ft rise per ft run, northward)
ROOF_REF_ELEV_FT = 7.5              # 7'6" roof underside elevation at F18-F1

# ── Page layout (points, 72 pt/inch) ───────────────────────────────────
PAGE_W, PAGE_H = 792, 612           # 11" x 8.5" landscape
M_LEFT = 108                        # 1.5" left margin (vertical dims)
M_RIGHT = 36                        # 0.5" right
M_TOP = 54                          # 0.75" top
M_BOTTOM = 90                       # 1.25" bottom (title block)

# ── Drawing style ──────────────────────────────────────────────────────
LW_OUTLINE = 1.6                    # building outline
LW_GRADE = 2.0                      # grade line
LW_OPENING = 0.75                   # window/door outlines
LW_DIM = 0.3                        # dimension lines
LW_ROOF = 1.0                       # roof edge lines
LW_SILL = 0.5                       # window sill line
LW_HATCH = 0.25                     # grade hatching
DIM_OFFSET_PT = 20                  # first dim line offset from drawing (pts)
DIM_GAP_PT = 14                     # gap between consecutive dim lines
DIM_TICK = 4                        # dim tick half-length (pts)
DIM_FONT = 7                        # dimension text size
TITLE_FONT = 14
SUBTITLE_FONT = 9
LABEL_FONT = 6.5

BLACK = (0, 0, 0)
GRAY = (0.5, 0.5, 0.5)
LIGHT_GRAY = (0.85, 0.85, 0.85)
WALL_FILL = (0.92, 0.92, 0.92)
OPENING_FILL = (1, 1, 1)            # white (cut through wall)


# ── Geometry helpers ───────────────────────────────────────────────────

def _sample_polygon(pts, segs, n_arc=80):
    """Sample outline/inner segments into a dense (E, N) point list."""
    result = []
    for seg in segs:
        p1 = pts[seg.start]
        if isinstance(seg, ArcSeg):
            c = pts[seg.center]
            R = seg.radius
            a1 = math.atan2(p1[1] - c[1], p1[0] - c[0])
            p2 = pts[seg.end]
            a2 = math.atan2(p2[1] - c[1], p2[0] - c[0])
            if seg.direction == "CW":
                sweep = (a1 - a2) % (2 * math.pi)
                for i in range(n_arc):
                    a = a1 - sweep * i / n_arc
                    result.append((c[0] + R * math.cos(a),
                                   c[1] + R * math.sin(a)))
            else:
                sweep = (a2 - a1) % (2 * math.pi)
                for i in range(n_arc):
                    a = a1 + sweep * i / n_arc
                    result.append((c[0] + R * math.cos(a),
                                   c[1] + R * math.sin(a)))
        else:
            result.append(p1)
    return result


def _scanline(polygon, h_func, n_samples=1000):
    """Scanline a polygon at evenly-spaced h values.

    Returns list of (h, [(E, N), ...]) for each h that intersects the
    polygon boundary.
    """
    h_all = [h_func(e, n) for e, n in polygon]
    h_min, h_max = min(h_all), max(h_all)
    if h_max - h_min < 1e-10:
        return []
    dh = (h_max - h_min) / (n_samples - 1)

    result = []
    n_poly = len(polygon)

    for i in range(n_samples):
        h = h_min + i * dh
        isects = []
        for j in range(n_poly):
            e1, n1 = polygon[j]
            e2, n2 = polygon[(j + 1) % n_poly]
            h1, h2 = h_func(e1, n1), h_func(e2, n2)

            if abs(h2 - h1) < 1e-10:
                if abs(h1 - h) < dh * 0.6:
                    isects.extend([(e1, n1), (e2, n2)])
                continue

            t = (h - h1) / (h2 - h1)
            if t < -0.001 or t > 1.001:
                continue
            isects.append((e1 + t * (e2 - e1), n1 + t * (n2 - n1)))

        if isects:
            result.append((h, isects))
    return result


def _seg_outward_normal(pts, seg):
    """Outward normal (exterior) unit vector at the midpoint of a segment.

    For CW outline traversal, the left normal points outward.
    """
    p1, p2 = pts[seg.start], pts[seg.end]
    if isinstance(seg, ArcSeg):
        c = pts[seg.center]
        a1 = math.atan2(p1[1] - c[1], p1[0] - c[0])
        a2 = math.atan2(p2[1] - c[1], p2[0] - c[0])
        if seg.direction == "CW":
            sweep = (a1 - a2) % (2 * math.pi)
            a_mid = a1 - sweep / 2
            # CW arc: outward = away from center
            return (math.cos(a_mid), math.sin(a_mid))
        else:
            sweep = (a2 - a1) % (2 * math.pi)
            a_mid = a1 + sweep / 2
            # CCW arc (concave): outward = toward center
            return (-math.cos(a_mid), -math.sin(a_mid))
    else:
        dx, dy = p2[0] - p1[0], p2[1] - p1[1]
        length = math.sqrt(dx * dx + dy * dy)
        if length < 1e-12:
            return (0.0, 0.0)
        # Left normal of (dx, dy) = (-dy, dx)
        return (-dy / length, dx / length)


# ── Formatting helpers ─────────────────────────────────────────────────

def _fmt_ft_in(ft):
    """Format feet as  X'-Y"  (architectural convention)."""
    total_in = abs(ft) * 12
    whole_ft = int(total_in) // 12
    remaining_in = total_in - whole_ft * 12
    # Round to nearest 1/2"
    remaining_in = round(remaining_in * 2) / 2
    if remaining_in >= 12:
        whole_ft += 1
        remaining_in -= 12
    if remaining_in == int(remaining_in):
        remaining_in = int(remaining_in)
    if remaining_in == 0:
        return f"{whole_ft}'-0\""
    return f"{whole_ft}'-{remaining_in}\""


def _fmt_slope():
    """Return roof slope notation."""
    return "2:12"


# ── Drawing helpers ────────────────────────────────────────────────────

def _draw_dim_h(shape, to_page, h1, h2, z_anchor, dim_y_pt, label=None):
    """Draw a horizontal dimension line between two h values.

    dim_y_pt is the page-y of the dimension line.
    """
    p1 = to_page(h1, z_anchor)
    p2 = to_page(h2, z_anchor)
    # Extension lines
    shape.draw_line(fitz.Point(p1.x, p1.y), fitz.Point(p1.x, dim_y_pt))
    shape.draw_line(fitz.Point(p2.x, p2.y), fitz.Point(p2.x, dim_y_pt))
    shape.finish(color=GRAY, width=LW_DIM)
    # Dimension line
    shape.draw_line(fitz.Point(p1.x, dim_y_pt), fitz.Point(p2.x, dim_y_pt))
    shape.finish(color=BLACK, width=LW_DIM)
    # Ticks
    for px in (p1.x, p2.x):
        shape.draw_line(fitz.Point(px, dim_y_pt - DIM_TICK),
                        fitz.Point(px, dim_y_pt + DIM_TICK))
    shape.finish(color=BLACK, width=LW_DIM + 0.2)
    # Label
    if label is None:
        label = _fmt_ft_in(abs(h2 - h1))
    mid_x = (p1.x + p2.x) / 2
    shape.insert_text(fitz.Point(mid_x - len(label) * 2, dim_y_pt - 3),
                      label, fontsize=DIM_FONT, color=BLACK)


def _draw_dim_v(shape, to_page, h_anchor, z1, z2, dim_x_pt, label=None):
    """Draw a vertical dimension line between two z values.

    dim_x_pt is the page-x of the dimension line.
    """
    p1 = to_page(h_anchor, z1)
    p2 = to_page(h_anchor, z2)
    # Extension lines
    shape.draw_line(fitz.Point(p1.x, p1.y), fitz.Point(dim_x_pt, p1.y))
    shape.draw_line(fitz.Point(p2.x, p2.y), fitz.Point(dim_x_pt, p2.y))
    shape.finish(color=GRAY, width=LW_DIM)
    # Dimension line
    shape.draw_line(fitz.Point(dim_x_pt, p1.y), fitz.Point(dim_x_pt, p2.y))
    shape.finish(color=BLACK, width=LW_DIM)
    # Ticks
    for py in (p1.y, p2.y):
        shape.draw_line(fitz.Point(dim_x_pt - DIM_TICK, py),
                        fitz.Point(dim_x_pt + DIM_TICK, py))
    shape.finish(color=BLACK, width=LW_DIM + 0.2)
    # Label
    if label is None:
        label = _fmt_ft_in(abs(z2 - z1))
    mid_y = (p1.y + p2.y) / 2
    # Rotate text 90 degrees for vertical dims
    shape.insert_text(fitz.Point(dim_x_pt - 3, mid_y + len(label) * 2),
                      label, fontsize=DIM_FONT, color=BLACK, rotate=90)


# ── Main elevation drawing ─────────────────────────────────────────────

def _interp_profile(profile, h):
    """Linearly interpolate z at a given h from a (h, z) profile."""
    if not profile:
        return WALL_HEIGHT_FT
    if h <= profile[0][0]:
        return profile[0][1]
    if h >= profile[-1][0]:
        return profile[-1][1]
    for i in range(len(profile) - 1):
        h0, z0 = profile[i]
        h1, z1 = profile[i + 1]
        if h0 <= h <= h1:
            if abs(h1 - h0) < 1e-10:
                return z0
            t = (h - h0) / (h1 - h0)
            return z0 + t * (z1 - z0)
    return profile[-1][1]


def _fmt_inches(ft):
    """Format feet as inches, e.g. 36\" — for opening dimensions."""
    inches = abs(ft) * 12
    # Round to nearest 0.5"
    inches = round(inches * 2) / 2
    if inches == int(inches):
        return f"{int(inches)}\""
    return f"{inches}\""


def _upper_envelope(pts_list, n_bins=300):
    """Compute the upper envelope of a set of (h, z) points.

    Groups points into h bins and takes the max z in each.
    Returns a clean list of (h, z).
    """
    if not pts_list:
        return []
    h_min = min(h for h, _ in pts_list)
    h_max = max(h for h, _ in pts_list)
    if h_max - h_min < 1e-10:
        return pts_list[:1]

    bin_w = (h_max - h_min) / n_bins
    bins = [None] * (n_bins + 1)

    # Bin the raw points
    for h, z in pts_list:
        idx = min(n_bins, max(0, int((h - h_min) / bin_w)))
        if bins[idx] is None or z > bins[idx]:
            bins[idx] = z

    # Also interpolate along segments between consecutive points
    for k in range(len(pts_list)):
        h1, z1 = pts_list[k]
        h2, z2 = pts_list[(k + 1) % len(pts_list)]
        i1 = min(n_bins, max(0, int((h1 - h_min) / bin_w)))
        i2 = min(n_bins, max(0, int((h2 - h_min) / bin_w)))
        lo, hi = min(i1, i2), max(i1, i2)
        for i in range(lo, hi + 1):
            h = h_min + i * bin_w
            dh = h2 - h1
            if abs(dh) > 1e-10:
                t = (h - h1) / dh
                if 0 <= t <= 1:
                    z = z1 + t * (z2 - z1)
                    if bins[i] is None or z > bins[i]:
                        bins[i] = z

    result = []
    for i in range(n_bins + 1):
        if bins[i] is not None:
            result.append((h_min + i * bin_w, bins[i]))
    return result


def _generate_elevation(name, h_func, depth_func, normal_test,
                        pts, outline_segs, outline_poly, roof_poly,
                        roof_z_offset, wall_openings, outer_openings,
                        filename):
    """Generate one elevation PDF."""

    is_south_north = name in ("South", "North")

    # ── 1. Project polygons to elevation coordinates ───────────────
    # Project roof polygon to (h, z_top) and (h, z_under)
    roof_top_proj = [(h_func(e, n),
                      ROOF_SLOPE * n + roof_z_offset + ROOF_THICK_FT)
                     for e, n in roof_poly]
    roof_under_proj = [(h_func(e, n),
                        ROOF_SLOPE * n + roof_z_offset)
                       for e, n in roof_poly]

    # Project outline polygon to get the base footprint h-extent
    outline_h = [h_func(e, n) for e, n in outline_poly]
    outline_n = [n for _, n in outline_poly]

    h_min_wall = min(outline_h)
    h_max_wall = max(outline_h)
    h_min_roof = min(h for h, _ in roof_top_proj)
    h_max_roof = max(h for h, _ in roof_top_proj)
    h_min = min(h_min_wall, h_min_roof)
    h_max = max(h_max_wall, h_max_roof)

    # Upper profile: upper envelope of the roof top projection
    upper_profile = _upper_envelope(roof_top_proj)

    # Make sure the profile extends to the wall edges (fill with wall height)
    if upper_profile:
        if upper_profile[0][0] > h_min_wall + 0.1:
            upper_profile.insert(0, (h_min_wall, WALL_HEIGHT_FT))
        if upper_profile[-1][0] < h_max_wall - 0.1:
            upper_profile.append((h_max_wall, WALL_HEIGHT_FT))

    z_max = max(z for _, z in upper_profile) if upper_profile else 10.0

    # Eave / ridge: scanline the roof polygon to find min/max N at each h.
    # min-N = south-closest edge (eave for south view),
    # max-N = north-farthest edge (ridge for south view).
    roof_scan = _scanline(roof_poly, h_func)
    eave_under = []   # (h, z) soffit line
    eave_top = []     # (h, z) fascia top
    ridge_half = []   # (h, z) ridge top (for label positioning)

    for h, isects in roof_scan:
        n_vals = [n for _, n in isects]
        n_min_r, n_max_r = min(n_vals), max(n_vals)

        if name == "South":
            n_near, n_far = n_min_r, n_max_r
        elif name == "North":
            n_near, n_far = n_max_r, n_min_r
        else:
            # East/West: near = south edge (low), far = north edge (high)
            n_near = n_min_r
            n_far = n_max_r

        z_eu = ROOF_SLOPE * n_near + roof_z_offset
        z_et = z_eu + ROOF_THICK_FT
        z_rt = ROOF_SLOPE * n_far + roof_z_offset + ROOF_THICK_FT

        eave_under.append((h, z_eu))
        eave_top.append((h, z_et))
        ridge_half.append((h, z_rt))

    # Representative eave Z for dimension annotation
    if eave_under:
        # Use midpoint of the eave profile
        mid_eu = eave_under[len(eave_under) // 2]
        z_eave_repr = mid_eu[1]
    else:
        z_eave_repr = WALL_HEIGHT_FT

    # ── 2. Compute visible openings ────────────────────────────────
    visible_openings = []
    for wo in wall_openings:
        seg = outline_segs[wo.seg_idx]
        nx, ny = _seg_outward_normal(pts, seg)
        if not normal_test(nx, ny):
            continue

        # Get opening horizontal extent
        p_start = pts[seg.start]
        p_end = pts[seg.end]
        if isinstance(seg, ArcSeg):
            # For arcs, sample the midpoint
            c = pts[seg.center]
            a1 = math.atan2(p_start[1] - c[1], p_start[0] - c[0])
            a2 = math.atan2(p_end[1] - c[1], p_end[0] - c[0])
            if seg.direction == "CW":
                sweep = (a1 - a2) % (2 * math.pi)
                a_s = a1 - sweep * wo.t_start
                a_e = a1 - sweep * wo.t_end
            else:
                sweep = (a2 - a1) % (2 * math.pi)
                a_s = a1 + sweep * wo.t_start
                a_e = a1 + sweep * wo.t_end
            R = seg.radius
            e_s = c[0] + R * math.cos(a_s)
            n_s = c[1] + R * math.sin(a_s)
            e_e = c[0] + R * math.cos(a_e)
            n_e = c[1] + R * math.sin(a_e)
        else:
            de = p_end[0] - p_start[0]
            dn = p_end[1] - p_start[1]
            e_s = p_start[0] + wo.t_start * de
            n_s = p_start[1] + wo.t_start * dn
            e_e = p_start[0] + wo.t_end * de
            n_e = p_start[1] + wo.t_end * dn

        h_s = h_func(e_s, n_s)
        h_e = h_func(e_e, n_e)
        h_left, h_right = min(h_s, h_e), max(h_s, h_e)

        # Vertical extent
        is_door = wo.name in ("O3", "O6")
        z_bot = 0.0 if is_door else LOWER_HEIGHT_FT
        z_top_op = WALL_HEIGHT_FT

        visible_openings.append((wo.name, h_left, h_right, z_bot, z_top_op,
                                 is_door))

    # ── 3. Create PDF ──────────────────────────────────────────────
    doc = fitz.open()
    page = doc.new_page(width=PAGE_W, height=PAGE_H)
    shape = page.new_shape()

    draw_w = PAGE_W - M_LEFT - M_RIGHT
    draw_h = PAGE_H - M_TOP - M_BOTTOM

    # Auto-scale: fit building + padding for dimension lines
    building_w = h_max - h_min
    building_h = z_max
    pad_h = max(3.0, building_w * 0.08)  # horizontal padding (ft)
    pad_top = max(1.5, building_h * 0.12)
    pad_bot = 1.5  # below grade

    total_w = building_w + 2 * pad_h
    total_h = building_h + pad_top + pad_bot

    sx = draw_w / total_w
    sy = draw_h / total_h
    scale = min(sx, sy)

    actual_w = total_w * scale
    actual_h = total_h * scale
    ox = M_LEFT + (draw_w - actual_w) / 2
    oy = M_TOP + (draw_h - actual_h) / 2

    def to_page(h, z):
        """Convert building (h, z) to page point."""
        px = ox + (h - h_min + pad_h) * scale
        py = oy + actual_h - (z + pad_bot) * scale
        return fitz.Point(px, py)

    # ── 4. Grade hatching ──────────────────────────────────────────
    gh_left = to_page(h_min - pad_h * 0.4, 0)
    gh_right = to_page(h_max + pad_h * 0.4, 0)
    spacing = 5
    n_hatches = int((gh_right.x - gh_left.x) / spacing)
    for i in range(n_hatches + 1):
        x = gh_left.x + i * spacing
        shape.draw_line(fitz.Point(x, gh_left.y),
                        fitz.Point(x - 7, gh_left.y + 7))
    shape.finish(color=LIGHT_GRAY, width=LW_HATCH)

    # Grade line
    shape.draw_line(to_page(h_min - pad_h * 0.4, 0),
                    to_page(h_max + pad_h * 0.4, 0))
    shape.finish(color=BLACK, width=LW_GRADE)

    # "GRADE" label
    gl_pt = to_page(h_max + pad_h * 0.15, 0)
    shape.insert_text(fitz.Point(gl_pt.x, gl_pt.y - 3), "GRADE",
                      fontsize=LABEL_FONT, color=GRAY)

    # ── 5. Building silhouette ─────────────────────────────────────
    # Wall polygon extends from grade up to the roof soffit (eave_under),
    # so walls seamlessly meet the bottom of the roof with no gap.

    wall_sil = [to_page(h_min_wall, 0), to_page(h_max_wall, 0)]
    if eave_under:
        # Right side up to soffit
        wall_sil.append(to_page(h_max_wall,
                                _interp_profile(eave_under, h_max_wall)))
        # Top follows soffit right-to-left (clipped to wall range)
        eu_clipped = [(h, z) for h, z in eave_under
                      if h_min_wall <= h <= h_max_wall]
        for h, z in reversed(eu_clipped):
            wall_sil.append(to_page(h, z))
        # Left side back down
        wall_sil.append(to_page(h_min_wall,
                                _interp_profile(eave_under, h_min_wall)))
    else:
        wall_sil.append(to_page(h_max_wall, WALL_HEIGHT_FT))
        wall_sil.append(to_page(h_min_wall, WALL_HEIGHT_FT))
    shape.draw_polyline(wall_sil)
    shape.finish(color=BLACK, width=LW_OUTLINE, fill=WALL_FILL,
                 closePath=True)

    # Roof slab polygon (from eave_under to roof top profile)
    # This sits on top of the wall and overhangs on both sides.
    if upper_profile and eave_under:
        roof_sil = []
        # Bottom edge of roof slab (soffit) left-to-right
        for h, z in eave_under:
            roof_sil.append(to_page(h, z))
        # Right fascia (up from soffit to roof top)
        roof_sil.append(to_page(upper_profile[-1][0], upper_profile[-1][1]))
        # Top edge of roof slab right-to-left
        for h, z in reversed(upper_profile):
            roof_sil.append(to_page(h, z))
        # Left fascia (down from roof top to soffit)
        roof_sil.append(to_page(eave_under[0][0], eave_under[0][1]))
        shape.draw_polyline(roof_sil)
        shape.finish(color=BLACK, width=LW_ROOF, fill=WALL_FILL,
                     closePath=True)

    # For South/North views, draw the near-side fascia top line (visible
    # feature: top edge of the fascia board, between soffit and ridge).
    # For East/West views the filled roof polygon silhouette is sufficient;
    # drawing additional lines would show hidden far-side edges.
    if is_south_north and eave_top:
        et_pts = [to_page(h, z) for h, z in eave_top]
        shape.draw_polyline(et_pts)
        shape.finish(color=BLACK, width=LW_ROOF)

    # Roof slope annotation for East/West (slope directly visible)
    # Standard architectural pitch symbol: compact right triangle with
    # hypotenuse ON the roof slope.  Right angle on the ridge (high) side.
    # Horizontal leg labelled "12", vertical leg labelled "2".
    if not is_south_north and len(upper_profile) > 4:
        # Fixed symbol size: horizontal leg ≈ 1" on paper → run_ft = 1/scale
        sym_page_w = 60.0   # ~0.8" on paper for horizontal leg
        run_ft = sym_page_w / scale
        rise_ft = run_ft * ROOF_SLOPE

        # Anchor at ~45% along the roof profile (from low side)
        idx_pos = int(len(upper_profile) * 0.45)
        h_anchor, z_anchor = upper_profile[idx_pos]

        if name == "East":
            # Roof slopes up to the right; right angle at upper-right
            h_lo = h_anchor
            z_lo = z_anchor
            h_hi = h_anchor + run_ft
            z_hi = z_anchor + rise_ft
            p_hyp_lo = to_page(h_lo, z_lo)          # low end of hypotenuse
            p_corner = to_page(h_hi, z_lo)           # right-angle corner
            p_hyp_hi = to_page(h_hi, z_hi)           # high end of hypotenuse
            # "2" to the right of vertical leg
            v_lbl = fitz.Point(p_corner.x + 3,
                               (p_corner.y + p_hyp_hi.y) / 2 + 3)
            # "12" below horizontal leg
            h_lbl = fitz.Point((p_hyp_lo.x + p_corner.x) / 2 - 5,
                               p_corner.y + 10)
        else:
            # West: roof slopes up to the left; right angle at upper-left
            h_lo = h_anchor
            z_lo = z_anchor
            h_hi = h_anchor - run_ft
            z_hi = z_anchor + rise_ft
            p_hyp_lo = to_page(h_lo, z_lo)          # low end (right)
            p_corner = to_page(h_hi, z_lo)           # right-angle corner (left)
            p_hyp_hi = to_page(h_hi, z_hi)           # high end (left)
            # "2" to the left of vertical leg
            v_lbl = fitz.Point(p_corner.x - 12,
                               (p_corner.y + p_hyp_hi.y) / 2 + 3)
            # "12" below horizontal leg
            h_lbl = fitz.Point((p_hyp_lo.x + p_corner.x) / 2 - 5,
                               p_corner.y + 10)

        # Draw: low → corner (horiz) → high (vert), close = hypotenuse
        tri = [p_hyp_lo, p_corner, p_hyp_hi]
        shape.draw_polyline(tri)
        shape.finish(color=BLACK, width=LW_DIM, closePath=True)
        shape.insert_text(v_lbl, "2", fontsize=LABEL_FONT, color=BLACK)
        shape.insert_text(h_lbl, "12", fontsize=LABEL_FONT, color=BLACK)

    # ── 6. Openings ────────────────────────────────────────────────
    for _, h_l, h_r, z_b, z_t, is_door in visible_openings:
        # Draw opening rectangle (white fill = cut through wall)
        corners = [
            to_page(h_l, z_b), to_page(h_r, z_b),
            to_page(h_r, z_t), to_page(h_l, z_t),
        ]
        shape.draw_polyline(corners)
        shape.finish(color=BLACK, width=LW_OPENING, fill=OPENING_FILL,
                     closePath=True)

        if not is_door:
            # Window sill line
            shape.draw_line(to_page(h_l, z_b), to_page(h_r, z_b))
            shape.finish(color=BLACK, width=LW_SILL)
            # Window cross (mullion indication)
            h_mid = (h_l + h_r) / 2
            z_mid = (z_b + z_t) / 2
            shape.draw_line(to_page(h_mid, z_b), to_page(h_mid, z_t))
            shape.finish(color=GRAY, width=0.3)
            shape.draw_line(to_page(h_l, z_mid), to_page(h_r, z_mid))
            shape.finish(color=GRAY, width=0.3)
        else:
            # Door threshold line
            shape.draw_line(to_page(h_l, 0), to_page(h_r, 0))
            shape.finish(color=BLACK, width=LW_OPENING)

        # Opening dimension label (width × height in inches)
        w_in = _fmt_inches(h_r - h_l)
        h_in = _fmt_inches(z_t - z_b)
        dim_label = f"{w_in} x {h_in}"
        h_ctr = (h_l + h_r) / 2
        z_ctr = (z_b + z_t) / 2
        lbl_pt = to_page(h_ctr, z_ctr)
        tw = len(dim_label) * LABEL_FONT * 0.3
        shape.insert_text(fitz.Point(lbl_pt.x - tw, lbl_pt.y + 2),
                          dim_label, fontsize=LABEL_FONT, color=BLACK)

    # ── 7. Dimension lines ─────────────────────────────────────────
    # --- Vertical dimensions (left side) ---
    dim_x = to_page(h_min, 0).x - DIM_OFFSET_PT

    # Grade to eave (roof underside at viewer-side)
    z_eave = z_eave_repr
    _draw_dim_v(shape, to_page, h_min, 0, z_eave, dim_x,
                label=_fmt_ft_in(z_eave))

    # Grade to ridge (overall height)
    dim_x2 = dim_x - DIM_GAP_PT
    max_z_entry = max(upper_profile, key=lambda x: x[1])
    _draw_dim_v(shape, to_page, h_min, 0, max_z_entry[1], dim_x2,
                label=_fmt_ft_in(max_z_entry[1]))

    # Eave to ridge (if South/North where they differ significantly)
    if is_south_north and abs(max_z_entry[1] - z_eave) > 0.5:
        dim_x3 = dim_x - 2 * DIM_GAP_PT
        _draw_dim_v(shape, to_page, h_min, z_eave, max_z_entry[1], dim_x3,
                    label=_fmt_ft_in(max_z_entry[1] - z_eave))

    # Wall height (grade to wall top = 6'8")
    dim_x_r = to_page(h_max, 0).x + DIM_OFFSET_PT
    _draw_dim_v(shape, to_page, h_max, 0, WALL_HEIGHT_FT, dim_x_r,
                label=_fmt_ft_in(WALL_HEIGHT_FT))

    # --- Horizontal dimensions (bottom) ---
    dim_y = to_page(h_min, 0).y + DIM_OFFSET_PT

    # Overall width
    _draw_dim_h(shape, to_page, h_min if not upper_profile else
                upper_profile[0][0],
                h_max if not upper_profile else upper_profile[-1][0],
                0, dim_y,
                label=_fmt_ft_in(abs(upper_profile[-1][0] -
                                     upper_profile[0][0])
                                 if upper_profile else 0))

    # ── 8. Roof slope annotation (South/North) ────────────────────
    if is_south_north:
        # Add "2:12 SLOPE" text near the ridge
        if ridge_half:
            mid_rt = ridge_half[len(ridge_half) // 2]
            rt_pt = to_page(mid_rt[0], mid_rt[1])
            shape.insert_text(fitz.Point(rt_pt.x - 20, rt_pt.y - 6),
                              f"ROOF SLOPE {_fmt_slope()}",
                              fontsize=LABEL_FONT, color=BLACK)

    # ── 9. Title block ─────────────────────────────────────────────
    # Border
    shape.draw_rect(fitz.Rect(18, 18, PAGE_W - 18, PAGE_H - 18))
    shape.finish(color=BLACK, width=0.75)

    # Title block area (bottom strip)
    tb_y = PAGE_H - M_BOTTOM + 18
    shape.draw_line(fitz.Point(18, tb_y), fitz.Point(PAGE_W - 18, tb_y))
    shape.finish(color=BLACK, width=0.5)

    # Title text
    title = f"{name.upper()} ELEVATION"
    shape.insert_text(fitz.Point(36, tb_y + 22), title,
                      fontsize=TITLE_FONT, color=BLACK)

    # Subtitle
    shape.insert_text(fitz.Point(36, tb_y + 36),
                      "Accessory Dwelling Unit — 2:12 Roof",
                      fontsize=SUBTITLE_FONT, color=GRAY)

    # Scale
    scale_ft_per_pt = 1.0 / scale
    # Find a nice scale ratio
    # At current scale, how many inches on paper = 1 foot real?
    in_per_ft = scale / 72.0
    # Express as fraction
    nice_scales = [
        (1/48, '1/4" = 1\'-0"'),
        (1/32, '3/8" = 1\'-0"'),
        (1/24, '1/2" = 1\'-0"'),
        (1/16, '3/4" = 1\'-0"'),
        (1/12, '1" = 1\'-0"'),
    ]
    # Pick closest standard scale
    best_label = f"Scale: 1\" = {1/in_per_ft:.1f}'"
    for ns, nl in nice_scales:
        if abs(in_per_ft - ns) / ns < 0.15:
            best_label = f"Scale: {nl}"
            break
    shape.insert_text(fitz.Point(400, tb_y + 22), best_label,
                      fontsize=SUBTITLE_FONT, color=BLACK)

    # Sheet number
    sheet_map = {"South": "A-201", "North": "A-202",
                 "East": "A-203", "West": "A-204"}
    shape.insert_text(fitz.Point(650, tb_y + 22),
                      f"Sheet {sheet_map.get(name, 'A-200')}",
                      fontsize=SUBTITLE_FONT, color=BLACK)

    # Date placeholder
    import datetime
    date_str = datetime.date.today().strftime("%Y-%m-%d")
    shape.insert_text(fitz.Point(400, tb_y + 36),
                      f"Date: {date_str}",
                      fontsize=SUBTITLE_FONT - 1, color=GRAY)

    # ── 10. Scale bar ──────────────────────────────────────────────
    sb_y = tb_y + 50
    sb_x = 400
    bar_ft = 5  # 5-foot scale bar
    bar_len = bar_ft * scale
    shape.draw_line(fitz.Point(sb_x, sb_y), fitz.Point(sb_x + bar_len, sb_y))
    shape.finish(color=BLACK, width=1.0)
    # Ticks at each foot
    for i in range(bar_ft + 1):
        x = sb_x + i * scale
        tick_h = 4 if i % 5 == 0 else 2
        shape.draw_line(fitz.Point(x, sb_y - tick_h),
                        fitz.Point(x, sb_y + tick_h))
        if i % 5 == 0 or i == bar_ft:
            shape.insert_text(fitz.Point(x - 3, sb_y + 12),
                              f"{i}'", fontsize=5, color=BLACK)
    shape.finish(color=BLACK, width=0.5)

    # ── Commit and save ────────────────────────────────────────────
    shape.commit()
    out_path = os.path.join(_DIR, filename)
    doc.save(out_path)
    doc.close()
    print(f"wrote {out_path}")


# ── Entry point ────────────────────────────────────────────────────────

def generate():
    """Compute geometry and generate all 4 elevation PDFs."""
    # Outline geometry
    geo = compute_outline_geometry()
    pts = dict(geo.fp_pts)
    outline_segs = geo.outline_segs
    radii = geo.radii

    # Inner wall points (W-series) — mutates pts in place
    inner_segs = compute_inner_walls(outline_segs, pts, WALL_OUTER, radii)

    # Roof geometry
    roof_geo = compute_roof_geometry(pts, radii)
    roof_poly = roof_polyline(roof_geo)

    # Roof Z offset
    ref_y = pts["F18"][1]
    roof_z_offset = ROOF_REF_ELEV_FT - ROOF_SLOPE * ref_y

    # Dense outline polygon
    outline_poly = _sample_polygon(pts, outline_segs)

    # Openings
    from floorplan.layout import compute_interior_layout
    from floorplan.openings import (compute_outer_openings,
                                    outer_to_wall_openings)
    inner_poly = path_polygon(inner_segs, pts)
    layout = compute_interior_layout(pts, inner_poly)
    outer_openings = compute_outer_openings(pts, layout)
    wall_openings = outer_to_wall_openings(outer_openings, outline_segs, pts)

    # Generate each elevation
    elevations = [
        ("South", lambda e, n: e,  lambda e, n: n,
         lambda nx, ny: ny < -0.5, "elevation_south.pdf"),
        ("North", lambda e, n: -e, lambda e, n: -n,
         lambda nx, ny: ny > 0.5,  "elevation_north.pdf"),
        ("East",  lambda e, n: n,  lambda e, n: -e,
         lambda nx, ny: nx > 0.5,  "elevation_east.pdf"),
        ("West",  lambda e, n: -n, lambda e, n: e,
         lambda nx, ny: nx < -0.5, "elevation_west.pdf"),
    ]

    for elev_name, h_func, depth_func, normal_test, filename in elevations:
        _generate_elevation(
            elev_name, h_func, depth_func, normal_test,
            pts, outline_segs, outline_poly, roof_poly,
            roof_z_offset, wall_openings, outer_openings, filename)


if __name__ == "__main__":
    generate()
