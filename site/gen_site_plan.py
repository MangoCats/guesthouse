"""Generate site_plan.pdf: building outline overlaid on site survey."""

import datetime
import math
import sys
import os
from collections import namedtuple

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

import fitz  # pymupdf
from shared.geometry import vert_isects, path_polygon, segment_polyline, seg_vecs
from shared.types import LineSeg, ArcSeg
from shared.svg import git_describe
from shared.survey import compute_pt1

# Survey coordinate calibration constants (PDF coords from least-squares fitting)
LINE_TOP = (698.9, 55.2)     # 251.53' meets 216.73' (upper-right corner)
LINE_BOT = (817.9, 557.8)    # 216.73' meets 275.08' (lower-right corner)
BOT_LEFT = (160.0, 561.9)    # 275.08' left endpoint
TL_251 = (108.0, 174.5)      # 251.53' upper-left corner of parcel

# Named corner aliases (intersections of the four labeled property lines).
#
# IMPORTANT — survey orientation: the recorded plat is NOT drawn north-up.
# True north points toward the LEFT of the page (and slightly down); a good
# working approximation is N≈left, E≈top (hence S≈right, W≈bottom).  The
# page-position constants above therefore map to these TRUE compass corners,
# which is what these aliases name:
#
#     property line   compass boundary   page edge
#     251.53'         EAST               top    (TL_251 .. LINE_TOP)
#     216.73'         SOUTH              right  (LINE_TOP .. LINE_BOT)
#     275.08'         WEST               bottom (LINE_BOT .. BOT_LEFT)
#     163.69'         NORTH              left   (BOT_LEFT .. TL_251)
#
# NOTE: the building/FC frame uses true Easting/Northing, so these compass
# names also match the parcel corners' positions in building coordinates.
CORNER_NE = TL_251     # page top-left:     East 251.53' meets North 163.69'
CORNER_SE = LINE_TOP   # page top-right:    East 251.53' meets South 216.73'
CORNER_SW = LINE_BOT   # page bottom-right: South 216.73' meets West 275.08'
CORNER_NW = BOT_LEFT   # page bottom-left:  West 275.08' meets North 163.69'

# Placement constraints: survey point distances from property lines (feet)
P45_DIST_216 = 11.0    # P4/P5 tangent-normal distance from 216.73' line
P3_DIST_275 = 25.5     # P3 tangent-normal distance from 275.08' line

# Existing FRAME & STONE RESIDENCE lower-right corner (14.4'/28.2' wall midline isect)
RESIDENCE_LR = (661.35, 380.80)

COLOR_PROPOSED = (0, 0, 0.6)  # dark blue for proposed/new construction

# Property line vectors (precomputed from corner coordinates)
_LINE_216_DX = LINE_BOT[0] - LINE_TOP[0]
_LINE_216_DY = LINE_BOT[1] - LINE_TOP[1]
_LINE_216_LEN = math.hypot(_LINE_216_DX, _LINE_216_DY)
_LINE_275_DX = LINE_BOT[0] - BOT_LEFT[0]
_LINE_275_DY = LINE_BOT[1] - BOT_LEFT[1]
_LINE_275_LEN = math.hypot(_LINE_275_DX, _LINE_275_DY)

# Rendering constants
DIM_LABEL_FS = 6.0
DIM_LABEL_OFFSET = -3.5
SETBACK_LABEL_FS = 9.0
OUTLINE_STROKE_W = 1.6
BLDG_LABEL_FS = 8.0

SitePlanData = namedtuple("SitePlanData", [
    "pts",              # building points dict
    "building_to_pdf",  # transform function (E,N) → (pdf_x, pdf_y)
    "rotation_deg",     # rotation angle in degrees
    "f15_pdf",          # F15 position in PDF coords
    "ew_dim_ft",        # E-W external dimension (_site_se_pt → _site_w_pt projected)
    "ns_dim_ft",        # N-S external dimension (F18→F6 along building N-S)
    "arc_dim_ft",       # arc midpoint of F02-F03 to F24b-F26 surface (along C02 line)
    "min_setback_216",  # min perpendicular dist of any F point from 216.73' line (ft)
    "min_setback_275",  # min perpendicular dist of any F point from 275.08' line (ft)
    "min275_pdf",       # PDF coords of the F-series point closest to 275.08' line
    "draw_poly",        # interpolated building outline (building coords)
    "inner_poly",       # inner wall polygon
    "span_s_pdf",       # N-S span south endpoint in PDF coords
    "span_n_pdf",       # N-S span north endpoint in PDF coords
    "f2_pdf",           # F2 position in PDF coords
    "SCALE",            # PDF pts per foot (2.4)
    "f_series_pdf",     # dict of F-series + FC points in PDF coords
    "residence_dist_ft", # distance from RESIDENCE_LR to closest F point (ft)
    "residence_closest",  # name of closest F point to RESIDENCE_LR
    "p_series_pdf",     # dict of angularly-corrected P-series points in PDF coords
])


def _find_site_refs_from_db(outline_segs, pts):
    """Find reference points from a DB chain outline, stored under internal alias names.

    All returned keys use the ``_site_`` prefix to avoid overwriting any real
    DB point coordinates in pts.  Callers must use these internal names (not
    legacy F-series names) when accessing the returned reference points.

    Keys returned:
      _site_sf_start  — start of the south-face LineSeg (most SW-directed)
      _site_sf_end    — end   of the south-face LineSeg
      _site_se_pt     — start of the CW arc preceding the south face (SE corner ref)
      _site_w_pt      — westernmost outline vertex (E-W dimension west ref)
      _site_n_pt      — northernmost outline vertex (N-S dimension north ref)
      _site_s_pt      — southernmost outline vertex (N-S dimension south ref)
      W9, W10         — virtual E-W direction points (centroid ± 1 ft E)
      W2, W5          — virtual N-S direction points (centroid ± 1 ft N)
    """
    # Target: seed south face goes SW, unit ≈ (-0.864, -0.504)
    _TGT_E, _TGT_N = -0.864, -0.504

    seg_from_end = {s.end: s for s in outline_segs}

    all_names = list(dict.fromkeys(
        n for s in outline_segs for n in (s.start, s.end) if n in pts))

    # Find south-face LineSeg: best dot-product with SW target,
    # length ≥ 3 ft, bonus if preceded by a CW arc.
    best_seg = None
    best_score = -1.0
    for seg in outline_segs:
        if not isinstance(seg, LineSeg):
            continue
        if seg.start not in pts or seg.end not in pts:
            continue
        dx = pts[seg.end][0] - pts[seg.start][0]
        dy = pts[seg.end][1] - pts[seg.start][1]
        length = math.hypot(dx, dy)
        if length < 3.0:
            continue
        dot = (dx * _TGT_E + dy * _TGT_N) / length
        if dot < 0.5:
            continue
        prev = seg_from_end.get(seg.start)
        bonus = 0.1 if (prev is not None and isinstance(prev, ArcSeg)
                        and prev.direction == "CW") else 0.0
        score = dot + bonus
        if score > best_score:
            best_score = score
            best_seg = seg

    if best_seg is None:
        # Fallback: unconstrained best-dot LineSeg
        for seg in outline_segs:
            if not isinstance(seg, LineSeg):
                continue
            if seg.start not in pts or seg.end not in pts:
                continue
            dx = pts[seg.end][0] - pts[seg.start][0]
            dy = pts[seg.end][1] - pts[seg.start][1]
            length = math.hypot(dx, dy) or 1.0
            dot = (dx * _TGT_E + dy * _TGT_N) / length
            if dot > best_score:
                best_score = dot
                best_seg = seg

    result = {}
    result["_site_sf_start"] = pts[best_seg.start]
    result["_site_sf_end"]   = pts[best_seg.end]

    # SE corner reference: start of the CW arc preceding the south face
    prev_arc = seg_from_end.get(best_seg.start)
    if prev_arc is not None and isinstance(prev_arc, ArcSeg):
        result["_site_se_pt"] = pts[prev_arc.start]
    else:
        result["_site_se_pt"] = pts[max(all_names, key=lambda n: pts[n][0])]

    # Axis-extreme reference points (building is approx. axis-aligned in FC coords)
    result["_site_w_pt"] = pts[min(all_names, key=lambda n: pts[n][0])]   # westernmost
    result["_site_n_pt"] = pts[max(all_names, key=lambda n: pts[n][1])]   # northernmost
    result["_site_s_pt"] = pts[min(all_names, key=lambda n: pts[n][1])]   # southernmost

    # Virtual W-series direction points for building-axis unit vectors
    _cx = sum(pts[n][0] for n in all_names) / len(all_names)
    _cy = sum(pts[n][1] for n in all_names) / len(all_names)
    result["W9"]  = (_cx - 1.0, _cy)
    result["W10"] = (_cx + 1.0, _cy)
    result["W2"]  = (_cx, _cy - 1.0)
    result["W5"]  = (_cx, _cy + 1.0)

    return result


def build_site_plan_data(gd=None):
    """Compute all site plan geometry — no PDF I/O.

    If gd (GeneratorData) is provided, uses it as the geometry source.
    Otherwise constructs one from the hardcoded procedural modules.
    """
    if gd is None:
        from floorplan.gen_floorplan import build_floorplan_data
        data = build_floorplan_data()
        pts = data.pts
        outer_poly = data.outer_poly
        inner_poly = data.inner_poly
        _outline_struct_names = (
            [f"F{i}" for i in range(1, 19) if i not in (3, 4)] + ["F11a", "F11b"])
        _f_pdf_names = _outline_struct_names + ["FC"]
        # Internal aliases for the original F-series reference points
        pts["_site_sf_start"] = pts["F16"]
        pts["_site_sf_end"]   = pts["F17"]
        pts["_site_se_pt"]    = pts["F15"]
        pts["_site_w_pt"]     = pts["F2"]
        pts["_site_n_pt"]     = pts["F6"]
        pts["_site_s_pt"]     = pts["F18"]
        p45_dist = P45_DIST_216
        p3_dist  = P3_DIST_275
    else:
        pts = dict(gd.pts)  # copy so we can add aliases without mutating gd
        outer_poly = gd.outline_poly
        inner_poly = gd.inner_poly
        # Add internal reference aliases (_site_*) from DB chain geometry.
        # Also add virtual axis-direction points (W9/W10/W2/W5) if not already
        # present as real inner-wall points (original non-zero-padded DB has them).
        _site_refs = _find_site_refs_from_db(gd.outline_segs, pts)
        pts.update({k: v for k, v in _site_refs.items() if k.startswith("_site_")})
        if "W9" not in pts:
            for k in ("W9", "W10", "W2", "W5"):
                if k in _site_refs:
                    pts[k] = _site_refs[k]
        _outline_struct_names = list(dict.fromkeys(
            n for s in gd.outline_segs for n in (s.start, s.end) if n in pts))
        _f_pdf_names = _outline_struct_names + (["FC"] if "FC" in pts else [])
        # Setback distances from DB constants (authoritative)
        p45_dist = gd.constants.get("SITE_P4_DIST_216", P45_DIST_216)
        p3_dist  = gd.constants.get("SITE_P3_DIST_275", P3_DIST_275)

    # --- Survey coordinate calibration ---
    # Scale: 1 inch = 30 ft on the survey; 1 inch = 72 PDF pts
    SCALE = 72.0 / 30.0  # 2.4 PDF pts per foot

    # Direction of 216.73' line in PDF coords (x-right, y-down)
    ldx, ldy, llen = _LINE_216_DX, _LINE_216_DY, _LINE_216_LEN

    # In real-world coords (E=x-right, N=y-up), the line direction is:
    # dE = ldx = +119, dN = -ldy = -502.6
    # Angle of property line from E-axis (real-world):
    prop_angle = math.atan2(-ldy, ldx)  # atan2(-502.6, 119) ≈ -76.7°

    # South-face direction in building coords
    f16 = pts["_site_sf_start"]
    f17 = pts["_site_sf_end"]
    f16f17_angle = math.atan2(f17[1] - f16[1], f17[0] - f16[0])

    # Rotation needed: rotate building so F16→F17 is parallel to property line
    rotation = prop_angle - f16f17_angle  # ≈ 73° CCW

    # --- Rotation and placement ---
    cos_r = math.cos(rotation)
    sin_r = math.sin(rotation)

    def rotate_pt(e, n, ce, cn):
        """Rotate point (e,n) by `rotation` around center (ce,cn)."""
        de = e - ce
        dn = n - cn
        return (ce + de * cos_r - dn * sin_r,
                cn + de * sin_r + dn * cos_r)

    # SE corner reference point for placement
    f15 = pts["_site_se_pt"]

    # --- Constraint-based placement (2×2 linear system) ---
    # After rotation, each building point has a fixed PDF offset from F15.
    # Unknowns: f15_pdf_x, f15_pdf_y (F15's position on the PDF page).
    # Constraints use P4/P5 and P3 distances from property lines.

    p4 = pts["P4"]
    p3 = pts["P3"]

    # PDF offset of P4 from F15 (fixed after rotation)
    off_p4_x = ((p4[0] - f15[0]) * cos_r - (p4[1] - f15[1]) * sin_r) * SCALE
    off_p4_y = -((p4[0] - f15[0]) * sin_r + (p4[1] - f15[1]) * cos_r) * SCALE

    # PDF offset of P3 from F15 (fixed after rotation)
    off_p3_x = ((p3[0] - f15[0]) * cos_r - (p3[1] - f15[1]) * sin_r) * SCALE
    off_p3_y = -((p3[0] - f15[0]) * sin_r + (p3[1] - f15[1]) * cos_r) * SCALE

    # Constraint A: P4 is P45_DIST_216 inside the 216.73' line.
    # Signed distance to left of LINE_TOP→LINE_BOT = property interior.
    a1 = -ldy / llen
    b1 = ldx / llen
    c1 = (p45_dist * SCALE
          + a1 * (LINE_TOP[0] - off_p4_x) + b1 * (LINE_TOP[1] - off_p4_y))

    # Constraint B: P3 is p3_dist inside the 275.08' line.
    # Direction BOT_LEFT→LINE_BOT; interior is to the right.
    bdx, bdy, blen = _LINE_275_DX, _LINE_275_DY, _LINE_275_LEN
    a2 = bdy / blen
    b2 = -bdx / blen
    c2 = (p3_dist * SCALE
          + a2 * (BOT_LEFT[0] - off_p3_x) + b2 * (BOT_LEFT[1] - off_p3_y))

    # Solve: a1*fx + b1*fy = c1,  a2*fx + b2*fy = c2
    det = a1 * b2 - a2 * b1
    f15_pdf_x = (c1 * b2 - c2 * b1) / det
    f15_pdf_y = (a1 * c2 - a2 * c1) / det

    def building_to_pdf(e, n):
        """Transform building coords (E,N) → PDF coords (x,y)."""
        re, rn = rotate_pt(e, n, f15[0], f15[1])
        pdf_x = f15_pdf_x + (re - f15[0]) * SCALE
        pdf_y = f15_pdf_y - (rn - f15[1]) * SCALE  # N-up → y-down
        return pdf_x, pdf_y

    # --- Build drawing polygon (outline offset inward by half stroke width) ---
    WALL_T = 8.0 / 12.0
    half_stroke_ft = (OUTLINE_STROKE_W / 2.0) / SCALE
    frac = half_stroke_ft / WALL_T

    n_out = len(outer_poly)
    n_inn = len(inner_poly)
    draw_poly = []
    for i, (oe, on) in enumerate(outer_poly):
        j = round(i * (n_inn - 1) / (n_out - 1))
        ie, inn = inner_poly[j]
        draw_poly.append((oe + frac * (ie - oe), on + frac * (inn - on)))

    # --- External dimensions (building-aligned, matching floorplan dim15/dim13) ---
    _bld_ew, _ = seg_vecs(pts["W9"], pts["W10"])   # building E-W direction
    _bld_ns, _ = seg_vecs(pts["W2"], pts["W5"])     # building N-S direction

    _df_ew = (pts["_site_se_pt"][0] - pts["_site_w_pt"][0],
              pts["_site_se_pt"][1] - pts["_site_w_pt"][1])
    ew_dim_ft = abs(_df_ew[0] * _bld_ew[0] + _df_ew[1] * _bld_ew[1])

    # NS dim: northernmost point in the western 33% of outline → shoot ray south.
    # Works for any outline shape (arc or line north face) without naming specific points.
    _segs_for_ns = gd.outline_segs if gd is not None else []
    arc_dim_ft = 0.0
    if _segs_for_ns:
        # E extent of all outline vertices
        _all_e = [pts[n][0] for s in _segs_for_ns
                  for n in (s.start, s.end) if n in pts]
        _e_lo, _e_hi = min(_all_e), max(_all_e)
        _e_west_limit = _e_lo + 0.33 * (_e_hi - _e_lo)  # western 33%

        def _angle_on_arc(a_test, a_start, a_end, direction):
            """Robust check: is a_test within the arc sweep from a_start to a_end?
            Uses modular arithmetic to avoid atan2 ±π sign ambiguity."""
            _2pi = 2 * math.pi
            if direction == "CW":
                sweep = (a_start - a_end) % _2pi
                to_test = (a_start - a_test) % _2pi
            else:
                sweep = (a_end - a_start) % _2pi
                to_test = (a_test - a_start) % _2pi
            return to_test <= sweep + 1e-9

        def _seg_north_extremum(seg):
            """Return (E, N) of the northernmost point of seg within the western 33%."""
            if seg.start not in pts or seg.end not in pts:
                return None
            p1, p2 = pts[seg.start], pts[seg.end]
            candidates = []
            for p in (p1, p2):
                if p[0] <= _e_west_limit:
                    candidates.append(p)
            if isinstance(seg, ArcSeg) and seg.center in pts:
                cx, cy = pts[seg.center]
                R = math.hypot(p1[0]-cx, p1[1]-cy)
                tip = (cx, cy + R)
                if tip[0] <= _e_west_limit:
                    a1 = math.atan2(p1[1]-cy, p1[0]-cx)
                    a2 = math.atan2(p2[1]-cy, p2[0]-cx)
                    if _angle_on_arc(math.pi / 2.0, a1, a2, seg.direction):
                        candidates.append(tip)
            return max(candidates, key=lambda p: p[1]) if candidates else None

        _north_candidates = [_seg_north_extremum(s) for s in _segs_for_ns]
        _north_candidates = [p for p in _north_candidates if p is not None]
        if _north_candidates:
            _E_n, _N_n = max(_north_candidates, key=lambda p: p[1])

            def _ray_south_isects(e_ray, n_start):
                """N values where the vertical ray E=e_ray crosses outline below n_start."""
                hits = []
                for seg in _segs_for_ns:
                    if seg.start not in pts or seg.end not in pts:
                        continue
                    p1, p2 = pts[seg.start], pts[seg.end]
                    if isinstance(seg, LineSeg):
                        dx = p2[0] - p1[0]
                        if abs(dx) < 1e-9:
                            continue
                        t = (e_ray - p1[0]) / dx
                        if -1e-9 <= t <= 1+1e-9:
                            n_hit = p1[1] + t * (p2[1] - p1[1])
                            if n_hit < n_start - 1e-6:
                                hits.append(n_hit)
                    elif isinstance(seg, ArcSeg):
                        if seg.center not in pts:
                            continue
                        cx, cy = pts[seg.center]
                        R = math.hypot(p1[0]-cx, p1[1]-cy)
                        disc = R*R - (e_ray - cx)**2
                        if disc < 0:
                            continue
                        for sign in (+1, -1):
                            n_hit = cy + sign * math.sqrt(disc)
                            if n_hit >= n_start - 1e-6:
                                continue
                            a_hit = math.atan2(n_hit-cy, e_ray-cx)
                            a1 = math.atan2(p1[1]-cy, p1[0]-cx)
                            a2 = math.atan2(p2[1]-cy, p2[0]-cx)
                            if not _angle_on_arc(a_hit, a1, a2, seg.direction):
                                continue
                            hits.append(n_hit)
                return hits

            _hits = _ray_south_isects(_E_n, _N_n)
            if _hits:
                _N_s = max(_hits)   # closest face south of the north extremum
                arc_dim_ft = _N_n - _N_s
                pts["_arc_dim_start"] = (_E_n, _N_n)
                pts["_arc_dim_end"]   = (_E_n, _N_s)

    _df_ns = (pts["_site_n_pt"][0] - pts["_site_s_pt"][0],
              pts["_site_n_pt"][1] - pts["_site_s_pt"][1])
    ns_dim_ft = abs(_df_ns[0] * _bld_ns[0] + _df_ns[1] * _bld_ns[1])

    # --- N-S Interior Max Span (dimension line position only) ---
    _inch = 1.0 / 12.0
    _e_min = min(p[0] for p in inner_poly)
    _e_max = max(p[0] for p in inner_poly)
    _best_span, _best_e, _best_s, _best_n = 0, 0, 0, 0
    _e = _e_min
    while _e <= _e_max + 1e-9:
        _ns = vert_isects(inner_poly, _e)
        if len(_ns) >= 2:
            _s, _n = min(_ns), max(_ns)
            if _n - _s > _best_span:
                _best_span, _best_e, _best_s, _best_n = _n - _s, _e, _s, _n
        _e += _inch

    span_s_pdf = building_to_pdf(_best_e, _best_s)
    span_n_pdf = building_to_pdf(_best_e, _best_n)

    f2_pdf = building_to_pdf(*pts["_site_w_pt"])
    f15_pdf = (f15_pdf_x, f15_pdf_y)

    # --- Outline PDF coordinates (F-series or DB chain names) ---
    f_series_pdf = {name: building_to_pdf(*pts[name]) for name in _f_pdf_names
                    if name in pts}

    # --- Min setback distances (outline structural pts, excluding centroid pts) ---
    min_setback_216 = min(
        ((pt[0] - LINE_TOP[0]) * (-ldy) + (pt[1] - LINE_TOP[1]) * ldx)
        / (llen * SCALE)
        for pt in (f_series_pdf[n] for n in _outline_struct_names if n in f_series_pdf))
    _min275_best_pt, min_setback_275 = None, float("inf")
    for _n in _outline_struct_names:
        if _n not in f_series_pdf:
            continue
        _pt = f_series_pdf[_n]
        _d = (((_pt[0] - BOT_LEFT[0]) * bdy - (_pt[1] - BOT_LEFT[1]) * bdx)
              / (blen * SCALE))
        if _d < min_setback_275:
            min_setback_275, _min275_best_pt = _d, _pt

    # --- Distance from existing residence corner to closest outline point ---
    _res_best_name, _res_best_dist = None, float("inf")
    for n in _outline_struct_names:
        if n not in f_series_pdf:
            continue
        pt = f_series_pdf[n]
        d = math.hypot(pt[0] - RESIDENCE_LR[0], pt[1] - RESIDENCE_LR[1])
        if d < _res_best_dist:
            _res_best_dist, _res_best_name = d, n
    residence_dist_ft = _res_best_dist / SCALE

    # --- P-series PDF coordinates (building coords → PDF, no angular correction) ---
    _P_NAMES = ["POB", "P2", "P3", "P4", "P5",
                "T1", "T2", "T3", "PA", "PX", "TC1", "TC2", "TC3"]
    p_series_pdf = {n: building_to_pdf(*pts[n]) for n in _P_NAMES}

    # PT1: tangency point of TC1 arc with P4-P5 line extension
    # (projected onto arc so path_polygon arc segments close exactly)
    _pt1 = compute_pt1(pts, 10.0)
    pts["PT1"] = _pt1
    p_series_pdf["PT1"] = building_to_pdf(*_pt1)

    return SitePlanData(
        pts=pts,
        building_to_pdf=building_to_pdf,
        rotation_deg=math.degrees(rotation),
        f15_pdf=f15_pdf,
        ew_dim_ft=ew_dim_ft,
        ns_dim_ft=ns_dim_ft,
        arc_dim_ft=arc_dim_ft,
        min_setback_216=min_setback_216,
        min_setback_275=min_setback_275,
        min275_pdf=_min275_best_pt,
        draw_poly=draw_poly,
        inner_poly=inner_poly,
        span_s_pdf=span_s_pdf,
        span_n_pdf=span_n_pdf,
        f2_pdf=f2_pdf,
        SCALE=SCALE,
        f_series_pdf=f_series_pdf,
        residence_dist_ft=residence_dist_ft,
        residence_closest=_res_best_name,
        p_series_pdf=p_series_pdf,
    )


def _label_aabb(cx, cy, tw, fs, angle_deg):
    """Axis-aligned bounding box of a rotated text label centered at (cx, cy).

    Returns (x_min, y_min, x_max, y_max).
    """
    th = fs * 1.2
    a = math.radians(angle_deg)
    ca, sa = math.cos(a), math.sin(a)
    hw, hh = tw / 2.0, th / 2.0
    corners = [(-hw, -hh), (hw, -hh), (hw, hh), (-hw, hh)]
    xs = [cx + dx * ca - dy * sa for dx, dy in corners]
    ys = [cy + dx * sa + dy * ca for dx, dy in corners]
    return min(xs), min(ys), max(xs), max(ys)


def _aabb_overlap(a, b):
    """Return True if two AABBs (x0,y0,x1,y1) overlap."""
    return a[0] < b[2] and a[2] > b[0] and a[1] < b[3] and a[3] > b[1]


def _draw_dim_line(shape, page, pt1, pt2, label, color,
                   fs=DIM_LABEL_FS, offset=DIM_LABEL_OFFSET,
                   along_frac=0.5, avoid_aabbs=()):
    """Draw a dimension line between pt1 and pt2 with a rotated label.

    along_frac: label position along the line (0=pt1, 0.5=midpoint, 1=pt2).
    avoid_aabbs: iterable of (x0,y0,x1,y1) AABBs to avoid with the label.
      If the default along_frac causes overlap, the label is shifted along the
      line in steps of 0.05 (trying both directions) until clear or exhausted.
    """
    shape.draw_line(fitz.Point(*pt1), fitz.Point(*pt2))
    shape.finish(color=color, width=0.3)
    dx = pt2[0] - pt1[0]
    dy = pt2[1] - pt1[1]
    length = math.hypot(dx, dy) or 1.0
    deg = math.degrees(math.atan2(dy, dx))
    tw = fitz.get_text_length(label, fontname="helv", fontsize=fs)

    def _label_center(frac):
        cx = pt1[0] + frac * dx + offset * dy / length
        cy = pt1[1] + frac * dy - offset * dx / length
        return cx, cy

    def _overlaps(frac):
        cx, cy = _label_center(frac)
        aabb = _label_aabb(cx, cy, tw, fs, -deg - 180)
        return any(_aabb_overlap(aabb, b) for b in avoid_aabbs)

    frac = along_frac
    if avoid_aabbs and _overlaps(frac):
        for step in (i * 0.05 for i in range(1, 18)):
            for candidate in (along_frac - step, along_frac + step):
                candidate = max(0.05, min(0.95, candidate))
                if not _overlaps(candidate):
                    frac = candidate
                    break
            else:
                continue
            break

    cx, cy = _label_center(frac)
    page.insert_text(
        fitz.Point(cx - tw / 2.0, cy + fs / 3.0),
        label, fontname="helv", fontsize=fs, color=color,
        morph=(fitz.Point(cx, cy), fitz.Matrix(-deg - 180)))


def _draw_setback_label(page, pt_pdf, line_p1, line_p2, value_ft, color,
                        fs=SETBACK_LABEL_FS):
    """Draw a label midway between pt and its perpendicular projection onto a line."""
    ldx = line_p2[0] - line_p1[0]
    ldy = line_p2[1] - line_p1[1]
    llen_sq = ldx * ldx + ldy * ldy
    t_proj = ((pt_pdf[0] - line_p1[0]) * ldx +
              (pt_pdf[1] - line_p1[1]) * ldy) / llen_sq
    proj_x = line_p1[0] + t_proj * ldx
    proj_y = line_p1[1] + t_proj * ldy
    cap_x = (pt_pdf[0] + proj_x) / 2.0
    cap_y = (pt_pdf[1] + proj_y) / 2.0
    perp_deg = math.degrees(math.atan2(proj_y - pt_pdf[1], proj_x - pt_pdf[0]))
    text = f"{value_ft:.1f}'"
    tw = fitz.get_text_length(text, fontname="helv", fontsize=fs)
    page.insert_text(
        fitz.Point(cap_x - tw / 2.0, cap_y + fs / 3.0),
        text, fontname="helv", fontsize=fs, color=color,
        morph=(fitz.Point(cap_x, cap_y), fitz.Matrix(-perp_deg)))


def render_site_plan(sp, corners=True):
    """Render base site plan PDF overlay. Returns unsaved fitz.Document."""
    src = fitz.open(os.path.join(os.path.dirname(__file__), "site_survey.pdf"))
    doc = fitz.open()
    doc.insert_pdf(src)
    src.close()
    page = doc[0]

    SCALE = sp.SCALE
    pts = sp.pts
    building_to_pdf = sp.building_to_pdf

    # --- Draw building outline ---
    shape = page.new_shape()
    x0, y0 = building_to_pdf(*sp.draw_poly[0])
    shape.draw_line(fitz.Point(x0, y0), fitz.Point(x0, y0))
    for pt in sp.draw_poly[1:]:
        x1, y1 = building_to_pdf(*pt)
        shape.draw_line(fitz.Point(x0, y0), fitz.Point(x1, y1))
        x0, y0 = x1, y1
    x1, y1 = building_to_pdf(*sp.draw_poly[0])
    shape.draw_line(fitz.Point(x0, y0), fitz.Point(x1, y1))
    shape.finish(color=COLOR_PROPOSED, width=OUTLINE_STROKE_W)

    # --- Parcel corner circles (2' radius) ---
    if corners:
        _corner_r = 2.0 * SCALE
        for cx, cy in (CORNER_NW, CORNER_NE, CORNER_SE, CORNER_SW):
            shape.draw_circle(fitz.Point(cx, cy), _corner_r)
        shape.finish(color=(1, 0, 0), width=0.5, fill=None,
                     stroke_opacity=0.4)

    # --- "PROPOSED 950SF MAX ADU" label (computed early for overlap avoidance) ---
    fc_pdf = building_to_pdf(*pts["FC"])
    _char_w = fitz.get_text_length("M", fontname="helv", fontsize=BLDG_LABEL_FS)
    label_pdf = (fc_pdf[0] + 0.25 * _char_w, fc_pdf[1] + 5.0 * SCALE - 2.3 * BLDG_LABEL_FS)
    label_lines = ["PROPOSED", "950SF MAX", "ADU"]
    label_lh = BLDG_LABEL_FS * 1.15
    block_h = label_lh * len(label_lines)
    start_y = label_pdf[1] - block_h / 2.0 + BLDG_LABEL_FS
    _max_lw = max(fitz.get_text_length(l, fontname="helv", fontsize=BLDG_LABEL_FS)
                  for l in label_lines)
    _bldg_label_aabb = (
        label_pdf[0] - _max_lw / 2.0, start_y - BLDG_LABEL_FS,
        label_pdf[0] + _max_lw / 2.0, start_y + (len(label_lines) - 1) * label_lh,
    )
    for i, line in enumerate(label_lines):
        lw = fitz.get_text_length(line, fontname="helv", fontsize=BLDG_LABEL_FS)
        page.insert_text(
            fitz.Point(label_pdf[0] - lw / 2.0, start_y + i * label_lh),
            line, fontname="helv", fontsize=BLDG_LABEL_FS, color=COLOR_PROPOSED)

    # --- E-W external dimension line (_site_se_pt → west ref) ---
    f15 = pts["_site_se_pt"]
    f15_pdf = building_to_pdf(*f15)
    foot_pdf = building_to_pdf(pts["_site_w_pt"][0], f15[1])
    _draw_dim_line(shape, page, f15_pdf, foot_pdf,
                   f"{sp.ew_dim_ft:.1f}'", COLOR_PROPOSED,
                   avoid_aabbs=(_bldg_label_aabb,))

    # --- N-S external dimension line (_site_s_pt → north ref, same easting) ---
    s_pt = pts["_site_s_pt"]
    s_pdf = building_to_pdf(*s_pt)
    n_pdf = building_to_pdf(s_pt[0], pts["_site_n_pt"][1])
    _draw_dim_line(shape, page, s_pdf, n_pdf,
                   f"{sp.ns_dim_ft:.1f}'", COLOR_PROPOSED,
                   avoid_aabbs=(_bldg_label_aabb,))

    # --- Arc dim: midpoint of F02-F03 arc to F24b-F26 surface ---
    if "_arc_dim_start" in pts and "_arc_dim_end" in pts:
        _arc_a_pdf = building_to_pdf(*pts["_arc_dim_end"])    # isect (south)
        _arc_b_pdf = building_to_pdf(*pts["_arc_dim_start"])  # arc_mid (north)
        _draw_dim_line(shape, page, _arc_a_pdf, _arc_b_pdf,
                       f"{sp.arc_dim_ft:.1f}'", COLOR_PROPOSED,
                       avoid_aabbs=(_bldg_label_aabb,))

    # --- Setback caption (from 216.73' line) ---
    f16_pdf = building_to_pdf(*pts["_site_sf_start"])
    f17_pdf = building_to_pdf(*pts["_site_sf_end"])
    mid_f16f17 = ((f16_pdf[0] + f17_pdf[0]) / 2.0,
                  (f16_pdf[1] + f17_pdf[1]) / 2.0)
    _draw_setback_label(page, mid_f16f17, LINE_TOP, LINE_BOT,
                        sp.min_setback_216, COLOR_PROPOSED)

    # --- Min setback from 275.08' line caption ---
    _draw_setback_label(page, sp.min275_pdf, BOT_LEFT, LINE_BOT,
                        sp.min_setback_275, COLOR_PROPOSED)

    # --- Distance from residence corner to closest F point ---
    _res_pt = sp.f_series_pdf[sp.residence_closest]
    _res_mid_x = (RESIDENCE_LR[0] + _res_pt[0]) / 2.0
    _res_mid_y = (RESIDENCE_LR[1] + _res_pt[1]) / 2.0
    _res_deg = math.degrees(math.atan2(
        _res_pt[1] - RESIDENCE_LR[1], _res_pt[0] - RESIDENCE_LR[0]))
    _res_text = f"{sp.residence_dist_ft:.1f}'"
    _res_fs = 9.0
    _res_tw = fitz.get_text_length(_res_text, fontname="helv", fontsize=_res_fs)
    page.insert_text(
        fitz.Point(_res_mid_x - _res_tw / 2.0, _res_mid_y + _res_fs / 3.0),
        _res_text, fontname="helv", fontsize=_res_fs, color=COLOR_PROPOSED,
        morph=(fitz.Point(_res_mid_x, _res_mid_y), fitz.Matrix(-_res_deg)))

    # --- "FRONT ↑" annotation above 251.53' line ---
    _mid_251_x = (TL_251[0] + LINE_TOP[0]) / 2.0
    _mid_251_y = (TL_251[1] + LINE_TOP[1]) / 2.0

    front_fs = 16.8
    front_gap = 4.0 + 4.0 * SCALE

    front_text = "FRONT"
    front_tw = fitz.get_text_length(front_text, fontname="helv", fontsize=front_fs)

    _arr_space = front_fs * 0.3
    _arr_h = front_fs * 0.65
    _arr_hw = front_fs * 0.15
    _arr_ah = front_fs * 0.2

    _total_w = front_tw + _arr_space

    _front_x = _mid_251_x - _total_w / 2.0
    _front_y = _mid_251_y - front_gap

    page.insert_text(
        fitz.Point(_front_x, _front_y),
        front_text, fontname="helv", fontsize=front_fs, color=(0, 0, 0))

    # Draw up-arrow next to "FRONT"
    _arr_cx = _front_x + front_tw + _arr_space
    _arr_bot = _front_y
    _arr_top = _front_y - _arr_h

    shape.draw_line(fitz.Point(_arr_cx, _arr_bot), fitz.Point(_arr_cx, _arr_top))
    shape.draw_line(fitz.Point(_arr_cx, _arr_top),
                    fitz.Point(_arr_cx - _arr_hw, _arr_top + _arr_ah))
    shape.draw_line(fitz.Point(_arr_cx, _arr_top),
                    fitz.Point(_arr_cx + _arr_hw, _arr_top + _arr_ah))
    shape.finish(color=(0, 0, 0), width=0.8)

    # --- Git describe / timestamp caption (right of GRAPHIC SCALE) ---
    _cap_x = 680.0
    _cap_y = 600.0
    _cap_fs = 5.5
    _cap_now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    _cap_git = git_describe()
    _cap_text = f"Generated {_cap_now}  {_cap_git}"
    _cap_tw = fitz.get_text_length(_cap_text, fontname="helv", fontsize=_cap_fs)
    page.insert_text(
        fitz.Point(_cap_x - _cap_tw / 2.0, _cap_y),
        _cap_text, fontname="helv", fontsize=_cap_fs, color=(0.4, 0.4, 0.4))

    shape.commit()
    return doc


def render_site_plan_df(doc, sp):
    """Add drainfield annotations to an existing site plan document."""
    page = doc[0]
    SCALE = sp.SCALE
    building_to_pdf = sp.building_to_pdf

    # Drainfield: dashed rounded rectangle, 30'W x 8'H, 2' corner radius
    _df_line_x = 661.5
    _df_line_y = 316.0

    _df_w = 25.0 * SCALE
    _df_h = 10.0 * SCALE
    _df_r = 2.0 * SCALE
    _df_left = _df_line_x + 8.0 * SCALE
    _df_top = _df_line_y - _df_h / 2.0
    _df_right = _df_left + _df_w
    _df_bot = _df_line_y + _df_h / 2.0

    _df_fs = 7.2
    _df_lh = _df_fs * 1.15
    _df_lines_existing = ["EXISTING", "DRAINFIELD"]
    _n_arc = 8

    def _draw_drainfield(left, top, right, bot, r, lines=None, angle_deg=0,
                         color=(0, 0, 0)):
        if lines is None:
            lines = _df_lines_existing
        cx = (left + right) / 2.0
        cy = (top + bot) / 2.0
        hw = (right - left) / 2.0
        hh = (bot - top) / 2.0

        rel = []
        rel.append((-hw + r, -hh))
        rel.append((hw - r, -hh))
        for i in range(1, _n_arc + 1):
            a = -math.pi / 2 + math.pi / 2 * i / _n_arc
            rel.append((hw - r + r * math.cos(a), -hh + r + r * math.sin(a)))
        rel.append((hw, hh - r))
        for i in range(1, _n_arc + 1):
            a = math.pi / 2 * i / _n_arc
            rel.append((hw - r + r * math.cos(a), hh - r + r * math.sin(a)))
        rel.append((-hw + r, hh))
        for i in range(1, _n_arc + 1):
            a = math.pi / 2 + math.pi / 2 * i / _n_arc
            rel.append((-hw + r + r * math.cos(a), hh - r + r * math.sin(a)))
        rel.append((-hw, -hh + r))
        for i in range(1, _n_arc + 1):
            a = math.pi + math.pi / 2 * i / _n_arc
            rel.append((-hw + r + r * math.cos(a), -hh + r + r * math.sin(a)))

        rad = math.radians(angle_deg)
        cos_r = math.cos(rad)
        sin_r = math.sin(rad)
        pts = []
        for rx, ry in rel:
            pts.append(fitz.Point(cx + rx * cos_r - ry * sin_r,
                                  cy + rx * sin_r + ry * cos_r))

        page.draw_polyline(pts, color=color, width=0.8,
                           dashes="[4 3] 0", closePath=True)

        block_h = _df_lh * len(lines)
        start_y = cy - block_h / 2.0 + _df_fs
        morph = None
        if angle_deg != 0:
            morph = (fitz.Point(cx, cy), fitz.Matrix(-angle_deg))
        for li, lt in enumerate(lines):
            tw = fitz.get_text_length(lt, fontname="helv", fontsize=_df_fs)
            page.insert_text(
                fitz.Point(cx - tw / 2.0, start_y + li * _df_lh),
                lt, fontname="helv", fontsize=_df_fs, color=color,
                morph=morph)

    # Right drainfield
    _draw_drainfield(_df_left, _df_top, _df_right, _df_bot, _df_r)

    # Left drainfield
    _res_left = 534.0
    _res_bl_y = 348.0
    _df2_right = _res_left - 8.0 * SCALE
    _df2_left = _df2_right - _df_w
    _df2_cy = _res_bl_y - 1.0 * SCALE
    _df2_top = _df2_cy - _df_h / 2.0
    _df2_bot = _df2_cy + _df_h / 2.0
    _draw_drainfield(_df2_left, _df2_top, _df2_right, _df2_bot, _df_r)

    # New drainfield: P3-near end is 8' from P3, aligned along P3→LINE_BOT
    p3_pdf = sp.p_series_pdf["P3"]
    _ndf_dx = LINE_BOT[0] - p3_pdf[0]
    _ndf_dy = LINE_BOT[1] - p3_pdf[1]
    _ndf_len = math.hypot(_ndf_dx, _ndf_dy)
    _ndf_ux = _ndf_dx / _ndf_len
    _ndf_uy = _ndf_dy / _ndf_len
    _ndf_cx = p3_pdf[0] + (6.0 * SCALE + _df_w / 2.0) * _ndf_ux
    _ndf_cy = p3_pdf[1] + (6.0 * SCALE + _df_w / 2.0) * _ndf_uy
    _ndf_angle = math.degrees(math.atan2(_ndf_dy, _ndf_dx))
    _ndf_left = _ndf_cx - _df_w / 2.0
    _ndf_top = _ndf_cy - _df_h / 2.0
    _ndf_right = _ndf_cx + _df_w / 2.0
    _ndf_bot = _ndf_cy + _df_h / 2.0
    _draw_drainfield(_ndf_left, _ndf_top, _ndf_right, _ndf_bot, _df_r,
                     lines=["NEW", "DRAINFIELD"], angle_deg=_ndf_angle,
                     color=COLOR_PROPOSED)

    return doc


COLOR_OUTER_PATH = (0, 0.4, 0)  # dark green for survey outer path



def render_site_plan_fs(doc, sp):
    """Add outer survey path (P-series traverse) overlay to an existing site plan."""
    page = doc[0]
    building_to_pdf = sp.building_to_pdf

    # Outer path segments (same definition as survey/gen_path_svg.py)
    R1, R2, R3 = 10.0, 12.5, 11.0
    outer_segs = [
        LineSeg("POB", "P2"), LineSeg("P2", "P3"), LineSeg("P3", "T3"),
        ArcSeg("T3", "PX", "TC3", R3, "CW", 60),
        LineSeg("PX", "P4"), LineSeg("P4", "P5"), LineSeg("P5", "PT1"),
        ArcSeg("PT1", "T1", "TC1", R1, "CW", 60),
        ArcSeg("T1", "PA", "TC1", R1, "CW", 60),
        ArcSeg("PA", "T2", "TC2", R2, "CW", 60),
        LineSeg("T2", "POB"),
    ]

    # Building coords → polygon → PDF coords (no angular correction)
    corrected_pts = dict(sp.pts)
    poly_bld = path_polygon(outer_segs, corrected_pts)
    poly_pdf = [fitz.Point(*building_to_pdf(*p)) for p in poly_bld]

    # Draw closed polyline: 40% opaque dark green
    page.draw_polyline(poly_pdf, color=COLOR_OUTER_PATH, width=1.2,
                       closePath=True, stroke_opacity=0.4)

    return doc


def main():
    sp = build_site_plan_data()

    # --- Base site_plan.pdf (includes corner circles) ---
    doc = render_site_plan(sp)
    out_path = os.path.join(os.path.dirname(__file__), "site_plan.pdf")
    doc.save(out_path)
    doc.close()
    print(f"Written to {out_path}")

    # --- site_plan_df.pdf (fresh render without corner circles) ---
    doc_df = render_site_plan(sp, corners=False)
    render_site_plan_df(doc_df, sp)
    df_path = os.path.join(os.path.dirname(__file__), "site_plan_df.pdf")
    doc_df.save(df_path)
    doc_df.close()
    print(f"Written to {df_path}")

    # --- site_plan_fs.pdf (df content + outer survey path) ---
    doc_fs = render_site_plan(sp, corners=False)
    render_site_plan_df(doc_fs, sp)
    render_site_plan_fs(doc_fs, sp)
    fs_path = os.path.join(os.path.dirname(__file__), "site_plan_fs.pdf")
    doc_fs.save(fs_path)
    doc_fs.close()
    print(f"Written to {fs_path}")

    print(f"Rotation: {sp.rotation_deg:.1f}\u00b0 CCW")
    print(f"F15 PDF position: ({sp.f15_pdf[0]:.1f}, {sp.f15_pdf[1]:.1f})")


if __name__ == "__main__":
    main()
