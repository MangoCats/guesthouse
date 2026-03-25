"""Generator data provider — bridge between DB-driven engine and SVG generators.

Provides GeneratorData, a typed container of native Python geometry objects
(Point tuples, LineSeg/ArcSeg namedtuples) that SVG generators can consume
directly.  This replaces direct imports from floorplan/ modules, allowing
generators to render from database-edited state.

Also contains the shared geometry pipeline (compute_native_geometry) used by
both compute_geometry() and GeneratorData, plus helpers for DB-driven traverse,
outline chain walking, and derived constant computation.
"""
import math

from shared.types import LineSeg, ArcSeg
from shared.geometry import (
    compute_inner_walls, path_polygon, poly_area,
    f8f9_corner_polyline, GEOM_EPS,
)
from shared.wall_shells import compute_inset_path, enumerate_wall_sections
from floorplan.openings import compute_outer_openings, outer_to_wall_openings
from floorplan.roof import compute_roof_geometry, roof_polyline


# ---------------------------------------------------------------------------
# Walk override chain (Phase 15½-B)
# ---------------------------------------------------------------------------

def walk_override_chain(chain, start_pt, start_bearing_deg):
    """Walk a parametric override chain, returning a polyline.

    Parameters:
        chain       — list of dicts with keys: seg_type, bearing, distance,
                       radius, sweep, n_pts
        start_pt    — (E, N) starting point
        start_bearing_deg — compass bearing in degrees at chain start

    Returns list of (E, N) points including start_pt and final endpoint.
    Bearings use compass convention (0=N, 90=E, 180=S, 270=W).
    """
    polyline = [start_pt]
    cur = start_pt
    cur_bearing = start_bearing_deg

    for sub in chain:
        seg_type = sub["seg_type"]
        if seg_type == "L":
            brg = sub["bearing"]
            dist = sub["distance"]
            brg_rad = math.radians(brg)
            dx = dist * math.sin(brg_rad)
            dy = dist * math.cos(brg_rad)
            cur = (cur[0] + dx, cur[1] + dy)
            polyline.append(cur)
            cur_bearing = brg
        else:
            # Arc: CW or CCW
            radius = sub["radius"]
            sweep_deg = sub["sweep"]
            n_pts = sub.get("n_pts", 20)
            sweep_rad = math.radians(sweep_deg)

            # Current direction vector from compass bearing
            dir_x = math.sin(math.radians(cur_bearing))
            dir_y = math.cos(math.radians(cur_bearing))

            if seg_type == "CCW":
                # Center is to the LEFT of travel direction
                cx = cur[0] - dir_y * radius
                cy = cur[1] + dir_x * radius
                entry_angle = math.atan2(cur[1] - cy, cur[0] - cx)
                for i in range(1, n_pts + 1):
                    a = entry_angle + i * sweep_rad / n_pts
                    polyline.append((cx + radius * math.cos(a),
                                     cy + radius * math.sin(a)))
            else:
                # CW: center is to the RIGHT of travel direction
                cx = cur[0] + dir_y * radius
                cy = cur[1] - dir_x * radius
                entry_angle = math.atan2(cur[1] - cy, cur[0] - cx)
                for i in range(1, n_pts + 1):
                    a = entry_angle - i * sweep_rad / n_pts
                    polyline.append((cx + radius * math.cos(a),
                                     cy + radius * math.sin(a)))

            cur = polyline[-1]
            # Update bearing: rotate by sweep
            if seg_type == "CCW":
                cur_bearing = (cur_bearing - sweep_deg) % 360
            else:
                cur_bearing = (cur_bearing + sweep_deg) % 360

    return polyline


# ---------------------------------------------------------------------------
# Inner poly splice helpers
# ---------------------------------------------------------------------------

def _splice_poly(poly, start_pt, end_pt, replacement):
    """Replace the sub-path from start_pt to end_pt in poly with replacement.

    Finds start_pt and end_pt by coordinate match (within GEOM_EPS),
    then replaces poly[start_idx:end_idx+1] with replacement.
    """
    start_idx = next(i for i, p in enumerate(poly)
                     if abs(p[0] - start_pt[0]) < GEOM_EPS
                     and abs(p[1] - start_pt[1]) < GEOM_EPS)
    end_idx = next(i for i, p in enumerate(poly)
                   if i > start_idx
                   and abs(p[0] - end_pt[0]) < GEOM_EPS
                   and abs(p[1] - end_pt[1]) < GEOM_EPS)
    poly[start_idx:end_idx + 1] = replacement


def _seg_start_bearing(seg, pts):
    """Compute compass bearing (degrees) at the start of an inner_seg.

    For LineSeg: bearing from start→end.
    For ArcSeg: CW tangent at start point (right normal of radius vector).
    """
    if isinstance(seg, LineSeg):
        s, e = pts[seg.start], pts[seg.end]
        dx, dy = e[0] - s[0], e[1] - s[1]
    else:
        # Arc: tangent at start = perpendicular to radius
        c = pts[seg.center]
        s = pts[seg.start]
        rx, ry = s[0] - c[0], s[1] - c[1]
        if seg.direction == "CW":
            dx, dy = ry, -rx  # CW tangent = right normal of radius
        else:
            dx, dy = -ry, rx  # CCW tangent = left normal of radius
    return math.degrees(math.atan2(dx, dy)) % 360


def compute_default_override(seg_index, inner_segs, pts, constants_dict):
    """Compute the default parametric override chain for an inner segment.

    For line segments: single L sub-segment with bearing and distance.
    For arc segments: single CW/CCW sub-segment with radius and sweep.
    For seg 5 (W8-W9 concave corner): L-CCW-L straight-arc-straight chain
    computed from the F-series geometry (same logic as _seed_inner_wall_overrides).

    Returns list of sub-segment dicts, or empty list if seg_index is invalid.
    """
    if seg_index < 0 or seg_index >= len(inner_segs):
        return []

    seg = inner_segs[seg_index]
    start_pt = pts[seg.start]
    end_pt = pts[seg.end]

    # Special case: W8-W9 concave corner (seg 5) — straight-arc-straight
    if seg.start == "W8" and seg.end == "W9":
        return _compute_w8w9_default_chain(pts, constants_dict)

    if isinstance(seg, LineSeg):
        dx = end_pt[0] - start_pt[0]
        dy = end_pt[1] - start_pt[1]
        dist = math.hypot(dx, dy)
        bearing = math.degrees(math.atan2(dx, dy)) % 360
        return [{"seg_type": "L", "bearing": round(bearing, 4),
                 "distance": round(dist, 6), "n_pts": 20}]
    else:
        # Arc segment
        center = pts[seg.center]
        r = seg.radius
        a_start = math.atan2(start_pt[1] - center[1], start_pt[0] - center[0])
        a_end = math.atan2(end_pt[1] - center[1], end_pt[0] - center[0])
        if seg.direction == "CW":
            sweep = (a_start - a_end) % (2 * math.pi)
        else:
            sweep = (a_end - a_start) % (2 * math.pi)
        return [{"seg_type": seg.direction,
                 "radius": round(r, 6),
                 "sweep": round(math.degrees(sweep), 4),
                 "n_pts": seg.n_pts}]


def _walk_chain_exit(chain, start_bearing_deg):
    """Compute the exit bearing after walking a parametric chain.

    Returns the compass bearing (degrees) at the end of the chain.
    """
    cur_bearing = start_bearing_deg
    for sub in chain:
        seg_type = sub["seg_type"]
        if seg_type == "L":
            cur_bearing = sub["bearing"]
        elif seg_type == "CCW":
            cur_bearing = (cur_bearing - sub["sweep"]) % 360
        else:  # CW
            cur_bearing = (cur_bearing + sub["sweep"]) % 360
    return cur_bearing


def validate_override_endpoint(chain, start_pt, start_bearing_deg,
                               expected_end_pt, expected_exit_bearing_deg,
                               pos_tol=0.01, brg_tol=0.5):
    """Validate that an override chain arrives at the expected endpoint.

    Returns dict with ``ok`` (bool) and ``warnings`` (list of strings).
    pos_tol is in feet, brg_tol is in degrees.
    """
    warnings = []
    poly = walk_override_chain(chain, start_pt, start_bearing_deg)
    actual_end = poly[-1]
    dist_err = math.hypot(actual_end[0] - expected_end_pt[0],
                          actual_end[1] - expected_end_pt[1])
    if dist_err > pos_tol:
        warnings.append(
            f"Endpoint error: {dist_err:.4f} ft from target "
            f"(tolerance {pos_tol} ft)")
    exit_brg = _walk_chain_exit(chain, start_bearing_deg)
    brg_err = abs((exit_brg - expected_exit_bearing_deg + 180) % 360 - 180)
    if brg_err > brg_tol:
        warnings.append(
            f"Exit bearing error: {brg_err:.2f}\u00b0 from target "
            f"(tolerance {brg_tol}\u00b0)")
    return {"ok": len(warnings) == 0, "warnings": warnings}


def compute_default_span_override(seg_index, span_end, inner_segs, pts,
                                   constants_dict):
    """Compute the default parametric chain for a multi-segment span.

    Concatenates the default chain for each inner_seg in the range
    [seg_index, span_end] inclusive.  Returns list of sub-segment dicts.
    """
    if seg_index < 0 or span_end >= len(inner_segs) or span_end < seg_index:
        return []
    chain = []
    for idx in range(seg_index, span_end + 1):
        chain.extend(compute_default_override(idx, inner_segs, pts,
                                              constants_dict))
    return chain


def _compute_w8w9_default_chain(pts, constants_dict):
    """Compute the W8-W9 straight-arc-straight chain from current geometry.

    Same logic as _seed_inner_wall_overrides in database.py, but runs
    from live geometry instead of import-time constants.
    """
    wall_outer = constants_dict.get("WALL_OUTER", 8.0 / 12.0)
    shell_t = constants_dict.get("SHELL_THICKNESS", 2.0 / 12.0)
    opening_r = constants_dict.get("OPENING_INSIDE_RADIUS", 1.5 / 12.0)
    R_turn = opening_r + shell_t

    F8, C7 = pts["F8"], pts["C7"]
    F9, F10 = pts["F9"], pts["F10"]

    # CW traversal directions
    r8x, r8y = F8[0] - C7[0], F8[1] - C7[1]
    r8_len = math.hypot(r8x, r8y)
    dir_f8 = (r8y / r8_len, -r8x / r8_len)

    d9x, d9y = F10[0] - F9[0], F10[1] - F9[1]
    d9_len = math.hypot(d9x, d9y)
    dir_f9 = (d9x / d9_len, d9y / d9_len)

    brg_f8 = math.degrees(math.atan2(dir_f8[0], dir_f8[1]))
    brg_f9 = math.degrees(math.atan2(dir_f9[0], dir_f9[1]))

    # Inset points
    ins_f8 = (dir_f8[1], -dir_f8[0])
    ins_f9 = (dir_f9[1], -dir_f9[0])
    W8 = (F8[0] + wall_outer * ins_f8[0], F8[1] + wall_outer * ins_f8[1])
    W9 = (F9[0] + wall_outer * ins_f9[0], F9[1] + wall_outer * ins_f9[1])

    # Arc center and tangent distances
    left_f8 = (-dir_f8[1], dir_f8[0])
    left_f9 = (-dir_f9[1], dir_f9[0])
    P1 = (W8[0] + R_turn * left_f8[0], W8[1] + R_turn * left_f8[1])
    P2 = (W9[0] + R_turn * left_f9[0], W9[1] + R_turn * left_f9[1])
    dp = (P2[0] - P1[0], P2[1] - P1[1])
    cross = dir_f8[0] * dir_f9[1] - dir_f8[1] * dir_f9[0]
    t = (dp[0] * dir_f9[1] - dp[1] * dir_f9[0]) / cross
    arc_tangent1 = (W8[0] + t * dir_f8[0], W8[1] + t * dir_f8[1])

    t2_dp = (P1[0] - P2[0], P1[1] - P2[1])
    t2_cross = dir_f9[0] * dir_f8[1] - dir_f9[1] * dir_f8[0]
    t2 = (t2_dp[0] * dir_f8[1] - t2_dp[1] * dir_f8[0]) / t2_cross
    arc_tangent2 = (W9[0] + t2 * dir_f9[0], W9[1] + t2 * dir_f9[1])

    dist_start = math.hypot(arc_tangent1[0] - W8[0], arc_tangent1[1] - W8[1])
    dist_end = math.hypot(W9[0] - arc_tangent2[0], W9[1] - arc_tangent2[1])

    entry = math.atan2(-left_f8[1], -left_f8[0])
    exit_ = math.atan2(-left_f9[1], -left_f9[0])
    sweep_deg = math.degrees((exit_ - entry) % (2 * math.pi))

    return [
        {"seg_type": "L", "bearing": round(brg_f8, 4),
         "distance": round(dist_start, 6), "n_pts": 20},
        {"seg_type": "CCW", "radius": round(R_turn, 6),
         "sweep": round(sweep_deg, 4), "n_pts": 20},
        {"seg_type": "L", "bearing": round(brg_f9, 4),
         "distance": round(dist_end, 6), "n_pts": 20},
    ]


def apply_overrides_to_poly(inner_poly, inner_segs, pts, overrides):
    """Apply DB-driven inner wall overrides to an inner_poly list (in-place).

    Processes overrides in descending seg_index order so that splice
    index positions remain valid for earlier segments.

    Supports multi-segment spans: when a chain's first entry has a non-None
    ``span_end``, the override replaces inner_segs[seg_index] through
    inner_segs[span_end] inclusive.

    Parameters:
        inner_poly  — mutable list of (E, N) tuples
        inner_segs  — list of LineSeg/ArcSeg (for bearing computation)
        pts         — point dict {name: (E, N)}
        overrides   — dict {seg_index: [sub-segment dicts]}
    """
    for seg_idx in sorted(overrides.keys(), reverse=True):
        if seg_idx >= len(inner_segs):
            continue  # override references a segment that no longer exists
        chain = overrides[seg_idx]
        seg = inner_segs[seg_idx]
        if seg.start not in pts or seg.end not in pts:
            continue  # referenced points absent (e.g. non-F-series chain)
        start_pt = pts[seg.start]
        # Determine span end
        span_end = chain[0].get("span_end") if chain else None
        end_seg_idx = span_end if span_end is not None else seg_idx
        if end_seg_idx >= len(inner_segs):
            continue
        end_seg = inner_segs[end_seg_idx]
        if end_seg.end not in pts:
            continue
        end_pt = pts[end_seg.end]
        start_bearing = _seg_start_bearing(seg, pts)
        poly = walk_override_chain(chain, start_pt, start_bearing)
        _splice_poly(inner_poly, start_pt, end_pt, poly)


# ---------------------------------------------------------------------------
# Helpers (moved from engine.py to avoid circular imports)
# ---------------------------------------------------------------------------

def _derive_constant(constants_dict, name):
    """Compute a derived constant from the constants dict.

    Handles derived constants that depend on other constant values:
      WALL_EXTRA = WALL_OUTER - 8/12
      CORNER_SW_R = 10/12 + WALL_EXTRA
    """
    wall_outer = constants_dict.get("WALL_OUTER", 8.0 / 12.0)
    if name == "WALL_EXTRA":
        return wall_outer - 8.0 / 12.0
    if name == "CORNER_SW_R":
        if "CORNER_SW_R" in constants_dict:
            return constants_dict["CORNER_SW_R"]
        return 10.0 / 12.0 + (wall_outer - 8.0 / 12.0)
    return constants_dict.get(name)


def _build_outline_segs_from_chain(chain):
    """Build outline_segs from solved chain (matching geometry.py rotation).

    Returns list of LineSeg/ArcSeg in outline convention (F1->F2 first).
    """
    point_names = [seg.end_name for seg in chain]
    chain_start = chain[-1].end_name   # closure arc returns to chain start
    start_names = [chain_start] + point_names[:-1]

    segs = []
    for entry, start, end in zip(chain, start_names, point_names):
        if entry.seg_type == "L":
            segs.append(LineSeg(start, end))
        else:
            segs.append(ArcSeg(start, end, entry.center_name,
                               entry.radius, entry.seg_type, entry.n_pts))

    # Rotate: F1->F2 first (last entry becomes first)
    return segs[-1:] + segs[:-1]


def _compute_traverse_from_db(db_path=None):
    """Compute traverse from DB survey data (Phase 14-B).

    Same math as shared/survey.py:compute_traverse() but reads legs and
    config from the database instead of hardcoded values.
    """
    from app.database import get_survey_legs, get_survey_config

    legs = get_survey_legs(db_path)
    config = get_survey_config(db_path)

    # Walk legs from origin
    trav = [(0.0, 0.0)]
    for leg in legs:
        brg = leg["bearing_deg"] + leg["bearing_min"] / 60.0 + leg["bearing_sec"] / 3600.0
        dist_in = leg["distance_ft"] * 12 + leg["distance_inch"]
        brg_rad = math.radians(brg)
        dE = dist_in * math.sin(brg_rad)
        dN = dist_in * math.cos(brg_rad)
        last = trav[-1]
        trav.append((last[0] + dE, last[1] + dN))

    # Convert to feet, take first 5 points
    trav_ft = [(e / 12, n / 12) for e, n in trav[:5]]

    # Apply manual corrections from config
    p3_e_override = config.get("P3_EASTING_OVERRIDE", -19.1177)
    p2_p3_n_offset = config.get("P2_P3_NORTHING_OFFSET", 29.0)
    trav_ft[2] = (p3_e_override, trav_ft[3][1])
    trav_ft[1] = (trav_ft[2][0], trav_ft[2][1] + p2_p3_n_offset)

    # Shift from P3 origin to FC origin
    fc_e = config.get("FC_IN_P3_E", 18.5141152720)
    fc_n = config.get("FC_IN_P3_N", 13.3968094375)
    p3 = trav_ft[2]
    pts = {}
    labels = ["POB", "P2", "P3", "P4", "P5"]
    for i, label in enumerate(labels):
        pts[label] = (trav_ft[i][0] - p3[0] - fc_e,
                      trav_ft[i][1] - p3[1] - fc_n)

    return pts


def _db_openings_to_wall_openings(db_openings, outline_segs, pts):
    """Convert DB-driven outer opening dicts to WallOpening objects.

    Each db_opening has name, seg_start, seg_end, poly.  We find the
    matching outline segment by (seg_start, seg_end) and compute
    parametric positions from the first two polygon vertices.
    """
    from floorplan.openings import WallOpening

    seg_map = {(seg.start, seg.end): i for i, seg in enumerate(outline_segs)}
    result = []
    for o in db_openings:
        key = (o["seg_start"], o["seg_end"])
        idx = seg_map.get(key)
        if idx is None:
            continue  # skip openings on segments not in current chain
        seg = outline_segs[idx]
        poly = o.get("poly", [])
        if len(poly) < 2:
            continue
        # Compute parametric position along the segment
        A = pts[seg.start]
        B = pts[seg.end]
        dx = B[0] - A[0]
        dy = B[1] - A[1]
        if abs(dx) < 1e-12 and abs(dy) < 1e-12:
            continue
        if abs(dx) > abs(dy):
            t1 = (poly[0][0] - A[0]) / dx
            t2 = (poly[1][0] - A[0]) / dx
        else:
            t1 = (poly[0][1] - A[1]) / dy
            t2 = (poly[1][1] - A[1]) / dy
        result.append(WallOpening(o["name"], idx, min(t1, t2), max(t1, t2)))
    return result


# ---------------------------------------------------------------------------
# GeneratorData
# ---------------------------------------------------------------------------

class GeneratorData:
    """All geometry a generator needs, sourced from the database.

    Attributes — Outline (Phase 15-A):
        pts           — dict of named points {str: (E, N)} (F/W/C/P/TC series)
        outline_segs  — list of LineSeg/ArcSeg (CW traversal)
        outline_poly  — list of (E, N) tuples (outline polygon vertices)
        inner_segs    — list of LineSeg/ArcSeg (inner wall boundary)
        inner_poly    — list of (E, N) tuples (inner wall polygon)
        radii         — dict of arc radii {str: float}
        constants     — dict of all constants {name: value}
        wall_t        — outer wall thickness in feet
        outer_area    — outline polygon area in sq ft
        inner_area    — inner polygon area in sq ft

    Attributes — Shell (Phase 15-B):
        s_pts         — S-series shell boundary points
        s_segs        — S-series segments (SHELL_THICKNESS inset)
        g_pts         — G-series shell boundary points
        g_segs        — G-series segments (SHELL_THICKNESS + AIR_GAP inset)
        w_f8f9_poly   — W-series F8-F9 corner polyline
        g_f8f9_poly   — G-series F8-F9 corner polyline
        openings      — WallOpening list (parametric outer wall openings)
        wall_sections — enumerated wall sections with openings
        layout        — InteriorLayout namedtuple (rooms, appliances, furniture)

    Attributes — Roof (Phase 15-C):
        roof          — RoofGeometry namedtuple (pts, radii, centers, area)
        roof_poly     — list of (E, N) tuples (roof outline polygon)
    """

    def __init__(self, pts, outline_segs, inner_segs, radii, constants_dict,
                 overrides=None, db_openings=None):
        self.pts = pts
        self.outline_segs = outline_segs
        self.inner_segs = inner_segs
        self.radii = radii
        self.constants = constants_dict

        self.wall_t = constants_dict.get("WALL_OUTER", 8.0 / 12.0)
        self.outline_poly = path_polygon(outline_segs, pts)
        self.inner_poly = path_polygon(inner_segs, pts)
        self.outer_area = poly_area(self.outline_poly)

        # Shell geometry — modifies inner_poly with override splices
        self._compute_shell_geometry(constants_dict, overrides=overrides,
                                     db_openings=db_openings)

        # Inner area computed after F8-F9 replacement (matches build_floorplan_data)
        self.inner_area = poly_area(self.inner_poly)

        # Roof geometry (Phase 15-C)
        self._compute_roof_geometry()

    def _compute_shell_geometry(self, constants_dict, overrides=None,
                                db_openings=None):
        """Compute S/G-series shell paths, F8-F9 polylines, and wall sections.

        overrides   — dict {seg_index: [sub-segment dicts]} from DB, or None.
                      Inner wall geometry comes exclusively from DB overrides;
                      there is no hardcoded fallback.
        db_openings — list of DB-driven outer opening dicts with keys
                      name, seg_start, seg_end, poly.  When provided, these
                      are used instead of hardcoded compute_outer_openings.
        """
        shell_t = constants_dict.get("SHELL_THICKNESS", 2.0 / 12.0)
        air_gap = constants_dict.get("AIR_GAP",
                                     self.wall_t - 2 * shell_t)
        opening_r = constants_dict.get("OPENING_INSIDE_RADIUS", 1.5 / 12.0)

        # Apply inner wall overrides (DB-only — no hardcoded fallback)
        _has_f8f9 = "F8" in self.pts and "W8" in self.pts and "W9" in self.pts
        if overrides:
            self._apply_inner_wall_overrides(overrides)

        # Find F8→W9 inner-seg index and its DB override chain
        f8f9_inner_idx = next((i for i, s in enumerate(self.inner_segs)
                               if s.start == "W8" and s.end == "W9"), None)
        _f8f9_chain = (overrides.get(f8f9_inner_idx)
                       if overrides and f8f9_inner_idx is not None else None)
        _render_f8f9 = bool(_f8f9_chain) and _has_f8f9

        # W/G-series F8-F9 polylines: only present when DB override exists
        if _render_f8f9:
            # W-series: walk the DB override chain from W8
            seg = self.inner_segs[f8f9_inner_idx]
            start_bearing = _seg_start_bearing(seg, self.pts)
            self.w_f8f9_poly = walk_override_chain(
                _f8f9_chain, self.pts["W8"], start_bearing)
            # G-series: parallel path at shell gap offset (from DB constants)
            self.g_f8f9_poly = f8f9_corner_polyline(
                self.pts, shell_t + air_gap, opening_r)
        else:
            self.w_f8f9_poly = []
            self.g_f8f9_poly = []

        # S-series (inner face of outer shell)
        self.s_pts, self.s_segs = compute_inset_path(
            self.outline_segs, self.pts, self.radii, shell_t, "S")
        self.pts.update(self.s_pts)

        # G-series (outer face of inner shell)
        self.g_pts, self.g_segs = compute_inset_path(
            self.outline_segs, self.pts, self.radii,
            shell_t + air_gap, "G")
        self.pts.update(self.g_pts)

        # Layout — hardcoded seed values; used only on the standalone
        # subprocess path.  DB-driven callers should use gd.iw_polys instead.
        try:
            from floorplan.layout import compute_interior_layout
            self.layout = compute_interior_layout(self.pts, self.inner_poly)
        except Exception:
            self.layout = None

        # DB-driven IW polygons — populated by build_generator_data() when
        # db_path is available.  Span generators prefer these over layout.
        self.iw_polys = None

        # Openings (parametric on outline segments) for wall section enumeration
        if db_openings is not None:
            # DB-driven: convert DB opening dicts to WallOpening objects
            self.openings = _db_openings_to_wall_openings(
                db_openings, self.outline_segs, self.pts)
        else:
            # Legacy fallback: hardcoded openings (requires layout)
            outer_openings = compute_outer_openings(self.pts, self.layout)
            self.openings = outer_to_wall_openings(
                outer_openings, self.outline_segs, self.pts)

        # Wall sections
        self.wall_sections = enumerate_wall_sections(
            self.openings, self.outline_segs)

    def _apply_inner_wall_overrides(self, overrides):
        """Apply DB-driven inner wall overrides to inner_poly."""
        apply_overrides_to_poly(
            self.inner_poly, self.inner_segs, self.pts, overrides)

    def _compute_roof_geometry(self):
        """Compute R-series roof geometry (F-series chain only)."""
        if "R_a5" not in self.radii:
            self.roof = None
            self.roof_poly = []
            return
        self.roof = compute_roof_geometry(self.pts, self.radii)
        self.roof_poly = roof_polyline(self.roof)
        # Merge roof points into main pts dict
        self.pts.update(self.roof.pts)


# ---------------------------------------------------------------------------
# Shared geometry pipeline
# ---------------------------------------------------------------------------

def compute_native_geometry(constants_dict, chain_rows=None, db_path=None):
    """Compute outline + inner wall geometry as native Python objects.

    This is the shared computation used by both compute_geometry() (which
    serializes to JSON) and GeneratorData (which keeps native objects).

    Returns (pts, outline_segs, inner_segs, radii) tuple where:
        pts          — dict {str: (E, N)} of named points
        outline_segs — list of LineSeg/ArcSeg
        inner_segs   — list of LineSeg/ArcSeg
        radii        — dict {str: float} of arc radii
    """
    from floorplan.geometry import compute_outline_geometry, align_pts_to_f_series
    from shared.survey import compute_traverse, compute_three_arc, compute_inset

    # 1. Survey traverse
    if db_path is not None:
        trav_pts = _compute_traverse_from_db(db_path)
    else:
        trav_pts = compute_traverse()
    three_arc = compute_three_arc(trav_pts)
    inset_res = compute_inset(
        trav_pts, three_arc["R1"], three_arc["R2"], three_arc["R3"],
        three_arc["nE"], three_arc["nN"],
    )
    trav_pts.update(inset_res.pts_update)
    align_pts_to_f_series(trav_pts)

    # 2. F-series outline
    if chain_rows is not None:
        from app.outline_solver import (db_rows_to_chain,
                                        solve_closure_general, walk_chain,
                                        flex_specs_from_chain_rows,
                                        point_name_to_seq)
        from app.database import get_outline_anchor_pos, get_outline_anchor_pivot

        chain = db_rows_to_chain(chain_rows)
        n = len(chain)

        anchor_pos = get_outline_anchor_pos(db_path) if db_path else None
        if anchor_pos is not None:
            # Anchor-relative rotated walk — DB values already solved by pivot-aware solver
            anchor_name, _pivot_name = get_outline_anchor_pivot(db_path)
            anchor_pt_seq = point_name_to_seq(chain, anchor_name)
            a_start = (anchor_pt_seq + 1) % n
            start_E, start_N, start_brg = anchor_pos
            rotated = [chain[(a_start + i) % n] for i in range(n)]
            walk_result = walk_chain(rotated, start_E, start_N, start_brg)
        else:
            # Fallback: standard DB-order walk (default anchor = chain wrap point, brg=0)
            R_a1 = _derive_constant(constants_dict, "CORNER_SW_R")
            flex_specs = flex_specs_from_chain_rows(chain_rows)
            solver_result = solve_closure_general(chain, flex_specs)
            if not solver_result.valid:
                raise ValueError(
                    f"Outline chain does not close: error={solver_result.closure_error:.6f}")
            chain = list(chain)
            for seq, (param, value) in solver_result.solved_values.items():
                chain[seq] = chain[seq]._replace(**{param: value})
            chain_start = chain[-1].end_name
            start_E = constants_dict.get(f"{chain_start}_EASTING",
                                         constants_dict.get("F2_EASTING", -18.5))
            start_N = constants_dict.get(f"{chain_start}_NORTHING",
                                         constants_dict.get("F2_NORTHING", -13.5)) + R_a1
            walk_result = walk_chain(chain, start_E, start_N)

        fp_pts = walk_result.points
        radii = walk_result.radii
        outline_segs = _build_outline_segs_from_chain(chain)

        pts = dict(fp_pts)
        pts.update(trav_pts)
    else:
        geom = compute_outline_geometry()
        pts = dict(geom.fp_pts)
        pts.update(trav_pts)
        outline_segs = geom.outline_segs
        radii = geom.radii

    # 3. Inner walls (W-series)
    wall_t = constants_dict.get("WALL_OUTER", 8.0 / 12.0)
    inner_segs = compute_inner_walls(outline_segs, pts, wall_t, radii)

    return pts, outline_segs, inner_segs, radii


def build_generator_data(constants_dict, chain_rows=None, db_path=None,
                         overrides=None):
    """Build a GeneratorData object from DB state.

    This is the main entry point for SVG generators.  It computes all
    geometry from the database and returns a GeneratorData with native
    Python objects ready for rendering.

    When db_path is provided, calls compute_geometry() once (with chain_rows
    so F-series points are DB-driven) to populate:
      - db_openings: DB-driven outer opening positions for wall sections
      - gd.iw_polys: DB-driven interior wall polygons for span generators

    overrides — dict {seg_index: [sub-segment dicts]} of inner wall
                overrides.  If None and db_path is set, loads from DB.
    """
    pts, outline_segs, inner_segs, radii = compute_native_geometry(
        constants_dict, chain_rows=chain_rows, db_path=db_path)

    # Load overrides from DB if not explicitly provided
    if overrides is None and db_path is not None:
        from app.database import get_inner_wall_overrides
        overrides = get_inner_wall_overrides(db_path)

    # Compute DB-driven geometry once for outer openings + IW polys
    db_openings = None
    iw_polys = None
    if db_path is not None:
        try:
            from app.engine import compute_geometry
            geom = compute_geometry(constants_dict, variant="standard",
                                    chain_rows=chain_rows, db_path=db_path)
            db_openings = geom.get("outer_openings", [])
            iw_polys = {
                name: data["poly"]
                for name, data in geom.get("interior_walls", {}).items()
                if "poly" in data
            }
        except Exception:
            pass

    gd = GeneratorData(pts, outline_segs, inner_segs, radii, constants_dict,
                       overrides=overrides, db_openings=db_openings)
    if iw_polys:
        gd.iw_polys = iw_polys
    return gd
