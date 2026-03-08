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


def apply_overrides_to_poly(inner_poly, inner_segs, pts, overrides):
    """Apply DB-driven inner wall overrides to an inner_poly list (in-place).

    Processes overrides in descending seg_index order so that splice
    index positions remain valid for earlier segments.

    Parameters:
        inner_poly  — mutable list of (E, N) tuples
        inner_segs  — list of LineSeg/ArcSeg (for bearing computation)
        pts         — point dict {name: (E, N)}
        overrides   — dict {seg_index: [sub-segment dicts]}
    """
    for seg_idx in sorted(overrides.keys(), reverse=True):
        chain = overrides[seg_idx]
        seg = inner_segs[seg_idx]
        start_pt = pts[seg.start]
        end_pt = pts[seg.end]
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
        return 10.0 / 12.0 + (wall_outer - 8.0 / 12.0)
    return constants_dict.get(name)


def _build_outline_segs_from_chain(chain):
    """Build outline_segs from solved chain (matching geometry.py rotation).

    Returns list of LineSeg/ArcSeg in outline convention (F1->F2 first).
    """
    point_names = [seg.end_name for seg in chain]
    start_names = ["F2"] + point_names[:-1]

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
                 overrides=None):
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
        self._compute_shell_geometry(constants_dict, overrides=overrides)

        # Inner area computed after F8-F9 replacement (matches build_floorplan_data)
        self.inner_area = poly_area(self.inner_poly)

        # Roof geometry (Phase 15-C)
        self._compute_roof_geometry()

    def _compute_shell_geometry(self, constants_dict, overrides=None):
        """Compute S/G-series shell paths, F8-F9 polylines, and wall sections.

        overrides — dict {seg_index: [sub-segment dicts]} from DB, or None
                    to use the hardcoded f8f9_corner_polyline fallback.
        """
        shell_t = constants_dict.get("SHELL_THICKNESS", 2.0 / 12.0)
        air_gap = constants_dict.get("AIR_GAP",
                                     self.wall_t - 2 * shell_t)
        opening_r = constants_dict.get("OPENING_INSIDE_RADIUS", 1.5 / 12.0)
        f8f9_inner_r = opening_r + shell_t

        # Apply inner wall overrides to inner_poly (W-series)
        if overrides:
            self._apply_inner_wall_overrides(overrides)
        else:
            # Legacy fallback: hardcoded F8-F9 splice
            self._apply_f8f9_hardcoded(f8f9_inner_r)

        # W-series F8-F9 polyline (for shell geometry, always from hardcoded)
        self.w_f8f9_poly = f8f9_corner_polyline(
            self.pts, self.wall_t, f8f9_inner_r)

        # S-series (inner face of outer shell)
        self.s_pts, self.s_segs = compute_inset_path(
            self.outline_segs, self.pts, self.radii, shell_t, "S")
        self.pts.update(self.s_pts)

        # G-series (outer face of inner shell)
        self.g_pts, self.g_segs = compute_inset_path(
            self.outline_segs, self.pts, self.radii,
            shell_t + air_gap, "G")
        self.pts.update(self.g_pts)

        # G-series F8-F9 polyline
        self.g_f8f9_poly = f8f9_corner_polyline(
            self.pts, shell_t + air_gap, opening_r)

        # Openings (parametric on outline segments) for wall section enumeration
        from floorplan.layout import compute_interior_layout
        self.layout = compute_interior_layout(self.pts, self.inner_poly)
        outer_openings = compute_outer_openings(self.pts, self.layout)
        self.openings = outer_to_wall_openings(
            outer_openings, self.outline_segs, self.pts)

        # Wall sections
        self.wall_sections = enumerate_wall_sections(
            self.openings, self.outline_segs)

    def _apply_f8f9_hardcoded(self, f8f9_inner_r):
        """Legacy hardcoded F8-F9 straight-arc-straight splice."""
        poly = f8f9_corner_polyline(self.pts, self.wall_t, f8f9_inner_r)
        _splice_poly(self.inner_poly, self.pts["W8"], self.pts["W9"], poly)

    def _apply_inner_wall_overrides(self, overrides):
        """Apply DB-driven inner wall overrides to inner_poly."""
        apply_overrides_to_poly(
            self.inner_poly, self.inner_segs, self.pts, overrides)

    def _compute_roof_geometry(self):
        """Compute R-series roof geometry."""
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
        from app.outline_solver import db_rows_to_chain, solve_closure, walk_chain
        chain = db_rows_to_chain(chain_rows)
        R_a1 = _derive_constant(constants_dict, "CORNER_SW_R")

        solver_result = solve_closure(chain, R_a1)
        if not solver_result.valid:
            raise ValueError(
                f"Outline chain does not close: error={solver_result.closure_error:.6f}")

        chain = list(chain)
        chain[0] = chain[0]._replace(distance=solver_result.d_F2_F5)
        chain[-2] = chain[-2]._replace(distance=solver_result.d_F18_F1)
        chain[-1] = chain[-1]._replace(sweep=solver_result.sweep_closure)

        F2_E = constants_dict.get("F2_EASTING", -18.5)
        F2_N = constants_dict.get("F2_NORTHING", -13.5) + R_a1
        walk_result = walk_chain(chain, F2_E, F2_N)

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

    overrides — dict {seg_index: [sub-segment dicts]} of inner wall
                overrides.  If None and db_path is set, loads from DB.
    """
    pts, outline_segs, inner_segs, radii = compute_native_geometry(
        constants_dict, chain_rows=chain_rows, db_path=db_path)

    # Load overrides from DB if not explicitly provided
    if overrides is None and db_path is not None:
        from app.database import get_inner_wall_overrides
        overrides = get_inner_wall_overrides(db_path)

    return GeneratorData(pts, outline_segs, inner_segs, radii, constants_dict,
                         overrides=overrides)
