"""Generator data provider — bridge between DB-driven engine and SVG generators.

Provides GeneratorData, a typed container of native Python geometry objects
(Point tuples, LineSeg/ArcSeg namedtuples) that SVG generators can consume
directly.  This replaces direct imports from floorplan/ modules, allowing
generators to render from database-edited state.

Phase 15-A: core outline + inner wall geometry
Phase 15-B: shell geometry (S/G-series), wall sections, U-turn polygons
Phase 15-C: roof geometry (R-series corners, roof polyline, roof area)
"""
import math

from shared.geometry import (
    compute_inner_walls, path_polygon, poly_area,
    f8f9_corner_polyline, GEOM_EPS,
)
from shared.wall_shells import (
    compute_inset_path, enumerate_wall_sections, uturn_polygon,
    build_section_outlines,
)
from floorplan.openings import (
    compute_outer_openings, compute_rough_openings, outer_to_wall_openings,
)
from floorplan.roof import compute_roof_geometry, roof_polyline


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

    Attributes — Roof (Phase 15-C):
        roof          — RoofGeometry namedtuple (pts, radii, centers, area)
        roof_poly     — list of (E, N) tuples (roof outline polygon)

    Attributes — Layout (from procedural modules, for identity with generators):
        layout        — InteriorLayout namedtuple (rooms, appliances, furniture)
    """

    def __init__(self, pts, outline_segs, inner_segs, radii, constants_dict,
                 db_path=None):
        self.pts = pts
        self.outline_segs = outline_segs
        self.inner_segs = inner_segs
        self.radii = radii
        self.constants = constants_dict

        self.wall_t = constants_dict.get("WALL_OUTER", 8.0 / 12.0)
        self.outline_poly = path_polygon(outline_segs, pts)
        self.inner_poly = path_polygon(inner_segs, pts)
        self.outer_area = poly_area(self.outline_poly)

        # Shell geometry (Phase 15-B) — modifies inner_poly with F8-F9 polyline
        self._compute_shell_geometry(constants_dict)

        # Inner area computed after F8-F9 replacement (matches build_floorplan_data)
        self.inner_area = poly_area(self.inner_poly)

        # Roof geometry (Phase 15-C)
        self._compute_roof_geometry(constants_dict)

        # Layout from procedural modules (for generator compatibility)
        self._compute_layout()

    def _compute_shell_geometry(self, constants_dict):
        """Compute S/G-series shell paths, F8-F9 polylines, and wall sections."""
        shell_t = constants_dict.get("SHELL_THICKNESS", 2.0 / 12.0)
        air_gap = constants_dict.get("AIR_GAP",
                                     self.wall_t - 2 * shell_t)
        opening_r = constants_dict.get("OPENING_INSIDE_RADIUS", 1.5 / 12.0)
        f8f9_inner_r = opening_r + shell_t

        # W-series F8-F9 straight-arc-straight polyline
        self.w_f8f9_poly = f8f9_corner_polyline(
            self.pts, self.wall_t, f8f9_inner_r)

        # Replace W8-W9 arc in inner_poly with straight-arc-straight path
        w8 = self.pts["W8"]
        w9 = self.pts["W9"]
        w8_idx = next(i for i, p in enumerate(self.inner_poly)
                      if abs(p[0] - w8[0]) < GEOM_EPS
                      and abs(p[1] - w8[1]) < GEOM_EPS)
        w9_idx = next(i for i, p in enumerate(self.inner_poly)
                      if i > w8_idx
                      and abs(p[0] - w9[0]) < GEOM_EPS
                      and abs(p[1] - w9[1]) < GEOM_EPS)
        self.inner_poly[w8_idx:w9_idx + 1] = self.w_f8f9_poly

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
        layout = compute_interior_layout(self.pts, self.inner_poly)
        outer_openings = compute_outer_openings(self.pts, layout)
        self.openings = outer_to_wall_openings(
            outer_openings, self.outline_segs, self.pts)

        # Wall sections
        self.wall_sections = enumerate_wall_sections(
            self.openings, self.outline_segs)

        # Store layout for later use
        self._layout = layout

    def _compute_roof_geometry(self, constants_dict):
        """Compute R-series roof geometry."""
        overhang = constants_dict.get("ROOF_OVERHANG", 1.5)
        self.roof = compute_roof_geometry(self.pts, self.radii)
        self.roof_poly = roof_polyline(self.roof)
        # Merge roof points into main pts dict
        self.pts.update(self.roof.pts)

    def _compute_layout(self):
        """Store the procedural layout for generator compatibility."""
        self.layout = self._layout


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
        from app.engine import _compute_traverse_from_db
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
        from app.engine import _build_outline_segs_from_chain, _derive_constant
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


def build_generator_data(constants_dict, chain_rows=None, db_path=None):
    """Build a GeneratorData object from DB state.

    This is the main entry point for SVG generators.  It computes all
    geometry from the database and returns a GeneratorData with native
    Python objects ready for rendering.
    """
    pts, outline_segs, inner_segs, radii = compute_native_geometry(
        constants_dict, chain_rows=chain_rows, db_path=db_path)
    return GeneratorData(pts, outline_segs, inner_segs, radii,
                         constants_dict, db_path=db_path)
