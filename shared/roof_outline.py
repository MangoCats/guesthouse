"""DB-driven roof outline geometry: external tangent lines from arc corner centers.

Algorithm
---------
Given N roof corners (arc centers) in CW order:
  1. For each adjacent pair (C_i, C_{i+1}), compute the *outside* external
     common tangent of the two roof-radius circles (outline radius + overhang).
     "Outside" = left of the CW C_i→C_{i+1} direction = building exterior.
  2. For each corner:
     - *Sharp*: polygon vertex = intersection of the two adjacent tangent lines.
     - *Radiused*: polygon arc (CW) from the tangent point on the incoming line
       to the tangent point on the outgoing line, around the corner center.
  3. The resulting polygon is the roof outline; area is computed via poly_area.
"""
import math
from typing import NamedTuple

from shared.geometry import line_isect, poly_area
from shared.types import Point


class DbRoofResult(NamedTuple):
    """Result of DB-driven roof outline computation."""
    poly: list         # sampled polygon [(E, N), ...]
    area: float        # sq ft
    corner_pts: dict   # center_name → (E, N) label position
    corner_radii: dict # center_name → roof radius in feet


def _ext_tangent_outside(ca: Point, ra: float,
                          cb: Point, rb: float) -> tuple:
    """Outside external tangent of circles (ca, ra) and (cb, rb).

    "Outside" is the tangent on the left side of the CW ca→cb direction,
    which is the building exterior for a CW outline traversal.

    Derivation
    ----------
    The outward unit normal n satisfies:
        n · (ca - cb) = ra - rb
    → n · e = (ra - rb) / d    where e = unit(ca→cb), d = |ca - cb|

    Normal decomposition in the (e, left_perp(e)) frame:
        n = cos_t * e + sin_t * left_perp(e)
        cos_t = (ra - rb) / d
        sin_t = sqrt(1 - cos_t²)   ← positive ↔ left (exterior) side

    Tangent point on each circle: C + R * n.

    Returns
    -------
    (pa, pb, nx, ny)
        pa  — tangent point on circle ca
        pb  — tangent point on circle cb
        nx, ny — components of outward unit normal
    """
    dx, dy = cb[0] - ca[0], cb[1] - ca[1]
    d = math.sqrt(dx * dx + dy * dy)
    if d < 1e-12:
        # Degenerate: coincident centers — fall back to perpendicular
        nx, ny = 0.0, 1.0
    else:
        ex, ey = dx / d, dy / d
        cos_t = (ra - rb) / d
        cos_t = max(-1.0, min(1.0, cos_t))          # clamp for float safety
        sin_t = math.sqrt(max(0.0, 1.0 - cos_t * cos_t))
        # left_perp(e) = (-ey, ex)
        nx = cos_t * ex + sin_t * (-ey)
        ny = cos_t * ey + sin_t * ex

    pa = (ca[0] + ra * nx, ca[1] + ra * ny)
    pb = (cb[0] + rb * nx, cb[1] + rb * ny)
    return pa, pb, nx, ny


def _sample_arc_cw(center: Point, radius: float,
                   start_pt: Point, end_pt: Point, n: int = 30) -> list:
    """Interior CW arc sample points (start and end endpoints excluded)."""
    cx, cy = center
    a_start = math.atan2(start_pt[1] - cy, start_pt[0] - cx)
    a_end   = math.atan2(end_pt[1]   - cy, end_pt[0]   - cx)
    if a_end > a_start:
        a_end -= 2.0 * math.pi          # ensure CW (decreasing angle)
    return [
        (cx + radius * math.cos(a_start + (a_end - a_start) * i / n),
         cy + radius * math.sin(a_start + (a_end - a_start) * i / n))
        for i in range(1, n)
    ]


def compute_db_roof_outline(corner_names: list,
                             corner_radiused: list,
                             pts: dict,
                             radii: dict,
                             overhang: float,
                             n_arc: int = 30) -> DbRoofResult:
    """Compute DB-driven roof outline polygon.

    Parameters
    ----------
    corner_names    CW-ordered list of arc center names (e.g. ['C01','C03',...])
    corner_radiused Parallel bool list; True = arc corner, False = sharp corner
    pts             Geometry points dict (must include all center names)
    radii           Radii dict from outline_solver: 'R_a{suffix}' → radius
    overhang        Overhang distance in feet
    n_arc           Arc discretisation point count for radiused corners

    Returns
    -------
    DbRoofResult with sampled polygon, area, label positions, and roof radii.
    """
    N = len(corner_names)
    if N < 3:
        raise ValueError("roof outline requires at least 3 corner centers")

    # Derive the arc center point name from the corner name.
    # R-series names (R01, R03, …) map to C-series centers (C01, C03, …).
    # C-series names are passed through unchanged (legacy / direct use).
    def _cname(name):
        return ("C" + name[1:]) if name.startswith("R") else name

    # Roof arc radius for each corner = outline radius + overhang
    r_roof = []
    for name in corner_names:
        ra_key = "R_a" + name[1:]      # 'R01'/'C01' → 'R_a01'
        r_roof.append(radii[ra_key] + overhang)

    # Tangent line for each consecutive pair i → (i+1)%N
    # Each entry: (pa, pb, nx, ny)
    #   pa  = tangent point on circle i
    #   pb  = tangent point on circle (i+1)%N
    #   n   = outward unit normal of the tangent line
    tangent_lines = []
    for i in range(N):
        j = (i + 1) % N
        ca, ra = pts[_cname(corner_names[i])], r_roof[i]
        cb, rb = pts[_cname(corner_names[j])], r_roof[j]
        tangent_lines.append(_ext_tangent_outside(ca, ra, cb, rb))

    # Build polygon and label positions
    poly = []
    corner_pts = {}

    for i in range(N):
        prev_i = (i - 1) % N
        pa_in,  pb_in,  nx_in,  ny_in  = tangent_lines[prev_i]
        pa_out, pb_out, nx_out, ny_out = tangent_lines[i]

        # pb_in  = tangent point of incoming line ON circle i (entry point)
        # pa_out = tangent point of outgoing line ON circle i (exit point)

        if not corner_radiused[i]:
            # Sharp: intersection of the two tangent lines
            d_in  = (-ny_in,  nx_in)   # line direction ⊥ to normal
            d_out = (-ny_out, nx_out)
            v = line_isect(pb_in, d_in, pa_out, d_out)
            poly.append(v)
            corner_pts[corner_names[i]] = v
        else:
            # Radiused: CW arc from pb_in to pa_out around corner center
            center = pts[_cname(corner_names[i])]
            r = r_roof[i]
            tp_in  = pb_in   # == center + r * (nx_in, ny_in)
            tp_out = pa_out  # == center + r * (nx_out, ny_out)
            poly.append(tp_in)
            arc_pts = _sample_arc_cw(center, r, tp_in, tp_out, n_arc)
            poly.extend(arc_pts)
            poly.append(tp_out)
            # Label at arc midpoint
            mid = arc_pts[len(arc_pts) // 2] if arc_pts else (
                (tp_in[0] + tp_out[0]) / 2, (tp_in[1] + tp_out[1]) / 2)
            corner_pts[corner_names[i]] = mid

    area = poly_area(poly)
    return DbRoofResult(
        poly=poly,
        area=area,
        corner_pts=corner_pts,
        corner_radii=dict(zip(corner_names, r_roof)),
    )


def db_roof_segments(corner_names: list,
                     corner_radiused: list,
                     pts: dict,
                     radii: dict,
                     overhang: float) -> list:
    """Compute DB-driven roof outline as T-path arc/line segments.

    Returns a list of ("line", x1, y1, x2, y2) and
    ("arc", cx, cy, r, a1_deg, a2_deg) tuples — same format as
    floorplan.roof.roof_segments().  Radiused corners produce CW arc elements;
    straight sections produce line elements.  Arc geometry is preserved exactly
    (not sampled to polylines).

    Parameters
    ----------
    corner_names    CW-ordered list of arc center names (e.g. ['R01','R03',...])
    corner_radiused Parallel bool list; True = arc corner, False = sharp corner
    pts             Geometry points dict (must include all C-series center names)
    radii           Radii dict: 'R_a{suffix}' → radius in feet
    overhang        Overhang distance in feet
    """
    N = len(corner_names)
    if N < 3:
        raise ValueError("roof outline requires at least 3 corner centers")

    def _cname(name):
        return ("C" + name[1:]) if name.startswith("R") else name

    # Roof arc radius for each corner = outline radius + overhang
    r_roof = []
    for name in corner_names:
        ra_key = "R_a" + name[1:]      # 'R01'/'C01' → 'R_a01'
        r_roof.append(radii[ra_key] + overhang)

    # Tangent line for each consecutive pair i → (i+1)%N
    tangent_lines = []
    for i in range(N):
        j = (i + 1) % N
        ca, ra = pts[_cname(corner_names[i])], r_roof[i]
        cb, rb = pts[_cname(corner_names[j])], r_roof[j]
        tangent_lines.append(_ext_tangent_outside(ca, ra, cb, rb))

    # Corner geometry: entry/exit tangent points and sharp vertices
    corner_entry = []   # tangent point where incoming line arrives at corner i
    corner_exit = []    # tangent point where outgoing line leaves corner i
    corner_vertex = []  # sharp vertex (None for radiused corners)

    for i in range(N):
        prev_i = (i - 1) % N
        _, pb_in, nx_in, ny_in = tangent_lines[prev_i]
        pa_out, _, nx_out, ny_out = tangent_lines[i]

        corner_entry.append(pb_in)
        corner_exit.append(pa_out)

        if not corner_radiused[i]:
            d_in  = (-ny_in,  nx_in)
            d_out = (-ny_out, nx_out)
            v = line_isect(pb_in, d_in, pa_out, d_out)
            corner_vertex.append(v)
        else:
            corner_vertex.append(None)

    def _pt_exit(i):
        return corner_exit[i] if corner_radiused[i] else corner_vertex[i]

    def _pt_entry(i):
        return corner_entry[i] if corner_radiused[i] else corner_vertex[i]

    def _arc_cw(center, radius, sp, ep):
        a1 = math.degrees(math.atan2(sp[1] - center[1], sp[0] - center[0]))
        a2 = math.degrees(math.atan2(ep[1] - center[1], ep[0] - center[0]))
        sweep = (a1 - a2) % 360
        return ("arc", center[0], center[1], radius, a1, a1 - sweep)

    segments = []
    for i in range(N):
        if corner_radiused[i]:
            center = pts[_cname(corner_names[i])]
            segments.append(_arc_cw(center, r_roof[i],
                                    corner_entry[i], corner_exit[i]))

        next_i = (i + 1) % N
        start = _pt_exit(i)
        end = _pt_entry(next_i)
        segments.append(("line", start[0], start[1], end[0], end[1]))

    return segments
