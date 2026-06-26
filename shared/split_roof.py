"""Split single-pitch roof geometry (the 'split2' roof style, MarkZ).

The roof is divided by a north-south seam line into two planar sections:

* **West section** — identical to the existing 2:12 shed roof (low eave at the
  south, sloping up 2:12 to the north).  Underside plane: ``z = slope*N + z_off``.

* **East section** — a single tilted plane that:
    1. coincides with the west plane along the N-S seam line (so the two roofs
       intersect on that line continuously from the south end up to the north,
       rising at 2:12), and
    2. meets the closet's diagonal south wall face at a *constant* elevation
       (the level "back wall" / low eave of the east plane).

The east plane ``z = A + B*E + C*N`` is fully determined by those two
conditions:

    C = slope                                  (N-S slope == west plane)
    B = -C * wy / wx                            (level along the closet S wall)
    A = z_off - B * seam_x                      (coincide along the seam)

where ``(wx, wy)`` is the plan direction of the closet south wall and
``seam_x`` is the easting of the seam (the centreline of the short N-S exterior
wall east of the office / half-bath).

Coordinates are FC-origin feet (E = +X east, N = +Y north, Z = up).
"""
from typing import NamedTuple


class EastPlane(NamedTuple):
    a: float          # z = a + b*E + c*N
    b: float
    c: float
    seam_x: float     # easting of the N-S seam line
    slope: float      # west (2:12) plane N-S slope
    z_off: float      # west plane: z = slope*N + z_off
    low_elev: float   # constant underside elevation along the closet S wall

    def z(self, x, y):
        return self.a + self.b * x + self.c * y

    def z_west(self, y):
        return self.slope * y + self.z_off


def compute_east_plane(pts, *, slope, eave_elev, ref_y, wall_outer,
                       seam_wall="F23ab", closet_wall=("F17", "F18"),
                       west_extension=0.0):
    """Solve the east roof plane for the split2 roof style.

    Parameters
    ----------
    pts        geometry points dict (must contain seam_wall + closet_wall names)
    slope      west 2:12 plane N-S slope (e.g. SHED_ROOF_SLOPE)
    eave_elev  west plane south-eave underside elevation (SHED_ROOF_EAVE_ELEV)
    ref_y      N coordinate the eave_elev is referenced to (min N of roof outline)
    wall_outer outer wall thickness (the seam is the wall centreline)
    seam_wall  point name on the outer face of the short N-S seam wall
    closet_wall (start, end) point names of the closet's diagonal south wall face
    west_extension  feet to extend the west plane further east (shifts the seam
                    east by this much; the east plane is re-solved to stay
                    coincident with the west plane along the new seam line)
    """
    z_off = eave_elev - slope * ref_y

    # Seam = centreline of the short N-S wall whose outer face is at seam_wall,
    # then shifted east by west_extension (extending the west plane eastward).
    # The wall's interior (office) is to the west, so the centreline is half a
    # wall-thickness west (more negative E) of the outer face.
    seam_x = pts[seam_wall][0] - wall_outer / 2.0 + west_extension

    pa, pb = pts[closet_wall[0]], pts[closet_wall[1]]
    wx, wy = pb[0] - pa[0], pb[1] - pa[1]

    c = slope
    b = -c * wy / wx
    a = z_off - b * seam_x
    low_elev = a + b * pa[0] + c * pa[1]   # constant along the closet S wall
    return EastPlane(a=a, b=b, c=c, seam_x=seam_x, slope=slope,
                     z_off=z_off, low_elev=low_elev)


def split_polygon_x(poly, seam_x):
    """Clip a 2D polygon by the vertical line E = seam_x.

    Returns (west_poly, east_poly): the parts with E <= seam_x and E >= seam_x.
    Either may be empty.  Standard Sutherland-Hodgman half-plane clipping.
    """
    def clip(keep_west):
        out = []
        n = len(poly)
        for i in range(n):
            cur = poly[i]
            nxt = poly[(i + 1) % n]
            cur_in = (cur[0] <= seam_x) if keep_west else (cur[0] >= seam_x)
            nxt_in = (nxt[0] <= seam_x) if keep_west else (nxt[0] >= seam_x)
            if cur_in:
                out.append(cur)
            if cur_in != nxt_in:
                # Edge crosses the seam — add the intersection point.
                dx = nxt[0] - cur[0]
                t = (seam_x - cur[0]) / dx if abs(dx) > 1e-12 else 0.0
                out.append((seam_x, cur[1] + t * (nxt[1] - cur[1])))
        return out

    return clip(True), clip(False)
