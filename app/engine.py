"""Geometry computation engine that reads from the database.

Patches floorplan.constants with DB values, runs existing computation
pipeline, and returns JSON-serialisable results.
"""
import math
import os
import subprocess
import sys

from app.apputil import point_to_list, bbox_from_poly, seg_to_dict
from app.database import get_variant_exclusions, get_room_label_offsets

_PROJECT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def patch_constants(constants_dict: dict):
    """Monkey-patch floorplan.constants with database values."""
    import floorplan.constants as mod
    for name, value in constants_dict.items():
        if hasattr(mod, name):
            setattr(mod, name, value)
    # Recompute derived constants that depend on others
    mod.WALL_EXTRA = mod.WALL_OUTER - 8.0 / 12.0
    mod.AIR_GAP = mod.WALL_OUTER - 2 * mod.SHELL_THICKNESS
    mod.DOOR_FLAT_FACE = mod.WALL_OUTER - 2 * (mod.OPENING_INSIDE_RADIUS + mod.SHELL_THICKNESS)
    mod.F8F9_INNER_TURN_R = mod.OPENING_INSIDE_RADIUS + mod.SHELL_THICKNESS
    mod.CORNER_SW_R = 10.0 / 12.0 + mod.WALL_EXTRA


def _wall_to_dict(wall):
    """Convert Wall namedtuple to JSON dict."""
    return {
        "poly": [point_to_list(p) for p in wall.poly],
        "bbox": {"w": wall.w, "s": wall.s, "e": wall.e, "n": wall.n},
    }


def _centroid(poly):
    """Area-weighted centroid of a simple polygon."""
    n = len(poly)
    area6 = 0.0
    cx = cy = 0.0
    for i in range(n):
        x0, y0 = poly[i]
        x1, y1 = poly[(i + 1) % n]
        cross = x0 * y1 - x1 * y0
        area6 += cross
        cx += (x0 + x1) * cross
        cy += (y0 + y1) * cross
    if abs(area6) < 1e-12:
        # Degenerate — fall back to vertex average
        return (sum(p[0] for p in poly) / n, sum(p[1] for p in poly) / n)
    cx /= 3.0 * area6
    cy /= 3.0 * area6
    return (cx, cy)


def _compute_room_labels(pts, layout, inner_segs, radii, variant):
    """Compute room label positions, areas, and SF partition lines.

    Builds the actual room-area polygons (matching compute_room_areas) so
    labels sit at the true centroid and SF variant can highlight on click.

    Returns dict with:
      room_labels: [{name, pos, area?, poly?}]
      sf_lines: [{start, end}]  (sf variant only)
    """
    from shared.geometry import seg_vecs, line_isect, segment_polyline
    from floorplan.openings import compute_outer_openings, compute_rough_openings

    ro_list = compute_rough_openings(pts, layout)
    outer_openings = compute_outer_openings(pts, layout)
    o6 = next((o for o in outer_openings if o.name == "O6"), None)
    if o6 is None:
        return []
    o6_w = o6.poly[0]
    ro1 = next((r for r in ro_list if r.name == "RO1"), None)
    if ro1 is None:
        return []
    ro1_bd = ro1.poly
    ro1_w_nf = ro1_bd[3]
    w9w10_al, _ = seg_vecs(pts["W9"], pts["W10"])
    w2w5_al, _ = seg_vecs(pts["W2"], pts["W5"])
    w18w1_al, _ = seg_vecs(pts["W18"], pts["W1"])
    iw3_nw = layout.iw3.poly[3]
    iw3_w2w5 = line_isect(iw3_nw, w18w1_al, pts["W2"], w2w5_al)
    iw2s_e_al, _ = seg_vecs(layout.iw2s.poly[1], layout.iw2s.poly[2])
    iw2s_at_w9 = line_isect(layout.iw2s.poly[1], iw2s_e_al, pts["W9"], w9w10_al)

    is_sf = (variant == "sf")

    # --- Build room polygons (matching compute_room_areas) ---
    rooms = {}

    # BEDROOM
    rooms["BEDROOM"] = [
        (layout.iw9.poly[2][0], layout.iw1.poly[0][1]),
        (layout.iw11.poly[3][0], layout.iw1.poly[0][1]),
        (layout.iw11.poly[3][0], pts["W1"][1]),
        (layout.iw9.poly[2][0], pts["W1"][1]),
    ]

    # UTIL_N
    rooms["UTIL_N"] = [
        iw3_w2w5,
        layout.iw8.poly[0], layout.iw8.poly[1],
        layout.iw1.poly[0], layout.iw9.poly[3],
        (layout.iw9.poly[0][0], layout.iw7.poly[2][1]),
        layout.iw7.poly[3], iw3_nw,
    ]

    # UTIL_S
    _util_s = [iw3_nw, (layout.iw3.poly[0][0], pts["W1"][1]), pts["W1"]]
    _util_s.extend(segment_polyline(inner_segs[0], pts)[1:])
    _util_s.append(iw3_w2w5)
    rooms["UTIL_S"] = _util_s

    # KITCHEN
    rooms["KITCHEN"] = [
        o6_w, iw2s_at_w9,
        layout.iw2s.poly[1], layout.iw2o.poly[3],
        layout.iw2o.poly[0], layout.iw2.poly[2],
        layout.iw2.poly[1], ro1_w_nf,
    ]

    # LIVING
    _living = [o6_w]
    _living.append(segment_polyline(inner_segs[6], pts)[-1])
    for si in range(7, 13):
        _living.extend(segment_polyline(inner_segs[si], pts)[1:])
    _living.append(layout.iw1.poly[2])
    _living.append(ro1_w_nf)
    rooms["LIVING"] = _living

    # BATH
    _seg2 = segment_polyline(inner_segs[2], pts)
    _seg3 = segment_polyline(inner_segs[3], pts)
    _bath = [
        layout.iw8.poly[3], layout.iw8.poly[2],
        layout.iw2.poly[3], layout.iw2o.poly[1],
        layout.iw2o.poly[2], layout.iw2s.poly[0],
        layout.iw2s.poly[3], _seg3[0],
    ]
    _bath.extend(reversed(_seg2[:-1]))
    rooms["BATH"] = _bath

    # OFFICE
    _office = [layout.iw5.poly[0]]
    _office.append((pts["W15"][0], layout.iw5.poly[0][1]))
    _office.append(pts["W15"])
    for si in [14, 15, 16]:
        _office.extend(segment_polyline(inner_segs[si], pts)[1:])
    _office.append(layout.iw4.poly[1])
    _office.append(layout.iw4.poly[2])
    _office.append(layout.iw12.poly[2])
    _office.append(layout.iw12.poly[3])
    rooms["OFFICE"] = _office

    # E CLOSET
    rooms["E CLOSET"] = [
        (layout.iw11.poly[1][0], layout.iw12.poly[0][1]),
        (layout.iw4.poly[0][0], layout.iw12.poly[0][1]),
        layout.iw4.poly[0],
        (layout.iw11.poly[1][0], pts["W1"][1]),
    ]

    # W CLOSET
    rooms["W CLOSET"] = [
        (layout.iw3.poly[1][0], layout.iw7.poly[0][1]),
        (layout.iw9.poly[0][0], layout.iw7.poly[0][1]),
        (layout.iw9.poly[0][0], pts["W1"][1]),
        (layout.iw3.poly[1][0], pts["W1"][1]),
    ]

    # STORAGE
    rooms["STORAGE"] = [
        (layout.iw11.poly[1][0], layout.iw5.poly[3][1]),
        (pts["W14"][0], layout.iw5.poly[3][1]),
        (pts["W14"][0], layout.iw1.poly[0][1]),
        (layout.iw11.poly[1][0], layout.iw1.poly[0][1]),
    ]

    # WH
    _seg3b = segment_polyline(inner_segs[3], pts)
    _seg4 = segment_polyline(inner_segs[4], pts)
    _seg5 = segment_polyline(inner_segs[5], pts)
    _wh = [layout.iw2s.poly[2], _seg3b[-1]]
    _wh.extend(_seg4[1:])
    _wh.extend(_seg5[1:])
    _wh.append((layout.iw2s.poly[2][0], pts["W9"][1]))
    rooms["WH"] = _wh

    # --- Get authoritative area values ---
    from floorplan.gen_floorplan import compute_room_areas
    from collections import namedtuple
    Data = namedtuple("Data", ["pts", "inner_segs", "radii"])
    data = Data(pts=pts, inner_segs=inner_segs, radii=radii)
    areas = compute_room_areas(data, layout)

    # --- Build label list from polygon centroids + DB offsets ---
    offsets = get_room_label_offsets()
    labels = []
    for name, poly in rooms.items():
        cx, cy = _centroid(poly)
        de, dn = offsets.get(name, (0.0, 0.0))
        lbl = {
            "name": name,
            "pos": point_to_list((cx + de, cy + dn)),
            "centroid": point_to_list((cx, cy)),
        }
        if is_sf:
            lbl["area"] = round(areas[name], 1)
            lbl["poly"] = [point_to_list(p) for p in poly]
        labels.append(lbl)

    result = {"room_labels": labels}

    # SF-specific: dashed partition lines
    if is_sf:
        result["sf_lines"] = [
            {"start": point_to_list(ro1_w_nf), "end": point_to_list(o6_w)},
            {"start": point_to_list(pts["W9"]), "end": point_to_list(iw2s_at_w9)},
            {"start": point_to_list(iw3_nw), "end": point_to_list(iw3_w2w5)},
        ]

    return result


def _build_outline_segs_from_chain(chain):
    """Build outline_segs from solved chain (matching geometry.py rotation).

    Returns list of LineSeg/ArcSeg in outline convention (F1->F2 first).
    """
    from shared.types import LineSeg, ArcSeg

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


def compute_geometry(constants_dict: dict, variant: str = "standard",
                     chain_rows: list[dict] | None = None) -> dict:
    """Compute all building geometry from constants and return JSON-serialisable dict.

    If chain_rows is provided (Phase 5+), the app solver computes the
    outline from DB chain data, bypassing the module-scope solver in
    floorplan/geometry.py.
    """
    patch_constants(constants_dict)

    # Reload all floorplan modules so module-scope code re-executes with
    # the patched constant values (e.g. OUTLINE_CHAIN closure solver).
    import importlib
    import floorplan.constants
    import floorplan.geometry
    import floorplan.layout
    import floorplan.openings
    importlib.reload(floorplan.constants)
    patch_constants(constants_dict)  # re-patch after reload
    importlib.reload(floorplan.geometry)
    importlib.reload(floorplan.layout)
    importlib.reload(floorplan.openings)

    from floorplan.geometry import compute_outline_geometry, align_pts_to_f_series
    from shared.geometry import compute_inner_walls, path_polygon
    from floorplan.layout import compute_interior_layout
    from floorplan.openings import compute_outer_openings, compute_rough_openings
    from shared.survey import compute_traverse, compute_three_arc, compute_inset

    # 1. Survey traverse
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
        # Phase 5+: use app solver with DB chain
        from app.outline_solver import db_rows_to_chain, solve_closure, walk_chain
        chain = db_rows_to_chain(chain_rows)
        import floorplan.constants as fc
        R_a1 = fc.CORNER_SW_R

        solver_result = solve_closure(chain, R_a1)
        if not solver_result.valid:
            raise ValueError(
                f"Outline chain does not close: error={solver_result.closure_error:.6f}")

        # Inject solved values (distances + closure arc sweep)
        chain = list(chain)
        chain[0] = chain[0]._replace(distance=solver_result.d_F2_F5)
        chain[-2] = chain[-2]._replace(distance=solver_result.d_F18_F1)
        chain[-1] = chain[-1]._replace(sweep=solver_result.sweep_closure)

        F2_E = -18.5
        F2_N = -13.5 + R_a1
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
    inner_segs = compute_inner_walls(outline_segs, pts, constants_dict.get("WALL_OUTER", 8.0/12.0), radii)
    inner_poly = path_polygon(inner_segs, pts)

    # 4. Interior layout
    layout = compute_interior_layout(pts, inner_poly)

    # 5. Outer openings
    outer_openings = compute_outer_openings(pts, layout)

    # 6. Rough openings
    rough_openings = compute_rough_openings(pts, layout)

    # Load variant exclusions from database
    exclusions = get_variant_exclusions(variant)

    # Build result
    result = {
        "points": {},
        "outline_segments": [],
        "inner_segments": [],
        "outline_poly": [],
        "inner_poly": [point_to_list(p) for p in inner_poly],
        "interior_walls": {},
        "outer_openings": [],
        "rough_openings": [],
        "appliances": {},
        "furniture": {},
    }

    # Points
    for name, pt in sorted(pts.items()):
        result["points"][name] = point_to_list(pt)

    # Outline segments
    result["outline_segments"] = [seg_to_dict(s) for s in outline_segs]

    # Inner segments
    result["inner_segments"] = [seg_to_dict(s) for s in inner_segs]

    # Outline polygon
    outline_poly = path_polygon(outline_segs, pts)
    result["outline_poly"] = [point_to_list(p) for p in outline_poly]

    # Interior walls
    excluded_walls = exclusions.get("wall", set())
    iw_names = ["iw1", "iw2", "iw2o", "iw2s", "iw3", "iw4", "iw5", "iw6",
                "iw7", "iw8", "iw9", "iw11", "iw12"]
    for name in iw_names:
        if name.upper() in excluded_walls:
            continue
        wall = getattr(layout, name, None)
        if wall:
            result["interior_walls"][name.upper()] = _wall_to_dict(wall)

    # Outer openings
    for op in outer_openings:
        result["outer_openings"].append({
            "name": op.name,
            "seg_start": op.seg_start,
            "seg_end": op.seg_end,
            "poly": [point_to_list(p) for p in op.poly],
        })

    # Rough openings
    excluded_openings = exclusions.get("rough_opening", set())
    for ro in rough_openings:
        if ro.name in excluded_openings:
            continue
        d = {
            "name": ro.name,
            "bbox": {"w": ro.bbox.w, "s": ro.bbox.s, "e": ro.bbox.e, "n": ro.bbox.n},
            "wall_name": ro.wall_name,
            "orientation": ro.orientation,
            "width": ro.width,
        }
        if ro.poly:
            d["poly"] = [point_to_list(p) for p in ro.poly]
        result["rough_openings"].append(d)

    # Appliances
    for name in ("dryer", "washer"):
        appl = getattr(layout, name, None)
        if appl:
            result["appliances"][name] = _wall_to_dict(appl)

    # Counter
    ctr = getattr(layout, "ctr", None)
    if ctr:
        result["appliances"]["counter"] = _wall_to_dict(ctr)
        clip = getattr(layout, "ctr_clip", None)
        if clip:
            result["appliances"]["counter"]["clip"] = [point_to_list(p) for p in clip]

    # Furniture
    for name in ("bed", "dresser", "shelves"):
        item = getattr(layout, name, None)
        if item:
            result["furniture"][name] = _wall_to_dict(item)

    # Bounding box
    all_pts = result["outline_poly"]
    if all_pts:
        result["bbox"] = bbox_from_poly(all_pts)

    # Variant items
    from app.variants import compute_variant_items, VARIANTS
    variant_items = compute_variant_items(pts, inner_poly, layout, radii, variant)
    result["variant_items"] = variant_items
    result["variant"] = variant
    result["available_variants"] = list(VARIANTS.keys())

    # Dimension lines
    from floorplan.gen_floorplan import compute_dimension_endpoints
    dim_endpoints = compute_dimension_endpoints(pts, layout, radii, bare=(variant == "bare"))
    ep_dict = {name: pt for name, pt in dim_endpoints}
    dim_names = sorted(set(k.rsplit("_", 1)[0] for k in ep_dict))
    dimensions = {}
    for dname in dim_names:
        a_key = f"{dname}_A"
        b_key = f"{dname}_B"
        if a_key in ep_dict and b_key in ep_dict:
            pa, pb = ep_dict[a_key], ep_dict[b_key]
            dist = math.sqrt((pb[0] - pa[0])**2 + (pb[1] - pa[1])**2)
            dimensions[dname] = {
                "A": point_to_list(pa),
                "B": point_to_list(pb),
                "dist": round(dist, 6),
            }
    result["dimensions"] = dimensions

    # Room labels and SF extras
    room_data = _compute_room_labels(pts, layout, inner_segs, radii, variant)
    result["room_labels"] = room_data["room_labels"]
    if "sf_lines" in room_data:
        result["sf_lines"] = room_data["sf_lines"]

    return result


def generate_svg(view_name: str, script_path: str) -> bool:
    """Run a generator script and return True on success."""
    full_path = os.path.join(_PROJECT, script_path)
    if not os.path.exists(full_path):
        return False
    try:
        subprocess.check_call(
            [sys.executable, full_path],
            cwd=_PROJECT,
            timeout=60,
        )
        return True
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired):
        return False


def generate_all_svgs(views: list[dict]) -> dict:
    """Regenerate all enabled views. Returns {name: success} dict."""
    results = {}
    # Group by script to avoid running the same generator twice
    seen_scripts = set()
    for v in views:
        script = v["script"]
        if script in seen_scripts:
            results[v["name"]] = True
            continue
        seen_scripts.add(script)
        results[v["name"]] = generate_svg(v["name"], script)
    return results


def get_svg_content(svg_path: str) -> str | None:
    """Read an SVG file and return its content."""
    full_path = os.path.join(_PROJECT, svg_path)
    if not os.path.exists(full_path):
        return None
    with open(full_path, "r", encoding="utf-8") as f:
        return f.read()
