"""Geometry computation engine that reads from the database.

Patches floorplan.constants with DB values, runs existing computation
pipeline, and returns JSON-serialisable results.
"""
import math
import os
import subprocess
import sys

from app.apputil import point_to_list, bbox_from_poly, seg_to_dict
from app.database import get_variant_exclusions

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


def _compute_room_labels(pts, layout, inner_segs, radii, variant):
    """Compute room label positions, areas, and SF partition lines.

    Returns dict with:
      room_labels: [{name, pos, area (sf variant only)}]
      sf_lines: [{start, end}]  (dashed partition lines, sf variant only)
    """
    from shared.geometry import seg_vecs, line_isect, offset_pt
    from floorplan.openings import compute_outer_openings, compute_rough_openings

    ro_list = compute_rough_openings(pts, layout)
    outer_openings = compute_outer_openings(pts, layout)
    o6 = next(o for o in outer_openings if o.name == "O6")
    o6_w = o6.poly[0]
    ro1_bd = next(r for r in ro_list if r.name == "RO1").poly
    ro1_w_nf = ro1_bd[3]
    w9w10_al, _ = seg_vecs(pts["W9"], pts["W10"])
    w2w5_al, _ = seg_vecs(pts["W2"], pts["W5"])
    w18w1_al, _ = seg_vecs(pts["W18"], pts["W1"])
    iw3_nw = layout.iw3.poly[3]
    iw3_w2w5 = line_isect(iw3_nw, w18w1_al, pts["W2"], w2w5_al)
    iw2s_e_al, _ = seg_vecs(layout.iw2s.poly[1], layout.iw2s.poly[2])
    iw2s_at_w9 = line_isect(layout.iw2s.poly[1], iw2s_e_al, pts["W9"], w9w10_al)

    labels = []

    # BEDROOM: between IW9 and IW11, IW1 south to W1 south
    bd_e = (layout.iw9.poly[2][0] + layout.iw11.poly[3][0]) / 2
    bd_n = (layout.iw1.poly[0][1] + pts["W1"][1]) / 2
    labels.append({"name": "BEDROOM", "pos": point_to_list((bd_e, bd_n))})

    # UTIL (north): between IW3-W2W5 line and IW8, horizontal center between IW8/IW9
    ut_e = (layout.iw8.poly[0][0] + layout.iw9.poly[3][0]) / 2
    ut_n = (iw3_w2w5[1] + layout.iw8.poly[0][1]) / 2
    labels.append({"name": "UTIL", "pos": point_to_list((ut_e, ut_n))})

    # KITCHEN: between RO1-O6 (E-W) and IW2-W9 (N-S)
    k_e = (ro1_w_nf[0] + o6_w[0]) / 2
    k_n = (layout.iw2.poly[1][1] + iw2s_at_w9[1]) / 2
    labels.append({"name": "KITCHEN", "pos": point_to_list((k_e, k_n))})

    # LIVING: between O6 and the E wall arc, N of IW1
    o6_cx = sum(p[0] for p in o6.poly) / len(o6.poly)
    lv_n = (layout.iw1.poly[2][1] + pts["W10"][1]) / 2
    labels.append({"name": "LIVING", "pos": point_to_list((o6_cx, lv_n))})

    # BATH: between IW2s west and W2-W5, vertically between IW8 and W6
    ba_e = (layout.iw2s.poly[0][0] + pts["W2"][0]) / 2
    ba_n = (layout.iw8.poly[3][1] + pts["W6"][1]) / 2
    labels.append({"name": "BATH", "pos": point_to_list((ba_e, ba_n))})

    # OFFICE: between IW4 east and W15, vertically between IW1 and IW5/counter
    of_e = (layout.iw4.poly[1][0] + pts["W15"][0]) / 2
    of_n = (layout.iw1.poly[0][1] + layout.iw5.poly[0][1]) / 2
    labels.append({"name": "OFFICE", "pos": point_to_list((of_e, of_n))})

    # WH: between IW2s east and W8, vertically between W7 and W9
    wh_e = (layout.iw2s.poly[2][0] + pts["W8"][0]) / 2
    wh_n = (pts["W7"][1] + pts["W9"][1]) / 2
    labels.append({"name": "WH", "pos": point_to_list((wh_e, wh_n))})

    # E CLOSET: between IW11 and IW4, south of IW12
    ec_e = (layout.iw11.poly[1][0] + layout.iw4.poly[0][0]) / 2
    ec_n = (layout.iw12.poly[0][1] + pts["W1"][1]) / 2
    labels.append({"name": "E CLOSET", "pos": point_to_list((ec_e, ec_n))})

    # W CLOSET: between IW3 and IW9, south of IW7
    wc_e = (layout.iw3.poly[1][0] + layout.iw9.poly[0][0]) / 2
    wc_n = (layout.iw7.poly[0][1] + pts["W1"][1]) / 2
    labels.append({"name": "W CLOSET", "pos": point_to_list((wc_e, wc_n))})

    # STORAGE: between IW11 and W14, south of IW5
    st_e = (layout.iw11.poly[1][0] + pts["W14"][0]) / 2
    st_n = (layout.iw5.poly[3][1] + layout.iw1.poly[0][1]) / 2
    labels.append({"name": "STORAGE", "pos": point_to_list((st_e, st_n))})

    result = {"room_labels": labels}

    # SF-specific: room areas and dashed partition lines
    if variant == "sf":
        from floorplan.gen_floorplan import compute_room_areas
        from collections import namedtuple
        Data = namedtuple("Data", ["pts", "inner_segs", "radii"])
        data = Data(pts=pts, inner_segs=inner_segs, radii=radii)
        areas = compute_room_areas(data, layout)
        for lbl in labels:
            key = lbl["name"]
            if key in areas:
                lbl["area"] = round(areas[key], 1)
            elif key == "UTIL":
                lbl["area"] = round(areas.get("UTIL_N", 0) + areas.get("UTIL_S", 0), 1)

        # Three dashed partition lines
        sf_lines = []
        # 1. RO1 west (IW1 north) → O6 west (W-surface)
        sf_lines.append({"start": point_to_list(ro1_w_nf), "end": point_to_list(o6_w)})
        # 2. W9 → IW2s east face at W9 northing
        sf_lines.append({"start": point_to_list(pts["W9"]), "end": point_to_list(iw2s_at_w9)})
        # 3. IW3 NW → W2-W5 face
        sf_lines.append({"start": point_to_list(iw3_nw), "end": point_to_list(iw3_w2w5)})
        result["sf_lines"] = sf_lines

    return result


def compute_geometry(constants_dict: dict, variant: str = "standard") -> dict:
    """Compute all building geometry from constants and return JSON-serialisable dict."""
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
    geom = compute_outline_geometry()
    pts = dict(geom.fp_pts)
    pts.update(trav_pts)

    # 3. Inner walls (W-series)
    inner_segs = compute_inner_walls(geom.outline_segs, pts, constants_dict.get("WALL_OUTER", 8.0/12.0), geom.radii)
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
    result["outline_segments"] = [seg_to_dict(s) for s in geom.outline_segs]

    # Inner segments
    result["inner_segments"] = [seg_to_dict(s) for s in inner_segs]

    # Outline polygon
    outline_poly = path_polygon(geom.outline_segs, pts)
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
    variant_items = compute_variant_items(pts, inner_poly, layout, geom.radii, variant)
    result["variant_items"] = variant_items
    result["variant"] = variant
    result["available_variants"] = list(VARIANTS.keys())

    # Dimension lines
    from floorplan.gen_floorplan import compute_dimension_endpoints
    dim_endpoints = compute_dimension_endpoints(pts, layout, geom.radii, bare=(variant == "bare"))
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
    room_data = _compute_room_labels(pts, layout, inner_segs, geom.radii, variant)
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
