"""Geometry computation engine that reads from the database.

Patches floorplan.constants with DB values, runs existing computation
pipeline, and returns JSON-serialisable results.
"""
import math
import os
import subprocess
import sys

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


def _point_to_list(pt):
    """Convert (e, n) tuple to [e, n] list for JSON."""
    return [round(pt[0], 6), round(pt[1], 6)]


def _wall_to_dict(wall):
    """Convert Wall namedtuple to JSON dict."""
    return {
        "poly": [_point_to_list(p) for p in wall.poly],
        "bbox": {"w": wall.w, "s": wall.s, "e": wall.e, "n": wall.n},
    }


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

    # Build result
    result = {
        "points": {},
        "outline_segments": [],
        "inner_segments": [],
        "outline_poly": [],
        "inner_poly": [_point_to_list(p) for p in inner_poly],
        "interior_walls": {},
        "outer_openings": [],
        "rough_openings": [],
        "appliances": {},
        "furniture": {},
    }

    # Points
    for name, pt in sorted(pts.items()):
        result["points"][name] = _point_to_list(pt)

    # Outline segments
    from shared.types import LineSeg, ArcSeg
    for seg in geom.outline_segs:
        if isinstance(seg, LineSeg):
            result["outline_segments"].append({
                "type": "line", "start": seg.start, "end": seg.end,
            })
        else:
            result["outline_segments"].append({
                "type": "arc", "start": seg.start, "end": seg.end,
                "center": seg.center, "radius": seg.radius,
                "direction": seg.direction, "n_pts": seg.n_pts,
            })

    # Inner segments
    for seg in inner_segs:
        if isinstance(seg, LineSeg):
            result["inner_segments"].append({
                "type": "line", "start": seg.start, "end": seg.end,
            })
        else:
            result["inner_segments"].append({
                "type": "arc", "start": seg.start, "end": seg.end,
                "center": seg.center, "radius": seg.radius,
                "direction": seg.direction, "n_pts": seg.n_pts,
            })

    # Outline polygon
    outline_poly = path_polygon(geom.outline_segs, pts)
    result["outline_poly"] = [_point_to_list(p) for p in outline_poly]

    # Interior walls
    iw_names = ["iw1", "iw2", "iw2o", "iw2s", "iw3", "iw4", "iw5", "iw6",
                "iw7", "iw8", "iw9", "iw11", "iw12"]
    for name in iw_names:
        wall = getattr(layout, name, None)
        if wall:
            result["interior_walls"][name.upper()] = _wall_to_dict(wall)

    # Outer openings
    for op in outer_openings:
        result["outer_openings"].append({
            "name": op.name,
            "seg_start": op.seg_start,
            "seg_end": op.seg_end,
            "poly": [_point_to_list(p) for p in op.poly],
        })

    # Rough openings
    for ro in rough_openings:
        d = {
            "name": ro.name,
            "bbox": {"w": ro.bbox.w, "s": ro.bbox.s, "e": ro.bbox.e, "n": ro.bbox.n},
            "wall_name": ro.wall_name,
            "orientation": ro.orientation,
            "width": ro.width,
        }
        if ro.poly:
            d["poly"] = [_point_to_list(p) for p in ro.poly]
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
            result["appliances"]["counter"]["clip"] = [_point_to_list(p) for p in clip]

    # Furniture
    for name in ("bed", "dresser", "shelves"):
        item = getattr(layout, name, None)
        if item:
            result["furniture"][name] = _wall_to_dict(item)

    # Bounding box
    all_pts = result["outline_poly"]
    if all_pts:
        es = [p[0] for p in all_pts]
        ns = [p[1] for p in all_pts]
        result["bbox"] = {
            "w": min(es), "s": min(ns), "e": max(es), "n": max(ns),
        }

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
                "A": _point_to_list(pa),
                "B": _point_to_list(pb),
                "dist": round(dist, 6),
            }
    result["dimensions"] = dimensions

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
