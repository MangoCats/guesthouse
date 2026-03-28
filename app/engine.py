"""Geometry computation engine — formula-driven architecture (Phase 12g).

All element geometry (walls, openings, furniture, appliances) is computed by
the FormulaEvaluator from database-stored JSON formulas.  The procedural
floorplan modules are used only for outline geometry and as a metadata source
for variant items (labels, door configs, clearance configs, product URLs).
Module-level constants are never patched or reloaded.
"""
import json
import math
import os
import subprocess
import sys

from app.apputil import point_to_list, bbox_from_poly, seg_to_dict
from app.database import get_variant_exclusions, get_room_label_offsets, get_all_elements

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


from shared.geometry import poly_centroid as _centroid


def _compute_door_arcs(outer_openings, rough_openings, doors_data, exclusions):
    """Compute door swing arc geometry for all openings with doors.

    Returns list of door_arc dicts, each with opening_name, door_type,
    and leaves list (each leaf has hinge, tip, arc_pts).
    """
    if not doors_data:
        return []

    # Build lookup: opening_name → door config dict
    door_lookup = {d["opening_name"]: d for d in doors_data}

    # Gather all openings with their polys (dict-based: {"name", "poly", ...})
    all_ops = []
    for op in outer_openings:
        name = op["name"] if isinstance(op, dict) else op.name
        poly = op.get("poly") if isinstance(op, dict) else op.poly
        if poly:
            all_ops.append({"name": name, "poly": poly})
    excluded_openings = exclusions.get("rough_opening", set())
    for ro in rough_openings:
        name = ro["name"] if isinstance(ro, dict) else ro.name
        if name in excluded_openings:
            continue
        poly = ro.get("poly") if isinstance(ro, dict) else ro.poly
        bbox = ro.get("bbox") if isinstance(ro, dict) else getattr(ro, "bbox", None)
        if poly:
            all_ops.append({"name": name, "poly": poly})
        elif bbox:
            if isinstance(bbox, dict):
                all_ops.append({"name": name, "poly": [
                    (bbox["w"], bbox["s"]), (bbox["e"], bbox["s"]),
                    (bbox["e"], bbox["n"]), (bbox["w"], bbox["n"])]})
            else:
                all_ops.append({"name": name, "poly": [
                    (bbox.w, bbox.s), (bbox.e, bbox.s),
                    (bbox.e, bbox.n), (bbox.w, bbox.n)]})

    result = []
    for op in all_ops:
        door = door_lookup.get(op["name"])
        if not door:
            continue

        poly = op["poly"]
        if len(poly) < 4:
            continue

        leaves = _compute_door_leaves(poly, door)
        if leaves:
            result.append({
                "opening_name": op["name"],
                "door_type": door.get("door_type", "single"),
                "leaves": leaves,
            })

    return result


def _compute_door_leaves(poly, door):
    """Compute door leaf geometry from opening polygon and door config.

    poly: 4 vertices of the opening rectangle.  Vertex winding varies by
    opening type — some have poly[0]→poly[1] along the opening width, others
    have it across the wall thickness.  We detect which edge pair is the
    opening width (longer) vs wall thickness (shorter) and orient accordingly.

    door: dict with opening_name, width, hinge_side, swing_direction, door_type.
    """
    p0, p1, p2, p3 = poly[0], poly[1], poly[2], poly[3]

    # Compute both edge lengths
    dx_01 = p1[0] - p0[0]
    dy_01 = p1[1] - p0[1]
    len_01 = math.sqrt(dx_01**2 + dy_01**2)

    dx_03 = p3[0] - p0[0]
    dy_03 = p3[1] - p0[1]
    len_03 = math.sqrt(dx_03**2 + dy_03**2)

    if len_01 < 1e-12 or len_03 < 1e-12:
        return []

    # "along" = opening width direction (longer edge pair)
    # "cross" = wall thickness direction (shorter edge pair)
    if len_01 >= len_03:
        # Normal: poly[0]→poly[1] is along the opening width
        along_len = len_01
        along = (dx_01 / len_01, dy_01 / len_01)
        cross = (dx_03 / len_03, dy_03 / len_03)
        mid_start = ((p0[0] + p3[0]) / 2, (p0[1] + p3[1]) / 2)
        mid_end = ((p1[0] + p2[0]) / 2, (p1[1] + p2[1]) / 2)
    else:
        # Swapped: poly[0]→poly[3] is along the opening width
        along_len = len_03
        along = (dx_03 / len_03, dy_03 / len_03)
        cross = (dx_01 / len_01, dy_01 / len_01)
        mid_start = ((p0[0] + p1[0]) / 2, (p0[1] + p1[1]) / 2)
        mid_end = ((p3[0] + p2[0]) / 2, (p3[1] + p2[1]) / 2)

    # Door width in feet
    door_width_ft = door.get("width", 36) / 12.0
    hinge_side = door.get("hinge_side", "east")
    swing_dir = door.get("swing_direction", "south")
    door_type = door.get("door_type", "single")

    # Determine which end is the hinge side by projecting cardinal directions
    # onto the opening's along vector
    cardinals = {
        "east": (1, 0), "west": (-1, 0),
        "north": (0, 1), "south": (0, -1),
    }

    hinge_vec = cardinals.get(hinge_side, (1, 0))
    # Dot product of hinge direction with along vector
    # Positive → hinge at "end" (poly[1] side), Negative → hinge at "start" (poly[0] side)
    dot_hinge = hinge_vec[0] * along[0] + hinge_vec[1] * along[1]

    # Swing direction vector
    swing_vec = cardinals.get(swing_dir, (0, -1))
    # Determine swing unit vector: project onto cross and along, pick the
    # perpendicular direction that best matches
    dot_swing_cross = swing_vec[0] * cross[0] + swing_vec[1] * cross[1]
    dot_swing_along = swing_vec[0] * along[0] + swing_vec[1] * along[1]

    # The swing direction in opening-local coords
    if abs(dot_swing_cross) >= abs(dot_swing_along):
        # Swing is perpendicular to the opening (through the wall)
        swing_unit = cross if dot_swing_cross > 0 else (-cross[0], -cross[1])
    else:
        # Swing is along the opening direction
        swing_unit = along if dot_swing_along > 0 else (-along[0], -along[1])

    gap = (along_len - door_width_ft) / 2.0
    if gap < 0:
        gap = 0

    if door_type == "double":
        # Double door: two leaves, hinged at opposite ends
        leaf_width = door_width_ft  # width per leaf (from DB, already per-leaf)
        total = 2 * leaf_width
        d_gap = (along_len - total) / 2.0
        if d_gap < 0:
            d_gap = 0

        # Leaf 1: hinge at start end
        h1 = (mid_start[0] + d_gap * along[0], mid_start[1] + d_gap * along[1])
        # Leaf 2: hinge at end end
        h2 = (mid_end[0] - d_gap * along[0], mid_end[1] - d_gap * along[1])

        leaves = []
        for hinge_pt, closed_dir in [(h1, along), (h2, (-along[0], -along[1]))]:
            tip = (hinge_pt[0] + leaf_width * swing_unit[0],
                   hinge_pt[1] + leaf_width * swing_unit[1])
            arc_pts = _swing_arc(hinge_pt, leaf_width, swing_unit, closed_dir)
            leaves.append({
                "hinge": list(hinge_pt),
                "tip": list(tip),
                "arc_pts": [list(p) for p in arc_pts],
            })
        return leaves
    else:
        # Single door: hinge at one end
        if dot_hinge > 0:
            # Hinge at end (poly[1] side)
            hinge_pt = (mid_end[0] - gap * along[0], mid_end[1] - gap * along[1])
            closed_dir = (-along[0], -along[1])
        else:
            # Hinge at start (poly[0] side)
            hinge_pt = (mid_start[0] + gap * along[0], mid_start[1] + gap * along[1])
            closed_dir = along

        tip = (hinge_pt[0] + door_width_ft * swing_unit[0],
               hinge_pt[1] + door_width_ft * swing_unit[1])
        arc_pts = _swing_arc(hinge_pt, door_width_ft, swing_unit, closed_dir)

        return [{
            "hinge": list(hinge_pt),
            "tip": list(tip),
            "arc_pts": [list(p) for p in arc_pts],
        }]


def _swing_arc(hinge, radius, dir_from, dir_to, n_pts=20):
    """Compute 90-degree door swing arc points.

    Replicates floorplan/gen_floorplan.py::_swing_arc_svg math.
    Rotates dir_from toward dir_to (perpendicular unit vectors).
    Returns list of (E, N) tuples.
    """
    cross = dir_from[0] * dir_to[1] - dir_from[1] * dir_to[0]
    sweep = math.pi / 2 if cross > 0 else -math.pi / 2
    pts = []
    for i in range(n_pts + 1):
        a = sweep * i / n_pts
        ca, sa = math.cos(a), math.sin(a)
        ae = hinge[0] + radius * (dir_from[0] * ca - dir_from[1] * sa)
        an = hinge[1] + radius * (dir_from[0] * sa + dir_from[1] * ca)
        pts.append((ae, an))
    return pts


def _compute_clearance_zones(variant, variant_items=None):
    """Compute clearance zone polygons for fixture/furniture items.

    Returns list of clearance zone dicts with name, poly, style.
    All clearance zones are driven by 'clearance' metadata in variant_items
    (face indices + distance), including dresser clearance.
    """
    zones = []
    # Bare and SF variants have no furniture
    if variant in ("bare", "sf"):
        return zones

    # Clearance zones from variant_items metadata
    if variant_items:
        for item_name, item in variant_items.items():
            cl = item.get("clearance")
            if not cl:
                continue
            # Support single clearance dict or list of clearance dicts
            cl_list = cl if isinstance(cl, list) else [cl]
            poly = item["poly"]  # [[e,n], ...]
            cx = sum(p[0] for p in poly) / len(poly)
            cy = sum(p[1] for p in poly) / len(poly)
            for cl_idx, cl_item in enumerate(cl_list):
                face = cl_item["face"]  # [i, j] vertex indices
                dist = cl_item["distance"]
                p_i = poly[face[0]]
                p_j = poly[face[1]]
                # Compute face direction and perpendicular
                dx = p_j[0] - p_i[0]
                dy = p_j[1] - p_i[1]
                face_len = math.sqrt(dx**2 + dy**2)
                if face_len < 1e-12:
                    continue
                # Right-hand perpendicular of face direction
                perp = (-dy / face_len, dx / face_len)
                # Check that perp points away from item centroid
                mid = ((p_i[0] + p_j[0]) / 2, (p_i[1] + p_j[1]) / 2)
                to_center = (cx - mid[0], cy - mid[1])
                if perp[0] * to_center[0] + perp[1] * to_center[1] > 0:
                    perp = (-perp[0], -perp[1])  # flip to point outward
                # Build extension polygon
                ext_i = (p_i[0] + dist * perp[0], p_i[1] + dist * perp[1])
                ext_j = (p_j[0] + dist * perp[0], p_j[1] + dist * perp[1])
                # Single-clearance items keep legacy name; multi-clearance use index suffix
                name = (f"{item_name}_clearance" if len(cl_list) == 1
                        else f"{item_name}_clearance_{cl_idx}")
                zones.append({
                    "name": name,
                    "parent": item_name,
                    "poly": [list(p_i), list(p_j), list(ext_j), list(ext_i)],
                    "style": "dashed",
                })

    return zones


def _compute_appliance_doors(variant_items, variant="standard"):
    """Compute door swing arcs for appliances with 'door' metadata.

    Door config uses intrinsic polygon indices:
      hinge_idx: polygon vertex index of the hinge point
      target_idx: polygon vertex index defining the closed-door direction
    open_dir, closed_dir, and width are computed from the polygon geometry.

    For variant-keyed door configs (fridge, microwave), the config for the
    current variant is selected.  Returns list of appliance_door dicts.
    """
    if not variant_items:
        return []

    result = []
    for item_name, item in variant_items.items():
        door_meta = item.get("door")
        if not door_meta:
            continue

        # Handle variant-keyed door configs: {"standard": {...}, "minik": {...}}
        first_val = next(iter(door_meta.values()), None)
        if isinstance(first_val, dict):
            door = door_meta.get(variant)
            if not door:
                continue
        else:
            door = door_meta

        poly = item["poly"]  # [[e,n], ...]
        hinge_idx = door["hinge_idx"]
        target_idx = door["target_idx"]

        hinge_pt = (poly[hinge_idx][0], poly[hinge_idx][1])
        target_pt = (poly[target_idx][0], poly[target_idx][1])

        # closed_dir: unit vector from hinge to target
        dx = target_pt[0] - hinge_pt[0]
        dy = target_pt[1] - hinge_pt[1]
        width = math.sqrt(dx * dx + dy * dy)
        if width < 1e-12:
            continue
        closed_dir = (dx / width, dy / width)

        # open_dir: perpendicular to closed_dir, pointing away from polygon
        perp = (-closed_dir[1], closed_dir[0])
        cx = sum(p[0] for p in poly) / len(poly)
        cy = sum(p[1] for p in poly) / len(poly)
        to_center = (cx - hinge_pt[0], cy - hinge_pt[1])
        if perp[0] * to_center[0] + perp[1] * to_center[1] > 0:
            perp = (-perp[0], -perp[1])
        open_dir = perp

        tip = (hinge_pt[0] + width * open_dir[0],
               hinge_pt[1] + width * open_dir[1])
        arc_pts = _swing_arc(hinge_pt, width, open_dir, closed_dir)

        entry = {
            "item_name": item_name,
            "hinge": list(hinge_pt),
            "tip": list(tip),
            "arc_pts": [list(p) for p in arc_pts],
        }
        if item.get("stacked"):
            entry["stacked"] = True
        result.append(entry)

    return result


_FACE_MAP = {"south": (0, 1), "east": (1, 2), "north": (2, 3), "west": (3, 0)}


def _face_midpoint(poly, face):
    """Midpoint of a named face of a 4-vertex polygon.

    Convention: south=poly[0]→poly[1], east=poly[1]→poly[2],
    north=poly[2]→poly[3], west=poly[3]→poly[0].
    poly vertices are [E, N] lists or tuples.
    """
    indices = _FACE_MAP.get(face)
    if not indices or len(poly) < 4:
        return None
    i, j = indices
    a, b = poly[i], poly[j]
    return [(a[0] + b[0]) / 2, (a[1] + b[1]) / 2]


def _face_point_at(poly, face, dist_in):
    """Point at dist_in inches from the near corner of a named face of a 4-vertex polygon."""
    indices = _FACE_MAP.get(face)
    if not indices or len(poly) < 4:
        return None
    i, j = indices
    a, b = poly[i], poly[j]
    dx, dy = b[0] - a[0], b[1] - a[1]
    length = math.hypot(dx, dy)
    if length < 1e-9:
        return None
    d = dist_in / 12.0
    return [a[0] + d * dx / length, a[1] + d * dy / length]


def _face_vertices(poly, face):
    """Return (vertex_i, vertex_j) for a named face of a 4-vertex polygon."""
    indices = _FACE_MAP.get(face)
    if not indices or len(poly) < 4:
        return None, None
    i, j = indices
    return poly[i], poly[j]


def _find_wall_poly(geom, name):
    """Look up a wall polygon by name from geometry result."""
    wall = geom.get("interior_walls", {}).get(name)
    return wall["poly"] if wall else None


def _find_opening_poly(geom, name):
    """Look up an opening polygon by name from geometry result."""
    for op in geom.get("outer_openings", []):
        if op["name"] == name and "poly" in op:
            return op["poly"]
    for ro in geom.get("rough_openings", []):
        if ro["name"] == name and "poly" in ro:
            return ro["poly"]
    return None


def _resolve_point_spec(spec, geom):
    """Resolve a point spec to [E, N] coordinates.

    Point specs:
      str            → geom["points"][spec]
      {face_mid: T, face: F}          → midpoint of wall face
      {opening_face_mid: T, face: F}  → midpoint of opening face
      {opening_centroid: T}            → centroid of opening poly vertices
      {midpoint: [A, B]}              → midpoint of two point specs
      {offset: base, dir: D, dist: N} → base + dist * direction
      {arc_point: {...}}              → point on arc at reference northing
    """
    if spec is None:
        return None
    if isinstance(spec, str):
        pt = geom.get("points", {}).get(spec)
        return list(pt) if pt else None

    if isinstance(spec, dict):
        if "face_mid" in spec:
            poly = _find_wall_poly(geom, spec["face_mid"])
            if poly:
                return _face_midpoint(poly, spec.get("face"))
            return None

        if "opening_face_mid" in spec:
            poly = _find_opening_poly(geom, spec["opening_face_mid"])
            if poly:
                return _face_midpoint(poly, spec.get("face"))
            return None

        if "opening_centroid" in spec:
            poly = _find_opening_poly(geom, spec["opening_centroid"])
            if poly and len(poly) >= 3:
                n = len(poly)
                return [sum(p[0] for p in poly) / n, sum(p[1] for p in poly) / n]
            return None

        if "midpoint" in spec:
            pair = spec["midpoint"]
            if len(pair) != 2:
                return None
            a = _resolve_point_spec(pair[0], geom)
            b = _resolve_point_spec(pair[1], geom)
            if a and b:
                return [(a[0] + b[0]) / 2, (a[1] + b[1]) / 2]
            return None

        if "offset" in spec:
            base = _resolve_point_spec(spec["offset"], geom)
            direction = _resolve_dir_spec(spec.get("dir"), geom)
            dist = spec.get("dist", 0)
            if base and direction:
                return [base[0] + dist * direction[0], base[1] + dist * direction[1]]
            return None

        if "arc_point" in spec:
            ap = spec["arc_point"]
            center = _resolve_point_spec(ap.get("center"), geom)
            radii = geom.get("radii", {})
            radius = radii.get(ap.get("radius_key"))
            ref = _resolve_point_spec(ap.get("reference"), geom)
            side = ap.get("side", "east")
            if center and radius and ref:
                # Find point on circle at reference's northing
                ns_dir = _resolve_dir_spec("north", geom)
                ew_dir = _resolve_dir_spec("east", geom)
                dn = (ref[0] - center[0]) * ns_dir[0] + (ref[1] - center[1]) * ns_dir[1]
                disc = radius**2 - dn**2
                if disc < 0:
                    return None
                de = math.sqrt(disc)
                if side == "west":
                    de = -de
                base = [center[0] + dn * ns_dir[0], center[1] + dn * ns_dir[1]]
                return [base[0] + de * ew_dir[0], base[1] + de * ew_dir[1]]
            return None

    return None


def _resolve_dir_spec(spec, geom):
    """Resolve a direction spec to [dE, dN] unit vector.

    Direction specs:
      "east"  → [1, 0]
      "north" → [0, 1]
      {face_along: T, face: F}   → along direction of wall face
      {face_perp: T, face: F}    → perpendicular of wall face (right-hand)
      {segment: [A, B]}          → unit vector from A to B
      {segment_perp: [A, B]}     → perpendicular of A→B (right-hand)
    """
    if spec is None:
        return None
    if isinstance(spec, str):
        if spec == "east":
            return [1.0, 0.0]
        if spec == "north":
            return [0.0, 1.0]
        return None

    if isinstance(spec, dict):
        if "face_along" in spec:
            poly = _find_wall_poly(geom, spec["face_along"])
            if not poly:
                return None
            a, b = _face_vertices(poly, spec.get("face"))
            if a is None:
                return None
            dx, dy = b[0] - a[0], b[1] - a[1]
            length = math.sqrt(dx * dx + dy * dy)
            if length < 1e-12:
                return None
            return [dx / length, dy / length]

        if "face_perp" in spec:
            poly = _find_wall_poly(geom, spec["face_perp"])
            if not poly:
                return None
            a, b = _face_vertices(poly, spec.get("face"))
            if a is None:
                return None
            dx, dy = b[0] - a[0], b[1] - a[1]
            length = math.sqrt(dx * dx + dy * dy)
            if length < 1e-12:
                return None
            # Right-hand perpendicular: (dy, -dx) / length
            return [dy / length, -dx / length]

        if "segment" in spec:
            pair = spec["segment"]
            if len(pair) != 2:
                return None
            a = _resolve_point_spec(pair[0], geom)
            b = _resolve_point_spec(pair[1], geom)
            if a and b:
                dx, dy = b[0] - a[0], b[1] - a[1]
                length = math.sqrt(dx * dx + dy * dy)
                if length < 1e-12:
                    return None
                return [dx / length, dy / length]
            return None

        if "segment_perp" in spec:
            pair = spec["segment_perp"]
            if len(pair) != 2:
                return None
            a = _resolve_point_spec(pair[0], geom)
            b = _resolve_point_spec(pair[1], geom)
            if a and b:
                dx, dy = b[0] - a[0], b[1] - a[1]
                length = math.sqrt(dx * dx + dy * dy)
                if length < 1e-12:
                    return None
                return [dy / length, -dx / length]
            return None

    return None


def _resolve_anchor(anchor, geometry_result):
    """Resolve a dimension anchor to [E, N] coordinates.

    Anchor types:
      point            → named geometry point
      wall_face        → midpoint of wall polygon face
      opening_face     → midpoint of opening polygon face
      line_intersection → intersection of two infinite lines
      computed         → arbitrary point via _resolve_point_spec
    Returns [E, N] or None if the target is not found.
    """
    if not anchor:
        return None
    atype = anchor.get("type")
    if not atype:
        return None

    if atype == "point":
        target = anchor.get("target")
        if not target:
            return None
        pt = geometry_result.get("points", {}).get(target)
        return list(pt) if pt else None

    if atype == "wall_face":
        target = anchor.get("target")
        face = anchor.get("face")
        if not target or not face:
            return None
        wall = geometry_result.get("interior_walls", {}).get(target)
        if wall:
            dist_in = anchor.get("distIn")
            if dist_in is not None:
                return _face_point_at(wall["poly"], face, dist_in)
            return _face_midpoint(wall["poly"], face)
        return None

    if atype == "opening_face":
        target = anchor.get("target")
        face = anchor.get("face")
        if not target or not face:
            return None
        for op in geometry_result.get("outer_openings", []):
            if op["name"] == target and "poly" in op:
                return _face_midpoint(op["poly"], face)
        for ro in geometry_result.get("rough_openings", []):
            if ro["name"] == target and "poly" in ro:
                return _face_midpoint(ro["poly"], face)
        return None

    if atype == "line_intersection":
        p1 = _resolve_point_spec(anchor.get("line1_point"), geometry_result)
        d1 = _resolve_dir_spec(anchor.get("line1_dir"), geometry_result)
        p2 = _resolve_point_spec(anchor.get("line2_point"), geometry_result)
        d2 = _resolve_dir_spec(anchor.get("line2_dir"), geometry_result)
        if p1 and d1 and p2 and d2:
            det = d1[0] * d2[1] - d1[1] * d2[0]
            if abs(det) < 1e-12:
                return None
            t = ((p2[0] - p1[0]) * d2[1] - (p2[1] - p1[1]) * d2[0]) / det
            return [p1[0] + t * d1[0], p1[1] + t * d1[1]]
        return None

    if atype == "computed":
        return _resolve_point_spec(anchor.get("spec"), geometry_result)

    return None


def _polygon_area(poly):
    """Signed area of a simple polygon via the shoelace formula."""
    n = len(poly)
    if n < 3:
        return 0.0
    area = 0.0
    for i in range(n):
        x0 = poly[i][0] if isinstance(poly[i], (list, tuple)) else poly[i][0]
        y0 = poly[i][1] if isinstance(poly[i], (list, tuple)) else poly[i][1]
        j = (i + 1) % n
        x1 = poly[j][0] if isinstance(poly[j], (list, tuple)) else poly[j][0]
        y1 = poly[j][1] if isinstance(poly[j], (list, tuple)) else poly[j][1]
        area += x0 * y1 - x1 * y0
    return abs(area) / 2.0


def _iwp(walls, name, idx):
    """Get interior wall poly vertex: walls["IW1"]["poly"][0] → (E, N) tuple."""
    w = walls.get(name)
    if not w:
        return (0.0, 0.0)
    p = w["poly"][idx]
    return (p[0], p[1]) if isinstance(p, (list, tuple)) else p


def _compute_room_labels(pts, walls, openings_result, inner_segs, variant,
                         db_path=None):
    """Compute room label positions, areas, and SF partition lines.

    Uses formula-evaluated wall positions (walls dict from
    result["interior_walls"]) and opening positions from the result dict,
    eliminating the dependency on the procedural layout object.

    Args:
        pts: dict of named points (W-series, F-series, etc.)
        walls: dict of wall name → {"poly": [...], "bbox": {...}}
        openings_result: dict with "outer_openings" and "rough_openings" lists
        inner_segs: list of inner wall segments
        variant: variant name
        db_path: database path

    Returns dict with:
      room_labels: [{name, pos, area?, poly?}]
      sf_lines: [{start, end}]  (sf variant only)
    """
    from shared.geometry import seg_vecs, line_isect, segment_polyline, offset_pt

    # Look up O6 and RO1 from the result openings
    o6_poly = None
    for op in openings_result.get("outer_openings", []):
        if op["name"] == "O6" and "poly" in op:
            o6_poly = op["poly"]
            break
    if o6_poly is None:
        return {"room_labels": []}
    o6_w = (o6_poly[0][0], o6_poly[0][1])

    ro1_poly = None
    for ro in openings_result.get("rough_openings", []):
        if ro["name"] == "RO1" and "poly" in ro:
            ro1_poly = ro["poly"]
            break
    if ro1_poly is None:
        return {"room_labels": []}
    ro1_w_nf = (ro1_poly[3][0], ro1_poly[3][1])

    w9w10_al, _ = seg_vecs(pts["W9"], pts["W10"])
    w2w5_al, _ = seg_vecs(pts["W2"], pts["W5"])
    w18w1_al, _ = seg_vecs(pts["W18"], pts["W1"])
    iw3_nw = _iwp(walls, "IW3", 3)
    iw3_w2w5 = line_isect(iw3_nw, w18w1_al, pts["W2"], w2w5_al)
    iw2s_1 = _iwp(walls, "IW2S", 1)
    iw2s_2 = _iwp(walls, "IW2S", 2)
    iw2s_e_al, _ = seg_vecs(iw2s_1, iw2s_2)
    iw2s_at_w9 = line_isect(iw2s_1, iw2s_e_al, pts["W9"], w9w10_al)

    is_sf = (variant == "sf")

    # --- Build room polygons from formula-evaluated wall positions ---
    rooms = {}

    # BEDROOM
    rooms["BEDROOM"] = [
        (_iwp(walls, "IW9", 2)[0], _iwp(walls, "IW1", 0)[1]),
        (_iwp(walls, "IW11", 3)[0], _iwp(walls, "IW1", 0)[1]),
        (_iwp(walls, "IW11", 3)[0], pts["W1"][1]),
        (_iwp(walls, "IW9", 2)[0], pts["W1"][1]),
    ]

    # UTIL_N
    rooms["UTIL_N"] = [
        iw3_w2w5,
        _iwp(walls, "IW8", 0), _iwp(walls, "IW8", 1),
        _iwp(walls, "IW1", 0), _iwp(walls, "IW9", 3),
        (_iwp(walls, "IW9", 0)[0], _iwp(walls, "IW7", 2)[1]),
        _iwp(walls, "IW7", 3), iw3_nw,
    ]

    # UTIL_S
    _util_s = [iw3_nw, (_iwp(walls, "IW3", 0)[0], pts["W1"][1]), pts["W1"]]
    _util_s.extend(segment_polyline(inner_segs[0], pts)[1:])
    _util_s.append(iw3_w2w5)
    rooms["UTIL_S"] = _util_s

    # KITCHEN
    rooms["KITCHEN"] = [
        o6_w, iw2s_at_w9,
        iw2s_1, _iwp(walls, "IW2O", 3),
        _iwp(walls, "IW2O", 0), _iwp(walls, "IW2", 2),
        _iwp(walls, "IW2", 1), ro1_w_nf,
    ]

    # LIVING
    _living = [o6_w]
    _living.append(segment_polyline(inner_segs[6], pts)[-1])
    for si in range(7, 13):
        _living.extend(segment_polyline(inner_segs[si], pts)[1:])
    _living.append(_iwp(walls, "IW1", 2))
    _living.append(ro1_w_nf)
    rooms["LIVING"] = _living

    # BATH
    _seg2 = segment_polyline(inner_segs[2], pts)
    _seg3 = segment_polyline(inner_segs[3], pts)
    _bath = [
        _iwp(walls, "IW8", 3), _iwp(walls, "IW8", 2),
        _iwp(walls, "IW2", 3), _iwp(walls, "IW2O", 1),
        _iwp(walls, "IW2O", 2), _iwp(walls, "IW2S", 0),
        _iwp(walls, "IW2S", 3), _seg3[0],
    ]
    _bath.extend(reversed(_seg2[:-1]))
    rooms["BATH"] = _bath

    # OFFICE
    _office = [_iwp(walls, "IW5", 0)]
    _office.append((pts["W15"][0], _iwp(walls, "IW5", 0)[1]))
    _office.append(pts["W15"])
    for si in [14, 15, 16]:
        _office.extend(segment_polyline(inner_segs[si], pts)[1:])
    _office.append(_iwp(walls, "IW4", 1))
    _office.append(_iwp(walls, "IW4", 2))
    _office.append(_iwp(walls, "IW12", 2))
    _office.append(_iwp(walls, "IW12", 3))
    rooms["OFFICE"] = _office

    # E CLOSET
    rooms["E CLOSET"] = [
        (_iwp(walls, "IW11", 1)[0], _iwp(walls, "IW12", 0)[1]),
        (_iwp(walls, "IW4", 0)[0], _iwp(walls, "IW12", 0)[1]),
        _iwp(walls, "IW4", 0),
        (_iwp(walls, "IW11", 1)[0], pts["W1"][1]),
    ]

    # W CLOSET
    rooms["W CLOSET"] = [
        (_iwp(walls, "IW3", 1)[0], _iwp(walls, "IW7", 0)[1]),
        (_iwp(walls, "IW9", 0)[0], _iwp(walls, "IW7", 0)[1]),
        (_iwp(walls, "IW9", 0)[0], pts["W1"][1]),
        (_iwp(walls, "IW3", 1)[0], pts["W1"][1]),
    ]

    # STORAGE
    rooms["STORAGE"] = [
        (_iwp(walls, "IW11", 1)[0], _iwp(walls, "IW5", 3)[1]),
        (pts["W14"][0], _iwp(walls, "IW5", 3)[1]),
        (pts["W14"][0], _iwp(walls, "IW1", 0)[1]),
        (_iwp(walls, "IW11", 1)[0], _iwp(walls, "IW1", 0)[1]),
    ]

    # WH
    _seg3b = segment_polyline(inner_segs[3], pts)
    _seg4 = segment_polyline(inner_segs[4], pts)
    _seg5 = segment_polyline(inner_segs[5], pts)
    _wh = [iw2s_2, _seg3b[-1]]
    _wh.extend(_seg4[1:])
    _wh.extend(_seg5[1:])
    _wh.append((iw2s_2[0], pts["W9"][1]))
    rooms["WH"] = _wh

    # --- Compute areas directly via shoelace (replaces compute_room_areas) ---
    areas = {}
    for name, poly in rooms.items():
        areas[name] = _polygon_area(poly)
    # BATH: subtract IW6 area (IW6 sits entirely inside the bath room)
    iw6 = walls.get("IW6")
    if iw6 and "BATH" in areas:
        areas["BATH"] -= _polygon_area(iw6["poly"])

    # --- Build label list with wall-relative positions (matching reference) ---
    # Reference uses specific wall/opening-relative positions, not centroids.
    # Compute those positions here, falling back to centroids for unknown rooms.
    label_positions = {}

    # BEDROOM: right edge at W end of RO1, top at N end of RO3
    ro3_poly = None
    for ro in openings_result.get("rough_openings", []):
        if ro["name"] == "RO3" and "poly" in ro:
            ro3_poly = ro["poly"]
            break
    ro1_w_mid = ((ro1_poly[0][0] + ro1_poly[3][0]) / 2,
                 (ro1_poly[0][1] + ro1_poly[3][1]) / 2)
    if ro3_poly:
        ro3_n_mid = ((ro3_poly[2][0] + ro3_poly[3][0]) / 2,
                     (ro3_poly[2][1] + ro3_poly[3][1]) / 2)
        label_positions["BEDROOM"] = (ro1_w_mid[0], ro3_n_mid[1])

        # UTIL_N: same northing as BEDROOM, easting = midpoint of
        # south toilet center and util sink center (matching reference)
        vi = openings_result.get("variant_items", {})
        toilet_item = vi.get("toilet_s", {})
        sink_item = vi.get("util_sink", {})
        if toilet_item.get("poly") and sink_item.get("poly"):
            tp = toilet_item["poly"]
            sp = sink_item["poly"]
            toilet_cx = sum(p[0] for p in tp) / len(tp)
            sink_cx = sum(p[0] for p in sp) / len(sp)
            util_e = (toilet_cx + sink_cx) / 2
        else:
            util_e = _centroid(rooms.get("UTIL_N", []))[0]
        label_positions["UTIL_N"] = (util_e, ro3_n_mid[1])

    # KITCHEN: centered beneath kitchen sink, just above dim02 line
    # dim02 line is at F12-F13 northing; label is 3" above that
    if "F12" in pts and "F13" in pts:
        from floorplan.constants import (NORTH_CTR_LENGTH, KITCHEN_APPL_GAP,
                                         STOVE_WIDTH, KITCHEN_SINK_WIDTH)
        dim02_n = (pts["F12"][1] + pts["F13"][1]) / 2
        iw2s_ne = _iwp(walls, "IW2S", 2)
        iw2_d_k = ((iw2s_ne[0] - pts["W9"][0]) * w9w10_al[0] +
                    (iw2s_ne[1] - pts["W9"][1]) * w9w10_al[1])
        st_d_k = iw2_d_k + NORTH_CTR_LENGTH + KITCHEN_APPL_GAP
        ks_d_k = st_d_k + STOVE_WIDTH + KITCHEN_APPL_GAP + 2.0 / 12.0
        sink_ctr = offset_pt(pts["W9"], ks_d_k + KITCHEN_SINK_WIDTH / 2, w9w10_al)
        label_positions["KITCHEN"] = (sink_ctr[0], dim02_n + 3.0 / 12.0)

        # LIVING: centered under O6 at KITCHEN label's northing
        o6_cx = sum(p[0] for p in o6_poly) / len(o6_poly)
        label_positions["LIVING"] = (o6_cx, dim02_n + 3.0 / 12.0)

    # BATH: centered between IW2s west face and W2-W5, at RO4 north end northing
    ro4_poly = None
    for ro in openings_result.get("rough_openings", []):
        if ro["name"] == "RO4" and "poly" in ro:
            ro4_poly = ro["poly"]
            break
    if ro4_poly:
        iw2s_0 = _iwp(walls, "IW2S", 0)
        bath_e = (iw2s_0[0] + pts["W2"][0]) / 2
        ro4_n_mid = ((ro4_poly[2][0] + ro4_poly[3][0]) / 2,
                     (ro4_poly[2][1] + ro4_poly[3][1]) / 2)
        label_positions["BATH"] = (bath_e, ro4_n_mid[1])

    # OFFICE: midpoint between IW4 east face and W15 (easting),
    # N-S: between ctr+5'+WALL_3IN and IW1, adjusted by -2'+26"
    from floorplan.constants import WALL_3IN
    iw4_e_mid = ((_iwp(walls, "IW4", 1)[0] + _iwp(walls, "IW4", 2)[0]) / 2,
                 (_iwp(walls, "IW4", 1)[1] + _iwp(walls, "IW4", 2)[1]) / 2)
    of_ew = ((iw4_e_mid[0] + pts["W15"][0]) / 2,
             (iw4_e_mid[1] + pts["W15"][1]) / 2)
    iw1_n_al, _iw1_n_cw = seg_vecs(
        _iwp(walls, "IW1", 3), _iwp(walls, "IW1", 2))
    iw1_n_out = (-_iw1_n_cw[0], -_iw1_n_cw[1])  # outward = north
    iw1_s_mid = ((_iwp(walls, "IW1", 0)[0] + _iwp(walls, "IW1", 1)[0]) / 2,
                 (_iwp(walls, "IW1", 0)[1] + _iwp(walls, "IW1", 1)[1]) / 2)
    # Counter south face midpoint (closet counter from variant_items)
    ctr_item_of = vi.get("counter", {}) if vi else {}
    ctr_of_poly = ctr_item_of.get("poly", [])
    if len(ctr_of_poly) >= 2:
        ctr_s_mid = ((ctr_of_poly[0][0] + ctr_of_poly[1][0]) / 2,
                     (ctr_of_poly[0][1] + ctr_of_poly[1][1]) / 2)
    else:
        ctr_s_mid = iw1_s_mid
    ctr_offset = offset_pt(ctr_s_mid, 5.0 + WALL_3IN, iw1_n_out)
    of_ns = ((ctr_offset[0] + iw1_s_mid[0]) / 2,
             (ctr_offset[1] + iw1_s_mid[1]) / 2)
    of_ns_adj = offset_pt(of_ns, -2.0 + 26.0 / 12.0, iw1_n_out)
    of_ns_d = ((of_ns_adj[0] - of_ew[0]) * iw1_n_out[0] +
               (of_ns_adj[1] - of_ew[1]) * iw1_n_out[1])
    of_cx = of_ew[0] + of_ns_d * iw1_n_out[0]
    of_cy = of_ew[1] + of_ns_d * iw1_n_out[1]
    label_positions["OFFICE"] = (of_cx, of_cy)

    # Apply DB-stored offsets on top of computed positions
    label_offsets = {}
    try:
        all_elems = get_all_elements(db_path)
        for e in all_elems:
            if e["type"] == "label":
                props = json.loads(e["properties"]) if isinstance(e["properties"], str) else e["properties"]
                if props.get("source") == "room":
                    label_offsets[e["name"]] = (
                        props.get("offset_e", 0.0),
                        props.get("offset_n", 0.0),
                    )
    except Exception:
        pass
    if not label_offsets:
        label_offsets = get_room_label_offsets()

    labels = []
    for name, poly in rooms.items():
        cx, cy = _centroid(poly)
        # Use wall-relative position if available, else centroid
        if name in label_positions:
            lx, ly = label_positions[name]
        else:
            lx, ly = cx, cy
        de, dn = label_offsets.get(name, (0.0, 0.0))
        lbl = {
            "name": name,
            "pos": point_to_list((lx + de, ly + dn)),
            "centroid": point_to_list((cx, cy)),
        }
        lbl["area"] = round(areas[name], 1)
        if is_sf:
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
    """Build outline_segs from solved chain — delegated to gen_provider."""
    from app.gen_provider import _build_outline_segs_from_chain as _impl
    return _impl(chain)


def _derive_constant(constants_dict, name):
    """Compute a derived constant — delegated to gen_provider."""
    from app.gen_provider import _derive_constant as _impl
    return _impl(constants_dict, name)


def _build_elements_from_formulas(ev, variant, exclusions, db_path):
    """Build result dicts from FormulaEvaluator output + DB element metadata.

    Queries the elements table for metadata (label, type, shape, door,
    clearance, etc.) and combines with computed geometry from ev.elements.

    Returns (interior_walls, outer_openings, rough_openings, variant_items,
             furniture, appliances) tuple.
    """
    from app.variants import get_variant_flags

    flags = get_variant_flags(variant, db_path)
    minik = flags.get("minik", False)
    db_flag = flags.get("db", False)

    # Load element metadata from DB
    all_elems = get_all_elements(db_path)
    meta_by_name = {}
    for e in all_elems:
        props = json.loads(e["properties"]) if isinstance(e["properties"], str) else (e["properties"] or {})
        meta_by_name[e["name"]] = {
            "type": e["type"],
            "variant": e.get("variant"),
            "props": props,
        }

    excluded_walls = exclusions.get("wall", set())
    excluded_openings = exclusions.get("rough_opening", set())

    interior_walls = {}
    outer_openings = []
    rough_openings = []
    variant_items = {}
    furniture = {}
    appliances = {}

    for elem_name, computed in ev.elements.items():
        if "poly" not in computed:
            continue

        poly = computed["poly"]
        bbox = computed.get("bbox", {})

        meta = meta_by_name.get(elem_name, {})
        props = meta.get("props", {})
        item_name = elem_name
        elem_type = meta.get("type", "")

        # Interior walls
        if elem_type == "wall":
            if elem_name in excluded_walls:
                continue
            interior_walls[elem_name] = {"poly": poly, "bbox": bbox}

        # Outer openings (O-prefixed, e.g. O1-O11 seed + user-created O12+)
        elif elem_type == "opening" and elem_name.startswith("O") and not elem_name.startswith("O_"):
            entry = {
                "name": elem_name,
                "seg_start": props.get("seg_start", ""),
                "seg_end": props.get("seg_end", ""),
                "poly": poly,
            }
            if props.get("opening_type"):
                entry["opening_type"] = props["opening_type"]
            if props.get("product_url"):
                entry["product_url"] = props["product_url"]
            outer_openings.append(entry)

        # Rough openings (RO1-RO7)
        elif elem_type == "opening" and elem_name.startswith("RO"):
            if elem_name in excluded_openings:
                continue
            rough_openings.append({
                "name": elem_name,
                "bbox": bbox,
                "wall_name": props.get("wall_name", ""),
                "orientation": props.get("orientation", ""),
                "width": bbox.get("e", 0) - bbox.get("w", 0)
                         if props.get("orientation") == "H"
                         else bbox.get("n", 0) - bbox.get("s", 0),
                "poly": poly,
            })

        # Variant items (furniture, appliances, fixtures)
        elif elem_type in ("appliance", "fixture", "furniture"):
            # Variant filtering: check if item belongs to current variant
            variants_list = props.get("variants")
            if variants_list and variant not in variants_list:
                continue

            entry = {
                "name": item_name,
                "type": elem_type,
                "poly": [[p[0], p[1]] for p in poly] if poly else [],
                "bbox": bbox,
                "label": props.get("label", item_name.upper()),
                "shape": props.get("shape", "rect"),
            }
            if computed.get("center"):
                entry["center"] = computed["center"]
            if computed.get("pos_origin"):
                entry["pos_origin"] = computed["pos_origin"]
            if computed.get("width"):
                entry["width"] = computed["width"]
            if computed.get("depth"):
                entry["depth"] = computed["depth"]
            if computed.get("rotation") is not None:
                entry["rotation"] = computed["rotation"]
            if computed.get("radius"):
                entry["radius"] = computed["radius"]
            if props.get("stacked"):
                entry["stacked"] = True
            if props.get("door"):
                door_cfg = props["door"]
                # Resolve variant-keyed door configs: {"standard": {...}, ...}
                first_val = next(iter(door_cfg.values()), None)
                if isinstance(first_val, dict):
                    door_cfg = door_cfg.get(variant, {})
                if door_cfg:
                    if props.get("door_flipped"):
                        door_cfg = dict(door_cfg)
                        door_cfg["hinge_idx"], door_cfg["target_idx"] = (
                            door_cfg["target_idx"], door_cfg["hinge_idx"]
                        )
                    entry["door"] = door_cfg
            if props.get("clearance"):
                entry["clearance"] = props["clearance"]
            if props.get("clip_to_inner"):
                entry["clip_to_inner"] = True

            # Product URL from DB properties
            raw_url = props.get("product_url")
            if raw_url:
                if isinstance(raw_url, dict):
                    # Variant-keyed: {"minik": "...", "default": "..."}
                    if minik and "minik" in raw_url:
                        entry["product_url"] = raw_url["minik"]
                    elif db_flag and "db" in raw_url:
                        entry["product_url"] = raw_url["db"]
                    elif "default" in raw_url:
                        entry["product_url"] = raw_url["default"]
                else:
                    entry["product_url"] = raw_url

            variant_items[item_name] = entry

            # Backward-compatible furniture/appliances dicts
            compat_entry = {"poly": poly, "bbox": bbox}
            if elem_type == "furniture":
                furniture[item_name] = compat_entry
            elif elem_type == "appliance":
                appliances[item_name] = compat_entry

    return (interior_walls, outer_openings, rough_openings,
            variant_items, furniture, appliances)


def _compute_traverse_from_db(db_path=None):
    """Compute traverse from DB survey data — delegated to gen_provider."""
    from app.gen_provider import _compute_traverse_from_db as _impl
    return _impl(db_path)


def compute_geometry(constants_dict: dict, variant: str = "standard",
                     chain_rows: list[dict] | None = None,
                     doors_data: list[dict] | None = None,
                     db_path: str | None = None) -> dict:
    """Compute all building geometry from constants and return JSON-serialisable dict.

    Phase 12h: fully formula-driven architecture.  The FormulaEvaluator is
    the sole source for all element geometry (interior walls, openings,
    furniture, appliances).  Element metadata (labels, door configs,
    clearance, URLs) comes from the elements table in the database.

    If chain_rows is provided, the app solver computes the outline from
    DB chain data (the normal path).

    If doors_data is provided, door swing arcs are computed from opening
    polygons and door configurations.
    """
    from shared.geometry import path_polygon
    from app.evaluator import FormulaEvaluator
    from app.database import (get_all_formulas, get_variants as _get_db_variants,
                              get_inner_wall_overrides)
    from app.gen_provider import compute_native_geometry, apply_overrides_to_poly

    # 1-3. Survey traverse, F-series outline, inner walls (shared with GeneratorData)
    pts, outline_segs, inner_segs, radii = compute_native_geometry(
        constants_dict, chain_rows=chain_rows, db_path=db_path)
    inner_poly = path_polygon(inner_segs, pts)

    # 3b. Apply inner wall overrides (W-series geometry exceptions)
    overrides = get_inner_wall_overrides(db_path)
    if overrides:
        apply_overrides_to_poly(inner_poly, inner_segs, pts, overrides)

    # 4. Load variant exclusions
    exclusions = get_variant_exclusions(variant, db_path)

    # 5. FormulaEvaluator — sole source for all element geometry
    base_points = {k: point_to_list(v) for k, v in pts.items()}
    inner = [(p[0], p[1]) for p in inner_poly]
    ev = FormulaEvaluator(constants_dict, base_points, inner, radii)
    ev.load_formulas_from_db(db_path=db_path, variant=variant)
    ev.evaluate_all()

    # Locked element names for canvas rendering
    formulas = get_all_formulas(variant=variant, db_path=db_path)
    locked = sorted({row["element_name"] for row in formulas if row.get("locked")})

    # 6. Build element result dicts from formula output + DB metadata
    (iw, oo, ro, vi, furn, appl) = _build_elements_from_formulas(
        ev, variant, exclusions, db_path)

    # 6a. Include placed walls from elements table (no formulas, geometry
    #     stored directly in properties.poly — e.g., user-drawn CW walls)
    excluded_walls = exclusions.get("wall", set())
    all_elems = get_all_elements(db_path)
    for e in all_elems:
        if e["type"] != "wall" or e["name"] in iw or e["name"] in excluded_walls:
            continue
        props = json.loads(e["properties"]) if isinstance(e["properties"], str) else (e["properties"] or {})
        poly = props.get("poly")
        if poly and len(poly) >= 3:
            bbox = bbox_from_poly([[p[0], p[1]] for p in poly])
            iw[e["name"]] = {"poly": poly, "bbox": bbox}

    # 6b. Site path (survey outer boundary as dense polyline)
    site_path = _compute_site_path(pts)

    # 6c. Roof outline polygon (DB-driven if roof corners configured)
    roof_poly_pts = []
    try:
        from app.database import get_roof_corners
        roof_data = get_roof_corners(db_path)
        if roof_data and roof_data.get("corners"):
            from shared.roof_outline import compute_db_roof_outline
            corner_names = [c["center"] for c in roof_data["corners"]]
            corner_radiused = [c.get("radiused", False) for c in roof_data["corners"]]
            overhang = float(roof_data.get("overhang", 0.5))
            result_r = compute_db_roof_outline(
                corner_names, corner_radiused, pts, radii, overhang)
            roof_poly_pts = [point_to_list(p) for p in result_r.poly]
    except Exception:
        pass

    # 7. Build result
    outline_poly_pts = path_polygon(outline_segs, pts)
    result = {
        "points": {name: point_to_list(pt) for name, pt in sorted(pts.items())},
        "outline_segments": [seg_to_dict(s) for s in outline_segs],
        "inner_segments": [seg_to_dict(s) for s in inner_segs],
        "outline_poly": [point_to_list(p) for p in outline_poly_pts],
        "inner_poly": [point_to_list(p) for p in inner_poly],
        "interior_walls": iw,
        "outer_openings": oo,
        "rough_openings": ro,
        "appliances": appl,
        "furniture": furn,
        "variant_items": vi,
        "variant": variant,
        "locked_elements": locked,
        "radii": dict(radii),
        "site_path": site_path,
        "roof_poly": roof_poly_pts,
    }

    # Bounding box
    if result["outline_poly"]:
        result["bbox"] = bbox_from_poly(result["outline_poly"])

    # Available variants
    try:
        result["available_variants"] = [v["name"] for v in _get_db_variants(db_path)]
    except Exception:
        from app.variants import VARIANTS
        result["available_variants"] = list(VARIANTS.keys())

    # Room labels and SF extras
    room_data = _compute_room_labels(pts, result["interior_walls"],
                                      result, inner_segs, variant, db_path)
    result["room_labels"] = room_data["room_labels"]
    if "sf_lines" in room_data:
        result["sf_lines"] = room_data["sf_lines"]

    # Door arcs — uses formula-evaluated opening polygons (now dict-based)
    result["door_arcs"] = _compute_door_arcs(
        oo, ro, doors_data or [], exclusions)

    # Clearance zones — all from variant_items metadata
    result["clearance_zones"] = _compute_clearance_zones(
        variant, result.get("variant_items"))

    # Appliance door arcs — intrinsic door configs from DB metadata
    result["appliance_doors"] = _compute_appliance_doors(
        result.get("variant_items"), variant)

    # Dimension and label elements
    all_elements = get_all_elements(db_path)
    excluded_dims = exclusions.get("dimension", set())
    user_dims = []
    label_elems = []
    for e in all_elements:
        props = json.loads(e["properties"]) if isinstance(e["properties"], str) else e["properties"]
        variants_list = props.get("variants")
        if variants_list is not None:
            if variant not in variants_list:
                continue
        else:
            elem_variant = e.get("variant")
            if elem_variant is not None and elem_variant != variant:
                continue
        if e["type"] == "dimension":
            if e["name"] in excluded_dims:
                continue
            for anchor_key, coord_key in [("start_anchor", "start"), ("end_anchor", "end")]:
                anchor = props.get(anchor_key)
                if anchor:
                    resolved = _resolve_anchor(anchor, result)
                    if resolved:
                        props[coord_key] = resolved
            user_dims.append({
                "id": e["id"], "name": e["name"], "properties": props,
            })
        elif e["type"] == "label":
            entry = {"id": e["id"], "name": e["name"], "properties": props}
            if props.get("source") == "room":
                rl = next((r for r in result["room_labels"] if r["name"] == e["name"]), None)
                if rl:
                    entry["centroid"] = rl["centroid"]
                    entry["pos"] = rl["pos"]
            label_elems.append(entry)
    result["user_dimensions"] = user_dims
    result["label_elements"] = label_elems

    return result


def _compute_site_path(pts: dict) -> list[list[float]]:
    """Compute site boundary (survey outer path) as a dense polyline.

    Uses the P-series traverse points already present in pts (aligned to
    F-series coordinates) to build the outer boundary segments, then
    returns a polygon approximation as [[E, N], ...].
    """
    from shared.types import LineSeg, ArcSeg
    from shared.geometry import path_polygon
    from shared.survey import compute_traverse, compute_three_arc, compute_pt1
    from floorplan.geometry import align_pts_to_f_series

    trav_pts = compute_traverse()
    three_arc = compute_three_arc(trav_pts)
    R1, R2, R3 = three_arc["R1"], three_arc["R2"], three_arc["R3"]
    align_pts_to_f_series(trav_pts)
    trav_pts["PT1"] = compute_pt1(trav_pts, R1)

    outer_segs = [
        LineSeg("POB", "P2"), LineSeg("P2", "P3"), LineSeg("P3", "T3"),
        ArcSeg("T3", "PX", "TC3", R3, "CW", 60),
        LineSeg("PX", "P4"), LineSeg("P4", "P5"), LineSeg("P5", "PT1"),
        ArcSeg("PT1", "T1", "TC1", R1, "CW", 60),
        ArcSeg("T1", "PA", "TC1", R1, "CW", 60),
        ArcSeg("PA", "T2", "TC2", R2, "CW", 60),
        LineSeg("T2", "POB"),
    ]
    poly = path_polygon(outer_segs, trav_pts)
    return [point_to_list(p) for p in poly]


def compute_survey_points(constants_dict: dict) -> dict:
    """Return P-series survey points and inter-point distances."""
    from shared.survey import compute_traverse, compute_three_arc, compute_pt1
    pts = compute_traverse()
    three_arc = compute_three_arc(pts)
    pts["PT1"] = compute_pt1(pts, three_arc["R1"])

    traverse_order = ["POB", "P2", "P3", "P4", "P5"]
    arc_point_names = ["T1", "T2", "T3", "PA", "PX", "TC1", "TC2", "TC3", "PT1"]

    points = {}
    for name in traverse_order + arc_point_names:
        if name in pts:
            points[name] = point_to_list(pts[name])

    distances = []
    for i in range(len(traverse_order)):
        a = traverse_order[i]
        b = traverse_order[(i + 1) % len(traverse_order)]
        if a in pts and b in pts:
            d = math.hypot(pts[b][0] - pts[a][0], pts[b][1] - pts[a][1])
            distances.append({"from": a, "to": b, "distance": round(d, 4)})

    return {
        "points": points,
        "distances": distances,
        "arc_radii": {
            "R1": three_arc["R1"],
            "R2": three_arc["R2"],
            "R3": three_arc["R3"],
        },
    }


def generate_svg(view_name: str, script_path: str) -> bool:
    """Run a generator script as a subprocess and return True on success.

    This is the standalone path — generators construct their own
    GeneratorData from hardcoded procedural modules.  For DB-driven
    generation, use generate_svg_db() instead.
    """
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


def _run_generator_inprocess(script_path: str, gd, db_path: str = None) -> bool:
    """Run a generator in-process using DB-sourced GeneratorData.

    Returns True on success, False on error, or None if no in-process
    handler exists for the given script (caller should fall back to
    subprocess).
    """
    _has_f_series = "F1" in gd.pts and "F2" in gd.pts

    if script_path == "floorplan/gen_floorplan.py":
        from floorplan.gen_floorplan import build_floorplan_data
        from plumbing.gen_plumbing import _compute_boundary_corners
        from app.db_render import render_floorplan_svg_db
        from app.database import get_all_doors
        from app.plumbing import get_plumbing_elements

        data = build_floorplan_data(gd)
        try:
            boundary = _compute_boundary_corners(data.pts)
        except (KeyError, TypeError):
            boundary = None  # F-series points absent (non-standard chain)
        base = os.path.join(_PROJECT, "floorplan")

        # Variant configs: (suffix, variant_name, room_title)
        _VARIANTS = [
            ("",          "standard",  "Parent Suite"),
            ("_minik",    "minik",     "Parent Suite w/Small Kitchen"),
            ("_db",       "daybed",    "Parent Suite with Daybed"),
            ("_bare",     "bare",      "Room Dimensions"),
            ("_sf",       "sf",        "Room Dimensions"),
            ("_plumbing", "plumbing",  "Plumbing Plan"),
        ]

        constants_dict = gd.constants
        doors_data = get_all_doors(db_path) if db_path else []

        # Load chain_rows from DB so compute_geometry uses DB-driven outline
        # geometry (F-series points) rather than hardcoded floorplan/constants.py.
        # Without this, dimension anchor resolution uses stale F-series positions.
        from app.database import get_outline_chain
        chain_rows = get_outline_chain(db_path) if db_path else None

        # Fetch fixture connections once for plumbing variant
        fc_elems = [e for e in get_plumbing_elements(db_path)
                    if e["type"] == "fixture_connection"]

        for suffix, variant, room_title in _VARIANTS:
            geom = compute_geometry(
                constants_dict, variant=variant,
                chain_rows=chain_rows,
                doors_data=doors_data, db_path=db_path)
            bdy = boundary if variant == "plumbing" else None
            if variant == "plumbing":
                geom["fixture_connections"] = fc_elems
            svg = render_floorplan_svg_db(
                geom, data, room_title=room_title, boundary=bdy)
            with open(os.path.join(base, f"floorplan{suffix}.svg"), "w", encoding="utf-8") as f:
                f.write(svg)
        return True

    if script_path == "roof/gen_roof.py":
        _has_db_roof = bool(getattr(gd, '_roof_corners_data', None) and
                            gd._roof_corners_data.get("corners"))
        if not _has_f_series and not _has_db_roof:
            return True  # roof outline requires F-series chain or DB corners
        from roof.gen_roof import render_roof_svg
        svg = render_roof_svg(gd)
        with open(os.path.join(_PROJECT, "roof", "roof.svg"), "w", encoding="utf-8") as f:
            f.write(svg)
        return True

    if script_path == "walls/gen_walls.py":
        if not _has_f_series:
            return True  # walls detail requires F-series chain
        from walls.gen_walls import build_wall_data, render_walls_svg
        data = build_wall_data(gd)
        base = os.path.join(_PROJECT, "walls")
        svg = render_walls_svg(data, title="Outer Walls")
        with open(os.path.join(base, "walls.svg"), "w", encoding="utf-8") as f:
            f.write(svg)
        svg_all = render_walls_svg(data, title="All Walls", include_interior=True)
        with open(os.path.join(base, "all_walls.svg"), "w", encoding="utf-8") as f:
            f.write(svg_all)
        return True

    if script_path == "span/gen_span.py":
        if not _has_f_series:
            return True  # span analysis requires F-series chain
        from span.gen_span import _generate_svg as _span_svg
        from span._common import build_geometry
        pts, _, outer_poly, inner_poly, layout, roof_poly = build_geometry(gd)
        svg = _span_svg(pts, outer_poly, inner_poly, layout, roof_poly)
        with open(os.path.join(_PROJECT, "span", "span.svg"), "w", encoding="utf-8") as f:
            f.write(svg)
        return True

    if script_path == "span/gen_span_minmax.py":
        if not _has_f_series:
            return True  # span analysis requires F-series chain
        from span.gen_span_minmax import _generate_svg as _span_mm_svg
        from span._common import build_geometry
        pts, _, outer_poly, inner_poly, layout, roof_poly = build_geometry(gd)
        svg = _span_mm_svg(pts, outer_poly, inner_poly, layout, roof_poly)
        with open(os.path.join(_PROJECT, "span", "span_minmax.svg"), "w", encoding="utf-8") as f:
            f.write(svg)
        return True

    if script_path == "span/gen_span_min.py":
        if not _has_f_series:
            return True  # span analysis requires F-series chain
        from span.gen_span_min import _generate_svg as _span_min_svg
        from span._common import build_geometry
        pts, _, outer_poly, inner_poly, layout, roof_poly = build_geometry(gd)
        svg = _span_min_svg(pts, outer_poly, inner_poly, layout, roof_poly)
        with open(os.path.join(_PROJECT, "span", "span_min.svg"), "w", encoding="utf-8") as f:
            f.write(svg)
        return True

    if script_path == "site/gen_site_plan.py":
        # 'site' shadows stdlib module — use importlib
        import importlib.util
        _sp = importlib.util.spec_from_file_location(
            "site_gen_site_plan",
            os.path.join(_PROJECT, "site", "gen_site_plan.py"))
        _mod = importlib.util.module_from_spec(_sp)
        _sp.loader.exec_module(_mod)
        build_site_plan_data = _mod.build_site_plan_data
        render_site_plan = _mod.render_site_plan
        render_site_plan_df = _mod.render_site_plan_df
        render_site_plan_fs = _mod.render_site_plan_fs
        sp = build_site_plan_data(gd)
        base = os.path.join(_PROJECT, "site")
        doc = render_site_plan(sp)
        doc.save(os.path.join(base, "site_plan.pdf"))
        doc.close()
        doc_df = render_site_plan(sp, corners=False)
        render_site_plan_df(doc_df, sp)
        doc_df.save(os.path.join(base, "site_plan_df.pdf"))
        doc_df.close()
        doc_fs = render_site_plan(sp, corners=False)
        render_site_plan_df(doc_fs, sp)
        render_site_plan_fs(doc_fs, sp)
        doc_fs.save(os.path.join(base, "site_plan_fs.pdf"))
        doc_fs.close()
        return True

    if script_path in ("survey/gen_path_svg.py", "survey\\gen_path_svg.py"):
        import importlib.util
        _sp = importlib.util.spec_from_file_location(
            "survey_gen_path_svg",
            os.path.join(_PROJECT, "survey", "gen_path_svg.py"))
        _mod = importlib.util.module_from_spec(_sp)
        _sp.loader.exec_module(_mod)
        _data = _mod.compute_all()
        # Override with DB-sourced geometry.
        # pts.update() adds PA/PB/WA/WB/PC/FC from DB while keeping survey pts.
        # outline_segs/inner_segs replaced with DB chain (PA/PB naming).
        # generate_svg() detects non-F-series segs and uses build_outline_cfg_db.
        _data["pts"].update(gd.pts)
        _data["outline_segs"] = gd.outline_segs
        _data["inner_segs"] = gd.inner_segs
        _data["outer_poly"] = gd.outline_poly
        _data["inner_poly"] = gd.inner_poly
        _data["outline_area"] = gd.outer_area
        _out = os.path.join(_PROJECT, "survey", "path_area.svg")
        _mod.generate_svg(_data, iw_polys=gd.iw_polys, out_path=_out)
        return True

    if script_path == "scad/gen_flat_roof.py":
        from scad.gen_flat_roof import generate as gen_flat
        gen_flat(gd)
        return True

    if script_path == "scad/gen_2in12.py":
        from scad.gen_2in12 import generate as gen_2in12
        gen_2in12(gd)
        return True

    # No in-process handler — caller should fall back to subprocess
    return None


def generate_svg_db(view_name: str, script_path: str, gd,
                    db_path: str = None) -> bool:
    """Generate SVG using DB-sourced GeneratorData (in-process).

    Falls back to subprocess for scripts without in-process dispatch
    (e.g., gen_views.py, gen_line_drawings.py, gen_3views.py).
    """
    try:
        result = _run_generator_inprocess(script_path, gd, db_path=db_path)
        if result is not None:
            return result
    except Exception:
        return False
    # Fallback to subprocess for unhandled scripts
    return generate_svg(view_name, script_path)


def build_generator_data_from_db(db_path: str):
    """Build a GeneratorData from the current database state.

    This is the bridge between the app database and the generator
    data provider.  Used by the regeneration API to pass DB-driven
    geometry to generators.
    """
    from app.database import get_constants_dict, get_outline_chain
    from app.gen_provider import build_generator_data
    constants = get_constants_dict(db_path)
    chain_rows = get_outline_chain(db_path)
    return build_generator_data(constants, chain_rows=chain_rows,
                                db_path=db_path)


def generate_all_svgs(views: list[dict], gd=None) -> dict:
    """Regenerate all enabled views. Returns {name: success} dict.

    If gd (GeneratorData) is provided, uses in-process DB-driven
    generation.  Otherwise falls back to subprocess execution.
    """
    results = {}
    # Group by script to avoid running the same generator twice
    seen_scripts = set()
    for v in views:
        script = v["script"]
        if script in seen_scripts:
            results[v["name"]] = True
            continue
        seen_scripts.add(script)
        if gd is not None:
            results[v["name"]] = generate_svg_db(v["name"], script, gd)
        else:
            results[v["name"]] = generate_svg(v["name"], script)
    return results


def get_svg_content(svg_path: str) -> str | None:
    """Read an SVG file and return its content."""
    full_path = os.path.join(_PROJECT, svg_path)
    if not os.path.exists(full_path):
        return None
    with open(full_path, "r", encoding="utf-8") as f:
        return f.read()


# ---------------------------------------------------------------------------
# Span analysis (ANALYSIS-1, ANALYSIS-2)
# ---------------------------------------------------------------------------

def _extract_iw_centerlines_from_geo(geo):
    """Extract IW1/IW2/IW8 centerlines from compute_geometry() result.

    Same output format as span/_common.py:extract_iw_centerlines(layout).
    """
    iw = geo.get("interior_walls", {})
    cls = []
    for name in ("IW1", "IW2", "IW8"):
        wall = iw.get(name)
        if not wall:
            continue
        poly = wall["poly"]
        bbox = wall["bbox"]
        if name in ("IW1", "IW8"):
            # Horizontal wall: midline at vertical center
            mid_n = (bbox["s"] + bbox["n"]) / 2
            cl = ((poly[0][0], mid_n),
                  ((poly[1][0] + poly[2][0]) / 2, mid_n))
        else:
            # Vertical wall (IW2): midline at horizontal center
            mid_e = (bbox["w"] + bbox["e"]) / 2
            cl = ((mid_e, bbox["s"]), (mid_e, bbox["n"]))
        cls.append(cl)
    return cls


def _compute_spans_from_geo(inner_poly, geo):
    """Compute N-S spans using IW polygons from geometry result.

    Same logic as span/gen_span.py:_compute_spans but uses geo result
    instead of layout namedtuple.
    """
    from shared.geometry import vert_isects

    iw = geo.get("interior_walls", {})
    iw1_poly = iw["IW1"]["poly"] if "IW1" in iw else []
    iw8_poly = iw["IW8"]["poly"] if "IW8" in iw else []

    e_min = min(p[0] for p in inner_poly)
    e_max = max(p[0] for p in inner_poly)
    inch = 1.0 / 12.0
    eastings, spans, south_spans, north_spans = [], [], [], []
    e = e_min
    while e <= e_max + 1e-9:
        ns = vert_isects(inner_poly, e)
        if len(ns) >= 2:
            south_n = min(ns)
            north_n = max(ns)
            span = north_n - south_n
        else:
            span = south_n = north_n = 0.0

        spans.append(span)

        mid_n = None
        if iw1_poly:
            iw1_ns = vert_isects(iw1_poly, e)
            if len(iw1_ns) >= 2:
                mid_n = (min(iw1_ns) + max(iw1_ns)) / 2
        if iw8_poly:
            iw8_ns = vert_isects(iw8_poly, e)
            if len(iw8_ns) >= 2:
                iw8_mid = (min(iw8_ns) + max(iw8_ns)) / 2
                if mid_n is None or iw8_mid < mid_n:
                    mid_n = iw8_mid

        if mid_n is not None and span > 0:
            south_spans.append(mid_n - south_n)
            north_spans.append(north_n - mid_n)
        else:
            south_spans.append(span)
            north_spans.append(span)

        eastings.append(e)
        e += inch
    return eastings, spans, south_spans, north_spans


def compute_span_data(constants, db_path=None):
    """Return N-S span profile data for the current geometry.

    Returns dict with eastings, spans, south_spans, north_spans arrays.
    Uses compute_geometry() result (Phase 14-C), no module patching.
    """
    from app.database import get_outline_chain
    chain_rows = get_outline_chain(db_path) if db_path else None
    geo = compute_geometry(constants, chain_rows=chain_rows, db_path=db_path)
    inner_poly = [(p[0], p[1]) for p in geo["inner_poly"]]
    eastings, spans, south_spans, north_spans = _compute_spans_from_geo(
        inner_poly, geo)
    return {
        "eastings": eastings,
        "spans": spans,
        "south_spans": south_spans,
        "north_spans": north_spans,
    }


def compute_span_rotation(constants, db_path=None):
    """Return span-vs-rotation analysis with min/max.

    Returns dict with min_angle, min_span, max_angle, max_span, and
    data array of [angle, max_span] pairs at 5-degree steps.
    Uses compute_geometry() result (Phase 14-C), no module patching.
    """
    from span._common import max_span_at_angle, find_min_span_angle
    from app.database import get_outline_chain
    chain_rows = get_outline_chain(db_path) if db_path else None
    geo = compute_geometry(constants, chain_rows=chain_rows, db_path=db_path)
    inner_poly = [(p[0], p[1]) for p in geo["inner_poly"]]
    iw_cls = _extract_iw_centerlines_from_geo(geo)

    # Centroid for rotation
    cx = sum(p[0] for p in inner_poly) / len(inner_poly)
    cy = sum(p[1] for p in inner_poly) / len(inner_poly)

    # Sweep 5-175 degrees in 5-degree steps
    data = []
    for angle in range(5, 176, 5):
        ms = max_span_at_angle(inner_poly, iw_cls, angle, cx, cy)
        data.append([angle, round(ms, 4)])

    # Find precise min
    min_angle, min_span = find_min_span_angle(inner_poly, iw_cls, cx, cy,
                                               normalize=False)

    # Find max from sweep data
    max_entry = max(data, key=lambda d: d[1])

    return {
        "min_angle": round(min_angle, 1),
        "min_span": round(min_span, 4),
        "max_angle": max_entry[0],
        "max_span": max_entry[1],
        "data": data,
    }
