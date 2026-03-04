"""Geometry computation engine that reads from the database.

Patches floorplan.constants with DB values, runs existing computation
pipeline, and returns JSON-serialisable results.
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


def _compute_door_arcs(outer_openings, rough_openings, doors_data, exclusions):
    """Compute door swing arc geometry for all openings with doors.

    Returns list of door_arc dicts, each with opening_name, door_type,
    and leaves list (each leaf has hinge, tip, arc_pts).
    """
    if not doors_data:
        return []

    # Build lookup: opening_name → door config dict
    door_lookup = {d["opening_name"]: d for d in doors_data}

    # Gather all openings with their polys
    all_ops = []
    for op in outer_openings:
        all_ops.append({"name": op.name, "poly": op.poly})
    excluded_openings = exclusions.get("rough_opening", set())
    for ro in rough_openings:
        if ro.name in excluded_openings:
            continue
        if ro.poly:
            all_ops.append({"name": ro.name, "poly": ro.poly})
        elif ro.bbox:
            b = ro.bbox
            all_ops.append({"name": ro.name, "poly": [
                (b.w, b.s), (b.e, b.s), (b.e, b.n), (b.w, b.n)]})

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


def _compute_clearance_zones(layout, variant, variant_items=None):
    """Compute clearance zone polygons for fixture/furniture items.

    Returns list of clearance zone dicts with name, poly, style.
    Reads 'clearance' metadata from variant_items and also includes the
    dresser's hardcoded 15" clearance (matching gen_floorplan.py).
    """
    from shared.geometry import seg_vecs, offset_pt

    zones = []
    # Bare and SF variants have no furniture
    if variant in ("bare", "sf"):
        return zones

    dresser = getattr(layout, "dresser", None)
    if dresser:
        # Dresser clearance: 15" south from south face (poly[0]→poly[1] = SW→SE)
        al, outward = seg_vecs(dresser.poly[0], dresser.poly[1])
        cl_sw = offset_pt(dresser.poly[0], 15.0 / 12.0, outward)
        cl_se = offset_pt(dresser.poly[1], 15.0 / 12.0, outward)
        zones.append({
            "name": "dresser_clearance",
            "poly": [point_to_list(dresser.poly[0]),
                     point_to_list(dresser.poly[1]),
                     point_to_list(cl_se),
                     point_to_list(cl_sw)],
            "style": "dashed",
        })

    # Clearance zones from variant_items metadata
    if variant_items:
        for item_name, item in variant_items.items():
            cl = item.get("clearance")
            if not cl:
                continue
            face = cl["face"]  # [i, j] vertex indices
            dist = cl["distance"]
            poly = item["poly"]  # [[e,n], ...]
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
            cx = sum(p[0] for p in poly) / len(poly)
            cy = sum(p[1] for p in poly) / len(poly)
            mid = ((p_i[0] + p_j[0]) / 2, (p_i[1] + p_j[1]) / 2)
            to_center = (cx - mid[0], cy - mid[1])
            if perp[0] * to_center[0] + perp[1] * to_center[1] > 0:
                perp = (-perp[0], -perp[1])  # flip to point outward
            # Build extension polygon
            ext_i = (p_i[0] + dist * perp[0], p_i[1] + dist * perp[1])
            ext_j = (p_j[0] + dist * perp[0], p_j[1] + dist * perp[1])
            zones.append({
                "name": f"{item_name}_clearance",
                "poly": [list(p_i), list(p_j), list(ext_j), list(ext_i)],
                "style": "dashed",
            })

    return zones


def _compute_appliance_doors(variant_items):
    """Compute door swing arcs for appliances with 'door' metadata.

    Reads door metadata from variant_items (set in variants.py) and computes
    arc geometry using _swing_arc.  Returns list of appliance_door dicts.
    """
    if not variant_items:
        return []

    result = []
    for item_name, item in variant_items.items():
        door = item.get("door")
        if not door:
            continue
        poly = item["poly"]  # [[e,n], ...]
        hinge_idx = door["hinge_idx"]
        width = door["width"]
        open_dir = tuple(door["open_dir"])
        closed_dir = tuple(door["closed_dir"])

        hinge_pt = (poly[hinge_idx][0], poly[hinge_idx][1])
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


def _face_midpoint(poly, face):
    """Midpoint of a named face of a 4-vertex polygon.

    Convention: south=poly[0]→poly[1], east=poly[1]→poly[2],
    north=poly[2]→poly[3], west=poly[3]→poly[0].
    poly vertices are [E, N] lists or tuples.
    """
    face_map = {"south": (0, 1), "east": (1, 2), "north": (2, 3), "west": (3, 0)}
    indices = face_map.get(face)
    if not indices or len(poly) < 4:
        return None
    i, j = indices
    a, b = poly[i], poly[j]
    return [(a[0] + b[0]) / 2, (a[1] + b[1]) / 2]


def _resolve_anchor(anchor, geometry_result):
    """Resolve a dimension anchor to [E, N] coordinates.

    Returns [E, N] or None if the target is not found.
    """
    if not anchor:
        return None
    atype = anchor.get("type")
    target = anchor.get("target")
    if not atype or not target:
        return None

    if atype == "point":
        pt = geometry_result.get("points", {}).get(target)
        return list(pt) if pt else None

    if atype == "wall_face":
        face = anchor.get("face")
        wall = geometry_result.get("interior_walls", {}).get(target)
        if wall and face:
            return _face_midpoint(wall["poly"], face)
        return None

    if atype == "opening_face":
        face = anchor.get("face")
        for op in geometry_result.get("outer_openings", []):
            if op["name"] == target and "poly" in op:
                return _face_midpoint(op["poly"], face)
        for ro in geometry_result.get("rough_openings", []):
            if ro["name"] == target and "poly" in ro:
                return _face_midpoint(ro["poly"], face)
        return None

    return None


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

    # --- Build label list from polygon centroids + element offsets ---
    # Try label elements first (Phase 8), fall back to room_label_offsets
    label_offsets = {}
    try:
        all_elems = get_all_elements()
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
        de, dn = label_offsets.get(name, (0.0, 0.0))
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
                     chain_rows: list[dict] | None = None,
                     doors_data: list[dict] | None = None) -> dict:
    """Compute all building geometry from constants and return JSON-serialisable dict.

    If chain_rows is provided (Phase 5+), the app solver computes the
    outline from DB chain data, bypassing the module-scope solver in
    floorplan/geometry.py.

    If doors_data is provided (Phase 6+), door swing arcs are computed
    from opening polygons and door configurations.
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

    # Door arcs (Phase 6)
    result["door_arcs"] = _compute_door_arcs(
        outer_openings, rough_openings, doors_data or [], exclusions)

    # Clearance zones (Phase 6)
    result["clearance_zones"] = _compute_clearance_zones(layout, variant, variant_items)

    # Appliance door arcs (Phase 6)
    result["appliance_doors"] = _compute_appliance_doors(variant_items)

    # User dimensions and label elements (Phase 8)
    all_elements = get_all_elements()
    user_dims = []
    label_elems = []
    for e in all_elements:
        props = json.loads(e["properties"]) if isinstance(e["properties"], str) else e["properties"]
        if e["type"] == "dimension":
            # Resolve anchors to absolute coordinates
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
                # Merge centroid position from room_labels
                rl = next((r for r in result["room_labels"] if r["name"] == e["name"]), None)
                if rl:
                    entry["centroid"] = rl["centroid"]
                    entry["pos"] = rl["pos"]
            label_elems.append(entry)
    result["user_dimensions"] = user_dims
    result["label_elements"] = label_elems

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
