"""Graphical diagnostic: rasterize room polygons onto W-series boundary,
check for overlaps, gaps, and boundary violations."""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))

import numpy as np
from PIL import Image, ImageDraw, ImageFont
from shared.geometry import (poly_area, segment_polyline, seg_vecs,
                             line_isect, offset_pt, path_polygon)
from floorplan.gen_floorplan import build_floorplan_data, compute_iw_area
from floorplan.openings import compute_outer_openings, compute_rough_openings
from floorplan.constants import (SINK_RY, NORTH_CTR_LENGTH, KITCHEN_APPL_GAP,
                                  STOVE_WIDTH, KITCHEN_SINK_WIDTH, WALL_3IN)

# ---------------------------------------------------------------------------
# Build geometry
# ---------------------------------------------------------------------------
data = build_floorplan_data()
pts = data.pts
layout = data.layout

# ---------------------------------------------------------------------------
# Raster parameters: 100 pixels per foot for high resolution
# ---------------------------------------------------------------------------
PPF = 100  # pixels per foot
# Bounding box of inner polygon (with margin)
inner_poly = path_polygon(data.inner_segs, pts)
all_e = [p[0] for p in inner_poly]
all_n = [p[1] for p in inner_poly]
margin = 1.0  # 1 foot margin
e_min, e_max = min(all_e) - margin, max(all_e) + margin
n_min, n_max = min(all_n) - margin, max(all_n) + margin
W = int((e_max - e_min) * PPF)
H = int((n_max - n_min) * PPF)

def survey_to_px(e, n):
    """Convert survey (E, N) to pixel (x, y) where y=0 is top (north)."""
    x = (e - e_min) * PPF
    y = (n_max - n) * PPF  # flip: north = top
    return (x, y)

def poly_to_px(poly):
    """Convert survey polygon to pixel coordinates."""
    return [survey_to_px(p[0], p[1]) for p in poly]

# ---------------------------------------------------------------------------
# Build room polygons (same logic as _render_sf_extras)
# ---------------------------------------------------------------------------
_ro_list = compute_rough_openings(pts, layout)
_ro1_bd = [r for r in _ro_list if r.name == "RO1"][0].poly
_ro1_w_nf = _ro1_bd[3]    # RO1 NW: IW1 north face at RO1 west end
outer_openings = compute_outer_openings(pts, layout)
o6 = [o for o in outer_openings if o.name == "O6"][0]
o6_w = o6.poly[0]          # O6 SW: W-surface at O6 west edge
_w9w10_al, _ = seg_vecs(pts["W9"], pts["W10"])
_iw2s_e_al, _ = seg_vecs(layout.iw2s.poly[1], layout.iw2s.poly[2])
_iw2s_e_at_w9 = line_isect(layout.iw2s.poly[1], _iw2s_e_al,
                            pts["W9"], _w9w10_al)

# BEDROOM
bedroom_poly = [
    (layout.iw9.poly[2][0], layout.iw1.poly[0][1]),
    (layout.iw11.poly[3][0], layout.iw1.poly[0][1]),
    (layout.iw11.poly[3][0], pts["W1"][1]),
    (layout.iw9.poly[2][0], pts["W1"][1]),
]

# UTIL
util_poly = [
    layout.iw8.poly[0],
    layout.iw8.poly[1],
    layout.iw1.poly[0],
    layout.iw9.poly[3],
    (layout.iw9.poly[0][0], layout.iw7.poly[2][1]),
    layout.iw7.poly[3],
    layout.iw3.poly[3],
    (layout.iw3.poly[0][0], pts["W1"][1]),
    pts["W1"],
]
util_poly.extend(segment_polyline(data.inner_segs[0], pts)[1:])

# KITCHEN
kitchen_poly = [
    o6_w,
    _iw2s_e_at_w9,
    layout.iw2s.poly[1],
    layout.iw2o.poly[3],
    layout.iw2o.poly[0],
    layout.iw2.poly[2],
    layout.iw2.poly[1],
    _ro1_w_nf,
]

# LIVING
_iw1_ne = layout.iw1.poly[2]
_iw1_n_n = _iw1_ne[1]
_dash_at_iw1_n = _ro1_w_nf
living_poly = [o6_w]
living_poly.append(segment_polyline(data.inner_segs[6], pts)[-1])
for _si in range(7, 13):
    _pl = segment_polyline(data.inner_segs[_si], pts)
    living_poly.extend(_pl[1:])
living_poly.append(_iw1_ne)
living_poly.append(_dash_at_iw1_n)

# BATH
_seg2_pl = segment_polyline(data.inner_segs[2], pts)
_seg3_pl = segment_polyline(data.inner_segs[3], pts)
bath_poly = [
    layout.iw8.poly[3],
    layout.iw8.poly[2],
    layout.iw2.poly[3],
    layout.iw2o.poly[1],
    layout.iw2o.poly[2],
    layout.iw2s.poly[0],
    layout.iw2s.poly[3],
    _seg3_pl[0],
]
bath_poly.extend(reversed(_seg2_pl[:-1]))

# OFFICE
office_poly = [layout.iw5.poly[0]]
office_poly.append((pts["W15"][0], layout.iw5.poly[0][1]))
office_poly.append(pts["W15"])
for _si in [14, 15, 16]:
    office_poly.extend(segment_polyline(data.inner_segs[_si], pts)[1:])
office_poly.append(layout.iw4.poly[1])
office_poly.append(layout.iw4.poly[2])
office_poly.append(layout.iw12.poly[2])
office_poly.append(layout.iw12.poly[3])

# E CLOSET
e_closet_poly = [
    (layout.iw11.poly[1][0], layout.iw12.poly[0][1]),
    (layout.iw4.poly[0][0], layout.iw12.poly[0][1]),
    layout.iw4.poly[0],
    (layout.iw11.poly[1][0], pts["W1"][1]),
]

# W CLOSET
w_closet_poly = [
    (layout.iw3.poly[1][0], layout.iw7.poly[0][1]),
    (layout.iw9.poly[0][0], layout.iw7.poly[0][1]),
    (layout.iw9.poly[0][0], pts["W1"][1]),
    (layout.iw3.poly[1][0], pts["W1"][1]),
]

# STORAGE
storage_poly = [
    (layout.iw11.poly[1][0], layout.iw5.poly[3][1]),
    (pts["W14"][0], layout.iw5.poly[3][1]),
    (pts["W14"][0], layout.iw1.poly[0][1]),
    (layout.iw11.poly[1][0], layout.iw1.poly[0][1]),
]

# WH
_seg3_pl2 = segment_polyline(data.inner_segs[3], pts)
_seg4_pl = segment_polyline(data.inner_segs[4], pts)
_seg5_pl = segment_polyline(data.inner_segs[5], pts)
wh_poly = [layout.iw2s.poly[2]]
wh_poly.append(_seg3_pl2[-1])
wh_poly.extend(_seg4_pl[1:])
wh_poly.extend(_seg5_pl[1:])
wh_poly.append((layout.iw2s.poly[2][0], pts["W9"][1]))

# ---------------------------------------------------------------------------
# IW polygons
# ---------------------------------------------------------------------------
iw_names = ["iw1","iw8","iw2","iw2o","iw2s","iw3","iw7","iw9","iw6",
            "iw4","iw11","iw12","iw5"]
iw_polys = {name: getattr(layout, name).poly for name in iw_names}

# ---------------------------------------------------------------------------
# Room definitions: name, polygon, color (RGBA)
# ---------------------------------------------------------------------------
rooms = [
    ("BEDROOM",  bedroom_poly,  (255, 180, 180)),  # pink
    ("UTIL",     util_poly,     (180, 255, 180)),  # green
    ("KITCHEN",  kitchen_poly,  (180, 180, 255)),  # blue
    ("LIVING",   living_poly,   (255, 255, 150)),  # yellow
    ("BATH",     bath_poly,     (200, 150, 255)),  # purple
    ("OFFICE",   office_poly,   (255, 200, 150)),  # orange
    ("E CLOSET", e_closet_poly, (150, 220, 220)),  # teal
    ("W CLOSET", w_closet_poly, (220, 220, 150)),  # olive
    ("STORAGE",  storage_poly,  (220, 180, 220)),  # mauve
    ("WH",       wh_poly,       (180, 220, 200)),  # sage
]

# ---------------------------------------------------------------------------
# Create images
# ---------------------------------------------------------------------------
print(f"Raster: {W}x{H} pixels at {PPF} px/ft")
print(f"Survey bounds: E=[{e_min:.1f}, {e_max:.1f}], N=[{n_min:.1f}, {n_max:.1f}]")

# 1) W-series boundary mask (white = inside)
boundary_img = Image.new("L", (W, H), 0)
boundary_draw = ImageDraw.Draw(boundary_img)
boundary_draw.polygon(poly_to_px(inner_poly), fill=255)
boundary_mask = np.array(boundary_img)
boundary_pixels = np.sum(boundary_mask > 0)
print(f"W-boundary: {boundary_pixels} pixels = {boundary_pixels / PPF**2:.1f} sf")

# 2) Room accumulator: count how many rooms claim each pixel;
#    also track which room index owns each pixel (last writer wins for the ID)
room_count = np.zeros((H, W), dtype=np.int32)
room_id = np.full((H, W), -1, dtype=np.int32)   # -1 = unclaimed
room_img = Image.new("RGB", (W, H), (40, 40, 40))  # dark gray background
room_draw = ImageDraw.Draw(room_img)

# Draw W-boundary outline first (thin white line)
bnd_px = poly_to_px(inner_poly)
for i in range(len(bnd_px)):
    p1 = bnd_px[i]
    p2 = bnd_px[(i+1) % len(bnd_px)]
    room_draw.line([p1, p2], fill=(100, 100, 100), width=1)

print("\nRoom flood fills:")
room_masks = []   # keep masks for IW per-room analysis
for idx, (name, poly, color) in enumerate(rooms):
    # Create mask for this room
    room_mask = Image.new("L", (W, H), 0)
    rd = ImageDraw.Draw(room_mask)
    px_poly = poly_to_px(poly)
    rd.polygon(px_poly, fill=255)
    mask_arr = np.array(room_mask)
    room_masks.append(mask_arr)

    # Count pixels
    room_pixels = np.sum(mask_arr > 0)
    room_sf = room_pixels / PPF**2

    # Check: pixels outside W-boundary
    outside = np.sum((mask_arr > 0) & (boundary_mask == 0))
    outside_sf = outside / PPF**2

    # Check: pixels already claimed by another room
    overlap = np.sum((mask_arr > 0) & (room_count > 0))
    overlap_sf = overlap / PPF**2

    # Accumulate
    room_count += (mask_arr > 0).astype(np.int32)
    room_id[mask_arr > 0] = idx

    status = ""
    if outside > 0:
        status += f"  ** {outside_sf:.2f} sf OUTSIDE boundary **"
    if overlap > 0:
        status += f"  ** {overlap_sf:.2f} sf OVERLAPS prior room **"

    print(f"  {name:12s} {room_sf:7.2f} sf  ({room_pixels:7d} px){status}")

    # Draw filled polygon on composite image
    room_draw.polygon(px_poly, fill=color)

# 3) IW analysis — outlines only (room colors show through)
print("\nInterior wall fills:")
iw_count = np.zeros((H, W), dtype=np.int32)
for name in iw_names:
    poly = iw_polys[name]
    iw_mask = Image.new("L", (W, H), 0)
    iwd = ImageDraw.Draw(iw_mask)
    px_poly = poly_to_px(poly)
    iwd.polygon(px_poly, fill=255)
    mask_arr = np.array(iw_mask)

    iw_pixels = np.sum(mask_arr > 0)
    iw_sf = iw_pixels / PPF**2

    # Check overlap with each room
    in_room = np.sum((mask_arr > 0) & (room_count > 0))
    in_room_sf = in_room / PPF**2

    # Check outside boundary
    outside = np.sum((mask_arr > 0) & (boundary_mask == 0))
    outside_sf = outside / PPF**2

    iw_count += (mask_arr > 0).astype(np.int32)

    # Per-room breakdown
    room_overlaps = []
    for ri, (rname, _, _) in enumerate(rooms):
        ri_overlap = np.sum((mask_arr > 0) & (room_masks[ri] > 0))
        if ri_overlap > 0:
            room_overlaps.append((rname, ri_overlap / PPF**2))

    status = ""
    if in_room > 0:
        pct = 100.0 * in_room / max(iw_pixels, 1)
        room_detail = ", ".join(f"{rn}={rv:.2f}" for rn, rv in room_overlaps)
        status += f"  ** {in_room_sf:.2f} sf ({pct:.0f}%) in rooms: {room_detail} **"
    if outside > 0:
        status += f"  ** {outside_sf:.2f} sf outside boundary **"

    print(f"  {name:6s} {iw_sf:6.2f} sf  ({iw_pixels:6d} px){status}")

    # Draw IW as outline only — room color remains visible inside
    px_list = list(px_poly)
    px_list.append(px_list[0])  # close the polygon
    room_draw.line(px_list, fill=(50, 50, 50), width=2)

# 4) Gap analysis: pixels inside W-boundary not claimed by any room or IW
total_claimed = room_count + iw_count
gap_pixels = np.sum((boundary_mask > 0) & (total_claimed == 0))
gap_sf = gap_pixels / PPF**2

multi_room = np.sum(room_count > 1)
multi_room_sf = multi_room / PPF**2

print(f"\n{'='*60}")
print(f"SUMMARY")
print(f"{'='*60}")
print(f"  W-boundary interior:      {boundary_pixels/PPF**2:7.1f} sf")
print(f"  Room pixels (unique):     {np.sum(room_count > 0)/PPF**2:7.1f} sf")
print(f"  Room-room overlap:        {multi_room_sf:7.2f} sf ({np.sum(room_count > 1)} px)")
print(f"  IW material inside rooms: {np.sum((iw_count > 0) & (room_count > 0))/PPF**2:7.2f} sf")
print(f"  Gap (unclaimed):          {gap_sf:7.2f} sf ({gap_pixels} px)")
print(f"  Sum(room areas):          {sum(poly_area(p) for _,p,_ in rooms):7.2f} sf")

# Draw room-room overlaps in red
overlap_mask = (room_count > 1)
for y in range(H):
    for x in range(W):
        if overlap_mask[y, x]:
            room_img.putpixel((x, y), (255, 0, 0))

# Draw gaps in bright white
gap_mask = (boundary_mask > 0) & (total_claimed == 0)
for y in range(H):
    for x in range(W):
        if gap_mask[y, x]:
            room_img.putpixel((x, y), (255, 255, 255))

# Draw W-boundary outline (white)
for i in range(len(bnd_px)):
    p1 = bnd_px[i]
    p2 = bnd_px[(i+1) % len(bnd_px)]
    room_draw.line([p1, p2], fill=(255, 255, 255), width=2)

# Add legend
try:
    font = ImageFont.truetype("arial.ttf", 16)
    small_font = ImageFont.truetype("arial.ttf", 12)
except:
    font = ImageFont.load_default()
    small_font = font

legend_x = 10
legend_y = H - len(rooms) * 22 - 80
room_draw.rectangle([legend_x-5, legend_y-5, legend_x+200, H-10],
                     fill=(30, 30, 30), outline=(150, 150, 150))
for i, (name, _, color) in enumerate(rooms):
    y = legend_y + i * 22
    room_draw.rectangle([legend_x, y, legend_x+16, y+16], fill=color)
    room_draw.text((legend_x+22, y), name, fill=(255, 255, 255), font=small_font)

# Legend for special items
y = legend_y + len(rooms) * 22 + 4
room_draw.rectangle([legend_x, y, legend_x+16, y+16], fill=(80, 80, 80))
room_draw.text((legend_x+22, y), "Interior wall", fill=(255, 255, 255), font=small_font)
y += 22
room_draw.rectangle([legend_x, y, legend_x+16, y+16], fill=(255, 0, 0))
room_draw.text((legend_x+22, y), "Room overlap", fill=(255, 255, 255), font=small_font)
y += 22
room_draw.rectangle([legend_x, y, legend_x+16, y+16], fill=(255, 255, 255))
room_draw.text((legend_x+22, y), "Gap (unclaimed)", fill=(255, 255, 255), font=small_font)

out_path = os.path.join(os.path.dirname(__file__), "diag_floodfill.png")
room_img.save(out_path)
print(f"\nImage saved to {out_path}")
