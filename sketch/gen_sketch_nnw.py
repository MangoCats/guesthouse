#!/usr/bin/env python3
"""Generate an architectural pen-and-ink sketch SVG of the ADU from NNW perspective.

Style: hand-drawn architectural illustration with siding texture,
window details, trees, and ground vegetation.

Output: sketch/sketch_nnw.svg
"""
import math
import os
import sys
import random

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

from floorplan.geometry import compute_outline_geometry
from floorplan.roof import compute_roof_geometry, roof_polyline

_DIR = os.path.dirname(os.path.abspath(__file__))

# ── Height constants (feet) ───────────────────────────────────────────
WALL_HT    = 80.0 / 12.0       # 6'8"
LOWER_HT   = 20.0 / 12.0       # 1'8" window sill
ROOF_THICK  = 18.0 / 12.0      # 1'6"
ROOF_REF    = 7.5               # 7'6" soffit at south edge
ROOF_SLOPE  = 2.0 / 12.0       # 2:12 northward
SOUTH_N     = -13.5

def rz_bot(n):
    return ROOF_REF + ROOF_SLOPE * (n - SOUTH_N)

def rz_top(n):
    return rz_bot(n) + ROOF_THICK

# ── Canvas ────────────────────────────────────────────────────────────
SVG_W, SVG_H = 1400, 1050
BG      = "#f5f0e0"
INK     = "#1a1208"
INK_L   = "#5a4d3a"
INK_M   = "#3a2d1a"
WALL_F  = "#e0d8c8"
ROOF_F  = "#d0c8b8"
FASC_F  = "#b8b0a0"
WIN_F   = "#4a5868"
DOOR_F  = "#2a2418"

random.seed(42)

# ── Vector helpers ────────────────────────────────────────────────────
def _norm(v):
    l = math.sqrt(sum(c*c for c in v))
    return tuple(c/l for c in v) if l > 1e-12 else v

def _cross(a, b):
    return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])

def _dot(a, b):
    return sum(x*y for x, y in zip(a, b))

# ── Camera (NNW, ~8' above grade) ────────────────────────────────────
AZ   = math.radians(332)
DIST = 55
LOOK = (-1.0, 0.0, 4.5)
_cam_e = DIST * math.sin(AZ)
_cam_n = DIST * math.cos(AZ)
_cam_z = 8.0                        # 8' above grade

_cam = (_cam_e, _cam_n, _cam_z)
_fwd = _norm((LOOK[0]-_cam[0], LOOK[1]-_cam[1], LOOK[2]-_cam[2]))
_rt  = _norm(_cross(_fwd, (0, 0, 1)))
_up  = _cross(_rt, _fwd)
FOC  = 1000

def proj(e, n, z):
    """Project world (E,N,Z) → screen (sx, sy)."""
    dx, dy, dz = e-_cam[0], n-_cam[1], z-_cam[2]
    x   = _dot((dx, dy, dz), _rt)
    y   = _dot((dx, dy, dz), _up)
    dep = _dot((dx, dy, dz), _fwd)
    if dep < 0.1:
        return None
    s = FOC / dep
    return (SVG_W/2 + x*s, SVG_H/2 - y*s)

# ── Geometry helpers ──────────────────────────────────────────────────
def arc_pts(center, r, start, end, cw=True, n=20):
    cx, cy = center
    a1 = math.atan2(start[1]-cy, start[0]-cx)
    a2 = math.atan2(end[1]-cy, end[0]-cx)
    if cw:
        if a2 > a1: a2 -= 2*math.pi
    else:
        if a2 < a1: a2 += 2*math.pi
    return [(cx+r*math.cos(a1+(a2-a1)*t/n),
             cy+r*math.sin(a1+(a2-a1)*t/n)) for t in range(n+1)]

def dedup(pts, tol=0.01):
    if not pts: return pts
    out = [pts[0]]
    for p in pts[1:]:
        if abs(p[0]-out[-1][0]) > tol or abs(p[1]-out[-1][1]) > tol:
            out.append(p)
    return out

# ── SVG path helpers ─────────────────────────────────────────────────
def sline(x1, y1, x2, y2, j=0.5):
    mx = (x1+x2)/2 + random.uniform(-j, j)
    my = (y1+y2)/2 + random.uniform(-j, j)
    return f"M{x1:.1f},{y1:.1f} Q{mx:.1f},{my:.1f} {x2:.1f},{y2:.1f}"

def polypath(pts):
    if not pts or any(p is None for p in pts): return ""
    d = f"M{pts[0][0]:.1f},{pts[0][1]:.1f}"
    for p in pts[1:]:
        d += f" L{p[0]:.1f},{p[1]:.1f}"
    return d + " Z"

def linepath(pts):
    if not pts or any(p is None for p in pts): return ""
    d = f"M{pts[0][0]:.1f},{pts[0][1]:.1f}"
    for p in pts[1:]:
        d += f" L{p[0]:.1f},{p[1]:.1f}"
    return d

# ── Tree / vegetation drawing ────────────────────────────────────────
def draw_tree(svg, bx, by, h, spread, bare=False):
    tx = bx + random.uniform(-2, 2)
    ty = by - h
    # Trunk
    for dx, sw in [(0, 2.5), (1.2, 1.2), (-0.8, 0.8)]:
        svg.append(f'<path d="{sline(bx+dx, by, tx+dx, ty+h*0.15, 2)}"'
                   f' fill="none" stroke="#3a2a15" stroke-width="{sw}"/>')
    # Branches
    tips = []
    for _ in range(random.randint(5, 8)):
        t = random.uniform(0.35, 0.95)
        sx = bx + t*(tx-bx); sy = by + t*(ty+h*0.15-by)
        angle = random.uniform(math.pi*0.15, math.pi*0.85)
        side = random.choice([-1, 1])
        bl = random.uniform(h*0.15, h*0.45)
        ex = sx + side*bl*math.cos(angle)
        ey = sy - bl*math.sin(angle)
        sw = random.uniform(0.8, 1.8)
        svg.append(f'<path d="{sline(sx, sy, ex, ey, 2.5)}"'
                   f' fill="none" stroke="#3a2a15" stroke-width="{sw}"/>')
        tips.append((ex, ey))
        n_sub = random.randint(2, 4) if bare else random.randint(1, 3)
        for _ in range(n_sub):
            sa = angle + random.uniform(-0.5, 0.5)
            sl = bl*random.uniform(0.2, 0.6)
            sex = ex + side*sl*math.cos(sa)
            sey = ey - sl*math.sin(sa)
            svg.append(f'<path d="{sline(ex, ey, sex, sey, 1.5)}"'
                       f' fill="none" stroke="#4a3a25" stroke-width="{random.uniform(0.5,1.0)}"/>')
            tips.append((sex, sey))
            # Tertiary twigs for bare trees
            if bare:
                for _ in range(random.randint(1, 2)):
                    ta = sa + random.uniform(-0.6, 0.6)
                    tl = sl * random.uniform(0.3, 0.5)
                    tex = sex + side*tl*math.cos(ta)
                    tey = sey - tl*math.sin(ta)
                    svg.append(f'<path d="{sline(sex, sey, tex, tey, 1)}"'
                               f' fill="none" stroke="#5a4a35" stroke-width="{random.uniform(0.3,0.6)}"/>')
    if bare:
        return
    # Foliage masses
    centers = [(tx, ty - h*0.1)] + [(t[0], t[1]-random.uniform(3, 10)) for t in tips]
    for fcx, fcy in centers:
        for _ in range(random.randint(20, 40)):
            r = random.uniform(3, spread*0.35)
            a = random.uniform(0, 2*math.pi)
            x1 = fcx + r*math.cos(a); y1 = fcy + r*math.sin(a)
            sl = random.uniform(3, 10)
            sa = random.uniform(0, 2*math.pi)
            x2 = x1+sl*math.cos(sa); y2 = y1+sl*math.sin(sa)
            op = random.uniform(0.3, 0.7)
            svg.append(f'<path d="{sline(x1, y1, x2, y2, 1)}"'
                       f' fill="none" stroke="#3a4a28" stroke-width="0.7" opacity="{op:.2f}"/>')
    # Canopy outline
    rx = spread*0.45; ry = h*0.35
    cpts = []
    for i in range(24):
        a = 2*math.pi*i/24
        rv = 1 + random.uniform(-0.15, 0.15)
        cpts.append((tx+rx*rv*math.cos(a), (ty-h*0.05)+ry*rv*math.sin(a)))
    cpts.append(cpts[0])
    d = f"M{cpts[0][0]:.1f},{cpts[0][1]:.1f}"
    for i in range(1, len(cpts)):
        mx = (cpts[i-1][0]+cpts[i][0])/2 + random.uniform(-2, 2)
        my = (cpts[i-1][1]+cpts[i][1])/2 + random.uniform(-2, 2)
        d += f" Q{mx:.1f},{my:.1f} {cpts[i][0]:.1f},{cpts[i][1]:.1f}"
    svg.append(f'<path d="{d}" fill="none" stroke="#3a4a28" stroke-width="1" opacity="0.5"/>')

def draw_ground_veg(svg, x1, x2, y, n_strokes=80):
    for _ in range(n_strokes):
        gx = random.uniform(x1, x2)
        gy = y + random.uniform(-3, 3)
        gh = random.uniform(4, 12)
        ga = random.uniform(math.pi*0.3, math.pi*0.7)
        gex = gx + gh*math.cos(ga)*random.choice([-1, 1])
        gey = gy - gh*math.sin(ga)
        op = random.uniform(0.3, 0.6)
        svg.append(f'<path d="{sline(gx, gy, gex, gey, 0.5)}"'
                   f' fill="none" stroke="#4a5a30" stroke-width="0.6" opacity="{op:.2f}"/>')

def draw_shrub(svg, cx, cy, w, h, density=30):
    for _ in range(density):
        r = random.uniform(2, w*0.4); a = random.uniform(0, 2*math.pi)
        x1 = cx+r*math.cos(a); y1 = cy+(h/w)*r*math.sin(a)
        sl = random.uniform(3, 8); sa = random.uniform(0, 2*math.pi)
        x2 = x1+sl*math.cos(sa); y2 = y1+sl*math.sin(sa)
        svg.append(f'<path d="{sline(x1, y1, x2, y2, 0.5)}"'
                   f' fill="none" stroke="#4a5838" stroke-width="0.6" opacity="0.5"/>')


# ══════════════════════════════════════════════════════════════════════
def main():
    geo = compute_outline_geometry()
    fp  = geo.fp_pts
    rg  = compute_roof_geometry(fp, geo.radii)
    rp  = rg.pts

    # ── Visible outline (F2 → ... → F12) ─────────────────────────────
    vis = []
    vis += [fp["F2"], fp["F5"]]
    vis += arc_pts(fp["C5"],  geo.radii["R_a5"],  fp["F5"],  fp["F6"],  cw=True,  n=15)
    vis += [fp["F6"], fp["F7"]]
    vis += arc_pts(fp["C7"],  geo.radii["R_a7"],  fp["F7"],  fp["F8"],  cw=True,  n=12)
    vis += arc_pts(fp["C8"],  geo.radii["R_a8"],  fp["F8"],  fp["F9"],  cw=False, n=8)
    vis += [fp["F9"], fp["F10"]]
    vis += arc_pts(fp["C10"], geo.radii["R_a10"], fp["F10"], fp["F11"], cw=False, n=10)
    vis += arc_pts(fp["C11a"],geo.radii["R_a11"], fp["F11"], fp["F11a"],cw=True,  n=10)
    vis += [fp["F11a"], fp["F11b"]]
    vis += arc_pts(fp["C11"], geo.radii["R_a11"], fp["F11b"],fp["F12"], cw=True,  n=10)
    vis = dedup(vis)

    # ── Roof north edge ──────────────────────────────────────────────
    ne = [rp["R1"], rp["R3s"]]
    ne += arc_pts(rg.r3_center, rg.r3_radius, rp["R3s"], rp["R3e"], cw=True, n=15)
    ne += [rp["R3e"], rp["R4s"]]
    ne += arc_pts(rg.r4_center, rg.r4_radius, rp["R4s"], rp["R4e"], cw=True, n=15)
    ne += [rp["R4e"]]
    ne = dedup(ne)

    rpoly = roof_polyline(rg)

    # ── Opening positions ────────────────────────────────────────────
    we = fp["F2"][0]                          # west wall E
    o3_ne = fp["F5"][1] - 8.0/12;            o3_ns = o3_ne - 32.0/12
    o2_ne = o3_ns - 48.0/12;                 o2_ns = o2_ne - 19.0/12
    o1_ne = o2_ns - 72.0/12;                 o1_ns = o1_ne - 19.0/12
    nn   = fp["F6"][1]                        # north wall N
    o4_ce = (fp["F6"][0]+fp["F7"][0])/2;      o4_hw = 4.5/12
    n910 = fp["F9"][1]                        # F9-F10 wall N
    o6_ee = fp["F10"][0]-6.0/12;              o6_es = o6_ee - 44.0/12
    o5_es = fp["F9"][0]+2.0;                  o5_ee = o5_es + 68.0/12

    # ══════════════════════════════════════════════════════════════════
    # Build SVG
    # ══════════════════════════════════════════════════════════════════
    svg = [f'<svg xmlns="http://www.w3.org/2000/svg" '
           f'width="{SVG_W}" height="{SVG_H}" viewBox="0 0 {SVG_W} {SVG_H}">']

    svg.append('''<defs>
  <filter id="sk" x="-2%" y="-2%" width="104%" height="104%">
    <feTurbulence type="turbulence" baseFrequency="0.025" numOctaves="3" seed="7" result="n"/>
    <feDisplacementMap in="SourceGraphic" in2="n" scale="1.0"
      xChannelSelector="R" yChannelSelector="G"/>
  </filter>
</defs>''')

    svg.append(f'<rect width="{SVG_W}" height="{SVG_H}" fill="{BG}"/>')
    svg.append('<g filter="url(#sk)">')

    # ── 1. Background trees ──────────────────────────────────────────
    # Trees placed well behind building so trunks don't cross the roof
    for (te, tn, th, ts) in [(25, 28, 160, 110), (-28, 30, 180, 130),
                              (35, 26, 130, 80), (-35, 25, 120, 75),
                              (0, 30, 150, 100)]:
        b = proj(te, tn, 0)
        if b:
            draw_tree(svg, b[0], b[1], th, ts)

    # ── 1b. Fence (behind building — drawn before wall fill) ────────
    # 4' two-rail fence parallel to F16-F17, 11' outward (exterior)
    f16e, f16n = fp["F16"]; f17e, f17n = fp["F17"]
    _fdx = f17e - f16e; _fdy = f17n - f16n
    _fln = math.sqrt(_fdx*_fdx + _fdy*_fdy)
    fence_dir = (_fdx/_fln, _fdy/_fln)      # parallel to F16-F17
    # Exterior normal (right of CW traversal = opposite of left normal)
    _ext_nx, _ext_ny = -_fdy/_fln, _fdx/_fln
    # Reference point: midpoint of F16-F17 offset 11' outward
    _fmid_e = (f16e + f17e) / 2 + 11.0 * _ext_nx
    _fmid_n = (f16n + f17n) / 2 + 11.0 * _ext_ny
    fence_ref = (_fmid_e, _fmid_n)
    FENCE_HT = 4.0
    RAIL1_Z = 1.5; RAIL2_Z = 3.2
    POST_SPACING = 8.0
    n_posts = 25                             # enough to span the page
    fence_pts = []
    for i in range(n_posts):
        t = (i - n_posts // 2) * POST_SPACING
        fe = fence_ref[0] + t * fence_dir[0]
        fn = fence_ref[1] + t * fence_dir[1]
        fence_pts.append((fe, fn))
    for fe, fn in fence_pts:
        pb = proj(fe, fn, 0)
        pt = proj(fe, fn, FENCE_HT)
        if pb and pt:
            svg.append(f'<path d="{sline(pb[0],pb[1],pt[0],pt[1],0.3)}"'
                       f' fill="none" stroke="{INK_L}" stroke-width="1.2" opacity="0.5"/>')
    for rz in [RAIL1_Z, RAIL2_Z]:
        rail_scr = []
        for fe, fn in fence_pts:
            p = proj(fe, fn, rz)
            if p:
                rail_scr.append(p)
        if len(rail_scr) >= 2:
            svg.append(f'<path d="{linepath(rail_scr)}"'
                       f' fill="none" stroke="{INK_L}" stroke-width="0.8" opacity="0.45"/>')

    # ── 1c. Background ground vegetation (behind building) ──────────
    bl = proj(-25, -14, 0)
    br = proj(20, -14, 0)
    if bl and br:
        draw_ground_veg(svg, bl[0]-80, br[0]+80, (bl[1]+br[1])/2+15, n_strokes=180)

    # ── 2. Building silhouette ─────────────────────────────────────────
    # Two overlapping opaque fills to ensure nothing behind bleeds through:
    #   A) Wall face: ground to roof-top along vis outline
    #   B) Full roof at soffit level: covers back/south overhang underside
    gnd  = [proj(e, n, 0)          for e, n in vis]
    rtop = [proj(e, n, rz_top(n))  for e, n in vis]
    gnd_ok  = [p for p in gnd if p]
    rtop_ok = [p for p in reversed(rtop) if p]
    wall_pts = gnd_ok + rtop_ok
    wp = polypath(wall_pts)
    if wp:
        svg.append(f'<path d="{wp}" fill="{WALL_F}" stroke="none"/>')
    # Full roof at soffit height — covers the back overhang underside
    rsoff = [proj(e, n, rz_bot(n)) for e, n in rpoly]
    rsoff_ok = [p for p in rsoff if p]
    rsp = polypath(rsoff_ok)
    if rsp:
        svg.append(f'<path d="{rsp}" fill="{WALL_F}" stroke="none"/>')

    # ── 3. Roof top fill ─────────────────────────────────────────────
    rt_proj = [proj(e, n, rz_top(n)) for e, n in rpoly]
    rtp = polypath([p for p in rt_proj if p])
    if rtp:
        svg.append(f'<path d="{rtp}" fill="{ROOF_F}" stroke="none"/>')

    # ── 4. Fascia strips (connect wall top to roof edge) ─────────────
    # The roof overhangs the wall by ROOF_OH. The fascia fills the gap
    # between the wall-face soffit line and the roof-edge soffit/top.
    # North fascia
    for i in range(len(ne)-1):
        e1, n1 = ne[i]; e2, n2 = ne[i+1]
        if abs(e1-e2) < 0.005 and abs(n1-n2) < 0.005:
            continue
        b1, b2 = proj(e1, n1, rz_bot(n1)), proj(e2, n2, rz_bot(n2))
        t1, t2 = proj(e1, n1, rz_top(n1)), proj(e2, n2, rz_top(n2))
        if all(x is not None for x in [b1, b2, t1, t2]):
            fd = polypath([b1, b2, t2, t1])
            if fd:
                svg.append(f'<path d="{fd}" fill="{FASC_F}" stroke="none"/>')

    # West fascia (R7→R1), sampled for smooth slope
    e1w, n1w = rp["R7"]; e2w, n2w = rp["R1"]
    for j in range(20):
        t0 = j/20; t1v = (j+1)/20
        se1 = e1w+t0*(e2w-e1w); sn1 = n1w+t0*(n2w-n1w)
        se2 = e1w+t1v*(e2w-e1w); sn2 = n1w+t1v*(n2w-n1w)
        b1, b2 = proj(se1, sn1, rz_bot(sn1)), proj(se2, sn2, rz_bot(sn2))
        t1, t2 = proj(se1, sn1, rz_top(sn1)), proj(se2, sn2, rz_top(sn2))
        if all(x is not None for x in [b1, b2, t1, t2]):
            fd = polypath([b1, b2, t2, t1])
            if fd:
                svg.append(f'<path d="{fd}" fill="{FASC_F}" stroke="none"/>')

    # Soffit underside strip (visible between wall face and roof edge)
    # Draw as a filled band connecting wall outline at soffit Z to roof edge at soffit Z
    for i in range(len(vis)-1):
        we1, wn1 = vis[i]; we2, wn2 = vis[i+1]
        if abs(we1-we2) < 0.005 and abs(wn1-wn2) < 0.005:
            continue
        # Find corresponding roof-edge points (offset outward by ROOF_OH)
        # Approximate: use the vis outline direction to compute outward offset
        dx, dy = we2-we1, wn2-wn1
        ln = math.sqrt(dx*dx+dy*dy)
        if ln < 0.01: continue
        # Left normal (outward for CW traversal)
        nx, ny = dy/ln, -dx/ln
        oh = 6.0/12.0
        re1, rn1 = we1+oh*nx, wn1+oh*ny
        re2, rn2 = we2+oh*nx, wn2+oh*ny
        wp1 = proj(we1, wn1, rz_bot(wn1)); wp2 = proj(we2, wn2, rz_bot(wn2))
        rp1 = proj(re1, rn1, rz_bot(rn1)); rp2 = proj(re2, rn2, rz_bot(rn2))
        if all(x is not None for x in [wp1, wp2, rp1, rp2]):
            fd = polypath([wp1, wp2, rp2, rp1])
            if fd:
                svg.append(f'<path d="{fd}" fill="{FASC_F}" stroke="none"/>')

    # ── 5. Adobe texture (subtle stipple/spatter) ────────────────────
    # Smooth adobe walls: sparse random short marks for surface grain
    for _ in range(400):
        idx = random.randint(0, len(vis)-2)
        e1, n1 = vis[idx]; e2, n2 = vis[idx+1]
        t = random.uniform(0, 1)
        pe, pn = e1+t*(e2-e1), n1+t*(n2-n1)
        z = random.uniform(0.3, min(WALL_HT, rz_bot(pn)-0.2))
        p = proj(pe, pn, z)
        if p:
            ml = random.uniform(1, 3)
            ma = random.uniform(0, 2*math.pi)
            x2 = p[0]+ml*math.cos(ma); y2 = p[1]+ml*math.sin(ma)
            op = random.uniform(0.06, 0.15)
            svg.append(f'<path d="M{p[0]:.1f},{p[1]:.1f} L{x2:.1f},{y2:.1f}"'
                       f' fill="none" stroke="{INK_L}" stroke-width="0.3" opacity="{op:.2f}"/>')

    # ── 6. Openings ──────────────────────────────────────────────────
    def draw_op(wall_val, is_west, s, e_val, zb, zt, door=False, pv=0, ph=0):
        if is_west:
            bl, br = proj(wall_val, s, zb), proj(wall_val, e_val, zb)
            tr, tl = proj(wall_val, e_val, zt), proj(wall_val, s, zt)
        else:
            bl, br = proj(s, wall_val, zb), proj(e_val, wall_val, zb)
            tr, tl = proj(e_val, wall_val, zt), proj(s, wall_val, zt)
        if not all(x is not None for x in [bl, br, tr, tl]):
            return
        fill = DOOR_F if door else WIN_F
        pp = polypath([bl, br, tr, tl])
        svg.append(f'<path d="{pp}" fill="{fill}" stroke="{INK_M}" stroke-width="1.2"/>')
        # Pane dividers
        for iv in range(1, max(pv, 1)):
            t = iv/pv
            p1 = (bl[0]+t*(br[0]-bl[0]), bl[1]+t*(br[1]-bl[1]))
            p2 = (tl[0]+t*(tr[0]-tl[0]), tl[1]+t*(tr[1]-tl[1]))
            svg.append(f'<path d="{sline(p1[0],p1[1],p2[0],p2[1],0.3)}"'
                       f' fill="none" stroke="#8a8070" stroke-width="0.8"/>')
        for ih in range(1, max(ph, 1)):
            t = ih/ph
            p1 = (bl[0]+t*(tl[0]-bl[0]), bl[1]+t*(tl[1]-bl[1]))
            p2 = (br[0]+t*(tr[0]-br[0]), br[1]+t*(tr[1]-br[1]))
            svg.append(f'<path d="{sline(p1[0],p1[1],p2[0],p2[1],0.3)}"'
                       f' fill="none" stroke="#8a8070" stroke-width="0.8"/>')

    # West wall: O1, O2 (windows), O3 (door)
    draw_op(we, True, o1_ns, o1_ne, LOWER_HT, WALL_HT)
    draw_op(we, True, o2_ns, o2_ne, LOWER_HT, WALL_HT)
    draw_op(we, True, o3_ns, o3_ne, 0, WALL_HT, door=True)
    # North wall: O4
    draw_op(nn, False, o4_ce-o4_hw, o4_ce+o4_hw, LOWER_HT, WALL_HT)
    # F9-F10 wall: O5 (window, 42" sill to 66" head), O6 (door)
    draw_op(n910, False, o5_es, o5_ee, 42.0/12.0, 66.0/12.0)
    draw_op(n910, False, o6_es, o6_ee, 0, WALL_HT, door=True)

    # ── 7. Heavy building edges ──────────────────────────────────────
    # Ground
    gp = [p for p in (proj(e, n, 0) for e, n in vis) if p]
    if gp:
        svg.append(f'<path d="{linepath(gp)}" fill="none" stroke="{INK}" stroke-width="2.5"/>')
    # Soffit line
    sp = [p for p in (proj(e, n, rz_bot(n)) for e, n in vis) if p]
    if sp:
        svg.append(f'<path d="{linepath(sp)}" fill="none" stroke="{INK}" stroke-width="1.8"/>')
    # Roof top north edge
    rtep = [p for p in (proj(e, n, rz_top(n)) for e, n in ne) if p]
    if rtep:
        svg.append(f'<path d="{linepath(rtep)}" fill="none" stroke="{INK}" stroke-width="2.2"/>')
    # Roof top west edge
    r7t = proj(rp["R7"][0], rp["R7"][1], rz_top(rp["R7"][1]))
    r1t = proj(rp["R1"][0], rp["R1"][1], rz_top(rp["R1"][1]))
    if r7t and r1t:
        svg.append(f'<path d="{sline(r7t[0],r7t[1],r1t[0],r1t[1],0.6)}"'
                   f' fill="none" stroke="{INK}" stroke-width="2.2"/>')
    # South/east roof edges (visible above building)
    for pa, pb in [(rp["R7"], rp["R6"]), (rp["R6"], rp["R5"]), (rp["R4e"], rp["R5"])]:
        p1 = proj(pa[0], pa[1], rz_top(pa[1]))
        p2 = proj(pb[0], pb[1], rz_top(pb[1]))
        if p1 and p2:
            svg.append(f'<path d="{sline(p1[0],p1[1],p2[0],p2[1],0.4)}"'
                       f' fill="none" stroke="{INK}" stroke-width="1.5" opacity="0.6"/>')
    # Vertical corners
    for name in ["F2", "F12"]:
        e, n = fp[name]
        b, t = proj(e, n, 0), proj(e, n, rz_bot(n))
        if b and t:
            svg.append(f'<path d="{sline(b[0],b[1],t[0],t[1],0.5)}"'
                       f' fill="none" stroke="{INK}" stroke-width="1.8"/>')
    # Roof fascia vertical corners
    for pt2 in [rp["R1"], rp["R7"]]:
        e, n = pt2
        b, t = proj(e, n, rz_bot(n)), proj(e, n, rz_top(n))
        if b and t:
            svg.append(f'<path d="{sline(b[0],b[1],t[0],t[1],0.3)}"'
                       f' fill="none" stroke="{INK}" stroke-width="1.5"/>')

    # ── Soffit shadow line (dark line just below roof overhang) ─────
    soff_pts = [proj(e, n, rz_bot(n) - 0.15) for e, n in vis]
    soff_ok2 = [p for p in soff_pts if p]
    if len(soff_ok2) >= 2:
        svg.append(f'<path d="{linepath(soff_ok2)}" fill="none" stroke="{INK}"'
                   f' stroke-width="1.2" opacity="0.3"/>')

    # ── 8. Roof texture (standing seam) ──────────────────────────────
    # Clip to roof top polygon so lines don't bleed onto wall face
    rt_clip = [p for p in rt_proj if p]
    _has_roofclip = False
    if rt_clip:
        cp = polypath(rt_clip)
        if cp:
            svg.append(f'<clipPath id="roofclip"><path d="{cp}"/></clipPath>')
            _has_roofclip = True
    if _has_roofclip:
        svg.append('<g clip-path="url(#roofclip)">')
    re_min = min(p[0] for p in rpoly)
    re_max = max(p[0] for p in rpoly)
    rn_min = min(p[1] for p in rpoly)
    rn_max = max(p[1] for p in rpoly)
    ev = re_min + 2.0
    while ev < re_max:
        p1 = proj(ev, rn_min, rz_top(rn_min))
        p2 = proj(ev, rn_max, rz_top(rn_max))
        if p1 and p2:
            svg.append(f'<path d="{sline(p1[0],p1[1],p2[0],p2[1],0.8)}"'
                       f' fill="none" stroke="{INK_L}" stroke-width="0.3" opacity="0.3"/>')
        ev += 2.0
    if _has_roofclip:
        svg.append('</g>')

    # ── 9. Ground and vegetation ─────────────────────────────────────
    # (Fence and background ground veg moved to before wall fill)
    # Foundation plantings — east of F11a and west of F7 only
    for se, sn in [(-16, 13.5), (-13, 13.5),               # west of F7
                    (11, 13.5), (14, 13), (17, 12)]:        # east of F11a
        sb = proj(se, sn, 0)
        if sb:
            draw_shrub(svg, sb[0], sb[1]-5, random.uniform(30, 50),
                       random.uniform(18, 30), density=50)
    # Foundation plantings along west wall (skip O3 zone: N ~ 8-11)
    for se, sn in [(-19, -10), (-19, -5), (-19, 0), (-19, 5)]:
        sb = proj(se, sn, 0)
        if sb:
            draw_shrub(svg, sb[0]-8, sb[1]-3, random.uniform(22, 35),
                       random.uniform(14, 22), density=35)

    # ── 10. Foreground trees ─────────────────────────────────────────
    # Large tree right side — moved further east to clear building
    fg = proj(36, -2, 0)
    if fg:
        draw_tree(svg, fg[0], fg[1], 350, 200, bare=False)
    # Bare/sparse tree left foreground
    fg2 = proj(-24, -10, 0)
    if fg2:
        draw_tree(svg, fg2[0], fg2[1], 300, 140, bare=True)

    # ── 11. Extra foreground vegetation ──────────────────────────────
    # Dense grass across full foreground (thicker strokes)
    for _ in range(200):
        gx = random.uniform(SVG_W*0.05, SVG_W*0.95)
        gy = random.uniform(SVG_H*0.75, SVG_H*0.96)
        gh = random.uniform(5, 20)
        svg.append(f'<path d="{sline(gx, gy, gx+random.uniform(-5,5), gy-gh, 0.5)}"'
                   f' fill="none" stroke="#4a5830" stroke-width="1.0" opacity="0.45"/>')
    # Extra dense patches at bottom corners
    for _ in range(100):
        gx = random.uniform(SVG_W*0.55, SVG_W*0.98)
        gy = random.uniform(SVG_H*0.85, SVG_H*0.98)
        gh = random.uniform(6, 22)
        svg.append(f'<path d="{sline(gx, gy, gx+random.uniform(-6,6), gy-gh, 0.5)}"'
                   f' fill="none" stroke="#4a5830" stroke-width="1.1" opacity="0.5"/>')
    for _ in range(100):
        gx = random.uniform(SVG_W*0.02, SVG_W*0.45)
        gy = random.uniform(SVG_H*0.85, SVG_H*0.98)
        gh = random.uniform(5, 18)
        svg.append(f'<path d="{sline(gx, gy, gx+random.uniform(-5,5), gy-gh, 0.5)}"'
                   f' fill="none" stroke="#4a5830" stroke-width="1.0" opacity="0.45"/>')

    # ── 12. Shadow hatching ──────────────────────────────────────────
    # Diagonal hatching on upper wall (above window height) for depth
    hatch_z = WALL_HT + 0.3
    while hatch_z < 13.0:
        pts_a = []
        pts_b = []
        for e, n in vis:
            if hatch_z < rz_bot(n) and hatch_z + 0.6 < rz_bot(n):
                pa = proj(e, n, hatch_z)
                pb = proj(e, n, hatch_z + 0.6)
                if pa and pb:
                    pts_a.append(pa)
                    pts_b.append(pb)
        for i in range(0, min(len(pts_a), len(pts_b)) - 1, 2):
            svg.append(f'<path d="{sline(pts_a[i][0],pts_a[i][1],pts_b[i+1][0],pts_b[i+1][1],0.5)}"'
                       f' fill="none" stroke="{INK_L}" stroke-width="0.3" opacity="0.2"/>')
        hatch_z += 0.8

    # Cross-hatching on north-facing wall sections for shade effect
    # North wall area: segments with N > 10 (the northern-facing walls)
    north_vis = [(e, n) for e, n in vis if n > 8.0]
    for z in [i * 0.5 for i in range(1, int(WALL_HT / 0.5))]:
        for k in range(0, len(north_vis) - 1, 3):
            e1, n1 = north_vis[k]
            e2_idx = min(k + 2, len(north_vis) - 1)
            e2, n2 = north_vis[e2_idx]
            pa = proj(e1, n1, z)
            pb = proj(e2, n2, z + 0.4)
            if pa and pb:
                svg.append(f'<path d="{sline(pa[0],pa[1],pb[0],pb[1],0.4)}"'
                           f' fill="none" stroke="{INK_L}" stroke-width="0.25" opacity="0.12"/>')

    svg.append('</g>')

    # Title
    svg.append(f'<text x="{SVG_W/2}" y="{SVG_H-20}" text-anchor="middle"'
               f' font-family="serif" font-size="14" fill="{INK_L}" opacity="0.6">'
               f'ADU \u2014 North-Northwest View</text>')

    svg.append('</svg>')

    out = os.path.join(_DIR, "sketch_nnw.svg")
    with open(out, "w") as f:
        f.write("\n".join(svg))
    print(f"Wrote {out}")

    # ── Render to PDF ──────────────────────────────────────────────────
    try:
        import fitz  # PyMuPDF
        svg_doc = fitz.open(out)                       # open SVG
        pix = svg_doc[0].get_pixmap(dpi=150)           # rasterize
        svg_doc.close()
        img_bytes = pix.tobytes("png")
        pdf_doc = fitz.open()
        page = pdf_doc.new_page(width=SVG_W, height=SVG_H)
        page.insert_image(fitz.Rect(0, 0, SVG_W, SVG_H), stream=img_bytes)
        pdf_out = os.path.join(_DIR, "sketch_nnw.pdf")
        pdf_doc.save(pdf_out)
        pdf_doc.close()
        print(f"Wrote {pdf_out}")
    except ImportError:
        print("PyMuPDF (fitz) not available — PDF output skipped")
    except Exception as ex:
        print(f"PDF rendering failed: {ex}")


if __name__ == "__main__":
    main()
