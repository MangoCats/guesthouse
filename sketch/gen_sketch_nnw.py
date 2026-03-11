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

# ── Camera (NNW, slightly elevated) ──────────────────────────────────
AZ   = math.radians(332)
EL   = math.radians(12)
DIST = 55
LOOK = (-1.0, 0.0, 4.5)

_cam = (DIST*math.sin(AZ), DIST*math.cos(AZ), LOOK[2]+DIST*math.sin(EL))
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

    # ── Visible outline (F2 → ... → F10) ─────────────────────────────
    # Skip F1-F2 arc (SW corner faces away) and F10+ (NE bump is edge-on)
    vis = []
    vis += [fp["F2"], fp["F5"]]
    vis += arc_pts(fp["C5"],  geo.radii["R_a5"],  fp["F5"],  fp["F6"],  cw=True,  n=15)
    vis += [fp["F6"], fp["F7"]]
    vis += arc_pts(fp["C7"],  geo.radii["R_a7"],  fp["F7"],  fp["F8"],  cw=True,  n=12)
    vis += arc_pts(fp["C8"],  geo.radii["R_a8"],  fp["F8"],  fp["F9"],  cw=False, n=8)
    vis += [fp["F9"], fp["F10"]]
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

    # ── 2. Wall fill ─────────────────────────────────────────────────
    gnd  = [proj(e, n, 0)          for e, n in vis]
    soff = [proj(e, n, rz_bot(n))  for e, n in vis]
    gnd_ok  = [p for p in gnd if p]
    soff_ok = [p for p in reversed(soff) if p]
    wall_pts = gnd_ok + soff_ok
    wp = polypath(wall_pts)
    if wp:
        svg.append(f'<path d="{wp}" fill="{WALL_F}" stroke="none"/>')
    # Mask below ground line to clean up any projection artifacts
    if len(gnd_ok) >= 2:
        base_y = max(p[1] for p in gnd_ok) + 3
        max_y = base_y + 60
        mask_pts = [(gnd_ok[0][0] - 40, base_y),
                    (gnd_ok[-1][0] + 40, base_y),
                    (gnd_ok[-1][0] + 40, max_y),
                    (gnd_ok[0][0] - 40, max_y)]
        mfp = polypath(mask_pts)
        if mfp:
            svg.append(f'<path d="{mfp}" fill="{BG}" stroke="none"/>')

    # ── 3. Roof top fill ─────────────────────────────────────────────
    rt_proj = [proj(e, n, rz_top(n)) for e, n in rpoly]
    rtp = polypath([p for p in rt_proj if p])
    if rtp:
        svg.append(f'<path d="{rtp}" fill="{ROOF_F}" stroke="none"/>')

    # ── 4. Fascia strips ─────────────────────────────────────────────
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

    # West fascia (R7→R1), sampled
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

    # ── 5. Siding lines ──────────────────────────────────────────────
    z = 4.0/12.0
    while z < 14.0:
        pts_z = []
        for e, n in vis:
            if z < rz_bot(n):
                p2 = proj(e, n, z)
                if p2:
                    pts_z.append(p2)
        if len(pts_z) >= 2:
            d = linepath(pts_z)
            if d:
                op = 0.35 if z < WALL_HT else 0.25
                svg.append(f'<path d="{d}" fill="none" stroke="{INK_L}"'
                           f' stroke-width="0.4" opacity="{op}"/>')
        z += 6.0/12.0

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
    draw_op(we, True, o1_ns, o1_ne, LOWER_HT, WALL_HT, pv=2, ph=2)
    draw_op(we, True, o2_ns, o2_ne, LOWER_HT, WALL_HT, pv=2, ph=2)
    draw_op(we, True, o3_ns, o3_ne, 0, WALL_HT, door=True, pv=2, ph=3)
    # North wall: O4
    draw_op(nn, False, o4_ce-o4_hw, o4_ce+o4_hw, LOWER_HT, WALL_HT)
    # F9-F10 wall: O5 (window), O6 (door)
    draw_op(n910, False, o5_es, o5_ee, LOWER_HT, WALL_HT, pv=4, ph=2)
    draw_op(n910, False, o6_es, o6_ee, 0, WALL_HT, door=True, pv=2, ph=3)

    # ── 7. Heavy building edges ──────────────────────────────────────
    # Ground
    gp = [p for p in (proj(e, n, 0) for e, n in vis) if p]
    if gp:
        svg.append(f'<path d="{linepath(gp)}" fill="none" stroke="{INK}" stroke-width="2.5"/>')
    # Soffit line
    sp = [p for p in (proj(e, n, rz_bot(n)) for e, n in vis) if p]
    if sp:
        svg.append(f'<path d="{linepath(sp)}" fill="none" stroke="{INK}" stroke-width="1.8"/>')
    # Wall-height line
    wtp = [p for p in (proj(e, n, WALL_HT) for e, n in vis) if p]
    if wtp:
        svg.append(f'<path d="{linepath(wtp)}" fill="none" stroke="{INK_L}"'
                   f' stroke-width="0.6" opacity="0.5"/>')
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
    for name in ["F2", "F10"]:
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

    # ── 9. Ground and vegetation ─────────────────────────────────────
    # Extended ground line
    for offset in [0, 0.5, -0.5]:
        gl1 = proj(-35, -14+offset, 0)
        gl2 = proj(30, -14+offset, 0)
        if gl1 and gl2:
            sw = 1.8 if offset == 0 else 0.6
            svg.append(f'<path d="{sline(gl1[0],gl1[1],gl2[0],gl2[1],1.5)}"'
                       f' fill="none" stroke="{INK}" stroke-width="{sw}" opacity="0.5"/>')
    bl = proj(-25, -14, 0)
    br = proj(20, -14, 0)
    if bl and br:
        draw_ground_veg(svg, bl[0]-80, br[0]+80, (bl[1]+br[1])/2+15, n_strokes=180)
    # Foundation plantings along north wall
    for se, sn in [(-16, 13.5), (-12, 13.5), (-8, 13), (-4, 12),
                    (0, 12), (3, 12), (6, 12)]:
        sb = proj(se, sn, 0)
        if sb:
            draw_shrub(svg, sb[0], sb[1]-5, random.uniform(30, 50),
                       random.uniform(18, 30), density=50)
    # Foundation plantings along west wall
    for se, sn in [(-19, -10), (-19, -5), (-19, 0), (-19, 5), (-19, 9)]:
        sb = proj(se, sn, 0)
        if sb:
            draw_shrub(svg, sb[0]-8, sb[1]-3, random.uniform(22, 35),
                       random.uniform(14, 22), density=35)

    # ── 10. Foreground trees ─────────────────────────────────────────
    # Large tree right side (partially overlapping building)
    fg = proj(20, -6, 0)
    if fg:
        draw_tree(svg, fg[0], fg[1], 350, 200, bare=False)
    # Bare/sparse tree left foreground
    fg2 = proj(-24, -10, 0)
    if fg2:
        draw_tree(svg, fg2[0], fg2[1], 300, 140, bare=True)

    # ── 11. Extra foreground vegetation ──────────────────────────────
    for _ in range(60):
        gx = random.uniform(SVG_W*0.55, SVG_W*0.98)
        gy = random.uniform(SVG_H*0.82, SVG_H*0.95)
        gh = random.uniform(5, 18)
        svg.append(f'<path d="{sline(gx, gy, gx+random.uniform(-5,5), gy-gh, 0.5)}"'
                   f' fill="none" stroke="#4a5830" stroke-width="0.6" opacity="0.4"/>')
    for _ in range(40):
        gx = random.uniform(SVG_W*0.02, SVG_W*0.35)
        gy = random.uniform(SVG_H*0.82, SVG_H*0.95)
        gh = random.uniform(5, 15)
        svg.append(f'<path d="{sline(gx, gy, gx+random.uniform(-5,5), gy-gh, 0.5)}"'
                   f' fill="none" stroke="#4a5830" stroke-width="0.5" opacity="0.35"/>')

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


if __name__ == "__main__":
    main()
