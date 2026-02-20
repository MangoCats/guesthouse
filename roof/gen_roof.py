"""Generate roof outline SVG showing R-series roof corners over building walls."""
import math
import os
from datetime import date

from shared.svg import W, H, git_describe
from shared.geometry import path_polygon, poly_area, fmt_dist
from floorplan.gen_floorplan import build_floorplan_data
from floorplan.roof import compute_roof_geometry, roof_polyline


def _svg_pts(pts_list, to_svg):
    """Convert a list of (E,N) points to SVG polygon points string."""
    return " ".join(f"{to_svg(*p)[0]:.2f},{to_svg(*p)[1]:.2f}" for p in pts_list)


def render_roof_svg(fp_data, roof):
    """Render roof SVG with building outline and roof outline overlay."""
    pts = fp_data.pts
    to_svg = fp_data.to_svg

    # Building outline polygon
    bldg_poly = path_polygon(fp_data.outline_segs, pts)

    # Roof outline polygon
    roof_poly = roof_polyline(roof)

    # --- Page layout ---
    all_pts = bldg_poly + roof_poly
    all_svg = [to_svg(*p) for p in all_pts]
    xmin = min(p[0] for p in all_svg) - 15
    xmax = max(p[0] for p in all_svg) + 15
    ymin = min(p[1] for p in all_svg) - 40
    ymax = max(p[1] for p in all_svg) + 15

    cx = (xmin + xmax) / 2
    cy = (ymin + ymax) / 2

    # Scale to fit US Letter landscape
    content_w = xmax - xmin
    content_h = ymax - ymin
    scale_x = W / content_w
    scale_y = H / content_h
    fit_scale = min(scale_x, scale_y) * 0.95  # 5% margin
    vb_w = W / fit_scale
    vb_h = H / fit_scale
    vb_x = cx - vb_w / 2
    vb_y = cy - vb_h / 2

    out = []
    out.append(f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}"'
               f' viewBox="{vb_x:.2f} {vb_y:.2f} {vb_w:.2f} {vb_h:.2f}">')
    out.append(f'<rect x="{vb_x:.2f}" y="{vb_y:.2f}"'
               f' width="{vb_w:.2f}" height="{vb_h:.2f}" fill="white"/>')

    # Title
    title_x = cx
    title_y = ymin + 8
    out.append(f'<text x="{title_x:.1f}" y="{title_y:.1f}"'
               f' text-anchor="middle" font-family="Arial" font-size="14"'
               f' font-weight="bold">Roof Plan</text>')

    # Building outline (solid gray)
    bldg_svg = _svg_pts(bldg_poly, to_svg)
    out.append(f'<polygon points="{bldg_svg}" fill="rgba(220,220,220,0.3)"'
               f' stroke="#999" stroke-width="0.5"/>')

    # Roof outline (dotted)
    roof_svg = _svg_pts(roof_poly, to_svg)
    out.append(f'<polygon points="{roof_svg}" fill="none"'
               f' stroke="#333" stroke-width="0.8" stroke-dasharray="3,2"/>')

    # R-series point labels
    for name in ["R1", "R2", "R5", "R6", "R7"]:
        sx, sy = to_svg(*roof.pts[name])
        out.append(f'<circle cx="{sx:.1f}" cy="{sy:.1f}" r="1.5"'
                   f' fill="#333" stroke="none"/>')
        out.append(f'<text x="{sx + 3:.1f}" y="{sy + 1:.1f}"'
                   f' font-family="Arial" font-size="5" fill="#333">{name}</text>')

    # R3, R4 labels at the midpoint of each fillet
    for name, start_name, end_name in [("R3", "R3s", "R3e"), ("R4", "R4s", "R4e")]:
        ps = roof.pts[start_name]
        pe = roof.pts[end_name]
        mid_e = (ps[0] + pe[0]) / 2
        mid_n = (ps[1] + pe[1]) / 2
        sx, sy = to_svg(mid_e, mid_n)
        out.append(f'<text x="{sx:.1f}" y="{sy - 3:.1f}"'
                   f' text-anchor="middle" font-family="Arial" font-size="5"'
                   f' fill="#333">{name}</text>')

    # Title block
    tb_w, tb_h = 100, 45
    tb_x = vb_x + vb_w - tb_w - 5
    tb_y = vb_y + vb_h - tb_h - 5
    out.append(f'<rect x="{tb_x:.1f}" y="{tb_y:.1f}"'
               f' width="{tb_w}" height="{tb_h}"'
               f' fill="white" stroke="#999" stroke-width="0.5"/>')

    _y = tb_y + 10
    out.append(f'<text x="{tb_x + tb_w / 2:.1f}" y="{_y:.1f}"'
               f' text-anchor="middle" font-family="Arial" font-size="7"'
               f' font-weight="bold">Roof Plan</text>')
    _y += 10
    out.append(f'<text x="{tb_x + tb_w / 2:.1f}" y="{_y:.1f}"'
               f' text-anchor="middle" font-family="Arial" font-size="6"'
               f' fill="#333">Area: {roof.area:.1f} sq ft</text>')
    _y += 8
    out.append(f'<text x="{tb_x + tb_w / 2:.1f}" y="{_y:.1f}"'
               f' text-anchor="middle" font-family="Arial" font-size="5"'
               f' fill="#666">{date.today().isoformat()}  {git_describe()}</text>')

    out.append('</svg>')
    return "\n".join(out)


if __name__ == "__main__":
    fp_data = build_floorplan_data()
    roof = compute_roof_geometry(fp_data.pts, fp_data.radii)

    svg_content = render_roof_svg(fp_data, roof)
    _dir = os.path.dirname(os.path.abspath(__file__))
    svg_path = os.path.join(_dir, "roof.svg")
    with open(svg_path, "w", encoding="utf-8") as f:
        f.write(svg_content)
    print(f"Wrote {svg_path}")
    print(f"Roof area: {roof.area:.2f} sq ft")
