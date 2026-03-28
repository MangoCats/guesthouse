"""Generate roof outline SVG showing roof corners over building walls."""
import math
import os
import sys
from datetime import date

# Ensure project root is on sys.path for package imports
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

from shared.svg import W, H, git_describe


def _svg_pts(pts_list, to_svg):
    """Convert a list of (E,N) points to SVG polygon points string."""
    return " ".join(f"{to_svg(*p)[0]:.2f},{to_svg(*p)[1]:.2f}" for p in pts_list)


def render_roof_svg(fp_data, roof=None):
    """Render roof SVG with building outline and roof outline overlay.

    fp_data may be a FloorplanData or a GeneratorData.  If roof is None,
    it is taken from fp_data.roof (GeneratorData path).

    Handles both legacy RoofGeometry (R-series labels) and DbRoofResult
    (corner-center labels).
    """
    from shared.svg import make_svg_transform
    from shared.roof_outline import DbRoofResult

    # Support both GeneratorData and FloorplanData
    pts = fp_data.pts
    to_svg = getattr(fp_data, 'to_svg', None) or make_svg_transform()
    bldg_poly = getattr(fp_data, 'outline_poly', None)
    if bldg_poly is None:
        from shared.geometry import path_polygon
        bldg_poly = path_polygon(fp_data.outline_segs, pts)

    if roof is None:
        roof = fp_data.roof
    roof_poly = getattr(fp_data, 'roof_poly', None)
    if roof_poly is None and roof is not None:
        if isinstance(roof, DbRoofResult):
            roof_poly = roof.poly
        else:
            from floorplan.roof import roof_polyline
            roof_poly = roof_polyline(roof)

    is_db_roof = isinstance(roof, DbRoofResult)

    # --- Page layout ---
    all_pts = list(bldg_poly)
    if roof_poly:
        all_pts += list(roof_poly)
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
    if roof_poly:
        roof_svg = _svg_pts(roof_poly, to_svg)
        out.append(f'<polygon points="{roof_svg}" fill="none"'
                   f' stroke="#333" stroke-width="0.8" stroke-dasharray="3,2"/>')

    # Corner labels
    if is_db_roof and roof is not None:
        for name, pt in roof.corner_pts.items():
            sx, sy = to_svg(*pt)
            out.append(f'<circle cx="{sx:.1f}" cy="{sy:.1f}" r="1.5"'
                       f' fill="#333" stroke="none"/>')
            out.append(f'<text x="{sx + 3:.1f}" y="{sy + 1:.1f}"'
                       f' font-family="Arial" font-size="5" fill="#333">{name}</text>')
    elif roof is not None and not is_db_roof:
        # Legacy R-series labels
        for name in ["R1", "R5", "R6", "R7"]:
            if name not in roof.pts:
                continue
            sx, sy = to_svg(*roof.pts[name])
            out.append(f'<circle cx="{sx:.1f}" cy="{sy:.1f}" r="1.5"'
                       f' fill="#333" stroke="none"/>')
            out.append(f'<text x="{sx + 3:.1f}" y="{sy + 1:.1f}"'
                       f' font-family="Arial" font-size="5" fill="#333">{name}</text>')
        for name, start_name, end_name in [("R3", "R3s", "R3e"), ("R4", "R4s", "R4e")]:
            if start_name not in roof.pts or end_name not in roof.pts:
                continue
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
    area = roof.area if roof is not None else 0.0
    out.append(f'<text x="{tb_x + tb_w / 2:.1f}" y="{_y:.1f}"'
               f' text-anchor="middle" font-family="Arial" font-size="6"'
               f' fill="#333">Area: {area:.1f} sq ft</text>')
    _y += 8
    out.append(f'<text x="{tb_x + tb_w / 2:.1f}" y="{_y:.1f}"'
               f' text-anchor="middle" font-family="Arial" font-size="5"'
               f' fill="#666">{date.today().isoformat()}  {git_describe()}</text>')

    out.append('</svg>')
    return "\n".join(out)


if __name__ == "__main__":
    from floorplan.gen_floorplan import build_floorplan_data

    fp_data = build_floorplan_data()
    from floorplan.roof import compute_roof_geometry
    roof = compute_roof_geometry(fp_data.pts, fp_data.radii)

    svg_content = render_roof_svg(fp_data, roof)
    _dir = os.path.dirname(os.path.abspath(__file__))
    svg_path = os.path.join(_dir, "roof.svg")
    with open(svg_path, "w", encoding="utf-8") as f:
        f.write(svg_content)
    print(f"Wrote {svg_path}")
    print(f"Roof area: {roof.area:.2f} sq ft")
