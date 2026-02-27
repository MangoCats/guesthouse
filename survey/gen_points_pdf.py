"""Generate points.pdf: coordinate table and key measurements."""

import datetime
import math
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

import fitz  # pymupdf
from shared.survey import compute_traverse
from shared.geometry import fmt_dist
from shared.svg import git_describe
from floorplan.geometry import compute_outline_geometry, align_pts_to_f_series
from survey.gen_path_svg_ks import compute_k_points

_DIR = os.path.dirname(os.path.abspath(__file__))

# Point names in display order
F_NAMES = [
    "F1", "F2", "F5", "F6", "F7", "F8", "F9", "F10",
    "F11", "F11a", "F11b", "F12", "F13", "F14", "F15",
    "F16", "F17", "F18",
]
P_NAMES = ["POB", "P2", "P3", "P4", "P5"]
K_NAMES = ["K1", "K2", "K3", "K4", "K5", "K6"]

# Layout constants (points)
PAGE_W, PAGE_H = 612, 792  # US Letter
MARGIN = 54
FONT = "helv"
TITLE_SIZE = 12
HEADER_SIZE = 9
CELL_SIZE = 9
ROW_H = 16
HEADER_ROW_H = 18


def _point_to_line_dist(pt, line_a, line_b):
    """Perpendicular distance from pt to the line through line_a and line_b."""
    dx = line_b[0] - line_a[0]
    dy = line_b[1] - line_a[1]
    length = math.sqrt(dx * dx + dy * dy)
    return abs(dy * pt[0] - dx * pt[1] + line_b[0] * line_a[1] - line_b[1] * line_a[0]) / length


def _draw_coord_table(page, title, sections, y_start):
    """Draw a titled three-column table with section separators.

    sections: list of (section_label, [(name, (E, N)), ...])
    """
    col_widths = [60, 120, 120]
    table_w = sum(col_widths)
    x0 = (PAGE_W - table_w) / 2

    # Title
    tw = fitz.get_text_length(title, fontname=FONT, fontsize=TITLE_SIZE)
    page.insert_text(fitz.Point((PAGE_W - tw) / 2, y_start + TITLE_SIZE),
                     title, fontname=FONT, fontsize=TITLE_SIZE)
    y = y_start + TITLE_SIZE + 8

    # Count total rows (header + section headers + data)
    total_data_rows = sum(len(rows) for _, rows in sections)
    n_separators = len(sections)  # thin separator line before each section
    table_h = HEADER_ROW_H + total_data_rows * ROW_H

    shape = page.new_shape()

    # Outer border
    shape.draw_rect(fitz.Rect(x0, y, x0 + table_w, y + table_h))
    # Header bottom line
    shape.draw_line(fitz.Point(x0, y + HEADER_ROW_H),
                    fitz.Point(x0 + table_w, y + HEADER_ROW_H))
    # Vertical column lines
    cx = x0
    for j in range(len(col_widths) - 1):
        cx += col_widths[j]
        shape.draw_line(fitz.Point(cx, y), fitz.Point(cx, y + table_h))

    shape.finish(color=(0, 0, 0), width=0.5)

    # Section separator lines (lighter weight)
    row_y = y + HEADER_ROW_H
    for sec_idx, (_, rows) in enumerate(sections):
        if sec_idx > 0:
            shape.draw_line(fitz.Point(x0, row_y), fitz.Point(x0 + table_w, row_y))
            shape.finish(color=(0.5, 0.5, 0.5), width=0.3)
        row_y += len(rows) * ROW_H

    shape.commit()

    # Header text
    headers = ["Point", "Easting (ft)", "Northing (ft)"]
    for j, hdr in enumerate(headers):
        cx = x0 + sum(col_widths[:j])
        tw = fitz.get_text_length(hdr, fontname="hebo", fontsize=HEADER_SIZE)
        tx = cx + (col_widths[j] - tw) / 2
        ty = y + HEADER_ROW_H - 5
        page.insert_text(fitz.Point(tx, ty), hdr, fontname="hebo",
                         fontsize=HEADER_SIZE)

    # Data rows
    row_y = y + HEADER_ROW_H
    for _, rows in sections:
        for name, (e, n) in rows:
            cells = [name, f"{e:+.4f}", f"{n:+.4f}"]
            for j, val in enumerate(cells):
                cx = x0 + sum(col_widths[:j])
                tw = fitz.get_text_length(val, fontname=FONT, fontsize=CELL_SIZE)
                tx = cx + (col_widths[j] - tw) / 2
                ty = row_y + ROW_H - 4
                page.insert_text(fitz.Point(tx, ty), val, fontname=FONT,
                                 fontsize=CELL_SIZE)
            row_y += ROW_H

    return y + table_h


def _draw_measurements_table(page, title, rows, y_start):
    """Draw a two-column key measurements table.

    rows: list of (description, value_string)
    """
    col_widths = [200, 100]
    table_w = sum(col_widths)
    x0 = (PAGE_W - table_w) / 2

    # Title
    tw = fitz.get_text_length(title, fontname=FONT, fontsize=TITLE_SIZE)
    page.insert_text(fitz.Point((PAGE_W - tw) / 2, y_start + TITLE_SIZE),
                     title, fontname=FONT, fontsize=TITLE_SIZE)
    y = y_start + TITLE_SIZE + 8

    table_h = HEADER_ROW_H + len(rows) * ROW_H
    shape = page.new_shape()

    # Outer border
    shape.draw_rect(fitz.Rect(x0, y, x0 + table_w, y + table_h))
    # Header bottom
    shape.draw_line(fitz.Point(x0, y + HEADER_ROW_H),
                    fitz.Point(x0 + table_w, y + HEADER_ROW_H))
    # Vertical separator
    shape.draw_line(fitz.Point(x0 + col_widths[0], y),
                    fitz.Point(x0 + col_widths[0], y + table_h))
    shape.finish(color=(0, 0, 0), width=0.5)
    shape.commit()

    # Headers
    headers = ["Measurement", "Distance"]
    for j, hdr in enumerate(headers):
        cx = x0 + sum(col_widths[:j])
        tw = fitz.get_text_length(hdr, fontname="hebo", fontsize=HEADER_SIZE)
        tx = cx + (col_widths[j] - tw) / 2
        ty = y + HEADER_ROW_H - 5
        page.insert_text(fitz.Point(tx, ty), hdr, fontname="hebo",
                         fontsize=HEADER_SIZE)

    # Data rows
    row_y = y + HEADER_ROW_H
    PAD = 6
    for desc, val in rows:
        # Description left-aligned with padding
        ty = row_y + ROW_H - 4
        page.insert_text(fitz.Point(x0 + PAD, ty), desc,
                         fontname=FONT, fontsize=CELL_SIZE)
        # Value centered
        tw = fitz.get_text_length(val, fontname=FONT, fontsize=CELL_SIZE)
        tx = x0 + col_widths[0] + (col_widths[1] - tw) / 2
        page.insert_text(fitz.Point(tx, ty), val,
                         fontname=FONT, fontsize=CELL_SIZE)
        row_y += ROW_H

    return y + table_h


def main():
    # Compute all points
    pts = compute_traverse()
    align_pts_to_f_series(pts)
    outline = compute_outline_geometry()
    pts.update(outline.fp_pts)
    k_pts = compute_k_points(pts)
    pts.update(k_pts)

    # Build coordinate table sections
    sections = [
        ("F-series", [(n, pts[n]) for n in F_NAMES]),
        ("P-series", [(n, pts[n]) for n in P_NAMES]),
        ("K-series", [(n, pts[n]) for n in K_NAMES]),
    ]

    # Compute key measurements
    p4, p5 = pts["P4"], pts["P5"]
    k5_dist = _point_to_line_dist(pts["K5"], p4, p5)
    k6_dist = _point_to_line_dist(pts["K6"], p4, p5)

    measurements = [
        ("K5 to P4-P5 line", fmt_dist(k5_dist)),
        ("K6 to P4-P5 line", fmt_dist(k6_dist)),
    ]

    doc = fitz.open()
    page = doc.new_page(width=PAGE_W, height=PAGE_H)

    y = MARGIN
    y = _draw_coord_table(page, "Point Coordinates (FC frame)", sections, y)
    y += 24
    _draw_measurements_table(page, "Key Measurements", measurements, y)

    # Footer: Generated timestamp + repo info
    cap_fs = 7
    cap_text = f"Generated {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}  {git_describe()}"
    cap_tw = fitz.get_text_length(cap_text, fontname=FONT, fontsize=cap_fs)
    page.insert_text(fitz.Point((PAGE_W - cap_tw) / 2, PAGE_H - MARGIN / 2),
                     cap_text, fontname=FONT, fontsize=cap_fs,
                     color=(0.4, 0.4, 0.4))

    out_path = os.path.join(_DIR, "points.pdf")
    doc.save(out_path)
    doc.close()
    print(f"Written to {out_path}")


if __name__ == "__main__":
    main()
