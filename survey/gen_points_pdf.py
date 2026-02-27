"""Generate points.pdf: F-series and P-series coordinate tables."""

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

import fitz  # pymupdf
from shared.survey import compute_traverse
from floorplan.geometry import compute_outline_geometry, align_pts_to_f_series

_DIR = os.path.dirname(os.path.abspath(__file__))

# F-series point names in outline traversal order (no F3/F4)
F_NAMES = [
    "F1", "F2", "F5", "F6", "F7", "F8", "F9", "F10",
    "F11", "F11a", "F11b", "F12", "F13", "F14", "F15",
    "F16", "F17", "F18",
]

# P-series point names
P_NAMES = ["POB", "P2", "P3", "P4", "P5"]

# Layout constants (points)
PAGE_W, PAGE_H = 612, 792  # US Letter
MARGIN = 72
FONT = "helv"
TITLE_SIZE = 14
HEADER_SIZE = 10
CELL_SIZE = 10
COL_WIDTHS = [60, 120, 120]  # Point, Easting, Northing
ROW_H = 18
HEADER_ROW_H = 20


def _draw_table(page, title, names, pts, y_start):
    """Draw a titled three-column table. Returns y after the table."""
    table_w = sum(COL_WIDTHS)
    x0 = (PAGE_W - table_w) / 2

    # Title
    tw = fitz.get_text_length(title, fontname=FONT, fontsize=TITLE_SIZE)
    page.insert_text(fitz.Point((PAGE_W - tw) / 2, y_start + TITLE_SIZE),
                     title, fontname=FONT, fontsize=TITLE_SIZE)
    y = y_start + TITLE_SIZE + 10

    headers = ["Point", "Easting (ft)", "Northing (ft)"]
    n_rows = len(names) + 1  # +1 for header
    shape = page.new_shape()

    # Horizontal lines
    for i in range(n_rows + 1):
        h = HEADER_ROW_H if i == 0 else ROW_H
        ly = y + HEADER_ROW_H + (i - 1) * ROW_H if i > 0 else y
        if i == n_rows:
            ly = y + HEADER_ROW_H + (n_rows - 1) * ROW_H
        shape.draw_line(fitz.Point(x0, ly), fitz.Point(x0 + table_w, ly))

    # Vertical lines
    table_h = HEADER_ROW_H + len(names) * ROW_H
    cx = x0
    for j in range(len(COL_WIDTHS) + 1):
        shape.draw_line(fitz.Point(cx, y), fitz.Point(cx, y + table_h))
        if j < len(COL_WIDTHS):
            cx += COL_WIDTHS[j]

    shape.finish(color=(0, 0, 0), width=0.5)
    shape.commit()

    # Header text (bold via "hebo")
    for j, hdr in enumerate(headers):
        cx = x0 + sum(COL_WIDTHS[:j])
        tw = fitz.get_text_length(hdr, fontname="hebo", fontsize=HEADER_SIZE)
        tx = cx + (COL_WIDTHS[j] - tw) / 2
        ty = y + HEADER_ROW_H - 6
        page.insert_text(fitz.Point(tx, ty), hdr, fontname="hebo",
                         fontsize=HEADER_SIZE)

    # Data rows
    for i, name in enumerate(names):
        ry = y + HEADER_ROW_H + i * ROW_H
        e, n = pts[name]
        cells = [name, f"{e:+.4f}", f"{n:+.4f}"]
        for j, val in enumerate(cells):
            cx = x0 + sum(COL_WIDTHS[:j])
            tw = fitz.get_text_length(val, fontname=FONT, fontsize=CELL_SIZE)
            tx = cx + (COL_WIDTHS[j] - tw) / 2
            ty = ry + ROW_H - 5
            page.insert_text(fitz.Point(tx, ty), val, fontname=FONT,
                             fontsize=CELL_SIZE)

    return y + table_h


def main():
    # Compute points
    pts = compute_traverse()
    align_pts_to_f_series(pts)
    outline = compute_outline_geometry()
    pts.update(outline.fp_pts)

    doc = fitz.open()
    page = doc.new_page(width=PAGE_W, height=PAGE_H)

    y = MARGIN
    y = _draw_table(page, "F-Series Points (FC frame)", F_NAMES, pts, y)
    y += 30
    _draw_table(page, "P-Series Points (FC frame)", P_NAMES, pts, y)

    out_path = os.path.join(_DIR, "points.pdf")
    doc.save(out_path)
    doc.close()
    print(f"Written to {out_path}")


if __name__ == "__main__":
    main()
