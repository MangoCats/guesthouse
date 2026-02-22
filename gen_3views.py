"""Generate 3views.pdf — 2×2 landscape layout of elevation line drawings
and the bare floorplan.

Top-left:     patio elevation    (scad/2in12_patio_ld.png)
Top-right:    corner elevation   (scad/2in12_corner_ld.png)
Bottom-left:  bumpout elevation  (scad/2in12_bumpout_ld.png)
Bottom-right: bare floorplan     (floorplan/floorplan_bare.svg)
"""

import os

import fitz  # pymupdf

_DIR = os.path.dirname(os.path.abspath(__file__))

# US Letter landscape (points at 72 DPI)
_PAGE_W, _PAGE_H = 792, 612
_MARGIN = 24
_GUTTER = 12

# Input paths (relative to project root)
_IMAGES = [
    os.path.join(_DIR, "scad", "2in12_patio_ld.png"),
    os.path.join(_DIR, "scad", "2in12_corner_ld.png"),
    os.path.join(_DIR, "scad", "2in12_bumpout_ld.png"),
    os.path.join(_DIR, "floorplan", "floorplan_bare.svg"),
]

_OUT = os.path.join(_DIR, "3views.pdf")


def _cell_rect(col: int, row: int) -> fitz.Rect:
    """Return the bounding rect for a cell in the 2×2 grid."""
    cell_w = (_PAGE_W - 2 * _MARGIN - _GUTTER) / 2
    cell_h = (_PAGE_H - 2 * _MARGIN - _GUTTER) / 2
    x0 = _MARGIN + col * (cell_w + _GUTTER)
    y0 = _MARGIN + row * (cell_h + _GUTTER)
    return fitz.Rect(x0, y0, x0 + cell_w, y0 + cell_h)


def _fit_rect(img_w: float, img_h: float, cell: fitz.Rect) -> fitz.Rect:
    """Scale image to fit inside cell, preserving aspect ratio, centered."""
    scale = min(cell.width / img_w, cell.height / img_h)
    w, h = img_w * scale, img_h * scale
    cx, cy = (cell.x0 + cell.x1) / 2, (cell.y0 + cell.y1) / 2
    return fitz.Rect(cx - w / 2, cy - h / 2, cx + w / 2, cy + h / 2)


def _svg_to_pixmap(svg_path: str) -> fitz.Pixmap:
    """Render an SVG to a high-resolution pixmap via PyMuPDF."""
    doc = fitz.open(svg_path)
    page = doc[0]
    # Render at 3x for crisp output
    pix = page.get_pixmap(matrix=fitz.Matrix(3, 3))
    doc.close()
    return pix


def generate():
    doc = fitz.open()
    page = doc.new_page(width=_PAGE_W, height=_PAGE_H)

    positions = [(0, 0), (1, 0), (0, 1), (1, 1)]

    for path, (col, row) in zip(_IMAGES, positions):
        cell = _cell_rect(col, row)

        if path.endswith(".svg"):
            pix = _svg_to_pixmap(path)
            target = _fit_rect(pix.width, pix.height, cell)
            page.insert_image(target, pixmap=pix)
        else:
            # Get image dimensions
            img_doc = fitz.open(path)
            img_w, img_h = img_doc[0].rect.width, img_doc[0].rect.height
            img_doc.close()
            target = _fit_rect(img_w, img_h, cell)
            page.insert_image(target, filename=path)

    doc.save(_OUT)
    doc.close()
    print(f"Written to {_OUT}")


if __name__ == "__main__":
    generate()
