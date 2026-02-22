"""Render OpenSCAD preview images for each .scad model.

Each render is post-processed into a black-on-white line drawing using
Sobel edge detection.
"""

import os
import subprocess

import numpy as np
from PIL import Image
from scipy.ndimage import convolve, gaussian_filter, minimum_filter

_DIR = os.path.dirname(os.path.abspath(__file__))

# OpenSCAD executable
_OPENSCAD = r"C:\Program Files\OpenSCAD\openscad.com"

# Camera views: (name_suffix, eye_x, eye_y, eye_z, center_x, center_y, center_z)
_VIEWS = [
    ("patio",   12,    95,   5,   18.5, 13, 3.5),  # north face from patio, 5'
    ("corner", -33,   -68,   5,   18.5, 13, 3.5),  # SW property corner, 5'
    ("bumpout", 101.9, 50.6, 12,  18.5, 13, 6),    # NE residence bump-out, 12'
]

_IMG_SIZE = "1200,900"

_SOBEL_X = np.array([[-1, 0, 1], [-2, 0, 2], [-1, 0, 1]], dtype=np.float64)
_SOBEL_Y = np.array([[-1, -2, -1], [0, 0, 0], [1, 2, 1]], dtype=np.float64)


def _to_line_drawing(png_path: str) -> None:
    """Convert a rendered PNG in-place to a black-on-white line drawing."""
    gray = np.array(Image.open(png_path).convert("L"), dtype=np.float64)
    gray = gaussian_filter(gray, sigma=0.8)
    gx = convolve(gray, _SOBEL_X)
    gy = convolve(gray, _SOBEL_Y)
    edges = np.sqrt(gx**2 + gy**2)
    edges = edges / edges.max() * 255
    # Binary threshold → black lines on white
    result = 255 - (edges > 12).astype(np.uint8) * 255
    # Dilate lines slightly for a bolder pen stroke
    result = minimum_filter(result, size=2)
    Image.fromarray(result, mode="L").save(png_path)


def generate():
    if not os.path.isfile(_OPENSCAD):
        print(f"OpenSCAD not found at {_OPENSCAD}, skipping views")
        return

    scad_files = sorted(f for f in os.listdir(_DIR) if f.endswith(".scad"))
    for scad in scad_files:
        scad_path = os.path.join(_DIR, scad)
        base = scad.removesuffix(".scad")
        for name, ex, ey, ez, cx, cy, cz in _VIEWS:
            out_png = os.path.join(_DIR, f"{base}_{name}.png")
            cam = f"{ex},{ey},{ez},{cx},{cy},{cz}"
            subprocess.run(
                [_OPENSCAD, "--preview", f"--imgsize={_IMG_SIZE}",
                 "--projection=p", f"--camera={cam}",
                 "-o", out_png, scad_path],
                cwd=_DIR, capture_output=True, text=True,
            )
            _to_line_drawing(out_png)
        print(f"rendered {len(_VIEWS)} views for {scad}")


if __name__ == "__main__":
    generate()
