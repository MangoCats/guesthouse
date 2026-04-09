"""Generate black-on-white line drawing variants of SCAD view PNGs.

For each *_<view>.png in the scad/ directory, produces a *_<view>_ld.png
using Sobel edge detection.
"""

import os

import numpy as np
from PIL import Image
from scipy.ndimage import convolve, gaussian_filter, minimum_filter

_DIR = os.path.dirname(os.path.abspath(__file__))

_SOBEL_X = np.array([[-1, 0, 1], [-2, 0, 2], [-1, 0, 1]], dtype=np.float64)
_SOBEL_Y = np.array([[-1, -2, -1], [0, 0, 0], [1, 2, 1]], dtype=np.float64)

# View suffixes that gen_views.py produces
_VIEW_SUFFIXES = ("_patio", "_corner", "_bumpout")


def _to_line_drawing(src_path: str, dst_path: str) -> None:
    """Convert a rendered PNG to a black-on-white line drawing."""
    gray = np.array(Image.open(src_path).convert("L"), dtype=np.float64)
    gray = gaussian_filter(gray, sigma=0.8)
    gx = convolve(gray, _SOBEL_X)
    gy = convolve(gray, _SOBEL_Y)
    edges = np.sqrt(gx**2 + gy**2)
    edges = edges / edges.max() * 255
    # Binary threshold -> black lines on white
    result = 255 - (edges > 12).astype(np.uint8) * 255
    # Dilate lines slightly for a bolder pen stroke
    result = minimum_filter(result, size=2)
    Image.fromarray(result, mode="L").save(dst_path)


def generate():
    count = 0
    for fname in sorted(os.listdir(_DIR)):
        if not fname.endswith(".png"):
            continue
        base = fname.removesuffix(".png")
        if not any(base.endswith(s) for s in _VIEW_SUFFIXES):
            continue
        src = os.path.join(_DIR, fname)
        dst = os.path.join(_DIR, f"{base}_ld.png")
        _to_line_drawing(src, dst)
        count += 1
    print(f"generated {count} line drawings")


if __name__ == "__main__":
    generate()
