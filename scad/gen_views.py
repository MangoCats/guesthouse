"""Render OpenSCAD preview images for each .scad model."""

import os
import subprocess

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
        print(f"rendered {len(_VIEWS)} views for {scad}")


if __name__ == "__main__":
    generate()
