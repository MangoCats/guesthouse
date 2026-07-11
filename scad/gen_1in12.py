"""Generate 1in12.scad — 1:12 sloped shed roof.

Identical to the 2:12 shed roof (scad/gen_2in12.py) but with a shallower 1:12
N-ward rising slope.  Delegates to gen_2in12.generate() with a slope override so
the two variants stay in lockstep.  All coordinates in feet.
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

from scad.gen_2in12 import generate as _generate_shed

_DIR = os.path.dirname(os.path.abspath(__file__))
_OUT = os.path.join(_DIR, "1in12.scad")

ROOF_SLOPE = 1.0 / 12.0  # 1:12 (for test introspection)


def generate(gd=None):
    _generate_shed(gd, slope_override=ROOF_SLOPE, out_path=_OUT)


if __name__ == "__main__":
    generate()
