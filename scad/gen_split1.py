"""Generate split1.scad — 1:12 split (2-panel) roof.

Like the split2 roof (scad/gen_split2.py) but:

* the west panel (over office + kitchen) is strictly 1:12,
* the level east (bedroom-closet) low eave is pinned to 7'6",
* the two panels still meet along the same N-S seam as the split2 model, so the
  east panel's N-S slope matches the west (1:12) along the seam while its E-W
  tilt keeps the closet eave level — giving an east-panel slope of ~1.20:12
  (perpendicular to the back wall).

The west plane's height floats to whatever keeps the panels meeting at the seam
given the pinned closet eave.  Delegates to gen_split2.generate() so the two
variants stay in lockstep.  All coordinates in feet.
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

from scad.gen_split2 import generate as _generate_split

_DIR = os.path.dirname(os.path.abspath(__file__))
_OUT = os.path.join(_DIR, "split1.scad")

ROOF_SLOPE = 1.0 / 12.0            # west panel slope (for test introspection)
CLOSET_EAVE_ELEV = 7.5            # 7'6" pinned back-of-closet wall height


def generate(gd=None):
    _generate_split(gd, slope_override=ROOF_SLOPE,
                    low_elev_target=CLOSET_EAVE_ELEV, out_path=_OUT)


if __name__ == "__main__":
    generate()
