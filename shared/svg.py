"""SVG transform factory and page constants."""
from typing import Callable
from .types import Point

# US Letter landscape at 72 dpi (11" x 8.5")
W, H = 792, 612

# SVG calibration: derived from known P3 survey coordinates on a reference PDF.
# P3 maps to x=368.79pt, POB maps to x=151.26pt; survey distance P3–POB = 18.66'.
_CALIB_X_P3 = 368.79
_CALIB_X_POB = 151.26
_CALIB_DIST = 18.66
_s = (_CALIB_X_P3 - _CALIB_X_POB) / _CALIB_DIST  # SVG points per survey foot

_CALIB_Y_P3 = 124.12  # P3 y-position in SVG points

def make_svg_transform(p3_trav: Point) -> Callable[[float, float], tuple[float, float]]:
    """Create to_svg closure from P3 traverse position."""
    px = _CALIB_X_P3 + p3_trav[0] * _s
    py = _CALIB_Y_P3 - p3_trav[1] * _s
    def to_svg(e: float, n: float) -> tuple[float, float]:
        return (px + e * _s, py - n * _s)
    return to_svg
