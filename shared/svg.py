"""SVG transform factory and page constants."""
import os, subprocess
from typing import Callable
from .types import Point
from .survey import FC_IN_P3

# Cache file written by gen_all.py so all SVGs embed the same git describe.
_GIT_DESCRIBE_CACHE = os.path.join(os.path.dirname(__file__), os.pardir, ".git_describe")


def git_describe() -> str:
    """Return git describe string, preferring a cached value from gen_all.py."""
    try:
        with open(_GIT_DESCRIBE_CACHE) as f:
            return f.read().strip()
    except FileNotFoundError:
        return subprocess.check_output(
            ["git", "describe", "--always", "--dirty=-DEV"], text=True
        ).strip()

# US Letter landscape at 72 dpi (11" x 8.5")
W, H = 792, 612

# SVG calibration: derived from known P3 survey coordinates on a reference PDF.
# P3 maps to x=368.79pt, POB maps to x=151.26pt; survey distance P3–POB = 18.66'.
_CALIB_X_P3 = 368.79
_CALIB_X_POB = 151.26
_CALIB_DIST = 18.66
_s = (_CALIB_X_P3 - _CALIB_X_POB) / _CALIB_DIST  # SVG points per survey foot

_CALIB_Y_P3 = 124.12  # P3 y-position in SVG points

def svg_polygon_pts(points, to_svg, prec=1) -> str:
    """Format polygon points string from (E,N) coords via to_svg transform."""
    fmt = f".{prec}f"
    return " ".join(f"{to_svg(*p)[0]:{fmt}},{to_svg(*p)[1]:{fmt}}" for p in points)

def normalize_svg_angle(deg: float) -> float:
    """Normalize angle to [-90, 90] range for readable SVG text rotation."""
    if deg > 90:
        deg -= 180
    elif deg < -90:
        deg += 180
    return deg

def make_svg_transform(p3_trav: Point) -> Callable[[float, float], tuple[float, float]]:
    """Create to_svg closure from P3 traverse position.

    Coordinates are FC-based (building center = origin).
    The FC offset compensates so to_svg(0,0) maps FC to the correct SVG pixel.
    """
    px = _CALIB_X_P3 + (p3_trav[0] + FC_IN_P3[0]) * _s
    py = _CALIB_Y_P3 - (p3_trav[1] + FC_IN_P3[1]) * _s
    def to_svg(e: float, n: float) -> tuple[float, float]:
        return (px + e * _s, py - n * _s)
    return to_svg
