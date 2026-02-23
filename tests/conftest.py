"""Shared test fixtures for hut2 geometry tests."""
import importlib.util
import io
import os
import pytest
from unittest.mock import patch
from shared.survey import compute_traverse, compute_three_arc, compute_inset, rotate_pts, COORD_ROTATION
from shared.geometry import compute_inner_walls, path_polygon, left_norm
from floorplan.geometry import compute_outline_geometry, OutlineAnchors
from floorplan.layout import compute_interior_layout
from floorplan.constants import WALL_OUTER, WALL_EXTRA

_PROJECT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _import_from(subdir, module_name):
    """Import a module from a non-package subdirectory."""
    path = os.path.join(_PROJECT, subdir, f"{module_name}.py")
    spec = importlib.util.spec_from_file_location(f"{subdir}.{module_name}", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def dist_sq(a, b):
    """Squared distance between two 2D points."""
    return (a[0] - b[0])**2 + (a[1] - b[1])**2


def generate_scad_content(generate_fn):
    """Capture SCAD generator output to string."""
    buf = io.StringIO()
    with patch("builtins.open", return_value=buf):
        with patch.object(buf, "close"):
            generate_fn()
    return buf.getvalue()


@pytest.fixture(scope="session")
def traverse():
    """Base traverse: (pts, p3_trav)."""
    return compute_traverse()


@pytest.fixture(scope="session")
def pts_base(traverse):
    """pts dict with P3/POB/P2/P4/P5."""
    return traverse[0]


@pytest.fixture(scope="session")
def arc_info(pts_base):
    """Three-arc result dict; also mutates pts_base with T*/TC* points."""
    return compute_three_arc(pts_base)


@pytest.fixture(scope="session")
def inset_result(pts_base, arc_info):
    """InsetResult from compute_inset (pure, does not mutate pts_base)."""
    return compute_inset(
        pts_base, arc_info["R1"], arc_info["R2"], arc_info["R3"],
        arc_info["nE"], arc_info["nN"],
    )


@pytest.fixture(scope="session")
def pts_full(pts_base, inset_result):
    """pts dict with all points including inset, rotated by COORD_ROTATION."""
    pts = dict(pts_base)
    pts.update(inset_result.pts_update)
    rotate_pts(pts, COORD_ROTATION)
    return pts


@pytest.fixture(scope="session")
def outline_geo(pts_full, inset_result):
    """OutlineGeometry from compute_outline_geometry."""
    # Shift anchors outward for 10" wall (2" beyond original 8")
    _ln = left_norm(pts_full["PiX"], pts_full["Pi5"])
    anchors = OutlineAnchors(
        Pi2=(pts_full["Pi2"][0] - WALL_EXTRA, pts_full["Pi2"][1]),
        Pi3=(pts_full["Pi3"][0] - WALL_EXTRA, pts_full["Pi3"][1] - WALL_EXTRA),
        Ti3=pts_full["Ti3"],
        PiX=(pts_full["PiX"][0] - WALL_EXTRA * _ln[0], pts_full["PiX"][1] - WALL_EXTRA * _ln[1]),
        Pi5=(pts_full["Pi5"][0] - WALL_EXTRA * _ln[0], pts_full["Pi5"][1] - WALL_EXTRA * _ln[1]),
        TC1=pts_full["TC1"], R1i=inset_result.R1i,
    )
    return compute_outline_geometry(anchors)


@pytest.fixture(scope="session")
def pts_with_outline(pts_full, outline_geo):
    """pts dict with F-series and W-series."""
    pts = dict(pts_full)
    pts.update(outline_geo.fp_pts)
    inner_segs = compute_inner_walls(
        outline_geo.outline_segs, pts, WALL_OUTER, outline_geo.radii,
    )
    return pts, inner_segs


@pytest.fixture(scope="session")
def inner_poly(pts_with_outline):
    """Inner wall polygon."""
    pts, inner_segs = pts_with_outline
    return path_polygon(inner_segs, pts)


@pytest.fixture(scope="session")
def layout(pts_with_outline, inner_poly):
    """InteriorLayout result."""
    pts, _ = pts_with_outline
    return compute_interior_layout(pts, inner_poly)


@pytest.fixture(scope="session")
def span_geometry():
    """Full geometry tuple from span _build_geometry()."""
    mod = _import_from("span", "gen_span")
    return mod._build_geometry()
