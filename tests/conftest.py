"""Shared test fixtures for ADU geometry tests."""
import importlib.util
import io
import os
import pytest
from unittest.mock import patch
from shared.survey import compute_traverse, compute_three_arc, compute_inset
from shared.geometry import compute_inner_walls, path_polygon
from floorplan.geometry import compute_outline_geometry, align_pts_to_f_series
from floorplan.layout import compute_interior_layout
from floorplan.constants import WALL_OUTER

from app.database import init_db as _init_db


@pytest.fixture(scope="session")
def _db_template(tmp_path_factory):
    """Seed a database once per session; fresh_db copies it per test.

    Session-scoped so init_db (~0.3s) runs only once regardless of how many
    test files import fresh_db from test_zapp_conftest.

    Also runs _seed_reference_plumbing_if_needed so the template is fully
    equivalent to a DB that has been through create_app() startup.
    """
    db_path = str(tmp_path_factory.mktemp("db_template") / "template.db")
    _init_db(db_path)
    return db_path

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
def pts_base():
    """pts dict with P3/POB/P2/P4/P5 (FC-based)."""
    return compute_traverse()


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
    """pts dict with all points including inset, aligned to F-series."""
    pts = dict(pts_base)
    pts.update(inset_result.pts_update)
    align_pts_to_f_series(pts)
    return pts


@pytest.fixture(scope="session")
def outline_geo():
    """OutlineGeometry from compute_outline_geometry."""
    return compute_outline_geometry()


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
    """Full geometry tuple from span._common.build_geometry()."""
    from span._common import build_geometry
    return build_geometry()
