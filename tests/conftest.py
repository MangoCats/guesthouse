"""Shared test fixtures for hut2 geometry tests."""
import pytest
from shared.survey import compute_traverse, compute_three_arc, compute_inset
from shared.geometry import compute_inner_walls, path_polygon
from floorplan.geometry import compute_outline_geometry, OutlineAnchors
from floorplan.layout import compute_interior_layout
from floorplan.constants import WALL_OUTER


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
    """pts dict with all points including inset and F-series."""
    pts = dict(pts_base)
    pts.update(inset_result.pts_update)
    return pts


@pytest.fixture(scope="session")
def outline_geo(pts_full, inset_result):
    """OutlineGeometry from compute_outline_geometry."""
    anchors = OutlineAnchors(
        Pi2=pts_full["Pi2"], Pi3=pts_full["Pi3"], Ti3=pts_full["Ti3"],
        PiX=pts_full["PiX"], Pi5=pts_full["Pi5"],
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
