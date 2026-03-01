"""Shared test fixtures for app-layer tests.

Provides isolated database, Flask test client, and geometry fixtures
that do not interfere with the project's production database.

The key challenge: app.engine.compute_geometry() monkey-patches
floorplan.constants then importlib.reload()s geometry, layout, and
openings.  Those modules have 'from floorplan.constants import ...'
at module scope, so reload gives them new copies of the patched
values.  Restoring floorplan.constants alone leaves the other three
modules stale.  We therefore snapshot and restore ALL four modules.
"""
import importlib
import os
import sys
import pytest

_PROJECT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _PROJECT not in sys.path:
    sys.path.insert(0, _PROJECT)

from app.database import init_db, get_constants_dict

_FLOORPLAN_MODULES = [
    "floorplan.constants",
    "floorplan.geometry",
    "floorplan.layout",
    "floorplan.openings",
]


def _snapshot_modules():
    """Snapshot every numeric attribute across all four floorplan modules."""
    snapshots = {}
    for mod_name in _FLOORPLAN_MODULES:
        mod = importlib.import_module(mod_name)
        snapshots[mod_name] = {
            a: getattr(mod, a) for a in dir(mod)
            if not a.startswith("__") and isinstance(getattr(mod, a), (int, float))
        }
    return snapshots


def _restore_modules(snapshots):
    """Restore every numeric attribute across all four floorplan modules."""
    for mod_name, values in snapshots.items():
        mod = importlib.import_module(mod_name)
        for name, value in values.items():
            setattr(mod, name, value)


@pytest.fixture
def fresh_db(tmp_path):
    """Create an isolated database seeded with default data."""
    snapshot = _snapshot_modules()
    db_path = str(tmp_path / "test.db")
    init_db(db_path)
    yield db_path
    _restore_modules(snapshot)


@pytest.fixture
def app_client(fresh_db):
    """Flask test client backed by a fresh isolated database."""
    from app.server import create_app
    app = create_app(db_path=fresh_db)
    app.config["TESTING"] = True
    with app.test_client() as client:
        yield client


@pytest.fixture
def geometry(fresh_db):
    """Pre-computed geometry from a fresh database."""
    from app.engine import compute_geometry
    constants = get_constants_dict(fresh_db)
    return compute_geometry(constants)
