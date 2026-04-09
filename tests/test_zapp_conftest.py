"""Shared test fixtures for app-layer tests.

Provides isolated database, Flask test client, and geometry fixtures
that do not interfere with the project's production database.

The key challenge: app.engine.compute_geometry() monkey-patches
floorplan.constants then importlib.reload()s geometry, layout, and
openings.  Those modules have 'from floorplan.constants import ...'
at module scope, so reload gives them new copies of the patched
values.  Restoring floorplan.constants alone leaves the other three
modules stale.  We therefore snapshot and restore ALL four modules.

Performance note: fresh_db uses a session-scoped template DB that is
seeded once and copied via shutil.copy2 (~0.003s) rather than calling
init_db() (~0.307s) per test.
"""
import importlib
import os
import shutil
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
def fresh_db(tmp_path, _db_template):
    """Create an isolated database by copying the session template (~0.003s vs ~0.307s for init_db)."""
    snapshot = _snapshot_modules()
    db_path = str(tmp_path / "test.db")
    shutil.copy2(_db_template, db_path)
    yield db_path
    _restore_modules(snapshot)


@pytest.fixture
def app_client(fresh_db):
    """Flask test client backed by a fresh isolated database.

    Passes skip_init=True since fresh_db is already fully seeded.
    """
    from app.server import create_app
    app = create_app(db_path=fresh_db, skip_init=True)
    app.config["TESTING"] = True
    with app.test_client() as client:
        yield client


@pytest.fixture
def geometry(fresh_db):
    """Pre-computed geometry from a fresh database."""
    from app.engine import compute_geometry
    constants = get_constants_dict(fresh_db)
    return compute_geometry(constants, db_path=fresh_db)
