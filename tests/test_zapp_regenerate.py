"""Tests for Phase 19: DB-driven in-process regeneration.

Verifies that the app's regeneration endpoints produce SVGs that reflect
database state (element deletions, constant changes) rather than hardcoded
procedural module output.
"""
import json
import os
import pytest

from app.database import init_db, get_constants_dict, get_outline_chain
from app.engine import (
    build_generator_data_from_db, generate_svg_db,
    _run_generator_inprocess,
)

_PROJECT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


# ── fixtures ──────────────────────────────────────────────────────────

@pytest.fixture
def fresh_db(tmp_path):
    db_path = str(tmp_path / "test.db")
    init_db(db_path)
    return db_path


@pytest.fixture
def gd(fresh_db):
    return build_generator_data_from_db(fresh_db)


# ── build_generator_data_from_db ──────────────────────────────────────

class TestBuildGeneratorData:
    def test_returns_generator_data(self, gd):
        assert hasattr(gd, "pts")
        assert hasattr(gd, "outline_segs")
        assert hasattr(gd, "outline_poly")
        assert hasattr(gd, "inner_poly")
        assert hasattr(gd, "roof_poly")
        assert hasattr(gd, "layout")

    def test_pts_contain_f_series(self, gd):
        assert "F1" in gd.pts
        assert "F2" in gd.pts

    def test_outline_poly_nonempty(self, gd):
        assert len(gd.outline_poly) > 10


# ── _run_generator_inprocess dispatch ─────────────────────────────────

class TestInprocessDispatch:
    def test_floorplan_generates_files(self, gd, fresh_db):
        # Ensure output dir exists
        fp_dir = os.path.join(_PROJECT, "floorplan")
        assert os.path.isdir(fp_dir)
        result = _run_generator_inprocess("floorplan/gen_floorplan.py", gd,
                                          db_path=fresh_db)
        assert result is True
        assert os.path.exists(os.path.join(fp_dir, "floorplan.svg"))

    def test_roof_generates_file(self, gd):
        result = _run_generator_inprocess("roof/gen_roof.py", gd)
        assert result is True
        assert os.path.exists(os.path.join(_PROJECT, "roof", "roof.svg"))

    def test_walls_generates_files(self, gd):
        result = _run_generator_inprocess("walls/gen_walls.py", gd)
        assert result is True
        assert os.path.exists(os.path.join(_PROJECT, "walls", "walls.svg"))
        assert os.path.exists(os.path.join(_PROJECT, "walls", "all_walls.svg"))

    def test_span_generates_file(self, gd):
        result = _run_generator_inprocess("span/gen_span.py", gd)
        assert result is True
        assert os.path.exists(os.path.join(_PROJECT, "span", "span.svg"))

    def test_span_minmax_generates_file(self, gd):
        result = _run_generator_inprocess("span/gen_span_minmax.py", gd)
        assert result is True
        assert os.path.exists(os.path.join(_PROJECT, "span", "span_minmax.svg"))

    def test_span_min_generates_file(self, gd):
        result = _run_generator_inprocess("span/gen_span_min.py", gd)
        assert result is True
        assert os.path.exists(os.path.join(_PROJECT, "span", "span_min.svg"))

    def test_plumbing_generates_file(self, gd):
        result = _run_generator_inprocess("plumbing/gen_plumbing.py", gd)
        assert result is True
        assert os.path.exists(os.path.join(_PROJECT, "plumbing", "plumbing.svg"))

    def test_unknown_script_returns_none(self, gd):
        result = _run_generator_inprocess("nonexistent/script.py", gd)
        assert result is None

    def test_subprocess_fallback_scripts_return_none(self, gd):
        assert _run_generator_inprocess("scad/gen_views.py", gd) is None
        assert _run_generator_inprocess("scad/gen_line_drawings.py", gd) is None
        assert _run_generator_inprocess("gen_3views.py", gd) is None
        assert _run_generator_inprocess("survey/gen_path_svg.py", gd) is None


# ── generate_svg_db fallback ──────────────────────────────────────────

class TestGenerateSvgDb:
    def test_known_script_succeeds(self, gd):
        ok = generate_svg_db("roof", "roof/gen_roof.py", gd)
        assert ok is True

    def test_unknown_script_falls_back_to_subprocess(self, gd):
        # survey/gen_path_svg.py has no in-process handler but exists
        ok = generate_svg_db("path_area", "survey/gen_path_svg.py", gd)
        assert ok is True


# ── API integration ───────────────────────────────────────────────────

@pytest.fixture
def app_client(fresh_db):
    from app.server import create_app
    app = create_app(db_path=fresh_db)
    app.config["TESTING"] = True
    with app.test_client() as client:
        yield client


class TestRegenerateAPI:
    def test_regenerate_single_view_uses_db(self, app_client):
        resp = app_client.post(
            "/api/regenerate",
            data=json.dumps({"view": "floorplan"}),
            content_type="application/json",
        )
        assert resp.status_code == 200
        data = resp.get_json()
        assert data["ok"] is True

    def test_regenerate_all_uses_db(self, app_client):
        resp = app_client.post("/api/regenerate")
        assert resp.status_code == 200
        data = resp.get_json()
        assert data["ok"] is True
        assert "results" in data

    def test_regenerate_floorplan_reflects_db_state(self, app_client):
        """Verify regenerated SVG comes from DB, not hardcoded modules.

        Delete an element from the DB, regenerate, and confirm the SVG
        content changes (element label absent from output).
        """
        # First, verify OTTO is present in seed-state SVG
        resp = app_client.post(
            "/api/regenerate",
            data=json.dumps({"view": "floorplan"}),
            content_type="application/json",
        )
        assert resp.status_code == 200
        svg_resp = app_client.get("/api/svg/floorplan")
        assert svg_resp.status_code == 200
        seed_svg = svg_resp.data.decode("utf-8")
        assert "OTTO" in seed_svg, "OTTO should appear in seed-state SVG"

        # Delete ottoman element from DB
        resp = app_client.get("/api/elements")
        elements = resp.get_json()
        ottoman = next((e for e in elements if e["name"] == "ottoman"), None)
        assert ottoman is not None, "ottoman element should exist in DB"
        app_client.delete(f"/api/elements/{ottoman['id']}")

        # Regenerate floorplan
        resp = app_client.post(
            "/api/regenerate",
            data=json.dumps({"view": "floorplan"}),
            content_type="application/json",
        )
        assert resp.status_code == 200

        # Verify OTTO is absent from the regenerated SVG
        svg_resp = app_client.get("/api/svg/floorplan")
        assert svg_resp.status_code == 200
        new_svg = svg_resp.data.decode("utf-8")
        assert "OTTO" not in new_svg, "OTTO should be absent after deletion"
