"""Flask web server for the ADU Editor application.

Provides REST API for constants/geometry/SVG, WebSocket for real-time
updates, and serves the single-page frontend.
"""
import datetime
import json
import os
import queue
import sqlite3
import subprocess
import threading
import time

from flask import Flask, jsonify, render_template, request, Response, send_file

from app.database import (
    DB_PATH, init_db, get_all_constants, get_constants_dict,
    get_constant_value, update_constant, update_constants_batch,
    get_categories, get_outline_chain, get_views, reset_constants,
    get_all_elements, get_element, get_element_by_name,
    create_element, update_element, delete_element,
    get_all_doors, get_door, create_door, update_door, delete_door,
    get_outline_chain_row, update_outline_segment, insert_outline_segment,
    delete_outline_segment, restore_outline_chain, reset_outline_chain,
    reset_elements, get_config, set_config,
    validate_db, reset_db,
    get_variants, get_variant, get_variant_by_id, update_variant,
    create_variant, delete_variant, create_variant_raw,
    clone_variant_exclusions, delete_variant_exclusions,
    set_variant_exclusion, get_variant_exclusions,
    clone_variant_elements, unclone_variant_elements,
    get_element_formulas, get_all_formulas, upsert_formula,
    delete_formula, set_formula_lock, get_formula_deps,
    get_dependents as db_get_dependents, rebuild_formula_deps,
    get_all_formula_deps,
    get_survey_legs, get_survey_config,
    update_survey_leg, update_survey_config, reset_survey,
    export_project, import_project,
    get_inner_wall_overrides, get_inner_wall_override,
    upsert_inner_wall_override, delete_inner_wall_override,
    snapshot_inner_wall_overrides, restore_inner_wall_overrides,
    reset_inner_wall_overrides, check_override_overlap,
)
from app.doors import validate_door
from app.elements import compute_constant_delta, IW_CONSTANT_MAP, IW_HOSTED_OPENINGS
from app.engine import (
    compute_geometry, generate_svg, generate_svg_db,
    build_generator_data_from_db, get_svg_content, patch_constants,
    compute_survey_points,
)
from app.plumbing import (
    get_plumbing_elements, get_plumbing_element,
    create_plumbing_element, create_plumbing_raw,
    update_plumbing_element, delete_plumbing_element,
    seed_reference_plumbing,
)
from app.outline_solver import db_rows_to_chain, solve_closure, validate_chain
from app.undo import UndoManager
import floorplan.constants as fc

_PROJECT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Build/version info — computed once at import time
_START_TIME = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
try:
    _GIT_DESCRIBE = subprocess.check_output(
        ["git", "describe", "--always", "--dirty"],
        cwd=_PROJECT, stderr=subprocess.DEVNULL, text=True,
    ).strip()
except Exception:
    _GIT_DESCRIBE = "unknown"

# SSE event bus — all connected clients get geometry update notifications
_sse_clients: list[queue.Queue] = []
_sse_lock = threading.Lock()


def _broadcast(event: str, data: dict | None = None):
    """Send an SSE event to all connected clients."""
    msg = f"event: {event}\ndata: {json.dumps(data or {})}\n\n"
    with _sse_lock:
        alive = []
        for q in _sse_clients:
            try:
                q.put_nowait(msg)
                alive.append(q)
            except queue.Full:
                pass  # drop dead/full clients
        _sse_clients[:] = alive


def _seed_reference_plumbing_if_needed(db):
    """Seed reference plumbing pipes if none exist yet."""
    existing = get_plumbing_elements(db)
    if any(e["type"] in ("supply_pipe", "drain_pipe") for e in existing):
        return
    try:
        c = get_constants_dict(db)
        ch = get_outline_chain(db)
        dr = get_all_doors(db)
        geom = compute_geometry(c, "standard", ch, doors_data=dr, db_path=db)
        seed_reference_plumbing(geom, c.get("WALL_OUTER", 10.0 / 12.0), db)
    except Exception:
        pass  # skip if geometry computation fails (e.g. corrupt DB)


def create_app(db_path=None):
    """Create and configure the Flask application."""
    db = db_path or DB_PATH
    init_db(db)

    # Validate database health on startup
    _db_ok, _db_issues = validate_db(db)
    if not _db_ok:
        print(f"WARNING: Database issues detected: {_db_issues}")

    _seed_reference_plumbing_if_needed(db)

    app = Flask(
        __name__,
        template_folder=os.path.join(os.path.dirname(__file__), "templates"),
        static_folder=os.path.join(os.path.dirname(__file__), "static"),
    )
    app.config["SEND_FILE_MAX_AGE_DEFAULT"] = 0  # no caching during dev

    # Undo/redo manager
    undo_mgr = UndoManager(db)

    # Cache for computed geometry — keyed by variant
    _geom_cache = {}  # {variant: {"data": ..., "dirty": True}}
    _geom_lock = threading.Lock()

    def _get_geometry(variant="standard"):
        """Get geometry for a variant, computing if dirty."""
        with _geom_lock:
            entry = _geom_cache.get(variant)
            if entry is None or entry["dirty"] or entry["data"] is None:
                constants = get_constants_dict(db)
                chain_rows = get_outline_chain(db)
                doors_data = get_all_doors(db)
                data = compute_geometry(constants, variant, chain_rows,
                                        doors_data=doors_data, db_path=db)
                _geom_cache[variant] = {"data": data, "dirty": False}
                return data
            return entry["data"]

    def _invalidate():
        """Mark all variant caches as dirty and notify clients."""
        with _geom_lock:
            for entry in _geom_cache.values():
                entry["dirty"] = True
        _broadcast("constants_changed")
        _broadcast("geometry_changed")
        _broadcast_undo_status()

    def _broadcast_undo_status():
        """Broadcast current undo/redo availability to all clients."""
        _broadcast("undo_status", {
            "can_undo": undo_mgr.can_undo,
            "can_redo": undo_mgr.can_redo,
        })

    def _evaluate_element(elem_name, variant="standard"):
        """Run the formula evaluator and return computed data for elem_name.

        Returns dict with 'poly' key, or None if the element has no formula
        or evaluation fails.
        """
        from app.evaluator import FormulaEvaluator
        constants = get_constants_dict(db)
        # Build a minimal geometry context for the evaluator
        geom = _get_geometry(variant)
        base_points = geom.get("points", {})
        inner = [(p[0], p[1]) for p in (geom.get("inner_poly") or [])]
        radii = geom.get("radii", {})

        ev = FormulaEvaluator(constants, base_points, inner, radii)
        ev.load_formulas_from_db(db_path=db, variant=variant)
        try:
            ev.evaluate_all()
        except Exception:
            return None
        return ev.elements.get(elem_name)

    def _solve_and_update_closure():
        """Re-solve closure and update solved values in DB.

        Returns (solver, chain_rows) on success.
        Returns (solver, None) if closure is invalid.
        """
        chain_rows = get_outline_chain(db)
        chain = db_rows_to_chain(chain_rows)
        patch_constants(get_constants_dict(db))
        R_a1 = fc.CORNER_SW_R
        solver = solve_closure(chain, R_a1)
        if not solver.valid:
            return solver, None
        update_outline_segment(0, {"distance": solver.d_F2_F5}, db)
        f18_seq = len(chain_rows) - 2
        update_outline_segment(f18_seq, {"distance": solver.d_F18_F1}, db)
        closure_seq = len(chain_rows) - 1
        update_outline_segment(closure_seq, {"sweep": solver.sweep_closure}, db)
        return solver, chain_rows

    # ------------------------------------------------------------------
    # Prevent browser caching of API JSON responses
    # ------------------------------------------------------------------

    @app.after_request
    def add_no_cache(response):
        if request.path.startswith("/api/") and "text/event-stream" not in (response.content_type or ""):
            response.headers["Cache-Control"] = "no-store"
        return response

    # ------------------------------------------------------------------
    # Routes
    # ------------------------------------------------------------------

    @app.route("/")
    def index():
        return render_template("index.html")

    # -- Constants API --

    @app.route("/api/constants")
    def api_constants():
        cat = request.args.get("category")
        constants = get_all_constants(db)
        if cat:
            constants = [c for c in constants if c["category"] == cat]
        return jsonify(constants)

    @app.route("/api/constants/categories")
    def api_categories():
        return jsonify(get_categories(db))

    @app.route("/api/constants/<name>", methods=["PUT"])
    def api_update_constant(name):
        body = request.get_json(force=True)
        value = body.get("value")
        if value is None:
            return jsonify({"error": "missing value"}), 400
        try:
            value = float(value)
        except (TypeError, ValueError):
            return jsonify({"error": "invalid value"}), 400
        old = get_constant_value(name, db)
        ok = update_constant(name, value, db)
        if ok:
            undo_mgr.record(
                "constant_update", {name: old}, {name: value},
                f"Change {name}",
            )
            _invalidate()
            return jsonify({"ok": True, "name": name, "value": value})
        return jsonify({"error": "not found"}), 404

    @app.route("/api/constants/batch", methods=["PUT"])
    def api_update_constants_batch():
        body = request.get_json(force=True)
        updates = body.get("updates", {})
        try:
            updates = {k: float(v) for k, v in updates.items()}
        except (TypeError, ValueError):
            return jsonify({"error": "invalid values"}), 400
        before = {k: v for k, v in get_constants_dict(db).items() if k in updates}
        n = update_constants_batch(updates, db)
        if n > 0:
            undo_mgr.record(
                "constant_batch", before, updates,
                f"Batch update {n} constants",
            )
            _invalidate()
        return jsonify({"ok": True, "changed": n})

    @app.route("/api/constants/reset", methods=["POST"])
    def api_reset_constants():
        before = {
            "constants": get_constants_dict(db),
            "outline": get_outline_chain(db),
            "elements": get_all_elements(db),
            "doors": get_all_doors(db),
            "survey_legs": get_survey_legs(db),
            "survey_config": get_survey_config(db),
            "inner_wall_overrides": snapshot_inner_wall_overrides(db),
        }
        reset_constants(db)
        reset_outline_chain(db)
        reset_elements(db)
        reset_survey(db)
        reset_inner_wall_overrides(db)
        after = {
            "constants": get_constants_dict(db),
            "outline": get_outline_chain(db),
            "elements": get_all_elements(db),
            "doors": get_all_doors(db),
            "survey_legs": get_survey_legs(db),
            "survey_config": get_survey_config(db),
            "inner_wall_overrides": snapshot_inner_wall_overrides(db),
        }
        undo_mgr.record("full_reset", before, after, "Reset to defaults")
        _invalidate()
        _broadcast("element_changed")
        _broadcast("outline_changed")
        return jsonify({"ok": True})

    # -- Undo/Redo API --

    @app.route("/api/undo", methods=["POST"])
    def api_undo():
        entry = undo_mgr.undo()
        if entry is None:
            return jsonify({"error": "nothing to undo"}), 400
        _invalidate()
        return jsonify({
            "ok": True,
            "action": entry["action_type"],
            "description": entry["description"],
            "can_undo": undo_mgr.can_undo,
            "can_redo": undo_mgr.can_redo,
        })

    @app.route("/api/redo", methods=["POST"])
    def api_redo():
        entry = undo_mgr.redo()
        if entry is None:
            return jsonify({"error": "nothing to redo"}), 400
        _invalidate()
        return jsonify({
            "ok": True,
            "action": entry["action_type"],
            "description": entry["description"],
            "can_undo": undo_mgr.can_undo,
            "can_redo": undo_mgr.can_redo,
        })

    # -- Elements API --

    @app.route("/api/elements")
    def api_elements():
        variant = request.args.get("variant")
        if variant:
            from app.elements import get_elements_for_variant
            return jsonify(get_elements_for_variant(variant, db))
        return jsonify(get_all_elements(db))

    @app.route("/api/elements", methods=["POST"])
    def api_create_element():
        body = request.get_json(force=True)
        type_ = body.get("type")
        name = body.get("name")
        if not type_ or not name:
            return jsonify({"error": "missing type or name"}), 400
        properties = body.get("properties", {})
        variant = body.get("variant")
        try:
            record = create_element(type_, name, properties, variant, db)
        except Exception as exc:
            return jsonify({"error": str(exc)}), 400
        undo_mgr.record(
            "element_create", record, record,
            f"Create {type_} {name}",
        )
        _broadcast("element_changed")
        _invalidate()
        return jsonify(record), 201

    @app.route("/api/elements/<int:element_id>", methods=["PUT"])
    def api_update_element(element_id):
        old = get_element(element_id, db)
        if not old:
            return jsonify({"error": "not found"}), 404
        body = request.get_json(force=True)
        updated = update_element(element_id, body, db)
        undo_mgr.record(
            "element_update", old, updated,
            f"Update {old['name']}",
        )
        _broadcast("element_changed")
        _invalidate()
        return jsonify(updated)

    @app.route("/api/elements/<int:element_id>", methods=["DELETE"])
    def api_delete_element(element_id):
        old = get_element(element_id, db)
        if not old:
            return jsonify({"error": "not found"}), 404
        # Collect all records that will be deleted (for undo)
        deleted_records = [old]
        # Check for cascade targets (openings hosted by this wall)
        if old["type"] == "wall":
            from app.database import get_db as _get_db
            with _get_db(db) as conn:
                hosted = conn.execute(
                    "SELECT id, type, name, properties, variant FROM elements "
                    "WHERE type = 'opening' AND "
                    "json_extract(properties, '$.host_wall') = ?",
                    (old["name"],),
                ).fetchall()
                for h in hosted:
                    deleted_records.append(dict(h))
        deleted_ids = delete_element(element_id, db)
        undo_mgr.record(
            "element_delete", deleted_records, {"ids": deleted_ids},
            f"Delete {old['name']}",
        )
        _broadcast("element_changed")
        _invalidate()
        return jsonify({"ok": True, "deleted": deleted_ids})

    @app.route("/api/elements/<int:element_id>/move", methods=["POST"])
    def api_move_element(element_id):
        el = get_element(element_id, db)
        if not el:
            return jsonify({"error": "not found"}), 404
        body = request.get_json(force=True)

        # Resolve anchor-based format to dx/dy
        if "anchor" in body:
            anchor_name = body["anchor"]
            face = body.get("face")
            offset_in = body.get("offset", 0)
            if not face or face not in ("north", "south", "east", "west"):
                return jsonify({"error": "invalid face"}), 400
            anchor_el = get_element_by_name(anchor_name, db)
            if not anchor_el:
                return jsonify({"error": f"anchor '{anchor_name}' not found"}), 404
            # Get anchor bbox from geometry
            variant = el.get("variant")
            geom = _get_geometry(variant)
            el_geom = geom.get("elements", {}).get(anchor_name)
            if not el_geom or "bbox" not in el_geom:
                return jsonify({"error": f"no geometry for anchor '{anchor_name}'"}), 400
            abbox = el_geom["bbox"]  # [min_x, min_y, max_x, max_y]
            # Current element geometry
            cur_geom = geom.get("elements", {}).get(el["name"])
            if not cur_geom or "bbox" not in cur_geom:
                return jsonify({"error": "no geometry for this element"}), 400
            cbbox = cur_geom["bbox"]
            offset_ft = offset_in / 12.0
            # Compute target position based on anchor face
            if face == "east":
                target_x = abbox[2] + offset_ft
                dx = target_x - cbbox[0]
                dy = 0
            elif face == "west":
                target_x = abbox[0] - offset_ft
                dx = target_x - cbbox[2]
                dy = 0
            elif face == "north":
                target_y = abbox[3] + offset_ft
                dx = 0
                dy = target_y - cbbox[1]
            else:  # south
                target_y = abbox[1] - offset_ft
                dx = 0
                dy = target_y - cbbox[3]
        else:
            dx = body.get("dx", 0)
            dy = body.get("dy", 0)
            if dx == 0 and dy == 0:
                return jsonify({"error": "no movement specified"}), 400

        name = el["name"]

        # Case 1: IW wall — update controlling constant
        if el["type"] == "wall" and name in IW_CONSTANT_MAP:
            result = compute_constant_delta(name, dx, dy)
            if result is None:
                return jsonify({"error": f"wall {name} is not movable"}), 400
            const_name, delta = result
            old_val = get_constant_value(const_name, db)
            if old_val is None:
                return jsonify({"error": f"constant {const_name} not found"}), 500
            new_val = old_val + delta
            update_constant(const_name, new_val, db)
            undo_mgr.record(
                "element_move",
                {"move_type": "constant", "constant": const_name, "value": old_val},
                {"move_type": "constant", "constant": const_name, "value": new_val},
                f"Move {name}",
            )
            _invalidate()
            return jsonify({
                "ok": True, "constant": const_name,
                "old_value": old_val, "new_value": new_val,
                "can_undo": undo_mgr.can_undo,
                "can_redo": undo_mgr.can_redo,
            })

        # Case 2: Custom/override element with offset properties
        props = el.get("properties") or {}
        if isinstance(props, str):
            props = json.loads(props)
        old_ox = props.get("offset_x", 0)
        old_oy = props.get("offset_y", 0)
        new_ox = old_ox + dx
        new_oy = old_oy + dy
        old_props = dict(props)
        new_props = dict(props, offset_x=new_ox, offset_y=new_oy)
        update_element(element_id, {"properties": new_props}, db)
        undo_mgr.record(
            "element_move",
            {"move_type": "position", "id": element_id, "properties": old_props},
            {"move_type": "position", "id": element_id, "properties": new_props},
            f"Move {name}",
        )
        _broadcast("element_changed")
        _invalidate()
        return jsonify({
            "ok": True,
            "offset_x": new_ox, "offset_y": new_oy,
            "can_undo": undo_mgr.can_undo,
            "can_redo": undo_mgr.can_redo,
        })

    # -- Openings API (stored as type='opening' elements) --

    @app.route("/api/openings", methods=["POST"])
    def api_create_opening():
        body = request.get_json(force=True)
        name = body.get("name")
        segment = body.get("segment")
        if not name or not segment:
            return jsonify({"error": "missing name or segment"}), 400
        properties = {
            "segment": segment,
            "width": body.get("width", 0),
            "offset": body.get("offset", 0),
            "host_wall": body.get("host_wall"),
        }
        variant = body.get("variant")
        try:
            record = create_element("opening", name, properties, variant, db)
        except Exception as exc:
            return jsonify({"error": str(exc)}), 400
        undo_mgr.record(
            "opening_create", record, record,
            f"Create opening {name}",
        )
        _broadcast("element_changed")
        _invalidate()
        return jsonify(record), 201

    @app.route("/api/openings/<name>", methods=["PUT"])
    def api_update_opening(name):
        old = get_element_by_name(name, db)
        if not old or old["type"] != "opening":
            return jsonify({"error": "not found"}), 404
        body = request.get_json(force=True)
        # Merge updates into properties
        props = json.loads(old["properties"]) if isinstance(old["properties"], str) else old["properties"]
        for k in ("width", "offset", "segment", "host_wall"):
            if k in body:
                props[k] = body[k]
        updated = update_element(old["id"], {"properties": props}, db)
        undo_mgr.record(
            "opening_update", old, updated,
            f"Update opening {name}",
        )
        _broadcast("element_changed")
        _invalidate()
        return jsonify(updated)

    @app.route("/api/openings/<name>", methods=["DELETE"])
    def api_delete_opening(name):
        old = get_element_by_name(name, db)
        if not old or old["type"] != "opening":
            return jsonify({"error": "not found"}), 404
        # Also collect door record if exists
        deleted_records = [old]
        door = get_door(name, db)
        if door:
            deleted_records.append({"_type": "door", **door})
            delete_door(name, db)
        deleted_ids = delete_element(old["id"], db)
        undo_mgr.record(
            "opening_delete", deleted_records, {"ids": deleted_ids},
            f"Delete opening {name}",
        )
        _broadcast("element_changed")
        _invalidate()
        return jsonify({"ok": True, "deleted": deleted_ids})

    # -- Doors API --

    @app.route("/api/doors")
    def api_doors():
        return jsonify(get_all_doors(db))

    @app.route("/api/doors", methods=["POST"])
    def api_create_door():
        body = request.get_json(force=True)
        opening = body.get("opening_name") or body.get("opening")
        if not opening:
            return jsonify({"error": "missing opening_name"}), 400
        hinge = body.get("hinge_side", "east")
        swing = body.get("swing_direction", "south")
        width = body.get("width", 36)
        dtype = body.get("door_type", "single")
        err = validate_door(hinge, swing, dtype)
        if err:
            return jsonify({"error": err}), 400
        try:
            record = create_door(opening, width, hinge, swing, dtype, db)
        except Exception as exc:
            return jsonify({"error": str(exc)}), 400
        undo_mgr.record(
            "door_create", record, record,
            f"Create door on {opening}",
        )
        _broadcast("element_changed")
        _invalidate()
        return jsonify(record), 201

    @app.route("/api/doors/<opening_name>", methods=["PUT"])
    def api_update_door(opening_name):
        old = get_door(opening_name, db)
        if not old:
            return jsonify({"error": "not found"}), 404
        body = request.get_json(force=True)
        if "hinge_side" in body or "swing_direction" in body or "door_type" in body:
            err = validate_door(
                body.get("hinge_side", old["hinge_side"]),
                body.get("swing_direction", old["swing_direction"]),
                body.get("door_type", old["door_type"]),
            )
            if err:
                return jsonify({"error": err}), 400
        updated = update_door(opening_name, body, db)
        undo_mgr.record(
            "door_update", old, updated,
            f"Update door on {opening_name}",
        )
        _broadcast("element_changed")
        _invalidate()
        return jsonify(updated)

    @app.route("/api/doors/<opening_name>", methods=["DELETE"])
    def api_delete_door(opening_name):
        old = get_door(opening_name, db)
        if not old:
            return jsonify({"error": "not found"}), 404
        delete_door(opening_name, db)
        undo_mgr.record(
            "door_delete", old, {"opening_name": opening_name},
            f"Delete door on {opening_name}",
        )
        _broadcast("element_changed")
        _invalidate()
        return jsonify({"ok": True, "opening_name": opening_name})

    # -- Geometry API --

    @app.route("/api/geometry")
    def api_geometry():
        try:
            variants = get_variants(db)
            valid_names = {v["name"] for v in variants}
            variant = request.args.get("variant", "standard")
            if variant not in valid_names:
                variant = "standard"
            geom = _get_geometry(variant)
            return jsonify(geom)
        except Exception as exc:
            msg = str(exc)
            is_db = (isinstance(exc, (sqlite3.OperationalError, IndexError))
                     or "no such table" in msg)
            resp = {"error": msg}
            if is_db:
                ok, issues = validate_db(db)
                resp["db_issue"] = True
                resp["db_issues"] = issues
                resp["hint"] = ("Database may be corrupt or incomplete. "
                                "Use File \u203a Reset Database to recreate it.")
            return jsonify(resp), 500

    @app.route("/api/variants")
    def api_variants():
        return jsonify(get_variants(db))

    @app.route("/api/variants/<int:variant_id>", methods=["PUT"])
    def api_update_variant(variant_id):
        data = request.get_json(force=True)
        old = get_variant_by_id(variant_id, db)
        if not old:
            return jsonify({"error": "variant not found"}), 404
        updated = update_variant(variant_id, data, db)
        # Record undo with only changed fields
        before = {"id": variant_id}
        after = {"id": variant_id}
        for k in data:
            if k in old:
                before[k] = old[k]
            if k in updated:
                after[k] = updated[k]
        undo_mgr.record("variant_update", before, after,
                        f"Update variant {old['name']}")
        _broadcast("variant_changed")
        return jsonify(updated)

    @app.route("/api/variants", methods=["POST"])
    def api_create_variant():
        body = request.get_json(force=True)
        name = body.get("name")
        if not name:
            return jsonify({"error": "name is required"}), 400
        label = body.get("label", name)
        source = body.get("source_variant", "standard")
        # Validate source exists
        source_v = get_variant(source, db)
        if not source_v:
            return jsonify({"error": f"source variant '{source}' not found"}), 400
        # Check for duplicate name
        if get_variant(name, db):
            return jsonify({"error": f"variant '{name}' already exists"}), 400
        flags = source_v["flags"]
        record = create_variant(name, label, source, flags, db)
        clone_variant_exclusions(source, name, db)
        clone_variant_elements(source, name, db)
        undo_mgr.record("variant_create", record, record,
                        f"Create variant {name}")
        _invalidate()
        _broadcast("variant_changed")
        return jsonify(record), 201

    @app.route("/api/variants/<int:variant_id>", methods=["DELETE"])
    def api_delete_variant(variant_id):
        v = get_variant_by_id(variant_id, db)
        if not v:
            return jsonify({"error": "not found"}), 404
        if v["is_builtin"]:
            return jsonify({"error": "cannot delete built-in variant"}), 400
        unclone_variant_elements(v["name"], db)
        delete_variant_exclusions(v["name"], db)
        delete_variant(variant_id, db)
        undo_mgr.record("variant_delete", v, v,
                        f"Delete variant {v['name']}")
        _invalidate()
        _broadcast("variant_changed")
        return jsonify({"deleted": variant_id})

    # -- Exclusions API --

    @app.route("/api/exclusions")
    def api_get_exclusions():
        variant = request.args.get("variant")
        if variant:
            excl = get_variant_exclusions(variant, db)
            return jsonify({k: sorted(v) for k, v in excl.items()})
        # No variant param → return all variants' exclusions
        all_excl = {}
        for v in get_variants(db):
            excl = get_variant_exclusions(v["name"], db)
            if excl:
                all_excl[v["name"]] = {k: sorted(vs) for k, vs in excl.items()}
        return jsonify(all_excl)

    @app.route("/api/exclusions", methods=["PUT"])
    def api_set_exclusion():
        body = request.get_json(force=True)
        variant = body.get("variant")
        element_type = body.get("element_type")
        element_name = body.get("element_name")
        excluded = body.get("excluded", True)
        if not variant or not element_type or not element_name:
            return jsonify({"error": "variant, element_type, element_name required"}), 400
        changed = set_variant_exclusion(
            variant, element_type, element_name, excluded, db)
        # Cascade: when a wall is excluded/included, also exclude/include
        # its hosted rough openings so they stay in sync.
        if element_type == "wall":
            for ro_name in IW_HOSTED_OPENINGS.get(element_name, []):
                if set_variant_exclusion(
                        variant, "rough_opening", ro_name, excluded, db):
                    changed = True
        if changed:
            _invalidate()
            _broadcast("element_changed")
        return jsonify({"ok": True, "changed": changed})

    # -- Shapes API --

    @app.route("/api/shapes")
    def api_shapes():
        from app.database import get_shapes
        return jsonify(get_shapes(db))

    @app.route("/api/shapes", methods=["POST"])
    def api_create_shape():
        from app.database import create_shape
        body = request.get_json(force=True)
        name = body.get("name")
        poly_json = body.get("poly_json")
        if not name or poly_json is None:
            return jsonify({"error": "name and poly_json required"}), 400
        try:
            shape = create_shape(
                name, poly_json,
                scale=body.get("scale", 1.0),
                origin=body.get("origin", "center"),
                width_key=body.get("width_key"),
                depth_key=body.get("depth_key"),
                description=body.get("description", ""),
                db_path=db,
            )
            return jsonify(shape), 201
        except Exception as exc:
            return jsonify({"error": str(exc)}), 400

    @app.route("/api/shapes/<name>", methods=["PUT"])
    def api_update_shape(name):
        from app.database import get_shape, update_shape
        if not get_shape(name, db):
            return jsonify({"error": "not found"}), 404
        body = request.get_json(force=True)
        update_shape(name, db_path=db, **{
            k: body[k] for k in
            ("poly_json", "scale", "origin", "width_key", "depth_key", "description")
            if k in body
        })
        return jsonify(get_shape(name, db))

    # -- Outline chain API --

    @app.route("/api/outline")
    def api_outline():
        return jsonify(get_outline_chain(db))

    @app.route("/api/outline/<int:seq>", methods=["PUT"])
    def api_update_outline(seq):
        """API-16: Update outline chain segment, re-solve closure."""


        old_row = get_outline_chain_row(seq, db)
        if not old_row:
            return jsonify({"error": f"segment {seq} not found"}), 404

        body = request.get_json(force=True)
        updates = {}

        # Map dist_or_radius to distance or radius based on seg_type
        if "dist_or_radius" in body:
            val = float(body["dist_or_radius"])
            if old_row["seg_type"] == "L":
                updates["distance"] = val
            else:
                updates["radius"] = val
        if "sweep" in body:
            updates["sweep"] = float(body["sweep"])
            updates["sweep_name"] = str(body["sweep"])

        if not updates:
            return jsonify({"error": "no valid fields to update"}), 400

        # Protect solved segments
        chain_rows = get_outline_chain(db)
        n = len(chain_rows)
        solved_seqs_dist = {0, n - 2}  # F2→F5 and F18→F1
        if seq in solved_seqs_dist and "distance" in updates:
            return jsonify({"error": "Cannot directly edit solved distance"}), 400
        if seq == n - 1 and "sweep" in updates:
            return jsonify({"error": "Cannot directly edit closure arc sweep"}), 400

        # Snapshot before state
        before_chain = get_outline_chain(db)

        # Apply update to DB
        update_outline_segment(seq, updates, db)

        # Re-solve closure
        solver, _ = _solve_and_update_closure()
        if not solver.valid:
            restore_outline_chain(before_chain, db)
            return jsonify({
                "error": "Closure failed after edit",
                "closure_error": solver.closure_error,
            }), 400

        # Record undo
        after_chain = get_outline_chain(db)
        undo_mgr.record(
            "outline_update", before_chain, after_chain,
            f"Edit outline seg {seq}",
        )

        _invalidate()
        _broadcast("outline_changed")

        updated = get_outline_chain_row(seq, db)
        return jsonify({
            "ok": True,
            "segment": updated,
            "closure_valid": True,
            "d_F2_F5": solver.d_F2_F5,
            "d_F18_F1": solver.d_F18_F1,
            "sweep_closure": solver.sweep_closure,
        })

    @app.route("/api/outline/validate", methods=["POST"])
    def api_validate_outline():
        """API-17: Dry-run validation without committing."""


        body = request.get_json(force=True)
        changes = body.get("changes", {})  # {seq_str: {field: value}}

        # Get current chain and apply proposed changes in memory
        chain_rows = get_outline_chain(db)
        for seq_str, upd in changes.items():
            seq_i = int(seq_str)
            for row in chain_rows:
                if row["seq"] == seq_i:
                    for k, v in upd.items():
                        if k == "dist_or_radius":
                            if row["seg_type"] == "L":
                                row["distance"] = float(v)
                            else:
                                row["radius"] = float(v)
                        elif k in ("sweep", "distance", "radius"):
                            row[k] = float(v)

        chain = db_rows_to_chain(chain_rows)

        patch_constants(get_constants_dict(db))
        R_a1 = fc.CORNER_SW_R
        result = validate_chain(chain, R_a1)

        return jsonify(result)

    @app.route("/api/outline/add-point", methods=["POST"])
    def api_add_outline_point():
        """API-18: Insert new F-point by splitting a segment."""


        body = request.get_json(force=True)
        after_seq = body.get("after_seq")
        end_name = body.get("end_name")
        seg_type = body.get("seg_type", "L")

        if after_seq is None or not end_name:
            return jsonify({"error": "missing after_seq or end_name"}), 400

        before_chain = get_outline_chain(db)

        old_seg = get_outline_chain_row(after_seq, db)
        if not old_seg:
            return jsonify({"error": f"segment {after_seq} not found"}), 404

        # Check for duplicate end_name
        for row in before_chain:
            if row["end_name"] == end_name:
                return jsonify({"error": f"point name {end_name} already exists"}), 400

        if old_seg["seg_type"] == "L":
            # Line split: halve distance
            half_dist = (old_seg["distance"] or 0) / 2.0
            update_outline_segment(after_seq, {
                "distance": half_dist, "end_name": end_name,
            }, db)
            new_row = {
                "seg_type": "L",
                "distance": half_dist,
                "end_name": old_seg["end_name"],
            }
            insert_outline_segment(after_seq + 1, new_row, db)
        else:
            # Arc split: halve sweep, keep radius
            half_sweep = (old_seg["sweep"] or 0) / 2.0
            center_name = body.get("center_name", old_seg["center_name"])
            update_outline_segment(after_seq, {
                "sweep": half_sweep, "end_name": end_name,
            }, db)
            new_row = {
                "seg_type": old_seg["seg_type"],
                "radius": old_seg["radius"],
                "sweep": half_sweep,
                "sweep_name": None,
                "center_name": center_name,
                "n_pts": old_seg["n_pts"],
                "end_name": old_seg["end_name"],
            }
            insert_outline_segment(after_seq + 1, new_row, db)

        # Re-solve closure
        solver, _ = _solve_and_update_closure()
        if not solver.valid:
            restore_outline_chain(before_chain, db)
            return jsonify({
                "error": "Closure failed after adding point",
                "closure_error": solver.closure_error,
            }), 400

        after_chain = get_outline_chain(db)
        undo_mgr.record(
            "outline_add_point", before_chain, after_chain,
            f"Add point {end_name}",
        )

        _invalidate()
        _broadcast("outline_changed")

        return jsonify({
            "ok": True,
            "chain": get_outline_chain(db),
            "closure_valid": True,
        })

    @app.route("/api/outline/<int:seq>", methods=["DELETE"])
    def api_delete_outline_point(seq):
        """API-19: Remove a point from the outline chain."""


        chain_rows = get_outline_chain(db)
        n = len(chain_rows)

        if seq < 0 or seq >= n:
            return jsonify({"error": "invalid seq"}), 400

        # Protect closure arc (last segment, F1->F2)
        if seq == n - 1:
            return jsonify({"error": "Cannot delete the closure arc"}), 400

        # Protect the chain from becoming too small (need at least 3 segments)
        if n <= 3:
            return jsonify({"error": "Chain too small to remove segments"}), 400

        before_chain = list(chain_rows)

        # For consecutive line merging
        deleted_row = chain_rows[seq]
        if (deleted_row["seg_type"] == "L" and seq > 0
                and chain_rows[seq - 1]["seg_type"] == "L"):
            merged_dist = ((chain_rows[seq - 1].get("distance") or 0)
                           + (deleted_row.get("distance") or 0))
            update_outline_segment(seq - 1, {
                "distance": merged_dist,
                "end_name": deleted_row["end_name"],
            }, db)

        delete_outline_segment(seq, db)

        # Re-solve closure
        solver, _ = _solve_and_update_closure()
        if not solver.valid:
            restore_outline_chain(before_chain, db)
            return jsonify({
                "error": "Closure failed after removing point",
                "closure_error": solver.closure_error,
            }), 400

        after_chain = get_outline_chain(db)
        undo_mgr.record(
            "outline_remove_point", before_chain, after_chain,
            f"Remove seg {seq}",
        )

        _invalidate()
        _broadcast("outline_changed")

        return jsonify({
            "ok": True,
            "chain": get_outline_chain(db),
            "closure_valid": True,
        })

    # -- Views & SVG API --

    @app.route("/api/views")
    def api_views():
        return jsonify(get_views(db))

    _FLOORPLAN_VARIANT_SUFFIX = {
        "standard": "", "minik": "_minik", "daybed": "_db",
        "bare": "_bare", "sf": "_sf",
    }

    @app.route("/api/svg/<view_name>")
    def api_svg(view_name):
        views = get_views(db)
        view = next((v for v in views if v["name"] == view_name), None)
        if not view:
            return jsonify({"error": "unknown view"}), 404
        svg_path = view["svg_path"]
        # Support variant-specific floorplan SVGs
        variant = request.args.get("variant", "standard")
        if view_name == "floorplan":
            suffix = _FLOORPLAN_VARIANT_SUFFIX.get(variant)
            if suffix is None:
                # User-defined variant: fall back to source variant's suffix
                v = get_variant(variant, db)
                source = v["source_variant"] if v else "standard"
                suffix = _FLOORPLAN_VARIANT_SUFFIX.get(source, "")
            svg_path = svg_path.replace(".svg", f"{suffix}.svg")
        content = get_svg_content(svg_path)
        if content is None:
            return jsonify({"error": "SVG not found — run regenerate first"}), 404
        return Response(content, mimetype="image/svg+xml")

    @app.route("/api/svg/<view_name>/file")
    def api_svg_file(view_name):
        views = get_views(db)
        view = next((v for v in views if v["name"] == view_name), None)
        if not view:
            return jsonify({"error": "unknown view"}), 404
        full_path = os.path.join(_PROJECT, view["svg_path"])
        if not os.path.exists(full_path):
            return jsonify({"error": "file not found"}), 404
        if full_path.endswith(".pdf"):
            mime = "application/pdf"
        elif full_path.endswith(".png"):
            mime = "image/png"
        else:
            mime = "image/svg+xml"
        return send_file(full_path, mimetype=mime)

    # -- Span Analysis API (ANALYSIS-1, ANALYSIS-2) --

    @app.route("/api/span-data")
    def api_span_data():
        from app.engine import compute_span_data
        constants = get_constants_dict(db)
        data = compute_span_data(constants, db_path=db)
        # Downsample to every 6th point (~2-inch resolution)
        step = 6
        return jsonify({
            "eastings": data["eastings"][::step],
            "spans": data["spans"][::step],
            "south_spans": data["south_spans"][::step],
            "north_spans": data["north_spans"][::step],
        })

    @app.route("/api/span-rotation")
    def api_span_rotation():
        from app.engine import compute_span_rotation
        constants = get_constants_dict(db)
        return jsonify(compute_span_rotation(constants, db_path=db))

    # -- Config API (SCAD-2) --

    @app.route("/api/config/<key>")
    def api_config_get(key):
        value = get_config(key, db)
        if value is None:
            return jsonify({"error": "unknown config key"}), 404
        return jsonify({"key": key, "value": value})

    @app.route("/api/config/<key>", methods=["PUT"])
    def api_config_set(key):
        body = request.get_json(force=True)
        value = body.get("value")
        if value is None:
            return jsonify({"error": "value required"}), 400
        set_config(key, value, db)
        return jsonify({"ok": True, "key": key, "value": str(value)})

    # -- Database Reset API --

    @app.route("/api/reset-database", methods=["POST"])
    def api_reset_database():
        reset_db(db)
        _seed_reference_plumbing_if_needed(db)
        ok, issues = validate_db(db)
        _invalidate()
        _broadcast("element_changed")
        _broadcast("outline_changed")
        _broadcast("plumbing_changed")
        return jsonify({"ok": ok, "issues": issues})

    # -- Project Export/Import API (Phase 14-D) --

    @app.route("/api/project/export")
    def api_project_export():
        data = export_project(db)
        return jsonify(data)

    @app.route("/api/project/import", methods=["POST"])
    def api_project_import():
        data = request.get_json(force=True)
        # Capture current state for undo
        before = export_project(db)
        try:
            import_project(data, db)
        except ValueError as e:
            return jsonify({"error": str(e)}), 400
        after = export_project(db)
        undo_mgr.record("project_import", before, after, "Import project")
        _invalidate()
        _broadcast("element_changed")
        _broadcast("outline_changed")
        return jsonify({"ok": True})

    # -- Survey Points API (SITE-4) --

    @app.route("/api/survey-points")
    def api_survey_points():
        constants = get_constants_dict(db)
        return jsonify(compute_survey_points(constants))

    # -- Survey Data API (Phase 14-B) --

    @app.route("/api/survey/legs")
    def api_survey_legs():
        return jsonify(get_survey_legs(db))

    @app.route("/api/survey/config")
    def api_survey_config():
        return jsonify(get_survey_config(db))

    @app.route("/api/survey/legs/<int:seq>", methods=["PUT"])
    def api_update_survey_leg(seq):
        data = request.get_json(force=True)
        before = {"legs": get_survey_legs(db)}
        update_survey_leg(seq, data, db)
        after = {"legs": get_survey_legs(db)}
        undo_mgr.record("survey_leg_update", before, after,
                        f"Update survey leg {seq}")
        _invalidate()
        _broadcast("geometry_changed")
        return jsonify({"ok": True})

    @app.route("/api/survey/config/<key>", methods=["PUT"])
    def api_update_survey_config(key):
        data = request.get_json(force=True)
        before = {"config": get_survey_config(db)}
        update_survey_config(key, data["value"], db)
        after = {"config": get_survey_config(db)}
        undo_mgr.record("survey_config_update", before, after,
                        f"Update survey config {key}")
        _invalidate()
        _broadcast("geometry_changed")
        return jsonify({"ok": True})

    @app.route("/api/survey/reset", methods=["POST"])
    def api_reset_survey():
        before = {"legs": get_survey_legs(db), "config": get_survey_config(db)}
        reset_survey(db)
        after = {"legs": get_survey_legs(db), "config": get_survey_config(db)}
        undo_mgr.record("survey_reset", before, after, "Reset survey data")
        _invalidate()
        _broadcast("geometry_changed")
        return jsonify({"ok": True})

    # -- Inner Wall Overrides API (Phase 15½-C) --

    @app.route("/api/inner-wall-overrides")
    def api_inner_wall_overrides():
        """Get all inner wall overrides."""
        return jsonify(get_inner_wall_overrides(db))

    @app.route("/api/inner-wall-overrides/<int:seg_index>")
    def api_inner_wall_override(seg_index):
        """Get override chain for a single segment."""
        chain = get_inner_wall_override(seg_index, db)
        if not chain:
            return jsonify({"error": f"No override for segment {seg_index}"}), 404
        return jsonify(chain)

    @app.route("/api/inner-wall-overrides/<int:seg_index>/compute-default")
    def api_compute_default_override(seg_index):
        """Compute the default parametric chain for a segment from geometry."""
        from app.gen_provider import (
            compute_native_geometry, compute_default_override,
            compute_default_span_override,
        )
        span_end = request.args.get("span_end", type=int)
        constants = get_constants_dict(db)
        chain_rows = get_outline_chain(db)
        pts, outline_segs, inner_segs, radii = compute_native_geometry(
            constants, chain_rows=chain_rows, db_path=db)
        if span_end is not None:
            chain = compute_default_span_override(
                seg_index, span_end, inner_segs, pts, constants)
        else:
            chain = compute_default_override(
                seg_index, inner_segs, pts, constants)
        if not chain:
            return jsonify({"error": f"No default for segment {seg_index}"}), 404
        return jsonify({"seg_index": seg_index, "span_end": span_end,
                        "chain": chain})

    @app.route("/api/inner-wall-overrides/<int:seg_index>", methods=["PUT"])
    def api_upsert_inner_wall_override(seg_index):
        """Create or update an inner wall override chain."""
        body = request.get_json(force=True)
        chain = body.get("chain")
        span_end = body.get("span_end")  # None for single-segment
        if not chain or not isinstance(chain, list):
            return jsonify({"error": "chain must be a non-empty list"}), 400

        # Validate span_end
        if span_end is not None:
            if not isinstance(span_end, int) or span_end < seg_index:
                return jsonify({"error": "span_end must be >= seg_index"}), 400

        # Check overlap with existing overrides
        conflicts = check_override_overlap(seg_index, span_end, db)
        if conflicts:
            return jsonify({
                "error": f"Span overlaps existing override(s) at seg {conflicts}"
            }), 400

        # Validate chain entries
        for i, sub in enumerate(chain):
            st = sub.get("seg_type")
            if st not in ("L", "CW", "CCW"):
                return jsonify({"error": f"sub[{i}]: invalid seg_type '{st}'"}), 400
            if st == "L":
                if sub.get("bearing") is None or sub.get("distance") is None:
                    return jsonify({"error": f"sub[{i}]: line requires bearing and distance"}), 400
            else:
                if sub.get("radius") is None or sub.get("sweep") is None:
                    return jsonify({"error": f"sub[{i}]: arc requires radius and sweep"}), 400

        before = snapshot_inner_wall_overrides(db)
        upsert_inner_wall_override(seg_index, chain, span_end=span_end,
                                   db_path=db)
        after = snapshot_inner_wall_overrides(db)

        label = (f"Set inner wall override seg {seg_index}-{span_end}"
                 if span_end is not None
                 else f"Set inner wall override seg {seg_index}")
        undo_mgr.record("inner_wall_override_upsert", before, after, label)

        # Endpoint validation (warnings, not blocking)
        warnings = []
        try:
            from app.gen_provider import (
                compute_native_geometry, validate_override_endpoint,
                _seg_start_bearing,
            )
            constants = get_constants_dict(db)
            chain_rows = get_outline_chain(db)
            pts, outline_segs, inner_segs, radii = compute_native_geometry(
                constants, chain_rows=chain_rows, db_path=db)
            seg = inner_segs[seg_index]
            end_idx = span_end if span_end is not None else seg_index
            end_seg = inner_segs[end_idx]
            start_pt = pts[seg.start]
            end_pt = pts[end_seg.end]
            start_brg = _seg_start_bearing(seg, pts)
            # Expected exit bearing: bearing at start of the NEXT segment
            next_idx = (end_idx + 1) % len(inner_segs)
            expected_exit = _seg_start_bearing(inner_segs[next_idx], pts)
            result = validate_override_endpoint(
                chain, start_pt, start_brg, end_pt, expected_exit)
            warnings = result.get("warnings", [])
        except Exception:
            pass  # Don't block save on validation errors

        _invalidate()
        _broadcast("geometry_changed")
        return jsonify({"ok": True, "seg_index": seg_index,
                        "span_end": span_end,
                        "warnings": warnings,
                        "chain": get_inner_wall_override(seg_index, db)})

    @app.route("/api/inner-wall-overrides/<int:seg_index>", methods=["DELETE"])
    def api_delete_inner_wall_override(seg_index):
        """Delete an inner wall override, reverting to default computation."""
        existing = get_inner_wall_override(seg_index, db)
        if not existing:
            return jsonify({"error": f"No override for segment {seg_index}"}), 404

        before = snapshot_inner_wall_overrides(db)
        delete_inner_wall_override(seg_index, db)
        after = snapshot_inner_wall_overrides(db)

        undo_mgr.record("inner_wall_override_delete", before, after,
                        f"Remove inner wall override seg {seg_index}")

        _invalidate()
        _broadcast("geometry_changed")
        return jsonify({"ok": True})

    # -- Plumbing Elements API (PLUMB-5) --

    @app.route("/api/plumbing")
    def api_plumbing():
        return jsonify(get_plumbing_elements(db))

    @app.route("/api/plumbing", methods=["POST"])
    def api_create_plumbing():
        data = request.get_json(force=True)
        elem = create_plumbing_element(
            data["type"], data["name"],
            path=data.get("path"),
            properties=data.get("properties"),
            fixture=data.get("fixture"),
            db_path=db,
        )
        undo_mgr.record("plumbing_create", {"id": elem["id"]}, elem,
                        f"Create plumbing {elem['name']}")
        _broadcast("plumbing_changed")
        return jsonify(elem), 201

    @app.route("/api/plumbing/<int:element_id>", methods=["PUT"])
    def api_update_plumbing(element_id):
        old = get_plumbing_element(element_id, db)
        if not old:
            return jsonify({"error": "not found"}), 404
        data = request.get_json(force=True)
        updated = update_plumbing_element(element_id, data, db)
        undo_mgr.record("plumbing_update",
                        {"id": element_id, **{k: old[k] for k in data if k in old}},
                        {"id": element_id, **{k: updated[k] for k in data if k in updated}},
                        f"Update plumbing {old['name']}")
        _broadcast("plumbing_changed")
        return jsonify(updated)

    @app.route("/api/plumbing/<int:element_id>", methods=["DELETE"])
    def api_delete_plumbing(element_id):
        old = get_plumbing_element(element_id, db)
        if not old:
            return jsonify({"error": "not found"}), 404
        delete_plumbing_element(element_id, db)
        undo_mgr.record("plumbing_delete", old, {"id": element_id},
                        f"Delete plumbing {old['name']}")
        _broadcast("plumbing_changed")
        return jsonify({"deleted": element_id})

    # -- Site Plan Generation API (SITE-1) --

    @app.route("/api/generate-site-plan", methods=["POST"])
    def api_generate_site_plan():
        gd = build_generator_data_from_db(db)
        ok = generate_svg_db("site_plan", "site/gen_site_plan.py", gd, db_path=db)
        return jsonify({
            "ok": ok,
            "setback_216": get_config("setback_216", db),
            "setback_275": get_config("setback_275", db),
        })

    # -- 3D Generation API (SCAD-1, SCAD-3) --

    _ROOF_SCRIPTS = {
        "flat": "scad/gen_flat_roof.py",
        "2in12": "scad/gen_2in12.py",
    }

    @app.route("/api/generate-3d", methods=["POST"])
    def api_generate_3d():
        roof_style = get_config("roof_style", db) or "flat"
        script = _ROOF_SCRIPTS.get(roof_style)
        if not script:
            return jsonify({"ok": False, "error": f"unknown roof style: {roof_style}"}), 400
        gd = build_generator_data_from_db(db)
        ok = generate_svg_db("3d_" + roof_style, script, gd, db_path=db)
        return jsonify({
            "ok": ok,
            "roof_style": roof_style,
            "output": f"scad/{'flat_roof' if roof_style == 'flat' else '2in12'}.scad",
        })

    @app.route("/api/generate-views", methods=["POST"])
    def api_generate_views():
        roof_style = get_config("roof_style", db) or "flat"
        gd = build_generator_data_from_db(db)
        # 1. Generate SCAD file for selected roof style
        scad_script = _ROOF_SCRIPTS.get(roof_style, "scad/gen_flat_roof.py")
        ok = generate_svg_db("scad", scad_script, gd, db_path=db)
        if not ok:
            return jsonify({"ok": False, "error": "SCAD generation failed"})
        # 2. Render views (requires OpenSCAD CLI) — subprocess only
        ok = generate_svg("views", "scad/gen_views.py")
        if not ok:
            return jsonify({"ok": False, "error": "View rendering failed (OpenSCAD required)"})
        # 3. Generate line drawings — subprocess only
        generate_svg("line_drawings", "scad/gen_line_drawings.py")
        # 4. Compose 3-view PDF — subprocess only
        ok = generate_svg("3views", "gen_3views.py")
        return jsonify({
            "ok": ok,
            "roof_style": roof_style,
            "output": "3views.pdf",
        })

    # -- Regeneration API --

    @app.route("/api/regenerate", methods=["POST"])
    def api_regenerate():
        body = request.get_json(force=True) if request.data else {}
        view_name = body.get("view")

        # Build GeneratorData from DB for in-process generation
        gd = build_generator_data_from_db(db)

        views = get_views(db)
        if view_name:
            view = next((v for v in views if v["name"] == view_name), None)
            if not view:
                return jsonify({"error": "unknown view"}), 404
            ok = generate_svg_db(view_name, view["script"], gd, db_path=db)
            _broadcast("svg_updated", {"view": view_name})
            return jsonify({"ok": ok, "view": view_name})
        else:
            results = {}
            seen = set()
            for v in views:
                if v["script"] not in seen:
                    seen.add(v["script"])
                    ok = generate_svg_db(v["name"], v["script"], gd, db_path=db)
                    results[v["name"]] = ok
                    _broadcast("svg_updated", {"view": v["name"]})
                else:
                    results[v["name"]] = True
            return jsonify({"ok": True, "results": results})

    # -- Formula endpoints (Phase 12b) --

    @app.route("/api/formulas")
    def api_formulas():
        variant = request.args.get("variant")
        rows = get_all_formulas(variant=variant, db_path=db)
        return jsonify([dict(r) for r in rows])

    @app.route("/api/formulas/<element_name>")
    def api_element_formulas(element_name):
        variant = request.args.get("variant")
        rows = get_element_formulas(element_name, variant=variant, db_path=db)
        return jsonify([dict(r) for r in rows])

    @app.route("/api/formulas/<element_name>/<param_name>", methods=["PUT"])
    def api_upsert_formula(element_name, param_name):
        body = request.get_json(force=True)
        formula_json = body.get("formula")
        if not formula_json:
            return jsonify({"error": "missing formula"}), 400
        variant = body.get("variant")
        try:
            result = upsert_formula(element_name, param_name, formula_json,
                                    variant=variant, db_path=db)
        except Exception as exc:
            return jsonify({"error": str(exc)}), 400
        # Rebuild dependency cache
        from app.evaluator import extract_deps
        deps = extract_deps(formula_json)
        rebuild_formula_deps(element_name, param_name, deps, db_path=db)
        _broadcast("element_changed")
        _invalidate()
        return jsonify(result)

    @app.route("/api/formulas/<element_name>/element", methods=["DELETE"])
    def api_delete_formula_element(element_name):
        """Delete an element by removing its formula and re-basing dependents.

        Dependents' formulas are structurally rebased: references to the
        deleted element are replaced with the element's own formula
        sub-expressions (so dependents reference the same parent elements).
        Falls back to literal coordinate inlining when structural rebasing
        is not possible.  Supports undo.
        """
        from app.evaluator import rebase_element_refs, extract_deps
        variant = request.args.get("variant", "standard")

        # 1. Get the element's own formula (for structural rebasing)
        own_formulas = get_element_formulas(element_name, variant=variant,
                                            db_path=db)
        if not own_formulas:
            return jsonify({"error": "no formula found"}), 404

        # Parse the element's position formula for corner extraction
        elem_formula = None
        for f in own_formulas:
            if f["param_name"] == "position":
                fj = f["formula_json"]
                elem_formula = json.loads(fj) if isinstance(fj, str) else fj
                break

        # 2. Get evaluated geometry as fallback for literal inlining
        elem_data = _evaluate_element(element_name, variant)

        # 3. Find dependent formulas
        dependents = db_get_dependents(element_name, db_path=db)
        dep_names = {d["element_name"] for d in dependents}

        # 4. Capture before state for undo
        before_state = {
            "element_name": element_name,
            "own_formulas": [
                {"element_name": f["element_name"],
                 "param_name": f["param_name"],
                 "formula_json": f["formula_json"],
                 "locked": f.get("locked", 0),
                 "locked_value": f.get("locked_value"),
                 "variant": f.get("variant")}
                for f in own_formulas
            ],
            "updated_deps": {},
        }

        # 5. Rebase element references in dependent formulas.
        #    If structural rebasing produces a broken formula (collapsed to
        #    a non-dict, lost its "type", or still references the deleted
        #    element), fall back to a four_corner literal formula that
        #    preserves the dependent's current evaluated position.
        rebased_deps = {}  # key → new_json (for redo)
        for dep_name in dep_names:
            dep_formulas = get_element_formulas(dep_name, variant=variant,
                                                db_path=db)
            for f in dep_formulas:
                old_json = f["formula_json"]
                formula = json.loads(old_json) if isinstance(old_json, str) else old_json
                new_formula = rebase_element_refs(
                    formula, element_name, elem_formula, elem_data)
                # Check if rebasing produced a broken formula
                broken = (not isinstance(new_formula, dict)
                          or "type" not in new_formula)
                if not broken:
                    new_deps = extract_deps(new_formula)
                    broken = any(dt == "element" and dn == element_name
                                 for dt, dn in new_deps)
                if broken:
                    # Fall back: convert to four_corner with literal coords
                    dep_data = _evaluate_element(dep_name, variant)
                    poly = dep_data.get("poly", []) if dep_data else []
                    if len(poly) >= 4:
                        new_formula = {
                            "type": "four_corner",
                            "sw": [round(poly[0][0], 10),
                                   round(poly[0][1], 10)],
                            "se": [round(poly[1][0], 10),
                                   round(poly[1][1], 10)],
                            "ne": [round(poly[2][0], 10),
                                   round(poly[2][1], 10)],
                            "nw": [round(poly[3][0], 10),
                                   round(poly[3][1], 10)],
                        }
                    else:
                        continue  # can't rebase, leave unchanged
                new_json = json.dumps(new_formula)
                if new_json != old_json:
                    key = f"{f['element_name']}/{f['param_name']}"
                    before_state["updated_deps"][key] = old_json
                    rebased_deps[key] = new_json
                    upsert_formula(f["element_name"], f["param_name"],
                                   new_formula, variant=f.get("variant"),
                                   db_path=db)
                    new_deps = extract_deps(new_formula)
                    rebuild_formula_deps(f["element_name"], f["param_name"],
                                         list(new_deps), db_path=db)

        # 6. Delete the element's own formula(s)
        for f in own_formulas:
            delete_formula(f["element_name"], f["param_name"],
                           variant=f.get("variant"), db_path=db)
            rebuild_formula_deps(f["element_name"], f["param_name"],
                                  [], db_path=db)

        # 7. Also delete the elements DB record if one exists
        deleted_element = None
        elem_rec = get_element_by_name(element_name, db)
        if elem_rec:
            deleted_element = dict(elem_rec)
            delete_element(elem_rec["id"], db)
            before_state["deleted_element"] = deleted_element

        # 8. Record undo
        after_state = {"element_name": element_name,
                       "variant": variant,
                       "rebased_deps": rebased_deps}
        undo_mgr.record("formula_delete_element", before_state,
                         after_state, f"Delete {element_name}")

        _broadcast("element_changed")
        _invalidate()
        return jsonify({"ok": True, "rebased": list(dep_names)})

    @app.route("/api/formulas/<element_name>/<param_name>", methods=["DELETE"])
    def api_delete_formula(element_name, param_name):
        variant = request.args.get("variant")
        ok = delete_formula(element_name, param_name, variant=variant,
                            db_path=db)
        if not ok:
            return jsonify({"error": "not found"}), 404
        # Clear dependency cache
        rebuild_formula_deps(element_name, param_name, [], db_path=db)
        _broadcast("element_changed")
        _invalidate()
        return jsonify({"ok": True})

    @app.route("/api/formulas/<element_name>/<param_name>/lock", methods=["POST"])
    def api_lock_formula(element_name, param_name):
        body = request.get_json(force=True)
        locked = body.get("locked", True)
        locked_value = body.get("locked_value")
        variant = body.get("variant")
        # Capture before state for undo
        formulas = get_element_formulas(element_name, db_path=db)
        old_formula = None
        for f in formulas:
            if f["param_name"] == param_name:
                if variant is None and f.get("variant") is None:
                    old_formula = f
                    break
                elif variant is not None and f.get("variant") == variant:
                    old_formula = f
                    break
        ok = set_formula_lock(element_name, param_name, locked,
                              locked_value=locked_value, variant=variant,
                              db_path=db)
        if not ok:
            return jsonify({"error": "formula not found"}), 404
        # Record undo
        old_lv = None
        if old_formula and old_formula.get("locked_value"):
            try:
                old_lv = json.loads(old_formula["locked_value"])
            except (TypeError, ValueError):
                old_lv = old_formula["locked_value"]
        before = {"element_name": element_name, "param_name": param_name,
                  "locked": bool(old_formula["locked"]) if old_formula else False,
                  "locked_value": old_lv, "variant": variant}
        after = {"element_name": element_name, "param_name": param_name,
                 "locked": locked, "locked_value": locked_value,
                 "variant": variant}
        desc = f"{'Lock' if locked else 'Unlock'} {element_name}.{param_name}"
        undo_mgr.record("formula_lock", before, after, desc)
        _broadcast("formula_locked")
        _broadcast("element_changed")
        _invalidate()
        return jsonify({"ok": True})

    @app.route("/api/formulas/<element_name>/deps")
    def api_formula_deps(element_name):
        param_name = request.args.get("param")
        rows = get_formula_deps(element_name, param_name=param_name,
                                db_path=db)
        return jsonify([dict(r) for r in rows])

    @app.route("/api/formulas/<element_name>/dependents")
    def api_formula_dependents(element_name):
        dep_type = request.args.get("type", "element")
        rows = db_get_dependents(element_name, dep_type=dep_type, db_path=db)
        return jsonify([dict(r) for r in rows])

    @app.route("/api/deps/graph")
    def api_deps_graph():
        """Return the full formula dependency DAG as nodes and edges."""
        all_formulas = get_all_formulas(db_path=db)
        all_deps = get_all_formula_deps(db_path=db)
        nodes = {}
        for f in all_formulas:
            name = f["element_name"]
            if name not in nodes:
                nodes[name] = {"name": name, "locked": False, "params": []}
            nodes[name]["params"].append(f["param_name"])
            if f.get("locked"):
                nodes[name]["locked"] = True
        edges = []
        seen = set()
        for d in all_deps:
            if d["dep_type"] == "element":
                edge = (d["dep_name"], d["element_name"])
                if edge not in seen:
                    edges.append({"from": d["dep_name"],
                                  "to": d["element_name"]})
                    seen.add(edge)
        return jsonify({"nodes": list(nodes.values()), "edges": edges})

    # -- SSE endpoint --

    @app.route("/api/events")
    def api_events():
        def stream():
            q = queue.Queue(maxsize=50)
            with _sse_lock:
                _sse_clients.append(q)
            try:
                yield "event: connected\ndata: {}\n\n"
                while True:
                    try:
                        msg = q.get(timeout=30)
                        yield msg
                    except queue.Empty:
                        yield ": keepalive\n\n"
            except GeneratorExit:
                pass
            finally:
                with _sse_lock:
                    if q in _sse_clients:
                        _sse_clients.remove(q)

        return Response(stream(), mimetype="text/event-stream",
                        headers={"Cache-Control": "no-cache", "X-Accel-Buffering": "no"})

    # -- Version endpoint --

    @app.route("/api/version")
    def api_version():
        return jsonify({
            "git": _GIT_DESCRIBE,
            "started": _START_TIME,
        })

    return app
