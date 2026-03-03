"""Flask web server for the ADU Editor application.

Provides REST API for constants/geometry/SVG, WebSocket for real-time
updates, and serves the single-page frontend.
"""
import json
import os
import queue
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
)
from app.doors import validate_door
from app.engine import compute_geometry, generate_svg, get_svg_content, patch_constants
from app.undo import UndoManager

_PROJECT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# SSE event bus — all connected clients get geometry update notifications
_sse_clients: list[queue.Queue] = []
_sse_lock = threading.Lock()


def _broadcast(event: str, data: dict | None = None):
    """Send an SSE event to all connected clients."""
    msg = f"event: {event}\ndata: {json.dumps(data or {})}\n\n"
    with _sse_lock:
        dead = []
        for q in _sse_clients:
            try:
                q.put_nowait(msg)
            except queue.Full:
                dead.append(q)
        for q in dead:
            _sse_clients.remove(q)


def create_app(db_path=None):
    """Create and configure the Flask application."""
    db = db_path or DB_PATH
    init_db(db)

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
                data = compute_geometry(constants, variant)
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
        before = get_constants_dict(db)
        reset_constants(db)
        after = get_constants_dict(db)
        undo_mgr.record("constant_reset", before, after, "Reset all constants")
        _invalidate()
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
        _broadcast_undo_status()
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
        _broadcast_undo_status()
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
        _broadcast_undo_status()
        return jsonify({"ok": True, "deleted": deleted_ids})

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
        _broadcast_undo_status()
        return jsonify(record), 201

    @app.route("/api/openings/<name>", methods=["PUT"])
    def api_update_opening(name):
        old = get_element_by_name(name, db)
        if not old or old["type"] != "opening":
            return jsonify({"error": "not found"}), 404
        body = request.get_json(force=True)
        # Merge updates into properties
        import json as _json
        props = _json.loads(old["properties"]) if isinstance(old["properties"], str) else old["properties"]
        for k in ("width", "offset", "segment", "host_wall"):
            if k in body:
                props[k] = body[k]
        updated = update_element(old["id"], {"properties": props}, db)
        undo_mgr.record(
            "opening_update", old, updated,
            f"Update opening {name}",
        )
        _broadcast("element_changed")
        _broadcast_undo_status()
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
        _broadcast_undo_status()
        return jsonify({"ok": True, "deleted": deleted_ids})

    # -- Doors API --

    @app.route("/api/doors")
    def api_doors():
        return jsonify(get_all_doors(db))

    @app.route("/api/doors", methods=["POST"])
    def api_create_door():
        body = request.get_json(force=True)
        opening = body.get("opening")
        if not opening:
            return jsonify({"error": "missing opening"}), 400
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
        _broadcast_undo_status()
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
        _broadcast_undo_status()
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
        _broadcast_undo_status()
        return jsonify({"ok": True, "opening_name": opening_name})

    # -- Geometry API --

    @app.route("/api/geometry")
    def api_geometry():
        try:
            from app.variants import VARIANTS
            variant = request.args.get("variant", "standard")
            if variant not in VARIANTS:
                variant = "standard"
            geom = _get_geometry(variant)
            return jsonify(geom)
        except Exception as exc:
            return jsonify({"error": str(exc)}), 500

    @app.route("/api/variants")
    def api_variants():
        from app.variants import VARIANTS
        return jsonify([{"name": k, "label": v["label"]} for k, v in VARIANTS.items()])

    # -- Outline chain API --

    @app.route("/api/outline")
    def api_outline():
        return jsonify(get_outline_chain(db))

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
        if view_name == "floorplan" and variant in _FLOORPLAN_VARIANT_SUFFIX:
            suffix = _FLOORPLAN_VARIANT_SUFFIX[variant]
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
        mime = "application/pdf" if full_path.endswith(".pdf") else "image/svg+xml"
        return send_file(full_path, mimetype=mime)

    # -- Regeneration API --

    @app.route("/api/regenerate", methods=["POST"])
    def api_regenerate():
        body = request.get_json(force=True) if request.data else {}
        view_name = body.get("view")

        # Apply constants to the module before regenerating
        constants = get_constants_dict(db)
        patch_constants(constants)

        views = get_views(db)
        if view_name:
            view = next((v for v in views if v["name"] == view_name), None)
            if not view:
                return jsonify({"error": "unknown view"}), 404
            ok = generate_svg(view_name, view["script"])
            _broadcast("svg_updated", {"view": view_name})
            return jsonify({"ok": ok, "view": view_name})
        else:
            results = {}
            seen = set()
            for v in views:
                if v["script"] not in seen:
                    seen.add(v["script"])
                    ok = generate_svg(v["name"], v["script"])
                    results[v["name"]] = ok
                    _broadcast("svg_updated", {"view": v["name"]})
                else:
                    results[v["name"]] = True
            return jsonify({"ok": True, "results": results})

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

    return app
