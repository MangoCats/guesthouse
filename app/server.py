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
)
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
        })

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
