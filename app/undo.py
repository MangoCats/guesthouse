"""Undo/Redo manager for the ADU Editor.

Command-pattern undo system with 50-level depth.  Maintains an in-memory
stack with a position pointer and persists entries to the undo_history table.
"""
import json
import time

from app.database import (
    get_db, update_constant, update_constants_batch,
    create_element_raw, update_element, delete_element,
    get_element, get_element_by_name,
    create_door_raw, update_door, delete_door,
    set_formula_lock, get_element_formulas, upsert_formula,
)


class UndoManager:
    """Manages undo/redo history for all mutation operations."""

    def __init__(self, db_path, max_depth=50):
        self._db_path = db_path
        self._max_depth = max_depth
        self._stack = []      # list of entry dicts
        self._position = -1   # index of last applied action
        self._load()

    def record(self, action_type, before_state, after_state, description=""):
        """Record a new action.  Discards any redo entries past position."""
        # Truncate redo entries
        self._stack = self._stack[:self._position + 1]
        entry = {
            "action_type": action_type,
            "before_state": before_state,
            "after_state": after_state,
            "timestamp": time.time(),
            "description": description,
        }
        self._stack.append(entry)
        self._position = len(self._stack) - 1
        # Enforce depth limit
        if len(self._stack) > self._max_depth:
            excess = len(self._stack) - self._max_depth
            self._stack = self._stack[excess:]
            self._position = len(self._stack) - 1
        self._persist()

    def undo(self):
        """Undo the action at current position.  Returns entry dict or None."""
        if not self.can_undo:
            return None
        entry = self._stack[self._position]
        self._apply(entry["action_type"], entry["before_state"])
        self._position -= 1
        return entry

    def redo(self):
        """Redo the next action.  Returns entry dict or None."""
        if not self.can_redo:
            return None
        self._position += 1
        entry = self._stack[self._position]
        self._apply(entry["action_type"], entry["after_state"])
        return entry

    @property
    def can_undo(self):
        """True if there is an action to undo."""
        return self._position >= 0

    @property
    def can_redo(self):
        """True if there is an action to redo."""
        return self._position < len(self._stack) - 1

    def _apply(self, action_type, state):
        """Apply a state snapshot to the database.

        Dispatches by action_type:
        - constant_*: batch-update constants from name→value dict
        - element_create undo: delete the created element (state = {"id": N})
        - element_delete undo: re-insert full record (state = full element dict)
        - element_update: restore old field values (state = {"id": N, ...fields})
        - door_create undo: delete the created door (state = {"opening_name": X})
        - door_delete undo: re-insert full record (state = full door dict)
        - door_update: restore old field values (state = {"opening_name": X, ...})
        - opening_*: same pattern as element_* (openings are elements)
        """
        if action_type in ("constant_update", "constant_batch", "constant_reset"):
            update_constants_batch(state, self._db_path)
        elif action_type == "element_create":
            # Undo: delete the element that was created
            delete_element(state["id"], self._db_path)
        elif action_type == "element_delete":
            # Undo: re-insert the deleted element(s)
            recs = state if isinstance(state, list) else [state]
            for rec in recs:
                create_element_raw(rec, self._db_path)
                # Re-create formula for placed items (absolute fallback)
                props = rec.get("properties", {})
                if isinstance(props, str):
                    try:
                        props = json.loads(props)
                    except (ValueError, TypeError):
                        props = {}
                if props.get("source") == "placed":
                    center = props.get("center", [0, 0])
                    shape = props.get("shape", "rect")
                    w = props.get("width", 1)
                    d = props.get("depth", 1)
                    rotation = props.get("rotation", 0)
                    if shape == "circle":
                        r = props.get("radius") or w / 2
                        formula = {"type": "item_circle",
                                   "center": center, "radius": r}
                    else:
                        formula = {"type": "item_rect",
                                   "anchor": center,
                                   "along": [1, 0], "across": [0, 1],
                                   "width": w, "depth": d,
                                   "anchor_corner": "center",
                                   "rotation_deg": rotation}
                    upsert_formula(rec["name"], "position", formula,
                                   variant=None, db_path=self._db_path)
        elif action_type == "element_update":
            # Undo: restore old field values
            eid = state["id"]
            fields = {k: v for k, v in state.items() if k != "id"}
            update_element(eid, fields, self._db_path)
        elif action_type == "door_create":
            # Undo: delete the door that was created
            delete_door(state["opening_name"], self._db_path)
        elif action_type == "door_delete":
            # Undo: re-insert the deleted door
            create_door_raw(state, self._db_path)
        elif action_type == "door_update":
            # Undo: restore old field values
            oname = state["opening_name"]
            fields = {k: v for k, v in state.items() if k != "opening_name" and k != "id"}
            update_door(oname, fields, self._db_path)
        elif action_type == "opening_create":
            delete_element(state["id"], self._db_path)
        elif action_type == "opening_delete":
            if isinstance(state, list):
                for rec in state:
                    create_element_raw(rec, self._db_path)
            else:
                create_element_raw(state, self._db_path)
        elif action_type == "opening_update":
            eid = state["id"]
            fields = {k: v for k, v in state.items() if k != "id"}
            update_element(eid, fields, self._db_path)
        elif action_type == "element_move":
            # Move undo/redo: state contains move_type + details
            if state.get("move_type") == "constant":
                update_constant(state["constant"], state["value"], self._db_path)
            elif state.get("move_type") == "formula":
                # Unified formula-based move: restore formula + properties
                if state.get("formula"):
                    upsert_formula(state["name"], "position",
                                   state["formula"], variant=None,
                                   db_path=self._db_path)
                if state.get("properties"):
                    update_element(state["id"],
                                   {"properties": state["properties"]},
                                   self._db_path)
            elif state.get("move_type") == "placed":
                # Legacy placed item move (for undo history compat)
                update_element(state["id"],
                               {"properties": state["properties"]}, self._db_path)
                props = state["properties"]
                el = get_element(state["id"], self._db_path)
                if el:
                    center = props.get("center", [0, 0])
                    shape = props.get("shape", "rect")
                    if shape == "circle":
                        r = props.get("radius") or props.get("width", 1) / 2
                        formula = {"type": "item_circle",
                                   "center": center, "radius": r}
                    else:
                        w = props.get("width", 1)
                        d = props.get("depth", 1)
                        sw = [center[0] - w / 2, center[1] - d / 2]
                        formula = {"type": "item_rect", "anchor": sw,
                                   "along": [1, 0], "across": [0, 1],
                                   "width": w, "depth": d,
                                   "anchor_corner": "sw"}
                    upsert_formula(el["name"], "position", formula,
                                   variant=None, db_path=self._db_path)
            else:
                update_element(state["id"],
                               {"properties": state["properties"]}, self._db_path)
        elif action_type == "formula_update":
            # Formula field update (width/depth/rotation): restore formula + props
            if state.get("formula"):
                upsert_formula(state["name"], "position",
                               state["formula"], variant=None,
                               db_path=self._db_path)
            if state.get("properties"):
                update_element(state["id"],
                               {"properties": state["properties"]},
                               self._db_path)
        elif action_type == "plumbing_create":
            from app.plumbing import delete_plumbing_element
            delete_plumbing_element(state["id"], self._db_path)
        elif action_type == "plumbing_delete":
            from app.plumbing import create_plumbing_raw
            create_plumbing_raw(state, self._db_path)
        elif action_type == "plumbing_update":
            from app.plumbing import update_plumbing_element
            eid = state["id"]
            fields = {k: v for k, v in state.items() if k != "id"}
            update_plumbing_element(eid, fields, self._db_path)
        elif action_type in ("outline_update", "outline_add_point",
                             "outline_remove_point"):
            # Outline undo/redo: state is full chain snapshot (list of dicts)
            from app.database import restore_outline_chain
            restore_outline_chain(state, self._db_path)
        elif action_type == "outline_parcel_reset":
            # Parcel reset undo/redo: restore chain + anchor position
            from app.database import (restore_outline_chain,
                                      set_outline_anchor_pivot,
                                      clear_outline_pivot)
            restore_outline_chain(state["chain"], self._db_path)
            apos = state.get("anchor_pos")
            if state.get("anchor_name") and state.get("pivot_name") and apos:
                set_outline_anchor_pivot(
                    state["anchor_name"], state["pivot_name"],
                    apos[0], apos[1], apos[2],
                    self._db_path, user_set=True,
                )
            else:
                clear_outline_pivot(self._db_path)
        elif action_type == "outline_pivot":
            # Pivot undo/redo: state = {"chain": [...], "anchor": str|None,
            #                           "pivot": str|None}
            from app.database import (restore_outline_chain,
                                      set_outline_anchor_pivot,
                                      clear_outline_pivot)
            restore_outline_chain(state["chain"], self._db_path)
            if state.get("anchor") and state.get("pivot"):
                set_outline_anchor_pivot(state["anchor"], state["pivot"],
                                         self._db_path)
            else:
                # Clear config keys (no pivot active)
                with get_db(self._db_path) as conn:
                    conn.execute("DELETE FROM config WHERE key IN "
                                 "('outline_anchor', 'outline_pivot')")
        elif action_type == "variant_update":
            from app.database import update_variant
            vid = state["id"]
            fields = {k: v for k, v in state.items() if k != "id"}
            update_variant(vid, fields, self._db_path)
        elif action_type == "variant_create":
            from app.database import (get_variant_by_id, create_variant_raw,
                                      clone_variant_exclusions,
                                      clone_variant_elements,
                                      unclone_variant_elements,
                                      delete_variant_exclusions,
                                      delete_variant)
            existing = get_variant_by_id(state["id"], self._db_path)
            if existing:
                # Undo: remove the variant and its cloned elements/exclusions
                unclone_variant_elements(state["name"], self._db_path)
                delete_variant_exclusions(state["name"], self._db_path)
                delete_variant(state["id"], self._db_path)
            else:
                # Redo: re-create the variant
                create_variant_raw(state, self._db_path)
                source = state.get("source_variant", "standard")
                clone_variant_exclusions(source, state["name"], self._db_path)
                clone_variant_elements(source, state["name"], self._db_path)
        elif action_type == "variant_delete":
            from app.database import (get_variant_by_id, create_variant_raw,
                                      clone_variant_exclusions,
                                      clone_variant_elements,
                                      unclone_variant_elements,
                                      delete_variant_exclusions,
                                      delete_variant)
            existing = get_variant_by_id(state["id"], self._db_path)
            if existing:
                # Redo: re-delete the variant
                unclone_variant_elements(state["name"], self._db_path)
                delete_variant_exclusions(state["name"], self._db_path)
                delete_variant(state["id"], self._db_path)
            else:
                # Undo: re-create the variant and restore cloned data
                create_variant_raw(state, self._db_path)
                source = state.get("source_variant", "standard")
                clone_variant_exclusions(source, state["name"], self._db_path)
                clone_variant_elements(source, state["name"], self._db_path)
        elif action_type == "formula_lock":
            set_formula_lock(
                state["element_name"], state["param_name"],
                state["locked"], locked_value=state.get("locked_value"),
                variant=state.get("variant"), db_path=self._db_path,
            )
        elif action_type == "formula_delete_element":
            self._apply_formula_delete(state)
        elif action_type in ("survey_leg_update", "survey_config_update",
                             "survey_reset"):
            from app.database import restore_survey
            restore_survey(state.get("legs", []), state.get("config", {}),
                           self._db_path)
        elif action_type in ("inner_wall_override_upsert",
                             "inner_wall_override_delete"):
            from app.database import restore_inner_wall_overrides
            restore_inner_wall_overrides(state, self._db_path)
        elif action_type == "project_import":
            from app.database import import_project
            import_project(state, self._db_path)
        elif action_type == "full_reset":
            # Combined constants + outline + elements + doors + survey + overrides
            from app.database import (restore_outline_chain, restore_elements,
                                      restore_survey,
                                      restore_inner_wall_overrides)
            update_constants_batch(state["constants"], self._db_path)
            restore_outline_chain(state["outline"], self._db_path)
            if "elements" in state:
                restore_elements(state["elements"], state.get("doors", []),
                                 self._db_path)
            if "survey_legs" in state:
                restore_survey(state["survey_legs"],
                               state.get("survey_config", {}), self._db_path)
            if "inner_wall_overrides" in state:
                restore_inner_wall_overrides(state["inner_wall_overrides"],
                                             self._db_path)
        else:
            raise ValueError(f"Unknown undo action type: {action_type}")

    def _apply_formula_delete(self, state):
        """Apply undo/redo for formula_delete_element action."""
        from app.database import (upsert_formula, delete_formula,
                                  rebuild_formula_deps)
        from app.evaluator import extract_deps
        if "own_formulas" in state:
            # Undo: re-insert the deleted formula(s)
            for f in state["own_formulas"]:
                formula = json.loads(f["formula_json"]) if isinstance(
                    f["formula_json"], str) else f["formula_json"]
                upsert_formula(f["element_name"], f["param_name"],
                               formula, variant=f.get("variant"),
                               db_path=self._db_path)
                deps = extract_deps(formula)
                rebuild_formula_deps(f["element_name"], f["param_name"],
                                     list(deps), db_path=self._db_path)
                if f.get("locked"):
                    set_formula_lock(f["element_name"], f["param_name"],
                                     True, locked_value=f.get("locked_value"),
                                     variant=f.get("variant"),
                                     db_path=self._db_path)
            # Restore dependent formulas to their pre-inlined state
            for key, old_json in state.get("updated_deps", {}).items():
                parts = key.split("/", 1)
                if len(parts) == 2:
                    formula = json.loads(old_json) if isinstance(
                        old_json, str) else old_json
                    upsert_formula(parts[0], parts[1], formula,
                                   db_path=self._db_path)
                    deps = extract_deps(formula)
                    rebuild_formula_deps(parts[0], parts[1],
                                          list(deps), db_path=self._db_path)
            # Re-insert the elements record(s) if deleted
            if "deleted_elements" in state:
                for de in state["deleted_elements"]:
                    create_element_raw(de, self._db_path)
            elif "deleted_element" in state:
                create_element_raw(state["deleted_element"], self._db_path)
        else:
            # Redo: re-inline dependents then re-delete
            elem_name = state["element_name"]
            redo_variant = state.get("variant")
            # Re-apply rebased formulas to dependents
            for key, new_json in state.get("rebased_deps",
                                           state.get("inlined_deps", {})).items():
                parts = key.split("/", 1)
                if len(parts) == 2:
                    formula = json.loads(new_json) if isinstance(
                        new_json, str) else new_json
                    upsert_formula(parts[0], parts[1], formula,
                                   db_path=self._db_path)
                    deps = extract_deps(formula)
                    rebuild_formula_deps(parts[0], parts[1],
                                          list(deps), db_path=self._db_path)
            # Delete the element's own formulas for this variant
            formulas = get_element_formulas(elem_name,
                                            variant=redo_variant,
                                            db_path=self._db_path)
            for f in formulas:
                delete_formula(f["element_name"], f["param_name"],
                               variant=f.get("variant"),
                               db_path=self._db_path)
                rebuild_formula_deps(f["element_name"], f["param_name"],
                                     [], db_path=self._db_path)
            # Delete element record only if no formulas remain
            with get_db(self._db_path) as conn:
                remaining = conn.execute(
                    "SELECT 1 FROM element_formulas "
                    "WHERE element_name = ? LIMIT 1", (elem_name,),
                ).fetchone()
            if not remaining:
                elem_rec = get_element_by_name(elem_name,
                                               self._db_path)
                if elem_rec:
                    delete_element(elem_rec["id"], self._db_path)

    def _load(self):
        """Rebuild the in-memory stack from the undo_history table."""
        with get_db(self._db_path) as conn:
            rows = conn.execute(
                "SELECT action_type, before_state, after_state, timestamp, description "
                "FROM undo_history ORDER BY id"
            ).fetchall()
        self._stack = []
        for r in rows:
            self._stack.append({
                "action_type": r["action_type"],
                "before_state": json.loads(r["before_state"]),
                "after_state": json.loads(r["after_state"]),
                "timestamp": r["timestamp"],
                "description": r["description"],
            })
        self._position = len(self._stack) - 1

    def _persist(self):
        """Sync the in-memory stack to the undo_history table."""
        with get_db(self._db_path) as conn:
            conn.execute("DELETE FROM undo_history")
            for entry in self._stack:
                conn.execute(
                    "INSERT INTO undo_history "
                    "(timestamp, action_type, before_state, after_state, description) "
                    "VALUES (?, ?, ?, ?, ?)",
                    (
                        entry["timestamp"],
                        entry["action_type"],
                        json.dumps(entry["before_state"]),
                        json.dumps(entry["after_state"]),
                        entry.get("description", ""),
                    ),
                )
