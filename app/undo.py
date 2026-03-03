"""Undo/Redo manager for the ADU Editor.

Command-pattern undo system with 50-level depth.  Maintains an in-memory
stack with a position pointer and persists entries to the undo_history table.
"""
import json
import time

from app.database import (
    get_db, update_constants_batch,
    create_element_raw, update_element, delete_element,
    create_door_raw, update_door, delete_door,
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
            if isinstance(state, list):
                for rec in state:
                    create_element_raw(rec, self._db_path)
            else:
                create_element_raw(state, self._db_path)
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
        else:
            raise ValueError(f"Unknown undo action type: {action_type}")

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
