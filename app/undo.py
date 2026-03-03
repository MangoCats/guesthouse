"""Undo/Redo manager for the ADU Editor.

Command-pattern undo system with 50-level depth.  Maintains an in-memory
stack with a position pointer and persists entries to the undo_history table.
"""
import json
import time

from app.database import get_db, update_constants_batch


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

        For Phase 2, all action types (constant_update, constant_batch,
        constant_reset) store name→value dicts and apply via batch update.
        Future phases extend this with additional action type dispatch.
        """
        if action_type in ("constant_update", "constant_batch", "constant_reset"):
            update_constants_batch(state, self._db_path)
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
