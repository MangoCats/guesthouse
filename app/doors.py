"""Door business logic for the ADU Editor.

Provides higher-level operations on the doors table and door-specific
validation.  Door geometry rendering is currently handled by the engine
(gen_floorplan.py); this module manages the database records.
"""
from app.database import get_all_doors, get_door, create_door, update_door, delete_door

VALID_SIDES = {"east", "west", "north", "south"}
VALID_TYPES = {"single", "double"}


def validate_door(hinge_side, swing_direction, door_type="single"):
    """Validate door parameters. Returns error string or None."""
    if hinge_side not in VALID_SIDES:
        return f"invalid hinge_side: {hinge_side}"
    if swing_direction not in VALID_SIDES:
        return f"invalid swing_direction: {swing_direction}"
    if door_type not in VALID_TYPES:
        return f"invalid door_type: {door_type}"
    return None
