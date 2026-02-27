"""Shared type definitions for the ADU project."""
from typing import Literal, NamedTuple

Point = tuple[float, float]

class LineSeg(NamedTuple):
    start: str; end: str

class ArcSeg(NamedTuple):
    start: str; end: str; center: str
    radius: float; direction: Literal["CW", "CCW"]
    n_pts: int

Segment = LineSeg | ArcSeg

class BBox(NamedTuple):
    """Axis-aligned bounding box: west, south, east, north edges."""
    w: float; s: float; e: float; n: float

class Wall(NamedTuple):
    """Rectangular element: polygon corners + axis-aligned bounding box."""
    poly: list[Point]  # [SW, SE, NE, NW] corners
    w: float           # min easting  (= BBox west)
    s: float           # min northing (= BBox south)
    e: float           # max easting  (= BBox east)
    n: float           # max northing (= BBox north)
