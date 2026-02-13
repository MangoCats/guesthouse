"""Shared type definitions for the hut2 project."""
from typing import Literal, NamedTuple

Point = tuple[float, float]

class LineSeg(NamedTuple):
    start: str; end: str

class ArcSeg(NamedTuple):
    start: str; end: str; center: str
    radius: float; direction: Literal["CW", "CCW"]
    n_pts: int

Segment = LineSeg | ArcSeg
