"""Generate plumbing plan SVG.

Based on the daybed floorplan variant, showing only plumbing-relevant
fixtures: washer, toilets, sinks, water heater, fridge, and dishwasher.
"""
import os, sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

from floorplan.gen_floorplan import build_floorplan_data, render_floorplan_svg

if __name__ == "__main__":
    data = build_floorplan_data()
    svg = render_floorplan_svg(data, room_title="Plumbing Plan",
                               db=True, plumbing=True)
    base_dir = os.path.dirname(os.path.abspath(__file__))
    path = os.path.join(base_dir, "plumbing.svg")
    with open(path, "w") as f:
        f.write(svg)
    print(f"Plumbing plan written to {path}")
