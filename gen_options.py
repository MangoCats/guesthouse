"""Generate Daybed-variant floorplan SVGs for each config DB into options/."""
import os, sys

_DIR = os.path.dirname(os.path.abspath(__file__))
if _DIR not in sys.path:
    sys.path.insert(0, _DIR)

CONFIGS = [
    ("Base",  os.path.join(_DIR, "app", "configs", "Base.db")),
    ("Mark6", os.path.join(_DIR, "app", "configs", "Mark6.db")),
    ("Mark8", os.path.join(_DIR, "app", "configs", "Mark8.db")),
    ("Mark9", os.path.join(_DIR, "app", "configs", "Mark9.db")),
    ("MarkX", os.path.join(_DIR, "app", "configs", "MarkX.db")),
]

OUT_DIR = os.path.join(_DIR, "options")
os.makedirs(OUT_DIR, exist_ok=True)

from app.engine import build_generator_data_from_db, compute_geometry
from app.database import get_constants_dict, get_outline_chain, get_all_doors
from floorplan.gen_floorplan import build_floorplan_data
from plumbing.gen_plumbing import _compute_boundary_corners
from app.db_render import render_floorplan_svg_db

for name, db_path in CONFIGS:
    print(f"Generating {name} daybed SVG...")
    gd = build_generator_data_from_db(db_path)
    data = build_floorplan_data(gd)
    constants_dict = get_constants_dict(db_path)
    chain_rows = get_outline_chain(db_path)
    doors_data = get_all_doors(db_path)

    geom = compute_geometry(
        constants_dict, variant="daybed",
        chain_rows=chain_rows,
        doors_data=doors_data, db_path=db_path)

    svg = render_floorplan_svg_db(
        geom, data, room_title=f"Parent Suite with Daybed — {name}")

    out_path = os.path.join(OUT_DIR, f"{name}_daybed.svg")
    with open(out_path, "w", encoding="utf-8") as f:
        f.write(svg)
    print(f"  -> {out_path}")

print("Done.")
