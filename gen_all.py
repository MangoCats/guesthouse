"""Regenerate all SVG files from the database with a consistent git describe stamp.

The database (app/adu.db) is the authoritative source for all building geometry.
Scripts with in-process handlers are run in-process using DB-sourced GeneratorData.
Scripts without in-process handlers (survey, scad/gen_views, scad/gen_line_drawings,
gen_3views, plumbing) fall back to subprocess — acceptable because those derive
geometry from survey traversal / 3D model parameters rather than DB formulas.

Captures `git describe --always --dirty=-DEV` once before any SVG is written,
so that all title blocks show the same version string even though writing the
first SVG makes the working tree dirty.
"""
import os, sys, subprocess
import fitz  # PyMuPDF

_DIR = os.path.dirname(os.path.abspath(__file__))
_DB  = os.path.join(_DIR, "app", "adu.db")
_CACHE = os.path.join(_DIR, ".git_describe")

# Scripts that have in-process handlers in engine._run_generator_inprocess().
# Stored as relative paths (same form engine.py uses for dispatch).
_INPROCESS_SCRIPTS = [
    "floorplan/gen_floorplan.py",
    "walls/gen_walls.py",
    "span/gen_span.py",
    "span/gen_span_minmax.py",
    "span/gen_span_min.py",
    "roof/gen_roof.py",
    "site/gen_site_plan.py",
    "scad/gen_flat_roof.py",
    "scad/gen_2in12.py",
    "survey/gen_path_svg.py",
]

# Scripts without in-process handlers — run as subprocesses.
# These use hardcoded procedural geometry (survey traversal, SCAD dimensions)
# which is not stored in the DB.
_SUBPROCESS_SCRIPTS = [
    os.path.join(_DIR, "survey", "gen_path_svg_wo.py"),
    os.path.join(_DIR, "survey", "gen_path_svg_ks.py"),
    os.path.join(_DIR, "survey", "gen_points_pdf.py"),
    os.path.join(_DIR, "scad", "gen_views.py"),
    os.path.join(_DIR, "scad", "gen_line_drawings.py"),
    os.path.join(_DIR, "gen_3views.py"),
    os.path.join(_DIR, "plumbing", "gen_plumbing.py"),
    os.path.join(_DIR, "scad", "gen_gltf.py"),
]


def main():
    # Ensure app/ is importable from the project root
    if _DIR not in sys.path:
        sys.path.insert(0, _DIR)

    # 1. Capture git describe before any SVG is written
    desc = subprocess.check_output(
        ["git", "describe", "--always", "--dirty=-DEV"],
        cwd=_DIR, text=True,
    ).strip()
    with open(_CACHE, "w") as f:
        f.write(desc)
    print(f"git describe: {desc}")

    try:
        # 2a. Build GeneratorData from DB (DB is authoritative)
        from app.engine import build_generator_data_from_db, generate_svg_db
        gd = build_generator_data_from_db(_DB)

        # 2b. Run in-process generators (DB-driven)
        for rel_path in _INPROCESS_SCRIPTS:
            print(f"  running {rel_path} ...")
            ok = generate_svg_db(rel_path, rel_path, gd, db_path=_DB)
            if not ok:
                print(f"    WARNING: {rel_path} failed")

        # 2c. Run subprocess generators (hardcoded geometry, no DB formulas)
        for script in _SUBPROCESS_SCRIPTS:
            rel = os.path.relpath(script, _DIR)
            print(f"  running {rel} ...")
            subprocess.check_call(["python", script], cwd=_DIR)

    finally:
        # 3. Always clean up the cache file
        if os.path.exists(_CACHE):
            os.remove(_CACHE)

    # 4. Render PDFs → PNGs at 1920px wide
    pdf_pngs = [
        os.path.join(_DIR, "site", "site_plan_df.pdf"),
        os.path.join(_DIR, "site", "site_plan_fs.pdf"),
        os.path.join(_DIR, "3views.pdf"),
    ]
    for pdf_path in pdf_pngs:
        png_path = pdf_path.replace(".pdf", ".png")
        doc = fitz.open(pdf_path)
        page = doc[0]
        scale = 1920 / page.rect.width
        pix = page.get_pixmap(matrix=fitz.Matrix(scale, scale))
        pix.save(png_path)
        doc.close()
        print(f"  rendered {os.path.relpath(png_path, _DIR)} ({pix.width}x{pix.height})")

    print("done.")


if __name__ == "__main__":
    main()
