"""Regenerate all SVG files with a consistent git describe stamp.

Captures `git describe --always --dirty=-DEV` once before any SVG is
written, so that all title blocks show the same version string even
though writing the first SVG makes the working tree dirty.
"""
import os, subprocess
import fitz  # PyMuPDF

_DIR = os.path.dirname(os.path.abspath(__file__))
_CACHE = os.path.join(_DIR, ".git_describe")

_SCRIPTS = [
    os.path.join(_DIR, "survey", "gen_path_svg.py"),
    os.path.join(_DIR, "floorplan", "gen_floorplan.py"),
    os.path.join(_DIR, "walls", "gen_walls.py"),
    os.path.join(_DIR, "span", "gen_span.py"),
    os.path.join(_DIR, "span", "gen_span_minmax.py"),
    os.path.join(_DIR, "span", "gen_span_min.py"),
    os.path.join(_DIR, "roof", "gen_roof.py"),
    os.path.join(_DIR, "site", "gen_site_plan.py"),
    os.path.join(_DIR, "scad", "gen_flat_roof.py"),
    os.path.join(_DIR, "scad", "gen_2in12.py"),
    os.path.join(_DIR, "scad", "gen_views.py"),
]


def main():
    # 1. Capture git describe before any SVG is written
    desc = subprocess.check_output(
        ["git", "describe", "--always", "--dirty=-DEV"],
        cwd=_DIR, text=True,
    ).strip()
    with open(_CACHE, "w") as f:
        f.write(desc)
    print(f"git describe: {desc}")

    try:
        # 2. Run each generator
        for script in _SCRIPTS:
            print(f"  running {os.path.relpath(script, _DIR)} ...")
            subprocess.check_call(["python", script], cwd=_DIR)
    finally:
        # 3. Always clean up the cache file
        if os.path.exists(_CACHE):
            os.remove(_CACHE)

    # 4. Render site_plan_df.pdf → site_plan_df.png at 1920px wide
    pdf_path = os.path.join(_DIR, "site", "site_plan_df.pdf")
    png_path = os.path.join(_DIR, "site", "site_plan_df.png")
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
