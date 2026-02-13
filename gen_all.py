"""Regenerate all SVG files with a consistent git describe stamp.

Captures `git describe --always --dirty=-DEV` once before any SVG is
written, so that all title blocks show the same version string even
though writing the first SVG makes the working tree dirty.
"""
import os, subprocess

_DIR = os.path.dirname(os.path.abspath(__file__))
_CACHE = os.path.join(_DIR, ".git_describe")

_SCRIPTS = [
    os.path.join(_DIR, "survey", "gen_path_svg.py"),
    os.path.join(_DIR, "floorplan", "gen_floorplan.py"),
    os.path.join(_DIR, "walls", "gen_walls.py"),
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

    print("done.")


if __name__ == "__main__":
    main()
