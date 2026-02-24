"""Update d² regression expected values in test files.

Usage:
    python tests/update_d2.py [--check]

Without flags: recompute all d² values and update test files in place.
With --check:  report which values changed, but don't modify files.

Handles three d² tables:
  1. test_d2_regression.py   EXPECTED       (layout points → F1, F6, F12, F15)
  2. test_gen_site_plan.py   EXPECTED_D2    (F-series PDF → parcel corners)
  3. test_gen_site_plan.py   EXPECTED_P_D2  (P-series PDF → parcel corners)
"""
import importlib.util
import math
import os
import re
import sys

_PROJECT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _PROJECT)
sys.path.insert(0, os.path.join(_PROJECT, "tests"))

from shared.survey import compute_traverse, compute_three_arc, compute_inset
from shared.geometry import compute_inner_walls, path_polygon
from floorplan.geometry import compute_outline_geometry, align_pts_to_f_series
from floorplan.layout import compute_interior_layout
from floorplan.constants import WALL_OUTER


def _dist_sq(a, b):
    return (a[0] - b[0])**2 + (a[1] - b[1])**2


# ============================================================
# 1. Layout d² (test_d2_regression.py)
# ============================================================

def _build_layout_fixtures():
    """Replicate the conftest fixture chain to get pts, layout, radii."""
    pts = compute_traverse()
    arc = compute_three_arc(pts)
    inset = compute_inset(pts, arc["R1"], arc["R2"], arc["R3"], arc["nE"], arc["nN"])
    pts.update(inset.pts_update)
    align_pts_to_f_series(pts)
    geo = compute_outline_geometry()
    pts.update(geo.fp_pts)
    inner_segs = compute_inner_walls(geo.outline_segs, pts, WALL_OUTER, geo.radii)
    inner_poly = path_polygon(inner_segs, pts)
    layout = compute_interior_layout(pts, inner_poly)
    return pts, layout, geo.radii


def _compute_layout_d2():
    """Compute current layout d² values using _collect_all_points from the test."""
    pts, layout, radii = _build_layout_fixtures()

    # Import _collect_all_points from the test module
    spec = importlib.util.spec_from_file_location(
        "test_d2_regression",
        os.path.join(_PROJECT, "tests", "test_d2_regression.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)

    all_pts = mod._collect_all_points(pts, layout, radii)
    ref_names = mod.REF_NAMES
    ref_pts = {n: pts[n] for n in ref_names}

    rows = []
    for name, pt in all_pts:
        d2s = tuple(_dist_sq(pt, ref_pts[rn]) for rn in ref_names)
        rows.append((name, *d2s))
    return rows


# ============================================================
# 2. Site plan d² (test_gen_site_plan.py)
# ============================================================

def _build_site_plan_data():
    """Import and call build_site_plan_data from gen_site_plan.py."""
    spec = importlib.util.spec_from_file_location(
        "gen_site_plan",
        os.path.join(_PROJECT, "site", "gen_site_plan.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.build_site_plan_data(), mod


def _compute_site_d2():
    """Compute current F-series site plan d² values."""
    sp, mod = _build_site_plan_data()
    ref_pts = {
        "CORNER_NW": mod.CORNER_NW,
        "CORNER_NE": mod.CORNER_NE,
        "CORNER_SE": mod.CORNER_SE,
        "CORNER_SW": mod.CORNER_SW,
    }
    ref_names = ["CORNER_NW", "CORNER_NE", "CORNER_SE", "CORNER_SW"]

    f_rows = []
    for name, pt in sp.f_series_pdf.items():
        d2s = tuple(_dist_sq(pt, ref_pts[rn]) for rn in ref_names)
        f_rows.append((name, *d2s))

    p_rows = []
    for name, pt in sp.p_series_pdf.items():
        d2s = tuple(_dist_sq(pt, ref_pts[rn]) for rn in ref_names)
        p_rows.append((name, *d2s))

    return f_rows, p_rows


# ============================================================
# Formatting
# ============================================================

def _format_row(row):
    """Format a d² row as a Python tuple string."""
    name = row[0]
    vals = ", ".join(f"{v:.6f}" for v in row[1:])
    return f'    ("{name}", {vals}),'


def _format_table(rows):
    """Format a complete d² table."""
    return "\n".join(_format_row(r) for r in rows)


# ============================================================
# Diffing
# ============================================================

def _parse_expected(text, var_name):
    """Parse an EXPECTED-style list from source text. Returns list of (name, *vals) tuples."""
    # Find the variable assignment
    pattern = rf'^{re.escape(var_name)}\s*=\s*\['
    match = re.search(pattern, text, re.MULTILINE)
    if not match:
        return []

    # Find the closing bracket
    start = match.end()
    depth = 1
    pos = start
    while depth > 0 and pos < len(text):
        if text[pos] == '[':
            depth += 1
        elif text[pos] == ']':
            depth -= 1
        pos += 1
    block = text[start:pos - 1]

    rows = []
    for m in re.finditer(r'\("([^"]+)",\s*([\d., e+\-]+)\)', block):
        name = m.group(1)
        vals = tuple(float(x.strip()) for x in m.group(2).split(","))
        rows.append((name, *vals))
    return rows


def _diff_tables(old_rows, new_rows, table_name, ref_names):
    """Compare old and new d² tables. Returns (changed_names, report_lines)."""
    old_dict = {r[0]: r[1:] for r in old_rows}
    new_dict = {r[0]: r[1:] for r in new_rows}
    tol = 1e-6

    changed = []
    added = []
    removed = []
    report = []

    # Check for name/count changes
    old_names = [r[0] for r in old_rows]
    new_names = [r[0] for r in new_rows]
    if old_names != new_names:
        for n in new_names:
            if n not in old_dict:
                added.append(n)
        for n in old_names:
            if n not in new_dict:
                removed.append(n)

    # Check value changes
    for name in new_names:
        if name not in old_dict:
            continue
        old_vals = old_dict[name]
        new_vals = new_dict[name]
        diffs = []
        for i, (ov, nv) in enumerate(zip(old_vals, new_vals)):
            if abs(ov - nv) > tol:
                diffs.append((ref_names[i], ov, nv, nv - ov))
        if diffs:
            changed.append(name)
            for rn, ov, nv, delta in diffs:
                report.append(
                    f"  {name} -> {rn}: {ov:.6f} -> {nv:.6f} (delta {delta:+.6f})")

    header = [f"\n--- {table_name} ---"]
    if not changed and not added and not removed:
        header.append("  No changes.")
    else:
        if changed:
            header.append(f"  {len(changed)} point(s) changed:")
        if added:
            header.append(f"  {len(added)} point(s) added: {', '.join(added)}")
        if removed:
            header.append(f"  {len(removed)} point(s) removed: {', '.join(removed)}")

    return changed + added + removed, header + report


# ============================================================
# File update
# ============================================================

def _replace_table(text, var_name, new_rows):
    """Replace an EXPECTED table in source text with new values."""
    pattern = rf'^({re.escape(var_name)}\s*=\s*\[)\n'
    match = re.search(pattern, text, re.MULTILINE)
    if not match:
        raise ValueError(f"Cannot find {var_name} in source text")

    # Find the closing bracket
    start = match.start()
    pos = match.end()
    depth = 1
    while depth > 0 and pos < len(text):
        if text[pos] == '[':
            depth += 1
        elif text[pos] == ']':
            depth -= 1
        pos += 1

    # Build replacement
    replacement = f"{var_name} = [\n{_format_table(new_rows)}\n]"
    return text[:start] + replacement + text[pos:]


# ============================================================
# Main
# ============================================================

def main():
    check_only = "--check" in sys.argv

    print("Computing layout d² values...")
    layout_rows = _compute_layout_d2()

    print("Computing site plan d² values...")
    site_f_rows, site_p_rows = _compute_site_d2()

    # Read current test files
    d2_path = os.path.join(_PROJECT, "tests", "test_d2_regression.py")
    sp_path = os.path.join(_PROJECT, "tests", "test_gen_site_plan.py")

    with open(d2_path, "r") as f:
        d2_text = f.read()
    with open(sp_path, "r") as f:
        sp_text = f.read()

    # Parse existing tables
    old_layout = _parse_expected(d2_text, "EXPECTED")
    old_site_f = _parse_expected(sp_text, "EXPECTED_D2")
    old_site_p = _parse_expected(sp_text, "EXPECTED_P_D2")

    # Diff
    all_report = []
    any_changes = False

    names1, lines1 = _diff_tables(
        old_layout, layout_rows, "EXPECTED (test_d2_regression.py)",
        ["F1", "F6", "F12", "F15"])
    all_report.extend(lines1)
    if names1:
        any_changes = True

    names2, lines2 = _diff_tables(
        old_site_f, site_f_rows, "EXPECTED_D2 (test_gen_site_plan.py)",
        ["CORNER_NW", "CORNER_NE", "CORNER_SE", "CORNER_SW"])
    all_report.extend(lines2)
    if names2:
        any_changes = True

    names3, lines3 = _diff_tables(
        old_site_p, site_p_rows, "EXPECTED_P_D2 (test_gen_site_plan.py)",
        ["CORNER_NW", "CORNER_NE", "CORNER_SE", "CORNER_SW"])
    all_report.extend(lines3)
    if names3:
        any_changes = True

    # Print report
    for line in all_report:
        print(line)

    if not any_changes:
        print("\nAll d² values are up to date.")
        return

    if check_only:
        print(f"\n--check: {len(names1) + len(names2) + len(names3)} total change(s) found. "
              "Run without --check to update.")
        sys.exit(1)

    # Update files
    print("\nUpdating test files...")

    new_d2_text = _replace_table(d2_text, "EXPECTED", layout_rows)
    with open(d2_path, "w") as f:
        f.write(new_d2_text)
    print(f"  Updated {d2_path}")

    new_sp_text = _replace_table(sp_text, "EXPECTED_D2", site_f_rows)
    new_sp_text = _replace_table(new_sp_text, "EXPECTED_P_D2", site_p_rows)
    with open(sp_path, "w") as f:
        f.write(new_sp_text)
    print(f"  Updated {sp_path}")

    print("\nDone. Run 'python -m pytest tests/ -x' to verify.")


if __name__ == "__main__":
    main()
