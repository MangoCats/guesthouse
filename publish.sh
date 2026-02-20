#!/usr/bin/env bash
set -euo pipefail

branch=$(git rev-parse --abbrev-ref HEAD)

files=(
    floorplan/floorplan.svg
    floorplan/floorplan_db.svg
    floorplan/floorplan_bare.svg
    floorplan/floorplan_minik.svg
         roof/roof.svg
         span/span.svg
         span/span_min.svg
         site/site_plan.pdf
         site/site_plan_df.pdf
         site/site_plan_df.png
       survey/path_area.svg
        walls/walls.svg
        walls/all_walls.svg
         scad/2in12_patio.png
         scad/2in12_corner.png
         scad/2in12_bumpout.png
         scad/flat_roof_patio.png
         scad/flat_roof_corner.png
         scad/flat_roof_bumpout.png
)

# Files that are never branch-suffixed
fixed_files=(
    index.html
)

dest_files=()
for src in "${fixed_files[@]}"; do
    base=$(basename "$src")
    dest="../www/adu/$base"
    cp "$src" "$dest"
    dest_files+=("adu/$base")
done

for src in "${files[@]}"; do
    base=$(basename "$src")
    if [[ "$branch" != "main" ]]; then
        name="${base%.*}"
        ext="${base##*.}"
        base="${name}_${branch}.${ext}"
    fi
    dest="../www/adu/$base"
    cp "$src" "$dest"
    dest_files+=("adu/$base")
done

git -C ../www add "${dest_files[@]}"
git -C ../www commit -m "update floorplan and related files"
git -C ../www push
