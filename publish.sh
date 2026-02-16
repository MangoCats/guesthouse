#!/usr/bin/env bash
set -euo pipefail

branch=$(git rev-parse --abbrev-ref HEAD)

files=(
    floorplan/floorplan.svg
    floorplan/floorplan_minik.svg
    walls/walls.svg
    walls/all_walls.svg
    survey/path_area.svg
)

dest_files=()
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
