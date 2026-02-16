#!/usr/bin/env bash
set -euo pipefail

cp floorplan/floorplan.svg       ../www/adu/floorplan.svg
cp floorplan/floorplan_minik.svg ../www/adu/floorplan_minik.svg
cp     walls/walls.svg           ../www/adu/walls.svg
cp     walls/all_walls.svg ../www/adu/all_walls.svg
cp    survey/path_area.svg ../www/adu/path_area.svg
git -C ../www add adu/floorplan.svg adu/floorplan_minik.svg adu/walls.svg adu/all_walls.svg adu/path_area.svg
git -C ../www commit -m "update floorplan and related files"
git -C ../www push
