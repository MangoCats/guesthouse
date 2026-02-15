#!/usr/bin/env bash
set -euo pipefail

cp floorplan/floorplan.svg ../www/adu/floorplan.svg
git -C ../www add adu/floorplan.svg
git -C ../www commit -m "update floorplan"
git -C ../www push
