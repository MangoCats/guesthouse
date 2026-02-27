#!/bin/bash
# Generate floorplan animation from git history
# Creates SVGs by checking out each commit in a worktree and running gen_floorplan.py
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
ANI_DIR="$REPO_ROOT/ani"
SVG_DIR="$ANI_DIR/svgs"
WORKTREE_DIR="$ANI_DIR/_worktree"

mkdir -p "$SVG_DIR"

# Clean up any leftover worktree
if [ -d "$WORKTREE_DIR" ]; then
    cd "$REPO_ROOT"
    git worktree remove --force "$WORKTREE_DIR" 2>/dev/null || rm -rf "$WORKTREE_DIR"
fi

# Get all commits from when gen_floorplan.py was introduced (22f34c0) to HEAD
COMMITS_FILE="$ANI_DIR/commits.txt"
cd "$REPO_ROOT"
git log --reverse --format="%H %s" 22f34c0^..HEAD > "$COMMITS_FILE"
TOTAL=$(wc -l < "$COMMITS_FILE")
echo "Found $TOTAL commits to process"

# Create a worktree for generation
git worktree add "$WORKTREE_DIR" HEAD --detach 2>/dev/null

COUNT=0
SKIPPED=0
GENERATED=0

while IFS=' ' read -r HASH MSG; do
    COUNT=$((COUNT + 1))
    PADDED=$(printf "%04d" $COUNT)
    SVG_OUT="$SVG_DIR/frame_${PADDED}.svg"

    # Skip if already generated
    if [ -f "$SVG_OUT" ]; then
        GENERATED=$((GENERATED + 1))
        continue
    fi

    # Checkout this commit in the worktree
    cd "$WORKTREE_DIR"
    git checkout --force "$HASH" 2>/dev/null || {
        SKIPPED=$((SKIPPED + 1))
        printf "\r[%d/%d] SKIP (checkout) %s" "$COUNT" "$TOTAL" "${HASH:0:8}"
        continue
    }
    git clean -fdx 2>/dev/null

    # Check if gen_floorplan.py exists
    if [ ! -f "floorplan/gen_floorplan.py" ]; then
        SKIPPED=$((SKIPPED + 1))
        printf "\r[%d/%d] SKIP (no script) %s" "$COUNT" "$TOTAL" "${HASH:0:8}"
        continue
    fi

    # Remove old output
    rm -f floorplan/floorplan.svg

    # Generate with PYTHONPATH set to worktree root so imports resolve correctly
    if PYTHONPATH="$WORKTREE_DIR" python floorplan/gen_floorplan.py >/dev/null 2>&1; then
        if [ -f "floorplan/floorplan.svg" ]; then
            cp "floorplan/floorplan.svg" "$SVG_OUT"
            GENERATED=$((GENERATED + 1))
            printf "\r[%d/%d] OK   %s %.50s" "$COUNT" "$TOTAL" "${HASH:0:8}" "$MSG"
        else
            SKIPPED=$((SKIPPED + 1))
            printf "\r[%d/%d] SKIP (no output) %s" "$COUNT" "$TOTAL" "${HASH:0:8}"
        fi
    else
        SKIPPED=$((SKIPPED + 1))
        printf "\r[%d/%d] FAIL %s %.50s" "$COUNT" "$TOTAL" "${HASH:0:8}" "$MSG"
    fi
done < "$COMMITS_FILE"

echo ""
echo "Generation complete: $GENERATED succeeded, $SKIPPED skipped out of $TOTAL"

# Clean up worktree
cd "$REPO_ROOT"
git worktree remove --force "$WORKTREE_DIR" 2>/dev/null || true

echo "SVGs saved in $SVG_DIR"
