# Floorplan History Animation

Generates an MP4 video showing the evolution of `floorplan.svg` across the
repository's commit history.

## Prerequisites

- Python 3 with `Pillow` (`pip install Pillow`)
- Node.js with `sharp-cli` installed globally (`npm install -g sharp-cli`)
- `ffmpeg` on PATH
- Run from a bash shell (Git Bash on Windows)

## Steps

### 1. Generate SVGs from git history

```bash
bash ani/gen_animation.sh
```

This creates a detached git worktree at `ani/_worktree`, iterates every commit
from `22f34c0` (where `floorplan/gen_floorplan.py` was introduced) to `HEAD`,
checks out each commit, runs `gen_floorplan.py`, and copies the resulting SVG to
`ani/svgs/frame_NNNN.svg`. The worktree is removed when the script finishes.

`PYTHONPATH` is set to the worktree root so that Python module imports resolve
against the checked-out code rather than the main working tree.

Previously generated SVGs are cached -- the script skips any frame whose output
file already exists. To regenerate everything, delete `ani/svgs/` first.

Typical result: ~767 SVGs out of ~776 commits (the earliest 8-9 predate the
current module structure and fail harmlessly).

### 2. Convert to PNG and assemble the video

```bash
python ani/make_video.py
```

This runs five phases:

1. **SVG to PNG** -- `sharp-cli` renders each SVG at 1920 px wide, preserving
   the native aspect ratio, with a white background.
2. **Find max height** -- scans all PNGs to determine the tallest frame.
3. **Pad frames** -- centers each PNG on a white 1920 x max_height canvas so
   every frame has identical dimensions (H.264 requires this).
4. **Build concat list** -- writes `ani/frames.txt` for ffmpeg's concat
   demuxer, with each frame held for 7/30 s (~0.233 s).
5. **Encode MP4** -- H.264, yuv420p, CRF 18, 30 fps.

Output: `ani/floorplan_history.mp4`

Intermediate directories (`pngs/`, `padded/`) are cached -- delete them to
force a full re-render.

## Disk usage

| Directory | Typical size |
|-----------|-------------|
| `svgs/`   | ~70 MB      |
| `pngs/`   | ~145 MB     |
| `padded/` | ~158 MB     |
| MP4       | ~6 MB       |

Once the MP4 is satisfactory, `svgs/`, `pngs/`, and `padded/` can be deleted
to reclaim ~373 MB.

## Tuning

Edit the constants at the top of `make_video.py`:

| Variable           | Default | Effect                              |
|--------------------|---------|-------------------------------------|
| `WIDTH`            | 1920    | Frame width in pixels               |
| `FRAMES_PER_IMAGE` | 7       | How many video frames per floorplan |
| `FPS`              | 30      | Video frame rate                    |

Hold time per floorplan = `FRAMES_PER_IMAGE / FPS` (default 0.233 s).
