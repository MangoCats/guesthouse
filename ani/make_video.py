"""Convert SVGs to PNGs and assemble into MP4 video."""
import os, sys, subprocess, glob
from PIL import Image

ANI_DIR = os.path.dirname(os.path.abspath(__file__))
SVG_DIR = os.path.join(ANI_DIR, "svgs")
PNG_DIR = os.path.join(ANI_DIR, "pngs")
PADDED_DIR = os.path.join(ANI_DIR, "padded")
OUTPUT = os.path.join(ANI_DIR, "floorplan_history.mp4")
WIDTH = 1920
FRAMES_PER_IMAGE = 7
FPS = 30

os.makedirs(PNG_DIR, exist_ok=True)
os.makedirs(PADDED_DIR, exist_ok=True)

# Get sorted list of SVG files
svg_files = sorted(glob.glob(os.path.join(SVG_DIR, "frame_*.svg")))
total = len(svg_files)
print(f"Found {total} SVG files to process")

# Phase 1: Convert SVGs to PNGs at 1920px wide, native aspect ratio
print("\n=== Phase 1: Converting SVGs to PNGs ===")
for i, svg in enumerate(svg_files, 1):
    base = os.path.splitext(os.path.basename(svg))[0]
    png = os.path.join(PNG_DIR, f"{base}.png")
    if os.path.exists(png):
        sys.stdout.write(f"\r[{i}/{total}] skip (cached) {base}   ")
        sys.stdout.flush()
        continue
    try:
        cmd = f'npx sharp-cli -i "{svg}" -o "{png}" -- resize {WIDTH} -- flatten "#ffffff"'
        subprocess.run(cmd, check=True, capture_output=True, timeout=30, shell=True)
        sys.stdout.write(f"\r[{i}/{total}] OK   {base}   ")
    except Exception as e:
        sys.stdout.write(f"\r[{i}/{total}] FAIL {base} ({e})   ")
    sys.stdout.flush()
print()

# Phase 2: Find max dimensions
print("\n=== Phase 2: Finding max dimensions ===")
png_files = sorted(glob.glob(os.path.join(PNG_DIR, "frame_*.png")))
max_h = 0
for png in png_files:
    im = Image.open(png)
    w, h = im.size
    if h > max_h:
        max_h = h
    im.close()

# Make height even (H.264 requirement)
max_h = (max_h + 1) // 2 * 2
print(f"Uniform frame size: {WIDTH}x{max_h}")

# Phase 3: Create padded frames (centered on white background)
print("\n=== Phase 3: Creating padded frames ===")
for i, png in enumerate(png_files, 1):
    base = os.path.basename(png)
    padded = os.path.join(PADDED_DIR, base)
    if os.path.exists(padded):
        sys.stdout.write(f"\r[{i}/{len(png_files)}] skip (cached) {base}   ")
        sys.stdout.flush()
        continue
    im = Image.open(png)
    w, h = im.size
    canvas = Image.new("RGB", (WIDTH, max_h), (255, 255, 255))
    x = (WIDTH - w) // 2
    y = (max_h - h) // 2
    # Handle RGBA -> RGB paste
    if im.mode == "RGBA":
        canvas.paste(im, (x, y), im)
    else:
        canvas.paste(im, (x, y))
    canvas.save(padded)
    im.close()
    canvas.close()
    sys.stdout.write(f"\r[{i}/{len(png_files)}] OK   {base}   ")
    sys.stdout.flush()
print()

# Phase 4: Build ffmpeg concat file
print("\n=== Phase 4: Building frame list ===")
padded_files = sorted(glob.glob(os.path.join(PADDED_DIR, "frame_*.png")))
concat_file = os.path.join(ANI_DIR, "frames.txt")
duration = FRAMES_PER_IMAGE / FPS

with open(concat_file, "w") as f:
    for png in padded_files:
        abs_path = os.path.abspath(png).replace("\\", "/")
        f.write(f"file '{abs_path}'\n")
        f.write(f"duration {duration:.6f}\n")
    # Last frame repeated (ffmpeg concat demuxer quirk)
    abs_path = os.path.abspath(padded_files[-1]).replace("\\", "/")
    f.write(f"file '{abs_path}'\n")

frame_count = len(padded_files)
video_duration = frame_count * FRAMES_PER_IMAGE / FPS
print(f"Total unique frames: {frame_count}")
print(f"Video duration: {video_duration:.1f} seconds ({video_duration/60:.1f} min)")

# Phase 5: Encode MP4
print("\n=== Phase 5: Encoding MP4 ===")
cmd = [
    "ffmpeg", "-y", "-f", "concat", "-safe", "0", "-i", concat_file,
    "-c:v", "libx264", "-pix_fmt", "yuv420p", "-crf", "18",
    "-r", str(FPS),
    OUTPUT
]
print(f"Running: {' '.join(cmd)}")
result = subprocess.run(cmd, capture_output=True, text=True)
if result.returncode != 0:
    print(f"ffmpeg stderr:\n{result.stderr[-2000:]}")
    sys.exit(1)

size_mb = os.path.getsize(OUTPUT) / (1024 * 1024)
print(f"\nDone! Output: {OUTPUT}")
print(f"File size: {size_mb:.1f} MB")
