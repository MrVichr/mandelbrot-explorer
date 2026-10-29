#!/usr/bin/env python3
"""check_images.py — regression checker for mandel_runner.

Invokes mandel_runner in the requested mode (text or images), then compares
the output against the "blessed" reference in ./blessed/.

Text mode  — compares the JSON file byte-for-byte against blessed/
Images mode — converts PPMs to PNG for inspection and compares raw PPM
               bytes against those in blessed/

Usage
-----
    python3 runner/check_images.py --fractal man|jul|lag --mode text   [options] [--bless]
    python3 runner/check_images.py --fractal man|jul|lag --mode images [options] [--bless]

Options
-------
    --fractal man|jul|lag    Fractal type (default: man)
    --mode text|images       Which artefacts to produce / compare (required)
    --bless                  Run with --bless so artefacts go to blessed/
    --runner PATH            Path to mandel_runner binary (default: ./mandel_runner)
    --cx NUM                 Centre real part              (default: -0.5)
    --cy NUM                 Centre imaginary part         (default: 0.0)
    --scale FLOAT            Scale                         (default: 0.125)
    --width INT              Image width in pixels         (default: 40)
    --height INT             Image height in pixels        (default: 20)
    --title STR              Label for the run             (default: Wide view)
    --c-re NUM               Julia/Laguerre c, real part   (default: runner's)
    --c-im NUM               Julia/Laguerre c, imag part   (default: runner's)
    --period INT             Julia/Laguerre period         (default: runner's)
    --precision NAME         double|float128|dd|qd|real642|ddheap (default: double);
                             non-double outputs are tagged, e.g. Man_dd_...

Exit codes
----------
    0 — all checks passed (or bless succeeded)
    1 — at least one artefact differs from its blessed copy
    2 — runner failed or some other error

Notes
-----
Image comparison is done on raw PPM bytes, independent of PNG encoding.
PNG files are generated in current/ solely for human inspection.
If Pillow is not installed, PNG generation is silently skipped.
"""

import argparse
import subprocess
import sys
from pathlib import Path

try:
    from PIL import Image
    HAS_PILLOW = True
except ImportError:
    HAS_PILLOW = False

# ---------------------------------------------------------------------------
# Paint-style names and filename prefixes — per fractal type.
# Must match the tables in main_runner.cpp.
# ---------------------------------------------------------------------------
STYLE_NAMES: dict[str, list[str]] = {
    "man": ["Kind", "Cls", "Exter", "ExterAngle", "Inter", "Near", "FZ", "FC"],
    "jul": ["Kind", "Cls", "Exter", "ExterAngle", "Inter",
            "Near0", "NearZ", "Phi", "Phi1", "Phi2", "F1", "FirstUnder1", "NearMR"],
    "lag": ["Cls", "NthFz", "NthFz1", "Roots"],
}

FRACTAL_PREFIX: dict[str, str] = {
    "man": "Man_",
    "jul": "Jul_",
    "lag": "Lag_",
}

HERE        = Path(__file__).parent
CURRENT_DIR = Path("current")
BLESSED_DIR = Path("blessed")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def sanitize_title(title: str) -> str:
    """Mirror the C++ sanitize_title: replace non-alphanumeric chars with '_'."""
    return "".join(c if c.isalnum() else "_" for c in title)


def file_prefix(args) -> str:
    """Mirror the C++ file_prefix: fractal prefix plus precision tag."""
    prefix = FRACTAL_PREFIX[args.fractal]
    if args.precision != "double":
        prefix += f"{args.precision}_"
    return prefix


def ppm_to_png(ppm_path: Path, png_path: Path) -> None:
    """Convert a PPM file to PNG using Pillow.  Skips if Pillow unavailable."""
    if not HAS_PILLOW:
        return
    Image.open(ppm_path).save(png_path)


def run_runner(runner: Path, args, mode: str, bless: bool) -> int:
    """Invoke mandel_runner with the given mode and bless flag."""
    cmd = [
        str(runner),
        "--fractal", args.fractal,
        "--mode",   mode,
        "--cx",     str(args.cx),
        "--cy",     str(args.cy),
        "--scale",  str(args.scale),
        "--width",  str(args.width),
        "--height", str(args.height),
        "--title",  args.title,
        "--precision", args.precision,
    ]
    if args.c_re is not None:
        cmd += ["--c-re", str(args.c_re)]
    if args.c_im is not None:
        cmd += ["--c-im", str(args.c_im)]
    if args.period is not None:
        cmd += ["--period", str(args.period)]
    if bless:
        cmd.append("--bless")
    print("Running:", " ".join(cmd))
    result = subprocess.run(cmd)
    return result.returncode


# ---------------------------------------------------------------------------
# Text mode
# ---------------------------------------------------------------------------

def bless_text(runner: Path, args) -> int:
    """Run text mode with --bless so the JSON lands in blessed/."""
    rc = run_runner(runner, args, "text", bless=True)
    if rc != 0:
        print(f"ERROR: mandel_runner exited with code {rc}", file=sys.stderr)
        return 2
    prefix = file_prefix(args)
    stem = sanitize_title(args.title)
    blessed_json = BLESSED_DIR / f"{prefix}{stem}_text.json"
    if blessed_json.exists():
        print(f"  Blessed: {blessed_json}")
    else:
        print(f"WARNING: {blessed_json} not produced by runner.", file=sys.stderr)
    return 0


def check_text(runner: Path, args) -> int:
    """Run text mode into current/ and compare against blessed/."""
    if not BLESSED_DIR.exists():
        print(
            f"ERROR: No blessed artefacts found ({BLESSED_DIR}).\n"
            "       Run with --bless first.",
            file=sys.stderr,
        )
        return 2

    rc = run_runner(runner, args, "text", bless=False)
    if rc != 0:
        print(f"ERROR: mandel_runner exited with code {rc}", file=sys.stderr)
        return 2

    prefix = file_prefix(args)
    stem        = sanitize_title(args.title)
    current_json = CURRENT_DIR / f"{prefix}{stem}_text.json"
    blessed_json = BLESSED_DIR / f"{prefix}{stem}_text.json"

    if not current_json.exists():
        print(f"  SKIP (runner did not produce): {current_json}", file=sys.stderr)
        return 2
    if not blessed_json.exists():
        print(f"  SKIP (no blessed JSON): {blessed_json}", file=sys.stderr)
        print("\nText check: SKIP")
        return 2

    if current_json.read_bytes() == blessed_json.read_bytes():
        print(f"  PASS: text JSON")
        print("\nText check: PASS")
        return 0
    else:
        print(f"  FAIL: {current_json} differs from {blessed_json}")
        print("\nText check: FAIL", file=sys.stderr)
        return 1


# ---------------------------------------------------------------------------
# Images mode
# ---------------------------------------------------------------------------

def bless_images(runner: Path, args) -> int:
    """Run images mode with --bless so PPMs land in blessed/."""
    rc = run_runner(runner, args, "images", bless=True)
    if rc != 0:
        print(f"ERROR: mandel_runner exited with code {rc}", file=sys.stderr)
        return 2

    blessed = 0
    prefix = file_prefix(args)
    stem = sanitize_title(args.title)
    for name in STYLE_NAMES[args.fractal]:
        ppm = BLESSED_DIR / f"{prefix}{stem}_{name}.ppm"
        if not ppm.exists():
            print(f"WARNING: {ppm} not produced by runner — skipping.", file=sys.stderr)
            continue
        png = ppm.with_suffix(".png")
        ppm_to_png(ppm, png)
        print(f"  Blessed PPM: {ppm}" + (f" → PNG: {png}" if HAS_PILLOW else ""))
        blessed += 1

    print(f"\nBlessed {blessed} images to {BLESSED_DIR}/")
    if not HAS_PILLOW:
        print("Note: Pillow not available — PNG files were not generated.\n"
              "      Install with: pip install Pillow")
    return 0


def check_images(runner: Path, args) -> int:
    """Run images mode into current/ and compare against blessed/."""
    if not BLESSED_DIR.exists():
        print(
            f"ERROR: No blessed artefacts found ({BLESSED_DIR}).\n"
            "       Run with --bless first.",
            file=sys.stderr,
        )
        return 2

    rc = run_runner(runner, args, "images", bless=False)
    if rc != 0:
        print(f"ERROR: mandel_runner exited with code {rc}", file=sys.stderr)
        return 2

    # Convert produced PPMs to PNG for easy inspection.
    prefix = file_prefix(args)
    stem = sanitize_title(args.title)
    for name in STYLE_NAMES[args.fractal]:
        ppm = CURRENT_DIR / f"{prefix}{stem}_{name}.ppm"
        if ppm.exists():
            ppm_to_png(ppm, ppm.with_suffix(".png"))

    all_pass = True
    all_exist = True
    for name in STYLE_NAMES[args.fractal]:
        src_ppm   = CURRENT_DIR / f"{prefix}{stem}_{name}.ppm"
        bless_ppm = BLESSED_DIR / f"{prefix}{stem}_{name}.ppm"

        if not src_ppm.exists():
            print(f"  SKIP (runner did not produce): {name}")
            all_exist = False
            continue
        if not bless_ppm.exists():
            print(f"  SKIP (no blessed):             {name}")
            all_exist = False
            continue

        if src_ppm.read_bytes() == bless_ppm.read_bytes():
            print(f"  PASS: {name}")
        else:
            print(f"  FAIL: {name}  (see {src_ppm.with_suffix('.png')})")
            all_pass = False

    if not all_exist:
        print("\nImage check: SKIP")
        return 0
    elif all_pass:
        print("\nImage check: PASS")
        return 0
    else:
        print("\nImage check: FAIL", file=sys.stderr)
        return 1


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> int:
    parser = argparse.ArgumentParser(
        description="Regression checker for mandel_runner (text and images modes).",
    )
    parser.add_argument("--fractal", default="man", choices=["man", "jul", "lag"],
                        help="Fractal type: man(delbrot), jul(ia), lag(uerre) (default: man)")
    parser.add_argument("--mode",   required=True, choices=["text", "images"],
                        help="Which artefacts to produce and compare")
    parser.add_argument("--bless",  action="store_true",
                        help="Record current output as the blessed reference")
    parser.add_argument("--runner", default=str(HERE.parent / "mandel_runner"),
                        help="Path to mandel_runner binary")
    # NUM values stay strings: the runner parses them in the selected precision
    parser.add_argument("--cx",     default="-0.5")
    parser.add_argument("--cy",     default="0.0")
    parser.add_argument("--scale",  type=float, default=0.125)
    parser.add_argument("--width",  type=int,   default=80)
    parser.add_argument("--height", type=int,   default=60)
    parser.add_argument("--title",  default="XWide view")
    parser.add_argument("--c-re",   dest="c_re", default=None)
    parser.add_argument("--c-im",   dest="c_im", default=None)
    parser.add_argument("--period", type=int, default=None)
    parser.add_argument("--precision", default="double",
                        choices=["double", "float128", "dd", "qd", "real642", "ddheap"])
    args = parser.parse_args()

    runner = Path(args.runner)
    if not runner.exists():
        print(f"ERROR: runner not found: {runner}\n"
              "       Run 'make' first.", file=sys.stderr)
        return 2

    if args.mode == "text":
        if args.bless:
            return bless_text(runner, args)
        else:
            return check_text(runner, args)
    else:  # images
        if args.bless:
            return bless_images(runner, args)
        else:
            return check_images(runner, args)


if __name__ == "__main__":
    sys.exit(main())
