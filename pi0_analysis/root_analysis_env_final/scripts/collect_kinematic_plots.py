#!/usr/bin/env python3
"""Collect existing per-run plots into one inspection PDF per kinematic setting."""

from __future__ import annotations

import argparse
import os
import re
import shutil
import subprocess
import tempfile
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

from PIL import Image, ImageChops
from reportlab.lib.pagesizes import TABLOID, landscape
from reportlab.lib.utils import ImageReader
from reportlab.pdfgen import canvas


SOURCE_DIR = Path(
    "/w/hallc-scshelf2102/nps/singhav/nps_analysis/pi0_analysis/root_analysis_env_main/output/"
)
OUTPUT_DIR = Path(
    "/w/hallc-scshelf2102/nps/singhav/nps_analysis/pi0_analysis/"
    "root_analysis_env_main/output/plots_misc"
)
PLOT_SUFFIXES = {".png", ".jpg", ".jpeg", ".pdf"}
RUN_RE = re.compile(r"run_?(\d+)", re.IGNORECASE)
COLS, ROWS = 4, 4
PLOTS_PER_PAGE = COLS * ROWS
PDF_DPI = 120
# Known producer canvases; unknown images retain one slot. Never split a bitmap
# into inferred subplots or remove either PNG/PDF representation.
MULTIPANEL_PREFIXES = ("cluster_E_T_", "mass_cut_", "combbg_")


def plot_span(label: str) -> int:
    if label.startswith(MULTIPANEL_PREFIXES):
        return 2
    if label.startswith("cut_debug_"):
        page = re.search(r"\[(\d+)/(\d+)\]$", label)
        # The producer's final page is the single dead-block map.
        return 1 if page and page[1] == page[2] else 2
    return 1


def paginate(plots: list[tuple[Path, str]]) -> list[list[tuple[Path, str, int, int, int]]]:
    pages = []
    page = []
    occupied: set[tuple[int, int]] = set()
    for path, label in plots:
        span = plot_span(label)
        while True:
            slot = next(((r, c) for r in range(ROWS - span + 1)
                         for c in range(COLS - span + 1)
                         if all((r + dr, c + dc) not in occupied
                                for dr in range(span) for dc in range(span))), None)
            if slot is not None:
                break
            pages.append(page)
            page, occupied = [], set()
        row, col = slot
        occupied.update((row + dr, col + dc) for dr in range(span) for dc in range(span))
        page.append((path, label, row, col, span))
    if page:
        pages.append(page)
    return pages


def plot_image(path: Path) -> Image.Image:
    """Remove only exterior pure-white margins; preserve every nonwhite pixel."""
    with Image.open(path) as source:
        rgba = source.convert("RGBA")
        image = Image.new("RGB", rgba.size, "white")
        image.paste(rgba, mask=rgba.getchannel("A"))
    bounds = ImageChops.difference(image, Image.new("RGB", image.size, "white")).getbbox()
    if bounds:
        x0, y0, x1, y1 = bounds
        pad = 10
        image = image.crop((max(0, x0 - pad), max(0, y0 - pad),
                            min(image.width, x1 + pad), min(image.height, y1 + pad)))
    return image


def natural_key(value: str) -> list[object]:
    return [int(part) if part.isdigit() else part.lower() for part in re.split(r"(\d+)", value)]


def group_plot_files(plots_dir: Path) -> dict[int | None, list[Path]]:
    groups: dict[int | None, list[Path]] = defaultdict(list)
    files = (path for path in plots_dir.rglob("*") if path.is_file())
    for path in sorted(files, key=lambda item: natural_key(str(item.relative_to(plots_dir)))):
        if path.suffix.lower() not in PLOT_SUFFIXES:
            continue
        runs = RUN_RE.findall(str(path.relative_to(plots_dir)))
        groups[int(runs[-1]) if runs else None].append(path)
    return groups


def render_pdf_pages(pdf_path: Path, temp_dir: Path) -> list[Path]:
    prefix = temp_dir / pdf_path.stem
    subprocess.run(
        [
            "pdftoppm",
            "-png",
            "-r",
            str(PDF_DPI),
            str(pdf_path),
            str(prefix),
        ],
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.PIPE,
        text=True,
    )
    pages = sorted(temp_dir.glob(f"{pdf_path.stem}-*.png"), key=lambda p: natural_key(p.name))
    if not pages:
        raise RuntimeError(f"No pages rendered from {pdf_path}")
    return pages


def expand_plots(plot_files: list[Path], temp_root: Path) -> list[tuple[Path, str]]:
    expanded: list[tuple[Path, str]] = []
    for index, path in enumerate(plot_files):
        if path.suffix.lower() != ".pdf":
            expanded.append((path, path.name))
            continue

        pdf_temp = temp_root / f"pdf_{index:03d}"
        pdf_temp.mkdir()
        pages = render_pdf_pages(path, pdf_temp)
        expanded.extend(
            (page, f"{path.name} [{page_number}/{len(pages)}]")
            for page_number, page in enumerate(pages, start=1)
        )
    return expanded


def draw_page(
    pdf: canvas.Canvas,
    title: str,
    plots: list[tuple[Path, str, int, int, int]],
    page_number: int,
    page_count: int,
) -> None:
    page_width, page_height = landscape(TABLOID)
    margin, gap, title_height, label_height = 18, 6, 22, 11
    cell_width = (page_width - 2 * margin - (COLS - 1) * gap) / COLS
    cell_height = (page_height - 2 * margin - title_height - (ROWS - 1) * gap) / ROWS

    pdf.setFont("Helvetica-Bold", 12)
    pdf.drawString(margin, page_height - margin - 10, f"{title}  |  page {page_number}/{page_count}")

    for image_path, label, row, col, span in plots:
        x = margin + col * (cell_width + gap)
        width = span * cell_width + (span - 1) * gap
        height = span * cell_height + (span - 1) * gap
        y = page_height - margin - title_height - row * (cell_height + gap) - height

        image = ImageReader(plot_image(image_path))
        image_width, image_height = image.getSize()
        available_height = height - label_height
        scale = min(width / image_width, available_height / image_height)
        draw_width, draw_height = image_width * scale, image_height * scale
        draw_x = x + (width - draw_width) / 2
        draw_y = y + label_height + (available_height - draw_height) / 2
        pdf.drawImage(image, draw_x, draw_y, draw_width, draw_height, mask="auto")

        pdf.setFont("Helvetica", 6.5)
        pdf.drawCentredString(x + width / 2, y + 2, label[:90])

    pdf.showPage()


def make_kinematic_pdf(kin_dir: Path, output_dir: Path) -> tuple[Path, int, int]:
    groups = group_plot_files(kin_dir / "plots")
    run_numbers = sorted(run for run in groups if run is not None)
    if not run_numbers:
        raise RuntimeError(f"No run-tagged plot files found under {kin_dir / 'plots'}")

    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / f"all_plots_{kin_dir.name}.pdf"
    temporary_output = output_dir / f".{output_path.name}.tmp"
    page_total = 0

    try:
        pdf = canvas.Canvas(str(temporary_output), pagesize=landscape(TABLOID), pageCompression=1)
        with tempfile.TemporaryDirectory(prefix=f"collect_{kin_dir.name}_") as temp_name:
            temp_root = Path(temp_name)
            for group_index, run_number in enumerate(run_numbers):
                group_temp = temp_root / f"group_{group_index:04d}"
                group_temp.mkdir()
                plots = expand_plots(groups[run_number], group_temp)
                pages = paginate(plots)
                page_count = len(pages)
                title = f"{kin_dir.name} - run {run_number}"
                for page_index in range(page_count):
                    draw_page(
                        pdf,
                        title,
                        pages[page_index],
                        page_index + 1,
                        page_count,
                    )
                    page_total += 1
        pdf.save()
        os.replace(temporary_output, output_path)
    except Exception:
        temporary_output.unlink(missing_ok=True)
        raise

    return output_path, len(run_numbers), page_total


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, default=SOURCE_DIR)
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR)
    parser.add_argument(
        "--kin",
        action="append",
        help="Generate only this kinematic setting (repeat for more than one)",
    )
    parser.add_argument("--jobs", type=int, default=1, help="Kinematic PDFs to build concurrently")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if shutil.which("pdftoppm") is None:
        raise SystemExit("pdftoppm is required to include source PDF pages")
    if not args.input_dir.is_dir():
        raise SystemExit(f"Input directory does not exist: {args.input_dir}")

    available = {
        path.name: path
        for path in args.input_dir.glob("KinC_*")
        if path.is_dir() and (path / "plots").is_dir()
    }
    selected = args.kin or sorted(available, key=natural_key)
    unknown = [kin for kin in selected if kin not in available]
    if unknown:
        raise SystemExit(f"Unknown kinematic setting(s): {', '.join(unknown)}")
    if args.jobs < 1:
        raise SystemExit("--jobs must be at least 1")

    with ProcessPoolExecutor(max_workers=args.jobs) as executor:
        futures = {
            executor.submit(make_kinematic_pdf, available[kin], args.output_dir): kin
            for kin in selected
        }
        for future in as_completed(futures):
            output_path, run_count, page_count = future.result()
            print(f"Wrote {output_path} ({run_count} runs, {page_count} pages)", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
