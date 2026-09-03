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

from reportlab.lib.pagesizes import TABLOID, landscape
from reportlab.lib.utils import ImageReader
from reportlab.pdfgen import canvas


SOURCE_DIR = Path(
    "/lustre24/expphy/volatile/hallc/nps/singhav/nps_analysis/"
    "nps_analysis_20260831_165046"
)
OUTPUT_DIR = Path(
    "/w/hallc-scshelf2102/nps/singhav/nps_analysis/pi0_analysis/"
    "root_analysis_env_main/output/plots_misc"
)
PLOT_SUFFIXES = {".png", ".jpg", ".jpeg", ".pdf"}
RUN_RE = re.compile(r"run_?(\d+)", re.IGNORECASE)
COLS, ROWS = 5, 3
PLOTS_PER_PAGE = COLS * ROWS
PDF_DPI = 120


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
    plots: list[tuple[Path, str]],
    page_number: int,
    page_count: int,
) -> None:
    page_width, page_height = landscape(TABLOID)
    margin, gap, title_height, label_height = 18, 6, 22, 11
    cell_width = (page_width - 2 * margin - (COLS - 1) * gap) / COLS
    cell_height = (page_height - 2 * margin - title_height - (ROWS - 1) * gap) / ROWS

    pdf.setFont("Helvetica-Bold", 12)
    pdf.drawString(margin, page_height - margin - 10, f"{title}  |  page {page_number}/{page_count}")

    for slot, (image_path, label) in enumerate(plots):
        row, col = divmod(slot, COLS)
        x = margin + col * (cell_width + gap)
        y = page_height - margin - title_height - (row + 1) * cell_height - row * gap

        image = ImageReader(str(image_path))
        image_width, image_height = image.getSize()
        available_height = cell_height - label_height
        scale = min(cell_width / image_width, available_height / image_height)
        draw_width, draw_height = image_width * scale, image_height * scale
        draw_x = x + (cell_width - draw_width) / 2
        draw_y = y + label_height + (available_height - draw_height) / 2
        pdf.drawImage(image, draw_x, draw_y, draw_width, draw_height, mask="auto")

        pdf.setFont("Helvetica", 6.5)
        pdf.drawCentredString(x + cell_width / 2, y + 2, label[:90])

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
                page_count = (len(plots) + PLOTS_PER_PAGE - 1) // PLOTS_PER_PAGE
                title = f"{kin_dir.name} - run {run_number}"
                for page_index in range(page_count):
                    start = page_index * PLOTS_PER_PAGE
                    draw_page(
                        pdf,
                        title,
                        plots[start : start + PLOTS_PER_PAGE],
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
