#!/usr/bin/env python3
"""
Render Filtered.vtu snapshots for every case that provides a view.pvsm file.

The script loads ParaView via the module system before dispatching pvbatch for
each case:

    python render_all_snapshots.py

To skip the module load step, pass --module "" (empty string).
"""

import argparse
import shlex
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import List, Optional, Sequence, Tuple

STATE_FILE = "view.pvsm"
DATASET_FILE = "Filtered.vtu"
OUTPUT_FILE = "Filtered.png"
SURFACE_FILE = "surface.off"
SURFACE_OUTPUT_FILE = "Surface.png"
DEFAULT_MODULE = "paraview/5.13.0/opengl2-renderer/gui"
A4_WIDTH = 2480
A4_HEIGHT = 3508


def _parse_arguments(script_dir: Path) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run render_filtered_snapshot.py for every case with a view.pvsm file."
    )
    parser.add_argument(
        "--cases-root",
        default=str(script_dir / "CASES"),
        help="Directory containing the case folders (default: CASES next to this script).",
    )
    parser.add_argument(
        "--render-script",
        default=str(script_dir / "render_filtered_snapshot.py"),
        help="Path to render_filtered_snapshot.py (default: alongside this driver).",
    )
    parser.add_argument(
        "--module",
        default=DEFAULT_MODULE,
        help=(
            "Environment module to load before running pvbatch "
            f"(default: {DEFAULT_MODULE!r}; pass an empty string to skip)."
        ),
    )
    parser.add_argument(
        "--pvbatch",
        default="pvbatch",
        help="pvbatch executable to use (default: pvbatch resolved through the loaded module).",
    )
    parser.add_argument(
        "--skip-existing",
        action="store_true",
        help="Skip cases where the output PNG already exists.",
    )
    parser.add_argument(
        "--legacy-layout",
        action="store_true",
        help="Render PDF pages using the original montage layout instead of A4 pages.",
    )
    return parser.parse_args()


def _candidate_case_dirs(cases_root: Path) -> List[Path]:
    """Return either the single case (if the root is one) or its subdirectories."""
    single_case_state = cases_root / STATE_FILE
    if single_case_state.is_file():
        return [cases_root]
    return [entry for entry in sorted(cases_root.iterdir()) if entry.is_dir()]


def _discover_cases(cases_root: Path, skip_existing: bool) -> List[Path]:
    cases: List[Path] = []
    for entry in _candidate_case_dirs(cases_root):
        state_path = entry / STATE_FILE
        dataset_path = entry / DATASET_FILE
        output_path = entry / OUTPUT_FILE
        surface_geom_path = entry / SURFACE_FILE
        surface_output_path = entry / SURFACE_OUTPUT_FILE
        if not state_path.is_file():
            continue
        if not dataset_path.is_file():
            print(f"Skipping {entry}: missing dataset {dataset_path}", file=sys.stderr)
            continue
        surface_required = surface_geom_path.is_file()
        if skip_existing and output_path.is_file():
            if not surface_required or surface_output_path.is_file():
                print(
                    f"Skipping {entry}: existing snapshots detected.",
                    file=sys.stderr,
                )
                continue
        cases.append(entry)
    return cases


def _build_render_arguments(case_dir: Path) -> str:
    state_path = case_dir / STATE_FILE
    dataset_path = case_dir / DATASET_FILE
    output_path = case_dir / OUTPUT_FILE
    surface_path = case_dir / SURFACE_FILE
    surface_output_path = case_dir / SURFACE_OUTPUT_FILE
    return " ".join(
        [
            f"--state {shlex.quote(str(state_path))}",
            f"--dataset {shlex.quote(str(dataset_path))}",
            f"--output {shlex.quote(str(output_path))}",
            f"--surface {shlex.quote(str(surface_path))}",
            f"--surface-output {shlex.quote(str(surface_output_path))}",
        ]
    )


def _build_batch_command(
    module_name: str,
    pvbatch_exe: str,
    render_script: Path,
    case_dir: Path,
) -> str:
    render_args = _build_render_arguments(case_dir)
    pvbatch_cmd = " ".join(
        [
            shlex.quote(pvbatch_exe),
            "--force-offscreen-rendering",
            shlex.quote(str(render_script)),
            render_args,
        ]
    )
    if module_name:
        return f"module load {shlex.quote(module_name)} && {pvbatch_cmd}"
    return pvbatch_cmd


def _trim_image(image_path: Path) -> bool:
    """Trim whitespace from the image using ImageMagick."""
    if not image_path.is_file():
        return False
    cmd = ["convert", str(image_path), "-trim", "+repage", str(image_path)]
    result = subprocess.run(cmd, check=False)
    if result.returncode != 0:
        print(f"  Failed to trim {image_path} (exit code {result.returncode}).", file=sys.stderr)
        return False
    return True


def _load_monitor_summary(case_dir: Path) -> Optional[str]:
    summary_path = case_dir / "monitor_summary.txt"
    if not summary_path.is_file():
        return None
    try:
        return summary_path.read_text(encoding="utf-8").strip()
    except OSError:
        return None


def _extract_nel(case_dir: Path) -> Optional[str]:
    header_path = case_dir / "Filtered.tri"
    if not header_path.is_file():
        return None
    try:
        for line in header_path.open("r", encoding="utf-8", errors="ignore"):
            if "NEL" in line:
                parts = line.strip().split()
                if not parts:
                    continue
                try:
                    int(parts[0])
                    return parts[0]
                except ValueError:
                    continue
    except OSError:
        return None
    return None


def _create_pdf_page(
    case_name: str,
    summary: Optional[str],
    nel: Optional[str],
    images: Sequence[Path],
    output_path: Path,
    legacy_layout: bool,
) -> bool:
    inputs = [str(image) for image in images if image.is_file()]
    if not inputs:
        return False
    tile = f"1x{len(inputs)}"
    cmd = [
        "montage",
        *inputs,
        "-tile",
        tile,
        "-geometry",
        "+20+20",
        "-background",
        "white",
        "-title",
        "\n".join(filter(None, [case_name, f"NEL={nel}" if nel else None, summary])) if (summary or nel) else case_name,
        str(output_path),
    ]
    if legacy_layout:
        result = subprocess.run(cmd, check=False)
        return result.returncode == 0

    tmp_image = output_path.with_suffix(".tmp.png")
    result = subprocess.run([*cmd[:-1], str(tmp_image)], check=False)
    if result.returncode != 0:
        return False
    image_info = subprocess.run(
        ["identify", "-format", "%w %h", str(tmp_image)], capture_output=True, text=True, check=False
    )
    if image_info.returncode != 0:
        tmp_image.unlink(missing_ok=True)
        return False
    width_str, height_str = image_info.stdout.strip().split()
    width = int(width_str)
    height = int(height_str)

    scale_factor = min(A4_WIDTH / width, A4_HEIGHT / height)
    resize_arg = f"{int(width * scale_factor)}x{int(height * scale_factor)}!"
    convert_cmd = [
        "convert",
        str(tmp_image),
        "-resize",
        resize_arg,
        "-background",
        "white",
        "-gravity",
        "center",
        "-extent",
        f"{A4_WIDTH}x{A4_HEIGHT}",
        str(output_path),
    ]
    convert_result = subprocess.run(convert_cmd, check=False)
    try:
        tmp_image.unlink()
    except OSError:
        pass
    return convert_result.returncode == 0


def _build_cases_pdf(
    successful_cases: Sequence[Tuple[Path, Path, Optional[Path], Optional[str], Optional[str]]],
    output_pdf: Path,
    legacy_layout: bool,
) -> bool:
    tmp_dir = Path(tempfile.mkdtemp(prefix="cases_pdf_"))
    page_paths: List[Path] = []
    try:
        for idx, (case_dir, filtered, surface, summary, nel) in enumerate(successful_cases):
            images = [filtered]
            if surface and surface.is_file():
                images.append(surface)
            page_path = tmp_dir / f"page_{idx:04d}.png"
            if _create_pdf_page(case_dir.name, summary, nel, images, page_path, legacy_layout):
                page_paths.append(page_path)
            else:
                print(f"Failed to build PDF page for {case_dir}", file=sys.stderr)
        if not page_paths:
            return False
        cmd = ["convert", *[str(page) for page in page_paths], str(output_pdf)]
        result = subprocess.run(cmd, check=False)
        return result.returncode == 0
    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)


def main() -> int:
    script_dir = Path(__file__).resolve().parent
    args = _parse_arguments(script_dir)

    cases_root = Path(args.cases_root).expanduser().resolve()
    render_script = Path(args.render_script).expanduser().resolve()

    if not cases_root.is_dir():
        print(f"Cases root {cases_root} does not exist or is not a directory.", file=sys.stderr)
        return 1
    if not render_script.is_file():
        print(f"Render script {render_script} does not exist.", file=sys.stderr)
        return 1

    cases = _discover_cases(cases_root, args.skip_existing)
    if not cases:
        print("No cases with both view.pvsm and Filtered.vtu were found.")
        return 0

    commands = [
        ["bash", "-lc", _build_batch_command(args.module, args.pvbatch, render_script, case)]
        for case in cases
    ]

    failures = 0
    successful_cases: List[Tuple[Path, Path, Optional[Path], Optional[str], Optional[str]]] = []
    for idx, (case, cmd) in enumerate(zip(cases, commands), start=1):
        print(f"[{idx}/{len(cases)}] Rendering {case} ...")
        result = subprocess.run(cmd, cwd=str(case), check=False)
        filtered_image = case / OUTPUT_FILE
        surface_image = case / SURFACE_OUTPUT_FILE
        if result.returncode != 0:
            failures += 1
            print(f"  Failed to render {case} (exit code {result.returncode}).", file=sys.stderr)
            continue
        if _trim_image(filtered_image):
            print(f"  Trimmed whitespace from {OUTPUT_FILE}.")
        surface_available = surface_image if _trim_image(surface_image) else None
        if surface_available:
            print(f"  Trimmed whitespace from {SURFACE_OUTPUT_FILE}.")
        summary_line = _load_monitor_summary(case)
        nel_value = _extract_nel(case)
        successful_cases.append((case, filtered_image, surface_available, summary_line, nel_value))

    if failures:
        print(f"Completed with {failures} failure(s).", file=sys.stderr)
        return 1

    single_case_root = (cases_root / STATE_FILE).is_file()
    if not single_case_root and successful_cases:
        output_pdf = cases_root / "cases.pdf"
        if _build_cases_pdf(successful_cases, output_pdf, args.legacy_layout):
            print(f"Created {output_pdf}.")
        else:
            print("Failed to create combined PDF.", file=sys.stderr)

    print("Snapshots generated for all cases.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
