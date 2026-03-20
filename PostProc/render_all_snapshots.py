#!/usr/bin/env python3
"""
Render Filtered.vtu snapshots for every case that provides a view.pvsm file and
assemble a PDF report with LaTeX.

Usage:
    python render_all_snapshots.py --cases-root CASES
"""

import argparse
import re
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
IMAGE_HEIGHT_SHARE = 0.78  # fraction of page height reserved for case images
LATEX_MONITOR_CAPTION = r"\textrm{Element monitor distribution (\% of kept elements)}"
CaseRecord = Tuple[Path, Path, Optional[Path], Optional[str], Optional[str], Optional[str]]


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
        "--skip-render",
        action="store_true",
        help="Reuse existing Filtered/Surface PNGs and skip pvbatch rendering.",
    )
    parser.add_argument(
        "--keep-temp",
        action="store_true",
        help="Preserve the temporary LaTeX workspace/images instead of deleting them.",
    )
    return parser.parse_args()


def _candidate_case_dirs(cases_root: Path) -> List[Path]:
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
        return f"module purge && module load {shlex.quote(module_name)} && {pvbatch_cmd}"
    return pvbatch_cmd


def _trim_image(image_path: Path) -> bool:
    if not image_path.is_file():
        return False
    cmd = ["convert", str(image_path), "-trim", "+repage", str(image_path)]
    result = subprocess.run(cmd, check=False)
    if result.returncode != 0:
        print(f"  Failed to trim {image_path} (exit code {result.returncode}).", file=sys.stderr)
        return False
    return True


def _load_monitor_summary(case_dir: Path) -> Tuple[Optional[str], Optional[str]]:
    summary_path = case_dir / "monitor_summary.txt"
    volume_path = case_dir / "monitor_summary_volumetric.txt"

    def _read(path: Path) -> Optional[str]:
        if not path.is_file():
            return None
        try:
            return path.read_text(encoding="utf-8").strip()
        except OSError:
            return None

    return _read(summary_path), _read(volume_path)


def _extract_nel(case_dir: Path) -> Optional[str]:
    header_path = case_dir / "Filtered.tri"
    if not header_path.is_file():
        return None
    try:
        for line in header_path.open("r", encoding="utf-8", errors="ignore"):
            if "NEL" in line:
                parts = line.strip().split()
                if parts:
                    try:
                        int(parts[0])
                        return parts[0]
                    except ValueError:
                        continue
    except OSError:
        return None
    return None


def _latex_escape(text: str) -> str:
    replacements = {
        "\\": r"\textbackslash{}",
        "&": r"\&",
        "%": r"\%",
        "$": r"\$",
        "#": r"\#",
        "_": r"\_",
        "{": r"\{",
        "}": r"\}",
        "~": r"\textasciitilde{}",
        "^": r"\textasciicircum{}",
    }
    return "".join(replacements.get(ch, ch) for ch in text)


def _parse_monitor_values(summary: Optional[str]) -> Optional[List[str]]:
    if not summary:
        return None
    matches = re.findall(r"([0-3])=\[\s*([0-9.\-]+)%\]", summary)
    if len(matches) != 4:
        return None
    values = [""] * 4
    for bin_id, pct_text in matches:
        try:
            pct_value = float(pct_text)
        except ValueError:
            return None
        values[int(bin_id)] = f"{pct_value:.1f}\\%"
    return values


def _latex_monitor_table(
    count_values: Optional[List[str]], volume_values: Optional[List[str]]
) -> List[str]:
    rows: List[Tuple[str, Optional[List[str]]]] = []
    if count_values:
        rows.append(("Count", count_values))
    if volume_values:
        rows.append(("Volume", volume_values))
    if not rows:
        return []

    def _prepare(values: Optional[List[str]]) -> List[str]:
        if not values:
            return ["--"] * 4
        return [val if val else "--" for val in values]

    lines = [
        r"\begin{tabular}{lcccc}",
        r"\textbf{Method\textbackslash Bin} & \textbf{0} & \textbf{1} & \textbf{2} & \textbf{3} \\ \hline",
    ]
    for label, vals in rows:
        prepared = _prepare(vals)
        lines.append(
            "%s & %s & %s & %s & %s \\\\" % (
                _latex_escape(label),
                prepared[0],
                prepared[1],
                prepared[2],
                prepared[3],
            )
        )
    lines.append(r"\end{tabular}")
    return lines


def _build_cases_pdf(
    successful_cases: Sequence[CaseRecord],
    output_pdf: Path,
    keep_temp: bool,
) -> bool:
    tex_root = Path(tempfile.mkdtemp(prefix="cases_pdf_tex_"))
    tex_dir = tex_root / "tex"
    tex_dir.mkdir(parents=True, exist_ok=True)
    entries: List[dict] = []
    try:
        for idx, (case_dir, filtered, surface, summary_count, summary_volume, nel) in enumerate(successful_cases):
            images: List[str] = []
            for image_idx, source in enumerate([filtered, surface] if surface else [filtered]):
                if source is None or not source.is_file():
                    continue
                dest = tex_dir / f"case_{idx:04d}_img{image_idx:02d}{source.suffix.lower()}"
                try:
                    shutil.copy2(source, dest)
                except OSError as exc:
                    print(f"Failed to copy {source} -> {dest}: {exc}", file=sys.stderr)
                    return False
                images.append(dest.name)
            entries.append(
                {
                    "case_name": case_dir.name,
                    "nel": nel,
                    "summary_count": summary_count,
                    "summary_volume": summary_volume,
                    "images": images,
                }
            )

        tex_content = _render_latex_document(entries)
        tex_path = tex_dir / "report.tex"
        tex_path.write_text(tex_content, encoding="utf-8")
        pdflatex = shutil.which("pdflatex")
        if pdflatex is None:
            print("pdflatex executable not found in PATH.", file=sys.stderr)
            return False
        cmd = [pdflatex, "-interaction=nonstopmode", tex_path.name]
        result = subprocess.run(cmd, cwd=tex_dir, check=False)
        if result.returncode != 0:
            print("Failed to build PDF via pdflatex.", file=sys.stderr)
            return False
        pdf_path = tex_dir / "report.pdf"
        if not pdf_path.is_file():
            print("pdflatex did not produce report.pdf.", file=sys.stderr)
            return False
        shutil.move(str(pdf_path), str(output_pdf))
        return True
    finally:
        if keep_temp and tex_root.exists():
            print(f"[debug] Keeping LaTeX workspace under {tex_root}")
        elif tex_root.exists():
            shutil.rmtree(tex_root, ignore_errors=True)


def _render_latex_document(entries: Sequence[dict]) -> str:
    lines = [
        r"\documentclass[a4paper]{article}",
        r"\usepackage[T1]{fontenc}",
        r"\usepackage{graphicx}",
        r"\usepackage{geometry}",
        r"\geometry{margin=1.5cm}",
        r"\setlength{\parindent}{0pt}",
        r"\begin{document}",
    ]
    for idx, entry in enumerate(entries):
        if idx > 0:
            lines.append(r"\newpage")
        lines.append(r"\section*{%s}" % _latex_escape(entry["case_name"]))
        if entry.get("nel"):
            lines.append(r"\textbf{NEL=%s}\\[0.5em]" % _latex_escape(entry["nel"]))
        table_lines = _latex_monitor_table(
            _parse_monitor_values(entry.get("summary_count")),
            _parse_monitor_values(entry.get("summary_volume")),
        )
        if table_lines:
            lines.append(r"\begin{center}")
            lines.append(LATEX_MONITOR_CAPTION + r"\\[0.5em]")
            lines.extend(table_lines)
            lines.append(r"\vspace{0.5em}")
            lines.append(r"\end{center}")
        images = entry.get("images", [])
        if images:
            max_images = max(1, len(images))
            height_fraction = min(IMAGE_HEIGHT_SHARE / max_images, IMAGE_HEIGHT_SHARE)
            height_spec = f"{height_fraction:.3f}\\textheight"
            lines.append(r"\begin{center}")
            for image_idx, image_name in enumerate(images):
                lines.append(
                    r"\includegraphics[width=\linewidth,height=%s,keepaspectratio]{%s}"
                    % (height_spec, image_name)
                )
                if image_idx != len(images) - 1:
                    lines.append(r"\\[1em]")
            lines.append(r"\end{center}")
        lines.append(r"\bigskip")
    lines.append(r"\end{document}")
    return "\n".join(lines)


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

    if args.skip_render:
        commands: List[Optional[List[str]]] = [None] * len(cases)
    else:
        commands = [
            ["bash", "-lc", _build_batch_command(args.module, args.pvbatch, render_script, case)]
            for case in cases
        ]

    failures = 0
    successful_cases: List[CaseRecord] = []
    for idx, (case, cmd) in enumerate(zip(cases, commands), start=1):
        filtered_image = case / OUTPUT_FILE
        surface_image = case / SURFACE_OUTPUT_FILE
        if cmd is None:
            print(f"[{idx}/{len(cases)}] Using existing outputs in {case} ...")
        else:
            print(f"[{idx}/{len(cases)}] Rendering {case} ...")
            result = subprocess.run(cmd, cwd=str(case), check=False)
            if result.returncode != 0:
                failures += 1
                print(f"  Failed to render {case} (exit code {result.returncode}).", file=sys.stderr)
                continue
        if not filtered_image.is_file():
            print(f"  Missing {OUTPUT_FILE} in {case}.", file=sys.stderr)
            failures += 1
            continue
        if not args.skip_render and _trim_image(filtered_image):
            print(f"  Trimmed whitespace from {OUTPUT_FILE}.")
        if surface_image.is_file():
            if not args.skip_render:
                trimmed = _trim_image(surface_image)
            else:
                trimmed = True
            surface_available = surface_image if trimmed else None
            if surface_available and not args.skip_render:
                print(f"  Trimmed whitespace from {SURFACE_OUTPUT_FILE}.")
        else:
            surface_available = None
        summary_count, summary_volume = _load_monitor_summary(case)
        nel_value = _extract_nel(case)
        successful_cases.append(
            (case, filtered_image, surface_available, summary_count, summary_volume, nel_value)
        )

    if failures:
        print(f"Completed with {failures} failure(s).", file=sys.stderr)
    single_case_root = (cases_root / STATE_FILE).is_file()
    if not single_case_root and successful_cases:
        output_pdf = cases_root / "cases.pdf"
        if _build_cases_pdf(successful_cases, output_pdf, args.keep_temp):
            print(f"Created {output_pdf}.")
        else:
            print("Failed to create combined PDF.", file=sys.stderr)

    print("Snapshots generated for all cases.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
