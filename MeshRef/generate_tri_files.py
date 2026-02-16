#!/usr/bin/env python3
"""Generate .tri files with custom KNPR sections based on a list of patterns."""

from __future__ import annotations

import argparse
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Create .tri files by copying a template and replacing the KNPR "
            "section with binary patterns listed in a text file."
        )
    )
    parser.add_argument(
        "--template",
        required=True,
        type=Path,
        help="Path to TEMPLATE.tri",
    )
    parser.add_argument(
        "--list",
        required=True,
        type=Path,
        help="Path to the file containing 8-character binary patterns (one per line).",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("."),
        help="Directory where the generated .tri files will be written.",
    )
    return parser.parse_args()


def load_template(template_path: Path) -> tuple[list[str], list[str]]:
    lines = template_path.read_text().splitlines()
    try:
        knpr_idx = lines.index("KNPR")
    except ValueError as exc:
        raise RuntimeError("Template file does not contain a KNPR section.") from exc

    after_knpr = lines[knpr_idx + 1 :]
    if len(after_knpr) < 8:
        raise RuntimeError("Template KNPR section must have at least 8 value lines.")

    prefix = lines[: knpr_idx + 1]
    suffix = after_knpr[8:]
    return prefix, suffix


def load_patterns(list_path: Path) -> list[str]:
    patterns: list[str] = []
    for raw_line in list_path.read_text().splitlines():
        pattern = raw_line.strip()
        if not pattern:
            continue
        if len(pattern) != 8 or any(ch not in {"0", "1"} for ch in pattern):
            raise ValueError(f"Invalid pattern '{raw_line}'. Patterns must be 8 binary digits.")
        patterns.append(pattern)
    if not patterns:
        raise RuntimeError("No valid patterns found in the list file.")
    return patterns


def write_files(prefix: list[str], suffix: list[str], patterns: list[str], output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    for pattern in patterns:
        knpr_lines = list(pattern)
        new_contents = prefix + knpr_lines + suffix
        out_path = output_dir / f"{pattern}.tri"
        out_path.write_text("\n".join(new_contents) + "\n")


def main() -> None:
    args = parse_args()
    prefix, suffix = load_template(args.template)
    patterns = load_patterns(args.list)
    write_files(prefix, suffix, patterns, args.output_dir)


if __name__ == "__main__":
    main()
