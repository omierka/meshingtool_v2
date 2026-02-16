#!/usr/bin/env python3
"""Batch driver for meshref."""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run meshref on *.tri inputs.")
    parser.add_argument(
        "--input-folder",
        required=True,
        type=Path,
        help="Folder containing *.tri input meshes",
    )
    parser.add_argument(
        "--output-folder",
        required=True,
        type=Path,
        help="Folder where *.vtu outputs are written",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()

    input_dir = args.input_folder.resolve()
    output_dir = args.output_folder.resolve()
    binary_dir = Path(__file__).resolve().parent
    meshref_binary = binary_dir / "meshref"

    if not input_dir.is_dir():
        print(f"Input folder not found: {input_dir}", file=sys.stderr)
        return 1

    if not meshref_binary.exists():
        print(f"meshref binary not found at {meshref_binary}", file=sys.stderr)
        return 1

    output_dir.mkdir(parents=True, exist_ok=True)

    tri_files = sorted(input_dir.glob("*.tri"))
    if not tri_files:
        print(f"No *.tri files found in {input_dir}", file=sys.stderr)
        return 1

    for tri_file in tri_files:
        output_file = output_dir / (tri_file.stem + ".vtu")
        print(f"Processing {tri_file} -> {output_file}")

        result = subprocess.run(
            [str(meshref_binary), "-i", str(tri_file), "-o", str(output_file)],
            cwd=binary_dir,
        )
        if result.returncode != 0:
            print(
                f"meshref failed for {tri_file} (exit code {result.returncode})",
                file=sys.stderr,
            )
            return result.returncode

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
