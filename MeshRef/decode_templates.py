#!/usr/bin/env python3
"""Decode 8-bit refinement patterns into canonical template IDs.

Usage:
  python decode_templates.py 11101000 10111110
  python decode_templates.py list
  python decode_templates.py --file patterns.txt

Patterns may optionally be wrapped in brackets (e.g. "[10101010]").
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable, Iterator, List, Sequence, Tuple

CANONICAL_TEMPLATES: Tuple[str, ...] = (
    "00000000",
    "10000000",
    "11000000",
    "10100000",
    "10000010",
    "11100000",
    "10101000",
    "10100100",
    "11110000",
    "11011000",
    "10111000",
    "10110100",
    "10101010",
    "10100101",
    "11111000",
    "10111100",
    "11011010",
    "11111100",
    "11111010",
    "10111110",
    "11111110",
    "11111111",
    "11101000",
)

RROT: Tuple[Tuple[int, ...], ...] = (
    (1, 2, 3, 4, 5, 6, 7, 8),
    (2, 3, 4, 1, 6, 7, 8, 5),
    (3, 4, 1, 2, 7, 8, 5, 6),
    (4, 1, 2, 3, 8, 5, 6, 7),
    (5, 8, 7, 6, 1, 4, 3, 2),
    (6, 5, 8, 7, 2, 1, 4, 3),
    (7, 6, 5, 8, 3, 2, 1, 4),
    (8, 7, 6, 5, 4, 3, 2, 1),
)

R2ROT: Tuple[int, ...] = (1, 5, 6, 2, 4, 8, 7, 3)
R3ROT: Tuple[int, ...] = (1, 4, 8, 5, 2, 3, 7, 6)


def parse_code(raw: str) -> Tuple[List[bool], str]:
    digits = [ch for ch in raw if ch in "01"]
    if len(digits) != 8:
        raise ValueError(f"expected 8 binary digits, got {raw!r}")
    return [ch == "1" for ch in digits], "".join(digits)


def bools_to_string(bits: Sequence[bool]) -> str:
    return "".join("1" if bit else "0" for bit in bits)


def rotate_pattern(bits: Sequence[bool]) -> List[bool]:
    jsum_init = 64 * 8
    jsum = jsum_init
    best = list(bits)

    for idx in range(8):
        if not bits[idx]:
            continue

        r1 = [bits[i - 1] for i in RROT[idx]]
        isum = sum((pos + 1) ** 2 for pos, val in enumerate(r1) if val)
        if isum < jsum:
            best = r1
            jsum = isum

        r2 = [r1[i - 1] for i in R2ROT]
        isum = sum((pos + 1) ** 2 for pos, val in enumerate(r2) if val)
        if isum < jsum:
            best = r2
            jsum = isum

        r3 = [r1[i - 1] for i in R3ROT]
        isum = sum((pos + 1) ** 2 for pos, val in enumerate(r3) if val)
        if isum < jsum:
            best = r3
            jsum = isum

    return best


def decode(bits: Sequence[bool]) -> Tuple[int, str, bool] | None:
    as_str = bools_to_string(bits)
    for idx, tmpl in enumerate(CANONICAL_TEMPLATES, start=1):
        if as_str == tmpl:
            return idx, tmpl, False

    rotated = rotate_pattern(bits)
    rotated_str = bools_to_string(rotated)
    for idx, tmpl in enumerate(CANONICAL_TEMPLATES, start=1):
        if rotated_str == tmpl:
            return idx, tmpl, True

    return None


def read_code_file(path: Path) -> Iterator[str]:
    for line in path.read_text().splitlines():
        raw = line.strip()
        if raw:
            yield raw


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Decode 8-bit refinement patterns into canonical template IDs."
    )
    parser.add_argument("codes", nargs="*", help="8-bit sequences or filenames")
    parser.add_argument(
        "-f",
        "--file",
        dest="files",
        action="append",
        help="Optional file containing one pattern per line (can be specified multiple times)",
    )
    args = parser.parse_args()

    if not args.codes and not args.files:
        parser.error("provide codes as arguments or via --file")

    inputs: List[str] = []
    if args.files:
        for filename in args.files:
            inputs.extend(read_code_file(Path(filename)))

    for token in args.codes:
        path = Path(token)
        if path.is_file():
            inputs.extend(read_code_file(path))
        else:
            inputs.append(token)

    for raw in inputs:
        try:
            bits, digits = parse_code(raw)
        except ValueError as exc:
            print(f"{raw!r}: {exc}")
            continue

        result = decode(bits)
        if result is None:
            print(f"{digits}: no canonical match (even after rotation)")
            continue

        idx, tmpl, rotated = result
        status = "rotated" if rotated else "direct"
        print(f"{digits} -> {tmpl} [{status}] (template {idx:02d})")


if __name__ == "__main__":
    main()
