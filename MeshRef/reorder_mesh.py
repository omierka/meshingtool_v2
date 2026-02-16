#!/usr/bin/env python3
"""Utility for reordering vertices in TRI3D hex mesh files."""

from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass
from typing import List, Sequence, Tuple

TARGET_VERTICES: Tuple[Tuple[float, float, float], ...] = (
    (-1.0, -1.0, -1.0),
    ( 1.0, -1.0, -1.0),
    ( 1.0, +1.0, -1.0),
    (-1.0, +1.0, -1.0),
    (-1.0, -1.0, +1.0),
    ( 1.0, -1.0, +1.0),
    ( 1.0, +1.0, +1.0),
    (-1.0, +1.0, +1.0),
)

COORD_TOL = 1e-7


@dataclass
class MeshData:
    header: List[str]
    counts_line: str
    vertices: List[Tuple[Tuple[float, float, float], str]]
    connectivity: List[List[int]]
    nodal_props: List[str]
    tail: List[str]


def parse_counts(line: str) -> Tuple[int, int]:
    tokens = line.split()
    ints: List[int] = []
    for token in tokens:
        try:
            ints.append(int(token))
        except ValueError:
            break
    if len(ints) < 2:
        raise ValueError("Count line must specify at least NEL and NVT.")
    return ints[0], ints[1]


def read_mesh(path: str) -> MeshData:
    with open(path, "r", encoding="utf-8") as handle:
        lines = handle.read().splitlines()

    if len(lines) < 5:
        raise ValueError("Mesh file is too short.")

    header = lines[:2]
    counts_line = lines[2]
    nel, nvt = parse_counts(counts_line)

    index = 3
    if lines[index].strip().upper() != "DCORVG":
        raise ValueError("DCORVG section missing.")
    index += 1

    vertices: List[Tuple[Tuple[float, float, float], str]] = []
    for _ in range(nvt):
        if index >= len(lines):
            raise ValueError("Unexpected end of file in DCORVG section.")
        line = lines[index]
        parts = line.split()
        if len(parts) < 3:
            raise ValueError("Malformed vertex line encountered.")
        coords = (float(parts[0]), float(parts[1]), float(parts[2]))
        vertices.append((coords, line))
        index += 1

    if index >= len(lines) or lines[index].strip().upper() != "KVERT":
        raise ValueError("KVERT section missing.")
    index += 1

    connectivity: List[List[int]] = []
    while index < len(lines):
        stripped = lines[index].strip()
        if stripped == "":
            index += 1
            continue
        if stripped.upper() == "KNPR":
            break
        numbers = [int(value) for value in stripped.split()]
        connectivity.append(numbers)
        index += 1

    if len(connectivity) != nel:
        raise ValueError(f"Expected {nel} connectivity rows, got {len(connectivity)}.")

    if index >= len(lines) or lines[index].strip().upper() != "KNPR":
        raise ValueError("KNPR section missing.")
    index += 1

    nodal_props: List[str] = []
    for _ in range(nvt):
        if index >= len(lines):
            raise ValueError("Unexpected end of file in KNPR section.")
        stripped = lines[index].strip()
        if not stripped:
            raise ValueError("Empty KNPR entry encountered.")
        nodal_props.append(stripped.split()[0])
        index += 1

    tail = lines[index:]

    return MeshData(
        header=header,
        counts_line=counts_line,
        vertices=vertices,
        connectivity=connectivity,
        nodal_props=nodal_props,
        tail=tail,
    )


def find_vertex_indices(
    vertices: Sequence[Tuple[Tuple[float, float, float], str]]
) -> List[int]:
    used = set()
    order: List[int] = []

    for target in TARGET_VERTICES:
        match_index = None
        for idx, (coords, _) in enumerate(vertices):
            if idx in used:
                continue
            if all(abs(c - t) <= COORD_TOL for c, t in zip(coords, target)):
                match_index = idx
                break
        if match_index is None:
            raise ValueError(
                "Could not find vertex with coordinates "
                f"{target} (±{COORD_TOL})."
            )
        used.add(match_index)
        order.append(match_index)

    order.extend(idx for idx in range(len(vertices)) if idx not in used)
    return order


def reorder_mesh(mesh: MeshData) -> MeshData:
    new_order = find_vertex_indices(mesh.vertices)
    old_to_new = {old_idx + 1: new_idx + 1 for new_idx, old_idx in enumerate(new_order)}

    reordered_vertices = [mesh.vertices[idx] for idx in new_order]
    reordered_props = [mesh.nodal_props[idx] for idx in new_order]

    remapped_connectivity: List[List[int]] = []
    for elem in mesh.connectivity:
        remapped = []
        for vertex_id in elem:
            if vertex_id not in old_to_new:
                raise ValueError(f"Vertex index {vertex_id} out of bounds in KVERT.")
            remapped.append(old_to_new[vertex_id])
        remapped_connectivity.append(remapped)

    return MeshData(
        header=mesh.header,
        counts_line=mesh.counts_line,
        vertices=reordered_vertices,
        connectivity=remapped_connectivity,
        nodal_props=reordered_props,
        tail=mesh.tail,
    )


def mesh_to_lines(mesh: MeshData) -> List[str]:
    lines: List[str] = []
    lines.extend(mesh.header)
    lines.append(mesh.counts_line)
    lines.append("DCORVG")
    lines.extend(line for _, line in mesh.vertices)
    lines.append("KVERT")
    lines.extend(" ".join(str(value) for value in row) for row in mesh.connectivity)
    lines.append("KNPR")
    lines.extend(mesh.nodal_props)
    lines.extend(mesh.tail)
    return lines


def write_output(lines: Sequence[str], path: str | None) -> None:
    text = "\n".join(lines) + "\n"
    if path:
        with open(path, "w", encoding="utf-8") as handle:
            handle.write(text)
    else:
        sys.stdout.write(text)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Ensure first 8 vertices follow a canonical ordering."
    )
    parser.add_argument("-i", "--input", required=True, help="Input TRI3D mesh path.")
    parser.add_argument(
        "-o",
        "--output",
        help="Output path. Defaults to stdout if not provided.",
    )
    args = parser.parse_args()

    mesh = read_mesh(args.input)
    reordered = reorder_mesh(mesh)
    lines = mesh_to_lines(reordered)
    write_output(lines, args.output)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:  # pragma: no cover - CLI entry point
        print(f"Error: {exc}", file=sys.stderr)
        sys.exit(1)
