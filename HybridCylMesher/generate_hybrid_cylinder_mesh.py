#!/usr/bin/env python3
"""Structured hexahedral mesh generator for hollow cylinders.

The script reads the unified `setup.e3d` preprocessing configuration file and
produces an output that mimics the `.tri` layout seen in `Mesh.tri`. Only the
fields needed for the provided mesh (geometry + tangential/radial/axial
resolutions) are parsed.
"""

from __future__ import annotations

import argparse
import collections
import configparser
import math
import shutil
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Set, Tuple, Union


MM_TO_CM = 0.1


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate a hollow-cylinder hexahedral mesh in TRI format."
    )
    parser.add_argument(
        "config_path",
        nargs="?",
        type=Path,
        default=None,
        help=(
            "Path to setup.e3d (Ini-style) configuration file. "
            "Defaults to ./setup.e3d."
        ),
    )
    parser.add_argument(
        "-i",
        "--input",
        dest="input_path",
        type=Path,
        help="Alternative way to specify the setup.e3d file.",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=None,
        help=(
            "Where to write the generated TRI mesh. Defaults to "
            "<input-folder>/meshDir/Mesh.tri."
        ),
    )
    parser.add_argument(
        "--theta-offset",
        type=float,
        default=None,
        help=(
            "Optional angular offset (radians) applied to the first "
            "tangential node. Defaults to pi/2 - dtheta to mimic the "
            "provided mesh."
        ),
    )
    parser.add_argument(
        "-s",
        "--grid-size",
        type=float,
        default=None,
        help=(
            "Target element size used to back-fill missing nEl_* values "
            "from the configuration file."
        ),
    )
    return parser.parse_args()


@dataclass
class AnnularMeshInput:
    outer_diameter: float
    inner_diameter: float
    length: float
    start_z: float
    n_t: int | None
    n_r: int | None
    n_z: int | None
    scale_t: float
    scale_r: float
    scale_z: float
    grid_size: float | None


@dataclass
class BoxMeshInput:
    start: Tuple[float, float, float]
    lengths: Tuple[float, float, float]
    n_x: int | None
    n_y: int | None
    n_z: int | None
    scale_x: float
    scale_y: float
    scale_z: float
    grid_size: float | None


@dataclass
class FullCylinderMeshInput:
    outer_diameter: float
    length: float
    start_z: float
    periodicity: int
    n_inner: int | None
    n_outer: int | None
    n_z: int | None
    scale_t: float
    scale_r: float
    scale_z: float
    grid_size: float | None


@dataclass(frozen=True)
class HybridCylinderSection:
    start_z: float
    end_z: float
    mode: str


@dataclass
class HybridCylinderMeshInput:
    outer_diameter: float
    inner_diameter: float
    length: float
    start_z: float
    periodicity: int
    n_inner: int | None
    n_outer: int | None
    n_z: int | None
    scale_t: float
    scale_r: float
    scale_z: float
    grid_size: float | None
    sections: List[HybridCylinderSection]


MeshInput = Union[
    AnnularMeshInput,
    BoxMeshInput,
    FullCylinderMeshInput,
    HybridCylinderMeshInput,
]


def read_parameters(
    path: Path,
) -> MeshInput:
    cfg = configparser.ConfigParser()
    if not cfg.read(path):
        raise FileNotFoundError(f"Unable to read configuration file {path!s}")

    if not cfg.has_section("E3DGeometryData/Preprocessing"):
        raise KeyError("Missing [E3DGeometryData/Preprocessing] section.")

    preprocess = cfg["E3DGeometryData/Preprocessing"]
    geom = preprocess
    sim = preprocess
    grid_size = _get_float_with_default(sim, "MinGap", default=None)
    hex_mesher = sim.get("HexMesher", "Axi").strip().lower()

    if hex_mesher == "box":
        start = _parse_float_list(geom["geometryStart"], expected=3)
        lengths = _parse_float_list(geom["geometryLength"], expected=3)
        n_x = _get_optional_int(sim, "nEl_x")
        n_y = _get_optional_int(sim, "nEl_y")
        n_z = _get_optional_int(sim, "nEl_z")
        scale_x = _get_float_with_default(sim, "sEl_x", default=1.0)
        scale_y = _get_float_with_default(sim, "sEl_y", default=1.0)
        scale_z = _get_float_with_default(sim, "sEl_z", default=1.0)
        return BoxMeshInput(
            start, lengths, n_x, n_y, n_z, scale_x, scale_y, scale_z, grid_size
        )
    if hex_mesher == "fullcylinder":
        outer_diameter = float(geom["BarrelDiameter"])
        length = float(geom["BarrelLength"])
        start_z = float(geom.get("AxialStartPosition", 0.0))
        periodicity = _get_optional_int(sim, "FullCylinderPeriodicity")
        if periodicity is None or periodicity <= 0:
            periodicity = 4
        n_inner = _get_optional_int(sim, "nEl_Tangential")
        n_outer = _get_optional_int(sim, "nEl_Radial")
        n_z = _get_optional_int(sim, "nEl_Axial")
        scale_t = _get_float_with_default(sim, "sEl_Tangential", default=1.0)
        scale_r = _get_float_with_default(sim, "sEl_Radial", default=1.0)
        scale_z = _get_float_with_default(sim, "sEl_Axial", default=1.0)
        return FullCylinderMeshInput(
            outer_diameter,
            length,
            start_z,
            periodicity,
            n_inner,
            n_outer,
            n_z,
            scale_t,
            scale_r,
            scale_z,
            grid_size,
        )
    if hex_mesher == "hybridcylinder":
        outer_diameter = float(geom["BarrelDiameter"])
        inner_diameter = float(geom["InnerDiameter"])
        length = float(geom["BarrelLength"])
        start_z = float(geom.get("AxialStartPosition", 0.0))
        periodicity = _get_optional_int(sim, "FullCylinderPeriodicity")
        if periodicity is None or periodicity <= 0:
            periodicity = 4
        n_inner = _get_optional_int(sim, "nEl_Tangential")
        n_outer = _get_optional_int(sim, "nEl_Radial")
        n_z = _get_optional_int(sim, "nEl_Axial")
        scale_t = _get_float_with_default(sim, "sEl_Tangential", default=1.0)
        scale_r = _get_float_with_default(sim, "sEl_Radial", default=1.0)
        scale_z = _get_float_with_default(sim, "sEl_Axial", default=1.0)
        raw_sections = sim.get("HybridCylinderSections", fallback="").strip()
        if not raw_sections:
            raise ValueError(
                "HybridCylinder requires HybridCylinderSections in "
                "[E3DGeometryData/Preprocessing]."
            )
        sections = _parse_hybrid_sections(
            raw_sections, start_z=start_z, length=length
        )
        return HybridCylinderMeshInput(
            outer_diameter,
            inner_diameter,
            length,
            start_z,
            periodicity,
            n_inner,
            n_outer,
            n_z,
            scale_t,
            scale_r,
            scale_z,
            grid_size,
            sections,
        )

    outer_diameter = float(geom["BarrelDiameter"])
    inner_diameter = float(geom["InnerDiameter"])
    length = float(geom["BarrelLength"])
    start_z = float(geom.get("AxialStartPosition", 0.0))
    n_t = _get_optional_int(sim, "nEl_Tangential")
    n_r = _get_optional_int(sim, "nEl_Radial")
    n_z = _get_optional_int(sim, "nEl_Axial")
    scale_t = _get_float_with_default(sim, "sEl_Tangential", default=1.0)
    scale_r = _get_float_with_default(sim, "sEl_Radial", default=1.0)
    scale_z = _get_float_with_default(sim, "sEl_Axial", default=1.0)
    return AnnularMeshInput(
        outer_diameter,
        inner_diameter,
        length,
        start_z,
        n_t,
        n_r,
        n_z,
        scale_t,
        scale_r,
        scale_z,
        grid_size,
    )


def _get_optional_int(section: configparser.SectionProxy, key: str) -> int | None:
    raw_value = section.get(key, fallback=None)
    if raw_value is None or not raw_value.strip():
        return None
    return int(raw_value)


def _parse_float_list(value: str, *, expected: int) -> Tuple[float, ...]:
    parts = [part.strip() for part in value.split(",")]
    if len(parts) != expected:
        raise ValueError(
            f"Expected {expected} comma-separated values, got {len(parts)}."
        )
    return tuple(float(part) for part in parts)


def _get_float_with_default(
    section: configparser.SectionProxy, key: str, *, default: float | None
) -> float | None:
    raw_value = section.get(key, fallback=None)
    if raw_value is None or not raw_value.strip():
        return default
    return float(raw_value)


def _parse_hybrid_sections(
    raw_value: str, *, start_z: float, length: float
) -> List[HybridCylinderSection]:
    end_z = start_z + length
    parsed_sections: List[HybridCylinderSection] = []

    for raw_entry in raw_value.split(","):
        entry = raw_entry.strip()
        if not entry:
            continue
        parts = [part.strip() for part in entry.split(":")]
        if len(parts) != 3:
            raise ValueError(
                "Each HybridCylinderSections entry must have the form "
                "'start:end:mode'."
            )
        sec_start = float(parts[0])
        sec_end = float(parts[1])
        mode = parts[2].lower()
        if mode not in {"full", "hollow"}:
            raise ValueError(
                "HybridCylinder section mode must be 'full' or 'hollow'."
            )
        if sec_end <= sec_start:
            raise ValueError(
                "HybridCylinder section end must be greater than start."
            )
        parsed_sections.append(HybridCylinderSection(sec_start, sec_end, mode))

    if not parsed_sections:
        raise ValueError("HybridCylinderSections must not be empty.")

    parsed_sections.sort(key=lambda section: section.start_z)
    tolerance = max(1e-9, length * 1e-9)
    cursor = start_z

    for section in parsed_sections:
        if abs(section.start_z - cursor) > tolerance:
            raise ValueError(
                "HybridCylinderSections must form a contiguous coverage from "
                f"{start_z} to {end_z}."
            )
        cursor = section.end_z

    if abs(cursor - end_z) > tolerance:
        raise ValueError(
            "HybridCylinderSections must end at BarrelLength + "
            "AxialStartPosition."
        )

    return parsed_sections


def resolve_annular_resolutions(
    outer_diameter: float,
    inner_diameter: float,
    length: float,
    n_t: int | None,
    n_r: int | None,
    n_z: int | None,
    grid_size: float | None,
    scale_t: float,
    scale_r: float,
    scale_z: float,
) -> Tuple[int, int, int]:
    """Return tangential/radial/axial counts, estimating missing ones."""
    for label, scale in (
        ("sEl_Tangential", scale_t),
        ("sEl_Radial", scale_r),
        ("sEl_Axial", scale_z),
    ):
        if scale <= 0.0:
            raise ValueError(f"{label} must be positive.")

    missing = [n_t is None, n_r is None, n_z is None]
    if grid_size is None:
        if any(missing):
            raise ValueError(
                "nEl_Tangential, nEl_Radial, and nEl_Axial must be provided "
                "when --grid-size is not specified."
            )
        assert n_t is not None and n_r is not None and n_z is not None
        return n_t, n_r, n_z

    if grid_size <= 0.0:
        raise ValueError("--grid-size must be a positive number.")

    outer_radius = outer_diameter / 2.0
    inner_radius = inner_diameter / 2.0
    if outer_radius <= inner_radius:
        raise ValueError("Outer diameter must be larger than inner diameter.")

    radial_span = outer_radius - inner_radius
    if radial_span <= 0.0:
        raise ValueError("Invalid geometry: zero radial thickness.")

    avg_radius = 0.5 * (outer_radius + inner_radius)
    circumference = 2.0 * math.pi * avg_radius

    n_t_final = (
        n_t
        if n_t is not None
        else _estimate_elements(
            circumference, grid_size * scale_t, minimum=3
        )
    )
    n_r_final = (
        n_r
        if n_r is not None
        else _estimate_elements(radial_span, grid_size * scale_r, minimum=1)
    )
    n_z_final = (
        n_z
        if n_z is not None
        else _estimate_elements(length, grid_size * scale_z, minimum=1)
    )

    for label, value in (
        ("nEl_Tangential", n_t_final),
        ("nEl_Radial", n_r_final),
        ("nEl_Axial", n_z_final),
    ):
        if value <= 0:
            raise ValueError(f"{label} must be positive.")

    return n_t_final, n_r_final, n_z_final


def _estimate_elements(
    span: float, grid_size: float, *, minimum: int
) -> int:
    if span <= 0.0:
        raise ValueError("Cannot estimate elements for non-positive span.")
    return max(minimum, int(math.ceil(span / grid_size)))


def resolve_box_resolutions(
    lengths: Tuple[float, float, float],
    n_x: int | None,
    n_y: int | None,
    n_z: int | None,
    grid_size: float | None,
    scale_x: float,
    scale_y: float,
    scale_z: float,
) -> Tuple[int, int, int]:
    if any(length <= 0.0 for length in lengths):
        raise ValueError("geometryLength values must be positive.")
    for label, scale in (
        ("sEl_x", scale_x),
        ("sEl_y", scale_y),
        ("sEl_z", scale_z),
    ):
        if scale <= 0.0:
            raise ValueError(f"{label} must be positive.")
    missing = [n_x is None, n_y is None, n_z is None]
    if grid_size is None:
        if any(missing):
            raise ValueError(
                "nEl_x, nEl_y, and nEl_z must be provided when --grid-size is not specified."
            )
        assert n_x is not None and n_y is not None and n_z is not None
    else:
        if grid_size <= 0.0:
            raise ValueError("--grid-size must be a positive number.")
        n_x = (
            n_x
            if n_x is not None
            else _estimate_elements(lengths[0], grid_size * scale_x, minimum=1)
        )
        n_y = (
            n_y
            if n_y is not None
            else _estimate_elements(lengths[1], grid_size * scale_y, minimum=1)
        )
        n_z = (
            n_z
            if n_z is not None
            else _estimate_elements(lengths[2], grid_size * scale_z, minimum=1)
        )
    assert n_x is not None and n_y is not None and n_z is not None
    for label, value in (("nEl_x", n_x), ("nEl_y", n_y), ("nEl_z", n_z)):
        if value <= 0:
            raise ValueError(f"{label} must be positive.")
    return n_x, n_y, n_z


def build_annular_coordinates(
    outer_diameter: float,
    inner_diameter: float,
    length: float,
    start_z: float,
    n_t: int,
    n_r: int,
    n_z: int,
    theta_offset: float | None,
) -> Tuple[List[Tuple[float, float, float]], List[int]]:
    outer_radius = outer_diameter / 2.0
    inner_radius = inner_diameter / 2.0
    if outer_radius <= inner_radius:
        raise ValueError("Outer diameter must be larger than inner diameter.")

    n_r_nodes = n_r + 1
    n_z_nodes = n_z + 1
    tangential_step = 2.0 * math.pi / n_t
    if theta_offset is None:
        theta_offset = math.pi / 2.0 - tangential_step

    radii = [
        inner_radius + (outer_radius - inner_radius) * j / (n_r_nodes - 1)
        for j in range(n_r_nodes)
    ]
    if n_z_nodes > 1:
        z_levels = [
            start_z + length * k / (n_z_nodes - 1) for k in range(n_z_nodes)
        ]
    else:
        z_levels = [start_z]

    coords: List[Tuple[float, float, float]] = []
    knpr: List[int] = []

    for k, z in enumerate(z_levels):
        for j, radius in enumerate(radii):
            for i in range(n_t):
                theta = theta_offset + i * tangential_step
                x = radius * math.cos(theta)
                y = radius * math.sin(theta)
                coords.append((x, y, z))

                is_boundary = (
                    k == 0
                    or k == n_z_nodes - 1
                    or j == 0
                    or j == n_r_nodes - 1
                )
                knpr.append(1 if is_boundary else 0)

    return coords, knpr


def annular_node_index(
    t: int, r: int, z: int, n_t: int, n_r_nodes: int
) -> int:
    """Return 1-based vertex index for (t, r, z)."""
    return (
        t
        + n_t * (r + n_r_nodes * z)  # flatten (r, z)
        + 1
    )


def build_annular_connectivity(
    n_t: int, n_r: int, n_z: int
) -> Iterable[Tuple[int, int, int, int, int, int, int, int]]:
    n_r_nodes = n_r + 1
    for z in range(n_z):
        for r in range(n_r):
            for t in range(n_t):
                t_next = (t + 1) % n_t
                bottom_left = annular_node_index(t, r, z, n_t, n_r_nodes)
                bottom_right = annular_node_index(
                    t_next, r, z, n_t, n_r_nodes
                )
                top_right = annular_node_index(
                    t_next, r + 1, z, n_t, n_r_nodes
                )
                top_left = annular_node_index(
                    t, r + 1, z, n_t, n_r_nodes
                )

                bottom_left_up = annular_node_index(
                    t, r, z + 1, n_t, n_r_nodes
                )
                bottom_right_up = annular_node_index(
                    t_next, r, z + 1, n_t, n_r_nodes
                )
                top_right_up = annular_node_index(
                    t_next, r + 1, z + 1, n_t, n_r_nodes
                )
                top_left_up = annular_node_index(
                    t, r + 1, z + 1, n_t, n_r_nodes
                )

                yield (
                    bottom_left,
                    bottom_right,
                    top_right,
                    top_left,
                    bottom_left_up,
                    bottom_right_up,
                    top_right_up,
                top_left_up,
            )


def resolve_full_cylinder_resolutions(
    outer_diameter: float,
    length: float,
    periodicity: int,
    n_inner: int | None,
    n_outer: int | None,
    n_z: int | None,
    grid_size: float | None,
    scale_t: float,
    scale_r: float,
    scale_z: float,
) -> Tuple[int, int, int]:
    if outer_diameter <= 0.0:
        raise ValueError("BarrelDiameter must be positive for FullCylinder.")
    if length <= 0.0:
        raise ValueError("BarrelLength must be positive for FullCylinder.")
    if periodicity <= 0:
        raise ValueError("FullCylinderPeriodicity must be positive.")
    for label, scale in (
        ("sEl_Tangential", scale_t),
        ("sEl_Radial", scale_r),
        ("sEl_Axial", scale_z),
    ):
        if scale <= 0.0:
            raise ValueError(f"{label} must be positive.")
    if grid_size is None:
        if n_inner is None or n_outer is None or n_z is None:
            raise ValueError(
                "nEl_Tangential, nEl_Radial, and nEl_Axial must be provided "
                "for FullCylinder meshes when --grid-size is not specified."
            )
        n_inner_final = n_inner
        n_outer_final = n_outer
        n_z_final = n_z
    else:
        if grid_size <= 0.0:
            raise ValueError("--grid-size must be a positive number.")
        outer_radius = outer_diameter / 2.0
        circumference = 2.0 * math.pi * outer_radius
        tangential_target = _estimate_elements(
            circumference, grid_size * scale_t, minimum=2 * periodicity
        )
        if n_inner is None:
            n_inner_final = max(
                1, math.ceil(tangential_target / (2 * periodicity))
            )
        else:
            n_inner_final = n_inner
        radial_target_total = _estimate_elements(
            outer_radius, grid_size * scale_r, minimum=1
        )
        if n_outer is None:
            n_outer_final = max(1, radial_target_total - n_inner_final)
        else:
            n_outer_final = n_outer
        if n_z is None:
            n_z_final = _estimate_elements(
                length, grid_size * scale_z, minimum=1
            )
        else:
            n_z_final = n_z
    for label, value in (
        ("nEl_Tangential", n_inner_final),
        ("nEl_Radial", n_outer_final),
        ("nEl_Axial", n_z_final),
    ):
        if value <= 0:
            raise ValueError(f"{label} must be positive for FullCylinder.")
    return n_inner_final, n_outer_final, n_z_final


def build_full_cylinder_coordinates(
    outer_diameter: float,
    length: float,
    start_z: float,
    periodicity: int,
    n_inner: int,
    n_outer: int,
    n_z: int,
) -> Tuple[List[Tuple[float, float, float]], List[int]]:
    if n_inner <= 0 or n_outer <= 0 or n_z <= 0:
        raise ValueError("FullCylinder element counts must be positive.")
    if periodicity <= 0:
        raise ValueError("FullCylinderPeriodicity must be positive.")
    outer_radius = outer_diameter / 2.0
    n_z_nodes = n_z + 1
    coords: List[Tuple[float, float, float]] = []
    knpr: List[int] = []

    inner_fraction = n_inner / (n_inner + n_outer)
    inner_radius_span = outer_radius * inner_fraction
    outer_radius_span = max(0.0, outer_radius - inner_radius_span)

    def add_node(x: float, y: float, z: float, is_boundary: bool) -> None:
        coords.append((x, y, z))
        knpr.append(1 if is_boundary else 0)

    if n_z == 0:
        z_levels = [start_z]
    else:
        z_levels = [
            start_z + length * k / n_z for k in range(n_z_nodes)
        ]

    for k, z in enumerate(z_levels):
        at_axial_boundary = (k == 0) or (k == n_z_nodes - 1)
        add_node(0.0, 0.0, z, at_axial_boundary)

        for j in range(1, n_inner + 1):
            ni = 2 * periodicity * j
            radius = inner_radius_span * j / n_inner
            for i in range(ni):
                theta = 2.0 * math.pi * i / ni
                x = radius * math.cos(theta)
                y = radius * math.sin(theta)
                is_boundary = at_axial_boundary or (
                    n_outer == 0 and j == n_inner
                )
                add_node(x, y, z, is_boundary)

        if n_outer <= 0:
            continue
        outer_nodes = 2 * periodicity * n_inner
        for j in range(1, n_outer + 1):
            radius = inner_radius_span
            radius += outer_radius_span * j / n_outer
            for i in range(outer_nodes):
                theta = 2.0 * math.pi * i / outer_nodes
                x = radius * math.cos(theta)
                y = radius * math.sin(theta)
                is_boundary = at_axial_boundary or (j == n_outer)
                add_node(x, y, z, is_boundary)

    return coords, knpr


def build_hybrid_cylinder_coordinates(
    outer_diameter: float,
    inner_diameter: float,
    length: float,
    start_z: float,
    periodicity: int,
    n_inner: int,
    n_outer: int,
    n_z: int,
) -> Tuple[List[Tuple[float, float, float]], List[int]]:
    if n_inner <= 0 or n_outer <= 0 or n_z <= 0:
        raise ValueError("HybridCylinder element counts must be positive.")
    if periodicity <= 0:
        raise ValueError("FullCylinderPeriodicity must be positive.")

    outer_radius = outer_diameter / 2.0
    inner_radius = inner_diameter / 2.0
    if inner_radius <= 0.0:
        raise ValueError("InnerDiameter must be positive for HybridCylinder.")
    if inner_radius >= outer_radius:
        raise ValueError(
            "InnerDiameter must be smaller than BarrelDiameter for HybridCylinder."
        )

    n_z_nodes = n_z + 1
    coords: List[Tuple[float, float, float]] = []
    knpr: List[int] = []

    def add_node(x: float, y: float, z: float, is_boundary: bool) -> None:
        coords.append((x, y, z))
        knpr.append(1 if is_boundary else 0)

    z_levels = [start_z + length * k / n_z for k in range(n_z_nodes)]

    for k, z in enumerate(z_levels):
        at_axial_boundary = (k == 0) or (k == n_z_nodes - 1)
        add_node(0.0, 0.0, z, at_axial_boundary)

        for j in range(1, n_inner + 1):
            node_count = 2 * periodicity * j
            radius = inner_radius * j / n_inner
            for i in range(node_count):
                theta = 2.0 * math.pi * i / node_count
                x = radius * math.cos(theta)
                y = radius * math.sin(theta)
                add_node(x, y, z, at_axial_boundary)

        shell_node_count = 2 * periodicity * n_inner
        shell_span = outer_radius - inner_radius
        for j in range(1, n_outer + 1):
            radius = inner_radius + shell_span * j / n_outer
            for i in range(shell_node_count):
                theta = 2.0 * math.pi * i / shell_node_count
                x = radius * math.cos(theta)
                y = radius * math.sin(theta)
                is_boundary = at_axial_boundary or (j == n_outer)
                add_node(x, y, z, is_boundary)

    return coords, knpr


def build_full_cylinder_connectivity(
    periodicity: int, n_inner: int, n_outer: int, n_z: int
) -> Iterable[Tuple[int, int, int, int, int, int, int, int]]:
    if n_inner <= 0 or n_outer <= 0 or n_z <= 0:
        raise ValueError("FullCylinder element counts must be positive.")
    if periodicity <= 0:
        raise ValueError("FullCylinderPeriodicity must be positive.")
    level_nodes = (
        periodicity * (n_inner * n_inner + n_inner)
        + n_outer * 2 * n_inner * periodicity
        + 1
    )

    for k in range(1, n_z + 1):
        mm1 = (k - 1) * level_nodes
        mm2 = k * level_nodes
        kk = periodicity * 2 + 1
        ll = 1

        for j in range(1, periodicity):
            i = 2 * (j - 1) + 1
            yield (
                mm1 + i + 1,
                mm1 + i + 2,
                mm1 + i + 3,
                mm1 + 1,
                mm2 + i + 1,
                mm2 + i + 2,
                mm2 + i + 3,
                mm2 + 1,
            )
        i = 2 * (periodicity - 1) + 1
        yield (
            mm1 + i + 1,
            mm1 + i + 2,
            mm1 + 2,
            mm1 + 1,
            mm2 + i + 1,
            mm2 + i + 2,
            mm2 + 2,
            mm2 + 1,
        )

        for j in range(2, n_inner + 1):
            for i in range(1, 2 * periodicity * j + 1):
                mod_value = i % (2 * j)
                if mod_value != j and mod_value != j + 1:
                    kk += 1
                    ll += 1
                    if i < 2 * periodicity * j:
                        yield (
                            mm1 + kk,
                            mm1 + kk + 1,
                            mm1 + ll + 1,
                            mm1 + ll,
                            mm2 + kk,
                            mm2 + kk + 1,
                            mm2 + ll + 1,
                            mm2 + ll,
                        )
                    else:
                        yield (
                            mm1 + kk,
                            mm1 + ll + 1,
                            mm1 + ll + 1 - 2 * periodicity * (j - 1),
                            mm1 + ll,
                            mm2 + kk,
                            mm2 + ll + 1,
                            mm2 + ll + 1 - 2 * periodicity * (j - 1),
                            mm2 + ll,
                        )
                else:
                    if mod_value == j:
                        kk += 1
                        yield (
                            mm1 + kk,
                            mm1 + kk + 1,
                            mm1 + kk + 2,
                            mm1 + ll + 1,
                            mm2 + kk,
                            mm2 + kk + 1,
                            mm2 + kk + 2,
                            mm2 + ll + 1,
                        )
                    else:
                        kk += 1

        for j in range(1, n_outer + 1):
            for _ in range(1, 2 * periodicity * n_inner):
                kk += 1
                ll += 1
                yield (
                    mm1 + kk,
                    mm1 + kk + 1,
                    mm1 + ll + 1,
                    mm1 + ll,
                    mm2 + kk,
                    mm2 + kk + 1,
                    mm2 + ll + 1,
                    mm2 + ll,
                )
            kk += 1
            ll += 1
            yield (
                mm1 + kk,
                mm1 + ll + 1,
                mm1 + ll + 1 - 2 * periodicity * n_inner,
                mm1 + ll,
                mm2 + kk,
                mm2 + ll + 1,
                mm2 + ll + 1 - 2 * periodicity * n_inner,
                mm2 + ll,
            )


def snap_hybrid_sections_to_layers(
    sections: Sequence[HybridCylinderSection],
    *,
    start_z: float,
    length: float,
    n_z: int,
) -> List[Tuple[int, int, str]]:
    if n_z <= 0:
        raise ValueError("HybridCylinder requires a positive axial resolution.")

    dz = length / n_z
    snapped_sections: List[Tuple[int, int, str]] = []

    for section in sections:
        start_index = round((section.start_z - start_z) / dz)
        end_index = round((section.end_z - start_z) / dz)
        start_index = max(0, min(n_z, start_index))
        end_index = max(0, min(n_z, end_index))
        if end_index <= start_index:
            raise ValueError(
                "A snapped HybridCylinder section collapsed to zero axial thickness. "
                "Adjust nEl_Axial or the section boundaries."
            )
        snapped_sections.append((start_index, end_index, section.mode))

    if snapped_sections[0][0] != 0 or snapped_sections[-1][1] != n_z:
        raise ValueError(
            "HybridCylinder sections do not cover the full axial range after snapping."
        )

    for previous, current in zip(snapped_sections, snapped_sections[1:]):
        if previous[1] != current[0]:
            raise ValueError(
                "HybridCylinder sections became non-contiguous after snapping. "
                "Adjust nEl_Axial or the section boundaries."
            )

    return snapped_sections


def build_hybrid_cylinder_connectivity(
    periodicity: int,
    n_inner: int,
    n_outer: int,
    n_z: int,
    hollow_layers: set[int],
) -> Iterable[Tuple[int, int, int, int, int, int, int, int]]:
    if n_inner <= 0 or n_outer <= 0 or n_z <= 0:
        raise ValueError("HybridCylinder element counts must be positive.")
    if periodicity <= 0:
        raise ValueError("FullCylinderPeriodicity must be positive.")

    level_nodes = (
        periodicity * (n_inner * n_inner + n_inner)
        + n_outer * 2 * n_inner * periodicity
        + 1
    )

    for k in range(1, n_z + 1):
        mm1 = (k - 1) * level_nodes
        mm2 = k * level_nodes
        kk = periodicity * 2 + 1
        ll = 1
        layer_is_hollow = (k - 1) in hollow_layers

        if not layer_is_hollow:
            for j in range(1, periodicity):
                i = 2 * (j - 1) + 1
                yield (
                    mm1 + i + 1,
                    mm1 + i + 2,
                    mm1 + i + 3,
                    mm1 + 1,
                    mm2 + i + 1,
                    mm2 + i + 2,
                    mm2 + i + 3,
                    mm2 + 1,
                )
            i = 2 * (periodicity - 1) + 1
            yield (
                mm1 + i + 1,
                mm1 + i + 2,
                mm1 + 2,
                mm1 + 1,
                mm2 + i + 1,
                mm2 + i + 2,
                mm2 + 2,
                mm2 + 1,
            )

            for j in range(2, n_inner + 1):
                for i in range(1, 2 * periodicity * j + 1):
                    mod_value = i % (2 * j)
                    if mod_value != j and mod_value != j + 1:
                        kk += 1
                        ll += 1
                        if i < 2 * periodicity * j:
                            yield (
                                mm1 + kk,
                                mm1 + kk + 1,
                                mm1 + ll + 1,
                                mm1 + ll,
                                mm2 + kk,
                                mm2 + kk + 1,
                                mm2 + ll + 1,
                                mm2 + ll,
                            )
                        else:
                            yield (
                                mm1 + kk,
                                mm1 + ll + 1,
                                mm1 + ll + 1
                                - 2 * periodicity * (j - 1),
                                mm1 + ll,
                                mm2 + kk,
                                mm2 + ll + 1,
                                mm2 + ll + 1
                                - 2 * periodicity * (j - 1),
                                mm2 + ll,
                            )
                    else:
                        if mod_value == j:
                            kk += 1
                            yield (
                                mm1 + kk,
                                mm1 + kk + 1,
                                mm1 + kk + 2,
                                mm1 + ll + 1,
                                mm2 + kk,
                                mm2 + kk + 1,
                                mm2 + kk + 2,
                                mm2 + ll + 1,
                            )
                        else:
                            kk += 1
        else:
            kk = periodicity * (n_inner * n_inner + n_inner) + 1
            ll = periodicity * (n_inner * n_inner - n_inner) + 1

        for j in range(1, n_outer + 1):
            for _ in range(1, 2 * periodicity * n_inner):
                kk += 1
                ll += 1
                yield (
                    mm1 + kk,
                    mm1 + kk + 1,
                    mm1 + ll + 1,
                    mm1 + ll,
                    mm2 + kk,
                    mm2 + kk + 1,
                    mm2 + ll + 1,
                    mm2 + ll,
                )
            kk += 1
            ll += 1
            yield (
                mm1 + kk,
                mm1 + ll + 1,
                mm1 + ll + 1 - 2 * periodicity * n_inner,
                mm1 + ll,
                mm2 + kk,
                mm2 + ll + 1,
                mm2 + ll + 1 - 2 * periodicity * n_inner,
                mm2 + ll,
            )


def mark_hybrid_boundaries(
    knpr: List[int],
    sections: Sequence[Tuple[int, int, str]],
    *,
    periodicity: int,
    n_inner: int,
    n_outer: int,
    n_z: int,
) -> None:
    level_nodes = (
        periodicity * (n_inner * n_inner + n_inner)
        + n_outer * 2 * n_inner * periodicity
        + 1
    )
    inner_ring_nodes = 2 * periodicity * n_inner
    inner_ring_offset = 2 + periodicity * (n_inner * n_inner - n_inner)

    hollow_layers: set[int] = set()
    for start_idx, end_idx, mode in sections:
        if mode == "hollow":
            hollow_layers.update(range(start_idx, end_idx))

    for z_level in range(n_z + 1):
        touches_hollow = (
            z_level < n_z and z_level in hollow_layers
        ) or (z_level > 0 and (z_level - 1) in hollow_layers)
        if not touches_hollow:
            continue
        base = z_level * level_nodes + inner_ring_offset
        for offset in range(inner_ring_nodes):
            knpr[base + offset - 1] = 1


def remove_unused_nodes(
    coords: Sequence[Tuple[float, float, float]],
    cells: Sequence[Tuple[int, int, int, int, int, int, int, int]],
    knpr: Sequence[int],
) -> Tuple[
    List[Tuple[float, float, float]],
    List[Tuple[int, int, int, int, int, int, int, int]],
    List[int],
]:
    used_nodes = sorted({node for cell in cells for node in cell})
    index_map = {old_index: new_index for new_index, old_index in enumerate(used_nodes, start=1)}

    pruned_coords = [coords[node - 1] for node in used_nodes]
    pruned_knpr = [knpr[node - 1] for node in used_nodes]
    pruned_cells = [
        tuple(index_map[node] for node in cell)
        for cell in cells
    ]
    return pruned_coords, pruned_cells, pruned_knpr


def build_box_coordinates(
    start: Tuple[float, float, float],
    lengths: Tuple[float, float, float],
    n_x: int,
    n_y: int,
    n_z: int,
) -> Tuple[List[Tuple[float, float, float]], List[int]]:
    for label, value in (("nEl_x", n_x), ("nEl_y", n_y), ("nEl_z", n_z)):
        if value <= 0:
            raise ValueError(f"{label} must be positive.")
    if any(length <= 0.0 for length in lengths):
        raise ValueError("geometryLength values must be positive.")

    x0, y0, z0 = start
    lx, ly, lz = lengths
    x_nodes = n_x + 1
    y_nodes = n_y + 1
    z_nodes = n_z + 1

    x_coords = [x0 + lx * i / n_x for i in range(x_nodes)]
    y_coords = [y0 + ly * j / n_y for j in range(y_nodes)]
    z_coords = [z0 + lz * k / n_z for k in range(z_nodes)]

    coords: List[Tuple[float, float, float]] = []
    knpr: List[int] = []

    for k, z in enumerate(z_coords):
        for j, y in enumerate(y_coords):
            for i, x in enumerate(x_coords):
                coords.append((x, y, z))
                is_boundary = (
                    i == 0
                    or i == x_nodes - 1
                    or j == 0
                    or j == y_nodes - 1
                    or k == 0
                    or k == z_nodes - 1
                )
                knpr.append(1 if is_boundary else 0)

    return coords, knpr


def box_node_index(
    x: int, y: int, z: int, n_x_nodes: int, n_y_nodes: int
) -> int:
    return x + n_x_nodes * (y + n_y_nodes * z) + 1


def build_box_connectivity(
    n_x: int, n_y: int, n_z: int
) -> Iterable[Tuple[int, int, int, int, int, int, int, int]]:
    n_x_nodes = n_x + 1
    n_y_nodes = n_y + 1
    for z in range(n_z):
        for y in range(n_y):
            for x in range(n_x):
                bottom_left = box_node_index(x, y, z, n_x_nodes, n_y_nodes)
                bottom_right = box_node_index(
                    x + 1, y, z, n_x_nodes, n_y_nodes
                )
                top_right = box_node_index(
                    x + 1, y + 1, z, n_x_nodes, n_y_nodes
                )
                top_left = box_node_index(x, y + 1, z, n_x_nodes, n_y_nodes)

                bottom_left_up = box_node_index(
                    x, y, z + 1, n_x_nodes, n_y_nodes
                )
                bottom_right_up = box_node_index(
                    x + 1, y, z + 1, n_x_nodes, n_y_nodes
                )
                top_right_up = box_node_index(
                    x + 1, y + 1, z + 1, n_x_nodes, n_y_nodes
                )
                top_left_up = box_node_index(
                    x, y + 1, z + 1, n_x_nodes, n_y_nodes
                )

                yield (
                    bottom_left,
                    bottom_right,
                    top_right,
                    top_left,
                    bottom_left_up,
                    bottom_right_up,
                    top_right_up,
                    top_left_up,
                )


def write_tri(
    coords: Sequence[Tuple[float, float, float]],
    cells: Iterable[Tuple[int, int, int, int, int, int, int, int]],
    knpr: Sequence[int],
    output: Path,
) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    cell_list = list(cells)
    nel = len(cell_list)
    nvt = len(coords)

    with output.open("w", encoding="utf-8") as fh:
        fh.write(" Coarse mesh exported by AxiMesher\n")
        fh.write(" Parametrisierung PARXC, PARYC, TMAXC\n")
        fh.write(f"  {nel:5d}  {nvt:5d} 1 8 12 6     NEL,NVT,NBCT,NVE\n")
        fh.write(" DCORVG\n")
        for x, y, z in coords:
            fh.write(
                f"  {x * MM_TO_CM:24.15f}        "
                f"{y * MM_TO_CM:24.15f}        "
                f"{z * MM_TO_CM:24.15f}\n"
            )

        fh.write(" KVERT\n")
        for cell in cell_list:
            fh.write("".join(f"{node:8d}" for node in cell) + "\n")

        fh.write(" KNPR\n")
        for flag in knpr:
            fh.write(f"  {flag:d}\n")


def hex_faces(
    cell: Tuple[int, int, int, int, int, int, int, int]
) -> List[Tuple[int, int, int, int]]:
    n1, n2, n3, n4, n5, n6, n7, n8 = cell
    return [
        (n1, n2, n3, n4),
        (n5, n6, n7, n8),
        (n1, n2, n6, n5),
        (n2, n3, n7, n6),
        (n3, n4, n8, n7),
        (n4, n1, n5, n8),
    ]


def hex_edges(
    cell: Tuple[int, int, int, int, int, int, int, int]
) -> List[Tuple[int, int]]:
    n1, n2, n3, n4, n5, n6, n7, n8 = cell
    return [
        (n1, n2),
        (n2, n3),
        (n3, n4),
        (n4, n1),
        (n5, n6),
        (n6, n7),
        (n7, n8),
        (n8, n5),
        (n1, n5),
        (n2, n6),
        (n3, n7),
        (n4, n8),
    ]


def compute_boundary_faces(
    cells: Sequence[Tuple[int, int, int, int, int, int, int, int]]
) -> List[Tuple[int, int, int, int]]:
    face_map: Dict[
        Tuple[int, int, int, int], List[Tuple[int, int, int, int]]
    ] = collections.defaultdict(list)

    for cell in cells:
        for face in hex_faces(cell):
            face_map[tuple(sorted(face))].append(face)

    return [
        faces[0]
        for faces in face_map.values()
        if len(faces) == 1
    ]


def compute_shortest_edge_length(
    coords: Sequence[Tuple[float, float, float]],
    cells: Sequence[Tuple[int, int, int, int, int, int, int, int]],
) -> float:
    seen_edges: Set[Tuple[int, int]] = set()
    shortest = math.inf

    for cell in cells:
        for edge in hex_edges(cell):
            normalized = tuple(sorted(edge))
            if normalized in seen_edges:
                continue
            seen_edges.add(normalized)
            p0 = coords[normalized[0] - 1]
            p1 = coords[normalized[1] - 1]
            length = math.dist(p0, p1)
            if length > 0.0:
                shortest = min(shortest, length)

    if not math.isfinite(shortest):
        raise ValueError("Unable to determine a positive hexa edge length.")
    return shortest


def classify_boundary_faces(
    coords: Sequence[Tuple[float, float, float]],
    boundary_faces: Sequence[Tuple[int, int, int, int]],
    tolerance: float,
) -> Dict[str, List[Tuple[int, int, int, int]]]:
    z_values = [coord[2] for coord in coords]
    z_min = min(z_values)
    z_max = max(z_values)
    r_values = [math.hypot(coord[0], coord[1]) for coord in coords]
    r_max = max(r_values)

    classified: Dict[str, List[Tuple[int, int, int, int]]] = {
        "inflow": [],
        "outflow": [],
        "barrel": [],
        "inner_transition": [],
        "inner_cylinder": [],
        "unclassified": [],
    }

    for face in boundary_faces:
        points = [coords[node - 1] for node in face]
        z_face = [point[2] for point in points]
        r_face = [math.hypot(point[0], point[1]) for point in points]
        z_span = max(z_face) - min(z_face)
        r_span = max(r_face) - min(r_face)

        if max(abs(z - z_min) for z in z_face) <= tolerance:
            classified["inflow"].append(face)
        elif max(abs(z - z_max) for z in z_face) <= tolerance:
            classified["outflow"].append(face)
        elif r_span <= tolerance and abs(sum(r_face) / 4.0 - r_max) <= tolerance:
            classified["barrel"].append(face)
        elif z_span <= tolerance:
            classified["inner_transition"].append(face)
        elif r_span <= tolerance:
            classified["inner_cylinder"].append(face)
        else:
            classified["unclassified"].append(face)

    return classified


def validate_boundary_classification(
    boundary_faces: Sequence[Tuple[int, int, int, int]],
    classified: Dict[str, List[Tuple[int, int, int, int]]],
) -> None:
    classified_faces = [
        face
        for name, faces in classified.items()
        if name != "unclassified"
        for face in faces
    ]

    boundary_keys = {tuple(sorted(face)) for face in boundary_faces}
    classified_keys = {tuple(sorted(face)) for face in classified_faces}
    unclassified_keys = {
        tuple(sorted(face)) for face in classified["unclassified"]
    }

    if unclassified_keys:
        raise ValueError(
            f"Boundary classification left {len(unclassified_keys)} faces "
            "unassigned."
        )
    if classified_keys != boundary_keys:
        missing = len(boundary_keys - classified_keys)
        extra = len(classified_keys - boundary_keys)
        raise ValueError(
            "Boundary classification mismatch detected: "
            f"missing={missing}, extra={extra}."
        )


def split_faces_by_connected_component(
    faces: Sequence[Tuple[int, int, int, int]]
) -> List[List[Tuple[int, int, int, int]]]:
    if not faces:
        return []

    node_to_faces: Dict[int, List[int]] = collections.defaultdict(list)
    for face_index, face in enumerate(faces):
        for node in face:
            node_to_faces[node].append(face_index)

    visited: Set[int] = set()
    components: List[List[Tuple[int, int, int, int]]] = []

    for start_index in range(len(faces)):
        if start_index in visited:
            continue
        queue = collections.deque([start_index])
        visited.add(start_index)
        component_indices: List[int] = []

        while queue:
            face_index = queue.popleft()
            component_indices.append(face_index)
            for node in faces[face_index]:
                for neighbor_index in node_to_faces[node]:
                    if neighbor_index not in visited:
                        visited.add(neighbor_index)
                        queue.append(neighbor_index)

        components.append([faces[index] for index in component_indices])

    return components


def face_nodes_sorted(
    faces: Sequence[Tuple[int, int, int, int]]
) -> List[int]:
    return sorted({node for face in faces for node in face})


def write_par_file(
    path: Path,
    node_ids: Sequence[int],
    boundary_kind: str,
    parameterization: str,
) -> None:
    with path.open("w", encoding="utf-8") as fh:
        fh.write(f"{len(node_ids)} {boundary_kind}\n")
        fh.write(f"'{parameterization}'\n")
        for node_id in node_ids:
            fh.write(f"{node_id}\n")


def axial_parameterization(
    coords: Sequence[Tuple[float, float, float]],
    faces: Sequence[Tuple[int, int, int, int]],
) -> str:
    node_ids = face_nodes_sorted(faces)
    z_values = [coords[node_id - 1][2] for node_id in node_ids]
    z_position = (sum(z_values) / len(z_values)) * MM_TO_CM
    return f"4 0.0 0.0 1.0 {-z_position:.15f}"


def cylindrical_parameterization(
    coords: Sequence[Tuple[float, float, float]],
    faces: Sequence[Tuple[int, int, int, int]],
) -> str:
    node_ids = face_nodes_sorted(faces)
    radii = [
        math.hypot(coords[node_id - 1][0], coords[node_id - 1][1])
        for node_id in node_ids
    ]
    radius = (sum(radii) / len(radii)) * MM_TO_CM
    return f"7 0.0 0.0 0.0 {radius:.15f} 1.0 1.0 0.0"


def write_parameterization_files(
    coords: Sequence[Tuple[float, float, float]],
    cells: Sequence[Tuple[int, int, int, int, int, int, int, int]],
    output: Path,
) -> List[str]:
    boundary_faces = compute_boundary_faces(cells)
    shortest_edge = compute_shortest_edge_length(coords, cells)
    tolerance = 0.25 * shortest_edge
    classified = classify_boundary_faces(coords, boundary_faces, tolerance)
    validate_boundary_classification(boundary_faces, classified)
    output_dir = output.parent
    written_files: List[str] = []

    boundary_labels = {
        "inflow": "Inflow10",
        "outflow": "Outflow",
        "barrel": "Wall",
    }

    for name in ("inflow", "outflow", "barrel"):
        faces = classified[name]
        if not faces:
            continue
        if name in {"inflow", "outflow"}:
            parameterization = axial_parameterization(coords, faces)
        else:
            parameterization = cylindrical_parameterization(coords, faces)
        write_par_file(
            output_dir / f"{name}.par",
            face_nodes_sorted(faces),
            boundary_labels[name],
            parameterization,
        )
        written_files.append(f"{name}.par")

    for name in ("inner_cylinder", "inner_transition"):
        components = split_faces_by_connected_component(classified[name])
        for index, component in enumerate(components, start=1):
            if name == "inner_transition":
                parameterization = axial_parameterization(coords, component)
            else:
                parameterization = cylindrical_parameterization(
                    coords, component
                )
            write_par_file(
                output_dir / f"{name}_{index}.par",
                face_nodes_sorted(component),
                "Inflow13",
                parameterization,
            )
            written_files.append(f"{name}_{index}.par")

    return written_files


def write_project_file(output: Path, parameterization_files: Sequence[str]) -> None:
    prj_path = output.with_name("file.prj")
    with prj_path.open("w", encoding="utf-8") as fh:
        fh.write(f"{output.name}\n")
        for par_file in parameterization_files:
            fh.write(f"{par_file}\n")


def prepare_output_directory(output: Path) -> None:
    output_dir = output.parent
    if output_dir.name == "meshDir" and output_dir.exists():
        shutil.rmtree(output_dir)


def main() -> None:
    args = parse_args()
    try:
        config_path = (
            args.input_path or args.config_path or Path("setup.e3d")
        )
        mesh_input = read_parameters(config_path)
        output_path = (
            args.output
            if args.output is not None
            else config_path.resolve().parent / "meshDir" / "Mesh.tri"
        )
        effective_grid_size = (
            args.grid_size
            if args.grid_size is not None
            else mesh_input.grid_size
        )

        if isinstance(mesh_input, BoxMeshInput):
            n_x, n_y, n_z = resolve_box_resolutions(
                mesh_input.lengths,
                mesh_input.n_x,
                mesh_input.n_y,
                mesh_input.n_z,
                effective_grid_size,
                mesh_input.scale_x,
                mesh_input.scale_y,
                mesh_input.scale_z,
            )
            total_elements = n_x * n_y * n_z
            print(
                "Resolved element counts (box): "
                f"x={n_x}, y={n_y}, z={n_z}, total={total_elements}"
            )
            coords, knpr = build_box_coordinates(
                mesh_input.start, mesh_input.lengths, n_x, n_y, n_z
            )
            cells = build_box_connectivity(n_x, n_y, n_z)
        elif isinstance(mesh_input, FullCylinderMeshInput):
            n_inner, n_outer, n_z = resolve_full_cylinder_resolutions(
                mesh_input.outer_diameter,
                mesh_input.length,
                mesh_input.periodicity,
                mesh_input.n_inner,
                mesh_input.n_outer,
                mesh_input.n_z,
                effective_grid_size,
                mesh_input.scale_t,
                mesh_input.scale_r,
                mesh_input.scale_z,
            )
            total_elements = n_z * mesh_input.periodicity * (
                n_inner * n_inner + 2 * n_inner * n_outer
            )
            print(
                "Resolved element counts (full cylinder): "
                f"inner={n_inner}, outer={n_outer}, axial={n_z}, "
                f"periodicity={mesh_input.periodicity}, total={total_elements}"
            )
            coords, knpr = build_full_cylinder_coordinates(
                mesh_input.outer_diameter,
                mesh_input.length,
                mesh_input.start_z,
                mesh_input.periodicity,
                n_inner,
                n_outer,
                n_z,
            )
            cells = build_full_cylinder_connectivity(
                mesh_input.periodicity, n_inner, n_outer, n_z
            )
        elif isinstance(mesh_input, HybridCylinderMeshInput):
            n_inner, n_outer, n_z = resolve_full_cylinder_resolutions(
                mesh_input.outer_diameter,
                mesh_input.length,
                mesh_input.periodicity,
                mesh_input.n_inner,
                mesh_input.n_outer,
                mesh_input.n_z,
                effective_grid_size,
                mesh_input.scale_t,
                mesh_input.scale_r,
                mesh_input.scale_z,
            )
            snapped_sections = snap_hybrid_sections_to_layers(
                mesh_input.sections,
                start_z=mesh_input.start_z,
                length=mesh_input.length,
                n_z=n_z,
            )
            hollow_layers = {
                layer
                for start_idx, end_idx, mode in snapped_sections
                if mode == "hollow"
                for layer in range(start_idx, end_idx)
            }
            full_layers = n_z - len(hollow_layers)
            total_elements = mesh_input.periodicity * (
                full_layers * (n_inner * n_inner + 2 * n_inner * n_outer)
                + len(hollow_layers) * (2 * n_inner * n_outer)
            )
            print(
                "Resolved element counts (hybrid cylinder): "
                f"inner={n_inner}, outer={n_outer}, axial={n_z}, "
                f"periodicity={mesh_input.periodicity}, total={total_elements}"
            )
            print("Snapped hybrid sections:")
            dz = mesh_input.length / n_z
            for start_idx, end_idx, mode in snapped_sections:
                snapped_start = mesh_input.start_z + start_idx * dz
                snapped_end = mesh_input.start_z + end_idx * dz
                print(
                    f"  {mode}: {snapped_start:.6f} -> {snapped_end:.6f} "
                    f"(layers {start_idx}:{end_idx})"
                )
            coords, knpr = build_hybrid_cylinder_coordinates(
                mesh_input.outer_diameter,
                mesh_input.inner_diameter,
                mesh_input.length,
                mesh_input.start_z,
                mesh_input.periodicity,
                n_inner,
                n_outer,
                n_z,
            )
            mark_hybrid_boundaries(
                knpr,
                snapped_sections,
                periodicity=mesh_input.periodicity,
                n_inner=n_inner,
                n_outer=n_outer,
                n_z=n_z,
            )
            cells = list(
                build_hybrid_cylinder_connectivity(
                    mesh_input.periodicity,
                    n_inner,
                    n_outer,
                    n_z,
                    hollow_layers,
                )
            )
            coords, cells, knpr = remove_unused_nodes(coords, cells, knpr)
        else:
            n_t, n_r, n_z = resolve_annular_resolutions(
                mesh_input.outer_diameter,
                mesh_input.inner_diameter,
                mesh_input.length,
                mesh_input.n_t,
                mesh_input.n_r,
                mesh_input.n_z,
                effective_grid_size,
                mesh_input.scale_t,
                mesh_input.scale_r,
                mesh_input.scale_z,
            )
            total_elements = n_t * n_r * n_z
            print(
                "Resolved element counts (annular): "
                f"tangential={n_t}, radial={n_r}, axial={n_z}, "
                f"total={total_elements}"
            )
            coords, knpr = build_annular_coordinates(
                mesh_input.outer_diameter,
                mesh_input.inner_diameter,
                mesh_input.length,
                mesh_input.start_z,
                n_t,
                n_r,
                n_z,
                args.theta_offset,
            )
            cells = build_annular_connectivity(n_t, n_r, n_z)
        cell_list = list(cells)
        prepare_output_directory(output_path)
        write_tri(coords, cell_list, knpr, output_path)
        parameterization_files = write_parameterization_files(
            coords, cell_list, output_path
        )
        write_project_file(output_path, parameterization_files)
    except (
        ValueError,
        FileNotFoundError,
        KeyError,
        configparser.Error,
    ) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        raise SystemExit(1) from exc


if __name__ == "__main__":
    main()
