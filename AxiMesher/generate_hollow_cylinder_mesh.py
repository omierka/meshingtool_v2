#!/usr/bin/env python3
"""Structured hexahedral mesh generator for hollow cylinders.

The script reads the unified `setup.e3d` preprocessing configuration file and
produces an output that mimics the `.tri` layout seen in `Mesh.tri`. Only the
fields needed for the provided mesh (geometry + tangential/radial/axial
resolutions) are parsed.
"""

from __future__ import annotations

import argparse
import configparser
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List, Tuple, Union


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
        default=Path("Mesh.tri"),
        help="Where to write the generated TRI mesh (default: ./Mesh.tri).",
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


MeshInput = Union[AnnularMeshInput, BoxMeshInput, FullCylinderMeshInput]


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
            start, lengths, n_x, n_y, n_z, scale_x, scale_y, scale_z
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
    section: configparser.SectionProxy, key: str, *, default: float
) -> float:
    raw_value = section.get(key, fallback=None)
    if raw_value is None or not raw_value.strip():
        return default
    return float(raw_value)


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
                f"  {x:24.15f}        {y:24.15f}        {z:24.15f}\n"
            )

        fh.write(" KVERT\n")
        for cell in cell_list:
            fh.write("".join(f"{node:8d}" for node in cell) + "\n")

        fh.write(" KNPR\n")
        for flag in knpr:
            fh.write(f"  {flag:d}\n")
    prj_path = output.with_name("file.prj")
    with prj_path.open("w", encoding="utf-8") as fh:
        fh.write(f"{output.name}\n")


def main() -> None:
    args = parse_args()
    try:
        config_path = (
            args.input_path or args.config_path or Path("setup.e3d")
        )
        mesh_input = read_parameters(config_path)

        if isinstance(mesh_input, BoxMeshInput):
            n_x, n_y, n_z = resolve_box_resolutions(
                mesh_input.lengths,
                mesh_input.n_x,
                mesh_input.n_y,
                mesh_input.n_z,
                args.grid_size,
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
                args.grid_size,
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
        else:
            n_t, n_r, n_z = resolve_annular_resolutions(
                mesh_input.outer_diameter,
                mesh_input.inner_diameter,
                mesh_input.length,
                mesh_input.n_t,
                mesh_input.n_r,
                mesh_input.n_z,
                args.grid_size,
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
        write_tri(coords, cells, knpr, args.output)
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
