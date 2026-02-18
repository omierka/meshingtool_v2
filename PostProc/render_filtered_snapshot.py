#!/usr/bin/env python3
"""
Render Filtered.vtu using a ParaView state file and save a PNG image.

Examples:
    pvbatch --force-offscreen-rendering render_filtered_snapshot.py
    pvbatch --force-offscreen-rendering render_filtered_snapshot.py \
        --state /path/to/view.pvsm \
        --dataset /path/to/Filtered.vtu \
        --output /path/to/Filtered.png
"""

import argparse
import math
import os
import re
import shutil
import sys
import tempfile

from paraview.simple import (
    FindSource,
    GetActiveView,
    GetAnimationScene,
    GetSources,
    Hide,
    LoadState,
    OpenDataFile,
    Render,
    SaveScreenshot,
    Show,
)


STATE_FILE = "view.pvsm"
DATASET_FILE = "Filtered.vtu"
OUTPUT_IMAGE = "Filtered.png"
SURFACE_FILE = "surface.off"
SURFACE_IMAGE = "Surface.png"
BASE_RESOLUTION = (1920, 1080)


def _parse_arguments(script_dir: str) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Render Filtered.vtu using a ParaView state file."
    )
    parser.add_argument(
        "--state",
        default=os.path.join(script_dir, STATE_FILE),
        help=f"Path to the ParaView state file (default: {STATE_FILE} next to this script).",
    )
    parser.add_argument(
        "--dataset",
        default=os.path.join(script_dir, DATASET_FILE),
        help=f"Path to the dataset file (default: {DATASET_FILE} next to this script).",
    )
    parser.add_argument(
        "--output",
        default=os.path.join(script_dir, OUTPUT_IMAGE),
        help=f"Path to the output image (default: {OUTPUT_IMAGE} next to this script).",
    )
    parser.add_argument(
        "--surface",
        default=os.path.join(script_dir, SURFACE_FILE),
        help=(
            f"Optional path to the surface geometry (default: {SURFACE_FILE} next to this script; "
            "pass an empty string to disable)."
        ),
    )
    parser.add_argument(
        "--surface-output",
        default=os.path.join(script_dir, SURFACE_IMAGE),
        help=f"Path to the surface snapshot (default: {SURFACE_IMAGE} next to this script).",
    )
    parser.add_argument(
        "--resolution-scale",
        type=float,
        default=1.5,
        help=(
            f"Scale factor applied to the default {BASE_RESOLUTION[0]}x{BASE_RESOLUTION[1]} "
            "resolution (default: 1.5)."
        ),
    )
    return parser.parse_args()


def _prepare_state_file(state_path: str, dataset_path: str) -> str:
    """Return the path of a temporary state file with absolute dataset paths."""
    state_text = open(state_path, "r", encoding="utf-8").read()
    dataset_path = os.path.abspath(dataset_path)

    # Replace any path that ends with Filtered.vtu (state stores absolute paths).
    pattern = re.compile(r'(?<=["\'])[^"\']*/Filtered\.vtu(?=["\'])')
    updated_state, replacements = pattern.subn(dataset_path, state_text)
    if replacements == 0:
        raise RuntimeError("Failed to update dataset path inside the state file.")

    tmp_dir = tempfile.mkdtemp(prefix="pv_state_")
    patched_state_path = os.path.join(tmp_dir, os.path.basename(state_path))
    with open(patched_state_path, "w", encoding="utf-8") as handle:
        handle.write(updated_state)
    return patched_state_path


def _find_dataset_source(dataset_path: str):
    dataset_name = os.path.basename(dataset_path)
    source = FindSource(dataset_name)
    if source is None:
        dataset_lower = dataset_name.lower()
        for (name, _), proxy in GetSources().items():
            if dataset_lower in name.lower():
                source = proxy
                break
    if source is None:
        raise RuntimeError(f"Could not locate a source matching {dataset_name} in the state.")
    return source


def _fit_camera_to_dataset(view, source) -> None:
    """Move the camera along its current direction so the dataset fills the frame."""

    bounds = source.GetDataInformation().GetBounds()
    if not bounds or bounds[0] >= bounds[1]:
        return

    center = [
        0.5 * (bounds[0] + bounds[1]),
        0.5 * (bounds[2] + bounds[3]),
        0.5 * (bounds[4] + bounds[5]),
    ]
    diag = math.sqrt(
        (bounds[1] - bounds[0]) ** 2
        + (bounds[3] - bounds[2]) ** 2
        + (bounds[5] - bounds[4]) ** 2
    )
    radius = max(diag * 0.5, 1e-6)

    camera_position = list(view.GetProperty("CameraPosition"))
    camera_focal_point = list(view.GetProperty("CameraFocalPoint"))
    direction = [
        camera_position[0] - camera_focal_point[0],
        camera_position[1] - camera_focal_point[1],
        camera_position[2] - camera_focal_point[2],
    ]
    distance = math.sqrt(direction[0] ** 2 + direction[1] ** 2 + direction[2] ** 2)
    if distance < 1e-6:
        direction = [0.0, 0.0, 1.0]
        distance = 1.0
    else:
        direction = [d / distance for d in direction]

    if bool(view.GetProperty("CameraParallelProjection")):
        scale = radius * 1.05
        view.CameraParallelScale = scale
        distance = max(distance, scale * 2.0)
    else:
        fov_vertical = math.radians(float(view.GetProperty("CameraViewAngle")))
        fov_vertical = max(fov_vertical, math.radians(5.0))
        aspect = max(view.ViewSize[0], 1) / max(view.ViewSize[1], 1)
        fov_horizontal = 2.0 * math.atan(math.tan(fov_vertical / 2.0) * aspect)
        max_vertical = radius / math.tan(fov_vertical / 2.0)
        max_horizontal = radius / math.tan(max(fov_horizontal, math.radians(5.0)) / 2.0)
        distance = max(max_vertical, max_horizontal) * 1.05

    view.CameraFocalPoint = center
    view.CameraPosition = [
        center[0] + direction[0] * distance,
        center[1] + direction[1] * distance,
        center[2] + direction[2] * distance,
    ]


def main() -> None:
    script_dir = os.path.abspath(os.path.dirname(__file__) or ".")
    args = _parse_arguments(script_dir)
    state_path = os.path.abspath(args.state)
    dataset_path = os.path.abspath(args.dataset)
    output_path = os.path.abspath(args.output)
    surface_path = os.path.abspath(args.surface) if args.surface else None
    surface_output_path = os.path.abspath(args.surface_output) if args.surface_output else None
    if args.resolution_scale <= 0:
        raise ValueError("resolution-scale must be positive.")
    resolution = [
        max(1, int(round(dim * args.resolution_scale))) for dim in BASE_RESOLUTION
    ]

    for path, label in [(state_path, "state"), (dataset_path, "dataset")]:
        if not os.path.isfile(path):
            raise FileNotFoundError(f"Missing {label} file: {path}")
    if surface_path and not os.path.isfile(surface_path):
        print(f"Surface file not found, skipping surface snapshot: {surface_path}", file=sys.stderr)
        surface_path = None

    patched_state_dir = None
    try:
        patched_state_path = _prepare_state_file(state_path, dataset_path)
        patched_state_dir = os.path.dirname(patched_state_path)

        LoadState(patched_state_path)
        animation_scene = GetAnimationScene()
        animation_scene.UpdateAnimationUsingDataTimeSteps()
        animation_scene.GoToLast()

        view = GetActiveView()
        if view is None:
            raise RuntimeError("State file did not create an active view.")
        dataset_source = _find_dataset_source(dataset_path)
        _fit_camera_to_dataset(view, dataset_source)
        Render(view)
        SaveScreenshot(output_path, viewOrLayout=view, ImageResolution=resolution)

        if surface_path and surface_output_path:
            surface_source = OpenDataFile(surface_path)
            Hide(dataset_source, view)
            Show(surface_source, view)
            Render(view)
            SaveScreenshot(surface_output_path, viewOrLayout=view, ImageResolution=resolution)
            Hide(surface_source, view)
            Show(dataset_source, view)
    finally:
        if patched_state_dir:
            shutil.rmtree(patched_state_dir, ignore_errors=True)


if __name__ == "__main__":
    main()
