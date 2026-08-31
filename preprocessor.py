#!/usr/bin/env python3
"""Unified driver that combines run_all_cases.sh and runner.sh functionality."""
from __future__ import annotations

import argparse
import configparser
import os
import re
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Dict, Iterable, List, Optional, Sequence


REQUIRED_BINARIES = (
    "meshhexer-cli",
    "generate_hollow_cylinder_mesh",
    "meshcleaner",
    "hex_VS_triangulation_intersection",
    "meshref",
)

MODULE_BOOTSTRAP = (
    "source /etc/profile >/dev/null 2>&1; "
    "module purge >/dev/null 2>&1; "
    "module load cmake/3.28.3 gcc/latest-v13 openmpi/4.1.6 python/3.13.5; "
    "env -0"
)

STAGE_LEGEND = """\
Stage legend:
[0] start
[1] preprocessing configuration
[2] mindist
[3] coarse mesher
[4] mesh filter
[5] monitor function
[6] mesh refinement
[7] fine mesh filter"""


@dataclass
class DriverConfig:
    monitor_zero_bin_threshold_percent: float = 1.0e-12


@dataclass
class CaseResult:
    folder: Path
    mindist: Optional[float]
    span: Optional[str]
    coarse_mesh_size: Optional[str]
    duration: float
    stage_events: List[str] = field(default_factory=list)
    iteration_marker_used: bool = False
    cleanup_performed: bool = False


class CaseRunner:
    """Runs the full meshing workflow for a single case folder."""

    def __init__(
        self,
        *,
        folder: Path,
        script_dir: Path,
        num_proc: int,
        use_srun: bool,
        cleanup: bool,
        silent: bool,
        base_env: Dict[str, str],
        driver_config: DriverConfig,
        stage_callback: Optional[Callable[[str, str], None]] = None,
        protocol_file: Optional[Path] = None,
    ) -> None:
        self.script_dir = script_dir
        self.folder = folder if folder.is_absolute() else script_dir / folder
        self.num_proc = num_proc
        self.use_srun = use_srun
        self.cleanup = cleanup
        self.silent = silent
        self.env = dict(base_env)
        self.stage_callback = stage_callback
        self.monitor_zero_bin_threshold_percent = driver_config.monitor_zero_bin_threshold_percent
        self.protocol_file = protocol_file

        self.monitor_summary_file = self.folder / "monitor_summary.txt"
        self.monitor_summary_vol_file = self.folder / "monitor_summary_volumetric.txt"

        self.stage_events: List[str] = []
        self.iteration_marker_used = False
        self.mindist_value: Optional[float] = None
        self.span_value: Optional[str] = None
        self.coarse_mesh_size_value: Optional[str] = None
        self._coarse_scaling_factor: Optional[float] = None

    def _cmd_path(self, name: str) -> str:
        return str((self.script_dir / name).resolve())

    def _mpi_command(self, executable: str, *args: str) -> List[str]:
        if self.use_srun:
            return ["srun", executable, *args]
        return ["mpirun", "-np", str(self.num_proc), executable, *args]

    def _record_stage(self, label: int) -> None:
        marker = f"[{label}]"
        self.stage_events.append(marker)
        self._append_protocol_line(marker)
        if self.stage_callback:
            self.stage_callback("stage", marker)
        elif self.silent:
            print(marker)

    def _record_iteration_marker(self) -> None:
        self.iteration_marker_used = True
        self.stage_events.append("*")
        self._append_protocol_line("*")
        if self.stage_callback:
            self.stage_callback("iteration", "*")
        elif self.silent:
            print("*")

    def _emit_info(self, value: str) -> None:
        self._append_protocol_line(value)
        if self.stage_callback:
            self.stage_callback("info", value)
            return
        print(value)

    def _append_protocol(self, text: str) -> None:
        if self.protocol_file is None or not text:
            return
        with self.protocol_file.open("a", encoding="utf-8") as handle:
            handle.write(text)

    def _append_protocol_line(self, line: str) -> None:
        if not line:
            return
        self._append_protocol(f"{line}\n")

    def _run_command(
        self,
        args: Sequence[str],
        *,
        capture_output: bool = False,
    ) -> Optional[str]:
        command = list(map(str, args))
        stdout = subprocess.PIPE if capture_output else None
        stderr = None
        if self.protocol_file is not None:
            stdout = subprocess.PIPE
            stderr = subprocess.PIPE
        elif self.silent and not capture_output:
            stdout = subprocess.DEVNULL

        self._append_protocol_line(f"$ {' '.join(command)}")
        try:
            result = subprocess.run(
                command,
                cwd=self.script_dir,
                env=self.env,
                stdout=stdout,
                stderr=stderr,
                check=True,
                text=True,
            )
        except subprocess.CalledProcessError as exc:
            self._append_protocol(exc.stdout or "")
            self._append_protocol(exc.stderr or "")
            raise
        self._append_protocol(result.stdout or "")
        self._append_protocol(result.stderr or "")
        return result.stdout if capture_output else None

    def _ensure_case_exists(self) -> None:
        if not self.folder.is_dir():
            raise FileNotFoundError(f"Case folder '{self.folder}' does not exist")

    def _cleanup_generated_files(self) -> None:
        print(f"Cleanup enabled. Removing generated files under/for: {self.folder}")
        patterns = ("*.vtu", "*.tri", "*meshDir*")
        for pattern in patterns:
            for path in self.folder.glob(pattern):
                if path.is_dir():
                    shutil.rmtree(path, ignore_errors=True)
                else:
                    path.unlink(missing_ok=True)
        extra_files = (
            "area.txt",
            "coarse_size_distribution_histogram.txt",
            "hex_intersection_intersections.pvtu",
            "hex_intersection_tets.pvtu",
            "hex_mesh.pvtu",
            "preprocessing_protocol.txt",
            "size_distribution_histogram.txt",
        )
        for file_name in extra_files:
            (self.folder / file_name).unlink(missing_ok=True)
        self.monitor_summary_file.unlink(missing_ok=True)
        self.monitor_summary_vol_file.unlink(missing_ok=True)

    def _remove_monitor_files(self) -> None:
        self.monitor_summary_file.unlink(missing_ok=True)
        self.monitor_summary_vol_file.unlink(missing_ok=True)

    def _configure_case_for_preprocessing(self) -> None:
        try:
            self._run_command(
                (
                    self._cmd_path("meshhexer-cli"),
                    "report",
                    "--configure-case-for-preprocessing",
                    str(self.folder),
                )
            )
        except subprocess.CalledProcessError as exc:
            raise RuntimeError(
                "meshhexer-cli report --configure-case-for-preprocessing failed."
            ) from exc

    def _compute_initial_parameters(self) -> None:
        try:
            output = self._run_command(
                (
                    self._cmd_path("meshhexer-cli"),
                    "--checkpoint-path",
                    str(self.folder / "MINGAP.vtu"),
                    "min-gap",
                    str(self.folder / "surface.off"),
                ),
                capture_output=True,
            )
        except subprocess.CalledProcessError as exc:
            raise RuntimeError(
                "meshhexer-cli min-gap failed."
            ) from exc
        if output is None:
            raise RuntimeError(
                "meshhexer-cli min-gap produced no output."
            )
        parts = output.split()
        if len(parts) < 3:
            raise RuntimeError(
                "Unable to parse meshhexer-cli min-gap output: "
                f"'{output.strip()}'."
            )
        mindist_raw, span_raw, coarse_raw = parts[:3]
        self.mindist_value = float(mindist_raw)
        self.span_value = span_raw
        self.coarse_mesh_size_value = coarse_raw
        self._coarse_scaling_factor = self._derive_coarse_scaling_factor(
            self.mindist_value, span_raw, coarse_raw
        )
        self._emit_info(f"mindist={mindist_raw}")
        self._emit_info(f"histogram_span={span_raw}")
        self._emit_info(f"CoarseMeshSize={coarse_raw}")
        manual_mingap = self._read_manual_mingap_override()
        if manual_mingap is not None:
            self.mindist_value = manual_mingap
            self._emit_info(f"mindist (manual override)={manual_mingap}")
            self._recompute_coarse_mesh_size()

    def _read_manual_mingap_override(self) -> Optional[float]:
        config_path = self.folder / "setup.e3d"
        if not config_path.is_file():
            return None
        parser = configparser.ConfigParser()
        try:
            if not parser.read(config_path):
                return None
        except configparser.Error:
            return None
        section = "E3DGeometryData/Machine"
        if not parser.has_section(section):
            return None
        raw_value = parser.get(section, "UserDefinedMinGap", fallback=None)
        if raw_value is None:
            return None
        raw_value = raw_value.strip()
        if not raw_value:
            return None
        try:
            return float(raw_value)
        except ValueError:
            self._emit_info(
                f"Warning: invalid UserDefinedMinGap value '{raw_value}' in setup.e3d; ignoring."
            )
            return None

    def _derive_coarse_scaling_factor(
        self, mindist: float, span_raw: str, coarse_raw: str
    ) -> Optional[float]:
        if mindist <= 0.0:
            return None
        try:
            coarse_value = float(coarse_raw)
        except ValueError:
            return None
        span_power = self._span_exponent_value(span_raw)
        if span_power is None:
            return None
        try:
            denom = mindist * (3 ** span_power)
        except OverflowError:
            return None
        if denom == 0.0:
            return None
        return coarse_value / denom

    def _span_exponent_value(self, span_raw: Optional[str]) -> Optional[int]:
        if span_raw is None:
            return None
        try:
            value = int(float(span_raw))
        except (TypeError, ValueError):
            return None
        if value < 0:
            return None
        return value

    def _recompute_coarse_mesh_size(self) -> None:
        if (
            self.mindist_value is None
            or self.span_value is None
            or self._coarse_scaling_factor is None
        ):
            return
        span_power = self._span_exponent_value(self.span_value)
        if span_power is None:
            return
        try:
            new_value = self.mindist_value * self._coarse_scaling_factor * (
                3 ** span_power
            )
        except OverflowError:
            return
        self.coarse_mesh_size_value = f"{new_value:.6f}"
        self._emit_info(
            f"CoarseMeshSize (manual override)={self.coarse_mesh_size_value}"
        )

    def _run_mesh_workflow(self, mindist_value: str) -> None:
        assert self.span_value is not None
        assert self.coarse_mesh_size_value is not None
        coarse_dir = self.folder / "Coarse_meshDir"
        coarse_dir.mkdir(parents=True, exist_ok=True)

        self._run_command(
            (
                self._cmd_path("generate_hollow_cylinder_mesh"),
                "-i",
                str(self.folder / "setup.e3d"),
                "-o",
                str(coarse_dir / "Mesh.tri"),
                "-s",
                self.coarse_mesh_size_value,
            )
        )

        self._record_stage(3)
        self._run_command(
            self._mpi_command(
                self._cmd_path("meshcleaner"),
                "-h",
                str(coarse_dir / "Mesh.tri"),
                "-t",
                str(self.folder / "surface.off"),
                "-s",
                "1.0",
                "-o",
                str(self.folder),
            )
        )

        coarse_mesh = coarse_dir / "Mesh.tri"
        filtered_mesh = self.folder / "Filtered.tri"
        backup_path = coarse_dir / "Mesh_BU.tri"
        if backup_path.exists():
            backup_path.unlink()
        if coarse_mesh.exists():
            coarse_mesh.rename(backup_path)
        if not filtered_mesh.exists():
            raise FileNotFoundError(f"Filtered mesh '{filtered_mesh}' missing")
        if coarse_mesh.exists():
            coarse_mesh.unlink()
        filtered_mesh.rename(coarse_dir / "Mesh.tri")

        self._record_stage(4)
        self._run_command(
            self._mpi_command(
                self._cmd_path("hex_VS_triangulation_intersection"),
                "-h",
                str(coarse_dir / "Mesh.tri"),
                "-t",
                str(self.folder / "MINGAP.vtu"),
                "-m",
                mindist_value,
                "-o",
                str(self.folder),
            )
        )

        mesh_names_file = self.script_dir / "mesh_names.offs"
        mesh_names_file.write_text(f"1\n{self.folder / 'surface.off'}\n", encoding="utf-8")

        self._record_stage(5)
        self._run_command(
            (
                self._cmd_path("meshref"),
                "-f",
                str(self.folder),
                "-d",
                self.span_value,
            )
        )

        self._record_stage(6)
        try:
            span_int = int(float(self.span_value))
        except ValueError as exc:
            raise RuntimeError(f"Invalid span value '{self.span_value}'") from exc
        span_int = max(span_int, 0)
        ispan = max(span_int - 1, 0)
        refined_mesh = (
            self.folder / f"RefinedCleanMesh_lvl{ispan}_refined_clean.vtu"
        )

        self._remove_monitor_files()
        self._run_command(
            self._mpi_command(
                self._cmd_path("meshcleaner"),
                "-h",
                str(refined_mesh),
                "-t",
                str(self.folder / "surface.off"),
                "-s",
                "1.0",
                "-o",
                str(self.folder),
            )
        )
        self._record_stage(7)

    def _read_monitor_first_bin_percent(self, path: Path) -> Optional[float]:
        if not path.exists():
            return None
        data = path.read_text(encoding="utf-8")
        match = re.search(r"0=\[\s*([0-9.]+)%\]", data)
        if not match:
            return None
        try:
            return float(match.group(1))
        except ValueError:
            return None

    def _maybe_apply_monitor_correction(self) -> None:
        volum_zero = self._read_monitor_first_bin_percent(self.monitor_summary_vol_file)
        threshold = max(self.monitor_zero_bin_threshold_percent, 0.0)
        if volum_zero is None or volum_zero > threshold:
            return
        assert self.coarse_mesh_size_value is not None
        assert self.span_value is not None

        corrected_coarse = float(self.coarse_mesh_size_value) / 3.0
        corrected_span = max(int(float(self.span_value)) - 1, 1)

        self.coarse_mesh_size_value = f"{corrected_coarse:.12g}"
        self.span_value = str(corrected_span)
        self._emit_info(f"histogram_span={self.span_value}")
        self._emit_info(f"CoarseMeshSize={self.coarse_mesh_size_value}")
        self._record_iteration_marker()
        self._run_mesh_workflow(str(self.mindist_value))

    def run(self) -> CaseResult:
        start = time.perf_counter()
        if self.cleanup:
            self._ensure_case_exists()
            self._cleanup_generated_files()
            duration = time.perf_counter() - start
            return CaseResult(
                folder=self.folder,
                mindist=None,
                span=None,
                coarse_mesh_size=None,
                duration=duration,
                stage_events=[],
                iteration_marker_used=False,
                cleanup_performed=True,
            )

        self._ensure_case_exists()
        self._remove_monitor_files()
        self._record_stage(0)
        self._configure_case_for_preprocessing()
        self._record_stage(1)
        self._compute_initial_parameters()
        self._record_stage(2)
        assert self.mindist_value is not None
        self._run_mesh_workflow(str(self.mindist_value))
        self._maybe_apply_monitor_correction()
        duration = time.perf_counter() - start
        return CaseResult(
            folder=self.folder,
            mindist=self.mindist_value,
            span=self.span_value,
            coarse_mesh_size=self.coarse_mesh_size_value,
            duration=duration,
            stage_events=list(self.stage_events),
            iteration_marker_used=self.iteration_marker_used,
            cleanup_performed=False,
        )


def load_driver_config(config_path: Path) -> DriverConfig:
    """Load optional driver-specific settings from preprocessor.cfg."""
    config = configparser.ConfigParser(
        interpolation=None,
        comment_prefixes=("#", ";"),
        inline_comment_prefixes=("#", ";"),
    )
    driver_settings = DriverConfig()
    read_files = config.read(config_path, encoding="utf-8")
    if not read_files or not config.has_section("Driver"):
        return driver_settings
    section = config["Driver"]
    raw_threshold = section.get("monitor_zero_bin_threshold_percent")
    if raw_threshold is not None:
        try:
            driver_settings.monitor_zero_bin_threshold_percent = float(raw_threshold)
        except ValueError as exc:
            raise RuntimeError(
                f"Invalid monitor_zero_bin_threshold_percent '{raw_threshold}' in {config_path}"
            ) from exc
    return driver_settings


def ensure_required_binaries(script_dir: Path) -> None:
    missing = []
    for binary in REQUIRED_BINARIES:
        path = script_dir / binary
        if not path.exists() or not os.access(path, os.X_OK):
            missing.append(binary)
    if missing:
        raise FileNotFoundError(
            "Required binaries missing: " + ", ".join(sorted(missing))
        )


def capture_module_environment() -> Dict[str, str]:
    result = subprocess.run(
        ["/bin/bash", "-lc", MODULE_BOOTSTRAP],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=False,
        check=True,
    )
    env: Dict[str, str] = {}
    for chunk in result.stdout.split(b"\0"):
        if not chunk:
            continue
        key, _, value = chunk.partition(b"=")
        if key:
            env[key.decode()] = value.decode()
    return env


def prepare_base_environment(script_dir: Path, skip_modules: bool) -> Dict[str, str]:
    env = os.environ.copy()
    if not skip_modules:
        env.update(capture_module_environment())
    env.setdefault("PREPROCESSOR_CONFIG", str(script_dir / "preprocessor.cfg"))
    return env


def gather_cases(cases_dir: Path) -> List[Path]:
    if not cases_dir.is_dir():
        raise FileNotFoundError(f"Cases directory '{cases_dir}' not found")
    return sorted(path for path in cases_dir.iterdir() if path.is_dir())


def format_duration(seconds: float) -> str:
    return f"{seconds:.3f} [s]"


def read_mesh_counts(filtered_tri: Path) -> tuple[int, int]:
    if not filtered_tri.exists():
        raise FileNotFoundError(f"Filtered mesh '{filtered_tri}' missing")
    with filtered_tri.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if "NEL,NVT" in line:
                parts = line.split()
                if len(parts) < 2:
                    break
                try:
                    return int(parts[0]), int(parts[1])
                except ValueError:
                    break
    raise RuntimeError(f"Unable to parse NEL/NVT from '{filtered_tri}'")


def read_monitor_summary(case_path: Path, filename: str = "monitor_summary.txt") -> str:
    summary = case_path / filename
    if not summary.exists():
        return ""
    # Only keep the portion starting at the 0-bin definition, e.g. "0=[...]".
    for line in summary.read_text(encoding="utf-8").splitlines():
        zero_idx = line.find("0=")
        if zero_idx != -1:
            return line[zero_idx:].strip()
    return ""


def cleanup_case_pngs(case_path: Path) -> None:
    for png in case_path.glob("*.png"):
        if png.is_file():
            png.unlink()
    extra_files = (
        "size_distribution_histogram.txt",
        "coarse_size_distribution_histogram.txt",
        "hex_intersection_intersections.pvtu",
        "hex_intersection_tets.pvtu",
        "hex_mesh.pvtu",
        "preprocessing_protocol.txt",
    )
    for file_name in extra_files:
        (case_path / file_name).unlink(missing_ok=True)


def run_case_command(args: argparse.Namespace, script_dir: Path) -> None:
    ensure_required_binaries(script_dir)
    env = prepare_base_environment(script_dir, args.skip_modules)
    config_location = env.get("PREPROCESSOR_CONFIG") or str(script_dir / "preprocessor.cfg")
    config_path = Path(config_location)
    driver_config = load_driver_config(config_path)
    raw_folder = Path(args.folder)
    case_folder = raw_folder if raw_folder.is_absolute() else Path.cwd() / raw_folder
    runner = CaseRunner(
        folder=case_folder,
        script_dir=script_dir,
        num_proc=args.num_proc,
        use_srun=args.use_srun,
        cleanup=args.cleanup,
        silent=args.silent,
        base_env=env,
        driver_config=driver_config,
    )
    try:
        result = runner.run()
    except subprocess.CalledProcessError as exc:
        raise SystemExit(exc.returncode) from exc
    except Exception as exc:
        raise SystemExit(f"Case execution failed: {exc}") from exc
    if result.cleanup_performed:
        print(f"Cleanup complete for '{result.folder}'.")


def run_all_command(args: argparse.Namespace, script_dir: Path) -> None:
    raw_cases_dir = Path(args.cases_dir)
    ensure_required_binaries(script_dir)
    resolved_cases_dir = (
        raw_cases_dir if raw_cases_dir.is_absolute() else Path.cwd() / raw_cases_dir
    )
    cases = gather_cases(resolved_cases_dir)
    if not cases:
        raise SystemExit(f"No cases found under '{resolved_cases_dir}'")

    env = prepare_base_environment(script_dir, args.skip_modules)
    config_location = env.get("PREPROCESSOR_CONFIG") or str(script_dir / "preprocessor.cfg")
    config_path = Path(config_location)
    driver_config = load_driver_config(config_path)
    print(STAGE_LEGEND)
    case_names = [case.name for case in cases]
    name_width = max(len(name) for name in case_names)
    time_column_width = len("9999.999 [s]")
    mindist_column_width = len("mindist=999.999")
    nel_column_width = len("99,999,999")
    iteration_pad_len = len("*[3][4][5][6][7]")

    for case_path, case_name in zip(cases, case_names):
        print(f"[{case_name:<{name_width}}]:", end="", flush=True)
        protocol_file = case_path / "preprocessing_protocol.txt"
        protocol_file.write_text("", encoding="utf-8")

        def stage_printer(event_type: str, payload: str) -> None:
            if event_type in {"stage", "iteration"}:
                print(payload, end="", flush=True)

        runner = CaseRunner(
            folder=case_path,
            script_dir=script_dir,
            num_proc=args.num_proc,
            use_srun=args.use_srun,
            cleanup=args.clean,
            silent=True,
            base_env=env,
            driver_config=driver_config,
            stage_callback=stage_printer,
            protocol_file=protocol_file,
        )
        try:
            result = runner.run()
        except subprocess.CalledProcessError as exc:
            print()
            raise SystemExit(
                f"Case '{case_name}' failed with exit code {exc.returncode}. "
                f"See '{protocol_file}' for the full protocol."
            ) from exc
        except Exception as exc:
            print()
            raise SystemExit(
                f"Case '{case_name}' failed: {exc}. "
                f"See '{protocol_file}' for the full protocol."
            ) from exc
        if args.clean:
            print()
            cleanup_case_pngs(case_path)
            continue

        filtered_tri = case_path / "Filtered.tri"
        nel, _ = read_mesh_counts(filtered_tri)
        monitor_summary = read_monitor_summary(case_path)
        monitor_summary_vol = read_monitor_summary(
            case_path, "monitor_summary_volumetric.txt"
        )

        elapsed = format_duration(result.duration)
        print(" :: ", end="")
        if not result.iteration_marker_used:
            print(" " * iteration_pad_len, end="")
        print(f"{elapsed:>{time_column_width}}", end="")
        if result.mindist is not None:
            mindist_str = f"mindist={result.mindist:.3f}"
            print(f"  {mindist_str:<{mindist_column_width}}", end="")
            if result.span:
                print(f" span={result.span}", end="")
        formatted_nel = f"{nel:,}"
        print(f"  :: NEL={formatted_nel:>{nel_column_width}}", end="")
        if monitor_summary:
            print(f"  [COUNT]: {monitor_summary}", end="")
        if monitor_summary_vol:
            print(f"  [VOLUME]: {monitor_summary_vol}", end="")
        print()


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Unified case/all-case driver for the meshing tool.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument(
        "--skip-modules",
        action="store_true",
        help="Skip invoking environment modules before running commands.",
    )
    common.add_argument(
        "-u",
        "--use-srun",
        action="store_true",
        help="Launch MPI stages with 'srun <executable>' instead of mpirun.",
    )

    subparsers = parser.add_subparsers(dest="command", required=True)

    case_parser = subparsers.add_parser(
        "case",
        parents=[common],
        help="Run the workflow for a single case (replacement for runner.sh).",
    )
    case_parser.add_argument(
        "-f",
        "--folder",
        default="PROFEX",
        help="Case folder to process, relative to the current working directory unless absolute.",
    )
    case_parser.add_argument(
        "-n",
        "--num-proc",
        type=int,
        default=4,
        help="Number of MPI ranks for mesh tools.",
    )
    case_parser.add_argument(
        "-c",
        "--cleanup",
        "--clean",
        action="store_true",
        dest="cleanup",
        help="Remove generated files instead of running the workflow.",
    )
    case_parser.add_argument(
        "-s",
        "--silent",
        action="store_true",
        help="Silence downstream tool output while keeping stage markers.",
    )
    case_parser.set_defaults(func=run_case_command)

    all_parser = subparsers.add_parser(
        "all",
        parents=[common],
        help="Run the workflow for every case in CASES (replacement for run_all_cases.sh).",
    )
    all_parser.add_argument(
        "-c",
        "--clean",
        action="store_true",
        help="Propagate cleanup to each case instead of running meshes.",
    )
    all_parser.add_argument(
        "-n",
        "--num-proc",
        type=int,
        default=64,
        help="Number of MPI ranks to forward to each case run.",
    )
    all_parser.add_argument(
        "--cases-dir",
        default="CASES",
        help="Directory that stores individual case folders, relative to the current working directory unless absolute.",
    )
    all_parser.set_defaults(func=run_all_command)
    return parser


def main(argv: Optional[Sequence[str]] = None) -> None:
    script_dir = Path(__file__).resolve().parent
    parser = build_parser()
    args = parser.parse_args(argv)
    args.func(args, script_dir)


if __name__ == "__main__":
    main()
