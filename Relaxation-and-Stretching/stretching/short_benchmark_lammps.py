#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import importlib.util
import os
import re
import shlex
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from types import ModuleType


MPIEXEC = "mpiexec"
LAMMPS_BIN = "/home/hkust/mylammps/build/lmp"
DEFAULT_MPI_MIN = 1
DEFAULT_MPI_MAX = 6
DEFAULT_OMP_MIN = 1
DEFAULT_OMP_MAX = 6
DEFAULT_RUN_TIME = 0.2
DEFAULT_MAX_PARALLEL = 1
DEFAULT_OUTPUT_INTERVAL = 10**9
LOOP_TIME_PATTERN = re.compile(r"Loop time of\s+([0-9.eE+-]+)\s+on\s+\d+\s+procs\s+for\s+(\d+)\s+steps")
VARIABLE_LINE = re.compile(r"^\s*variable\s+(\w+)\s+index\s+(.+?)\s*$")
REQUIRED_VARIABLES = ("restart_file", "target_T", "force_mag", "protocol_tag")


@dataclass(frozen=True)
class RunConfig:
    mpi_ranks: int
    omp_threads: int

    @property
    def label(self) -> str:
        return f"np{self.mpi_ranks}_omp{self.omp_threads}"


@dataclass
class InputConfig:
    input_file: Path
    restart_file: str
    target_t: str
    force_mag: str
    protocol_tag: str


@dataclass
class RunResult:
    config: RunConfig
    return_code: int
    total_steps: int
    total_loop_time: float
    timesteps_per_second: float
    stdout_log: Path
    lammps_log: Path


def load_output_root_helper(helper_root: Path) -> ModuleType:
    helper_path = helper_root / "build_output_root.py"
    spec = importlib.util.spec_from_file_location("build_output_root", helper_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load helper module: {helper_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run a very short benchmark for one packaged stretching .lmp across an MPI x OMP grid. "
            "The benchmark runs from the chosen workdir so restart_file and outputs resolve there."
        )
    )
    parser.add_argument(
        "--input",
        required=True,
        type=Path,
        help="Path to one packaged stretching .lmp file.",
    )
    parser.add_argument(
        "--workdir",
        type=Path,
        default=Path.cwd(),
        help="Directory containing the restart file and receiving benchmark outputs (default: current directory).",
    )
    parser.add_argument(
        "--helper-root",
        type=Path,
        default=Path(__file__).resolve().parent.parent / "helpers",
        help="Directory containing build_output_root.py (default: sibling ../helpers relative to this script).",
    )
    parser.add_argument(
        "--restart-file",
        default=None,
        help="Override restart_file for the benchmark run.",
    )
    parser.add_argument(
        "--run-time",
        type=float,
        default=DEFAULT_RUN_TIME,
        help=f"Short benchmark runtime in LJ time units (default: {DEFAULT_RUN_TIME}).",
    )
    parser.add_argument(
        "--sample-every",
        type=int,
        default=DEFAULT_OUTPUT_INTERVAL,
        help=f"Thermo/print interval override used during the benchmark (default: {DEFAULT_OUTPUT_INTERVAL}).",
    )
    parser.add_argument(
        "--dump-every",
        type=int,
        default=DEFAULT_OUTPUT_INTERVAL,
        help=f"Dump interval override used during the benchmark (default: {DEFAULT_OUTPUT_INTERVAL}).",
    )
    parser.add_argument(
        "--restart-every",
        type=int,
        default=DEFAULT_OUTPUT_INTERVAL,
        help=f"Restart interval override used during the benchmark (default: {DEFAULT_OUTPUT_INTERVAL}).",
    )
    parser.add_argument("--mpi-min", type=int, default=DEFAULT_MPI_MIN, help=f"Minimum MPI ranks (default: {DEFAULT_MPI_MIN}).")
    parser.add_argument("--mpi-max", type=int, default=DEFAULT_MPI_MAX, help=f"Maximum MPI ranks (default: {DEFAULT_MPI_MAX}).")
    parser.add_argument("--omp-min", type=int, default=DEFAULT_OMP_MIN, help=f"Minimum OMP threads (default: {DEFAULT_OMP_MIN}).")
    parser.add_argument("--omp-max", type=int, default=DEFAULT_OMP_MAX, help=f"Maximum OMP threads (default: {DEFAULT_OMP_MAX}).")
    parser.add_argument(
        "--max-parallel",
        type=int,
        default=DEFAULT_MAX_PARALLEL,
        help=f"Maximum benchmark configurations to run concurrently (default: {DEFAULT_MAX_PARALLEL}).",
    )
    parser.add_argument("--mpiexec", default=MPIEXEC, help=f"MPI launcher executable (default: {MPIEXEC}).")
    parser.add_argument("--lammps-bin", default=LAMMPS_BIN, help=f"LAMMPS executable (default: {LAMMPS_BIN}).")
    parser.add_argument("--dry-run", action="store_true", help="Print the benchmark matrix and commands without running LAMMPS.")
    return parser.parse_args()


def validate_environment(*, dry_run: bool, mpiexec: str, lammps_bin: str) -> None:
    if dry_run:
        return

    missing = []
    for binary in (mpiexec, lammps_bin):
        resolved = binary if os.path.sep in binary else shutil.which(binary)
        if resolved is None or (os.path.sep in binary and not Path(binary).exists()):
            missing.append(binary)

    if missing:
        raise SystemExit(f"Cannot find required executable(s): {', '.join(missing)}")


def read_index_variables(input_file: Path) -> dict[str, str]:
    variables: dict[str, str] = {}
    for raw_line in input_file.read_text(encoding="utf-8").splitlines():
        line = raw_line.split("#", 1)[0].rstrip()
        if not line:
            continue
        match = VARIABLE_LINE.match(line)
        if match:
            variables[match.group(1)] = match.group(2).strip()
    return variables


def load_input_config(input_file: Path, restart_override: str | None) -> InputConfig:
    variables = read_index_variables(input_file)
    missing = [name for name in REQUIRED_VARIABLES if name not in variables]
    if missing:
        raise SystemExit(f"{input_file} is missing required variable(s): {', '.join(missing)}")

    return InputConfig(
        input_file=input_file,
        restart_file=restart_override or variables["restart_file"],
        target_t=variables["target_T"],
        force_mag=variables["force_mag"],
        protocol_tag=variables["protocol_tag"],
    )


def build_output_root_name(helper: ModuleType, config: InputConfig) -> str:
    restart_tag = helper.extract_restart_tag(config.restart_file)
    extra_tag = helper.sanitize_optional_tag(config.protocol_tag)
    output_root = f"output_restart{restart_tag}_Tstar_{config.target_t}_F{config.force_mag}"
    if extra_tag:
        output_root = f"{output_root}_{extra_tag}"
    return output_root


def resolve_restart_file(restart_file: str, workdir: Path) -> str:
    path = Path(restart_file)
    if path.is_absolute():
        return str(path)
    return str((workdir / path).resolve())


def build_run_configs(mpi_min: int, mpi_max: int, omp_min: int, omp_max: int) -> list[RunConfig]:
    if mpi_min <= 0 or omp_min <= 0 or mpi_max < mpi_min or omp_max < omp_min:
        raise SystemExit("Invalid MPI/OMP range. Require 1 <= min <= max.")
    return [
        RunConfig(mpi_ranks=mpi_ranks, omp_threads=omp_threads)
        for mpi_ranks in range(mpi_min, mpi_max + 1)
        for omp_threads in range(omp_min, omp_max + 1)
    ]


def build_command(
    *,
    run_config: RunConfig,
    input_file: Path,
    workdir: Path,
    helper_root: Path,
    restart_file: str,
    run_time: float,
    sample_every: int,
    dump_every: int,
    restart_every: int,
    mpiexec: str,
    lammps_bin: str,
) -> list[str]:
    return [
        mpiexec,
        "-np",
        str(run_config.mpi_ranks),
        lammps_bin,
        "-sf",
        "omp",
        "-pk",
        "omp",
        str(run_config.omp_threads),
        "-var",
        "helper_root",
        str(helper_root.resolve()),
        "-var",
        "output_base",
        str(workdir.resolve()),
        "-var",
        "restart_file",
        restart_file,
        "-var",
        "run_time",
        f"{run_time:g}",
        "-var",
        "sample_every",
        str(sample_every),
        "-var",
        "dump_every",
        str(dump_every),
        "-var",
        "restart_every",
        str(restart_every),
        "-in",
        str(input_file.resolve()),
    ]


def parse_log_metrics(log_path: Path) -> tuple[int, float] | None:
    if not log_path.exists():
        return None
    text = log_path.read_text(encoding="utf-8", errors="replace")
    matches = LOOP_TIME_PATTERN.findall(text)
    if not matches:
        return None
    total_loop_time = sum(float(loop_time) for loop_time, _ in matches)
    total_steps = sum(int(steps) for _, steps in matches)
    return total_steps, total_loop_time


def write_results_csv(results: list[RunResult], output_path: Path) -> None:
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "config",
                "mpi_ranks",
                "omp_threads",
                "return_code",
                "total_steps",
                "loop_time_s",
                "timesteps_per_second",
                "stdout_log",
                "lammps_log",
            ]
        )
        for result in results:
            writer.writerow(
                [
                    result.config.label,
                    result.config.mpi_ranks,
                    result.config.omp_threads,
                    result.return_code,
                    result.total_steps,
                    f"{result.total_loop_time:.6f}",
                    f"{result.timesteps_per_second:.6f}",
                    str(result.stdout_log),
                    str(result.lammps_log),
                ]
            )


def run_single_config(
    *,
    run_config: RunConfig,
    benchmark_root: Path,
    helper: ModuleType,
    input_config: InputConfig,
    input_file: Path,
    helper_root: Path,
    source_workdir: Path,
    run_time: float,
    sample_every: int,
    dump_every: int,
    restart_every: int,
    mpiexec: str,
    lammps_bin: str,
) -> RunResult:
    run_dir = benchmark_root / run_config.label
    run_dir.mkdir(parents=True, exist_ok=True)
    output_root = run_dir / build_output_root_name(helper, input_config)
    stdout_log = run_dir / "stdout.log"
    lammps_log = output_root / "log.force_clamp_aligned_box.lammps"
    command = build_command(
        run_config=run_config,
        input_file=input_file,
        workdir=run_dir,
        helper_root=helper_root,
        restart_file=resolve_restart_file(input_config.restart_file, source_workdir),
        run_time=run_time,
        sample_every=sample_every,
        dump_every=dump_every,
        restart_every=restart_every,
        mpiexec=mpiexec,
        lammps_bin=lammps_bin,
    )

    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(run_config.omp_threads)

    with stdout_log.open("w", encoding="utf-8") as handle:
        handle.write(f"Working directory: {run_dir}\n")
        handle.write(f"Command: {' '.join(shlex.quote(part) for part in command)}\n")
        handle.write(f"OMP_NUM_THREADS={run_config.omp_threads}\n\n")
        handle.flush()
        process = subprocess.run(
            command,
            cwd=run_dir,
            env=env,
            stdout=handle,
            stderr=subprocess.STDOUT,
            text=True,
            check=False,
        )

    metrics = parse_log_metrics(lammps_log)
    if metrics is None:
        metrics = parse_log_metrics(stdout_log)

    if metrics is None:
        total_steps = 0
        total_loop_time = 0.0
        timesteps_per_second = 0.0
    else:
        total_steps, total_loop_time = metrics
        timesteps_per_second = total_steps / total_loop_time if total_loop_time > 0 else 0.0

    return RunResult(
        config=run_config,
        return_code=process.returncode,
        total_steps=total_steps,
        total_loop_time=total_loop_time,
        timesteps_per_second=timesteps_per_second,
        stdout_log=stdout_log,
        lammps_log=lammps_log,
    )


def main() -> int:
    args = parse_args()
    validate_environment(dry_run=args.dry_run, mpiexec=args.mpiexec, lammps_bin=args.lammps_bin)

    input_file = args.input.expanduser().resolve()
    if not input_file.exists():
        raise SystemExit(f"Input file not found: {input_file}")

    workdir = args.workdir.expanduser().resolve()
    helper_root = args.helper_root.expanduser().resolve()
    helper = load_output_root_helper(helper_root)
    input_config = load_input_config(input_file, args.restart_file)
    run_configs = build_run_configs(args.mpi_min, args.mpi_max, args.omp_min, args.omp_max)
    benchmark_root = workdir / f"short_benchmark_{input_file.stem}"
    benchmark_root.mkdir(parents=True, exist_ok=True)

    print(f"Input file        : {input_file}")
    print(f"Working directory : {workdir}")
    print(f"Benchmark root    : {benchmark_root}")
    print(f"Helper root       : {helper_root}")
    print(f"Restart file      : {input_config.restart_file}")
    print(f"Run time          : {args.run_time:g}")
    print(f"Sample every      : {args.sample_every}")
    print(f"Dump every        : {args.dump_every}")
    print(f"Restart every     : {args.restart_every}")
    print(f"MPI range         : {args.mpi_min}..{args.mpi_max}")
    print(f"OMP range         : {args.omp_min}..{args.omp_max}")
    print(f"Max parallel      : {args.max_parallel}")
    print()

    if args.dry_run:
        for run_config in run_configs:
            run_dir = benchmark_root / run_config.label
            output_root = run_dir / build_output_root_name(helper, input_config)
            command = " ".join(
                shlex.quote(part)
                for part in build_command(
                    run_config=run_config,
                    input_file=input_file,
                    workdir=run_dir,
                    helper_root=helper_root,
                    restart_file=resolve_restart_file(input_config.restart_file, workdir),
                    run_time=args.run_time,
                    sample_every=args.sample_every,
                    dump_every=args.dump_every,
                    restart_every=args.restart_every,
                    mpiexec=args.mpiexec,
                    lammps_bin=args.lammps_bin,
                )
            )
            print(f"[PLAN] {run_config.label}")
            print(f"       run_dir    : {run_dir}")
            print(f"       output_root: {output_root}")
            print(f"       command    : {command}")
        return 0

    results: list[RunResult] = []
    max_parallel = max(1, args.max_parallel)
    running: list[subprocess.Popen[str]] = []
    del running  # explicit: benchmark is intentionally sequential unless max_parallel is expanded later
    if max_parallel != 1:
        print("Warning: current implementation runs sequentially; max_parallel is accepted for forward compatibility only.")

    for run_config in run_configs:
        print(f"[RUN ] {run_config.label}")
        result = run_single_config(
            run_config=run_config,
            benchmark_root=benchmark_root,
            helper=helper,
            input_config=input_config,
            input_file=input_file,
            helper_root=helper_root,
            source_workdir=workdir,
            run_time=args.run_time,
            sample_every=args.sample_every,
            dump_every=args.dump_every,
            restart_every=args.restart_every,
            mpiexec=args.mpiexec,
            lammps_bin=args.lammps_bin,
        )
        results.append(result)
        status = "OK" if result.return_code == 0 else f"FAIL({result.return_code})"
        print(
            f"[DONE] {run_config.label}: {status}, "
            f"{result.timesteps_per_second:.2f} steps/s, log={result.stdout_log}"
        )

    results.sort(key=lambda item: item.timesteps_per_second, reverse=True)
    csv_path = benchmark_root / "benchmark_results.csv"
    write_results_csv(results, csv_path)

    print("\nResults")
    print("-------")
    for result in results:
        status = "OK" if result.return_code == 0 else f"FAIL({result.return_code})"
        print(
            f"{result.config.label:>10}  {status:>8}  "
            f"{result.timesteps_per_second:10.2f} steps/s  "
            f"steps={result.total_steps:>8}  "
            f"loop={result.total_loop_time:9.6f}s"
        )

    print(f"\nCSV summary: {csv_path}")
    if any(result.return_code != 0 for result in results):
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
