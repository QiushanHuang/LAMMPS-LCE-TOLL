#!/usr/bin/env python3

from __future__ import annotations

import argparse
import importlib.util
import re
import shlex
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from types import ModuleType


TMUX_BIN = "tmux"
MPIEXEC = "mpiexec"
LAMMPS_BIN = "/home/hkust/mylammps/build/lmp"
MPI_RANKS = 2
OMP_THREADS = 2
INPUT_PATTERN = "*.lmp"
TMUX_PREFIX = "qiushan"
SKIP_EXISTING_RESULT_DIR = True
SKIP_EXISTING_TMUX_SESSION = True

REQUIRED_VARIABLES = ("restart_file", "target_T", "force_mag", "protocol_tag")
VARIABLE_LINE = re.compile(r"^\s*variable\s+(\w+)\s+index\s+(.+?)\s*$")


@dataclass
class InputConfig:
    input_file: Path
    restart_file: str
    target_t: str
    force_mag: str
    protocol_tag: str
    output_root: Path


@dataclass
class LaunchPlan:
    config: InputConfig
    workdir: Path
    session_name: str
    window_name: str
    command: str
    launch_script: Path
    log_file: Path


def load_output_root_helper(helper_root: Path) -> ModuleType:
    helper_path = helper_root / "build_output_root.py"
    spec = importlib.util.spec_from_file_location("build_output_root", helper_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load helper module: {helper_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def natural_sort_key(path: Path) -> list[object]:
    parts = re.split(r"(\d+)", path.stem)
    key: list[object] = []
    for part in parts:
        if part.isdigit():
            key.append(int(part))
        else:
            key.append(part.lower())
    return key


def discover_inputs(script_dir: Path, pattern: str) -> list[Path]:
    return sorted(
        (path for path in script_dir.glob(pattern) if path.is_file()),
        key=natural_sort_key,
    )


def read_index_variables(input_file: Path) -> dict[str, str]:
    variables: dict[str, str] = {}
    for raw_line in input_file.read_text(encoding="utf-8").splitlines():
        line = raw_line.split("#", 1)[0].rstrip()
        if not line:
            continue
        match = VARIABLE_LINE.match(line)
        if not match:
            continue
        name = match.group(1)
        value = match.group(2).strip()
        variables[name] = value
    return variables


def resolve_output_base(raw_output_base: str, workdir: Path) -> Path:
    output_base = Path(raw_output_base)
    if output_base.is_absolute():
        return output_base
    return (workdir / output_base).resolve()


def build_output_root_name(
    helper: ModuleType,
    restart_file: str,
    target_t: str,
    force_mag: str,
    protocol_tag: str,
) -> str:
    restart_tag = helper.extract_restart_tag(restart_file)
    extra_tag = helper.sanitize_optional_tag(protocol_tag)
    output_root = f"output_restart{restart_tag}_Tstar_{target_t}_F{force_mag}"
    if extra_tag:
        output_root = f"{output_root}_{extra_tag}"
    return output_root


def build_input_config(
    input_file: Path,
    workdir: Path,
    helper: ModuleType,
    restart_override: str | None,
) -> InputConfig:
    variables = read_index_variables(input_file)
    missing = [name for name in REQUIRED_VARIABLES if name not in variables]
    if missing:
        raise ValueError(f"{input_file.name} is missing required variable(s): {', '.join(missing)}")

    restart_file = restart_override or variables["restart_file"]
    target_t = variables["target_T"]
    force_mag = variables["force_mag"]
    protocol_tag = variables["protocol_tag"]
    output_base = resolve_output_base(variables.get("output_base", "."), workdir)
    output_root = output_base / build_output_root_name(
        helper=helper,
        restart_file=restart_file,
        target_t=target_t,
        force_mag=force_mag,
        protocol_tag=protocol_tag,
    )

    return InputConfig(
        input_file=input_file,
        restart_file=restart_file,
        target_t=target_t,
        force_mag=force_mag,
        protocol_tag=protocol_tag,
        output_root=output_root,
    )


def session_exists(tmux_bin: str, session_name: str) -> bool:
    result = subprocess.run(
        [tmux_bin, "has-session", "-t", session_name],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        check=False,
    )
    return result.returncode == 0


def make_tmux_name(prefix: str, input_file: Path) -> str:
    stem = input_file.stem
    stem = re.sub(r"^in\.single_chain_constant_force_aligned_box-?", "", stem)
    sanitized = re.sub(r"[^A-Za-z0-9_-]+", "_", stem).strip("_")
    if not sanitized:
        sanitized = "job"
    max_suffix_len = 64
    if len(sanitized) > max_suffix_len:
        sanitized = sanitized[:max_suffix_len].rstrip("_")
    return f"{prefix}_{sanitized}"


def ensure_unique_session_name(tmux_bin: str, base_name: str) -> str:
    if not session_exists(tmux_bin, base_name):
        return base_name

    suffix = 2
    while True:
        candidate = f"{base_name}_{suffix}"
        if not session_exists(tmux_bin, candidate):
            return candidate
        suffix += 1


def build_inner_command(
    plan: LaunchPlan,
    mpiexec: str,
    mpi_ranks: int,
    lammps_bin: str,
    omp_threads: int,
    restart_override: str | None,
    helper_root: Path,
) -> str:
    command_parts = [
        mpiexec,
        "-np",
        str(mpi_ranks),
        lammps_bin,
        "-sf",
        "omp",
        "-pk",
        "omp",
        str(omp_threads),
    ]
    command_parts.extend(["-var", "helper_root", str(helper_root.resolve())])
    command_parts.extend(["-var", "output_base", str(plan.workdir.resolve())])
    if restart_override:
        command_parts.extend(["-var", "restart_file", restart_override])
    command_parts.extend(["-in", str(plan.config.input_file.resolve())])
    command_text = " ".join(shlex.quote(part) for part in command_parts)
    return (
        "set -euo pipefail\n"
        f"cd {shlex.quote(str(plan.workdir))}\n"
        f"export OMP_NUM_THREADS={omp_threads}\n"
        f"{command_text} 2>&1 | tee {shlex.quote(str(plan.log_file))}\n"
    )


def validate_environment(
    *,
    skip_runtime_checks: bool,
    tmux_bin: str,
    mpiexec: str,
    lammps_bin: str,
) -> None:
    if skip_runtime_checks:
        return

    missing = []
    for binary in (tmux_bin, mpiexec, lammps_bin):
        if shutil.which(binary) is None:
            missing.append(binary)

    if missing:
        raise SystemExit(f"Cannot find required executable(s): {', '.join(missing)}")


def write_launch_script(
    plan: LaunchPlan,
    mpiexec: str,
    mpi_ranks: int,
    lammps_bin: str,
    omp_threads: int,
    helper_root: Path,
    restart_override: str | None,
) -> None:
    plan.config.output_root.mkdir(parents=True, exist_ok=True)
    plan.launch_script.write_text(plan.command, encoding="utf-8")
    plan.launch_script.chmod(0o755)
    command_txt = (
        f"cd {plan.workdir}\n"
        f"OMP_NUM_THREADS={omp_threads} "
        f"{mpiexec} -np {mpi_ranks} {lammps_bin} -sf omp -pk omp {omp_threads} "
        f"-var helper_root {helper_root.resolve()} "
        f"-var output_base {plan.workdir.resolve()} "
        f"{f'-var restart_file {restart_override} ' if restart_override else ''}"
        f"-in {plan.config.input_file.resolve()}\n"
    )
    (plan.config.output_root / "command.txt").write_text(command_txt, encoding="utf-8")


def launch_tmux_session(tmux_bin: str, plan: LaunchPlan) -> None:
    tmux_command = f"bash -lc {shlex.quote(f'bash {plan.launch_script.name}')}"
    subprocess.run(
        [
            tmux_bin,
            "new-session",
            "-d",
            "-s",
            plan.session_name,
            "-n",
            plan.window_name,
            tmux_command,
        ],
        cwd=plan.config.output_root,
        check=True,
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Launch each packaged stretching LAMMPS input in its own detached tmux session. "
            "The script discovers .lmp files next to itself, but runs them from the current working "
            "directory so restart_file and output paths resolve there."
        )
    )
    parser.add_argument("--dry-run", action="store_true", help="Only print the planned jobs and commands.")
    parser.add_argument(
        "--workdir",
        type=Path,
        default=Path.cwd(),
        help="Working directory that contains restart files and will receive outputs (default: current directory).",
    )
    parser.add_argument(
        "--script-dir",
        type=Path,
        default=Path(__file__).resolve().parent,
        help="Directory containing the packaged stretching .lmp files (default: this script's directory).",
    )
    parser.add_argument("--pattern", default=INPUT_PATTERN, help=f"Glob used to discover input files (default: {INPUT_PATTERN}).")
    parser.add_argument(
        "--restart-file",
        default=None,
        help="Override restart_file for every launched job. Useful when all scripts should read the same relaxed restart.",
    )
    parser.add_argument("--mpi-ranks", type=int, default=MPI_RANKS, help=f"MPI ranks per job (default: {MPI_RANKS}).")
    parser.add_argument("--omp-threads", type=int, default=OMP_THREADS, help=f"OpenMP threads per job (default: {OMP_THREADS}).")
    parser.add_argument("--tmux-bin", default=TMUX_BIN, help=f"tmux executable name/path (default: {TMUX_BIN}).")
    parser.add_argument("--mpiexec", default=MPIEXEC, help=f"MPI launcher executable (default: {MPIEXEC}).")
    parser.add_argument("--lammps-bin", default=LAMMPS_BIN, help=f"LAMMPS executable (default: {LAMMPS_BIN}).")
    parser.add_argument(
        "--helper-root",
        type=Path,
        default=Path(__file__).resolve().parent.parent / "helpers",
        help="Directory containing build_output_root.py (default: sibling ../helpers relative to this script).",
    )
    parser.add_argument("--tmux-prefix", default=TMUX_PREFIX, help=f"Prefix used in tmux session names (default: {TMUX_PREFIX}).")
    parser.add_argument(
        "--rerun-existing-output",
        action="store_true",
        help="Do not skip jobs whose output_root directory already exists.",
    )
    parser.add_argument(
        "--reuse-session-names",
        action="store_true",
        help="Do not skip jobs when the default tmux session name already exists; append a numeric suffix instead.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    workdir = args.workdir.resolve()
    script_dir = args.script_dir.resolve()
    helper_root = args.helper_root.resolve()
    helper = load_output_root_helper(helper_root)
    validate_environment(
        skip_runtime_checks=args.dry_run,
        tmux_bin=args.tmux_bin,
        mpiexec=args.mpiexec,
        lammps_bin=args.lammps_bin,
    )

    inputs = discover_inputs(script_dir, args.pattern)
    if not inputs:
        print(f"No input files matched {args.pattern!r} in {script_dir}")
        return 1

    print(f"Script directory : {script_dir}")
    print(f"Working directory: {workdir}")
    print(f"Found input files: {len(inputs)}")
    print(f"MPI ranks/job    : {args.mpi_ranks}")
    print(f"OMP threads/job  : {args.omp_threads}")
    if args.restart_file:
        print(f"Restart override : {args.restart_file}")
    print(f"tmux prefix      : {args.tmux_prefix}")
    print()

    launched = 0
    skipped = 0
    failed = 0
    seen_output_roots: dict[Path, Path] = {}

    for input_file in inputs:
        try:
            config = build_input_config(
                input_file=input_file,
                workdir=workdir,
                helper=helper,
                restart_override=args.restart_file,
            )
        except Exception as exc:
            print(f"[ERR ] {input_file.name}: {exc}")
            failed += 1
            continue

        previous_owner = seen_output_roots.get(config.output_root)
        if previous_owner is not None:
            print(
                f"[ERR ] {input_file.name}: output_root collision -> {config.output_root.name} "
                f"(already claimed by {previous_owner.name}); adjust protocol_tag or launch separately"
            )
            failed += 1
            continue
        seen_output_roots[config.output_root] = input_file

        session_name = make_tmux_name(args.tmux_prefix, input_file)
        if not args.dry_run:
            if session_exists(args.tmux_bin, session_name) and not args.reuse_session_names:
                print(f"[SKIP] {input_file.name}: tmux session exists -> {session_name}")
                skipped += 1
                continue
            if session_exists(args.tmux_bin, session_name) and args.reuse_session_names:
                session_name = ensure_unique_session_name(args.tmux_bin, session_name)

        if config.output_root.exists() and not args.rerun_existing_output:
            print(f"[SKIP] {input_file.name}: output_root exists -> {config.output_root.name}")
            skipped += 1
            continue

        plan = LaunchPlan(
            config=config,
            workdir=workdir,
            session_name=session_name,
            window_name=session_name,
            command="",
            launch_script=config.output_root / "launch_in_tmux.sh",
            log_file=config.output_root / "tmux.log",
        )
        plan.command = build_inner_command(
            plan=plan,
            mpiexec=args.mpiexec,
            mpi_ranks=args.mpi_ranks,
            lammps_bin=args.lammps_bin,
            omp_threads=args.omp_threads,
            restart_override=args.restart_file,
            helper_root=helper_root,
        )

        if args.dry_run:
            print(f"[PLAN] {input_file.name}")
            print(f"       tmux session: {plan.session_name}")
            print(f"       workdir     : {workdir}")
            print(f"       helper_root : {helper_root}")
            print(f"       output_root : {config.output_root}")
            print(f"       restart_file: {config.restart_file}")
            print(
                f"       command     : {args.mpiexec} -np {args.mpi_ranks} {args.lammps_bin} "
                f"-sf omp -pk omp {args.omp_threads} "
                f"-var helper_root {helper_root} -var output_base {workdir}"
                f"{f' -var restart_file {args.restart_file}' if args.restart_file else ''} "
                f"-in {input_file.resolve()}"
            )
            continue

        write_launch_script(
            plan=plan,
            mpiexec=args.mpiexec,
            mpi_ranks=args.mpi_ranks,
            lammps_bin=args.lammps_bin,
            omp_threads=args.omp_threads,
            helper_root=helper_root,
            restart_override=args.restart_file,
        )
        launch_tmux_session(args.tmux_bin, plan)
        launched += 1
        print(f"[ OK ] {input_file.name}: tmux={plan.session_name}, output={config.output_root.name}")

    print()
    print(f"Summary: launched={launched}, skipped={skipped}, failed={failed}")
    if failed:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
