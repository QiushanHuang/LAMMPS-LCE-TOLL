#!/usr/bin/env python3

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import os
from pathlib import Path
import re
import subprocess
import sys

from analyze_stretch_single import (
    ANALYSIS_CACHE_DIRNAME,
    WAIT_LIST_DIRNAME,
    dataset_filename_from_output_dir,
    ensure_analysis_root,
)


DEFAULT_DATA_FILENAME = "force_clamp_response.dat"
OUTPUT_PREFIX = "output_"
SCRIPT_BY_VARIANT = {
    "base": "analyze_stretch_single.py",
    "loess": "analyze_stretch_single_loess.py",
    "savgol": "analyze_stretch_single_savgol.py",
}


def natural_sort_key(text: str) -> list[object]:
    parts = re.split(r"(\d+)", text)
    key: list[object] = []
    for part in parts:
        if part.isdigit():
            key.append(int(part))
        else:
            key.append(part.lower())
    return key


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Batch-run cached stretch/Rg analysis across output_* directories."
    )
    parser.add_argument(
        "--root",
        type=Path,
        default=Path.cwd(),
        help="Root directory containing output_* subdirectories (default: current working directory)",
    )
    parser.add_argument(
        "--variant",
        choices=tuple(SCRIPT_BY_VARIANT),
        default="base",
        help="Which stretch analyzer entrypoint to run",
    )
    parser.add_argument(
        "--show",
        action="store_true",
        help="Pass --show to each per-directory stretch analyzer; windows open sequentially",
    )
    parser.add_argument(
        "--python-bin",
        default=sys.executable,
        help="Python interpreter used to launch the per-directory analyzer",
    )
    parser.add_argument(
        "--jobs",
        type=int,
        default=None,
        help="Parallel directory count for batch analysis (default: auto; forced to 1 with --show)",
    )
    parser.add_argument(
        "--analysis-jobs",
        type=int,
        default=None,
        help="Parallel dump-parser worker count forwarded to each stretch analyzer (default: auto when batch jobs=1, else 1)",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print planned commands without executing them",
    )
    return parser.parse_args()


def discover_output_dirs(root: Path) -> tuple[list[Path], list[Path]]:
    matched: list[Path] = []
    skipped: list[Path] = []
    for child in sorted(root.iterdir(), key=lambda path: natural_sort_key(path.name)):
        if not child.is_dir() or not child.name.startswith(OUTPUT_PREFIX):
            continue
        if (child / DEFAULT_DATA_FILENAME).exists():
            matched.append(child)
        else:
            skipped.append(child)
    return matched, skipped


def is_waitlisted(output_dir: Path) -> bool:
    wait_list_path = output_dir.parent / ANALYSIS_CACHE_DIRNAME / WAIT_LIST_DIRNAME / dataset_filename_from_output_dir(output_dir)
    return wait_list_path.exists()


def resolve_max_workers(requested_workers: int | None, task_count: int) -> int:
    if task_count <= 0:
        return 1
    if requested_workers is not None:
        return max(1, min(int(requested_workers), task_count))
    return max(1, min(len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else (os.cpu_count() or 1), task_count))


def build_command(script_path: Path, python_bin: str, show: bool, analysis_jobs: int | None) -> list[str]:
    command = [python_bin, str(script_path)]
    if show:
        command.append("--show")
    if analysis_jobs is not None:
        command.extend(["--jobs", str(analysis_jobs)])
    return command


def run_output_dir_task(task: tuple[list[str], str]) -> dict[str, str | int]:
    command, cwd = task
    completed = subprocess.run(
        command,
        cwd=cwd,
        capture_output=True,
        text=True,
    )
    return {
        "cwd": cwd,
        "returncode": completed.returncode,
        "stdout": completed.stdout,
        "stderr": completed.stderr,
    }


def main() -> int:
    args = parse_args()
    root = args.root.expanduser().resolve()
    if not root.is_dir():
        raise SystemExit(f"Root directory does not exist: {root}")
    analysis_root = ensure_analysis_root(root)

    script_path = Path(__file__).resolve().parent / SCRIPT_BY_VARIANT[args.variant]
    if not script_path.is_file():
        raise SystemExit(f"Analyzer script not found: {script_path}")

    output_dirs, skipped_dirs = discover_output_dirs(root)
    if not output_dirs:
        raise SystemExit(f"No output_* directories with {DEFAULT_DATA_FILENAME} found under {root}")

    print(f"Root directory : {root}")
    print(f"Analysis root  : {analysis_root}")
    print(f"Analyzer       : {script_path.name}")
    print(f"Directories    : {len(output_dirs)} runnable, {len(skipped_dirs)} skipped")
    if skipped_dirs:
        print("Skipped:")
        for skipped_dir in skipped_dirs:
            print(f"  - {skipped_dir}")

    batch_worker_count = resolve_max_workers(args.jobs, len(output_dirs))
    if args.show and batch_worker_count != 1:
        print("[INFO] --show requested; forcing serial batch execution.")
        batch_worker_count = 1

    if args.analysis_jobs is not None:
        analysis_jobs = max(1, args.analysis_jobs)
    elif batch_worker_count == 1:
        analysis_jobs = None
    else:
        analysis_jobs = 1

    success_count = 0
    failure_count = 0
    planned_tasks: list[tuple[Path, list[str]]] = []
    waitlisted_dirs: list[Path] = []
    for output_dir in output_dirs:
        if is_waitlisted(output_dir):
            waitlisted_dirs.append(output_dir)
            continue
        command = build_command(script_path, args.python_bin, args.show, analysis_jobs)
        printable = " ".join(command)
        print(f"[RUN ] {output_dir.name}")
        print(f"       cwd={output_dir}")
        print(f"       cmd={printable}")
        planned_tasks.append((output_dir, command))
        if args.dry_run:
            continue

    if waitlisted_dirs:
        print("Wait-listed:")
        for output_dir in waitlisted_dirs:
            print(f"  - {output_dir.name}")

    if args.dry_run:
        print(f"[SUMMARY] Planned {len(planned_tasks)} directories")
        return 0

    if batch_worker_count == 1:
        for output_dir, command in planned_tasks:
            completed = subprocess.run(command, cwd=output_dir)
            if completed.returncode == 0:
                success_count += 1
                print(f"[DONE] {output_dir.name}: PASS")
            else:
                failure_count += 1
                print(f"[DONE] {output_dir.name}: FAIL({completed.returncode})")
    else:
        task_specs = [(command, str(output_dir)) for output_dir, command in planned_tasks]
        with ThreadPoolExecutor(max_workers=batch_worker_count) as executor:
            future_to_name = {
                executor.submit(run_output_dir_task, task_spec): Path(task_spec[1]).name
                for task_spec in task_specs
            }
            for future in as_completed(future_to_name):
                output_name = future_to_name[future]
                result = future.result()
                if result["stdout"]:
                    print(result["stdout"], end="" if str(result["stdout"]).endswith("\n") else "\n")
                if result["stderr"]:
                    print(result["stderr"], file=sys.stderr, end="" if str(result["stderr"]).endswith("\n") else "\n")
                if int(result["returncode"]) == 0:
                    success_count += 1
                    print(f"[DONE] {output_name}: PASS")
                else:
                    failure_count += 1
                    print(f"[DONE] {output_name}: FAIL({result['returncode']})")

    print(
        f"[SUMMARY] success={success_count} failed={failure_count} skipped={len(skipped_dirs) + len(waitlisted_dirs)} total={len(output_dirs)}"
    )
    return 1 if failure_count else 0


if __name__ == "__main__":
    raise SystemExit(main())
