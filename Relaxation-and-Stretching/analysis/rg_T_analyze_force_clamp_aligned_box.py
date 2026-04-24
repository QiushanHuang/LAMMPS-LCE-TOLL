#!/usr/bin/env python3

from __future__ import annotations

import argparse
from array import array
import json
import math
import multiprocessing as mp
import os
from pathlib import Path
import re
from concurrent.futures import ProcessPoolExecutor
import sys

from PIL import Image, ImageColor, ImageDraw, ImageFont


BASE_COLUMN_NAMES = [
    "step",
    "time",
    "L",
    "dL",
    "strain",
    "contour_fraction",
    "Fext",
    "Fchain",
    "Ftotal",
]

TAIL_DIAGNOSTIC_COLUMN_NAMES = [
    "x_tail",
    "vx_tail",
    "fx_total_tail",
    "fx_chain_tail",
]

SUPPORTED_COLUMN_LAYOUTS = {
    len(BASE_COLUMN_NAMES): BASE_COLUMN_NAMES,
    len(BASE_COLUMN_NAMES) + len(TAIL_DIAGNOSTIC_COLUMN_NAMES): BASE_COLUMN_NAMES
    + TAIL_DIAGNOSTIC_COLUMN_NAMES,
}

DEFAULT_DATA_FILENAME = "force_clamp_response.dat"
DUMP_GLOB = "traj.force_clamp_aligned.*.dump"
PLOT_FILE_PREFIX = "rg_T_"
ANALYSIS_CACHE_DIRNAME = "0_stretch_analysis"
WAIT_LIST_DIRNAME = "wait-list-file"
STRETCH_DATASET_PREFIX = "stretch_dataset_"
STRETCH_DATASET_COLUMNS = [
    "source_tag",
    "timestep",
    "time",
    "L",
    "dL",
    "strain",
    "contour_fraction",
    "Fext",
    "Fchain",
    "Ftotal",
    "Lx",
    "Rg",
    "Rg2",
    "Rg_x",
    "count",
    "total_mass",
    "coord_source",
    "dump_file",
]
STRETCH_DATASET_NUMERIC_COLUMNS = [
    "timestep",
    "time",
    "L",
    "dL",
    "strain",
    "contour_fraction",
    "Fext",
    "Fchain",
    "Ftotal",
    "Lx",
    "Rg",
    "Rg2",
    "Rg_x",
    "count",
    "total_mass",
]


def choose_marker_indices(num_points: int, max_markers: int = 1200) -> list[int]:
    if num_points <= 0:
        return []
    if num_points <= max_markers:
        return list(range(num_points))

    # Sample evenly across the full range so the tail of long trajectories
    # does not collapse into a single endpoint.
    indices = [
        round(marker_idx * (num_points - 1) / (max_markers - 1))
        for marker_idx in range(max_markers)
    ]
    deduped: list[int] = []
    for idx in indices:
        if not deduped or idx != deduped[-1]:
            deduped.append(idx)
    return deduped


def generate_ticks(min_value: float, max_value: float, num_ticks: int = 6) -> list[float]:
    if num_ticks < 2 or max_value == min_value:
        return [min_value, max_value]
    step = (max_value - min_value) / (num_ticks - 1)
    return [min_value + idx * step for idx in range(num_ticks)]


def median(values: list[float]) -> float:
    ordered = sorted(values)
    n = len(ordered)
    mid = n // 2
    if n % 2 == 1:
        return ordered[mid]
    return 0.5 * (ordered[mid - 1] + ordered[mid])


def maybe_import_matplotlib():
    try:
        import matplotlib.pyplot as plt  # type: ignore
    except ModuleNotFoundError:
        return None
    return plt


def _decode_force_clamp_line(raw_line: bytes, path: Path) -> str:
    sanitized = raw_line.replace(b"\x00", b"")
    try:
        return sanitized.decode("utf-8").strip()
    except UnicodeDecodeError as exc:
        raise ValueError(f"Could not decode {path} as UTF-8 text") from exc


def load_force_clamp_table(path: str | Path) -> dict[str, array]:
    path = Path(path)
    rows: list[list[float]] = []
    active_columns: list[str] | None = None

    with path.open("rb") as handle:
        for raw_line in handle:
            line = _decode_force_clamp_line(raw_line, path)
            if not line or line.startswith("#"):
                continue
            fields = [float(value) for value in line.split()]
            if active_columns is None:
                active_columns = SUPPORTED_COLUMN_LAYOUTS.get(len(fields))
            if active_columns is None or len(fields) != len(active_columns):
                supported = ", ".join(str(width) for width in sorted(SUPPORTED_COLUMN_LAYOUTS))
                raise ValueError(
                    f"Expected one of [{supported}] columns in {path}, got {len(fields)}: {line}"
                )
            rows.append(fields)

    if not rows:
        raise ValueError(f"No data rows found in {path}")

    return {
        name: array("d", (row[idx] for row in rows))
        for idx, name in enumerate(active_columns or BASE_COLUMN_NAMES)
    }


def has_tail_diagnostics(table: dict[str, array]) -> bool:
    return all(name in table for name in TAIL_DIAGNOSTIC_COLUMN_NAMES)


def _natural_sort_key(text: str) -> list[object]:
    parts = re.split(r"(\d+)", text)
    key: list[object] = []
    for part in parts:
        if part.isdigit():
            key.append(int(part))
        else:
            key.append(part.lower())
    return key


def resolve_max_workers(requested_workers: int | None, task_count: int) -> int:
    if task_count <= 0:
        return 1
    if requested_workers is not None:
        return max(1, min(int(requested_workers), task_count))
    cpu_count = os.cpu_count() or 1
    return max(1, min(cpu_count, task_count))


def create_process_pool(max_workers: int) -> ProcessPoolExecutor | None:
    try:
        if sys.platform != "win32":
            available_methods = mp.get_all_start_methods()
            if "fork" in available_methods:
                return ProcessPoolExecutor(
                    max_workers=max_workers,
                    mp_context=mp.get_context("fork"),
                )
        return ProcessPoolExecutor(max_workers=max_workers)
    except (OSError, PermissionError):
        return None


def analysis_root_for_parent(root_dir: Path) -> Path:
    return root_dir / ANALYSIS_CACHE_DIRNAME


def ensure_analysis_root(root_dir: Path) -> Path:
    analysis_root = analysis_root_for_parent(root_dir)
    analysis_root.mkdir(parents=True, exist_ok=True)
    (analysis_root / WAIT_LIST_DIRNAME).mkdir(parents=True, exist_ok=True)
    return analysis_root


def dataset_tag_from_output_dir(output_dir: Path) -> str:
    if output_dir.name.startswith("output_"):
        return output_dir.name[len("output_") :]
    return output_dir.name


def dataset_filename_from_output_dir(output_dir: Path) -> str:
    return f"{STRETCH_DATASET_PREFIX}{dataset_tag_from_output_dir(output_dir)}.dat"


def dataset_path_for_output_dir(output_dir: Path) -> Path:
    return analysis_root_for_parent(output_dir.parent) / dataset_filename_from_output_dir(output_dir)


def wait_list_path_for_output_dir(output_dir: Path) -> Path:
    return analysis_root_for_parent(output_dir.parent) / WAIT_LIST_DIRNAME / dataset_filename_from_output_dir(output_dir)


def dataset_cache_status(output_dir: Path) -> str:
    dataset_path = dataset_path_for_output_dir(output_dir)
    if dataset_path.exists():
        return "cached"
    if wait_list_path_for_output_dir(output_dir).exists():
        return "wait-listed"
    return "missing"


def discover_cached_dataset_paths(root_dir: Path) -> list[Path]:
    analysis_root = ensure_analysis_root(root_dir)
    dataset_paths = [
        path
        for path in analysis_root.iterdir()
        if path.is_file()
        and path.name.startswith(STRETCH_DATASET_PREFIX)
        and path.suffix == ".dat"
        and stretch_dataset_schema_matches(path)
    ]
    return sorted(dataset_paths, key=lambda path: _natural_sort_key(path.name))


def stretch_dataset_header_columns(path: Path) -> list[str] | None:
    try:
        with path.open("r", encoding="utf-8", errors="replace") as handle:
            for raw_line in handle:
                line = raw_line.strip()
                if not line:
                    continue
                if not line.startswith("#"):
                    return None
                return line[1:].strip().split()
    except FileNotFoundError:
        return None
    return None


def stretch_dataset_schema_matches(path: Path) -> bool:
    return stretch_dataset_header_columns(path) == STRETCH_DATASET_COLUMNS


def compute_mass_weighted_rg_components(
    positions: list[tuple[float, float, float]],
    masses: list[float],
) -> dict[str, float]:
    if len(positions) != len(masses):
        raise ValueError("positions and masses must have the same length")
    total_mass = sum(masses)
    if total_mass <= 0.0:
        raise ValueError("Total mass must be positive to compute Rg")

    x_cm = sum(mass * position[0] for position, mass in zip(positions, masses)) / total_mass
    y_cm = sum(mass * position[1] for position, mass in zip(positions, masses)) / total_mass
    z_cm = sum(mass * position[2] for position, mass in zip(positions, masses)) / total_mass

    rg2 = sum(
        mass
        * (
            (position[0] - x_cm) ** 2
            + (position[1] - y_cm) ** 2
            + (position[2] - z_cm) ** 2
        )
        for position, mass in zip(positions, masses)
    ) / total_mass
    rg_x2 = sum(mass * (position[0] - x_cm) ** 2 for position, mass in zip(positions, masses)) / total_mass
    x_values = [position[0] for position in positions]
    return {
        "Rg": math.sqrt(max(rg2, 0.0)),
        "Rg2": max(rg2, 0.0),
        "Rg_x": math.sqrt(max(rg_x2, 0.0)),
        "Rg_x2": max(rg_x2, 0.0),
        "Lx": max(x_values) - min(x_values),
    }


def compute_mass_weighted_rg(positions: list[tuple[float, float, float]], masses: list[float]) -> float:
    return compute_mass_weighted_rg_components(positions, masses)["Rg"]


def parse_rg_dump_frames(
    dump_path: Path,
    selected_types: set[int] | None = None,
) -> list[dict[str, float]]:
    if selected_types is None:
        selected_types = {1, 2}

    frames: list[dict[str, float]] = []
    with dump_path.open("r", encoding="utf-8", errors="replace") as handle:
        while True:
            line = handle.readline()
            if not line:
                break
            if line.strip() == "":
                continue
            if line.strip() != "ITEM: TIMESTEP":
                raise ValueError(f"{dump_path} has an unexpected line: {line.strip()}")

            timestep = float(handle.readline().strip())

            if handle.readline().strip() != "ITEM: NUMBER OF ATOMS":
                raise ValueError(f"{dump_path} is missing ITEM: NUMBER OF ATOMS")
            natoms = int(float(handle.readline().strip()))

            box_line = handle.readline().strip()
            if not box_line.startswith("ITEM: BOX BOUNDS"):
                raise ValueError(f"{dump_path} is missing ITEM: BOX BOUNDS")
            for _ in range(3):
                handle.readline()

            atom_line = handle.readline().strip()
            if not atom_line.startswith("ITEM: ATOMS"):
                raise ValueError(f"{dump_path} is missing ITEM: ATOMS")

            header = atom_line.split()[2:]
            col_index = {name: idx for idx, name in enumerate(header)}
            for required in ("type", "mass"):
                if required not in col_index:
                    raise ValueError(f"{dump_path} is missing required dump column: {required}")

            if all(column in col_index for column in ("xu", "yu", "zu")):
                x_col, y_col, z_col = "xu", "yu", "zu"
                coord_source = "unwrap"
            elif all(column in col_index for column in ("x", "y", "z")):
                x_col, y_col, z_col = "x", "y", "z"
                coord_source = "wrapped"
            else:
                raise ValueError(f"{dump_path} must contain either xu/yu/zu or x/y/z")

            positions: list[tuple[float, float, float]] = []
            masses: list[float] = []
            for _ in range(natoms):
                fields = handle.readline().split()
                particle_type = int(float(fields[col_index["type"]]))
                if particle_type not in selected_types:
                    continue
                positions.append(
                    (
                        float(fields[col_index[x_col]]),
                        float(fields[col_index[y_col]]),
                        float(fields[col_index[z_col]]),
                    )
                )
                masses.append(float(fields[col_index["mass"]]))

            if positions:
                metrics = compute_mass_weighted_rg_components(positions, masses)
                frames.append(
                    {
                        "step": timestep,
                        "rg": metrics["Rg"],
                        "rg2": metrics["Rg2"],
                        "rg_x": metrics["Rg_x"],
                        "rg_x2": metrics["Rg_x2"],
                        "Lx": metrics["Lx"],
                        "count": float(len(positions)),
                        "total_mass": float(sum(masses)),
                        "coord_source": coord_source,
                        "dump_file": str(dump_path),
                    }
                )

    return frames


def parse_rg_dump_file_task(task: tuple[str, tuple[int, ...]]) -> list[dict[str, float]]:
    dump_path, selected_types = task
    return parse_rg_dump_frames(Path(dump_path), selected_types=set(selected_types))


def interpolate_series(
    x_source: array | list[float],
    y_source: array | list[float],
    x_target: array | list[float],
) -> array:
    ordered = sorted((float(x_value), float(y_value)) for x_value, y_value in zip(x_source, y_source))
    if not ordered:
        raise ValueError("Cannot interpolate an empty source series")
    x_sorted = [item[0] for item in ordered]
    y_sorted = [item[1] for item in ordered]
    return array("d", _interpolate_sorted_series(x_sorted, y_sorted, [float(value) for value in x_target]))


def collect_rg_frame_rows(data_dir: Path, max_workers: int | None = None) -> list[dict[str, float | str]]:
    dump_files = sorted(data_dir.glob(DUMP_GLOB), key=lambda path: _natural_sort_key(path.name))
    if not dump_files:
        raise SystemExit(f"No dump files matched {DUMP_GLOB!r} in {data_dir}")

    step_to_frame: dict[float, dict[str, float | str]] = {}
    tasks = [(str(dump_path), (1, 2)) for dump_path in dump_files]
    worker_count = resolve_max_workers(max_workers, len(tasks))

    if worker_count == 1:
        results = map(parse_rg_dump_file_task, tasks)
        executor: ProcessPoolExecutor | None = None
    else:
        chunksize = max(1, len(tasks) // (worker_count * 4))
        executor = create_process_pool(worker_count)
        if executor is None:
            results = map(parse_rg_dump_file_task, tasks)
        else:
            results = executor.map(parse_rg_dump_file_task, tasks, chunksize=chunksize)

    try:
        for frames in results:
            for frame in frames:
                step_to_frame[float(frame["step"])] = frame
    finally:
        if executor is not None:
            executor.shutdown(wait=True)

    if not step_to_frame:
        raise SystemExit(f"No valid type-1/type-2 frames were found in dump files under {data_dir}")

    return [step_to_frame[step] for step in sorted(step_to_frame)]


def load_rg_kinematics(
    data_dir: Path,
    table: dict[str, array],
    max_workers: int | None = None,
) -> dict[str, array]:
    ordered_frames = collect_rg_frame_rows(data_dir, max_workers=max_workers)
    rg_step = array("d", (float(item["step"]) for item in ordered_frames))
    rg = array("d", (float(item["rg"]) for item in ordered_frames))
    rg_time = interpolate_series(table["step"], table["time"], rg_step)
    rg_length = interpolate_series(table["step"], table["L"], rg_step)
    rg_derivative_time = finite_difference(rg_time, rg)
    return {
        "step": rg_step,
        "time": rg_time,
        "length": rg_length,
        "length_x": array("d", (float(item["Lx"]) for item in ordered_frames)),
        "rg": rg,
        "rg_x": array("d", (float(item["rg_x"]) for item in ordered_frames)),
        "rg2": array("d", (float(item["rg2"]) for item in ordered_frames)),
        "count": array("d", (float(item["count"]) for item in ordered_frames)),
        "total_mass": array("d", (float(item["total_mass"]) for item in ordered_frames)),
        "coord_source": [str(item["coord_source"]) for item in ordered_frames],
        "dump_file": [str(item["dump_file"]) for item in ordered_frames],
        "rg_derivative_time": rg_derivative_time,
    }


def build_stretch_dataset_rows(
    output_dir: Path,
    force_table: dict[str, array],
    rg_frames: list[dict[str, float | str]],
) -> list[dict[str, float | str]]:
    source_tag = dataset_tag_from_output_dir(output_dir)
    rg_step = array("d", (float(item["step"]) for item in rg_frames))
    rows: list[dict[str, float | str]] = []

    interpolated = {
        column: interpolate_series(force_table["step"], force_table[column], rg_step)
        for column in BASE_COLUMN_NAMES[1:]
    }

    for index, rg_frame in enumerate(rg_frames):
        rows.append(
            {
                "source_tag": source_tag,
                "timestep": float(rg_frame["step"]),
                "time": float(interpolated["time"][index]),
                "L": float(interpolated["L"][index]),
                "dL": float(interpolated["dL"][index]),
                "strain": float(interpolated["strain"][index]),
                "contour_fraction": float(interpolated["contour_fraction"][index]),
                "Fext": float(interpolated["Fext"][index]),
                "Fchain": float(interpolated["Fchain"][index]),
                "Ftotal": float(interpolated["Ftotal"][index]),
                "Lx": float(rg_frame["Lx"]),
                "Rg": float(rg_frame["rg"]),
                "Rg2": float(rg_frame["rg2"]),
                "Rg_x": float(rg_frame["rg_x"]),
                "count": float(rg_frame["count"]),
                "total_mass": float(rg_frame["total_mass"]),
                "coord_source": str(rg_frame["coord_source"]),
                "dump_file": str(rg_frame["dump_file"]),
            }
        )
    return rows


def write_stretch_dataset(path: Path, rows: list[dict[str, float | str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        handle.write("# " + " ".join(STRETCH_DATASET_COLUMNS) + "\n")
        for row in rows:
            handle.write(
                (
                    f"{row['source_tag']} "
                    f"{int(round(float(row['timestep'])))} "
                    f"{float(row['time']):.8f} "
                    f"{float(row['L']):.8f} "
                    f"{float(row['dL']):.8f} "
                    f"{float(row['strain']):.8f} "
                    f"{float(row['contour_fraction']):.8f} "
                    f"{float(row['Fext']):.8f} "
                    f"{float(row['Fchain']):.8f} "
                    f"{float(row['Ftotal']):.8f} "
                    f"{float(row['Lx']):.8f} "
                    f"{float(row['Rg']):.8f} "
                    f"{float(row['Rg2']):.8f} "
                    f"{float(row['Rg_x']):.8f} "
                    f"{int(round(float(row['count'])))} "
                    f"{float(row['total_mass']):.8f} "
                    f"{row['coord_source']} "
                    f"{row['dump_file']}\n"
                )
            )


def load_stretch_dataset_rows(path: Path) -> list[dict[str, float | str]]:
    rows: list[dict[str, float | str]] = []
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) != len(STRETCH_DATASET_COLUMNS):
                raise ValueError(
                    f"Expected {len(STRETCH_DATASET_COLUMNS)} columns in {path}, got {len(fields)}: {line}"
                )
            row: dict[str, float | str] = {}
            for column, value in zip(STRETCH_DATASET_COLUMNS, fields):
                if column in STRETCH_DATASET_NUMERIC_COLUMNS:
                    row[column] = float(value)
                else:
                    row[column] = value
            rows.append(row)

    if not rows:
        raise ValueError(f"No data rows found in {path}")
    return rows


def dataset_rows_to_analysis(rows: list[dict[str, float | str]]) -> dict[str, object]:
    table: dict[str, array] = {}
    for column in BASE_COLUMN_NAMES:
        source_column = "timestep" if column == "step" else column
        table[column] = array("d", (float(row[source_column]) for row in rows))
    table["Lx"] = array("d", (float(row["Lx"]) for row in rows))
    rg = array("d", (float(row["Rg"]) for row in rows))
    rg_time = array("d", (float(row["time"]) for row in rows))
    rg_kinematics = {
        "step": table["step"],
        "time": rg_time,
        "length": table["L"],
        "length_x": table["Lx"],
        "rg": rg,
        "rg_x": array("d", (float(row["Rg_x"]) for row in rows)),
        "rg2": array("d", (float(row["Rg2"]) for row in rows)),
        "count": array("d", (float(row["count"]) for row in rows)),
        "total_mass": array("d", (float(row["total_mass"]) for row in rows)),
        "coord_source": [str(row["coord_source"]) for row in rows],
        "dump_file": [str(row["dump_file"]) for row in rows],
        "rg_derivative_time": finite_difference(rg_time, rg),
    }
    return {
        "source_tag": str(rows[0]["source_tag"]),
        "table": table,
        "rg_kinematics": rg_kinematics,
    }


def ensure_stretch_dataset_file(
    data_file: Path,
    max_workers: int | None = None,
) -> tuple[Path | None, str]:
    output_dir = data_file.parent
    analysis_root = ensure_analysis_root(output_dir.parent)
    dataset_path = analysis_root / dataset_filename_from_output_dir(output_dir)
    wait_list_path = wait_list_path_for_output_dir(output_dir)
    if dataset_path.exists() and stretch_dataset_schema_matches(dataset_path):
        return dataset_path, "reused"
    if wait_list_path.exists():
        return None, "wait-listed"
    if dataset_path.exists():
        dataset_path.unlink()

    force_table = load_force_clamp_table(data_file)
    rg_frames = collect_rg_frame_rows(output_dir, max_workers=max_workers)
    rows = build_stretch_dataset_rows(output_dir, force_table, rg_frames)
    write_stretch_dataset(dataset_path, rows)
    return dataset_path, "created"


def prepare_stretch_dataset(
    data_file: Path,
    max_workers: int | None = None,
) -> tuple[Path | None, dict[str, object] | None, str]:
    dataset_path, status = ensure_stretch_dataset_file(data_file, max_workers=max_workers)
    if status == "wait-listed" or dataset_path is None:
        return None, None, status
    return dataset_path, dataset_rows_to_analysis(load_stretch_dataset_rows(dataset_path)), status


def _tricube_weight(distance: float, bandwidth: float) -> float:
    if bandwidth <= 0.0:
        return 1.0 if distance == 0.0 else 0.0
    ratio = abs(distance) / bandwidth
    if ratio >= 1.0:
        return 0.0
    return (1.0 - ratio**3) ** 3


def _lowess_core(
    x_values: list[float], y_values: list[float], frac: float, robust_iters: int
) -> list[float]:
    n = len(x_values)
    if n <= 2:
        return list(y_values)

    window = max(2, min(n, math.ceil(frac * n)))
    robust_weights = [1.0] * n
    y_hat = [0.0] * n

    for iteration in range(robust_iters + 1):
        for idx, x0 in enumerate(x_values):
            distances = [abs(xj - x0) for xj in x_values]
            bandwidth = sorted(distances)[window - 1]

            weights = [
                _tricube_weight(distance, bandwidth) * robust_weights[j]
                for j, distance in enumerate(distances)
            ]

            sum_w = sum(weights)
            if sum_w == 0.0:
                y_hat[idx] = y_values[idx]
                continue

            sum_wx = sum(w * xv for w, xv in zip(weights, x_values))
            sum_wy = sum(w * yv for w, yv in zip(weights, y_values))
            sum_wxx = sum(w * xv * xv for w, xv in zip(weights, x_values))
            sum_wxy = sum(w * xv * yv for w, xv, yv in zip(weights, x_values, y_values))
            denominator = sum_w * sum_wxx - sum_wx * sum_wx

            if abs(denominator) < 1e-12:
                y_hat[idx] = sum_wy / sum_w
                continue

            slope = (sum_w * sum_wxy - sum_wx * sum_wy) / denominator
            intercept = (sum_wy - slope * sum_wx) / sum_w
            y_hat[idx] = intercept + slope * x0

        if iteration == robust_iters:
            break

        residuals = [abs(yv - yfit) for yv, yfit in zip(y_values, y_hat)]
        scale = median(residuals)
        if scale <= 1e-12:
            break

        cutoff = 6.0 * scale
        new_weights = []
        for residual in residuals:
            if residual >= cutoff:
                new_weights.append(0.0)
                continue
            ratio = residual / cutoff
            new_weights.append((1.0 - ratio * ratio) ** 2)
        robust_weights = new_weights

    return y_hat


def _interpolate_sorted_series(
    x_source: list[float], y_source: list[float], x_target: list[float]
) -> list[float]:
    if len(x_source) == 1:
        return [y_source[0]] * len(x_target)

    result: list[float] = []
    cursor = 0
    for x_value in x_target:
        while cursor < len(x_source) - 2 and x_source[cursor + 1] < x_value:
            cursor += 1

        x_left = x_source[cursor]
        x_right = x_source[cursor + 1]
        y_left = y_source[cursor]
        y_right = y_source[cursor + 1]

        if x_right == x_left:
            result.append(y_right)
            continue

        alpha = (x_value - x_left) / (x_right - x_left)
        result.append(y_left + alpha * (y_right - y_left))
    return result


def lowess_smooth(
    x: array | list[float],
    y: array | list[float],
    frac: float = 0.08,
    robust_iters: int = 2,
    max_points: int = 2000,
) -> array:
    if len(x) != len(y):
        raise ValueError("x and y must have the same length")
    if len(x) == 0:
        raise ValueError("x and y must not be empty")
    if not 0.0 < frac <= 1.0:
        raise ValueError("frac must satisfy 0 < frac <= 1")

    original = [(float(xv), float(yv), idx) for idx, (xv, yv) in enumerate(zip(x, y))]
    original.sort(key=lambda item: item[0])

    x_sorted = [item[0] for item in original]
    y_sorted = [item[1] for item in original]

    if len(original) > max_points:
        sample_indices = choose_marker_indices(len(original), max_markers=max_points)
        x_work = [x_sorted[idx] for idx in sample_indices]
        y_work = [y_sorted[idx] for idx in sample_indices]
        y_work_smooth = _lowess_core(x_work, y_work, frac=frac, robust_iters=robust_iters)
        y_sorted_smooth = _interpolate_sorted_series(x_work, y_work_smooth, x_sorted)
    else:
        y_sorted_smooth = _lowess_core(x_sorted, y_sorted, frac=frac, robust_iters=robust_iters)

    reordered = [0.0] * len(original)
    for (_, _, original_idx), y_value in zip(original, y_sorted_smooth):
        reordered[original_idx] = y_value
    return array("d", reordered)


def sorted_pairs(x: array | list[float], y: array | list[float]) -> tuple[array, array]:
    ordered = sorted((float(xv), float(yv)) for xv, yv in zip(x, y))
    return array("d", (pair[0] for pair in ordered)), array("d", (pair[1] for pair in ordered))


def make_ratio(numerator: array, denominator: array) -> array:
    values = []
    for num, den in zip(numerator, denominator):
        if abs(den) < 1.0e-12:
            values.append(0.0)
        else:
            values.append(float(num) / float(den))
    return array("d", values)


def negate(values: array | list[float]) -> array:
    return array("d", (-float(value) for value in values))


def scale_array(values: array | list[float], factor: float) -> array:
    return array("d", (float(value) * factor for value in values))


def infer_contour_length(table: dict[str, array]) -> float:
    inferred = [
        float(length) / float(contour_fraction)
        for length, contour_fraction in zip(table["L"], table["contour_fraction"])
        if abs(contour_fraction) > 1.0e-12
    ]
    if not inferred:
        raise ValueError("Cannot infer contour length because contour_fraction is zero everywhere")
    return median(inferred)


def derive_mechanics(table: dict[str, array]) -> dict[str, float | array]:
    initial_length = float(table["L"][0])
    contour_length = infer_contour_length(table)
    denominator = contour_length - initial_length
    if abs(denominator) < 1.0e-12:
        raise ValueError("Contour length is too close to initial length; cannot define strain")

    dL = array("d", (float(length) - initial_length for length in table["L"]))
    strain = array("d", (float(delta) / denominator for delta in dL))

    dL_mismatch = max(abs(float(a) - float(b)) for a, b in zip(table["dL"], dL))
    strain_mismatch = max(abs(float(a) - float(b)) for a, b in zip(table["strain"], strain))

    return {
        "initial_length": initial_length,
        "contour_length": contour_length,
        "dL": dL,
        "strain": strain,
        "max_dL_column_mismatch": dL_mismatch,
        "max_strain_column_mismatch": strain_mismatch,
    }


def finite_difference(x: array | list[float], y: array | list[float]) -> array:
    x_values = [float(value) for value in x]
    y_values = [float(value) for value in y]
    n = len(x_values)
    if n == 0:
        return array("d")
    if n == 1:
        return array("d", [0.0])

    derivatives = [0.0] * n
    for idx in range(n):
        if idx == 0:
            dx = x_values[1] - x_values[0]
            derivatives[idx] = 0.0 if abs(dx) < 1.0e-12 else (y_values[1] - y_values[0]) / dx
        elif idx == n - 1:
            dx = x_values[-1] - x_values[-2]
            derivatives[idx] = 0.0 if abs(dx) < 1.0e-12 else (y_values[-1] - y_values[-2]) / dx
        else:
            dx = x_values[idx + 1] - x_values[idx - 1]
            derivatives[idx] = (
                0.0 if abs(dx) < 1.0e-12 else (y_values[idx + 1] - y_values[idx - 1]) / dx
            )
    return array("d", derivatives)


def compute_length_kinematics(table: dict[str, array]) -> tuple[array, array]:
    length_derivative = finite_difference(table["time"], table["L"])
    length_second_derivative = finite_difference(table["time"], length_derivative)
    return length_derivative, length_second_derivative


def _solve_linear_system(matrix: list[list[float]], rhs: list[float]) -> list[float]:
    size = len(rhs)
    aug = [row[:] + [rhs[row_idx]] for row_idx, row in enumerate(matrix)]

    for col in range(size):
        pivot = max(range(col, size), key=lambda row_idx: abs(aug[row_idx][col]))
        if abs(aug[pivot][col]) < 1.0e-12:
            return [0.0] * size
        if pivot != col:
            aug[col], aug[pivot] = aug[pivot], aug[col]

        scale = aug[col][col]
        for idx in range(col, size + 1):
            aug[col][idx] /= scale

        for row_idx in range(size):
            if row_idx == col:
                continue
            factor = aug[row_idx][col]
            if factor == 0.0:
                continue
            for idx in range(col, size + 1):
                aug[row_idx][idx] -= factor * aug[col][idx]

    return [aug[row_idx][-1] for row_idx in range(size)]


def savitzky_golay_series(
    x: array | list[float],
    y: array | list[float],
    derivative_order: int = 0,
    window_size: int = 31,
    poly_order: int = 3,
) -> array:
    x_values = [float(value) for value in x]
    y_values = [float(value) for value in y]
    n = len(x_values)
    if n == 0:
        return array("d")
    if n == 1:
        return array("d", [0.0 if derivative_order > 0 else y_values[0]])

    window = min(max(window_size, 3), n)
    if window < derivative_order + 1:
        window = derivative_order + 1
    if window > n:
        window = n

    values: list[float] = []
    for idx, x0 in enumerate(x_values):
        half = window // 2
        start = max(0, idx - half)
        end = min(n, start + window)
        start = max(0, end - window)

        x_window = x_values[start:end]
        y_window = y_values[start:end]
        effective_order = min(poly_order, len(x_window) - 1)
        if effective_order < derivative_order:
            values.append(0.0)
            continue

        offsets = [xv - x0 for xv in x_window]
        degree_count = effective_order + 1
        normal = [[0.0] * degree_count for _ in range(degree_count)]
        rhs = [0.0] * degree_count

        for offset, yv in zip(offsets, y_window):
            powers = [1.0]
            for _ in range(effective_order):
                powers.append(powers[-1] * offset)
            for row in range(degree_count):
                rhs[row] += yv * powers[row]
                for col in range(degree_count):
                    normal[row][col] += powers[row] * powers[col]

        coeffs = _solve_linear_system(normal, rhs)
        if derivative_order == 0:
            values.append(coeffs[0])
        else:
            values.append(coeffs[derivative_order] * math.factorial(derivative_order))

    return array("d", values)


def build_summary(table: dict[str, array]) -> dict[str, float]:
    derived = derive_mechanics(table)
    resisting_force = negate(table["Fchain"])
    net_drag_force = array(
        "d", (float(fext) + float(fresist) for fext, fresist in zip(table["Fext"], resisting_force))
    )
    length_derivative, _ = compute_length_kinematics(table)
    extension_rate = scale_array(
        length_derivative,
        1.0 / (float(derived["contour_length"]) - float(derived["initial_length"])),
    )
    summary = {
        "force_basis": "x_projected",
        "prescribed_force": float(table["Fext"][0]),
        "initial_length": float(derived["initial_length"]),
        "inferred_contour_length": float(derived["contour_length"]),
        "final_length": float(table["L"][-1]),
        "final_extension": float(derived["dL"][-1]),
        "final_strain": float(derived["strain"][-1]),
        "final_contour_fraction": float(table["contour_fraction"][-1]),
        "final_resisting_force": float(resisting_force[-1]),
        "max_resisting_force": float(max(resisting_force)),
        "min_resisting_force": float(min(resisting_force)),
        "max_chain_force": float(max(resisting_force)),
        "min_chain_force": float(min(resisting_force)),
        "mean_chain_force": float(sum(resisting_force) / len(resisting_force)),
        "max_abs_total_force": float(max(abs(value) for value in table["Ftotal"])),
        "mean_total_force": float(sum(table["Ftotal"]) / len(table["Ftotal"])),
        "final_net_drag_force": float(net_drag_force[-1]),
        "max_net_drag_force": float(max(net_drag_force)),
        "min_net_drag_force": float(min(net_drag_force)),
        "max_length_derivative": float(max(length_derivative)),
        "min_length_derivative": float(min(length_derivative)),
        "max_extension_rate": float(max(extension_rate)),
        "min_extension_rate": float(min(extension_rate)),
        "max_dL_column_mismatch": float(derived["max_dL_column_mismatch"]),
        "max_strain_column_mismatch": float(derived["max_strain_column_mismatch"]),
    }
    if has_tail_diagnostics(table):
        tail_acceleration = finite_difference(table["time"], table["vx_tail"])
        summary.update(
            {
                "final_x_tail": float(table["x_tail"][-1]),
                "final_vx_tail": float(table["vx_tail"][-1]),
                "final_fx_total_tail": float(table["fx_total_tail"][-1]),
                "final_fx_chain_tail": float(table["fx_chain_tail"][-1]),
                "max_abs_vx_tail": float(max(abs(value) for value in table["vx_tail"])),
                "max_abs_fx_total_tail": float(max(abs(value) for value in table["fx_total_tail"])),
                "max_abs_fx_chain_tail": float(max(abs(value) for value in table["fx_chain_tail"])),
                "max_abs_tail_acceleration": float(max(abs(value) for value in tail_acceleration)),
            }
        )
    return summary


def _load_font(size: int) -> ImageFont.ImageFont:
    for candidate in ("DejaVuSans.ttf", "Arial.ttf", "Helvetica.ttf"):
        try:
            return ImageFont.truetype(candidate, size=size)
        except OSError:
            continue
    return ImageFont.load_default()


def _color_with_alpha(color: str, opacity: float = 1.0) -> tuple[int, int, int, int]:
    rgb = ImageColor.getrgb(color)
    alpha = max(0, min(255, int(round(255 * opacity))))
    return rgb[0], rgb[1], rgb[2], alpha


def _text_size(draw: ImageDraw.ImageDraw, text: str, font: ImageFont.ImageFont) -> tuple[int, int]:
    left, top, right, bottom = draw.textbbox((0, 0), text, font=font)
    return right - left, bottom - top


def _draw_text(
    draw: ImageDraw.ImageDraw,
    x: float,
    y: float,
    text: str,
    font: ImageFont.ImageFont,
    fill: tuple[int, int, int, int],
    anchor: str = "lt",
    angle: float = 0.0,
) -> None:
    text_width, text_height = _text_size(draw, text, font)
    x0 = x
    y0 = y
    if "m" in anchor:
        x0 -= text_width / 2.0
    elif "r" in anchor:
        x0 -= text_width
    if anchor.startswith("m"):
        y0 -= text_height / 2.0
    elif anchor.startswith("b"):
        y0 -= text_height
    if angle == 0.0:
        draw.text((x0, y0), text, font=font, fill=fill)
        return

    text_image = Image.new("RGBA", (max(1, text_width + 8), max(1, text_height + 8)), (255, 255, 255, 0))
    text_draw = ImageDraw.Draw(text_image, "RGBA")
    text_draw.text((4, 4), text, font=font, fill=fill)
    rotated = text_image.rotate(angle, expand=True)
    draw._image.alpha_composite(rotated, (int(round(x0)), int(round(y0))))  # type: ignore[attr-defined]


def _draw_line_segment(
    draw: ImageDraw.ImageDraw,
    start: tuple[float, float],
    end: tuple[float, float],
    color: tuple[int, int, int, int],
    width: int,
) -> None:
    draw.line((start[0], start[1], end[0], end[1]), fill=color, width=max(1, width))


def _draw_polyline(
    draw: ImageDraw.ImageDraw,
    points: list[tuple[float, float]],
    color: str,
    width: float = 2.0,
    opacity: float = 0.9,
    dasharray: str | None = None,
) -> None:
    if len(points) < 2:
        return
    rgba = _color_with_alpha(color, opacity)
    pixel_width = max(1, int(round(width)))

    if not dasharray:
        draw.line(points, fill=rgba, width=pixel_width)
        return

    dash_values = [max(1.0, float(value)) for value in dasharray.split()]
    if not dash_values:
        draw.line(points, fill=rgba, width=pixel_width)
        return

    for start, end in zip(points[:-1], points[1:]):
        dx = end[0] - start[0]
        dy = end[1] - start[1]
        segment_length = math.hypot(dx, dy)
        if segment_length <= 1.0e-12:
            continue
        ux = dx / segment_length
        uy = dy / segment_length
        cursor = 0.0
        dash_idx = 0
        draw_on = True
        while cursor < segment_length:
            piece = dash_values[dash_idx % len(dash_values)]
            next_cursor = min(segment_length, cursor + piece)
            if draw_on:
                p0 = (start[0] + ux * cursor, start[1] + uy * cursor)
                p1 = (start[0] + ux * next_cursor, start[1] + uy * next_cursor)
                _draw_line_segment(draw, p0, p1, rgba, pixel_width)
            cursor = next_cursor
            dash_idx += 1
            draw_on = not draw_on


def _draw_scatter(
    draw: ImageDraw.ImageDraw,
    points: list[tuple[float, float]],
    color: str,
    radius: float = 2.5,
    opacity: float = 0.9,
) -> None:
    rgba = _color_with_alpha(color, opacity)
    r = max(1.0, radius)
    for x_value, y_value in points:
        draw.ellipse((x_value - r, y_value - r, x_value + r, y_value + r), fill=rgba)


def _save_png(image: Image.Image, output_path: Path) -> None:
    image.convert("RGB").save(output_path, format="PNG")


def save_series_plot(
    x: array | list[float],
    series: list[dict[str, object]],
    xlabel: str,
    ylabel: str,
    title: str,
    output_path: Path,
    subtitle: str | None = None,
) -> None:
    width = 1600
    height = 980
    pad_left = 150
    pad_right = 60
    pad_top = 95
    pad_bottom = 120

    image = Image.new("RGBA", (width, height), (255, 255, 255, 255))
    draw = ImageDraw.Draw(image, "RGBA")
    font_title = _load_font(30)
    font_subtitle = _load_font(16)
    font_axis = _load_font(22)
    font_tick = _load_font(16)
    font_legend = _load_font(16)

    x_values = [float(value) for value in x]
    y_values: list[float] = []
    for item in series:
        y_values.extend(float(value) for value in item["y"])  # type: ignore[index]

    x_min = min(x_values)
    x_max = max(x_values)
    y_min = min(y_values)
    y_max = max(y_values)
    if x_max == x_min:
        x_max = x_min + 1.0
    if y_max == y_min:
        y_max = y_min + 1.0

    padding_y = 0.06 * (y_max - y_min) or 1.0
    y_min -= padding_y
    y_max += padding_y

    plot_width = width - pad_left - pad_right
    plot_height = height - pad_top - pad_bottom

    def sx(value: float) -> float:
        return pad_left + (value - x_min) / (x_max - x_min) * plot_width

    def sy(value: float) -> float:
        return height - pad_bottom - (value - y_min) / (y_max - y_min) * plot_height

    _draw_text(draw, width / 2.0, 34, title, font_title, _color_with_alpha("#000000"), anchor="mt")
    if subtitle:
        _draw_text(draw, width / 2.0, 63, subtitle, font_subtitle, _color_with_alpha("#555555"), anchor="mt")

    _draw_line_segment(draw, (pad_left, height - pad_bottom), (width - pad_right, height - pad_bottom), _color_with_alpha("#000000"), 2)
    _draw_line_segment(draw, (pad_left, pad_top), (pad_left, height - pad_bottom), _color_with_alpha("#000000"), 2)
    _draw_text(draw, width / 2.0, height - 28, xlabel, font_axis, _color_with_alpha("#000000"), anchor="mt")
    _draw_text(draw, 34, height / 2.0, ylabel, font_axis, _color_with_alpha("#000000"), anchor="mm", angle=90.0)

    for tick in generate_ticks(x_min, x_max):
        tx = sx(tick)
        _draw_line_segment(draw, (tx, pad_top), (tx, height - pad_bottom), _color_with_alpha("#d9d9d9"), 1)
        _draw_text(draw, tx, height - pad_bottom + 30, f"{tick:.4g}", font_tick, _color_with_alpha("#000000"), anchor="mt")
    for tick in generate_ticks(y_min, y_max):
        ty = sy(tick)
        _draw_line_segment(draw, (pad_left, ty), (width - pad_right, ty), _color_with_alpha("#d9d9d9"), 1)
        _draw_text(draw, pad_left - 14, ty, f"{tick:.4g}", font_tick, _color_with_alpha("#000000"), anchor="rm")

    for item in series:
        y_series = [float(value) for value in item["y"]]  # type: ignore[index]
        x_series = [float(value) for value in item.get("x", x)]  # type: ignore[arg-type]
        points = [(sx(xv), sy(yv)) for xv, yv in zip(x_series, y_series)]
        mode = str(item.get("mode", "line"))
        if mode == "scatter":
            marker_indices = choose_marker_indices(len(points), max_markers=1200)
            _draw_scatter(
                draw,
                [points[idx] for idx in marker_indices],
                color=str(item.get("color", "#1f77b4")),
                radius=2.3,
                opacity=float(item.get("opacity", 0.9)),
            )
        else:
            _draw_polyline(
                draw,
                points,
                color=str(item.get("color", "#1f77b4")),
                width=float(item.get("width", 2.5)),
                opacity=float(item.get("opacity", 0.9)),
                dasharray=item.get("dasharray") if isinstance(item.get("dasharray"), str) else None,
            )

    legend_x0 = width - pad_right - 330
    legend_y = pad_top + 16
    for item in series:
        color = str(item.get("color", "#1f77b4"))
        width_value = float(item.get("width", 2.5))
        dasharray = item.get("dasharray") if isinstance(item.get("dasharray"), str) else None
        mode = str(item.get("mode", "line"))
        if mode == "scatter":
            _draw_scatter(draw, [(legend_x0 + 24, legend_y)], color=color, radius=4.0)
        else:
            _draw_polyline(
                draw,
                [(legend_x0, legend_y), (legend_x0 + 48, legend_y)],
                color=color,
                width=width_value,
                dasharray=dasharray,
                opacity=1.0,
            )
        _draw_text(draw, legend_x0 + 60, legend_y + 1, str(item["label"]), font_legend, _color_with_alpha("#000000"), anchor="lm")
        legend_y += 26

    _save_png(image, output_path)


def save_dual_axis_plot(
    x: array | list[float],
    left_series: list[dict[str, object]],
    right_series: list[dict[str, object]],
    xlabel: str,
    left_ylabel: str,
    right_ylabel: str,
    title: str,
    output_path: Path,
    subtitle: str | None = None,
) -> None:
    width = 1650
    height = 980
    pad_left = 150
    pad_right = 150
    pad_top = 95
    pad_bottom = 120

    image = Image.new("RGBA", (width, height), (255, 255, 255, 255))
    draw = ImageDraw.Draw(image, "RGBA")
    font_title = _load_font(30)
    font_subtitle = _load_font(16)
    font_axis = _load_font(22)
    font_tick = _load_font(16)
    font_legend = _load_font(16)

    x_values = [float(value) for value in x]
    left_values: list[float] = []
    right_values: list[float] = []
    for item in left_series:
        left_values.extend(float(value) for value in item["y"])  # type: ignore[index]
    for item in right_series:
        right_values.extend(float(value) for value in item["y"])  # type: ignore[index]

    x_min = min(x_values)
    x_max = max(x_values)
    left_min = min(left_values)
    left_max = max(left_values)
    right_min = min(right_values)
    right_max = max(right_values)
    if x_max == x_min:
        x_max = x_min + 1.0
    if left_max == left_min:
        left_max = left_min + 1.0
    if right_max == right_min:
        right_max = right_min + 1.0

    left_pad = 0.06 * (left_max - left_min) or 1.0
    right_pad = 0.06 * (right_max - right_min) or 1.0
    left_min -= left_pad
    left_max += left_pad
    right_min -= right_pad
    right_max += right_pad

    plot_width = width - pad_left - pad_right
    plot_height = height - pad_top - pad_bottom

    def sx(value: float) -> float:
        return pad_left + (value - x_min) / (x_max - x_min) * plot_width

    def sy_left(value: float) -> float:
        return height - pad_bottom - (value - left_min) / (left_max - left_min) * plot_height

    def sy_right(value: float) -> float:
        return height - pad_bottom - (value - right_min) / (right_max - right_min) * plot_height

    _draw_text(draw, width / 2.0, 34, title, font_title, _color_with_alpha("#000000"), anchor="mt")
    if subtitle:
        _draw_text(draw, width / 2.0, 63, subtitle, font_subtitle, _color_with_alpha("#555555"), anchor="mt")

    _draw_line_segment(draw, (pad_left, height - pad_bottom), (width - pad_right, height - pad_bottom), _color_with_alpha("#000000"), 2)
    _draw_line_segment(draw, (pad_left, pad_top), (pad_left, height - pad_bottom), _color_with_alpha("#000000"), 2)
    _draw_line_segment(draw, (width - pad_right, pad_top), (width - pad_right, height - pad_bottom), _color_with_alpha("#000000"), 2)
    _draw_text(draw, width / 2.0, height - 28, xlabel, font_axis, _color_with_alpha("#000000"), anchor="mt")
    _draw_text(draw, 34, height / 2.0, left_ylabel, font_axis, _color_with_alpha("#000000"), anchor="mm", angle=90.0)
    _draw_text(draw, width - 34, height / 2.0, right_ylabel, font_axis, _color_with_alpha("#000000"), anchor="mm", angle=-90.0)

    for tick in generate_ticks(x_min, x_max):
        tx = sx(tick)
        _draw_line_segment(draw, (tx, pad_top), (tx, height - pad_bottom), _color_with_alpha("#d9d9d9"), 1)
        _draw_text(draw, tx, height - pad_bottom + 30, f"{tick:.4g}", font_tick, _color_with_alpha("#000000"), anchor="mt")
    for tick in generate_ticks(left_min, left_max):
        ty = sy_left(tick)
        _draw_line_segment(draw, (pad_left, ty), (width - pad_right, ty), _color_with_alpha("#ececec"), 1)
        _draw_text(draw, pad_left - 14, ty, f"{tick:.4g}", font_tick, _color_with_alpha("#1f77b4"), anchor="rm")
    for tick in generate_ticks(right_min, right_max):
        ty = sy_right(tick)
        _draw_text(draw, width - pad_right + 14, ty, f"{tick:.4g}", font_tick, _color_with_alpha("#d62728"), anchor="lm")

    for item in left_series:
        points = [(sx(xv), sy_left(float(yv))) for xv, yv in zip(x_values, item["y"])]  # type: ignore[index]
        _draw_polyline(draw, points, color=str(item.get("color", "#1f77b4")), width=float(item.get("width", 2.5)), opacity=1.0)
    for item in right_series:
        points = [(sx(xv), sy_right(float(yv))) for xv, yv in zip(x_values, item["y"])]  # type: ignore[index]
        _draw_polyline(draw, points, color=str(item.get("color", "#d62728")), width=float(item.get("width", 2.5)), opacity=1.0)

    legend_x0 = width - pad_right - 360
    legend_y = pad_top + 16
    for item in left_series + right_series:
        _draw_polyline(
            draw,
            [(legend_x0, legend_y), (legend_x0 + 48, legend_y)],
            color=str(item.get("color", "#1f77b4")),
            width=float(item.get("width", 2.5)),
            opacity=1.0,
        )
        _draw_text(draw, legend_x0 + 60, legend_y + 1, str(item["label"]), font_legend, _color_with_alpha("#000000"), anchor="lm")
        legend_y += 26

    _save_png(image, output_path)


def _append_panel_series_plot(
    image: Image.Image,
    draw: ImageDraw.ImageDraw,
    x: array | list[float],
    series: list[dict[str, object]],
    xlabel: str,
    ylabel: str,
    title: str,
    box_left: float,
    box_top: float,
    box_width: float,
    box_height: float,
) -> None:
    pad_left = 70
    pad_right = 20
    pad_top = 40
    pad_bottom = 55
    font_title = _load_font(18)
    font_axis = _load_font(14)
    font_tick = _load_font(11)
    font_legend = _load_font(11)

    x_values = [float(value) for value in x]
    y_values: list[float] = []
    for item in series:
        y_values.extend(float(value) for value in item["y"])  # type: ignore[index]

    x_min = min(x_values)
    x_max = max(x_values)
    y_min = min(y_values)
    y_max = max(y_values)
    if x_max == x_min:
        x_max = x_min + 1.0
    if y_max == y_min:
        y_max = y_min + 1.0

    y_pad = 0.06 * (y_max - y_min) or 1.0
    y_min -= y_pad
    y_max += y_pad

    plot_left = box_left + pad_left
    plot_right = box_left + box_width - pad_right
    plot_top = box_top + pad_top
    plot_bottom = box_top + box_height - pad_bottom
    plot_width = plot_right - plot_left
    plot_height = plot_bottom - plot_top

    def sx(value: float) -> float:
        return plot_left + (value - x_min) / (x_max - x_min) * plot_width

    def sy(value: float) -> float:
        return plot_bottom - (value - y_min) / (y_max - y_min) * plot_height

    draw.rectangle((box_left, box_top, box_left + box_width, box_top + box_height), outline=_color_with_alpha("#cccccc"), fill=_color_with_alpha("#ffffff"))
    _draw_text(draw, box_left + box_width / 2.0, box_top + 22, title, font_title, _color_with_alpha("#000000"), anchor="mt")
    _draw_line_segment(draw, (plot_left, plot_bottom), (plot_right, plot_bottom), _color_with_alpha("#000000"), 1)
    _draw_line_segment(draw, (plot_left, plot_top), (plot_left, plot_bottom), _color_with_alpha("#000000"), 1)
    _draw_text(draw, box_left + box_width / 2.0, box_top + box_height - 12, xlabel, font_axis, _color_with_alpha("#000000"), anchor="mt")
    _draw_text(draw, box_left + 18, box_top + box_height / 2.0, ylabel, font_axis, _color_with_alpha("#000000"), anchor="mm", angle=90.0)

    for tick in generate_ticks(x_min, x_max):
        tx = sx(tick)
        _draw_line_segment(draw, (tx, plot_top), (tx, plot_bottom), _color_with_alpha("#e6e6e6"), 1)
        _draw_text(draw, tx, plot_bottom + 18, f"{tick:.4g}", font_tick, _color_with_alpha("#000000"), anchor="mt")
    for tick in generate_ticks(y_min, y_max):
        ty = sy(tick)
        _draw_line_segment(draw, (plot_left, ty), (plot_right, ty), _color_with_alpha("#e6e6e6"), 1)
        _draw_text(draw, plot_left - 8, ty, f"{tick:.4g}", font_tick, _color_with_alpha("#000000"), anchor="rm")

    legend_x = plot_right - 165
    legend_y = plot_top + 16
    for item in series:
        x_series = [float(value) for value in item.get("x", x)]  # type: ignore[arg-type]
        y_series = [float(value) for value in item["y"]]  # type: ignore[index]
        points = [(sx(xv), sy(yv)) for xv, yv in zip(x_series, y_series)]
        _draw_polyline(
            draw,
            points,
            color=str(item.get("color", "#1f77b4")),
            width=float(item.get("width", 2.0)),
            opacity=float(item.get("opacity", 0.9)),
            dasharray=item.get("dasharray") if isinstance(item.get("dasharray"), str) else None,
        )
        _draw_polyline(
            draw,
            [(legend_x, legend_y), (legend_x + 28, legend_y)],
            color=str(item.get("color", "#1f77b4")),
            width=float(item.get("width", 2.0)),
            opacity=1.0,
            dasharray=item.get("dasharray") if isinstance(item.get("dasharray"), str) else None,
        )
        _draw_text(draw, legend_x + 36, legend_y + 1, str(item["label"]), font_legend, _color_with_alpha("#000000"), anchor="lm")
        legend_y += 16


def save_length_kinematics_canvas(
    output_path: Path,
    canvas_title: str,
    derivative_time_series: list[dict[str, object]],
    derivative_length_series: list[dict[str, object]],
    second_time_series: list[dict[str, object]],
    second_length_series: list[dict[str, object]],
    time: array | list[float],
    length: array | list[float],
) -> None:
    width = 1800
    height = 1220
    outer_left = 40
    outer_top = 70
    outer_right = 40
    outer_bottom = 40
    gap_x = 28
    gap_y = 28

    image = Image.new("RGBA", (width, height), (255, 255, 255, 255))
    draw = ImageDraw.Draw(image, "RGBA")
    font_title = _load_font(30)
    font_subtitle = _load_font(16)

    panel_width = (width - outer_left - outer_right - gap_x) / 2.0
    panel_height = (height - outer_top - outer_bottom - gap_y - 40) / 2.0

    _draw_text(draw, width / 2.0, 34, canvas_title, font_title, _color_with_alpha("#000000"), anchor="mt")
    _draw_text(draw, width / 2.0, 58, "Top row: first derivative. Bottom row: second derivative.", font_subtitle, _color_with_alpha("#555555"), anchor="mt")

    _append_panel_series_plot(
        image,
        draw,
        time,
        derivative_time_series,
        xlabel="Time",
        ylabel="dL/dt",
        title="Length Derivative vs Time",
        box_left=outer_left,
        box_top=outer_top,
        box_width=panel_width,
        box_height=panel_height,
    )
    _append_panel_series_plot(
        image,
        draw,
        length,
        derivative_length_series,
        xlabel="Projected length L",
        ylabel="dL/dt",
        title="Length Derivative vs Length",
        box_left=outer_left + panel_width + gap_x,
        box_top=outer_top,
        box_width=panel_width,
        box_height=panel_height,
    )
    _append_panel_series_plot(
        image,
        draw,
        time,
        second_time_series,
        xlabel="Time",
        ylabel="d²L/dt²",
        title="Length Second Derivative vs Time",
        box_left=outer_left,
        box_top=outer_top + panel_height + gap_y,
        box_width=panel_width,
        box_height=panel_height,
    )
    _append_panel_series_plot(
        image,
        draw,
        length,
        second_length_series,
        xlabel="Projected length L",
        ylabel="d²L/dt²",
        title="Length Second Derivative vs Length",
        box_left=outer_left + panel_width + gap_x,
        box_top=outer_top + panel_height + gap_y,
        box_width=panel_width,
        box_height=panel_height,
    )
    _save_png(image, output_path)


def smooth_series(
    x: array | list[float],
    y: array | list[float],
    trend_method: str,
    trend_frac: float,
    trend_max_points: int,
    sg_window: int,
    sg_order: int,
) -> array:
    if trend_method == "savgol":
        return savitzky_golay_series(
            x,
            y,
            derivative_order=0,
            window_size=sg_window,
            poly_order=sg_order,
        )
    return lowess_smooth(x, y, frac=trend_frac, max_points=trend_max_points)


def assemble_dashboard_png(
    output_path: Path,
    top_plot_path: Path,
    paired_plot_paths: list[tuple[Path, Path]],
    title: str,
) -> None:
    top_plot = Image.open(top_plot_path).convert("RGBA")
    pair_images = [
        (Image.open(left_path).convert("RGBA"), Image.open(right_path).convert("RGBA"))
        for left_path, right_path in paired_plot_paths
    ]

    left_width = max(left_image.width for left_image, _ in pair_images)
    right_width = max(right_image.width for _, right_image in pair_images)
    panel_height = max(max(left_image.height, right_image.height) for left_image, right_image in pair_images)
    total_width = left_width + right_width + 48
    header_height = 70
    total_height = header_height + top_plot.height + 28 + len(pair_images) * (panel_height + 28) + 28

    dashboard = Image.new("RGBA", (total_width + 80, total_height), (255, 255, 255, 255))
    draw = ImageDraw.Draw(dashboard, "RGBA")
    font_title = _load_font(30)

    _draw_text(draw, dashboard.width / 2.0, 18, title, font_title, _color_with_alpha("#000000"), anchor="mt")

    y_cursor = header_height
    top_resized = top_plot.resize((total_width, int(round(top_plot.height * total_width / top_plot.width))), Image.Resampling.LANCZOS)
    dashboard.alpha_composite(top_resized, (40, y_cursor))
    y_cursor += top_resized.height + 28

    for left_image, right_image in pair_images:
        left_resized = left_image.resize((left_width, int(round(left_image.height * left_width / left_image.width))), Image.Resampling.LANCZOS)
        right_resized = right_image.resize((right_width, int(round(right_image.height * right_width / right_image.width))), Image.Resampling.LANCZOS)
        left_y = y_cursor + max(0, (panel_height - left_resized.height) // 2)
        right_y = y_cursor + max(0, (panel_height - right_resized.height) // 2)
        dashboard.alpha_composite(left_resized, (40, left_y))
        dashboard.alpha_composite(right_resized, (40 + left_width + 48, right_y))
        y_cursor += panel_height + 28

    _save_png(dashboard, output_path)


def maybe_show_with_matplotlib(
    table: dict[str, array],
    length_derivative_time: array,
    length_derivative_trend_time: array,
    length_derivative_trend_length: tuple[array, array],
    length_second_derivative_time: array,
    length_second_derivative_trend_time: array,
    length_second_derivative_trend_length: tuple[array, array],
    rg_time: array,
    rg_length: array,
    rg_values: array,
    rg_trend_time: array,
    rg_trend_length: tuple[array, array],
    rg_derivative_time: array,
    rg_derivative_trend_time: array,
    rg_derivative_trend_length: tuple[array, array],
    trend_label: str,
    length_second_trend_label: str,
    trend_color: str,
) -> bool:
    plt = maybe_import_matplotlib()
    if plt is None:
        return False

    row_count = 5
    fig = plt.figure(figsize=(18, 4.6 * row_count))
    grid = fig.add_gridspec(row_count, 2)

    ax = fig.add_subplot(grid[0, :])
    ax_right = ax.twinx()
    ax.plot(table["time"], table["L"], color="#1f77b4", linewidth=1.8, label="L")
    ax_right.plot(table["time"], table["contour_fraction"], color="#d62728", linewidth=1.8, label="L/Lcontour")
    ax.set_title("Length and Contour Fraction vs Time")
    ax.set_xlabel("Time")
    ax.set_ylabel("Projected length L")
    ax_right.set_ylabel("Contour fraction")
    ax.grid(alpha=0.25)
    handles_left, labels_left = ax.get_legend_handles_labels()
    handles_right, labels_right = ax_right.get_legend_handles_labels()
    ax.legend(handles_left + handles_right, labels_left + labels_right, loc="best")

    ax = fig.add_subplot(grid[1, 0])
    ax.plot(table["time"], length_derivative_time, color="#1f77b4", linewidth=1.8, label="dL/dt")
    ax.plot(table["time"], length_derivative_trend_time, color=trend_color, linewidth=2.8, label=trend_label)
    ax.axhline(0.0, color="#666", linewidth=1.2, linestyle=":")
    ax.set_title("Length Derivative vs Time")
    ax.set_xlabel("Time")
    ax.set_ylabel("dL/dt")
    ax.grid(alpha=0.25)
    ax.legend(loc="best")

    ax = fig.add_subplot(grid[1, 1])
    ax.plot(table["L"], length_derivative_time, color="#1f77b4", linewidth=1.8, label="dL/dt")
    ax.plot(
        length_derivative_trend_length[0],
        length_derivative_trend_length[1],
        color=trend_color,
        linewidth=2.8,
        label=trend_label,
    )
    ax.axhline(0.0, color="#666", linewidth=1.2, linestyle=":")
    ax.set_title("Length Derivative vs Length")
    ax.set_xlabel("Projected length L")
    ax.set_ylabel("dL/dt")
    ax.grid(alpha=0.25)
    ax.legend(loc="best")

    ax = fig.add_subplot(grid[2, 0])
    ax.plot(table["time"], length_second_derivative_time, color="#1f77b4", linewidth=1.8, label="d²L/dt²")
    ax.plot(table["time"], length_second_derivative_trend_time, color=trend_color, linewidth=2.8, label=length_second_trend_label)
    ax.axhline(0.0, color="#666", linewidth=1.2, linestyle=":")
    ax.set_title("Length Second Derivative vs Time")
    ax.set_xlabel("Time")
    ax.set_ylabel("d²L/dt²")
    ax.grid(alpha=0.25)
    ax.legend(loc="best")

    ax = fig.add_subplot(grid[2, 1])
    ax.plot(table["L"], length_second_derivative_time, color="#1f77b4", linewidth=1.8, label="d²L/dt²")
    ax.plot(
        length_second_derivative_trend_length[0],
        length_second_derivative_trend_length[1],
        color=trend_color,
        linewidth=2.8,
        label=length_second_trend_label,
    )
    ax.axhline(0.0, color="#666", linewidth=1.2, linestyle=":")
    ax.set_title("Length Second Derivative vs Length")
    ax.set_xlabel("Projected length L")
    ax.set_ylabel("d²L/dt²")
    ax.grid(alpha=0.25)
    ax.legend(loc="best")

    ax = fig.add_subplot(grid[3, 0])
    ax.plot(rg_time, rg_values, color="#1f77b4", linewidth=1.8, alpha=0.9, label="Rg")
    ax.plot(rg_time, rg_trend_time, color=trend_color, linewidth=2.8, label=trend_label)
    ax.set_title("Radius of Gyration vs Time")
    ax.set_xlabel("Time")
    ax.set_ylabel("Radius of gyration Rg")
    ax.grid(alpha=0.25)
    ax.legend(loc="best")

    ax = fig.add_subplot(grid[3, 1])
    ax.plot(rg_length, rg_values, color="#1f77b4", linewidth=1.8, alpha=0.9, label="Rg")
    ax.plot(
        rg_trend_length[0],
        rg_trend_length[1],
        color=trend_color,
        linewidth=2.8,
        label=trend_label,
    )
    ax.set_title("Radius of Gyration vs Length")
    ax.set_xlabel("Projected length L")
    ax.set_ylabel("Radius of gyration Rg")
    ax.grid(alpha=0.25)
    ax.legend(loc="best")

    ax = fig.add_subplot(grid[4, 0])
    ax.plot(rg_time, rg_derivative_time, color="#1f77b4", linewidth=1.8, alpha=0.9, label="dRg/dt")
    ax.plot(rg_time, rg_derivative_trend_time, color=trend_color, linewidth=2.8, label=trend_label)
    ax.axhline(0.0, color="#666", linewidth=1.2, linestyle=":")
    ax.set_title("Rg Derivative vs Time")
    ax.set_xlabel("Time")
    ax.set_ylabel("dRg/dt")
    ax.grid(alpha=0.25)
    ax.legend(loc="best")

    ax = fig.add_subplot(grid[4, 1])
    ax.plot(rg_length, rg_derivative_time, color="#1f77b4", linewidth=1.8, alpha=0.9, label="dRg/dt")
    ax.plot(
        rg_derivative_trend_length[0],
        rg_derivative_trend_length[1],
        color=trend_color,
        linewidth=2.8,
        label=trend_label,
    )
    ax.axhline(0.0, color="#666", linewidth=1.2, linestyle=":")
    ax.set_title("Rg Derivative vs Length")
    ax.set_xlabel("Projected length L")
    ax.set_ylabel("dRg/dt")
    ax.grid(alpha=0.25)
    ax.legend(loc="best")

    fig.tight_layout()
    plt.show()
    return True


def resolve_data_file(data_file: Path | None) -> Path:
    if data_file is not None:
        resolved = data_file.resolve()
        if not resolved.exists():
            raise SystemExit(f"Data file not found: {resolved}")
        return resolved

    candidate = (Path.cwd() / DEFAULT_DATA_FILENAME).resolve()
    if candidate.exists():
        return candidate

    raise SystemExit(
        "No data_file was provided and "
        f"{DEFAULT_DATA_FILENAME!r} was not found in the current working directory: {Path.cwd()}"
    )


def parse_args(default_trend_method: str = "loess", default_sg_window: int = 31) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Analyze aligned-box single-chain constant-force output from LAMMPS."
    )
    parser.add_argument(
        "data_file",
        nargs="?",
        type=Path,
        default=None,
        help="Path to force_clamp_response.dat (default: ./force_clamp_response.dat)",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Directory for plots and summary.json (default: sibling analysis/ directory)",
    )
    parser.add_argument(
        "--show",
        action="store_true",
        help="Open an interactive matplotlib window if matplotlib is available in the current Python environment",
    )
    parser.add_argument(
        "--trend-frac",
        type=float,
        default=0.08,
        help="LOWESS neighborhood fraction for orange trend lines",
    )
    parser.add_argument(
        "--trend-max-points",
        type=int,
        default=2000,
        help="Maximum number of points used internally by LOWESS before downsampling+interpolation",
    )
    parser.add_argument(
        "--sg-window",
        type=int,
        default=default_sg_window,
        help="Window size for the Savitzky-Golay-style local polynomial derivative estimate",
    )
    parser.add_argument(
        "--sg-order",
        type=int,
        default=3,
        help="Polynomial order for the Savitzky-Golay-style local polynomial derivative estimate",
    )
    parser.add_argument(
        "--trend-method",
        choices=("loess", "savgol"),
        default=default_trend_method,
        help="Which smoothing method to use for trend lines and the aggregate dashboard",
    )
    parser.add_argument(
        "--jobs",
        type=int,
        default=None,
        help="Parallel worker count for parsing traj.force_clamp_aligned.*.dump files (default: auto)",
    )
    return parser.parse_args()


def main(
    default_trend_method: str = "loess",
    default_output_dir_name: str = "analysis_rg_T",
    default_sg_window: int = 31,
) -> None:
    args = parse_args(default_trend_method=default_trend_method, default_sg_window=default_sg_window)
    data_file = resolve_data_file(args.data_file)
    dataset_path, dataset_bundle, dataset_status = prepare_stretch_dataset(data_file, max_workers=args.jobs)
    if dataset_status == "wait-listed":
        print(
            "Skipped analysis because the dataset is parked in "
            f"{wait_list_path_for_output_dir(data_file.parent)}"
        )
        return

    if dataset_path is None or dataset_bundle is None:
        raise SystemExit("Stretch dataset preparation did not return analysis data.")

    output_dir = args.output_dir or data_file.parent / default_output_dir_name
    output_dir.mkdir(parents=True, exist_ok=True)
    obsolete_plot_names = (
        "length_contour_vs_time",
        "length_derivative_vs_time",
        "length_derivative_vs_length",
        "length_second_derivative_vs_time",
        "length_second_derivative_vs_length",
        "length_kinematics_loess_canvas",
        "length_kinematics_savgol_canvas",
        "mechanics_dashboard_loess",
        "mechanics_dashboard_savgol",
        "resisting_force_vs_time",
        "net_drag_force_vs_time",
        "resisting_force_vs_length",
        "net_drag_force_vs_length",
        "tail_velocity_vs_time",
        "tail_velocity_vs_length",
        "tail_acceleration_vs_time",
        "tail_acceleration_vs_length",
        "tail_total_force_vs_time",
        "tail_total_force_vs_length",
        "tail_chain_force_vs_time",
        "tail_chain_force_vs_length",
        "rg_vs_time",
        "rg_vs_length",
        "rg_derivative_vs_time",
        "rg_derivative_vs_length",
    )
    cleanup_names = {f"{PLOT_FILE_PREFIX}{name}" for name in obsolete_plot_names}
    for obsolete_name in cleanup_names:
        for extension in ("png", "svg"):
            obsolete_path = output_dir / f"{obsolete_name}.{extension}"
            if obsolete_path.exists():
                obsolete_path.unlink()

    table = dataset_bundle["table"]  # type: ignore[assignment]
    derived = derive_mechanics(table)
    length_derivative_time, length_second_derivative_time = compute_length_kinematics(table)
    rg_kinematics = dataset_bundle["rg_kinematics"]  # type: ignore[assignment]
    trend_method = args.trend_method
    if trend_method == "savgol":
        trend_label = "Savitzky-Golay"
        length_second_trend_label = "Savitzky-Golay"
        trend_color = "#ff6a3d"
        length_derivative_trend_time = savitzky_golay_series(
            table["time"],
            table["L"],
            derivative_order=1,
            window_size=args.sg_window,
            poly_order=args.sg_order,
        )
        length_second_derivative_trend_time = savitzky_golay_series(
            table["time"],
            table["L"],
            derivative_order=2,
            window_size=args.sg_window,
            poly_order=args.sg_order,
        )
    else:
        trend_label = "LOESS mean"
        length_second_trend_label = "d/dt(LOESS dL/dt)"
        trend_color = "#ff7f0e"
        length_derivative_trend_time = lowess_smooth(
            table["time"], length_derivative_time, frac=args.trend_frac, max_points=args.trend_max_points
        )
        length_second_derivative_trend_time = finite_difference(table["time"], length_derivative_trend_time)

    length_derivative_trend_length = sorted_pairs(
        table["L"],
        smooth_series(
            table["L"],
            length_derivative_time,
            trend_method=trend_method,
            trend_frac=args.trend_frac,
            trend_max_points=args.trend_max_points,
            sg_window=args.sg_window,
            sg_order=args.sg_order,
        ),
    )
    length_second_derivative_trend_length = sorted_pairs(
        table["L"], length_second_derivative_trend_time
    )

    rg_trend_time = smooth_series(
        rg_kinematics["time"],
        rg_kinematics["rg"],
        trend_method=trend_method,
        trend_frac=args.trend_frac,
        trend_max_points=args.trend_max_points,
        sg_window=args.sg_window,
        sg_order=args.sg_order,
    )
    rg_trend_length = sorted_pairs(
        rg_kinematics["length"],
        smooth_series(
            rg_kinematics["length"],
            rg_kinematics["rg"],
            trend_method=trend_method,
            trend_frac=args.trend_frac,
            trend_max_points=args.trend_max_points,
            sg_window=args.sg_window,
            sg_order=args.sg_order,
        ),
    )
    rg_derivative_trend_time = smooth_series(
        rg_kinematics["time"],
        rg_kinematics["rg_derivative_time"],
        trend_method=trend_method,
        trend_frac=args.trend_frac,
        trend_max_points=args.trend_max_points,
        sg_window=args.sg_window,
        sg_order=args.sg_order,
    )
    rg_derivative_trend_length = sorted_pairs(
        rg_kinematics["length"],
        smooth_series(
            rg_kinematics["length"],
            rg_kinematics["rg_derivative_time"],
            trend_method=trend_method,
            trend_frac=args.trend_frac,
            trend_max_points=args.trend_max_points,
            sg_window=args.sg_window,
            sg_order=args.sg_order,
        ),
    )

    length_contour_path = output_dir / f"{PLOT_FILE_PREFIX}length_contour_vs_time.png"
    length_derivative_time_path = output_dir / f"{PLOT_FILE_PREFIX}length_derivative_vs_time.png"
    length_derivative_length_path = output_dir / f"{PLOT_FILE_PREFIX}length_derivative_vs_length.png"
    length_second_time_path = output_dir / f"{PLOT_FILE_PREFIX}length_second_derivative_vs_time.png"
    length_second_length_path = output_dir / f"{PLOT_FILE_PREFIX}length_second_derivative_vs_length.png"
    rg_time_path = output_dir / f"{PLOT_FILE_PREFIX}rg_vs_time.png"
    rg_length_path = output_dir / f"{PLOT_FILE_PREFIX}rg_vs_length.png"
    rg_derivative_time_path = output_dir / f"{PLOT_FILE_PREFIX}rg_derivative_vs_time.png"
    rg_derivative_length_path = output_dir / f"{PLOT_FILE_PREFIX}rg_derivative_vs_length.png"

    save_dual_axis_plot(
        table["time"],
        left_series=[
            {"y": table["L"], "label": "L", "color": "#1f77b4", "width": 2.2},
        ],
        right_series=[
            {"y": table["contour_fraction"], "label": "L/Lcontour", "color": "#d62728", "width": 2.2},
        ],
        xlabel="Time",
        left_ylabel="Projected length L",
        right_ylabel="Contour fraction",
        title="Length and Contour Fraction vs Time",
        output_path=length_contour_path,
        subtitle="The only default time-axis plot. The remaining panels compare length and radius-of-gyration responses against time or projected length.",
    )
    save_series_plot(
        table["time"],
        [
            {"y": length_derivative_time, "label": "dL/dt", "color": "#1f77b4", "width": 2.0},
            {"y": length_derivative_trend_time, "label": trend_label, "color": trend_color, "width": 4.0},
            {"y": array("d", [0.0] * len(table["time"])), "label": "zero", "color": "#666666", "width": 1.3, "dasharray": "4 6"},
        ],
        xlabel="Time",
        ylabel="Length derivative dL/dt",
        title="Length Derivative vs Time",
        output_path=length_derivative_time_path,
        subtitle="Direct finite-difference derivative of the sampled L(t) series. Positive means the projected end-to-end length is increasing.",
    )
    save_series_plot(
        table["L"],
        [
            {"y": length_derivative_time, "label": "dL/dt", "color": "#1f77b4", "width": 1.8, "opacity": 0.85},
            {"x": length_derivative_trend_length[0], "y": length_derivative_trend_length[1], "label": trend_label, "color": trend_color, "width": 4.0},
            {"y": array("d", [0.0] * len(table["L"])), "label": "zero", "color": "#666666", "width": 1.3, "dasharray": "4 6"},
        ],
        xlabel="Projected length L",
        ylabel="Length derivative dL/dt",
        title="Length Derivative vs Length",
        output_path=length_derivative_length_path,
        subtitle="Time-parametrized curve: x is current projected length, y is the simultaneous dL/dt.",
    )
    save_series_plot(
        table["time"],
        [
            {"y": length_second_derivative_time, "label": "d²L/dt²", "color": "#1f77b4", "width": 2.0},
            {"y": length_second_derivative_trend_time, "label": length_second_trend_label, "color": trend_color, "width": 4.0},
            {"y": array("d", [0.0] * len(table["time"])), "label": "zero", "color": "#666666", "width": 1.3, "dasharray": "4 6"},
        ],
        xlabel="Time",
        ylabel="Length second derivative d²L/dt²",
        title="Length Second Derivative vs Time",
        output_path=length_second_time_path,
        subtitle="Second finite-difference derivative of the sampled L(t) series. Positive means acceleration of lengthening; negative means deceleration.",
    )
    save_series_plot(
        table["L"],
        [
            {"y": length_second_derivative_time, "label": "d²L/dt²", "color": "#1f77b4", "width": 1.8, "opacity": 0.85},
            {"x": length_second_derivative_trend_length[0], "y": length_second_derivative_trend_length[1], "label": length_second_trend_label, "color": trend_color, "width": 4.0},
            {"y": array("d", [0.0] * len(table["L"])), "label": "zero", "color": "#666666", "width": 1.3, "dasharray": "4 6"},
        ],
        xlabel="Projected length L",
        ylabel="Length second derivative d²L/dt²",
        title="Length Second Derivative vs Length",
        output_path=length_second_length_path,
        subtitle="Time-parametrized curve: x is current projected length, y is the simultaneous d²L/dt².",
    )
    save_series_plot(
        rg_kinematics["time"],
        [
            {"y": rg_kinematics["rg"], "label": "Rg", "color": "#1f77b4", "width": 1.8, "opacity": 0.9},
            {"y": rg_trend_time, "label": trend_label, "color": trend_color, "width": 4.0},
        ],
        xlabel="Time",
        ylabel="Radius of gyration Rg",
        title="Radius of Gyration vs Time",
        output_path=rg_time_path,
        subtitle="Mass-weighted Rg computed from type-1 and type-2 particles in traj.force_clamp_aligned.*.dump using unwrapped coordinates when available.",
    )
    save_series_plot(
        rg_kinematics["length"],
        [
            {"y": rg_kinematics["rg"], "label": "Rg", "color": "#1f77b4", "width": 1.8, "opacity": 0.9},
            {"x": rg_trend_length[0], "y": rg_trend_length[1], "label": trend_label, "color": trend_color, "width": 4.0},
        ],
        xlabel="Projected length L",
        ylabel="Radius of gyration Rg",
        title="Radius of Gyration vs Length",
        output_path=rg_length_path,
        subtitle="Time-parametrized curve: x is projected chain length from force_clamp_response.dat, y is the simultaneous mass-weighted Rg from the dump trajectory.",
    )
    save_series_plot(
        rg_kinematics["time"],
        [
            {"y": rg_kinematics["rg_derivative_time"], "label": "dRg/dt", "color": "#1f77b4", "width": 1.8, "opacity": 0.9},
            {"y": rg_derivative_trend_time, "label": trend_label, "color": trend_color, "width": 4.0},
            {"y": array("d", [0.0] * len(rg_kinematics["time"])), "label": "zero", "color": "#666666", "width": 1.3, "dasharray": "4 6"},
        ],
        xlabel="Time",
        ylabel="Rg derivative dRg/dt",
        title="Rg Derivative vs Time",
        output_path=rg_derivative_time_path,
        subtitle="Finite-difference derivative of the mass-weighted Rg(t) series measured from the dump trajectory.",
    )
    save_series_plot(
        rg_kinematics["length"],
        [
            {"y": rg_kinematics["rg_derivative_time"], "label": "dRg/dt", "color": "#1f77b4", "width": 1.8, "opacity": 0.9},
            {"x": rg_derivative_trend_length[0], "y": rg_derivative_trend_length[1], "label": trend_label, "color": trend_color, "width": 4.0},
            {"y": array("d", [0.0] * len(rg_kinematics["length"])), "label": "zero", "color": "#666666", "width": 1.3, "dasharray": "4 6"},
        ],
        xlabel="Projected length L",
        ylabel="Rg derivative dRg/dt",
        title="Rg Derivative vs Length",
        output_path=rg_derivative_length_path,
        subtitle="Time-parametrized curve: x is projected length from the response table, y is the simultaneous dRg/dt from the dump-derived Rg series.",
    )

    summary = build_summary(table)
    summary.update(
        {
            "final_rg": float(rg_kinematics["rg"][-1]),
            "max_rg": float(max(rg_kinematics["rg"])),
            "min_rg": float(min(rg_kinematics["rg"])),
            "max_abs_rg_derivative": float(max(abs(value) for value in rg_kinematics["rg_derivative_time"])),
        }
    )
    with (output_dir / "summary.json").open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)

    if args.show and not maybe_show_with_matplotlib(
        table,
        length_derivative_time=length_derivative_time,
        length_derivative_trend_time=length_derivative_trend_time,
        length_derivative_trend_length=length_derivative_trend_length,
        length_second_derivative_time=length_second_derivative_time,
        length_second_derivative_trend_time=length_second_derivative_trend_time,
        length_second_derivative_trend_length=length_second_derivative_trend_length,
        rg_time=rg_kinematics["time"],
        rg_length=rg_kinematics["length"],
        rg_values=rg_kinematics["rg"],
        rg_trend_time=rg_trend_time,
        rg_trend_length=rg_trend_length,
        rg_derivative_time=rg_kinematics["rg_derivative_time"],
        rg_derivative_trend_time=rg_derivative_trend_time,
        rg_derivative_trend_length=rg_derivative_trend_length,
        trend_label=trend_label,
        length_second_trend_label=length_second_trend_label,
        trend_color=trend_color,
    ):
        print("Matplotlib is not available in the current Python environment; kept PNG output only.")

    if float(derived["max_dL_column_mismatch"]) > 1.0e-6 or float(derived["max_strain_column_mismatch"]) > 1.0e-6:
        print(
            "Warning: input file dL/strain columns are not self-consistent with L and contour_fraction; "
            "analysis recomputed dL and strain from L instead."
        )

    if dataset_status == "created":
        print(f"Wrote stretch dataset to {dataset_path}")
    else:
        print(f"Reused stretch dataset {dataset_path}")
    print(f"Wrote analysis outputs to {output_dir}")


if __name__ == "__main__":
    main()
