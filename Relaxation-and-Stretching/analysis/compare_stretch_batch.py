#!/usr/bin/env python3

from __future__ import annotations

import argparse
from array import array
import math
from pathlib import Path
import re
from tempfile import TemporaryDirectory

from PIL import Image, ImageDraw

from analyze_stretch_single import (
    ANALYSIS_CACHE_DIRNAME,
    DEFAULT_DATA_FILENAME,
    _color_with_alpha,
    _draw_line_segment,
    _draw_polyline,
    _draw_scatter,
    _draw_text,
    _load_font,
    _save_png,
    create_process_pool,
    dataset_rows_to_analysis,
    discover_cached_dataset_paths,
    ensure_analysis_root,
    ensure_stretch_dataset_file,
    generate_ticks,
    load_stretch_dataset_rows,
    maybe_import_matplotlib,
    resolve_max_workers,
    save_dual_axis_plot,
    save_series_plot,
)


OUTPUT_PREFIX = "output_"
PALETTE = [
    "#1f77b4",
    "#ff7f0e",
    "#2ca02c",
    "#d62728",
    "#9467bd",
    "#8c564b",
    "#e377c2",
    "#7f7f7f",
    "#bcbd22",
    "#17becf",
    "#4c78a8",
    "#f58518",
    "#54a24b",
    "#e45756",
    "#72b7b2",
    "#b279a2",
    "#ff9da6",
    "#9d755d",
    "#bab0ab",
    "#6f4c9b",
]
FORCE_PATTERN = re.compile(r"(?:^|_)F(?P<force>-?\d+(?:\.\d+)?)")


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
        description="Build or reuse stretch datasets, then compare cached length, contour, Lx, Rg, and Rg_x trajectories."
    )
    parser.add_argument(
        "--root",
        type=Path,
        default=Path.cwd(),
        help="Root directory containing output_* folders and the sibling central cache (default: current working directory)",
    )
    parser.add_argument(
        "--output-file",
        type=Path,
        default=None,
        help=f"Output PNG path for the time-series canvas (default: <root>/{ANALYSIS_CACHE_DIRNAME}/stretch_batch_compare.png)",
    )
    parser.add_argument(
        "--show",
        action="store_true",
        help="Open the saved comparison figures sequentially if matplotlib is available",
    )
    parser.add_argument(
        "--jobs",
        type=int,
        default=None,
        help="Parallel output-directory count for dataset completion / dataset loading (default: auto)",
    )
    parser.add_argument(
        "--analysis-jobs",
        type=int,
        default=None,
        help="Parallel dump-parser worker count used inside each dataset build (default: auto when jobs=1, else 1)",
    )
    parser.add_argument(
        "--window-points",
        type=int,
        default=3000,
        help="Number of cached tail samples used for mean/error-bar summaries vs F (default: 3000)",
    )
    return parser.parse_args()


def discover_output_data_files(root: Path) -> tuple[list[Path], list[Path]]:
    matched: list[Path] = []
    skipped: list[Path] = []
    for child in sorted(root.iterdir(), key=lambda path: natural_sort_key(path.name)):
        if not child.is_dir() or not child.name.startswith(OUTPUT_PREFIX):
            continue
        data_file = child / DEFAULT_DATA_FILENAME
        if data_file.exists():
            matched.append(data_file)
        else:
            skipped.append(child)
    return matched, skipped


def ensure_dataset_task(task: tuple[str, int | None]) -> dict[str, str]:
    data_file_str, analysis_jobs = task
    dataset_path, status = ensure_stretch_dataset_file(Path(data_file_str), max_workers=analysis_jobs)
    return {
        "data_file": data_file_str,
        "dataset_path": "" if dataset_path is None else str(dataset_path),
        "status": status,
    }


def load_dataset_bundle_task(dataset_path_str: str) -> dict[str, object]:
    dataset_path = Path(dataset_path_str)
    return {
        "dataset_path": str(dataset_path),
        "bundle": dataset_rows_to_analysis(load_stretch_dataset_rows(dataset_path)),
    }


def paste_centered(
    canvas: Image.Image,
    image: Image.Image,
    x0: int,
    y0: int,
    box_width: int,
    box_height: int,
) -> None:
    resized = image.resize(
        (box_width, int(round(image.height * box_width / image.width))),
        Image.Resampling.LANCZOS,
    )
    y_shift = max(0, (box_height - resized.height) // 2)
    canvas.alpha_composite(resized, (x0, y0 + y_shift))


def parse_force_from_source_tag(source_tag: str) -> float:
    match = FORCE_PATTERN.search(source_tag)
    if not match:
        raise ValueError(f"Could not parse force from source tag: {source_tag}")
    return float(match.group("force"))


def compute_mean_std(values: list[float]) -> tuple[float, float]:
    if not values:
        raise ValueError("Cannot summarize an empty value list")
    mean_value = sum(values) / len(values)
    variance = sum((value - mean_value) ** 2 for value in values) / len(values)
    return mean_value, math.sqrt(max(variance, 0.0))


def tail_window(values: array | list[float], window_points: int) -> list[float]:
    numeric = [float(value) for value in values]
    if not numeric:
        raise ValueError("Cannot take a tail window of an empty series")
    if window_points <= 0:
        raise ValueError("window_points must be positive")
    return numeric[-min(window_points, len(numeric)) :]


def build_force_summary_rows(
    payloads: list[dict[str, object]],
    window_points: int,
) -> list[dict[str, float]]:
    grouped: dict[float, dict[str, list[float]]] = {}

    for payload in payloads:
        force = parse_force_from_source_tag(str(payload["source_tag"]))
        table = payload["table"]  # type: ignore[assignment]
        rg_kinematics = payload["rg_kinematics"]  # type: ignore[assignment]
        bucket = grouped.setdefault(
            force,
            {
                "L": [],
                "contour_fraction": [],
                "Lx": [],
                "Rg": [],
                "Rg_x": [],
            },
        )
        bucket["L"].extend(tail_window(table["L"], window_points))
        bucket["contour_fraction"].extend(tail_window(table["contour_fraction"], window_points))
        bucket["Lx"].extend(tail_window(table.get("Lx", table["L"]), window_points))
        bucket["Rg"].extend(tail_window(rg_kinematics["rg"], window_points))
        bucket["Rg_x"].extend(tail_window(rg_kinematics.get("rg_x", rg_kinematics["rg"]), window_points))

    rows: list[dict[str, float]] = []
    for force in sorted(grouped):
        metrics = grouped[force]
        mean_l, std_l = compute_mean_std(metrics["L"])
        mean_contour, std_contour = compute_mean_std(metrics["contour_fraction"])
        mean_lx, std_lx = compute_mean_std(metrics["Lx"])
        mean_rg, std_rg = compute_mean_std(metrics["Rg"])
        mean_rg_x, std_rg_x = compute_mean_std(metrics["Rg_x"])
        rows.append(
            {
                "force": force,
                "mean_L": mean_l,
                "std_L": std_l,
                "mean_contour_fraction": mean_contour,
                "std_contour_fraction": std_contour,
                "mean_Lx": mean_lx,
                "std_Lx": std_lx,
                "mean_Rg": mean_rg,
                "std_Rg": std_rg,
                "mean_Rg_x": mean_rg_x,
                "std_Rg_x": std_rg_x,
            }
        )
    return rows


def _draw_error_bars(
    draw: ImageDraw.ImageDraw,
    x_points: list[float],
    low_points: list[float],
    high_points: list[float],
    color: str,
    width: int = 2,
    cap_half_width: int = 8,
    opacity: float = 0.95,
) -> None:
    rgba = _color_with_alpha(color, opacity)
    for x_value, low_value, high_value in zip(x_points, low_points, high_points):
        y0 = min(low_value, high_value)
        y1 = max(low_value, high_value)
        draw.line((x_value, y0, x_value, y1), fill=rgba, width=width)
        draw.line(
            (x_value - cap_half_width, low_value, x_value + cap_half_width, low_value),
            fill=rgba,
            width=width,
        )
        draw.line(
            (x_value - cap_half_width, high_value, x_value + cap_half_width, high_value),
            fill=rgba,
            width=width,
        )


def save_errorbar_series_plot(
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
        for y_value, yerr_value in zip(item["y"], item["yerr"]):  # type: ignore[index]
            y_values.append(float(y_value) - float(yerr_value))
            y_values.append(float(y_value) + float(yerr_value))

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

    _draw_line_segment(
        draw,
        (pad_left, height - pad_bottom),
        (width - pad_right, height - pad_bottom),
        _color_with_alpha("#000000"),
        2,
    )
    _draw_line_segment(
        draw,
        (pad_left, pad_top),
        (pad_left, height - pad_bottom),
        _color_with_alpha("#000000"),
        2,
    )
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
        x_series = [float(value) for value in item.get("x", x)]  # type: ignore[arg-type]
        y_series = [float(value) for value in item["y"]]  # type: ignore[index]
        yerr_series = [float(value) for value in item["yerr"]]  # type: ignore[index]
        x_points = [sx(value) for value in x_series]
        y_points = [sy(value) for value in y_series]
        low_points = [sy(y_value - yerr_value) for y_value, yerr_value in zip(y_series, yerr_series)]
        high_points = [sy(y_value + yerr_value) for y_value, yerr_value in zip(y_series, yerr_series)]
        points = list(zip(x_points, y_points))
        color = str(item.get("color", "#1f77b4"))
        _draw_error_bars(draw, x_points, low_points, high_points, color=color, width=2)
        _draw_polyline(
            draw,
            points,
            color=color,
            width=float(item.get("width", 2.5)),
            opacity=float(item.get("opacity", 0.95)),
        )
        _draw_scatter(
            draw,
            points,
            color=color,
            radius=3.3,
            opacity=float(item.get("opacity", 0.95)),
        )

    legend_x0 = width - pad_right - 330
    legend_y = pad_top + 16
    for item in series:
        color = str(item.get("color", "#1f77b4"))
        _draw_polyline(draw, [(legend_x0, legend_y), (legend_x0 + 48, legend_y)], color=color, width=3.0, opacity=1.0)
        _draw_scatter(draw, [(legend_x0 + 24, legend_y)], color=color, radius=4.0, opacity=1.0)
        _draw_text(draw, legend_x0 + 60, legend_y + 1, str(item["label"]), font_legend, _color_with_alpha("#000000"), anchor="lm")
        legend_y += 26

    _save_png(image, output_path)


def save_dual_axis_errorbar_plot(
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
        for y_value, yerr_value in zip(item["y"], item["yerr"]):  # type: ignore[index]
            left_values.append(float(y_value) - float(yerr_value))
            left_values.append(float(y_value) + float(yerr_value))
    for item in right_series:
        for y_value, yerr_value in zip(item["y"], item["yerr"]):  # type: ignore[index]
            right_values.append(float(y_value) - float(yerr_value))
            right_values.append(float(y_value) + float(yerr_value))

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
        x_series = [float(value) for value in item.get("x", x)]  # type: ignore[arg-type]
        y_series = [float(value) for value in item["y"]]  # type: ignore[index]
        yerr_series = [float(value) for value in item["yerr"]]  # type: ignore[index]
        x_points = [sx(value) for value in x_series]
        y_points = [sy_left(value) for value in y_series]
        low_points = [sy_left(y_value - yerr_value) for y_value, yerr_value in zip(y_series, yerr_series)]
        high_points = [sy_left(y_value + yerr_value) for y_value, yerr_value in zip(y_series, yerr_series)]
        points = list(zip(x_points, y_points))
        color = str(item.get("color", "#1f77b4"))
        _draw_error_bars(draw, x_points, low_points, high_points, color=color, width=2)
        _draw_polyline(draw, points, color=color, width=float(item.get("width", 2.4)), opacity=0.95)
        _draw_scatter(draw, points, color=color, radius=3.3, opacity=0.95)

    for item in right_series:
        x_series = [float(value) for value in item.get("x", x)]  # type: ignore[arg-type]
        y_series = [float(value) for value in item["y"]]  # type: ignore[index]
        yerr_series = [float(value) for value in item["yerr"]]  # type: ignore[index]
        x_points = [sx(value) for value in x_series]
        y_points = [sy_right(value) for value in y_series]
        low_points = [sy_right(y_value - yerr_value) for y_value, yerr_value in zip(y_series, yerr_series)]
        high_points = [sy_right(y_value + yerr_value) for y_value, yerr_value in zip(y_series, yerr_series)]
        points = list(zip(x_points, y_points))
        color = str(item.get("color", "#d62728"))
        _draw_error_bars(draw, x_points, low_points, high_points, color=color, width=2)
        _draw_polyline(draw, points, color=color, width=float(item.get("width", 2.4)), opacity=0.95)
        _draw_scatter(draw, points, color=color, radius=3.3, opacity=0.95)

    legend_x0 = width - pad_right - 330
    legend_y = pad_top + 16
    for item in left_series + right_series:
        color = str(item.get("color", "#1f77b4"))
        _draw_polyline(draw, [(legend_x0, legend_y), (legend_x0 + 48, legend_y)], color=color, width=3.0, opacity=1.0)
        _draw_scatter(draw, [(legend_x0 + 24, legend_y)], color=color, radius=4.0, opacity=1.0)
        _draw_text(draw, legend_x0 + 60, legend_y + 1, str(item["label"]), font_legend, _color_with_alpha("#000000"), anchor="lm")
        legend_y += 26

    _save_png(image, output_path)


def save_comparison_canvas(
    image_paths: list[Path],
    output_file: Path,
    title: str,
    subtitle: str,
) -> None:
    opened = [Image.open(path).convert("RGBA") for path in image_paths]
    try:
        panel_width = max(image.width for image in opened)
        panel_height = max(image.height for image in opened)
        gap_x = 36
        gap_y = 36
        header_height = 92
        margin = 36

        canvas = Image.new(
            "RGBA",
            (
                margin * 2 + panel_width * 2 + gap_x,
                margin * 2 + header_height + panel_height * 2 + gap_y,
            ),
            (255, 255, 255, 255),
        )
        draw = ImageDraw.Draw(canvas, "RGBA")
        font_title = _load_font(30)
        font_subtitle = _load_font(16)
        _draw_text(draw, canvas.width / 2.0, 24, title, font_title, _color_with_alpha("#000000"), anchor="mt")
        _draw_text(draw, canvas.width / 2.0, 56, subtitle, font_subtitle, _color_with_alpha("#555555"), anchor="mt")

        top_y = margin + header_height
        bottom_y = top_y + panel_height + gap_y
        left_x = margin
        right_x = margin + panel_width + gap_x

        paste_centered(canvas, opened[0], left_x, top_y, panel_width, panel_height)
        paste_centered(canvas, opened[1], right_x, top_y, panel_width, panel_height)
        paste_centered(canvas, opened[2], left_x, bottom_y, panel_width, panel_height)
        paste_centered(canvas, opened[3], right_x, bottom_y, panel_width, panel_height)
        _save_png(canvas, output_file)
    finally:
        for image in opened:
            image.close()


def save_time_canvas(payloads: list[dict[str, object]], output_file: Path) -> None:
    output_file.parent.mkdir(parents=True, exist_ok=True)

    combined_time = array(
        "d",
        (value for payload in payloads for value in payload["table"]["time"]),  # type: ignore[index]
    )
    combined_rg_time = array(
        "d",
        (value for payload in payloads for value in payload["rg_kinematics"]["time"]),  # type: ignore[index]
    )

    length_series: list[dict[str, object]] = []
    contour_series: list[dict[str, object]] = []
    lx_series: list[dict[str, object]] = []
    rg_series: list[dict[str, object]] = []
    rg_x_series: list[dict[str, object]] = []

    for index, payload in enumerate(payloads):
        color = PALETTE[index % len(PALETTE)]
        label = str(payload["source_tag"])
        table = payload["table"]  # type: ignore[assignment]
        rg_kinematics = payload["rg_kinematics"]  # type: ignore[assignment]
        time_values = [float(value) for value in table["time"]]
        rg_time_values = [float(value) for value in rg_kinematics["time"]]

        length_series.append({"x": time_values, "y": table["L"], "label": f"{label} | L", "color": color, "width": 2.0})
        contour_series.append(
            {
                "x": time_values,
                "y": table["contour_fraction"],
                "label": f"{label} | contour",
                "color": color,
                "width": 1.6,
                "dasharray": "6 4",
            }
        )
        lx_series.append({"x": time_values, "y": table["Lx"], "label": label, "color": color, "width": 2.0})
        rg_series.append({"x": rg_time_values, "y": rg_kinematics["rg"], "label": label, "color": color, "width": 2.0})
        rg_x_series.append({"x": rg_time_values, "y": rg_kinematics["rg_x"], "label": label, "color": color, "width": 2.0})

    with TemporaryDirectory(prefix="rg_t_compare_time_") as temp_dir_str:
        temp_dir = Path(temp_dir_str)
        length_path = temp_dir / "length_contour.png"
        lx_path = temp_dir / "lx.png"
        rg_path = temp_dir / "rg.png"
        rg_x_path = temp_dir / "rg_x.png"

        save_dual_axis_plot(
            combined_time,
            left_series=length_series,
            right_series=contour_series,
            xlabel="Time",
            left_ylabel="Projected length L",
            right_ylabel="Contour fraction",
            title="Length and Contour Fraction vs Time",
            output_path=length_path,
            subtitle="Solid lines = L from force_clamp_response.dat. Dashed lines = contour fraction from the same cached stretch dataset.",
        )
        save_series_plot(
            combined_time,
            lx_series,
            xlabel="Time",
            ylabel="X-projected span Lx",
            title="Lx vs Time",
            output_path=lx_path,
            subtitle="Lx is computed from the dump trajectory as max(x) - min(x) using type-1/type-2 particles and unwrapped coordinates when available.",
        )
        save_series_plot(
            combined_rg_time,
            rg_series,
            xlabel="Time",
            ylabel="Radius of gyration Rg",
            title="Rg vs Time",
            output_path=rg_path,
            subtitle="Mass-weighted Rg computed from type-1/type-2 particles in traj.force_clamp_aligned.*.dump.",
        )
        save_series_plot(
            combined_rg_time,
            rg_x_series,
            xlabel="Time",
            ylabel="X-projected radius of gyration Rg_x",
            title="Rg_x vs Time",
            output_path=rg_x_path,
            subtitle="Rg_x is the mass-weighted x-axis projection of Rg from the same dump trajectory.",
        )

        save_comparison_canvas(
            [length_path, lx_path, rg_path, rg_x_path],
            output_file=output_file,
            title="Batch Comparison of Stretch Trajectories",
            subtitle=f"Data source: {ANALYSIS_CACHE_DIRNAME}/stretch_dataset_*.dat | panels = L/contour, Lx, Rg, Rg_x",
        )


def save_force_summary_canvas(
    payloads: list[dict[str, object]],
    output_file: Path,
    window_points: int,
) -> None:
    output_file.parent.mkdir(parents=True, exist_ok=True)
    summary_rows = build_force_summary_rows(payloads, window_points=window_points)
    if not summary_rows:
        raise ValueError("No summary rows were produced for the force comparison canvas")

    force_values = array("d", (float(row["force"]) for row in summary_rows))
    with TemporaryDirectory(prefix="rg_t_compare_force_") as temp_dir_str:
        temp_dir = Path(temp_dir_str)
        length_path = temp_dir / "length_contour_force.png"
        lx_path = temp_dir / "lx_force.png"
        rg_path = temp_dir / "rg_force.png"
        rg_x_path = temp_dir / "rg_x_force.png"

        save_dual_axis_errorbar_plot(
            force_values,
            left_series=[
                {
                    "x": force_values,
                    "y": array("d", (row["mean_L"] for row in summary_rows)),
                    "yerr": array("d", (row["std_L"] for row in summary_rows)),
                    "label": "mean L ± σ",
                    "color": "#1f77b4",
                }
            ],
            right_series=[
                {
                    "x": force_values,
                    "y": array("d", (row["mean_contour_fraction"] for row in summary_rows)),
                    "yerr": array("d", (row["std_contour_fraction"] for row in summary_rows)),
                    "label": "mean contour ± σ",
                    "color": "#d62728",
                }
            ],
            xlabel="Applied force F",
            left_ylabel="Projected length L",
            right_ylabel="Contour fraction",
            title="Tail-window Mean Length and Contour Fraction vs F",
            output_path=length_path,
            subtitle=f"Each point uses the last {window_points} cached samples (or fewer if the series is shorter). Error bars are population standard deviations.",
        )
        save_errorbar_series_plot(
            force_values,
            [
                {
                    "x": force_values,
                    "y": array("d", (row["mean_Lx"] for row in summary_rows)),
                    "yerr": array("d", (row["std_Lx"] for row in summary_rows)),
                    "label": "mean Lx ± σ",
                    "color": "#ff7f0e",
                }
            ],
            xlabel="Applied force F",
            ylabel="X-projected span Lx",
            title="Tail-window Mean Lx vs F",
            output_path=lx_path,
            subtitle=f"Each point uses the last {window_points} cached samples of dump-derived Lx. Error bars are population standard deviations.",
        )
        save_errorbar_series_plot(
            force_values,
            [
                {
                    "x": force_values,
                    "y": array("d", (row["mean_Rg"] for row in summary_rows)),
                    "yerr": array("d", (row["std_Rg"] for row in summary_rows)),
                    "label": "mean Rg ± σ",
                    "color": "#2ca02c",
                }
            ],
            xlabel="Applied force F",
            ylabel="Radius of gyration Rg",
            title="Tail-window Mean Rg vs F",
            output_path=rg_path,
            subtitle=f"Each point uses the last {window_points} cached samples of dump-derived Rg. Error bars are population standard deviations.",
        )
        save_errorbar_series_plot(
            force_values,
            [
                {
                    "x": force_values,
                    "y": array("d", (row["mean_Rg_x"] for row in summary_rows)),
                    "yerr": array("d", (row["std_Rg_x"] for row in summary_rows)),
                    "label": "mean Rg_x ± σ",
                    "color": "#9467bd",
                }
            ],
            xlabel="Applied force F",
            ylabel="X-projected radius of gyration Rg_x",
            title="Tail-window Mean Rg_x vs F",
            output_path=rg_x_path,
            subtitle=f"Each point uses the last {window_points} cached samples of dump-derived Rg_x. Error bars are population standard deviations.",
        )

        save_comparison_canvas(
            [length_path, lx_path, rg_path, rg_x_path],
            output_file=output_file,
            title="Tail-window Force Response Summary",
            subtitle=f"Data source: {ANALYSIS_CACHE_DIRNAME}/stretch_dataset_*.dat | each point = last {window_points} cached samples per dataset",
        )


def show_png_sequence(image_paths: list[Path]) -> None:
    plt = maybe_import_matplotlib()
    if plt is None:
        print("Matplotlib is not available in the current Python environment; kept PNG output only.")
        return

    for image_path in image_paths:
        image = Image.open(image_path).convert("RGBA")
        try:
            figure = plt.figure(figsize=(16, 10))
            axis = figure.add_subplot(111)
            axis.imshow(image)
            axis.axis("off")
            plt.tight_layout()
            plt.show()
        finally:
            image.close()


def collect_dataset_completion_results(
    data_files: list[Path],
    directory_jobs: int,
    analysis_jobs: int | None,
) -> list[dict[str, str]]:
    tasks = [(str(data_file), analysis_jobs) for data_file in data_files]
    if directory_jobs == 1:
        return [ensure_dataset_task(task) for task in tasks]

    chunksize = max(1, len(tasks) // (directory_jobs * 4))
    executor = create_process_pool(directory_jobs)
    if executor is None:
        return [ensure_dataset_task(task) for task in tasks]
    try:
        return list(executor.map(ensure_dataset_task, tasks, chunksize=chunksize))
    finally:
        executor.shutdown(wait=True)


def load_dataset_payloads(dataset_paths: list[Path], directory_jobs: int) -> list[dict[str, object]]:
    tasks = [str(path) for path in dataset_paths]
    if directory_jobs == 1:
        loaded = [load_dataset_bundle_task(task) for task in tasks]
    else:
        chunksize = max(1, len(tasks) // (directory_jobs * 4))
        executor = create_process_pool(directory_jobs)
        if executor is None:
            loaded = [load_dataset_bundle_task(task) for task in tasks]
        else:
            try:
                loaded = list(executor.map(load_dataset_bundle_task, tasks, chunksize=chunksize))
            finally:
                executor.shutdown(wait=True)

    payloads = [item["bundle"] for item in loaded]  # type: ignore[index]
    payloads.sort(key=lambda payload: natural_sort_key(str(payload["source_tag"])))
    return payloads  # type: ignore[return-value]


def main() -> int:
    args = parse_args()
    root = args.root.expanduser().resolve()
    if not root.is_dir():
        raise SystemExit(f"Root directory does not exist: {root}")

    analysis_root = ensure_analysis_root(root)
    output_file = (args.output_file or (analysis_root / "stretch_batch_compare.png")).expanduser().resolve()
    summary_output_file = output_file.with_name(f"{output_file.stem}_force_window{output_file.suffix}")

    data_files, skipped_dirs = discover_output_data_files(root)
    directory_jobs = resolve_max_workers(args.jobs, max(len(data_files), 1))
    if args.analysis_jobs is not None:
        analysis_jobs = max(1, args.analysis_jobs)
    elif directory_jobs == 1:
        analysis_jobs = None
    else:
        analysis_jobs = 1

    print(f"Root directory : {root}")
    print(f"Analysis root  : {analysis_root}")
    print(f"Time canvas    : {output_file}")
    print(f"Force canvas   : {summary_output_file}")
    print(f"Directory jobs : {directory_jobs}")
    print(f"Analysis jobs  : {'auto' if analysis_jobs is None else analysis_jobs}")
    print(f"Window points  : {args.window_points}")
    print(f"Output dirs    : {len(data_files)} runnable, {len(skipped_dirs)} skipped")
    if skipped_dirs:
        print("Skipped output dirs:")
        for skipped_dir in skipped_dirs:
            print(f"  - {skipped_dir}")

    completion_results = collect_dataset_completion_results(
        data_files=data_files,
        directory_jobs=directory_jobs,
        analysis_jobs=analysis_jobs,
    )
    created_count = sum(1 for result in completion_results if result["status"] == "created")
    reused_count = sum(1 for result in completion_results if result["status"] == "reused")
    wait_listed_count = sum(1 for result in completion_results if result["status"] == "wait-listed")
    if wait_listed_count:
        print("Wait-listed datasets:")
        for result in completion_results:
            if result["status"] == "wait-listed":
                print(f"  - {Path(result['data_file']).parent.name}")

    dataset_paths = discover_cached_dataset_paths(root)
    if not dataset_paths:
        raise SystemExit(
            f"No cached stretch datasets were found in {analysis_root}. "
            "All candidates may be in wait-list-file or still missing."
        )

    print(
        f"Cached datasets: {len(dataset_paths)} top-level files "
        f"({created_count} created, {reused_count} reused, {wait_listed_count} wait-listed)"
    )
    payloads = load_dataset_payloads(dataset_paths, directory_jobs=directory_jobs)
    save_time_canvas(payloads=payloads, output_file=output_file)
    save_force_summary_canvas(payloads=payloads, output_file=summary_output_file, window_points=args.window_points)

    print(f"Wrote time-series comparison canvas to {output_file}")
    print(f"Wrote force-summary comparison canvas to {summary_output_file}")

    if args.show:
        show_png_sequence([output_file, summary_output_file])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
