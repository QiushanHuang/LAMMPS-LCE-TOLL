import re
import sys
import argparse
import textwrap
import os
import multiprocessing as mp
from pathlib import Path
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor

try:
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.ticker import ScalarFormatter
    from matplotlib.lines import Line2D
except ModuleNotFoundError as e:
    raise SystemExit(
        f"缺少依赖模块: {e.name}\n"
        f"当前解释器: {sys.executable}\n"
        f"请运行: {sys.executable} -m pip install numpy matplotlib"
    ) from e


PANEL_TITLE_FONTSIZE = 11
SINGLE_FIGURE_SUPTITLE_FONTSIZE = 13
COMBINED_FIGURE_SUPTITLE_FONTSIZE = 11
DEFAULT_PARTICLE_AXIAL_LENGTH = {1: 3.0, 2: 1.0}
AXIAL_LENGTH_RATIO_OVERRIDE_THRESHOLD = 20.0
S_COMPUTE_AXIAL_LENGTH_DEFAULT = 1.0
PERCENT_AXIS_LABEL = "Relative $R_g$ (% of straight-chain max)"
PERCENT_MEAN_AXIS_LABEL = "Relative mean tail $R_g$ (% of straight-chain max)"
ABSOLUTE_AXIS_LABEL = "Radius of gyration $R_g$"
ABSOLUTE_MEAN_AXIS_LABEL = "Mean tail $R_g$"
ROW_KIND_FRAME = "FRAME"
ROW_KIND_TAIL = "TAIL_SUMMARY"
SUMMARY_SUFFIX_DROP_TOKENS = {"result", "md", "0406", "changedof"}


# =========================
# 1. 基础计算：质量加权 Rg
# =========================
def compute_mass_weighted_rg(positions, masses):
    """
    计算质量加权回转半径 Rg

    Parameters
    ----------
    positions : np.ndarray, shape (N, 3)
        粒子坐标
    masses : np.ndarray, shape (N,)
        粒子质量

    Returns
    -------
    rg : float
        回转半径
    rg2 : float
        回转半径平方
    r_cm : np.ndarray, shape (3,)
        质量加权质心
    """
    total_mass = np.sum(masses)
    if total_mass <= 0:
        raise ValueError("总质量 <= 0，无法计算 Rg。")

    r_cm = np.sum(positions * masses[:, None], axis=0) / total_mass
    dr = positions - r_cm
    rg2 = np.sum(masses * np.sum(dr**2, axis=1)) / total_mass
    rg = np.sqrt(rg2)

    return rg, rg2, r_cm


def compute_linear_chain_max_rg(axial_lengths, masses):
    axial_lengths = np.asarray(axial_lengths, dtype=float)
    masses = np.asarray(masses, dtype=float)

    if axial_lengths.size == 0:
        raise ValueError("没有可用于构建拉直链参考的粒子。")
    if axial_lengths.shape != masses.shape:
        raise ValueError("axial_lengths 与 masses 维度不一致。")
    if np.any(axial_lengths <= 0):
        raise ValueError("粒子轴向长度必须为正数。")

    starts = np.concatenate(([0.0], np.cumsum(axial_lengths[:-1])))
    centers = starts + 0.5 * axial_lengths
    positions = np.column_stack((centers, np.zeros_like(centers), np.zeros_like(centers)))
    rg, _, r_cm = compute_mass_weighted_rg(positions, masses)
    return float(rg), float(np.sum(axial_lengths)), float(r_cm[0])


def resolve_compute_lengths(representative_lengths):
    compute_lengths = dict(representative_lengths)
    raw_e_length = representative_lengths.get(1)
    raw_s_length = representative_lengths.get(2)
    raw_ratio = np.nan
    use_s_override = False

    if raw_e_length is not None and raw_s_length is not None and raw_s_length > 0:
        raw_ratio = float(raw_e_length) / float(raw_s_length)
        if raw_ratio > AXIAL_LENGTH_RATIO_OVERRIDE_THRESHOLD:
            compute_lengths[2] = float(S_COMPUTE_AXIAL_LENGTH_DEFAULT)
            use_s_override = True

    return {
        "compute_lengths": compute_lengths,
        "raw_ratio_e_over_s": raw_ratio,
        "use_s_override": use_s_override,
    }


def read_first_frame_chain_reference(filename, selected_types={1, 2}):
    filename = str(filename)
    selected_types = set(selected_types)

    with open(filename, "r", encoding="utf-8") as f:
        line = f.readline()
        if not line:
            raise ValueError(f"{filename} 是空文件。")
        if line.strip() != "ITEM: TIMESTEP":
            raise ValueError(f"{filename} 文件格式错误，预期 ITEM: TIMESTEP。")

        timestep = int(f.readline().strip())

        line = f.readline().strip()
        if line != "ITEM: NUMBER OF ATOMS":
            raise ValueError(f"{filename} 文件格式错误，预期 ITEM: NUMBER OF ATOMS，实际得到: {line}")
        natoms = int(f.readline().strip())

        line = f.readline().strip()
        if not line.startswith("ITEM: BOX BOUNDS"):
            raise ValueError(f"{filename} 文件格式错误，预期 ITEM: BOX BOUNDS，实际得到: {line}")

        for _ in range(3):
            f.readline()

        line = f.readline().strip()
        if not line.startswith("ITEM: ATOMS"):
            raise ValueError(f"{filename} 文件格式错误，预期 ITEM: ATOMS，实际得到: {line}")

        header = line.split()[2:]
        col_index = {name: idx for idx, name in enumerate(header)}

        for col in ("id", "type", "mass"):
            if col not in col_index:
                raise ValueError(f"{filename} 缺少必要列: {col}")

        has_shape = all(col in col_index for col in ("shapex", "shapey", "shapez"))
        particles = []

        for _ in range(natoms):
            atom_line = f.readline().split()
            ptype = int(float(atom_line[col_index["type"]]))
            if ptype not in selected_types:
                continue

            pid = int(float(atom_line[col_index["id"]]))
            mass = float(atom_line[col_index["mass"]])

            if has_shape:
                axial_length = max(
                    float(atom_line[col_index["shapex"]]),
                    float(atom_line[col_index["shapey"]]),
                    float(atom_line[col_index["shapez"]]),
                )
            else:
                axial_length = DEFAULT_PARTICLE_AXIAL_LENGTH.get(ptype, 1.0)

            particles.append(
                {
                    "id": pid,
                    "type": ptype,
                    "mass": mass,
                    "axial_length": axial_length,
                }
            )

    if len(particles) == 0:
        raise ValueError(f"{filename} 第一帧中没有匹配 selected_types={sorted(selected_types)} 的粒子。")

    particles.sort(key=lambda item: item["id"])
    raw_axial_lengths = np.array([item["axial_length"] for item in particles], dtype=float)
    masses = np.array([item["mass"] for item in particles], dtype=float)
    raw_max_rg, raw_total_length, raw_reference_cm_x = compute_linear_chain_max_rg(raw_axial_lengths, masses)

    type_counts = defaultdict(int)
    type_lengths = defaultdict(list)
    for item in particles:
        type_counts[item["type"]] += 1
        type_lengths[item["type"]].append(item["axial_length"])

    representative_lengths = {
        ptype: float(np.median(lengths))
        for ptype, lengths in type_lengths.items()
        if len(lengths) > 0
    }

    compute_length_info = resolve_compute_lengths(representative_lengths)
    compute_lengths = compute_length_info["compute_lengths"]
    compute_axial_lengths = np.array(
        [float(compute_lengths.get(item["type"], item["axial_length"])) for item in particles],
        dtype=float,
    )
    max_rg, total_length, reference_cm_x = compute_linear_chain_max_rg(compute_axial_lengths, masses)

    return {
        "source_file": filename,
        "timestep": timestep,
        "particle_count": len(particles),
        "total_mass": float(np.sum(masses)),
        "raw_total_length": raw_total_length,
        "total_length": total_length,
        "raw_max_rg": raw_max_rg,
        "max_rg": max_rg,
        "raw_reference_cm_x": raw_reference_cm_x,
        "reference_cm_x": reference_cm_x,
        "type_counts": dict(type_counts),
        "representative_lengths": representative_lengths,
        "compute_lengths": compute_lengths,
        "raw_ratio_e_over_s": compute_length_info["raw_ratio_e_over_s"],
        "use_s_compute_length_override": compute_length_info["use_s_override"],
    }


# ==========================================
# 2. 解析单个 dump 文件中的所有 frame（若多帧）
# ==========================================
def parse_dump_frames(filename, selected_types={1, 2}):
    """
    解析一个 LAMMPS dump 文件中的所有 frame，
    对每一帧计算 type1/type2 的质量加权 Rg。

    Parameters
    ----------
    filename : str or Path
        dump 文件路径
    selected_types : set
        粒子类型，例如 {1, 2}

    Returns
    -------
    frames : list of dict
        每个元素包含:
        {
            "timestep": int,
            "rg": float,
            "rg2": float,
            "count": int,
            "total_mass": float,
            "coord_source": str,   # "unwrap" 或 "wrapped"
            "file": str
        }
    """
    frames = []
    filename = str(filename)

    with open(filename, "r", encoding="utf-8") as f:
        while True:
            line = f.readline()
            if not line:
                break  # EOF

            line = line.strip()
            if line == "":
                continue

            if line != "ITEM: TIMESTEP":
                raise ValueError(f"{filename} 文件格式错误，预期 ITEM: TIMESTEP，实际得到: {line}")

            timestep = int(f.readline().strip())

            line = f.readline().strip()
            if line != "ITEM: NUMBER OF ATOMS":
                raise ValueError(f"{filename} 文件格式错误，预期 ITEM: NUMBER OF ATOMS，实际得到: {line}")
            natoms = int(f.readline().strip())

            line = f.readline().strip()
            if not line.startswith("ITEM: BOX BOUNDS"):
                raise ValueError(f"{filename} 文件格式错误，预期 ITEM: BOX BOUNDS，实际得到: {line}")

            # 跳过 box bounds 三行
            for _ in range(3):
                f.readline()

            line = f.readline().strip()
            if not line.startswith("ITEM: ATOMS"):
                raise ValueError(f"{filename} 文件格式错误，预期 ITEM: ATOMS，实际得到: {line}")

            header = line.split()[2:]
            col_index = {name: idx for idx, name in enumerate(header)}

            required_cols = ["type", "mass"]
            for col in required_cols:
                if col not in col_index:
                    raise ValueError(f"{filename} 缺少必要列: {col}")

            # 坐标优先级：先用 unwrap 坐标 xu/yu/zu；若不存在再回退到 x/y/z
            if all(col in col_index for col in ("xu", "yu", "zu")):
                x_col, y_col, z_col = "xu", "yu", "zu"
                coord_source = "unwrap"
            elif all(col in col_index for col in ("x", "y", "z")):
                x_col, y_col, z_col = "x", "y", "z"
                coord_source = "wrapped"
            else:
                raise ValueError(
                    f"{filename} 缺少坐标列，需至少包含 (xu,yu,zu) 或 (x,y,z)。"
                )

            selected_positions = []
            selected_masses = []
            selected_type_counts = defaultdict(int)

            for _ in range(natoms):
                atom_line = f.readline().split()
                ptype = int(float(atom_line[col_index["type"]]))

                if ptype not in selected_types:
                    continue

                x = float(atom_line[col_index[x_col]])
                y = float(atom_line[col_index[y_col]])
                z = float(atom_line[col_index[z_col]])
                mass = float(atom_line[col_index["mass"]])

                selected_positions.append([x, y, z])
                selected_masses.append(mass)
                selected_type_counts[ptype] += 1

            if len(selected_positions) == 0:
                frames.append({
                    "timestep": timestep,
                    "rg": np.nan,
                    "rg2": np.nan,
                    "count": 0,
                    "total_mass": 0.0,
                    "coord_source": coord_source,
                    "file": filename,
                    "r_cm_x": np.nan,
                    "r_cm_y": np.nan,
                    "r_cm_z": np.nan,
                    "type1_count": 0,
                    "type2_count": 0,
                })
                continue

            positions = np.array(selected_positions, dtype=float)
            masses = np.array(selected_masses, dtype=float)

            rg, rg2, r_cm = compute_mass_weighted_rg(positions, masses)

            frames.append({
                "timestep": timestep,
                "rg": rg,
                "rg2": rg2,
                "count": len(selected_positions),
                "total_mass": float(np.sum(masses)),
                "coord_source": coord_source,
                "file": filename,
                "r_cm_x": float(r_cm[0]),
                "r_cm_y": float(r_cm[1]),
                "r_cm_z": float(r_cm[2]),
                "type1_count": int(selected_type_counts.get(1, 0)),
                "type2_count": int(selected_type_counts.get(2, 0)),
            })

    return frames


# =====================================
# 3. 从路径中提取 Tstar_xxx 文件夹名
# =====================================
def extract_tstar_from_path(path_obj):
    """
    从路径各级目录中提取形如 Tstar_* 的文件夹名
    """
    for part in path_obj.parts:
        if part.startswith("Tstar_"):
            return part
    return None


# ============================================
# 4. 递归搜索整个根目录下的所有目标 dump 文件
# ============================================
def find_all_dump_files(root_dir):
    """
    在 root_dir 下递归搜索:
    */Tstar_*/**/rho030/HEAV.rho030.*.dump
    更稳妥起见，只要路径中有 Tstar_* 且父路径中有 rho030，
    文件名匹配 HEAV.rho030.*.dump 即纳入。

    回退规则：
    若 root_dir 下不存在任何 Tstar_* 相关子文件夹，且 root_dir 当前层
    直接存在 *.dump 文件，则直接处理这些 dump（不再要求 rho030 与文件名模式）。

    Returns
    -------
    grouped_files : dict
        {
            "Tstar_0.60": [Path(...), Path(...), ...],
            ...
        }
    """
    root = Path(root_dir)
    if not root.exists():
        raise FileNotFoundError(f"根目录不存在: {root_dir}")

    grouped_files = defaultdict(list)

    pattern = re.compile(r"^HEAV\.rho030\..*\.dump$")

    for file_path in root.rglob("*.dump"):
        if not file_path.is_file():
            continue

        if not pattern.match(file_path.name):
            continue

        # 路径中必须含有 rho030
        if "rho030" not in file_path.parts:
            continue

        tstar = extract_tstar_from_path(file_path)
        if tstar is None:
            continue

        grouped_files[tstar].append(file_path)

    if len(grouped_files) > 0:
        return grouped_files

    # 没有找到 Tstar 分组时，判断是否可启用“当前目录直读 dump”模式
    has_tstar_subdir = any(p.is_dir() for p in root.rglob("Tstar_*"))
    if not has_tstar_subdir:
        direct_dumps = sorted([p for p in root.glob("*.dump") if p.is_file()])
        if len(direct_dumps) > 0:
            group_name = extract_tstar_from_path(root)
            if group_name is None:
                group_name = f"DIRECT_{root.name}" if root.name else "DIRECT_CURRENT_DIR"
            grouped_files[group_name] = direct_dumps
            print(
                f"[INFO] 未检测到 Tstar_* 子目录，启用当前目录直读模式，"
                f"共 {len(direct_dumps)} 个 dump。"
            )

    return grouped_files


def resolve_max_workers(requested_workers, task_count):
    if task_count <= 0:
        return 1

    if requested_workers is not None:
        return max(1, min(int(requested_workers), task_count))

    cpu_count = os.cpu_count() or 1
    return max(1, min(cpu_count, task_count))


def build_parallel_tasks(grouped_files, selected_types):
    tasks = []
    selected_types_tuple = tuple(sorted(selected_types))

    for tstar in sorted(grouped_files.keys()):
        for file_path in sorted(grouped_files[tstar]):
            tasks.append((tstar, str(file_path), selected_types_tuple))

    return tasks


def create_process_pool(max_workers):
    if sys.platform != "win32":
        available_methods = mp.get_all_start_methods()
        if "fork" in available_methods:
            return ProcessPoolExecutor(
                max_workers=max_workers,
                mp_context=mp.get_context("fork"),
            )
    return ProcessPoolExecutor(max_workers=max_workers)


def parse_dump_file_task(task):
    tstar, file_path, selected_types = task
    try:
        frames = parse_dump_frames(file_path, selected_types=set(selected_types))
        return {
            "tstar": tstar,
            "frames": frames,
            "file_path": file_path,
            "error": None,
        }
    except Exception as exc:
        return {
            "tstar": tstar,
            "frames": None,
            "file_path": file_path,
            "error": str(exc),
        }


def collect_chain_reference(grouped_files, selected_types={1, 2}):
    for tstar in sorted(grouped_files.keys()):
        files = sorted(grouped_files[tstar])
        for file_path in files:
            try:
                return read_first_frame_chain_reference(file_path, selected_types=selected_types)
            except Exception as exc:
                print(f"[WARNING] 参考链信息读取失败 {file_path}，原因: {exc}")
    raise RuntimeError("无法从任何 dump 文件中构建拉直链参考信息。")


# ============================================
# 5. 聚合所有 Tstar 的数据，并按 timestep 排序
# ============================================
def collect_all_rg_data(root_dir, selected_types={1, 2}, max_workers=None, grouped_files=None):
    """
    搜索整个目录并计算所有 Tstar 下的 Rg 数据

    Returns
    -------
    all_data : dict
        {
            "Tstar_0.60": [
                {"timestep":..., "rg":..., ...},
                ...
            ],
            ...
        }
    """
    if grouped_files is None:
        grouped_files = find_all_dump_files(root_dir)

    if len(grouped_files) == 0:
        raise RuntimeError(f"在 {root_dir} 下没有找到符合条件的 dump 文件。")

    for tstar, files in grouped_files.items():
        print(f"\n[INFO] 正在处理 {tstar}，文件数 = {len(files)}")

    tasks = build_parallel_tasks(grouped_files, selected_types)
    worker_count = resolve_max_workers(max_workers, len(tasks))
    mode_label = "串行模式" if worker_count == 1 else "并行模式"
    print(
        f"\n[INFO] {mode_label}：使用 {worker_count} 个进程处理 "
        f"{len(tasks)} 个 dump 文件。"
    )

    all_data = {tstar: [] for tstar in sorted(grouped_files.keys())}

    if worker_count == 1:
        results = map(parse_dump_file_task, tasks)
    else:
        chunksize = max(1, len(tasks) // (worker_count * 4))
        executor = create_process_pool(worker_count)
        results = executor.map(parse_dump_file_task, tasks, chunksize=chunksize)

    try:
        for result in results:
            if result["error"] is not None:
                print(f"[WARNING] 跳过文件 {result['file_path']}，原因: {result['error']}")
                continue
            all_data[result["tstar"]].extend(result["frames"])
    finally:
        if worker_count != 1:
            executor.shutdown(wait=True)

    for tstar in sorted(all_data.keys()):
        all_data[tstar].sort(key=lambda d: d["timestep"])

    return all_data


# ===============================
# 6.1 汇总绘图使用的坐标来源
# ===============================
def summarize_coord_source(all_data):
    sources = set()
    for items in all_data.values():
        for item in items:
            source = item.get("coord_source")
            if source in {"unwrap", "wrapped"}:
                sources.add(source)

    if sources == {"unwrap"}:
        return "unwrap (xu, yu, zu)"
    if sources == {"wrapped"}:
        return "wrapped (x, y, z)"
    if sources == {"unwrap", "wrapped"}:
        return "mixed (prefer unwrap, fallback wrapped)"
    return "unknown"


# ===============================
# 6.2 绘图辅助函数
# ===============================
def extract_tstar_numeric(label):
    match = re.search(r"Tstar_([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)", label)
    if match is None:
        return None
    return float(match.group(1))


def build_curve_xy(items, use_real_time=False, dt=None):
    x = np.array([item["timestep"] for item in items], dtype=float)
    y = np.array([item["rg"] for item in items], dtype=float)

    if use_real_time:
        if dt is None:
            raise ValueError("若 use_real_time=True，则必须提供 dt。")
        x = x * dt

    return x, y


def get_sg_window_length(n_points, desired_window, polyorder):
    if n_points <= polyorder:
        return None

    min_window = polyorder + 1
    if min_window % 2 == 0:
        min_window += 1

    max_window = n_points if n_points % 2 == 1 else n_points - 1
    if max_window < min_window:
        return None

    window = min(desired_window, max_window)
    if window % 2 == 0:
        window -= 1
    if window < min_window:
        window = min_window
    if window > max_window:
        return None

    return window


def resolve_sg_window_length(
    n_points,
    requested_window,
    polyorder,
    auto_fraction=0.25,
    auto_min=201,
    auto_max=1001,
):
    if requested_window is None:
        desired_window = int(np.ceil(n_points * auto_fraction))
        if desired_window % 2 == 0:
            desired_window += 1
        desired_window = max(auto_min, desired_window)
        desired_window = min(auto_max, desired_window)
    else:
        desired_window = requested_window

    return get_sg_window_length(n_points, desired_window, polyorder)


def format_sg_label(window_length, polyorder):
    if window_length is None:
        return f"SG(window=auto~25%, polyorder={polyorder})"
    return f"SG(window={window_length}, polyorder={polyorder})"


def savitzky_golay_smooth(y, window_length=None, polyorder=2):
    y = np.asarray(y, dtype=float)
    if y.size == 0:
        return y.copy()

    finite_mask = np.isfinite(y)
    if np.count_nonzero(finite_mask) <= polyorder:
        return y.copy()

    y_filled = y.copy()
    valid_idx = np.flatnonzero(finite_mask)
    all_idx = np.arange(y.size)
    y_filled[~finite_mask] = np.interp(all_idx[~finite_mask], valid_idx, y[finite_mask])

    actual_window = resolve_sg_window_length(y.size, window_length, polyorder)
    if actual_window is None:
        return y.copy()

    half_window = actual_window // 2
    offsets = np.arange(-half_window, half_window + 1, dtype=float)
    vandermonde = np.vander(offsets, N=polyorder + 1, increasing=True)
    coeffs = np.linalg.pinv(vandermonde)[0]

    y_padded = np.pad(y_filled, (half_window, half_window), mode="edge")
    y_smooth = np.convolve(y_padded, coeffs[::-1], mode="valid")

    result = y_smooth.astype(float)
    result[~finite_mask] = np.nan
    return result


def collect_tail_statistics(all_data, tail_points=3000):
    stats = []

    for tstar, items in sorted(all_data.items()):
        temperature = extract_tstar_numeric(tstar)
        if temperature is None:
            continue

        y = np.array([item["rg"] for item in items], dtype=float)
        y = y[np.isfinite(y)]
        if y.size == 0:
            continue

        tail = y[-tail_points:] if y.size >= tail_points else y
        mean_rg = float(np.mean(tail))
        std_rg = float(np.std(tail, ddof=1)) if tail.size > 1 else 0.0

        stats.append(
            {
                "label": tstar,
                "temperature": temperature,
                "mean_rg": mean_rg,
                "std_rg": std_rg,
                "tail_count": int(tail.size),
            }
        )

    stats.sort(key=lambda item: item["temperature"])
    return stats


def style_curve_axis(ax, use_real_time=False):
    ax.grid(True, which="major", alpha=0.25)
    ax.minorticks_on()
    ax.grid(True, which="minor", alpha=0.12)

    if use_real_time:
        ax.set_xlabel("Time")
    else:
        ax.set_xlabel("Timestep")
        formatter = ScalarFormatter(useMathText=True)
        formatter.set_scientific(True)
        formatter.set_powerlimits((0, 0))
        ax.xaxis.set_major_formatter(formatter)


def format_analysis_title(analysis_name):
    wrapped_name = textwrap.fill(analysis_name, width=60)
    return f"Analyzed folder:\n{wrapped_name}"


def format_combined_title_line(analysis_name):
    if analysis_name:
        return f"All Rg Views | Folder: {analysis_name}"
    return "All Rg Views"


def build_analysis_suffix(analysis_name):
    tokens = re.split(r"[^A-Za-z0-9]+", analysis_name)
    kept = [token.lower() for token in tokens if token and token.lower() not in SUMMARY_SUFFIX_DROP_TOKENS]
    if not kept:
        return "current_dir"
    return "_".join(kept)


def resolve_analysis_name(root_dir):
    root_dir = Path(root_dir).resolve()
    informative_patterns = (
        "result_",
        "md_",
        "changedof",
        "unwrap-length",
        "mass1over",
        "onlyellipsoid",
        "ball_ellipsoid",
        "multibond",
        "anchor",
    )

    for candidate in (root_dir, *root_dir.parents):
        name = candidate.name
        if not name:
            continue
        lower_name = name.lower()
        if lower_name.startswith("tstar_") or lower_name.startswith("rho"):
            continue
        if any(pattern in lower_name for pattern in informative_patterns):
            return name

    return root_dir.name or str(root_dir)


def compute_percent_value(value, reference_max_rg):
    if reference_max_rg is None or not np.isfinite(reference_max_rg) or reference_max_rg <= 0:
        return np.nan
    if value is None or not np.isfinite(value):
        return np.nan
    return float(value) * 100.0 / float(reference_max_rg)


def format_chain_reference_line(chain_reference):
    type_counts = chain_reference.get("type_counts", {})
    representative_lengths = chain_reference.get("representative_lengths", {})
    compute_lengths = chain_reference.get("compute_lengths", {})
    max_rg = chain_reference.get("max_rg", np.nan)
    total_length = chain_reference.get("total_length", np.nan)
    s_count = int(type_counts.get(2, 0))
    e_count = int(type_counts.get(1, 0))
    s_length = representative_lengths.get(2, np.nan)
    e_length = representative_lengths.get(1, np.nan)
    s_compute_length = compute_lengths.get(2, s_length)
    use_s_override = bool(chain_reference.get("use_s_compute_length_override", False))
    s_part = f"S={s_count} x {s_length:.6f}"
    if use_s_override:
        s_part += f" [comp {s_compute_length:.6f}]"
    return (
        f"100% $R_g$={max_rg:.6f} | "
        f"$L_{{chain}}$={total_length:.6f} | "
        f"{s_part} | "
        f"E={e_count} x {e_length:.6f}"
    )


def compute_curve_series(items, chain_reference, use_real_time=False, dt=None, sg_window=None, sg_polyorder=2):
    x, y_raw = build_curve_xy(items, use_real_time=use_real_time, dt=dt)
    y_sg = savitzky_golay_smooth(y_raw, window_length=sg_window, polyorder=sg_polyorder)
    max_rg = chain_reference.get("max_rg")
    raw_percent = np.array([compute_percent_value(value, max_rg) for value in y_raw], dtype=float)
    sg_percent = np.array([compute_percent_value(value, max_rg) for value in y_sg], dtype=float)
    return {
        "x": x,
        "rg_raw": y_raw,
        "rg_sg": y_sg,
        "rg_percent": raw_percent,
        "rg_sg_percent": sg_percent,
    }


def build_tail_summary_rows(all_data, chain_reference, tail_points=3000):
    rows = []
    stats = collect_tail_statistics(all_data, tail_points=tail_points)
    max_rg = chain_reference.get("max_rg")

    for item in stats:
        mean_rg_percent = compute_percent_value(item["mean_rg"], max_rg)
        std_rg_percent = compute_percent_value(item["std_rg"], max_rg)
        rows.append(
            {
                "row_kind": ROW_KIND_TAIL,
                "tstar": item["label"],
                "temperature": item["temperature"],
                "tail_count": item["tail_count"],
                "tail_points": tail_points,
                "tail_mean_rg": item["mean_rg"],
                "tail_mean_rg_percent": mean_rg_percent,
                "tail_std_rg": item["std_rg"],
                "tail_std_rg_percent": std_rg_percent,
            }
        )

    return rows


def attach_dual_y_axis(ax, reference_max_rg, absolute_label, percent_label):
    ax.set_ylabel(absolute_label)
    ax.yaxis.set_label_position("right")
    ax.yaxis.tick_right()
    ax.spines["right"].set_visible(True)

    if reference_max_rg is None or not np.isfinite(reference_max_rg) or reference_max_rg <= 0:
        return None

    def absolute_to_percent(values):
        return np.asarray(values, dtype=float) * 100.0 / reference_max_rg

    def percent_to_absolute(values):
        return np.asarray(values, dtype=float) * reference_max_rg / 100.0

    secax = ax.secondary_yaxis("left", functions=(absolute_to_percent, percent_to_absolute))
    secax.set_ylabel(percent_label)
    return secax


def plot_curve_panel(
    ax,
    all_data,
    coord_desc,
    chain_reference,
    use_real_time=False,
    dt=None,
    apply_sg=False,
    overlay_sg=False,
    sg_window=None,
    sg_polyorder=2,
):
    tstar_names = sorted(all_data.keys())
    cmap = plt.get_cmap("tab20", len(tstar_names))

    for i, tstar in enumerate(tstar_names):
        items = all_data[tstar]
        if len(items) == 0:
            continue

        curve = compute_curve_series(
            items,
            chain_reference,
            use_real_time=use_real_time,
            dt=dt,
            sg_window=sg_window,
            sg_polyorder=sg_polyorder,
        )
        x = curve["x"]
        y_raw = curve["rg_raw"]
        y_sg = curve["rg_sg"]
        y_plot = y_raw

        color = cmap(i)
        if overlay_sg:
            ax.plot(
                x,
                y_raw,
                linewidth=1.0,
                linestyle="-",
                color=color,
                alpha=0.16,
            )
            y_plot = y_sg
            ax.plot(
                x,
                y_plot,
                linewidth=2.4,
                linestyle="--",
                color=color,
                label=tstar,
            )
        else:
            if apply_sg:
                y_plot = y_sg
            ax.plot(x, y_plot, linewidth=1.5, linestyle="-", color=color, label=tstar)

        finite_mask = np.isfinite(y_plot)
        if np.any(finite_mask):
            x_valid = x[finite_mask]
            y_valid = y_plot[finite_mask]
            ax.text(
                x_valid[-1],
                y_valid[-1],
                f"  {tstar}",
                color=color,
                fontsize=9,
                va="center",
            )

    if overlay_sg:
        title = (
            f"Figure 2: Original Rg with SG curve overlay\n"
            f"Coordinates: {coord_desc} | raw=faint, SG=bold | "
            f"{format_sg_label(sg_window, sg_polyorder)}"
        )
    elif apply_sg:
        title = (
            f"Figure 2: SG-smoothed Rg curves\n"
            f"Coordinates: {coord_desc} | {format_sg_label(sg_window, sg_polyorder)}\n"
            f"{format_chain_reference_line(chain_reference)}"
        )
    else:
        title = (
            f"Figure 1: Original Rg curves\n"
            f"Coordinates: {coord_desc}\n"
            f"{format_chain_reference_line(chain_reference)}"
        )

    if overlay_sg:
        title = (
            f"Figure 2: Original Rg with SG curve overlay\n"
            f"Coordinates: {coord_desc} | raw=faint, SG=bold | "
            f"{format_sg_label(sg_window, sg_polyorder)}\n"
            f"{format_chain_reference_line(chain_reference)}"
        )

    ax.set_title(title, fontsize=PANEL_TITLE_FONTSIZE)
    style_curve_axis(ax, use_real_time=use_real_time)
    attach_dual_y_axis(
        ax,
        chain_reference.get("max_rg"),
        absolute_label=ABSOLUTE_AXIS_LABEL,
        percent_label=PERCENT_AXIS_LABEL,
    )
    if overlay_sg:
        style_handles = [
            Line2D([0], [0], color="0.5", linewidth=1.0, linestyle="-", alpha=0.35, label="raw"),
            Line2D([0], [0], color="0.2", linewidth=2.2, linestyle="--", label="SG"),
        ]
        ax.legend(handles=style_handles, title="Style", loc="upper right")


def plot_temperature_panel(ax, all_data, chain_reference, tail_points=3000):
    stats = collect_tail_statistics(all_data, tail_points=tail_points)

    if len(stats) == 0:
        ax.text(
            0.5,
            0.5,
            "No numeric Tstar values found.\nFigure 3 is unavailable for this dataset.",
            ha="center",
            va="center",
            fontsize=11,
            transform=ax.transAxes,
        )
        ax.set_axis_off()
        return

    temperatures = np.array([item["temperature"] for item in stats], dtype=float)
    mean_rg = np.array([item["mean_rg"] for item in stats], dtype=float)
    std_rg = np.array([item["std_rg"] for item in stats], dtype=float)

    ax.errorbar(
        temperatures,
        mean_rg,
        yerr=std_rg,
        fmt="o-",
        linewidth=1.5,
        markersize=5,
        capsize=4,
        color="tab:blue",
        ecolor="tab:gray",
    )

    for item in stats:
        ax.text(
            item["temperature"],
            item["mean_rg"],
            f"  {item['label']}",
            fontsize=8,
            va="center",
        )

    ax.set_xlabel("Temperature $T^*$")
    ax.set_title(
        f"Figure 3: Mean $R_g$ of last {tail_points} points vs $T^*$\n"
        f"Error bar = standard deviation of the tail window\n"
        f"{format_chain_reference_line(chain_reference)}",
        fontsize=PANEL_TITLE_FONTSIZE,
    )
    ax.grid(True, which="major", alpha=0.25)
    ax.minorticks_on()
    ax.grid(True, which="minor", alpha=0.12)
    attach_dual_y_axis(
        ax,
        chain_reference.get("max_rg"),
        absolute_label=ABSOLUTE_MEAN_AXIS_LABEL,
        percent_label=PERCENT_MEAN_AXIS_LABEL,
    )


def save_three_single_plots(
    all_data,
    out_dir,
    coord_desc,
    chain_reference,
    use_real_time=False,
    dt=None,
    sg_window=None,
    sg_polyorder=2,
    tail_points=3000,
    analysis_name="",
):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if use_real_time:
        fig1_name = "figure1_rg_vs_time_all_Tstar.png"
        fig2_name = "figure2_rg_vs_time_all_Tstar_sg.png"
    else:
        fig1_name = "figure1_rg_vs_timestep_all_Tstar.png"
        fig2_name = "figure2_rg_vs_timestep_all_Tstar_sg.png"
    fig3_name = "figure3_rg_tail_mean_vs_temperature.png"

    fig1, ax1 = plt.subplots(figsize=(16, 7))
    plot_curve_panel(
        ax1,
        all_data,
        coord_desc,
        chain_reference,
        use_real_time=use_real_time,
        dt=dt,
        apply_sg=False,
        sg_window=sg_window,
        sg_polyorder=sg_polyorder,
    )
    fig1.suptitle(format_analysis_title(analysis_name), fontsize=SINGLE_FIGURE_SUPTITLE_FONTSIZE)
    fig1.tight_layout(rect=(0, 0, 1, 0.90))
    fig1.savefig(out_dir / fig1_name, dpi=200)
    plt.close(fig1)

    fig2, ax2 = plt.subplots(figsize=(16, 7))
    plot_curve_panel(
        ax2,
        all_data,
        coord_desc,
        chain_reference,
        use_real_time=use_real_time,
        dt=dt,
        overlay_sg=True,
        sg_window=sg_window,
        sg_polyorder=sg_polyorder,
    )
    fig2.suptitle(format_analysis_title(analysis_name), fontsize=SINGLE_FIGURE_SUPTITLE_FONTSIZE)
    fig2.tight_layout(rect=(0, 0, 1, 0.90))
    fig2.savefig(out_dir / fig2_name, dpi=200)
    plt.close(fig2)

    fig3, ax3 = plt.subplots(figsize=(16, 7))
    plot_temperature_panel(ax3, all_data, chain_reference, tail_points=tail_points)
    fig3.suptitle(format_analysis_title(analysis_name), fontsize=SINGLE_FIGURE_SUPTITLE_FONTSIZE)
    fig3.tight_layout(rect=(0, 0, 1, 0.90))
    fig3.savefig(out_dir / fig3_name, dpi=200)
    plt.close(fig3)

    print(f"[INFO] 图片已保存: {out_dir / fig1_name}")
    print(f"[INFO] 图片已保存: {out_dir / fig2_name}")
    print(f"[INFO] 图片已保存: {out_dir / fig3_name}")


# ===============================
# 6. 输出每个 Tstar 的数据文件
# ===============================
def build_frame_export_rows(
    all_data,
    chain_reference,
    sg_window=None,
    sg_polyorder=2,
    dt=None,
):
    rows = []
    max_rg = chain_reference.get("max_rg")
    type_counts = chain_reference.get("type_counts", {})
    representative_lengths = chain_reference.get("representative_lengths", {})
    compute_lengths = chain_reference.get("compute_lengths", {})

    for tstar in sorted(all_data.keys()):
        items = all_data[tstar]
        curve = compute_curve_series(
            items,
            chain_reference,
            use_real_time=dt is not None,
            dt=dt,
            sg_window=sg_window,
            sg_polyorder=sg_polyorder,
        )

        time_values = curve["x"] if dt is not None else np.full(len(items), np.nan, dtype=float)
        temperature = extract_tstar_numeric(tstar)

        for idx, item in enumerate(items):
            rows.append(
                {
                    "row_kind": ROW_KIND_FRAME,
                    "tstar": tstar,
                    "temperature": temperature,
                    "timestep": int(item["timestep"]),
                    "time": float(time_values[idx]) if idx < len(time_values) else np.nan,
                    "rg": float(item["rg"]) if np.isfinite(item["rg"]) else np.nan,
                    "rg_percent": float(curve["rg_percent"][idx]) if idx < len(curve["rg_percent"]) else np.nan,
                    "rg2": float(item["rg2"]) if np.isfinite(item["rg2"]) else np.nan,
                    "sg_rg": float(curve["rg_sg"][idx]) if idx < len(curve["rg_sg"]) and np.isfinite(curve["rg_sg"][idx]) else np.nan,
                    "sg_rg_percent": float(curve["rg_sg_percent"][idx]) if idx < len(curve["rg_sg_percent"]) else np.nan,
                    "count": int(item["count"]),
                    "type1_count": int(item.get("type1_count", 0)),
                    "type2_count": int(item.get("type2_count", 0)),
                    "frame_total_mass": float(item["total_mass"]),
                    "coord_source": item.get("coord_source", "unknown"),
                    "r_cm_x": float(item.get("r_cm_x", np.nan)),
                    "r_cm_y": float(item.get("r_cm_y", np.nan)),
                    "r_cm_z": float(item.get("r_cm_z", np.nan)),
                    "raw_straight_chain_max_rg": float(chain_reference.get("raw_max_rg", np.nan)),
                    "straight_chain_max_rg": float(max_rg),
                    "raw_straight_chain_cm_x": float(chain_reference.get("raw_reference_cm_x", np.nan)),
                    "straight_chain_cm_x": float(chain_reference.get("reference_cm_x", np.nan)),
                    "raw_contour_length": float(chain_reference.get("raw_total_length", np.nan)),
                    "contour_length": float(chain_reference.get("total_length", np.nan)),
                    "reference_particle_count": int(chain_reference.get("particle_count", 0)),
                    "reference_total_mass": float(chain_reference.get("total_mass", np.nan)),
                    "raw_ratio_e_over_s": float(chain_reference.get("raw_ratio_e_over_s", np.nan)),
                    "use_s_compute_length_override": int(bool(chain_reference.get("use_s_compute_length_override", False))),
                    "e_count": int(type_counts.get(1, 0)),
                    "s_count": int(type_counts.get(2, 0)),
                    "e_length": float(representative_lengths.get(1, np.nan)),
                    "s_length": float(representative_lengths.get(2, np.nan)),
                    "e_compute_length": float(compute_lengths.get(1, representative_lengths.get(1, np.nan))),
                    "s_compute_length": float(compute_lengths.get(2, representative_lengths.get(2, np.nan))),
                    "tail_points": np.nan,
                    "tail_count": np.nan,
                    "tail_mean_rg": np.nan,
                    "tail_mean_rg_percent": np.nan,
                    "tail_std_rg": np.nan,
                    "tail_std_rg_percent": np.nan,
                    "sg_window": sg_window if sg_window is not None else np.nan,
                    "sg_polyorder": int(sg_polyorder),
                    "dt": float(dt) if dt is not None else np.nan,
                    "reference_file": chain_reference.get("source_file", ""),
                    "source_file": item.get("file", ""),
                }
            )

    return rows


def enrich_tail_export_rows(
    tail_rows,
    chain_reference,
    sg_window=None,
    sg_polyorder=2,
    dt=None,
):
    rows = []
    type_counts = chain_reference.get("type_counts", {})
    representative_lengths = chain_reference.get("representative_lengths", {})
    compute_lengths = chain_reference.get("compute_lengths", {})

    for item in tail_rows:
        rows.append(
            {
                "row_kind": ROW_KIND_TAIL,
                "tstar": item["tstar"],
                "temperature": item["temperature"],
                "timestep": np.nan,
                "time": np.nan,
                "rg": np.nan,
                "rg_percent": np.nan,
                "rg2": np.nan,
                "sg_rg": np.nan,
                "sg_rg_percent": np.nan,
                "count": np.nan,
                "type1_count": np.nan,
                "type2_count": np.nan,
                "frame_total_mass": np.nan,
                "coord_source": "summary",
                "r_cm_x": np.nan,
                "r_cm_y": np.nan,
                "r_cm_z": np.nan,
                "raw_straight_chain_max_rg": float(chain_reference.get("raw_max_rg", np.nan)),
                "straight_chain_max_rg": float(chain_reference.get("max_rg", np.nan)),
                "raw_straight_chain_cm_x": float(chain_reference.get("raw_reference_cm_x", np.nan)),
                "straight_chain_cm_x": float(chain_reference.get("reference_cm_x", np.nan)),
                "raw_contour_length": float(chain_reference.get("raw_total_length", np.nan)),
                "contour_length": float(chain_reference.get("total_length", np.nan)),
                "reference_particle_count": int(chain_reference.get("particle_count", 0)),
                "reference_total_mass": float(chain_reference.get("total_mass", np.nan)),
                "raw_ratio_e_over_s": float(chain_reference.get("raw_ratio_e_over_s", np.nan)),
                "use_s_compute_length_override": int(bool(chain_reference.get("use_s_compute_length_override", False))),
                "e_count": int(type_counts.get(1, 0)),
                "s_count": int(type_counts.get(2, 0)),
                "e_length": float(representative_lengths.get(1, np.nan)),
                "s_length": float(representative_lengths.get(2, np.nan)),
                "e_compute_length": float(compute_lengths.get(1, representative_lengths.get(1, np.nan))),
                "s_compute_length": float(compute_lengths.get(2, representative_lengths.get(2, np.nan))),
                "tail_points": int(item["tail_points"]),
                "tail_count": int(item["tail_count"]),
                "tail_mean_rg": float(item["tail_mean_rg"]),
                "tail_mean_rg_percent": float(item["tail_mean_rg_percent"]),
                "tail_std_rg": float(item["tail_std_rg"]),
                "tail_std_rg_percent": float(item["tail_std_rg_percent"]),
                "sg_window": sg_window if sg_window is not None else np.nan,
                "sg_polyorder": int(sg_polyorder),
                "dt": float(dt) if dt is not None else np.nan,
                "reference_file": chain_reference.get("source_file", ""),
                "source_file": "",
            }
        )

    return rows


def write_export_table(outfile, rows, chain_reference, selected_types, sg_window, sg_polyorder, tail_points, dt, analysis_name=""):
    columns = [
        "row_kind",
        "tstar",
        "temperature",
        "timestep",
        "time",
        "rg",
        "rg_percent",
        "rg2",
        "sg_rg",
        "sg_rg_percent",
        "count",
        "type1_count",
        "type2_count",
        "frame_total_mass",
        "coord_source",
        "r_cm_x",
        "r_cm_y",
        "r_cm_z",
        "raw_straight_chain_max_rg",
        "straight_chain_max_rg",
        "raw_straight_chain_cm_x",
        "straight_chain_cm_x",
        "raw_contour_length",
        "contour_length",
        "reference_particle_count",
        "reference_total_mass",
        "raw_ratio_e_over_s",
        "use_s_compute_length_override",
        "e_count",
        "s_count",
        "e_length",
        "s_length",
        "e_compute_length",
        "s_compute_length",
        "tail_points",
        "tail_count",
        "tail_mean_rg",
        "tail_mean_rg_percent",
        "tail_std_rg",
        "tail_std_rg_percent",
        "sg_window",
        "sg_polyorder",
        "dt",
        "reference_file",
        "source_file",
    ]

    with open(outfile, "w", encoding="utf-8") as f:
        f.write("# Complete export for reconstructing Figure 1, Figure 2, and Figure 3 directly from this file.\n")
        f.write(f"# analysis_name = {analysis_name}\n")
        f.write(f"# selected_types = {sorted(selected_types)}\n")
        f.write(f"# raw_straight_chain_max_rg = {chain_reference.get('raw_max_rg', np.nan):.8f}\n")
        f.write(f"# straight_chain_max_rg = {chain_reference.get('max_rg', np.nan):.8f}\n")
        f.write(f"# raw_straight_chain_cm_x = {chain_reference.get('raw_reference_cm_x', np.nan):.8f}\n")
        f.write(f"# straight_chain_cm_x = {chain_reference.get('reference_cm_x', np.nan):.8f}\n")
        f.write(f"# raw_contour_length = {chain_reference.get('raw_total_length', np.nan):.8f}\n")
        f.write(f"# contour_length = {chain_reference.get('total_length', np.nan):.8f}\n")
        f.write(f"# raw_ratio_e_over_s = {chain_reference.get('raw_ratio_e_over_s', np.nan):.8f}\n")
        f.write(f"# use_s_compute_length_override = {int(bool(chain_reference.get('use_s_compute_length_override', False)))}\n")
        f.write(f"# e_count = {int(chain_reference.get('type_counts', {}).get(1, 0))}\n")
        f.write(f"# s_count = {int(chain_reference.get('type_counts', {}).get(2, 0))}\n")
        f.write(f"# e_length = {float(chain_reference.get('representative_lengths', {}).get(1, np.nan)):.8f}\n")
        f.write(f"# s_length = {float(chain_reference.get('representative_lengths', {}).get(2, np.nan)):.8f}\n")
        f.write(f"# e_compute_length = {float(chain_reference.get('compute_lengths', {}).get(1, chain_reference.get('representative_lengths', {}).get(1, np.nan))):.8f}\n")
        f.write(f"# s_compute_length = {float(chain_reference.get('compute_lengths', {}).get(2, chain_reference.get('representative_lengths', {}).get(2, np.nan))):.8f}\n")
        f.write(f"# sg_window = {sg_window if sg_window is not None else 'auto'}\n")
        f.write(f"# sg_polyorder = {sg_polyorder}\n")
        f.write(f"# tail_points = {tail_points}\n")
        f.write(f"# dt = {dt if dt is not None else 'nan'}\n")
        f.write("# " + " ".join(columns) + "\n")

        for row in rows:
            values = []
            for col in columns:
                value = row.get(col, np.nan)
                if isinstance(value, str):
                    values.append(value)
                elif isinstance(value, (int, np.integer)):
                    values.append(str(int(value)))
                elif value is None or (isinstance(value, float) and not np.isfinite(value)):
                    values.append("nan")
                elif isinstance(value, (float, np.floating)):
                    values.append(f"{float(value):.8f}")
                else:
                    values.append(str(value))
            f.write(" ".join(values) + "\n")


def save_all_data(all_data, out_dir, chain_reference, selected_types, sg_window=None, sg_polyorder=2, tail_points=3000, dt=None, analysis_name=""):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    frame_rows = build_frame_export_rows(
        all_data,
        chain_reference,
        sg_window=sg_window,
        sg_polyorder=sg_polyorder,
        dt=dt,
    )
    tail_rows = enrich_tail_export_rows(
        build_tail_summary_rows(all_data, chain_reference, tail_points=tail_points),
        chain_reference,
        sg_window=sg_window,
        sg_polyorder=sg_polyorder,
        dt=dt,
    )
    all_rows = frame_rows + tail_rows

    summary_file = out_dir / "Rg_all_Tstar.dat"
    summary_file_with_suffix = out_dir / f"Rg_all_Tstar_{build_analysis_suffix(analysis_name)}.dat"
    write_export_table(
        summary_file,
        all_rows,
        chain_reference,
        selected_types,
        sg_window,
        sg_polyorder,
        tail_points,
        dt,
        analysis_name=analysis_name,
    )
    write_export_table(
        summary_file_with_suffix,
        all_rows,
        chain_reference,
        selected_types,
        sg_window,
        sg_polyorder,
        tail_points,
        dt,
        analysis_name=analysis_name,
    )

    for tstar in sorted(all_data.keys()):
        outfile = out_dir / f"{tstar}_Rg.dat"
        tstar_rows = [row for row in frame_rows if row["tstar"] == tstar]
        tstar_rows.extend(row for row in tail_rows if row["tstar"] == tstar)
        write_export_table(
            outfile,
            tstar_rows,
            chain_reference,
            selected_types,
            sg_window,
            sg_polyorder,
            tail_points,
            dt,
            analysis_name=analysis_name,
        )

    print(f"\n[INFO] 数据已输出到: {out_dir}")
    print(f"[INFO] 总表文件: {summary_file}")
    print(f"[INFO] 带后缀总表文件: {summary_file_with_suffix}")


# ============================================
# 7. 绘图：输出三张单图和一张聚合图
# ============================================
def plot_all_tstar(
    all_data,
    out_dir,
    chain_reference,
    use_real_time=False,
    dt=None,
    show_plot=False,
    sg_window=None,
    sg_polyorder=2,
    tail_points=3000,
    analysis_name="",
):
    """
    输出：
    1. 原始 Rg 曲线图
    2. SG 滤波后的 Rg 曲线图
    3. 各 Tstar 末 1000 点平均值与温度关系图
    4. 包含上述三个视图的聚合图
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    coord_desc = summarize_coord_source(all_data)

    save_three_single_plots(
        all_data,
        out_dir,
        coord_desc,
        chain_reference,
        use_real_time=use_real_time,
        dt=dt,
        sg_window=sg_window,
        sg_polyorder=sg_polyorder,
        tail_points=tail_points,
        analysis_name=analysis_name,
    )

    fig, axes = plt.subplots(3, 1, figsize=(16, 19.5))
    plot_curve_panel(
        axes[0],
        all_data,
        coord_desc,
        chain_reference,
        use_real_time=use_real_time,
        dt=dt,
        apply_sg=False,
        sg_window=sg_window,
        sg_polyorder=sg_polyorder,
    )
    plot_curve_panel(
        axes[1],
        all_data,
        coord_desc,
        chain_reference,
        use_real_time=use_real_time,
        dt=dt,
        overlay_sg=True,
        sg_window=sg_window,
        sg_polyorder=sg_polyorder,
    )
    plot_temperature_panel(axes[2], all_data, chain_reference, tail_points=tail_points)

    fig.suptitle(
        format_combined_title_line(analysis_name),
        fontsize=COMBINED_FIGURE_SUPTITLE_FONTSIZE,
        y=0.985,
    )
    fig.tight_layout(rect=(0.03, 0.03, 0.98, 0.93), h_pad=3.8)
    fig.subplots_adjust(top=0.90, hspace=0.42)
    combined_name = "figure_all_rg_views.png"
    fig.savefig(out_dir / combined_name, dpi=200)
    if show_plot:
        plt.show()
    plt.close(fig)

    print(f"[INFO] 图片已保存: {out_dir / combined_name}")


def print_chain_reference_summary(chain_reference):
    type_counts = chain_reference.get("type_counts", {})
    representative_lengths = chain_reference.get("representative_lengths", {})
    compute_lengths = chain_reference.get("compute_lengths", {})

    e_count = type_counts.get(1, 0)
    s_count = type_counts.get(2, 0)
    e_length = representative_lengths.get(1)
    s_length = representative_lengths.get(2)
    e_compute_length = compute_lengths.get(1, e_length)
    s_compute_length = compute_lengths.get(2, s_length)

    print("\n[INFO] 拉直链参考:")
    print(f"[INFO]   source_file = {chain_reference['source_file']}")
    print(f"[INFO]   particle_count = {chain_reference['particle_count']}")
    print(f"[INFO]   E_count(type 1) = {e_count}")
    print(f"[INFO]   S_count(type 2) = {s_count}")
    if e_length is not None:
        print(f"[INFO]   E axial length = {e_length:.6f}")
    if s_length is not None:
        print(f"[INFO]   S axial length = {s_length:.6f}")
    if e_compute_length is not None and (e_length is None or abs(e_compute_length - e_length) > 1e-12):
        print(f"[INFO]   E compute axial length = {e_compute_length:.6f}")
    if s_compute_length is not None and (s_length is None or abs(s_compute_length - s_length) > 1e-12):
        print(f"[INFO]   S compute axial length = {s_compute_length:.6f}")
    if np.isfinite(chain_reference.get("raw_ratio_e_over_s", np.nan)):
        print(f"[INFO]   raw E/S axial ratio = {chain_reference['raw_ratio_e_over_s']:.6f}")
    if "raw_total_length" in chain_reference:
        print(f"[INFO]   raw_contour_length = {chain_reference['raw_total_length']:.6f}")
    print(f"[INFO]   contour_length = {chain_reference['total_length']:.6f}")
    if "raw_max_rg" in chain_reference:
        print(f"[INFO]   raw_straight_chain_max_rg = {chain_reference['raw_max_rg']:.6f}")
    print(f"[INFO]   straight_chain_max_rg = {chain_reference['max_rg']:.6f}")


# ==========================
# 8. 主程序
# ==========================
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            "计算所有可发现 dump 的质量加权 Rg。默认处理当前工作目录，"
            "优先按 Tstar_* 规则分组；若无 Tstar_* 子目录且当前目录有 dump，"
            "则直接处理当前目录 dump。"
        )
    )
    parser.add_argument(
        "--root-dir",
        default=".",
        help="待处理根目录，默认当前工作目录。",
    )
    parser.add_argument(
        "--out-dir",
        default="Rg_results",
        help="输出目录。若给相对路径，则相对于 root-dir。",
    )
    parser.add_argument(
        "--types",
        type=int,
        nargs="+",
        default=[1, 2],
        help="参与 Rg 计算的粒子类型，默认: 1 2",
    )
    parser.add_argument(
        "--use-real-time",
        action="store_true",
        help="开启后，横坐标用时间而非 timestep。",
    )
    parser.add_argument(
        "--dt",
        type=float,
        default=None,
        help="每个 timestep 对应的物理时间。仅在 --use-real-time 时生效。",
    )
    parser.add_argument(
        "--show",
        action="store_true",
        help="显示绘图窗口（默认不弹窗，仅保存图片）。",
    )
    parser.add_argument(
        "--sg-window",
        type=int,
        default=2001,
        help="SG 滤波窗口长度，默认 2001。程序会自动调整为合法奇数并受轨迹长度限制。",
    )
    parser.add_argument(
        "--sg-polyorder",
        type=int,
        default=2,
        help="SG 滤波多项式阶数，默认 2。",
    )
    parser.add_argument(
        "--tail-points",
        type=int,
        default=3000,
        help="图 3 末尾平均窗口大小，默认 3000。",
    )
    parser.add_argument(
        "--jobs",
        type=int,
        default=None,
        help="并行进程数。默认自动使用可用 CPU 核心数；设为 1 可禁用并行。",
    )
    args = parser.parse_args()

    root_dir = Path(args.root_dir).resolve()
    out_dir = Path(args.out_dir)
    if not out_dir.is_absolute():
        out_dir = root_dir / out_dir

    selected_types = set(args.types)
    use_real_time = args.use_real_time
    dt = args.dt
    analysis_name = resolve_analysis_name(root_dir)

    print(f"[INFO] root_dir = {root_dir}")
    print(f"[INFO] out_dir  = {out_dir}")
    print(f"[INFO] types    = {sorted(selected_types)}")
    print(f"[INFO] jobs     = {'auto' if args.jobs is None else args.jobs}")

    grouped_files = find_all_dump_files(root_dir)
    chain_reference = collect_chain_reference(grouped_files, selected_types=selected_types)
    print_chain_reference_summary(chain_reference)

    all_data = collect_all_rg_data(
        root_dir,
        selected_types=selected_types,
        max_workers=args.jobs,
        grouped_files=grouped_files,
    )
    save_all_data(
        all_data,
        out_dir,
        chain_reference,
        selected_types=selected_types,
        sg_window=args.sg_window,
        sg_polyorder=args.sg_polyorder,
        tail_points=args.tail_points,
        dt=dt,
        analysis_name=analysis_name,
    )
    plot_all_tstar(
        all_data,
        out_dir,
        chain_reference,
        use_real_time=use_real_time,
        dt=dt,
        show_plot=args.show,
        sg_window=args.sg_window,
        sg_polyorder=args.sg_polyorder,
        tail_points=args.tail_points,
        analysis_name=analysis_name,
    )

    print("\n===== 全部完成 =====")
    print(f"共处理 Tstar 数量: {len(all_data)}")
    for tstar in sorted(all_data.keys()):
        print(f"{tstar}: {len(all_data[tstar])} frames")
