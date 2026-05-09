import re
import sys
import argparse
from pathlib import Path
from collections import defaultdict

try:
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.ticker import ScalarFormatter
except ModuleNotFoundError as e:
    raise SystemExit(
        f"缺少依赖模块: {e.name}\n"
        f"当前解释器: {sys.executable}\n"
        f"请运行: {sys.executable} -m pip install numpy matplotlib"
    ) from e


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

            if len(selected_positions) == 0:
                frames.append({
                    "timestep": timestep,
                    "rg": np.nan,
                    "rg2": np.nan,
                    "count": 0,
                    "total_mass": 0.0,
                    "coord_source": coord_source,
                    "file": filename
                })
                continue

            positions = np.array(selected_positions, dtype=float)
            masses = np.array(selected_masses, dtype=float)

            rg, rg2, _ = compute_mass_weighted_rg(positions, masses)

            frames.append({
                "timestep": timestep,
                "rg": rg,
                "rg2": rg2,
                "count": len(selected_positions),
                "total_mass": float(np.sum(masses)),
                "coord_source": coord_source,
                "file": filename
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


# ============================================
# 5. 聚合所有 Tstar 的数据，并按 timestep 排序
# ============================================
def collect_all_rg_data(root_dir, selected_types={1, 2}):
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
    grouped_files = find_all_dump_files(root_dir)

    if len(grouped_files) == 0:
        raise RuntimeError(f"在 {root_dir} 下没有找到符合条件的 dump 文件。")

    all_data = {}

    for tstar, files in grouped_files.items():
        print(f"\n[INFO] 正在处理 {tstar}，文件数 = {len(files)}")

        frames_all = []
        for fp in sorted(files):
            try:
                frames = parse_dump_frames(fp, selected_types=selected_types)
                frames_all.extend(frames)
            except Exception as e:
                print(f"[WARNING] 跳过文件 {fp}，原因: {e}")

        # 按 timestep 排序
        frames_all.sort(key=lambda d: d["timestep"])

        all_data[tstar] = frames_all

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


def savitzky_golay_smooth(y, window_length=51, polyorder=3):
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

    actual_window = get_sg_window_length(y.size, window_length, polyorder)
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


def collect_tail_statistics(all_data, tail_points=1000):
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
    ax.set_ylabel("Radius of gyration $R_g$")
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


def plot_curve_panel(
    ax,
    all_data,
    coord_desc,
    use_real_time=False,
    dt=None,
    apply_sg=False,
    overlay_sg=False,
    sg_window=51,
    sg_polyorder=3,
):
    tstar_names = sorted(all_data.keys())
    cmap = plt.get_cmap("tab20", len(tstar_names))

    for i, tstar in enumerate(tstar_names):
        items = all_data[tstar]
        if len(items) == 0:
            continue

        x, y_raw = build_curve_xy(items, use_real_time=use_real_time, dt=dt)
        y_plot = y_raw
        if apply_sg or overlay_sg:
            y_sg = savitzky_golay_smooth(y_raw, window_length=sg_window, polyorder=sg_polyorder)
        else:
            y_sg = None

        color = cmap(i)
        if overlay_sg:
            ax.plot(
                x,
                y_raw,
                linewidth=1.0,
                linestyle="-",
                color=color,
                alpha=0.25,
            )
            y_plot = y_sg
            ax.plot(
                x,
                y_plot,
                linewidth=2.0,
                linestyle="-",
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
            f"SG(window={sg_window}, polyorder={sg_polyorder})"
        )
    elif apply_sg:
        title = (
            f"Figure 2: SG-smoothed Rg curves\n"
            f"Coordinates: {coord_desc} | SG(window={sg_window}, polyorder={sg_polyorder})"
        )
    else:
        title = f"Figure 1: Original Rg curves\nCoordinates: {coord_desc}"

    ax.set_title(title)
    style_curve_axis(ax, use_real_time=use_real_time)
    ax.legend(title="$T^*$", loc="best")


def plot_temperature_panel(ax, all_data, tail_points=1000):
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
    ax.set_ylabel("Mean tail $R_g$")
    ax.set_title(
        f"Figure 3: Mean $R_g$ of last {tail_points} points vs $T^*$\n"
        f"Error bar = standard deviation of the tail window"
    )
    ax.grid(True, which="major", alpha=0.25)
    ax.minorticks_on()
    ax.grid(True, which="minor", alpha=0.12)


def save_three_single_plots(
    all_data,
    out_dir,
    coord_desc,
    use_real_time=False,
    dt=None,
    sg_window=51,
    sg_polyorder=3,
    tail_points=1000,
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

    fig1, ax1 = plt.subplots(figsize=(10, 6))
    plot_curve_panel(
        ax1,
        all_data,
        coord_desc,
        use_real_time=use_real_time,
        dt=dt,
        apply_sg=False,
        sg_window=sg_window,
        sg_polyorder=sg_polyorder,
    )
    fig1.tight_layout()
    fig1.savefig(out_dir / fig1_name, dpi=200)
    plt.close(fig1)

    fig2, ax2 = plt.subplots(figsize=(10, 6))
    plot_curve_panel(
        ax2,
        all_data,
        coord_desc,
        use_real_time=use_real_time,
        dt=dt,
        overlay_sg=True,
        sg_window=sg_window,
        sg_polyorder=sg_polyorder,
    )
    fig2.tight_layout()
    fig2.savefig(out_dir / fig2_name, dpi=200)
    plt.close(fig2)

    fig3, ax3 = plt.subplots(figsize=(10, 6))
    plot_temperature_panel(ax3, all_data, tail_points=tail_points)
    fig3.tight_layout()
    fig3.savefig(out_dir / fig3_name, dpi=200)
    plt.close(fig3)

    print(f"[INFO] 图片已保存: {out_dir / fig1_name}")
    print(f"[INFO] 图片已保存: {out_dir / fig2_name}")
    print(f"[INFO] 图片已保存: {out_dir / fig3_name}")


# ===============================
# 6. 输出每个 Tstar 的数据文件
# ===============================
def save_all_data(all_data, out_dir):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # 总表
    summary_file = out_dir / "Rg_all_Tstar.dat"
    with open(summary_file, "w", encoding="utf-8") as f:
        f.write("# Tstar    timestep    Rg    Rg2    count    total_mass    coord_source    file\n")
        for tstar in sorted(all_data.keys()):
            for item in all_data[tstar]:
                f.write(
                    f"{tstar:15s} "
                    f"{item['timestep']:12d} "
                    f"{item['rg']:15.8f} "
                    f"{item['rg2']:15.8f} "
                    f"{item['count']:8d} "
                    f"{item['total_mass']:15.8f} "
                    f"{item.get('coord_source', 'unknown'):12s} "
                    f"{item['file']}\n"
                )

    # 每个 Tstar 单独一个文件
    for tstar in sorted(all_data.keys()):
        outfile = out_dir / f"{tstar}_Rg.dat"
        with open(outfile, "w", encoding="utf-8") as f:
            f.write("# timestep    Rg    Rg2    count    total_mass    coord_source    file\n")
            for item in all_data[tstar]:
                f.write(
                    f"{item['timestep']:12d} "
                    f"{item['rg']:15.8f} "
                    f"{item['rg2']:15.8f} "
                    f"{item['count']:8d} "
                    f"{item['total_mass']:15.8f} "
                    f"{item.get('coord_source', 'unknown'):12s} "
                    f"{item['file']}\n"
                )

    print(f"\n[INFO] 数据已输出到: {out_dir}")


# ============================================
# 7. 绘图：输出三张单图和一张聚合图
# ============================================
def plot_all_tstar(
    all_data,
    out_dir,
    use_real_time=False,
    dt=None,
    show_plot=False,
    sg_window=51,
    sg_polyorder=3,
    tail_points=1000,
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
        use_real_time=use_real_time,
        dt=dt,
        sg_window=sg_window,
        sg_polyorder=sg_polyorder,
        tail_points=tail_points,
    )

    fig, axes = plt.subplots(3, 1, figsize=(12, 18))
    plot_curve_panel(
        axes[0],
        all_data,
        coord_desc,
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
        use_real_time=use_real_time,
        dt=dt,
        overlay_sg=True,
        sg_window=sg_window,
        sg_polyorder=sg_polyorder,
    )
    plot_temperature_panel(axes[2], all_data, tail_points=tail_points)

    fig.suptitle("All Rg Views", fontsize=14)
    fig.tight_layout()
    fig.subplots_adjust(top=0.96)
    combined_name = "figure_all_rg_views.png"
    fig.savefig(out_dir / combined_name, dpi=200)
    if show_plot:
        plt.show()
    plt.close(fig)

    print(f"[INFO] 图片已保存: {out_dir / combined_name}")


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
        default=51,
        help="SG 滤波窗口长度，默认 51。程序会自动调整为合法奇数。",
    )
    parser.add_argument(
        "--sg-polyorder",
        type=int,
        default=3,
        help="SG 滤波多项式阶数，默认 3。",
    )
    parser.add_argument(
        "--tail-points",
        type=int,
        default=1000,
        help="图 3 末尾平均窗口大小，默认 1000。",
    )
    args = parser.parse_args()

    root_dir = Path(args.root_dir).resolve()
    out_dir = Path(args.out_dir)
    if not out_dir.is_absolute():
        out_dir = root_dir / out_dir

    selected_types = set(args.types)
    use_real_time = args.use_real_time
    dt = args.dt

    print(f"[INFO] root_dir = {root_dir}")
    print(f"[INFO] out_dir  = {out_dir}")
    print(f"[INFO] types    = {sorted(selected_types)}")

    all_data = collect_all_rg_data(root_dir, selected_types=selected_types)
    save_all_data(all_data, out_dir)
    plot_all_tstar(
        all_data,
        out_dir,
        use_real_time=use_real_time,
        dt=dt,
        show_plot=args.show,
        sg_window=args.sg_window,
        sg_polyorder=args.sg_polyorder,
        tail_points=args.tail_points,
    )

    print("\n===== 全部完成 =====")
    print(f"共处理 Tstar 数量: {len(all_data)}")
    for tstar in sorted(all_data.keys()):
        print(f"{tstar}: {len(all_data[tstar])} frames")
