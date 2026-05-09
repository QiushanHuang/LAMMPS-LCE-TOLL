# Radius of Gyration Analysis

This directory contains Python scripts for mass-weighted radius-of-gyration (`Rg`) analysis from LAMMPS dump files.

## Recommended Script

Use:

- `compute_rg_all_tstar_unwrap_v11.py`

This version:

- prefers unwrapped coordinates `xu/yu/zu`, and falls back to `x/y/z` if needed
- can be run from the directory you want to analyze
- supports two directory modes:
  - `Tstar_*` grouped recursive search under the target root
  - direct processing of `*.dump` files in the current directory if no `Tstar_*` subdirectory exists
- supports multi-process acceleration with `--jobs`
- keeps the current plotting layout:
  - Figure 1: raw `Rg` curves
  - Figure 2: raw `Rg` curves with an SG-smoothed overlay
  - Figure 3: tail-average `Rg` vs `T*` with error bars
  - combined figure: three stacked views in one PNG and one matplotlib window

## Plot Outputs

The script writes:

- `figure1_rg_vs_timestep_all_Tstar.png` or `figure1_rg_vs_time_all_Tstar.png`
- `figure2_rg_vs_timestep_all_Tstar_sg.png` or `figure2_rg_vs_time_all_Tstar_sg.png`
- `figure3_rg_tail_mean_vs_temperature.png`
- `figure_all_rg_views.png`

## Usage

Run from the dump directory:

```bash
cd /path/to/dump-directory
python /path/to/LAMMPS-LCE-TOLL/Figure-Analyze/Radius-of-Gyration/compute_rg_all_tstar_unwrap_v11.py --show
```

Or specify an explicit root:

```bash
python /path/to/LAMMPS-LCE-TOLL/Figure-Analyze/Radius-of-Gyration/compute_rg_all_tstar_unwrap_v11.py \
  --root-dir /path/to/data \
  --out-dir Rg_results \
  --show
```

Use multiple CPU cores explicitly:

```bash
python /path/to/LAMMPS-LCE-TOLL/Figure-Analyze/Radius-of-Gyration/compute_rg_all_tstar_unwrap_v11.py \
  --root-dir /path/to/data \
  --jobs 8 \
  --show
```

## Main Options

- `--root-dir`: root directory to analyze, default `.`.
- `--out-dir`: output directory, default `Rg_results`.
- `--types`: atom types included in the `Rg` calculation, default `1 2`.
- `--use-real-time`: use physical time instead of timestep on the x-axis.
- `--dt`: timestep-to-time conversion factor when `--use-real-time` is enabled.
- `--show`: open the combined matplotlib window.
- `--sg-window`: Savitzky-Golay window length, default `2001`.
- `--sg-polyorder`: Savitzky-Golay polynomial order, default `2`.
- `--tail-points`: number of tail points used in Figure 3, default `3000`.
- `--jobs`: number of worker processes. Default is automatic. Use `--jobs 1` to disable parallel execution.

## Notes

- Run `v11` as a normal Python script. This is the intended way to use the multi-process mode.
- On large datasets, start with automatic parallelism. If the machine becomes memory-bound or the environment has trouble spawning worker processes, retry with a smaller value such as `--jobs 4`, or disable parallelism with `--jobs 1`.
- If you `cd` into a directory that already contains dump files and there are no `Tstar_*` subdirectories below it, the script switches to direct-current-directory mode automatically.
- Relative `--out-dir` paths are resolved relative to `--root-dir`.
- Figure 3 uses the last `3000` valid points of each `Tstar` by default.
- Figure 1 and Figure 2 keep the end-of-line `Tstar` text annotations; the old per-curve legend box is removed.

## Data Outputs

The script also writes:

- `Rg_all_Tstar.dat`
- `<Tstar>_Rg.dat`

Each row records timestep, `Rg`, `Rg^2`, particle count, total mass, coordinate source, and source file path.
