# Radius of Gyration Analysis

This directory contains Python scripts for mass-weighted radius-of-gyration (`Rg`) analysis from LAMMPS dump files.

## Recommended Script

Use:

- `compute_rg_all_tstar_unwrap_v4.py`

This version:

- prefers unwrapped coordinates `xu/yu/zu`, and falls back to `x/y/z` if needed
- can be run from the directory you want to analyze
- supports two directory modes:
  - `Tstar_*` grouped recursive search under the target root
  - direct processing of `*.dump` files in the current directory if no `Tstar_*` subdirectory exists

## Plot Outputs

The script writes:

- `figure1_rg_vs_timestep_all_Tstar.png` or `figure1_rg_vs_time_all_Tstar.png`
- `figure2_rg_vs_timestep_all_Tstar_sg.png` or `figure2_rg_vs_time_all_Tstar_sg.png`
- `figure3_rg_tail_mean_vs_temperature.png`
- `figure_all_rg_views.png`

Figure behavior:

- Figure 1: original `Rg` curves
- Figure 2: original `Rg` curves with SG-smoothed curves overlaid
- Figure 3: mean `Rg` from the last `N` points for each `Tstar`, with error bars from the tail-window standard deviation
- Combined figure: all three views stacked in one matplotlib window and one PNG

## Usage

Run from the dump directory:

```bash
cd /path/to/dump-directory
python /path/to/LAMMPS-LCE-TOLL/Figure-Analyze/Radius-of-Gyration/compute_rg_all_tstar_unwrap_v4.py --show
```

Or specify an explicit root:

```bash
python /path/to/LAMMPS-LCE-TOLL/Figure-Analyze/Radius-of-Gyration/compute_rg_all_tstar_unwrap_v4.py \
  --root-dir /path/to/data \
  --out-dir Rg_results \
  --show
```

## Main Options

- `--root-dir`: root directory to analyze, default `.`.
- `--out-dir`: output directory, default `Rg_results`.
- `--types`: atom types included in the `Rg` calculation, default `1 2`.
- `--use-real-time`: use physical time instead of timestep on the x-axis.
- `--dt`: timestep-to-time conversion factor when `--use-real-time` is enabled.
- `--show`: open the combined matplotlib window.
- `--sg-window`: Savitzky-Golay window length, default `51`.
- `--sg-polyorder`: Savitzky-Golay polynomial order, default `3`.
- `--tail-points`: number of tail points used in Figure 3, default `1000`.

## Data Outputs

The script also writes:

- `Rg_all_Tstar.dat`
- `<Tstar>_Rg.dat`

Each row records timestep, `Rg`, `Rg^2`, particle count, total mass, coordinate source, and source file path.
