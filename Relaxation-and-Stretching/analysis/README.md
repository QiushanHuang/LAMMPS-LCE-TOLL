# Relaxation-and-Stretching Analysis

This directory is the packaged analysis tool group for the `Relaxation-and-Stretching` workflow.

Contributor:

- Qiushan Huang (GitHub: `QiushanHuang`)

It contains two analysis families:

- legacy `analyze_force_clamp_aligned_box*.py`
- new `rg_T_*` and batch comparison tools

The new `rg_T_*` tools are the ones to use when you want:

- cached dump-derived `Rg` data
- batch processing over many `output_*` folders
- comparison of chain `Rg` versus applied force `F`
- two sequential `--show` canvases for trajectory view and force-summary view

---

## File Map

### Legacy single-output analyzers

- `analyze_force_clamp_aligned_box.py`
  - original single-directory analyzer
- `analyze_force_clamp_aligned_box_loess.py`
  - legacy LOESS wrapper
- `analyze_force_clamp_aligned_box_savgol.py`
  - legacy Savitzky-Golay wrapper

### New single-output and batch analyzers

- `rg_T_analyze_force_clamp_aligned_box.py`
  - new single-output analyzer
  - builds or reuses a central cached dataset
  - computes dump-derived `Rg`, `Rg_x`, and `Lx`
- `rg_T_analyze_force_clamp_aligned_box_loess.py`
  - new LOESS wrapper for the `rg_T` single-output analyzer
- `rg_T_analyze_force_clamp_aligned_box_savgol.py`
  - new Savitzky-Golay wrapper for the `rg_T` single-output analyzer
- `batch_rg_T_analyze.py`
  - runs the new `rg_T` single-output analyzer over every `output_*` folder under one root
- `compare_batch_rg_T.py`
  - builds/reuses all cached datasets under one root
  - generates the batch comparison figures
  - this is the main script for analyzing chain `Rg` and `F` relationships

---

## Central Cache Layout

When you analyze a root directory that contains many `output_*` folders, the new tools use one sibling cache directory:

```text
<root>/
├── output_*/
└── 0_stretch_analysis/
    ├── stretch_dataset_*.dat
    ├── rg_T_batch_compare.png
    ├── rg_T_batch_compare_force_window.png
    └── wait-list-file/
```

Rules:

- every `output_*` folder maps to one cached dataset:
  - `stretch_dataset_<output-tag>.dat`
- cached datasets are stored only in `0_stretch_analysis/`
- `wait-list-file/` is ignored by scans
- if you move one cached dataset into `wait-list-file/`, later batch runs skip that sample
- if an old cached dataset is missing the new `Lx` / `Rg_x` columns, the script deletes and rebuilds it automatically

---

## Cached Dataset Content

Each `stretch_dataset_*.dat` is a reusable text table with one row per cached sampled timestep.

Key columns:

- `time`
- `L`
- `contour_fraction`
- `Fext`
- `Fchain`
- `Ftotal`
- `Lx`
- `Rg`
- `Rg_x`

Definitions:

- `L`
  - projected chain length from `force_clamp_response.dat`
- `contour_fraction`
  - projected length divided by contour length
- `Lx`
  - x-span from dump data: `max(x) - min(x)` over type 1 and type 2 particles
- `Rg`
  - mass-weighted radius of gyration from dump data
- `Rg_x`
  - x-projected mass-weighted radius of gyration from dump data

Dump coordinate rule:

- use `xu/yu/zu` if present
- otherwise fall back to `x/y/z`

---

## Single-Output Analysis

If you are already inside one `output_*` directory, run:

```bash
cd /path/to/output_...
python /Users/joshua/Desktop/MD/LAMMPS-LCE-TOLL/Relaxation-and-Stretching/analysis/rg_T_analyze_force_clamp_aligned_box.py --jobs 4
```

This does two things:

- creates or reuses the cached dataset in the sibling `0_stretch_analysis/`
- writes per-output figures into `analysis_rg_T/`

To open the interactive figures:

```bash
cd /path/to/output_...
python /Users/joshua/Desktop/MD/LAMMPS-LCE-TOLL/Relaxation-and-Stretching/analysis/rg_T_analyze_force_clamp_aligned_box.py --jobs 4 --show
```

Shortcuts:

```bash
python /Users/joshua/Desktop/MD/LAMMPS-LCE-TOLL/Relaxation-and-Stretching/analysis/rg_T_analyze_force_clamp_aligned_box_loess.py --jobs 4
python /Users/joshua/Desktop/MD/LAMMPS-LCE-TOLL/Relaxation-and-Stretching/analysis/rg_T_analyze_force_clamp_aligned_box_savgol.py --jobs 4
```

---

## Batch Per-Output Analysis

If a root directory contains many `output_*` folders, run:

```bash
python /Users/joshua/Desktop/MD/LAMMPS-LCE-TOLL/Relaxation-and-Stretching/analysis/batch_rg_T_analyze.py \
  --root /path/to/root \
  --jobs 4 \
  --analysis-jobs 1
```

Meaning:

- `--jobs`
  - parallelizes across `output_*` directories
- `--analysis-jobs`
  - parallelizes dump parsing inside each single-output analysis

If you also want interactive windows, use:

```bash
python /Users/joshua/Desktop/MD/LAMMPS-LCE-TOLL/Relaxation-and-Stretching/analysis/batch_rg_T_analyze.py \
  --root /path/to/root \
  --show \
  --analysis-jobs 4
```

With `--show`, directories are handled serially so windows open one after another.

---

## Analyze Chain `Rg` vs Force `F`

This is the main workflow for comparing chain response across force.

Run:

```bash
python /Users/joshua/Desktop/MD/LAMMPS-LCE-TOLL/Relaxation-and-Stretching/analysis/compare_batch_rg_T.py \
  --root /path/to/root \
  --jobs 4 \
  --analysis-jobs 1
```

This script:

1. scans the root for `output_*`
2. creates missing `stretch_dataset_*.dat` files in `0_stretch_analysis/`
3. reuses existing valid cached datasets
4. ignores anything placed in `0_stretch_analysis/wait-list-file/`
5. generates two comparison canvases

### What the two canvases mean

#### Canvas 1: trajectory comparison

Saved as:

- `0_stretch_analysis/rg_T_batch_compare.png`

It contains 4 panels:

- `Length and Contour Fraction vs Time`
- `Lx vs Time`
- `Rg vs Time`
- `Rg_x vs Time`

#### Canvas 2: force-summary comparison

Saved as:

- `0_stretch_analysis/rg_T_batch_compare_force_window.png`

It contains 4 panels:

- mean `L` vs `F` with error bars
- mean `Lx` vs `F` with error bars
- mean `Rg` vs `F` with error bars
- mean `Rg_x` vs `F` with error bars

The x-axis force `F` is parsed from the cached dataset name, for example:

- `stretch_dataset_restart860000_Tstar_0.50_F20_tail_v0_no_thermal_xuyu_zu.dat`

### What the error bars mean

The second canvas uses the final tail window of each cached dataset.

Default:

- last `3000` cached samples

Error bars:

- population standard deviation over that tail window

To change the window size:

```bash
python /Users/joshua/Desktop/MD/LAMMPS-LCE-TOLL/Relaxation-and-Stretching/analysis/compare_batch_rg_T.py \
  --root /path/to/root \
  --window-points 1000
```

### Interactive display

To show both canvases:

```bash
python /Users/joshua/Desktop/MD/LAMMPS-LCE-TOLL/Relaxation-and-Stretching/analysis/compare_batch_rg_T.py \
  --root /path/to/root \
  --show
```

Behavior:

- first window opens for the trajectory canvas
- after you close it, the second window opens for the force-summary canvas

---

## Recommended Command for Your `Rg`-`F` Question

If your root directory is:

```text
/Users/joshua/Desktop/MD/2026_04/0408/02_Youngs_Modulus/length7-Tstar0.7
```

run:

```bash
python /Users/joshua/Desktop/MD/LAMMPS-LCE-TOLL/Relaxation-and-Stretching/analysis/compare_batch_rg_T.py \
  --root /Users/joshua/Desktop/MD/2026_04/0408/02_Youngs_Modulus/length7-Tstar0.7 \
  --jobs 4 \
  --analysis-jobs 1 \
  --show
```

Use this result as follows:

- if you want the full time evolution of chain compaction/extension:
  - read Canvas 1
- if you want the final-state relationship between chain size and applied force:
  - read Canvas 2
- if you want just the `Rg`-`F` relationship:
  - read the lower-left panel of Canvas 2
- if you want the x-projected `Rg_x`-`F` relationship:
  - read the lower-right panel of Canvas 2

---

## Practical Notes

- batch comparison reads only top-level cached datasets in `0_stretch_analysis/`
- anything under `wait-list-file/` is skipped
- if you add new `output_*` directories later, rerun the same command; only missing datasets are built
- if you regenerate one trajectory and want to force a fresh cache, delete its corresponding `stretch_dataset_*.dat` from `0_stretch_analysis/`
