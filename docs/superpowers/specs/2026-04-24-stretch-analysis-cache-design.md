# Stretch Analysis Cache Design

**Date:** 2026-04-24

**Goal**

Refactor the stretching analysis workflow so dump-derived stretch datasets are cached in one central directory per run root, reused across repeated analysis runs, and used as the sole data source for downstream comparison plots.

**Scope**

- Applies to the `rg_T_*` single-directory analyzer, batch analyzer, and batch comparison workflow under `/Users/joshua/Desktop/MD/LAMMPS-LCE-TOLL/Relaxation-and-Stretching/analysis`.
- Uses a central cache directory named `0_stretch_analysis` created directly under the analysis root that currently contains `output_*` directories.
- Introduces a sibling ignore directory `0_stretch_analysis/wait-list-file` that is never scanned for datasets or plots.

**Directory Layout**

Given a root such as:

- `/path/to/run-root/output_*`

the workflow creates:

- `/path/to/run-root/0_stretch_analysis/`
- `/path/to/run-root/0_stretch_analysis/wait-list-file/`

The cache directory stores:

- one processed dataset per `output_*` directory
- batch-level comparison figures generated from those cached datasets

The batch scanners must only read cached datasets from `0_stretch_analysis` itself. They must not recurse into `wait-list-file` or re-scan raw `output_*` directories once datasets are available.

**Dataset Naming**

For each `output_*` directory, generate:

- `stretch_dataset_<output-tag>.dat`

where `<output-tag>` is the original directory name with the leading `output_` removed.

Example:

- `output_restart16000000_Tstar_0.7_F0.12_tail_v0_no_thermal_xuyu_zu`
- `stretch_dataset_restart16000000_Tstar_0.7_F0.12_tail_v0_no_thermal_xuyu_zu.dat`

This naming is stable and directly reversible back to the source run tag.

**Dataset Format**

The dataset is a plain-text table similar in style to `Rg_all_Tstar.dat`: one header line followed by one row per sampled timestep.

Required columns:

- `source_tag`
- `timestep`
- `time`
- `L`
- `dL`
- `strain`
- `contour_fraction`
- `Fext`
- `Fchain`
- `Ftotal`
- `Rg`
- `Rg2`
- `count`
- `total_mass`
- `coord_source`
- `dump_file`

The dataset must contain enough information to regenerate:

- `Length and Contour Fraction vs Time`
- `Rg vs Time`
- any smoothed variants of those figures

without reopening raw dump files.

**Processing Rules**

1. Discover root-level `output_*` directories.
2. Map each `output_*` to its target cached dataset inside `0_stretch_analysis`.
3. If the cached dataset already exists, skip dump processing for that run.
4. If the cached dataset is missing, build it from:
   - `force_clamp_response.dat`
   - `traj.force_clamp_aligned.*.dump`
5. After all missing datasets are filled, generate plots from cached datasets only.

This guarantees resumable runs and avoids recomputing dump-derived `Rg` when only plots need updating.

**Parallelism**

- Single-directory dataset construction retains dump-file parallel parsing via `--jobs`.
- Batch processing may parallelize across `output_*` directories via `--jobs`.
- Batch plotting reads cached datasets, so it should be materially cheaper than raw dump parsing.
- `--show` still forces serial plot display order at the batch level.

**Plot Outputs**

The batch comparison workflow writes comparison figures into `0_stretch_analysis` rather than the raw run root.

The existing comparison layout remains:

- one total figure
- top-left raw `Length and Contour Fraction vs Time`
- top-right smoothed `Length and Contour Fraction vs Time`
- bottom-left raw `Rg vs Time`
- bottom-right smoothed `Rg vs Time`

The default smoothing method remains Savitzky-Golay.

**Non-Goals**

- No recursive scan of nested directories below `0_stretch_analysis`
- No use of `wait-list-file` as an automatic input source
- No overwrite of existing cached datasets unless an explicit future refresh flag is added

**Error Handling**

- If `force_clamp_response.dat` is missing for an `output_*`, skip that run and report it.
- If no valid dump frames exist for a run, fail that dataset build and report it without deleting any already-built datasets.
- If plot generation encounters no cached datasets, fail with a clear error naming `0_stretch_analysis`.

**Implementation Impact**

- Update the `rg_T_*` analyzer to support writing and reading the new cached dataset format.
- Update batch analysis to populate `0_stretch_analysis` first, then plot from cache.
- Update batch comparison to scan only `0_stretch_analysis/*.dat` and ignore `wait-list-file`.
