# LAMMPS-LCE-TOLL

This repository collects LAMMPS workflows and post-processing utilities for LCE-related simulations.

## Directory Layout

- `Relaxation-and-Stretching/`: restart-based relaxation and constant-force stretching inputs, plus helper scripts
- `Figure-Analyze/Radius-of-Gyration/`: radius-of-gyration analysis scripts for dump trajectories

## Current Analysis Utility

The latest radius-of-gyration analysis entry point is:

- `Figure-Analyze/Radius-of-Gyration/compute_rg_all_tstar_unwrap_v11.py`

That script:

- prefers `xu/yu/zu` and falls back to `x/y/z`
- can be launched from the directory that contains the target dumps
- supports both `Tstar_*` grouped layouts and direct-current-directory dump processing
- supports multi-process acceleration through `--jobs`
- exports raw and smoothed `Rg` plots plus a tail-averaged `Rg` vs `T*` figure
- draws Figure 2 as a raw-curve background with a visually distinct SG overlay
- uses `--sg-window 2001`, `--sg-polyorder 2`, and `--tail-points 3000` by default

See `Figure-Analyze/Radius-of-Gyration/README.md` for usage details.
