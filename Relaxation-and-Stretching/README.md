# Relaxation and Stretching Workflow

This package collects restart-based LAMMPS input scripts for two stages:

- `relaxation/`: rebuild the chain geometry from restart `xu/yu/zu`, expand the box, rotate/center the cluster, then run free thermal relaxation in either `fff` or `ppp`
- `stretching/`: take a relaxed restart, align the chain into a slender box, pin the head, and apply a constant force to the tail

The scripts share a helper in `helpers/build_output_root.py` that builds the output directory name from the restart file, target temperature, force magnitude, and protocol tag.

## Directory Layout

File names now encode the main distinguishing parameters directly:

- relaxation: boundary + `T*` + cubic `box_length`
- stretching: force + `T*` + `x_head_margin` + protocol variant

- `analysis/analyze_force_clamp_aligned_box.py`: main analysis entry point for `force_clamp_response.dat`
- `analysis/analyze_force_clamp_aligned_box_loess.py`: wrapper around the main analyzer with default `trend_method=loess` and output directory `analysis_loess`
- `analysis/analyze_force_clamp_aligned_box_savgol.py`: wrapper around the main analyzer with default `trend_method=savgol`, output directory `analysis_savgol`, and default `sg_window=11`
- `helpers/build_output_root.py`: generates `output_root` names
- `relaxation/in.restart_relaxation_fff_T050_box200_xuyu-zu.lmp`: baseline `fff`, `T*=0.50`, cubic `box_length=200`
- `relaxation/in.restart_relaxation_fff_T080_box400_xuyu-zu.lmp`: `fff`, `T*=0.80`, cubic `box_length=400`
- `relaxation/in.restart_relaxation_ppp_T040_box200_xuyu-zu.lmp` ... `relaxation/in.restart_relaxation_ppp_T100_box200_xuyu-zu.lmp`: `ppp` temperature grid at `T*=0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00`, all with cubic `box_length=200`
- `stretching/in.single_chain_constant_force_aligned_box-F5.00-T050-margin50-xuyu-zu.lmp`: standard `F=5`, `T*=0.50`, `x_head_margin=50`
- `stretching/in.single_chain_constant_force_aligned_box-F5.00-T050-margin80-xuyu-zu.lmp`: standard `F=5`, `T*=0.50`, `x_head_margin=80`
- `stretching/in.single_chain_constant_force_aligned_box-F10.00-T050-margin80-xuyu-zu.lmp`: standard `F=10`, `T*=0.50`, `x_head_margin=80`
- `stretching/in.single_chain_constant_force_aligned_box-F10.00-T080-margin60-tail_v0_no_thermal-xuyu-zu.lmp`: force-clamp stretching, `F=10`, `T*=0.80`, `x_head_margin=60`
- `stretching/in.single_chain_constant_force_aligned_box-F20.00-T050-margin150-tail_v0_no_thermal-xuyu-zu.lmp`: force-clamp stretching, `F=20`, `T*=0.50`, `x_head_margin=150`, tail velocity zeroed before loading
- `stretching/in.single_chain_constant_force_aligned_box-F20.00-T040-margin100-xuyu-zu.lmp` ... `stretching/in.single_chain_constant_force_aligned_box-F20.00-T100-margin100-xuyu-zu.lmp`: standard `F=20` temperature grid at `T*=0.40, 0.45, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00`, all with `x_head_margin=100`
- `stretching/run_all_lammps_tmux_linux.py`: Linux/tmux batch launcher for the packaged stretching inputs; it discovers `.lmp` files in `stretching/`, but runs them from your current working directory so restart paths and output directories stay local to the run directory
- `stretching/short_benchmark_lammps.py`: short benchmark runner for one specified stretching `.lmp`; it sweeps an MPI × OMP grid, overrides runtime/output intervals via `-var`, and writes a CSV summary under a dedicated benchmark directory

## Endpoint Selection

All packaged scripts now support automatic endpoint detection:

- `head_id_override index auto`
- `tail_id_override index auto`

Default behavior:

- the smallest atom ID is used as `head_id`
- the largest atom ID is used as `tail_id`

If your restart does **not** place the physical chain ends at the minimum and maximum atom IDs, override them manually by editing the script or from the command line, for example:

```bash
lmp_mpi -var head_id_override 1 -var tail_id_override 8006 -in in.single_chain_constant_force_aligned_box-F20.00-T050-margin150-tail_v0_no_thermal-xuyu-zu.lmp
```

Use the **same** endpoint override pair in both relaxation and stretching when you want a consistent workflow.

## How to Use

### 1. Run relaxation

Change into the chosen relaxation directory and run one of the scripts:

```bash
cd relaxation
mpirun -np 4 lmp_mpi -in in.restart_relaxation_ppp_T050_box200_xuyu-zu.lmp
```

What it does:

- reads `restart_file`
- snapshots the restart's `xu/yu/zu`
- expands to a temporary large `fff` box
- restores coordinates from the restart snapshot
- rotates the chain using the head→tail direction
- centers the cluster
- shrinks to the final cubic box
- runs free relaxation

For the `ppp` variants, only the **final** relaxation box is periodic; the transform stage remains nonperiodic to avoid premature wrapping during rotation and translation. All packaged `ppp` variants in this directory now use `box_length=200`.

### 2. Feed the relaxed restart into stretching

Choose a relaxation restart or `Final.relaxation.bin`, then update the stretching script:

- set `restart_file`
- keep `target_T` consistent if that is part of your protocol
- keep `head_id_override` and `tail_id_override` consistent
- packaged stretching defaults currently use `restart_file = Restart.relaxation.860000` and `run_time = 5000`

Then run, for example:

```bash
cd ../stretching
mpirun -np 4 lmp_mpi -in in.single_chain_constant_force_aligned_box-F20.00-T080-margin100-xuyu-zu.lmp
```

## Recommended Command Templates

The most robust workflow is:

1. keep the `.lmp` scripts in this repository
2. `cd` into the working directory that contains your restart files
3. launch LAMMPS with the **absolute path** to the script in this repository

In this mode:

- `restart_file` is resolved relative to your current working directory unless you override it with `-var restart_file ...`
- `output_base` defaults to `.`, so outputs are written into your current working directory
- the packaged `helper_root` already points to this local repository clone

### Relaxation Template

If `Restart.GB.rho030.500000` is in your current working directory:

```bash
cd /path/to/working-directory
/opt/homebrew/bin/mpirun -np 4 /opt/homebrew/bin/lmp_mpi \
  -in /Users/joshua/Desktop/MD/LAMMPS-LCE-TOLL/Relaxation-and-Stretching/relaxation/in.restart_relaxation_ppp_T050_box200_xuyu-zu.lmp
```

If the restart name is different, override it explicitly:

```bash
cd /path/to/working-directory
/opt/homebrew/bin/mpirun -np 4 /opt/homebrew/bin/lmp_mpi \
  -var restart_file Restart.GB.rho030.900000 \
  -in /Users/joshua/Desktop/MD/LAMMPS-LCE-TOLL/Relaxation-and-Stretching/relaxation/in.restart_relaxation_ppp_T050_box200_xuyu-zu.lmp
```

### Stretching Template

If `Restart.relaxation.860000` is in your current working directory:

```bash
cd /path/to/working-directory
/opt/homebrew/bin/mpirun -np 4 /opt/homebrew/bin/lmp_mpi \
  -in /Users/joshua/Desktop/MD/LAMMPS-LCE-TOLL/Relaxation-and-Stretching/stretching/in.single_chain_constant_force_aligned_box-F20.00-T080-margin100-xuyu-zu.lmp
```

If you want to select a different relaxed restart:

```bash
cd /path/to/working-directory
/opt/homebrew/bin/mpirun -np 4 /opt/homebrew/bin/lmp_mpi \
  -var restart_file Restart.relaxation.920000 \
  -in /Users/joshua/Desktop/MD/LAMMPS-LCE-TOLL/Relaxation-and-Stretching/stretching/in.single_chain_constant_force_aligned_box-F20.00-T080-margin100-xuyu-zu.lmp
```

### Batch Launch Template

If you want one detached tmux session per packaged stretching script, stay in the working directory that contains your relaxed restart and run:

```bash
cd /path/to/working-directory
python /Users/joshua/Desktop/MD/LAMMPS-LCE-TOLL/Relaxation-and-Stretching/stretching/run_all_lammps_tmux_linux.py
```

What the launcher does:

- discovers `.lmp` files in `/Users/joshua/Desktop/MD/LAMMPS-LCE-TOLL/Relaxation-and-Stretching/stretching`
- keeps your current directory as the runtime working directory
- computes the same `output_root` that the LAMMPS script will generate
- creates `launch_in_tmux.sh`, `command.txt`, and `tmux.log` inside each `output_root`
- starts one detached tmux session per job

Useful options:

```bash
python /Users/joshua/Desktop/MD/LAMMPS-LCE-TOLL/Relaxation-and-Stretching/stretching/run_all_lammps_tmux_linux.py \
  --dry-run \
  --mpi-ranks 4 \
  --omp-threads 2 \
  --restart-file Restart.relaxation.920000
```

- `--dry-run`: print the launch plan without requiring `tmux`, `mpiexec`, or `lmp`
- `--restart-file`: override `restart_file` for every discovered stretching script
- `--workdir`: run against a directory other than the current one
- `--pattern`: restrict discovery to a subset such as `in.single_chain_constant_force_aligned_box-F20*.lmp`
- `--rerun-existing-output`: allow launching a job even if its `output_root` already exists
- `--reuse-session-names`: append a numeric suffix instead of skipping when the default tmux session name already exists

Important constraint:

- the launcher rejects two scripts that resolve to the same `output_root`
- in practice this happens when two scripts differ in filename but share the same `restart_file`, `target_T`, `force_mag`, and `protocol_tag`
- if you want both scripts in the same batch, give them distinct `protocol_tag` values or launch them separately

### Short Benchmark Template

If you want to benchmark one specific stretching input with a very short run, stay in the working directory that contains the restart file and run:

```bash
cd /path/to/working-directory
python /Users/joshua/Desktop/MD/LAMMPS-LCE-TOLL/Relaxation-and-Stretching/stretching/short_benchmark_lammps.py \
  --input /Users/joshua/Desktop/MD/LAMMPS-LCE-TOLL/Relaxation-and-Stretching/stretching/in.single_chain_constant_force_aligned_box-F20.00-T080-margin100-xuyu-zu.lmp
```

Default benchmark behavior:

- benchmarks exactly one `.lmp`
- sweeps `MPI ranks = 1..6` and `OMP threads = 1..6`
- overrides `run_time` to a short value (`0.2` by default)
- overrides `sample_every`, `dump_every`, and `restart_every` to very large values so I/O does not dominate the timing
- writes results to `short_benchmark_<input-stem>/`
- writes one subdirectory per `(mpi, omp)` pair plus `benchmark_results.csv`

Useful options:

```bash
python /Users/joshua/Desktop/MD/LAMMPS-LCE-TOLL/Relaxation-and-Stretching/stretching/short_benchmark_lammps.py \
  --input /path/to/in.single_chain_constant_force_aligned_box-F20.00-T080-margin100-xuyu-zu.lmp \
  --dry-run \
  --mpi-min 1 --mpi-max 6 \
  --omp-min 1 --omp-max 6 \
  --run-time 0.05 \
  --restart-file Restart.relaxation.860000
```

- `--dry-run`: print all benchmark commands without running LAMMPS
- `--run-time`: set the short test duration in LJ time units
- `--mpi-min/--mpi-max`, `--omp-min/--omp-max`: define the benchmark matrix
- `--restart-file`: override the restart used by the selected `.lmp`
- `--sample-every`, `--dump-every`, `--restart-every`: adjust output suppression if you want a more realistic but noisier timing
- `--mpiexec`, `--lammps-bin`: override the server defaults when testing on a different machine

### Server Deployment Layout

If the goal is "upload a small folder set to Linux and then batch-launch from the run directory", the minimum layout is:

```text
Relaxation-and-Stretching/
├── helpers/
│   └── build_output_root.py
└── stretching/
    ├── run_all_lammps_tmux_linux.py
    └── in.single_chain_constant_force_aligned_box-*.lmp
```

If you also want packaged relaxation and analysis on the server, upload the whole `Relaxation-and-Stretching/` directory. For stretching-only batch jobs, `helpers/` + `stretching/` is sufficient.

Recommended server workflow:

1. copy `Relaxation-and-Stretching/helpers` and `Relaxation-and-Stretching/stretching` to the server without changing their relative positions
2. `cd` into the directory that contains `Restart.relaxation.860000` or your chosen restart file
3. run:

```bash
python /path/to/Relaxation-and-Stretching/stretching/run_all_lammps_tmux_linux.py --dry-run
python /path/to/Relaxation-and-Stretching/stretching/run_all_lammps_tmux_linux.py
```

The batch launcher now always passes these values explicitly to LAMMPS:

- `-var helper_root /path/to/Relaxation-and-Stretching/helpers`
- `-var output_base <your current working directory>`

That is the reason this script is more reliable than calling the `.lmp` files by hand on a new machine: it removes the ambiguity around local absolute paths and current working directory.

### Why Helper Paths Seemed Confusing

There are two different path-resolution rules in play:

- Python scripts usually resolve helper paths relative to `__file__`
- LAMMPS `shell` / `include` commands resolve paths relative to the LAMMPS process working directory unless you pass explicit absolute paths or `-var` overrides

So if a `.lmp` file contains a local Mac path such as `/Users/.../helpers`, or a relative helper path that only worked from one directory, moving the script to another machine breaks it immediately.

The packaged batch launcher avoids that by discovering the helper folder from its own location and then injecting it into each LAMMPS job with `-var helper_root ...`. On Linux the default LAMMPS executable in the launcher is:

```text
/home/hkust/mylammps/build/lmp
```

Override it only if the server uses a different binary:

```bash
python /path/to/Relaxation-and-Stretching/stretching/run_all_lammps_tmux_linux.py \
  --lammps-bin /some/other/lmp
```

## Analysis Scripts

For the packaged analysis tool group, including the new central-cache `rg_T_*` workflow, batch processing, and `Rg`-vs-`F` comparison commands, see:

- `Relaxation-and-Stretching/analysis/README.md`

### Script Relationship

- `helpers/build_output_root.py` is used by the LAMMPS input scripts to construct `output_root`
- `analysis/analyze_force_clamp_aligned_box.py` is the main analysis program
- `analysis/analyze_force_clamp_aligned_box_loess.py` imports the main program and runs it with LOESS defaults
- `analysis/analyze_force_clamp_aligned_box_savgol.py` imports the main program and runs it with Savitzky-Golay defaults

### Default Current-Directory Workflow

If you `cd` into a simulation output directory that contains `force_clamp_response.dat`, the analysis scripts now default to that file automatically. You do **not** need to pass the data path explicitly.

Main analyzer:

```bash
cd /path/to/output-directory
python /Users/joshua/Desktop/MD/LAMMPS-LCE-TOLL/Relaxation-and-Stretching/analysis/analyze_force_clamp_aligned_box.py
```

LOESS shortcut:

```bash
cd /path/to/output-directory
python /Users/joshua/Desktop/MD/LAMMPS-LCE-TOLL/Relaxation-and-Stretching/analysis/analyze_force_clamp_aligned_box_loess.py
```

Savitzky-Golay shortcut:

```bash
cd /path/to/output-directory
python /Users/joshua/Desktop/MD/LAMMPS-LCE-TOLL/Relaxation-and-Stretching/analysis/analyze_force_clamp_aligned_box_savgol.py
```

If `force_clamp_response.dat` is not in the current working directory, pass it explicitly:

```bash
python /Users/joshua/Desktop/MD/LAMMPS-LCE-TOLL/Relaxation-and-Stretching/analysis/analyze_force_clamp_aligned_box.py \
  /path/to/force_clamp_response.dat
```

### Analysis Outputs

- main analyzer default output directory: `analysis/`
- LOESS wrapper default output directory: `analysis_loess/`
- Savitzky-Golay wrapper default output directory: `analysis_savgol/`
- all three write PNG figures plus `summary.json`
- the analyzer tolerates `force_clamp_response.dat` files that contain NUL-padded gaps from interrupted or resumed runs; pure NUL blocks are ignored during parsing

### Important Analysis Parameters

- `--trend-method loess|savgol`
  - choose the smoothing backend explicitly when using the main analyzer
- `--trend-frac`
  - LOESS neighborhood fraction
  - default is `0.08`
  - larger values smooth more aggressively and suppress local fluctuations
  - smaller values preserve local structure but track noise more easily
- `--trend-max-points`
  - maximum number of points used internally by LOESS before downsampling + interpolation
  - default is `2000`
- `--sg-window`
  - controls Savitzky-Golay smoothing strength
  - wrapper default is `11`
  - larger windows smooth more strongly but can wash out short-lived force or velocity features
  - smaller windows preserve sharper transients but are more sensitive to noise
  - in practice: start with `11`, try `21` or `31` only for longer and noisier trajectories
- `--sg-order`
  - local polynomial order for the Savitzky-Golay fit
  - default is `3`
  - keep this low unless you have a specific reason to fit higher-order local curvature
- `--output-dir`
  - override the default analysis output directory if you want results somewhere else
- `--show`
  - opens an interactive matplotlib window if matplotlib is available in the current Python environment

## Input / Output Naming

Output directories are generated as:

```text
output_restart<restart-tag>_Tstar_<target_T>_F<force_mag>_<protocol_tag>
```

Examples:

- relaxation: `output_restart500000_Tstar_0.50_F20_relaxation_ppp_xuyu_zu`
- stretching: `output_restart860000_Tstar_0.5_F20_tail_v0_no_thermal_xuyu_zu`

Relaxation outputs:

- `traj.relaxation.*.dump`
- `Restart.relaxation.*`
- `Final.relaxation.bin`
- `relaxation_state.dat`

Stretching outputs:

- `traj.force_clamp_aligned.*.dump`
- `Restart.force_clamp_aligned.*`
- `Final.force_clamp_aligned.bin`
- `force_clamp_response.dat`

## Parameters You Should Review Before Every Run

### Shared / consistency-sensitive

- `restart_file`: must point to the correct upstream restart
- `target_T`: keep consistent if you want the same thermal state across stages
- `head_id_override`, `tail_id_override`: keep consistent across stages
- `ts`, `Tdamp`, `Tchain`: changing these changes the thermostat/integration protocol
- `sample_every`, `dump_every`, `restart_every`: review together when changing output density

### Relaxation-specific

- `box_length`: final cubic box size
- `protocol_tag`: affects output naming
- `relax_time`: total relaxation duration
- `fff` vs `ppp`: choose based on the ensemble you want after centering

### Stretching-specific

- `force_mag`: applied load
- `run_time`: total stretching duration
- `x_head_margin`: left margin for placing the head before loading
- `x_length_factor`: final slender-box length scaling
- `dry_run`: use `yes` to build geometry and outputs without time integration

## Important Notes

- These scripts dump `xu yu zu`, but frame 0 in OVITO is the **post-transform** geometry, not the raw restart geometry.
- The packaged scripts can be invoked via absolute `.lmp` paths from another working directory. `restart_file` and `output_base` are still interpreted relative to the directory where you launch LAMMPS. The packaged default `helper_root` points to this local clone at `/Users/joshua/Desktop/MD/LAMMPS-LCE-TOLL/Relaxation-and-Stretching/helpers`; if you move the repo, update `helper_root` or override it with `-var helper_root /new/path/to/helpers`.
- The `ppp` relaxation variants can have noticeably different thermodynamics from `fff`. Re-check `box_length` after switching boundary conditions.
- The warning `Temperature for fix modify is not for group all` comes from `fix_modify thermostat_sph temp T_sph`. It is expected in this workflow.
- The relaxation scripts compute `restart_every` from `dump_every`; the stretching scripts use a fixed `restart_every index 50000`.

## Dependency Assumptions

- LAMMPS compatible with the `22 Jul 2025` restart format used here
- packages needed by these scripts, including Gay-Berne and rigid-body integration
- Python 3 available for `helpers/build_output_root.py`
- restart files are **not** included in this repository; supply them locally and update `restart_file` as needed
