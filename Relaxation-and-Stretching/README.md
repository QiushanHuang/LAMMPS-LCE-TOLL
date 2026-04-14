# Relaxation and Stretching Workflow

This package collects restart-based LAMMPS input scripts for two stages:

- `relaxation/`: rebuild the chain geometry from restart `xu/yu/zu`, expand the box, rotate/center the cluster, then run free thermal relaxation in either `fff` or `ppp`
- `stretching/`: take a relaxed restart, align the chain into a slender box, pin the head, and apply a constant force to the tail

The scripts share a helper in `helpers/build_output_root.py` that builds the output directory name from the restart file, target temperature, force magnitude, and protocol tag.

## Directory Layout

- `helpers/build_output_root.py`: generates `output_root` names
- `relaxation/in.restart_relaxation_fff_xuyu-zu.lmp`: baseline `fff`, `T*=0.50`, cubic `box_length=200`
- `relaxation/in.restart_relaxation_fff_xuyu-zu-08.lmp`: `fff`, `T*=0.80`, cubic `box_length=400`
- `relaxation/in.restart_relaxation_ppp_xuyu-zu-05.lmp`: `ppp`, `T*=0.50`, cubic `box_length=200`
- `relaxation/in.restart_relaxation_ppp_xuyu-zu-07.lmp`: `ppp`, `T*=0.70`, cubic `box_length=300`
- `relaxation/in.restart_relaxation_ppp_xuyu-zu-08.lmp`: `ppp`, `T*=0.80`, cubic `box_length=300`
- `stretching/in.single_chain_constant_force_aligned_box-F10.00-tail_v0_no_thermal-xuyu-zu.lmp`: force-clamp stretching, `F=10`
- `stretching/in.single_chain_constant_force_aligned_box-F20.00-tail_v0_no_thermal-xuyu-zu.lmp`: force-clamp stretching, `F=20`, tail velocity zeroed before loading
- `stretching/in.single_chain_constant_force_aligned_box-F20.00-xuyu-zu.lmp`: force-clamp stretching, `F=20`

## Endpoint Selection

All packaged scripts now support automatic endpoint detection:

- `head_id_override index auto`
- `tail_id_override index auto`

Default behavior:

- the smallest atom ID is used as `head_id`
- the largest atom ID is used as `tail_id`

If your restart does **not** place the physical chain ends at the minimum and maximum atom IDs, override them manually by editing the script or from the command line, for example:

```bash
lmp_mpi -var head_id_override 1 -var tail_id_override 8006 -in in.single_chain_constant_force_aligned_box-F20.00-tail_v0_no_thermal-xuyu-zu.lmp
```

Use the **same** endpoint override pair in both relaxation and stretching when you want a consistent workflow.

## How to Use

### 1. Run relaxation

Change into the chosen relaxation directory and run one of the scripts:

```bash
cd relaxation
mpirun -np 4 lmp_mpi -in in.restart_relaxation_ppp_xuyu-zu-05.lmp
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

For the `ppp` variants, only the **final** relaxation box is periodic; the transform stage remains nonperiodic to avoid premature wrapping during rotation and translation.

### 2. Feed the relaxed restart into stretching

Choose a relaxation restart or `Final.relaxation.bin`, then update the stretching script:

- set `restart_file`
- keep `target_T` consistent if that is part of your protocol
- keep `head_id_override` and `tail_id_override` consistent

Then run, for example:

```bash
cd ../stretching
mpirun -np 4 lmp_mpi -in in.single_chain_constant_force_aligned_box-F20.00-tail_v0_no_thermal-xuyu-zu.lmp
```

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
- The helper path is relative. Run scripts from inside their own `relaxation/` or `stretching/` directory unless you intentionally want outputs in a different working directory.
- The `ppp` relaxation variants can have noticeably different thermodynamics from `fff`. Re-check `box_length` after switching boundary conditions.
- The warning `Temperature for fix modify is not for group all` comes from `fix_modify thermostat_sph temp T_sph`. It is expected in this workflow.
- The relaxation scripts compute `restart_every` from `dump_every`; the stretching scripts use a fixed `restart_every index 50000`.

## Dependency Assumptions

- LAMMPS compatible with the `22 Jul 2025` restart format used here
- packages needed by these scripts, including Gay-Berne and rigid-body integration
- Python 3 available for `helpers/build_output_root.py`
- restart files are **not** included in this repository; supply them locally and update `restart_file` as needed
