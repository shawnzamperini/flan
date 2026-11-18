# Input Options

Input options are set in your input file through the `inpts` map, e.g.:

```cpp
inpts["case_name"] = "my_case";
inpts["imp_num"] = 10000;
inpts["imp_mass_amu"] = 183.84;
```

Any option you don't set uses the default value listed below. The options
are grouped by what they control.

!!! note "Mistakes in options stop the run"
    Flan stops with an error before the simulation starts if an option name
    isn't recognized (e.g., a typo), a value has the wrong type, or a string
    option has a value that isn't one of its valid values. The type must
    match the table exactly. For example, `inpts["imp_num"] = 1e4` is a
    double, so use `inpts["imp_num"] = 10000` for an int option.

## Simulation

| Option | Type | Default | Description |
|---|---|---|---|
| `case_name` | string | `"undefined"` | Name of the simulation case. The results are saved to `<case_name>.nc`. |
| `use_gpu` | string | `"off"` | Whether to run the simulation on GPUs. Valid values: `"off"`, `"cuda"`. |
| `seed` | int | `-1` | Base seed for the random number generators. A negative value (the default) seeds them from the clock so each run is different. Set it to 0 or above to make a run repeatable. The MPI rank is added to the seed so each rank gets a different sequence. |
| `slot_cap` | int | `1048576` | Maximum number of particle slots, i.e., the most particles followed at once. The number of slots used is the smaller of this and `imp_num`. |
| `print_interval` | int | `10` | Number of progress updates printed while particles are followed. For example, 10 prints an update every 10% of particles. |

## Background plasma

| Option | Type | Default | Description |
|---|---|---|---|
| `bkg_source` | string | `"gkeyll"` | Source of the background plasma: the slab background (a box of plasma in simple fields, set entirely by the `slab_*` options and used by the regression cases and unit tests), a Gkeyll simulation, a SOLEDGE3X simulation or a HERMES (JAX-Hermes3) simulation. Valid values: `"slab"`, `"gkeyll"`, `"soledge3x"`, `"hermes"`. |
| `slab_ne` | float | `1e19` | Electron (and ion) density of the slab background (m^-3). Only used with `bkg_source = "slab"`. |
| `slab_te` | float | `10.0` | Electron temperature of the slab background (eV). Only used with `bkg_source = "slab"`. |
| `slab_ti` | float | `10.0` | Ion temperature of the slab background (eV). Only used with `bkg_source = "slab"`. |
| `slab_flow` | float | `0.0` | Plasma flow of the slab background along the magnetic field, which is in +X (m/s). Only used with `bkg_source = "slab"`. |
| `slab_heat_flux_i` | float | `0.0` | Ion parallel conductive heat flux of the slab background, in units of n T_i v_th,i with v_th,i = sqrt(T_i / m_i). Positive is along the field (+X). This is the size of the heat-flux skew that gives the ion thermal force; about 0.1 is a mild skew. Only used with `bkg_source = "slab"`. |
| `slab_heat_flux_e` | float | `0.0` | Electron parallel conductive heat flux of the slab background, in units of n T_e v_th,e with v_th,e = sqrt(T_e / m_e). Positive is along the field (+X). Only used with `bkg_source = "slab"`. |
| `slab_bmag` | float | `1.0` | Magnetic field strength of the slab background (T). The field is uniform and in +X. Only used with `bkg_source = "slab"`. |
| `slab_tmax` | float | `50e-6` | End of the slab background's time range, which starts at 0 (s). The output has 128 time frames evenly spaced over it. Only used with `bkg_source = "slab"`. |
| `slab_efield_y` | float | `0.0` | Electric field of the slab background in +Y at t = 0 (V/m), uniform in space. It is across the magnetic field, so it gives an ExB drift in -Z of E / B. Only for a straight slab (`slab_major_radius = 0`). Only used with `bkg_source = "slab"`. |
| `slab_efield_y_rate` | float | `0.0` | Rate at which the slab background's electric field in +Y changes (V/m/s): E_Y = `slab_efield_y` + `slab_efield_y_rate` * t. A changing field gives a polarization drift along it. Only for a straight slab. Only used with `bkg_source = "slab"`. |
| `slab_bgrad_y` | float | `0.0` | Relative gradient of the slab background's magnetic field strength in +Y (1/m): B_X = `slab_bmag` * (1 + `slab_bgrad_y` * Y). It gives a grad-B drift. Don't make the y boundaries periodic with it. Only for a straight slab. Only used with `bkg_source = "slab"`. |
| `slab_major_radius` | float | `0.0` | If above 0, the slab background's box is bent into a quarter of a torus of this major radius (m), to give the magnetic field curvature: x becomes the major radius R (from 2.5 cm inside to 2.5 cm outside this value), y the height Z and z the toroidal angle (0 to pi/2). The field is then toroidal with strength `slab_bmag` * `slab_major_radius` / R, and the flow and heat fluxes are along it. 0 is the straight slab, with the field along +X. Only used with `bkg_source = "slab"`. |
| `min_ne` | float | `1e15` | Lower limit on the electron density (m^-3). |
| `min_te` | float | `0.1` | Lower limit on the electron temperature (eV). |
| `min_ti` | float | `0.1` | Lower limit on the ion temperature (eV), applied when calculating collisions with the background ions and when sampling starting velocities with `imp_temp_start_opt = "main_ion"`. |

## Testing

| Option | Type | Default | Description |
|---|---|---|---|
| `save_track` | string | `"off"` | Whether to save particle tracks to the output file. Valid values: `"off"`, `"on"`. |

## Gkeyll

| Option | Type | Default | Description |
|---|---|---|---|
| `gkyl_dir` | string | `"undefined"` | Directory containing the Gkeyll simulation files. |
| `gkyl_casename` | string | `"undefined"` | Name of the Gkeyll case, i.e., the prefix of its output files. |
| `gkyl_frame_start` | int | `0` | First Gkeyll frame to load. |
| `gkyl_frame_end` | int | `1` | Last Gkeyll frame to load. The range from `gkyl_frame_start` to `gkyl_frame_end` is inclusive. |
| `gkyl_elec_name` | string | `"elc"` | Name of the electron species in the Gkeyll simulation. |
| `gkyl_ion_name` | string | `"ion"` | Name of the main ion species in the Gkeyll simulation. |
| `gkyl_elec_mass_amu` | float | `0.000548` | Electron mass used in the Gkeyll simulation (amu). |
| `gkyl_ion_mass_amu` | float | `2.014` | Main ion mass used in the Gkeyll simulation (amu). |
| `gkyl_file_type` | string | `"binary"` | File format of the Gkeyll output. Valid values: `"binary"`. |
| `gkyl_moment_type` | string | `"bimaxwellian"` | Type of moment files in the Gkeyll output. Collisions need `"bimaxwellian"`, which provides the ion parallel flow and perpendicular velocity used to calculate the background flow. Valid values: `"bimaxwellian"`, `"maxwellian"`, `"m0m1m2"`. |
| `calc_grad_elec` | string | `"off"` | Whether to calculate the electric field gradient. This is not implemented yet, so turning it on only prints a warning. Valid values: `"off"`, `"on"`. |

## SOLEDGE3X

Used with `bkg_source = "soledge3x"`. See [SOLEDGE3X Backgrounds](running/soledge3x.md).

| Option | Type | Default | Description |
|---|---|---|---|
| `s3x_dir` | string | `"undefined"` | Full path to the directory containing the SOLEDGE3X files: the mesh (`s3x_mesh_file`), the reference file (`s3x_ref_file`) and the plasma files `plasma_00000.h5`, `plasma_00001.h5`, etc. Must be set when `bkg_source = "soledge3x"`. |
| `s3x_mesh_file` | string | `"mesh_S3X.h5"` | Name of the SOLEDGE3X mesh file within `s3x_dir`. It must be the mesh the plasma was run on. |
| `s3x_ref_file` | string | `"refParam_raptorX.h5"` | Name of the SOLEDGE3X reference file within `s3x_dir`, which holds the values SOLEDGE3X's output is normalized by (n0, T0, c0, tau0, phi0, ...). Its R0, a0 and B0 must agree with the mesh file. |
| `s3x_frame_start` | int | `0` | First SOLEDGE3X plasma frame to load, i.e., the number in `plasma_%05d.h5`. |
| `s3x_frame_end` | int | `1` | Last SOLEDGE3X plasma frame to load. The frames from `s3x_frame_start` to `s3x_frame_end` are loaded, which must be at least two, since the background is interpolated in time between frames. |
| `s3x_chi_absorb` | string | `"on"` | SOLEDGE3X marks the cells inside the wall with chi = 1 (its penalization mask). With `"on"`, a particle that ends a time step in a cell with chi = 1 is absorbed, as if it hit the wall. With `"off"`, particles move through the wall cells, which hold SOLEDGE3X's penalized plasma. Valid values: `"on"`, `"off"`. |

## HERMES

Used with `bkg_source = "hermes"`. The background is read from a
single-null hypnotoad grid and the run's snapshots, and split into 6
zones: 1 inner PFR, 2 core, 3 outer PFR, 4 inner leg, 5 main SOL and
6 outer leg. In every zone, x and y are the grid's global cell indices
(x counting the 2 guard cells on each side, so the separatrix is at
x = `ixseps1`) and z is the toroidal angle in radians. The normalization
must be the run's: for JAX-Hermes3 it's in the case's input JSON
(`pars.norm`, `average_atomic_mass`, `zperiod`). The plasma is
interpolated, and differentiated in y, along the field line using the
grid's `zShift`. See the [HERMES Coupling Plan](developer/hermes-coupling.md)
for what isn't included yet.

| Option | Type | Default | Description |
|---|---|---|---|
| `hermes_format` | string | `"jax"` | Output format of the HERMES run: `"jax"` for JAX-Hermes3 (one `.npz` snapshot per time) or `"bout"` for BOUT++ Hermes-3 dump files, which isn't implemented yet. Valid values: `"jax"`, `"bout"`. |
| `hermes_grid_file` | string | `"undefined"` | Full path to the run's hypnotoad grid file (netCDF). Must be set. |
| `hermes_snap_dir` | string | `"undefined"` | Full path to the directory with the snapshots (`snap*.npz` for JAX-Hermes3). Must be set. |
| `hermes_snap_start` | int | `0` | First snapshot to load, as an index (from 0) into the snapshots in `hermes_snap_dir` sorted by name, which for JAX-Hermes3 is time order. |
| `hermes_snap_end` | int | `1` | Last snapshot to load (inclusive), as an index like `hermes_snap_start`. Negative values count from the end: -1 is the last snapshot. At least two snapshots must be loaded, since the background is interpolated in time between them. |
| `hermes_snap_stride` | int | `1` | Load every `hermes_snap_stride`-th snapshot from `hermes_snap_start` to `hermes_snap_end`. Each loaded snapshot holds many 3D fields, so this keeps long time windows within memory. |
| `hermes_nnorm` | float | `1e19` | Density (m^-3) the run is normalized by (Nnorm). |
| `hermes_tnorm` | float | `100.0` | Temperature (eV) the run is normalized by (Tnorm). Potentials are in units of Tnorm too, and velocities in units of sqrt(e Tnorm / m_p). |
| `hermes_bnorm` | float | `1.0` | Magnetic field (T) the run is normalized by (Bnorm). Times are in units of m_p / (e Bnorm). |
| `hermes_ion_mass_amu` | float | `2.0` | Mass number (AA) of the main ions. The parallel flow is NVi / (AA Ne), and this is also the ion mass the impurities collide with. |
| `hermes_zperiod` | int | `1` | The run covers 1 / `hermes_zperiod` of the torus in the toroidal direction (zperiod in the HERMES input). |

## Geometry and boundaries

| Option | Type | Default | Description |
|---|---|---|---|
| `tbound_type` | string | `"absorbing"` | Boundary condition at the start and end of the background time range. Valid values: `"absorbing"`, `"periodic"`. |
| `min_xbound_type` | string | `"absorbing"` | Boundary condition at the minimum x boundary. Valid values: `"absorbing"`, `"periodic"`, `"core"`. |
| `max_xbound_type` | string | `"absorbing"` | Boundary condition at the maximum x boundary. Valid values: `"absorbing"`, `"periodic"`, `"core"`. |
| `min_ybound_type` | string | `"absorbing"` | Boundary condition at the minimum y boundary. Valid values: `"absorbing"`, `"periodic"`, `"core"`. |
| `max_ybound_type` | string | `"absorbing"` | Boundary condition at the maximum y boundary. Valid values: `"absorbing"`, `"periodic"`, `"core"`. |
| `min_zbound_type` | string | `"absorbing"` | Boundary condition at the minimum z boundary. Valid values: `"absorbing"`, `"periodic"`, `"core"`, `"lcfs"`. `"lcfs"` is periodic for particles inside the last closed flux surface (`x < lcfs_x`) and absorbing outside it (assumes x increases outward). It must be used on both z sides and is CPU-only for now. |
| `max_zbound_type` | string | `"absorbing"` | Boundary condition at the maximum z boundary. Valid values: `"absorbing"`, `"periodic"`, `"core"`, `"lcfs"`. See `min_zbound_type`. |
| `lcfs_x` | float | `0.0` | x coordinate of the last closed flux surface. Required (must be inside the background x range) when the z boundaries are `"lcfs"`. |
| `lcfs_y` | float | `0.0` | y coordinate of the last closed flux surface. |
| `lcfs_z` | float | `0.0` | z coordinate of the last closed flux surface. |
| `imp_xbound_buffer` | float | `0.0` | Buffer distance at the x boundaries. With an absorbing boundary, a particle is absorbed once it is within this distance of the boundary. |
| `imp_ybound_buffer` | float | `0.0` | Buffer distance at the y boundaries. With an absorbing boundary, a particle is absorbed once it is within this distance of the boundary. |
| `imp_zbound_buffer` | float | `0.0` | Buffer distance at the z boundaries. With an absorbing boundary, a particle is absorbed once it is within this distance of the boundary. |
| `sep_x_bc_xp_z1` | float | `0.0` | First z value for the separatrix x boundary condition near the X-point. |
| `sep_x_bc_xp_z2` | float | `0.0` | Second z value for the separatrix x boundary condition near the X-point. |

!!! note "Periodic boundaries must be periodic on both sides"
    A periodic boundary wraps a particle around to the other side of the
    domain, so it only makes sense if both the minimum and maximum
    boundaries in that direction are periodic. If only one side of x, y or z
    is set to `"periodic"`, Flan prints a warning and sets both sides to
    `"periodic"`.

## Impurity properties

| Option | Type | Default | Description |
|---|---|---|---|
| `imp_atom_num` | int | `74` | Atomic number of the impurity. |
| `imp_mass_amu` | float | `183.84` | Mass of the impurity (amu). |
| `imp_init_charge` | int | `1` | Charge state the impurities start with. |

## Impurity transport

| Option | Type | Default | Description |
|---|---|---|---|
| `imp_num` | int | `1` | Total number of impurity particles to follow. |
| `imp_tstart_opt` | string | `"single_value"` | How each impurity's starting time is chosen: a single value (`imp_tstart_val`), uniformly between `imp_trange_min` and `imp_trange_max`, or uniformly over the full background time range. Valid values: `"single_value"`, `"range"`, `"full_range"`. |
| `imp_tstart_val` | float | `0.0` | Starting time of every impurity when `imp_tstart_opt` is `"single_value"` (s). |
| `imp_trange_min` | float | `0.0` | Lower limit of the starting time range when `imp_tstart_opt` is `"range"` (s). |
| `imp_trange_max` | float | `0.0` | Upper limit of the starting time range when `imp_tstart_opt` is `"range"` (s). |
| `imp_xstart_opt` | string | `"single_value"` | How each impurity's starting x coordinate is chosen: a single value (`imp_xstart_val`), uniformly between `imp_xrange_min` and `imp_xrange_max`, or uniformly over the full x range of the background. Valid values: `"single_value"`, `"range"`, `"full_range"`. |
| `imp_xstart_val` | float | `0.0` | Starting x coordinate of every impurity when `imp_xstart_opt` is `"single_value"`. |
| `imp_xrange_min` | float | `0.0` | Lower limit of the starting x range when `imp_xstart_opt` is `"range"`. |
| `imp_xrange_max` | float | `0.0` | Upper limit of the starting x range when `imp_xstart_opt` is `"range"`. |
| `imp_ystart_opt` | string | `"single_value"` | How each impurity's starting y coordinate is chosen: a single value (`imp_ystart_val`), uniformly between `imp_yrange_min` and `imp_yrange_max`, or uniformly over the full y range of the background. Valid values: `"single_value"`, `"range"`, `"full_range"`. |
| `imp_ystart_val` | float | `0.0` | Starting y coordinate of every impurity when `imp_ystart_opt` is `"single_value"`. |
| `imp_yrange_min` | float | `0.0` | Lower limit of the starting y range when `imp_ystart_opt` is `"range"`. |
| `imp_yrange_max` | float | `0.0` | Upper limit of the starting y range when `imp_ystart_opt` is `"range"`. |
| `imp_zstart_opt` | string | `"single_value"` | How each impurity's starting z coordinate is chosen: a single value (`imp_zstart_val`), uniformly between `imp_zrange_min` and `imp_zrange_max`, or uniformly over the full z range of the background. Valid values: `"single_value"`, `"range"`, `"full_range"`. |
| `imp_zstart_val` | float | `0.0` | Starting z coordinate of every impurity when `imp_zstart_opt` is `"single_value"`. |
| `imp_zrange_min` | float | `0.0` | Lower limit of the starting z range when `imp_zstart_opt` is `"range"`. |
| `imp_zrange_max` | float | `0.0` | Upper limit of the starting z range when `imp_zstart_opt` is `"range"`. |
| `imp_zone_start_val` | int | `1` | Zone particles start in, numbered from 1 like the SOLEDGE3X zones. The starting x, y and z (`imp_xstart_opt`, etc.) are in that zone's coordinates. Gkeyll and slab backgrounds have a single zone, so this must be 1 for them. |
| `imp_temp_start_opt` | string | `"main_ion"` | How each impurity's starting temperature is chosen: a single value (`imp_temp_start_val`) or the local main ion temperature. Starting velocities are sampled from an isotropic Maxwellian at this temperature. Valid values: `"single_value"`, `"main_ion"`. |
| `imp_temp_start_val` | float | `0.0` | Starting temperature of every impurity when `imp_temp_start_opt` is `"single_value"` (eV). |
| `imp_vX_start_val` | float | `0.0` | Velocity in Cartesian X that every impurity starts with (m/s), added to the velocity sampled at the starting temperature. With `imp_temp_start_val = 0` every impurity starts with exactly this velocity. Not used for sputtered impurities. |
| `imp_vY_start_val` | float | `0.0` | Same as `imp_vX_start_val`, in Cartesian Y (m/s). |
| `imp_vZ_start_val` | float | `0.0` | Same as `imp_vX_start_val`, in Cartesian Z (m/s). |
| `imp_collisions` | string | `"nanbu"` | Collision model for the impurities, either off or the Nanbu model. Valid values: `"off"`, `"nanbu"`. |
| `imp_thermal_force` | string | `"on"` | Whether the Nanbu collisions include the electron and ion thermal forces, by colliding impurities with background distributions skewed by their parallel heat fluxes. The heat fluxes come from SOLEDGE3X plasma files with the fluxes on the phi faces; other backgrounds have none, so a warning says the thermal force isn't included. With heat fluxes, it runs on the CPU only. Valid values: `"off"`, `"on"`. |
| `imp_coll_subcycle_target` | float | `0.03` | Collision sub-cycling: when a time step is too collisional for the Nanbu model, the collision is split into sub-steps, each with a newly sampled background particle, so that each has $s_0$ (the Nanbu $s$ for a background particle at its thermal speed) at or below this value. Smaller values bring friction and the thermal force closer to their fluid values and cost more sub-steps. Must be above 0. See [The time step](physics/thermal-force.md#the-time-step). |
| `imp_coll_subcycle_cap` | int | `128` | Largest number of sub-steps one collision is split into. `1` turns collision sub-cycling off. Must be at least 1. |
| `imp_time_step_opt` | string | `"constant"` | Whether the impurity time step is constant or variable. Valid values: `"constant"`, `"variable"`. |
| `imp_time_step` | float | `1e-8` | Impurity time step (s). |
| `imp_time_step_min` | float | `1e-12` | Minimum impurity time step when the time step is variable (s). |
| `imp_source_scale_fact` | float | `1.0` | How many impurities the source releases, shared between the `imp_num` particles, so the density comes out in m^-3. For a source spread over time (`imp_tstart_opt = "full_range"` or `"range"`) it's a rate (atoms/s); for a pulse (`"single_value"`) it's a number of atoms. With the default of 1, the density is per atom/s (or per atom) of source. Each time frame's density is the density near that frame's time. Not used with `imp_source_opt = "sputter"`, which uses the sputtered rate instead. |
| `imp_iz_recomb` | string | `"on"` | Whether impurity ionization and recombination are on. Valid values: `"off"`, `"on"`. |
| `imp_max_time` | float | `1e10` | Maximum time a particle is followed for (s). A particle that exceeds it is killed and counted in a warning at the end of the run. |

## Sputtering

Impurities sputtered off the walls of a SOLEDGE3X background by the
background plasma. See [Sputtering](physics/sputtering.md).

| Option | Type | Default | Description |
|---|---|---|---|
| `imp_source_opt` | string | `"manual"` | Where impurities come from. `"manual"` starts them where the `imp_*start_opt` options say. `"sputter"` makes them the atoms the background plasma sputters off the wall of a SOLEDGE3X background (CPU only): the rate off each part of the wall is calculated and saved to the output file, particles start from it neutral, and its total rate replaces `imp_source_scale_fact` (each particle stands for that rate / `imp_num` atoms/s). The `imp_*start_opt`, `imp_zone_start_val`, `imp_temp_start_opt` and `imp_init_charge` options aren't used then. Valid values: `"manual"`, `"sputter"`. |
| `sput_table_file` | string | `"undefined"` | Path to the sputtering yield table, e.g., `data/sputtering/D_on_W.h5`. Needed with `imp_source_opt = "sputter"`. |
| `sput_flux_opt` | string | `"parallel_exb"` | Which part of the background ion flow carries ions onto the wall. `"parallel_exb"` is the parallel flow plus the E x B drift, `"parallel"` the parallel flow only. Valid values: `"parallel_exb"`, `"parallel"`. |
| `sput_bohm_floor` | string | `"off"` | Whether the parallel flow onto the wall is at least the ion sound speed (the Bohm criterion), for backgrounds whose flow next to the wall is too slow. Valid values: `"off"`, `"on"`. |
| `sput_angle_val` | float | `60` | Angle at which the background ions hit the wall, in degrees from the surface normal. Used to look up the sputtering yield. |
| `sput_max_wall_dist` | float | `2` | How far, in cell widths, the wall outline can be from a wall face before the face uses its own normal and area instead of the outline's. |

## Variance reduction (temporarily removed)

| Option | Type | Default | Description |
|---|---|---|---|
| `var_red_split` | string | `"off"` | When particles are split for variance reduction: never, at ionization/recombination events, or at collisions. Valid values: `"off"`, `"iz_rec"`, `"coll"`. |
| `var_red_import` | string | `"median"` | Method used to calculate particle importance for variance reduction. Valid values: `"median"`, `"exp_dist"`, `"exp_time"`. |
| `var_red_freq` | float | `0.1` | How often variance reduction is applied. |
| `var_red_min_weight` | float | `0.1` | Minimum particle weight for variance reduction. |
| `var_red_med_mod` | float | `1.0` | Modifier applied to the median importance for variance reduction. |
| `var_red_rusrol` | string | `"off"` | Whether to use Russian roulette for variance reduction. Valid values: `"off"`, `"on"`. |
| `var_red_rusrol_prob` | float | `0.5` | Probability used for Russian roulette. |

## OpenADAS

| Option | Type | Default | Description |
|---|---|---|---|
| `openadas_root` | string | `"undefined"` | Directory containing the OpenADAS rate coefficient files. |
| `openadas_year` | int | `50` | Year of the OpenADAS data files to use. |
