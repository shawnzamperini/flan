# Input Options

Input options are set in your input file through the `inpts` map, e.g.:

```cpp
inpts["case_name"] = "my_case";
inpts["imp_num"] = 10000;
inpts["imp_mass_amu"] = 183.84;
```

Any option you don't set uses the default value listed below. The options
are grouped by what they control.

## Simulation

| Option | Type | Default | Description |
|---|---|---|---|
| `case_name` | string | `"undefined"` | Name of the simulation case. The results are saved to `<case_name>.nc`. |
| `use_gpu` | string | `"off"` | Whether to run the simulation on GPUs. Valid values: `"off"`, `"cuda"`. |
| `seed` | int | `4` | Base seed for the random number generators. A negative value seeds them from the clock so each run is different. The MPI rank is added to the seed so each rank gets a different sequence. |
| `slot_cap` | int | `1048576` | Maximum number of particle slots, i.e., the most particles followed at once. The number of slots used is the smaller of this and `imp_num`. |
| `print_interval` | int | `10` | Number of progress updates printed while particles are followed. For example, 10 prints an update every 10% of particles. |

## Background plasma

| Option | Type | Default | Description |
|---|---|---|---|
| `bkg_source` | string | `"gkeyll"` | Source of the background plasma, either a built-in test background or a Gkeyll simulation. Valid values: `"test"`, `"gkeyll"`. |
| `min_ne` | float | `1e15` | Lower limit on the electron density (m^-3). |
| `min_te` | float | `0.1` | Lower limit on the electron temperature (eV). |
| `min_ti` | float | `0.1` | Lower limit on the ion temperature (eV), applied when calculating collisions with the background ions. |

## Testing

| Option | Type | Default | Description |
|---|---|---|---|
| `test_opt` | string | `"gyrate"` | Which built-in test case to run when using the test background. Valid values: `"gyrate"`, `"exb"`, `"gradb"`, `"polarization"`, `"curvature"`, `"friction_force"`. |
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

## Geometry and boundaries

| Option | Type | Default | Description |
|---|---|---|---|
| `tbound_type` | string | `"absorbing"` | Boundary condition at the start and end of the background time range. Valid values: `"absorbing"`, `"periodic"`. |
| `min_xbound_type` | string | `"absorbing"` | Boundary condition at the minimum x boundary. Valid values: `"absorbing"`, `"periodic"`, `"core"`. |
| `max_xbound_type` | string | `"absorbing"` | Boundary condition at the maximum x boundary. Valid values: `"absorbing"`, `"periodic"`, `"core"`. |
| `min_ybound_type` | string | `"absorbing"` | Boundary condition at the minimum y boundary. Valid values: `"absorbing"`, `"periodic"`, `"core"`. |
| `max_ybound_type` | string | `"absorbing"` | Boundary condition at the maximum y boundary. Valid values: `"absorbing"`, `"periodic"`, `"core"`. |
| `min_zbound_type` | string | `"absorbing"` | Boundary condition at the minimum z boundary. Valid values: `"absorbing"`, `"periodic"`, `"core"`. |
| `max_zbound_type` | string | `"absorbing"` | Boundary condition at the maximum z boundary. Valid values: `"absorbing"`, `"periodic"`, `"core"`. |
| `lcfs_x` | float | `0.0` | x coordinate of the last closed flux surface. |
| `lcfs_y` | float | `0.0` | y coordinate of the last closed flux surface. |
| `lcfs_z` | float | `0.0` | z coordinate of the last closed flux surface. |
| `imp_xbound_buffer` | float | `0.0` | Buffer distance at the x boundaries. With an absorbing boundary, a particle is absorbed once it is within this distance of the boundary. |
| `imp_ybound_buffer` | float | `0.0` | Buffer distance at the y boundaries. With an absorbing boundary, a particle is absorbed once it is within this distance of the boundary. |
| `imp_zbound_buffer` | float | `0.0` | Buffer distance at the z boundaries. With an absorbing boundary, a particle is absorbed once it is within this distance of the boundary. |
| `sep_x_bc_xp_z1` | float | `0.0` | First z value for the separatrix x boundary condition near the X-point. |
| `sep_x_bc_xp_z2` | float | `0.0` | Second z value for the separatrix x boundary condition near the X-point. |

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
| `imp_temp_start_opt` | string | `"main_ion"` | How each impurity's starting temperature is chosen: a single value (`imp_temp_start_val`) or the local main ion temperature. Starting velocities are sampled from an isotropic Maxwellian at this temperature. Valid values: `"single_value"`, `"main_ion"`. |
| `imp_temp_start_val` | float | `0.0` | Starting temperature of every impurity when `imp_temp_start_opt` is `"single_value"` (eV). |
| `imp_collisions` | string | `"nanbu"` | Collision model for the impurities, either off or the Nanbu model. Valid values: `"off"`, `"nanbu"`. |
| `imp_time_step_opt` | string | `"constant"` | Whether the impurity time step is constant or variable. Valid values: `"constant"`, `"variable"`. |
| `imp_time_step` | float | `1e-8` | Impurity time step (s). |
| `imp_time_step_min` | float | `1e-12` | Minimum impurity time step when the time step is variable (s). |
| `imp_source_scale_fact` | float | `1.0` | Impurity source scaling factor (particles/s), i.e., how many real particles a single simulated particle represents. Used to convert the Monte Carlo counts into a density. |
| `imp_iz_recomb` | string | `"on"` | Whether impurity ionization and recombination are on. Valid values: `"off"`, `"on"`. |
| `imp_max_time` | float | `1e10` | Maximum time a particle is followed for (s). A particle that exceeds it is killed and counted in a warning at the end of the run. |

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
