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
| `case_name` | string | `"undefined"` | TODO |
| `use_gpu` | string | `"off"` | TODO |
| `seed` | int | `4` | TODO |
| `slot_cap` | int | `1048576` | TODO |
| `print_interval` | int | `10` | TODO |

## Background plasma

| Option | Type | Default | Description |
|---|---|---|---|
| `bkg_source` | string | `"gkeyll"` | TODO |
| `min_ne` | float | `1e15` | TODO |
| `min_te` | float | `0.1` | TODO |
| `min_ti` | float | `0.1` | TODO |

## Testing

| Option | Type | Default | Description |
|---|---|---|---|
| `test_opt` | string | `"gyrate"` | TODO |
| `save_track` | string | `"off"` | TODO |

## Gkeyll

| Option | Type | Default | Description |
|---|---|---|---|
| `gkyl_dir` | string | `"undefined"` | TODO |
| `gkyl_casename` | string | `"undefined"` | TODO |
| `gkyl_frame_start` | int | `0` | TODO |
| `gkyl_frame_end` | int | `1` | TODO |
| `gkyl_elec_name` | string | `"elc"` | TODO |
| `gkyl_ion_name` | string | `"ion"` | TODO |
| `gkyl_elec_mass_amu` | float | `0.000548` | TODO |
| `gkyl_ion_mass_amu` | float | `2.014` | TODO |
| `gkyl_file_type` | string | `"binary"` | TODO |
| `gkyl_moment_type` | string | `"bimaxwellian"` | TODO |
| `calc_grad_elec` | string | `"off"` | TODO |

## Geometry and boundaries

| Option | Type | Default | Description |
|---|---|---|---|
| `tbound_type` | string | `"absorbing"` | TODO |
| `min_xbound_type` | string | `"absorbing"` | TODO |
| `max_xbound_type` | string | `"absorbing"` | TODO |
| `min_ybound_type` | string | `"absorbing"` | TODO |
| `max_ybound_type` | string | `"absorbing"` | TODO |
| `min_zbound_type` | string | `"absorbing"` | TODO |
| `max_zbound_type` | string | `"absorbing"` | TODO |
| `lcfs_x` | float | `0.0` | TODO |
| `lcfs_y` | float | `0.0` | TODO |
| `lcfs_z` | float | `0.0` | TODO |
| `imp_xbound_buffer` | float | `0.0` | TODO |
| `imp_ybound_buffer` | float | `0.0` | TODO |
| `imp_zbound_buffer` | float | `0.0` | TODO |
| `sep_x_bc_xp_z1` | float | `0.0` | TODO |
| `sep_x_bc_xp_z2` | float | `0.0` | TODO |

## Impurity properties

| Option | Type | Default | Description |
|---|---|---|---|
| `imp_atom_num` | int | `74` | TODO |
| `imp_mass_amu` | float | `183.84` | TODO |
| `imp_init_charge` | int | `1` | TODO |

## Impurity transport

| Option | Type | Default | Description |
|---|---|---|---|
| `imp_num` | int | `1` | TODO |
| `imp_tstart_opt` | string | `"single_value"` | TODO |
| `imp_tstart_val` | float | `0.0` | TODO |
| `imp_trange_min` | float | `0.0` | TODO |
| `imp_trange_max` | float | `0.0` | TODO |
| `imp_xstart_opt` | string | `"single_value"` | TODO |
| `imp_xstart_val` | float | `0.0` | TODO |
| `imp_xrange_min` | float | `0.0` | TODO |
| `imp_xrange_max` | float | `0.0` | TODO |
| `imp_ystart_opt` | string | `"single_value"` | TODO |
| `imp_ystart_val` | float | `0.0` | TODO |
| `imp_yrange_min` | float | `0.0` | TODO |
| `imp_yrange_max` | float | `0.0` | TODO |
| `imp_zstart_opt` | string | `"single_value"` | TODO |
| `imp_zstart_val` | float | `0.0` | TODO |
| `imp_zrange_min` | float | `0.0` | TODO |
| `imp_zrange_max` | float | `0.0` | TODO |
| `imp_temp_start_opt` | string | `"main_ion"` | TODO |
| `imp_temp_start_val` | float | `0.0` | TODO |
| `imp_collisions` | string | `"nanbu"` | TODO |
| `imp_time_step_opt` | string | `"constant"` | TODO |
| `imp_time_step` | float | `1e-8` | TODO |
| `imp_time_step_min` | float | `1e-12` | TODO |
| `imp_source_scale_fact` | float | `1.0` | TODO |
| `imp_iz_recomb` | string | `"on"` | TODO |
| `imp_max_time` | float | `1e10` | TODO |

## Variance reduction

| Option | Type | Default | Description |
|---|---|---|---|
| `var_red_split` | string | `"off"` | TODO |
| `var_red_import` | string | `"median"` | TODO |
| `var_red_freq` | float | `0.1` | TODO |
| `var_red_min_weight` | float | `0.1` | TODO |
| `var_red_med_mod` | float | `1.0` | TODO |
| `var_red_rusrol` | string | `"off"` | TODO |
| `var_red_rusrol_prob` | float | `0.5` | TODO |

## OpenADAS

| Option | Type | Default | Description |
|---|---|---|---|
| `openadas_root` | string | `"undefined"` | TODO |
| `openadas_year` | int | `50` | TODO |
