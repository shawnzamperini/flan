======================================================================================================
Input options
======================================================================================================

This page contains all the input options available to the user, their definition/allowed values, and finally the default values that each variable takes on. 

Last Update: 8/22/25 (commit XXX)

.. list-table:: Input Options
   :header-rows: 1
   :widths: 30 50 20

   * - Variable
     - Description
     - Default

   * - :ref:`case_name <case_name>`
     - Name of the simulation
     - "undefined"

   * - :ref:`gkyl_dir <gkyl_dir>`
     - Full path to the directory containing all the Gkeyll files
     - "undefined"

   * - :ref:`gkyl_casename <gkyl_casename>`
     - Gkeyll simulation name within :ref:`gkyl_dir <gkyl_dir>`
     - "undefined"

   * - :ref:`gkyl_frame_start <gkyl_frame_start>`
     - Starting frame index
     - 0

   * - :ref:`gkyl_frame_end <gkyl_frame_end>`
     - Ending frame index
     - 1

   * - :ref:`gkyl_elec_name <gkyl_elec_name>`
     - Electron species name in Gkeyll simulation
     - "elc"

   * - :ref:`gkyl_ion_name <gkyl_ion_name>`
     - Ion species name in Gkeyll simulation
     - "ion"

   * - :ref:`gkyl_elec_mass_amu <gkyl_elec_mass_amu>`
     - Electron mass (amu)
     - 0.000548

   * - :ref:`gkyl_ion_mass_amu <gkyl_ion_mass_amu>`
     - Ion mass (amu)
     - 2.014

   * - :ref:`gkyl_file_type <gkyl_file_type>`
     - Gkeyll file format
     - "binary"

   * - :ref:`gkyl_moment_type <gkyl_moment_type>`
     - Moment type used
     - "bimaxwellian"

   * - :ref:`lcfs_x <lcfs_x>`
     - LCFS x-location in Gkeyll simulation (if applicable)
     - 0.0

   * - :ref:`imp_xbound_buffer <imp_xbound_buffer>`
     - Move absorbing x boundary condition off boundary by this much
     - 0.0

   * - :ref:`min_xbound_type <min_xbound_type>`
     - Min x-bound condition
     - "absorbing"

   * - :ref:`imp_atom_num <imp_atom_num>`
     - Impurity atomic number
     - 74

   * - :ref:`imp_mass_amu <imp_mass_amu>`
     - Impurity mass (amu)
     - 183.84

   * - :ref:`imp_init_charge <imp_init_charge>`
     - Initial impurity charge
     - 1

   * - :ref:`imp_num <imp_num>`
     - Number of primary impurities to follow
     - 1

   * - :ref:`imp_tstart_opt <imp_tstart_opt>`
     - Time start option
     - "single_value"

   * - :ref:`imp_tstart_val <imp_tstart_val>`
     - Time start value when :ref:`imp_tstart_opt <imp_tstart_opt>` = "single_value"
     - 0.0

   * - :ref:`imp_trange_min <imp_trange_min>`
     - Time range minimum when :ref:`imp_tstart_opt <imp_tstart_opt>` = "range"
     - 0.0

   * - :ref:`imp_trange_max <imp_trange_max>`
     - Time range maximum when :ref:`imp_tstart_opt <imp_tstart_opt>` = "range"
     - 0.0

   * - :ref:`imp_xstart_opt <imp_xstart_opt>`
     - x start option
     - "single_value"

   * - :ref:`imp_xstart_val <imp_xstart_val>`
     - X start value when :ref:`imp_xstart_opt <imp_xstart_opt>` = "single_value"
     - 0.0

   * - :ref:`imp_xrange_min <imp_xrange_min>`
     - x range minimum when :ref:`imp_xstart_opt <imp_xstart_opt>` = "range"
     - 0.0

   * - :ref:`imp_xrange_max <imp_xrange_max>`
     - x range maximum when :ref:`imp_xstart_opt <imp_xstart_opt>` = "range"
     - 0.0

   * - :ref:`imp_ystart_opt <imp_ystart_opt>`
     - y start option
     - "single_value"

   * - :ref:`imp_ystart_val <imp_ystart_val>`
     - y start value when :ref:`imp_ystart_opt <imp_ystart_opt>` = "single_value"
     - 0.0

   * - :ref:`imp_yrange_min <imp_yrange_min>`
     - y range minimum when :ref:`imp_ystart_opt <imp_ystart_opt>` = "range"
     - 0.0

   * - :ref:`imp_yrange_max <imp_yrange_max>`
     - y range maximum when :ref:`imp_ystart_opt <imp_ystart_opt>` = "range"
     - 0.0

   * - :ref:`imp_zstart_opt <imp_zstart_opt>`
     - Z start option
     - "single_value"

   * - :ref:`imp_zstart_val <imp_zstart_val>`
     - Z start value when :ref:`imp_zstart_opt <imp_zstart_opt>` = "single_value"
     - 0.0

   * - :ref:`imp_zrange_min <imp_zrange_min>`
     - Z range minimum when :ref:`imp_zstart_opt <imp_zstart_opt>` = "range"
     - 0.0

   * - :ref:`imp_zrange_max <imp_zrange_max>`
     - Z range maximum when :ref:`imp_zstart_opt <imp_zstart_opt>` = "range"
     - 0.0

   * - :ref:`imp_collisions <imp_collisions>`
     - Impurity collisions toggle
     - "off"

   * - :ref:`imp_time_step_opt <imp_time_step_opt>`
     - Time step option
     - "variable"

   * - :ref:`imp_time_step <imp_time_step>`
     - Time step value
     - 1e-7

   * - :ref:`imp_time_step_min <imp_time_step_min>`
     - Minimum time step
     - 1e-12

   * - :ref:`imp_source_scale_fact <imp_source_scale_fact>`
     - Source scaling factor
     - 1.0

   * - :ref:`imp_vel_stats <imp_vel_stats>`
     - Velocity statistics toggle
     - "off"

   * - :ref:`imp_iz_recomb <imp_iz_recomb>`
     - Ionization/recombination toggle
     - "on"

   * - :ref:`print_interval <print_interval>`
     - Print interval
     - 10

   * - :ref:`var_red_split <var_red_split>`
     - Variance reduction split toggle
     - "off"

   * - :ref:`var_red_import <var_red_import>`
     - Import method
     - "median"

   * - :ref:`var_red_freq <var_red_freq>`
     - Frequency threshold
     - 0.1

   * - :ref:`var_red_min_weight <var_red_min_weight>`
     - Minimum weight
     - 0.1

   * - :ref:`var_red_med_mod <var_red_med_mod>`
     - Median modifier
     - 1.0

   * - :ref:`var_red_rusrol <var_red_rusrol>`
     - Russian roulette toggle
     - "off"

   * - :ref:`var_red_rusrol_prob <var_red_rusrol_prob>`
     - Russian roulette probability
     - 0.5

   * - :ref:`openadas_root <openadas_root>`
     - OpenADAS root path
     - "undefined"

   * - :ref:`openadas_year <openadas_year>`
     - OpenADAS year
     - 50



case_name
---------


gkyl_dir
--------



gkyl_casename
-------------



gkyl_frame_start
----------------



gkyl_frame_end
--------------



gkyl_elec_name
--------------



gkyl_ion_name
-------------



gkyl_elec_mass_amu
------------------



gkyl_ion_mass_amu
-----------------



gkyl_file_type
--------------


gkyl_moment_type
----------------



lcfs_x
------



imp_xbound_buffer
-----------------



min_xbound_type
---------------



imp_atom_num
------------



imp_mass_amu
------------



imp_init_charge
---------------



imp_num
-------



imp_tstart_opt
--------------



imp_tstart_val
--------------



imp_trange_min
--------------



imp_trange_max
--------------



imp_xstart_opt
--------------



imp_xstart_val
--------------



imp_xrange_min
--------------



imp_xrange_max
--------------



imp_ystart_opt
--------------



imp_ystart_val
--------------



imp_yrange_min
--------------



imp_yrange_max
--------------



imp_zstart_opt
--------------



imp_zstart_val
--------------



imp_zrange_min
--------------



imp_zrange_max
--------------



imp_collisions
--------------



imp_time_step_opt
-----------------



imp_time_step
-------------



imp_time_step_min
-----------------



imp_source_scale_fact
---------------------



imp_vel_stats
-------------



imp_iz_recomb
-------------



print_interval
--------------



var_red_split
-------------



var_red_import
--------------



var_red_freq
------------



var_red_min_weight
------------------



var_red_med_mod
---------------



var_red_rusrol
--------------



var_red_rusrol_prob
-------------------



openadas_root
-------------



openadas_year
-------------
