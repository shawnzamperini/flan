======================================================================================================
Input options
======================================================================================================

This page contains all the input options available to the user, their definition/allowed values, and finally the default values that each variable takes on. 

Last Update: 8/22/25 (commit XXX)

.. list-table:: Input Options
   :header-rows: 1
   :widths: 30 50 20

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


.. _gkyl_dir:

**gkyl_dir**  
  Describe std::string gkyl_dir (default: "undefined")

.. _gkyl_casename:

**gkyl_casename**  
  Describe std::string gkyl_casename (default: "undefined")

.. _gkyl_frame_start:

**gkyl_frame_start**  
  Describe int gkyl_frame_start (default: 0)

.. _gkyl_frame_end:

**gkyl_frame_end**  
  Describe int gkyl_frame_end (default: 1)

.. _gkyl_elec_name:

**gkyl_elec_name**  
  Describe std::string gkyl_elec_name (default: "elc")

.. _gkyl_ion_name:

**gkyl_ion_name**  
  Describe std::string gkyl_ion_name (default: "ion")

.. _A1:

**gkyl_elec_mass_amu**  
  Describe double gkyl_elec_mass_amu (default: 0.000548)

.. _A2:

**gkyl_ion_mass_amu**  
  Describe double gkyl_ion_mass_amu (default: 2.014)

.. _gkyl_file_type:

**gkyl_file_type**  
  Describe std::string gkyl_file_type (default: "binary")

.. _lcfs_x:

**lcfs_x**  
  Describe double lcfs_x (default: 0.0)

.. _A3:

**imp_xbound_buffer**  
  Describe double imp_xbound_buffer (default: 0.0)

.. _min_xbound_type:

**min_xbound_type**  
  Describe std::string min_xbound_type (default: "absorbing")

.. _imp_atonum:

**imp_atonum**  
  Describe int imp_atonum (default: 74)

.. _imp_mass_amu:

**imp_mass_amu**  
  Describe double imp_mass_amu (default: 183.84)

.. _imp_init_charge:

**imp_init_charge**  
  Describe int imp_init_charge (default: 1)

.. _imp_num:

**imp_num**  
  Describe int imp_num (default: 1)

.. _imp_xmin:

**imp_xmin**  
  Describe double imp_xmin (default: 0.0)

.. _imp_xmax:

**imp_xmax**  
  Describe double imp_xmax (default: 0.0)

.. _imp_ystart_opt:

**imp_ystart_opt**  
  Describe std::string imp_ystart_opt (default: "single_value")

.. _imp_ystart_val:

**imp_ystart_val**  
  Describe double imp_ystart_val (default: 0.0)

.. _imp_zstart_opt:

**imp_zstart_opt**  
  Describe std::string imp_zstart_opt (default: "single_value")

.. _imp_zstart_val:

**imp_zstart_val**  
  Describe double imp_zstart_val (default: 0.0)

.. _imp_collisions:

**imp_collisions**  
  Describe std::string imp_collisions (default: "off")

.. _A4:

**imp_time_step_opt**  
  Describe std::string imp_time_step_opt (default: "variable")

.. _imp_time_step:

**imp_time_step**  
  Describe double imp_time_step (default: 1e-07)

.. _A5:

**_imp_time_step_min**  
  Describe double imp_time_step_min (default: 1e-12)

.. _A6:

**imp_source_scale_fact**  
  Describe double imp_source_scale_fact (default: 1.0)

.. _imp_vel_stats:

**imp_vel_stats**  
  Describe std::string imp_vel_stats (default: "off")

.. _imp_iz_recomb:

**imp_iz_recomb**  
  Describe std::string imp_iz_recomb (default: "on")

.. _print_interval:

**print_interval**  
  Describe int print_interval (default: 10)

.. _var_red:

**var_red**  
  Describe std::string var_red (default: "off")

.. _var_red_mode:

**var_red_mode**  
  Describe std::string var_red_mode (default: "median")

.. _var_red_freq:

**var_red_freq**  
  Describe double var_red_freq (default: 0.1)

.. _A7:

**var_red_min_weight**  
  Describe double var_red_min_weight (default: 0.1)

.. _var_red_med_mod:

**var_red_med_mod**  
  Describe double var_red_med_mod (default: 1.0)

.. _openadas_root:

**openadas_root**  
  Describe std::string openadas_root (default: "undefined")

.. _openadas_year:

**openadas_year**  
  Describe int openadas_year (default: 50)
