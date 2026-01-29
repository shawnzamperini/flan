======================================================================================================
Input Options
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
     - x start value when :ref:`imp_xstart_opt <imp_xstart_opt>` = "single_value"
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
     - z start option
     - "single_value"

   * - :ref:`imp_zstart_val <imp_zstart_val>`
     - z start value when :ref:`imp_zstart_opt <imp_zstart_opt>` = "single_value"
     - 0.0

   * - :ref:`imp_zrange_min <imp_zrange_min>`
     - z range minimum when :ref:`imp_zstart_opt <imp_zstart_opt>` = "range"
     - 0.0

   * - :ref:`imp_zrange_max <imp_zrange_max>`
     - z range maximum when :ref:`imp_zstart_opt <imp_zstart_opt>` = "range"
     - 0.0

   * - :ref:`imp_collisions <imp_collisions>`
     - Impurity collisions with background toggle
     - "on"

   * - :ref:`imp_time_step_opt <imp_time_step_opt>`
     - Time step option
     - "constant"

   * - :ref:`imp_time_step <imp_time_step>`
     - Time step value
     - 1e-8

   * - :ref:`imp_time_step_min <imp_time_step_min>`
     - Minimum time step
     - 1e-12

   * - :ref:`imp_source_scale_fact <imp_source_scale_fact>`
     - Source scaling factor
     - 1.0

   * - :ref:`imp_vel_stats <imp_vel_stats>`
     - Velocity statistics toggle
     - "on"

   * - :ref:`imp_iz_recomb <imp_iz_recomb>`
     - Ionization/recombination toggle
     - "on"

   * - :ref:`print_interval <print_interval>`
     - Print interval
     - 10

   * - :ref:`var_red_split <var_red_split>`
     - Variance reduction split particles toggle
     - "off"

   * - :ref:`var_red_import <var_red_import>`
     - Importance method for variance reduction
     - "median"

   * - :ref:`var_red_freq <var_red_freq>`
     - Frequency at which variance reduction is recalculated
     - 0.1

   * - :ref:`var_red_min_weight <var_red_min_weight>`
     - Minimum allowed weight in variance reduction scheme
     - 0.1

   * - :ref:`var_red_med_mod <var_red_med_mod>`
     - Modifier to consider a fraction of median below which variance reduction occurs
     - 1.0

   * - :ref:`var_red_rusrol <var_red_rusrol>`
     - Variance reduction Russian roulette toggle
     - "off"

   * - :ref:`var_red_rusrol_prob <var_red_rusrol_prob>`
     - Russian roulette probability in low-statistic regions
     - 0.5

   * - :ref:`openadas_root <openadas_root>`
     - OpenADAS root path
     - "undefined"

   * - :ref:`openadas_year <openadas_year>`
     - OpenADAS year
     - 50



case_name
---------

This is only used when naming the NetCDf file. Useful when doing parameter scans, since it allows the user to reuse the same background files for different iterations of the same simulation.

gkyl_dir
--------

Must be the full path to the directory containing the Gkeyll background files, e.g., "/home/zamp/gkyldir/my_gkyl_sim".

gkyl_casename
-------------

Name of the Gkeyll simulation contained within :ref:`gkyl_dir <gkyl_dir>`. All files must share the same case name.

gkyl_frame_start
----------------

Gkeyll frame at which the Flan simulation starts from.

gkyl_frame_end
--------------

Gkeyll frame at which the Flan simulation ends at (inclusive).

gkyl_elec_name
--------------

Name of the electron species in the Gkeyll simulation, e.g., "elc".

gkyl_ion_name
-------------

Name of the ion species in the Gkeyll simulation, e.g., "ion".

gkyl_elec_mass_amu
------------------

Mass of electron in Gkeyll simulation in atomic mass units, e.g., 0.000548.

gkyl_ion_mass_amu
-----------------

Mass of ion in Gkeyll simulation in atomic mass units, e.g., 2.014.

gkyl_file_type
--------------

(Slated for removal)
The format the Gkeyll files are saved in. Currently, only "binary" (.gkyl) is supported and there is no expectation to extend support beyond this. This is the format output by gkylzero.

gkyl_moment_type
----------------

(Slated for removal)
The type of moments output by Gkeyll. Right now, only "bimaxwellian" is really supported. You can run with "maxwellian", but the collision model will be incomplete.

lcfs_x
------

The x coordinate that corresponds to the last closed flux surface, when applicable. This has implications on the boundary conditions, since boundary conditions in the core are treated differently from those in the SOL. When using a geometry that is only in the SOL, leave as 0.0 (or more technically, something below the minimum x bound). Conversely, if a simulation only takes place in the core, set this to a value above the maximum x bound. This assumes the x coordinate is the radial coordinate, which is traditional but not necessarily required in Gkeyll.

Note: This introduces potential pitfalls for the user if they are not aware of this option. Could be likely be improved.

imp_xbound_buffer
-----------------

(Slated for removal)
This is mostly an experimental option. It is possible that the Gkeyll plasma behaves strangely at the x bounds due to artificial boundary conditions. If this input option is a value greater than zero, then Flan will trigger the x boundary conditions when the impurity is within that much distance of either x boundary. This is a bit of a hack option, and hasn't really been necessary so it will probably be removed.

min_xbound_type
---------------

Type of boundary condition to impose at the minimum x boundary. 

  - **"absorbing"**: Particle following is ended when encountering minimum x boundary.
  - **"core"**: Particle is teleported to a random y,z cell along the minimum x boundary. This mimics a particle entering the not-simulated core region and coming out at some other location instantaneously. This is fine for steady-state simulations, but for time-varying simulations (e.g., simulating a pellet injection) one should consider that in the real world some amount of time will pass before the particle is ejected from the not-simulated core region.

imp_atom_num
------------

Atomic number of impurity.

imp_mass_amu
------------

Atomic mass of the impurity.

imp_init_charge
---------------

Initial charge of the impurity.

imp_num
-------

Number of primary impurities to follow. If a variance reduction scheme is on, then secondary impurities could be generated and thus the number of followed particles could be significantly larger than the value entered here. 

imp_tstart_opt
--------------

Option for where the particle starts in time. t=0 corresponds to the start of the Flan simulation (NOT the start of the Gkeyll simulation).

  - **"single_value"**: Particle starts at a specific time designated by :ref:`imp_tstart_opt <imp_tstart_opt>`. Typical for time-varying simulations.
  - **"range"**: Particle is uniformily distributed between :ref:`imp_trange_min <imp_trange_min>` and :ref:`imp_trange_max <imp_trange_max>`. Typical for time-varying simulations.
  - **"full_range"**: Particle is uniformily distributed between the full time range of the simulation. Typical for steady-state simulations.

imp_tstart_val
--------------

Starting time value in seconds when :ref:`imp_tstart_opt <imp_tstart_opt>` = "single_value".

imp_trange_min
--------------

Minimum time value in seconds when :ref:`imp_tstart_opt <imp_tstart_opt>` = "range".

imp_trange_max
--------------

Maximum time value in seconds when :ref:`imp_tstart_opt <imp_tstart_opt>` = "range".

imp_xstart_opt
--------------

Option for what x the particle starts.

  - **"single_value"**: Particle starts at a specific x designated by :ref:`imp_xstart_opt <imp_xstart_opt>`.
  - **"range"**: Particle is uniformily distributed between :ref:`imp_xrange_min <imp_xrange_min>` and :ref:`imp_xrange_max <imp_xrange_max>`.
  - **"full_range"**: Particle is uniformily distributed between the full x range of the simulation.

imp_xstart_val
--------------

Starting x value when :ref:`imp_xstart_opt <imp_xstart_opt>` = "single_value".

imp_xrange_min
--------------

Minimum x value when :ref:`imp_xstart_opt <imp_xstart_opt>` = "range".

imp_xrange_max
--------------

Maximum x value when :ref:`imp_xstart_opt <imp_xstart_opt>` = "range".

imp_ystart_opt
--------------

Option for what y the particle starts.

  - **"single_value"**: Particle starts at a specific y designated by :ref:`imp_ystart_opt <imp_ystart_opt>`.
  - **"range"**: Particle is uniformily distributed between :ref:`imp_yrange_min <imp_yrange_min>` and :ref:`imp_yrange_max <imp_yrange_max>`.
  - **"full_range"**: Particle is uniformily distributed between the full y range of the simulation.

imp_ystart_val
--------------

Starting y value when :ref:`imp_ystart_opt <imp_ystart_opt>` = "single_value".

imp_yrange_min
--------------

Minimum y value when :ref:`imp_ystart_opt <imp_ystart_opt>` = "range".

imp_yrange_max
--------------

Maximum y value when :ref:`imp_ystart_opt <imp_ystart_opt>` = "range".

imp_zstart_opt
--------------

Option for what z the particle starts.

  - **"single_value"**: Particle starts at a specific z designated by :ref:`imp_zstart_opt <imp_zstart_opt>`.
  - **"range"**: Particle is uniformily distributed between :ref:`imp_zrange_min <imp_zrange_min>` and :ref:`imp_zrange_max <imp_zrange_max>`.
  - **"full_range"**: Particle is uniformily distributed between the full z range of the simulation.

imp_zstart_val
--------------

Starting z value when :ref:`imp_zstart_opt <imp_zstart_opt>` = "single_value".

imp_zrange_min
--------------

Minimum z value when :ref:`imp_zstart_opt <imp_zstart_opt>` = "range".

imp_zrange_max
--------------

Maximum z value when :ref:`imp_zstart_opt <imp_zstart_opt>` = "range".

imp_collisions
--------------

Toggle for if collisions with the background plasma should be included. The collision model is an implementation of the Nanbu collision model for cumulative small angle collisions. It is a very general collision model that should be applicable to all regions of the plasma. The only approximation within the model is assuming that the impurity ions are in thermal equlibrium with the ions but *only when calculating the Coloumb logarithm*, which is likely has a negligible impact. 

Nanbu, K. Theory of cumulative small-angle collisions in plasmas. Phys. Rev. E 55, 4642–4652 (1997).

  - **"off"**: Collision model turned off.
  - **"on"**: Collision model turned on.

imp_time_step_opt
-----------------

(Slated for removal) Leave as "constant". 

imp_time_step
-------------

Time step for the simulation in seconds. The time spent following impurities is inversely proportional to the time step. The correct time step depends on the simulation, but it is unlikely one would need to go below :math:`10^{-10}` seconds.

imp_time_step_min
-----------------

(Slated for removal) Leave as :math:`10^{-12}` seconds.

imp_source_scale_fact
---------------------

The meaning of this variable depends on what type of situation the simulation represents, but it is used to convert the Monte Carlo weight density in each cell to actual density in (:math:`m^{-3}`). Internally, every density value is just multiplied by this value. One could equally run a simulation with this set to 1.0 and multiply the impurity density by this value in post-processing if they wanted to.

  - Steady-state simulations: This value represents a particle source rate in units of particles/s. Generally used when :ref:`imp_tstart_opt <imp_tstart_opt>` = "full_range".
  - Time-dependent: This value represents the total number of particles in the real world. For example, if the source represents the instantaneous source from a pellet, one would calculate the number of atoms in the pellet and set this option equal to that. Generally used when :ref:`imp_tstart_opt <imp_tstart_opt>` = "single_value" or "range".

imp_vel_stats
-------------

Toggle to include the average Cartesian impurity velocity components in each cell. This exists just because each component is a 4D array that can take up substantial memory and isn't always needed. If you aren't memory-bound, you can just set this to "on" and forget about it.

  - **"off"**: Do not track average Cartesian velocity components. 
  - **"on"**: The average Cartesian velocity components are tracked and saved in the NetCDF file as imp_vX, imp_vY and imp_vZ. 

imp_iz_recomb
-------------

Toggle to turn on/off impurity ionization and recombination (for whatever reason, if you want to do that you can). Ionization/recombination probabilities are calculated within the code by loading the corresponding rate coefficients at a given ne, Te and converting it to a probability of ionizing/recombining via ``prob = rate_coef * ne * imp_time_step``. If :ref:`imp_time_step <imp_time_step>` is too large, probabilities above 1.0 could occur and the results should be treated with caution (a warning will be output if this happens). The solution is to just use a smaller time step. 

  - **"off"**: Particles remain at their initial charge state determined by :ref:`imp_init_charge <imp_init_charge>`.
  - **"on"**: Particles are free to ionize and recombine according to probabilities calculated from ADAS.

print_interval
--------------

A number telling Flan how often to print how many impurities have been followed. 10 would be every 10%, 100 would be every 1%, etc.

var_red_split
-------------

Turn on the particle splitting variance reduction scheme. If a particle is deemed in a "high importance" (low-count) region, the particle is split into two. The weight of the split particles are assigned proportional to the probability of ionizing or recombing, whichever is greater. For example, say a particle's intial weight is 0.4, and it's probability of ionizing is 0.25 and recombining is 0.1. If it is in a high-importance region, the original particle's weight will be reduced by 0.4*0.25=0.1 and "given" to the split particle, which will also be one charge state higher (because ionizing had a greater probability). The end result is a 0.3 weight particle and a 0.1 weight particle with one higher charge state. This process repeats until the starting particle weight is below :ref:`var_red_min_weight <var_red_min_weight>`.

  - **"off"**: Variance reduction is turned off.
  - **"iz_rec"**: Described above. :ref:`imp_iz_recomb <imp_iz_recomb>` must be turned on for this to work (an error will be issued if not). 
  - **"coll"**: (Not implemented yet)

var_red_import
--------------

This option determines what counts as a high-importance" region. Right now only median is implemented. It only matters when :ref:`var_red_split <var_red_split>` or :ref:`var_red_rusrol <var_red_rusrol>` is in use.

  - **"median"**: High-importance region is anywhere with counts less than :ref:`var_red_med_mod <var_red_med_mod>` * (median counts). Median counts is just the median number of counts across every cell in the simulation volume for each frame.
  - **"exp_dist"**: (Not implemented yet) 
  - **"exp_time"**: (Not implemented yet)

var_red_freq
------------

The frequency at which to update the high-importance regions. What this means is how often you want Flan to recalculate the criteria in :ref:`var_red_import <var_red_import>`. For instance, 10 recalculates it every 10% of particles, 100 would do it every 1% of particles, etc.

var_red_min_weight
------------------

Minimum allowed weight of a particle below which variance reduction is not performed. This is needed, otherwise the particle would continue to be split forever.

var_red_med_mod
---------------

Scalar that the median number of counts for each frame is multiplied by to determine the threshold for what is considered a high-importance (low count) region.

var_red_rusrol
--------------

Switch to turn on Russian roulette variance reduction. Russian roulette variance reduction encourages a simulation to spend more computational time in high-importance regions by "killing" off particles in low importance regions. If a particle is in a low-importance region, it will have a probability of :ref:`var_red_rusrol_prob <var_red_rusrol_prob>` of being killed. If it survives, then its weight is increased by weight / :ref:`var_red_rusrol_prob <var_red_rusrol_prob>`. This is a traditional Monte Carlo scheme, and the implementation in Flan is rather straightforward. 

The criteria for a high-importance region is determined with :ref:`var_red_import <var_red_import>`.

  - **"off"**: Russian roulette is turned off.
  - **"on"**: Russian roulette is turned on.

var_red_rusrol_prob
-------------------

Probability used in Russian roulette scheme when :ref:`var_red_rusrol <var_red_rusrol>` is "on". 

openadas_root
-------------

Full path to the directory containing ADAS data. In this directory should be subdirectories of the acd and scd files. E.g., if :ref:`openadas_year <openadas_year>` = 89, then the directories ``acd89`` and ``scd89`` should exist within openadas_root. 

openadas_year
-------------

The year of the ADAS data. Specifically, the scd and acd data from the ADF11 data type. The files should be saved in the standard naming convention. E.g., if simulating tungsten and using ADAS data from year 50, then the files ``[openadas_root]/acd50/acd50_w.dat`` and ``[openadas_root]/scd50/scd50_w.dat`` should exist.
