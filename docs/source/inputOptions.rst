======================================================================================================
Input options
======================================================================================================

This page contains all the input options available to the user, their definition/allowed values, and finally the default values that each variable takes on. 

Last Update: 8/22/25 (commit XXX)

.. list-table:: Input Options
  :header-rows: 1
  :widths: 25 50 25

  * - Variable
    - Description
    - Default Value
  * - :ref:`gkyl_casename` <gkyl_casename>
    - The name of the Gkeyll simulation
    - "undefined"


+-----------------------------------------------+--------------------------------------------------------------+---------------------+
| Variable                                      | Description                                                  | Default Value       |
+===============================================+==============================================================+=====================+
| :ref:`gkyl_dir <gkyl_dir>`                    |                                                              | "undefined"         |
+-----------------------------------------------+--------------------------------------------------------------+---------------------+
| :ref:`gkyl_casename <gkyl_casename>`      |                                                              | "undefined"         |
+-----------------------------------------------+--------------------------------------------------------------+---------------------+
| :ref:`gkyl_frame_start <gkyl_frame_start>`|                                                              | 0                   |
+-----------------------------------------------+--------------------------------------------------------------+---------------------+
| :ref:`gkyl_frame_end <gkyl_frame_end>`    |                                                              | 1                   |
+-----------------------------------------------+--------------------------------------------------------------+---------------------+
| :ref:`gkyl_elec_name <gkyl_elec_name>`    |                                                              | "elc"               |
+-----------------------------------------------+--------------------------------------------------------------+---------------------+
| :ref:`gkyl_ion_name <gkyl_ion_name>`      |                                                              | "ion"               |
+-----------------------------------------------+--------------------------------------------------------------+---------------------+
| :ref:`gkyl_elec_mass_amu <A1>`              |                                                              | 0.000548            |
+-----------------------------------------------+--------------------------------------------------------------+---------------------+
| :ref:`gkyl_ion_mass_amu <A2>`               |                                                              | 2.014               |
+-----------------------------------------------+--------------------------------------------------------------+---------------------+
| :ref:`gkyl_file_type <gkyl_file_type>`    |                                                              | "binary"            |
+-----------------------------------------------+--------------------------------------------------------------+---------------------+
| :ref:`lcfs_x <lcfs_x>`                    |                                                              | 0.0                 |
+-----------------------------------------------+--------------------------------------------------------------+---------------------+
| :ref:`imp_xbound_buffer <A3>`               |                                                              | 0.0                 |
+-----------------------------------------------+--------------------------------------------------------------+---------------------+
| :ref:`min_xbound_type <min_xbound_type>`  |                                                              | "absorbing"         |
+-----------------------------------------------+--------------------------------------------------------------+---------------------+
| :ref:`imp_atonum <imp_atonum>`        |                                                              | 74                  |
+-----------------------------------------------+--------------------------------------------------------------+---------------------+
| :ref:`imp_mass_amu <imp_mass_amu>`        |                                                              | 183.84              |
+-----------------------------------------------+--------------------------------------------------------------+---------------------+
| :ref:`imp_init_charge <imp_init_charge>`  |                                                              | 1                   |
+-----------------------------------------------+--------------------------------------------------------------+---------------------+
| :ref:`imp_num <imp_num>`                  |                                                              | 1                   |
+-----------------------------------------------+--------------------------------------------------------------+---------------------+
| :ref:`imp_xmin <imp_xmin>`                |                                                              | 0.0                 |
+-----------------------------------------------+--------------------------------------------------------------+---------------------+
| :ref:`imp_xmax <imp_xmax>`                |                                                              | 0.0                 |
+-----------------------------------------------+--------------------------------------------------------------+---------------------+
| :ref:`imp_ystart_opt <imp_ystart_opt>`    |                                                              | "single_value"      |
+-----------------------------------------------+--------------------------------------------------------------+---------------------+
| :ref:`imp_ystart_val <imp_ystart_val>`    |                                                              | 0.0                 |
+-----------------------------------------------+--------------------------------------------------------------+---------------------+
| :ref:`imp_zstart_opt <imp_zstart_opt>`    |                                                              | "single_value"      |
+-----------------------------------------------+--------------------------------------------------------------+---------------------+
| :ref:`imp_zstart_val <imp_zstart_val>`    |                                                              | 0.0                 |
+-----------------------------------------------+--------------------------------------------------------------+---------------------+
| :ref:`imp_collisions <imp_collisions>`    |                                                              | "off"               |
+-----------------------------------------------+--------------------------------------------------------------+---------------------+
| :ref:`imp_time_step_opt <A4>`               |                                                              | "variable"          |
+-----------------------------------------------+--------------------------------------------------------------+---------------------+
| :ref:`imp_time_step <imp_time_step>`      |                                                              | 1e-07               |
+-----------------------------------------------+--------------------------------------------------------------+---------------------+
| :ref:`imp_time_step_min <A5>`               |                                                              | 1e-12               |
+-----------------------------------------------+--------------------------------------------------------------+---------------------+
| :ref:`imp_source_scale_fact <A6>`           |                                                              | 1.0                 |
+-----------------------------------------------+--------------------------------------------------------------+---------------------+
| :ref:`imp_vel_stats <imp_vel_stats>`      |                                                              | "off"               |
+-----------------------------------------------+--------------------------------------------------------------+---------------------+
| :ref:`imp_iz_recomb <imp_iz_recomb>`      |                                                              | "on"                |
+-----------------------------------------------+--------------------------------------------------------------+---------------------+
| :ref:`print_interval <print_interval>`    |                                                              | 10                  |
+-----------------------------------------------+--------------------------------------------------------------+---------------------+
| :ref:`var_red <var_red>`                  |                                                              | "off"               |
+-----------------------------------------------+--------------------------------------------------------------+---------------------+
| :ref:`var_red_mode <var_red_mode>`        |                                                              | "median"            |
+-----------------------------------------------+--------------------------------------------------------------+---------------------+
| :ref:`var_red_freq <var_red_freq>`        |                                                              | 0.1                 |
+-----------------------------------------------+--------------------------------------------------------------+---------------------+
| :ref:`var_red_min_weight <A7>`              |                                                              | 0.1                 |
+-----------------------------------------------+--------------------------------------------------------------+---------------------+
| :ref:`var_red_med_mod <var_red_med_mod>`  |                                                              | 1.0                 |
+-----------------------------------------------+--------------------------------------------------------------+---------------------+
| :ref:`openadas_root <openadas_root>`      |                                                              | "undefined"         |
+-----------------------------------------------+--------------------------------------------------------------+---------------------+
| :ref:`openadas_year <openadas_year>`      |                                                              | 50                  |
+-----------------------------------------------+--------------------------------------------------------------+---------------------+




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
