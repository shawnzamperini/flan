======================================================================================================
Input options
======================================================================================================

This page covers the names of various variables within the code, their definition/allowed values, and finally the default values that each variable takes on. 

+-----------------------------------------------+--------------------------------------------------------------+---------------------+
| Variable                                      | Description                                                  | Default Value       |
+===============================================+==============================================================+=====================+
| :ref:`m_gkyl_dir <m_gkyl_dir>`                |                                                              | "undefined"         |
+-----------------------------------------------+--------------------------------------------------------------+---------------------+
| :ref:`m_gkyl_casename <m_gkyl_casename>`      |                                                              | "undefined"         |
+-----------------------------------------------+--------------------------------------------------------------+---------------------+
| :ref:`m_gkyl_frame_start <m_gkyl_frame_start>`|                                                              | 0                   |
+-----------------------------------------------+--------------------------------------------------------------+---------------------+
| :ref:`m_gkyl_frame_end <m_gkyl_frame_end>`    |                                                              | 1                   |
+-----------------------------------------------+--------------------------------------------------------------+---------------------+
| :ref:`m_gkyl_elec_name <m_gkyl_elec_name>`    |                                                              | "elc"               |
+-----------------------------------------------+--------------------------------------------------------------+---------------------+
| :ref:`m_gkyl_ion_name <m_gkyl_ion_name>`      |                                                              | "ion"               |
+-----------------------------------------------+--------------------------------------------------------------+---------------------+
| :ref:`m_gkyl_elec_mass_amu <A1>`              |                                                              | 0.000548            |
+-----------------------------------------------+--------------------------------------------------------------+---------------------+
| :ref:`m_gkyl_ion_mass_amu <A2>`               |                                                              | 2.014               |
+-----------------------------------------------+--------------------------------------------------------------+---------------------+
| :ref:`m_gkyl_file_type <m_gkyl_file_type>`    |                                                              | "binary"            |
+-----------------------------------------------+--------------------------------------------------------------+---------------------+
| :ref:`m_lcfs_x <m_lcfs_x>`                    |                                                              | 0.0                 |
+-----------------------------------------------+--------------------------------------------------------------+---------------------+
| :ref:`m_imp_xbound_buffer <A3>`               |                                                              | 0.0                 |
+-----------------------------------------------+--------------------------------------------------------------+---------------------+
| :ref:`m_min_xbound_type <m_min_xbound_type>`  |                                                              | "absorbing"         |
+-----------------------------------------------+--------------------------------------------------------------+---------------------+
| :ref:`m_imp_atom_num <m_imp_atom_num>`        |                                                              | 74                  |
+-----------------------------------------------+--------------------------------------------------------------+---------------------+
| :ref:`m_imp_mass_amu <m_imp_mass_amu>`        |                                                              | 183.84              |
+-----------------------------------------------+--------------------------------------------------------------+---------------------+
| :ref:`m_imp_init_charge <m_imp_init_charge>`  |                                                              | 1                   |
+-----------------------------------------------+--------------------------------------------------------------+---------------------+
| :ref:`m_imp_num <m_imp_num>`                  |                                                              | 1                   |
+-----------------------------------------------+--------------------------------------------------------------+---------------------+
| :ref:`m_imp_xmin <m_imp_xmin>`                |                                                              | 0.0                 |
+-----------------------------------------------+--------------------------------------------------------------+---------------------+
| :ref:`m_imp_xmax <m_imp_xmax>`                |                                                              | 0.0                 |
+-----------------------------------------------+--------------------------------------------------------------+---------------------+
| :ref:`m_imp_ystart_opt <m_imp_ystart_opt>`    |                                                              | "single_value"      |
+-----------------------------------------------+--------------------------------------------------------------+---------------------+
| :ref:`m_imp_ystart_val <m_imp_ystart_val>`    |                                                              | 0.0                 |
+-----------------------------------------------+--------------------------------------------------------------+---------------------+
| :ref:`m_imp_zstart_opt <m_imp_zstart_opt>`    |                                                              | "single_value"      |
+-----------------------------------------------+--------------------------------------------------------------+---------------------+
| :ref:`m_imp_zstart_val <m_imp_zstart_val>`    |                                                              | 0.0                 |
+-----------------------------------------------+--------------------------------------------------------------+---------------------+
| :ref:`m_imp_collisions <m_imp_collisions>`    |                                                              | "off"               |
+-----------------------------------------------+--------------------------------------------------------------+---------------------+
| :ref:`m_imp_time_step_opt <A4>`               |                                                              | "variable"          |
+-----------------------------------------------+--------------------------------------------------------------+---------------------+
| :ref:`m_imp_time_step <m_imp_time_step>`      |                                                              | 1e-07               |
+-----------------------------------------------+--------------------------------------------------------------+---------------------+
| :ref:`m_imp_time_step_min <A5>`               |                                                              | 1e-12               |
+-----------------------------------------------+--------------------------------------------------------------+---------------------+
| :ref:`m_imp_source_scale_fact <A6>`           |                                                              | 1.0                 |
+-----------------------------------------------+--------------------------------------------------------------+---------------------+
| :ref:`m_imp_vel_stats <m_imp_vel_stats>`      |                                                              | "off"               |
+-----------------------------------------------+--------------------------------------------------------------+---------------------+
| :ref:`m_imp_iz_recomb <m_imp_iz_recomb>`      |                                                              | "on"                |
+-----------------------------------------------+--------------------------------------------------------------+---------------------+
| :ref:`m_print_interval <m_print_interval>`    |                                                              | 10                  |
+-----------------------------------------------+--------------------------------------------------------------+---------------------+
| :ref:`m_var_red <m_var_red>`                  |                                                              | "off"               |
+-----------------------------------------------+--------------------------------------------------------------+---------------------+
| :ref:`m_var_red_mode <m_var_red_mode>`        |                                                              | "median"            |
+-----------------------------------------------+--------------------------------------------------------------+---------------------+
| :ref:`m_var_red_freq <m_var_red_freq>`        |                                                              | 0.1                 |
+-----------------------------------------------+--------------------------------------------------------------+---------------------+
| :ref:`m_var_red_min_weight <A7>`              |                                                              | 0.1                 |
+-----------------------------------------------+--------------------------------------------------------------+---------------------+
| :ref:`m_var_red_med_mod <m_var_red_med_mod>`  |                                                              | 1.0                 |
+-----------------------------------------------+--------------------------------------------------------------+---------------------+
| :ref:`m_openadas_root <m_openadas_root>`      |                                                              | "undefined"         |
+-----------------------------------------------+--------------------------------------------------------------+---------------------+
| :ref:`m_openadas_year <m_openadas_year>`      |                                                              | 50                  |
+-----------------------------------------------+--------------------------------------------------------------+---------------------+




.. _m_gkyl_dir:

**m_gkyl_dir**  
  Describe std::string m_gkyl_dir (default: "undefined")

.. _m_gkyl_casename:

**m_gkyl_casename**  
  Describe std::string m_gkyl_casename (default: "undefined")

.. _m_gkyl_frame_start:

**m_gkyl_frame_start**  
  Describe int m_gkyl_frame_start (default: 0)

.. _m_gkyl_frame_end:

**m_gkyl_frame_end**  
  Describe int m_gkyl_frame_end (default: 1)

.. _m_gkyl_elec_name:

**m_gkyl_elec_name**  
  Describe std::string m_gkyl_elec_name (default: "elc")

.. _m_gkyl_ion_name:

**m_gkyl_ion_name**  
  Describe std::string m_gkyl_ion_name (default: "ion")

.. A1:

**m_gkyl_elec_mass_amu**  
  Describe double m_gkyl_elec_mass_amu (default: 0.000548)

.. A2:

**m_gkyl_ion_mass_amu**  
  Describe double m_gkyl_ion_mass_amu (default: 2.014)

.. _m_gkyl_file_type:

**m_gkyl_file_type**  
  Describe std::string m_gkyl_file_type (default: "binary")

.. _m_lcfs_x:

**m_lcfs_x**  
  Describe double m_lcfs_x (default: 0.0)

.. A3:

**m_imp_xbound_buffer**  
  Describe double m_imp_xbound_buffer (default: 0.0)

.. _m_min_xbound_type:

**m_min_xbound_type**  
  Describe std::string m_min_xbound_type (default: "absorbing")

.. _m_imp_atom_num:

**m_imp_atom_num**  
  Describe int m_imp_atom_num (default: 74)

.. _m_imp_mass_amu:

**m_imp_mass_amu**  
  Describe double m_imp_mass_amu (default: 183.84)

.. _m_imp_init_charge:

**m_imp_init_charge**  
  Describe int m_imp_init_charge (default: 1)

.. _m_imp_num:

**m_imp_num**  
  Describe int m_imp_num (default: 1)

.. _m_imp_xmin:

**m_imp_xmin**  
  Describe double m_imp_xmin (default: 0.0)

.. _m_imp_xmax:

**m_imp_xmax**  
  Describe double m_imp_xmax (default: 0.0)

.. _m_imp_ystart_opt:

**m_imp_ystart_opt**  
  Describe std::string m_imp_ystart_opt (default: "single_value")

.. _m_imp_ystart_val:

**m_imp_ystart_val**  
  Describe double m_imp_ystart_val (default: 0.0)

.. _m_imp_zstart_opt:

**m_imp_zstart_opt**  
  Describe std::string m_imp_zstart_opt (default: "single_value")

.. _m_imp_zstart_val:

**m_imp_zstart_val**  
  Describe double m_imp_zstart_val (default: 0.0)

.. _m_imp_collisions:

**m_imp_collisions**  
  Describe std::string m_imp_collisions (default: "off")

.. A4:

**m_imp_time_step_opt**  
  Describe std::string m_imp_time_step_opt (default: "variable")

.. _m_imp_time_step:

**m_imp_time_step**  
  Describe double m_imp_time_step (default: 1e-07)

.. A5:

**_m_imp_time_step_min**  
  Describe double m_imp_time_step_min (default: 1e-12)

.. A6:

**m_imp_source_scale_fact**  
  Describe double m_imp_source_scale_fact (default: 1.0)

.. _m_imp_vel_stats:

**m_imp_vel_stats**  
  Describe std::string m_imp_vel_stats (default: "off")

.. _m_imp_iz_recomb:

**m_imp_iz_recomb**  
  Describe std::string m_imp_iz_recomb (default: "on")

.. _m_print_interval:

**m_print_interval**  
  Describe int m_print_interval (default: 10)

.. _m_var_red:

**m_var_red**  
  Describe std::string m_var_red (default: "off")

.. _m_var_red_mode:

**m_var_red_mode**  
  Describe std::string m_var_red_mode (default: "median")

.. _m_var_red_freq:

**m_var_red_freq**  
  Describe double m_var_red_freq (default: 0.1)

.. A7:

**m_var_red_min_weight**  
  Describe double m_var_red_min_weight (default: 0.1)

.. _m_var_red_med_mod:

**m_var_red_med_mod**  
  Describe double m_var_red_med_mod (default: 1.0)

.. _m_openadas_root:

**m_openadas_root**  
  Describe std::string m_openadas_root (default: "undefined")

.. _m_openadas_year:

**m_openadas_year**  
  Describe int m_openadas_year (default: 50)
