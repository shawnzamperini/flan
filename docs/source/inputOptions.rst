======================================================================================================
Input options
======================================================================================================

This page covers the names of various variables within the code, their definition/allowed values, and finally the default values that each variable takes on. 


Merged Simulation Options Table
===============================

+--------------------------+-----------------------------+
| Tag                      |        Description          |
+==========================+=============================+
| :ref:`A1 <A1>`           |                             |
+--------------------------+-----------------------------+
| :ref:`A2 <A2>`           |                             |
+--------------------------+-----------------------------+
| :ref:`A3 <A3>`           |                             |
+--------------------------+-----------------------------+
| :ref:`A4 <A4>`           |                             |
+--------------------------+-----------------------------+
| :ref:`A5 <A5>`           |                             |
+--------------------------+-----------------------------+
| :ref:`A6 <A6>`           |                             |
+--------------------------+-----------------------------+
| :ref:`A7 <A7>`           |                             |
+--------------------------+-----------------------------+
| :ref:`A8 <A8>`           |                             |
+--------------------------+-----------------------------+
| :ref:`A9 <A9>`           |                             |
+--------------------------+-----------------------------+
| :ref:`A10 <A10>`         |                             |
+--------------------------+-----------------------------+
| :ref:`A11 <A11>`         |                             |
+--------------------------+-----------------------------+
| :ref:`A12 <A12>`         |                             |
+--------------------------+-----------------------------+
| :ref:`A13 <A13>`         |                             |
+--------------------------+-----------------------------+
| :ref:`A14 <A14>`         |                             |
+--------------------------+-----------------------------+
| :ref:`A15 <A15>`         |                             |
+--------------------------+-----------------------------+
| :ref:`A16 <A16>`         |                             |
+--------------------------+-----------------------------+
| :ref:`A17 <A17>`         |                             |
+--------------------------+-----------------------------+
| :ref:`A18 <A18>`         |                             |
+--------------------------+-----------------------------+
| :ref:`A19 <A19>`         |                             |
+--------------------------+-----------------------------+
| :ref:`A20 <A20>`         |                             |
+--------------------------+-----------------------------+
| :ref:`A21 <A21>`         |                             |
+--------------------------+-----------------------------+
| :ref:`A22 <A22>`         |                             |
+--------------------------+-----------------------------+
| :ref:`A23 <A23>`         |                             |
+--------------------------+-----------------------------+
| :ref:`A24 <A24>`         |                             |
+--------------------------+-----------------------------+
| :ref:`A25 <A25>`         |                             |
+--------------------------+-----------------------------+
| :ref:`A26 <A26>`         |                             |
+--------------------------+-----------------------------+
| :ref:`A27 <A27>`         |                             |
+--------------------------+-----------------------------+
| :ref:`A28 <A28>`         |                             |
+--------------------------+-----------------------------+
| :ref:`A29 <A29>`         |                             |
+--------------------------+-----------------------------+
| :ref:`A30 <A30>`         |                             |
+--------------------------+-----------------------------+
| :ref:`A31 <A31>`         |                             |
+--------------------------+-----------------------------+
| :ref:`A32 <A32>`         |                             |
+--------------------------+-----------------------------+
| :ref:`A33 <A33>`         |                             |
+--------------------------+-----------------------------+
| :ref:`A34 <A34>`         |                             |
+--------------------------+-----------------------------+
| :ref:`A35 <A35>`         |                             |
+--------------------------+-----------------------------+
| :ref:`A36 <A36>`         |                             |
+--------------------------+-----------------------------+
| :ref:`A37 <A37>`         |                             |
+--------------------------+-----------------------------+
| :ref:`A38 <A38>`         |                             |
+--------------------------+-----------------------------+
| :ref:`A39 <A39>`         |                             |
+--------------------------+-----------------------------+
| :ref:`A40 <A40>`         |                             |
+--------------------------+-----------------------------+
| :ref:`A41 <A41>`         |                             |
+--------------------------+-----------------------------+
| :ref:`A42 <A42>`         |                             |
+--------------------------+-----------------------------+
| :ref:`A43 <A43>`         |                             |
+--------------------------+-----------------------------+
| :ref:`A44 <A44>`         |                             |
+--------------------------+-----------------------------+
| :ref:`A45 <A45>`         |                             |
+--------------------------+-----------------------------+



.. _A1:

**A1**  
  Describe std::string m_gkyl_dir (default: "undefined")

.. _A2:

**A2**  
  Describe std::string m_gkyl_casename (default: "undefined")

.. _A3:

**A3**  
  Describe int m_gkyl_frame_start (default: 0)

.. _A4:

**A4**  
  Describe int m_gkyl_frame_end (default: 1)

.. _A5:

**A5**  
  Describe std::string m_gkyl_elec_name (default: "elc")

.. _A6:

**A6**  
  Describe std::string m_gkyl_ion_name (default: "ion")

.. _A7:

**A7**  
  Describe double m_gkyl_elec_mass_amu (default: 0.000548)

.. _A8:

**A8**  
  Describe double m_gkyl_ion_mass_amu (default: 2.014)

.. _A9:

**A9**  
  Describe std::string m_gkyl_file_type (default: "binary")

.. _A10:

**A10**  
  Describe double m_lcfs_x (default: 0.0)

.. _A11:

**A11**  
  Describe double m_imp_xbound_buffer (default: 0.0)

.. _A12:

**A12**  
  Describe std::string m_min_xbound_type (default: "absorbing")

.. _A13:

**A13**  
  Describe int m_imp_atom_num (default: 74)

.. _A14:

**A14**  
  Describe double m_imp_mass_amu (default: 183.84)

.. _A15:

**A15**  
  Describe int m_imp_init_charge (default: 1)

.. _A16:

**A16**  
  Describe int m_imp_num (default: 1)

.. _A17:

**A17**  
  Describe double m_imp_xmin (default: 0.0)

.. _A18:

**A18**  
  Describe double m_imp_xmax (default: 0.0)

.. _A19:

**A19**  
  Describe std::string m_imp_ystart_opt (default: "single_value")

.. _A20:

**A20**  
  Describe double m_imp_ystart_val (default: 0.0)

.. _A21:

**A21**  
  Describe std::string m_imp_zstart_opt (default: "single_value")

.. _A22:

**A22**  
  Describe double m_imp_zstart_val (default: 0.0)

.. _A23:

**A23**  
  Describe std::string m_imp_collisions (default: "off")

.. _A24:

**A24**  
  Describe std::string m_imp_time_step_opt (default: "variable")

.. _A25:

**A25**  
  Describe double m_imp_time_step (default: 1e-07)

.. _A26:

**A26**  
  Describe double m_imp_time_step_min (default: 1e-12)

.. _A27:

**A27**  
  Describe double m_imp_source_scale_fact (default: 1.0)

.. _A28:

**A28**  
  Describe std::string m_imp_vel_stats (default: "off")

.. _A29:

**A29**  
  Describe std::string m_imp_iz_recomb (default: "on")

.. _A30:

**A30**  
  Describe int m_print_interval (default: 10)

.. _A31:

**A31**  
  Describe std::string m_var_red (default: "off")

.. _A32:

**A32**  
  Describe std::string m_var_red_mode (default: "median")

.. _A33:

**A33**  
  Describe double m_var_red_freq (default: 0.1)

.. _A34:

**A34**  
  Describe double m_var_red_min_weight (default: 0.1)

.. _A35:

**A35**  
  Describe double m_var_red_med_mod (default: 1.0)

.. _A36:

**A36**  
  Describe std::string m_openadas_root (default: "undefined")

.. _A37:

**A37**  
  Describe int m_openadas_year (default: 50)

.. _A38:

**A38**  
  Describe int m_imp_ystart_opt_int (default: 0)

.. _A39:

**A39**  
  Describe int m_imp_zstart_opt_int (default: 0)

.. _A40:

**A40**  
  Describe int m_imp_collisions_int (default: 0)

.. _A41:

**A41**  
  Describe int m_var_red_int (default: 0)

.. _A42:

**A42**  
  Describe int m_var_red_mode_int (default: 0)

.. _A43:

**A43**  
  Describe int m_imp_time_step_opt_int (default: 0)

.. _A44:

**A44**  
  Describe int m_imp_vel_stats_int (default: 0)

.. _A45:

**A45**  
  Describe int m_imp_iz_recomb_int (default: 0)
