/**
* @file options.cpp
* @brief Function definitoins for Options class.
*/
#include <iostream>
#include <string>
#include <tuple>
#include <initializer_list>
#include <sstream>

#include "errors.h"
#include "options.h"
#include "options_device.h"


namespace Options
{

	/**
	* Constructs an Options object with every option at its default value.
	*
	* The defaults are set where each member variable is declared in
	* options.h. The control integers are then set to match the default
	* string options by calling initialize_control_ints.
	*/
	Options::Options()
	{
		// Set internal control integers based on (at this point, default)
		// options. They'll get changed later if the option changes.
		initialize_control_ints();
	}
	
	/**
	* Sets every control integer to match its string option.
	*
	* Control integers are integer versions of string options that are
	* faster
	* to check in the code than comparing strings. The logic that assigns
	* each one lives in the option's setter, so this calls every such setter
	* with the option's current value. Any new option with a control integer
	* must have its setter called here too. Otherwise its control integer
	* stays at -1, which get_control_int reports as an error.
	*/
	void Options::initialize_control_ints()
	{
		// We are setting each value to the value it is already at, why?
		// Because the logic for assigning the respective control integers
		// lives within each set_ function. This is not the best, since if
		// one were to add an additional control integer for a new option and
		// then forgot to call the corresponding set_ function within this
		// function, then they may introduce a silent bug! This can be silent
		// because if the user does not have that option in their input file,
		// the default string value will be used but there is no guarantee that
		// the default control integer will be the correct value! Luckily,
		// this bug is avoided by accessing the control integers via
		// get_control_int, which will tell you if a control integer hasn't
		// been correctly initalized (by not including it here).
		set_use_gpu(m_use_gpu);
		set_bkg_source(m_bkg_source);
		set_save_track(m_save_track);
		set_imp_tstart_opt(m_imp_tstart_opt);
		set_imp_xstart_opt(m_imp_xstart_opt);
		set_imp_ystart_opt(m_imp_ystart_opt);
		set_imp_zstart_opt(m_imp_zstart_opt);
		set_imp_temp_start_opt(m_imp_temp_start_opt);
		set_imp_collisions(m_imp_collisions);
		set_imp_thermal_force(m_imp_thermal_force);
		set_var_red_split(m_var_red_split);
		set_var_red_import(m_var_red_import);
		set_var_red_rusrol(m_var_red_rusrol);
		set_imp_time_step_opt(m_imp_time_step_opt);
		set_imp_iz_recomb(m_imp_iz_recomb);
		set_tbound_type(m_tbound_type);
		set_min_xbound_type(m_min_xbound_type);
		set_max_xbound_type(m_max_xbound_type);
		set_min_ybound_type(m_min_ybound_type);
		set_max_ybound_type(m_max_ybound_type);
		set_min_zbound_type(m_min_zbound_type);
		set_max_zbound_type(m_max_zbound_type);
		set_calc_grad_elec(m_calc_grad_elec);
		set_s3x_chi_absorb(m_s3x_chi_absorb);
		set_hermes_format(m_hermes_format);
		set_hermes_heat_flux(m_hermes_heat_flux);
		set_hermes_field_aligned(m_hermes_field_aligned);
		set_hermes_inductive_epar(m_hermes_inductive_epar);
		set_imp_source_opt(m_imp_source_opt);
		set_sput_flux_opt(m_sput_flux_opt);
		set_sput_bohm_floor(m_sput_bohm_floor);
		set_sheath_model(m_sheath_model);
		set_sput_self(m_sput_self);
		set_sput_impact_opt(m_sput_impact_opt);
	}

	// Setter definitions. They're all the same, so just throwing together
	// as a mass of code. Can organize if necessary, but it would probably
	// just take up three times as much space.

	/**
	* Sets the case_name option.
	*
	* Name of the simulation case. The results are saved to <case_name>.nc.
	*
	* @param case_name New value of the case_name option.
	*/
	void Options::set_case_name(std::string case_name) 
		{m_case_name = case_name;}

	/**
	* Sets the use_gpu option.
	*
	* Whether to run the simulation on GPUs.
	*
	* Valid values are "off" and "cuda". An invalid value prints an error
	* and leaves the option unchanged. Also sets the control integer
	* returned by use_gpu_int() to 0 for "off", 1 for "cuda". An invalid
	* value leaves the control integer unchanged.
	*
	* @param use_gpu New value of the use_gpu option.
	*/
	void Options::set_use_gpu(std::string use_gpu) 
	{
		if (check_input<std::string>("use_gpu", use_gpu,
			{"off", "cuda"}))
		{
			m_use_gpu = use_gpu;
		}

		// Assign control integers
		if (use_gpu == "off") m_use_gpu_int = 0;
		else if (use_gpu == "cuda") m_use_gpu_int = 1;
	}

	/**
	* Sets the seed option.
	*
	* Base seed for the random number generators. A negative value seeds
	* them from the clock so each run is different. The MPI rank is added to
	* the seed so each rank gets a different sequence.
	*
	* @param seed New value of the seed option.
	*/
	void Options::set_seed(int seed) 
		{m_seed = seed;}

	/**
	* Sets the slot_cap option.
	*
	* Maximum number of particle slots, i.e., the most particles followed at
	* once. The number of slots used is the smaller of this and imp_num.
	*
	* @param slot_cap New value of the slot_cap option.
	*/
	void Options::set_slot_cap(int slot_cap) 
		{m_slot_cap = slot_cap;}

	/**
	* Sets the bkg_source option.
	*
	* Source of the background plasma: the slab background (a box of
	* plasma in simple fields, set entirely by the slab_* options), a
	* Gkeyll simulation, a SOLEDGE3X simulation or a HERMES (Hermes-3 or
	* JAX-Hermes3) simulation.
	*
	* Valid values are "slab", "gkeyll", "soledge3x" and "hermes". An
	* invalid value prints an error and leaves the option unchanged. Also
	* sets the control integer returned by bkg_source_int() to 0 for
	* "slab", 1 for "gkeyll", 2 for "soledge3x", 3 for "hermes". An invalid
	* value leaves the control integer unchanged.
	*
	* @param bkg_source New value of the bkg_source option.
	*/
	void Options::set_bkg_source(std::string bkg_source) 
	{
		if (check_input<std::string>("bkg_source", bkg_source,
			{"slab", "gkeyll", "soledge3x", "hermes"}))
		{
			m_bkg_source = bkg_source;
		}

		// Assign control integers
		if (bkg_source == "gkeyll") m_bkg_source_int = 1;
		else if (bkg_source == "soledge3x") m_bkg_source_int = 2;
		else if (bkg_source == "hermes") m_bkg_source_int = 3;
		else if (bkg_source == "slab") m_bkg_source_int = 0;
	}

	/**
	* Sets the save_track option.
	*
	* Whether to save particle tracks to the output file.
	*
	* Valid values are "off" and "on". An invalid value prints an error and
	* leaves the option unchanged. Also sets the control integer returned by
	* save_track_int() to 0 for "off", 1 for "on". An invalid value leaves
	* the control integer unchanged.
	*
	* @param save_track New value of the save_track option.
	*/
	void Options::set_save_track(std::string save_track) 
	{
		if (check_input<std::string>("save_track", save_track,
			{"off", "on"}))
		{
			m_save_track = save_track;
		}

		// Assign control integers
		if (save_track == "off") m_save_track_int = 0;
		else if (save_track == "on") m_save_track_int = 1;
	}
	
	/**
	* Sets the gkyl_dir option.
	*
	* Directory containing the Gkeyll simulation files.
	*
	* @param gkyl_dir New value of the gkyl_dir option.
	*/
	void Options::set_gkyl_dir(std::string gkyl_dir) 
		{m_gkyl_dir = gkyl_dir;}
	
	/**
	* Sets the gkyl_casename option.
	*
	* Name of the Gkeyll case, i.e., the prefix of its output files.
	*
	* @param gkyl_casename New value of the gkyl_casename option.
	*/
	void Options::set_gkyl_casename(std::string gkyl_casename) 
		{m_gkyl_casename = gkyl_casename;}
	
	/**
	* Sets the gkyl_frame_start option.
	*
	* First Gkeyll frame to load.
	*
	* @param gkyl_frame_start New value of the gkyl_frame_start option.
	*/
	void Options::set_gkyl_frame_start(int gkyl_frame_start) 
		{m_gkyl_frame_start = gkyl_frame_start;}
	
	/**
	* Sets the gkyl_frame_end option.
	*
	* Last Gkeyll frame to load. The range from gkyl_frame_start to
	* gkyl_frame_end is inclusive.
	*
	* @param gkyl_frame_end New value of the gkyl_frame_end option.
	*/
	void Options::set_gkyl_frame_end(int gkyl_frame_end) 
		{m_gkyl_frame_end = gkyl_frame_end;}

	/**
	* Sets the gkyl_elec_name option.
	*
	* Name of the electron species in the Gkeyll simulation.
	*
	* @param gkyl_elec_name New value of the gkyl_elec_name option.
	*/
	void Options::set_gkyl_elec_name(std::string gkyl_elec_name) 
		{m_gkyl_elec_name = gkyl_elec_name;}
	
	/**
	* Sets the gkyl_ion_name option.
	*
	* Name of the main ion species in the Gkeyll simulation.
	*
	* @param gkyl_ion_name New value of the gkyl_ion_name option.
	*/
	void Options::set_gkyl_ion_name(std::string gkyl_ion_name) 
		{m_gkyl_ion_name = gkyl_ion_name;}
	
	/**
	* Sets the gkyl_elec_mass_amu option.
	*
	* Electron mass used in the Gkeyll simulation (amu).
	*
	* @param gkyl_elec_mass_amu New value of the gkyl_elec_mass_amu option.
	*/
	void Options::set_gkyl_elec_mass_amu(double gkyl_elec_mass_amu) 
		{m_gkyl_elec_mass_amu = gkyl_elec_mass_amu;}
	
	/**
	* Sets the gkyl_ion_mass_amu option.
	*
	* Main ion mass used in the Gkeyll simulation (amu).
	*
	* @param gkyl_ion_mass_amu New value of the gkyl_ion_mass_amu option.
	*/
	void Options::set_gkyl_ion_mass_amu(double gkyl_ion_mass_amu) 
		{m_gkyl_ion_mass_amu = gkyl_ion_mass_amu;}

	/**
	* Sets the gkyl_file_type option.
	*
	* File format of the Gkeyll output.
	*
	* Valid values are "binary". An invalid value prints an error and leaves
	* the option unchanged.
	*
	* @param gkyl_file_type New value of the gkyl_file_type option.
	*/
	void Options::set_gkyl_file_type(std::string gkyl_file_type) 
	{
		if (check_input<std::string>("gkyl_file_type", gkyl_file_type,
			{"binary"}))
		{
			m_gkyl_file_type = gkyl_file_type;
		}
	}

	/**
	* Sets the gkyl_moment_type option.
	*
	* Type of moment files in the Gkeyll output. Collisions need
	* "bimaxwellian", which provides the ion parallel flow and perpendicular
	* velocity used to calculate the background flow.
	*
	* Valid values are "bimaxwellian", "maxwellian" and "m0m1m2". An invalid
	* value prints an error and leaves the option unchanged.
	*
	* @param gkyl_moment_type New value of the gkyl_moment_type option.
	*/
	void Options::set_gkyl_moment_type(std::string gkyl_moment_type) 
	{
		if (check_input<std::string>("gkyl_moment_type", gkyl_moment_type,
			{"bimaxwellian", "maxwellian", "m0m1m2"}))
		{
			m_gkyl_moment_type = gkyl_moment_type;
		}
	}

	/**
	* Sets the calc_grad_elec option.
	*
	* Whether to calculate the electric field gradient. This is not
	* implemented yet, so turning it on only prints a warning.
	*
	* Valid values are "off" and "on". An invalid value prints an error and
	* leaves the option unchanged. Also sets the control integer returned by
	* calc_grad_elec_int() to 0 for "off", 1 for "on". An invalid value
	* leaves the control integer unchanged.
	*
	* @param calc_grad_elec New value of the calc_grad_elec option.
	*/
	void Options::set_calc_grad_elec(std::string calc_grad_elec) 
		{
			// Only assign if valid option, leaving as default if not
			if (check_input<std::string>("calc_grad_elec", calc_grad_elec, 
				{"off", "on"}))
			{
				m_calc_grad_elec = calc_grad_elec;
			}

			// Assign control integers
			if (calc_grad_elec == "off") m_calc_grad_elec_int = 0;
			else if (calc_grad_elec == "on") m_calc_grad_elec_int = 1;
		}
	
	/**
	* Sets the min_ne option.
	*
	* Lower limit on the electron density (m^-3).
	*
	* @param min_ne New value of the min_ne option.
	*/
	void Options::set_min_ne(double min_ne) 
		{m_min_ne = min_ne;}

	/**
	* Sets the min_te option.
	*
	* Lower limit on the electron temperature (eV).
	*
	* @param min_te New value of the min_te option.
	*/
	void Options::set_min_te(double min_te) 
		{m_min_te = min_te;}

	/**
	* Sets the min_ti option.
	*
	* Lower limit on the ion temperature (eV), applied when calculating
	* collisions with the background ions.
	*
	* @param min_ti New value of the min_ti option.
	*/
	void Options::set_min_ti(double min_ti) 
		{m_min_ti = min_ti;}

	/**
	* Sets the slab_ne option.
	*
	* Electron density of the slab background, uniform in space and time
	* (m^-3). The ion density is the same. Only used with bkg_source =
	* "slab".
	*
	* @param slab_ne New value of the slab_ne option.
	*/
	void Options::set_slab_ne(double slab_ne)
		{m_slab_ne = slab_ne;}

	/**
	* Sets the slab_te option.
	*
	* Electron temperature of the slab background, uniform in space and time
	* (eV). Only used with bkg_source = "slab".
	*
	* @param slab_te New value of the slab_te option.
	*/
	void Options::set_slab_te(double slab_te)
		{m_slab_te = slab_te;}

	/**
	* Sets the slab_ti option.
	*
	* Ion temperature of the slab background, uniform in space and time
	* (eV). Only used with bkg_source = "slab".
	*
	* @param slab_ti New value of the slab_ti option.
	*/
	void Options::set_slab_ti(double slab_ti)
		{m_slab_ti = slab_ti;}

	/**
	* Sets the slab_flow option.
	*
	* Flow velocity of the slab background's plasma along the magnetic
	* field, which is in +X (m/s). Impurities are dragged toward it by
	* friction. Only used with bkg_source = "slab".
	*
	* @param slab_flow New value of the slab_flow option.
	*/
	void Options::set_slab_flow(double slab_flow)
		{m_slab_flow = slab_flow;}

	/**
	* Sets the slab_heat_flux_i option.
	*
	* Parallel conductive heat flux of the slab background's ions, in units
	* of n T_i v_th,i with v_th,i = sqrt(T_i / m_i). Positive is along the
	* magnetic field (+X). This is the size of the heat-flux skew (eps in
	* Collisions::skew_eps) that gives the ion thermal force; about 0.1 is a
	* mild skew. Only used with bkg_source = "slab".
	*
	* @param slab_heat_flux_i New value of the slab_heat_flux_i option.
	*/
	void Options::set_slab_heat_flux_i(double slab_heat_flux_i)
		{m_slab_heat_flux_i = slab_heat_flux_i;}

	/**
	* Sets the slab_heat_flux_e option.
	*
	* Parallel conductive heat flux of the slab background's electrons, in
	* units of n T_e v_th,e with v_th,e = sqrt(T_e / m_e). Positive is along
	* the magnetic field (+X). This is the size of the heat-flux skew that
	* gives the electron thermal force. Only used with bkg_source = "slab".
	*
	* @param slab_heat_flux_e New value of the slab_heat_flux_e option.
	*/
	void Options::set_slab_heat_flux_e(double slab_heat_flux_e)
		{m_slab_heat_flux_e = slab_heat_flux_e;}

	/**
	* Sets the slab_bmag option.
	*
	* Strength of the slab background's magnetic field, which is uniform and
	* in +X (T). Only used with bkg_source = "slab".
	*
	* @param slab_bmag New value of the slab_bmag option.
	*/
	void Options::set_slab_bmag(double slab_bmag)
		{m_slab_bmag = slab_bmag;}

	/**
	* Sets the slab_tmax option.
	*
	* End of the slab background's time range, which starts at 0 (s). The
	* output has 128 time frames evenly spaced over it. Only used with
	* bkg_source = "slab".
	*
	* @param slab_tmax New value of the slab_tmax option.
	*/
	void Options::set_slab_tmax(double slab_tmax)
		{m_slab_tmax = slab_tmax;}

	/**
	* Sets the slab_efield_y option.
	*
	* Electric field of the slab background in +Y at t = 0 (V/m), uniform in
	* space. It is across the magnetic field (+X), so it gives an ExB drift
	* in -Z of E / B. Only for a straight slab (slab_major_radius = 0). Only
	* used with bkg_source = "slab".
	*
	* @param slab_efield_y New value of the slab_efield_y option.
	*/
	void Options::set_slab_efield_y(double slab_efield_y)
		{m_slab_efield_y = slab_efield_y;}

	/**
	* Sets the slab_efield_y_rate option.
	*
	* Rate at which the slab background's electric field in +Y changes
	* (V/m/s): E_Y = slab_efield_y + slab_efield_y_rate * t. A changing
	* field gives a polarization drift along it. Only for a straight slab.
	* Only used with bkg_source = "slab".
	*
	* @param slab_efield_y_rate New value of the slab_efield_y_rate option.
	*/
	void Options::set_slab_efield_y_rate(double slab_efield_y_rate)
		{m_slab_efield_y_rate = slab_efield_y_rate;}

	/**
	* Sets the slab_bgrad_y option.
	*
	* Relative gradient of the slab background's magnetic field strength in
	* +Y (1/m): B_X = slab_bmag * (1 + slab_bgrad_y * Y). It is across the
	* field, so it gives a grad-B drift. Don't make the y boundaries
	* periodic with it. Only for a straight slab. Only used with bkg_source
	* = "slab".
	*
	* @param slab_bgrad_y New value of the slab_bgrad_y option.
	*/
	void Options::set_slab_bgrad_y(double slab_bgrad_y)
		{m_slab_bgrad_y = slab_bgrad_y;}

	/**
	* Sets the slab_major_radius option.
	*
	* If above 0, the slab background's box is bent into a quarter of a
	* torus of this major radius (m), to give the magnetic field curvature:
	* x becomes the major radius R (from 2.5 cm inside to 2.5 cm outside
	* this value), y the height Z and z the toroidal angle (0 to pi/2). The
	* field is then toroidal with strength slab_bmag * slab_major_radius /
	* R, and the flow and heat fluxes are along it. 0 is the straight slab,
	* with the field along +X. Only used with bkg_source = "slab".
	*
	* @param slab_major_radius New value of the slab_major_radius option.
	*/
	void Options::set_slab_major_radius(double slab_major_radius)
		{m_slab_major_radius = slab_major_radius;}

	/**
	* Sets the slab_wall option.
	*
	* Which face of the slab background's box is a material wall, such as
	* a divertor target: particles that reach it are absorbed and counted
	* as hitting the wall, and with imp_source_opt = "sputter" the plasma
	* sputters it. The face replaces whatever its *bound_type gives. Only
	* for a straight slab (slab_major_radius = 0).
	*
	* Valid values are "none", "xmin", "xmax", "ymin" and "ymax". An
	* invalid value prints an error and leaves the option unchanged.
	*
	* @param slab_wall New value of the slab_wall option.
	*/
	void Options::set_slab_wall(std::string slab_wall)
	{
		if (check_input<std::string>("slab_wall", slab_wall,
			{"none", "xmin", "xmax", "ymin", "ymax"}))
		{
			m_slab_wall = slab_wall;
		}
	}

	/**
	* Sets the slab_b_angle option.
	*
	* Angle (degrees) the slab background's magnetic field is turned from
	* +X toward -Y, in the X-Y plane: the unit vector along B is (cos a,
	* -sin a, 0). With the wall at the y min face (slab_wall = "ymin"), it
	* is the angle between the field and the wall, and the field runs into
	* the wall going along +B. 0 (the default) is the field along +X. Must
	* be from -90 to 90. Only for a straight slab, without slab_efield_y,
	* slab_efield_y_rate or slab_bgrad_y.
	*
	* @param slab_b_angle New value of the slab_b_angle option.
	*/
	void Options::set_slab_b_angle(double slab_b_angle)
		{m_slab_b_angle = slab_b_angle;}

	/**
	* Sets the tbound_type option.
	*
	* Boundary condition at the start and end of the background time range.
	*
	* Valid values are "absorbing" and "periodic". An invalid value prints
	* an error and leaves the option unchanged. Also sets the control
	* integer returned by tbound_type_int() to 0 for "absorbing", 1 for
	* "periodic". An invalid value leaves the control integer unchanged.
	*
	* @param tbound_type New value of the tbound_type option.
	*/
	void Options::set_tbound_type(std::string tbound_type) 
	{
		// Only assign if valid option, leaving as default if not
		if (check_input<std::string>("tbound_type", tbound_type,
			{"absorbing", "periodic"}))
		{
			m_tbound_type = tbound_type;
		}

		// Assign control integers
		if (tbound_type == "absorbing") m_tbound_type_int = 0;
		else if (tbound_type == "periodic") m_tbound_type_int = 1;
	}
	/**
	* Sets the s3x_dir option.
	*
	* Directory containing the SOLEDGE3X mesh and plasma files.
	*
	* @param s3x_dir New value of the s3x_dir option.
	*/
	void Options::set_s3x_dir(std::string s3x_dir) 
		{m_s3x_dir = s3x_dir;}

	/**
	* Sets the s3x_mesh_file option.
	*
	* Name of the SOLEDGE3X mesh file in s3x_dir.
	*
	* @param s3x_mesh_file New value of the s3x_mesh_file option.
	*/
	void Options::set_s3x_mesh_file(std::string s3x_mesh_file) 
		{m_s3x_mesh_file = s3x_mesh_file;}

	/**
	* Sets the s3x_frame_start option.
	*
	* First SOLEDGE3X plasma frame to load, i.e., the number in
	* plasma_%05d.h5.
	*
	* @param s3x_frame_start New value of the s3x_frame_start option.
	*/
	void Options::set_s3x_frame_start(int s3x_frame_start) 
		{m_s3x_frame_start = s3x_frame_start;}

	/**
	* Sets the s3x_frame_end option.
	*
	* Last SOLEDGE3X plasma frame to load. The range from
	* s3x_frame_start to s3x_frame_end is inclusive, and must have at least
	* two frames.
	*
	* @param s3x_frame_end New value of the s3x_frame_end option.
	*/
	void Options::set_s3x_frame_end(int s3x_frame_end) 
		{m_s3x_frame_end = s3x_frame_end;}

	/**
	* Sets the s3x_ref_file option.
	*
	* Name of the SOLEDGE3X file in s3x_dir with the reference values the
	* output is normalized by (n0, T0, c0, tau0, phi0, ...).
	*
	* @param s3x_ref_file New value of the s3x_ref_file option.
	*/
	void Options::set_s3x_ref_file(std::string s3x_ref_file) 
		{m_s3x_ref_file = s3x_ref_file;}

	/**
	* Sets the s3x_chi_absorb option.
	*
	* Whether particles are absorbed when they enter a cell inside the
	* SOLEDGE3X wall (chi = 1).
	*
	* Valid values are "off" and "on". An invalid value prints an error and
	* leaves the option unchanged. Also sets the control integer returned by
	* s3x_chi_absorb_int() to 0 for "off", 1 for "on". An invalid value
	* leaves the control integer unchanged.
	*
	* @param s3x_chi_absorb New value of the s3x_chi_absorb option.
	*/
	void Options::set_s3x_chi_absorb(std::string s3x_chi_absorb) 
	{
		// Only assign if valid option, leaving as default if not
		if (check_input<std::string>("s3x_chi_absorb", s3x_chi_absorb, 
			{"off", "on"}))
		{
			m_s3x_chi_absorb = s3x_chi_absorb;
		}

		// Assign control integers
		if (s3x_chi_absorb == "off") m_s3x_chi_absorb_int = 0;
		else if (s3x_chi_absorb == "on") m_s3x_chi_absorb_int = 1;
	}

	/**
	* Sets the hermes_format option.
	*
	* Output format of the HERMES run: "jax" for JAX-Hermes3 (a hypnotoad grid
	* plus one .npz snapshot per time) or "bout" for BOUT++ Hermes-3 dump
	* files (one per MPI process). Also sets the control integer returned by
	* hermes_format_int() to 0 for "jax", 1 for "bout".
	*
	* @param hermes_format New value of the hermes_format option.
	*/
	void Options::set_hermes_format(std::string hermes_format)
	{
		if (check_input<std::string>("hermes_format", hermes_format,
			{"jax", "bout"}))
		{
			m_hermes_format = hermes_format;
		}

		// Assign control integers
		if (hermes_format == "jax") m_hermes_format_int = 0;
		else if (hermes_format == "bout") m_hermes_format_int = 1;
	}

	/**
	* Sets the hermes_grid_file option.
	*
	* Path to the hypnotoad grid file (netCDF) of the HERMES run.
	*
	* @param hermes_grid_file New value of the hermes_grid_file option.
	*/
	void Options::set_hermes_grid_file(std::string hermes_grid_file)
		{m_hermes_grid_file = hermes_grid_file;}

	/**
	* Sets the hermes_snap_dir option.
	*
	* Directory with the HERMES snapshots (snap*.npz for JAX-Hermes3), or
	* the BOUT++ run directory with BOUT.dmp.*.nc (BOUT++ Hermes-3).
	*
	* @param hermes_snap_dir New value of the hermes_snap_dir option.
	*/
	void Options::set_hermes_snap_dir(std::string hermes_snap_dir)
		{m_hermes_snap_dir = hermes_snap_dir;}

	/**
	* Sets the hermes_ion_species option.
	*
	* Name of the main ion species in the Hermes-3 input, which names the
	* ions' variables in BOUT++ dump files (P<name>, V<name>, NV<name>).
	* Only used with hermes_format = "bout".
	*
	* @param hermes_ion_species New value of the hermes_ion_species option.
	*/
	void Options::set_hermes_ion_species(std::string hermes_ion_species)
		{m_hermes_ion_species = hermes_ion_species;}

	/**
	* Sets the hermes_snap_start option.
	*
	* First snapshot to load, as an index into the snapshots in
	* hermes_snap_dir sorted by name (time), from 0. Negative values count
	* from the end, like hermes_snap_end.
	*
	* @param hermes_snap_start New value of the hermes_snap_start option.
	*/
	void Options::set_hermes_snap_start(int hermes_snap_start)
		{m_hermes_snap_start = hermes_snap_start;}

	/**
	* Sets the hermes_snap_end option.
	*
	* Last snapshot to load (inclusive), as an index like hermes_snap_start.
	* Negative counts from the end (-1 is the last).
	*
	* @param hermes_snap_end New value of the hermes_snap_end option.
	*/
	void Options::set_hermes_snap_end(int hermes_snap_end)
		{m_hermes_snap_end = hermes_snap_end;}

	/**
	* Sets the hermes_snap_stride option.
	*
	* Load every hermes_snap_stride-th snapshot from hermes_snap_start to
	* hermes_snap_end.
	*
	* @param hermes_snap_stride New value of the hermes_snap_stride option.
	*/
	void Options::set_hermes_snap_stride(int hermes_snap_stride)
		{m_hermes_snap_stride = hermes_snap_stride;}

	/**
	* Sets the hermes_nnorm option.
	*
	* Density (m^-3) the HERMES output is normalized by (Nnorm).
	*
	* @param hermes_nnorm New value of the hermes_nnorm option.
	*/
	void Options::set_hermes_nnorm(double hermes_nnorm)
		{m_hermes_nnorm = hermes_nnorm;}

	/**
	* Sets the hermes_tnorm option.
	*
	* Temperature (eV) the HERMES output is normalized by (Tnorm).
	*
	* @param hermes_tnorm New value of the hermes_tnorm option.
	*/
	void Options::set_hermes_tnorm(double hermes_tnorm)
		{m_hermes_tnorm = hermes_tnorm;}

	/**
	* Sets the hermes_bnorm option.
	*
	* Magnetic field (T) the HERMES output is normalized by (Bnorm), which
	* sets the time unit 1 / Omega_ci = m_p / (e Bnorm).
	*
	* @param hermes_bnorm New value of the hermes_bnorm option.
	*/
	void Options::set_hermes_bnorm(double hermes_bnorm)
		{m_hermes_bnorm = hermes_bnorm;}

	/**
	* Sets the hermes_ion_mass_amu option.
	*
	* Mass number (AA) of the HERMES main ions, used to turn NVi into the
	* parallel flow and as the background ion mass.
	*
	* @param hermes_ion_mass_amu New value of the hermes_ion_mass_amu option.
	*/
	void Options::set_hermes_ion_mass_amu(double hermes_ion_mass_amu)
		{m_hermes_ion_mass_amu = hermes_ion_mass_amu;}

	/**
	* Sets the hermes_zperiod option.
	*
	* The HERMES domain covers 1 / hermes_zperiod of the torus (zperiod in the
	* HERMES input).
	*
	* @param hermes_zperiod New value of the hermes_zperiod option.
	*/
	void Options::set_hermes_zperiod(int hermes_zperiod)
		{m_hermes_zperiod = hermes_zperiod;}

	/**
	* Sets the hermes_heat_flux option.
	*
	* Where a HERMES background's parallel heat fluxes (for the thermal
	* force) come from: "conduction" calculates them as -kappa_par grad_par
	* T with Braginskii's kappa_par, as Hermes-3 does, and "off" leaves them
	* out (no thermal force). Also sets the control integer returned by
	* hermes_heat_flux_int() to 0 for "off", 1 for "conduction".
	*
	* @param hermes_heat_flux New value of the hermes_heat_flux option.
	*/
	void Options::set_hermes_heat_flux(std::string hermes_heat_flux)
	{
		if (check_input<std::string>("hermes_heat_flux", hermes_heat_flux,
			{"conduction", "off"}))
		{
			m_hermes_heat_flux = hermes_heat_flux;
		}

		// Assign control integers
		if (hermes_heat_flux == "off") m_hermes_heat_flux_int = 0;
		else if (hermes_heat_flux == "conduction") m_hermes_heat_flux_int = 1;
	}

	/**
	* Sets the hermes_inductive_epar option.
	*
	* "on" (the default) adds the inductive parallel electric field
	* -dApar/dt to E_par whenever the run has Apar, i.e. a BOUT++ Hermes-3
	* background from an electromagnetic run (Apar is read at the output
	* times on either side of each loaded one). Runs without Apar
	* (electrostatic, or JAX-Hermes3) go on without it. "off" leaves
	* E = -grad phi. Also sets the control integer
	* returned by hermes_inductive_epar_int() to 0 for "off", 1 for "on".
	*
	* @param hermes_inductive_epar New value of the hermes_inductive_epar
	*        option.
	*/
	void Options::set_hermes_inductive_epar(std::string hermes_inductive_epar)
	{
		if (check_input<std::string>("hermes_inductive_epar",
			hermes_inductive_epar, {"off", "on"}))
		{
			m_hermes_inductive_epar = hermes_inductive_epar;
		}

		// Assign control integers
		if (hermes_inductive_epar == "off") m_hermes_inductive_epar_int = 0;
		else if (hermes_inductive_epar == "on")
			m_hermes_inductive_epar_int = 1;
	}

	/**
	* Sets the hermes_field_aligned option.
	*
	* "on" (the default) interpolates and differentiates a HERMES
	* background in y along the field line, using the grid's zShift, as
	* the shifted metric intends. "off" does it at fixed z like any other
	* background, and takes parallel gradients as b . grad. "off" is only
	* for comparison (e.g., to see what field alignment changes): on a
	* turbulent background it mixes unrelated structures between y planes.
	* Also sets the control integer returned by hermes_field_aligned_int()
	* to 0 for "off", 1 for "on".
	*
	* @param hermes_field_aligned New value of the hermes_field_aligned
	*        option.
	*/
	void Options::set_hermes_field_aligned(std::string hermes_field_aligned)
	{
		if (check_input<std::string>("hermes_field_aligned",
			hermes_field_aligned, {"on", "off"}))
		{
			m_hermes_field_aligned = hermes_field_aligned;
		}

		// Assign control integers
		if (hermes_field_aligned == "off") m_hermes_field_aligned_int = 0;
		else if (hermes_field_aligned == "on") m_hermes_field_aligned_int = 1;
	}

	/**
	* Sets the hermes_kappa_coef_e option.
	*
	* Electron parallel conduction coefficient (kappa_coefficient of the
	* electrons in the HERMES input). The default, 3.16 / sqrt(2), is
	* Hermes-3's.
	*
	* @param hermes_kappa_coef_e New value of the hermes_kappa_coef_e option.
	*/
	void Options::set_hermes_kappa_coef_e(double hermes_kappa_coef_e)
		{m_hermes_kappa_coef_e = hermes_kappa_coef_e;}

	/**
	* Sets the hermes_kappa_coef_i option.
	*
	* Ion parallel conduction coefficient (kappa_coefficient of the ions in
	* the HERMES input). The default, 3.9, is Hermes-3's.
	*
	* @param hermes_kappa_coef_i New value of the hermes_kappa_coef_i option.
	*/
	void Options::set_hermes_kappa_coef_i(double hermes_kappa_coef_i)
		{m_hermes_kappa_coef_i = hermes_kappa_coef_i;}

	/**
	* Sets the hermes_kappa_limit_alpha option.
	*
	* Flux limiter of the parallel conduction (kappa_limit_alpha in the
	* HERMES input). Zero or negative, the default, means no limit.
	*
	* @param hermes_kappa_limit_alpha New value of the
	*        hermes_kappa_limit_alpha option.
	*/
	void Options::set_hermes_kappa_limit_alpha(
		double hermes_kappa_limit_alpha)
		{m_hermes_kappa_limit_alpha = hermes_kappa_limit_alpha;}

	/**
	* Sets the imp_source_opt option.
	*
	* Where impurities come from: "manual" starts them where the
	* imp_*start_opt options say, and "sputter" makes them the atoms the
	* background plasma sputters off the walls of a SOLEDGE3X background or
	* the targets of a HERMES background (see Sputtering::Source).
	*
	* Valid values are "manual" and "sputter". An invalid value prints an
	* error and leaves the option unchanged. Also sets the control integer
	* returned by imp_source_opt_int() to 0 for "manual", 1 for "sputter".
	* An invalid value leaves the control integer unchanged.
	*
	* @param imp_source_opt New value of the imp_source_opt option.
	*/
	void Options::set_imp_source_opt(std::string imp_source_opt)
	{
		if (check_input<std::string>("imp_source_opt", imp_source_opt,
			{"manual", "sputter"}))
		{
			m_imp_source_opt = imp_source_opt;
		}

		if (imp_source_opt == "manual") m_imp_source_opt_int = 0;
		else if (imp_source_opt == "sputter") m_imp_source_opt_int = 1;
	}

	/**
	* Sets the sput_table_file option.
	*
	* Path to the sputtering yield table, made by python/rustbca_yields.py
	* (e.g., data/sputtering/D_on_W.h5). Needed when imp_source_opt is
	* "sputter".
	*
	* @param sput_table_file New value of the sput_table_file option.
	*/
	void Options::set_sput_table_file(std::string sput_table_file)
		{m_sput_table_file = sput_table_file;}

	/**
	* Sets the sput_flux_opt option.
	*
	* Which part of the background ion flow carries ions onto the wall:
	* "parallel_exb" is the parallel flow plus the E x B drift, "parallel"
	* is the parallel flow only.
	*
	* Valid values are "parallel_exb" and "parallel". An invalid value
	* prints an error and leaves the option unchanged. Also sets the
	* control integer returned by sput_flux_opt_int() to 0 for
	* "parallel_exb", 1 for "parallel". An invalid value leaves the control
	* integer unchanged.
	*
	* @param sput_flux_opt New value of the sput_flux_opt option.
	*/
	void Options::set_sput_flux_opt(std::string sput_flux_opt)
	{
		if (check_input<std::string>("sput_flux_opt", sput_flux_opt,
			{"parallel_exb", "parallel"}))
		{
			m_sput_flux_opt = sput_flux_opt;
		}

		if (sput_flux_opt == "parallel_exb") m_sput_flux_opt_int = 0;
		else if (sput_flux_opt == "parallel") m_sput_flux_opt_int = 1;
	}

	/**
	* Sets the sput_bohm_floor option.
	*
	* Whether the parallel flow onto the wall is at least the ion sound
	* speed (the Bohm criterion), for backgrounds whose flow next to the
	* wall is too slow.
	*
	* Valid values are "off" and "on". An invalid value prints an error and
	* leaves the option unchanged. Also sets the control integer returned
	* by sput_bohm_floor_int() to 0 for "off", 1 for "on". An invalid value
	* leaves the control integer unchanged.
	*
	* @param sput_bohm_floor New value of the sput_bohm_floor option.
	*/
	void Options::set_sput_bohm_floor(std::string sput_bohm_floor)
	{
		if (check_input<std::string>("sput_bohm_floor", sput_bohm_floor,
			{"off", "on"}))
		{
			m_sput_bohm_floor = sput_bohm_floor;
		}

		if (sput_bohm_floor == "off") m_sput_bohm_floor_int = 0;
		else if (sput_bohm_floor == "on") m_sput_bohm_floor_int = 1;
	}

	/**
	* Sets the sheath_model option.
	*
	* The electric field of the sheath in front of the wall (the magnetic
	* presheath and the Debye sheath), added to the push in the plasma
	* cells next to the wall surfaces (see sheath.cpp). The layers are far
	* thinner than a cell, so the background doesn't hold them. "brooks" is
	* Brooks' two-exponential profile, as in ERO and GITR, and "borodkina"
	* Borodkina and Tsvetkov's; both split the drop the same way. The
	* default is "brooks". It only does anything on a background with wall
	* surfaces (as for imp_source_opt = "sputter"), and is left out on
	* other backgrounds.
	*
	* Valid values are "off", "brooks" and "borodkina". An invalid value
	* prints an error and leaves the option unchanged. Also sets the
	* control integer returned by sheath_model_int() to 0 for "off", 1 for
	* "brooks", 2 for "borodkina" (the values of Sheath::Model). An invalid
	* value leaves the control integer unchanged.
	*
	* @param sheath_model New value of the sheath_model option.
	*/
	void Options::set_sheath_model(std::string sheath_model)
	{
		if (check_input<std::string>("sheath_model", sheath_model,
			{"off", "brooks", "borodkina"}))
		{
			m_sheath_model = sheath_model;
		}

		if (sheath_model == "off") m_sheath_model_int = 0;
		else if (sheath_model == "brooks") m_sheath_model_int = 1;
		else if (sheath_model == "borodkina") m_sheath_model_int = 2;
	}

	/**
	* Sets the sput_angle_val option.
	*
	* Angle (degrees from the surface normal) at which background ions hit
	* the wall, used to look up the sputtering yield.
	*
	* @param sput_angle_val New value of the sput_angle_val option.
	*/
	void Options::set_sput_angle_val(double sput_angle_val)
		{m_sput_angle_val = sput_angle_val;}

	/**
	* Sets the sput_max_wall_dist option.
	*
	* How far, in cell widths, the nearest segment of the wall outline can
	* be from a wall face before that face falls back to its own normal and
	* area (see Sputtering::Surface). Only used with a wall mask
	* (SOLEDGE3X): HERMES targets always use their faces' own.
	*
	* @param sput_max_wall_dist New value of the sput_max_wall_dist option.
	*/
	void Options::set_sput_max_wall_dist(double sput_max_wall_dist)
		{m_sput_max_wall_dist = sput_max_wall_dist;}


	/**
	* Sets the min_xbound_type option.
	*
	* Boundary condition at the minimum x boundary.
	*
	* Valid values are "absorbing", "periodic" and "core". An invalid value
	* prints an error and leaves the option unchanged. Also sets the control
	* integer returned by min_xbound_type_int() to 0 for "absorbing", 1 for
	* "periodic", 2 for "core". An invalid value leaves the control integer
	* unchanged.
	*
	* @param min_xbound_type New value of the min_xbound_type option.
	*/
	void Options::set_min_xbound_type(std::string min_xbound_type) 
	{
		// Only assign if valid option, leaving as default if not
		if (check_input<std::string>("min_xbound_type", min_xbound_type,
			{"absorbing", "periodic", "core"}))
		{
			m_min_xbound_type = min_xbound_type;
		}

		// Assign control integers
		if (min_xbound_type == "absorbing") m_min_xbound_type_int = 0;
		else if (min_xbound_type == "periodic") m_min_xbound_type_int = 1;
		else if (min_xbound_type == "core") m_min_xbound_type_int = 2;
	}

	/**
	* Sets the max_xbound_type option.
	*
	* Boundary condition at the maximum x boundary.
	*
	* Valid values are "absorbing", "periodic" and "core". An invalid value
	* prints an error and leaves the option unchanged. Also sets the control
	* integer returned by max_xbound_type_int() to 0 for "absorbing", 1 for
	* "periodic", 2 for "core". An invalid value leaves the control integer
	* unchanged.
	*
	* @param max_xbound_type New value of the max_xbound_type option.
	*/
	void Options::set_max_xbound_type(std::string max_xbound_type) 
	{
		// Only assign if valid option, leaving as default if not
		if (check_input<std::string>("max_xbound_type", max_xbound_type,
			{"absorbing", "periodic", "core"}))
		{
			m_max_xbound_type = max_xbound_type;
		}

		// Assign control integers
		if (max_xbound_type == "absorbing") m_max_xbound_type_int = 0;
		else if (max_xbound_type == "periodic") m_max_xbound_type_int = 1;
		else if (max_xbound_type == "core") m_max_xbound_type_int = 2;
	}

	/**
	* Sets the min_ybound_type option.
	*
	* Boundary condition at the minimum y boundary.
	*
	* Valid values are "absorbing", "periodic" and "core". An invalid value
	* prints an error and leaves the option unchanged. Also sets the control
	* integer returned by min_ybound_type_int() to 0 for "absorbing", 1 for
	* "periodic", 2 for "core". An invalid value leaves the control integer
	* unchanged.
	*
	* @param min_ybound_type New value of the min_ybound_type option.
	*/
	void Options::set_min_ybound_type(std::string min_ybound_type) 
	{
		// Only assign if valid option, leaving as default if not
		if (check_input<std::string>("min_ybound_type", min_ybound_type,
			{"absorbing", "periodic", "core"}))
		{
			m_min_ybound_type = min_ybound_type;
		}

		// Assign control integers
		if (min_ybound_type == "absorbing") m_min_ybound_type_int = 0;
		else if (min_ybound_type == "periodic") m_min_ybound_type_int = 1;
		else if (min_ybound_type == "core") m_min_ybound_type_int = 2;
	}

	/**
	* Sets the max_ybound_type option.
	*
	* Boundary condition at the maximum y boundary.
	*
	* Valid values are "absorbing", "periodic" and "core". An invalid value
	* prints an error and leaves the option unchanged. Also sets the control
	* integer returned by max_ybound_type_int() to 0 for "absorbing", 1 for
	* "periodic", 2 for "core". An invalid value leaves the control integer
	* unchanged.
	*
	* @param max_ybound_type New value of the max_ybound_type option.
	*/
	void Options::set_max_ybound_type(std::string max_ybound_type) 
	{
		// Only assign if valid option, leaving as default if not
		if (check_input<std::string>("max_ybound_type", max_ybound_type,
			{"absorbing", "periodic", "core"}))
		{
			m_max_ybound_type = max_ybound_type;
		}

		// Assign control integers
		if (max_ybound_type == "absorbing") m_max_ybound_type_int = 0;
		else if (max_ybound_type == "periodic") m_max_ybound_type_int = 1;
		else if (max_ybound_type == "core") m_max_ybound_type_int = 2;
	}

	/**
	* Sets the min_zbound_type option.
	*
	* Boundary condition at the minimum z boundary.
	*
	* Valid values are "absorbing", "periodic", "core" and "lcfs". "lcfs"
	* is periodic for particles inside the last closed flux surface
	* (x < lcfs_x) and absorbing outside it. An invalid value prints an
	* error and leaves the option unchanged. Also sets the control integer
	* returned by min_zbound_type_int() to 0 for "absorbing", 1 for
	* "periodic", 2 for "core", 3 for "lcfs". An invalid value leaves the
	* control integer unchanged.
	*
	* @param min_zbound_type New value of the min_zbound_type option.
	*/
	void Options::set_min_zbound_type(std::string min_zbound_type) 
	{
		// Only assign if valid option, leaving as default if not
		if (check_input<std::string>("min_zbound_type", min_zbound_type,
			{"absorbing", "periodic", "core", "lcfs"}))
		{
			m_min_zbound_type = min_zbound_type;
		}

		// Assign control integers
		if (min_zbound_type == "absorbing") m_min_zbound_type_int = 0;
		else if (min_zbound_type == "periodic") m_min_zbound_type_int = 1;
		else if (min_zbound_type == "core") m_min_zbound_type_int = 2;
		else if (min_zbound_type == "lcfs") m_min_zbound_type_int = 3;
	}

	/**
	* Sets the max_zbound_type option.
	*
	* Boundary condition at the maximum z boundary.
	*
	* Valid values are "absorbing", "periodic", "core" and "lcfs". "lcfs"
	* is periodic for particles inside the last closed flux surface
	* (x < lcfs_x) and absorbing outside it. An invalid value prints an
	* error and leaves the option unchanged. Also sets the control integer
	* returned by max_zbound_type_int() to 0 for "absorbing", 1 for
	* "periodic", 2 for "core", 3 for "lcfs". An invalid value leaves the
	* control integer unchanged.
	*
	* @param max_zbound_type New value of the max_zbound_type option.
	*/
	void Options::set_max_zbound_type(std::string max_zbound_type) 
	{
		// Only assign if valid option, leaving as default if not
		if (check_input<std::string>("max_zbound_type", max_zbound_type,
			{"absorbing", "periodic", "core", "lcfs"}))
		{
			m_max_zbound_type = max_zbound_type;
		}

		// Assign control integers
		if (max_zbound_type == "absorbing") m_max_zbound_type_int = 0;
		else if (max_zbound_type == "periodic") m_max_zbound_type_int = 1;
		else if (max_zbound_type == "core") m_max_zbound_type_int = 2;
		else if (max_zbound_type == "lcfs") m_max_zbound_type_int = 3;
	}

	/**
	* Sets the lcfs_x option.
	*
	* x coordinate of the last closed flux surface.
	*
	* @param lcfs_x New value of the lcfs_x option.
	*/
	void Options::set_lcfs_x(double lcfs_x) 
		{m_lcfs_x = lcfs_x;}

	/**
	* Sets the lcfs_y option.
	*
	* y coordinate of the last closed flux surface.
	*
	* @param lcfs_y New value of the lcfs_y option.
	*/
	void Options::set_lcfs_y(double lcfs_y) 
		{m_lcfs_y = lcfs_y;}

	/**
	* Sets the lcfs_z option.
	*
	* z coordinate of the last closed flux surface.
	*
	* @param lcfs_z New value of the lcfs_z option.
	*/
	void Options::set_lcfs_z(double lcfs_z) 
		{m_lcfs_z = lcfs_z;}

	/**
	* Sets the sep_x_bc_xp_z1 option.
	*
	* First z value for the separatrix x boundary condition near the
	* X-point.
	*
	* @param sep_x_bc_xp_z1 New value of the sep_x_bc_xp_z1 option.
	*/
	void Options::set_sep_x_bc_xp_z1(double sep_x_bc_xp_z1) 
		{m_sep_x_bc_xp_z1 = sep_x_bc_xp_z1;}

	/**
	* Sets the sep_x_bc_xp_z2 option.
	*
	* Second z value for the separatrix x boundary condition near the
	* X-point.
	*
	* @param sep_x_bc_xp_z2 New value of the sep_x_bc_xp_z2 option.
	*/
	void Options::set_sep_x_bc_xp_z2(double sep_x_bc_xp_z2) 
		{m_sep_x_bc_xp_z2 = sep_x_bc_xp_z2;}

	/**
	* Sets the imp_atom_num option.
	*
	* Atomic number of the impurity.
	*
	* @param imp_atom_num New value of the imp_atom_num option.
	*/
	void Options::set_imp_atom_num(int imp_atom_num) 
		{m_imp_atom_num = imp_atom_num;}

	/**
	* Sets the imp_iz_max_charge option.
	*
	* Highest charge the impurity can ionize to. It is for regression
	* cases that compare with models leaving out further ionization, e.g.,
	* prompt redeposition after ionizing once (Fussmann's model), not for
	* physics runs. Below 0 (the default) there is no limit.
	*
	* @param imp_iz_max_charge New value of the imp_iz_max_charge option.
	*/
	void Options::set_imp_iz_max_charge(int imp_iz_max_charge)
		{m_imp_iz_max_charge = imp_iz_max_charge;}

	/**
	* Sets the imp_mass_amu option.
	*
	* Mass of the impurity (amu).
	*
	* @param imp_mass_amu New value of the imp_mass_amu option.
	*/
	void Options::set_imp_mass_amu(double imp_mass_amu) 
		{m_imp_mass_amu = imp_mass_amu;}
	
	/**
	* Sets the imp_init_charge option.
	*
	* Charge state the impurities start with.
	*
	* @param imp_init_charge New value of the imp_init_charge option.
	*/
	void Options::set_imp_init_charge(int imp_init_charge) 
		{m_imp_init_charge = imp_init_charge;}
	
	/**
	* Sets the imp_num option.
	*
	* Total number of impurity particles to follow.
	*
	* @param imp_num New value of the imp_num option.
	*/
	void Options::set_imp_num(int imp_num) 
		{m_imp_num = imp_num;}
	
	/**
	* Sets the imp_tstart_opt option.
	*
	* How each impurity's starting time is chosen: a single value
	* (imp_tstart_val), uniformly between imp_trange_min and imp_trange_max,
	* or uniformly over the full background time range.
	*
	* Valid values are "single_value", "range" and "full_range". An invalid
	* value prints an error and leaves the option unchanged. Also sets the
	* control integer returned by imp_tstart_opt_int() to 0 for
	* "single_value", 1 for "range", 2 for "full_range". An invalid value
	* leaves the control integer unchanged.
	*
	* @param imp_tstart_opt New value of the imp_tstart_opt option.
	*/
	void Options::set_imp_tstart_opt(std::string imp_tstart_opt) 
		{
			// Only assign if valid option, leaving as default if not
			if (check_input<std::string>("imp_tstart_opt", imp_tstart_opt, 
				{"single_value", "range", "full_range"}))
			{
				m_imp_tstart_opt = imp_tstart_opt;
			}

			// Assign control integers
			if (imp_tstart_opt == "single_value") m_imp_tstart_opt_int = 0;
			else if (imp_tstart_opt == "range") m_imp_tstart_opt_int = 1;
			else if (imp_tstart_opt == "full_range") m_imp_tstart_opt_int = 2;
		}

	/**
	* Sets the imp_tstart_val option.
	*
	* Starting time of every impurity when imp_tstart_opt is "single_value"
	* (s).
	*
	* @param imp_tstart_val New value of the imp_tstart_val option.
	*/
	void Options::set_imp_tstart_val(double imp_tstart_val) 
		{m_imp_tstart_val = imp_tstart_val;}
		
	/**
	* Sets the imp_trange_min option.
	*
	* Lower limit of the starting time range when imp_tstart_opt is "range"
	* (s).
	*
	* @param imp_trange_min New value of the imp_trange_min option.
	*/
	void Options::set_imp_trange_min(double imp_trange_min) 
		{m_imp_trange_min = imp_trange_min;}
	
	/**
	* Sets the imp_trange_max option.
	*
	* Upper limit of the starting time range when imp_tstart_opt is "range"
	* (s).
	*
	* @param imp_trange_max New value of the imp_trange_max option.
	*/
	void Options::set_imp_trange_max(double imp_trange_max) 
		{m_imp_trange_max = imp_trange_max;}

	/**
	* Sets the imp_xstart_opt option.
	*
	* How each impurity's starting x coordinate is chosen: a single value
	* (imp_xstart_val), uniformly between imp_xrange_min and imp_xrange_max,
	* or uniformly over the full x range of the background.
	*
	* Valid values are "single_value", "range" and "full_range". An invalid
	* value prints an error and leaves the option unchanged. Also sets the
	* control integer returned by imp_xstart_opt_int() to 0 for
	* "single_value", 1 for "range", 2 for "full_range". An invalid value
	* leaves the control integer unchanged.
	*
	* @param imp_xstart_opt New value of the imp_xstart_opt option.
	*/
	void Options::set_imp_xstart_opt(std::string imp_xstart_opt) 
		{
			// Only assign if valid option, leaving as default if not
			if (check_input<std::string>("imp_xstart_opt", imp_xstart_opt, 
				{"single_value", "range", "full_range"}))
			{
				m_imp_xstart_opt = imp_xstart_opt;
			}

			// Assign control integers
			if (imp_xstart_opt == "single_value") m_imp_xstart_opt_int = 0;
			else if (imp_xstart_opt == "range") m_imp_xstart_opt_int = 1;
			else if (imp_xstart_opt == "full_range") m_imp_xstart_opt_int = 2;
		}

	/**
	* Sets the imp_xstart_val option.
	*
	* Starting x coordinate of every impurity when imp_xstart_opt is
	* "single_value".
	*
	* @param imp_xstart_val New value of the imp_xstart_val option.
	*/
	void Options::set_imp_xstart_val(double imp_xstart_val) 
		{m_imp_xstart_val = imp_xstart_val;}
		
	/**
	* Sets the imp_xrange_min option.
	*
	* Lower limit of the starting x range when imp_xstart_opt is "range".
	*
	* @param imp_xrange_min New value of the imp_xrange_min option.
	*/
	void Options::set_imp_xrange_min(double imp_xrange_min) 
		{m_imp_xrange_min = imp_xrange_min;}
	
	/**
	* Sets the imp_xrange_max option.
	*
	* Upper limit of the starting x range when imp_xstart_opt is "range".
	*
	* @param imp_xrange_max New value of the imp_xrange_max option.
	*/
	void Options::set_imp_xrange_max(double imp_xrange_max) 
		{m_imp_xrange_max = imp_xrange_max;}
	
	/**
	* Sets the imp_ystart_opt option.
	*
	* How each impurity's starting y coordinate is chosen: a single value
	* (imp_ystart_val), uniformly between imp_yrange_min and imp_yrange_max,
	* or uniformly over the full y range of the background.
	*
	* Valid values are "single_value", "range" and "full_range". An invalid
	* value prints an error and leaves the option unchanged. Also sets the
	* control integer returned by imp_ystart_opt_int() to 0 for
	* "single_value", 1 for "range", 2 for "full_range". An invalid value
	* leaves the control integer unchanged.
	*
	* @param imp_ystart_opt New value of the imp_ystart_opt option.
	*/
	void Options::set_imp_ystart_opt(std::string imp_ystart_opt) 
		{
			// Only assign if valid option, leaving as default if not
			if (check_input<std::string>("imp_ystart_opt", imp_ystart_opt, 
				{"single_value", "range", "full_range"}))
			{
				m_imp_ystart_opt = imp_ystart_opt;
			}

			// Assign control integers
			if (imp_ystart_opt == "single_value") m_imp_ystart_opt_int = 0;
			else if (imp_ystart_opt == "range") m_imp_ystart_opt_int = 1;
			else if (imp_ystart_opt == "full_range") m_imp_ystart_opt_int = 2;
		}
	
	/**
	* Sets the imp_ystart_val option.
	*
	* Starting y coordinate of every impurity when imp_ystart_opt is
	* "single_value".
	*
	* @param imp_ystart_val New value of the imp_ystart_val option.
	*/
	void Options::set_imp_ystart_val(double imp_ystart_val) 
		{m_imp_ystart_val = imp_ystart_val;}

	/**
	* Sets the imp_yrange_min option.
	*
	* Lower limit of the starting y range when imp_ystart_opt is "range".
	*
	* @param imp_yrange_min New value of the imp_yrange_min option.
	*/
	void Options::set_imp_yrange_min(double imp_yrange_min) 
		{m_imp_yrange_min = imp_yrange_min;}
	
	/**
	* Sets the imp_yrange_max option.
	*
	* Upper limit of the starting y range when imp_ystart_opt is "range".
	*
	* @param imp_yrange_max New value of the imp_yrange_max option.
	*/
	void Options::set_imp_yrange_max(double imp_yrange_max) 
		{m_imp_yrange_max = imp_yrange_max;}
	
	/**
	* Sets the imp_zstart_opt option.
	*
	* How each impurity's starting z coordinate is chosen: a single value
	* (imp_zstart_val), uniformly between imp_zrange_min and imp_zrange_max,
	* or uniformly over the full z range of the background.
	*
	* Valid values are "single_value", "range" and "full_range". An invalid
	* value prints an error and leaves the option unchanged. Also sets the
	* control integer returned by imp_zstart_opt_int() to 0 for
	* "single_value", 1 for "range", 2 for "full_range". An invalid value
	* leaves the control integer unchanged.
	*
	* @param imp_zstart_opt New value of the imp_zstart_opt option.
	*/
	void Options::set_imp_zstart_opt(std::string imp_zstart_opt) 
		{
			// Only assign if valid option, leaving as default if not
			if (check_input<std::string>("imp_zstart_opt", imp_zstart_opt, 
				{"single_value", "range", "full_range"}))
			{
				m_imp_zstart_opt = imp_zstart_opt;
			}

			// Assign control integers
			if (imp_zstart_opt == "single_value") m_imp_zstart_opt_int = 0;
			else if (imp_zstart_opt == "range") m_imp_zstart_opt_int = 1;
			else if (imp_zstart_opt == "full_range") m_imp_zstart_opt_int = 2;
		}
	
	/**
	* Sets the imp_zstart_val option.
	*
	* Starting z coordinate of every impurity when imp_zstart_opt is
	* "single_value".
	*
	* @param imp_zstart_val New value of the imp_zstart_val option.
	*/
	void Options::set_imp_zstart_val(double imp_zstart_val) 
		{m_imp_zstart_val = imp_zstart_val;}

	/**
	* Sets the imp_zone_start_val option.
	*
	* Zone every impurity starts in, numbered from 1 like the SOLEDGE3X
	* zones. The starting x, y and z are in that zone's coordinates. Only
	* used for backgrounds with more than one zone.
	*
	* @param imp_zone_start_val New value of the imp_zone_start_val option.
	*/
	void Options::set_imp_zone_start_val(int imp_zone_start_val) 
		{m_imp_zone_start_val = imp_zone_start_val;}

	/**
	* Sets the imp_zrange_min option.
	*
	* Lower limit of the starting z range when imp_zstart_opt is "range".
	*
	* @param imp_zrange_min New value of the imp_zrange_min option.
	*/
	void Options::set_imp_zrange_min(double imp_zrange_min) 
		{m_imp_zrange_min = imp_zrange_min;}
	
	/**
	* Sets the imp_zrange_max option.
	*
	* Upper limit of the starting z range when imp_zstart_opt is "range".
	*
	* @param imp_zrange_max New value of the imp_zrange_max option.
	*/
	void Options::set_imp_zrange_max(double imp_zrange_max) 
		{m_imp_zrange_max = imp_zrange_max;}
	
	/**
	* Sets the imp_temp_start_opt option.
	*
	* How each impurity's starting temperature is chosen: a single value
	* (imp_temp_start_val) or the local main ion temperature. Starting
	* velocities are sampled from an isotropic Maxwellian at this
	* temperature.
	*
	* Valid values are "single_value" and "main_ion". An invalid value
	* prints an error and leaves the option unchanged. Also sets the control
	* integer returned by imp_temp_start_opt_int() to 0 for "single_value",
	* 1 for "main_ion". An invalid value leaves the control integer
	* unchanged.
	*
	* @param imp_temp_start_opt New value of the imp_temp_start_opt option.
	*/
	void Options::set_imp_temp_start_opt(std::string imp_temp_start_opt) 
	{
		// Only assign if valid option, leaving as default if not
		if (check_input<std::string>("imp_temp_start_opt", 
			imp_temp_start_opt, {"single_value", "main_ion"}))
		{
			m_imp_temp_start_opt = imp_temp_start_opt;
		}

		// Assign control integers
		if (imp_temp_start_opt == "single_value") 
			m_imp_temp_start_opt_int = 0;
		else if (imp_temp_start_opt == "main_ion") 
			m_imp_temp_start_opt_int = 1;
	}
	
	/**
	* Sets the imp_temp_start_val option.
	*
	* Starting temperature of every impurity when imp_temp_start_opt is
	* "single_value" (eV).
	*
	* @param imp_temp_start_val New value of the imp_temp_start_val option.
	*/
	void Options::set_imp_temp_start_val(double imp_temp_start_val) 
		{m_imp_temp_start_val = imp_temp_start_val;}

	/**
	* Sets the imp_vX_start_val option.
	*
	* Velocity in Cartesian X that every impurity starts with (m/s), added
	* to the velocity sampled at the starting temperature. With
	* imp_temp_start_val = 0 every impurity starts with exactly this
	* velocity. Not used for sputtered impurities.
	*
	* @param imp_vX_start_val New value of the imp_vX_start_val option.
	*/
	void Options::set_imp_vX_start_val(double imp_vX_start_val)
		{m_imp_vX_start_val = imp_vX_start_val;}

	/**
	* Sets the imp_vY_start_val option.
	*
	* Velocity in Cartesian Y that every impurity starts with (m/s), added
	* to the velocity sampled at the starting temperature. See
	* imp_vX_start_val.
	*
	* @param imp_vY_start_val New value of the imp_vY_start_val option.
	*/
	void Options::set_imp_vY_start_val(double imp_vY_start_val)
		{m_imp_vY_start_val = imp_vY_start_val;}

	/**
	* Sets the imp_vZ_start_val option.
	*
	* Velocity in Cartesian Z that every impurity starts with (m/s), added
	* to the velocity sampled at the starting temperature. See
	* imp_vX_start_val.
	*
	* @param imp_vZ_start_val New value of the imp_vZ_start_val option.
	*/
	void Options::set_imp_vZ_start_val(double imp_vZ_start_val)
		{m_imp_vZ_start_val = imp_vZ_start_val;}

	/**
	* Sets the imp_collisions option.
	*
	* Collision model for the impurities, either off or the Nanbu model.
	*
	* Valid values are "off" and "nanbu". An invalid value prints an error
	* and leaves the option unchanged. Also sets the control integer
	* returned by imp_collisions_int() to 0 for "off", 1 for "nanbu". An
	* invalid value leaves the control integer unchanged.
	*
	* @param imp_collisions New value of the imp_collisions option.
	*/
	void Options::set_imp_collisions(std::string imp_collisions) 
		{
			// Only assign if valid option, leaving as default if not
			if (check_input<std::string>("imp_collision", imp_collisions, 
				{"off", "nanbu"}))
			{
				m_imp_collisions = imp_collisions;
			}

			// Assign control integers
			if (imp_collisions == "off") m_imp_collisions_int = 0;
			else if (imp_collisions == "nanbu") m_imp_collisions_int = 1;
		}

	/**
	* Sets the imp_thermal_force option.
	*
	* Whether the Nanbu collisions include the thermal force: "on" samples
	* the colliding background electrons and ions from distributions skewed
	* by their parallel heat fluxes (see Collisions::sample_bkg_velocity),
	* "off" from Maxwellians. "on" (the default) only has an effect with a
	* background that has heat fluxes (SOLEDGE3X files with face fluxes),
	* and then needs the CPU (see Collisions::check_thermal_force).
	*
	* Valid values are "off" and "on". An invalid value prints an error and
	* leaves the option unchanged. Also sets the control integer returned
	* by imp_thermal_force_int() to 0 for "off", 1 for "on". An invalid
	* value leaves the control integer unchanged.
	*
	* @param imp_thermal_force New value of the imp_thermal_force option.
	*/
	void Options::set_imp_thermal_force(std::string imp_thermal_force)
	{
		if (check_input<std::string>("imp_thermal_force", imp_thermal_force,
			{"off", "on"}))
		{
			m_imp_thermal_force = imp_thermal_force;
		}

		if (imp_thermal_force == "off") m_imp_thermal_force_int = 0;
		else if (imp_thermal_force == "on") m_imp_thermal_force_int = 1;
	}

	/**
	* Sets the imp_coll_subcycle_target option.
	*
	* The Nanbu collision of one time step is split into sub-steps, each
	* with a newly sampled background particle, so that each sub-step is
	* at most this collisional. It is measured by s0, the Nanbu s term for
	* a background particle at its thermal speed (see
	* Collisions::nanbu_num_substeps for why this is needed). Smaller
	* values bring the friction and thermal force closer to their fluid
	* values and cost more sub-steps. Must be above 0, or the run is
	* stopped.
	*
	* @param imp_coll_subcycle_target New value of the
	*        imp_coll_subcycle_target option.
	*/
	void Options::set_imp_coll_subcycle_target(
		double imp_coll_subcycle_target)
	{
		if (!(imp_coll_subcycle_target > 0.0))
		{
			Errors::fatal_error("imp_coll_subcycle_target must be above 0, "
				"got " + std::to_string(imp_coll_subcycle_target) + ". To "
				"turn collision sub-cycling off, set imp_coll_subcycle_cap "
				"= 1.");
		}
		m_imp_coll_subcycle_target = imp_coll_subcycle_target;
	}

	/**
	* Sets the imp_coll_subcycle_cap option.
	*
	* The largest number of sub-steps one Nanbu collision is split into
	* (see set_imp_coll_subcycle_target). It limits the cost in very
	* collisional plasma. 1 turns sub-cycling off. Must be at least 1, or
	* the run is stopped.
	*
	* @param imp_coll_subcycle_cap New value of the imp_coll_subcycle_cap
	*        option.
	*/
	void Options::set_imp_coll_subcycle_cap(int imp_coll_subcycle_cap)
	{
		if (imp_coll_subcycle_cap < 1)
		{
			Errors::fatal_error("imp_coll_subcycle_cap must be at least 1 "
				"(1 turns collision sub-cycling off), got "
				+ std::to_string(imp_coll_subcycle_cap) + ".");
		}
		m_imp_coll_subcycle_cap = imp_coll_subcycle_cap;
	}

	/**
	* Sets the var_red_split option.
	*
	* When particles are split for variance reduction: never, at
	* ionization/recombination events, or at collisions.
	*
	* Valid values are "off", "iz_rec" and "coll". An invalid value prints
	* an error and leaves the option unchanged. Also sets the control
	* integer returned by var_red_split_int() to 0 for "off", 1 for
	* "iz_rec", 2 for "coll". An invalid value leaves the control integer
	* unchanged.
	*
	* @param var_red_split New value of the var_red_split option.
	*/
	void Options::set_var_red_split(std::string var_red_split) 
		{
			// Only assign if valid option, leaving as default if not
			if (check_input<std::string>("var_red_split", var_red_split, 
				{"off", "iz_rec", "coll"}))
			{
				m_var_red_split = var_red_split;
			}

			// Assign control integers
			if (var_red_split == "off") m_var_red_split_int = 0;
			else if (var_red_split == "iz_rec") m_var_red_split_int = 1;
			else if (var_red_split == "coll") m_var_red_split_int = 2;
		}

	/**
	* Sets the var_red_import option.
	*
	* Method used to calculate particle importance for variance reduction.
	*
	* Valid values are "median", "exp_dist" and "exp_time". An invalid value
	* prints an error and leaves the option unchanged. Also sets the control
	* integer returned by var_red_import_int() to 0 for "median", 1 for
	* "exp_dist", 2 for "exp_time". An invalid value leaves the control
	* integer unchanged.
	*
	* @param var_red_import New value of the var_red_import option.
	*/
	void Options::set_var_red_import(std::string var_red_import) 
		{
			// Only assign if valid option, leaving as default if not
			if (check_input<std::string>("var_red_import", var_red_import, 
				{"median", "exp_dist", "exp_time"}))
			{
				m_var_red_import = var_red_import;
			}

			// Assign control integers
			if (var_red_import == "median") m_var_red_import_int = 0;
			else if (var_red_import == "exp_dist") m_var_red_import_int = 1;
			else if (var_red_import == "exp_time") m_var_red_import_int = 2;
		}

	/**
	* Sets the var_red_freq option.
	*
	* How often variance reduction is applied.
	*
	* @param var_red_freq New value of the var_red_freq option.
	*/
	void Options::set_var_red_freq(double var_red_freq) 
		{m_var_red_freq = var_red_freq;}
	
	/**
	* Sets the var_red_min_weight option.
	*
	* Minimum particle weight for variance reduction.
	*
	* @param var_red_min_weight New value of the var_red_min_weight option.
	*/
	void Options::set_var_red_min_weight(double var_red_min_weight) 
		{m_var_red_min_weight = var_red_min_weight;}
	
	/**
	* Sets the var_red_med_mod option.
	*
	* Modifier applied to the median importance for variance reduction.
	*
	* @param var_red_med_mod New value of the var_red_med_mod option.
	*/
	void Options::set_var_red_med_mod(double var_red_med_mod) 
		{m_var_red_med_mod = var_red_med_mod;}
	
	/**
	* Sets the var_red_rusrol option.
	*
	* Whether to use Russian roulette for variance reduction.
	*
	* Valid values are "off" and "on". An invalid value prints an error and
	* leaves the option unchanged. Also sets the control integer returned by
	* var_red_rusrol_int() to 0 for "off", 1 for "on". An invalid value
	* leaves the control integer unchanged.
	*
	* @param var_red_rusrol New value of the var_red_rusrol option.
	*/
	void Options::set_var_red_rusrol(std::string var_red_rusrol) 
		{
			// Only assign if valid option, leaving as default if not
			if (check_input<std::string>("var_red_rusrol", var_red_rusrol, 
				{"off", "on"}))
			{
				m_var_red_rusrol = var_red_rusrol;
			}

			// Assign control integers
			if (var_red_rusrol == "off") m_var_red_rusrol_int = 0;
			else if (var_red_rusrol == "on") m_var_red_rusrol_int = 1;
		}

	/**
	* Sets the var_red_rusrol_prob option.
	*
	* Probability used for Russian roulette.
	*
	* @param var_red_rusrol_prob New value of the var_red_rusrol_prob
	*                            option.
	*/
	void Options::set_var_red_rusrol_prob(double var_red_rusrol_prob) 
		{m_var_red_rusrol_prob = var_red_rusrol_prob;}

	/**
	* Sets the imp_time_step_opt option.
	*
	* Whether the impurity time step is constant or variable.
	*
	* Valid values are "constant" and "variable". An invalid value prints an
	* error and leaves the option unchanged. Also sets the control integer
	* returned by imp_time_step_opt_int() to 0 for "constant", 1 for
	* "variable". An invalid value leaves the control integer unchanged.
	*
	* @param imp_time_step_opt New value of the imp_time_step_opt option.
	*/
	void Options::set_imp_time_step_opt(std::string imp_time_step_opt) 
		{
			// Only assign if valid option, leaving as default if not
			if (check_input<std::string>("imp_time_step_opt", imp_time_step_opt, 
				{"constant", "variable"}))
			{
				m_imp_time_step_opt = imp_time_step_opt;
			}

			// Assign control integers
			if (imp_time_step_opt == "constant") m_imp_time_step_opt_int = 0;
			else if (imp_time_step_opt == "variable") 
				m_imp_time_step_opt_int = 1;
		}
	
	/**
	* Sets the imp_time_step option.
	*
	* Impurity time step (s).
	*
	* @param imp_time_step New value of the imp_time_step option.
	*/
	void Options::set_imp_time_step(double imp_time_step) 
		{m_imp_time_step = imp_time_step;}
	
	/**
	* Sets the imp_time_step_min option.
	*
	* Minimum impurity time step when the time step is variable (s).
	*
	* @param imp_time_step_min New value of the imp_time_step_min option.
	*/
	void Options::set_imp_time_step_min(double imp_time_step_min) 
		{m_imp_time_step_min = imp_time_step_min;}
	
	/**
	* Sets the imp_source_scale_fact option.
	*
	* Impurity source scaling factor (particles/s), i.e., how many real
	* particles a single simulated particle represents. Used to convert the
	* Monte Carlo counts into a density.
	*
	* @param imp_source_scale_fact New value of the imp_source_scale_fact
	*                              option.
	*/
	void Options::set_imp_source_scale_fact(double imp_source_scale_fact) 
		{m_imp_source_scale_fact = imp_source_scale_fact;}
	
	/**
	* Sets the imp_xbound_buffer option.
	*
	* Buffer distance at the x boundaries. With an absorbing boundary, a
	* particle is absorbed once it is within this distance of the boundary.
	*
	* @param imp_xbound_buffer New value of the imp_xbound_buffer option.
	*/
	void Options::set_imp_xbound_buffer(double imp_xbound_buffer) 
		{m_imp_xbound_buffer = imp_xbound_buffer;}

	/**
	* Sets the imp_ybound_buffer option.
	*
	* Buffer distance at the y boundaries. With an absorbing boundary, a
	* particle is absorbed once it is within this distance of the boundary.
	*
	* @param imp_ybound_buffer New value of the imp_ybound_buffer option.
	*/
	void Options::set_imp_ybound_buffer(double imp_ybound_buffer) 
		{m_imp_ybound_buffer = imp_ybound_buffer;}

	/**
	* Sets the imp_zbound_buffer option.
	*
	* Buffer distance at the z boundaries. With an absorbing boundary, a
	* particle is absorbed once it is within this distance of the boundary.
	*
	* @param imp_zbound_buffer New value of the imp_zbound_buffer option.
	*/
	void Options::set_imp_zbound_buffer(double imp_zbound_buffer) 
		{m_imp_zbound_buffer = imp_zbound_buffer;}
	
	/**
	* Sets the imp_iz_recomb option.
	*
	* Whether impurity ionization and recombination are on.
	*
	* Valid values are "off" and "on". An invalid value prints an error and
	* leaves the option unchanged. Also sets the control integer returned by
	* imp_iz_recomb_int() to 0 for "off", 1 for "on". An invalid value
	* leaves the control integer unchanged.
	*
	* @param imp_iz_recomb New value of the imp_iz_recomb option.
	*/
	void Options::set_imp_iz_recomb(std::string imp_iz_recomb) 
		{
			// Only assign if valid option, leaving as default if not
			if (check_input<std::string>("imp_iz_recomb", imp_iz_recomb, 
				{"off", "on"}))
			{
				m_imp_iz_recomb = imp_iz_recomb;
			}

			// Set control integers
			if (imp_iz_recomb == "off") m_imp_iz_recomb_int = 0;
			else if (imp_iz_recomb == "on") m_imp_iz_recomb_int = 1;
		}

	/**
	* Sets the imp_max_time option.
	*
	* Maximum time a particle is followed for (s). A particle that exceeds
	* it is killed and counted in a warning at the end of the run.
	*
	* @param imp_max_time New value of the imp_max_time option.
	*/
	void Options::set_imp_max_time(double imp_max_time) 
		{m_imp_max_time = imp_max_time;}
	
	/**
	* Sets the print_interval option.
	*
	* Number of progress updates printed while particles are followed. For
	* example, 10 prints an update every 10% of particles.
	*
	* @param print_interval New value of the print_interval option.
	*/
	void Options::set_print_interval(int print_interval) 
		{m_print_interval = print_interval;}
	
	/**
	* Sets the openadas_root option.
	*
	* Directory containing the OpenADAS rate coefficient files.
	*
	* @param openadas_root New value of the openadas_root option.
	*/
	void Options::set_openadas_root(std::string openadas_root) 
		{m_openadas_root = openadas_root;}
	
	/**
	* Sets the openadas_year option.
	*
	* Year of the OpenADAS data files to use.
	*
	* @param openadas_year New value of the openadas_year option.
	*/
	void Options::set_openadas_year(int openadas_year) 
		{m_openadas_year = openadas_year;}
	
	// Accessor definitions

	/**
	* Returns the case_name option.
	*
	* Name of the simulation case. The results are saved to <case_name>.nc.
	*
	* @return Value of the case_name option.
	*/
	const std::string& Options::case_name() const 
		{return m_case_name;}

	/**
	* Returns the use_gpu option.
	*
	* Whether to run the simulation on GPUs.
	*
	* @return Value of the use_gpu option.
	*/
	const std::string& Options::use_gpu() const 
		{return m_use_gpu;}

	/**
	* Returns the seed option.
	*
	* Base seed for the random number generators. A negative value seeds
	* them from the clock so each run is different. The MPI rank is added to
	* the seed so each rank gets a different sequence.
	*
	* @return Value of the seed option.
	*/
	int Options::seed() const 
		{return m_seed;}

	/**
	* Returns the slot_cap option.
	*
	* Maximum number of particle slots, i.e., the most particles followed at
	* once. The number of slots used is the smaller of this and imp_num.
	*
	* @return Value of the slot_cap option.
	*/
	int Options::slot_cap() const 
		{return m_slot_cap;}

	/**
	* Returns the bkg_source option.
	*
	* Source of the background plasma: "slab", "gkeyll" or "soledge3x" (see
	* set_bkg_source).
	*
	* @return Value of the bkg_source option.
	*/
	const std::string& Options::bkg_source() const 
		{return m_bkg_source;}


	/**
	* Returns the save_track option.
	*
	* Whether to save particle tracks to the output file.
	*
	* @return Value of the save_track option.
	*/
	const std::string& Options::save_track() const 
		{return m_save_track;}

	/**
	* Returns the gkyl_dir option.
	*
	* Directory containing the Gkeyll simulation files.
	*
	* @return Value of the gkyl_dir option.
	*/
	const std::string& Options::gkyl_dir() const 
		{return m_gkyl_dir;}

	/**
	* Returns the gkyl_casename option.
	*
	* Name of the Gkeyll case, i.e., the prefix of its output files.
	*
	* @return Value of the gkyl_casename option.
	*/
	const std::string& Options::gkyl_casename() const 
		{return m_gkyl_casename;}

	/**
	* Returns the gkyl_frame_start option.
	*
	* First Gkeyll frame to load.
	*
	* @return Value of the gkyl_frame_start option.
	*/
	int Options::gkyl_frame_start() const 
		{return m_gkyl_frame_start;}

	/**
	* Returns the gkyl_frame_end option.
	*
	* Last Gkeyll frame to load. The range from gkyl_frame_start to
	* gkyl_frame_end is inclusive.
	*
	* @return Value of the gkyl_frame_end option.
	*/
	int Options::gkyl_frame_end() const 
		{return m_gkyl_frame_end;}

	/**
	* Returns the gkyl_elec_name option.
	*
	* Name of the electron species in the Gkeyll simulation.
	*
	* @return Value of the gkyl_elec_name option.
	*/
	const std::string& Options::gkyl_elec_name() const 
		{return m_gkyl_elec_name;}

	/**
	* Returns the gkyl_ion_name option.
	*
	* Name of the main ion species in the Gkeyll simulation.
	*
	* @return Value of the gkyl_ion_name option.
	*/
	const std::string& Options::gkyl_ion_name() const 
		{return m_gkyl_ion_name;}

	/**
	* Returns the gkyl_elec_mass_amu option.
	*
	* Electron mass used in the Gkeyll simulation (amu).
	*
	* @return Value of the gkyl_elec_mass_amu option.
	*/
	double Options::gkyl_elec_mass_amu() const 
		{return m_gkyl_elec_mass_amu;}

	/**
	* Returns the gkyl_ion_mass_amu option.
	*
	* Main ion mass used in the Gkeyll simulation (amu).
	*
	* @return Value of the gkyl_ion_mass_amu option.
	*/
	double Options::gkyl_ion_mass_amu() const 
		{return m_gkyl_ion_mass_amu;}

	/**
	* Returns the gkyl_file_type option.
	*
	* File format of the Gkeyll output.
	*
	* @return Value of the gkyl_file_type option.
	*/
	const std::string& Options::gkyl_file_type() const 
		{return m_gkyl_file_type;}

	/**
	* Returns the gkyl_moment_type option.
	*
	* Type of moment files in the Gkeyll output. Collisions need
	* "bimaxwellian", which provides the ion parallel flow and perpendicular
	* velocity used to calculate the background flow.
	*
	* @return Value of the gkyl_moment_type option.
	*/
	const std::string& Options::gkyl_moment_type() const 
		{return m_gkyl_moment_type;}

	/**
	* Returns the calc_grad_elec option.
	*
	* Whether to calculate the electric field gradient.
	*
	* @return Value of the calc_grad_elec option.
	*/
	const std::string& Options::calc_grad_elec() const 
		{return m_calc_grad_elec;}

	/**
	* Returns the s3x_dir option.
	*
	* Directory containing the SOLEDGE3X mesh and plasma files.
	*
	* @return Value of the s3x_dir option.
	*/
	const std::string& Options::s3x_dir() const 
		{return m_s3x_dir;}

	/**
	* Returns the s3x_mesh_file option.
	*
	* Name of the SOLEDGE3X mesh file in s3x_dir.
	*
	* @return Value of the s3x_mesh_file option.
	*/
	const std::string& Options::s3x_mesh_file() const 
		{return m_s3x_mesh_file;}

	/**
	* Returns the s3x_frame_start option.
	*
	* First SOLEDGE3X plasma frame to load, i.e., the number in
	* plasma_%05d.h5.
	*
	* @return Value of the s3x_frame_start option.
	*/
	int Options::s3x_frame_start() const 
		{return m_s3x_frame_start;}

	/**
	* Returns the s3x_frame_end option.
	*
	* Last SOLEDGE3X plasma frame to load. The range from
	* s3x_frame_start to s3x_frame_end is inclusive, and must have at least
	* two frames.
	*
	* @return Value of the s3x_frame_end option.
	*/
	int Options::s3x_frame_end() const 
		{return m_s3x_frame_end;}

	/**
	* Returns the s3x_ref_file option.
	*
	* Name of the SOLEDGE3X file in s3x_dir with the reference values the
	* output is normalized by (n0, T0, c0, tau0, phi0, ...).
	*
	* @return Value of the s3x_ref_file option.
	*/
	const std::string& Options::s3x_ref_file() const 
		{return m_s3x_ref_file;}

	/**
	* Returns the s3x_chi_absorb option.
	*
	* Whether particles are absorbed when they enter a cell inside the
	* SOLEDGE3X wall (chi = 1).
	*
	* @return Value of the s3x_chi_absorb option.
	*/
	const std::string& Options::s3x_chi_absorb() const 
		{return m_s3x_chi_absorb;}

	/**
	* Returns the hermes_format option.
	*
	* Output format of the HERMES run: "jax" for JAX-Hermes3 (a hypnotoad grid
	* plus one .npz snapshot per time) or "bout" for BOUT++ Hermes-3 dump
	* files (one per MPI process).
	*
	* @return Value of the hermes_format option.
	*/
	const std::string& Options::hermes_format() const
		{return m_hermes_format;}

	/**
	* Returns the hermes_grid_file option.
	*
	* Path to the hypnotoad grid file (netCDF) of the HERMES run.
	*
	* @return Value of the hermes_grid_file option.
	*/
	const std::string& Options::hermes_grid_file() const
		{return m_hermes_grid_file;}

	/**
	* Returns the hermes_snap_dir option.
	*
	* Directory with the HERMES snapshots (snap*.npz for JAX-Hermes3), or
	* the BOUT++ run directory with BOUT.dmp.*.nc (BOUT++ Hermes-3).
	*
	* @return Value of the hermes_snap_dir option.
	*/
	const std::string& Options::hermes_snap_dir() const
		{return m_hermes_snap_dir;}

	/**
	* Returns the hermes_ion_species option.
	*
	* Name of the main ion species in the Hermes-3 input (BOUT++ dumps).
	*
	* @return Value of the hermes_ion_species option.
	*/
	const std::string& Options::hermes_ion_species() const
		{return m_hermes_ion_species;}

	/**
	* Returns the hermes_snap_start option.
	*
	* First snapshot to load, as an index into the snapshots in
	* hermes_snap_dir sorted by name (time), from 0.
	*
	* @return Value of the hermes_snap_start option.
	*/
	int Options::hermes_snap_start() const
		{return m_hermes_snap_start;}

	/**
	* Returns the hermes_snap_end option.
	*
	* Last snapshot to load (inclusive), as an index like hermes_snap_start.
	* Negative counts from the end (-1 is the last).
	*
	* @return Value of the hermes_snap_end option.
	*/
	int Options::hermes_snap_end() const
		{return m_hermes_snap_end;}

	/**
	* Returns the hermes_snap_stride option.
	*
	* Load every hermes_snap_stride-th snapshot from hermes_snap_start to
	* hermes_snap_end.
	*
	* @return Value of the hermes_snap_stride option.
	*/
	int Options::hermes_snap_stride() const
		{return m_hermes_snap_stride;}

	/**
	* Returns the hermes_nnorm option.
	*
	* Density (m^-3) the HERMES output is normalized by (Nnorm).
	*
	* @return Value of the hermes_nnorm option.
	*/
	double Options::hermes_nnorm() const
		{return m_hermes_nnorm;}

	/**
	* Returns the hermes_tnorm option.
	*
	* Temperature (eV) the HERMES output is normalized by (Tnorm).
	*
	* @return Value of the hermes_tnorm option.
	*/
	double Options::hermes_tnorm() const
		{return m_hermes_tnorm;}

	/**
	* Returns the hermes_bnorm option.
	*
	* Magnetic field (T) the HERMES output is normalized by (Bnorm), which
	* sets the time unit 1 / Omega_ci = m_p / (e Bnorm).
	*
	* @return Value of the hermes_bnorm option.
	*/
	double Options::hermes_bnorm() const
		{return m_hermes_bnorm;}

	/**
	* Returns the hermes_ion_mass_amu option.
	*
	* Mass number (AA) of the HERMES main ions, used to turn NVi into the
	* parallel flow and as the background ion mass.
	*
	* @return Value of the hermes_ion_mass_amu option.
	*/
	double Options::hermes_ion_mass_amu() const
		{return m_hermes_ion_mass_amu;}

	/**
	* Returns the hermes_zperiod option.
	*
	* The HERMES domain covers 1 / hermes_zperiod of the torus (zperiod in the
	* HERMES input).
	*
	* @return Value of the hermes_zperiod option.
	*/
	int Options::hermes_zperiod() const
		{return m_hermes_zperiod;}

	/**
	* Returns the hermes_heat_flux option ("conduction" or "off").
	*
	* @return Value of the hermes_heat_flux option.
	*/
	const std::string& Options::hermes_heat_flux() const
		{return m_hermes_heat_flux;}

	/**
	* Returns the hermes_inductive_epar option ("off" or "on").
	*
	* @return Value of the hermes_inductive_epar option.
	*/
	const std::string& Options::hermes_inductive_epar() const
		{return m_hermes_inductive_epar;}

	/**
	* Returns the hermes_field_aligned option ("on" or "off").
	*
	* @return Value of the hermes_field_aligned option.
	*/
	const std::string& Options::hermes_field_aligned() const
		{return m_hermes_field_aligned;}

	/**
	* Returns the hermes_kappa_coef_e option.
	*
	* @return Electron parallel conduction coefficient.
	*/
	double Options::hermes_kappa_coef_e() const
		{return m_hermes_kappa_coef_e;}

	/**
	* Returns the hermes_kappa_coef_i option.
	*
	* @return Ion parallel conduction coefficient.
	*/
	double Options::hermes_kappa_coef_i() const
		{return m_hermes_kappa_coef_i;}

	/**
	* Returns the hermes_kappa_limit_alpha option.
	*
	* @return Conduction flux limiter, zero or negative for none.
	*/
	double Options::hermes_kappa_limit_alpha() const
		{return m_hermes_kappa_limit_alpha;}

	/**
	* @name Sputtering option getters
	* See the setters for what each option means.
	*/
	///@{
	const std::string& Options::imp_source_opt() const
		{return m_imp_source_opt;}
	const std::string& Options::sput_table_file() const
		{return m_sput_table_file;}
	const std::string& Options::sput_flux_opt() const
		{return m_sput_flux_opt;}
	const std::string& Options::sput_bohm_floor() const
		{return m_sput_bohm_floor;}

	/**
	* Sets the sheath_fd option.
	*
	* The share of the sheath drop in the Debye sheath, f_D, in place of
	* the one from the field's angle to the wall (1 - ln cos(theta) /
	* lam_w). For comparing with models that fix it, such as Brooks'
	* reference case (0.25), and for sensitivity studies. -1 (the
	* default) works it out from the angle.
	*
	* Valid values are -1 and 0 to 1. An invalid value prints an error and
	* leaves the option unchanged.
	*
	* @param sheath_fd New value of the sheath_fd option.
	*/
	void Options::set_sheath_fd(double sheath_fd)
	{
		if (sheath_fd == -1.0 || (sheath_fd >= 0.0 && sheath_fd <= 1.0))
			m_sheath_fd = sheath_fd;
		else
		{
			std::cerr << "Error! sheath_fd must be -1 or from 0 to 1, not "
				<< sheath_fd << ". Leaving it at " << m_sheath_fd << ".\n";
		}
	}

	/**
	* Sets the sheath_drop option.
	*
	* The sheath's potential drop in units of Te (e |phi_wall| / Te), in
	* place of the floating drop (-0.5 ln[2 pi (m_e / m_i)(1 + Ti / Te)],
	* about 2.8 for D at Ti = Te). For comparing with models that fix it,
	* such as Brooks' reference case (3), and for checks with a vanishing
	* drop. -1 (the default) uses the floating drop.
	*
	* Valid values are -1 and above 0. An invalid value prints an error and
	* leaves the option unchanged.
	*
	* @param sheath_drop New value of the sheath_drop option.
	*/
	void Options::set_sheath_drop(double sheath_drop)
	{
		if (sheath_drop == -1.0 || sheath_drop > 0.0)
			m_sheath_drop = sheath_drop;
		else
		{
			std::cerr << "Error! sheath_drop must be -1 or above 0, not "
				<< sheath_drop << ". Leaving it at " << m_sheath_drop
				<< ".\n";
		}
	}

	/**
	* Returns the sheath_fd option.
	*
	* @return Value of the sheath_fd option.
	*/
	double Options::sheath_fd() const {return m_sheath_fd;}

	/**
	* Returns the sheath_drop option.
	*
	* @return Value of the sheath_drop option.
	*/
	double Options::sheath_drop() const {return m_sheath_drop;}

	/**
	* Sets the sput_self option.
	*
	* Whether impurities that land on a wall element sputter more atoms
	* (self-sputtering). The particle carries on as one of the atoms it
	* sputtered, with its weight multiplied by the yield from the
	* sput_self_table_file table (see Sputtering::SelfSource). Needs
	* imp_source_opt = "sputter".
	*
	* Valid values are "off" and "on"; anything else stops the run
	* (check_input). Also sets the control integer returned by
	* sput_self_int() to 0 for "off", 1 for "on".
	*
	* @param sput_self New value of the sput_self option.
	*/
	void Options::set_sput_self(std::string sput_self)
	{
		if (check_input<std::string>("sput_self", sput_self, {"off", "on"}))
			m_sput_self = sput_self;

		if (sput_self == "off") m_sput_self_int = 0;
		else if (sput_self == "on") m_sput_self_int = 1;
	}

	/**
	* Sets the sput_impact_opt option.
	*
	* Where the energies and angles the background's ions hit the wall
	* with, for the sputtered source, come from. "formula": the sheath and
	* presheath energy (Sputtering::sheath_energy) plus a thermal energy
	* at Ti, all at sput_angle_val. "sheath": D+ test ions pushed through
	* the sheath_model field (Sheath::ImpactTable), which gives the
	* energies and angles together. "sheath" is the default and needs
	* sheath_model to be on, so a run with the sheath off must set
	* "formula". The test ions are pushed on the CPU, also for a GPU run.
	*
	* Valid values are "formula" and "sheath"; anything else stops the run
	* (check_input). Also sets the control integer returned by
	* sput_impact_opt_int() to 0 for "formula", 1 for "sheath".
	*
	* @param sput_impact_opt New value of the sput_impact_opt option.
	*/
	void Options::set_sput_impact_opt(std::string sput_impact_opt)
	{
		if (check_input<std::string>("sput_impact_opt", sput_impact_opt,
			{"formula", "sheath"}))
		{
			m_sput_impact_opt = sput_impact_opt;
		}

		if (sput_impact_opt == "formula") m_sput_impact_opt_int = 0;
		else if (sput_impact_opt == "sheath") m_sput_impact_opt_int = 1;
	}

	/**
	* Sets the sput_self_table_file option.
	*
	* Path to the self-sputtering yield table (the impurity on the wall
	* material, e.g. W on W), made by python/rustbca_yields.py. Needed
	* with sput_self = "on".
	*
	* @param sput_self_table_file New value of the sput_self_table_file
	*        option.
	*/
	void Options::set_sput_self_table_file(std::string sput_self_table_file)
		{m_sput_self_table_file = sput_self_table_file;}

	/**
	* Sets the sput_self_min_weight option.
	*
	* The weight, as a fraction of a particle's weight from the D+
	* sputtered source, below which a self-sputtered particle plays Russian
	* roulette: it carries on with this weight with probability (its new
	* weight) / (this weight), and ends otherwise. That keeps the expected
	* weight, while particles of little weight don't cost as much as the
	* others.
	*
	* Valid values are 0 to 1. 0 turns roulette off. An invalid value
	* prints an error and leaves the option unchanged.
	*
	* @param sput_self_min_weight New value of the sput_self_min_weight
	*        option.
	*/
	void Options::set_sput_self_min_weight(double sput_self_min_weight)
	{
		if (sput_self_min_weight >= 0.0 && sput_self_min_weight <= 1.0)
			m_sput_self_min_weight = sput_self_min_weight;
		else
		{
			std::cerr << "Error! sput_self_min_weight must be from 0 to 1, "
				"not " << sput_self_min_weight << ". Leaving it at "
				<< m_sput_self_min_weight << ".\n";
		}
	}

	/**
	* Sets the sput_self_max_gen option.
	*
	* The most self-sputtering generations a particle can go through. A
	* particle reaching it stops the run: with roulette, chains only get
	* that long when the self-sputtering multiplication (the share coming
	* back times the mean yield) is close to or above 1, the
	* self-sputtering avalanche, where the self-sputtered rate has no
	* steady value.
	*
	* Valid values are 1 and above. An invalid value prints an error and
	* leaves the option unchanged.
	*
	* @param sput_self_max_gen New value of the sput_self_max_gen option.
	*/
	void Options::set_sput_self_max_gen(int sput_self_max_gen)
	{
		if (sput_self_max_gen >= 1) m_sput_self_max_gen = sput_self_max_gen;
		else
		{
			std::cerr << "Error! sput_self_max_gen must be 1 or above, not "
				<< sput_self_max_gen << ". Leaving it at "
				<< m_sput_self_max_gen << ".\n";
		}
	}

	/**
	* @name Self-sputtering option getters
	* See the setters for what each option means.
	*/
	///@{
	const std::string& Options::sput_self() const {return m_sput_self;}

	/**
	* Returns the sput_impact_opt option.
	*
	* @return Value of the sput_impact_opt option.
	*/
	const std::string& Options::sput_impact_opt() const
		{return m_sput_impact_opt;}
	const std::string& Options::sput_self_table_file() const
		{return m_sput_self_table_file;}
	double Options::sput_self_min_weight() const
		{return m_sput_self_min_weight;}
	int Options::sput_self_max_gen() const {return m_sput_self_max_gen;}
	///@}

	/**
	* Returns the sheath_model option.
	*
	* @return Value of the sheath_model option.
	*/
	const std::string& Options::sheath_model() const
		{return m_sheath_model;}
	double Options::sput_angle_val() const {return m_sput_angle_val;}
	double Options::sput_max_wall_dist() const
		{return m_sput_max_wall_dist;}
	///@}

	/**
	* Returns the min_ne option.
	*
	* Lower limit on the electron density (m^-3).
	*
	* @return Value of the min_ne option.
	*/
	double Options::min_ne() const 
		{return m_min_ne;}

	/**
	* Returns the min_te option.
	*
	* Lower limit on the electron temperature (eV).
	*
	* @return Value of the min_te option.
	*/
	double Options::min_te() const 
		{return m_min_te;}

	/**
	* Returns the min_ti option.
	*
	* Lower limit on the ion temperature (eV), applied when calculating
	* collisions with the background ions.
	*
	* @return Value of the min_ti option.
	*/
	double Options::min_ti() const 
		{return m_min_ti;}

	/**
	* Returns the slab_ne option (slab background: electron (and ion) density).
	*
	* @return Value of the slab_ne option.
	*/
	double Options::slab_ne() const
		{return m_slab_ne;}

	/**
	* Returns the slab_te option (slab background: electron temperature).
	*
	* @return Value of the slab_te option.
	*/
	double Options::slab_te() const
		{return m_slab_te;}

	/**
	* Returns the slab_ti option (slab background: ion temperature).
	*
	* @return Value of the slab_ti option.
	*/
	double Options::slab_ti() const
		{return m_slab_ti;}

	/**
	* Returns the slab_flow option (slab background: plasma flow along B (+X)).
	*
	* @return Value of the slab_flow option.
	*/
	double Options::slab_flow() const
		{return m_slab_flow;}

	/**
	* Returns the slab_heat_flux_i option (slab background: ion parallel heat flux / (n Ti v_th,i)).
	*
	* @return Value of the slab_heat_flux_i option.
	*/
	double Options::slab_heat_flux_i() const
		{return m_slab_heat_flux_i;}

	/**
	* Returns the slab_heat_flux_e option (slab background: electron parallel heat flux / (n Te v_th,e)).
	*
	* @return Value of the slab_heat_flux_e option.
	*/
	double Options::slab_heat_flux_e() const
		{return m_slab_heat_flux_e;}

	/**
	* Returns the slab_bmag option (slab background: magnetic field strength (along +X)).
	*
	* @return Value of the slab_bmag option.
	*/
	double Options::slab_bmag() const
		{return m_slab_bmag;}

	/**
	* Returns the slab_tmax option (slab background: end of the time range).
	*
	* @return Value of the slab_tmax option.
	*/
	double Options::slab_tmax() const
		{return m_slab_tmax;}

	/**
	* Returns the slab_efield_y option (slab background: electric field in +Y at t = 0).
	*
	* @return Value of the slab_efield_y option.
	*/
	double Options::slab_efield_y() const
		{return m_slab_efield_y;}

	/**
	* Returns the slab_efield_y_rate option (slab background: rate of change of the electric field in +Y).
	*
	* @return Value of the slab_efield_y_rate option.
	*/
	double Options::slab_efield_y_rate() const
		{return m_slab_efield_y_rate;}

	/**
	* Returns the slab_bgrad_y option (slab background: relative gradient of the field strength in +Y).
	*
	* @return Value of the slab_bgrad_y option.
	*/
	double Options::slab_bgrad_y() const
		{return m_slab_bgrad_y;}

	/**
	* Returns the slab_major_radius option (slab background: major radius the box is bent to (0 = straight)).
	*
	* @return Value of the slab_major_radius option.
	*/
	double Options::slab_major_radius() const
		{return m_slab_major_radius;}

	/**
	* Returns the slab_wall option (slab background: which face is a wall).
	*
	* @return Value of the slab_wall option.
	*/
	const std::string& Options::slab_wall() const
		{return m_slab_wall;}

	/**
	* Returns the slab_b_angle option (slab background: angle of the field
	* from +X toward -Y).
	*
	* @return Value of the slab_b_angle option.
	*/
	double Options::slab_b_angle() const
		{return m_slab_b_angle;}

	/**
	* Returns the tbound_type option.
	*
	* Boundary condition at the start and end of the background time range.
	*
	* @return Value of the tbound_type option.
	*/
	const std::string& Options::tbound_type() const 
		{return m_tbound_type;}

	/**
	* Returns the min_xbound_type option.
	*
	* Boundary condition at the minimum x boundary.
	*
	* @return Value of the min_xbound_type option.
	*/
	const std::string& Options::min_xbound_type() const 
		{return m_min_xbound_type;}

	/**
	* Returns the max_xbound_type option.
	*
	* Boundary condition at the maximum x boundary.
	*
	* @return Value of the max_xbound_type option.
	*/
	const std::string& Options::max_xbound_type() const 
		{return m_max_xbound_type;}

	/**
	* Returns the min_ybound_type option.
	*
	* Boundary condition at the minimum y boundary.
	*
	* @return Value of the min_ybound_type option.
	*/
	const std::string& Options::min_ybound_type() const 
		{return m_min_ybound_type;}

	/**
	* Returns the max_ybound_type option.
	*
	* Boundary condition at the maximum y boundary.
	*
	* @return Value of the max_ybound_type option.
	*/
	const std::string& Options::max_ybound_type() const 
		{return m_max_ybound_type;}

	/**
	* Returns the min_zbound_type option.
	*
	* Boundary condition at the minimum z boundary.
	*
	* @return Value of the min_zbound_type option.
	*/
	const std::string& Options::min_zbound_type() const 
		{return m_min_zbound_type;}

	/**
	* Returns the max_zbound_type option.
	*
	* Boundary condition at the maximum z boundary.
	*
	* @return Value of the max_zbound_type option.
	*/
	const std::string& Options::max_zbound_type() const 
		{return m_max_zbound_type;}

	/**
	* Returns the lcfs_x option.
	*
	* x coordinate of the last closed flux surface.
	*
	* @return Value of the lcfs_x option.
	*/
	double Options::lcfs_x() const 
		{return m_lcfs_x;}

	/**
	* Returns the lcfs_y option.
	*
	* y coordinate of the last closed flux surface.
	*
	* @return Value of the lcfs_y option.
	*/
	double Options::lcfs_y() const 
		{return m_lcfs_y;}

	/**
	* Returns the lcfs_z option.
	*
	* z coordinate of the last closed flux surface.
	*
	* @return Value of the lcfs_z option.
	*/
	double Options::lcfs_z() const 
		{return m_lcfs_z;}

	/**
	* Returns the sep_x_bc_xp_z1 option.
	*
	* First z value for the separatrix x boundary condition near the
	* X-point.
	*
	* @return Value of the sep_x_bc_xp_z1 option.
	*/
	double Options::sep_x_bc_xp_z1() const 
		{return m_sep_x_bc_xp_z1;}

	/**
	* Returns the sep_x_bc_xp_z2 option.
	*
	* Second z value for the separatrix x boundary condition near the
	* X-point.
	*
	* @return Value of the sep_x_bc_xp_z2 option.
	*/
	double Options::sep_x_bc_xp_z2() const 
		{return m_sep_x_bc_xp_z2;}

	/**
	* Returns the imp_xbound_buffer option.
	*
	* Buffer distance at the x boundaries. With an absorbing boundary, a
	* particle is absorbed once it is within this distance of the boundary.
	*
	* @return Value of the imp_xbound_buffer option.
	*/
	double Options::imp_xbound_buffer() const 
		{return m_imp_xbound_buffer;}

	/**
	* Returns the imp_ybound_buffer option.
	*
	* Buffer distance at the y boundaries. With an absorbing boundary, a
	* particle is absorbed once it is within this distance of the boundary.
	*
	* @return Value of the imp_ybound_buffer option.
	*/
	double Options::imp_ybound_buffer() const 
		{return m_imp_ybound_buffer;}

	/**
	* Returns the imp_zbound_buffer option.
	*
	* Buffer distance at the z boundaries. With an absorbing boundary, a
	* particle is absorbed once it is within this distance of the boundary.
	*
	* @return Value of the imp_zbound_buffer option.
	*/
	double Options::imp_zbound_buffer() const 
		{return m_imp_zbound_buffer;}

	/**
	* Returns the imp_atom_num option.
	*
	* Atomic number of the impurity.
	*
	* @return Value of the imp_atom_num option.
	*/
	int Options::imp_atom_num() const 
		{return m_imp_atom_num;}

	/**
	* Returns the imp_iz_max_charge option.
	*
	* @return Value of the imp_iz_max_charge option (below 0: no limit).
	*/
	int Options::imp_iz_max_charge() const
		{return m_imp_iz_max_charge;}

	/**
	* Returns the imp_mass_amu option.
	*
	* Mass of the impurity (amu).
	*
	* @return Value of the imp_mass_amu option.
	*/
	double Options::imp_mass_amu() const 
		{return m_imp_mass_amu;}

	/**
	* Returns the imp_init_charge option.
	*
	* Charge state the impurities start with.
	*
	* @return Value of the imp_init_charge option.
	*/
	int Options::imp_init_charge() const 
		{return m_imp_init_charge;}

	/**
	* Returns the imp_num option.
	*
	* Total number of impurity particles to follow.
	*
	* @return Value of the imp_num option.
	*/
	int Options::imp_num() const
		{return m_imp_num;}

	/**
	* Returns the imp_tstart_opt option.
	*
	* How each impurity's starting time is chosen: a single value
	* (imp_tstart_val), uniformly between imp_trange_min and imp_trange_max,
	* or uniformly over the full background time range.
	*
	* @return Value of the imp_tstart_opt option.
	*/
	const std::string& Options::imp_tstart_opt() const 
		{return m_imp_tstart_opt;}

	/**
	* Returns the imp_tstart_val option.
	*
	* Starting time of every impurity when imp_tstart_opt is "single_value"
	* (s).
	*
	* @return Value of the imp_tstart_val option.
	*/
	double Options::imp_tstart_val() const 
		{return m_imp_tstart_val;}

	/**
	* Returns the imp_trange_min option.
	*
	* Lower limit of the starting time range when imp_tstart_opt is "range"
	* (s).
	*
	* @return Value of the imp_trange_min option.
	*/
	double Options::imp_trange_min() const 
		{return m_imp_trange_min;}

	/**
	* Returns the imp_trange_max option.
	*
	* Upper limit of the starting time range when imp_tstart_opt is "range"
	* (s).
	*
	* @return Value of the imp_trange_max option.
	*/
	double Options::imp_trange_max() const 
		{return m_imp_trange_max;}

	/**
	* Returns the imp_xstart_opt option.
	*
	* How each impurity's starting x coordinate is chosen: a single value
	* (imp_xstart_val), uniformly between imp_xrange_min and imp_xrange_max,
	* or uniformly over the full x range of the background.
	*
	* @return Value of the imp_xstart_opt option.
	*/
	const std::string& Options::imp_xstart_opt() const 
		{return m_imp_xstart_opt;}

	/**
	* Returns the imp_xstart_val option.
	*
	* Starting x coordinate of every impurity when imp_xstart_opt is
	* "single_value".
	*
	* @return Value of the imp_xstart_val option.
	*/
	double Options::imp_xstart_val() const 
		{return m_imp_xstart_val;}

	/**
	* Returns the imp_xrange_min option.
	*
	* Lower limit of the starting x range when imp_xstart_opt is "range".
	*
	* @return Value of the imp_xrange_min option.
	*/
	double Options::imp_xrange_min() const 
		{return m_imp_xrange_min;}

	/**
	* Returns the imp_xrange_max option.
	*
	* Upper limit of the starting x range when imp_xstart_opt is "range".
	*
	* @return Value of the imp_xrange_max option.
	*/
	double Options::imp_xrange_max() const 
		{return m_imp_xrange_max;}

	/**
	* Returns the imp_ystart_opt option.
	*
	* How each impurity's starting y coordinate is chosen: a single value
	* (imp_ystart_val), uniformly between imp_yrange_min and imp_yrange_max,
	* or uniformly over the full y range of the background.
	*
	* @return Value of the imp_ystart_opt option.
	*/
	const std::string& Options::imp_ystart_opt() const 
		{return m_imp_ystart_opt;}

	/**
	* Returns the imp_ystart_val option.
	*
	* Starting y coordinate of every impurity when imp_ystart_opt is
	* "single_value".
	*
	* @return Value of the imp_ystart_val option.
	*/
	double Options::imp_ystart_val() const 
		{return m_imp_ystart_val;}

	/**
	* Returns the imp_yrange_min option.
	*
	* Lower limit of the starting y range when imp_ystart_opt is "range".
	*
	* @return Value of the imp_yrange_min option.
	*/
	double Options::imp_yrange_min() const 
		{return m_imp_yrange_min;}

	/**
	* Returns the imp_yrange_max option.
	*
	* Upper limit of the starting y range when imp_ystart_opt is "range".
	*
	* @return Value of the imp_yrange_max option.
	*/
	double Options::imp_yrange_max() const 
		{return m_imp_yrange_max;}

	/**
	* Returns the imp_zstart_opt option.
	*
	* How each impurity's starting z coordinate is chosen: a single value
	* (imp_zstart_val), uniformly between imp_zrange_min and imp_zrange_max,
	* or uniformly over the full z range of the background.
	*
	* @return Value of the imp_zstart_opt option.
	*/
	const std::string& Options::imp_zstart_opt() const 
		{return m_imp_zstart_opt;}

	/**
	* Returns the imp_zstart_val option.
	*
	* Starting z coordinate of every impurity when imp_zstart_opt is
	* "single_value".
	*
	* @return Value of the imp_zstart_val option.
	*/
	double Options::imp_zstart_val() const 
		{return m_imp_zstart_val;}

	/**
	* Returns the imp_zone_start_val option.
	*
	* Zone every impurity starts in, numbered from 1 like the SOLEDGE3X
	* zones.
	*
	* @return Value of the imp_zone_start_val option.
	*/
	int Options::imp_zone_start_val() const 
		{return m_imp_zone_start_val;}

	/**
	* Returns the imp_zrange_min option.
	*
	* Lower limit of the starting z range when imp_zstart_opt is "range".
	*
	* @return Value of the imp_zrange_min option.
	*/
	double Options::imp_zrange_min() const 
		{return m_imp_zrange_min;}

	/**
	* Returns the imp_zrange_max option.
	*
	* Upper limit of the starting z range when imp_zstart_opt is "range".
	*
	* @return Value of the imp_zrange_max option.
	*/
	double Options::imp_zrange_max() const 
		{return m_imp_zrange_max;}

	/**
	* Returns the imp_temp_start_opt option.
	*
	* How each impurity's starting temperature is chosen: a single value
	* (imp_temp_start_val) or the local main ion temperature. Starting
	* velocities are sampled from an isotropic Maxwellian at this
	* temperature.
	*
	* @return Value of the imp_temp_start_opt option.
	*/
	const std::string& Options::imp_temp_start_opt() const 
		{return m_imp_temp_start_opt;}

	/**
	* Returns the imp_temp_start_val option.
	*
	* Starting temperature of every impurity when imp_temp_start_opt is
	* "single_value" (eV).
	*
	* @return Value of the imp_temp_start_val option.
	*/
	double Options::imp_temp_start_val() const 
		{return m_imp_temp_start_val;}

	/**
	* Returns the imp_vX_start_val option (impurity starting velocity in X, added to the thermal one).
	*
	* @return Value of the imp_vX_start_val option.
	*/
	double Options::imp_vX_start_val() const
		{return m_imp_vX_start_val;}

	/**
	* Returns the imp_vY_start_val option (impurity starting velocity in Y, added to the thermal one).
	*
	* @return Value of the imp_vY_start_val option.
	*/
	double Options::imp_vY_start_val() const
		{return m_imp_vY_start_val;}

	/**
	* Returns the imp_vZ_start_val option (impurity starting velocity in Z, added to the thermal one).
	*
	* @return Value of the imp_vZ_start_val option.
	*/
	double Options::imp_vZ_start_val() const
		{return m_imp_vZ_start_val;}

	/**
	* Returns the imp_collisions option.
	*
	* Collision model for the impurities, either off or the Nanbu model.
	*
	* @return Value of the imp_collisions option.
	*/
	const std::string& Options::imp_collisions() const 
		{return m_imp_collisions;}

	/**
	* Returns the imp_thermal_force option, whether the Nanbu collisions
	* include the thermal force ("on") or not ("off").
	*
	* @return Value of the imp_thermal_force option.
	*/
	const std::string& Options::imp_thermal_force() const
		{return m_imp_thermal_force;}

	/**
	* Returns the imp_coll_subcycle_target option, the largest s0 a
	* sub-step of a Nanbu collision is meant to have.
	*
	* @return Value of the imp_coll_subcycle_target option.
	*/
	double Options::imp_coll_subcycle_target() const
		{return m_imp_coll_subcycle_target;}

	/**
	* Returns the imp_coll_subcycle_cap option, the largest number of
	* sub-steps one Nanbu collision is split into (1 = no sub-cycling).
	*
	* @return Value of the imp_coll_subcycle_cap option.
	*/
	int Options::imp_coll_subcycle_cap() const
		{return m_imp_coll_subcycle_cap;}

	/**
	* Returns the var_red_split option.
	*
	* When particles are split for variance reduction: never, at
	* ionization/recombination events, or at collisions.
	*
	* @return Value of the var_red_split option.
	*/
	const std::string& Options::var_red_split() const 
		{return m_var_red_split;}

	/**
	* Returns the var_red_import option.
	*
	* Method used to calculate particle importance for variance reduction.
	*
	* @return Value of the var_red_import option.
	*/
	const std::string& Options::var_red_import() const 
		{return m_var_red_import;}

	/**
	* Returns the var_red_freq option.
	*
	* How often variance reduction is applied.
	*
	* @return Value of the var_red_freq option.
	*/
	double Options::var_red_freq() const 
		{return m_var_red_freq;}

	/**
	* Returns the var_red_min_weight option.
	*
	* Minimum particle weight for variance reduction.
	*
	* @return Value of the var_red_min_weight option.
	*/
	double Options::var_red_min_weight() const 
		{return m_var_red_min_weight;}

	/**
	* Returns the var_red_med_mod option.
	*
	* Modifier applied to the median importance for variance reduction.
	*
	* @return Value of the var_red_med_mod option.
	*/
	double Options::var_red_med_mod() const 
		{return m_var_red_med_mod;}

	/**
	* Returns the var_red_rusrol option.
	*
	* Whether to use Russian roulette for variance reduction.
	*
	* @return Value of the var_red_rusrol option.
	*/
	const std::string& Options::var_red_rusrol() const 
		{return m_var_red_rusrol;}

	/**
	* Returns the var_red_rusrol_prob option.
	*
	* Probability used for Russian roulette.
	*
	* @return Value of the var_red_rusrol_prob option.
	*/
	double Options::var_red_rusrol_prob() const 
		{return m_var_red_rusrol_prob;}

	/**
	* Returns the imp_time_step_opt option.
	*
	* Whether the impurity time step is constant or variable.
	*
	* @return Value of the imp_time_step_opt option.
	*/
	const std::string& Options::imp_time_step_opt() const 
		{return m_imp_time_step_opt;}

	/**
	* Returns the imp_time_step option.
	*
	* Impurity time step (s).
	*
	* @return Value of the imp_time_step option.
	*/
	double Options::imp_time_step() const 
		{return m_imp_time_step;}

	/**
	* Returns the imp_time_step_min option.
	*
	* Minimum impurity time step when the time step is variable (s).
	*
	* @return Value of the imp_time_step_min option.
	*/
	double Options::imp_time_step_min() const 
		{return m_imp_time_step_min;}

	/**
	* Returns the imp_source_scale_fact option.
	*
	* Impurity source scaling factor (particles/s), i.e., how many real
	* particles a single simulated particle represents. Used to convert the
	* Monte Carlo counts into a density.
	*
	* @return Value of the imp_source_scale_fact option.
	*/
	double Options::imp_source_scale_fact() const 
		{return m_imp_source_scale_fact;}

	/**
	* Returns the imp_iz_recomb option.
	*
	* Whether impurity ionization and recombination are on.
	*
	* @return Value of the imp_iz_recomb option.
	*/
	const std::string& Options::imp_iz_recomb() const 
		{return m_imp_iz_recomb;}

	/**
	* Returns the imp_max_time option.
	*
	* Maximum time a particle is followed for (s). A particle that exceeds
	* it is killed and counted in a warning at the end of the run.
	*
	* @return Value of the imp_max_time option.
	*/
	double Options::imp_max_time() const 
		{return m_imp_max_time;}

	/**
	* Returns the print_interval option.
	*
	* Number of progress updates printed while particles are followed. For
	* example, 10 prints an update every 10% of particles.
	*
	* @return Value of the print_interval option.
	*/
	int Options::print_interval() const 
		{return m_print_interval;}

	/**
	* Returns the openadas_root option.
	*
	* Directory containing the OpenADAS rate coefficient files.
	*
	* @return Value of the openadas_root option.
	*/
	const std::string& Options::openadas_root() const 
		{return m_openadas_root;}

	/**
	* Returns the openadas_year option.
	*
	* Year of the OpenADAS data files to use.
	*
	* @return Value of the openadas_year option.
	*/
	int Options::openadas_year() const 
		{return m_openadas_year;}
	
	/**
	* Returns a control integer after checking that it was set.
	*
	* Control integers start at -1, meaning not set. If control_int is still
	* negative, the run is stopped with Errors::fatal_error. This usually
	* means the option's setter is missing from initialize_control_ints.
	*
	* @param var_name Name of the option the control integer belongs to (for
	*                 the error message).
	* @param control_int Control integer to check.
	* @return The control integer, unchanged.
	*/
	int Options::get_control_int(const std::string& var_name, 
		const int control_int) const
	{
		if (control_int < 0)
		{
			// If you encounter this error, then it probably means the control
			// integer was not being set at the beginning within 
			// initialize_control_ints. Easy mistake to make.
			Errors::fatal_error("Control integer for " + var_name 
				+ " was never set! This is a programming error, contact a "
				"developer. See initialize_control_ints for details.");
		}
		return control_int;
	}

	// Accessors for internal control variables

	/**
	* Returns the control integer for the use_gpu option.
	*
	* The control integer is an integer version of use_gpu that is faster to
	* check than the string: 0 = "off", 1 = "cuda". It is fetched through
	* get_control_int, which prints an error if it was never set.
	*
	* @return Control integer for use_gpu, or -1 if it was never set.
	*/
	int Options::use_gpu_int() const
		{return get_control_int("use_gpu", m_use_gpu_int);}

	/**
	* Returns the control integer for the bkg_source option.
	*
	* The control integer is an integer version of bkg_source that is faster
	* to check than the string: 0 = "test", 1 = "gkeyll". It is fetched
	* through get_control_int, which prints an error if it was never set.
	*
	* @return Control integer for bkg_source, or -1 if it was never set.
	*/
	int Options::bkg_source_int() const
		{return get_control_int("bkg_source", m_bkg_source_int);}


	/**
	* Returns the control integer for the save_track option.
	*
	* The control integer is an integer version of save_track that is faster
	* to check than the string: 0 = "off", 1 = "on". It is fetched through
	* get_control_int, which prints an error if it was never set.
	*
	* @return Control integer for save_track, or -1 if it was never set.
	*/
	int Options::save_track_int() const
		{return get_control_int("save_track", m_save_track_int);}

	/**
	* Returns the control integer for the imp_tstart_opt option.
	*
	* The control integer is an integer version of imp_tstart_opt that is
	* faster to check than the string: 0 = "single_value", 1 = "range", 2 =
	* "full_range". It is fetched through get_control_int, which prints an
	* error if it was never set.
	*
	* @return Control integer for imp_tstart_opt, or -1 if it was never set.
	*/
	int Options::imp_tstart_opt_int() const
		{return get_control_int("imp_tstart_opt", m_imp_tstart_opt_int);}

	/**
	* Returns the control integer for the imp_xstart_opt option.
	*
	* The control integer is an integer version of imp_xstart_opt that is
	* faster to check than the string: 0 = "single_value", 1 = "range", 2 =
	* "full_range". It is fetched through get_control_int, which prints an
	* error if it was never set.
	*
	* @return Control integer for imp_xstart_opt, or -1 if it was never set.
	*/
	int Options::imp_xstart_opt_int() const
		{return get_control_int("imp_xstart_opt", m_imp_xstart_opt_int);}

	/**
	* Returns the control integer for the imp_ystart_opt option.
	*
	* The control integer is an integer version of imp_ystart_opt that is
	* faster to check than the string: 0 = "single_value", 1 = "range", 2 =
	* "full_range". It is fetched through get_control_int, which prints an
	* error if it was never set.
	*
	* @return Control integer for imp_ystart_opt, or -1 if it was never set.
	*/
	int Options::imp_ystart_opt_int() const
		{return get_control_int("imp_ystart_opt", m_imp_ystart_opt_int);}

	/**
	* Returns the control integer for the imp_zstart_opt option.
	*
	* The control integer is an integer version of imp_zstart_opt that is
	* faster to check than the string: 0 = "single_value", 1 = "range", 2 =
	* "full_range". It is fetched through get_control_int, which prints an
	* error if it was never set.
	*
	* @return Control integer for imp_zstart_opt, or -1 if it was never set.
	*/
	int Options::imp_zstart_opt_int() const
		{return get_control_int("imp_zstart_opt", m_imp_zstart_opt_int);}

	/**
	* Returns the control integer for the imp_temp_start_opt option.
	*
	* The control integer is an integer version of imp_temp_start_opt that
	* is faster to check than the string: 0 = "single_value", 1 =
	* "main_ion". It is fetched through get_control_int, which prints an
	* error if it was never set.
	*
	* @return Control integer for imp_temp_start_opt, or -1 if it was never
	*         set.
	*/
	int Options::imp_temp_start_opt_int() const
		{return get_control_int("imp_temp_start_opt", m_imp_temp_start_opt_int);}

	/**
	* Returns the control integer for the imp_collisions option.
	*
	* The control integer is an integer version of imp_collisions that is
	* faster to check than the string: 0 = "off", 1 = "nanbu". It is fetched
	* through get_control_int, which prints an error if it was never set.
	*
	* @return Control integer for imp_collisions, or -1 if it was never set.
	*/
	int Options::imp_collisions_int() const
		{return get_control_int("imp_collisions", m_imp_collisions_int);}

	/**
	* Returns the control integer for the imp_thermal_force option: 0 =
	* "off", 1 = "on". It is fetched through get_control_int, which prints
	* an error if it was never set.
	*
	* @return Control integer for imp_thermal_force, or -1 if it was never
	*         set.
	*/
	int Options::imp_thermal_force_int() const
	{
		return get_control_int("imp_thermal_force",
			m_imp_thermal_force_int);
	}

	/**
	* Returns the control integer for the var_red_split option.
	*
	* The control integer is an integer version of var_red_split that is
	* faster to check than the string: 0 = "off", 1 = "iz_rec", 2 = "coll".
	* It is fetched through get_control_int, which prints an error if it was
	* never set.
	*
	* @return Control integer for var_red_split, or -1 if it was never set.
	*/
	int Options::var_red_split_int() const
		{return get_control_int("var_red_split", m_var_red_split_int);}

	/**
	* Returns the control integer for the var_red_import option.
	*
	* The control integer is an integer version of var_red_import that is
	* faster to check than the string: 0 = "median", 1 = "exp_dist", 2 =
	* "exp_time". It is fetched through get_control_int, which prints an
	* error if it was never set.
	*
	* @return Control integer for var_red_import, or -1 if it was never set.
	*/
	int Options::var_red_import_int() const
		{return get_control_int("var_red_import", m_var_red_import_int);}

	/**
	* Returns the control integer for the var_red_rusrol option.
	*
	* The control integer is an integer version of var_red_rusrol that is
	* faster to check than the string: 0 = "off", 1 = "on". It is fetched
	* through get_control_int, which prints an error if it was never set.
	*
	* @return Control integer for var_red_rusrol, or -1 if it was never set.
	*/
	int Options::var_red_rusrol_int() const
		{return get_control_int("var_red_rusrol", m_var_red_rusrol_int);}

	/**
	* Returns the control integer for the imp_time_step_opt option.
	*
	* The control integer is an integer version of imp_time_step_opt that is
	* faster to check than the string: 0 = "constant", 1 = "variable". It is
	* fetched through get_control_int, which prints an error if it was never
	* set.
	*
	* @return Control integer for imp_time_step_opt, or -1 if it was never
	*         set.
	*/
	int Options::imp_time_step_opt_int() const
		{return get_control_int("imp_time_step_opt", m_imp_time_step_opt_int);}

	/**
	* Returns the control integer for the imp_iz_recomb option.
	*
	* The control integer is an integer version of imp_iz_recomb that is
	* faster to check than the string: 0 = "off", 1 = "on". It is fetched
	* through get_control_int, which prints an error if it was never set.
	*
	* @return Control integer for imp_iz_recomb, or -1 if it was never set.
	*/
	int Options::imp_iz_recomb_int() const
		{return get_control_int("imp_iz_recomb", m_imp_iz_recomb_int);}

	/**
	* Returns the control integer for the tbound_type option.
	*
	* The control integer is an integer version of tbound_type that is
	* faster to check than the string: 0 = "absorbing", 1 = "periodic". It
	* is fetched through get_control_int, which prints an error if it was
	* never set.
	*
	* @return Control integer for tbound_type, or -1 if it was never set.
	*/
	int Options::tbound_type_int() const
		{return get_control_int("tbound_type", m_tbound_type_int);}

	/**
	* Returns the control integer for the min_xbound_type option.
	*
	* The control integer is an integer version of min_xbound_type that is
	* faster to check than the string: 0 = "absorbing", 1 = "periodic", 2 =
	* "core". It is fetched through get_control_int, which prints an error
	* if it was never set.
	*
	* @return Control integer for min_xbound_type, or -1 if it was never
	*         set.
	*/
	int Options::min_xbound_type_int() const
		{return get_control_int("min_xbound_type", m_min_xbound_type_int);}

	/**
	* Returns the control integer for the max_xbound_type option.
	*
	* The control integer is an integer version of max_xbound_type that is
	* faster to check than the string: 0 = "absorbing", 1 = "periodic", 2 =
	* "core". It is fetched through get_control_int, which prints an error
	* if it was never set.
	*
	* @return Control integer for max_xbound_type, or -1 if it was never
	*         set.
	*/
	int Options::max_xbound_type_int() const
		{return get_control_int("max_xbound_type", m_max_xbound_type_int);}

	/**
	* Returns the control integer for the min_ybound_type option.
	*
	* The control integer is an integer version of min_ybound_type that is
	* faster to check than the string: 0 = "absorbing", 1 = "periodic", 2 =
	* "core". It is fetched through get_control_int, which prints an error
	* if it was never set.
	*
	* @return Control integer for min_ybound_type, or -1 if it was never
	*         set.
	*/
	int Options::min_ybound_type_int() const
		{return get_control_int("min_ybound_type", m_min_ybound_type_int);}

	/**
	* Returns the control integer for the max_ybound_type option.
	*
	* The control integer is an integer version of max_ybound_type that is
	* faster to check than the string: 0 = "absorbing", 1 = "periodic", 2 =
	* "core". It is fetched through get_control_int, which prints an error
	* if it was never set.
	*
	* @return Control integer for max_ybound_type, or -1 if it was never
	*         set.
	*/
	int Options::max_ybound_type_int() const
		{return get_control_int("max_ybound_type", m_max_ybound_type_int);}

	/**
	* Returns the control integer for the min_zbound_type option.
	*
	* The control integer is an integer version of min_zbound_type that is
	* faster to check than the string: 0 = "absorbing", 1 = "periodic", 2 =
	* "core". It is fetched through get_control_int, which prints an error
	* if it was never set.
	*
	* @return Control integer for min_zbound_type, or -1 if it was never
	*         set.
	*/
	int Options::min_zbound_type_int() const
		{return get_control_int("min_zbound_type", m_min_zbound_type_int);}

	/**
	* Returns the control integer for the max_zbound_type option.
	*
	* The control integer is an integer version of max_zbound_type that is
	* faster to check than the string: 0 = "absorbing", 1 = "periodic", 2 =
	* "core". It is fetched through get_control_int, which prints an error
	* if it was never set.
	*
	* @return Control integer for max_zbound_type, or -1 if it was never
	*         set.
	*/
	int Options::max_zbound_type_int() const
		{return get_control_int("max_zbound_type", m_max_zbound_type_int);}

	/**
	* Returns the control integer for the calc_grad_elec option.
	*
	* The control integer is an integer version of calc_grad_elec that is
	* faster to check than the string: 0 = "off", 1 = "on". It is fetched
	* through get_control_int, which prints an error if it was never set.
	*
	* @return Control integer for calc_grad_elec, or -1 if it was never set.
	*/
	int Options::calc_grad_elec_int() const
		{return get_control_int("calc_grad_elec", m_calc_grad_elec_int);}

	int Options::s3x_chi_absorb_int() const
		{return get_control_int("s3x_chi_absorb", m_s3x_chi_absorb_int);}

	int Options::hermes_format_int() const
		{return get_control_int("hermes_format", m_hermes_format_int);}

	int Options::hermes_heat_flux_int() const
		{return get_control_int("hermes_heat_flux", m_hermes_heat_flux_int);}

	int Options::hermes_inductive_epar_int() const
	{
		return get_control_int("hermes_inductive_epar",
			m_hermes_inductive_epar_int);
	}

	int Options::hermes_field_aligned_int() const
	{
		return get_control_int("hermes_field_aligned",
			m_hermes_field_aligned_int);
	}

	int Options::imp_source_opt_int() const
		{return get_control_int("imp_source_opt", m_imp_source_opt_int);}

	int Options::sput_flux_opt_int() const
		{return get_control_int("sput_flux_opt", m_sput_flux_opt_int);}

	int Options::sput_bohm_floor_int() const
		{return get_control_int("sput_bohm_floor", m_sput_bohm_floor_int);}

	int Options::sheath_model_int() const
		{return get_control_int("sheath_model", m_sheath_model_int);}

	/**
	* Returns the control integer for sput_self: 0 for "off", 1 for "on".
	*
	* @return Control integer of the sput_self option.
	*/
	int Options::sput_self_int() const
		{return get_control_int("sput_self", m_sput_self_int);}

	/**
	* Returns the control integer for sput_impact_opt: 0 for "formula", 1
	* for "sheath".
	*
	* @return Control integer of the sput_impact_opt option.
	*/
	int Options::sput_impact_opt_int() const
		{return get_control_int("sput_impact_opt", m_sput_impact_opt_int);}

	// Setter declarations for internal control variables

	/**
	* Sets the control integer for the var_red_split option directly.
	*
	* Unlike the other control integers, this one can be changed without
	* changing the var_red_split string option.
	*
	* @param var_red_split_int New value of the control integer.
	*/
	void Options::set_var_red_split_int(int var_red_split_int)
		{m_var_red_split_int = var_red_split_int;}
	
#ifndef USE_CUDA

	// If these definitions gets called, it means you're in the CPU-only version
	// of the code trying to call GPU-specific code.
	// Defined in src/cuda/options.cu since it calls CUDA code

	/**
	* Placeholder for copying the options to a GPU in CPU-only builds.
	*
	* The real definition is in src/cuda/options.cu. This version only
	* prints an error saying GPU support was not compiled in.
	*
	* @param device_id ID of the GPU to copy to (unused).
	* @return nullptr, since there is no device copy. Without a return
	* the call would be undefined behaviour, not just an error message.
	*/
	OptionsDevice* Options::to_device(int device_id)
	{
		std::cerr << "Error! Options::to_device was called but GPU support"
				  << " was not compiled in.\n";
		return nullptr;
	}

	/**
	* Placeholder for freeing device-side options in CPU-only builds.
	*
	* The real definition is in src/cuda/options.cu. This version only
	* prints an error saying GPU support was not compiled in.
	*
	* @param opts_d Device-side options to free (unused).
	* @param device_id ID of the GPU the options are on (unused).
	*/
	void free_opts(OptionsDevice* opts_d, int device_id)
	{
		std::cerr << "Error! Options::free_opts was called but GPU support"
				  << " was not compiled in.\n";
	}

#endif

	/**
	* Checks a value against a list of valid values for an option.
	*
	* If the value is not in the list, the run is stopped with
	* Errors::fatal_error, naming the option, the value and the valid
	* values. Running on with the default instead would silently ignore
	* what the user asked for.
	*
	* @tparam T Type of the option value.
	* @param var Name of the option (for the error message).
	* @param value Value to check.
	* @param valid_values List of valid values, e.g., {"on", "off"}.
	* @return True if the value is valid. Doesn't return otherwise.
	*/
	template <typename T>
	bool Options::check_input(const std::string& var, const T& value, 
		std::initializer_list<T> valid_values)
	{
		// Check input value against allowed values
		for (const auto& v : valid_values)
		{
			if (value == v) 
			{
				return true;
			}
		}

		// If allowed value not found, stop the run, listing the valid values
		std::ostringstream msg {};
		msg << "Value for " << var << " not recognized: " << value 
			<< ". Valid values are:";
		for (const auto& v : valid_values) msg << ' ' << v;
		Errors::fatal_error(msg.str());
	}
}
