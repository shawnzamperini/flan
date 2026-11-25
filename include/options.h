#pragma once

/**
* @file options.h
* @brief Header file for options.cpp.
*/
#include <string>
#include <tuple>

#include "flan_types.h"
#include "options_device.h"


namespace Options
{

	/**
	* @brief Class for storing all the options that control a simulation.
	*
	* Default values for every option are included in the definitions of each
	* member variable. These can then be overwritten as in @ref 
	* load_options.cpp
	*/
	class Options
	{
	private:

		// Input options and their defaults are defined here

		// Meta options for the simulation
		std::string m_case_name              {"undefined"};
		std::string m_use_gpu					   {"off"};
		int m_seed                                    {-1};
		int m_slot_cap                           {1048576};  
		int m_print_interval                          {10};

		// General background plasma options
		std::string m_bkg_source                {"gkeyll"};

		// Internal test control options
		std::string m_save_track                   {"off"};

		// Options related to reading in Gkeyll files
		std::string m_gkyl_dir               {"undefined"};
		std::string m_gkyl_casename          {"undefined"};
		int m_gkyl_frame_start                         {0};
		int m_gkyl_frame_end                           {1};
		std::string m_gkyl_elec_name               {"elc"};
		std::string m_gkyl_ion_name                {"ion"};
		double m_gkyl_elec_mass_amu             {0.000548};
		double m_gkyl_ion_mass_amu                 {2.014};
		std::string m_gkyl_file_type            {"binary"};
		std::string m_gkyl_moment_type    {"bimaxwellian"};
		std::string m_calc_grad_elec               {"off"};

		// Options related to reading in SOLEDGE3X files
		std::string m_s3x_dir                {"undefined"};
		std::string m_s3x_mesh_file        {"mesh_S3X.h5"};
		int m_s3x_frame_start                          {0};
		int m_s3x_frame_end                            {1};
		std::string m_s3x_ref_file {"refParam_raptorX.h5"};
		std::string m_s3x_chi_absorb                {"on"};

		// Options related to reading in HERMES (Hermes-3, JAX-Hermes3) runs
		std::string m_hermes_format {"jax"};
		std::string m_hermes_grid_file {"undefined"};
		std::string m_hermes_snap_dir {"undefined"};
		std::string m_hermes_ion_species {"i"};
		int m_hermes_snap_start {0};
		int m_hermes_snap_end {1};
		int m_hermes_snap_stride {1};
		double m_hermes_nnorm {1e19};
		double m_hermes_tnorm {100.0};
		double m_hermes_bnorm {1.0};
		double m_hermes_ion_mass_amu {2.0};
		int m_hermes_zperiod {1};
		std::string m_hermes_heat_flux {"conduction"};
		std::string m_hermes_field_aligned {"on"};
		std::string m_hermes_inductive_epar {"on"};
		double m_hermes_kappa_coef_e {3.16 / 1.4142135623730951};
		double m_hermes_kappa_coef_i {3.9};
		double m_hermes_kappa_limit_alpha {-1.0};

		// Options for sputtering of the walls by the background plasma
		std::string m_imp_source_opt            {"manual"};
		std::string m_sput_table_file        {"undefined"};
		std::string m_sput_flux_opt       {"parallel_exb"};
		std::string m_sput_bohm_floor              {"off"};
		std::string m_sheath_model                 {"brooks"};
		double m_sheath_fd                          {-1.0};
		double m_sheath_drop                        {-1.0};
		std::string m_sput_self                     {"off"};
		std::string m_sput_self_table_file    {"undefined"};
		double m_sput_self_min_weight               {0.1};
		int m_sput_self_max_gen                     {1000};
		double m_sput_angle_val                     {60.0};
		std::string m_sput_impact_opt               {"sheath"};
		double m_sput_max_wall_dist                  {2.0};
		
		// General background options
		double m_min_ne                             {1e15};
		double m_min_te                              {0.1};
		double m_min_ti                              {0.1};
		double m_slab_ne                              {1e19};
		double m_slab_te                              {10.0};
		double m_slab_ti                              {10.0};
		double m_slab_flow                            {0.0};
		double m_slab_heat_flux_i                     {0.0};
		double m_slab_heat_flux_e                     {0.0};
		double m_slab_bmag                            {1.0};
		double m_slab_tmax                            {50e-6};
		double m_slab_efield_y                        {0.0};
		double m_slab_efield_y_rate                   {0.0};
		double m_slab_bgrad_y                         {0.0};
		double m_slab_major_radius                    {0.0};
		std::string m_slab_wall                    {"none"};
		double m_slab_b_angle                         {0.0};

		// Geometry options
		std::string m_tbound_type            {"absorbing"};
		std::string m_min_xbound_type        {"absorbing"};
		std::string m_max_xbound_type        {"absorbing"};
		std::string m_min_ybound_type        {"absorbing"};
		std::string m_max_ybound_type        {"absorbing"};
		std::string m_min_zbound_type        {"absorbing"};
		std::string m_max_zbound_type        {"absorbing"};
		double m_lcfs_x                              {0.0};
		double m_lcfs_y                              {0.0};
		double m_lcfs_z                              {0.0};
		double m_imp_xbound_buffer                   {0.0};
		double m_imp_ybound_buffer                   {0.0};
		double m_imp_zbound_buffer                   {0.0};
		double m_sep_x_bc_xp_z1                      {0.0};
		double m_sep_x_bc_xp_z2                      {0.0};

		// Impurity characteristics
		int m_imp_atom_num                            {74};
		int m_imp_iz_max_charge                       {-1};
		double m_imp_mass_amu                     {183.84};
		int m_imp_init_charge                          {1};

		// Impurity transport options
		int m_imp_num                                  {1};
		std::string m_imp_tstart_opt      {"single_value"};
		double m_imp_tstart_val                      {0.0};
		double m_imp_trange_min                      {0.0};
		double m_imp_trange_max                      {0.0};
		std::string m_imp_xstart_opt      {"single_value"};
		double m_imp_xstart_val                      {0.0};
		double m_imp_xrange_min                      {0.0};
		double m_imp_xrange_max                      {0.0};
		std::string m_imp_ystart_opt      {"single_value"};
		double m_imp_ystart_val                      {0.0};
		double m_imp_yrange_min                      {0.0};
		double m_imp_yrange_max                      {0.0};
		std::string m_imp_zstart_opt      {"single_value"};
		double m_imp_zstart_val                      {0.0};
		int m_imp_zone_start_val                       {1};
		double m_imp_zrange_min                      {0.0};
		double m_imp_zrange_max                      {0.0};
		std::string m_imp_temp_start_opt      {"main_ion"};
		double m_imp_temp_start_val                  {0.0};
		double m_imp_vX_start_val                     {0.0};
		double m_imp_vY_start_val                     {0.0};
		double m_imp_vZ_start_val                     {0.0};
		std::string m_imp_collisions             {"nanbu"};
		std::string m_imp_thermal_force             {"on"};
		double m_imp_coll_subcycle_target           {0.03};
		int m_imp_coll_subcycle_cap                  {128};
		std::string m_imp_time_step_opt       {"constant"};
		double m_imp_time_step                      {1e-8};
		double m_imp_time_step_min                 {1e-12};
		double m_imp_source_scale_fact               {1.0};
		std::string m_imp_iz_recomb                 {"on"};
		double m_imp_max_time                       {1e10};

		// Variance reduction options
		std::string m_var_red_split                {"off"};
		std::string m_var_red_import            {"median"};
		double m_var_red_freq                        {0.1};
		double m_var_red_min_weight                  {0.1};
		double m_var_red_med_mod                     {1.0};
		std::string m_var_red_rusrol               {"off"};
		double m_var_red_rusrol_prob				 {0.5};

		// OpenADAS options
		std::string m_openadas_root          {"undefined"};
		int m_openadas_year                           {50};

		// Internal control variables for string options. The string options
		// are just there to be user-friendly, but it is more efficient to
		// use integers instead of comparing strings all the time. These are
		// all set automatically within the corresponding setters.
		// IMPORTANT: All control integers should be initialized below as -1,
		// which the code recognizes as "not assigned yet".
		int m_use_gpu_int	        	{-1};
		int m_bkg_source_int        	{-1};
		int m_save_track_int          	{-1};
		int m_imp_tstart_opt_int		{-1};
		int m_imp_xstart_opt_int		{-1};
		int m_imp_ystart_opt_int		{-1};
		int m_imp_zstart_opt_int		{-1};
		int m_imp_temp_start_opt_int	{-1};
		int m_imp_collisions_int		{-1};
		int m_imp_thermal_force_int		{-1};
		int m_var_red_split_int			{-1};
		int m_var_red_import_int		{-1};
		int m_var_red_rusrol_int		{-1};
		int m_imp_time_step_opt_int 	{-1};
		int m_imp_iz_recomb_int			{-1};
		int m_tbound_type_int			{-1};
		int m_min_xbound_type_int		{-1};
		int m_max_xbound_type_int		{-1};
		int m_min_ybound_type_int		{-1};
		int m_max_ybound_type_int		{-1};
		int m_min_zbound_type_int		{-1};
		int m_max_zbound_type_int		{-1};
		int m_calc_grad_elec_int		{-1};
		int m_s3x_chi_absorb_int                      {-1};
		int m_hermes_format_int                       {-1};
		int m_hermes_heat_flux_int                    {-1};
		int m_hermes_field_aligned_int                {-1};
		int m_hermes_inductive_epar_int               {-1};
		int m_imp_source_opt_int                      {-1};
		int m_sput_flux_opt_int                       {-1};
		int m_sput_bohm_floor_int                     {-1};
		int m_sheath_model_int                        {-1};
		int m_sput_self_int                           {-1};
		int m_sput_impact_opt_int                     {-1};

		template <typename T>
		bool check_input(const std::string& var, const T& value, 
			std::initializer_list<T> valid_values);

	public:

		// Constructor
		Options();

		// Initialize internal control integers (called during contructor)
		void initialize_control_ints();

		// Setter declarations
		void set_use_gpu(std::string use_gpu);
		void set_case_name(std::string case_name);
		void set_seed(int seed);
		void set_slot_cap(int slot_cap);
		void set_bkg_source(std::string bkg_source);
		void set_save_track(std::string save_track);
		void set_gkyl_dir(std::string gkyl_dir);
		void set_gkyl_casename(std::string gkyl_casename);
		void set_gkyl_frame_start(int gkyl_frame_start);
		void set_gkyl_frame_end(int gkyl_frame_end);
		void set_gkyl_elec_name(std::string gkyl_elec_name);
		void set_gkyl_ion_name(std::string gkyl_ion_name);
		void set_gkyl_elec_mass_amu(double gkyl_elec_mass_amu);
		void set_gkyl_ion_mass_amu(double gkyl_ion_mass_amu);
		void set_gkyl_file_type(std::string gkyl_file_type);
		void set_gkyl_moment_type(std::string gkyl_moment_type);
		void set_calc_grad_elec(std::string calc_grad_elec);
		void set_s3x_dir(std::string s3x_dir);
		void set_s3x_mesh_file(std::string s3x_mesh_file);
		void set_s3x_frame_start(int s3x_frame_start);
		void set_s3x_frame_end(int s3x_frame_end);
		void set_s3x_ref_file(std::string s3x_ref_file);
		void set_s3x_chi_absorb(std::string s3x_chi_absorb);
		void set_hermes_format(std::string hermes_format);
		void set_hermes_grid_file(std::string hermes_grid_file);
		void set_hermes_snap_dir(std::string hermes_snap_dir);
		void set_hermes_ion_species(std::string hermes_ion_species);
		void set_hermes_snap_start(int hermes_snap_start);
		void set_hermes_snap_end(int hermes_snap_end);
		void set_hermes_snap_stride(int hermes_snap_stride);
		void set_hermes_nnorm(double hermes_nnorm);
		void set_hermes_tnorm(double hermes_tnorm);
		void set_hermes_bnorm(double hermes_bnorm);
		void set_hermes_ion_mass_amu(double hermes_ion_mass_amu);
		void set_hermes_zperiod(int hermes_zperiod);
		void set_hermes_heat_flux(std::string hermes_heat_flux);
		void set_hermes_field_aligned(std::string hermes_field_aligned);
		void set_hermes_inductive_epar(std::string hermes_inductive_epar);
		void set_hermes_kappa_coef_e(double hermes_kappa_coef_e);
		void set_hermes_kappa_coef_i(double hermes_kappa_coef_i);
		void set_hermes_kappa_limit_alpha(double hermes_kappa_limit_alpha);
		void set_imp_source_opt(std::string imp_source_opt);
		void set_sput_table_file(std::string sput_table_file);
		void set_sput_flux_opt(std::string sput_flux_opt);
		void set_sput_bohm_floor(std::string sput_bohm_floor);
		void set_sheath_model(std::string sheath_model);
		void set_sheath_fd(double sheath_fd);
		void set_sheath_drop(double sheath_drop);
		void set_sput_self(std::string sput_self);
		void set_sput_impact_opt(std::string sput_impact_opt);
		void set_sput_self_table_file(std::string sput_self_table_file);
		void set_sput_self_min_weight(double sput_self_min_weight);
		void set_sput_self_max_gen(int sput_self_max_gen);
		void set_sput_angle_val(double sput_angle_val);
		void set_sput_max_wall_dist(double sput_max_wall_dist);
		void set_min_ne(double min_ne);
		void set_min_te(double min_te);
		void set_min_ti(double min_ti);
		void set_slab_ne(double slab_ne);
		void set_slab_te(double slab_te);
		void set_slab_ti(double slab_ti);
		void set_slab_flow(double slab_flow);
		void set_slab_heat_flux_i(double slab_heat_flux_i);
		void set_slab_heat_flux_e(double slab_heat_flux_e);
		void set_slab_bmag(double slab_bmag);
		void set_slab_tmax(double slab_tmax);
		void set_slab_efield_y(double slab_efield_y);
		void set_slab_efield_y_rate(double slab_efield_y_rate);
		void set_slab_bgrad_y(double slab_bgrad_y);
		void set_slab_major_radius(double slab_major_radius);
		void set_slab_wall(std::string slab_wall);
		void set_slab_b_angle(double slab_b_angle);
		void set_tbound_type(std::string tbound_type);
		void set_min_xbound_type(std::string min_xbound_type);
		void set_max_xbound_type(std::string max_xbound_type);
		void set_min_ybound_type(std::string min_ybound_type);
		void set_max_ybound_type(std::string max_ybound_type);
		void set_min_zbound_type(std::string min_zbound_type);
		void set_max_zbound_type(std::string max_zbound_type);
		void set_lcfs_x(double lcfs_x);
		void set_lcfs_y(double lcfs_y);
		void set_lcfs_z(double lcfs_z);
		void set_sep_x_bc_xp_z1(double sep_x_bc_xp_z1);
		void set_sep_x_bc_xp_z2(double sep_x_bc_xp_z2);
		void set_imp_xbound_buffer(double imp_xbound_buffer);
		void set_imp_ybound_buffer(double imp_ybound_buffer);
		void set_imp_zbound_buffer(double imp_zbound_buffer);
		void set_imp_atom_num(int imp_atom_num);
		void set_imp_iz_max_charge(int imp_iz_max_charge);
		void set_imp_mass_amu(double imp_mass_amu);
		void set_imp_init_charge(int imp_init_charge);
		void set_imp_num(int imp_num);
		void set_imp_tstart_opt(std::string imp_tstart_opt);
		void set_imp_tstart_val(double imp_tstart_val);
		void set_imp_trange_min(double imp_trange_min);
		void set_imp_trange_max(double imp_trange_max);
		void set_imp_xstart_opt(std::string imp_xstart_opt);
		void set_imp_xstart_val(double imp_xstart_val);
		void set_imp_xrange_min(double imp_xrange_min);
		void set_imp_xrange_max(double imp_xrange_max);
		void set_imp_ystart_opt(std::string imp_ystart_opt);
		void set_imp_ystart_val(double imp_ystart_val);
		void set_imp_yrange_min(double imp_yrange_min);
		void set_imp_yrange_max(double imp_yrange_max);
		void set_imp_zstart_opt(std::string imp_zstart_opt);
		void set_imp_zstart_val(double imp_zstart_val);
		void set_imp_zone_start_val(int imp_zone_start_val);
		void set_imp_zrange_min(double imp_zrange_min);
		void set_imp_zrange_max(double imp_zrange_max);
		void set_imp_temp_start_opt(std::string imp_temp_start_opt);
		void set_imp_temp_start_val(double imp_temp_start_val);
		void set_imp_vX_start_val(double imp_vX_start_val);
		void set_imp_vY_start_val(double imp_vY_start_val);
		void set_imp_vZ_start_val(double imp_vZ_start_val);
		void set_imp_collisions(std::string imp_collisions);
		void set_imp_thermal_force(std::string imp_thermal_force);
		void set_imp_coll_subcycle_target(double imp_coll_subcycle_target);
		void set_imp_coll_subcycle_cap(int imp_coll_subcycle_cap);
		void set_var_red_split(std::string var_red_split);
		void set_var_red_import(std::string var_red_import);
		void set_var_red_freq(double var_red_freq);
		void set_var_red_min_weight(double var_red_min_weight);
		void set_var_red_med_mod(double var_red_med_mod);
		void set_var_red_rusrol(std::string var_red_rusrol);
		void set_var_red_rusrol_prob(double var_red_rusrol_prob);
		void set_imp_time_step_opt(std::string imp_time_step_opt);
		void set_imp_time_step(double imp_time_step);
		void set_imp_time_step_min(double imp_time_step_min);
		void set_imp_source_scale_fact(double imp_source_scale_fact);
		void set_imp_iz_recomb(std::string imp_iz_recomb);
		void set_imp_max_time(double imp_max_time);
		void set_print_interval(int print_interval);
		void set_openadas_root(std::string openadas_root);
		void set_openadas_year(int openadas_year);

		// Accessor declarations
		const std::string& case_name() const;
		const std::string& use_gpu() const;
		int seed() const;
		int slot_cap() const;
		const std::string& bkg_source() const;
		const std::string& save_track() const;
		const std::string& gkyl_dir() const;
		const std::string& gkyl_casename() const;
		int gkyl_frame_start() const;
		int gkyl_frame_end() const;
		const std::string& gkyl_elec_name() const;
		const std::string& gkyl_ion_name() const;
		double gkyl_elec_mass_amu() const;
		double gkyl_ion_mass_amu() const;
		const std::string& gkyl_file_type() const;
		const std::string& gkyl_moment_type() const;
		const std::string& calc_grad_elec() const;
		const std::string& s3x_dir() const;
		const std::string& s3x_mesh_file() const;
		int s3x_frame_start() const;
		int s3x_frame_end() const;
		const std::string& s3x_ref_file() const;
		const std::string& s3x_chi_absorb() const;
		const std::string& hermes_format() const;
		const std::string& hermes_grid_file() const;
		const std::string& hermes_snap_dir() const;
		const std::string& hermes_ion_species() const;
		int hermes_snap_start() const;
		int hermes_snap_end() const;
		int hermes_snap_stride() const;
		double hermes_nnorm() const;
		double hermes_tnorm() const;
		double hermes_bnorm() const;
		double hermes_ion_mass_amu() const;
		int hermes_zperiod() const;
		const std::string& hermes_heat_flux() const;
		const std::string& hermes_field_aligned() const;
		const std::string& hermes_inductive_epar() const;
		double hermes_kappa_coef_e() const;
		double hermes_kappa_coef_i() const;
		double hermes_kappa_limit_alpha() const;
		const std::string& imp_source_opt() const;
		const std::string& sput_table_file() const;
		const std::string& sput_flux_opt() const;
		const std::string& sput_bohm_floor() const;
		const std::string& sheath_model() const;
		double sheath_fd() const;
		double sheath_drop() const;
		const std::string& sput_self() const;
		const std::string& sput_impact_opt() const;
		const std::string& sput_self_table_file() const;
		double sput_self_min_weight() const;
		int sput_self_max_gen() const;
		double sput_angle_val() const;
		double sput_max_wall_dist() const;
		double min_ne() const;
		double min_te() const;
		double min_ti() const;
		double slab_ne() const;
		double slab_te() const;
		double slab_ti() const;
		double slab_flow() const;
		double slab_heat_flux_i() const;
		double slab_heat_flux_e() const;
		double slab_bmag() const;
		double slab_tmax() const;
		double slab_efield_y() const;
		double slab_efield_y_rate() const;
		double slab_bgrad_y() const;
		double slab_major_radius() const;
		const std::string& slab_wall() const;
		double slab_b_angle() const;
		const std::string& tbound_type() const;
		const std::string& min_xbound_type() const;
		const std::string& max_xbound_type() const;
		const std::string& min_ybound_type() const;
		const std::string& max_ybound_type() const;
		const std::string& min_zbound_type() const;
		const std::string& max_zbound_type() const;
		double lcfs_x() const;
		double lcfs_y() const;
		double lcfs_z() const;
		double sep_x_bc_xp_z1() const;
		double sep_x_bc_xp_z2() const;
		double imp_xbound_buffer() const;
		double imp_ybound_buffer() const;
		double imp_zbound_buffer() const;
		int imp_atom_num() const;
		int imp_iz_max_charge() const;
		double imp_mass_amu() const;
		int imp_init_charge() const;
		int imp_num() const;
		const std::string& imp_tstart_opt() const;
		double imp_tstart_val() const;
		double imp_trange_min() const;
		double imp_trange_max() const;
		const std::string& imp_xstart_opt() const;
		double imp_xstart_val() const;
		double imp_xrange_min() const;
		double imp_xrange_max() const;
		const std::string& imp_ystart_opt() const;
		double imp_ystart_val() const;
		double imp_yrange_min() const;
		double imp_yrange_max() const;
		const std::string& imp_zstart_opt() const;
		double imp_zstart_val() const;
		int imp_zone_start_val() const;
		double imp_zrange_min() const;
		double imp_zrange_max() const;
		const std::string& imp_temp_start_opt() const;
		double imp_temp_start_val() const;
		double imp_vX_start_val() const;
		double imp_vY_start_val() const;
		double imp_vZ_start_val() const;
		const std::string& imp_collisions() const;
		const std::string& imp_thermal_force() const;
		double imp_coll_subcycle_target() const;
		int imp_coll_subcycle_cap() const;
		const std::string& var_red_split() const;
		const std::string& var_red_import() const;
		double var_red_freq() const;
		double var_red_min_weight() const;
		double var_red_med_mod() const;
		const std::string& var_red_rusrol() const;
		double var_red_rusrol_prob() const;
		const std::string& imp_time_step_opt() const;
		double imp_time_step() const;
		double imp_time_step_min() const;
		double imp_source_scale_fact() const;
		const std::string& imp_iz_recomb() const;
		double imp_max_time() const;
		int print_interval() const;
		const std::string& openadas_root() const;
		int openadas_year() const;

		// Wrapper function for retrieving control integers that ensures they
		// have been correctly set.
		int get_control_int(const std::string& var_name, 
			const int control_int) const;

		// Accessor declarations for internal control variables
		int use_gpu_int() const;
		int bkg_source_int() const;
		int save_track_int() const;
		int imp_tstart_opt_int() const;
		int imp_xstart_opt_int() const;
		int imp_ystart_opt_int() const;
		int imp_zstart_opt_int() const;
		int imp_temp_start_opt_int() const;
		int imp_collisions_int() const;
		int imp_thermal_force_int() const;
		int var_red_split_int() const;
		int var_red_import_int() const;
		int var_red_rusrol_int() const;
		int imp_time_step_opt_int() const;
		int imp_iz_recomb_int() const;
		int tbound_type_int() const;
		int min_xbound_type_int() const;
		int max_xbound_type_int() const;
		int min_ybound_type_int() const;
		int max_ybound_type_int() const;
		int min_zbound_type_int() const;
		int max_zbound_type_int() const;
		int calc_grad_elec_int() const;
		int s3x_chi_absorb_int() const;
		int hermes_format_int() const;
		int hermes_heat_flux_int() const;
		int hermes_field_aligned_int() const;
		int hermes_inductive_epar_int() const;
		int imp_source_opt_int() const;
		int sput_flux_opt_int() const;
		int sput_bohm_floor_int() const;
		int sheath_model_int() const;
		int sput_self_int() const;
		int sput_impact_opt_int() const;

		// Setter declarations for internal control variables. Added as needed.
		void set_var_red_split_int(int var_red_split_int);

		// Copy a select number of options to device
		OptionsDevice* to_device(int device_id = 0);
	};
	
	// Free options on device
	void free_opts(OptionsDevice* opts_d, int device_id = 0);
}
