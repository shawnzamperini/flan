/**
* @file read_gkyl.h
*
* @brief Header file for read_gkyl.cpp
*/

#ifndef READ_GKYL_H
#define READ_GKYL_H

#include <string>
#include <vector>

#include "background.h"
#include "flan_types.h"
#include "options.h"
#include "vectors.h"

namespace Gkyl
{
	// Alias for tuple to hold all the grid data to be passed around
	using grid_data_t = std::tuple<std::vector<BkgFPType>&, 
		std::vector<BkgFPType>&, std::vector<BkgFPType>&, 
		std::vector<BkgFPType>&, std::vector<BkgFPType>&, 
		std::vector<BkgFPType>&, std::vector<BkgFPType>&,
		Vectors::Vector3D<BkgFPType>&, Vectors::Vector3D<BkgFPType>&,
		Vectors::Vector3D<BkgFPType>&, Vectors::Vector3D<BkgFPType>&, 
		Vectors::Vector3D<BkgFPType>&, Vectors::Vector3D<BkgFPType>&>;

	Background::Background read_gkyl(const Options::Options& opts);
	
	std::string get_read_gkyl_py();

	std::vector<BkgFPType> load_times();

	std::tuple<std::vector<BkgFPType>, std::vector<BkgFPType>, 
		std::vector<BkgFPType>> load_grid();

	std::vector<std::string> load_interp_settings();

	template <typename T>
	void write_values(const Vectors::Vector4D<T>& vec4d, 
		const std::string& data_type);

	std::vector<BkgFPType> cell_centers(std::vector<BkgFPType>& grid);

	template <typename T>
	void read_data_pgkyl(const std::string& species, 
		const std::string& data_type, grid_data_t& grid_data, 
		Vectors::Vector4D<T>& gkyl_data, const Options::Options& opts, 
		const double species_mass_amu=0, const bool force_load=false);

	// Calcuate the electric field using the gradient of gkyl_vp. This must
	// be run after read_potential.
	void calc_elec_field();

	template <typename T>
	void read_elec_density(grid_data_t& grid_data, 
		Vectors::Vector4D<T>& gkyl_ne, const Options::Options& opts);

	template <typename T>
	void read_elec_temperature(grid_data_t& grid_data, 
		Vectors::Vector4D<T>& gkyl_te, const Options::Options& opts);

	template <typename T>
	void read_ion_temperature(grid_data_t& grid_data, 
		Vectors::Vector4D<T>& gkyl_ti, const Options::Options& opts);

	template <typename T>
	void read_potential(grid_data_t& grid_data, 
		Vectors::Vector4D<T>& gkyl_vp, const Options::Options& opts);

	template <typename T>
	void read_magnetic_magnitude(grid_data_t& grid_data, 
		Vectors::Vector4D<T>& gkyl_bmag, const Options::Options& opts);

	template <typename T>
	void read_covariant_comp_b(grid_data_t& grid_data, 
		Vectors::Vector4D<T>& gkyl_b_x, Vectors::Vector4D<T>& gkyl_b_y,
		Vectors::Vector4D<T>& gkyl_b_z, const Options::Options& opts);

	template <typename T>
	void read_magnetic_field(grid_data_t& grid_data, 
		Vectors::Vector4D<T>& gkyl_bX, Vectors::Vector4D<T>& gkyl_bY, 
		Vectors::Vector4D<T>& gkyl_bZ, Vectors::Vector4D<T>& gkyl_bR,
		Vectors::Vector4D<T>& gkyl_bmag, const Options::Options& opts);

	template <typename T>
	void calc_magnetic_field(grid_data_t& grid_data, 
		Vectors::Vector4D<T>& gkyl_bX, Vectors::Vector4D<T>& gkyl_bY, 
		Vectors::Vector4D<T>& gkyl_bZ, Vectors::Vector4D<T>& gkyl_bR,
		Vectors::Vector4D<T>& gkyl_b_x, Vectors::Vector4D<T>& gkyl_b_y, 
		Vectors::Vector4D<T>& gkyl_b_z, Vectors::Vector4D<T>& gkyl_dxdX, 
		Vectors::Vector4D<T>& gkyl_dxdY, Vectors::Vector4D<T>& gkyl_dxdZ, 
		Vectors::Vector4D<T>& gkyl_dydX, Vectors::Vector4D<T>& gkyl_dydY, 
		Vectors::Vector4D<T>& gkyl_dydZ, Vectors::Vector4D<T>& gkyl_dzdX, 
		Vectors::Vector4D<T>& gkyl_dzdY, Vectors::Vector4D<T>& gkyl_dzdZ, 
		const Options::Options& opts);

	template <typename T>
	void read_jacobian(grid_data_t& grid_data, 
		Vectors::Vector4D<T>& gkyl_J, const Options::Options& opts);

	template <typename T>
	void read_gij(grid_data_t& grid_data, 
		Vectors::Vector4D<T>& gkyl_gij, const Options::Options& opts,
		const std::string& idx);

	template <typename T>
	void read_tangent_basis(grid_data_t& grid_data, 
		Vectors::Vector4D<T>& gkyl_dxdz, const Options::Options& opts,
		const std::string& idx);

	template <typename T>
	void read_reciprocal_basis(grid_data_t& grid_data, 
		Vectors::Vector4D<T>& gkyl_dzdx, const Options::Options& opts,
		const std::string& idx);

	std::string get_calc_elec_field_py();

	template <typename T>
	void calc_gradients(grid_data_t& grid_data, 
		const Vectors::Vector4D<BkgFPType>& gkyl_vp, 
		Vectors::Vector4D<BkgFPType>& gkyl_eX, 
		Vectors::Vector4D<BkgFPType>& gkyl_eY, 
		Vectors::Vector4D<BkgFPType>& gkyl_eZ, 
		Vectors::Vector4D<BkgFPType>& gkyl_emag, 
		const Vectors::Vector4D<BkgFPType>& gkyl_bmag, 
		Vectors::Vector4D<BkgFPType>& gkyl_dbdX, 
		Vectors::Vector4D<BkgFPType>& gkyl_dbdY, 
		Vectors::Vector4D<BkgFPType>& gkyl_dbdZ, 
		const Vectors::Vector4D<T>& gkyl_dxdX, 
		const Vectors::Vector4D<T>& gkyl_dxdY, 
		const Vectors::Vector4D<T>& gkyl_dxdZ, 
		const Vectors::Vector4D<T>& gkyl_dydX, 
		const Vectors::Vector4D<T>& gkyl_dydY, 
		const Vectors::Vector4D<T>& gkyl_dydZ, 
		const Vectors::Vector4D<T>& gkyl_dzdX, 
		const Vectors::Vector4D<T>& gkyl_dzdY, 
		const Vectors::Vector4D<T>& gkyl_dzdZ,
		const Options::Options& opts);

	void calc_cell_XYZ_centers(grid_data_t& grid_data, 
		const Options::Options& opts);

	void calc_cell_XYZ_edges(grid_data_t& grid_data, 
		const Options::Options& opts);

	void read_gradient_data(grid_data_t& grid_data, 
		Vectors::Vector4D<BkgFPType>& dataX, 
		Vectors::Vector4D<BkgFPType>& dataY, 
		Vectors::Vector4D<BkgFPType>& dataZ,
		Vectors::Vector4D<BkgFPType>& gkyl_dataX, 
		Vectors::Vector4D<BkgFPType>& gkyl_dataY, 
		Vectors::Vector4D<BkgFPType>& gkyl_dataZ);

	void gradient_ybound_fix(grid_data_t& grid_data, 
		Vectors::Vector4D<BkgFPType>& ecomp);

	void read_elec_field(grid_data_t& grid_data, 
		Vectors::Vector4D<BkgFPType>& gkyl_eX, 
		Vectors::Vector4D<BkgFPType>& gkyl_eY, 
		Vectors::Vector4D<BkgFPType>& gkyl_eZ);

	void read_grad_elec_field(grid_data_t& grid_data, 
		Vectors::Vector4D<BkgFPType>& gkyl_gradeX, 
		Vectors::Vector4D<BkgFPType>& gkyl_gradeY, 
		Vectors::Vector4D<BkgFPType>& gkyl_gradeZ);

	void read_gradb(grid_data_t& grid_data, 
		Vectors::Vector4D<BkgFPType>& gkyl_gradbX, 
		Vectors::Vector4D<BkgFPType>& gkyl_gradbY, 
		Vectors::Vector4D<BkgFPType>& gkyl_gradbZ);

	void write_XYZ(grid_data_t& grid_data, const Options::Options& opts);

	template <typename T>
	void read_par_flow(grid_data_t& grid_data, 
		Vectors::Vector4D<T>& gkyl_uiz, const Options::Options& opts);
	
	template <typename T>
	void read_ion_vperp_sq(grid_data_t& grid_data, 
		Vectors::Vector4D<T>& gkyl_viperp_sq, const Options::Options& opts);

	void calc_bkg_flow(grid_data_t& grid_data, 
		Vectors::Vector4D<BkgFPType>& gkyl_uiz, 
		Vectors::Vector4D<BkgFPType>& gkyl_uX, 
		Vectors::Vector4D<BkgFPType>& gkyl_uY, 
		Vectors::Vector4D<BkgFPType>& gkyl_uZ,
		const Vectors::Vector4D<BkgFPType>& gkyl_bX, 
		const Vectors::Vector4D<BkgFPType>& gkyl_bY,
		const Vectors::Vector4D<BkgFPType>& gkyl_bZ, 
		Vectors::Vector4D<BkgFPType>& gkyl_eX, 
		Vectors::Vector4D<BkgFPType>& gkyl_eY, 
		Vectors::Vector4D<BkgFPType>& gkyl_eZ, 
		Vectors::Vector4D<BkgFPType>& gkyl_gradbX, 
		Vectors::Vector4D<BkgFPType>& gkyl_gradbY,
		Vectors::Vector4D<BkgFPType>& gkyl_gradbZ, 
		Vectors::Vector4D<BkgFPType>& gkyl_viperp_sq, 
		const Options::Options& opts);

	void check_tdim(Vectors::Vector4D<BkgFPType>& gkyl_arr, 
		const Options::Options& opts);
}

#endif
