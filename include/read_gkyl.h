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
		std::vector<BkgFPType>&, std::vector<BkgFPType>&, std::vector<BkgFPType>&, 
		std::vector<BkgFPType>&, std::vector<BkgFPType>&, std::vector<BkgFPType>&,
		Vectors::Vector3D<BkgFPType>&, Vectors::Vector3D<BkgFPType>&,
		Vectors::Vector3D<BkgFPType>&, Vectors::Vector3D<BkgFPType>&, 
		Vectors::Vector3D<BkgFPType>&, Vectors::Vector3D<BkgFPType>&>;

	/**
	* @brief Entry point for reading Gkeyll data into Flan. 
	*
	* @param opts Options object containing the controlling setting of the
	* simulation.
	*
	* @return Returns a filled in Background object
	*/
	Background::Background read_gkyl(const Options::Options& opts);
	
	// Simple helper function to return a string of the full path to a Gkeyll
	// file using it's file name convention.
	//std::string assemble_path(const std::string& species, 
	//	const std::string& ftype, int frame, const std::string& extension,
	//	const Options::Options& opts);

	/**
	* @brief Get path to python interface script to postgkyl, read_gkyl\.py.
	*
	* Function to return the full path to read_gkyl.py, which is used to
	* interface with postgkyl and create files that are easily read in by
	* Flan. This path is stored as an environment variable that is set
	* as part of the flan conda environment.
	*
	* @return Returns string containing the full path to read_gkyl.py.
	*/
	std::string get_read_gkyl_py();

	/**
	* @brief Reads in a 1D vector of times corresponding to each frame.
	*
	* @return Return a vector of the times.
	*/
	std::vector<BkgFPType> load_times();

	/**
	* @brief Reads in the x, y, z grid nodes using pgkyl.
	*
	* @return Returns a tuple of 3 vectors containing the grid edges - (x,y,z)
	*/
	std::tuple<std::vector<BkgFPType>, std::vector<BkgFPType>, std::vector<BkgFPType>> 
		load_grid();

	/**
	* @brief Read in data values using pgkyl, returning as a Vector4D.
	*
	* @param data_type String identifying which data file to load. This is one
	* of density, temperature, magnetic_field, potential or times.
	*
	* @return Returns a Vector4D object of the data (t,x,y,z).
	*/
	template <typename T>
	Vectors::Vector4D<T> load_values(const std::string& data_type);

	/**
	* @brief Load file with the interpolation setting used by Gkeyll
	*
	* This file is made automatically so it is available to other files which
	* do not have this data included with them (like the magnetic_field data).
	*
	* @returns Returns a vector of two elements: the basis type and the order.
	*/
	std::vector<std::string> load_interp_settings();

	/**
	* @brief Calculate cell centers for the given grid points.
	*
	* Calculate the cell centers for the given grid points. This returns a
	* vector one shorter than the grid point. This only works for uniform
	* grids, just in case that isn't obvious throughout the whole code.
	*
	* @param grid The grid edges for a given dimension.
	* @return Returns a vector of the cell centers that is one less in length
	* than grid.
	*/
	std::vector<BkgFPType> cell_centers(std::vector<BkgFPType>& grid);

	/**
	* @brief Function to read in Gkeyll data using a python interface to postgkyl
	* via read_gkyl.py. 
	*
	* This produces the following csv files:
	*   bkg_from_pgkyl_times.csv : The time for each frame
	*   bkg_from_pgkyl_grid.csv : Nodes for grid
	*   bkg_from_pgkyl_density.csv : Density arrays for all frames
	*   bkg_from_pgkyl_temperature.csv : Temperature arrays for all frames 
	*   bkg_from_pgkyl_potential.csv : Plasma potential arrays for all frames 
	*   bkg_from_pgkyl_magnetic_field.csv : Magnetic field arrays for all frames 
	* The data is loaded and placed into gkyl_data accordingly.
	*
	* @param species Name of the species for the Gkeyll run
	* @param data_type Name of the data to load.\ Can be one of temperature,
	* density, potential or magnetic field. 
	* @param gkyl_data One of the global (to the Gkyl namespace) vectors that
	* corredponds to the data_type being loaded.
	* @param opts Options object that contains all the controlling options
	* of the simulation.
	* @param species_mass_amu The mass of the species in amu.
	*/
	template <typename T>
	void read_data_pgkyl(const std::string& species, 
		const std::string& data_type, grid_data_t& grid_data, 
		Vectors::Vector4D<T>& gkyl_data, const Options::Options& opts, 
		const double species_mass_amu=0, const bool force_load=false);

	// Calcuate the electric field using the gradient of gkyl_vp. This must
	// be run after read_potential.
	void calc_elec_field();

	/**
	* @brief Read electron density into gkyl_ne.
	*/
	template <typename T>
	void read_elec_density(grid_data_t& grid_data, 
		Vectors::Vector4D<T>& gkyl_ne, const Options::Options& opts);

	/**
	* @brief Read electron temperature into gkyl_te.
	*/
	template <typename T>
	void read_elec_temperature(grid_data_t& grid_data, 
		Vectors::Vector4D<T>& gkyl_te, const Options::Options& opts);

	/**
	* @brief Read ion temperature into gkyl_ti.
	*/
	template <typename T>
	void read_ion_temperature(grid_data_t& grid_data, 
		Vectors::Vector4D<T>& gkyl_ti, const Options::Options& opts);

	/**
	* @brief Read plasma potential into gkyl_vp.
	*/
	template <typename T>
	void read_potential(grid_data_t& grid_data, 
		Vectors::Vector4D<T>& gkyl_vp, const Options::Options& opts);

	/**
	* @brief Read magnetic field into gkyl_bX, gkyl_bY and gkyl_bZ.
	*/
	template <typename T>
	void read_magnetic_field(grid_data_t& grid_data, 
		Vectors::Vector4D<T>& gkyl_bX, Vectors::Vector4D<T>& gkyl_bY, 
		Vectors::Vector4D<T>& gkyl_bZ, const Options::Options& opts);

	/**
	* @brief Read Jacobian into gkyl_J.
	*/
	template <typename T>
	void read_jacobian(grid_data_t& grid_data, 
		Vectors::Vector4D<T>& gkyl_J, const Options::Options& opts);

	/*
	* @brief Read metric coefficient gij component into gkyl_gijXX.
	*/
	template <typename T>
	void read_gij(grid_data_t& grid_data, 
		Vectors::Vector4D<T>& gkyl_gij, const Options::Options& opts,
		const std::string& idx);

	/**
	* @brief Get path to python script for calculating the electric field, 
	* calc_elec_field.py.
	*
	* Function to return the full path to calc_elec_field.py, which uses
	* the files cell_center_XYZ.csv and bkg_from_pgkyl_potential.csv to
	* calculate the electric field components from the potential gradient.
	* This path is stored as an environment variable that is set
	* as part of the flan conda environment.
	*
	* @return Returns string containing the full path to calc_elec_field.py.
	*/
	std::string get_calc_elec_field_py();

	/**
	* @brief Calculate the electric field components as the gradient of the 
	* potential. This function actually calls the python script calc_elec.py.
	*/
	void calc_elec_field();

	/**
	* @brief Calculate cell X,Y,Z coordinates for center of cells
	*/
	void calc_cell_XYZ_centers(grid_data_t& grid_data, 
		const Options::Options& opts);

	/**
	* @brief Calculate cell X,Y,Z coordinates for grid edges
	*/
	void calc_cell_XYZ_edges(grid_data_t& grid_data, 
		const Options::Options& opts);
	/**
	* @brief Read in electric field components (X,Y,Z) into arrays
	*/
	void read_elec_field(grid_data_t& grid_data, 
		Vectors::Vector4D<BkgFPType>& gkyl_eX, 
		Vectors::Vector4D<BkgFPType>& gkyl_eY, 
		Vectors::Vector4D<BkgFPType>& gkyl_eZ);

	/**
	* @brief Apply a fix to the y boundaries of the electric field where the
	* gradient generally messes up due to not as many points around.
	*
	* This is generally a hack, since it is not clear what to do in this 
	* situation. Ideally, the gradient calculation would handle this, but I
	* havent been able to figure it out.
	*/
	void elec_field_ybound_fix(grid_data_t& grid_data, 
		Vectors::Vector4D<BkgFPType>& ecomp);

	/**
	* @brief Write out Cartesian (X,Y,Z) coordinates of cell centers
	*
	* This function writes to a file of X,Y,Z coordinates,
	* so that we can read them into a python script and take advantage of
	* the scipy library to easily calculate a gradient on an irregular grid.
	*/
	void write_XYZ(grid_data_t& grid_data, const Options::Options& opts);

	/**
	* @brief Read in mean flow components (X,Y,Z) into arrays.
	*/
	template <typename T>
	void read_mean_flow(grid_data_t& grid_data, 
		Vectors::Vector4D<T>& gkyl_uX, Vectors::Vector4D<T>& gkyl_uY, 
		Vectors::Vector4D<T>& gkyl_uZ, const Options::Options& opts);
}

#endif
