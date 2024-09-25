/**
* @file read_gkyl.h
*
* @brief Header file for read_gkyl.cpp
*/

#ifndef READ_GKYL_H
#define READ_GKYL_H

#include <string>
#include <vector>

#include "vectors.h"
#include "background.h"

namespace Gkyl
{

	// Entry point for reading Gkeyll data into Flan. 
	Background::Background read_gkyl();
	
	// Simple helper function to return a string of the full path to a Gkeyll
	// file using it's file name convention.
	std::string assemble_path(const std::string& species, 
		const std::string& ftype, int frame, const std::string& extension);

	// Function to return the full path to read_gkyl.py, which is used to
	// interface with postgkyl and create files that are easily read in by
	// Flan. This path is stored as an environment variable that is set
	// before Flan is run.
	std::string get_read_gkyl_py();

	// Reads in a 1D vector of times corresponding to each frame.
	std::vector<double> load_times();

	// Reads in the x, y, z grid nodes using pgkyl.
	std::tuple<std::vector<double>, std::vector<double>, std::vector<double>> 
		load_grid();

	// Read in data values using pgkyl, returning as a Vector4D.
	template <typename T>
	Vectors::Vector4D<T> load_values(const std::string& data_type);

	// Load file with the interpolation setting used by Gkeyll
	std::vector<std::string> load_interp_settings();

	// Calculate cell centers for the given grid points.
	std::vector<double> cell_centers(std::vector<double>& grid);

	// Function to read in Gkeyll data using a python interface to postgkyl
	// via read_gkyl.py. This produces the following csv files:
	//   bkg_from_pgkyl_times.csv : The time for each frame
	//   bkg_from_pgkyl_grid.csv : Nodes for grid
	//   bkg_from_pgkyl_density.csv : Density arrays for all frames
	//   bkg_from_pgkyl_temperature.csv : Temperature arrays for all frames 
	// The data is loaded and placed into gkyl_data accordingly.
	template <typename T>
	void read_data_pgkyl(const std::string& species, 
		const std::string& data_type,
		Vectors::Vector4D<T>& gkyl_data);

	// Calcuate the electric field using the gradient of gkyl_vp. This must
	// be run after read_potential.
	void calc_elec_field();

	// Read in the corresponding Gkeyll data into vectors.
	void read_elec_density();
	void read_elec_temperature();
	void read_ion_temperature();
	void read_potential();
	void read_magnetic_field();

	// Implementation of gradient as used by numpy.gradient. This is a second
	// order approximation of the derivative.
	double calc_gradient(const double hd, const double hs, const double fd, 
		const double fs, const double f);

	// Calculate the electric field components as the gradient of the potential.
	void calc_elec_field();

	// Move all the data loaded into global arrays into a Background object and
	// return it.
	Background::Background create_bkg();
}

#endif
