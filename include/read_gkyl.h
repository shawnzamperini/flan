#ifndef READ_GKYL_H
#define READ_GKYL_H

#include <string>
#include <vector>

#include "vectors.h"

namespace Gkyl
{
	// Encapsulating class containing the Gkeyll background.
	class Background
	{
	private:
		std::vector<double> m_times {};
		Vectors::Vector4D m_ne {};

	public:

		// Accessors
		Vectors::Vector4D& get_ne() {return m_ne;}

		// Setters
		void move_into_ne(Vectors::Vector4D& ne) 
		{
			m_ne.move_into_data(ne);
		}

		// TBD
		void resize_vectors(int nframes, int dim1, int dim2, int dim3)
		{
			
		}
	};

	// Entry point for reading Gkeyll data into Flan. 
	Background read_gkyl();
	
	// Resize vectors so we aren't constantly doing this in the loop.
	//void resize_gkyl_data();

	// Simple helper function to return a string of the full path to a Gkeyll
	// file using it's file name convention.
	std::string assemble_path(const std::string& species, 
		const std::string& ftype, int frame, const std::string& extension);

	// Function to return the full path to read_gkyl.py, which is used to
	// interface with postgkyl and create files that are easily read in by
	// Flan. This path is stored as an environment variable that is set
	// before Flan is run.
	std::string get_read_gkyl_py();

	// Function to read in Gkeyll data using a python interface to postgkyl
	// via read_gkyl.py. This produces the following csv files:
	//   bkg_from_pgkyl_times.csv : The time for each frame
	//   bkg_from_pgkyl_grid.csv : Nodes for grid
	//   bkg_from_pgkyl_density.csv : Density arrays for all frames
	//   bkg_from_pgkyl_temperature.csv : Temperature arrays for all frames 
	// The data is loaded and placed into gkyl_data accordingly.
	void read_data_pgkyl(const std::string& species, 
		const std::string& data_type,
		std::vector<Vectors::Vector3D>& gkyl_data);

	// General function to read in data from Gkeyll into the relevant
	// vector specified by gkyl_data.
	void read_data(const std::string& species, const std::string& ftype,
		std::vector<Vectors::Vector3D>& gkyl_data, int comp);

	// Read in the corresponding Gkeyll data into vectors.
	void read_elec_density();
	void read_elec_temperature();
	void read_ion_temperature();
	void read_potential();
	void read_magnetic_field();

	Background create_bkg();
}

#endif
