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
		std::vector<Vectors::Vector3D> m_ne {};
	public:
		std::vector<Vectors::Vector3D>& get_ne() {return m_ne;}

		void resize_vectors(int nframes, int dim1, int dim2, int dim3)
		{
			
		}
	};

	// Entry point for reading Gkeyll data into Flan. 
	Background read_gkyl();
	
	// Resize vectors so we aren't constantly doing this in the loop.
	void resize_gkyl_data();

	// Simple helper function to return a string of the full path to a Gkeyll
	// file using it's file name convention.
	std::string assemble_path(const std::string& species, 
		const std::string& ftype, int frame, const std::string& extension);

	// General function to read in data from Gkeyll into the relevant
	// vector specified by gkyl_data.
	void read_data(const std::string& species, const std::string& ftype,
		std::vector<Vectors::Vector3D>& gkyl_data, int comp);

	// Read electron density into gkyl_ne.
	void read_elec_density();

	Background create_bkg();
}

#endif
