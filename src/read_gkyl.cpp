#include <iostream>
#include <string>
#include <tuple>
/*
 * Opting out of ADIOS2 since I think it's going to get removed...
	#include <adios2.h>
*/
#include "read_input.h"
#include "read_gkyl.h"
#include "read_gkyl_binary.h"
#include "vectors.h"

namespace Gkyl
{
	// Vectors to hold density, temperature, potential and magnetic field
	// for each frame (assuming electrostatic so only one magnetic field
	// entry is needed. Dimensions are (time, x, y, z). 
	std::vector<Vectors::Vector3D> gkyl_ne {};
	std::vector<Vectors::Vector3D> gkyl_te {};
	std::vector<Vectors::Vector3D> gkyl_vp {};
	std::vector<Vectors::Vector3D> gkyl_b {};

	// Entry point for reading Gkeyll data into Flan. 
	Background read_gkyl()
	{	
		// Resize vectors so we aren't constantly doing this in the loop.
		resize_gkyl_data();

		// Load each needed dataset from Gkeyll
		read_elec_density();

		// With all our arrays assembled, encapsulate them into a Background
		// class object with move semantics then return.
		Background bkg {create_bkg()};
		return bkg;
	}

	// Resize vectors so we aren't constantly doing this in the loop.
	void resize_gkyl_data()
	{
		// Load into local variables so code is easier to read.
		int gkyl_frame_start {Input::get_opt_int(Input::gkyl_frame_start)};
		int gkyl_frame_end {Input::get_opt_int(Input::gkyl_frame_end)};

		gkyl_ne.resize(gkyl_frame_end - gkyl_frame_start + 1);
	}
	
	// Simple helper function to return a string of the full path to a Gkeyll
	// file using it's file name convention.
	std::string assemble_path(const std::string& species, 
		const std::string& ftype, int frame, const std::string& extension) 
	{
		// Load into local variables so code is easier to read.
		std::string gkyl_dir {Input::get_opt_str(Input::gkyl_dir)};
		std::string gkyl_casename {Input::get_opt_str(Input::gkyl_casename)};

		return gkyl_dir + "/" + gkyl_casename + "-" + species + "_" + ftype + 
			"_" + std::to_string(frame) + extension;
	}

	// General function to read in data from Gkeyll into the relevant
	// vector specified by gkyl_data.
	void read_data(const std::string& species, const std::string& ftype,
		std::vector<Vectors::Vector3D>& gkyl_data, int comp)
	{
		// Load into local variables so code is easier to read.
		int gkyl_frame_start {Input::get_opt_int(Input::gkyl_frame_start)};
		int gkyl_frame_end {Input::get_opt_int(Input::gkyl_frame_end)};
		std::string gkyl_file_type {Input::get_opt_str(Input::gkyl_file_type)};
		
		// Loop through the frames.
		for (int f {gkyl_frame_start}; f <= gkyl_frame_end; ++f)
		{
			// Binary (.gkyl) type of file read in.
			if (gkyl_file_type == "binary")
			{
				// Assemble string of path to file.
				std::string path {assemble_path(species, ftype, f, ".gkyl")};
	
				// Read in a binary file.
				double time {};
				Vectors::Vector4D data {};
				std::tie(time, data) = GkylBinary::load_frame(path);
				std::cout << "Time = " << time << '\n';
				
				// The density is the first component in the last dimension.
				gkyl_data[f - gkyl_frame_start] = data.slice_dim4(comp);
			}
		}
	}

	// Read electron density into gkyl_ne.
	void read_elec_density()
	{
		// Load into local variables so code is easier to read.
		std::string gkyl_elec_name {Input::get_opt_str(Input::gkyl_elec_name)};

		// Electron density is component 0 of M0.
		read_data(gkyl_elec_name, "M0", gkyl_ne, 0);
	}

	Background create_bkg()
	{
		Background bkg {};

		// Move each of the vectors we've created into our bkg object
		bkg.get_ne() = std::move(gkyl_ne);
		
		// Okay to return by value since C++11 uses move semantics here.
		return bkg;
	}
}
