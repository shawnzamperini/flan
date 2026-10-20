/**
* @file  read_gkyl.cpp
*
* @brief Routines handling reading in a plasma background from Gkeyll.
*
* The Gkeyll data is read through read_gkyl.py, a python interface to
* postgkyl, which saves each dataset to a binary file named
* bkg_from_pgkyl_*.bin along with a few .csv files holding the times, grid
* and interpolation settings. These files are then read in here. If a file
* already exists from a previous run it is loaded directly instead of calling
* read_gkyl.py again, so delete the files to regenerate the background.
*/

#include <cmath>
#include <filesystem>
#include <fstream>
#include <functional>  // For std::multiplies
#include <iostream>
#include <numeric>     // For std::accumulate
#include <omp.h>
#include <string>
#include <sstream>
#include <type_traits> // For std::is_same
#include <tuple>

#include "background.h"
#include "constants.h"
#include "errors.h"
#include "flan_types.h"
#include "options.h"
#include "read_gkyl.h"
#include "vectors.h"

namespace Gkyl
{
	/**
	* Entry point for reading Gkeyll data into Flan.
	*
	* Each needed dataset is loaded from the Gkeyll case, either by calling
	* read_gkyl.py through postgkyl or, if they already exist, from the files
	* it saved on a previous run. The electric field and magnetic field
	* gradient are then calculated, along with the background flow if the
	* moments are BiMaxwellian. Finally every array is moved (not copied)
	* into a Background object.
	*
	* Time-dependent data is 4D with dimensions (t, x, y, z). Time-
	* independent data (the Jacobian, metric coefficients, basis vectors and
	* covariant components of the magnetic field) is read in as 4D, and a
	* single time slice is stored as 3D (x, y, z). In both cases the arrays
	* are indexed according to the Gkeyll computational coordinate indices.
	*
	* This file uses the convention that all variable names are lowercase,
	* except x,y,z always refers to computational (Gkeyll, curvilinear)
	* coordinates and X,Y,Z always refers to physical (Cartesian)
	* coordinates.
	*
	* @param opts Simulation options.
	* @return Background object holding the Gkeyll data.
	*/
	Background::Background read_gkyl(const Options::Options& opts)
	{
		using namespace std::string_literals;

		// Vectors to hold all the various arrays loaded from a Gkeyll case.
		// I feel like the logic could've shaken out better here, but in the
		// end the data in these arrays are moved into a Background object
		// (not copied) and then go out of scope, so end of the day it's fine.

		// Time of each frame, and the computational grid. x, y, z are the
		// cell centers and grid_x, grid_y, grid_z are the grid edges.
		std::vector<BkgFPType> gkyl_times {};
		std::vector<BkgFPType> gkyl_x {};
		std::vector<BkgFPType> gkyl_y {};
		std::vector<BkgFPType> gkyl_z {};
		std::vector<BkgFPType> gkyl_grid_x {};
		std::vector<BkgFPType> gkyl_grid_y {};
		std::vector<BkgFPType> gkyl_grid_z {};

		// Cell centers and grid edges in physical space
		Vectors::Vector3D<BkgFPType> gkyl_X {};
		Vectors::Vector3D<BkgFPType> gkyl_Y {};
		Vectors::Vector3D<BkgFPType> gkyl_Z {};
		Vectors::Vector3D<BkgFPType> gkyl_grid_X {};
		Vectors::Vector3D<BkgFPType> gkyl_grid_Y {};
		Vectors::Vector3D<BkgFPType> gkyl_grid_Z {};

		// Metric coefficients and Jacobian
		Vectors::Vector4D<BkgFPType> gkyl_gij_00 {};
		Vectors::Vector4D<BkgFPType> gkyl_gij_01 {};
		Vectors::Vector4D<BkgFPType> gkyl_gij_02 {};
		Vectors::Vector4D<BkgFPType> gkyl_gij_11 {};
		Vectors::Vector4D<BkgFPType> gkyl_gij_12 {};
		Vectors::Vector4D<BkgFPType> gkyl_gij_22 {};
		Vectors::Vector4D<BkgFPType> gkyl_J {};

		// Plasma density, temperatures, ion perpendicular thermal velocity
		// squared and plasma potential
		Vectors::Vector4D<BkgFPType> gkyl_ne {};
		Vectors::Vector4D<BkgFPType> gkyl_te {};
		Vectors::Vector4D<BkgFPType> gkyl_ti {};
		Vectors::Vector4D<BkgFPType> gkyl_viperp_sq {};
		Vectors::Vector4D<BkgFPType> gkyl_vp {};

		// Magnetic field magnitude and covariant components
		Vectors::Vector4D<BkgFPType> gkyl_bmag {};
		Vectors::Vector4D<BkgFPType> gkyl_b_x {};
		Vectors::Vector4D<BkgFPType> gkyl_b_y {};
		Vectors::Vector4D<BkgFPType> gkyl_b_z {};

		// Magnetic field components and its gradient in physical space
		Vectors::Vector4D<BkgFPType> gkyl_bX {};
		Vectors::Vector4D<BkgFPType> gkyl_bY {};
		Vectors::Vector4D<BkgFPType> gkyl_bZ {};
		Vectors::Vector4D<BkgFPType> gkyl_bR {};
		Vectors::Vector4D<BkgFPType> gkyl_dbdX {};
		Vectors::Vector4D<BkgFPType> gkyl_dbdY {};
		Vectors::Vector4D<BkgFPType> gkyl_dbdZ {};

		// Electric field components, magnitude and gradient components in
		// physical space
		Vectors::Vector4D<BkgFPType> gkyl_eX {};
		Vectors::Vector4D<BkgFPType> gkyl_eY {};
		Vectors::Vector4D<BkgFPType> gkyl_eZ {};
		Vectors::Vector4D<BkgFPType> gkyl_emag {};
		Vectors::Vector4D<BkgFPType> gkyl_gradeX {};
		Vectors::Vector4D<BkgFPType> gkyl_gradeY {};
		Vectors::Vector4D<BkgFPType> gkyl_gradeZ {};

		// Parallel (to z) ion flow, and the flow components in physical
		// space
		Vectors::Vector4D<BkgFPType> gkyl_uz {};
		Vectors::Vector4D<BkgFPType> gkyl_uX {};
		Vectors::Vector4D<BkgFPType> gkyl_uY {};
		Vectors::Vector4D<BkgFPType> gkyl_uZ {};

		// Tangent (covariant) basis vectors
		Vectors::Vector4D<BkgFPType> gkyl_dXdx {};
		Vectors::Vector4D<BkgFPType> gkyl_dYdx {};
		Vectors::Vector4D<BkgFPType> gkyl_dZdx {};
		Vectors::Vector4D<BkgFPType> gkyl_dXdy {};
		Vectors::Vector4D<BkgFPType> gkyl_dYdy {};
		Vectors::Vector4D<BkgFPType> gkyl_dZdy {};
		Vectors::Vector4D<BkgFPType> gkyl_dXdz {};
		Vectors::Vector4D<BkgFPType> gkyl_dYdz {};
		Vectors::Vector4D<BkgFPType> gkyl_dZdz {};

		// Reciprocal (contravariant) basis vectors
		Vectors::Vector4D<BkgFPType> gkyl_dxdX {};
		Vectors::Vector4D<BkgFPType> gkyl_dxdY {};
		Vectors::Vector4D<BkgFPType> gkyl_dxdZ {};
		Vectors::Vector4D<BkgFPType> gkyl_dydX {};
		Vectors::Vector4D<BkgFPType> gkyl_dydY {};
		Vectors::Vector4D<BkgFPType> gkyl_dydZ {};
		Vectors::Vector4D<BkgFPType> gkyl_dzdX {};
		Vectors::Vector4D<BkgFPType> gkyl_dzdY {};
		Vectors::Vector4D<BkgFPType> gkyl_dzdZ {};

		// These are passed to every function call, so save some time/space by
		// zipping them up into a tuple and passing together.
		//   grid_data[0]  = gkyl_times
		//   grid_data[1]  = gkyl_x
		//   grid_data[2]  = gkyl_y
		//   grid_data[3]  = gkyl_z
		//   grid_data[4]  = gkyl_grid_x
		//   grid_data[5]  = gkyl_grid_y
		//   grid_data[6]  = gkyl_grid_z
		//   grid_data[7]  = gkyl_X
		//   grid_data[8]  = gkyl_Y
		//   grid_data[9]  = gkyl_Z
		//   grid_data[10] = gkyl_grid_X
		//   grid_data[11] = gkyl_grid_Y
		//   grid_data[12] = gkyl_grid_Z
		grid_data_t grid_data {gkyl_times, gkyl_x, gkyl_y, gkyl_z,
			gkyl_grid_x, gkyl_grid_y, gkyl_grid_z, gkyl_X, gkyl_Y, gkyl_Z,
			gkyl_grid_X, gkyl_grid_Y, gkyl_grid_Z};

		// --------------------
		// Plasma
		// --------------------

		std::cout << "Loading Gkeyll data...\n";
		std::cout << "  - Electron density\n";
		read_elec_density(grid_data, gkyl_ne, opts);

		// After the first array, confirm that the number of times read
		// matches that from the input file. This can happen if the user
		// changes fstart or fend without deleting the existing background
		// files first. It is safe to assume the x,y,z dimensions don't
		// change, since that would mean a completely different simulation
		// (which would mean the user has no idea what they're doing, and we
		// should fail!).
		check_tdim(gkyl_ne, opts);

		std::cout << "  - Electron temperature\n";
		read_elec_temperature(grid_data, gkyl_te, opts);

		std::cout << "  - Ion temperature\n";
		read_ion_temperature(grid_data, gkyl_ti, opts);

		std::cout << "  - Plasma potential\n";
		read_potential(grid_data, gkyl_vp, opts);

		//std::cout << "  - Magnetic field magnitude\n";
		//read_magnetic_magnitude(grid_data, gkyl_bmag, opts);

		// Calculate the Cartesian (X,Y,Z) coordinate of each cell center.
		// These aren't actually used, but left here just in case they become
		// useful in the future.
		//std::cout << "  - X,Y,Z coordinates\n";
		//calc_cell_XYZ_centers(grid_data, opts);
		//calc_cell_XYZ_edges(grid_data, opts);

		// Write out the (X,Y,Z) coordinates so that we can load them in
		// python to take advantage of scipy library in calculating gradients
		// on irregular grids (i.e., for the electric field).
		//write_XYZ(grid_data, opts);

		// --------------------
		// Geometry
		// --------------------

		std::cout << "  - Jacobian\n";
		read_jacobian(grid_data, gkyl_J, opts);

		// Metric coefficients. Each has its own vector. Not actually used in
		// Flan but it can be useful in post-analysis.
		std::cout << "  - Metric coefficients\n";
		read_gij(grid_data, gkyl_gij_00, opts, "00");
		read_gij(grid_data, gkyl_gij_01, opts, "01");
		read_gij(grid_data, gkyl_gij_02, opts, "02");
		read_gij(grid_data, gkyl_gij_11, opts, "11");
		read_gij(grid_data, gkyl_gij_12, opts, "12");
		read_gij(grid_data, gkyl_gij_22, opts, "22");

		// Tangent (covariant) basis vectors (dX/dx). We don't actually need
		// these, but they're included here for completeness just in case
		// they become needed.
		std::cout << "  - Tangent basis vectors\n";
		read_tangent_basis(grid_data, gkyl_dXdx, opts, "Xx");
		read_tangent_basis(grid_data, gkyl_dYdx, opts, "Yx");
		read_tangent_basis(grid_data, gkyl_dZdx, opts, "Zx");
		read_tangent_basis(grid_data, gkyl_dXdy, opts, "Xy");
		read_tangent_basis(grid_data, gkyl_dYdy, opts, "Yy");
		read_tangent_basis(grid_data, gkyl_dZdy, opts, "Zy");
		read_tangent_basis(grid_data, gkyl_dXdz, opts, "Xz");
		read_tangent_basis(grid_data, gkyl_dYdz, opts, "Yz");
		read_tangent_basis(grid_data, gkyl_dZdz, opts, "Zz");

		// Reciprocal (contravariant) basis vectors (dx/dX)
		std::cout << "  - Reciprocal basis vectors\n";
		read_reciprocal_basis(grid_data, gkyl_dxdX, opts, "xX");
		read_reciprocal_basis(grid_data, gkyl_dxdY, opts, "xY");
		read_reciprocal_basis(grid_data, gkyl_dxdZ, opts, "xZ");
		read_reciprocal_basis(grid_data, gkyl_dydX, opts, "yX");
		read_reciprocal_basis(grid_data, gkyl_dydY, opts, "yY");
		read_reciprocal_basis(grid_data, gkyl_dydZ, opts, "yZ");
		read_reciprocal_basis(grid_data, gkyl_dzdX, opts, "zX");
		read_reciprocal_basis(grid_data, gkyl_dzdY, opts, "zY");
		read_reciprocal_basis(grid_data, gkyl_dzdZ, opts, "zZ");

		// Covariant components of the magnetic field
		std::cout << "  - Covariant components\n";
		read_covariant_comp_b(grid_data, gkyl_b_x, gkyl_b_y, gkyl_b_z, opts);

		// --------------------
		// Fields
		// --------------------

		// Calculate Cartesian components of the magnetic field. The function
		// calc_magnetic_field uses the covariant components of B and the
		// reciprocal basis vectors to calculate the Cartesian components of
		// B, but it is problematic in that it is not clear if Bmag is baked
		// into the covariant components or not. So instead we use
		// read_magnetic_field, where it loads the Cartesian unit vector,
		// rescales it to 1.0 magnitude and then multiplies each component
		// by Bmag. It seems to work fine, but keep an eye on it.
		std::cout << "  - Magnetic field\n";
		//calc_magnetic_field(grid_data, gkyl_bX, gkyl_bY, gkyl_bZ, gkyl_bR,
		//	gkyl_b_x, gkyl_b_y, gkyl_b_z, gkyl_dxdX, gkyl_dxdY, gkyl_dxdZ,
		//	gkyl_dydX, gkyl_dydY, gkyl_dydZ, gkyl_dzdX, gkyl_dzdY, gkyl_dzdZ,
		//	opts);
		read_magnetic_field(grid_data, gkyl_bX, gkyl_bY, gkyl_bZ, gkyl_bR,
			gkyl_bmag, opts);

		// At this point we can calculate gradient fields (E, gradB)
		std::cout << "  - Calculating gradient fields\n";
		calc_gradients(grid_data, gkyl_vp, gkyl_eX, gkyl_eY, gkyl_eZ,
			gkyl_emag, gkyl_bmag, gkyl_dbdX, gkyl_dbdY, gkyl_dbdZ,
			gkyl_dxdX, gkyl_dxdY, gkyl_dxdZ,
			gkyl_dydX, gkyl_dydY, gkyl_dydZ,
			gkyl_dzdX, gkyl_dzdY, gkyl_dzdZ, opts);

		// Electric field gradient (optional)
		if (opts.calc_grad_elec_int() == 1)
		{
			std::cout << "  - Warning: Electric field gradient requested but "
				<< "not implemented! It's easy though, ask Shawn to do it.\n";
		}

		// --------------------
		// Background flow
		// --------------------

		// These arrays are only needed if the gkyl moments are BiMaxwellian
		// (that is, it has vperp and vpar). We use vperp and gradB in the
		// background flow calculation, which are used in the collision
		// model. If you don't need collisions, then the other moment file
		// types are fine and this section is safely skipped.
		if (opts.gkyl_moment_type() == "bimaxwellian")
		{
			std::cout << "  - Ion parallel flow\n";
			read_par_flow(grid_data, gkyl_uz, opts);

			std::cout << "  - Ion perpendicular thermal velocity\n";
			read_ion_vperp_sq(grid_data, gkyl_viperp_sq, opts);

			// Calculate the background flows in each Cartesian direction
			std::cout << "  - Calculating background flow\n";
			calc_bkg_flow(grid_data, gkyl_uz, gkyl_uX, gkyl_uY, gkyl_uZ,
				gkyl_bX, gkyl_bY, gkyl_bZ, gkyl_eX, gkyl_eY, gkyl_eZ,
				gkyl_dbdX, gkyl_dbdY, gkyl_dbdZ, gkyl_viperp_sq, opts);
		}

		// --------------------
		// Assemble Background
		// --------------------

		// With all our arrays assembled, encapsulate them into a Background
		// class object with move semantics then return.
		Background::Background bkg {};

		// Move each of the vectors we've created into our bkg object
		bkg.move_into_times(gkyl_times);
		bkg.move_into_x(gkyl_x);
		bkg.move_into_y(gkyl_y);
		bkg.move_into_z(gkyl_z);
		bkg.move_into_grid_x(gkyl_grid_x);
		bkg.move_into_grid_y(gkyl_grid_y);
		bkg.move_into_grid_z(gkyl_grid_z);
		bkg.move_into_ne(gkyl_ne);
		bkg.move_into_te(gkyl_te);
		bkg.move_into_ti(gkyl_ti);
		bkg.move_into_vp(gkyl_vp);
		bkg.move_into_bX(gkyl_bX);
		bkg.move_into_bY(gkyl_bY);
		bkg.move_into_bZ(gkyl_bZ);
		bkg.move_into_bR(gkyl_bR);  // Just for debugging, not actually used
		bkg.move_into_bmag(gkyl_bmag);
		bkg.move_into_eX(gkyl_eX);
		bkg.move_into_eY(gkyl_eY);
		bkg.move_into_eZ(gkyl_eZ);
		bkg.move_into_emag(gkyl_emag);
		bkg.move_into_uX(gkyl_uX);
		bkg.move_into_uY(gkyl_uY);
		bkg.move_into_uZ(gkyl_uZ);
		bkg.move_into_X(gkyl_X);
		bkg.move_into_Y(gkyl_Y);
		bkg.move_into_Z(gkyl_Z);
		bkg.move_into_grid_X(gkyl_grid_X);
		bkg.move_into_grid_Y(gkyl_grid_Y);
		bkg.move_into_grid_Z(gkyl_grid_Z);

		// Only move these if they were calculated
		if (opts.calc_grad_elec_int() == 1)
		{
			bkg.move_into_gradeX(gkyl_gradeX);
			bkg.move_into_gradeY(gkyl_gradeY);
			bkg.move_into_gradeZ(gkyl_gradeZ);
		}

		// Only needed for debugging purposes, delete if not debugging
		bkg.move_into_gradbX(gkyl_dbdX);
		bkg.move_into_gradbY(gkyl_dbdY);
		bkg.move_into_gradbZ(gkyl_dbdZ);

		// Special case. gkyl_J is a 4D vector so we didn't have to write a
		// whole new function just for a 3D vector - less maintainable code.
		// So we just need to index any random timeslice as a Vector3D, pass
		// it in.
		Vectors::Vector3D gkyl_J_3D {gkyl_J.slice_dim1(0)};
		bkg.move_into_J(gkyl_J_3D);

		// Likewise for the 6 metric coefficient arrays
		Vectors::Vector3D gkyl_gij_00_3D {gkyl_gij_00.slice_dim1(0)};
		bkg.move_into_gij_00(gkyl_gij_00_3D);
		Vectors::Vector3D gkyl_gij_01_3D {gkyl_gij_01.slice_dim1(0)};
		bkg.move_into_gij_01(gkyl_gij_01_3D);
		Vectors::Vector3D gkyl_gij_02_3D {gkyl_gij_02.slice_dim1(0)};
		bkg.move_into_gij_02(gkyl_gij_02_3D);
		Vectors::Vector3D gkyl_gij_11_3D {gkyl_gij_11.slice_dim1(0)};
		bkg.move_into_gij_11(gkyl_gij_11_3D);
		Vectors::Vector3D gkyl_gij_12_3D {gkyl_gij_12.slice_dim1(0)};
		bkg.move_into_gij_12(gkyl_gij_12_3D);
		Vectors::Vector3D gkyl_gij_22_3D {gkyl_gij_22.slice_dim1(0)};
		bkg.move_into_gij_22(gkyl_gij_22_3D);

		// Likewise for the reciprocal basis vectors
		Vectors::Vector3D gkyl_dxdX_3D {gkyl_dxdX.slice_dim1(0)};
		bkg.move_into_dxdX(gkyl_dxdX_3D);
		Vectors::Vector3D gkyl_dxdY_3D {gkyl_dxdY.slice_dim1(0)};
		bkg.move_into_dxdY(gkyl_dxdY_3D);
		Vectors::Vector3D gkyl_dxdZ_3D {gkyl_dxdZ.slice_dim1(0)};
		bkg.move_into_dxdZ(gkyl_dxdZ_3D);
		Vectors::Vector3D gkyl_dydX_3D {gkyl_dydX.slice_dim1(0)};
		bkg.move_into_dydX(gkyl_dydX_3D);
		Vectors::Vector3D gkyl_dydY_3D {gkyl_dydY.slice_dim1(0)};
		bkg.move_into_dydY(gkyl_dydY_3D);
		Vectors::Vector3D gkyl_dydZ_3D {gkyl_dydZ.slice_dim1(0)};
		bkg.move_into_dydZ(gkyl_dydZ_3D);
		Vectors::Vector3D gkyl_dzdX_3D {gkyl_dzdX.slice_dim1(0)};
		bkg.move_into_dzdX(gkyl_dzdX_3D);
		Vectors::Vector3D gkyl_dzdY_3D {gkyl_dzdY.slice_dim1(0)};
		bkg.move_into_dzdY(gkyl_dzdY_3D);
		Vectors::Vector3D gkyl_dzdZ_3D {gkyl_dzdZ.slice_dim1(0)};
		bkg.move_into_dzdZ(gkyl_dzdZ_3D);

		// Likewise for the tangent basis vectors
		Vectors::Vector3D gkyl_dXdx_3D {gkyl_dXdx.slice_dim1(0)};
		bkg.move_into_dXdx(gkyl_dXdx_3D);
		Vectors::Vector3D gkyl_dYdx_3D {gkyl_dYdx.slice_dim1(0)};
		bkg.move_into_dYdx(gkyl_dYdx_3D);
		Vectors::Vector3D gkyl_dZdx_3D {gkyl_dZdx.slice_dim1(0)};
		bkg.move_into_dZdx(gkyl_dZdx_3D);
		Vectors::Vector3D gkyl_dXdy_3D {gkyl_dXdy.slice_dim1(0)};
		bkg.move_into_dXdy(gkyl_dXdy_3D);
		Vectors::Vector3D gkyl_dYdy_3D {gkyl_dYdy.slice_dim1(0)};
		bkg.move_into_dYdy(gkyl_dYdy_3D);
		Vectors::Vector3D gkyl_dZdy_3D {gkyl_dZdy.slice_dim1(0)};
		bkg.move_into_dZdy(gkyl_dZdy_3D);
		Vectors::Vector3D gkyl_dXdz_3D {gkyl_dXdz.slice_dim1(0)};
		bkg.move_into_dXdz(gkyl_dXdz_3D);
		Vectors::Vector3D gkyl_dYdz_3D {gkyl_dYdz.slice_dim1(0)};
		bkg.move_into_dYdz(gkyl_dYdz_3D);
		Vectors::Vector3D gkyl_dZdz_3D {gkyl_dZdz.slice_dim1(0)};
		bkg.move_into_dZdz(gkyl_dZdz_3D);

		// Likewise for the covariant components of magnetic field
		Vectors::Vector3D gkyl_b_x_3D {gkyl_b_x.slice_dim1(0)};
		bkg.move_into_b_x(gkyl_b_x_3D);
		Vectors::Vector3D gkyl_b_y_3D {gkyl_b_y.slice_dim1(0)};
		bkg.move_into_b_y(gkyl_b_y_3D);
		Vectors::Vector3D gkyl_b_z_3D {gkyl_b_z.slice_dim1(0)};
		bkg.move_into_b_z(gkyl_b_z_3D);

		// Okay to return by value since C++11 uses move semantics here.
		return bkg;
	}

	/**
	* Returns the full path to read_gkyl.py, which is used to interface with
	* postgkyl and create files that are easily read in by Flan.
	*
	* The path is stored in the READ_GKYL_PY environment variable, which is
	* set before Flan is run. If it is not set, the run is stopped with
	* Errors::fatal_error, since the background can't be read without it.
	*
	* @return Full path to read_gkyl.py.
	*/
	std::string get_read_gkyl_py()
	{
		if (const char* read_gkyl_py {std::getenv("READ_GKYL_PY")})
		{
			return static_cast<std::string>(read_gkyl_py);
		}
		else
		{
			Errors::fatal_error("READ_GKYL_PY is not set. Please set this "
				"to the full path to the read_gkyl.py script. E.g.:\n"
				"  READ_GKYL_PY=/path/to/flan/python/read_gkyl.py\n"
				"This error probably means you forgot to activate the flan "
				"conda environment.");
		}
	}

	/**
	* Reads in the time of each frame from bkg_from_pgkyl_times.csv.
	*
	* Lines starting with # are comments. The first value is the number of
	* frames, followed by one time per line.
	*
	* @return Time of each frame (s).
	*/
	std::vector<BkgFPType> load_times()
	{
		// Filename is treated as a constant.
		std::string filename {"bkg_from_pgkyl_times.csv"};
		std::ifstream fstream {filename};

		// Go through one line at a time
		std::string line {};
		std::vector<BkgFPType> times {};
		int ntimes {};
		bool ntimes_read {false};
		while(std::getline(fstream, line))
		{
			// If line starts with #, then it's a comment and ignore
			if (!line.starts_with("#"))
			{
				// The first line is a single integer telling us how many
				// times (frames) there are.
				if (!ntimes_read)
				{
					ntimes = std::stoi(line);
					ntimes_read = true;
					//std::cout << "Reading " << ntimes << " times\n";

					// Set the capacity of the vector now that we know what
					// it will be.
					times.reserve(ntimes);
				}
				else
				{
					// Append time to our vector. Note we already set
					// the capacity of the vector, so this should not be too
					// expensive.
					times.push_back(std::stod(line));
				}
			}
		}

		// Close and return (uses move semantics)
		fstream.close();
		return times;
	}

	/**
	* Reads in the x, y, z grid edges from bkg_from_pgkyl_grid.csv.
	*
	* Lines starting with # are comments. The first line holds the number
	* of grid edges in x, y and z. Then all the x values are listed, then the
	* y values, then the z values, one per line. The run is stopped with
	* Errors::fatal_error if there are not 3 dimensions or there are too
	* many values.
	*
	* @return Tuple of the grid edges in x, y and z.
	*/
	std::tuple<std::vector<BkgFPType>, std::vector<BkgFPType>,
		std::vector<BkgFPType>> load_grid()
	{
		// Filename is treated as a constant.
		std::string filename {"bkg_from_pgkyl_grid.csv"};
		std::ifstream fstream {filename};

		// Define variables needed for the following loop
		std::string line {};
		std::vector<BkgFPType> grid_x {};
		std::vector<BkgFPType> grid_y {};
		std::vector<BkgFPType> grid_z {};
		std::vector<int> dims {};
		[[maybe_unused]] int ngrid {};
		bool ngrid_read {false};
		int count {};

		// Go through one line at a time
		while(std::getline(fstream, line))
		{
			// If line starts with #, then it's a comment and ignore
			if (!line.starts_with("#"))
			{
				// The first line is three integers telling us the x, y, z
				// dimensions. These are the edges of each grid cell, so
				// there is one extra value in each dimension when compared
				// to the actual value (e.g., density) dimensions.
				if (!ngrid_read)
				{
					// Use string stream to pull the three integers out.
					std::stringstream ss {line};
					int tmp_int {};
					while (ss >> tmp_int)
					{
						dims.push_back(tmp_int);
					}

					// If we have anything other than 3 values something
					// went wrong.
					if (dims.size() != 3)
					{
						std::ostringstream msg {};
						msg << "Expected 3 grid dimensions in " << filename 
							<< " but read " << dims.size() << ". Check that "
							<< "this file is being generated correctly.\n"
							<< "  dims:";
						for (auto i : dims) msg << ' ' << i;
						Errors::fatal_error(msg.str());
					}

					// ngrid is the total number of points in the file. The
					// x values are all printed, then the y then the z.
					ngrid = dims[0] + dims[1] + dims[2];
					ngrid_read = true;
					//std::cout << "Reading " << ngrid << " grid values\n";

					// Set the capacity of the vectors now that we know them.
					grid_x.reserve(dims[0]);
					grid_y.reserve(dims[1]);
					grid_z.reserve(dims[2]);
				}
				else
				{
					// Append to each vector accordingly
					BkgFPType grid_value {static_cast<BkgFPType>(
						std::stod(line))};
					if (count < dims[0])
						grid_x.push_back(grid_value);
					else if (count >= dims[0] && count < dims[0] + dims[1])
						grid_y.push_back(grid_value);
					else if (count >= dims[0] + dims[1])
						grid_z.push_back(grid_value);

					// If count ends up being greater than all the values
					// we need to read in, something went wrong.
					if (count >= dims[0] + dims[1] + dims[2])
					{
						std::ostringstream msg {};
						msg << "There appear to be too many values in " 
							<< filename << ". Check that this file is being "
							<< "generated correctly.\n"
							<< "  count = " << count << ", dims = " << dims[0]
							<< ", " << dims[1] << ", " << dims[2];
						Errors::fatal_error(msg.str());
					}

					++count;
				}
			}
		}
		fstream.close();

		// Pack up into tuple and return
		return std::make_tuple(grid_x, grid_y, grid_z);
	}


	/**
	* Reads in the values of a dataset from bkg_from_pgkyl_[data_type].bin,
	* which was saved by read_gkyl.py.
	*
	* The file starts with the four dimensions (t, x, y, z) as 32-bit ints,
	* followed by the values as 64-bit doubles flattened using C-style
	* indexing.
	*
	* @tparam T Floating-point type of the returned values.
	* @param data_type Name of the dataset, which sets the file name. For
	*        species-specific data this is [species]_[data_type].
	* @return Values of the dataset with dimensions (t, x, y, z).
	*/
	template <typename T>
	Vectors::Vector4D<T> load_values(const std::string& data_type)
	{
		// Filename is treated as a constant.
		std::string filename {"bkg_from_pgkyl_" + data_type + ".bin"};
		std::ifstream fstream (filename, std::ios::binary);

		// First four values in file are 32-bit ints consisting of the
		// dimensions.
		int32_t tdim, xdim, ydim, zdim;
		fstream.read(reinterpret_cast<char*>(&tdim), sizeof(int32_t));
		fstream.read(reinterpret_cast<char*>(&xdim), sizeof(int32_t));
		fstream.read(reinterpret_cast<char*>(&ydim), sizeof(int32_t));
		fstream.read(reinterpret_cast<char*>(&zdim), sizeof(int32_t));

		// Next are the array values stored as 64-bit doubles. Calculate as
		// size_t to avoid potentially overflowing 32-bit int range.
		size_t num_vals = size_t(tdim) * xdim * ydim * zdim;

		// Load bytes into 1D vector.
		std::vector<double> values (num_vals);
		fstream.read(reinterpret_cast<char*>(values.data()),
			num_vals * sizeof(double));

		// If BkgFPType is not double, we need another vector of that type
		// to return. If it is already double, well I hope the compiler
		// can optimize this out. Not a huge deal.
		std::vector<BkgFPType> values_bkgtype (values.size());

		#pragma omp parallel for simd
		for (std::size_t i=0; i < values_bkgtype.size(); ++i)
			values_bkgtype[i] = static_cast<BkgFPType>(values[i]);

		// Return as a Vector4D
		return {values_bkgtype, tdim, xdim, ydim, zdim};
	}


	/**
	* Reads in the interpolation settings used by Gkeyll from
	* bkg_from_pgkyl_interp_settings.csv.
	*
	* This file is made automatically by read_gkyl.py when loading datasets
	* that include these settings, so they are available for files which do
	* not include them (like the magnetic field data).
	*
	* @return Vector of two elements: the basis type and the polynomial
	*         order.
	*/
	std::vector<std::string> load_interp_settings()
	{
		// Filename is treated as a constant.
		std::string filename {"bkg_from_pgkyl_interp_settings.csv"};
		std::ifstream fstream {filename};

		// The interpolation settings in a vector. First element is
		// basis type, second is the poly order.
		std::vector<std::string> interp_settings {};

		// Go through one line at a time
		std::string line {};
		while(std::getline(fstream, line))
		{
			// If line starts with #, then it's a comment and ignore
			if (!line.starts_with("#"))
			{
				interp_settings.push_back(line);
			}
		}

		return interp_settings;
	}


	/**
	* Writes a Vector4D to bkg_from_pgkyl_[data_type].csv in the same format
	* as write_csv in read_gkyl.py.
	*
	* @tparam T Data type of the values.
	* @param vec4d Values to write, with dimensions (t, x, y, z).
	* @param data_type Name of the dataset, which sets the file name.
	*/
	template <typename T>
	void write_values(const Vectors::Vector4D<T>& vec4d,
		const std::string& data_type)
	{
		// Filename is treated as a constant.
		std::string filename {"bkg_from_pgkyl_" + data_type + ".csv"};
		std::ofstream fstream {filename};

		// Write out the header info (it needs to match that in write_csv in
		// read_gkyl.py).
		fstream << "# The values for the specified type of data requested"
			<< " from\n";
		fstream << "# Gkeyll. The first row are integers: the number of"
			<< " frames,\n";
		fstream << "# and then the shape of each array (x, y, z dimensions)"
			<< ".\n";
		fstream << "# Then each array for each frame is printed out flattened"
			<< "\n";
		fstream << "# using C-style indexing.\n";

		// The 4 dimensions are expected next
		fstream << vec4d.get_dim1() << " " << vec4d.get_dim2() << " "
			<< vec4d.get_dim3() << " " << vec4d.get_dim4() << '\n';

		// Write out one value at at time. This assumes C-style flattening,
		// which is built into all parts of the code.
		for (const auto& val : vec4d.get_data())
		{
			fstream << val << '\n';
		}

		// Close file
		fstream.close();
	}


	/**
	* Calculates the cell centers from the grid edges.
	*
	* Each center is halfway between two neighboring edges.
	*
	* @param grid Grid edges for one dimension.
	* @return Cell centers, one shorter than grid.
	*/
	std::vector<BkgFPType> cell_centers(std::vector<BkgFPType>& grid)
	{
		std::vector<BkgFPType> centers (std::ssize(grid) - 1);
		for (std::size_t i {}; i < grid.size() - 1; ++i)
		{
			centers[i] = grid[i] + (grid[i+1] - grid[i]) / 2.0;
		}
		return centers;
	}


	/**
	* Replaces any values below a floor with the floor value.
	*
	* Useful for things like replacing negative temperatures with a
	* reasonable minimum. The number of replaced values is printed if any
	* were replaced.
	*
	* @tparam T Floating-point type of the values.
	* @param vec Values to modify in place.
	* @param floor_val Minimum allowed value.
	*/
	template <typename T>
	void set_vec_floor(Vectors::Vector4D<T>& vec, const double floor_val)
	{
		// Variable to keep track of number of replaced values
		int num_replaced {};

		// Loop through every value, replacing with floor_val if below it
		#pragma omp parallel for reduction(+:num_replaced)
		for (int i = 0; i < vec.get_dim1(); i++)
		{
		for (int j = 0; j < vec.get_dim2(); j++)
		{
		for (int k = 0; k < vec.get_dim3(); k++)
		{
		for (int l = 0; l < vec.get_dim4(); l++)
		{
			// If value is less than floor value, replace with floor value.
			if (vec(i,j,k,l) < floor_val)
			{
				vec(i,j,k,l) = floor_val;
				num_replaced++;
			}
		}
		}
		}
		}

		// Notify user of how many values were replaced
		if (num_replaced > 0) std::cout << "  " << num_replaced
			<< " values were below and replaced by " << floor_val << '\n';
	}


	/**
	* Reads in a Gkeyll dataset using a python interface to postgkyl,
	* read_gkyl.py.
	*
	* Species-specific data (density, temperature, par_flow and vperp_sq)
	* is saved to bkg_from_pgkyl_[species]_[data_type].bin, and everything
	* else to bkg_from_pgkyl_[data_type].bin. If that file already exists it
	* is loaded directly, unless force_load is true. Otherwise read_gkyl.py
	* is called to create it. The magnetic field and geometry files don't
	* include the DG interpolation settings, so for those the settings saved
	* when an earlier dataset was loaded are passed to read_gkyl.py.
	*
	* Either way, the times, grid edges and cell centers in grid_data are
	* also loaded.
	*
	* @tparam T Floating-point type of the values.
	* @param species Name of the Gkeyll species. Not used by data that isn't
	*        species-specific.
	* @param data_type Name of the dataset to load, e.g., density,
	*        temperature, potential or magnetic_magnitude.
	* @param grid_data Tuple of the time and grid vectors (see read_gkyl).
	* @param gkyl_data Array the values are moved into, with dimensions
	*        (t, x, y, z).
	* @param opts Simulation options.
	* @param species_mass_amu Mass of the species (amu), needed to calculate
	*        the temperature.
	* @param force_load If true, call read_gkyl.py even if the file already
	*        exists.
	*/
	template <typename T>
	void read_data_pgkyl(const std::string& species,
		const std::string& data_type, grid_data_t& grid_data,
		Vectors::Vector4D<T>& gkyl_data, const Options::Options& opts,
		const double species_mass_amu, const bool force_load)
	{
		// Check if file already exists and load it to save time. Can force
		// loading if we want to circumvent this. Account for the fact that
		// some files have the species name in them.
		std::string filename {};
		if (data_type == "density" || data_type == "temperature"
			|| data_type == "par_flow" || data_type == "vperp_sq")
		{
			//filename = "bkg_from_pgkyl_" + species + "_" + data_type + ".csv";
			filename = "bkg_from_pgkyl_" + species + "_" + data_type + ".bin";
		}
		else
		{
			//filename = "bkg_from_pgkyl_" + data_type + ".csv";
			filename = "bkg_from_pgkyl_" + data_type + ".bin";
		}

		if (std::filesystem::exists(filename) && !force_load)
		{
			//std::cout << "Previous file located\n";
			if (data_type == "density" || data_type == "temperature"
				|| data_type == "par_flow" || data_type == "vperp_sq")
			{
				Vectors::Vector4D<T> tmp_data {
					load_values<T>(species + "_" + data_type)};
				gkyl_data.move_into_data(tmp_data);
			}
			else
			{
				Vectors::Vector4D<T> tmp_data {load_values<T>(data_type)};
				gkyl_data.move_into_data(tmp_data);
			}

			// Load in grid data. Code copied from below in this function,
			// see it for the comments. I left them out for brevity.
			std::get<0>(grid_data) = {load_times()};  // gkyl_times
			std::tie(std::get<4>(grid_data), std::get<5>(grid_data),
				std::get<6>(grid_data)) = load_grid();
			std::get<1>(grid_data) = cell_centers(std::get<4>(grid_data));
			std::get<2>(grid_data) = cell_centers(std::get<5>(grid_data));
			std::get<3>(grid_data) = cell_centers(std::get<6>(grid_data));

			return;
		}

		// Ensure READ_GKYL_PY is defined, returning it if so.
		std::string read_gkyl_py {get_read_gkyl_py()};

		// Assemble command to call read_gkyl.py with stringstream
		std::stringstream read_gkyl_cmd_ss {};
		read_gkyl_cmd_ss << "python " << read_gkyl_py
			<< " --gkyl_dir=" << opts.gkyl_dir()
			<< " --gkyl_case_name=" << opts.gkyl_casename()
			<< " --gkyl_file_type=" << opts.gkyl_file_type()
			<< " --gkyl_frame_start=" << opts.gkyl_frame_start()
			<< " --gkyl_frame_end=" << opts.gkyl_frame_end()
			<< " --gkyl_species=" << species
			<< " --gkyl_species_mass_amu=" << species_mass_amu
			<< " --gkyl_data_type=" << data_type
			<< " --gkyl_moment_file_type=" << opts.gkyl_moment_type();

		// The bmag files don't include the needed DG interpolation data,
		// so we need to pass those in manually. To avoid user-error, we
		// will just rely on any of the other data file being created first,
		// where they will create a _interp_settings.csv file with the
		// basis type and poly order that we need to pass in.
		if (data_type == "magnetic_unit_X" ||
			data_type == "magnetic_unit_Y" ||
			data_type == "magnetic_unit_Z" ||
			data_type == "magnetic_magnitude" ||
			data_type == "jacobian" ||
			data_type == "metric_coeff00" ||
			data_type == "metric_coeff01" ||
			data_type == "metric_coeff02" ||
			data_type == "metric_coeff11" ||
			data_type == "metric_coeff12" ||
			data_type == "metric_coeff22" ||
			data_type == "tangent_basis_Xx" ||
			data_type == "tangent_basis_Yx" ||
			data_type == "tangent_basis_Zx" ||
			data_type == "tangent_basis_Xy" ||
			data_type == "tangent_basis_Yy" ||
			data_type == "tangent_basis_Zy" ||
			data_type == "tangent_basis_Xz" ||
			data_type == "tangent_basis_Yz" ||
			data_type == "tangent_basis_Zz" ||
			data_type == "reciprocal_basis_xX" ||
			data_type == "reciprocal_basis_xY" ||
			data_type == "reciprocal_basis_xZ" ||
			data_type == "reciprocal_basis_yX" ||
			data_type == "reciprocal_basis_yY" ||
			data_type == "reciprocal_basis_yZ" ||
			data_type == "reciprocal_basis_zX" ||
			data_type == "reciprocal_basis_zY" ||
			data_type == "reciprocal_basis_zZ" ||
			data_type == "covariant_comp_b_x" ||
			data_type == "covariant_comp_b_y" ||
			data_type == "covariant_comp_b_z")
		{
			std::vector<std::string> interp_settings = load_interp_settings();
			read_gkyl_cmd_ss << " --gkyl_basis_type=" << interp_settings[0]
				<< " --gkyl_poly_order=" << interp_settings[1];
		}

		// Execute the command to save the files
		std::string tmp_str {read_gkyl_cmd_ss.str()};
		const char* read_gkyl_cmd {tmp_str.c_str()};
		//std::cout << read_gkyl_cmd << '\n';
		[[maybe_unused]] int sys_result {system(read_gkyl_cmd)};

		// Load the times for each frame in a vector<double>
		//gkyl_times = {load_times()};
		std::get<0>(grid_data) = {load_times()};  // gkyl_times

		// Load the x, y, z grids each as their own vector<double>
		//std::tie(gkyl_grid_x, gkyl_grid_y, gkyl_grid_z) = load_grid();
		std::tie(std::get<4>(grid_data), std::get<5>(grid_data),
			std::get<6>(grid_data)) = load_grid();

		// Calculate the cell centers.
		//gkyl_x = cell_centers(gkyl_grid_x);
		//gkyl_y = cell_centers(gkyl_grid_y);
		//gkyl_z = cell_centers(gkyl_grid_z);
		std::get<1>(grid_data) = cell_centers(std::get<4>(grid_data));
		std::get<2>(grid_data) = cell_centers(std::get<5>(grid_data));
		std::get<3>(grid_data) = cell_centers(std::get<6>(grid_data));

		// Load the values at each t, x, y, z into a Vector4D. We use
		// a tmp variable here instead of passing into move_into_data
		// directly because move_into_data uses move semantics, which
		// cannot be used on a const lvalue reference. So we create
		// the tmp rvalue and then pass it in, where the data is then
		// moved from tmp into gkyl_data.
		//std::cout << "Loading " << data_type << "...\n";
		if (data_type == "density" || data_type == "temperature"
			|| data_type == "par_flow" || data_type == "vperp_sq")
		{
			// Some data is for specific species, which is going to be saved
			// as [species]_[data_type].
			Vectors::Vector4D<T> tmp_data {
				load_values<T>(species + "_" + data_type)};
			gkyl_data.move_into_data(tmp_data);
		}
		else
		{
			Vectors::Vector4D<T> tmp_data {load_values<T>(data_type)};
			gkyl_data.move_into_data(tmp_data);
		}
	}


	/**
	* Reads in the electron density (m^-3), replacing any values below
	* min_ne with min_ne.
	*
	* @tparam T Floating-point type of the values.
	* @param grid_data Tuple of the time and grid vectors (see read_gkyl).
	* @param gkyl_ne Array the density is moved into.
	* @param opts Simulation options.
	*/
	template <typename T>
	void read_elec_density(grid_data_t& grid_data,
		Vectors::Vector4D<T>& gkyl_ne, const Options::Options& opts)
	{
		// Call pgkyl to load data into gkyl_ne. We force load the density
		// since it gets called first, and we need at least one load to
		// load all the grid data.
		//read_data_pgkyl(opts.gkyl_elec_name(), "density", grid_data, gkyl_ne,
		//	opts, 0, true);

		// No longer force loading, as it can waste a large amount of time
		// for longer simulations. Delete above commented code when ready.
		read_data_pgkyl(opts.gkyl_elec_name(), "density", grid_data, gkyl_ne,
			opts, 0);
		set_vec_floor(gkyl_ne, opts.min_ne());
	}


	/**
	* Reads in the electron temperature (eV), replacing any values below
	* min_te with min_te.
	*
	* @tparam T Floating-point type of the values.
	* @param grid_data Tuple of the time and grid vectors (see read_gkyl).
	* @param gkyl_te Array the temperature is moved into.
	* @param opts Simulation options.
	*/
	template <typename T>
	void read_elec_temperature(grid_data_t& grid_data,
		Vectors::Vector4D<T>& gkyl_te, const Options::Options& opts)
	{
		// Call pgkyl to load data into gkyl_te. Need to pass in the species
		// mass so it can calculate temperature correctly.
		read_data_pgkyl(opts.gkyl_elec_name(), "temperature", grid_data,
			gkyl_te, opts, opts.gkyl_elec_mass_amu());
		set_vec_floor(gkyl_te, opts.min_te());
	}


	/**
	* Reads in the ion temperature (eV), replacing any values below min_ti
	* with min_ti.
	*
	* @tparam T Floating-point type of the values.
	* @param grid_data Tuple of the time and grid vectors (see read_gkyl).
	* @param gkyl_ti Array the temperature is moved into.
	* @param opts Simulation options.
	*/
	template <typename T>
	void read_ion_temperature(grid_data_t& grid_data,
		Vectors::Vector4D<T>& gkyl_ti, const Options::Options& opts)
	{
		// Call pgkyl to load data into gkyl_ti. Need to pass in the species
		// mass so it can calculate temperature correctly.
		read_data_pgkyl(opts.gkyl_ion_name(), "temperature", grid_data,
			gkyl_ti, opts, opts.gkyl_ion_mass_amu());
		set_vec_floor(gkyl_ti, opts.min_ti());
	}


	/**
	* Reads in the plasma potential (V).
	*
	* @tparam T Floating-point type of the values.
	* @param grid_data Tuple of the time and grid vectors (see read_gkyl).
	* @param gkyl_vp Array the potential is moved into.
	* @param opts Simulation options.
	*/
	template <typename T>
	void read_potential(grid_data_t& grid_data,
		Vectors::Vector4D<T>& gkyl_vp, const Options::Options& opts)
	{
		using namespace std::string_literals;

		// Call pgkyl to load data into gkyl_vp. There is no species
		// name with the potential file so just putting a placeholder
		// string. The script handles things.
		read_data_pgkyl("null"s, "potential", grid_data, gkyl_vp, opts);
	}


	/**
	* Reads in the magnetic field magnitude (T). Currently unused, since
	* read_magnetic_field loads it.
	*
	* @tparam T Floating-point type of the values.
	* @param grid_data Tuple of the time and grid vectors (see read_gkyl).
	* @param gkyl_bmag Array the magnitude is moved into.
	* @param opts Simulation options.
	*/
	template <typename T>
	void read_magnetic_magnitude(grid_data_t& grid_data,
		Vectors::Vector4D<T>& gkyl_bmag, const Options::Options& opts)
	{
		using namespace std::string_literals;

		// Call pgkyl to load data into gkyl_bmag. There is no species
		// name with the magnetic field magnitude file so just putting a
		// placeholder string. The script handles things.
		read_data_pgkyl("null"s, "magnetic_magnitude", grid_data, gkyl_bmag,
			opts);
	}


	/**
	* Reads in the covariant components of the magnetic field.
	*
	* @tparam T Floating-point type of the values.
	* @param grid_data Tuple of the time and grid vectors (see read_gkyl).
	* @param gkyl_b_x Array the x component is moved into.
	* @param gkyl_b_y Array the y component is moved into.
	* @param gkyl_b_z Array the z component is moved into.
	* @param opts Simulation options.
	*/
	template <typename T>
	void read_covariant_comp_b(grid_data_t& grid_data,
		Vectors::Vector4D<T>& gkyl_b_x, Vectors::Vector4D<T>& gkyl_b_y,
		Vectors::Vector4D<T>& gkyl_b_z, const Options::Options& opts)
	{
		using namespace std::string_literals;

		read_data_pgkyl("null"s, "covariant_comp_b_x", grid_data, gkyl_b_x,
			opts);
		read_data_pgkyl("null"s, "covariant_comp_b_y", grid_data, gkyl_b_y,
			opts);
		read_data_pgkyl("null"s, "covariant_comp_b_z", grid_data, gkyl_b_z,
			opts);
	}


	/**
	* Reads in the magnetic field and calculates its Cartesian components.
	*
	* Gkeyll provides the Cartesian components of the magnetic field unit
	* vector and the magnitude. The unit vector is not guaranteed to have a
	* magnitude of 1.0 at the cell centers, so each one is renormalized to
	* 1.0 and then multiplied by the magnitude. The number of unit vectors
	* with a magnitude outside [0.95, 1.05] is printed. bR is the magnitude
	* of the X-Y components.
	*
	* @tparam T Floating-point type of the values.
	* @param grid_data Tuple of the time and grid vectors (see read_gkyl).
	* @param gkyl_bX Array the X component (T) is moved into.
	* @param gkyl_bY Array the Y component (T) is moved into.
	* @param gkyl_bZ Array the Z component (T) is moved into.
	* @param gkyl_bR Array sqrt(bX^2 + bY^2) (T) is stored in.
	* @param gkyl_bmag Array the magnitude (T) is moved into.
	* @param opts Simulation options.
	*/
	template <typename T>
	void read_magnetic_field(grid_data_t& grid_data,
		Vectors::Vector4D<T>& gkyl_bX, Vectors::Vector4D<T>& gkyl_bY,
		Vectors::Vector4D<T>& gkyl_bZ, Vectors::Vector4D<T>& gkyl_bR,
		Vectors::Vector4D<T>& gkyl_bmag, const Options::Options& opts)
	{
		using namespace std::string_literals;

		// There is no species name with the magnetic field files so just
		// putting a placeholder string. The script handles things.
		// These are actually the unit vector components of the magnetic field
		// so call once for each component. We multiply by the magnitude
		// later.
		read_data_pgkyl("null"s, "magnetic_unit_X", grid_data, gkyl_bX, opts);
		read_data_pgkyl("null"s, "magnetic_unit_Y", grid_data, gkyl_bY, opts);
		read_data_pgkyl("null"s, "magnetic_unit_Z", grid_data, gkyl_bZ, opts);

		// Resize space for R component of magnetic field
		read_data_pgkyl("null"s, "magnetic_unit_X", grid_data, gkyl_bR, opts);
		//gkyl_bR = Vectors::Vector4D<T>(gkyl_bX.get_dim1(), gkyl_bX.get_dim2(),
		//	gkyl_bX.get_dim3(), gkyl_bX.get_dim4());

		// Load the magnetic field magnitude at each location
		//Vectors::Vector4D<BkgFPType> gkyl_bmag {};
		read_data_pgkyl("null"s, "magnetic_magnitude", grid_data, gkyl_bmag,
			opts);

		constexpr double unit_tol_min {0.95};
		constexpr double unit_tol_max {1.05};
		int outside_tol {};

		int tdim {static_cast<int>(std::get<0>(grid_data).size())};
		int xdim {static_cast<int>(std::get<1>(grid_data).size())};
		int ydim {static_cast<int>(std::get<2>(grid_data).size())};
		int zdim {static_cast<int>(std::get<3>(grid_data).size())};

		// The unit vector is actually a bit finnicky on the gkyl side of
		// things. It's not actually guaranteed to be 1.0 magnitude at the
		// cell centers which I guess is where it's calculated. So right now we
		// just hope the relative sizes of each component are mostly correct,
		// renormalize to 1.0 unit vector magnitude before multiplying by bmag.
		#pragma omp parallel for reduction(+:outside_tol)
		for (int i = 0; i < tdim; ++i)  // t
		{
		for (int j {}; j < xdim; ++j)  // x
		{
		for (int k {}; k < ydim; ++k)  // y
		{
		for (int l {}; l < zdim; ++l)  // z
		{
			// Unit vector magnitude. Ideally 1.0, but life is not always so
			// gentle. If the unit vector is just straight up wrong, then this
			// throws off the entire simulation. Need to ask a Gkeyll
			// developer if it's going to be correct or not.
			double unit_mag {std::sqrt(gkyl_bX(i,j,k,l) * gkyl_bX(i,j,k,l)
				+ gkyl_bY(i,j,k,l) * gkyl_bY(i,j,k,l)
				+ gkyl_bZ(i,j,k,l) * gkyl_bZ(i,j,k,l))};

			// Record if outside tolerance range.
			if (unit_mag < unit_tol_min || unit_mag > unit_tol_max)
				outside_tol += 1;

			// Renormalize each component, and multiply by the magnetic
			// field magnitude for each component.
			gkyl_bX(i,j,k,l) = gkyl_bX(i,j,k,l) * gkyl_bmag(i,j,k,l)
				/ unit_mag;
			gkyl_bY(i,j,k,l) = gkyl_bY(i,j,k,l) * gkyl_bmag(i,j,k,l)
				/ unit_mag;
			gkyl_bZ(i,j,k,l) = gkyl_bZ(i,j,k,l) * gkyl_bmag(i,j,k,l)
				/ unit_mag;

			// R component of magnetic field
			gkyl_bR(i,j,k,l) = std::sqrt(gkyl_bX(i,j,k,l) * gkyl_bX(i,j,k,l)
				+ gkyl_bY(i,j,k,l) * gkyl_bY(i,j,k,l));
		}
		}
		}
		}

		// Let user know how many were out of bounds. TBD if this is useful.
		if (outside_tol)
		{
			double outside_tol_perc {static_cast<double>(outside_tol)
				/ static_cast<double>(tdim* xdim * ydim * zdim) * 100};
			std::cout << "  " << outside_tol << " values ("
				<< outside_tol_perc << "%)  outside of [" << unit_tol_min
				<< ", " << unit_tol_max << "]\n";
			std::cout << "  This may be fine...\n";
		}
	}


	/**
	* Calculates the Cartesian components of the magnetic field from its
	* covariant components and the reciprocal basis vectors.
	*
	* Currently unused, since it is not clear if the magnitude is baked into
	* the covariant components or not (read_magnetic_field is used instead).
	* The geometry is taken from the first time slice, and the result is
	* multiplied by the magnitude, which is loaded here.
	*
	* @tparam T Floating-point type of the values.
	* @param grid_data Tuple of the time and grid vectors (see read_gkyl).
	* @param gkyl_bX Array the X component (T) is stored in.
	* @param gkyl_bY Array the Y component (T) is stored in.
	* @param gkyl_bZ Array the Z component (T) is stored in.
	* @param gkyl_bR Array sqrt(bX^2 + bY^2) (T) is stored in.
	* @param gkyl_b_x Covariant x component of the magnetic field.
	* @param gkyl_b_y Covariant y component of the magnetic field.
	* @param gkyl_b_z Covariant z component of the magnetic field.
	* @param gkyl_dxdX Reciprocal basis vector component dx/dX.
	* @param gkyl_dxdY Reciprocal basis vector component dx/dY.
	* @param gkyl_dxdZ Reciprocal basis vector component dx/dZ.
	* @param gkyl_dydX Reciprocal basis vector component dy/dX.
	* @param gkyl_dydY Reciprocal basis vector component dy/dY.
	* @param gkyl_dydZ Reciprocal basis vector component dy/dZ.
	* @param gkyl_dzdX Reciprocal basis vector component dz/dX.
	* @param gkyl_dzdY Reciprocal basis vector component dz/dY.
	* @param gkyl_dzdZ Reciprocal basis vector component dz/dZ.
	* @param opts Simulation options.
	*/
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
		const Options::Options& opts)
	{
		using namespace std::string_literals;

		// Load the magnetic field magnitude. Only used internally to scale
		// the Cartesian unit vector components we calculate
		Vectors::Vector4D<BkgFPType> gkyl_bmag {};
		read_data_pgkyl("null"s, "magnetic_magnitude", grid_data, gkyl_bmag,
			opts);

		// There is an inconsistency in the b_i.gkyl file. Sometimes it is
		// normalized in curvilinear space, such that it needs to be multiplied
		// by bmag for the actual B value. Other times it already includes
		// bmag and is not a unit vector. We try to detect that here so we
		// know if we need to multiply by bmag or not in the loop.

		// Create zero-filled Vector4Ds for each magnetic field component
		int tdim {static_cast<int>(std::get<0>(grid_data).size())};
		int xdim {static_cast<int>(std::get<1>(grid_data).size())};
		int ydim {static_cast<int>(std::get<2>(grid_data).size())};
		int zdim {static_cast<int>(std::get<3>(grid_data).size())};
		gkyl_bX = Vectors::Vector4D<BkgFPType>(tdim, xdim, ydim, zdim);
		gkyl_bY = Vectors::Vector4D<BkgFPType>(tdim, xdim, ydim, zdim);
		gkyl_bZ = Vectors::Vector4D<BkgFPType>(tdim, xdim, ydim, zdim);
		gkyl_bR = Vectors::Vector4D<BkgFPType>(tdim, xdim, ydim, zdim);

		// Note we are baking in the assumption that the grid does not change,
		// but the magnetic field magnitude is still free to change.
		#pragma omp parallel for
		for (int i = 0; i < tdim; ++i)  // t
		{
		for (int j {}; j < xdim; ++j)  // x
		{
		for (int k {}; k < ydim; ++k)  // y
		{
		for (int l {}; l < zdim; ++l)  // z
		{
			// Magnitude. I am still unsure if the magnetic field needs to be
			// multiplied by bmag or not! At some point it seemed like the
			// right thing to do, but now it seems like it isn't.
			const double bmag {gkyl_bmag(i,j,k,l)};
			//constexpr double bmag {1.0};

			// Calculate Cartesian components of magnetic field. We do this
			// by calculating the dot product between the reciprocal basis
			// vector (dz/dX) and the magnetic field covariant components (b_i).
			// BX
			gkyl_bX(i,j,k,l) = (gkyl_b_x(0,j,k,l) * gkyl_dxdX(0,j,k,l)
				+ gkyl_b_y(0,j,k,l) * gkyl_dydX(0,j,k,l)
				+ gkyl_b_z(0,j,k,l) * gkyl_dzdX(0,j,k,l)) * bmag;

			// BY
			gkyl_bY(i,j,k,l) = (gkyl_b_x(0,j,k,l) * gkyl_dxdY(0,j,k,l)
				+ gkyl_b_y(0,j,k,l) * gkyl_dydY(0,j,k,l)
				+ gkyl_b_z(0,j,k,l) * gkyl_dzdY(0,j,k,l)) * bmag;

			// BZ
			gkyl_bZ(i,j,k,l) = (gkyl_b_x(0,j,k,l) * gkyl_dxdZ(0,j,k,l)
				+ gkyl_b_y(0,j,k,l) * gkyl_dydZ(0,j,k,l)
				+ gkyl_b_z(0,j,k,l) * gkyl_dzdZ(0,j,k,l)) * bmag;

			// BR
			gkyl_bR(i,j,k,l) = std::sqrt(gkyl_bX(i,j,k,l) * gkyl_bX(i,j,k,l)
				+ gkyl_bY(i,j,k,l) * gkyl_bY(i,j,k,l));
		}
		}
		}
		}
	}


	/**
	* Reads in the Jacobian of the coordinate system.
	*
	* @tparam T Floating-point type of the values.
	* @param grid_data Tuple of the time and grid vectors (see read_gkyl).
	* @param gkyl_J Array the Jacobian is moved into.
	* @param opts Simulation options.
	*/
	template <typename T>
	void read_jacobian(grid_data_t& grid_data,
		Vectors::Vector4D<T>& gkyl_J, const Options::Options& opts)
	{
		using namespace std::string_literals;

		// Call pgkyl to load data into gkyl_J. There is no species
		// name with the Jacobian file so just putting a placeholder
		// string. The script handles things.
		read_data_pgkyl("null"s, "jacobian", grid_data, gkyl_J, opts);
	}


	/**
	* Reads in one component of the metric coefficients, g_ij.
	*
	* @tparam T Floating-point type of the values.
	* @param grid_data Tuple of the time and grid vectors (see read_gkyl).
	* @param gkyl_gij Array the component is moved into.
	* @param opts Simulation options.
	* @param idx Which component to read: 00, 01, 02, 11, 12 or 22.
	*/
	template <typename T>
	void read_gij(grid_data_t& grid_data,
		Vectors::Vector4D<T>& gkyl_gij, const Options::Options& opts,
		const std::string& idx)
	{
		using namespace std::string_literals;

		// Assign the correct string for the index needed. idx is one of
		// 00, 01, 02, 11, 12, 22.
		std::string metric_coeff_str {"metric_coeff"};
		metric_coeff_str += idx;

		read_data_pgkyl("null"s, metric_coeff_str, grid_data, gkyl_gij, opts);
	}


	/**
	* Reads in one component of the tangent (covariant) basis vectors,
	* e.g., dX/dx.
	*
	* @tparam T Floating-point type of the values.
	* @param grid_data Tuple of the time and grid vectors (see read_gkyl).
	* @param gkyl_dxdz Array the component is moved into.
	* @param opts Simulation options.
	* @param idx Which component to read, as [physical][computational],
	*        e.g., Xx for dX/dx.
	*/
	template <typename T>
	void read_tangent_basis(grid_data_t& grid_data,
		Vectors::Vector4D<T>& gkyl_dxdz, const Options::Options& opts,
		const std::string& idx)
	{
		using namespace std::string_literals;

		// Assign the correct string for the vector needed. idx is one of
		// Xx, Yx, Zx, Xy, Yy, Zy, Xz, Yz, Zz.
		std::string tangent_basis_str {"tangent_basis_"};
		tangent_basis_str += idx;

		read_data_pgkyl("null"s, tangent_basis_str, grid_data,
			gkyl_dxdz, opts);
	}


	/**
	* Reads in one component of the reciprocal (contravariant) basis
	* vectors, e.g., dx/dX.
	*
	* @tparam T Floating-point type of the values.
	* @param grid_data Tuple of the time and grid vectors (see read_gkyl).
	* @param gkyl_dzdx Array the component is moved into.
	* @param opts Simulation options.
	* @param idx Which component to read, as
	*        [computational][physical], e.g., xX for dx/dX.
	*/
	template <typename T>
	void read_reciprocal_basis(grid_data_t& grid_data,
		Vectors::Vector4D<T>& gkyl_dzdx, const Options::Options& opts,
		const std::string& idx)
	{
		using namespace std::string_literals;

		// Assign the correct string for the vector needed. idx is one of
		// xX, xY, xZ, yX, yY, yZ, zX, zY, zZ.
		std::string reciprocal_basis_str {"reciprocal_basis_"};
		reciprocal_basis_str += idx;

		read_data_pgkyl("null"s, reciprocal_basis_str, grid_data,
			gkyl_dzdx, opts);
	}


	/**
	* Returns the full path to calc_elec_field.py, which calculates the
	* electric field from the potential gradient.
	*
	* The path is stored in the CALC_ELEC_FIELD_PY environment variable. If
	* it is not set, the run is stopped with Errors::fatal_error. Currently
	* unused, since calc_gradients calculates the electric field directly.
	*
	* @return Full path to calc_elec_field.py.
	*/
	std::string get_calc_elec_field_py()
	{
		if (const char* calc_elec_field_py {std::getenv("CALC_ELEC_FIELD_PY")})
		{
			return static_cast<std::string>(calc_elec_field_py);
		}
		else
		{
			Errors::fatal_error("CALC_ELEC_FIELD_PY is not set. Please set "
				"this to the full path to the calc_elec_field.py script. "
				"E.g.:\n"
				"  CALC_ELEC_FIELD_PY=/path/to/flan/python/"
				"calc_elec_field.py");
		}
	}


	/**
	* Calculates the electric field and the magnetic field gradient in
	* Cartesian coordinates.
	*
	* The gradients of the potential and the magnetic field magnitude are
	* first calculated in computational coordinates with central differences
	* (forward or backward differences at the edges), assuming a uniform
	* grid. They are then converted to Cartesian components using the
	* reciprocal basis vectors (Eq. 2.6.35 in D'haeseleer), with
	* E = -grad(Vp). The basis vectors are taken from the first time slice.
	*
	* @tparam T Floating-point type of the basis vectors.
	* @param grid_data Tuple of the time and grid vectors (see read_gkyl).
	* @param gkyl_vp Plasma potential (V).
	* @param gkyl_eX Array the electric field X component (V/m) is stored
	*        in.
	* @param gkyl_eY Array the electric field Y component (V/m) is stored
	*        in.
	* @param gkyl_eZ Array the electric field Z component (V/m) is stored
	*        in.
	* @param gkyl_emag Array the electric field magnitude (V/m) is stored in.
	* @param gkyl_bmag Magnetic field magnitude (T).
	* @param gkyl_dbdX Array the X component of grad(B) (T/m) is stored in.
	* @param gkyl_dbdY Array the Y component of grad(B) (T/m) is stored in.
	* @param gkyl_dbdZ Array the Z component of grad(B) (T/m) is stored in.
	* @param gkyl_dxdX Reciprocal basis vector component dx/dX.
	* @param gkyl_dxdY Reciprocal basis vector component dx/dY.
	* @param gkyl_dxdZ Reciprocal basis vector component dx/dZ.
	* @param gkyl_dydX Reciprocal basis vector component dy/dX.
	* @param gkyl_dydY Reciprocal basis vector component dy/dY.
	* @param gkyl_dydZ Reciprocal basis vector component dy/dZ.
	* @param gkyl_dzdX Reciprocal basis vector component dz/dX.
	* @param gkyl_dzdY Reciprocal basis vector component dz/dY.
	* @param gkyl_dzdZ Reciprocal basis vector component dz/dZ.
	* @param opts Simulation options.
	*/
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
		const Options::Options& opts)
	{

		// Calculate electric field and magnetic field gradient using the
		// gradient in computational (curvilinear) coordinates and the
		// reciprocal (contravariant) vectors

		// Resize each Vector4D to the dimensions of Vp
		gkyl_eX.resize(gkyl_vp.get_dim1(), gkyl_vp.get_dim2(),
			gkyl_vp.get_dim3(), gkyl_vp.get_dim4());
		gkyl_eY.resize(gkyl_vp.get_dim1(), gkyl_vp.get_dim2(),
			gkyl_vp.get_dim3(), gkyl_vp.get_dim4());
		gkyl_eZ.resize(gkyl_vp.get_dim1(), gkyl_vp.get_dim2(),
			gkyl_vp.get_dim3(), gkyl_vp.get_dim4());
		gkyl_emag.resize(gkyl_vp.get_dim1(), gkyl_vp.get_dim2(),
			gkyl_vp.get_dim3(), gkyl_vp.get_dim4());
		gkyl_dbdX.resize(gkyl_vp.get_dim1(), gkyl_vp.get_dim2(),
			gkyl_vp.get_dim3(), gkyl_vp.get_dim4());
		gkyl_dbdY.resize(gkyl_vp.get_dim1(), gkyl_vp.get_dim2(),
			gkyl_vp.get_dim3(), gkyl_vp.get_dim4());
		gkyl_dbdZ.resize(gkyl_vp.get_dim1(), gkyl_vp.get_dim2(),
			gkyl_vp.get_dim3(), gkyl_vp.get_dim4());

		// Calculate dx, dy, dz, hence assuming a uniform grid, but this
		// can easily be generalized into arbitrary spacing if needed.
		double dx {std::get<1>(grid_data)[1] - std::get<1>(grid_data)[0]};
		double dy {std::get<2>(grid_data)[1] - std::get<2>(grid_data)[0]};
		double dz {std::get<3>(grid_data)[1] - std::get<3>(grid_data)[0]};

		// Loop through each index
		int tdim {static_cast<int>(std::get<0>(grid_data).size())};
		int xdim {static_cast<int>(std::get<1>(grid_data).size())};
		int ydim {static_cast<int>(std::get<2>(grid_data).size())};
		int zdim {static_cast<int>(std::get<3>(grid_data).size())};

		#pragma omp parallel for
		for (int i = 0; i < tdim; ++i)  // t
		{
		for (int j {}; j < xdim; ++j)  // x
		{
		for (int k {}; k < ydim; ++k)  // y
		{
		for (int l {}; l < zdim; ++l)  // z
		{
			// Calculate x,y,z partial derivative for Vp and B. Using
			// central difference formula, reverting to forward or
			// backward for the edge cases.
			//std::cout << i << " " << j << " " << k << " " << l << ": "
			//	<< "Vp = " << gkyl_vp(i,j,k,l) << '\n';

			// x
			double dVpdx {};
			double dBdx {};
			if (j == 0)
			{
				dVpdx = (gkyl_vp(i,j+1,k,l) - gkyl_vp(i,j,k,l)) / dx;
				dBdx = (gkyl_bmag(i,j+1,k,l) - gkyl_bmag(i,j,k,l)) / dx;
			}
			else if (j == (xdim-1))
			{
				dVpdx = (gkyl_vp(i,j,k,l) - gkyl_vp(i,j-1,k,l)) / dx;
				dBdx = (gkyl_bmag(i,j,k,l) - gkyl_bmag(i,j-1,k,l)) / dx;
			}
			else
			{
				dVpdx = (gkyl_vp(i,j+1,k,l) - gkyl_vp(i,j-1,k,l))
					/ (2 * dx);
				dBdx = (gkyl_bmag(i,j+1,k,l) - gkyl_bmag(i,j-1,k,l))
					/ (2 * dx);
			}

			// y
			double dVpdy {};
			double dBdy {};
			if (k == 0)
			{
				dVpdy = (gkyl_vp(i,j,k+1,l) - gkyl_vp(i,j,k,l)) / dy;
				dBdy = (gkyl_bmag(i,j,k+1,l) - gkyl_bmag(i,j,k,l)) / dy;
			}
			else if (k == (ydim-1))
			{
				dVpdy = (gkyl_vp(i,j,k,l) - gkyl_vp(i,j,k-1,l)) / dy;
				dBdy = (gkyl_bmag(i,j,k,l) - gkyl_bmag(i,j,k-1,l)) / dy;
			}
			else
			{
				dVpdy = (gkyl_vp(i,j,k+1,l) - gkyl_vp(i,j,k-1,l))
					/ (2 * dy);
				dBdy = (gkyl_bmag(i,j,k+1,l) - gkyl_bmag(i,j,k-1,l))
					/ (2 * dy);
			}

			// z
			double dVpdz {};
			double dBdz {};
			if (l == 0)
			{
				dVpdz = (gkyl_vp(i,j,k,l+1) - gkyl_vp(i,j,k,l)) / dz;
				dBdz = (gkyl_bmag(i,j,k,l+1) - gkyl_bmag(i,j,k,l)) / dz;
			}
			else if (l == (zdim-1))
			{
				dVpdz = (gkyl_vp(i,j,k,l) - gkyl_vp(i,j,k,l-1)) / dz;
				dBdz = (gkyl_bmag(i,j,k,l) - gkyl_bmag(i,j,k,l-1)) / dz;
			}
			else
			{
				dVpdz = (gkyl_vp(i,j,k,l+1) - gkyl_vp(i,j,k,l-1))
					/ (2 * dz);
				dBdz = (gkyl_bmag(i,j,k,l+1) - gkyl_bmag(i,j,k,l-1))
					/ (2 * dz);
			}

			// Assemble Cartesian X,Y,Z gradient. This is Eq. 2.6.35 in
			// D'haeseleer.
			// Warning! I've found a landmine that you need to be careful
			// with, a result of something slipping through the cracks in how
			// data is read in. The basis vectors are actually 3D data that's
			// stored in 4D arrays. This ends up meaning only i=0 has the data
			// in it, and it's undefined everywhere else! I need to convert
			// these to 3D arrays.
			double dVpdX {dVpdx * gkyl_dxdX(0,j,k,l)
				+ dVpdy * gkyl_dydX(0,j,k,l) + dVpdz * gkyl_dzdX(0,j,k,l)};
			double dVpdY {dVpdx * gkyl_dxdY(0,j,k,l)
				+ dVpdy * gkyl_dydY(0,j,k,l) + dVpdz * gkyl_dzdY(0,j,k,l)};
			double dVpdZ {dVpdx * gkyl_dxdZ(0,j,k,l)
				+ dVpdy * gkyl_dydZ(0,j,k,l) + dVpdz * gkyl_dzdZ(0,j,k,l)};
			double dBdX {dBdx * gkyl_dxdX(0,j,k,l)
				+ dBdy * gkyl_dydX(0,j,k,l) + dBdz * gkyl_dzdX(0,j,k,l)};
			double dBdY {dBdx * gkyl_dxdY(0,j,k,l)
				+ dBdy * gkyl_dydY(0,j,k,l) + dBdz * gkyl_dzdY(0,j,k,l)};
			double dBdZ {dBdx * gkyl_dxdZ(0,j,k,l)
				+ dBdy * gkyl_dydZ(0,j,k,l) + dBdz * gkyl_dzdZ(0,j,k,l)};

			// Store in Vector4D, not forgetting the minus sign!
			gkyl_eX(i,j,k,l) = -dVpdX;
			gkyl_eY(i,j,k,l) = -dVpdY;
			gkyl_eZ(i,j,k,l) = -dVpdZ;
			gkyl_dbdX(i,j,k,l) = dBdX;
			gkyl_dbdY(i,j,k,l) = dBdY;
			gkyl_dbdZ(i,j,k,l) = dBdZ;

			// Also store electric field magnitude
			gkyl_emag(i,j,k,l) = std::sqrt(dVpdX * dVpdX + dVpdY * dVpdY
				+ dVpdZ * dVpdZ);
		}
		}
		}
		}
	}


	/**
	* Calculates the Cartesian X, Y, Z coordinates of each cell center.
	*
	* Currently this only resizes gkyl_X, gkyl_Y and gkyl_Z in grid_data.
	* The coordinates were calculated with the mapc2p function, which is
	* planned to be removed, so that part is commented out.
	*
	* @param grid_data Tuple of the time and grid vectors (see read_gkyl).
	* @param opts Simulation options. Currently unused.
	*/
	void calc_cell_XYZ_centers(grid_data_t& grid_data,
		const Options::Options& opts)
	{
		//  grid_data[0]  = gkyl_times
		//  grid_data[1]  = gkyl_x
		//  grid_data[2]  = gkyl_y
		//  grid_data[3]  = gkyl_z
		//  grid_data[4]  = gkyl_grid_x
		//  grid_data[5]  = gkyl_grid_y
		//  grid_data[6]  = gkyl_grid_z
		//  grid_data[7]  = gkyl_X
		//  grid_data[8]  = gkyl_Y
		//  grid_data[9]  = gkyl_Z

		// Resize the arrays since we know by now how big they need to be
		int dim1 {static_cast<int>(std::ssize(std::get<1>(grid_data)))};
		int dim2 {static_cast<int>(std::ssize(std::get<2>(grid_data)))};
		int dim3 {static_cast<int>(std::ssize(std::get<3>(grid_data)))};
		std::get<7>(grid_data).resize(dim1, dim2, dim3); // gkyl_X
		std::get<8>(grid_data).resize(dim1, dim2, dim3); // gkyl_Y
		std::get<9>(grid_data).resize(dim1, dim2, dim3); // gkyl_Z

		// Loop through every x, y, z value to calculate each X, Y, Z
		#pragma omp parallel for
		for (int i = 0; i < dim1; ++i)
		{
		for (int j {}; j < dim2; ++j)
		{
		for (int k {}; k < dim3; ++k)
		{
			// Planning to remove mapc2p so likely will delete

			// Get the three Cartesian coordinates and store within
			// our 3D vector of 3-tuples
			//auto [X, Y, Z] = opts.mapc2p()(std::get<1>(grid_data)[i],
			//	std::get<2>(grid_data)[j], std::get<3>(grid_data)[k]);

			//std::get<7>(grid_data)(i,j,k) = X;
			//std::get<8>(grid_data)(i,j,k) = Y;
			//std::get<9>(grid_data)(i,j,k) = Z;
		}
		}
		}
	}


	/**
	* Calculates the Cartesian X, Y, Z coordinates of each grid edge.
	*
	* Currently this only resizes gkyl_grid_X, gkyl_grid_Y and gkyl_grid_Z
	* in grid_data. The coordinates were calculated with the mapc2p
	* function, which is planned to be removed, so that part is commented
	* out.
	*
	* @param grid_data Tuple of the time and grid vectors (see read_gkyl).
	* @param opts Simulation options. Currently unused.
	*/
	void calc_cell_XYZ_edges(grid_data_t& grid_data,
		const Options::Options& opts)
	{
		//  grid_data[0]  = gkyl_times
		//  grid_data[1]  = gkyl_x
		//  grid_data[2]  = gkyl_y
		//  grid_data[3]  = gkyl_z
		//  grid_data[4]  = gkyl_grid_x
		//  grid_data[5]  = gkyl_grid_y
		//  grid_data[6]  = gkyl_grid_z
		//  grid_data[7]  = gkyl_X
		//  grid_data[8]  = gkyl_Y
		//  grid_data[9]  = gkyl_Z
		//  grid_data[10] = gkyl_grid_X
		//  grid_data[11] = gkyl_grid_Y
		//  grid_data[12] = gkyl_grid_Z

		// Resize the arrays since we know by now how big they need to be
		int dim1 {static_cast<int>(std::ssize(std::get<4>(grid_data)))};
		int dim2 {static_cast<int>(std::ssize(std::get<5>(grid_data)))};
		int dim3 {static_cast<int>(std::ssize(std::get<6>(grid_data)))};
		std::get<10>(grid_data).resize(dim1, dim2, dim3); // gkyl_grid_X
		std::get<11>(grid_data).resize(dim1, dim2, dim3); // gkyl_grid_Y
		std::get<12>(grid_data).resize(dim1, dim2, dim3); // gkyl_grid_Z

		// Loop through every x, y, z value to calculate each X, Y, Z
		#pragma omp parallel for
		for (int i = 0; i < dim1; ++i)
		{
		//std::cout << i+1 << "/" << dim1 << '\n';
		for (int j {}; j < dim2; ++j)
		{
		for (int k {}; k < dim3; ++k)
		{
			// Get the three Cartesian coordinates and store within
			// our 3D vector of 3-tuples
			//auto [X, Y, Z] = opts.mapc2p()(std::get<4>(grid_data)[i],
			//	std::get<5>(grid_data)[j], std::get<6>(grid_data)[k]);
			//std::get<10>(grid_data)(i,j,k) = X;
			//std::get<11>(grid_data)(i,j,k) = Y;
			//std::get<12>(grid_data)(i,j,k) = Z;

			//std::cout << "X, Y, Z, R: " << X << ", " << Y << ", "
			//	<< Z << ", " << std::sqrt(X*X+Y*Y) << '\n';
		}
		}
		}
	}


	/**
	* Moves gradient-calculated data into the background arrays.
	*
	* The fix to the y boundaries in gradient_ybound_fix used to be applied
	* here, but it is currently commented out.
	*
	* @param grid_data Tuple of the time and grid vectors (see read_gkyl).
	*        Currently unused.
	* @param dataX X component of the calculated data.
	* @param dataY Y component of the calculated data.
	* @param dataZ Z component of the calculated data.
	* @param gkyl_dataX Array the X component is moved into.
	* @param gkyl_dataY Array the Y component is moved into.
	* @param gkyl_dataZ Array the Z component is moved into.
	*/
	void read_gradient_data(grid_data_t& grid_data,
		Vectors::Vector4D<BkgFPType>& dataX,
		Vectors::Vector4D<BkgFPType>& dataY,
		Vectors::Vector4D<BkgFPType>& dataZ,
		Vectors::Vector4D<BkgFPType>& gkyl_dataX,
		Vectors::Vector4D<BkgFPType>& gkyl_dataY,
		Vectors::Vector4D<BkgFPType>& gkyl_dataZ)
	{
		// Apply fix to y boundaries where the gradient generally messes up
		//gradient_ybound_fix(grid_data, dataX);
		//gradient_ybound_fix(grid_data, dataY);
		//gradient_ybound_fix(grid_data, dataZ);

		// Move into arrays
		gkyl_dataX.move_into_data(dataX);
		gkyl_dataY.move_into_data(dataY);
		gkyl_dataZ.move_into_data(dataZ);
	}


	/**
	* Applies a fix to the y boundaries of gradient-calculated data, where
	* the gradient generally messes up due to not as many points around.
	*
	* The y boundaries are treated as periodic, and the first and last cells
	* in y are replaced with weighted averages of the second and second to
	* last cells. This is generally a hack, since it is not clear what to do
	* in this situation. Ideally, the gradient calculation would handle this,
	* but I haven't been able to figure it out. Currently unused.
	*
	* @param grid_data Tuple of the time and grid vectors (see read_gkyl).
	* @param ecomp Gradient-calculated data to fix in place.
	*/
	void gradient_ybound_fix(grid_data_t& grid_data,
		Vectors::Vector4D<BkgFPType>& ecomp)
	{
		// The following comment is for the electric field, but the same
		// logic applies to any other gradient-calculated data (e.g., the
		// magnetic field gradient).
		//
		// We apply a fix here. The y boundaries are treated as periodic (in
		// Gkeyll this need not always be true!). At the edge of the boundaries
		// the gradient calculation for the electric field can get messed up
		// because it does not have as many cells to calculate a gradient from.
		// So we mitigate this by assigning the edge cells a weighted average
		// from the "bounding" values in the y direction.
		// Consider if y has 6 cells:
		//
		//    |-----|-----|-----|
		// iy 4     5     0     1
		//             ^
		//          end of grid here, periodic
		//          so we restart at iy=0
		//
		// The problem cells are 0 and 5. We (arbitrarily) assign those as:
		//  iy=0   Ex = (1/3) * Ex_4 + (2/3) * Ex_1
		//  iy=5   Ex = (2/3) * Ex_4 + (1/3) * Ex_1
		// So we are taking the bounding two "reliable" E values at iy=1,4 and
		// reassigning the values at iy=0,5 with the weighted averages since
		// it is not clear what else we can do in this situation.
		int leny {static_cast<int>(std::ssize(std::get<2>(grid_data)))};
		#pragma omp parallel for
		for (int i = 0; i < std::ssize(std::get<0>(grid_data)); i++)  // t
		{
		for (int j {}; j < std::ssize(std::get<1>(grid_data)); j++)  // x
		{
		for (int l {}; l < std::ssize(std::get<3>(grid_data)); l++)  // z
		{
			// Just naming these variables according to the example above.
			BkgFPType e1 {ecomp(i,j,1,l)};
			BkgFPType e4 {ecomp(i,j,leny-2,l)};
			ecomp(i,j,0,l) = (1.0/3.0) * e4 + (2.0/3.0) * e1;
			ecomp(i,j,leny-1,l) = (2.0/3.0) * e4 + (1.0/3.0) * e1;
		}
		}
		}
	}


	/**
	* Reads in the Cartesian components of the electric field (V/m)
	* calculated by calc_elec_field.py, from the
	* bkg_from_pgkyl_elec_field_X/Y/Z.bin files.
	*
	* Currently unused, since calc_gradients calculates these directly.
	*
	* @param grid_data Tuple of the time and grid vectors (see read_gkyl).
	* @param gkyl_eX Array the X component is moved into.
	* @param gkyl_eY Array the Y component is moved into.
	* @param gkyl_eZ Array the Z component is moved into.
	*/
	void read_elec_field(grid_data_t& grid_data,
		Vectors::Vector4D<BkgFPType>& gkyl_eX,
		Vectors::Vector4D<BkgFPType>& gkyl_eY,
		Vectors::Vector4D<BkgFPType>& gkyl_eZ)
	{
		// Ex
		Vectors::Vector4D<BkgFPType> tmp_dataX
			{load_values<BkgFPType>("elec_field_X")};

		// Ey
		Vectors::Vector4D<BkgFPType> tmp_dataY
			{load_values<BkgFPType>("elec_field_Y")};

		// Ez
		Vectors::Vector4D<BkgFPType> tmp_dataZ
			{load_values<BkgFPType>("elec_field_Z")};

		// Call general function to handle calculated-gradient data.
		read_gradient_data(grid_data, tmp_dataX, tmp_dataY, tmp_dataZ,
			gkyl_eX, gkyl_eY, gkyl_eZ);

		// Apply fix to y boundaries where the gradient generally messes up
		/*
		elec_field_ybound_fix(grid_data, tmp_dataX);
		elec_field_ybound_fix(grid_data, tmp_dataY);
		elec_field_ybound_fix(grid_data, tmp_dataZ);

		// Move into arrays
		gkyl_eX.move_into_data(tmp_dataX);
		gkyl_eY.move_into_data(tmp_dataY);
		gkyl_eZ.move_into_data(tmp_dataZ);
		*/
	}


	/**
	* Reads in the Cartesian components of the electric field gradient
	* calculated by calc_elec_field.py, from the
	* bkg_from_pgkyl_grad_elec_field_X/Y/Z.bin files.
	*
	* Currently unused, since calc_gradients calculates these directly.
	*
	* @param grid_data Tuple of the time and grid vectors (see read_gkyl).
	* @param gkyl_gradeX Array the X component is moved into.
	* @param gkyl_gradeY Array the Y component is moved into.
	* @param gkyl_gradeZ Array the Z component is moved into.
	*/
	void read_grad_elec_field(grid_data_t& grid_data,
		Vectors::Vector4D<BkgFPType>& gkyl_gradeX,
		Vectors::Vector4D<BkgFPType>& gkyl_gradeY,
		Vectors::Vector4D<BkgFPType>& gkyl_gradeZ)
	{
		// grad Ex
		Vectors::Vector4D<BkgFPType> tmp_dataX
			{load_values<BkgFPType>("grad_elec_field_X")};

		// grad Ey
		Vectors::Vector4D<BkgFPType> tmp_dataY
			{load_values<BkgFPType>("grad_elec_field_Y")};

		// grad Ez
		Vectors::Vector4D<BkgFPType> tmp_dataZ
			{load_values<BkgFPType>("grad_elec_field_Z")};

		// Call general function to handle calculated-gradient data.
		read_gradient_data(grid_data, tmp_dataX, tmp_dataY, tmp_dataZ,
			gkyl_gradeX, gkyl_gradeY, gkyl_gradeZ);
	}


	/**
	* Reads in the Cartesian components of the magnetic field gradient (T/m)
	* calculated by calc_elec_field.py, from the
	* bkg_from_pgkyl_gradb_X/Y/Z.bin files.
	*
	* Currently unused, since calc_gradients calculates these directly.
	*
	* @param grid_data Tuple of the time and grid vectors (see read_gkyl).
	* @param gkyl_dbdX Array the X component is moved into.
	* @param gkyl_dbdY Array the Y component is moved into.
	* @param gkyl_dbdZ Array the Z component is moved into.
	*/
	void read_gradb(grid_data_t& grid_data,
		Vectors::Vector4D<BkgFPType>& gkyl_dbdX,
		Vectors::Vector4D<BkgFPType>& gkyl_dbdY,
		Vectors::Vector4D<BkgFPType>& gkyl_dbdZ)
	{
		// gradB_X
		Vectors::Vector4D<BkgFPType> tmp_dataX
			{load_values<BkgFPType>("gradb_X")};

		// gradB_Y
		Vectors::Vector4D<BkgFPType> tmp_dataY
			{load_values<BkgFPType>("gradb_Y")};

		// gradB_Z
		Vectors::Vector4D<BkgFPType> tmp_dataZ
			{load_values<BkgFPType>("gradb_Z")};

		// Call general function to handle calculated-gradient data.
		read_gradient_data(grid_data, tmp_dataX, tmp_dataY, tmp_dataZ,
			gkyl_dbdX, gkyl_dbdY, gkyl_dbdZ);
	}


	/**
	* Writes the Cartesian X, Y, Z coordinates of the cell centers to
	* bkg_from_pgkyl_cell_XYZ.csv.
	*
	* This was used to read them into a python script and take advantage of
	* the scipy library to calculate a gradient on an irregular grid.
	* Currently unused.
	*
	* @param grid_data Tuple of the time and grid vectors (see read_gkyl).
	* @param opts Simulation options. Currently unused.
	*/
	void write_XYZ(grid_data_t& grid_data, const Options::Options& opts)
	{
		//  grid_data[0]  = gkyl_times
		//  grid_data[1]  = gkyl_x
		//  grid_data[2]  = gkyl_y
		//  grid_data[3]  = gkyl_z
		//  grid_data[4]  = gkyl_grid_x
		//  grid_data[5]  = gkyl_grid_y
		//  grid_data[6]  = gkyl_grid_z
		//  grid_data[7]  = gkyl_X
		//  grid_data[8]  = gkyl_Y
		//  grid_data[9]  = gkyl_Z

		// Open up file for writing
		std::ofstream XYZ_file {"bkg_from_pgkyl_cell_XYZ.csv"};
		if (!XYZ_file)
		{
			Errors::fatal_error("Could not create "
				"bkg_from_pgkyl_cell_XYZ.csv.");
		}

		// Print header info
		XYZ_file << "# This file contains the Cartesian X,Y,Z coordinates of\n"
			<< "# the Gkeyll grid at each cell center. They are\n"
			<< "# calculated using the mapc2p function passed in from\n"
			<< "# the input file. The first line in this file are the\n"
			<< "# dimensions of the data. The coordinates\n"
			<< "# are then printed one at a time with X Y Z on each \n"
			<< "# line.\n";

		// Print out shape of data
		XYZ_file << std::get<7>(grid_data).get_dim1() << " "
			<< std::get<7>(grid_data).get_dim2() << " "
			<< std::get<7>(grid_data).get_dim3() << '\n';

		// Print out the coordinates
		for (int i {}; i < std::ssize(std::get<1>(grid_data)); ++i)
		{
		for (int j {}; j < std::ssize(std::get<2>(grid_data)); ++j)
		{
		for (int k {}; k < std::ssize(std::get<3>(grid_data)); ++k)
		{
			XYZ_file << std::get<7>(grid_data)(i,j,k) << " "
				<< std::get<8>(grid_data)(i,j,k) << " "
				<< std::get<9>(grid_data)(i,j,k) << '\n';
		}
		}
		}
	}


	/**
	* Reads in the ion flow parallel to z (m/s). The cross-field flows are
	* drifts and are added in calc_bkg_flow.
	*
	* @tparam T Floating-point type of the values.
	* @param grid_data Tuple of the time and grid vectors (see read_gkyl).
	* @param gkyl_uz Array the flow is moved into.
	* @param opts Simulation options.
	*/
	template <typename T>
	void read_par_flow(grid_data_t& grid_data,
		Vectors::Vector4D<T>& gkyl_uz, const Options::Options& opts)
	{
		// Call pgkyl to load data into gkyl_uz. This loads just the parallel
		// to z flow. The cross-field flows are drifts and added elsewhere.
		read_data_pgkyl(opts.gkyl_ion_name(), "par_flow", grid_data,
			gkyl_uz, opts);
	}


	/**
	* Reads in the ion perpendicular thermal velocity squared (m^2/s^2).
	* This requires BiMaxwellian moments.
	*
	* @tparam T Floating-point type of the values.
	* @param grid_data Tuple of the time and grid vectors (see read_gkyl).
	* @param gkyl_viperp_sq Array the values are moved into.
	* @param opts Simulation options.
	*/
	template <typename T>
	void read_ion_vperp_sq(grid_data_t& grid_data,
		Vectors::Vector4D<T>& gkyl_viperp_sq, const Options::Options& opts)
	{
		// Call pgkyl to load data into gkyl_viperp2. This requires
		// BiMaxwellian moments and thus will fail if they are not there.
		read_data_pgkyl(opts.gkyl_ion_name(), "vperp_sq", grid_data,
			gkyl_viperp_sq, opts);
	}


	/**
	* Calculates the Cartesian components of the background ion flow.
	*
	* The flow is the projection of the parallel-to-z flow onto the
	* magnetic field direction, plus the ExB drift, E x B / B^2, plus the
	* grad-B drift, (m v_perp^2 / (2 q B^3)) B x grad(B). The ions are
	* assumed to be singly charged.
	*
	* @param grid_data Tuple of the time and grid vectors (see read_gkyl).
	* @param gkyl_uz Ion flow parallel to z (m/s).
	* @param gkyl_uX Array the X component of the flow (m/s) is stored in.
	* @param gkyl_uY Array the Y component of the flow (m/s) is stored in.
	* @param gkyl_uZ Array the Z component of the flow (m/s) is stored in.
	* @param gkyl_bX Magnetic field X component (T).
	* @param gkyl_bY Magnetic field Y component (T).
	* @param gkyl_bZ Magnetic field Z component (T).
	* @param gkyl_eX Electric field X component (V/m).
	* @param gkyl_eY Electric field Y component (V/m).
	* @param gkyl_eZ Electric field Z component (V/m).
	* @param gkyl_dbdX X component of grad(B) (T/m).
	* @param gkyl_dbdY Y component of grad(B) (T/m).
	* @param gkyl_dbdZ Z component of grad(B) (T/m).
	* @param gkyl_viperp_sq Ion perpendicular thermal velocity squared
	*        (m^2/s^2).
	* @param opts Simulation options, used for the ion mass.
	*/
	void calc_bkg_flow(grid_data_t& grid_data,
		Vectors::Vector4D<BkgFPType>& gkyl_uz,
		Vectors::Vector4D<BkgFPType>& gkyl_uX,
		Vectors::Vector4D<BkgFPType>& gkyl_uY,
		Vectors::Vector4D<BkgFPType>& gkyl_uZ,
		const Vectors::Vector4D<BkgFPType>& gkyl_bX,
		const Vectors::Vector4D<BkgFPType>& gkyl_bY,
		const Vectors::Vector4D<BkgFPType>& gkyl_bZ,
		Vectors::Vector4D<BkgFPType>& gkyl_eX,
		Vectors::Vector4D<BkgFPType>& gkyl_eY,
		Vectors::Vector4D<BkgFPType>& gkyl_eZ,
		Vectors::Vector4D<BkgFPType>& gkyl_dbdX,
		Vectors::Vector4D<BkgFPType>& gkyl_dbdY,
		Vectors::Vector4D<BkgFPType>& gkyl_dbdZ,
		Vectors::Vector4D<BkgFPType>& gkyl_viperp_sq,
		const Options::Options& opts)
	{
		// Resize the uX, uY, uZ vectors
		gkyl_uX.resize(gkyl_uz.get_dim1(), gkyl_uz.get_dim2(),
			gkyl_uz.get_dim3(), gkyl_uz.get_dim4());
		gkyl_uY.resize(gkyl_uz.get_dim1(), gkyl_uz.get_dim2(),
			gkyl_uz.get_dim3(), gkyl_uz.get_dim4());
		gkyl_uZ.resize(gkyl_uz.get_dim1(), gkyl_uz.get_dim2(),
			gkyl_uz.get_dim3(), gkyl_uz.get_dim4());

		// We are calculating uX = (parallel flow projection) +
		// (ExB drift) + (gradB drift) = uz * bX + ExB_X + gradB_X
		#pragma omp parallel for
		for (int i = 0; i < std::ssize(std::get<0>(grid_data)); ++i)  // t
		{
		for (int j {}; j < std::ssize(std::get<1>(grid_data)); ++j)  // x
		{
		for (int k {}; k < std::ssize(std::get<2>(grid_data)); ++k)  // y
		{
		for (int l {}; l < std::ssize(std::get<3>(grid_data)); ++l)  // z
		{
			// Some variables for cleaner code
			double bX {gkyl_bX(i,j,k,l)};
			double bY {gkyl_bY(i,j,k,l)};
			double bZ {gkyl_bZ(i,j,k,l)};
			double gradbX {gkyl_dbdX(i,j,k,l)};
			double gradbY {gkyl_dbdY(i,j,k,l)};
			double gradbZ {gkyl_dbdZ(i,j,k,l)};
			double eX {gkyl_eX(i,j,k,l)};
			double eY {gkyl_eY(i,j,k,l)};
			double eZ {gkyl_eZ(i,j,k,l)};
			double viperp_sq {gkyl_viperp_sq(i,j,k,l)};

			// Magnitude of magnetic field, used a couple times here
			double bmag_sq {bX * bX + bY * bY + bZ * bZ};
			double bmag {std::sqrt(bmag_sq)};

			// Projection of uz onto Cartesian magnetic field components
			gkyl_uX(i,j,k,l) = gkyl_uz(i,j,k,l) * bX / bmag;
			gkyl_uY(i,j,k,l) = gkyl_uz(i,j,k,l) * bY / bmag;
			gkyl_uZ(i,j,k,l) = gkyl_uz(i,j,k,l) * bZ / bmag;

			// ExB drift component
			gkyl_uX(i,j,k,l) = gkyl_uX(i,j,k,l) + (eY * bZ - eZ * bY)
				/ bmag_sq;
			gkyl_uY(i,j,k,l) = gkyl_uY(i,j,k,l) + (eZ * bX - eX * bZ)
				/ bmag_sq;
			gkyl_uZ(i,j,k,l) = gkyl_uZ(i,j,k,l) + (eX * bY - eY * bX)
				/ bmag_sq;

			// gradB drift component. Group terms out front into coefficient.
			// Assume singly charged ions (q = e).
			double ion_charge {-Constants::charge_e};
			double coef {opts.gkyl_ion_mass_amu() * Constants::amu_to_kg
				* viperp_sq / (2.0 * ion_charge * bmag_sq * bmag)};
			gkyl_uX(i,j,k,l) = gkyl_uX(i,j,k,l) + coef * (bY * gradbZ
				- bZ * gradbY);
			gkyl_uY(i,j,k,l) = gkyl_uY(i,j,k,l) + coef * (bZ * gradbX
				- bX * gradbZ);
			gkyl_uZ(i,j,k,l) = gkyl_uZ(i,j,k,l) + coef * (bX * gradbY
				- bY * gradbX);
		}
		}
		}
		}
	}


	/**
	* Checks that the number of frames in the loaded background matches that
	* requested for the simulation, stopping the run with
	* Errors::fatal_error if not.
	*
	* This only needs to be run once. It is useful if the user is reusing an
	* already calculated background, but changed the frames in the input
	* file. In that case the background files should be deleted so that they
	* can be regenerated, since the default is to read any existing files.
	* Note you can circumvent this error checking if the range of frames
	* stays the same, but that's a fringe case and not worth the effort to
	* protect against at the moment.
	*
	* @param gkyl_arr Any loaded array with dimensions (t, x, y, z).
	* @param opts Simulation options, used for gkyl_frame_start and
	*        gkyl_frame_end.
	*/
	void check_tdim(Vectors::Vector4D<BkgFPType>& gkyl_arr,
		const Options::Options& opts)
	{
		// Load trange from user input.
		const int tstart {opts.gkyl_frame_start()};
		const int tend {opts.gkyl_frame_end()};
		const int trange {tend - tstart + 1}; // Inclusive

		if (trange != gkyl_arr.get_dim1())
		{
			Errors::fatal_error("Number of time frames in loaded background "
				"(" + std::to_string(gkyl_arr.get_dim1()) + ") is different "
				"from the requested range (" + std::to_string(trange) + "). "
				"Delete all background files and rerun.");
		}
	}
}
