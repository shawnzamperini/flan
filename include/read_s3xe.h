/**
* @file read_s3xe.h
*
* @brief Header file for read_s3xe.cpp
*/

#pragma once

#include <array>
#include <string>
#include <vector>

#include "background.h"
#include "flan_types.h"
#include "vectors.h"
#include "zones.h"

namespace S3X
{
	/**
	* Sign relating the SOLEDGE3X toroidal angle phi to the Cartesian
	* frame, assuming (R, Z, phi) is right-handed. Then phi increases
	* clockwise seen from above, the opposite of the usual (R, phi, Z), so
	* a point is at (X, Y, Z) = (R cos phi, phi_sign R sin phi, Z) with
	* phi_sign = -1. This still needs to be confirmed (see
	* print_assumptions).
	*/
	inline constexpr double phi_sign {-1.0};

	/// Ghost cells on each side of every zone, the same as the plasma files
	inline constexpr int ghost {2};

	/**
	* Reference values SOLEDGE3X output is normalized by, from its
	* reference file (e.g., refParam_raptorX.h5). The reference mass is
	* 1 amu, so c0 = sqrt(e T0 / 1 amu) and tau0 = 1 amu / (e B0), not
	* values for the main ion.
	*/
	struct Reference
	{
		double n0 {};	///< Density (m^-3)
		double T0 {};	///< Temperature (eV)
		double c0 {};	///< Velocity (m/s)
		double tau0 {};	///< Time (s)
		double phi0 {};	///< Potential (V)
		double R0 {};	///< Major radius (m) R is offset by
		double a0 {};	///< Length (m) R and Z are normalized by
		double B0 {};	///< Magnetic field (T)
	};

	/**
	* The parts of a SOLEDGE3X mesh file Flan uses, per zone, with the zones
	* built from them. Lengths and fields are still normalized as in the
	* file. Per-zone arrays are flattened in the file's C order: corners are
	* [psi][theta] with Npsi + 1 by Ntheta + 1 values, and cell values are
	* [psi][theta] with Npsi by Ntheta values.
	*/
	struct Mesh
	{
		double R0 {};		///< Major radius (m) R is offset by
		double a0 {};		///< Length (m) R and Z are normalized by
		double B0 {};		///< Magnetic field (T) B is normalized by
		Zones::Zones zones {};	///< Zones with grids and faces
		std::vector<std::vector<double>> R_corners {};
		std::vector<std::vector<double>> Z_corners {};
		std::vector<std::vector<double>> Br {};
		std::vector<std::vector<double>> Bz {};
		std::vector<std::vector<double>> Bphi {};
		std::vector<std::vector<double>> chi {};
	};

	void print_assumptions();

	Reference read_reference(const std::string& ref_path);

	void print_reference(const Reference& ref);

	void check_reference(const Reference& ref, const Mesh& mesh);

	int read_nphi(const std::string& plasma_path, const int nzones);

	Mesh read_mesh(const std::string& mesh_path, const int nphi);

	std::array<double, 3> derivatives(const Zones::Zones& zones,
		const Vectors::ZonedVector<BkgFPType>& f, const int zone,
		const int t, const int i, const int j, const int k);

	std::array<double, 3> gradient(const Background::Background& bkg,
		const Vectors::ZonedVector<BkgFPType>& f, const int zone,
		const int t, const int i, const int j, const int k);

	void calc_geometry(const Mesh& mesh, Background::Background& bkg);

	void calc_bfield(const Mesh& mesh, const int nframes,
		Background::Background& bkg);

	void load_chi(const Mesh& mesh, Background::Background& bkg);
}
