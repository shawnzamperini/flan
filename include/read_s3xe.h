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
#include "hdf5_utils.h"
#include "options.h"
#include "vectors.h"
#include "zones.h"

namespace S3X
{
	/**
	* Sign relating the SOLEDGE3X toroidal angle phi to the Cartesian
	* frame. SOLEDGE3X's (R, Z, phi) is right-handed, so phi increases
	* clockwise seen from above, the opposite of the usual (R, phi, Z), and
	* a point is at (X, Y, Z) = (R cos phi, phi_sign R sin phi, Z) with
	* phi_sign = -1.
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

		/// Wall outlines (R, Z) in meters, from /walls (empty if the file
		/// has none). Unlike the corners, these aren't normalized.
		std::vector<Background::WallPolygon> walls {};
	};

	Reference read_reference(const std::string& ref_path);

	void print_reference(const Reference& ref);

	void check_reference(const Reference& ref, const Mesh& mesh);

	int read_nphi(const std::string& plasma_path);

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

	void report_ghosts(const Zones::Zones& zones,
		const Vectors::ZonedVector<BkgFPType>& v, const std::string& name);

	bool has_phi_face_fluxes(const Hdf5::File& file,
		const std::string& species);

	double calc_q_par(const Hdf5::File& file, const Reference& ref,
		const Mesh& mesh, const Zones::Zones& zones,
		const std::string& species, const double mass_amu, const int t,
		Vectors::ZonedVector<BkgFPType>& q_par, const bool strict = false);

	void read_frames(const std::vector<std::string>& frame_paths,
		const Reference& ref, const Mesh& mesh,
		Background::Background& bkg, const bool strict = false);

	void calc_efield_and_flows(Background::Background& bkg,
		const Vectors::ZonedVector<BkgFPType>& vp,
		const Vectors::ZonedVector<BkgFPType>& upar);

	Background::Background read_s3xe(const Options::Options& opts);
}
