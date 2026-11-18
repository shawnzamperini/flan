/**
* @file axisym_mesh.h
*
* @brief Header file for axisym_mesh.cpp: the geometry, magnetic field,
* electric field and flows of a multi-zone background on an axisymmetric
* mesh, shared by the SOLEDGE3X and HERMES readers.
*/

#pragma once

#include <array>
#include <cmath>
#include <string>
#include <utility>
#include <vector>

#include "background.h"
#include "braginskii.h"
#include "flan_types.h"
#include "vectors.h"
#include "zones.h"

namespace Axisym
{
	using Vec3 = std::array<double, 3>;

	/// a x b
	inline Vec3 cross(const Vec3& a, const Vec3& b)
	{
		return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2],
			a[0] * b[1] - a[1] * b[0]};
	}

	/// a . b
	inline double dot(const Vec3& a, const Vec3& b)
	{
		return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
	}

	/// a * u + b * v
	inline Vec3 combine(const double a, const Vec3& u, const double b,
		const Vec3& v)
	{
		return {a * u[0] + b * v[0], a * u[1] + b * v[1],
			a * u[2] + b * v[2]};
	}

	/**
	* Unit vector in R at toroidal angle phi in the Cartesian frame, where a
	* point is at (R cos phi, phi_sign R sin phi, Z).
	*
	* @param phi Toroidal angle (rad).
	* @param phi_sign +1 if (R, phi, Z) is right handed, -1 if not.
	* @return R_hat.
	*/
	inline Vec3 R_hat(const double phi, const double phi_sign)
		{return {std::cos(phi), phi_sign * std::sin(phi), 0.0};}

	/**
	* Unit vector in phi (the direction phi increases in) at toroidal angle
	* phi in the Cartesian frame. See R_hat.
	*
	* @param phi Toroidal angle (rad).
	* @param phi_sign +1 if (R, phi, Z) is right handed, -1 if not.
	* @return phi_hat.
	*/
	inline Vec3 phi_hat(const double phi, const double phi_sign)
		{return {-std::sin(phi), phi_sign * std::cos(phi), 0.0};}

	inline constexpr Vec3 Z_hat {0.0, 0.0, 1.0};

	/**
	* A poloidal-plane mesh, the same in every toroidal plane, in SI units.
	* Each zone's arrays are flattened in C order over (x, y): corners are
	* [i][j] with (nx + 1) by (ny + 1) values, where corner (i, j) is at the
	* low-x, low-y corner of cell (i, j), and cell values are [i][j] with nx
	* by ny values. The zones themselves (sizes, z grid and faces) are kept
	* in the background.
	*/
	struct Mesh
	{
		std::vector<std::vector<double>> R_corners {};	///< R (m)
		std::vector<std::vector<double>> Z_corners {};	///< Z (m)
		std::vector<std::vector<double>> Br {};		///< B_R (T)
		std::vector<std::vector<double>> Bz {};		///< B_Z (T)
		std::vector<std::vector<double>> Bphi {};	///< B_phi (T)

		/// A point at (R, phi, Z) is at (R cos phi, phi_sign R sin phi, Z):
		/// +1 if (R, phi, Z) is right handed, -1 if it's left handed.
		double phi_sign {1.0};
	};

	Vectors::ZonedVector<BkgFPType> zoned(const Zones::Zones& zones,
		const int dim1);

	std::pair<std::vector<double>, std::vector<double>> uniform_grid(
		const int n, const double lo, const double hi);

	std::pair<int, int> y_stencil(const Zones::Zones& zones, const int zone,
		const int j);

	std::array<double, 3> derivatives(const Zones::Zones& zones,
		const Vectors::ZonedVector<BkgFPType>& f, const int zone,
		const int t, const int i, const int j, const int k,
		const Vectors::ZonedVector<BkgFPType>* dfdy = nullptr);

	std::array<double, 3> gradient(const Background::Background& bkg,
		const Vectors::ZonedVector<BkgFPType>& f, const int zone,
		const int t, const int i, const int j, const int k,
		const Vectors::ZonedVector<BkgFPType>* dfdy = nullptr);

	double grad_par(const Background::Background& bkg,
		const Vectors::ZonedVector<BkgFPType>& dfdy, const int zone,
		const int t, const int i, const int j, const int k);

	std::vector<double> shift_periodic(const std::vector<double>& row,
		const double shift);

	Vectors::ZonedVector<BkgFPType> aligned_dfdy(const Zones::Zones& zones,
		const Vectors::ZonedVector<BkgFPType>& f);

	void calc_geometry(const Mesh& mesh, Background::Background& bkg);

	void calc_bfield(const Mesh& mesh, const int nframes,
		Background::Background& bkg);

	void report_ghosts(const Zones::Zones& zones,
		const Vectors::ZonedVector<BkgFPType>& v, const std::string& name);

	void calc_efield_and_flows(Background::Background& bkg,
		const Vectors::ZonedVector<BkgFPType>& vp,
		const Vectors::ZonedVector<BkgFPType>& upar);

	/// Settings of the conductive heat flux (see calc_conductive_heat_fluxes)
	struct Conduction
	{
		/// Electron and ion kappa coefficients
		double coef_e {Braginskii::kappa_coef_e_default};
		double coef_i {Braginskii::kappa_coef_i_default};
		double limit_alpha {-1.0};	///< Flux limiter, <= 0 for none
		double ion_z {1.0};			///< Main ion charge
		double ion_aa {2.0};		///< Main ion mass (amu)
	};

	void calc_conductive_heat_fluxes(Background::Background& bkg,
		const Conduction& c);
}
