/**
* @file read_slab.cpp
*
* @brief Routines that make the slab background, a uniform plasma set by
*        input options
*
* The slab background (bkg_source = "slab") is a box of uniform plasma in
* a simple magnetic field, with everything about the plasma and the fields
* given by the slab_* input options. The field is uniform by default; an
* electric field, a field gradient and curvature (by bending the box into
* part of a torus) can be added for the drifts. It exists so the physics in Flan can be checked
* against analytic results: a case is fully defined by its input file,
* with nothing about it hard-coded here, and the checks can be repeated at
* other densities, temperatures and so on. The cases in regression/ that
* start with reg_ use it.
*/
#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include "background.h"
#include "constants.h"
#include "errors.h"
#include "flan_types.h"
#include "options.h"
#include "read_bkg.h"
#include "read_slab.h"
#include "utilities.h"
#include "vectors.h"
#include "zones.h"

namespace Slab
{
	/**
	* Creates the slab background from the slab_* options.
	*
	* Straight slab (slab_major_radius = 0, the default)
	* --------------------------------------------------
	* A 5 cm cube, x from 0 to 0.05 m and y and z from -0.025 to 0.025 m,
	* with x = X, y = Y and z = Z (Cartesian, so the Jacobian is 1), in
	* 8 x 8 x 8 cells. The fields are uniform or linear, so the grid only
	* needs to be fine enough to interpolate on. For particles to stay in
	* the box, make the boundaries periodic in the input file.
	*
	* - Magnetic field along +X, B_X = slab_bmag (1 + slab_bgrad_y Y). The
	*   gradient is across the field and gives a grad-B drift. With a
	*   gradient the field isn't periodic in y, so leave those boundaries
	*   absorbing and keep the particles away from them.
	* - Electric field along +Y, E_Y = slab_efield_y + slab_efield_y_rate
	*   t, uniform in space. It gives an ExB drift, and a polarization
	*   drift when it changes in time.
	* - The field can be turned in the X-Y plane by slab_b_angle (a) from
	*   +X toward -Y: B = slab_bmag (cos a, -sin a, 0). Then there can be
	*   no electric field or gradient, since E_Y would have a part along B
	*   and the gradient would make div B nonzero.
	* - One face can be a material wall (slab_wall), such as a divertor
	*   target, in place of what its *bound_type gives. With the wall at
	*   y min, slab_b_angle is the angle between the field and the wall.
	*   The face across from it can't be periodic, or particles would wrap
	*   onto the wall.
	*
	* The zone is marked Cartesian (Zones::Zones::set_cartesian), so the
	* wall can be sputtered.
	*
	* Bent slab (slab_major_radius = R0 > 0)
	* --------------------------------------
	* The box is bent into a quarter of a torus, to give the magnetic field
	* curvature: cylindrical coordinates with x = R from R0 - 0.025 to
	* R0 + 0.025 m, y = Z from -0.025 to 0.025 m and z = phi (the toroidal
	* angle) from 0 to pi/2, in 8 x 8 x 32 cells. Cartesian X = R cos(phi),
	* Y = R sin(phi), and the Jacobian is R.
	*
	* - Magnetic field in the toroidal direction (+phi), B = slab_bmag R0
	*   / R, as for the vacuum field of a tokamak. Its curvature and
	*   gradient give the curvature and grad-B drifts.
	* - No electric field or extra gradient: the run is stopped if
	*   slab_efield_y, slab_efield_y_rate or slab_bgrad_y isn't 0.
	*
	* Both
	* ----
	* There are 128 time frames from 0 to slab_tmax. The frames set the
	* time resolution of the output.
	*
	* - Density slab_ne (electrons and ions), temperatures slab_te and
	*   slab_ti.
	* - Flow slab_flow along the magnetic field.
	* - Parallel conductive heat fluxes q_par = slab_heat_flux_i n T_i
	*   v_th,i for the ions and slab_heat_flux_e n T_e v_th,e for the
	*   electrons, with v_th = sqrt(T / m) of that species. So the
	*   options are the size of the heat-flux skew directly (see
	*   Collisions::skew_eps), which is what gives the thermal force. If
	*   both are 0 the background has no heat fluxes at all, as for a
	*   Gkeyll background, and can be run on the GPU.
	*
	* The ion and electron masses are gkyl_ion_mass_amu and
	* gkyl_elec_mass_amu. The run is stopped if a density, temperature,
	* field strength or slab_tmax isn't above 0, or if the field strength
	* would reach 0 inside the box because of slab_bgrad_y.
	*
	* @param opts Simulation options, used for the slab_* options and the
	*        masses.
	* @return Background holding the slab plasma.
	*/
	Background::Background read_slab(const Options::Options& opts)
	{
		// Values that must be above 0
		auto need_positive = [](const double val, const std::string& name)
		{
			if (!(val > 0.0))
			{
				Errors::fatal_error(name + " must be above 0, got "
					+ std::to_string(val) + ".");
			}
		};
		need_positive(opts.slab_ne(), "slab_ne");
		need_positive(opts.slab_te(), "slab_te");
		need_positive(opts.slab_ti(), "slab_ti");
		need_positive(opts.slab_bmag(), "slab_bmag");
		need_positive(opts.slab_tmax(), "slab_tmax");

		// Straight, or bent to this major radius
		const double R0 {opts.slab_major_radius()};
		const bool bent {R0 > 0.0};
		constexpr double half_width {0.025};
		if (R0 < 0.0 || (bent && R0 <= half_width))
		{
			Errors::fatal_error("slab_major_radius must be 0 (straight "
				"slab) or above 0.025 m, got " + std::to_string(R0) + ".");
		}
		if (bent && (opts.slab_efield_y() != 0.0
			|| opts.slab_efield_y_rate() != 0.0
			|| opts.slab_bgrad_y() != 0.0))
		{
			Errors::fatal_error("slab_efield_y, slab_efield_y_rate and "
				"slab_bgrad_y are only for a straight slab "
				"(slab_major_radius = 0).");
		}
		// Field angle and wall
		const double b_angle {opts.slab_b_angle() * Constants::pi / 180.0};
		if (!(std::abs(opts.slab_b_angle()) <= 90.0))
		{
			Errors::fatal_error("slab_b_angle must be from -90 to 90 "
				"degrees, got " + std::to_string(opts.slab_b_angle()) + ".");
		}
		if (opts.slab_b_angle() != 0.0 && (bent
			|| opts.slab_efield_y() != 0.0
			|| opts.slab_efield_y_rate() != 0.0
			|| opts.slab_bgrad_y() != 0.0))
		{
			Errors::fatal_error("slab_b_angle is only for a straight slab "
				"(slab_major_radius = 0) without slab_efield_y, "
				"slab_efield_y_rate or slab_bgrad_y.");
		}
		const bool has_wall {opts.slab_wall() != "none"};
		if (has_wall && bent)
		{
			Errors::fatal_error("slab_wall is only for a straight slab "
				"(slab_major_radius = 0).");
		}
		if (has_wall && opts.use_gpu_int() > 0)
		{
			Errors::fatal_error("slab_wall only runs on the CPU for now. Set "
				"use_gpu = \"off\".");
		}

		if (std::abs(opts.slab_bgrad_y()) * half_width >= 1.0)
		{
			Errors::fatal_error("slab_bgrad_y is too large: the magnetic "
				"field would reach 0 inside the box (|slab_bgrad_y| must "
				"be under 40 1/m).");
		}

		// The bent slab needs more cells along the toroidal angle, since
		// the Cartesian components of the field turn with it
		constexpr int tdim {128};
		constexpr int xdim {8};
		constexpr int ydim {8};
		const int zdim {bent ? 32 : 8};

		// Cell edges, and the time frames. Straight: x, y, z are X, Y, Z
		// (m). Bent: x is R (m), y is Z (m) and z is phi (rad).
		auto grid_x = bent
			? Utilities::linspace(R0 - half_width, R0 + half_width, xdim + 1)
			: Utilities::linspace(0.0, 2.0 * half_width, xdim + 1);
		auto grid_y = Utilities::linspace(-half_width, half_width, ydim + 1);
		auto grid_z = bent
			? Utilities::linspace(0.0, Constants::pi / 2.0, zdim + 1)
			: Utilities::linspace(-half_width, half_width, zdim + 1);
		auto times = Utilities::linspace(0.0, opts.slab_tmax(), tdim);

		// Cell centers: evenly spaced, starting half a cell in
		const double xwidth {grid_x[1] - grid_x[0]};
		const double ywidth {grid_y[1] - grid_y[0]};
		const double zwidth {grid_z[1] - grid_z[0]};
		auto x = Utilities::linspace(grid_x.front() + xwidth / 2.0,
			grid_x.back() - xwidth / 2.0, xdim);
		auto y = Utilities::linspace(grid_y.front() + ywidth / 2.0,
			grid_y.back() - ywidth / 2.0, ydim);
		auto z = Utilities::linspace(grid_z.front() + zwidth / 2.0,
			grid_z.back() - zwidth / 2.0, zdim);

		// Geometry: the reciprocal basis vectors (gradients of x, y, z in
		// Cartesian space), the Jacobian and the Cartesian coordinates of
		// the cell centers.
		Vectors::Vector3D<BkgFPType> dxdX {xdim, ydim, zdim};
		Vectors::Vector3D<BkgFPType> dxdY {xdim, ydim, zdim};
		Vectors::Vector3D<BkgFPType> dxdZ {xdim, ydim, zdim};
		Vectors::Vector3D<BkgFPType> dydX {xdim, ydim, zdim};
		Vectors::Vector3D<BkgFPType> dydY {xdim, ydim, zdim};
		Vectors::Vector3D<BkgFPType> dydZ {xdim, ydim, zdim};
		Vectors::Vector3D<BkgFPType> dzdX {xdim, ydim, zdim};
		Vectors::Vector3D<BkgFPType> dzdY {xdim, ydim, zdim};
		Vectors::Vector3D<BkgFPType> dzdZ {xdim, ydim, zdim};
		Vectors::Vector3D<BkgFPType> J {xdim, ydim, zdim};
		Vectors::Vector3D<BkgFPType> X {xdim, ydim, zdim};
		Vectors::Vector3D<BkgFPType> Y {xdim, ydim, zdim};
		Vectors::Vector3D<BkgFPType> Z {xdim, ydim, zdim};
		for (int i {}; i < xdim; ++i)
		{
		for (int j {}; j < ydim; ++j)
		{
		for (int k {}; k < zdim; ++k)
		{
			const int idx {dxdX.calc_index(i, j, k)};
			if (bent)
			{
				// Cylindrical: x = R, y = Z, z = phi. grad R = (cos phi,
				// sin phi, 0), grad Z = (0, 0, 1) and grad phi = (-sin
				// phi, cos phi, 0) / R.
				const double R {x[i]};
				const double phi {z[k]};
				dxdX.get_data()[idx] = std::cos(phi);
				dxdY.get_data()[idx] = std::sin(phi);
				dydZ.get_data()[idx] = 1.0;
				dzdX.get_data()[idx] = -std::sin(phi) / R;
				dzdY.get_data()[idx] = std::cos(phi) / R;
				J.get_data()[idx] = R;
				X.get_data()[idx] = R * std::cos(phi);
				Y.get_data()[idx] = R * std::sin(phi);
				Z.get_data()[idx] = y[j];
			}
			else
			{
				// x = X, y = Y, z = Z
				dxdX.get_data()[idx] = 1.0;
				dydY.get_data()[idx] = 1.0;
				dzdZ.get_data()[idx] = 1.0;
				J.get_data()[idx] = 1.0;
				X.get_data()[idx] = x[i];
				Y.get_data()[idx] = y[j];
				Z.get_data()[idx] = z[k];
			}
		}
		}
		}

		// Makes a field from a function of the time frame and cell
		// indices
		auto make_field = [&](auto&& func)
		{
			Vectors::Vector4D<BkgFPType> v {tdim, xdim, ydim, zdim};
			for (int it {}; it < tdim; ++it)
			{
			for (int i {}; i < xdim; ++i)
			{
			for (int j {}; j < ydim; ++j)
			{
			for (int k {}; k < zdim; ++k)
			{
				v.get_data()[v.calc_index(it, i, j, k)]
					= static_cast<BkgFPType>(func(it, i, j, k));
			}
			}
			}
			}
			return v;
		};
		auto uniform = [&](const double val)
		{
			return make_field([val](int, int, int, int) {return val;});
		};

		// Field strength, and the Cartesian components of the unit vector
		// along the field. Straight: along +X turned toward -Y by
		// slab_b_angle, with the optional gradient in Y. Bent: toroidal,
		// (-sin phi, cos phi, 0), falling as 1 / R.
		auto b_strength = [&](int, int i, int j, int)
		{
			return bent ? opts.slab_bmag() * R0 / x[i]
				: opts.slab_bmag() * (1.0 + opts.slab_bgrad_y() * y[j]);
		};
		auto dir_X = [&](int, int, int, int k)
			{return bent ? -std::sin(z[k]) : std::cos(b_angle);};
		auto dir_Y = [&](int, int, int, int k)
			{return bent ? std::cos(z[k]) : -std::sin(b_angle);};

		// Plasma. The flow is along the field. The electric field is
		// along +Y and changes linearly in time (straight slab only).
		auto ne = uniform(opts.slab_ne());
		auto te = uniform(opts.slab_te());
		auto ti = uniform(opts.slab_ti());
		auto bX = make_field([&](int it, int i, int j, int k)
			{return b_strength(it, i, j, k) * dir_X(it, i, j, k);});
		auto bY = make_field([&](int it, int i, int j, int k)
			{return b_strength(it, i, j, k) * dir_Y(it, i, j, k);});
		auto bZ = uniform(0.0);
		auto bmag = make_field(b_strength);
		auto e_field = [&](int it, int, int, int)
		{
			return opts.slab_efield_y()
				+ opts.slab_efield_y_rate() * times[it];
		};
		auto eX = uniform(0.0);
		auto eY = make_field(e_field);
		auto eZ = uniform(0.0);
		auto emag = make_field([&](int it, int i, int j, int k)
			{return std::abs(e_field(it, i, j, k));});
		auto uX = make_field([&](int it, int i, int j, int k)
			{return opts.slab_flow() * dir_X(it, i, j, k);});
		auto uY = make_field([&](int it, int i, int j, int k)
			{return opts.slab_flow() * dir_Y(it, i, j, k);});
		auto uZ = uniform(0.0);

		// Put it all into a Background
		Background::Background bkg {};
		bkg.move_into_times(times);
		bkg.set_elec_mass_amu(opts.gkyl_elec_mass_amu());
		bkg.set_ion_mass_amu(opts.gkyl_ion_mass_amu());
		Zones::Zones zones {Background::single_zone(x, y, z, grid_x,
			grid_y, grid_z, opts)};
		if (!bent) zones.set_cartesian(0);
		if (has_wall)
		{
			// The wall face, and the one across from it
			const std::string& w {opts.slab_wall()};
			const Zones::Face face {w == "xmin" ? Zones::xmin
				: w == "xmax" ? Zones::xmax
				: w == "ymin" ? Zones::ymin : Zones::ymax};
			const auto across {static_cast<Zones::Face>(
				static_cast<int>(face) ^ 1)};
			if (zones.get_face(0, across) == 0)
			{
				Errors::fatal_error("With slab_wall = \"" + w + "\", the "
					"face across from the wall can't be periodic, or "
					"particles would wrap onto the wall. Make it absorbing.");
			}
			zones.set_face(0, face, Zones::face_wall);
		}
		bkg.move_into_zones(zones);
		bkg.move_into_ne(ne);
		bkg.move_into_te(te);
		bkg.move_into_ti(ti);
		bkg.move_into_bX(bX);
		bkg.move_into_bY(bY);
		bkg.move_into_bZ(bZ);
		bkg.move_into_bmag(bmag);
		bkg.move_into_eX(eX);
		bkg.move_into_eY(eY);
		bkg.move_into_eZ(eZ);
		bkg.move_into_emag(emag);
		bkg.move_into_uX(uX);
		bkg.move_into_uY(uY);
		bkg.move_into_uZ(uZ);

		// Parallel conductive heat fluxes (W/m^2), positive along B:
		// q_par = eps n T v_th, with eps the option and n T v_th that
		// species' thermal energy moving at its thermal speed v_th =
		// sqrt(T / m). Only made if one of them is nonzero, so that a slab
		// without heat fluxes looks like any background without them.
		if (opts.slab_heat_flux_i() != 0.0 || opts.slab_heat_flux_e() != 0.0)
		{
			const double n {opts.slab_ne()};
			const double Ti_J {opts.slab_ti() * Constants::ev_to_j};
			const double Te_J {opts.slab_te() * Constants::ev_to_j};
			const double mi {opts.gkyl_ion_mass_amu()
				* Constants::amu_to_kg};
			const double me {opts.gkyl_elec_mass_amu()
				* Constants::amu_to_kg};
			auto q_par_i = uniform(opts.slab_heat_flux_i() * n * Ti_J
				* std::sqrt(Ti_J / mi));
			auto q_par_e = uniform(opts.slab_heat_flux_e() * n * Te_J
				* std::sqrt(Te_J / me));
			bkg.move_into_q_par_i(q_par_i);
			bkg.move_into_q_par_e(q_par_e);
		}

		bkg.move_into_dxdX(dxdX);
		bkg.move_into_dxdY(dxdY);
		bkg.move_into_dxdZ(dxdZ);
		bkg.move_into_dydX(dydX);
		bkg.move_into_dydY(dydY);
		bkg.move_into_dydZ(dydZ);
		bkg.move_into_dzdX(dzdX);
		bkg.move_into_dzdY(dzdY);
		bkg.move_into_dzdZ(dzdZ);
		bkg.move_into_J(J);
		bkg.move_into_X(X);
		bkg.move_into_Y(Y);
		bkg.move_into_Z(Z);

		return bkg;
	}
}
