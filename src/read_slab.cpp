/**
* @file read_slab.cpp
*
* @brief Routines that make the slab background, a uniform plasma set by
*        input options
*
* The slab background (bkg_source = "slab") is a box of uniform plasma in
* a uniform magnetic field, with everything about the plasma given by the
* slab_* input options. It exists so the physics in Flan can be checked
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
	* Geometry: a 5 cm cube, x from 0 to 0.05 m and y and z from -0.025 to
	* 0.025 m, with x = X, y = Y and z = Z (Cartesian, so the Jacobian is
	* 1), in 8 x 8 x 8 cells. The plasma is uniform, so the grid only needs
	* to be fine enough to interpolate on. For particles to stay in the
	* box, make the boundaries periodic in the input file. There are 128
	* time frames from 0 to slab_tmax, and the plasma doesn't change in
	* time. The frames set the time resolution of the output.
	*
	* Plasma:
	* - Magnetic field slab_bmag along +X, no electric field.
	* - Density slab_ne (electrons and ions), temperatures slab_te and
	*   slab_ti.
	* - Flow slab_flow along +X, i.e., along the field.
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
	* field strength or slab_tmax isn't above 0.
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

		constexpr int tdim {128};
		constexpr int xdim {8};
		constexpr int ydim {8};
		constexpr int zdim {8};

		// Cell edges of the 5 cm cube, and the time frames
		auto grid_x = Utilities::linspace(0.00, 0.05, xdim + 1);
		auto grid_y = Utilities::linspace(-0.025, 0.025, ydim + 1);
		auto grid_z = Utilities::linspace(-0.025, 0.025, zdim + 1);
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

		// Geometry. x = X, y = Y, z = Z: the reciprocal basis vectors are
		// the unit vectors and the Jacobian is 1.
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

		// A field with the same value everywhere and at all times
		auto uniform = [](const double val)
		{
			Vectors::Vector4D<BkgFPType> v {tdim, xdim, ydim, zdim};
			std::fill(v.get_data().begin(), v.get_data().end(),
				static_cast<BkgFPType>(val));
			return v;
		};

		// Plasma: B and the flow along +X, no electric field
		auto ne = uniform(opts.slab_ne());
		auto te = uniform(opts.slab_te());
		auto ti = uniform(opts.slab_ti());
		auto bX = uniform(opts.slab_bmag());
		auto bY = uniform(0.0);
		auto bZ = uniform(0.0);
		auto bmag = uniform(opts.slab_bmag());
		auto eX = uniform(0.0);
		auto eY = uniform(0.0);
		auto eZ = uniform(0.0);
		auto emag = uniform(0.0);
		auto uX = uniform(opts.slab_flow());
		auto uY = uniform(0.0);
		auto uZ = uniform(0.0);

		// Put it all into a Background
		Background::Background bkg {};
		bkg.move_into_times(times);
		bkg.set_elec_mass_amu(opts.gkyl_elec_mass_amu());
		bkg.set_ion_mass_amu(opts.gkyl_ion_mass_amu());
		Zones::Zones zones {Background::single_zone(x, y, z, grid_x,
			grid_y, grid_z, opts)};
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

		// Parallel conductive heat fluxes (W/m^2), positive along B (+X):
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
