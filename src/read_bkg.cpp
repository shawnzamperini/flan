/**
* @file  read_bkg.cpp
*
* @brief Entry point for loading the plasma background, which dispatches to
*        the reader for the requested background source.
*/
#include <string>
#include <utility>
#include <vector>

#include "background.h"
#include "flan_types.h"
#include "errors.h"
#include "options.h"
#include "read_gkyl.h"
#include "read_test.h"
#include "zones.h"

namespace Background
{
	/**
	* Reads in the plasma background from the source chosen by the
	* bkg_source option.
	*
	* "test" loads a hardcoded test background with Test::read_test, and
	* "gkeyll" loads a Gkeyll simulation with Gkyl::read_gkyl. SOLEDGE3X
	* backgrounds (control integer 2) are not implemented yet. If the source
	* is not implemented or not recognized, the run is stopped with
	* Errors::fatal_error. Options::set_bkg_source already rejects any other
	* value, so this is only a safeguard.
	*
	* @param opts Simulation options.
	* @return Background object holding the plasma background.
	*/
	Background read_bkg(const Options::Options& opts)
	{
		// Background object to fill in
		Background bkg {};

		// Load in hardcoded plasma background for testing
		if (opts.bkg_source_int() == 0)
		{
			bkg = Test::read_test(opts);
		}

		// Load in a Gkeyll background
		else if (opts.bkg_source_int() == 1)
		{
			bkg = Gkyl::read_gkyl(opts);
		}

		// Load in a SOLEDGE3X background
		else if (opts.bkg_source_int() == 2)
		{
			Errors::fatal_error("SOLEDGE3X backgrounds are not implemented "
				"yet.");
		}

		// Unrecognized option
		else
		{
			Errors::fatal_error("bkg_source = " + opts.bkg_source() 
				+ " not recognized.");
		}

		return bkg;
	}  // read_bkg

	/**
	* Builds the zones for a background made of a single zone, such as a
	* Gkeyll or test background.
	*
	* The grid vectors are moved in, so they are left empty. Each face of
	* the zone comes from the matching *bound_type option with
	* Zones::face_from_bound_type, e.g., a periodic min_ybound_type connects
	* the ymin face to the zone itself.
	*
	* @param x x cell centers.
	* @param y y cell centers.
	* @param z z cell centers.
	* @param grid_x x cell edges.
	* @param grid_y y cell edges.
	* @param grid_z z cell edges.
	* @param opts Simulation options.
	* @return Zones holding the one zone.
	*/
	Zones::Zones single_zone(std::vector<BkgFPType>& x,
		std::vector<BkgFPType>& y, std::vector<BkgFPType>& z,
		std::vector<BkgFPType>& grid_x, std::vector<BkgFPType>& grid_y,
		std::vector<BkgFPType>& grid_z, const Options::Options& opts)
	{
		Zones::Zones zones {};
		const int zone {zones.add_zone(std::move(x), std::move(y),
			std::move(z), std::move(grid_x), std::move(grid_y),
			std::move(grid_z))};

		zones.set_face(zone, Zones::xmin, Zones::face_from_bound_type(
			opts.min_xbound_type_int(), zone));
		zones.set_face(zone, Zones::xmax, Zones::face_from_bound_type(
			opts.max_xbound_type_int(), zone));
		zones.set_face(zone, Zones::ymin, Zones::face_from_bound_type(
			opts.min_ybound_type_int(), zone));
		zones.set_face(zone, Zones::ymax, Zones::face_from_bound_type(
			opts.max_ybound_type_int(), zone));
		zones.set_face(zone, Zones::zmin, Zones::face_from_bound_type(
			opts.min_zbound_type_int(), zone));
		zones.set_face(zone, Zones::zmax, Zones::face_from_bound_type(
			opts.max_zbound_type_int(), zone));

		return zones;
	}
}  // namespace Background
