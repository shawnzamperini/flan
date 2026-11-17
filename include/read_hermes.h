/**
* @file read_hermes.h
*
* @brief Header file for read_hermes.cpp
*/

#pragma once

#include <array>
#include <string>
#include <vector>

#include "axisym_mesh.h"
#include "background.h"
#include "options.h"
#include "zones.h"

namespace Hermes
{
	/// Ghost cells on each side of every zone. The same as BOUT++'s x
	/// guard cells (MXG), which the grid and snapshots include.
	inline constexpr int ghost {2};

	/// Constants BOUT++ normalizes with (SI::qe, SI::Mp), so unnormalizing
	/// gives back exactly what the run used
	inline constexpr double qe {1.602176634e-19};
	inline constexpr double mp {1.67262192369e-27};

	/**
	* The parts of a single-null hypnotoad grid Flan uses. 2D arrays are
	* (nx, ny) in C order, including the x guard cells; there are no y
	* guard cells. Corners follow hypnotoad's names: "corners" is the
	* cell's low-x, low-y corner, and in the others "right" is +x and
	* "upper" is +y.
	*/
	struct Grid
	{
		int nx {};			///< Cells in x, guard cells included
		int ny {};			///< Cells in y
		int ixseps1 {};		///< First x index outside the separatrix
		int jyseps1_1 {};	///< Last y index of the inner leg
		int jyseps2_2 {};	///< Last y index of the core
		std::vector<double> R_ll {}, R_lr {}, R_ul {}, R_ur {};	///< (m)
		std::vector<double> Z_ll {}, Z_lr {}, Z_ul {}, Z_ur {};	///< (m)
		std::vector<double> Br {}, Bz {}, Bt {};	///< At centers (T)
		std::vector<double> zShift {};	///< Toroidal shift (rad)
		std::vector<double> ShiftAngle {};	///< Twist over the core (nx)
	};

	/// A zone's cells as global grid indices: x0 <= x < x1, y0 <= y < y1
	struct ZoneRange
	{
		int x0 {};
		int x1 {};
		int y0 {};
		int y1 {};
		std::string name {};
	};

	Grid read_grid(const std::string& path);

	std::array<ZoneRange, 6> zone_ranges(const Grid& grid);

	Zones::Zones make_zones(const Grid& grid, const int nz,
		const int zperiod);

	Axisym::Mesh make_mesh(const Grid& grid);

	double parallel_sign(const Grid& grid);

	std::vector<std::string> select_snapshots(const std::string& dir,
		int start, int end, const int stride);

	void read_snapshots(const std::vector<std::string>& paths,
		const Grid& grid, const Options::Options& opts,
		Background::Background& bkg);

	Background::Background read_hermes(const Options::Options& opts);
}
