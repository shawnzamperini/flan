/**
* @file read_hermes.cpp
*
* @brief Reads a HERMES (JAX-Hermes3 or BOUT++ Hermes-3) simulation into a
* Background.
*
* A HERMES run is on a BOUT++ grid made by hypnotoad: a logically
* rectangular (x, y) grid in the poloidal plane, x across the flux surfaces
* and y along them, with z the toroidal angle (shifted metric, so z is the
* lab-frame angle). The grid is split at the separatrix (x = ixseps1) and,
* for a single null, at the ends of the divertor legs in y (jyseps1_1 and
* jyseps2_2), into 6 Flan zones:
*
*   1 inner PFR    x < ixseps1,   y <= jyseps1_1
*   2 core         x < ixseps1,   jyseps1_1 < y <= jyseps2_2
*   3 outer PFR    x < ixseps1,   y > jyseps2_2
*   4 inner leg    x >= ixseps1,  y <= jyseps1_1
*   5 main SOL     x >= ixseps1,  jyseps1_1 < y <= jyseps2_2
*   6 outer leg    x >= ixseps1,  y > jyseps2_2
*
* Each zone's x and y are the global grid cell indices (x counting the
* guard cells, so the physical cells are 2 to nx - 3 and the separatrix is
* at x = ixseps1), and z is the toroidal angle in radians, with cell k
* centered at k dz as in BOUT++. The y ends of the legs are the targets
* (wall faces), the core's inner x face is the core boundary and the PFR's
* inner and the SOL's outer x faces absorb particles. The core is periodic
* in y and the PFR connects its two legs.
*
* The grid is a shifted metric, so a field line moves toroidally by the
* change in zShift between y planes. Each zone gets a zShift table (see
* set_z_shift_tables), and plasma quantities are interpolated and
* differentiated in y along the field line.
*
* JAX-Hermes3 writes one uncompressed .npz file per time with the
* normalized Ne, Pe, Pi, NVi, phi and t on the (x, y, z) cell centers.
* BOUT++ Hermes-3 dump files are read by read_hermes_bout.cpp. Both go
* through Frames and finish_background here. The geometry, B and E come
* from the code shared with other multi-zone backgrounds (axisym_mesh.cpp).
*/

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <netcdf>
#include <string>
#include <vector>

#include "axisym_mesh.h"
#include "background.h"
#include "constants.h"
#include "errors.h"
#include "flan_types.h"
#include "npz.h"
#include "options.h"
#include "read_hermes.h"
#include "vectors.h"
#include "zones.h"

namespace Hermes
{
	namespace
	{
		// Reads a variable of a netCDF file as doubles, stopping the run if
		// it's missing or doesn't have the expected number of values
		std::vector<double> read_nc(const netCDF::NcFile& file,
			const std::string& name, const std::size_t want,
			const std::string& path)
		{
			const netCDF::NcVar var {file.getVar(name)};
			if (var.isNull())
			{
				Errors::fatal_error(path + " has no variable " + name + ".");
			}
			std::size_t n {1};
			for (const auto& d : var.getDims()) n *= d.getSize();
			if (n != want)
			{
				Errors::fatal_error(name + " in " + path + " has "
					+ std::to_string(n) + " values, but should have "
					+ std::to_string(want) + ".");
			}
			std::vector<double> out(n);
			var.getVar(out.data());
			return out;
		}

		// Reads an integer scalar of a netCDF file
		int read_nc_int(const netCDF::NcFile& file, const std::string& name,
			const std::string& path)
		{
			return static_cast<int>(std::lround(read_nc(file, name, 1,
				path)[0]));
		}
	}

	/**
	* Reads a single-null hypnotoad grid.
	*
	* The grid must be a single null (jyseps2_1 = jyseps1_2, the upper
	* legs empty), without y guard cells, and its separatrix and leg
	* indices must leave every zone at least one cell. Anything missing or
	* inconsistent stops the run.
	*
	* @param path Path to the grid file (netCDF).
	* @return The grid.
	*/
	Grid read_grid(const std::string& path)
	{
		if (!std::filesystem::exists(path))
			Errors::fatal_error("The HERMES grid file " + path
				+ " doesn't exist.");

		Grid g {};
		try
		{
			const netCDF::NcFile file(path, netCDF::NcFile::read);
			g.nx = read_nc_int(file, "nx", path);
			g.ny = read_nc_int(file, "ny", path);
			g.ixseps1 = read_nc_int(file, "ixseps1", path);
			g.jyseps1_1 = read_nc_int(file, "jyseps1_1", path);
			g.jyseps2_2 = read_nc_int(file, "jyseps2_2", path);
			const int jyseps2_1 {read_nc_int(file, "jyseps2_1", path)};
			const int jyseps1_2 {read_nc_int(file, "jyseps1_2", path)};
			if (jyseps2_1 != jyseps1_2)
			{
				Errors::fatal_error(path + " isn't a single null "
					"(jyseps2_1 = " + std::to_string(jyseps2_1)
					+ ", jyseps1_2 = " + std::to_string(jyseps1_2)
					+ "). Only single-null grids are handled.");
			}
			if (!file.getVar("y_boundary_guards").isNull()
				&& read_nc_int(file, "y_boundary_guards", path) != 0)
			{
				Errors::fatal_error(path + " has y guard cells, which "
					"aren't handled.");
			}

			const std::size_t n {static_cast<std::size_t>(g.nx) * g.ny};
			g.R_ll = read_nc(file, "Rxy_corners", n, path);
			g.R_lr = read_nc(file, "Rxy_lower_right_corners", n, path);
			g.R_ul = read_nc(file, "Rxy_upper_left_corners", n, path);
			g.R_ur = read_nc(file, "Rxy_upper_right_corners", n, path);
			g.Z_ll = read_nc(file, "Zxy_corners", n, path);
			g.Z_lr = read_nc(file, "Zxy_lower_right_corners", n, path);
			g.Z_ul = read_nc(file, "Zxy_upper_left_corners", n, path);
			g.Z_ur = read_nc(file, "Zxy_upper_right_corners", n, path);
			g.Br = read_nc(file, "Brxy", n, path);
			g.Bz = read_nc(file, "Bzxy", n, path);
			g.Bt = read_nc(file, "Btxy", n, path);
			g.zShift = read_nc(file, "zShift", n, path);
			g.ShiftAngle = read_nc(file, "ShiftAngle", g.nx, path);
		}
		catch (const netCDF::exceptions::NcException& e)
		{
			Errors::fatal_error("Couldn't read the HERMES grid file " + path
				+ ": " + e.what());
		}

		if (!(g.ixseps1 > ghost && g.ixseps1 < g.nx - ghost
			&& g.jyseps1_1 >= 0 && g.jyseps2_2 > g.jyseps1_1
			&& g.jyseps2_2 < g.ny - 1))
		{
			Errors::fatal_error(path + " has nx = " + std::to_string(g.nx)
				+ ", ny = " + std::to_string(g.ny) + ", ixseps1 = "
				+ std::to_string(g.ixseps1) + ", jyseps1_1 = "
				+ std::to_string(g.jyseps1_1) + " and jyseps2_2 = "
				+ std::to_string(g.jyseps2_2) + ", which don't make six "
				"zones with at least one cell each.");
		}
		return g;
	}

	/**
	* The six zones' cells as global grid indices (see the file comment).
	*
	* @param grid The grid.
	* @return x and y ranges of zones 0 to 5.
	*/
	std::array<ZoneRange, 6> zone_ranges(const Grid& grid)
	{
		const int xa {ghost};
		const int xs {grid.ixseps1};
		const int xb {grid.nx - ghost};
		const int ya {0};
		const int y1 {grid.jyseps1_1 + 1};
		const int y2 {grid.jyseps2_2 + 1};
		const int yb {grid.ny};
		return {{
			{xa, xs, ya, y1, "inner PFR"},
			{xa, xs, y1, y2, "core"},
			{xa, xs, y2, yb, "outer PFR"},
			{xs, xb, ya, y1, "inner leg"},
			{xs, xb, y1, y2, "main SOL"},
			{xs, xb, y2, yb, "outer leg"}}};
	}

	/**
	* Builds the six zones and their faces (see the file comment).
	*
	* The toroidal domain is 2 pi / zperiod wide with nz cells, centered at
	* k dz (dz = 2 pi / (zperiod nz)) as in BOUT++, and periodic.
	*
	* @param grid The grid.
	* @param nz Number of toroidal cells (from the snapshots).
	* @param zperiod The domain is 1 / zperiod of the torus.
	* @param aligned If true, give the zones their zShift tables, so the
	*        plasma is interpolated and differentiated along the field
	*        line (hermes_field_aligned = "on").
	* @return The zones, validated, with their zShift tables if aligned.
	*/
	Zones::Zones make_zones(const Grid& grid, const int nz,
		const int zperiod, const bool aligned)
	{
		const double dz {2.0 * Constants::pi / (zperiod * nz)};
		Zones::Zones zones(ghost);
		for (const auto& r : zone_ranges(grid))
		{
			auto [x, gx] = Axisym::uniform_grid(r.x1 - r.x0, r.x0, r.x1);
			auto [y, gy] = Axisym::uniform_grid(r.y1 - r.y0, r.y0, r.y1);
			auto [z, gz] = Axisym::uniform_grid(nz, -0.5 * dz,
				(nz - 0.5) * dz);
			const int zone {zones.add_zone(x, y, z, gx, gy, gz)};

			// z is the toroidal angle of a wedge 2 pi / zperiod wide, so
			// Cartesian vectors turn across its z faces
			zones.set_toroidal(zone);
		}

		using namespace Zones;
		constexpr int in_pfr {0}, core {1}, out_pfr {2}, in_leg {3},
			sol {4}, out_leg {5};
		auto set = [&](int zone, std::array<int, 4> f)
		{
			zones.set_face(zone, xmin, f[0]);
			zones.set_face(zone, xmax, f[1]);
			zones.set_face(zone, ymin, f[2]);
			zones.set_face(zone, ymax, f[3]);
			zones.set_face(zone, zmin, zone);
			zones.set_face(zone, zmax, zone);
		};
		set(in_pfr, {face_absorbing, in_leg, face_wall, out_pfr});
		set(core, {face_core, sol, core, core});
		set(out_pfr, {face_absorbing, out_leg, in_pfr, face_wall});
		set(in_leg, {in_pfr, face_absorbing, face_wall, sol});
		set(sol, {core, face_absorbing, in_leg, out_leg});
		set(out_leg, {out_pfr, face_absorbing, sol, face_wall});
		zones.validate();
		if (aligned) set_z_shift_tables(grid, zones);
		return zones;
	}

	/**
	* Gives every zone its zShift table (see Zones::set_z_shift), on its
	* cells and ghost cells.
	*
	* A cell's value is the grid's zShift at that cell. A ghost cell takes
	* the zShift of the cell its plasma values come from
	* (Zones::resolve_cell), plus the core's twist where the field line
	* crosses the branch cut at the X-point: zShift isn't continuous
	* around the core, but jumps by ShiftAngle(x) there. hypnotoad's zShift
	* is otherwise continuous along y, including from one leg of the PFR
	* to the other. A ghost source is across the cut when the source is in
	* the core rows (jyseps1_1 < y <= jyseps2_2) but the ghost's own row,
	* continued from its zone, is below them (then -ShiftAngle) or above
	* them (+ShiftAngle).
	*
	* A ghost at a corner of a zone can have several sources (at the
	* X-point), and gets their average, as its plasma values do. Sources
	* across the cut in the SOL, where ShiftAngle isn't defined (hypnotoad
	* writes NaN), only come from such corners, and are left out when
	* another source exists.
	*
	* @param grid The grid.
	* @param zones The zones from make_zones, with their faces set.
	*/
	void set_z_shift_tables(const Grid& grid, Zones::Zones& zones)
	{
		const auto ranges {zone_ranges(grid)};
		const int y1 {grid.jyseps1_1 + 1};
		const int y2 {grid.jyseps2_2 + 1};
		const int g {zones.get_ghost()};
		auto at = [&](int i, int j) {return static_cast<std::size_t>(i)
			* grid.ny + j;};

		// ShiftAngle of the last core column, for a corner whose only
		// source is across the cut in the SOL
		const double edge_twist {grid.ShiftAngle[grid.ixseps1 - 1]};
		for (int zone {}; zone < zones.get_nzones(); ++zone)
		{
			const auto& r {ranges[zone]};
			std::vector<double> table {};
			for (int i {-g}; i < zones.get_nx(zone) + g; ++i)
			for (int j {-g}; j < zones.get_ny(zone) + g; ++j)
			{
				const int row {r.y0 + j};
				double sum {}, sum_fallback {};
				int n {}, n_fallback {};
				for (const auto& s : zones.resolve_cell({zone, i, j, 0}))
				{
					const int gx {ranges[s.zone].x0 + s.i};
					const int gy {ranges[s.zone].y0 + s.j};
					const double zs {grid.zShift[at(gx, gy)]};
					const bool core_row {gy >= y1 && gy < y2};
					const int cut {core_row && row < y1 ? -1
						: core_row && row >= y2 ? 1 : 0};
					if (cut == 0)
					{
						sum += zs;
						++n;
						continue;
					}
					const double twist {grid.ShiftAngle[gx]};
					if (std::isfinite(twist))
					{
						sum += zs + cut * twist;
						++n;
					}
					else
					{
						sum_fallback += zs + cut * edge_twist;
						++n_fallback;
					}
				}
				table.push_back(n > 0 ? sum / n : sum_fallback / n_fallback);
			}
			zones.set_z_shift(zone, std::move(table));
		}
	}

	/**
	* The mesh for the shared geometry code: each zone's corners from the
	* grid's cell corners, B at the cell centers, and phi_sign = +1, since
	* (R, z, Z) is right handed on hypnotoad grids (checked by
	* python/check_hermes_conventions.py).
	*
	* Corner (i, j) of a zone is the low-x, low-y corner of its cell (i, j),
	* except on the zone's high-x and high-y edges, which come from the
	* other corners of the last cells.
	*
	* @param grid The grid.
	* @return The mesh.
	*/
	Axisym::Mesh make_mesh(const Grid& grid)
	{
		Axisym::Mesh mesh {};
		mesh.phi_sign = 1.0;
		const int ny {grid.ny};
		auto at = [ny](int i, int j) {return static_cast<std::size_t>(i)
			* ny + j;};
		for (const auto& r : zone_ranges(grid))
		{
			std::vector<double> R {}, Z {}, Br {}, Bz {}, Bphi {};
			for (int i {r.x0}; i <= r.x1; ++i)
			for (int j {r.y0}; j <= r.y1; ++j)
			{
				const bool hi_x {i == r.x1};
				const bool hi_y {j == r.y1};
				const std::size_t c {at(hi_x ? i - 1 : i, hi_y ? j - 1 : j)};
				const auto& gR {hi_x ? (hi_y ? grid.R_ur : grid.R_lr)
					: (hi_y ? grid.R_ul : grid.R_ll)};
				const auto& gZ {hi_x ? (hi_y ? grid.Z_ur : grid.Z_lr)
					: (hi_y ? grid.Z_ul : grid.Z_ll)};
				R.push_back(gR[c]);
				Z.push_back(gZ[c]);
			}
			for (int i {r.x0}; i < r.x1; ++i)
			for (int j {r.y0}; j < r.y1; ++j)
			{
				Br.push_back(grid.Br[at(i, j)]);
				Bz.push_back(grid.Bz[at(i, j)]);
				Bphi.push_back(grid.Bt[at(i, j)]);
			}
			mesh.R_corners.push_back(std::move(R));
			mesh.Z_corners.push_back(std::move(Z));
			mesh.Br.push_back(std::move(Br));
			mesh.Bz.push_back(std::move(Bz));
			mesh.Bphi.push_back(std::move(Bphi));
		}
		return mesh;
	}

	/**
	* The sign relating HERMES's parallel flow to the flow along B.
	*
	* HERMES's V is positive toward +y. In every physical cell, the
	* poloidal field (Br, Bz) is compared with the direction of +y (from
	* the cell's corners). hypnotoad makes y run along the poloidal field,
	* so they should agree everywhere; if they disagree everywhere the
	* sign is -1, and a mix stops the run.
	*
	* @param grid The grid.
	* @return +1 if +y is along B, -1 if against it.
	*/
	double parallel_sign(const Grid& grid)
	{
		int along {}, against {};
		for (int i {ghost}; i < grid.nx - ghost; ++i)
		for (int j {}; j < grid.ny; ++j)
		{
			const std::size_t c {static_cast<std::size_t>(i) * grid.ny + j};
			const double dRdy {0.5 * (grid.R_ul[c] + grid.R_ur[c]
				- grid.R_ll[c] - grid.R_lr[c])};
			const double dZdy {0.5 * (grid.Z_ul[c] + grid.Z_ur[c]
				- grid.Z_ll[c] - grid.Z_lr[c])};
			if (grid.Br[c] * dRdy + grid.Bz[c] * dZdy > 0.0) ++along;
			else ++against;
		}
		if (along > 0 && against > 0)
		{
			Errors::fatal_error("The poloidal field points along +y in "
				+ std::to_string(along) + " cells of the HERMES grid and "
				"against it in " + std::to_string(against) + ", so the "
				"direction of the parallel flow is ambiguous.");
		}
		return along > 0 ? 1.0 : -1.0;
	}

	/**
	* Picks which of n times to load: indices start to end inclusive in
	* steps of stride. A negative start or end counts from the end (-1 is
	* the last), so start = -50 and end = -1 picks the last 50. At least
	* two must be picked.
	*
	* @param n Number of times there are.
	* @param start Index of the first, from 0.
	* @param end Index of the last.
	* @param stride Step between them (>= 1).
	* @param what What the times are, for messages, e.g. "dir has 3
	*        snapshots (snap*.npz)".
	* @return The indices, increasing.
	*/
	std::vector<int> select_indices(const int n, int start, int end,
		const int stride, const std::string& what)
	{
		if (start < 0) start += n;
		if (end < 0) end += n;
		if (stride < 1)
		{
			Errors::fatal_error("hermes_snap_stride = "
				+ std::to_string(stride) + ", but must be at least 1.");
		}
		if (start < 0 || end >= n || start > end)
		{
			Errors::fatal_error(what + ", so hermes_snap_start and "
				"hermes_snap_end must be within 0 to " + std::to_string(n - 1)
				+ " with start <= end, but they're " + std::to_string(start)
				+ " and " + std::to_string(end) + ".");
		}
		std::vector<int> out {};
		for (int s {start}; s <= end; s += stride) out.push_back(s);
		if (out.size() < 2)
		{
			Errors::fatal_error("A HERMES background needs at least two "
				"snapshots, but hermes_snap_start = " + std::to_string(start)
				+ ", hermes_snap_end = " + std::to_string(end) + " and "
				"hermes_snap_stride = " + std::to_string(stride)
				+ " pick " + std::to_string(out.size()) + ".");
		}
		return out;
	}

	/**
	* Picks the JAX-Hermes3 snapshots to load: the snap*.npz files in dir,
	* sorted by name (JAX-Hermes3 names them by time, zero padded), picked
	* by index with select_indices.
	*
	* @param dir Snapshot directory.
	* @param start Index of the first snapshot, from 0.
	* @param end Index of the last snapshot.
	* @param stride Step between snapshots (>= 1).
	* @return Paths of the snapshots, in time order.
	*/
	std::vector<std::string> select_snapshots(const std::string& dir,
		int start, int end, const int stride)
	{
		if (!std::filesystem::is_directory(dir))
		{
			Errors::fatal_error("hermes_snap_dir = " + dir + " isn't a "
				"directory.");
		}
		std::vector<std::string> all {};
		for (const auto& f : std::filesystem::directory_iterator(dir))
		{
			const std::string name {f.path().filename().string()};
			if (name.starts_with("snap") && name.ends_with(".npz"))
				all.push_back(f.path().string());
		}
		std::sort(all.begin(), all.end());
		const int n {static_cast<int>(all.size())};
		std::vector<std::string> out {};
		for (const int i : select_indices(n, start, end, stride, dir
			+ " has " + std::to_string(n) + " snapshots (snap*.npz)"))
			out.push_back(all[i]);
		return out;
	}

	/**
	* The normalization of a JAX-Hermes3 run, which isn't in its snapshots:
	* hermes_nnorm, hermes_tnorm, hermes_bnorm, hermes_ion_mass_amu and
	* hermes_zperiod.
	*
	* @param opts Simulation options.
	* @return The normalization.
	*/
	Norm norm_from_options(const Options::Options& opts)
	{
		return {opts.hermes_nnorm(), opts.hermes_tnorm(), opts.hermes_bnorm(),
			opts.hermes_ion_mass_amu(), opts.hermes_zperiod()};
	}

	/**
	* Stops the run if the normalization can't be right: zperiod below 1,
	* or a normalization or ion mass that isn't positive.
	*
	* @param norm The normalization.
	*/
	void check_norm(const Norm& norm)
	{
		if (norm.zperiod < 1)
		{
			Errors::fatal_error("hermes_zperiod = "
				+ std::to_string(norm.zperiod) + ", but must be at least 1.");
		}
		for (const double v : {norm.nnorm, norm.tnorm, norm.bnorm,
			norm.ion_aa})
		{
			if (!(v > 0.0))
			{
				Errors::fatal_error("hermes_nnorm, hermes_tnorm, "
					"hermes_bnorm and hermes_ion_mass_amu must all be "
					"positive.");
			}
		}
	}

	/**
	* Sets up the zoned arrays for nframes times on the background's zones,
	* and the table from global (x, y) cells to zone cells.
	*
	* @param grid The grid.
	* @param zones The background's zones (from make_zones).
	* @param norm The run's normalization.
	* @param nframes Number of times that will be stored.
	*/
	Frames::Frames(const Grid& grid, const Zones::Zones& zones,
		const Norm& norm, const int nframes)
		: m_norm {norm}, m_nframes {nframes}, m_ny {grid.ny},
		m_cs0 {std::sqrt(qe * norm.tnorm / mp)},
		m_t_unit {mp / (qe * norm.bnorm)}, m_sign {parallel_sign(grid)},
		m_cells(static_cast<std::size_t>(grid.nx) * grid.ny),
		m_ne {Axisym::zoned(zones, nframes)},
		m_te {Axisym::zoned(zones, nframes)},
		m_ti {Axisym::zoned(zones, nframes)},
		m_vp {Axisym::zoned(zones, nframes)},
		m_upar {Axisym::zoned(zones, nframes)}, m_zones {&zones},
		m_times(nframes), m_labels(nframes)
	{
		const auto ranges {zone_ranges(grid)};
		for (int zone {}; zone < static_cast<int>(ranges.size()); ++zone)
		{
			const ZoneRange& r {ranges[zone]};
			for (int x {r.x0}; x < r.x1; ++x)
			for (int y {r.y0}; y < r.y1; ++y)
			{
				m_cells[static_cast<std::size_t>(x) * grid.ny + y]
					= {zone, x - r.x0, y - r.y0};
			}
		}
	}

	/**
	* Stores one cell of one time from the normalized plasma, as BOUT++
	* unnormalizes it: ne = Ne Nnorm, Te = Pe / Ne Tnorm, Ti = Pi / Ne
	* Tnorm (one ion species of charge 1, so ni = ne), the potential
	* phi Tnorm, and the parallel flow along B, V c_s0 times parallel_sign.
	* A cell with Ne <= 0 gets zero temperatures and flow (V isn't used),
	* and is counted.
	*
	* The reader works out V, the ion flow along +y in units of c_s0, from
	* what its format has (NVi / (AA Ne) for JAX-Hermes3).
	*
	* @param t Index of the time (0 to nframes - 1).
	* @param x, y Global cell indices (x counting the guard cells).
	* @param k Toroidal cell index.
	* @param Ne, Pe, Pi, V, phi The normalized plasma at the cell.
	*/
	void Frames::store(const int t, const int x, const int y, const int k,
		const double Ne, const double Pe, const double Pi, const double V,
		const double phi)
	{
		const Cell& c {m_cells[static_cast<std::size_t>(x) * m_ny + y]};
		if (c.zone < 0) return;
		const double tnorm {m_norm.tnorm};
		m_ne(c.zone, t, c.i, c.j, k) = Ne * m_norm.nnorm;
		m_vp(c.zone, t, c.i, c.j, k) = phi * tnorm;
		if (Ne > 0.0)
		{
			m_te(c.zone, t, c.i, c.j, k) = Pe / Ne * tnorm;
			m_ti(c.zone, t, c.i, c.j, k) = Pi / Ne * tnorm;
			m_upar(c.zone, t, c.i, c.j, k) = m_sign * V * m_cs0;
		}
		else ++m_empty_cells;
	}

	/**
	* Makes room for the inductive parallel electric field, -dApar/dt,
	* which the reader then fills with store_dapar_dt. finish_background
	* adds it to E_par.
	*/
	void Frames::enable_inductive()
	{
		m_epar = Axisym::zoned(*m_zones, m_nframes);
	}

	/**
	* Stores one cell of one time of the inductive parallel electric field
	* from dApar/dt, as BOUT++ unnormalizes it: Apar is b . A in units of
	* Bnorm rho_s0 and time in units of 1 / Omega_ci, so E = -dApar/dt
	* Bnorm c_s0 (V/m), along b. b is along +y in HERMES, so as for the
	* flow, E along B is that times parallel_sign.
	*
	* @param t Index of the time (0 to nframes - 1).
	* @param x, y Global cell indices (x counting the guard cells).
	* @param k Toroidal cell index.
	* @param dapar_dt dApar/dt, normalized.
	*/
	void Frames::store_dapar_dt(const int t, const int x, const int y,
		const int k, const double dapar_dt)
	{
		const Cell& c {m_cells[static_cast<std::size_t>(x) * m_ny + y]};
		if (c.zone < 0) return;
		m_epar(c.zone, t, c.i, c.j, k) = static_cast<BkgFPType>(-m_sign
			* dapar_dt * m_norm.bnorm * m_cs0);
	}

	/**
	* Sets the time of one of the stored times.
	*
	* @param t Index of the time (0 to nframes - 1).
	* @param t_norm The time, normalized (units of m_p / (e Bnorm)).
	* @param label What the time came from (e.g. a snapshot path), for
	*        messages.
	*/
	void Frames::set_time(const int t, const double t_norm,
		const std::string& label)
	{
		m_times[t] = t_norm * m_t_unit;
		m_labels[t] = label;
	}

	/**
	* Moves the stored plasma into the background and calculates what
	* follows from it, the same for every format.
	*
	* The background must already have its zones, geometry and magnetic
	* field. The times must increase. Ghost cells are filled with
	* Zones::fill_ghosts, and E and the flows come from
	* Axisym::calc_efield_and_flows, with the inductive E_par added if the
	* reader stored it (hermes_inductive_epar = "on"). The runs' outputs have no heat
	* fluxes, so with hermes_heat_flux = "conduction" they're calculated
	* from the conduction (Axisym::calc_conductive_heat_fluxes), with the
	* hermes_kappa_* options and one ion species of charge 1.
	*
	* @param frames The stored plasma, moved out of (left empty).
	* @param opts Simulation options (hermes_heat_flux, hermes_kappa_*).
	* @param bkg Background to move the plasma into.
	*/
	void finish_background(Frames& frames, const Options::Options& opts,
		Background::Background& bkg)
	{
		const Zones::Zones& zones {bkg.get_zones()};
		const int nframes {frames.m_nframes};
		const auto& times {frames.m_times};
		std::cout << "Read " << nframes << " HERMES snapshots, t = "
			<< times.front() << " to " << times.back() << " s\n";
		if (frames.m_empty_cells > 0)
		{
			std::cout << "  " << frames.m_empty_cells << " cells (over all "
				"snapshots) have Ne <= 0, and are given zero temperatures "
				"and flow.\n";
		}

		for (int t {1}; t < nframes; ++t)
		{
			if (!(times[t] > times[t - 1]))
			{
				Errors::fatal_error("The HERMES snapshots' times don't "
					"increase: " + frames.m_labels[t] + " is at "
					+ std::to_string(times[t]) + " s, but the snapshot before "
					"it is at " + std::to_string(times[t - 1]) + " s.");
			}
		}

		for (auto* v : {&frames.m_ne, &frames.m_te, &frames.m_ti,
			&frames.m_vp, &frames.m_upar}) zones.fill_ghosts(*v);
		bkg.move_into_times(frames.m_times);
		bkg.move_into_ne(frames.m_ne);
		bkg.move_into_te(frames.m_te);
		bkg.move_into_ti(frames.m_ti);
		const bool inductive {frames.inductive()};
		if (inductive)
		{
			zones.fill_ghosts(frames.m_epar);
			double sum2 {};
			std::size_t n {};
			for (const double e : frames.m_epar.get_data())
			{
				sum2 += e * e;
				++n;
			}
			std::cout << "  Inductive E_par = -dApar/dt added to E (rms "
				<< std::sqrt(sum2 / std::max<std::size_t>(n, 1))
				<< " V/m)\n";
		}
		Axisym::calc_efield_and_flows(bkg, frames.m_vp, frames.m_upar,
			inductive ? &frames.m_epar : nullptr);
		frames.m_epar = {};
		bkg.move_into_vp(frames.m_vp);

		// The runs don't write the heat fluxes, so they're calculated from
		// the conduction, as the run calculated them
		if (opts.hermes_heat_flux_int() == 1)
		{
			Axisym::Conduction c {};
			c.coef_e = opts.hermes_kappa_coef_e();
			c.coef_i = opts.hermes_kappa_coef_i();
			c.limit_alpha = opts.hermes_kappa_limit_alpha();
			c.ion_z = 1.0;
			c.ion_aa = frames.m_norm.ion_aa;
			Axisym::calc_conductive_heat_fluxes(bkg, c);
			std::cout << "  Parallel heat fluxes from Braginskii conduction: "
				"kappa coefficients " << c.coef_e << " (electrons), "
				<< c.coef_i << " (ions), ";
			if (c.limit_alpha > 0.0)
				std::cout << "flux limiter " << c.limit_alpha << '\n';
			else std::cout << "no flux limiter\n";
		}
		else
		{
			std::cout << "  hermes_heat_flux = \"off\": no parallel heat "
				"fluxes, so no thermal force.\n";
		}
	}

	/**
	* Reads JAX-Hermes3 snapshots into frames: each snapshot's Ne, Pe, Pi,
	* NVi, phi and t, with the ion flow V = NVi / (AA Ne).
	*
	* Every snapshot must have the grid's nx and ny and the same nz as the
	* frames' zones.
	*
	* @param paths Snapshot paths, in time order (frames.nframes() of them).
	* @param grid The grid.
	* @param frames Where the plasma goes.
	*/
	void read_snapshots(const std::vector<std::string>& paths,
		const Grid& grid, Frames& frames)
	{
		const int nframes {static_cast<int>(paths.size())};
		if (nframes != frames.nframes())
		{
			Errors::fatal_error("read_snapshots got "
				+ std::to_string(nframes) + " snapshots for "
				+ std::to_string(frames.nframes()) + " frames.");
		}
		const double aa {frames.norm().ion_aa};
		const std::vector<std::size_t> want_shape {
			static_cast<std::size_t>(grid.nx),
			static_cast<std::size_t>(grid.ny), 0};
		std::size_t nz {};

		for (int t {}; t < nframes; ++t)
		{
			const std::string& path {paths[t]};
			const Npz::File file(path);
			auto read = [&](const std::string& name)
			{
				auto a {file.read(name)};
				if (nz == 0 && a.shape.size() == 3) nz = a.shape[2];
				if (a.shape.size() != 3 || a.shape[0] != want_shape[0]
					|| a.shape[1] != want_shape[1] || a.shape[2] != nz)
				{
					Errors::fatal_error(name + " in " + path + " isn't ("
						+ std::to_string(grid.nx) + ", "
						+ std::to_string(grid.ny) + ", " + std::to_string(nz)
						+ "), the grid's nx and ny and the first snapshot's "
						"nz.");
				}
				return std::move(a.data);
			};
			const auto Ne {read("Ne")};
			const auto Pe {read("Pe")};
			const auto Pi {read("Pi")};
			const auto NVi {read("NVi")};
			const auto phi {read("phi")};
			frames.set_time(t, file.read("t").data.at(0), path);

			std::size_t c {};
			for (int x {}; x < grid.nx; ++x)
			for (int y {}; y < grid.ny; ++y)
			for (int k {}; k < static_cast<int>(nz); ++k, ++c)
			{
				const double n {Ne[c]};
				frames.store(t, x, y, k, n, Pe[c], Pi[c],
					n > 0.0 ? NVi[c] / (aa * n) : 0.0, phi[c]);
			}
		}
	}

	/**
	* Reads a HERMES simulation into a Background.
	*
	* The grid is hermes_grid_file. For JAX-Hermes3 (hermes_format = "jax")
	* the snapshots are picked from hermes_snap_dir with select_snapshots,
	* the normalization comes from the options and the number of toroidal
	* cells from the first snapshot. For BOUT++ Hermes-3 ("bout")
	* hermes_snap_dir is the run directory, the output times are picked
	* with select_indices, and the normalization, zperiod and nz come from
	* the dumps (read_bout_layout, norm_from_bout).
	*
	* @param opts Simulation options.
	* @return The background.
	*/
	Background::Background read_hermes(const Options::Options& opts)
	{
		if (opts.hermes_grid_file() == "undefined"
			|| opts.hermes_snap_dir() == "undefined")
		{
			Errors::fatal_error("hermes_grid_file and hermes_snap_dir must "
				"be set for a HERMES background.");
		}
		const bool bout {opts.hermes_format_int() == 1};

		// What to read, the normalization and nz, per format
		std::vector<std::string> paths {};
		std::vector<int> times {};
		BoutLayout layout {};
		Norm norm {};
		int nz {};
		if (bout)
		{
			layout = read_bout_layout(opts.hermes_snap_dir());
			norm = norm_from_bout(layout, opts);
			check_norm(norm);
			times = select_indices(static_cast<int>(layout.times.size()),
				opts.hermes_snap_start(), opts.hermes_snap_end(),
				opts.hermes_snap_stride(), opts.hermes_snap_dir() + " has "
				+ std::to_string(layout.times.size()) + " output times "
				"(t_array)");
			nz = layout.nz;
		}
		else
		{
			norm = norm_from_options(opts);
			check_norm(norm);
			paths = select_snapshots(opts.hermes_snap_dir(),
				opts.hermes_snap_start(), opts.hermes_snap_end(),
				opts.hermes_snap_stride());
		}
		const int nframes {static_cast<int>(bout ? times.size()
			: paths.size())};
		const int zperiod {norm.zperiod};
		const Grid grid {read_grid(opts.hermes_grid_file())};
		if (bout) check_bout_layout(layout, grid);
		else
		{
			// nz from the first snapshot's density
			const auto first {Npz::File(paths[0]).read("Ne")};
			if (first.shape.size() != 3)
			{
				Errors::fatal_error("Ne in " + paths[0]
					+ " isn't 3D (x, y, z).");
			}
			nz = static_cast<int>(first.shape[2]);
		}

		std::cout << "HERMES background: grid " << opts.hermes_grid_file()
			<< " (" << grid.nx << " x " << grid.ny << ", ixseps1 "
			<< grid.ixseps1 << ", jyseps1_1 " << grid.jyseps1_1
			<< ", jyseps2_2 " << grid.jyseps2_2 << "), nz = " << nz
			<< ", 1/" << zperiod << " of the torus\n"
			<< "  Nnorm = " << norm.nnorm << " m-3, Tnorm = "
			<< norm.tnorm << " eV, Bnorm = " << norm.bnorm
			<< " T, ion AA = " << norm.ion_aa
			<< (bout ? " (normalization from the dumps)" : "") << '\n';
		if (bout)
		{
			const double t_unit {mp / (qe * norm.bnorm)};
			const double gb {5.0 * nframes * grid.nx * grid.ny * nz
				* sizeof(BkgFPType) / 1e9};
			std::cout << "  BOUT++ dumps in " << layout.dir << ": "
				<< layout.paths.size() << " files (" << layout.nxpe << " x "
				<< layout.nype << " processes), " << layout.times.size()
				<< " output times; reading " << nframes << " of them, "
				<< times.front() << " to " << times.back() << " (t = "
				<< layout.times[times.front()] * t_unit << " to "
				<< layout.times[times.back()] * t_unit << " s), "
				<< gb << " GB of plasma before derived fields\n";
		}

		Background::Background bkg {};
		const bool aligned {opts.hermes_field_aligned_int() == 1};
		if (!aligned)
		{
			std::cout << "  hermes_field_aligned = \"off\": the plasma is "
				"interpolated and differentiated at fixed z, not along the "
				"field line. This is only meant for comparison.\n";
		}
		auto zones {make_zones(grid, nz, zperiod, aligned)};
		bkg.move_into_zones(zones);
		const Axisym::Mesh mesh {make_mesh(grid)};
		Axisym::calc_geometry(mesh, bkg);
		Axisym::calc_bfield(mesh, nframes, bkg);
		bkg.set_elec_mass_amu(Constants::mass_e / Constants::amu_to_kg);
		bkg.set_ion_mass_amu(norm.ion_aa);

		Frames frames(grid, bkg.get_zones(), norm, nframes);
		// The inductive E_par, -dApar/dt, needs Apar, which only BOUT++
		// dumps of an electromagnetic run have
		if (opts.hermes_inductive_epar_int() == 1)
		{
			if (bout && layout.has("Apar")) frames.enable_inductive();
			else
			{
				std::cout << "  No Apar in the " << (bout ? "dumps "
					"(electrostatic run)" : "JAX-Hermes3 snapshots")
					<< ", so E_par is -grad_par phi only.\n";
			}
		}
		if (bout)
		{
			read_bout_dumps(layout, times, opts.hermes_ion_species(),
				frames);
		}
		else read_snapshots(paths, grid, frames);
		finish_background(frames, opts, bkg);
		return bkg;
	}
}
