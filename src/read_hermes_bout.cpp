/**
* @file read_hermes_bout.cpp
*
* @brief Reads the plasma of a BOUT++ Hermes-3 run (hermes_format = "bout").
*
* A BOUT++ run writes one netCDF file per MPI process, BOUT.dmp.N.nc for
* N = 0 to NXPE * NYPE - 1. Each holds its process's block of MXSUB x MYSUB
* cells in (x, y), with MXG x and MYG y guard cells on each side, and all
* MZ toroidal cells, at every output time: variables are (t, x, y, z). The
* block's place in the global grid is PE_XIND and PE_YIND. The global grid
* is NXPE MXSUB + 2 MXG by NYPE MYSUB, the grid file's nx (guards
* included) and ny (no y guards).
*
* A global field is assembled from each block's interior cells, plus its x
* guard cells where the block is at the domain's x boundary (they are the
* grid's x guard cells, which the zones include). y guard cells are always
* dropped: the grid has none, and at the targets they hold sheath boundary
* values, not plasma.
*
* The normalization (Nnorm, Tnorm, Bnorm), zperiod and the output times
* (t_array, in units of 1 / Omega_ci) are in every file. The ion mass isn't,
* so it still comes from hermes_ion_mass_amu.
*
* Variables are named after the species in the Hermes-3 input: Ne, Pe for
* the electrons ("e") and P<ion>, V<ion>, NV<ion> for the ions
* (hermes_ion_species). In an electromagnetic run Hermes-3 evolves the
* canonical momentum, so the dumped NV<ion> = AA N V + Z N Apar. The ion
* flow is therefore taken from V<ion> if it was dumped, and otherwise from
* (NV<ion> - N Apar) / (AA N) (or NV<ion> / (AA N) without Apar).
*
* Each file is opened once and every selected time of each variable is read
* from it as one strided hyperslab: opening a file costs much more than
* reading from it, so looping over the times outside the files would be
* far slower.
*/

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <iostream>
#include <netcdf>
#include <string>
#include <vector>

#include "errors.h"
#include "options.h"
#include "read_hermes.h"

namespace Hermes
{
	namespace
	{
		/**
		* Reads a scalar variable of a netCDF file as a double, stopping the
		* run if it's missing.
		*
		* @param file Open netCDF file.
		* @param name Variable name.
		* @param path The file's path, for messages.
		* @return The value.
		*/
		double read_scalar(const netCDF::NcFile& file,
			const std::string& name, const std::string& path)
		{
			const netCDF::NcVar var {file.getVar(name)};
			if (var.isNull())
				Errors::fatal_error(path + " has no variable " + name + ".");
			if (var.getDimCount() != 0)
				Errors::fatal_error(name + " in " + path + " isn't a scalar.");
			double v {};
			var.getVar(&v);
			return v;
		}

		/// read_scalar, rounded to an int
		int read_int(const netCDF::NcFile& file, const std::string& name,
			const std::string& path)
		{
			return static_cast<int>(std::lround(read_scalar(file, name,
				path)));
		}

		/**
		* Reads a file's output times, t_array.
		*
		* @param file Open netCDF file.
		* @param path The file's path, for messages.
		* @return The times (units of 1 / Omega_ci).
		*/
		std::vector<double> read_times(const netCDF::NcFile& file,
			const std::string& path)
		{
			const netCDF::NcVar var {file.getVar("t_array")};
			if (var.isNull() || var.getDimCount() != 1)
			{
				Errors::fatal_error(path + " has no 1D t_array (the output "
					"times).");
			}
			std::vector<double> t(var.getDim(0).getSize());
			if (!t.empty()) var.getVar(t.data());
			return t;
		}

		/// How the ion flow is found from the dumps
		enum class FlowFrom {V, NV, NV_Apar};
	}

	/**
	* @param name Variable name.
	* @return Whether the dumps have the variable.
	*/
	bool BoutLayout::has(const std::string& name) const
	{
		return std::find(vars.begin(), vars.end(), name) != vars.end();
	}

	/**
	* Reads the layout, normalization and output times of a BOUT++ run from
	* its first dump file, and checks that all the dump files are there.
	*
	* Stops the run if dir has no BOUT.dmp.0.nc, a dump file is missing,
	* z is split over processes (NZPE > 1, not handled), or the dumps'
	* global size doesn't follow from the process layout.
	*
	* @param dir The run directory (hermes_snap_dir).
	* @return The layout.
	*/
	BoutLayout read_bout_layout(const std::string& dir)
	{
		if (!std::filesystem::is_directory(dir))
		{
			Errors::fatal_error("hermes_snap_dir = " + dir + " isn't a "
				"directory.");
		}
		const std::string first {(std::filesystem::path(dir)
			/ "BOUT.dmp.0.nc").string()};
		if (!std::filesystem::exists(first))
		{
			Errors::fatal_error("hermes_snap_dir = " + dir + " has no "
				"BOUT.dmp.0.nc. With hermes_format = \"bout\" it must be the "
				"directory of a BOUT++ Hermes-3 run's dump files.");
		}

		BoutLayout l {};
		l.dir = dir;
		try
		{
			const netCDF::NcFile file(first, netCDF::NcFile::read);
			l.nxpe = read_int(file, "NXPE", first);
			l.nype = read_int(file, "NYPE", first);
			if (!file.getVar("NZPE").isNull()
				&& read_int(file, "NZPE", first) != 1)
			{
				Errors::fatal_error("The BOUT++ run in " + dir + " splits z "
					"over processes (NZPE = "
					+ std::to_string(read_int(file, "NZPE", first))
					+ "), which isn't handled.");
			}
			l.mxsub = read_int(file, "MXSUB", first);
			l.mysub = read_int(file, "MYSUB", first);
			l.mxg = read_int(file, "MXG", first);
			l.myg = read_int(file, "MYG", first);
			l.nz = read_int(file, "MZ", first);
			l.zperiod = read_int(file, "zperiod", first);
			l.nnorm = read_scalar(file, "Nnorm", first);
			l.tnorm = read_scalar(file, "Tnorm", first);
			l.bnorm = read_scalar(file, "Bnorm", first);
			l.ixseps1 = read_int(file, "ixseps1", first);
			l.jyseps1_1 = read_int(file, "jyseps1_1", first);
			l.jyseps2_2 = read_int(file, "jyseps2_2", first);
			l.times = read_times(file, first);
			for (const auto& [name, var] : file.getVars())
				l.vars.push_back(name);

			l.nx = l.nxpe * l.mxsub + 2 * l.mxg;
			l.ny = l.nype * l.mysub;
			const int nx {read_int(file, "nx", first)};
			const int ny {read_int(file, "ny", first)};
			if (nx != l.nx || ny != l.ny)
			{
				Errors::fatal_error("The dumps in " + dir + " say nx = "
					+ std::to_string(nx) + ", ny = " + std::to_string(ny)
					+ ", but their process layout gives "
					+ std::to_string(l.nx) + " x " + std::to_string(l.ny)
					+ " (NXPE MXSUB + 2 MXG by NYPE MYSUB).");
			}
		}
		catch (const netCDF::exceptions::NcException& e)
		{
			Errors::fatal_error("Couldn't read the BOUT++ dump file " + first
				+ ": " + e.what());
		}

		const int nfiles {l.nxpe * l.nype};
		int missing {};
		std::string first_missing {};
		for (int n {}; n < nfiles; ++n)
		{
			const std::string p {(std::filesystem::path(dir)
				/ ("BOUT.dmp." + std::to_string(n) + ".nc")).string()};
			if (!std::filesystem::exists(p))
			{
				if (missing++ == 0) first_missing = p;
			}
			l.paths.push_back(p);
		}
		if (missing > 0)
		{
			Errors::fatal_error(std::to_string(missing) + " of the "
				+ std::to_string(nfiles) + " dump files (NXPE x NYPE) are "
				"missing, e.g. " + first_missing + ".");
		}
		return l;
	}

	/**
	* Stops the run if the dumps aren't of the grid: a different size or
	* different separatrix and X-point indices.
	*
	* @param layout The dumps' layout.
	* @param grid The grid (hermes_grid_file).
	*/
	void check_bout_layout(const BoutLayout& layout, const Grid& grid)
	{
		std::string bad {};
		auto cmp = [&](const std::string& name, const int dumps,
			const int g)
		{
			if (dumps != g)
			{
				bad += (bad.empty() ? "" : ", ") + name + " "
					+ std::to_string(dumps) + " vs " + std::to_string(g);
			}
		};
		cmp("nx", layout.nx, grid.nx);
		cmp("ny", layout.ny, grid.ny);
		cmp("ixseps1", layout.ixseps1, grid.ixseps1);
		cmp("jyseps1_1", layout.jyseps1_1, grid.jyseps1_1);
		cmp("jyseps2_2", layout.jyseps2_2, grid.jyseps2_2);
		if (!bad.empty())
		{
			Errors::fatal_error("The dumps in " + layout.dir + " aren't on "
				"the grid hermes_grid_file (dumps vs grid): " + bad + ".");
		}
	}

	/**
	* The normalization of a BOUT++ run: Nnorm, Tnorm, Bnorm and zperiod
	* from the dumps, and the ion mass from hermes_ion_mass_amu (it isn't
	* in the dumps).
	*
	* hermes_nnorm, hermes_tnorm, hermes_bnorm and hermes_zperiod aren't
	* used. If one was changed from its default to something other than the
	* dumps' value, a note says so, since the input file is then probably
	* meant for another run.
	*
	* @param layout The dumps' layout.
	* @param opts Simulation options.
	* @return The normalization.
	*/
	Norm norm_from_bout(const BoutLayout& layout,
		const Options::Options& opts)
	{
		const Options::Options defaults {};
		auto note = [](const std::string& opt, const double value,
			const double def, const std::string& dumped, const double used)
		{
			const bool changed {value != def};
			const bool differs {std::abs(value - used)
				> 1e-12 * std::abs(used)};
			if (changed && differs)
			{
				std::cout << "  Note: " << opt << " = " << value << " isn't "
					"used with hermes_format = \"bout\"; the dumps' " << dumped
					<< " = " << used << " is.\n";
			}
		};
		note("hermes_nnorm", opts.hermes_nnorm(), defaults.hermes_nnorm(),
			"Nnorm", layout.nnorm);
		note("hermes_tnorm", opts.hermes_tnorm(), defaults.hermes_tnorm(),
			"Tnorm", layout.tnorm);
		note("hermes_bnorm", opts.hermes_bnorm(), defaults.hermes_bnorm(),
			"Bnorm", layout.bnorm);
		note("hermes_zperiod", opts.hermes_zperiod(),
			defaults.hermes_zperiod(), "zperiod", layout.zperiod);
		return {layout.nnorm, layout.tnorm, layout.bnorm,
			opts.hermes_ion_mass_amu(), layout.zperiod};
	}

	/**
	* Reads the selected output times of a BOUT++ run into frames.
	*
	* Reads Ne, Pe, P<ion>, phi and the ion flow (see the file comment) of
	* every selected time from each dump file, opening each file once, and
	* stores them cell by cell. Stops the run if a variable is missing or
	* has the wrong shape, a file's times or place differ from what the
	* layout says, or the files don't fill every global cell exactly once.
	*
	* @param layout The dumps' layout (read_bout_layout).
	* @param times Indices into layout.times to read, evenly spaced and
	*        increasing (from select_indices).
	* @param ion Name of the ion species in the Hermes-3 input
	*        (hermes_ion_species).
	* @param frames Where the plasma goes (frames.nframes() == times.size()).
	*/
	void read_bout_dumps(const BoutLayout& layout,
		const std::vector<int>& times, const std::string& ion,
		Frames& frames)
	{
		const int nt {static_cast<int>(times.size())};
		if (nt != frames.nframes() || nt < 1)
		{
			Errors::fatal_error("read_bout_dumps got "
				+ std::to_string(nt) + " times for "
				+ std::to_string(frames.nframes()) + " frames.");
		}
		const int t_step {nt > 1 ? times[1] - times[0] : 1};
		for (int t {1}; t < nt; ++t)
		{
			if (times[t] - times[t - 1] != t_step || t_step < 1)
			{
				Errors::fatal_error("read_bout_dumps needs evenly spaced, "
					"increasing times.");
			}
		}

		// The variables to read, and where the ion flow comes from
		const std::string p_ion {"P" + ion};
		for (const std::string& name : {std::string("Ne"),
			std::string("Pe"), p_ion, std::string("phi")})
		{
			if (!layout.has(name))
			{
				Errors::fatal_error("The dumps in " + layout.dir + " have no "
					+ name + (name == p_ion ? " (is hermes_ion_species = \""
					+ ion + "\" the ion species' name in the Hermes-3 "
					"input?)" : "") + ".");
			}
		}
		FlowFrom flow {FlowFrom::V};
		std::vector<std::string> names {"Ne", "Pe", p_ion, "phi"};
		if (layout.has("V" + ion))
		{
			names.push_back("V" + ion);
			std::cout << "  Ion flow: V" << ion << '\n';
		}
		else if (layout.has("NV" + ion))
		{
			names.push_back("NV" + ion);
			if (layout.has("Apar"))
			{
				flow = FlowFrom::NV_Apar;
				names.push_back("Apar");
				std::cout << "  Ion flow: (NV" << ion << " - N Apar) / (AA N), "
					"since NV" << ion << " is the canonical momentum in an "
					"electromagnetic run\n";
			}
			else
			{
				flow = FlowFrom::NV;
				std::cout << "  Ion flow: NV" << ion << " / (AA N)\n";
			}
		}
		else
		{
			Errors::fatal_error("The dumps in " + layout.dir + " have neither "
				"V" + ion + " nor NV" + ion + ", so the ion flow can't be "
				"found.");
		}

		for (int t {}; t < nt; ++t)
		{
			frames.set_time(t, layout.times[times[t]], layout.dir
				+ " output time " + std::to_string(times[t]));
		}

		const double aa {frames.norm().ion_aa};
		const std::size_t nz {static_cast<std::size_t>(layout.nz)};
		const std::size_t mysub {static_cast<std::size_t>(layout.mysub)};
		const std::size_t nt_all {layout.times.size()};
		std::vector<int> count(static_cast<std::size_t>(layout.nx)
			* layout.ny);
		const int nfiles {static_cast<int>(layout.paths.size())};
		const int every {std::max(1, nfiles / 10)};
		const auto start_time {std::chrono::steady_clock::now()};
		std::vector<std::vector<double>> data(names.size());

		for (int n {}; n < nfiles; ++n)
		{
			const std::string& path {layout.paths[n]};
			int px {}, py {};
			int x0 {}, nxl {};
			try
			{
				const netCDF::NcFile file(path, netCDF::NcFile::read);
				px = read_int(file, "PE_XIND", path);
				py = read_int(file, "PE_YIND", path);
				if (px < 0 || px >= layout.nxpe || py < 0
					|| py >= layout.nype)
				{
					Errors::fatal_error(path + " is at PE_XIND = "
						+ std::to_string(px) + ", PE_YIND = "
						+ std::to_string(py) + ", outside the NXPE x NYPE = "
						+ std::to_string(layout.nxpe) + " x "
						+ std::to_string(layout.nype) + " processes.");
				}
				if (read_times(file, path) != layout.times)
				{
					Errors::fatal_error(path + " has different output times "
						"(t_array) from BOUT.dmp.0.nc.");
				}

				// Interior cells, plus the x guards at the domain's edges
				x0 = px == 0 ? 0 : layout.mxg;
				const int x1 {layout.mxsub + layout.mxg
					+ (px == layout.nxpe - 1 ? layout.mxg : 0)};
				nxl = x1 - x0;
				const std::vector<std::size_t> shape {nt_all,
					static_cast<std::size_t>(layout.mxsub + 2 * layout.mxg),
					mysub + 2 * layout.myg, nz};
				const std::vector<std::size_t> first {
					static_cast<std::size_t>(times[0]),
					static_cast<std::size_t>(x0),
					static_cast<std::size_t>(layout.myg), 0};
				const std::vector<std::size_t> cnt {
					static_cast<std::size_t>(nt),
					static_cast<std::size_t>(nxl), mysub, nz};
				const std::vector<std::ptrdiff_t> stride {t_step, 1, 1, 1};

				for (std::size_t v {}; v < names.size(); ++v)
				{
					const netCDF::NcVar var {file.getVar(names[v])};
					if (var.isNull())
					{
						Errors::fatal_error(path + " has no variable "
							+ names[v] + ".");
					}
					std::vector<std::size_t> dims {};
					for (const auto& d : var.getDims())
						dims.push_back(d.getSize());
					if (dims != shape)
					{
						Errors::fatal_error(names[v] + " in " + path
							+ " isn't (t, x, y, z) = ("
							+ std::to_string(shape[0]) + ", "
							+ std::to_string(shape[1]) + ", "
							+ std::to_string(shape[2]) + ", "
							+ std::to_string(shape[3]) + "), as the layout "
							"says.");
					}
					data[v].resize(static_cast<std::size_t>(nt) * nxl
						* mysub * nz);
					var.getVar(first, cnt, stride, data[v].data());
				}
			}
			catch (const netCDF::exceptions::NcException& e)
			{
				Errors::fatal_error("Couldn't read the BOUT++ dump file "
					+ path + ": " + e.what());
			}

			const int gx0 {px * layout.mxsub + x0};
			const int gy0 {py * layout.mysub};
			for (int i {}; i < nxl; ++i)
			for (std::size_t j {}; j < mysub; ++j)
				++count[static_cast<std::size_t>(gx0 + i) * layout.ny + gy0 + j];

			const auto& Ne {data[0]};
			const auto& Pe {data[1]};
			const auto& Pi {data[2]};
			const auto& phi {data[3]};
			const auto& F {data[4]};
			const std::vector<double>* apar {flow == FlowFrom::NV_Apar
				? &data[5] : nullptr};
			std::size_t c {};
			for (int t {}; t < nt; ++t)
			for (int i {}; i < nxl; ++i)
			for (std::size_t j {}; j < mysub; ++j)
			for (std::size_t k {}; k < nz; ++k, ++c)
			{
				const double ne {Ne[c]};
				double V {};
				if (flow == FlowFrom::V) V = F[c];
				else if (ne > 0.0)
				{
					// Ions of charge 1, so Z N = Ne
					const double nv {flow == FlowFrom::NV_Apar
						? F[c] - ne * (*apar)[c] : F[c]};
					V = nv / (aa * ne);
				}
				frames.store(t, gx0 + i, gy0 + static_cast<int>(j),
					static_cast<int>(k), ne, Pe[c], Pi[c], V, phi[c]);
			}

			if ((n + 1) % every == 0 || n + 1 == nfiles)
			{
				const double s {std::chrono::duration<double>(
					std::chrono::steady_clock::now() - start_time).count()};
				std::cout << "  Read " << n + 1 << " of " << nfiles
					<< " dump files (" << std::lround(s) << " s)\n"
					<< std::flush;
			}
		}

		const long unfilled {std::count(count.begin(), count.end(), 0)};
		const long overfilled {std::count_if(count.begin(), count.end(),
			[](const int c) {return c > 1;})};
		if (unfilled > 0 || overfilled > 0)
		{
			Errors::fatal_error("The dump files in " + layout.dir
				+ " left " + std::to_string(unfilled) + " global cells "
				"unfilled and filled " + std::to_string(overfilled)
				+ " more than once (check PE_XIND, PE_YIND and the layout).");
		}
	}
}
