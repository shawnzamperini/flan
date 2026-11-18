/*
* Writes small JAX-Hermes3-style files for tests: .npz snapshots laid out
* the way numpy.savez writes them, and a single-null hypnotoad-style grid,
* so the tests don't depend on real HERMES output.
*
* The test grid is a rectangle in (R, Z): cell corner (i, j) of the global
* index grid is at R = R0 + i dR, Z = Z0 + j dZ, with B = (Br, Bt, Bz)
* the same everywhere and Bz > 0, so the poloidal field points along +y.
* Only the index topology (where the separatrix and legs are) makes it a
* single null; the geometry doesn't have to look like a tokamak for the
* reader.
*
* zShift is what the field gives on this grid: along a field line,
* dz/dy = (B . grad z) / (B . grad y) = Bt dZ / (R Bz) (see twist_per_cell),
* integrated in y from 0 through the inner leg, core and outer leg in the
* SOL, and from one PFR leg straight into the other (skipping the core
* rows, as hypnotoad does). The core's ShiftAngle is its twist over the
* core rows; it's NaN in the SOL as in hypnotoad files.
*
* The snapshot fields are simple functions of the global indices (see
* field_value), so a test can check which grid cell each Flan cell came
* from, or any function given to write_snapshot.
*
* write_bout_dumps writes the same plasma as BOUT++ Hermes-3 dump files,
* split over processes with guard cells, so the two readers can be
* compared on the same plasma.
*/

#pragma once

#include <cstdint>
#include <cstring>
#include <fstream>
#include <functional>
#include <limits>
#include <map>
#include <netcdf>
#include <string>
#include <vector>

namespace HermesTest
{
	// ---------------------------------------------------------------------
	// .npy and .npz writing
	// ---------------------------------------------------------------------

	/// An array to write: its .npy dtype, shape and raw bytes
	struct NpyArray
	{
		std::string descr {"<f8"};
		std::vector<std::size_t> shape {};
		std::vector<char> bytes {};
		bool fortran {false};
	};

	/// A float64 array from values
	inline NpyArray f8(std::vector<std::size_t> shape,
		const std::vector<double>& v)
	{
		NpyArray a {"<f8", std::move(shape), {}};
		a.bytes.resize(v.size() * 8);
		std::memcpy(a.bytes.data(), v.data(), a.bytes.size());
		return a;
	}

	/// The bytes of a .npy file (version 1.0), header padded to 64 bytes
	/// as numpy does
	inline std::vector<char> npy_bytes(const NpyArray& a)
	{
		std::string shape {"("};
		for (std::size_t n {}; n < a.shape.size(); ++n)
			shape += (n ? ", " : "") + std::to_string(a.shape[n]);
		if (a.shape.size() == 1) shape += ",";
		shape += ")";
		std::string header {"{'descr': '" + a.descr + "', 'fortran_order': "
			+ (a.fortran ? "True" : "False") + ", 'shape': " + shape
			+ ", }"};
		while ((10 + header.size() + 1) % 64 != 0) header += ' ';
		header += '\n';

		std::vector<char> out {'\x93', 'N', 'U', 'M', 'P', 'Y', 1, 0};
		const std::uint16_t len {static_cast<std::uint16_t>(header.size())};
		out.push_back(static_cast<char>(len & 0xFF));
		out.push_back(static_cast<char>(len >> 8));
		out.insert(out.end(), header.begin(), header.end());
		out.insert(out.end(), a.bytes.begin(), a.bytes.end());
		return out;
	}

	namespace detail
	{
		template <typename T>
		void put(std::vector<char>& b, T v)
		{
			char c[sizeof(T)];
			std::memcpy(c, &v, sizeof(T));
			b.insert(b.end(), c, c + sizeof(T));
		}
	}

	/// How to lay out the zip archive
	struct ZipStyle
	{
		/// Local headers with 0xFFFFFFFF sizes and a zip64 extra field,
		/// as numpy.savez writes them
		bool zip64_local {true};
		/// Central directory offsets in a zip64 extra field instead
		bool zip64_central {false};
		/// Zip method to claim (8 = deflate, which the reader refuses)
		std::uint16_t method {0};
	};

	/**
	* Writes a .npz file holding the arrays (each as name.npy), stored
	* without compression, with the layout given by style.
	*/
	inline void write_npz(const std::string& path,
		const std::map<std::string, NpyArray>& arrays,
		const ZipStyle style = {})
	{
		using detail::put;
		std::vector<char> file {}, central {};
		std::uint16_t count {};
		for (const auto& [name, a] : arrays)
		{
			const std::string fname {name + ".npy"};
			const auto data {npy_bytes(a)};
			const std::uint64_t size {data.size()};
			const std::uint64_t offset {file.size()};

			// Local header
			put<std::uint32_t>(file, 0x04034b50);
			put<std::uint16_t>(file, 45);
			put<std::uint16_t>(file, 0);
			put<std::uint16_t>(file, style.method);
			put<std::uint32_t>(file, 0);	// Time and date
			put<std::uint32_t>(file, 0);	// CRC (not checked)
			put<std::uint32_t>(file, style.zip64_local ? 0xFFFFFFFF
				: static_cast<std::uint32_t>(size));
			put<std::uint32_t>(file, style.zip64_local ? 0xFFFFFFFF
				: static_cast<std::uint32_t>(size));
			put<std::uint16_t>(file, static_cast<std::uint16_t>(fname.size()));
			put<std::uint16_t>(file, style.zip64_local ? 20 : 0);
			file.insert(file.end(), fname.begin(), fname.end());
			if (style.zip64_local)
			{
				put<std::uint16_t>(file, 1);
				put<std::uint16_t>(file, 16);
				put<std::uint64_t>(file, size);
				put<std::uint64_t>(file, size);
			}
			file.insert(file.end(), data.begin(), data.end());

			// Central directory entry
			put<std::uint32_t>(central, 0x02014b50);
			put<std::uint16_t>(central, 45);
			put<std::uint16_t>(central, 45);
			put<std::uint16_t>(central, 0);
			put<std::uint16_t>(central, style.method);
			put<std::uint32_t>(central, 0);
			put<std::uint32_t>(central, 0);
			put<std::uint32_t>(central, static_cast<std::uint32_t>(size));
			put<std::uint32_t>(central, static_cast<std::uint32_t>(size));
			put<std::uint16_t>(central,
				static_cast<std::uint16_t>(fname.size()));
			put<std::uint16_t>(central, style.zip64_central ? 12 : 0);
			put<std::uint16_t>(central, 0);	// Comment
			put<std::uint16_t>(central, 0);	// Disk
			put<std::uint16_t>(central, 0);	// Internal attributes
			put<std::uint32_t>(central, 0);	// External attributes
			put<std::uint32_t>(central, style.zip64_central ? 0xFFFFFFFF
				: static_cast<std::uint32_t>(offset));
			central.insert(central.end(), fname.begin(), fname.end());
			if (style.zip64_central)
			{
				put<std::uint16_t>(central, 1);
				put<std::uint16_t>(central, 8);
				put<std::uint64_t>(central, offset);
			}
			++count;
		}
		const std::uint32_t cd_offset {static_cast<std::uint32_t>(
			file.size())};
		file.insert(file.end(), central.begin(), central.end());
		put<std::uint32_t>(file, 0x06054b50);
		put<std::uint16_t>(file, 0);
		put<std::uint16_t>(file, 0);
		put<std::uint16_t>(file, count);
		put<std::uint16_t>(file, count);
		put<std::uint32_t>(file, static_cast<std::uint32_t>(central.size()));
		put<std::uint32_t>(file, cd_offset);
		put<std::uint16_t>(file, 0);
		std::ofstream(path, std::ios::binary).write(file.data(),
			static_cast<std::streamsize>(file.size()));
	}

	// ---------------------------------------------------------------------
	// Grid and snapshots
	// ---------------------------------------------------------------------

	/// Layout of the test grid. Physical x cells are 2 to nx - 3.
	struct TestGrid
	{
		int nx {10};
		int ny {8};
		int ixseps1 {5};
		int jyseps1_1 {1};
		int jyseps2_2 {5};
		int jyseps2_1 {4};
		int jyseps1_2 {4};
		double R0 {1.0};
		double Z0 {-0.3};
		double dR {0.05};
		double dZ {0.08};
		double Br {0.02};
		double Bz {0.3};	///< > 0, so +y is along the poloidal field
		double Bt {-1.5};
		int bz_flip_row {-1};	///< y row with Bz reversed (-1: none)

		/// If nonzero, Bt is set in each column so the field line twists
		/// by this much per y cell everywhere (Bt then grows with R)
		double twist {0.0};

		/// Center R of x column i
		double R_center(int i) const {return R0 + (i + 0.5) * dR;}

		/// Bt in column i
		double Bt_at(int i) const
			{return twist != 0.0 ? twist * Bz * R_center(i) / dZ : Bt;}

		/// Toroidal twist of the field line per y cell in column i
		double twist_per_cell(int i) const
			{return Bt_at(i) * dZ / (R_center(i) * Bz);}

		/// zShift of cell (i, j) (see the file comment)
		double z_shift(int i, int j) const
		{
			const bool pfr_out {i < ixseps1 && j > jyseps2_2};
			return twist_per_cell(i) * (pfr_out ? j - (jyseps2_2
				- jyseps1_1) : j);
		}

		/// ShiftAngle of column i
		double shift_angle(int i) const
		{
			return i < ixseps1 ? twist_per_cell(i) * (jyseps2_2 - jyseps1_1)
				: std::numeric_limits<double>::quiet_NaN();
		}
	};

	/**
	* Writes the test grid as a hypnotoad-style netCDF file: the scalars,
	* the four corner arrays, the field and zShift on (x, y), and
	* ShiftAngle on x.
	*/
	inline void write_grid(const std::string& path, const TestGrid& g)
	{
		netCDF::NcFile f(path, netCDF::NcFile::replace);
		const auto dx {f.addDim("x", g.nx)};
		const auto dy {f.addDim("y", g.ny)};
		auto scalar = [&](const std::string& name, int v)
			{f.addVar(name, netCDF::ncInt).putVar(&v);};
		scalar("nx", g.nx);
		scalar("ny", g.ny);
		scalar("ixseps1", g.ixseps1);
		scalar("ixseps2", g.nx);
		scalar("jyseps1_1", g.jyseps1_1);
		scalar("jyseps2_1", g.jyseps2_1);
		scalar("jyseps1_2", g.jyseps1_2);
		scalar("jyseps2_2", g.jyseps2_2);
		scalar("y_boundary_guards", 0);

		auto field = [&](const std::string& name, auto fn)
		{
			std::vector<double> v {};
			for (int i {}; i < g.nx; ++i)
				for (int j {}; j < g.ny; ++j) v.push_back(fn(i, j));
			f.addVar(name, netCDF::ncDouble, {dx, dy}).putVar(v.data());
		};
		auto R = [&](double i) {return g.R0 + i * g.dR;};
		auto Z = [&](double j) {return g.Z0 + j * g.dZ;};
		field("Rxy_corners", [&](int i, int) {return R(i);});
		field("Rxy_lower_right_corners", [&](int i, int) {return R(i + 1);});
		field("Rxy_upper_left_corners", [&](int i, int) {return R(i);});
		field("Rxy_upper_right_corners", [&](int i, int) {return R(i + 1);});
		field("Zxy_corners", [&](int, int j) {return Z(j);});
		field("Zxy_lower_right_corners", [&](int, int j) {return Z(j);});
		field("Zxy_upper_left_corners", [&](int, int j) {return Z(j + 1);});
		field("Zxy_upper_right_corners", [&](int, int j) {return Z(j + 1);});
		field("Brxy", [&](int, int) {return g.Br;});
		field("Bzxy", [&](int, int j)
			{return j == g.bz_flip_row ? -g.Bz : g.Bz;});
		field("Btxy", [&](int i, int) {return g.Bt_at(i);});
		field("zShift", [&](int i, int j) {return g.z_shift(i, j);});
		std::vector<double> shift(g.nx);
		for (int i {}; i < g.nx; ++i) shift[i] = g.shift_angle(i);
		f.addVar("ShiftAngle", netCDF::ncDouble, {dx}).putVar(shift.data());
	}

	/// Value of a normalized snapshot field at global cell (i, j, k) at
	/// time index t, distinct for every field (f = 0 to 4), cell and time
	inline double field_value(int f, int t, int i, int j, int k)
	{
		return 1.0 + f + 0.1 * i + 0.01 * j + 0.001 * k + 0.0001 * t;
	}

	/// A normalized snapshot field's value from (field, i, j, k), with
	/// fields 0 to 4 being Ne, Pe, Pi, NVi and phi
	using FieldFn = std::function<double(int, int, int, int)>;

	/**
	* Writes snapshot snap<number>.npz in dir: Ne, Pe, Pi, NVi and phi
	* (fields 0 to 4 of field_value, or of fn if given) on (nx, ny, nz),
	* Vort (unused), and the normalized time t.
	*/
	inline void write_snapshot(const std::string& dir, const TestGrid& g,
		const int nz, const int number, const int t_index, const double t,
		const FieldFn& fn = {})
	{
		const std::vector<std::string> names {"Ne", "Pe", "Pi", "NVi",
			"phi"};
		std::map<std::string, NpyArray> arrays {};
		const std::vector<std::size_t> shape {
			static_cast<std::size_t>(g.nx), static_cast<std::size_t>(g.ny),
			static_cast<std::size_t>(nz)};
		for (std::size_t f {}; f < names.size(); ++f)
		{
			std::vector<double> v {};
			for (int i {}; i < g.nx; ++i)
			for (int j {}; j < g.ny; ++j)
			for (int k {}; k < nz; ++k)
				v.push_back(fn ? fn(static_cast<int>(f), i, j, k)
					: field_value(static_cast<int>(f), t_index, i, j, k));
			arrays[names[f]] = f8(shape, v);
		}
		arrays["Vort"] = arrays["phi"];
		arrays["t"] = f8({}, {t});
		char name[32];
		std::snprintf(name, sizeof(name), "/snap%07d.npz", number);
		write_npz(dir + name, arrays);
	}

	// ---------------------------------------------------------------------
	// BOUT++ Hermes-3 dump files
	// ---------------------------------------------------------------------

	/// How write_bout_dumps lays out and fills the dump files
	struct DumpSpec
	{
		int nxpe {2};		///< Processes in x; nx - 2 mxg must divide
		int nype {4};		///< Processes in y; ny must divide
		int mxg {2};
		int myg {2};
		int nz {4};
		int ntimes {3};		///< Output times, t = 100 (n + 1)
		int zperiod {2};
		double nnorm {2e19};
		double tnorm {40.0};
		double bnorm {1.5};
		double aa {2.0};		///< Ion mass number, for V and NV
		std::string ion {"i"};	///< Ion species name
		bool write_v {true};	///< Write V<ion>
		bool write_nv {true};	///< Write NV<ion>
		bool write_apar {true};	///< Write Apar (NV is then canonical)
		int nzpe {1};			///< NZPE written (the reader needs 1)
		int ixseps1_shift {0};	///< Added to the ixseps1 written
		int bad_times_file {-1};	///< File whose t_array is changed
		bool duplicate_place {false};	///< File 1 says it's at file 0's place
	};

	/// Apar (normalized) at global cell (i, j, k), time index t. Not
	/// linear in t, and its rate of change depends on the cell, so
	/// dApar/dt shows which output times and cells were used.
	inline double apar_value(int t, int i, int j, int k)
	{
		return 0.3 + 0.02 * i - 0.01 * j + 0.005 * k + 0.001 * t
			+ 0.0002 * t * t + 0.0001 * t * (i + 2 * j + 3 * k);
	}

	/**
	* Writes BOUT.dmp.N.nc files in dir with the plasma of field_value:
	* Ne, Pe, P<ion> and phi are fields 0, 1, 2 and 4, and the ion flow is
	* V = field 3 / (AA Ne), written as V<ion> and/or as NV<ion> = AA Ne V
	* + Ne Apar (the canonical momentum, or AA Ne V without Apar). So the
	* plasma is that of write_snapshot, whose NVi is field 3.
	*
	* File N is the process at PE_XIND = N % nxpe, PE_YIND = N / nxpe, as
	* in BOUT++. Each block has its guard cells; the y guards, and the x
	* guards between blocks, are NaN, so a reader that uses them gives NaN.
	* The x guards at the domain's edges are the grid's guard cells and get
	* plasma.
	*/
	inline void write_bout_dumps(const std::string& dir, const TestGrid& g,
		const DumpSpec& s)
	{
		const int mxsub {(g.nx - 2 * s.mxg) / s.nxpe};
		const int mysub {g.ny / s.nype};
		const int lx {mxsub + 2 * s.mxg};
		const int ly {mysub + 2 * s.myg};
		const double nan {std::numeric_limits<double>::quiet_NaN()};
		for (int n {}; n < s.nxpe * s.nype; ++n)
		{
			const int px {n % s.nxpe};
			const int py {n / s.nxpe};
			netCDF::NcFile f(dir + "/BOUT.dmp." + std::to_string(n) + ".nc",
				netCDF::NcFile::replace);
			const auto dt {f.addDim("t", s.ntimes)};
			const auto dx {f.addDim("x", lx)};
			const auto dy {f.addDim("y", ly)};
			const auto dz {f.addDim("z", s.nz)};
			auto iscalar = [&](const std::string& name, int v)
				{f.addVar(name, netCDF::ncInt).putVar(&v);};
			auto dscalar = [&](const std::string& name, double v)
				{f.addVar(name, netCDF::ncDouble).putVar(&v);};
			iscalar("NXPE", s.nxpe);
			iscalar("NYPE", s.nype);
			iscalar("NZPE", s.nzpe);
			iscalar("MXSUB", mxsub);
			iscalar("MYSUB", mysub);
			iscalar("MXG", s.mxg);
			iscalar("MYG", s.myg);
			iscalar("MZ", s.nz);
			iscalar("zperiod", s.zperiod);
			iscalar("nx", g.nx);
			iscalar("ny", g.ny);
			iscalar("ixseps1", g.ixseps1 + s.ixseps1_shift);
			iscalar("jyseps1_1", g.jyseps1_1);
			iscalar("jyseps2_2", g.jyseps2_2);
			iscalar("PE_XIND", s.duplicate_place && n == 1 ? 0 : px);
			iscalar("PE_YIND", s.duplicate_place && n == 1 ? 0 : py);
			dscalar("Nnorm", s.nnorm);
			dscalar("Tnorm", s.tnorm);
			dscalar("Bnorm", s.bnorm);

			std::vector<double> t(s.ntimes);
			for (int i {}; i < s.ntimes; ++i)
				t[i] = 100.0 * (i + 1) + (n == s.bad_times_file ? 1.0 : 0.0);
			f.addVar("t_array", netCDF::ncDouble, {dt}).putVar(t.data());

			// Values of one variable from the global-cell plasma
			auto var = [&](const std::string& name, auto fn)
			{
				std::vector<double> v {};
				for (int it {}; it < s.ntimes; ++it)
				for (int x {}; x < lx; ++x)
				for (int y {}; y < ly; ++y)
				for (int k {}; k < s.nz; ++k)
				{
					const bool yguard {y < s.myg || y >= s.myg + mysub};
					const bool inner_xguard {(px > 0 && x < s.mxg)
						|| (px < s.nxpe - 1 && x >= s.mxg + mxsub)};
					if (yguard || inner_xguard) v.push_back(nan);
					else
					{
						v.push_back(fn(it, px * mxsub + x,
							py * mysub + y - s.myg, k));
					}
				}
				f.addVar(name, netCDF::ncDouble, {dt, dx, dy, dz})
					.putVar(v.data());
			};
			auto ne = [](int it, int i, int j, int k)
				{return field_value(0, it, i, j, k);};
			auto v_par = [&](int it, int i, int j, int k)
				{return field_value(3, it, i, j, k) / (s.aa * ne(it, i, j, k));};
			auto apar = [&](int it, int i, int j, int k)
				{return s.write_apar ? apar_value(it, i, j, k) : 0.0;};
			var("Ne", ne);
			var("Pe", [](int it, int i, int j, int k)
				{return field_value(1, it, i, j, k);});
			var("P" + s.ion, [](int it, int i, int j, int k)
				{return field_value(2, it, i, j, k);});
			var("phi", [](int it, int i, int j, int k)
				{return field_value(4, it, i, j, k);});
			if (s.write_v) var("V" + s.ion, v_par);
			if (s.write_nv)
			{
				var("NV" + s.ion, [&](int it, int i, int j, int k)
				{
					const double n0 {ne(it, i, j, k)};
					return s.aa * n0 * v_par(it, i, j, k)
						+ n0 * apar(it, i, j, k);
				});
			}
			if (s.write_apar) var("Apar", apar);
		}
	}
}

