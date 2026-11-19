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

	/// Elementary charge and proton mass for unnormalizing (CODATA 2018,
	/// BOUT++'s current SI::qe and SI::Mp). Runs built with an older
	/// BOUT++ may have used an older proton mass, which changes c_s0 and
	/// the time unit by about 1e-8 relative: far below anything that
	/// matters here.
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
		const int zperiod, const bool aligned = true);

	void set_z_shift_tables(const Grid& grid, Zones::Zones& zones);

	Axisym::Mesh make_mesh(const Grid& grid);

	double parallel_sign(const Grid& grid);

	/// The run's normalization and main ions. Each output format fills it
	/// its own way (from the options for JAX-Hermes3), so the code shared
	/// by the formats doesn't read the options for them.
	struct Norm
	{
		double nnorm {};	///< Density Nnorm (m^-3)
		double tnorm {};	///< Temperature Tnorm (eV)
		double bnorm {};	///< Magnetic field Bnorm (T)
		double ion_aa {};	///< Mass number AA of the main ions
		int zperiod {};		///< The run covers 1 / zperiod of the torus
	};

	Norm norm_from_options(const Options::Options& opts);

	void check_norm(const Norm& norm);

	/**
	* The plasma of the loaded times, in SI units in the zoned arrays,
	* filled one cell at a time by a format's reader with store and
	* set_time, then moved into the background by finish_background.
	*
	* This keeps the unnormalization, the direction of the parallel flow
	* and the treatment of empty cells in one place for every format.
	*/
	class Frames
	{
	public:
		Frames(const Grid& grid, const Zones::Zones& zones,
			const Norm& norm, const int nframes);

		void store(const int t, const int x, const int y, const int k,
			const double Ne, const double Pe, const double Pi,
			const double V, const double phi);

		void set_time(const int t, const double t_norm,
			const std::string& label);

		void enable_inductive();
		bool inductive() const {return !m_epar.get_data().empty();}
		void store_dapar_dt(const int t, const int x, const int y,
			const int k, const double dapar_dt);

		int nframes() const {return m_nframes;}
		const Norm& norm() const {return m_norm;}

		friend void finish_background(Frames& frames,
			const Options::Options& opts, Background::Background& bkg);

	private:
		/// Where a global (x, y) cell is: zone (-1 if none), and the cell's
		/// index within the zone
		struct Cell
		{
			int zone {-1};
			int i {};
			int j {};
		};

		Norm m_norm {};
		int m_nframes {};
		int m_ny {};
		double m_cs0 {};	///< sqrt(e Tnorm / m_p) (m/s)
		double m_t_unit {};	///< m_p / (e Bnorm) (s)
		double m_sign {};	///< parallel_sign of the grid
		std::vector<Cell> m_cells {};	///< By x * ny + y
		Vectors::ZonedVector<BkgFPType> m_ne {};
		Vectors::ZonedVector<BkgFPType> m_te {};
		Vectors::ZonedVector<BkgFPType> m_ti {};
		Vectors::ZonedVector<BkgFPType> m_vp {};
		Vectors::ZonedVector<BkgFPType> m_upar {};

		/// Inductive E_par = -dApar/dt along B (V/m), only allocated if
		/// enable_inductive was called
		Vectors::ZonedVector<BkgFPType> m_epar {};
		const Zones::Zones* m_zones {};
		std::vector<BkgFPType> m_times {};
		std::vector<std::string> m_labels {};
		long m_empty_cells {};
	};

	void finish_background(Frames& frames, const Options::Options& opts,
		Background::Background& bkg);

	std::vector<int> select_indices(const int n, int start, int end,
		const int stride, const std::string& what);

	std::vector<std::string> select_snapshots(const std::string& dir,
		int start, int end, const int stride);

	/// What a BOUT++ Hermes-3 run's dump files say about the run, from
	/// BOUT.dmp.0.nc
	struct BoutLayout
	{
		std::string dir {};		///< Run directory
		std::vector<std::string> paths {};	///< BOUT.dmp.N.nc, N from 0
		int nxpe {};	///< Processes in x
		int nype {};	///< Processes in y
		int mxsub {};	///< Interior cells per process in x
		int mysub {};	///< Interior cells per process in y
		int mxg {};		///< x guard cells on each side
		int myg {};		///< y guard cells on each side
		int nz {};		///< Toroidal cells (MZ)
		int nx {};		///< Global x cells, x guards included
		int ny {};		///< Global y cells
		int ixseps1 {};		///< As in the grid
		int jyseps1_1 {};	///< As in the grid
		int jyseps2_2 {};	///< As in the grid
		int zperiod {};		///< The run covers 1 / zperiod of the torus
		double nnorm {};	///< Nnorm (m^-3)
		double tnorm {};	///< Tnorm (eV)
		double bnorm {};	///< Bnorm (T)
		std::vector<double> times {};	///< t_array (units of 1 / Omega_ci)
		std::vector<std::string> vars {};	///< Variables in the dumps

		bool has(const std::string& name) const;
	};

	BoutLayout read_bout_layout(const std::string& dir);

	void check_bout_layout(const BoutLayout& layout, const Grid& grid);

	Norm norm_from_bout(const BoutLayout& layout,
		const Options::Options& opts);

	void read_bout_dumps(const BoutLayout& layout,
		const std::vector<int>& times, const std::string& ion,
		Frames& frames);

	void read_snapshots(const std::vector<std::string>& paths,
		const Grid& grid, Frames& frames);

	Background::Background read_hermes(const Options::Options& opts);
}
