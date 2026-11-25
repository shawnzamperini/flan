/*
* Multi-zone backgrounds for the GPU tests, and particles spread over
* their zones, so the GPU kernels can be compared with the CPU routines on
* the same particles.
*
* - hermes_background: the HERMES test grid of hermes_test_files.h (six
*   zones with ghost cells and zShift tables, toroidal), with a plasma
*   that varies along and across the field lines and has a flow, so the
*   field-aligned interpolation and the rotation of vectors matter.
* - s3x_background: the SOLEDGE3X test mesh of s3x_test_files.h (three
*   zones with ghost cells, no zShift).
*
* Host only (the files are written with netCDF and HDF5).
*/

#pragma once

#include <cmath>
#include <filesystem>
#include <numbers>
#include <string>
#include <vector>

#include "background.h"
#include "hermes_test_files.h"
#include "indices.h"
#include "options.h"
#include "pcg32.h"
#include "read_hermes.h"
#include "read_s3xe.h"
#include "s3x_test_files.h"
#include "slots.h"

namespace MultiZone
{
	/**
	* Writes the HERMES test grid and two snapshots in dir and reads them.
	* The plasma is a function of the field-line label z - zShift, x, y
	* and time, so it changes along the field lines (y) as well as across
	* them, and has a parallel flow (NVi).
	*/
	inline Background::Background hermes_background(
		const std::filesystem::path& dir)
	{
		constexpr int nz {16};
		constexpr int zperiod {4};
		constexpr double period {2.0 * std::numbers::pi / zperiod};
		std::filesystem::create_directories(dir / "snaps");
		HermesTest::TestGrid g {};
		HermesTest::write_grid((dir / "grid.nc").string(), g);
		const double dz {period / nz};
		for (int t {}; t < 2; ++t)
		{
			auto fn = [&g, dz, t](int f, int i, int j, int k)
			{
				const double zeta {k * dz - g.z_shift(i, j)};
				const double p {std::cos(2.0 * std::numbers::pi * zeta
					/ period)};
				const double s {1.0 + 0.05 * i + 0.03 * j + 0.1 * t};
				switch (f)
				{
					case 0: return s * (1.0 + 0.5 * p);
					case 1: case 2: return s * (1.2 + 0.4 * p);
					case 3: return 0.2 * s * p;
					default: return 0.3 * s * p;
				}
			};
			HermesTest::write_snapshot((dir / "snaps").string(), g, nz, t,
				t, 10.0 * t, fn);
		}
		Options::Options opts {};
		opts.set_bkg_source("hermes");
		opts.set_hermes_grid_file((dir / "grid.nc").string());
		opts.set_hermes_snap_dir((dir / "snaps").string());
		opts.set_hermes_snap_start(0);
		opts.set_hermes_snap_end(-1);
		opts.set_hermes_tnorm(20.0);
		opts.set_hermes_zperiod(zperiod);
		return Hermes::read_hermes(opts);
	}

	/// Writes the SOLEDGE3X test mesh and two frames in dir and reads them
	inline Background::Background s3x_background(
		const std::filesystem::path& dir)
	{
		std::filesystem::create_directories(dir);
		const S3XTest::TestMesh m {};
		const int nphi {m.nphi};
		S3XTest::write_case(dir.string(), m, [nphi](const std::string& path)
			{return S3X::read_mesh(path, nphi).zones;});
		Options::Options opts {};
		opts.set_bkg_source("soledge3x");
		opts.set_s3x_dir(dir.string());
		opts.set_s3x_frame_start(0);
		opts.set_s3x_frame_end(1);
		return S3X::read_s3xe(opts);
	}

	/**
	* Live tungsten particles at random times and positions, n in each
	* zone, with their cell indices set, random charges 1 to 10 and
	* velocities of 1e4 m/s in random directions.
	*
	* @param bkg Background to put them in.
	* @param n Particles per zone.
	* @param seed Seed for the positions.
	* @return The particles.
	*/
	inline Slots::Slots spread_particles(const Background::Background& bkg,
		const int n, const int seed)
	{
		const int nzones {bkg.get_zones().get_nzones()};
		Slots::Slots slots(n * nzones);
		slots.set_mass(183.84 * 1.66053906660e-27);
		slots.set_Z(74);
		pcg32 rng(seed, 1);
		const auto& times {bkg.get_times()};
		for (int zone {}; zone < nzones; ++zone)
		{
			const auto& gx {bkg.get_grid_x(zone)};
			const auto& gy {bkg.get_grid_y(zone)};
			const auto& gz {bkg.get_grid_z(zone)};
			for (int p {}; p < n; ++p)
			{
				const int i {zone * n + p};
				const double t {times.front() + (times.back()
					- times.front()) * rng.next_double()};
				const double x {gx.front() + (gx.back() - gx.front())
					* rng.next_double()};
				const double y {gy.front() + (gy.back() - gy.front())
					* rng.next_double()};
				const double z {gz.front() + (gz.back() - gz.front())
					* rng.next_double()};
				slots.set_state(i, 0);
				slots.set_zone(i, zone);
				slots.set_t(i, t);
				slots.set_x(i, x);
				slots.set_y(i, y);
				slots.set_z(i, z);
				slots.set_tidx(i, Indices::get_nearest_index_cpu(times, t));
				slots.set_xidx(i, Indices::get_nearest_cell_index_cpu(gx, x));
				slots.set_yidx(i, Indices::get_nearest_cell_index_cpu(gy, y));
				slots.set_zidx(i, Indices::get_nearest_cell_index_cpu(gz, z));
				slots.set_q(i, 1 + static_cast<int>(10 * rng.next_double())
					% 10);
				slots.set_weight(i, 1.0 + rng.next_double());
				const double ct {2.0 * rng.next_double() - 1.0};
				const double st {std::sqrt(1.0 - ct * ct)};
				const double ph {2.0 * std::numbers::pi * rng.next_double()};
				slots.set_vX(i, 1e4 * st * std::cos(ph));
				slots.set_vY(i, 1e4 * st * std::sin(ph));
				slots.set_vZ(i, 1e4 * ct);
			}
		}
		return slots;
	}
}
