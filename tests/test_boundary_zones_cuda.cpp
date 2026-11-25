/*
* The GPU boundary check (check_bounds_gpu) against the CPU's
* (check_bounds_cpu), particle by particle, on:
*
* - the SOLEDGE3X test mesh: crossings between zones, the core faces split
*   over two zones (one picked at random, weighted by cells), the target
*   wall faces, absorbing faces and the wall mask chi;
* - the HERMES test grid: crossings between six zones, the core and the
*   periodic z faces of a toroidal wedge, where the velocity turns with
*   the particle;
* - the slab: one zone with a core x face, a periodic y and an lcfs z
*   boundary.
*
* The particles are spread over each zone and up to 40% of its size past
* its faces in every direction, so they cross one or several faces, and
* some are past the last background time. Slot i's generator on the GPU is
* pcg32(seed + i, i) (init_rngs_cuda), so the CPU is run one particle at a
* time with that generator, and the two must agree exactly where a core
* face draws random numbers too. The fates counted on the device must
* match the CPU's.
*/

#include <algorithm>
#include <cmath>
#include <cuda_runtime.h>
#include <filesystem>
#include <gtest/gtest.h>
#include <omp.h>
#include <string>
#include <unistd.h>
#include <vector>

#include "background.h"
#include "background_device.h"
#include "boundary.h"
#include "boundary_hd.h"
#include "multizone_cases.h"
#include "options.h"
#include "pcg32.h"
#include "read_bkg.h"
#include "slots.h"
#include "slots_device.h"

#include "boundary.cuh"
#include "init_rngs.cuh"

namespace
{
	constexpr std::uint64_t seed {77};

	/// Moves the particles past their zones' faces: each coordinate is
	/// uniform over the zone extended by 40% on both sides, and one in
	/// twenty is past the last background time
	void spread_past_faces(Slots::Slots& slots,
		const Background::Background& bkg)
	{
		pcg32 rng(5, 9);
		for (int i {}; i < slots.N(); ++i)
		{
			const int zone {slots.zone()[i]};
			const auto& zones {bkg.get_zones()};
			double pos[3] {};
			for (int d {}; d < 3; ++d)
			{
				const double lo {zones.get_min(zone, d)};
				const double hi {zones.get_max(zone, d)};
				pos[d] = lo + (hi - lo) * (1.8 * rng.next_double() - 0.4);
			}
			slots.set_x(i, pos[0]);
			slots.set_y(i, pos[1]);
			slots.set_z(i, pos[2]);
			if (rng.next_double() < 0.05)
				slots.set_t(i, bkg.get_t_max() * 1.01 + 1e-9);
		}
	}

	/// Number of the given fate over the per-particle CPU runs
	struct CpuResult
	{
		Slots::Slots slots;
		Boundary::Fates fates;
	};

	/// check_bounds_cpu on one particle at a time, each with slot i's GPU
	/// generator (every OpenMP thread gets a copy, since any may run it)
	CpuResult cpu_one_by_one(const Slots::Slots& in,
		const Background::Background& bkg, const Options::Options& opts)
	{
		CpuResult r {in, {}};
		for (int i {}; i < in.N(); ++i)
		{
			Slots::Slots one(1);
			one.set_mass(in.mass());
			one.set_Z(in.Z());
			one.set_state(0, in.state()[i]);
			one.set_zone(0, in.zone()[i]);
			one.set_t(0, in.t()[i]);
			one.set_x(0, in.x()[i]);
			one.set_y(0, in.y()[i]);
			one.set_z(0, in.z()[i]);
			one.set_vX(0, in.vX()[i]);
			one.set_vY(0, in.vY()[i]);
			one.set_vZ(0, in.vZ()[i]);
			std::vector<pcg32> rngs(omp_get_max_threads(),
				pcg32(seed + i, i));
			Boundary::check_bounds_cpu(one, bkg, opts, rngs, r.fates);
			r.slots.set_state(i, one.state()[0]);
			r.slots.set_zone(i, one.zone()[0]);
			r.slots.set_t(i, one.t()[0]);
			r.slots.set_x(i, one.x()[0]);
			r.slots.set_y(i, one.y()[0]);
			r.slots.set_z(i, one.z()[0]);
			r.slots.set_vX(i, one.vX()[0]);
			r.slots.set_vY(i, one.vY()[0]);
		}
		return r;
	}

	/**
	* Runs both and compares them.
	*
	* @param bkg_in Background.
	* @param opts Options (boundary types, buffers, lcfs_x, chi).
	* @param per_zone Particles per zone.
	* @return The CPU's fates, for checks of which ones happened.
	*/
	Boundary::Fates bounds_match_cpu(const Background::Background& bkg_in,
		const Options::Options& opts, const int per_zone)
	{
		Background::Background bkg {bkg_in};
		Slots::Slots slots {MultiZone::spread_particles(bkg, per_zone, 13)};
		spread_past_faces(slots, bkg);

		const CpuResult cpu {cpu_one_by_one(slots, bkg, opts)};

		Background::BackgroundDevice bkg_d {bkg.to_device()};
		Slots::Slots gpu {slots};
		Slots::SlotsDevice slots_d {gpu.to_device()};
		pcg32* rngs_d {init_rngs_cuda(slots_d, seed)};
		unsigned long long* fates_d {};
		const std::size_t fates_size {Boundary::num_fate_indices
			* sizeof(unsigned long long)};
		cudaMalloc(&fates_d, fates_size);
		cudaMemset(fates_d, 0, fates_size);
		Boundary::check_bounds_gpu(slots_d, bkg_d,
			Boundary::make_params(bkg, opts, bkg_d.chi), rngs_d, fates_d);
		gpu = gpu.to_host(slots_d);
		Boundary::Fates gpu_fates {};
		Boundary::add_fates_from_device(fates_d, 0, gpu_fates);
		cudaFree(fates_d);
		free_rngs(rngs_d, slots_d.device_id);
		Slots::free_slots(slots_d);
		Background::free_bkg(bkg_d);

		int bad {};
		for (int i {}; i < slots.N(); ++i)
		{
			const auto& c {cpu.slots};
			// Positions and velocities to rounding: a crossing or wrap is
			// an addition, and the turn a rotation
			auto near = [](double a, double b)
				{return std::abs(a - b) <= 1e-12 * (1.0 + std::abs(b));};
			const bool ok {gpu.state()[i] == c.state()[i]
				&& gpu.zone()[i] == c.zone()[i]
				&& near(gpu.t()[i], c.t()[i]) && near(gpu.x()[i], c.x()[i])
				&& near(gpu.y()[i], c.y()[i]) && near(gpu.z()[i], c.z()[i])
				&& near(gpu.vX()[i], c.vX()[i])
				&& near(gpu.vY()[i], c.vY()[i])};
			if (!ok && bad++ < 5)
			{
				ADD_FAILURE() << "particle " << i << " from zone "
					<< slots.zone()[i] << " at (" << slots.x()[i] << ", "
					<< slots.y()[i] << ", " << slots.z()[i] << "): GPU state "
					<< gpu.state()[i] << " zone " << gpu.zone()[i] << " ("
					<< gpu.x()[i] << ", " << gpu.y()[i] << ", " << gpu.z()[i]
					<< "), CPU state " << c.state()[i] << " zone "
					<< c.zone()[i] << " (" << c.x()[i] << ", " << c.y()[i]
					<< ", " << c.z()[i] << ")";
			}
		}
		EXPECT_EQ(bad, 0);

		const Boundary::Fates& f {cpu.fates};
		EXPECT_EQ(gpu_fates.crossings, f.crossings);
		EXPECT_EQ(gpu_fates.absorbed, f.absorbed);
		EXPECT_EQ(gpu_fates.wall, f.wall);
		EXPECT_EQ(gpu_fates.chi, f.chi);
		EXPECT_EQ(gpu_fates.core, f.core);
		EXPECT_EQ(gpu_fates.time_end, f.time_end);
		EXPECT_EQ(gpu_fates.runaway, f.runaway);
		return f;
	}

	std::filesystem::path test_dir(const std::string& name)
	{
		return std::filesystem::temp_directory_path()
			/ ("flan_test_boundary_zones_" + name + "_"
			+ std::to_string(getpid()));
	}
}

TEST(BoundaryZonesCUDA, MatchesCPUOnS3XGrid)
{
	const auto dir {test_dir("s3x")};
	const Background::Background bkg {MultiZone::s3x_background(dir)};
	std::filesystem::remove_all(dir);
	ASSERT_FALSE(bkg.get_chi().get_data().empty());

	Options::Options opts {};
	opts.set_s3x_chi_absorb("on");
	opts.set_imp_xbound_buffer(1e-4);
	const Boundary::Fates f {bounds_match_cpu(bkg, opts, 1500)};

	// Every kind of fate the mesh has happened
	EXPECT_GT(f.crossings, 0);
	EXPECT_GT(f.absorbed, 0);
	EXPECT_GT(f.wall, 0);
	EXPECT_GT(f.chi, 0);
	EXPECT_GT(f.core, 0);
	EXPECT_GT(f.time_end, 0);
}

TEST(BoundaryZonesCUDA, MatchesCPUOnHermesWedge)
{
	const auto dir {test_dir("hermes")};
	const Background::Background bkg {MultiZone::hermes_background(dir)};
	std::filesystem::remove_all(dir);
	ASSERT_TRUE(bkg.get_zones().any_wedge());

	const Options::Options opts {};
	const Boundary::Fates f {bounds_match_cpu(bkg, opts, 800)};
	EXPECT_GT(f.crossings, 0);
	EXPECT_GT(f.absorbed + f.wall, 0);
	EXPECT_GT(f.core, 0);
}

TEST(BoundaryZonesCUDA, MatchesCPUOnSlabWithLcfs)
{
	Options::Options opts {};
	opts.set_bkg_source("slab");
	opts.set_min_xbound_type("core");
	opts.set_max_xbound_type("absorbing");
	opts.set_min_ybound_type("periodic");
	opts.set_max_ybound_type("periodic");
	opts.set_min_zbound_type("lcfs");
	opts.set_max_zbound_type("lcfs");
	opts.set_lcfs_x(0.025);
	opts.set_imp_zbound_buffer(1e-4);
	const Background::Background bkg {Background::read_bkg(opts)};
	ASSERT_EQ(bkg.get_zones().get_face(0, Zones::zmin), Zones::face_lcfs);

	const Boundary::Fates f {bounds_match_cpu(bkg, opts, 3000)};
	EXPECT_EQ(f.crossings, 0);
	EXPECT_GT(f.absorbed, 0);
	EXPECT_GT(f.core, 0);
}
