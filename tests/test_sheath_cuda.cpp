/*
* The sheath field and the near-wall sub-steps on the GPU against the CPU:
*
* - the slab with a wall at y min (Brooks' profile, B at a few degrees to
*   the wall), a Cartesian zone;
* - the HERMES test grid's targets (Borodkina's profile), on a toroidal
*   wedge with field-aligned zones.
*
* Ions (and some neutrals) are put in the plasma cells next to the wall,
* from a few hundredths of a cell to the wall itself, moving in random
* directions. Two things are compared:
*
* - one Boris push with the sheath field added (update_velocity_gpu
*   against update_velocity_cpu), no sub-steps;
* - the full near-wall step: the sub-step counts (count_substeps_kernel
*   against Sheath::substeps), then the push of the rest and the
*   sub-stepped push and move of these (substep_kernel against
*   Boris::substep_cpu).
*
* Nothing is random, so the two must agree to rounding: velocities and
* positions to 1e-8 relative after up to hundreds of sub-steps, counts
* exactly.
*/

#include <algorithm>
#include <cmath>
#include <cuda_runtime.h>
#include <filesystem>
#include <gtest/gtest.h>
#include <string>
#include <unistd.h>
#include <vector>

#include "background.h"
#include "background_device.h"
#include "boris.h"
#include "constants.h"
#include "indices.h"
#include "multizone_cases.h"
#include "options.h"
#include "pcg32.h"
#include "read_bkg.h"
#include "sheath.h"
#include "sheath_hd.h"
#include "slots.h"
#include "slots_device.h"
#include "sputtering.h"

#include "boris.cuh"

namespace
{
	/// Ions in the plasma cells next to the wall elements, from the face
	/// to a fraction reach of the cell across it
	Slots::Slots near_wall(const Background::Background& bkg,
		const Sputtering::Surface& surface, const int n, const double reach,
		const double speed)
	{
		Slots::Slots slots(n);
		slots.set_mass(183.84 * 1.66053906660e-27);
		slots.set_Z(74);
		pcg32 rng(21, 4);
		const auto& els {surface.elements()};
		const auto& times {bkg.get_times()};
		const auto& zones {bkg.get_zones()};
		for (int p {}; p < n; ++p)
		{
			const auto& el {els[static_cast<std::size_t>(rng.next_double()
				* els.size()) % els.size()]};
			const int zone {el.zone};
			const auto& gx {zones.get_grid_x(zone)};
			const auto& gy {zones.get_grid_y(zone)};
			const auto& gz {zones.get_grid_z(zone)};
			double x {gx[el.i] + rng.next_double() * (gx[el.i + 1]
				- gx[el.i])};
			double y {gy[el.j] + rng.next_double() * (gy[el.j + 1]
				- gy[el.j])};
			const double z {gz[el.k] + rng.next_double() * (gz[el.k + 1]
				- gz[el.k])};
			const double f {reach * rng.next_double()};
			switch (el.face)
			{
				case Zones::xmin: x = gx[el.i] + f * (gx[el.i + 1]
					- gx[el.i]); break;
				case Zones::xmax: x = gx[el.i + 1] - f * (gx[el.i + 1]
					- gx[el.i]); break;
				case Zones::ymin: y = gy[el.j] + f * (gy[el.j + 1]
					- gy[el.j]); break;
				default: y = gy[el.j + 1] - f * (gy[el.j + 1]
					- gy[el.j]); break;
			}
			const double t {times.front() + (times.back() - times.front())
				* rng.next_double()};
			slots.set_state(p, 0);
			slots.set_zone(p, zone);
			slots.set_t(p, t);
			slots.set_x(p, x);
			slots.set_y(p, y);
			slots.set_z(p, z);
			slots.set_tidx(p, Indices::get_nearest_index_cpu(times, t));
			slots.set_xidx(p, Indices::get_nearest_cell_index_cpu(gx, x));
			slots.set_yidx(p, Indices::get_nearest_cell_index_cpu(gy, y));
			slots.set_zidx(p, Indices::get_nearest_cell_index_cpu(gz, z));
			slots.set_q(p, p % 7 == 0 ? 0 : 1 + p % 5);
			const double cos_th {2.0 * rng.next_double() - 1.0};
			const double sin_th {std::sqrt(1.0 - cos_th * cos_th)};
			const double az {2.0 * Constants::pi * rng.next_double()};
			const double v {speed * (0.2 + rng.next_double())};
			slots.set_vX(p, v * sin_th * std::cos(az));
			slots.set_vY(p, v * sin_th * std::sin(az));
			slots.set_vZ(p, v * cos_th);
		}
		return slots;
	}

	double rel(const double a, const double b, const double scale)
	{
		return std::abs(a - b) / std::max({std::abs(a), std::abs(b), scale});
	}

	/**
	* Runs both and compares them.
	*
	* @param bkg Background with walls.
	* @param sheath Its sheath.
	* @param slots The particles.
	* @param dt Time step (s).
	* @param substep Sub-step near the wall (else one push with the field).
	* @return Number of particles sub-stepped, and the most sub-steps.
	*/
	std::pair<int, int> push_matches_cpu(Background::Background& bkg,
		const Sheath::Sheath& sheath, const Slots::Slots& slots,
		const double dt, const bool substep)
	{
		// CPU, as boris_wrapper does it
		Slots::Slots cpu {slots};
		const Options::Options opts {};
		std::vector<int> nsub(slots.N(), 1);
		if (substep)
		{
			for (int i {}; i < cpu.N(); ++i)
			{
				const double q_m {cpu.q()[i] * -Constants::charge_e
					/ cpu.mass()};
				const double speed {std::hypot(cpu.vX()[i], cpu.vY()[i],
					cpu.vZ()[i])};
				nsub[i] = sheath.substeps(cpu.zone()[i], cpu.tidx()[i],
					cpu.xidx()[i], cpu.yidx()[i], cpu.zidx()[i],
					{cpu.x()[i], cpu.y()[i], cpu.z()[i]}, speed, q_m, dt);
			}
			Boris::update_velocity_cpu(cpu, bkg, opts, dt, &sheath, &nsub);
			Boris::substep_cpu(cpu, bkg, dt, sheath, nsub);
		}
		else Boris::update_velocity_cpu(cpu, bkg, opts, dt, &sheath);

		// GPU
		Background::BackgroundDevice bkg_d {bkg.to_device()};
		Slots::Slots gpu {slots};
		Slots::SlotsDevice slots_d {gpu.to_device()};
		Sheath::SheathDevice sheath_d {sheath.to_device(0)};
		EXPECT_EQ(sheath_d.device_id, 0);
		Boris::update_velocity_gpu(slots_d, bkg_d, dt, &sheath_d, substep);
		EXPECT_EQ(cudaDeviceSynchronize(), cudaSuccess);
		gpu = gpu.to_host(slots_d);
		std::vector<int> nsub_gpu(slots.N(), 1);
		if (substep)
		{
			cudaMemcpy(nsub_gpu.data(), slots_d.nsub, slots.N()
				* sizeof(int), cudaMemcpyDeviceToHost);
		}
		Sheath::free_sheath(sheath_d);
		Slots::free_slots(slots_d);
		Background::free_bkg(bkg_d);

		int bad {};
		int stepped {};
		int most {};
		for (int i {}; i < slots.N(); ++i)
		{
			const double v {std::max(std::hypot(cpu.vX()[i], cpu.vY()[i],
				cpu.vZ()[i]), std::hypot(slots.vX()[i], slots.vY()[i],
				slots.vZ()[i]))};
			bool ok {nsub_gpu[i] == nsub[i]};
			ok &= rel(gpu.vX()[i], cpu.vX()[i], v) < 1e-8;
			ok &= rel(gpu.vY()[i], cpu.vY()[i], v) < 1e-8;
			ok &= rel(gpu.vZ()[i], cpu.vZ()[i], v) < 1e-8;
			ok &= rel(gpu.x()[i], cpu.x()[i], 1e-6) < 1e-8;
			ok &= rel(gpu.y()[i], cpu.y()[i], 1e-6) < 1e-8;
			ok &= rel(gpu.z()[i], cpu.z()[i], 1e-6) < 1e-8;
			ok &= gpu.t()[i] == cpu.t()[i];
			if (!ok && bad++ < 5)
			{
				ADD_FAILURE() << "particle " << i << " (q " << slots.q()[i]
					<< ", nsub CPU " << nsub[i] << " GPU " << nsub_gpu[i]
					<< "): GPU v (" << gpu.vX()[i] << ", " << gpu.vY()[i]
					<< ", " << gpu.vZ()[i] << ") x (" << gpu.x()[i] << ", "
					<< gpu.y()[i] << ", " << gpu.z()[i] << "), CPU v ("
					<< cpu.vX()[i] << ", " << cpu.vY()[i] << ", "
					<< cpu.vZ()[i] << ") x (" << cpu.x()[i] << ", "
					<< cpu.y()[i] << ", " << cpu.z()[i] << ")";
			}
			if (nsub[i] > 1) ++stepped;
			most = std::max(most, nsub[i]);
		}
		EXPECT_EQ(bad, 0);
		return {stepped, most};
	}

	Options::Options slab_wall_options()
	{
		Options::Options opts {};
		opts.set_bkg_source("slab");
		opts.set_slab_ne(1e19);
		opts.set_slab_te(20.0);
		opts.set_slab_ti(20.0);
		opts.set_slab_bmag(2.0);
		opts.set_slab_tmax(1e-5);
		opts.set_slab_wall("ymin");
		opts.set_slab_b_angle(5.0);
		return opts;
	}
}

// Slab wall, Brooks: the field alone, then the sub-stepped step
TEST(SheathCUDA, SlabBrooksMatchesCPU)
{
	Background::Background bkg {Background::read_bkg(slab_wall_options())};
	const Sputtering::Surface surface(bkg);
	const Sheath::Sheath sheath(surface, bkg, Sheath::Model::brooks);
	const Slots::Slots slots {near_wall(bkg, surface, 3000, 0.1, 2e4)};

	push_matches_cpu(bkg, sheath, slots, 1e-8, false);
	const auto [stepped, most] {push_matches_cpu(bkg, sheath, slots, 1e-8,
		true)};
	EXPECT_GT(stepped, 500);
	EXPECT_GT(most, 100);
}

// HERMES targets, Borodkina, on a wedge with field-aligned zones
TEST(SheathCUDA, HermesBorodkinaMatchesCPU)
{
	const auto dir {std::filesystem::temp_directory_path()
		/ ("flan_test_sheath_cuda_" + std::to_string(getpid()))};
	Background::Background bkg {MultiZone::hermes_background(dir)};
	std::filesystem::remove_all(dir);
	const Sputtering::Surface surface(bkg);
	const Sheath::Sheath sheath(surface, bkg, Sheath::Model::borodkina);
	const Slots::Slots slots {near_wall(bkg, surface, 3000, 0.02, 2e4)};

	push_matches_cpu(bkg, sheath, slots, 1e-8, false);
	const auto [stepped, most] {push_matches_cpu(bkg, sheath, slots, 1e-8,
		true)};
	EXPECT_GT(stepped, 500);
	EXPECT_GT(most, 10);
}
