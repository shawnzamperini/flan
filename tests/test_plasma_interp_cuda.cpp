/*
* The plasma interpolation of the GPU kernels (Interpolate::interp_plasma
* and interp_plasma_vec3, used for ne, Te, Ti and the flow in the
* collision, ionization and birth kernels, and for E in the push) against
* the CPU's (Background::interp_ne, interp_te, interp_ti, interp_4d with
* aligned = true and interp_u), at points spread over every zone of the
* HERMES and SOLEDGE3X test backgrounds.
*
* On the HERMES grid every zone has a zShift table, so both interpolate
* along the field lines, with the shared Utilities::interp_aligned_flat;
* the flow is also rotated to each point's toroidal angle. On the
* SOLEDGE3X mesh there is no zShift: the GPU uses its stencils and the CPU
* Utilities::interp_vec4d, which pick the same cells (ghost cells at the
* zone faces) but differ in one place: the GPU doesn't clamp the time to
* the background's frames. The points are inside them here.
*/

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <gtest/gtest.h>
#include <string>
#include <unistd.h>
#include <vector>

#include "background.h"
#include "background_device.h"
#include "multizone_cases.h"
#include "plasma_interp_kernel.h"

namespace
{
	void interp_matches_cpu(const Background::Background& bkg_in)
	{
		Background::Background bkg {bkg_in};
		const Slots::Slots slots {MultiZone::spread_particles(bkg, 500, 3)};
		std::vector<PlasmaInterpTest::Point> points {};
		for (int i {}; i < slots.N(); ++i)
		{
			points.push_back({slots.zone()[i], slots.tidx()[i],
				slots.xidx()[i], slots.yidx()[i], slots.zidx()[i],
				slots.t()[i], slots.x()[i], slots.y()[i], slots.z()[i]});
		}

		Background::BackgroundDevice bkg_d {bkg.to_device()};
		const std::vector<double> gpu {PlasmaInterpTest::interp_on_gpu(bkg_d,
			points)};
		Background::free_bkg(bkg_d);

		constexpr int nv {PlasmaInterpTest::num_values};
		std::vector<double> cpu (gpu.size());
		for (std::size_t i {}; i < points.size(); ++i)
		{
			const auto& p {points[i]};
			double* c {cpu.data() + i * nv};
			c[0] = bkg.interp_ne(p.zone, p.t, p.x, p.y, p.z);
			c[1] = bkg.interp_te(p.zone, p.t, p.x, p.y, p.z);
			c[2] = bkg.interp_ti(p.zone, p.t, p.x, p.y, p.z);
			c[3] = bkg.interp_4d(bkg.get_emag(), p.zone, p.t, p.x, p.y, p.z,
				true);
			const auto u {bkg.interp_u(p.zone, p.t, p.x, p.y, p.z)};
			c[4] = u[0];
			c[5] = u[1];
			c[6] = u[2];
		}

		// Each quantity to rounding, relative to its largest value (the
		// three flow components together)
		const char* names[nv] {"ne", "Te", "Ti", "|E|", "uX", "uY", "uZ"};
		std::vector<double> scale (nv);
		for (std::size_t i {}; i < points.size(); ++i)
		{
			for (int v {}; v < nv; ++v)
			{
				const int s {v < 4 ? v : 4};
				scale[s] = std::max(scale[s], std::abs(cpu[i * nv + v]));
			}
		}
		for (int v {}; v < 5; ++v) ASSERT_GT(scale[v], 0.0) << names[v];
		int bad {};
		for (std::size_t i {}; i < points.size(); ++i)
		{
			for (int v {}; v < nv; ++v)
			{
				const double tol {1e-11 * scale[v < 4 ? v : 4]};
				const double d {std::abs(gpu[i * nv + v] - cpu[i * nv + v])};
				if (!(d <= tol) && bad++ < 5)
				{
					ADD_FAILURE() << names[v] << " at point " << i
						<< " (zone " << points[i].zone << "): GPU "
						<< gpu[i * nv + v] << ", CPU " << cpu[i * nv + v];
				}
			}
		}
		EXPECT_EQ(bad, 0);
	}
}

TEST(PlasmaInterpCUDA, MatchesCPUOnHermesGrid)
{
	const auto dir {std::filesystem::temp_directory_path()
		/ ("flan_test_plasma_interp_hermes_" + std::to_string(getpid()))};
	const Background::Background bkg {MultiZone::hermes_background(dir)};
	std::filesystem::remove_all(dir);
	interp_matches_cpu(bkg);
}

TEST(PlasmaInterpCUDA, MatchesCPUOnS3XGrid)
{
	const auto dir {std::filesystem::temp_directory_path()
		/ ("flan_test_plasma_interp_s3x_" + std::to_string(getpid()))};
	const Background::Background bkg {MultiZone::s3x_background(dir)};
	std::filesystem::remove_all(dir);
	interp_matches_cpu(bkg);
}
