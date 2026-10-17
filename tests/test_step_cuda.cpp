/*
* Tests for ImpurityTransport::step_gpu: that it moves particles, and that
* it kills and counts particles that have been followed for longer than
* imp_max_time, the same as step_cpu.
*/

#include <gtest/gtest.h>
#include <cuda_runtime.h>

#include "openadas_device.h"
#include "options.h"
#include "options_device.h"
#include "slots.h"
#include "slots_device.h"

#include "impurity_transport.cuh"
#include "openadas.cuh"


namespace
{
	// Copies slots and opts to the GPU, calls step_gpu nsteps times with
	// time step dt, copies the slots back and returns how many particles
	// were counted as exceeding imp_max_time. The counter is the one used
	// in a simulation, from OpenADAS::alloc_warnings.
	int step_on_device(Slots::Slots& slots, Options::Options& opts,
		const double dt, const int nsteps)
	{
		Slots::SlotsDevice slots_d {slots.to_device()};
		Options::OptionsDevice* opts_d {opts.to_device()};
		OpenADAS::IonizRecombWarningsDevice warnings_d {
			OpenADAS::alloc_warnings(0)};

		for (int n {}; n < nsteps; ++n)
		{
			ImpurityTransport::step_gpu(slots_d, dt, opts_d,
				warnings_d.max_time);
		}
		EXPECT_EQ(cudaDeviceSynchronize(), cudaSuccess);

		int ioniz_warnings {};
		int recomb_warnings {};
		int max_time_ctr {};
		OpenADAS::add_warnings_to_host(warnings_d, ioniz_warnings,
			recomb_warnings, max_time_ctr);

		// step_gpu never touches the ionization/recombination counters
		EXPECT_EQ(ioniz_warnings, 0);
		EXPECT_EQ(recomb_warnings, 0);

		slots = slots.to_host(slots_d);
		Slots::free_slots(slots_d);
		Options::free_opts(opts_d);
		OpenADAS::free_warnings(warnings_d);

		return max_time_ctr;
	}
}


// One step with particles at different times since launch. The ones that
// end up past imp_max_time are killed and counted, the rest are moved, and
// slots that were already dead are left alone.
TEST(StepCUDA, KillsParticlesPastMaxTime)
{
	const double dt {1e-8};
	const int N {1000};
	const int N_alive {900};

	// After one step, alive slot i has pt = (i % 10 + 1) * dt, so the ones
	// with i % 10 >= 5 are past 5.5 * dt
	Options::Options opts {};
	opts.set_imp_max_time(5.5 * dt);

	Slots::Slots slots(N);
	for (int i {}; i < N; ++i)
	{
		slots.set_x(i, 1.0);
		slots.set_y(i, 2.0);
		slots.set_z(i, 3.0);
		slots.set_vx(i, 100.0);
		slots.set_vy(i, 200.0);
		slots.set_vz(i, 300.0);

		// The last 100 slots are already dead, with a pt that would be past
		// imp_max_time if they were stepped
		if (i < N_alive)
		{
			slots.set_pt(i, (i % 10) * dt);
		}
		else
		{
			slots.set_pt(i, 1.0);
			slots.set_state(i, 1);
		}
	}

	int max_time_ctr {step_on_device(slots, opts, dt, 1)};

	// Half of the alive slots are past imp_max_time
	EXPECT_EQ(max_time_ctr, N_alive / 2);

	for (int i {}; i < N; ++i)
	{
		if (i >= N_alive)
		{
			// Already dead, so not stepped
			EXPECT_EQ(slots.state()[i], 1) << "i = " << i;
			EXPECT_DOUBLE_EQ(slots.pt()[i], 1.0) << "i = " << i;
			EXPECT_DOUBLE_EQ(slots.x()[i], 1.0) << "i = " << i;
			continue;
		}

		// Every alive particle is stepped, including the ones killed for
		// exceeding imp_max_time on this step
		EXPECT_DOUBLE_EQ(slots.pt()[i], (i % 10 + 1) * dt) << "i = " << i;
		EXPECT_DOUBLE_EQ(slots.t()[i], dt) << "i = " << i;
		EXPECT_DOUBLE_EQ(slots.x()[i], 1.0 + 100.0 * dt) << "i = " << i;
		EXPECT_DOUBLE_EQ(slots.y()[i], 2.0 + 200.0 * dt) << "i = " << i;
		EXPECT_DOUBLE_EQ(slots.z()[i], 3.0 + 300.0 * dt) << "i = " << i;

		const int expected_state {(i % 10 >= 5) ? 1 : 0};
		EXPECT_EQ(slots.state()[i], expected_state) << "i = " << i;
	}
}


// Several steps. Every particle dies on the third step, and is counted
// exactly once even though step_gpu is called after it's dead.
TEST(StepCUDA, CountsEachParticleOnce)
{
	const double dt {1e-8};
	const int N {300};

	Options::Options opts {};
	opts.set_imp_max_time(2.5 * dt);

	Slots::Slots slots(N);

	int max_time_ctr {step_on_device(slots, opts, dt, 5)};

	EXPECT_EQ(max_time_ctr, N);
	for (int i {}; i < N; ++i)
	{
		EXPECT_EQ(slots.state()[i], 1) << "i = " << i;

		// Stepped 3 times, then skipped because it's dead
		EXPECT_DOUBLE_EQ(slots.pt()[i], 3 * dt) << "i = " << i;
	}
}


// With the default imp_max_time (1e10 s), nothing is killed
TEST(StepCUDA, DefaultMaxTimeKillsNothing)
{
	const double dt {1e-8};
	const int N {300};

	Options::Options opts {};
	Slots::Slots slots(N);

	int max_time_ctr {step_on_device(slots, opts, dt, 5)};

	EXPECT_EQ(max_time_ctr, 0);
	for (int i {}; i < N; ++i)
	{
		EXPECT_EQ(slots.state()[i], 0) << "i = " << i;
		EXPECT_DOUBLE_EQ(slots.pt()[i], 5 * dt) << "i = " << i;
	}
}
