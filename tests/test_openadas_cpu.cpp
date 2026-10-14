#include <array>
#include <gtest/gtest.h>
#include <omp.h>
#include <vector>

#include "background.h"
#include "openadas.h"
#include "pcg32.h"
#include "slots.h"

#include "openadas_test_utils.h"

using namespace OpenADASTest;

namespace
{
	std::vector<pcg32> get_rngs(uint64_t base_seed = 4)
	{
		std::vector<pcg32> rngs;
		int nthreads = omp_get_max_threads();
		rngs.reserve(nthreads);

		for (int tid = 0; tid < nthreads; ++tid)
		{
			uint64_t seed = base_seed + 0;  // Assigned MPI rank = 0
			uint64_t stream = tid;          // unique per thread
			rngs.push_back(pcg32(seed, stream));
		}

		return rngs;
	}

	// Loads the fake SCD and ACD files written to root
	std::pair<OpenADAS::OpenADAS, OpenADAS::OpenADAS> load_fake(
		const fs::path& root)
	{
		return {OpenADAS::OpenADAS(root.string(), fake_year, fake_Z, "scd"),
			OpenADAS::OpenADAS(root.string(), fake_year, fake_Z, "acd")};
	}

	// Expects rate to equal the analytic value to within rounding of the
	// 12 decimal places written to the fake file.
	void expect_rate(double rate, double expected)
	{
		EXPECT_NEAR(rate, expected, 1e-9 * expected);
	}
}


// -------------------------------------------------------------------------
// Rate coefficient lookup
// -------------------------------------------------------------------------

// Every tabulated point for every charge block must come back exactly. This
// checks the file parsing: unit conversions, charge block order, and that
// the table is stored as (charge, Te, ne).
TEST(RateCoeff, MatchesTableAtGridPoints)
{
	const fs::path root {write_fake_adas_dir("grid", scd_rate, acd_rate)};
	auto [oa_ioniz, oa_recomb] = load_fake(root);

	for (int c {}; c < fake_Z; ++c)
	for (double te : fake_te)
	for (double ne : fake_ne)
	{
		SCOPED_TRACE("c = " + std::to_string(c) + ", te = "
			+ std::to_string(te) + ", ne = " + std::to_string(ne));
		expect_rate(oa_ioniz.get_rate_coeff(c, ne, te), scd_rate(c, te, ne));
		expect_rate(oa_recomb.get_rate_coeff(c, ne, te), acd_rate(c, te, ne));
	}
}


// The rates are bilinear in (Te, ne), so bilinear interpolation between grid
// points must also be exact. Points are chosen in different cells of the
// non-square table, including the edge cells.
TEST(RateCoeff, InterpolatesBetweenGridPoints)
{
	const fs::path root {write_fake_adas_dir("interp", scd_rate, acd_rate)};
	auto [oa_ioniz, oa_recomb] = load_fake(root);

	const std::vector<std::pair<double, double>> points {
		{3.0, 4e19}, {300.0, 2e18}, {1.5, 2e17}, {900.0, 8e21},
		{55.0, 7e20}, {10.0, 2e19}};

	for (int c {}; c < fake_Z; ++c)
	for (auto [te, ne] : points)
	{
		SCOPED_TRACE("c = " + std::to_string(c) + ", te = "
			+ std::to_string(te) + ", ne = " + std::to_string(ne));
		expect_rate(oa_ioniz.get_rate_coeff(c, ne, te), scd_rate(c, te, ne));
		expect_rate(oa_recomb.get_rate_coeff(c, ne, te), acd_rate(c, te, ne));
	}
}


// Outside the table, ne and Te are clamped to the nearest edge rather than
// extrapolated.
TEST(RateCoeff, ClampsOutsideTable)
{
	const fs::path root {write_fake_adas_dir("clamp", scd_rate, acd_rate)};
	auto [oa_ioniz, oa_recomb] = load_fake(root);

	const double te_min {fake_te.front()};
	const double te_max {fake_te.back()};
	const double ne_min {fake_ne.front()};
	const double ne_max {fake_ne.back()};

	for (int c {}; c < fake_Z; ++c)
	{
		SCOPED_TRACE("c = " + std::to_string(c));

		// Te below/above the table
		expect_rate(oa_ioniz.get_rate_coeff(c, 4e19, 0.1),
			scd_rate(c, te_min, 4e19));
		expect_rate(oa_ioniz.get_rate_coeff(c, 4e19, 1e5),
			scd_rate(c, te_max, 4e19));

		// ne below/above the table
		expect_rate(oa_ioniz.get_rate_coeff(c, 1e14, 30.0),
			scd_rate(c, 30.0, ne_min));
		expect_rate(oa_ioniz.get_rate_coeff(c, 1e25, 30.0),
			scd_rate(c, 30.0, ne_max));

		// Both outside at once (corners)
		expect_rate(oa_recomb.get_rate_coeff(c, 1e14, 0.1),
			acd_rate(c, te_min, ne_min));
		expect_rate(oa_recomb.get_rate_coeff(c, 1e25, 1e5),
			acd_rate(c, te_max, ne_max));
	}
}


// -------------------------------------------------------------------------
// Charge-state update (deterministic cases)
// -------------------------------------------------------------------------

// With certain ionization and no recombination, every particle that can
// ionize goes up by one. The fully stripped particles can't ionize.
TEST(IonizRecomb, IonizationOnly_IncreasesCharge)
{
	const fs::path root {write_fake_adas_dir("ioniz_only", huge_rate,
		tiny_rate)};
	auto [oa_ioniz, oa_recomb] = load_fake(root);
	Background::Background bkg {make_uniform_bkg(1e19, 10.0)};
	std::vector<pcg32> rngs {get_rngs()};

	constexpr int per_charge {16};
	Slots::Slots slots {make_slots(bkg, per_charge)};

	int ioniz_warnings {};
	int recomb_warnings {};
	OpenADAS::ioniz_recomb(slots, bkg, oa_ioniz, oa_recomb, 1e-6,
		ioniz_warnings, recomb_warnings, rngs);

	for (int i {}; i < slots.N(); ++i)
	{
		const int q0 {i / per_charge};
		EXPECT_EQ(slots.q()[i], std::min(q0 + 1, fake_Z)) << "q0 = " << q0;
	}

	// Probability > 1 for every particle that can ionize
	EXPECT_EQ(ioniz_warnings, fake_Z * per_charge);
	EXPECT_EQ(recomb_warnings, 0);
}


// With certain recombination and no ionization, every particle that can
// recombine goes down by one. Neutrals can't recombine.
TEST(IonizRecomb, RecombinationOnly_DecreasesCharge)
{
	const fs::path root {write_fake_adas_dir("recomb_only", tiny_rate,
		huge_rate)};
	auto [oa_ioniz, oa_recomb] = load_fake(root);
	Background::Background bkg {make_uniform_bkg(1e19, 10.0)};
	std::vector<pcg32> rngs {get_rngs()};

	constexpr int per_charge {16};
	Slots::Slots slots {make_slots(bkg, per_charge)};

	int ioniz_warnings {};
	int recomb_warnings {};
	OpenADAS::ioniz_recomb(slots, bkg, oa_ioniz, oa_recomb, 1e-6,
		ioniz_warnings, recomb_warnings, rngs);

	for (int i {}; i < slots.N(); ++i)
	{
		const int q0 {i / per_charge};
		EXPECT_EQ(slots.q()[i], std::max(q0 - 1, 0)) << "q0 = " << q0;
	}

	EXPECT_EQ(ioniz_warnings, 0);
	EXPECT_EQ(recomb_warnings, fake_Z * per_charge);
}


// With both certain, they cancel for partially ionized particles. Neutrals
// can only ionize and fully stripped particles can only recombine.
TEST(IonizRecomb, Cancellation_NoNetChange)
{
	const fs::path root {write_fake_adas_dir("cancel", huge_rate,
		huge_rate)};
	auto [oa_ioniz, oa_recomb] = load_fake(root);
	Background::Background bkg {make_uniform_bkg(1e19, 10.0)};
	std::vector<pcg32> rngs {get_rngs()};

	constexpr int per_charge {16};
	Slots::Slots slots {make_slots(bkg, per_charge)};

	int ioniz_warnings {};
	int recomb_warnings {};
	OpenADAS::ioniz_recomb(slots, bkg, oa_ioniz, oa_recomb, 1e-6,
		ioniz_warnings, recomb_warnings, rngs);

	for (int i {}; i < slots.N(); ++i)
	{
		const int q0 {i / per_charge};
		int expected {q0};
		if (q0 == 0) expected = 1;
		else if (q0 == fake_Z) expected = fake_Z - 1;
		EXPECT_EQ(slots.q()[i], expected) << "q0 = " << q0;
	}

	EXPECT_EQ(ioniz_warnings, fake_Z * per_charge);
	EXPECT_EQ(recomb_warnings, fake_Z * per_charge);
}


// Dead particles are skipped entirely, even with certain ionization and
// recombination.
TEST(IonizRecomb, DeadParticles_Unchanged)
{
	const fs::path root {write_fake_adas_dir("dead", huge_rate, huge_rate)};
	auto [oa_ioniz, oa_recomb] = load_fake(root);
	Background::Background bkg {make_uniform_bkg(1e19, 10.0)};
	std::vector<pcg32> rngs {get_rngs()};

	constexpr int per_charge {16};
	Slots::Slots slots {make_slots(bkg, per_charge)};
	for (int i {}; i < slots.N(); ++i) slots.set_state(i, 1);

	int ioniz_warnings {};
	int recomb_warnings {};
	OpenADAS::ioniz_recomb(slots, bkg, oa_ioniz, oa_recomb, 1e-6,
		ioniz_warnings, recomb_warnings, rngs);

	for (int i {}; i < slots.N(); ++i)
		EXPECT_EQ(slots.q()[i], i / per_charge);

	EXPECT_EQ(ioniz_warnings, 0);
	EXPECT_EQ(recomb_warnings, 0);
}


// -------------------------------------------------------------------------
// Charge-state update (statistical)
// -------------------------------------------------------------------------

// Many particles in each charge state take one step at a given (ne, Te).
// The fraction that ionize and recombine must match the analytic rates.
// This checks that the right charge block is used for each process
// (q for ionization, q-1 for recombination) and that the lookup uses the
// local ne and Te. The two plasma conditions are off-grid points where
// ionization dominates and where recombination is significant.
TEST(IonizRecomb, FractionsMatchRates)
{
	const fs::path root {write_fake_adas_dir("fractions", scd_rate,
		acd_rate)};
	auto [oa_ioniz, oa_recomb] = load_fake(root);

	// ne (m^-3), Te (eV), dt (s)
	const std::vector<std::array<double, 3>> cases {
		{4e19, 30.0, 1e-8}, {5e20, 3.0, 5e-10}};

	constexpr int per_charge {100000};
	uint64_t seed {11};
	for (auto [ne, te, dt] : cases)
	{
		Background::Background bkg {make_uniform_bkg(ne, te)};
		std::vector<pcg32> rngs {get_rngs(seed++)};
		Slots::Slots slots {make_slots(bkg, per_charge)};

		int ioniz_warnings {};
		int recomb_warnings {};
		OpenADAS::ioniz_recomb(slots, bkg, oa_ioniz, oa_recomb, dt,
			ioniz_warnings, recomb_warnings, rngs);

		expect_fractions_match(slots, per_charge, ne, te, dt);
		EXPECT_EQ(ioniz_warnings, 0);
		EXPECT_EQ(recomb_warnings, 0);
	}
}
