/**
* @file openadas_test_utils.h
* @brief Helpers shared by the CPU and CUDA OpenADAS tests.
*
* The fake ADAS files written here mimic the real ADF11 layout: the header
* lists ndens before ntemp, the density grid (log10 cm^-3) comes before the
* temperature grid (log10 eV), charge blocks are labeled starting at Z1 = 1
* (so charge_low = 1), and each block holds one row of ndens values per
* temperature, wrapped at 8 values per line.
*
* The rates are generated from analytic functions that are bilinear in
* (Te, ne), so bilinear interpolation of the table reproduces them exactly at
* any point inside the grid. The functions depend differently on Te and ne,
* and on the charge block, so swapping ne/Te or using the wrong charge block
* gives a clearly different answer.
*/
#pragma once

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <functional>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "background.h"
#include "options.h"
#include "read_bkg.h"
#include "slots.h"

namespace OpenADASTest
{
	namespace fs = std::filesystem;

	// Rate coefficient (m^3/s) as a function of the charge block (0-based),
	// Te (eV) and ne (m^-3).
	using RateFunc = std::function<double(int, double, double)>;

	// Fake impurity: lithium, so there are three charge blocks (Z1 = 1-3)
	constexpr int fake_Z {3};
	constexpr int fake_year {99};

	// Non-square table (10 densities x 4 temperatures). Ten densities means
	// each temperature row wraps onto a second line, like the real files.
	const std::vector<double> fake_ne {1e17, 1e18, 5e18, 1e19, 3e19, 1e20,
		3e20, 1e21, 3e21, 1e22};
	const std::vector<double> fake_te {1.0, 10.0, 100.0, 1000.0};

	// Ionization (SCD) rate for charge block c, i.e., c -> c+1
	inline double scd_rate(int c, double te, double ne)
	{
		return 1e-14 * (c + 1) * (1.0 + te / 10.0) * (1.0 + ne / 1e19);
	}

	// Recombination (ACD) rate for charge block c, i.e., c+1 -> c
	inline double acd_rate(int c, double te, double ne)
	{
		return 1e-14 * (c + 2) * (2.0 + te / 50.0) * (1.0 + 3.0 * ne / 1e20);
	}

	// Writes values in the ADF11 style, 8 per line
	inline void write_wrapped(std::ofstream& f, const std::vector<double>& v)
	{
		char buf[32];
		for (std::size_t i {}; i < v.size(); ++i)
		{
			std::snprintf(buf, sizeof(buf), " %.12f", v[i]);
			f << buf;
			if ((i + 1) % 8 == 0 || i + 1 == v.size()) f << '\n';
		}
	}

	// Writes a fake ADF11 file for the fake impurity using the rate function
	inline void write_fake_adas(const fs::path& path, const RateFunc& rate)
	{
		std::ofstream f(path);

		const int ndens {static_cast<int>(fake_ne.size())};
		const int ntemp {static_cast<int>(fake_te.size())};

		// atomic number, ndens, ntemp, charge_low, charge_high
		f << "    " << fake_Z << "   " << ndens << "    " << ntemp
			<< "    1    " << fake_Z << "     /LITHIUM            /TEST\n";
		f << " ---------------------------------------------------------------\n";

		// Density grid first (log10 cm^-3), then temperature (log10 eV)
		std::vector<double> log_ne {};
		for (double ne : fake_ne) log_ne.push_back(std::log10(ne * 1e-6));
		write_wrapped(f, log_ne);

		std::vector<double> log_te {};
		for (double te : fake_te) log_te.push_back(std::log10(te));
		write_wrapped(f, log_te);

		// One block per charge, each with one row of densities per Te
		for (int c {}; c < fake_Z; ++c)
		{
			f << " --------------------/ IPRT= 1  / IGRD= 1  /--------/ Z1= "
				<< c + 1 << "   / DATE= 01.01.26\n";
			for (double te : fake_te)
			{
				std::vector<double> row {};
				for (double ne : fake_ne)
					row.push_back(std::log10(rate(c, te, ne) * 1e6));
				write_wrapped(f, row);
			}
		}
	}

	// Writes fake SCD and ACD files into a fresh directory named after the
	// test, laid out the way OpenADAS::get_openadas_path expects. Returns
	// the root directory to pass to the OpenADAS constructor.
	inline fs::path write_fake_adas_dir(const std::string& name,
		const RateFunc& scd, const RateFunc& acd)
	{
		const fs::path root {fs::temp_directory_path()
			/ ("flan_fake_adas_" + name)};
		const std::string year {std::to_string(fake_year)};
		fs::remove_all(root);
		fs::create_directories(root / ("scd" + year));
		fs::create_directories(root / ("acd" + year));
		write_fake_adas(root / ("scd" + year) / ("scd" + year + "_li.dat"),
			scd);
		write_fake_adas(root / ("acd" + year) / ("acd" + year + "_li.dat"),
			acd);
		return root;
	}

	// Constant rate functions for the deterministic tests. A huge rate makes
	// the probability > 1 (always happens), a tiny one makes it ~0.
	inline double huge_rate(int, double, double) { return 1e-4; }
	inline double tiny_rate(int, double, double) { return 1e-90; }

	// Test background with ne and Te set to uniform values
	inline Background::Background make_uniform_bkg(double ne, double te)
	{
		Options::Options opts {};
		opts.set_bkg_source("test");
		opts.set_test_opt("gyrate");
		Background::Background bkg {Background::read_bkg(opts)};

		auto& ne_v = bkg.get_ne().get_data();
		auto& te_v = bkg.get_te().get_data();
		std::fill(ne_v.begin(), ne_v.end(), ne);
		std::fill(te_v.begin(), te_v.end(), te);
		return bkg;
	}

	// Slots with per_charge alive particles in each charge state 0..fake_Z,
	// ordered by charge, all inside the test grid.
	inline Slots::Slots make_slots(const Background::Background& bkg,
		int per_charge)
	{
		Slots::Slots slots((fake_Z + 1) * per_charge);
		slots.set_Z(fake_Z);
		for (int i {}; i < slots.N(); ++i)
		{
			slots.set_q(i, i / per_charge);
			slots.set_state(i, 0);
			slots.set_t(i, bkg.get_times()[0]);
			slots.set_x(i, 0.01);
			slots.set_y(i, 0.0);
			slots.set_z(i, 0.0);
			slots.set_tidx(i, 0);
			slots.set_xidx(i, 10);
			slots.set_yidx(i, 10);
			slots.set_zidx(i, 4);
		}
		return slots;
	}

	// Checks the outcome of one ionization/recombination step on slots made
	// by make_slots against the analytic rates. A particle with charge q
	// ionizes with P_iz = S_q ne dt (q < Z) and recombines with
	// P_rec = a_(q-1) ne dt (q > 0). Both can happen and cancel, so
	// P(q+1) = P_iz (1 - P_rec) and P(q-1) = P_rec (1 - P_iz). Observed
	// fractions must be within 5 binomial standard errors.
	inline void expect_fractions_match(const Slots::Slots& slots,
		int per_charge, double ne, double te, double dt)
	{
		for (int q {}; q <= fake_Z; ++q)
		{
			const double p_iz {q < fake_Z ? scd_rate(q, te, ne) * ne * dt
				: 0.0};
			const double p_rec {q > 0 ? acd_rate(q - 1, te, ne) * ne * dt
				: 0.0};
			const double exp_up {p_iz * (1.0 - p_rec)};
			const double exp_down {p_rec * (1.0 - p_iz)};

			int up {};
			int down {};
			int other {};
			for (int i {q * per_charge}; i < (q + 1) * per_charge; ++i)
			{
				const int dq {slots.q()[i] - q};
				up += (dq == 1);
				down += (dq == -1);
				other += (dq != 0 && dq != 1 && dq != -1);
			}

			const double n {static_cast<double>(per_charge)};
			const double tol_up {5.0 * std::sqrt(exp_up * (1.0 - exp_up) / n)};
			const double tol_down {5.0 * std::sqrt(exp_down * (1.0 - exp_down)
				/ n)};

			EXPECT_NEAR(up / n, exp_up, tol_up) << "q = " << q
				<< ", ne = " << ne << ", Te = " << te;
			EXPECT_NEAR(down / n, exp_down, tol_down) << "q = " << q
				<< ", ne = " << ne << ", Te = " << te;

			// Charge can only change by one per step
			EXPECT_EQ(other, 0) << "q = " << q;
		}
	}
}  // namespace OpenADASTest
