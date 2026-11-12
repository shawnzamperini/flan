/*
* Tests for reading and interpolating sputtering yield tables
* (Sputtering::YieldTable), on small tables written by the test itself and
* on the table in data/sputtering/.
*/

#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <gtest/gtest.h>
#include <string>
#include <unistd.h>
#include <vector>

#include "hdf5.h"
#include "hdf5_utils.h"
#include "sputtering.h"


namespace
{
	// Writes a dataset of doubles with the given shape (empty for a
	// scalar), creating groups along the way
	void write_doubles(hid_t file, const std::string& name,
		const std::vector<hsize_t>& shape, const std::vector<double>& data)
	{
		const hid_t lcpl {H5Pcreate(H5P_LINK_CREATE)};
		H5Pset_create_intermediate_group(lcpl, 1);
		const hid_t space {shape.empty() ? H5Screate(H5S_SCALAR)
			: H5Screate_simple(static_cast<int>(shape.size()), shape.data(),
			nullptr)};
		const hid_t dset {H5Dcreate2(file, name.c_str(), H5T_NATIVE_DOUBLE,
			space, lcpl, H5P_DEFAULT, H5P_DEFAULT)};
		H5Dwrite(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
			data.data());
		H5Dclose(dset);
		H5Sclose(space);
		H5Pclose(lcpl);
	}

	// The test table's yield above its lowest energy. Linear in log10(E)
	// and angle, so bilinear interpolation gives it exactly.
	double linear_yield(double energy, double angle)
	{
		return 0.01 + 0.02 * std::log10(energy) + 0.001 * angle;
	}

	// What goes in a test table. The defaults are a good table: 4 energies
	// and 4 angles, with the lowest energy below threshold (Y = 0) and
	// linear_yield above it.
	struct Table
	{
		std::vector<double> energy {10.0, 100.0, 1000.0, 10000.0};
		std::vector<double> angle {0.0, 30.0, 60.0, 89.0};
		std::vector<double> yield {};
		std::vector<hsize_t> yield_shape {};	// Default (nE, nA)
		bool target_Es {true};					// Write /target/Es

		Table()
		{
			for (const double e : energy)
			{
				for (const double a : angle)
				{
					yield.push_back(e == 10.0 ? 0.0 : linear_yield(e, a));
				}
			}
		}
	};

	// Test tables are written here and removed when the test ends
	class SputteringTable : public testing::Test
	{
	protected:
		std::string m_path {(std::filesystem::temp_directory_path()
			/ ("flan_test_sputtering_" + std::to_string(getpid())
			+ ".h5")).string()};

		void write(const Table& t)
		{
			const hid_t file {H5Fcreate(m_path.c_str(), H5F_ACC_TRUNC,
				H5P_DEFAULT, H5P_DEFAULT)};
			write_doubles(file, "/energy", {t.energy.size()}, t.energy);
			write_doubles(file, "/angle", {t.angle.size()}, t.angle);
			write_doubles(file, "/yield", t.yield_shape.empty()
				? std::vector<hsize_t> {t.energy.size(), t.angle.size()}
				: t.yield_shape, t.yield);
			write_doubles(file, "/projectile/Z", {}, {1.0});
			write_doubles(file, "/projectile/mass_amu", {}, {2.014});
			write_doubles(file, "/target/Z", {}, {74.0});
			write_doubles(file, "/target/mass_amu", {}, {183.84});
			if (t.target_Es) write_doubles(file, "/target/Es", {}, {8.79});
			H5Fclose(file);
		}

		void TearDown() override {std::filesystem::remove(m_path);}
	};
}

// The table's values come back exactly at the grid points, and the
// species and energy range are read
TEST_F(SputteringTable, Read)
{
	write(Table {});
	const Sputtering::YieldTable table(m_path);

	EXPECT_EQ(table.path(), m_path);
	EXPECT_DOUBLE_EQ(table.energy_min(), 10.0);
	EXPECT_DOUBLE_EQ(table.energy_max(), 10000.0);
	EXPECT_DOUBLE_EQ(table.projectile_Z(), 1.0);
	EXPECT_DOUBLE_EQ(table.projectile_mass_amu(), 2.014);
	EXPECT_DOUBLE_EQ(table.target_Z(), 74.0);
	EXPECT_DOUBLE_EQ(table.target_mass_amu(), 183.84);
	EXPECT_DOUBLE_EQ(table.target_Es(), 8.79);

	const Table t {};
	for (std::size_t i {}; i < t.energy.size(); ++i)
	{
		for (std::size_t j {}; j < t.angle.size(); ++j)
		{
			EXPECT_NEAR(table.get_yield(t.energy[i], t.angle[j]),
				t.yield[i * t.angle.size() + j], 1e-15);
		}
	}
}

// Between grid points the yield is linear in log10(E) and angle
TEST_F(SputteringTable, Interpolates)
{
	write(Table {});
	const Sputtering::YieldTable table(m_path);

	for (const double e : {150.0, 316.2, 999.0, 2500.0, 7777.0})
	{
		for (const double a : {0.0, 12.5, 45.0, 61.0, 88.0})
		{
			EXPECT_NEAR(table.get_yield(e, a), linear_yield(e, a), 1e-14)
				<< "E = " << e << ", angle = " << a;
		}
	}

	// Between the zero row and the next one, the yield goes linearly in
	// log10(E) from 0, e.g., half way at the geometric mean
	const double e_mid {std::sqrt(10.0 * 100.0)};
	EXPECT_NEAR(table.get_yield(e_mid, 30.0),
		0.5 * linear_yield(100.0, 30.0), 1e-14);
}

// Below the table the yield is 0, above it the highest energy is used,
// and angles are clamped to the table's
TEST_F(SputteringTable, OutsideTable)
{
	write(Table {});
	const Sputtering::YieldTable table(m_path);

	EXPECT_EQ(table.get_yield(9.99, 30.0), 0.0);
	EXPECT_EQ(table.get_yield(0.0, 30.0), 0.0);
	EXPECT_EQ(table.get_yield(-5.0, 30.0), 0.0);
	EXPECT_EQ(table.get_yield(NAN, 30.0), 0.0);

	EXPECT_NEAR(table.get_yield(1e6, 30.0), linear_yield(10000.0, 30.0),
		1e-14);
	EXPECT_NEAR(table.get_yield(1000.0, 90.0), linear_yield(1000.0, 89.0),
		1e-14);
	EXPECT_NEAR(table.get_yield(1000.0, -10.0), linear_yield(1000.0, 0.0),
		1e-14);
}

// Bad tables stop the run with a message saying what's wrong
TEST_F(SputteringTable, BadTables)
{
	EXPECT_EXIT(Sputtering::YieldTable("/no/such/table.h5"),
		testing::ExitedWithCode(EXIT_FAILURE), "table.h5");

	Table missing {};
	missing.target_Es = false;
	write(missing);
	EXPECT_EXIT(Sputtering::YieldTable {m_path},
		testing::ExitedWithCode(EXIT_FAILURE), "/target/Es");

	Table shape {};
	shape.yield_shape = {4, 2, 2};
	write(shape);
	EXPECT_EXIT(Sputtering::YieldTable {m_path},
		testing::ExitedWithCode(EXIT_FAILURE), "wrong shape");

	Table one {};
	one.energy = {100.0};
	one.yield = {0.1, 0.1, 0.1, 0.1};
	write(one);
	EXPECT_EXIT(Sputtering::YieldTable {m_path},
		testing::ExitedWithCode(EXIT_FAILURE), "At least 2");

	Table order {};
	order.energy = {10.0, 1000.0, 100.0, 10000.0};
	write(order);
	EXPECT_EXIT(Sputtering::YieldTable {m_path},
		testing::ExitedWithCode(EXIT_FAILURE), "not strictly increasing");

	Table zero {};
	zero.energy[0] = 0.0;
	write(zero);
	EXPECT_EXIT(Sputtering::YieldTable {m_path},
		testing::ExitedWithCode(EXIT_FAILURE), "energies <= 0");

	Table angle {};
	angle.angle.back() = 95.0;
	write(angle);
	EXPECT_EXIT(Sputtering::YieldTable {m_path},
		testing::ExitedWithCode(EXIT_FAILURE), "outside 0 to 90");

	Table negative {};
	negative.yield[5] = -0.1;
	write(negative);
	EXPECT_EXIT(Sputtering::YieldTable {m_path},
		testing::ExitedWithCode(EXIT_FAILURE), "negative or non-finite");
}

// The D on W table in the repository reads, and gives its own values back
// at its grid points
TEST(SputteringData, DOnW)
{
	const std::string path {std::string {FLAN_SOURCE_DIR}
		+ "/data/sputtering/D_on_W.h5"};
	if (!std::filesystem::exists(path))
	{
		GTEST_SKIP() << path << " doesn't exist yet.";
	}

	const Sputtering::YieldTable table(path);
	EXPECT_DOUBLE_EQ(table.projectile_Z(), 1.0);
	EXPECT_DOUBLE_EQ(table.target_Z(), 74.0);

	const Hdf5::File file(path);
	const auto energy {Hdf5::read_dataset<double>(file, "/energy").data};
	const auto angle {Hdf5::read_dataset<double>(file, "/angle").data};
	const auto yield {Hdf5::read_dataset<double>(file, "/yield").data};
	for (std::size_t i {}; i < energy.size(); ++i)
	{
		for (std::size_t j {}; j < angle.size(); ++j)
		{
			EXPECT_NEAR(table.get_yield(energy[i], angle[j]),
				yield[i * angle.size() + j],
				1e-12 * yield[i * angle.size() + j])
				<< "E = " << energy[i] << ", angle = " << angle[j];
		}
	}
}
