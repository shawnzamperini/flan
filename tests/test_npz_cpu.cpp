/*
* Tests for reading .npz files (npz.cpp): the zip layouts numpy writes,
* the .npy header, the data types, and the errors for files it can't
* read. The files are written by hermes_test_files.h.
*/

#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <gtest/gtest.h>
#include <map>
#include <string>
#include <unistd.h>
#include <vector>

#include "hermes_test_files.h"
#include "npz.h"

namespace
{
	class NpzTest : public testing::Test
	{
	protected:
		std::string m_path {(std::filesystem::temp_directory_path()
			/ ("flan_test_npz_" + std::to_string(getpid()) + ".npz"))
			.string()};

		void TearDown() override {std::filesystem::remove(m_path);}

		// A small set of arrays: a 3D float64, a scalar and a 1D one
		std::map<std::string, HermesTest::NpyArray> arrays() const
		{
			std::vector<double> v(2 * 3 * 4);
			for (std::size_t n {}; n < v.size(); ++n) v[n] = 0.5 * n - 3.0;
			return {{"field", HermesTest::f8({2, 3, 4}, v)},
				{"t", HermesTest::f8({}, {123.25})},
				{"line", HermesTest::f8({5}, {1, 2, 3, 4, 5})}};
		}
	};

	// Reads every array back, checking names, shapes and values
	void check_arrays(const std::string& path)
	{
		const Npz::File f(path);
		EXPECT_EQ(f.names(), (std::vector<std::string> {"field", "line",
			"t"}));
		EXPECT_TRUE(f.has("field"));
		EXPECT_FALSE(f.has("missing"));

		const auto a {f.read("field")};
		EXPECT_EQ(a.shape, (std::vector<std::size_t> {2, 3, 4}));
		ASSERT_EQ(a.data.size(), 24u);
		for (std::size_t n {}; n < a.data.size(); ++n)
			EXPECT_EQ(a.data[n], 0.5 * n - 3.0);

		const auto t {f.read("t")};
		EXPECT_TRUE(t.shape.empty());
		ASSERT_EQ(t.data.size(), 1u);
		EXPECT_EQ(t.data[0], 123.25);

		const auto line {f.read("line")};
		EXPECT_EQ(line.shape, (std::vector<std::size_t> {5}));
		EXPECT_EQ(line.data, (std::vector<double> {1, 2, 3, 4, 5}));
	}
}

// numpy.savez's layout: zip64 local headers with no sizes
TEST_F(NpzTest, Zip64LocalHeaders)
{
	HermesTest::write_npz(m_path, arrays());
	check_arrays(m_path);
}

// Plain local headers with the sizes in them
TEST_F(NpzTest, PlainLocalHeaders)
{
	HermesTest::write_npz(m_path, arrays(), {.zip64_local = false});
	check_arrays(m_path);
}

// Offsets in the central directory's zip64 extra field
TEST_F(NpzTest, Zip64CentralDirectory)
{
	HermesTest::write_npz(m_path, arrays(), {.zip64_central = true});
	check_arrays(m_path);
}

// float32, int64 and int32 are converted to double
TEST_F(NpzTest, OtherTypes)
{
	auto raw = [](const std::string& descr, auto values)
	{
		HermesTest::NpyArray a {descr, {values.size()}, {}};
		a.bytes.resize(values.size() * sizeof(values[0]));
		std::memcpy(a.bytes.data(), values.data(), a.bytes.size());
		return a;
	};
	HermesTest::write_npz(m_path, {
		{"f4", raw("<f4", std::vector<float> {1.5f, -2.25f})},
		{"i8", raw("<i8", std::vector<std::int64_t> {-7, 1LL << 40})},
		{"i4", raw("<i4", std::vector<std::int32_t> {3, -4})}});
	const Npz::File f(m_path);
	EXPECT_EQ(f.read("f4").data, (std::vector<double> {1.5, -2.25}));
	EXPECT_EQ(f.read("i8").data, (std::vector<double> {-7.0,
		static_cast<double>(1LL << 40)}));
	EXPECT_EQ(f.read("i4").data, (std::vector<double> {3.0, -4.0}));
}

// Files the reader can't handle stop the run with a message saying why
TEST_F(NpzTest, Errors)
{
	HermesTest::write_npz(m_path, arrays(), {.method = 8});
	EXPECT_EXIT(Npz::File(m_path).read("field"),
		testing::ExitedWithCode(EXIT_FAILURE), "compressed");

	HermesTest::write_npz(m_path, arrays());
	EXPECT_EXIT(Npz::File(m_path).read("missing"),
		testing::ExitedWithCode(EXIT_FAILURE), "has no array missing");

	auto fortran {arrays()};
	fortran["field"].fortran = true;
	HermesTest::write_npz(m_path, fortran);
	EXPECT_EXIT(Npz::File(m_path).read("field"),
		testing::ExitedWithCode(EXIT_FAILURE), "Fortran order");

	auto bad_type {arrays()};
	bad_type["line"].descr = ">f8";
	HermesTest::write_npz(m_path, bad_type);
	EXPECT_EXIT(Npz::File(m_path).read("line"),
		testing::ExitedWithCode(EXIT_FAILURE), "only little-endian");

	auto short_data {arrays()};
	short_data["line"].bytes.resize(16);
	HermesTest::write_npz(m_path, short_data);
	EXPECT_EXIT(Npz::File(m_path).read("line"),
		testing::ExitedWithCode(EXIT_FAILURE), "bytes of data");

	{
		std::ofstream out(m_path, std::ios::binary);
		out << "not a zip file at all, just some text";
	}
	EXPECT_EXIT(Npz::File {m_path}, testing::ExitedWithCode(EXIT_FAILURE),
		"isn't a zip");
}
