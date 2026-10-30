/*
* Tests for the HDF5 helpers in hdf5_utils.h, on a small file written by the
* test itself.
*/

#include <filesystem>
#include <gtest/gtest.h>
#include <string>
#include <unistd.h>
#include <vector>

#include "hdf5.h"
#include "hdf5_utils.h"


namespace
{
	// Writes a dataset of doubles with the given shape, creating groups
	// along the way
	void write_doubles(hid_t file, const std::string& name,
		const std::vector<hsize_t>& shape, const std::vector<double>& data)
	{
		const hid_t lcpl {H5Pcreate(H5P_LINK_CREATE)};
		H5Pset_create_intermediate_group(lcpl, 1);
		const hid_t space {H5Screate_simple(static_cast<int>(shape.size()),
			shape.data(), nullptr)};
		const hid_t dset {H5Dcreate2(file, name.c_str(), H5T_NATIVE_DOUBLE,
			space, lcpl, H5P_DEFAULT, H5P_DEFAULT)};
		H5Dwrite(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
			data.data());
		H5Dclose(dset);
		H5Sclose(space);
		H5Pclose(lcpl);
	}

	// Test file with a 2x3 dataset, a one-element dataset and an int
	// scalar, removed when the test ends
	class Hdf5Utils : public testing::Test
	{
	protected:
		std::string m_path {(std::filesystem::temp_directory_path()
			/ ("flan_test_hdf5_utils_" + std::to_string(getpid())
			+ ".h5")).string()};

		void SetUp() override
		{
			const hid_t file {H5Fcreate(m_path.c_str(), H5F_ACC_TRUNC,
				H5P_DEFAULT, H5P_DEFAULT)};
			write_doubles(file, "/zone1/spec0/n", {2, 3},
				{0.0, 1.0, 2.0, 10.0, 11.0, 12.0});
			write_doubles(file, "/B0", {1}, {1.9781});

			const int nzones {13};
			const hid_t space {H5Screate(H5S_SCALAR)};
			const hid_t dset {H5Dcreate2(file, "/NZones", H5T_NATIVE_INT,
				space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT)};
			H5Dwrite(dset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
				&nzones);
			H5Dclose(dset);
			H5Sclose(space);
			H5Fclose(file);
		}

		void TearDown() override {std::filesystem::remove(m_path);}
	};
}

// A dataset is read flattened in C order with its shape
TEST_F(Hdf5Utils, ReadDataset)
{
	const Hdf5::File file(m_path);
	const auto ds {Hdf5::read_dataset<double>(file, "/zone1/spec0/n")};
	EXPECT_EQ(ds.shape, (std::vector<hsize_t> {2, 3}));
	EXPECT_EQ(ds.data, (std::vector<double> {0, 1, 2, 10, 11, 12}));

	// Read as int, converted by HDF5
	const auto ints {Hdf5::read_dataset<int>(file, "/zone1/spec0/n")};
	EXPECT_EQ(ints.data, (std::vector<int> {0, 1, 2, 10, 11, 12}));
}

// One-element and scalar datasets read as a single value
TEST_F(Hdf5Utils, ReadScalar)
{
	const Hdf5::File file(m_path);
	EXPECT_DOUBLE_EQ(Hdf5::read_scalar<double>(file, "/B0"), 1.9781);
	EXPECT_EQ(Hdf5::read_scalar<int>(file, "/NZones"), 13);
	EXPECT_EXIT(Hdf5::read_scalar<double>(file, "/zone1/spec0/n"),
		testing::ExitedWithCode(EXIT_FAILURE), "should hold one value");
}

TEST_F(Hdf5Utils, Exists)
{
	const Hdf5::File file(m_path);
	EXPECT_TRUE(file.exists("/zone1"));
	EXPECT_TRUE(file.exists("/zone1/spec0/n"));
	EXPECT_FALSE(file.exists("/zone1/spec1/n"));
	EXPECT_FALSE(file.exists("/zone2/spec0/n"));
}

// Missing files and datasets stop the run with a message naming them
TEST_F(Hdf5Utils, Missing)
{
	EXPECT_EXIT(Hdf5::File("/no/such/file.h5"),
		testing::ExitedWithCode(EXIT_FAILURE), "/no/such/file.h5");

	const Hdf5::File file(m_path);
	EXPECT_EXIT(Hdf5::read_dataset<double>(file, "/zone2/PHI"),
		testing::ExitedWithCode(EXIT_FAILURE), "/zone2/PHI in .*");
}
