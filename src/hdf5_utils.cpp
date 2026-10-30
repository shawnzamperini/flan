/**
* @file hdf5_utils.cpp
*
* @brief Small helpers for reading HDF5 files with the HDF5 C API, used to
*        read SOLEDGE3X output.
*
* Any problem (a missing file or dataset, or a failed read) stops the run
* with Errors::fatal_error and a message naming the file and dataset.
*/

#include <cstddef>
#include <filesystem>
#include <string>
#include <type_traits>
#include <vector>

#include "errors.h"
#include "hdf5.h"
#include "hdf5_utils.h"

namespace Hdf5
{
	/**
	* Opens an HDF5 file read-only.
	*
	* HDF5's own error printing is turned off, since every failure is
	* reported with Errors::fatal_error instead. A file that doesn't exist
	* or can't be opened stops the run.
	*
	* @param path Path to the file.
	*/
	File::File(const std::string& path)
		: m_path {path}
	{
		H5Eset_auto2(H5E_DEFAULT, nullptr, nullptr);

		if (!std::filesystem::exists(path))
		{
			Errors::fatal_error("HDF5 file " + path + " does not exist.");
		}

		m_id = H5Fopen(path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
		if (m_id < 0)
		{
			Errors::fatal_error("Could not open HDF5 file " + path + ".");
		}
	}

	/**
	* Closes the file.
	*/
	File::~File()
	{
		if (m_id >= 0) H5Fclose(m_id);
	}

	/**
	* @name Accessors
	*
	* Return the HDF5 id of the file or the path it was opened from.
	*/
	///@{
	hid_t File::id() const {return m_id;}

	const std::string& File::path() const {return m_path;}
	///@}

	/**
	* Checks whether a group or dataset exists in the file.
	*
	* Each part of the path is checked in turn, since HDF5 can't check a
	* path whose parent groups don't exist.
	*
	* @param name Absolute path inside the file, e.g., "/zone1/spec0/n".
	* @return True if it exists.
	*/
	bool File::exists(const std::string& name) const
	{
		std::string partial {};
		std::size_t start {name.front() == '/' ? 1u : 0u};
		while (start <= name.size())
		{
			const std::size_t end {std::min(name.find('/', start),
				name.size())};
			partial += "/" + name.substr(start, end - start);
			if (H5Lexists(m_id, partial.c_str(), H5P_DEFAULT) <= 0)
				return false;
			start = end + 1;
		}
		return true;
	}

	/**
	* Reads a whole dataset, converting its values to T.
	*
	* The values are returned flattened in C (row-major) order, the order
	* HDF5 stores them in, along with the dataset's shape. A dataset that
	* doesn't exist or can't be read stops the run.
	*
	* @tparam T Type to read the values as: double or int.
	* @param file File to read from.
	* @param name Absolute path of the dataset in the file.
	* @return The values and the shape.
	*/
	template <typename T>
	Dataset<T> read_dataset(const File& file, const std::string& name)
	{
		static_assert(std::is_same_v<T, double> || std::is_same_v<T, int>,
			"Hdf5::read_dataset only reads double or int.");
		const std::string where {name + " in " + file.path()};

		if (!file.exists(name))
		{
			Errors::fatal_error("Dataset " + where + " does not exist.");
		}

		const hid_t dset {H5Dopen2(file.id(), name.c_str(), H5P_DEFAULT)};
		if (dset < 0) Errors::fatal_error("Could not open " + where + ".");

		// Shape
		const hid_t space {H5Dget_space(dset)};
		const int ndims {H5Sget_simple_extent_ndims(space)};
		Dataset<T> out {};
		out.shape.resize(std::max(ndims, 0));
		H5Sget_simple_extent_dims(space, out.shape.data(), nullptr);
		H5Sclose(space);

		hsize_t size {1};
		for (const hsize_t n : out.shape) size *= n;
		out.data.resize(size);

		// Read, converting to T
		const hid_t mem_type {std::is_same_v<T, double> ? H5T_NATIVE_DOUBLE
			: H5T_NATIVE_INT};
		const herr_t status {H5Dread(dset, mem_type, H5S_ALL, H5S_ALL,
			H5P_DEFAULT, out.data.data())};
		H5Dclose(dset);
		if (status < 0) Errors::fatal_error("Could not read " + where + ".");

		return out;
	}

	/**
	* Reads a dataset holding a single value.
	*
	* The dataset may be a scalar or have one element. Any other size stops
	* the run.
	*
	* @tparam T Type to read the value as: double or int.
	* @param file File to read from.
	* @param name Absolute path of the dataset in the file.
	* @return The value.
	*/
	template <typename T>
	T read_scalar(const File& file, const std::string& name)
	{
		const Dataset<T> ds {read_dataset<T>(file, name)};
		if (ds.data.size() != 1)
		{
			Errors::fatal_error("Dataset " + name + " in " + file.path()
				+ " should hold one value, but it has "
				+ std::to_string(ds.data.size()) + ".");
		}
		return ds.data[0];
	}
}

// Instantiate the templates since the definitions are in this file
template Hdf5::Dataset<double> Hdf5::read_dataset<double>(
	const Hdf5::File& file, const std::string& name);
template Hdf5::Dataset<int> Hdf5::read_dataset<int>(
	const Hdf5::File& file, const std::string& name);
template double Hdf5::read_scalar<double>(const Hdf5::File& file,
	const std::string& name);
template int Hdf5::read_scalar<int>(const Hdf5::File& file,
	const std::string& name);
