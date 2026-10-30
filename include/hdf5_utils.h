/**
* @file hdf5_utils.h
*
* @brief Header file for hdf5_utils.cpp
*/

#pragma once

#include <string>
#include <vector>

#include "hdf5.h"

namespace Hdf5
{
	/**
	* An HDF5 file opened read-only, closed when it goes out of scope. See
	* hdf5_utils.cpp for the member function docs.
	*/
	class File
	{
	private:
		hid_t m_id {H5I_INVALID_HID};	///< HDF5 file id
		std::string m_path {};			///< Path the file was opened from

	public:

		explicit File(const std::string& path);

		~File();

		File(const File&) = delete;

		File& operator=(const File&) = delete;

		hid_t id() const;

		const std::string& path() const;

		bool exists(const std::string& name) const;
	};

	/// Values of a dataset and its shape
	template <typename T>
	struct Dataset
	{
		std::vector<T> data {};
		std::vector<hsize_t> shape {};
	};

	template <typename T>
	Dataset<T> read_dataset(const File& file, const std::string& name);

	template <typename T>
	T read_scalar(const File& file, const std::string& name);
}
