/**
* @file npz.h
*
* @brief Header file for npz.cpp: reading arrays from uncompressed NumPy
* .npz files.
*/

#pragma once

#include <cstddef>
#include <cstdint>
#include <map>
#include <string>
#include <vector>

namespace Npz
{
	/**
	* An array read from a .npz file, converted to double. data is in C
	* (row-major) order, and shape is empty for a scalar.
	*/
	struct Array
	{
		std::vector<std::size_t> shape {};	///< Size of each dimension
		std::vector<double> data {};		///< Values, C order
	};

	/**
	* An open .npz file: a zip archive of .npy files, one per array, as
	* written by numpy.savez. Only stored (uncompressed) entries can be read,
	* so files from numpy.savez_compressed stop the run. The archive's
	* directory is read when the file is opened, and each array when it's
	* asked for.
	*/
	class File
	{
	private:
		/// Where an entry's data is and how big it is
		struct Entry
		{
			std::uint64_t header_offset {};	///< Local header in the file
			std::uint64_t size {};			///< Bytes of data
			std::uint16_t method {};		///< 0 = stored
		};

		std::string m_path {};
		std::map<std::string, Entry> m_entries {};	///< By array name

	public:
		explicit File(const std::string& path);

		const std::string& path() const;
		std::vector<std::string> names() const;
		bool has(const std::string& name) const;
		Array read(const std::string& name) const;
	};

	Array parse_npy(const std::vector<char>& bytes, const std::string& what);
}
