/**
* @file npz.cpp
*
* @brief Reads arrays from uncompressed NumPy .npz files.
*
* A .npz file (numpy.savez) is a zip archive holding one .npy file per
* array. numpy stores the entries without compression but writes them as
* zip64 entries, whose local headers give no sizes, so the sizes and
* offsets are taken from the archive's central directory instead. Each
* .npy file is a short text header (dtype, memory order and shape)
* followed by the raw values. This is enough for the plasma snapshots of
* codes like JAX-Hermes3 without depending on a zip library.
*/

#include <algorithm>
#include <bit>
#include <cctype>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <map>
#include <string>
#include <vector>

#include "errors.h"
#include "npz.h"

namespace Npz
{
	static_assert(std::endian::native == std::endian::little,
		"The .npz reader assumes a little-endian machine.");

	namespace
	{
		// Zip record signatures
		constexpr std::uint32_t sig_local {0x04034b50};
		constexpr std::uint32_t sig_central {0x02014b50};
		constexpr std::uint32_t sig_eocd {0x06054b50};
		constexpr std::uint32_t sig_zip64_eocd {0x06064b50};
		constexpr std::uint32_t sig_zip64_locator {0x07064b50};

		// Little-endian unsigned integer of type T at bytes[pos]
		template <typename T>
		T get(const std::vector<char>& bytes, const std::size_t pos)
		{
			T v {};
			std::memcpy(&v, bytes.data() + pos, sizeof(T));
			return v;
		}

		// Reads n bytes starting at offset, stopping the run if the file is
		// too short
		std::vector<char> read_bytes(std::ifstream& in,
			const std::uint64_t offset, const std::uint64_t n,
			const std::string& path)
		{
			std::vector<char> out(n);
			in.seekg(static_cast<std::streamoff>(offset));
			in.read(out.data(), static_cast<std::streamsize>(n));
			if (!in)
			{
				Errors::fatal_error("Couldn't read " + std::to_string(n)
					+ " bytes at offset " + std::to_string(offset) + " of "
					+ path + ". Is it a complete .npz file?");
			}
			return out;
		}

		// The value of key in a .npy header dict, e.g. "'<f8'" for 'descr',
		// up to the next comma at depth 0 (so a shape tuple stays whole)
		std::string header_value(const std::string& header,
			const std::string& key, const std::string& what)
		{
			const std::string k {"'" + key + "'"};
			auto pos {header.find(k)};
			if (pos == std::string::npos)
			{
				Errors::fatal_error("The .npy header of " + what + " has no "
					+ k + ": " + header);
			}
			pos = header.find(':', pos + k.size());
			std::size_t end {pos + 1};
			int depth {};
			for (; end < header.size(); ++end)
			{
				const char c {header[end]};
				if (c == '(') ++depth;
				else if (c == ')') --depth;
				else if ((c == ',' || c == '}') && depth == 0) break;
			}
			std::string v {header.substr(pos + 1, end - pos - 1)};
			const auto first {v.find_first_not_of(" ")};
			const auto last {v.find_last_not_of(" ")};
			return first == std::string::npos ? ""
				: v.substr(first, last - first + 1);
		}
	}

	/**
	* Opens a .npz file and reads its directory.
	*
	* The end of central directory record (or its zip64 version) gives
	* where the central directory is, and each entry there gives the
	* offset of its local header and its size, from the zip64 extra field
	* when the 32-bit fields don't hold them. Entry names lose their .npy
	* extension, so they're the names given to numpy.savez.
	*
	* @param path Path to the .npz file.
	*/
	File::File(const std::string& path)
		: m_path {path}
	{
		std::ifstream in(path, std::ios::binary | std::ios::ate);
		if (!in) Errors::fatal_error("Couldn't open " + path + ".");
		const std::uint64_t size {static_cast<std::uint64_t>(in.tellg())};

		// The end of central directory record is in the last 22 bytes plus
		// a comment of up to 65535 bytes
		const std::uint64_t tail_len {std::min<std::uint64_t>(size,
			22 + 65535)};
		const std::uint64_t tail_start {size - tail_len};
		const auto tail {read_bytes(in, tail_start, tail_len, path)};
		std::int64_t eocd {-1};
		for (std::int64_t p {static_cast<std::int64_t>(tail_len) - 22};
			p >= 0; --p)
		{
			if (get<std::uint32_t>(tail, p) == sig_eocd) {eocd = p; break;}
		}
		if (eocd < 0)
		{
			Errors::fatal_error(path + " isn't a zip (.npz) file: no end "
				"of central directory record.");
		}

		std::uint64_t n_entries {get<std::uint16_t>(tail, eocd + 10)};
		std::uint64_t cd_size {get<std::uint32_t>(tail, eocd + 12)};
		std::uint64_t cd_offset {get<std::uint32_t>(tail, eocd + 16)};

		// zip64: the locator is just before the end of central directory
		// record and points to the zip64 one, which has the real values
		if (n_entries == 0xFFFF || cd_size == 0xFFFFFFFF
			|| cd_offset == 0xFFFFFFFF)
		{
			const std::int64_t loc {eocd - 20};
			if (loc < 0 || get<std::uint32_t>(tail, loc) != sig_zip64_locator)
			{
				Errors::fatal_error(path + " needs a zip64 end of central "
					"directory record, but has no locator for it.");
			}
			const std::uint64_t z64 {get<std::uint64_t>(tail, loc + 8)};
			const auto rec {read_bytes(in, z64, 56, path)};
			if (get<std::uint32_t>(rec, 0) != sig_zip64_eocd)
			{
				Errors::fatal_error(path + " has a bad zip64 end of "
					"central directory record.");
			}
			n_entries = get<std::uint64_t>(rec, 32);
			cd_size = get<std::uint64_t>(rec, 40);
			cd_offset = get<std::uint64_t>(rec, 48);
		}

		const auto cd {read_bytes(in, cd_offset, cd_size, path)};
		std::size_t p {};
		for (std::uint64_t e {}; e < n_entries; ++e)
		{
			if (p + 46 > cd.size() || get<std::uint32_t>(cd, p) != sig_central)
			{
				Errors::fatal_error(path + " has a bad central directory "
					"entry.");
			}
			Entry entry {};
			entry.method = get<std::uint16_t>(cd, p + 10);
			std::uint64_t comp_size {get<std::uint32_t>(cd, p + 20)};
			std::uint64_t uncomp_size {get<std::uint32_t>(cd, p + 24)};
			const std::uint16_t name_len {get<std::uint16_t>(cd, p + 28)};
			const std::uint16_t extra_len {get<std::uint16_t>(cd, p + 30)};
			const std::uint16_t comment_len {get<std::uint16_t>(cd, p + 32)};
			std::uint64_t offset {get<std::uint32_t>(cd, p + 42)};
			std::string name(cd.data() + p + 46, name_len);

			// zip64 extra field (id 1): the 64-bit values, in this order,
			// of only those 32-bit fields that are 0xFFFFFFFF
			std::size_t x {p + 46 + name_len};
			const std::size_t x_end {x + extra_len};
			while (x + 4 <= x_end)
			{
				const std::uint16_t id {get<std::uint16_t>(cd, x)};
				const std::uint16_t len {get<std::uint16_t>(cd, x + 2)};
				if (id == 1)
				{
					std::size_t v {x + 4};
					if (uncomp_size == 0xFFFFFFFF)
						{uncomp_size = get<std::uint64_t>(cd, v); v += 8;}
					if (comp_size == 0xFFFFFFFF)
						{comp_size = get<std::uint64_t>(cd, v); v += 8;}
					if (offset == 0xFFFFFFFF)
						offset = get<std::uint64_t>(cd, v);
				}
				x += 4 + len;
			}

			entry.size = comp_size;
			entry.header_offset = offset;
			if (name.size() > 4 && name.ends_with(".npy"))
				name.resize(name.size() - 4);
			m_entries[name] = entry;
			p += 46 + name_len + extra_len + comment_len;
		}
	}

	/// The path the file was opened with.
	const std::string& File::path() const {return m_path;}

	/// The names of the arrays in the file, in alphabetical order.
	std::vector<std::string> File::names() const
	{
		std::vector<std::string> out {};
		for (const auto& [name, entry] : m_entries) out.push_back(name);
		return out;
	}

	/// Whether the file has an array with this name.
	bool File::has(const std::string& name) const
		{return m_entries.contains(name);}

	/**
	* Reads one array.
	*
	* @param name Name of the array (as given to numpy.savez).
	* @return The array, converted to double. A missing array or a
	*         compressed entry stops the run.
	*/
	Array File::read(const std::string& name) const
	{
		const auto it {m_entries.find(name)};
		if (it == m_entries.end())
		{
			Errors::fatal_error(m_path + " has no array " + name + ".");
		}
		const Entry& e {it->second};
		if (e.method != 0)
		{
			Errors::fatal_error("Array " + name + " in " + m_path + " is "
				"compressed (zip method " + std::to_string(e.method) + "). "
				"Only uncompressed .npz files (numpy.savez, not "
				"savez_compressed) can be read.");
		}

		std::ifstream in(m_path, std::ios::binary);
		if (!in) Errors::fatal_error("Couldn't open " + m_path + ".");
		const auto local {read_bytes(in, e.header_offset, 30, m_path)};
		if (get<std::uint32_t>(local, 0) != sig_local)
		{
			Errors::fatal_error(m_path + " has a bad local header for "
				+ name + ".");
		}
		const std::uint64_t data_offset {e.header_offset + 30
			+ get<std::uint16_t>(local, 26) + get<std::uint16_t>(local, 28)};
		return parse_npy(read_bytes(in, data_offset, e.size, m_path),
			name + " in " + m_path);
	}

	/**
	* Parses the bytes of a .npy file: the magic string and version, the
	* header (a Python dict with descr, fortran_order and shape), and the
	* values.
	*
	* Little-endian float64, float32, int64 and int32 are handled, in C
	* order. Anything else, or a size that doesn't match the shape, stops
	* the run.
	*
	* @param bytes The .npy file's bytes.
	* @param what Description for error messages (e.g., array and file).
	* @return The array, converted to double.
	*/
	Array parse_npy(const std::vector<char>& bytes, const std::string& what)
	{
		if (bytes.size() < 10 || std::memcmp(bytes.data(), "\x93NUMPY", 6))
			Errors::fatal_error(what + " isn't a .npy array.");
		const int major {static_cast<unsigned char>(bytes[6])};
		std::size_t header_len {}, header_start {};
		if (major == 1)
		{
			header_len = get<std::uint16_t>(bytes, 8);
			header_start = 10;
		}
		else if (major == 2 || major == 3)
		{
			header_len = get<std::uint32_t>(bytes, 8);
			header_start = 12;
		}
		else
		{
			Errors::fatal_error(what + " has .npy version "
				+ std::to_string(major) + ", which isn't handled.");
		}
		if (header_start + header_len > bytes.size())
			Errors::fatal_error(what + " has a truncated .npy header.");
		const std::string header(bytes.data() + header_start, header_len);

		const std::string descr {header_value(header, "descr", what)};
		const std::string order {header_value(header, "fortran_order", what)};
		const std::string shape_str {header_value(header, "shape", what)};
		if (order != "False")
		{
			Errors::fatal_error(what + " is in Fortran order, but only C "
				"order is handled.");
		}

		Array out {};
		std::size_t count {1};
		for (std::size_t p {shape_str.find('(') + 1};
			p < shape_str.size(); )
		{
			while (p < shape_str.size()
				&& !std::isdigit(static_cast<unsigned char>(shape_str[p]))) ++p;
			if (p >= shape_str.size()) break;
			std::size_t q {p};
			while (q < shape_str.size()
				&& std::isdigit(static_cast<unsigned char>(shape_str[q]))) ++q;
			out.shape.push_back(std::stoull(shape_str.substr(p, q - p)));
			count *= out.shape.back();
			p = q;
		}

		const std::size_t data_start {header_start + header_len};
		const std::size_t nbytes {bytes.size() - data_start};
		const char* data {bytes.data() + data_start};
		auto convert = [&]<typename T>(T)
		{
			if (nbytes != count * sizeof(T))
			{
				Errors::fatal_error(what + " has " + std::to_string(nbytes)
					+ " bytes of data, but its shape " + shape_str
					+ " and type " + descr + " need "
					+ std::to_string(count * sizeof(T)) + ".");
			}
			out.data.resize(count);
			for (std::size_t n {}; n < count; ++n)
			{
				T v {};
				std::memcpy(&v, data + n * sizeof(T), sizeof(T));
				out.data[n] = static_cast<double>(v);
			}
		};
		if (descr == "'<f8'") convert(double {});
		else if (descr == "'<f4'") convert(float {});
		else if (descr == "'<i8'") convert(std::int64_t {});
		else if (descr == "'<i4'") convert(std::int32_t {});
		else
		{
			Errors::fatal_error(what + " has type " + descr + ", but only "
				"little-endian float64, float32, int64 and int32 are "
				"handled.");
		}
		return out;
	}
}
