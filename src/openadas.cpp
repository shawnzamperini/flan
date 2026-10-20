/**
* @file  openadas.cpp
*
* @brief Definitions for the OpenADAS class, which reads in and stores the
*        ionization and recombination rate coefficients from OpenADAS. Also
*        contains the routines for updating an impurity's charge state.
*
* The rate coefficients are read from the OpenADAS SCD (ionization) and ACD
* (recombination) files and bilinearly interpolated in Te and ne. Each time
* step, each particle can ionize or recombine with a probability of
* rate * ne * dt.
*/
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <omp.h>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

#include "background.h"
#include "errors.h"
#include "openadas.h"
#include "openadas_device.h"
#include "pcg32.h"
#include "slots.h"
#include "utilities.h"
#include "vectors.h"

#ifdef USE_CUDA
#include <cuda_runtime.h>
#endif

namespace OpenADAS
{
	/**
	* Creates an empty OpenADAS object with no rate coefficients loaded.
	*/
	OpenADAS::OpenADAS() {};


	/**
	* Creates an OpenADAS object and loads the rate coefficients from the
	* matching OpenADAS file.
	*
	* The run is stopped with Errors::fatal_error if the file can't be
	* opened, e.g., because openadas_root or openadas_year is wrong.
	*
	* @param openadas_root Full path to the directory containing the
	*        OpenADAS data.
	* @param openadas_year Year of the OpenADAS data to load, e.g., 50.
	* @param imp_atom_num Atomic number of the impurity to load data for.
	* @param rate_type Either "scd" (ionization) or "acd" (recombination).
	*/
	OpenADAS::OpenADAS(const std::string_view openadas_root,
		const int openadas_year, const int imp_atom_num,
		const std::string_view rate_type)
	{
		// Create full path to recombination or ionization rate files.
		// These are the ACD and SCD files, respectively.
		std::string openadas_path {get_openadas_path(openadas_root,
			openadas_year, imp_atom_num, rate_type)};
		std::cout << "ADAS: Loading " << openadas_path << '\n' << std::flush;

		// Open file
		std::ifstream openadas_stream {openadas_path};
		if (!openadas_stream)
		{
			Errors::fatal_error("Could not open the OpenADAS file " 
				+ openadas_path + ". Check openadas_root and "
				"openadas_year.");
		}

		// Read in the rate coefficients along with all the header info.
		// Everything is stored internally in the class.
		read_rate_coefficients(openadas_stream);
	}


	/**
	* Assembles the full path to an OpenADAS data file, e.g.,
	* [openadas_root]/scd50/scd50_w.dat.
	*
	* The run is stopped with Errors::fatal_error if rate_type is not "scd"
	* or "acd", or if there is no OpenADAS data for imp_atom_num. Supported
	* elements are H through Ne, Al, Si, Cl, Ar, Fe, Kr, Mo, Xe and W.
	*
	* @param openadas_root Full path to the directory containing the
	*        OpenADAS data.
	* @param openadas_year Year of the OpenADAS data to load, e.g., 50.
	* @param imp_atom_num Atomic number of the impurity to load data for.
	* @param rate_type Either "scd" (ionization) or "acd" (recombination).
	* @return Full path to the OpenADAS file.
	*/
	std::string OpenADAS::get_openadas_path(
		const std::string_view openadas_root, const int openadas_year,
		const int imp_atom_num, const std::string_view rate_type)
	{
		// Check for valid entry for rate type
		if (rate_type != "scd" && rate_type != "acd")
		{
			Errors::fatal_error("rate_type = " + std::string(rate_type) 
				+ ". Only \"scd\" and \"acd\" are allowed.");
		}

		// Select element name from atomic number.
		std::string imp_name {};
		if (imp_atom_num == 1) imp_name = "h";
		else if (imp_atom_num == 2) imp_name = "he";
		else if (imp_atom_num == 3) imp_name = "li";
		else if (imp_atom_num == 4) imp_name = "be";
		else if (imp_atom_num == 5) imp_name = "b";
		else if (imp_atom_num == 6) imp_name = "c";
		else if (imp_atom_num == 7) imp_name = "n";
		else if (imp_atom_num == 8) imp_name = "o";
		else if (imp_atom_num == 9) imp_name = "f";
		else if (imp_atom_num == 10) imp_name = "ne";
		else if (imp_atom_num == 13) imp_name = "al";
		else if (imp_atom_num == 14) imp_name = "si";
		else if (imp_atom_num == 17) imp_name = "cl";
		else if (imp_atom_num == 18) imp_name = "ar";
		else if (imp_atom_num == 26) imp_name = "fe";
		else if (imp_atom_num == 36) imp_name = "kr";
		else if (imp_atom_num == 42) imp_name = "mo";
		else if (imp_atom_num == 54) imp_name = "xe";
		else if (imp_atom_num == 74) imp_name = "w";
		else
		{
			Errors::fatal_error("imp_atom_num = " 
				+ std::to_string(imp_atom_num) + " has no OpenADAS data in "
				"Flan. Supported atomic numbers are 1-10, 13, 14, 17, 18, 26, "
				"36, 42, 54 and 74.");
		}

		// Assemble the full path. Example:
		// /home/zamp/openadas/scd50/scd50_w.dat
		std::string openadas_path {openadas_root};
		openadas_path.append("/");
		openadas_path.append(rate_type);
		openadas_path.append(std::to_string(openadas_year));
		openadas_path.append("/");
		openadas_path.append(rate_type);
		openadas_path.append(std::to_string(openadas_year));
		openadas_path.append("_");
		openadas_path.append(imp_name);
		openadas_path.append(".dat");

		return openadas_path;
	}


	/**
	* Reads in the header, densities, temperatures and rate coefficients
	* from an OpenADAS file, storing them in the class.
	*
	* The densities and temperatures are stored in the file as log10(ne
	* [cm^-3]) and log10(Te [eV]), and the rate coefficients as log10(rate
	* [cm^3/s]). They are converted to ne (m^-3), Te (eV) and rate (m^3/s).
	* The rates are stored with dimensions (charge, Te, ne). The run is
	* stopped with Errors::fatal_error if the file doesn't cover every
	* charge state (1 to Z), since the rates would be indexed incorrectly,
	* or if there are fewer than two densities or temperatures, since the
	* bilinear interpolation needs at least two of each.
	*
	* @param openadas_stream Open file stream of the OpenADAS file.
	*/
	void OpenADAS::read_rate_coefficients(std::ifstream& openadas_stream)
	{
		// Get header info line as a string. Header contains (all on one line):
		// [atomic number] [# densities] [# temperatures] [minimum charge]
		// [maximum charge] [element name] ... other stuff we don't need.
		std::string line {};
		std::getline(openadas_stream, line);

		// Convert header to vector of strings for each element.
		std::vector<std::string> header {Utilities::split_str_at_spaces(line)};

		// Now pull out each into variables, casting as appropriate. Don't
		// really care about the rest of the header elements, the first 5 are
		// all we need here.
		m_atomic_number = Utilities::str_as_int(header[0]);
		m_ndens = Utilities::str_as_int(header[1]);
		m_ntemp = Utilities::str_as_int(header[2]);
		m_charge_low = Utilities::str_as_int(header[3]);
		m_charge_high = Utilities::str_as_int(header[4]);

		// The rate coefficients are indexed assuming the file covers every
		// charge state, i.e., IZ1 runs from 1 to Z. Slice 0 is then
		// W0 --> W1+ for ionization and W1+ --> W0 for recombination, which
		// is why get_rate_coeff and get_charge_slice index ionization by q
		// and recombination by q-1. A partial file would shift every rate.
		if (m_charge_low != 1 || m_charge_high != m_atomic_number)
		{
			Errors::fatal_error("The OpenADAS file covers charges " 
				+ std::to_string(m_charge_low) + " to " 
				+ std::to_string(m_charge_high) + ", but Flan expects 1 to "
				+ std::to_string(m_atomic_number) + ". The rate "
				"coefficients would be indexed incorrectly.");
		}

		// We need at least two ne and te values in order to do a bilinear
		// interpolation of the rate coefficients later in the program.
		if (m_ndens < 2)
		{
			Errors::fatal_error("The OpenADAS file has " 
				+ std::to_string(m_ndens) + " density value(s), but at least "
				"2 are needed to interpolate the rate coefficients.");
		}
		if (m_ntemp < 2)
		{
			Errors::fatal_error("The OpenADAS file has " 
				+ std::to_string(m_ntemp) + " temperature value(s), but at "
				"least 2 are needed to interpolate the rate coefficients.");
		}

		// Resize m_ne and m_te since we know how big they are now.
		m_ne.resize(m_ndens);
		m_te.resize(m_ntemp);

		// Next line is just a bunch of '----', so throw away.
		std::getline(openadas_stream, line);

		// First read in the electron density values used in this file (these
		// are actually log10(ne) values in cm-3). We do this by reading a
		// line at a time, and then putting those values into m_ne until we
		// have read ndens values.
		int count {};
		while (count < m_ndens)
		{
			// Read line, put into vector of strings.
			std::getline(openadas_stream, line);
			std::vector<std::string> words {
				Utilities::split_str_at_spaces(line)};

			// Convert each to a double value, undoing the log10(ne).
			for (auto word : words)
			{
				// Reads as log10(ne [cm-3]). 10 to the power of this to put
				// into cm-3, then * 1e6 to go to m-3.
				double ne {Utilities::str_as_dbl(word)};
				ne = std::pow(10, ne) * 1e6;

				// Store in m_ne.
				m_ne[count] = ne;
				++count;
			}
		}

		// Repeat for temperatures. These are also in log10(te [eV]).
		count = 0;
		while (count < m_ntemp)
		{
			// Read line, put into vector of strings.
			std::getline(openadas_stream, line);
			std::vector<std::string> words {
				Utilities::split_str_at_spaces(line)};

			// Convert each to a double value.
			for (auto word : words)
			{
				// Units of eV
				double te {Utilities::str_as_dbl(word)};

				// Store in m_te.
				m_te[count] = std::pow(10, te);
				++count;
			}
		}

		// Now read in the rate coefficients into a Vector3D. There are a Z
		// number of tables, where Z is the number of charge states for this
		// ion. The meaning of the tables depends on if it's for ionization
		// (scd) or recombination (acd). Using tungsten (W) as an example:
		//   scd: First table is rate coeff. for W0 --> W1+, second table is
		//        for W1+ --> W2+, up until the 74th table for W73+ --> W74+.
		//   acd: First table is rate coeff. for W1+ --> W0, second table is
		//        for W2+ -> W1+, up until the 74th table for W74+ --> W73+.
		// The tables contain values grouped by the same values of Te. So if
		// ndens = 5 and ntemp = 35, then the first 5 values are rate coeff.
		// at the first Te value for each density. Then the next 5 values are
		// at the second Te value for each density, and so on. This repeats 35
		// times for each Te value.
		// There is no guarantee that those 5 values for each Te will all be on
		// one line, they may span multiple lines, but each time values are
		// printed for a new Te it will start on a new line. This is the logic
		// I have picked up at least.

		// m_rates will be dimensions of (ncharges, ntemp, ndens).
		int ncharges {m_charge_high - m_charge_low + 1};
		m_rates = Vectors::Vector3D<double>(ncharges, m_ntemp, m_ndens);

		for (int i {}; i < ncharges; ++i)
		{
			// First line has some identifying info, but not needed
			std::getline(openadas_stream, line);

			for (int j {}; j < m_ntemp; ++j)
			{
				// Using a counter (k) with a while loop instead of a for loop
				// let's us handle the fact that entries at this charge and te
				// may span multiple lines and may not have the same amount of
				// entries per line. E.g., if ndens=9 there would be
				// alternating lines of 5 and then 4 values. We read values
				// until we've read ndens, irrespective of how many lines it
				// takes to get there.
				int k {};
				while (k < m_ndens)
				{
					// Grab rate coefficients as a vector of strings.
					std::getline(openadas_stream, line);
					std::vector<std::string> words {
						Utilities::split_str_at_spaces(line)};

					// Go through one at a time, adding to m_rates as a double.
					for (auto word : words)
					{
						// Units of log10(rate [cm3/s])
						double rate {Utilities::str_as_dbl(word)};

						// Store in m_rates, converting to m3/s.
						m_rates(i, j, k) = pow(10, rate) * 1e-6;
						++k;
					}
				}
			}
		}
	}


	/**
	* Finds the two neighboring values in a sorted vector that bound a value,
	* for use in a bilinear interpolation.
	*
	* Values outside the range of vec are clamped to its minimum or maximum,
	* so the rates are never extrapolated.
	*
	* @param vec Sorted vector of values, e.g., the OpenADAS densities.
	* @param value Value to find the bounding values of.
	* @return {value, value0, value1, index0, index1}, where value is
	*         clamped to the range of vec, value0 and value1 are the bounding
	*         values and index0 and index1 are their indices.
	*/
	std::tuple<double, double, double, int, int>
		OpenADAS::get_bilinear_interp_vals(const std::vector<double>& vec,
		double value) const
	{
		// Get the iterator of the first index in vec that is larger
		// than value. The bounding values for the bilinear interpolation
		// will then be lower and lower-1, depending on the conditions below.
		auto lower = std::lower_bound(vec.begin(), vec.end(), value);
		int index0 {};
		int index1 {};

		// In this case, value is at/below the minimum value in vec and thus
		// out of bounds. It is not a good idea to linearly extrapolate past
		// the bounds, so we instead treat it as value = min(vec) was passed in
		// (i.e. floor value) by reassigning value to the minimum.
		if (lower == vec.begin())
		{
			value = vec[0];
			index0 = 0;
			index1 = 1;
		}

		// Similar to the previous, except in this case value is above the
		// maximum value in vec.
		else if (lower == vec.end())
		{
			value = vec.back();
			index0 = vec.size() - 2;
			index1 = vec.size() - 1;
		}

		// If value is not out of bounds, then the bounding values for linear
		// interpolation are between lower-1 and lower. We can dereference an
		// iterator for its value.
		else
		{
			index0 = lower - vec.begin() - 1;
			index1 = lower - vec.begin();
		}

		double value0 {vec[index0]};
		double value1 {vec[index1]};

		// Return the values needed to perform the linear interpolation.
		return std::make_tuple(value, value0, value1, index0, index1);
	}


	/**
	* Returns the rate coefficient for a charge state, bilinearly
	* interpolated in Te and ne.
	*
	* ne and Te outside the range of the OpenADAS data are clamped to the
	* minimum or maximum value, so the rates are never extrapolated.
	*
	* @param charge Index of the rate coefficient table. For ionization this
	*        is the charge, and for recombination the charge - 1.
	* @param ne Electron density (m^-3).
	* @param te Electron temperature (eV).
	* @return Rate coefficient (m^3/s).
	*/
	double OpenADAS::get_rate_coeff(int charge, double ne, double te) const
	{
		// Load the values needed to perform a bilinear interpolation. Note
		// that if ne or te are out of range, it is reassigned to the
		// corresponding min/max value in m_ne and m_te, respectively. In other
		// words, we do not extrapolate beyond the values provided by OpenADAS.
		double ne0 {};
		double ne1 {};
		int ne0_index {};
		int ne1_index {};
		std::tie(ne, ne0, ne1, ne0_index, ne1_index) =
			get_bilinear_interp_vals(m_ne, ne);

		double te0 {};
		double te1 {};
		int te0_index {};
		int te1_index {};
		std::tie(te, te0, te1, te0_index, te1_index) =
			get_bilinear_interp_vals(m_te, te);

		// Get the four corner values of the interpolation rectangle.
		// te0, ne0 -> rate00
		// te1, ne0 -> rate10
		// te0, ne1 -> rate01
		// te1, ne1 -> rate11
		double rate00 = m_rates(charge, te0_index, ne0_index);
		double rate10 = m_rates(charge, te1_index, ne0_index);
		double rate01 = m_rates(charge, te0_index, ne1_index);
		double rate11 = m_rates(charge, te1_index, ne1_index);

		// Perform bilinear interpolation. In this case, x = Te and y = ne.
		return Utilities::bilinear_interpolate(te0, ne0, te1, ne1,
			rate00, rate10, rate01, rate11, te, ne);
	}

#ifndef USE_CUDA

	/**
	* CPU-only stub of OpenADAS::to_device, which is defined in
	* src/cuda/openadas.cu since it calls CUDA code. If it gets called, it
	* means you're in the CPU-only version of the code trying to call
	* GPU-specific code, so an error is printed.
	*
	* @param device_id ID of the GPU. Unused.
	* @return Zero-initialized OpenADASDevice struct.
	*/
	OpenADASDevice OpenADAS::to_device(int device_id)
	{
		std::cerr << "Error! OpenADAS::to_device was called but GPU support"
			<< " was not compiled in.\n";
		return OpenADASDevice {};
	}

	/**
	* CPU-only stub of free_oa, which is defined in src/cuda/openadas.cu
	* since it calls CUDA code. If it gets called, an error is printed.
	*
	* @param oa_d Device-side OpenADAS struct. Unused.
	* @param device_id ID of the GPU. Unused.
	*/
	void free_oa(OpenADASDevice& oa_d, int device_id)
	{
		std::cerr << "Error! OpenADAS::free_oa was called but GPU support"
			<< " was not compiled in.\n";
	}

#endif

	/**
	* Calculates the probabilities that a particle ionizes or recombines
	* during one time step.
	*
	* ne and Te are interpolated at the particle location, and each
	* probability is rate * ne * dt. A fully ionized particle can't ionize
	* and a neutral can't recombine, so those probabilities are 0.
	*
	* @param bkg Background plasma.
	* @param oa_ioniz OpenADAS ionization (scd) rate coefficients.
	* @param oa_recomb OpenADAS recombination (acd) rate coefficients.
	* @param dt Time step (s).
	* @param t Time of the particle (s).
	* @param x Computational x coordinate of the particle.
	* @param y Computational y coordinate of the particle.
	* @param z Computational z coordinate of the particle.
	* @param q Charge state of the particle.
	* @param Z Atomic number of the impurity.
	* @return {ionization probability, recombination probability}. Either
	*         can be greater than 1 if the time step is too large.
	*/
	std::pair<double, double> calc_ioniz_recomb_probs(
		const Background::Background& bkg,
		const OpenADAS& oa_ioniz, const OpenADAS& oa_recomb,
		const double dt, const double t, const double x,
		const double y, const double z, const double q, const double Z)
	{
		// Interpolate at the particle location
		double local_ne {bkg.interp_ne(t, x, y, z)};
		double local_te {bkg.interp_te(t, x, y, z)};

		// Branchless masks
		const bool can_ionize {q < Z};
		const bool can_recomb {q > 0};

		// Ionization rate coefficients are indexed by charge. This is
		// because the zeroeth charge index in the underlying rate data is
		// for neutral ionization (charge = 0), W0 --> W1+.
		const double ioniz_rate {can_ionize ?
			oa_ioniz.get_rate_coeff(q, local_ne, local_te) : 0.0};

		// Recombination rate coefficients are indexed by charge-1. This is
		// because this zeroeth entry is for W1+ --> W0. So if we want that
		// rate coefficient for, say, W1+, we need to pass it charge-1 so it
		// chooses that zeroeth index.
		const double recomb_rate {can_recomb ?
			oa_recomb.get_rate_coeff(q-1, local_ne, local_te) : 0.0};

		// The probability of either ionization or recombination occurring is:
		//   prob = rate [m3/s] * ne [m-3] * dt [s]
		const double ne_dt {local_ne * dt};
		double ioniz_prob {ioniz_rate * ne_dt};
		double recomb_prob {recomb_rate * ne_dt};

		return std::make_pair(ioniz_prob, recomb_prob);
	}


	/**
	* Updates the charge state of every particle due to ionization and
	* recombination.
	*
	* For each particle, a random number is drawn for each process and the
	* process occurs if it is less than its probability from
	* calc_ioniz_recomb_probs. If both occur they cancel out and the charge
	* is unchanged. Probabilities greater than 1 mean a smaller time step
	* should be used if they happen often, so they are counted in
	* ioniz_warnings and recomb_warnings.
	*
	* @param slots Particle slots. The charges are updated in place.
	* @param bkg Background plasma.
	* @param oa_ioniz OpenADAS ionization (scd) rate coefficients.
	* @param oa_recomb OpenADAS recombination (acd) rate coefficients.
	* @param dt Time step (s).
	* @param ioniz_warnings Incremented by the number of particles with an
	*        ionization probability greater than 1.
	* @param recomb_warnings Incremented by the number of particles with a
	*        recombination probability greater than 1.
	* @param rngs Random number generators, one per OpenMP thread.
	*/
	void ioniz_recomb(Slots::Slots& slots,
		const Background::Background& bkg, const OpenADAS& oa_ioniz,
		const OpenADAS& oa_recomb, const double dt,
		int& ioniz_warnings, int& recomb_warnings, std::vector<pcg32>& rngs)
	{
		// Atomic number
		int Z {slots.Z()};

		#pragma omp parallel
		{
			// Grab our RNG for this thread, each thread has it own (and they're
			// seeded uniquely).
			int tid = omp_get_thread_num();
			pcg32& rng = rngs[tid];

			#pragma omp for reduction(+:ioniz_warnings, recomb_warnings)
			for (int i=0; i < slots.N(); ++i)
			{
				// Skip dead particles
				if (slots.state()[i] > 0) continue;

				// Pull into local variables
				int q {slots.q()[i]};
				double t {slots.t()[i]};
				double x {slots.x()[i]};
				double y {slots.y()[i]};
				double z {slots.z()[i]};

				// Get ionization/recombination probabilities.
				auto [ioniz_prob, recomb_prob] = calc_ioniz_recomb_probs(bkg,
					oa_ioniz, oa_recomb, dt, t, x, y, z, q, Z);

				// Track number of times the probabilities are greater than 1.
				// This indicates that a smaller timestep should be used if
				// there are a significant number of warnings.
				ioniz_warnings += (ioniz_prob  > 1.0);
				recomb_warnings += (recomb_prob > 1.0);

				// For each process, pull a random number. If that number is
				// less than prob, then that event occurs. If both events occur,
				// then they just cancel each other out and there's no change.
				const double r1 = rng.next_double();
				const double r2 = rng.next_double();
				const int ionize = (r1 < ioniz_prob);
				const int recomb = (r2 < recomb_prob);

				// Clever branchless way to adjust the charge
				const int dq = ionize - recomb;
				slots.set_q(i, q + dq);
			}  // i loop
		}  // omp parallel
	}  // ioniz_recomb
}  // namespace OpenADAS

