/**
* @file openadas.cpp
* @brief Definitions for OpenADAS class. Contains routines related to reading
* in a storing OpenADAS data. Also contains two helper routines outside the 
* class for updating an impurity's ionization state.
*/
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <tuple>

#include "openadas.h"
#include "vectors.h"
#include "utilities.h"
#include "impurity.h"
#include "background.h"
#include "random.h"


namespace OpenADAS
{

	OpenADAS::OpenADAS(const std::string_view openadas_root, 
		const int openadas_year, const int imp_atom_num, 
		const std::string_view rate_type)
	{
		// Create full path to recombination or ionization rate files. 
		// These are the ACD and SCD files, respectively.
		std::string openadas_path {get_openadas_path(openadas_root, 
			openadas_year, imp_atom_num, rate_type)};
		std::cout << "ADAS: Loading " << openadas_path << '\n';

		// Open file
		std::ifstream openadas_stream {openadas_path};

		// Read in the rate coefficients along with all the header info. 
		// Everything is stored internally in the class. 
		read_rate_coefficients(openadas_stream);

		// Read in rate coefficients into Vector3D.
	}

	std::string OpenADAS::get_openadas_path(
		const std::string_view openadas_root, const int openadas_year, 
		const int imp_atom_num, const std::string_view rate_type)
	{
		// Check for valid entry for rate type
		if (rate_type != "scd" && rate_type != "acd")
		{
			std::cerr << "Error! rate_type = " << rate_type << ". Only " 
				<< "\"scd\" and \"acd\" are allowed.\n";
			// throw exception...
		}

		// Select element name from atomic number.
		std::string imp_name {};
		if (imp_atom_num == 2) imp_name = "he";
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
			std::cerr << "Error! imp_atom_num = " << imp_atom_num << " is not"
				<< " a valid option.\n";
			// throw exception...
		}

		// Assemble the full path. Example: /home/zamp/openadas/scd50_w.dat
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

	void OpenADAS::read_rate_coefficients(std::ifstream& openadas_stream)
	{
		// Get header info line as a string. Header contains (all on one line):
		// [atomic number] [# densities] [# temperatures] [minimum charge]
		// [maximimum charge] [element name] ... other stuff we don't need.
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

		// We need at least two ne and te values in order to do a bilinear
		// interpolation of the rate coefficients later in the program.
		if (m_ndens < 2)
		{
			std::cerr << "Error! Only " << m_ndens << " density values are in "
				<< "the OpenADAS file. The program will not work!\n";
		}
		if (m_ntemp < 2)
		{
			std::cerr << "Error! Only " << m_ndens << " temperature values " 
				<< "are in the OpenADAS file. The program will not work!\n";
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
		// There is no gaurantee that those 5 value for each Te will all be on
		// one line, they may span multiple lines, but each time values are
		// printed for a new Te it will start on a new line. This is the logic
		// I have picked up at least.

		// m_rates will be dimensions of (ncharges, ntemp, ndens). 
		int ncharges {m_charge_high - m_charge_low + 1};
		m_rates = Vectors::Vector3D<double>(ncharges, m_ntemp, m_ndens);

		for (int i {}; i < ncharges; ++i)
		{
			// First line has some indentifying info, but not needed
			std::getline(openadas_stream, line);

			for (int j {}; j < m_ntemp; ++j)
			{
				// Using a counter (k) with a while loop instead of a for loop
				// let's us handle the fact that entries at this charge and te
				// may span multiple lines and may not have the same amount of
				// entries per line. E.g., if ndens=9 there would be 
				// there would be alternating lines of 5 and then 4 values. We
				// read values until we've read ndens, irrespective of how many
				// line it takes to get there.
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
		// out of bounds. It is not a good idea to livalarly extrapolate past
		// the bounds, so we instead treat it as value = min(vec) was passed in
		// (i.e. floor value) by reassigning value to the minimum.
		if (lower == vec.begin())
		{
			value = vec[0];
			index0 = 0;
			index1 = 1;
		}

		// Similar to the previous, except in this case value is at/above the
		// maximimum value in vec.
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

	double OpenADAS::get_rate_coeff(int charge, double ne, double te)
		const
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

		// Get the value of the rate coefficients that we will be 
		// interpolating between.
		double rate0 = m_rates(charge, te0_index, ne0_index);
		double rate1 = m_rates(charge, te1_index, ne1_index);

		/*
		std::cout << "te0_index = " << te0_index << '\n';
		std::cout << "te1_index = " << te1_index << '\n';
		std::cout << "ne0_index = " << ne0_index << '\n';
		std::cout << "ne1_index = " << ne1_index << '\n';
		std::cout << "log10(te0) = " << std::log10(te0) << '\n';
		std::cout << "log10(te1) = " << std::log10(te1) << '\n';
		std::cout << "log10(ne0) = " << std::log10(ne0 * 1e-6) << '\n';
		std::cout << "log10(ne1) = " << std::log10(ne1 * 1e-6) << '\n';
		std::cout << "log10(rate0) = " << std::log10(rate0 * 1e-6) << '\n';
		std::cout << "log10(rate1) = " << std::log10(rate1 * 1e-6) << '\n';
		*/

		// Do a bilinear interpolation. In this case, x=Te, y=ne
		return Utilities::bilinear_interpolate(te0, ne0, rate0, te1, ne1, 
			rate1, te, ne);
			
	}

	std::pair<double, double> calc_ioniz_recomb_probs(
		Impurity::Impurity& imp, const Background::Background& bkg,
		const OpenADAS& oa_ioniz, const OpenADAS& oa_recomb, 
		const double imp_time_step, const int tidx, const int xidx, 
		const int yidx, const int zidx)
	{
		// We use ne and te more than once here, so to avoid indexing
		// multiple times it is cheaper to just do it once and save it in
		// a local variable.
		double local_ne = bkg.get_ne()(tidx, xidx, yidx, zidx);
		double local_te = bkg.get_te()(tidx, xidx, yidx, zidx);
		
		// Ionization rate coefficients are indexed by charge. This is 
		// because the zeroeth charge index in the underlying rate data is
		// for neutral ionization (charge = 0), W0 --> W1+.
		double ioniz_rate {};
		if (imp.get_charge() < imp.get_atom_num())
		{
			ioniz_rate = oa_ioniz.get_rate_coeff(imp.get_charge(), 
				local_ne, local_te);
		}

		// Recombination rate coefficients are indexed by charge-1. This is
		// because this zeroeth entry is for W1+ --> W0. So if we want that
		// rate coefficient for, say, W1+, we need to pass it charge-1 so it
		// chooses that zeroeth index.
		double recomb_rate {};
		if (imp.get_charge() > 0)
		{
			recomb_rate = oa_recomb.get_rate_coeff(imp.get_charge()-1, 
				local_ne, local_te);
		}

		/*
		std::cout << "-------------------------------\n";
		std::cout << "te = " << bkg.get_te()(tidx, xidx, yidx, zidx) << '\n';
		std::cout << "ne = " << bkg.get_ne()(tidx, xidx, yidx, zidx) << '\n';
		std::cout << "charge = " << imp.get_charge() << '\n';
		std::cout << "ioniz_rate =  " << ioniz_rate << '\n';
		std::cout << "recomb_rate = " << recomb_rate << '\n';
		std::cout << "-------------------------------\n";
		*/

		// The probability of either ionization or recombination occuring is:
		//   prob = rate [m3/s] * ne [m-3] * dt [s]
		double ioniz_prob {ioniz_rate * local_ne * imp_time_step};
		double recomb_prob {recomb_rate * local_ne * imp_time_step};

		return std::make_pair(ioniz_prob, recomb_prob);
	}

	void ioniz_recomb(Impurity::Impurity& imp, 
		const Background::Background& bkg, const OpenADAS& oa_ioniz, 
		const OpenADAS& oa_recomb, const double imp_time_step, 
		const int tidx, const int xidx, const int yidx, const int zidx, 
		int& ioniz_warnings, int& recomb_warnings)
	{
		// Get ionization/recombination probabilities.
		auto [ioniz_prob, recomb_prob] = calc_ioniz_recomb_probs(imp, bkg, 
			oa_ioniz, oa_recomb, imp_time_step, tidx, xidx, yidx, zidx);

		// Track number of times the probabilities are greater than 1. This
		// indicates that a smaller timestep should be used if there are a
		// significant number of warnings. 
		if (ioniz_prob > 1.0) ioniz_warnings += 1;
		if (recomb_prob > 1.0) recomb_warnings += 1;

		// For each process, pull a random number. If that number is less than
		// prob, then that event occurs. If both events occur, then they just 
		// cancel each other out and there's no change.
		if (Random::get(0.0, 1.0) < ioniz_prob)
		{
			imp.set_charge(imp.get_charge() + 1);
		}
		if (Random::get(0.0, 1.0) < recomb_prob)
		{
			imp.set_charge(imp.get_charge() - 1);
		}
	}
}

