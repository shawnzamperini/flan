#include <iostream>
#include <vector>
#include <algorithm>

#include "variance_reduction.h"
#include "impurity.h"
#include "impurity_stats.h"
#include "openadas.h"


namespace VarianceReduction
{

	void split_particle_main(Impurity::Impurity& imp, const int tidx, 
		const int xidx, const int yidx, const int zidx, 
		Impurity::Statistics& imp_stats, std::vector<Impurity::Impurity>& imps,
		const double imp_var_reduct_min_weight, 
		const std::vector<int>& imp_var_reduct_counts,
		const Background::Background& bkg, const OpenADAS::OpenADAS& oa_ioniz, 
		const OpenADAS::OpenADAS& oa_recomb, const double imp_time_step)
	{
		// First check if particle's weight is higher than the user-defined
		// minimum weight. Without this we'd go on splitting forever.
		if (imp.get_weight() < imp_var_reduct_min_weight) return;

		// Variance reduction (particle splitting) is done if the particle is
		// in a low-count region, defined by whatever is passed in for
		// imp_var_reduct_counts.
		if (imp_stats.get_counts()(tidx, xidx, yidx, zidx) < 
			imp_var_reduct_counts[tidx])
		{
			//std::cout << "Low-count region detected\n";
			// To split the particle, we need an outcome that has a 
			// probability of occuring. We use ionization/recombination.
			auto [ioniz_prob, recomb_prob] = calc_ioniz_recomb_probs(imp, bkg, 
				oa_ioniz, oa_recomb, imp_time_step, tidx, xidx, yidx, zidx);

			// The split particle is called a "secondary" particle. We use the
			// largest of the two possible outcomes and assign the secondary
			// weight to the probability. If particle is neutral, we just use
			// ionization, and if it's fully ionized use recombination. Care is
			// taken that we don't accidentally decrease the primary Impurity
			// weight below zero by requiring its weight is greater than the
			// weight (probability) that will be subtracted from it.
			if ((ioniz_prob > recomb_prob || imp.get_charge() == 0) &&
				(imp.get_charge() < imp.get_atom_num()) &&
				(imp.get_weight() > ioniz_prob) &&
				ioniz_prob > 0.0)
			{
				create_secondary(imp, imps, ioniz_prob, true);
			}
			else if (imp.get_weight() > recomb_prob &&
				imp.get_charge() > 0 && recomb_prob > 0.0)
			{
				create_secondary(imp, imps, recomb_prob, false);
			}
		}

	}

	std::vector<int> get_counts(Impurity::Statistics& imp_stats, 
		double modifier)
	{
		// Additional logic can be added, but for now consider low-count as
		// anything equal to or less than the median cell count at this point
		// in time (not including zeros).

		// Begin by creating a new vector of all the nonzero elements for each
		// frame.
		std::vector<int> counts (imp_stats.get_counts().get_dim1(), 0);
		for (int tidx {}; tidx < imp_stats.get_counts().get_dim1(); tidx++)
		{

		// For this frame, assemble all the nonzero counts.
		std::vector<int> frame_counts {};
		for (int xidx {}; xidx < imp_stats.get_counts().get_dim2(); xidx++)
		{
		for (int yidx {}; yidx < imp_stats.get_counts().get_dim3(); yidx++)
		{
		for (int zidx {}; zidx < imp_stats.get_counts().get_dim4(); zidx++)
		{
			int val {imp_stats.get_counts()(tidx, xidx, yidx, zidx)};		
			if (val > 0) frame_counts.push_back(val);
		}
		}
		}

		// If all zeros, array will be empty and just continue to next frame.
		std::size_t n {frame_counts.size()};
		if (n == 0) continue;

		// Sort vector
		std::sort(frame_counts.begin(), frame_counts.end());

		// Calculate and save median. Some potential integer divsion here
		// but it's fine.
		if (n % 2 == 0) 
		{
			counts[static_cast<std::size_t>(tidx)] = 
				static_cast<int>((frame_counts[n / 2 - 1] 
					+ frame_counts[n / 2]) / 2 * modifier); 
		} 
		else counts[static_cast<std::size_t>(tidx)] = static_cast<int>(
			frame_counts[n / 2] * modifier);
		}

		return counts;
	}

	void create_secondary(Impurity::Impurity& imp, 
		std::vector<Impurity::Impurity>& imps, const double secondary_weight,
		const bool ioniz)
	{
		// The primary/parent Impurity's weight must be reduced by the
		// amount that will be split off of it. If you accidentally call this
		// function with too large a secondary_weight, you've made a
		// programming error and need to fix it.
		if (imp.get_weight() > secondary_weight)
		{
			imp.set_weight(imp.get_weight() - secondary_weight);
		}
		else
		{
			std::cerr << "Error! Splitting Impurity would lead to negative "
				<< "weight.\n" << "  imp.get_weight() = " << imp.get_weight()
				<< '\n' << "  secondary_weight = " << secondary_weight << '\n';
			return;
		}

		// Secondary particle will have different charge depending on if it
		// represents a potential ionization or recombination.
		int secondary_charge {};
		if (ioniz) secondary_charge = imp.get_charge() + 1;
		else secondary_charge = imp.get_charge() - 1;

		// Create secondary Impurity with new weight and charge.
		Impurity::Impurity secondary_imp {imp.get_t(), imp.get_x(), 
			imp.get_y(), imp.get_z(), imp.get_X(), imp.get_Y(), imp.get_Z(), 
			imp.get_vX(), imp.get_vY(), imp.get_vZ(), secondary_weight, 
			secondary_charge, imp.get_mass(), imp.get_atom_num()};

		// Add to list of impurities to be followed.
		imps.push_back(secondary_imp);
		
		//std::cout << "Split: " << imp.get_weight() << " + " 
		//	<< secondary_imp.get_weight() << '\n';

	}
}
