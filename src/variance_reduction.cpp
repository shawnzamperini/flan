#include <iostream>
#include <vector>
#include <algorithm>

#include "background.h"
#include "impurity.h"
#include "impurity_stats.h"
#include "random.h"
#include "openadas.h"
#include "options.h"
#include "variance_reduction.h"


namespace VarianceReduction
{
	bool important_region(const int tidx, const int xidx, const int yidx, 
		const int zidx, const std::vector<int>& var_red_counts, 
		Impurity::Statistics& imp_stats, const Options::Options& opts)
	{
		// High importance region defined as a region where the counts are
		// less than that specified in var_red_counts.
		if (opts.var_red_import_int() == 0)
		{
			if (imp_stats.get_counts()(tidx, xidx, yidx, zidx) <= 
				var_red_counts[tidx]) 
			{
				return true;
			}
			else return false;
		}

		// Important region is randomly chosen based on how far the particle
		// has traveled from its source. Can optionally limit it to distance
		// in a specific gkyl coordinate (x,y,z). 
		else if (opts.var_red_import_int() == 1)
		{
			std::cerr << "Error! var_red_import = \"exp_dist\" not yet"
				<< " supported.\n";
			return false;
		}

		return false;

	}

	bool check_split_particle(Impurity::Impurity& imp, const int tidx, 
		const int xidx, const int yidx, const int zidx, 
		const Options::Options& opts, const std::vector<int>& var_red_counts, 
		Impurity::Statistics& imp_stats)
	{
		// First check if particle's weight is higher than the user-defined
		// minimum weight. Without this we'd go on splitting forever.
		if (imp.get_weight() < opts.var_red_min_weight()) return false;

		// Variance reduction (particle splitting) is done if the particle is
		// in a high importance region.
		return important_region(tidx, xidx, yidx, zidx, var_red_counts, 
			imp_stats, opts);
	}

	void split_iz_rec(Impurity::Impurity& imp, const int tidx, const int xidx,
		const int yidx, const int zidx, const Options::Options& opts, 
		const std::vector<int>& var_red_counts, const double imp_time_step,  
		Impurity::Statistics& imp_stats, const Background::Background& bkg,
		std::vector<Impurity::Impurity>& imps, 
		const OpenADAS::OpenADAS& oa_ioniz, 
		const OpenADAS::OpenADAS& oa_recomb)
	{
		// Check if in a region where splitting should occur
		bool split_particle {check_split_particle(
			imp, tidx, xidx, yidx, zidx, opts, var_red_counts, 
			imp_stats)};

		// If yes...
		if (split_particle)
		{ 
			// ... then calculate the ionization/recombination probs
			auto [ioniz_prob, recomb_prob] = 
				OpenADAS::calc_ioniz_recomb_probs(
				imp, bkg, oa_ioniz, oa_recomb, imp_time_step, tidx, 
				xidx, yidx, zidx);

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
				// Secondary is one charge state higher
				create_secondary(imp, imps, ioniz_prob, 1);
			}
			else if (imp.get_weight() > recomb_prob &&
				imp.get_charge() > 0 && recomb_prob > 0.0)
			{
				// Secondary is one charge state lower
				create_secondary(imp, imps, recomb_prob, -1);
			}
		}
	}

	bool russian_roulette(Impurity::Impurity& imp, 
		const Options::Options& opts, const int tidx, const int xidx, 
		const int yidx, const int zidx, const std::vector<int>& var_red_counts, 
		Impurity::Statistics& imp_stats)
	{
		// First check if in low-importance region since we only want to kill
		// particles in that type of region to let CPU time be used in high
		// importance regions.
		if (important_region(tidx, xidx, yidx, zidx, var_red_counts, 
			imp_stats, opts)) return false;

		// Sample random number between 0-1. If random number is less than the
		// survival probability, particle survives and weight increases by 1/p.
		// Otherwise it is killed and false is returned.
		double ran {Random::get(0.0, 1.0)};
		if (ran < opts.var_red_rusrol_prob())
		{
			imp.set_weight(imp.get_weight() / opts.var_red_rusrol_prob());
			return true;
		}
		else 
		{
			//std::cout << "rusrol: goodbye cruel world\n";
			return false;
		}

	}

	void split_particle_main(Impurity::Impurity& imp, const int tidx, 
		const int xidx, const int yidx, const int zidx, 
		Impurity::Statistics& imp_stats, std::vector<Impurity::Impurity>& imps,
		const double var_red_min_weight, 
		const std::vector<int>& var_red_counts,
		const Background::Background& bkg, const OpenADAS::OpenADAS& oa_ioniz, 
		const OpenADAS::OpenADAS& oa_recomb, const double imp_time_step)
	{
		// First check if particle's weight is higher than the user-defined
		// minimum weight. Without this we'd go on splitting forever.
		if (imp.get_weight() < var_red_min_weight) return;

		// Variance reduction (particle splitting) is done if the particle is
		// in a low-count region, defined by whatever is passed in for
		// var_red_counts.
		if (imp_stats.get_counts()(tidx, xidx, yidx, zidx) <= 
			var_red_counts[tidx])
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
			std::cerr << "ERROR THIS FUNCTION DOESNT WORK\n";
			if ((ioniz_prob > recomb_prob || imp.get_charge() == 0) &&
				(imp.get_charge() < imp.get_atom_num()) &&
				(imp.get_weight() > ioniz_prob) &&
				ioniz_prob > 0.0)
			{
				//create_secondary(imp, imps, ioniz_prob, true);
			}
			else if (imp.get_weight() > recomb_prob &&
				imp.get_charge() > 0 && recomb_prob > 0.0)
			{
				//create_secondary(imp, imps, recomb_prob, false);
			}
		}

	}

	std::vector<int> get_counts(Impurity::Statistics& imp_stats, 
		double modifier)
	{
		// Additional logic can be added, but for now consider low-count as
		// anything equal to or less than the median cell count at this point
		// in time (not including zeros) times some modifier value.

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

		// Calculate and save median (times some user-defined modifier). Some 
		// potential integer division here but it's fine.
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
		const int charge_increase)
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

		// Create secondary Impurity with new weight and charge.
		Impurity::Impurity secondary_imp {imp.get_t(), imp.get_x(), 
			imp.get_y(), imp.get_z(), imp.get_X(), imp.get_Y(), imp.get_Z(), 
			imp.get_vX(), imp.get_vY(), imp.get_vZ(), secondary_weight, 
			imp.get_charge() + charge_increase, imp.get_mass(), 
			imp.get_atom_num()};

		// Add to list of impurities to be followed.
		imps.push_back(secondary_imp);
		
		//std::cout << "Split: " << imp.get_weight() << " + " 
		//	<< secondary_imp.get_weight() << '\n';

	}
}
