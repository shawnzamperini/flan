#include <vector>

#include "background.h"
#include "impurity.h"
#include "impurity_stats.h"
#include "openadas.h"
#include "options.h"

namespace VarianceReduction
{
	/**
	*
	*/
	bool check_split_particle(Impurity::Impurity& imp, const int tidx, 
		const int xidx, const int yidx, const int zidx, 
		const Options::Options& opts, const std::vector<int>& var_red_counts, 
		Impurity::Statistics& imp_stats);

	/**
	* @brief Check if particle should be split or not to help improve 
	*	statistics in low-count areas.
	*
	* @param imp Impurity object reference
	* @param tidx Time index
	* @param xidx x index
	* @param yidx y index
	* @param zidx z index
	* @param imp_stats Statistics object reference
	* @param imps Vector of impurities that still need to be followed
	* @param imp_var_reduct_min_weight Impurity weight must be larger than this
	*	to be split
	* @param imp_var_reduct_counts Vector of ints containing a number of counts
	*	at each frame below which splitting occurs (i.e., defines the threshold
	*	for a low-count region)
	* @param bkg The Background object for the simulation
	* @param oa_ioniz The OpenADAS object containing ionization rates
	* @param oa_recomb The OpenADAS object containing recombination rates
	* @param imp_time_step The transport time step
	*/
	void split_particle_main(Impurity::Impurity& imp, const int tidx, 
		const int xidx, const int yidx, const int zidx, 
		Impurity::Statistics& imp_stats, 
		std::vector<Impurity::Impurity>& imps,
		const double imp_var_reduct_min_weight, 
		const std::vector<int>& imp_var_reduct_counts, 
		const Background::Background& bkg, const OpenADAS::OpenADAS& oa_ioniz, 
		const OpenADAS::OpenADAS& oa_recomb, const double imp_time_step);

	/**
	* @brief Calculate the number of counts in a cell below which causes a 
	*	particle to be split in that cell
	*
	* @param imp_stats The Statistics object for the simulation
	* @return Returns a vector of ints containing a number of counts
	*	at each frame below which splitting occurs (i.e., defines the threshold
	*	for a low-count region). 
	* @see Returned vector used in split_particle_main() as the 
	*	imp_var_reduct_counts parameter
	*/
	std::vector<int> get_counts(Impurity::Statistics& imp_stats, 
		double modifier=1.0);

	/**
	* @brief Split and create a secondary Impurity, append it to imps
	* @param imp Impurity that is being split
	* @param imps Vector of Impuritys that the secondary Impurity is appened to
	* @param secondary_weight Weight of the secondary Impurity
	* @param ioniz Boolean controlling if the new particle is one charge state
	*	higher (true) or lower (false)
	*/
	void create_secondary(Impurity::Impurity& imp, 
		std::vector<Impurity::Impurity>& imps, const double secondary_weight);
	//	const bool ioniz);
}
