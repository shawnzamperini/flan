#include "background_device.h"
#include "pcg32.h"
#include "slots_device.h"

#include "impurity_stats.cuh"


namespace Collisions
{
	void collision_gpu(Slots::SlotsDevice& slots_d, 
		const Background::BackgroundDevice& bkg_d, const bool elec, 
		const double dt, ImpurityStats::StatisticsDevice& imp_stats_d,
		pcg32* rngs_d, const double elec_mass_amu, const double ion_mass_amu);

}  // namespace Collisions
