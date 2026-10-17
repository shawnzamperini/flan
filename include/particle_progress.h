#pragma once

/**
* @file  particle_progress.h
*
* @brief Progress messages printed while the particles are followed.
*/
#include <cmath>
#include <iostream>

/**
* Prints a progress message each time another 1 / num_prints of the
* particles has been finished.
*
* Set up with the total number of particles and print_interval, then call
* update with the number left as the main loop runs. With num_prints = 10,
* a message is printed when 90%, 80%, ..., 0% of the particles are left.
*/
struct ParticleProgress 
{
    int tot_parts;		///< Total number of particles
    int num_prints;		///< Number of messages to print over the run
    int next_index;		///< Which checkpoint we are waiting for
    double step;		///< Percentage step per print (e.g., 10% if
						///< num_prints = 10)

	/**
	* @param total Total number of particles.
	* @param prints Number of progress messages to print over the run.
	*/
    ParticleProgress(int total, int prints)
        : tot_parts(total),
          num_prints(prints),
          next_index(1),
          step(100.0 / prints)
    {}

	/**
	* Prints a message for every checkpoint passed since the last call.
	*
	* Several messages are printed at once if rem_parts dropped past more
	* than one checkpoint. Does nothing if tot_parts or num_prints isn't
	* positive.
	*
	* @param rem_parts Number of particles not finished yet, including the
	*        ones being followed.
	*/
    void update(int rem_parts)
    {
        if (tot_parts <= 0 || num_prints <= 0) return;

        double pct_left = 100.0 * rem_parts / tot_parts;

        // Print while we've crossed one or more checkpoints
        while (pct_left <= 100.0 - next_index * step && next_index <= num_prints) 
		{
            std::cout << "Progress: " << std::lround(pct_left) 
				<< "% remaining (" << rem_parts << " particles left)\n";
            next_index++;
        }
    }
};

