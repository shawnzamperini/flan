#include <cmath>
#include <iostream>

struct ParticleProgress 
{
    int tot_parts;
    int num_prints;
    int next_index;        // which checkpoint we are waiting for
    double step;           // percentage step per print (e.g., 10% if num_prints=10)

    ParticleProgress(int total, int prints)
        : tot_parts(total),
          num_prints(prints),
          next_index(1),
          step(100.0 / prints)
    {}

    // Call this whenever rem_parts changes
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

