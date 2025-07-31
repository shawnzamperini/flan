#include <chrono>
#include <iomanip>
#include <iostream>
#include <omp.h>
#include <sstream>
#include <string>

#include "timer.h"

namespace Timer
{

	// Start the multi-purpose timer, checking that it isn't already in use.
	void Timer::start_timer()
	{
		if (m_timer_on) 
		{
			std::cerr << "Error! Timer is already in use. This is a "
				<< "programming error.\n";
		}
		else
		{
			m_timer = std::chrono::high_resolution_clock::now();
			m_timer_on = true;
		}
	}

	// Calculate duration, "turn off" timer, return duration
	double Timer::calc_duration(timer_t& timer)
	{
		// Calculate duration in seconds, store in a duration object
		auto duration {std::chrono::duration_cast<std::chrono::microseconds>
			(std::chrono::high_resolution_clock::now() - timer)};

		// With the duration calculated, the timer is effectively off and
		// free to be used elsewhere.
		m_timer_on = false;

		// Return the duration as a double
		return duration.count();
	}

	// Start/end timer for total simulation time
	void Timer::start_total_timer()
	{
		m_total_timer = std::chrono::high_resolution_clock::now();
	}
	void Timer::end_total_timer()
	{
		m_total_time += calc_duration(m_total_timer);
	}

	// Start/end timer for total time spent following impurities
	void Timer::start_imp_timer()
	{
		m_imp_timer = std::chrono::high_resolution_clock::now();
	}
	void Timer::end_imp_timer()
	{
		m_imp_time += calc_duration(m_imp_timer);
	}

	// Start/end timer for time spent reading in background files
	void Timer::start_read_timer()
	{
		start_timer();
	}
	void Timer::end_read_timer()
	{
		m_read_time += calc_duration(m_timer);
	}

	// Start/end timer for time spent in collisions
	void Timer::start_coll_timer()
	{
		start_timer();
	}
	void Timer::end_coll_timer()
	{
		m_coll_time += calc_duration(m_timer);
	}

	// Start/end timer for time spent saving data
	void Timer::start_save_timer()
	{
		start_timer();
	}
	void Timer::end_save_timer()
	{
		m_save_time += calc_duration(m_timer);
	}


	// Convert double to string with indicated precision
	std::string Timer::string_precision(double value, int precision)
	{
		std::ostringstream oss;
		oss << std::fixed << std::setprecision(precision) << value;
		return oss.str();
	}

	// Format time as a string
	std::string Timer::format_time(double time)
	{

		// Essentially choose the largest unit that is greater than 1
		// Hours
		if (time / (3600 * 1e6) > 1.0)
		{
			return string_precision(time / (3600 * 1e6), 2) + " hrs";
		}
		else if (time / (60 * 1e6) > 1.0)
		{
			return string_precision(time / (60 * 1e6), 2) + " min";
		}
		else if (time / 1e6 > 1.0)
		{
			return string_precision(time / 1e6, 2) + " s";
		}
		else if (time / 1e3 > 1.0)
		{
			return string_precision(time / 1e6, 2) + " ms";
		}
		else
		{
			return string_precision(time, 2) + " us";
		}

	}

	// Print summary in nice format
	void Timer::print_summary()
	{
		// Get how many threads are available since some tiem variables are
		// being "double counted". Example: If there are 10 threads, then
		// the actual coll_time, on average, will be coll_time / 10.
		int omp_num_threads {};
		#pragma omp parallel
		{
			if (omp_get_thread_num() == 0) 
				omp_num_threads = omp_get_num_threads();
		}

		std::cout << "Time Summary\n";
		std::cout << "  Total:  " << format_time(m_total_time) << '\n';
		std::cout << "    Calculating/reading background: " 
			<< format_time(m_read_time) << '\n';
		std::cout << "    Following impurities: " << format_time(m_imp_time)
			<< '\n';
		std::cout << "      Calculating collisions:  " 
			<< format_time(m_coll_time / omp_num_threads) << '\n';
		std::cout << "    Saving data: " << format_time(m_save_time) << '\n';
	}

	Timer Timer::operator+ (const Timer& other) const
	{
		Timer ret_timer {};

		// We assume that total timer is still running, so copy it over
		ret_timer.m_total_timer = m_total_timer;
		ret_timer.m_imp_timer = m_imp_timer;

		// These times are not from parallelized regions, so just copy over
		ret_timer.m_read_time = m_read_time;
		ret_timer.m_imp_time = m_imp_time;

		// Timers that are started/ended in parallel regions need to be 
		// aggregated. See print_summary for how we account for the fact that 
		// the time in parallel regions needs to be divided by the number of 
		// threads.
		ret_timer.m_coll_time = m_coll_time + other.m_coll_time;

		// Return combined Timer
		return ret_timer;
	}
}
