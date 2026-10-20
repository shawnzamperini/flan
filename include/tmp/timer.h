#ifndef TIMER_H
#define TIMER_H

#include <chrono>
#include <iostream>
#include <string>

namespace Timer
{
	// For conveinence
	using timer_t = std::chrono::time_point<std::chrono::high_resolution_clock>;

	/**
	* @brief Class object to store timers and time spent in various sections
	* of the code.
	*/
	class Timer
	{
	private:
		
		// Two timers: one to track the total time of the simulation, and 
		// another that is used to update the duration of various, 
		// non-overlapping subprocesses.
		timer_t m_total_timer {};
		timer_t m_imp_timer {};
		timer_t m_timer {};

		// Variables to accumulate time spent in various sections of the code
		double m_total_time {};
		double m_read_time {};
		double m_imp_time {};
		double m_boris_time {};
		double m_find_cell_time {};
		double m_bounds_time {};
		double m_fill_slots_time {};
		double m_record_time {};
		double m_all_dead_time {};
		double m_coll_time {};
		double m_save_time {};

		// Boolean that makes sure that m_timer is only being used in one 
		// section of code at a time
		bool m_timer_on {};

		/**
		* @brief Start the multi-purpose timer	
		*/
		void start_timer();

		/**
		* @brief Calculates how long a timer has been running
		*/
		double calc_duration(timer_t& timer);

		/**
		* @brief Convert a double to a string with indicated precision
		*/
		std::string string_precision(double value, int precision);

		/**
		* @brief Convert a time in us to a more reasonable unit
		*/
		std::string format_time(double time);

	public:
		
		/**
		* @brief Start timer for total simulation time
		*/
		void start_total_timer();

		/**
		* @brief End timer for total simulation time
		*/
		void end_total_timer();

		/**
		* @brief Start timer for time spent reading background files
		*/
		void start_read_timer();

		/**
		* @brief End timer for time spent reading background files
		*/
		void end_read_timer();

		/**
		* @brief Start timer for time spent following impurities
		*/
		void start_imp_timer();

		/**
		* @brief End timer for time spent following impurities
		*/
		void end_imp_timer();

		/**
		* @brief Start timer for time spent in the boris algorithm
		*/
		void start_boris_timer();

		/**
		* @brief End timer for time spent in the boris algorithm
		*/
		void end_boris_timer();

		/**
		* @brief Start timer for time spent checking boundary conditions
		*/
		void start_bounds_timer();

		/**
		* @brief End timer for time spent checking boundary conditions
		*/
		void end_bounds_timer();

		/**
		* @brief Start timer for time spent filling slots
		*/
		void start_fill_slots_timer();

		/**
		* @brief End timer for time spent filling slots
		*/
		void end_fill_slots_timer();

		/**
		* @brief Start timer for time spent in collisions
		*/
		void start_coll_timer();

		/**
		* @brief End timer for time spent in collisions
		*/
		void end_coll_timer();

		/**
		* @brief Start timer for time spent in finding cell indices
		*/
		void start_find_cell_timer();

		/**
		* @brief End timer for time spent in finding cell indices
		*/
		void end_find_cell_timer();

		/**
		* @brief Start timer for time spent recording statistics
		*/
		void start_record_timer();

		/**
		* @brief End timer for time spent in recording statistics
		*/
		void end_record_timer();

		/**
		* @brief Start timer for time spent in all_dead
		*/
		void start_all_dead_timer();

		/**
		* @brief End timer for time spent in all_dead
		*/
		void end_all_dead_timer();

		/**
		* @brief Start timer for time spent saving data
		*/
		void start_save_timer();

		/**
		* @brief End timer for time spent saving data
		*/
		void end_save_timer();

		/**
		* @brief Print a summary of how long simulation took in each tracked
		* section.
		*/
		void print_summary();

		/**
		* @brief operator+ overload to sum the various times spent in different
		* regions, but ONLY in parallelized regions!
		*/
		Timer operator+ (const Timer& other) const;
	};

}

#endif
