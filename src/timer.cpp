#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "timer.h"

namespace Timer
{
	// ---------------------------------------------------------------
	// ScopedTimer
	// ---------------------------------------------------------------

	ScopedTimer::ScopedTimer(Accumulator& acc)
		: m_acc(acc)
	{
		if (m_acc.running)
		{
			std::cerr << "Error! Timer is already in use. This is a "
				<< "programming error.\n";
			return;
		}
		m_acc.running = true;
		m_acc.start = std::chrono::high_resolution_clock::now();
	}

	ScopedTimer::~ScopedTimer()
	{
		// If start() above bailed out due to a re-entrancy error, don't
		// record a bogus (or negative) duration.
		if (!m_acc.running) return;

		auto dur {std::chrono::high_resolution_clock::now() - m_acc.start};
		m_acc.total_us += std::chrono::duration<double, std::micro>(dur).count();
		m_acc.calls++;
		m_acc.running = false;
	}

	// ---------------------------------------------------------------
	// Timer
	// ---------------------------------------------------------------

	Accumulator& Timer::acc(Section s)
	{
		return m_acc[static_cast<std::size_t>(s)];
	}

	const Accumulator& Timer::acc(Section s) const
	{
		return m_acc[static_cast<std::size_t>(s)];
	}

	void Timer::start(Section s)
	{
		Accumulator& a {acc(s)};
		if (a.running)
		{
			std::cerr << "Error! Timer for section '" << section_name(s)
				<< "' is already in use. This is a programming error.\n";
			return;
		}
		a.running = true;
		a.start = std::chrono::high_resolution_clock::now();
	}

	void Timer::stop(Section s)
	{
		Accumulator& a {acc(s)};
		if (!a.running)
		{
			std::cerr << "Error! Timer for section '" << section_name(s)
				<< "' was not started. This is a programming error.\n";
			return;
		}
		auto dur {std::chrono::high_resolution_clock::now() - a.start};
		a.total_us += std::chrono::duration<double, std::micro>(dur).count();
		a.calls++;
		a.running = false;
	}

	void Timer::merge(const Timer& other)
	{
		// Sum every field for every section. Because this loops over the
		// full array instead of naming fields individually, adding a new
		// Section can never be silently forgotten here (this was the bug
		// in the old hand-written operator+).
		for (std::size_t i = 0; i < m_acc.size(); ++i)
		{
			m_acc[i].total_us += other.m_acc[i].total_us;
			m_acc[i].calls += other.m_acc[i].calls;
			// Deliberately don't merge `running`/`start` -- those describe
			// live, in-progress timing state and don't make sense summed
			// across Timer instances.
		}
	}

	std::string Timer::string_precision(double value, int precision)
	{
		std::ostringstream oss;
		oss << std::fixed << std::setprecision(precision) << value;
		return oss.str();
	}

	std::string Timer::format_time(double time_us)
	{
		// Choose the largest unit that gives a value >= 1.0
		if (time_us / (3600.0 * 1e6) >= 1.0)
		{
			return string_precision(time_us / (3600.0 * 1e6), 2) + " hrs";
		}
		else if (time_us / (60.0 * 1e6) >= 1.0)
		{
			return string_precision(time_us / (60.0 * 1e6), 2) + " min";
		}
		else if (time_us / 1e6 >= 1.0)
		{
			return string_precision(time_us / 1e6, 2) + " s";
		}
		else if (time_us / 1e3 >= 1.0)
		{
			return string_precision(time_us / 1e3, 2) + " ms";
		}
		else
		{
			return string_precision(time_us, 4) + " us";
		}
	}

	const char* Timer::section_name(Section s)
	{
		switch (s)
		{
			case Section::Total:     return "Total";
			case Section::Read:      return "Read background";
			case Section::Imp:       return "Following impurities";
			case Section::Boris:     return "Boris algorithm";
			case Section::Step:      return "Position update";
			case Section::Coll:      return "Collisions";
			case Section::Bounds:    return "Boundary conditions";
			case Section::FillSlots: return "Filling slots";
			case Section::FindCell:  return "Updating cell indices";
			case Section::Record:    return "Recording statistics";
			case Section::AllDead:   return "Checking particle state";
			case Section::Save:      return "Saving data";
			default:                 return "Unknown";
		}
	}

	void Timer::print_summary() const
	{
		// (section, indent) pairs in the order/nesting we want to print.
		const std::vector<std::pair<Section, std::string>> rows
		{
			{Section::Total,     "  "},
			{Section::Read,      "    "},
			{Section::Imp,       "    "},
			{Section::Boris,     "      "},
			{Section::Coll,      "      "},
			{Section::FindCell,  "      "},
			{Section::FillSlots, "      "},
			{Section::Bounds,    "      "},
			{Section::AllDead,   "      "},
			{Section::Record,    "      "},
			{Section::Save,      "    "},
		};
 
		// Pass 1: build the label ("<indent><name>:") and time strings for
		// every row, and track the widest of each so column 2 (time) lines
		// up regardless of how long the section name or its indent is.
		std::vector<std::string> labels;
		std::vector<std::string> times;
		std::size_t label_width {0};
		std::size_t time_width {0};
 
		for (const auto& [section, indent] : rows)
		{
			const Accumulator& a {acc(section)};
 
			std::string label {indent + section_name(section) + ":"};
			label_width = std::max(label_width, label.size());
			labels.push_back(std::move(label));
 
			std::string time_str {format_time(a.total_us)};
			time_width = std::max(time_width, time_str.size());
			times.push_back(std::move(time_str));
		}
 
		// Pass 2: print with the label left-aligned and the time
		// right-aligned, so the numbers form a straight column.
		std::cout << "Time Summary (MPI Rank 0)\n";
		for (std::size_t i = 0; i < rows.size(); ++i)
		{
			const auto& [section, indent] {rows[i]};
			const Accumulator& a {acc(section)};
 
			std::cout << std::left << std::setw(static_cast<int>(label_width))
				<< labels[i] << "  "
				<< std::right << std::setw(static_cast<int>(time_width))
				<< times[i];
 
			std::cout << "  (" << a.calls << " call"
				<< (a.calls == 1 ? "" : "s");
			if (a.calls > 0)
			{
				std::cout << ", " << format_time(a.total_us / a.calls)
					<< "/call";
			}
			std::cout << ")\n";
		}
	}

}
