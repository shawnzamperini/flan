#ifndef TIMER_H
#define TIMER_H

#include <array>
#include <chrono>
#include <string>

namespace Timer
{
	// One entry per code section you want to time. Add new sections here
	// only -- nothing else needs to change to support a new section.
	enum class Section : std::size_t
	{
		Total,
		Read,
		Imp,
		Boris,
		Step,
		Coll,
		Bounds,
		FillSlots,
		FindCell,
		Record,
		AllDead,
		Save,
		Count  // sentinel -- must stay last, gives array size
	};

	// Per-section timing state. Each section owns its own start time and
	// "running" flag, so ending one section can never affect another.
	struct Accumulator
	{
		std::chrono::high_resolution_clock::time_point start{};
		double total_us {0.0};
		long long calls {0};
		bool running {false};
	};

	class Timer;

	// RAII scoped timer. Starts timing a section on construction and stops
	// it on destruction, so it can't be forgotten and is exception/early
	// return safe. Typical usage:
	//
	//   {
	//       ScopedTimer t(timer.acc(Section::Boris));
	//       push_boris(p);
	//   }
	//
	class ScopedTimer
	{
	public:
		explicit ScopedTimer(Accumulator& acc);
		~ScopedTimer();

		// Not copyable or movable -- a scoped timer represents exactly one
		// live timing interval.
		ScopedTimer(const ScopedTimer&) = delete;
		ScopedTimer& operator=(const ScopedTimer&) = delete;
		ScopedTimer(ScopedTimer&&) = delete;
		ScopedTimer& operator=(ScopedTimer&&) = delete;

	private:
		Accumulator& m_acc;
	};

	class Timer
	{
	public:
		// Access the accumulator for a given section. Use with ScopedTimer,
		// or call start()/stop() on the returned reference's section
		// manually via the Timer methods below.
		Accumulator& acc(Section s);
		const Accumulator& acc(Section s) const;

		// Manual start/stop, for cases where RAII scoping isn't convenient.
		// Prefer ScopedTimer where possible.
		void start(Section s);
		void stop(Section s);

		// Merge another Timer's accumulators into this one (e.g. combining
		// per-thread or per-rank timers). Sums every field for every
		// section -- unlike a hand-picked operator+, adding a new Section
		// automatically gets merged correctly.
		void merge(const Timer& other);

		// Print a formatted summary of all sections to stdout.
		void print_summary() const;

		// Format a duration given in microseconds as a human readable
		// string, choosing the largest sensible unit.
		static std::string format_time(double time_us);

	private:
		static std::string string_precision(double value, int precision);
		static const char* section_name(Section s);

		std::array<Accumulator, static_cast<std::size_t>(Section::Count)> m_acc{};
	};
}

#endif
