/**
* @file random.h
*
* @brief Header-only Mersenne Twister random number generator
*
* This file is taken from the useful website https://www.learncpp.com/cpp-tutorial/global-random-numbers-random-h/.
* It is freely redistributable, and saves a lot of headache trying to code 
* something like into Flan. The basic C++ random number libraries are 
* apparently notoriously bad, so this uses a more complex, robust RNG. Some
* small changes have been made to expand its flexibility and to adopt the 
* commenting format used in Flan.
*/
#ifndef RANDOM_MT_H
#define RANDOM_MT_H

#include <chrono>
#include <random>
#include <cmath>
#include <utility>

#include "constants.h"

/**
* @brief Namespace for Mersenne Twister random number generator.
*
* This header-only Random namespace implements a self-seeding Mersenne Twister.
* Requires C++17 or newer. It can be included into as many code files as 
* needed (The inline keyword avoids ODR violations). Freely redistributable, 
* courtesy of learncpp.com (https://www.learncpp.com/cpp-tutorial/global-random-numbers-random-h/)
*/
namespace Random
{
	/**
	* @brief Returns a seeded Mersenne Twister
	*
	* Note: we'd prefer to return a std::seed_seq (to initialize a 
	* std::mt19937), but std::seed can't be copied, so it can't be returned by 
	* value. Instead, we'll create a std::mt19937, seed it, and then return 
	* the std::mt19937 (which can be copied).
	*/
	inline std::mt19937 generate()
	{
		std::random_device rd{};

		// Create seed_seq with clock and 7 random numbers from 
		// std::random_device
		std::seed_seq ss{static_cast<std::seed_seq::result_type>(
			std::chrono::steady_clock::now().time_since_epoch().count()),
			rd(), rd(), rd(), rd(), rd(), rd(), rd()};

		return std::mt19937{ss};
	}

	/**
	* @brief Here's our global std::mt19937 object.
	*
	* Generates a seeded std::mt19937 and copies it into our global object. The
	* inline keyword means we only have one global instance for our whole 
	* program.
	*/
	inline std::mt19937 mt{ generate() };

	/**
	* @brief Generate a random int between [min, max] (inclusive)
	*/
	inline int get(int min, int max)
	{
		return std::uniform_int_distribution{min, max}(mt);
	}

	/**
	* @brief Generate a random double between [min, max] (inclusive)
	*
	* This function was added and is not part of the original learncpp example.
	*/
	inline double get(double min, double max)
	{
		return std::uniform_real_distribution{min, max}(mt);
	}

	/**
	* @brief Generate random double from a normal distribution
	*
	* This function was added and is not part of the original learncpp example.
	* It uses the Box-Muller transform, which outputs normally distributed
	* numbers two at a time.
	*
	* @param mu The mean of the distribution
	* @param sigma The standard deviation of the distribution
	*/
	inline std::pair<double, double> get_two_norm(double mu, double sigma)
	{

		constexpr double two_pi {2.0 * Constants::pi};
			
		// Generate two uniformally distributed numbers. Make sure ran1 is
		// greater than zero since we use it in cosine below. 
		double ran1 {};
		double ran2 {get(0.0, 1.0)};
		while (ran1 == 0.0)
		{
			ran1 = get(0.0, 1.0);
		}

		// Box-Muller transform. Calculate coefficient in front and then apply
		// the formulas.
		double coef {sigma * std::sqrt(-2.0 * std::log(ran1))};
		double ran1_norm {coef * std::cos(two_pi * ran2) + mu};
		double ran2_norm {coef * std::sin(two_pi * ran2) + mu};

		return std::make_pair(ran1_norm, ran2_norm);
		
	}

	// The following function templates can be used to generate random numbers
	// when min and/or max are not type int
	// See https://www.learncpp.com/cpp-tutorial/function-template-instantiation/
	// You can ignore these if you don't understand them

	/**
	* @brief Generate a random value between [min, max] (inclusive)
	*
	* * min and max have same type
	* * Return value has same type as min and max
	* * Supported types:
	* *    short, int, long, long long
	* *    unsigned short, unsigned int, unsigned long, or unsigned long long
	* Sample call: Random::get(1L, 6L);             // returns long
	* Sample call: Random::get(1u, 6u);             // returns unsigned int
	*/
	template <typename T>
	T get(T min, T max)
	{
		return std::uniform_int_distribution<T>{min, max}(mt);
	}

	/*
	* @brief Generate a random value between [min, max] (inclusive)
	*
	* * min and max can have different types
	* * Must explicitly specify return type as template type argument
	* * min and max will be converted to the return type
	* Sample call: Random::get<std::size_t>(0, 6);  // returns std::size_t
	* Sample call: Random::get<std::size_t>(0, 6u); // returns std::size_t
	* Sample call: Random::get<std::int>(0, 6u);    // returns int
	*/
	template <typename R, typename S, typename T>
	R get(S min, T max)
	{
		return get<R>(static_cast<R>(min), static_cast<R>(max));
	}

}

#endif
