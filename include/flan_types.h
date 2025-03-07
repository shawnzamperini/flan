#ifndef FLAN_TYPES_H
#define FLAN_TYPES_H
/**
* @file flan_types.h
*
* @brief Contains aliases for some of the messier types used in Flan.
*/

#include <functional>
#include <map>
#include <string>
#include <tuple>
#include <variant>


/**
* @brief Alias for function pointer to mapc2p function.
*
* Alias for a function pointer that returns a tuple of three doubles and 
* accepts three doubles as input. Mainly used for mapc2p.
*/
using Mapc2p_ptr = std::function<std::tuple<double, double, double>
	(double, double, double)>;

/**
* @brief Alias for input options map
*
* Define type alias for input option map. Options are either an int, double
* or string. This is essentially a python dictionary if you are familiar with
* those. 
*/
using Inputs = std::map<std::string, std::variant<int, double, std::string, 
	Mapc2p_ptr>>;

#endif
