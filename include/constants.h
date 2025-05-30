/**
* @file constants.h
* @brief Some physical constants
*/
#ifndef CONSTANTS_H
#define CONSTANTS_H
namespace Constants
{
	// Mass of electron in kg
	inline constexpr double mass_e {9.109e-31};

	// Charge of electron in C
	inline constexpr double charge_e {-1.602e-19};

	// Go from eV to J
	inline constexpr double ev_to_j {1.602e-19};

	// kg per amu
	inline constexpr double amu_to_kg {1.661e-27};

	// Vacuum permittivity in F*m^-1 (i.e., C^2 * kg^-1 * m^-3 *s^2)
	inline constexpr double eps0 {8.85e-12};

	// Pi
	inline constexpr double pi {3.141592653589793};

	// A numerically small value used for checking zero equivalence
	inline constexpr double small {1.0e-50};
}
#endif
