/**
* @file boris.cpp
*
* @brief Implementation of the Boris algorithm for particle stepping in an
* electromagentic field. See useful introduction and sample code that this
* was inspired by here: https://www.particleincell.com/2011/vxb-rotation/
*/
#include <array>

#include "background.h"
#include "constants.h"
#include "impurity.h"

namespace Boris
{

	// Calculate cross product and return array
	std::array<double, 3> cross_product(const std::array<double, 3>& a, 
		const std::array<double, 3>& b)
	{
		return 
		{ 
			a[1] * b[2] - a[2] * b[1],
			a[2] * b[0] - a[0] * b[2],
			a[0] * b[1] - a[1] * b[0] 
		};
	}

	void update_velocity(Impurity::Impurity& imp, 
		const Background::Background& bkg, const double dt, 
		const int tidx, const int xidx, const int yidx, const int zidx)
	{
		// Pull out electric and magnetic field components, and particle
		// velocity components and put into arrays for cleaner code
		std::array<double, 3> B {bkg.get_bX()(tidx, xidx, yidx, zidx), 
			bkg.get_bY()(tidx, xidx, yidx, zidx),
			bkg.get_bZ()(tidx, xidx, yidx, zidx)};
		std::array<double, 3> E {bkg.get_eX()(tidx, xidx, yidx, zidx), 
			bkg.get_eY()(tidx, xidx, yidx, zidx),
			bkg.get_eZ()(tidx, xidx, yidx, zidx)};
		std::array<double, 3> v {imp.get_vX(), imp.get_vY(), imp.get_vZ()};

		// Store in local variables for conveinence
		double q_m {imp.get_charge() * -Constants::charge_e / imp.get_mass()};

		// t vector
		std::array<double, 3> t {};
		for (int i {}; i < 3; ++i) t[i] = q_m * B[i] * 0.5 * dt;

		// Magnitude of t, squared
		double tmag2 {t[0]*t[0] + t[1]*t[1] + t[2]*t[2]};

		// s vector
		std::array<double, 3> s {};
		for (int i {}; i < 3; ++i) s[i] = 2 * t[i] / (1.0 + tmag2);

		// v minus
		std::array<double, 3> vminus {};
		for (int i {}; i < 3; ++i) vminus[i] = v[i] + q_m * E[i] * 0.5 * dt;

		// v prime
		std::array<double, 3> vprime {};
		std::array<double, 3> vminus_cross_t {cross_product(vminus, t)};
		for (int i {}; i < 3; ++i) vprime[i] = vminus[i] + vminus_cross_t[i];

		// v plus
		std::array<double, 3> vplus {};
		std::array<double ,3> vprime_cross_s {cross_product(vprime, s)};
		for (int i {}; i < 3; ++i) vplus[i] = vminus[i] + vprime_cross_s[i];

		// v n+1/2
		// Note we are storing particle velocity at the half time steps before
		// the position. I.e., xi = x(ti), vi = v(ti - dt/2). 
		imp.set_vX(vplus[0] + q_m * E[0] * 0.5 * dt);
		imp.set_vY(vplus[1] + q_m * E[1] * 0.5 * dt);
		imp.set_vZ(vplus[2] + q_m * E[2] * 0.5 * dt);
	}
}


