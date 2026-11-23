/**
* @file reg_sheath_angles.cpp
*
* @brief Regression cases: the angles D+ ions hit a wall at after falling
*        through the sheath field (Borodkina and Tsvetkov's model)
*
* D+ test ions start 10 ion sound gyroradii (rho_s) from a wall, with a
* Maxwellian at Ti shifted by the sound speed c_s along B toward the wall
* (only those moving toward the wall along B, which are the ones entering
* the sheath; the rest are drawn again), and fall through the magnetic presheath and Debye sheath in a uniform
* field at an angle alpha to the wall's normal. Borodkina and Tsvetkov
* (Physics Procedia 71, 25, 2015, Fig. 2a) give the distribution of the
* angle they hit the wall at, from the normal, from their sheath model and
* from PIC codes. check_sheath_angles.py compares Flan's with theirs.
*
* This checks the "borodkina" sheath field and the push near the wall,
* sub-stepping included: the ions take the same steps as in Flan's main
* loop (Sheath::Sheath::substeps, Boris::update_velocity_cpu,
* Boris::substep_cpu, then the move), with no collisions or atomic
* physics, which the comparison leaves out.
*
* Usage: reg_sheath_angles <variant>, where <variant> is
*
*   a45   alpha = 45 degrees, Te = Ti = 30 eV, n = 1e18 m^-3, B = 5 T
*   a60   alpha = 60 degrees, Te = Ti = 30 eV, n = 1e20 m^-3, B = 3 T
*
* (Borodkina and Tsvetkov's two cases). Their paper gives the start only
* as "a shifted Maxwell velocity distribution" from 10 rho_s. The background is the slab with a
* wall at y min (slab_wall), set entirely through the slab_* options. The
* output is reg_sheath_angles_<variant>.txt in the working directory: the
* inputs, how many ions hit the wall, their mean impact angle and energy,
* the largest error in energy conservation (each ion should hit with its
* starting energy plus e (phi(start) - phi_wall)), as a fraction of the
* drop, and the number in each 1 degree bin of angle from the normal. The
* run stops when 1 in 1000 ions are left (those moving slowly along B),
* which are counted as lost.
*/
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <omp.h>

#include "background.h"
#include "boris.h"
#include "constants.h"
#include "indices.h"
#include "mpi.h"
#include "options.h"
#include "pcg32.h"
#include "read_bkg.h"
#include "sheath.h"
#include "slots.h"
#include "sputtering.h"

namespace
{
	/// One of Borodkina and Tsvetkov's cases
	struct Case
	{
		double alpha {};	///< Field angle from the wall's normal (degrees)
		double ne {};		///< Density (m^-3)
		double bmag {};		///< Field strength (T)
	};

	constexpr double te {30.0};			///< Te = Ti (eV)
	constexpr double mi_amu {2.014};	///< D
	constexpr int num_ions {20000};
	constexpr double dt {5e-10};		///< Time step (s)
	constexpr int max_steps {20000};		///< 10 us
}

int main(int argc, char* argv[])
{
	MPI_Init(&argc, &argv);

	const std::string variant {argc > 1 ? argv[1] : ""};
	Case c {};
	if (variant == "a45") c = {45.0, 1e18, 5.0};
	else if (variant == "a60") c = {60.0, 1e20, 3.0};
	else
	{
		std::cerr << "Usage: reg_sheath_angles <variant>, with <variant> "
			"a45 or a60\n";
		MPI_Finalize();
		return 1;
	}

	// The slab with a wall at y min, B at 90 - alpha degrees to the wall
	Options::Options opts {};
	opts.set_bkg_source("slab");
	opts.set_gkyl_ion_mass_amu(mi_amu);
	opts.set_gkyl_elec_mass_amu(0.000548);
	opts.set_slab_ne(c.ne);
	opts.set_slab_te(te);
	opts.set_slab_ti(te);
	opts.set_slab_bmag(c.bmag);
	opts.set_slab_efield_y(0.0);
	opts.set_slab_efield_y_rate(0.0);
	opts.set_slab_bgrad_y(0.0);
	opts.set_slab_major_radius(0.0);
	opts.set_slab_tmax(1e-3);
	opts.set_slab_wall("ymin");
	opts.set_slab_b_angle(90.0 - c.alpha);
	opts.set_imp_atom_num(1);
	opts.set_imp_mass_amu(mi_amu);
	const auto bkg {Background::read_bkg(opts)};
	const Sputtering::Surface surface(bkg);
	const Sheath::Sheath sheath(surface, bkg, Sheath::Model::borodkina);
	sheath.print_summary(std::cout);
	const auto& p {sheath.params(0, 0)};

	// The ions: 10 rho_s from the wall, Maxwellian at Ti plus c_s along
	// b, which points toward the wall (its Y part is negative)
	const double mass {mi_amu * Constants::amu_to_kg};
	const double c_s {std::sqrt(Constants::ev_to_j * 2.0 * te / mass)};
	const double v_th {std::sqrt(Constants::ev_to_j * te / mass)};
	const double a {(90.0 - c.alpha) * Constants::pi / 180.0};
	const std::array<double, 3> b {std::cos(a), -std::sin(a), 0.0};
	const auto& gx {bkg.get_grid_x(0)};
	const auto& gy {bkg.get_grid_y(0)};
	const auto& gz {bkg.get_grid_z(0)};
	const double y0 {gy.front() + 10.0 * p.rho_s};

	// Each ion's energy at the start (eV), and the potential it starts at
	std::vector<double> e_start(num_ions, 0.0);
	const double phi_start {Sheath::potential(p, Sheath::Model::borodkina,
		y0 - gy.front())};

	Slots::Slots slots(num_ions);
	slots.set_mass(mass);
	pcg32 rng(12345, 1);
	auto set_cell = [&](const int i)
	{
		slots.set_tidx(i, 0);
		slots.set_xidx(i, Indices::get_nearest_cell_index_cpu(gx,
			slots.x()[i]));
		slots.set_yidx(i, Indices::get_nearest_cell_index_cpu(gy,
			slots.y()[i]));
		slots.set_zidx(i, Indices::get_nearest_cell_index_cpu(gz,
			slots.z()[i]));
	};
	for (int i {}; i < num_ions; ++i)
	{
		slots.set_state(i, 0);
		slots.set_q(i, 1);
		slots.set_weight(i, 1.0);
		slots.set_t(i, bkg.get_times()[0]);
		slots.set_x(i, gx.front() + (0.2 + 0.6 * rng.next_double())
			* (gx.back() - gx.front()));
		slots.set_y(i, y0);
		slots.set_z(i, gz.front() + (0.2 + 0.6 * rng.next_double())
			* (gz.back() - gz.front()));
		// Only ions moving toward the wall along B enter the sheath;
		// the others are drawn again
		std::array<double, 3> v {};
		do
		{
			for (int k {}; k < 3; ++k)
				v[k] = rng.normal(0.0, v_th) + c_s * b[k];
		} while (v[0] * b[0] + v[1] * b[1] + v[2] * b[2] <= 0.0);
		slots.set_vX(i, v[0]);
		slots.set_vY(i, v[1]);
		slots.set_vZ(i, v[2]);
		e_start[i] = 0.5 * mass * (v[0] * v[0] + v[1] * v[1] + v[2] * v[2])
			/ Constants::ev_to_j;
		set_cell(i);
	}

	// Offset the velocities by half a step, as the main loop does
	const Options::Options push_opts {};
	Boris::update_velocity_cpu(slots, bkg, push_opts, -0.5 * dt, &sheath);

	std::vector<double> hits(90, 0.0);
	double angle_sum {};
	double energy_sum {};
	int n_hit {};
	double max_energy_error {};	// Most |E - E_start - (phi_start -
								// phi_wall)| over the drop
	int n_lost {};
	std::vector<int> nsub(num_ions, 1);
	int alive {num_ions};
	int step {};
	for (; step < max_steps && alive > num_ions / 1000; ++step)
	{
		// Sub-steps near the wall, the push, then the move, as in the
		// main loop
		#pragma omp parallel for
		for (int i = 0; i < num_ions; ++i)
		{
			nsub[i] = 1;
			if (slots.state()[i] > 0) continue;
			const double speed {std::sqrt(slots.vX()[i] * slots.vX()[i]
				+ slots.vY()[i] * slots.vY()[i] + slots.vZ()[i]
				* slots.vZ()[i])};
			nsub[i] = sheath.substeps(0, 0, slots.xidx()[i],
				slots.yidx()[i], slots.zidx()[i], {slots.x()[i],
				slots.y()[i], slots.z()[i]}, speed,
				Constants::ev_to_j / mass, dt);
		}
		Boris::update_velocity_cpu(slots, bkg, push_opts, dt, &sheath,
			&nsub);
		Boris::substep_cpu(slots, bkg, dt, sheath, nsub);

		for (int i = 0; i < num_ions; ++i)
		{
			if (slots.state()[i] > 0) continue;
			if (nsub[i] <= 1)
			{
				slots.set_x(i, slots.x()[i] + slots.vx()[i] * dt);
				slots.set_y(i, slots.y()[i] + slots.vy()[i] * dt);
				slots.set_z(i, slots.z()[i] + slots.vz()[i] * dt);
			}

			// Into the wall: its angle from the normal (+Y) and energy
			if (slots.y()[i] < gy.front())
			{
				const double vx {slots.vX()[i]};
				const double vy {slots.vY()[i]};
				const double vz {slots.vZ()[i]};
				const double v {std::sqrt(vx * vx + vy * vy + vz * vz)};
				const double angle {std::acos(std::min(std::abs(vy) / v,
					1.0)) * 180.0 / Constants::pi};
				const double energy {0.5 * mass * v * v / Constants::ev_to_j};
				hits[std::min(static_cast<int>(angle), 89)] += 1.0;
				angle_sum += angle;
				energy_sum += energy;
				max_energy_error = std::max(max_energy_error,
					std::abs(energy - e_start[i] - (phi_start - p.phi_wall))
					/ -p.phi_wall);
				++n_hit;
				slots.set_state(i, 1);
				--alive;
				continue;
			}

			// Back into the plasma: lost
			if (slots.y()[i] > gy.back())
			{
				++n_lost;
				slots.set_state(i, 1);
				--alive;
				continue;
			}

			// Periodic along the wall
			const double lx {gx.back() - gx.front()};
			const double lz {gz.back() - gz.front()};
			if (slots.x()[i] < gx.front()) slots.set_x(i, slots.x()[i] + lx);
			if (slots.x()[i] > gx.back()) slots.set_x(i, slots.x()[i] - lx);
			if (slots.z()[i] < gz.front()) slots.set_z(i, slots.z()[i] + lz);
			if (slots.z()[i] > gz.back()) slots.set_z(i, slots.z()[i] - lz);
			set_cell(i);
		}
	}
	n_lost += alive;

	const std::string out {"reg_sheath_angles_" + variant + ".txt"};
	std::ofstream f(out);
	f << "# D+ impact angles on a wall, Borodkina and Tsvetkov's sheath\n";
	f << "alpha " << c.alpha << "\n";
	f << "te " << te << "\n";
	f << "ti " << te << "\n";
	f << "ne " << c.ne << "\n";
	f << "bmag " << c.bmag << "\n";
	f << "rho_s " << p.rho_s << "\n";
	f << "lambda_d " << p.lambda_d << "\n";
	f << "drop " << -p.phi_wall << "\n";
	f << "num_ions " << num_ions << "\n";
	f << "hit " << n_hit << "\n";
	f << "lost " << n_lost << "\n";
	f << "steps " << step << "\n";
	f << "mean_angle " << (n_hit > 0 ? angle_sum / n_hit : 0.0) << "\n";
	f << "mean_energy " << (n_hit > 0 ? energy_sum / n_hit : 0.0) << "\n";
	f << "max_energy_error " << max_energy_error << "\n";
	f << "# angle_lo count\n";
	for (int k {}; k < 90; ++k) f << k << " " << hits[k] << "\n";
	std::cout << out << ": " << n_hit << " of " << num_ions << " hit the "
		"wall, mean angle " << (n_hit > 0 ? angle_sum / n_hit : 0.0)
		<< " degrees, mean energy " << (n_hit > 0 ? energy_sum / n_hit
		: 0.0) << " eV, energy conserved to " << max_energy_error
		<< " of the drop\n";

	MPI_Finalize();
	return 0;
}
