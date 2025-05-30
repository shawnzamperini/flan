#include <iostream>
#include <tuple>
#include <cmath>

#include "flan.h"
#include "flan_types.h"
#include "constants.h"


// Namespace Mapc2p to hold all the accompanying functions that go into mapc2p.
// Copied out of the input file and modified to use the stdlib.
namespace Mapc2p
{

	// Minor radius
	double r_x(double x, double a_mid, double x_inner)
	{
		return x + a_mid - x_inner;
	}

	// Magnetic safety factor as a function of minor radius r.
	double qprofile(double r, double R_axis) 
	{
		double a[4] = {154.51071835546747,  -921.8584472748003, 
			1842.1077075366113, -1231.619813170522};
		return a[0] * std::pow(r + R_axis, 3.0) + a[1] 
			* std::pow(r + R_axis, 2.0) + a[2] * (r + R_axis) + a[3];
	}

	// Major radius as a function of minor radius r and poloidal angle theta.
	double R_rtheta(double r, double theta, double a_shift, double R_axis, 
		double delta)
	{
		return R_axis - a_shift * r * r / (2.0 * R_axis) + r 
			* std::cos(theta + std::asin(delta) * std::sin(theta));
	}

	// Z (height) as a function of minor radius r and poloidal angle theta.
	double Z_rtheta(double r, double theta, double Z_axis, double kappa)
	{
		return Z_axis + kappa * r * std::sin(theta);
	}

	// Partial derivatives of R(r,theta) and Z(r,theta)
	double dRdr(double r, double theta, double a_shift, double R_axis, 
		double delta)
	{
		return - a_shift * r / R_axis + std::cos(theta + std::asin(delta)
			* std::sin(theta));
	}
	double dRdtheta(double r, double theta, double delta)
	{
		return -r * std::sin(theta + std::asin(delta) * std::sin(theta))
			* (1.0 + std::asin(delta) * std::cos(theta));
	}
	double dZdr(double r, double theta, double kappa)
	{
		return kappa * std::sin(theta);
	}
	double dZdtheta(double r, double theta, double kappa)
	{
		return kappa * r * std::cos(theta);
	}

	// Looks like the Jacobian
	double Jr(double r, double theta, double a_shift, double R_axis, 
		double delta, double kappa)
	{
		return R_rtheta(r, theta, a_shift, R_axis, delta) 
			* (dRdr(r, theta, a_shift, R_axis, delta) 
			* dZdtheta(r, theta, kappa) - dRdtheta(r, theta, delta) 
			* dZdr(r, theta, kappa));
	}
	
	// Result from numerical quadrature. Used in below function gkyl_dbl_exp.
	struct gkyl_qr_res {
	  double res; // result of quadrature
	  double error; // error estimate
	  int status; // 0 - success.
	  int nevals; // number of function evaluations
	  int nlevels; // number of levels
	};

	// This is I think something used to efficiently compute an integral. I
	// don't really know, but it is used in the geometry equations defined 
	// here by Tess. I've modified it to use the C++ stdlib.
	struct gkyl_qr_res gkyl_dbl_exp(double (*func)(double, double, double, 
		double, double, double),
		double a, double b, int n, double eps, double r, 
		double a_shift, double R_axis, double delta, double kappa)
	{
		int nev = 0;  
		double thr = 10*std::sqrt(eps); // too generous for larger eps, e.g. eps=1e-9
		//double thr = eps; // too generous for larger eps, e.g. eps=1e-9
		double c = (a+b)/2; // center (mean)
		double d = (b-a)/2; // half distance
		double s = func(c, r, a_shift, R_axis, delta, kappa); nev += 1;
		double fp = 0, fm = 0;
		double p, e, v, h = 2;
		double tmax = std::log(2/Constants::pi * std::log((d < 1 ? 2*d : 2) / eps));
		int k = 0; // level
		do 
		{
			double q, t;
			int j = 1;
			v = s*d*Constants::pi/2*h; // last sum
			p = 0;
			h /= 2;
			t = h;
			do
			{
				double ch = std::cosh(t);
				double ecs = std::cosh(Constants::pi/2 * std::sqrt(ch*ch - 1)); // = cosh(pi/2*sinh(t))
				double w = 1/(ecs*ecs);
				double r = std::sqrt(ecs*ecs - 1)/ecs;
				double x = d*r;
				if (c+x > a) 
				{
					double y = func(c+x, r, a_shift, R_axis, delta, kappa); nev += 1;
					if (std::isfinite(y))
					fp = y;
				}
				if (c-x < b) 
				{
					double y = func(c-x, r, a_shift, R_axis, delta, kappa); nev += 1;
					if (std::isfinite(y)) 
					fm = y;
				}
				q = ch*w*(fp+fm);
				p += q;
				j += 1+(k>0);
				t = j*h;
			} 
			while (t <= tmax && std::abs(q) > eps*std::abs(p));
			s += p;
			++k;
		} 
		while (s && std::abs(2*std::abs(p) - std::abs(s)) >= std::abs(thr*s) && k <= n);
		s *= d*Constants::pi/2*h;
		e = std::abs(v-s);
		if (10*e >= std::abs(s)) 
		{
			e += std::abs(s);
			s = 0;
		}
	  
		return (struct gkyl_qr_res) 
		{
			.res = s,
			.error = e,
			.status = k>n ? 1 : 0,
			.nevals = nev,
			.nlevels = k
		};
	}

	// Integrand passed to gkyl_dbl_exp
	double integrand(double theta, double r, double a_shift, double R_axis, 
		double delta, double kappa)
	{
		return Jr(r, theta, a_shift, R_axis, delta, kappa) / 
			std::pow(R_rtheta(r, theta, a_shift, R_axis, delta), 2);
	}

	// I guess the radial derivative of the magnetic flux function
	double dPsidr(double r, double theta, double a_shift, 
		double R_axis, double delta, double kappa, double B0, double a_mid)
	{
		struct gkyl_qr_res integral;
		integral = gkyl_dbl_exp(integrand, 0.0, 2.0 * Constants::pi, 7, 1e-10,
			r, a_shift, R_axis, delta, kappa);

		return (B0 * R_axis / (2.0 * Constants::pi * qprofile(r, R_axis))) 
			* integral.res;
	}

	// Some geometric function needed for mapc2p
	double alpha(double r, double theta, double phi, double a_shift, 
		double R_axis, double delta, double kappa, double B0, double a_mid)
	{
		double twrap = theta;
		while (twrap < -Constants::pi) twrap = twrap + 2.0 *Constants::pi;
		while (Constants::pi < twrap) twrap = twrap - 2.0 * Constants::pi;

		struct gkyl_qr_res integral;
		if (0.0 < twrap) 
		{
			integral = gkyl_dbl_exp(integrand, 0.0, twrap, 7, 1e-10, r, 
				a_shift, R_axis, delta, kappa);
		} 
		else 
		{
			integral = gkyl_dbl_exp(integrand, twrap, 0.0, 7, 1e-10, r, 
				a_shift, R_axis, delta, kappa);
			integral.res = -integral.res;
		}

		return phi - B0 * R_axis * integral.res / dPsidr(r, theta, a_shift, 
			R_axis, delta, kappa, B0, a_mid);
	}
	

}


std::tuple<double, double, double> mapc2p(double xc, double yc, double zc)
{
	// Geometry and magnetic field.
	double a_shift   = 1.0;                // Parameter in Shafranov shift.
	double Z_axis    = -0.0014645315;      // Magnetic axis height [m].
	double R_axis    = 1.7074685;          // Magnetic axis major radius [m].
	double B_axis    = 2.0;                // Magnetic field at the magnetic axis [T].
	double R_LCFSmid = 2.17;               // Major radius of the LCFS at the outboard midplane [m].
	double x_inner   = 0.10;               // Radial extent inside LCFS    
	double x_outer   = 0.05;               // Radial extent outside LCFS
	double Rmid_min  = R_LCFSmid - x_inner;      // Minimum midplane major radius of simulation box [m].
	double Rmid_max  = R_LCFSmid + x_outer;      // Maximum midplane major radius of simulation box [m].
	double R0        = 0.5 * (Rmid_min + Rmid_max);  // Major radius of the simulation box [m].
	double a_mid     = R_LCFSmid - R_axis;   // Minor radius at outboard midplane [m].

	// Redefine a_mid with Shafranov shift, to ensure LCFS radial location.
	a_mid = R_axis / a_shift - std::sqrt(R_axis * (R_axis - 2.0 * a_shift
		* R_LCFSmid + 2.0 * a_shift * R_axis)) / a_shift;

	double r0        = R0 - R_axis;           // Minor radius of the simulation box [m].
	double B0        = B_axis * (R_axis / R0);  // Magnetic field magnitude in the simulation box [T].
	double kappa     = 1.35;                // Elongation (=1 for no elongation).
	double delta     = -0.4;                // Triangularity (=0 for no triangularity).
	
	// Configuration domain parameters 
	double Lx        = Rmid_max - Rmid_min;   // Domain size along x.
	double x_min     = 0.;
	double x_max     = Lx;

	// Magnetic safety factor in the center of domain.
	double q0        = Mapc2p::qprofile(Mapc2p::r_x(0.5 * (x_min + x_max), a_mid, 
		x_inner), R_axis);   

	// Minor radius
	double r = Mapc2p::r_x(xc, a_mid, x_inner);

	// Map to cylindrical (R, Z, phi) coordinates.
	double R   = Mapc2p::R_rtheta(r, zc, a_shift, R_axis, delta);
	double Z   = Mapc2p::Z_rtheta(r, zc, Z_axis, kappa);
	double phi = -q0 / r0 * yc - Mapc2p::alpha(r, zc, 0.0, a_shift, R_axis, 
		delta, kappa, B0, a_mid);
	
	// Map to Cartesian (X, Y, Z) coordinates.
	double X = R * std::cos(phi);
	double Y = R * std::sin(phi);

	//std::cout << "mapc2p: (" << xc << ", " << yc << ", " << zc << ") --> ("
	//	<< X << ", " << Y << ", " << Z << ")\n";

	// Return as three-tuple of doubles.
	return std::make_tuple(X, Y, Z);
}

Inputs create_inputs()
{
	Inputs inpts;
	
	// Int, double and string input options go here
	inpts["case_name"] = "iwl_test";
	inpts["imp_num"] = 250000;
	inpts["gkyl_dir"] = "/home/zamp/gkyldir/NT-current";
	inpts["gkyl_casename"] = "gk_d3d_negD_iwl_3x2v";
	inpts["gkyl_frame_start"] = 433;
	inpts["gkyl_frame_end"] = 533;
	inpts["gkyl_elec_name"] = "elc";
	inpts["gkyl_ion_name"] = "ion";
	inpts["imp_mass_amu"] = 183.84;
	inpts["imp_init_charge"] = 15;
	inpts["imp_time_step_opt"] = "constant";
	inpts["imp_time_step"] = 1e-9;
	inpts["imp_time_step_min"] = 1e-12;
	inpts["imp_collisions"] = "off";
	inpts["imp_vel_stats"] = "on";
	inpts["imp_xmin"] = 0.04;
	inpts["imp_xmax"] = 0.04;
	inpts["imp_xbound_buffer"] = 0.000;
	inpts["imp_ystart_opt"] = "single_value";
	inpts["imp_ystart_val"] = 0.0;
	inpts["imp_zstart_opt"] = "range";
	inpts["imp_zstart_val"] = 0.0;
	inpts["imp_iz_recomb"] = "on";
	inpts["openadas_root"] = "/home/zamp/flandir/openadas";
	inpts["openadas_year"] = 50;
	inpts["imp_var_reduct"] = "off";
	inpts["imp_var_reduct_freq"] = 0.01;
	inpts["imp_var_reduct_min_weight"] = 0.99;

	// Pointers to functions go here
	Mapc2p_ptr mapc2p_ptr {&mapc2p};

	// Function pointer input options go here
	inpts["mapc2p"] = mapc2p_ptr;

	return inpts;
}


int main()
{
	// Run Flan
	Inputs inpts {create_inputs()};
	flan(inpts);
}
