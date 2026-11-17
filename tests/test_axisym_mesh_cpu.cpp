/*
* Tests for the geometry, magnetic field, electric field and flows shared
* by the axisymmetric multi-zone backgrounds (axisym_mesh.cpp), on a mesh
* built here: a rectangle in (R, Z) with uniform cells, one zone, periodic
* over the whole torus in z. On it the cell centers, Jacobian and the
* gradients of fields linear in R, Z and phi are known exactly, for either
* handedness of the toroidal angle.
*/

#include <array>
#include <cmath>
#include <gtest/gtest.h>
#include <vector>

#include "axisym_mesh.h"
#include "background.h"
#include "constants.h"
#include "flan_types.h"
#include "vectors.h"
#include "zones.h"

namespace
{
	// Mesh: nx by ny cells from (R0, Z0) in steps of dR and dZ, nz
	// toroidal cells over 2 pi
	constexpr int nx {6};
	constexpr int ny {5};
	constexpr int nz {12};
	constexpr int ghost {2};
	constexpr double R0 {1.5};
	constexpr double Z0 {-0.2};
	constexpr double dR {0.02};
	constexpr double dZ {0.03};

	// Toroidal field B_phi = Bt0 R0 / R, with a small B_R and B_Z
	constexpr double Bt0 {2.0};
	constexpr double BR0 {0.05};
	constexpr double BZ0 {-0.1};

	double R_corner(int i) {return R0 + i * dR;}
	double Z_corner(int j) {return Z0 + j * dZ;}
	double R_center(int i) {return R0 + (i + 0.5) * dR;}
	double Z_center(int j) {return Z0 + (j + 0.5) * dZ;}

	/**
	* The test mesh with toroidal handedness phi_sign.
	*/
	Axisym::Mesh make_mesh(double phi_sign)
	{
		Axisym::Mesh m {};
		m.phi_sign = phi_sign;
		std::vector<double> R, Z, Br, Bz, Bphi;
		for (int i {}; i <= nx; ++i)
		for (int j {}; j <= ny; ++j)
		{
			R.push_back(R_corner(i));
			Z.push_back(Z_corner(j));
		}
		for (int i {}; i < nx; ++i)
		for (int j {}; j < ny; ++j)
		{
			Br.push_back(BR0);
			Bz.push_back(BZ0);
			Bphi.push_back(Bt0 * R0 / R_center(i));
		}
		m.R_corners = {R};
		m.Z_corners = {Z};
		m.Br = {Br};
		m.Bz = {Bz};
		m.Bphi = {Bphi};
		return m;
	}

	/**
	* A background with the test mesh's zone (absorbing in x and y,
	* periodic in z), nframes time frames, and its geometry and B field.
	*/
	Background::Background make_background(double phi_sign,
		int nframes = 2)
	{
		Zones::Zones zones(ghost);
		auto [x, gx] = Axisym::uniform_grid(nx, 0.0, nx);
		auto [y, gy] = Axisym::uniform_grid(ny, 0.0, ny);
		auto [z, gz] = Axisym::uniform_grid(nz, 0.0, 2.0 * Constants::pi);
		zones.add_zone(x, y, z, gx, gy, gz);
		for (auto f : {Zones::xmin, Zones::xmax, Zones::ymin, Zones::ymax})
			zones.set_face(0, f, Zones::face_absorbing);
		zones.set_face(0, Zones::zmin, 0);
		zones.set_face(0, Zones::zmax, 0);
		zones.validate();

		Background::Background bkg {};
		std::vector<BkgFPType> times(nframes);
		for (int t {}; t < nframes; ++t) times[t] = t;
		bkg.move_into_times(times);
		bkg.move_into_zones(zones);
		const auto mesh {make_mesh(phi_sign)};
		Axisym::calc_geometry(mesh, bkg);
		Axisym::calc_bfield(mesh, nframes, bkg);
		return bkg;
	}

	/**
	* A field f(R, Z, phi) on the background's zone, with ghost cells
	* filled.
	*/
	template <typename F>
	Vectors::ZonedVector<BkgFPType> field(const Background::Background& bkg,
		int nframes, F f)
	{
		const auto& zones {bkg.get_zones()};
		auto v {Axisym::zoned(zones, nframes)};
		for (int t {}; t < nframes; ++t)
		for (int i {}; i < nx; ++i)
		for (int j {}; j < ny; ++j)
		for (int k {}; k < nz; ++k)
		{
			v(0, t, i, j, k) = f(R_center(i), Z_center(j),
				zones.get_z(0)[k]);
		}
		zones.fill_ghosts(v);
		return v;
	}

	class AxisymMeshTest : public testing::TestWithParam<double> {};
}

// Cell centers are at (R cos phi, phi_sign R sin phi, Z), the corners
// likewise, and |J| is the cell volume per unit x, y and z: dR dZ R. The
// background keeps the handedness.
TEST_P(AxisymMeshTest, Geometry)
{
	const double s {GetParam()};
	const auto bkg {make_background(s)};
	EXPECT_EQ(bkg.get_phi_sign(), s);
	const auto& zones {bkg.get_zones()};
	for (int i {}; i < nx; ++i)
	for (int j {}; j < ny; ++j)
	for (int k {}; k < nz; ++k)
	{
		const double phi {zones.get_z(0)[k]};
		const double R {R_center(i)};
		EXPECT_NEAR(bkg.get_X()(0, i, j, k), R * std::cos(phi), 1e-12);
		EXPECT_NEAR(bkg.get_Y()(0, i, j, k), s * R * std::sin(phi), 1e-12);
		EXPECT_NEAR(bkg.get_Z()(0, i, j, k), Z_center(j), 1e-12);
		EXPECT_NEAR(bkg.get_J()(0, i, j, k), dR * dZ * R, 1e-14);
	}
	for (int i {}; i <= nx; ++i)
	for (int k {}; k <= nz; ++k)
	{
		const double phi {zones.get_grid_z(0)[k]};
		EXPECT_NEAR(bkg.get_grid_Y()(0, i, 0, k),
			s * R_corner(i) * std::sin(phi), 1e-12);
	}
}

// The gradient of f = a R + b Z + c phi is a R_hat + b Z_hat
// + (c / R) phi_hat, with phi_hat = (-sin phi, phi_sign cos phi, 0).
// Central differences are exact for it. phi jumps by 2 pi across the
// periodic z face, so the cells next to it are skipped.
TEST_P(AxisymMeshTest, Gradient)
{
	const double s {GetParam()};
	const auto bkg {make_background(s)};
	constexpr double a {3.0}, b {-2.0}, c {0.7};
	const auto f {field(bkg, 1, [&](double R, double Z, double phi)
		{return a * R + b * Z + c * phi;})};
	const auto& zones {bkg.get_zones()};
	for (int i {}; i < nx; ++i)
	for (int j {}; j < ny; ++j)
	for (int k {1}; k < nz - 1; ++k)
	{
		const double phi {zones.get_z(0)[k]};
		const auto Rh {Axisym::R_hat(phi, s)};
		const auto Ph {Axisym::phi_hat(phi, s)};
		const auto g {Axisym::gradient(bkg, f, 0, 0, i, j, k)};
		const double R {R_center(i)};
		for (int d {}; d < 3; ++d)
		{
			const double want {a * Rh[d] + b * Axisym::Z_hat[d]
				+ c / R * Ph[d]};
			EXPECT_NEAR(g[d], want, 1e-9) << "i " << i << " j " << j
				<< " k " << k << " d " << d;
		}
	}
}

// B is (B_R, B_phi, B_Z) rotated to Cartesian, the same in every frame,
// and grad |B| matches the analytic gradient to second order in dR.
TEST_P(AxisymMeshTest, MagneticField)
{
	const double s {GetParam()};
	constexpr int nframes {3};
	const auto bkg {make_background(s, nframes)};
	const auto& zones {bkg.get_zones()};
	for (int t {}; t < nframes; ++t)
	for (int i {1}; i < nx - 1; ++i)
	for (int j {}; j < ny; ++j)
	for (int k {}; k < nz; ++k)
	{
		const double phi {zones.get_z(0)[k]};
		const double R {R_center(i)};
		const double Bphi {Bt0 * R0 / R};
		const auto B {Axisym::combine(BR0, Axisym::R_hat(phi, s), Bphi,
			Axisym::phi_hat(phi, s))};
		EXPECT_NEAR(bkg.get_bX()(0, t, i, j, k), B[0], 1e-12);
		EXPECT_NEAR(bkg.get_bY()(0, t, i, j, k), B[1], 1e-12);
		EXPECT_NEAR(bkg.get_bZ()(0, t, i, j, k), BZ0, 1e-12);

		// |B| = sqrt(BR0^2 + BZ0^2 + (Bt0 R0 / R)^2), so
		// d|B|/dR = -(Bt0 R0)^2 / (R^3 |B|)
		const double Bmag {std::sqrt(BR0 * BR0 + BZ0 * BZ0 + Bphi * Bphi)};
		const double dBdR {-Bphi * Bphi / (R * Bmag)};
		const auto Rh {Axisym::R_hat(phi, s)};
		EXPECT_NEAR(bkg.get_bmag()(0, t, i, j, k), Bmag, 1e-12);
		EXPECT_NEAR(bkg.get_gradbX()(0, t, i, j, k), dBdR * Rh[0],
			1e-3 * std::abs(dBdR));
		EXPECT_NEAR(bkg.get_gradbY()(0, t, i, j, k), dBdR * Rh[1],
			1e-3 * std::abs(dBdR));
		EXPECT_NEAR(bkg.get_gradbZ()(0, t, i, j, k), 0.0, 1e-12);
	}
}

// With vp = -E0 Z, E = E0 Z_hat. With no parallel flow and no ion
// temperature, the flow is E x B / B^2, which for a purely toroidal B
// points along -phi_sign R_hat: B_phi > 0 is the opposite physical
// direction for the two handednesses, so the drift flips with it.
TEST_P(AxisymMeshTest, ElectricFieldAndDrift)
{
	const double s {GetParam()};
	auto bkg {make_background(s)};
	constexpr double E0 {500.0};
	const auto vp {field(bkg, 2, [&](double, double Z, double)
		{return -E0 * Z;})};
	const auto upar {field(bkg, 2, [](double, double, double)
		{return 0.0;})};
	Axisym::calc_efield_and_flows(bkg, vp, upar);
	const auto& zones {bkg.get_zones()};
	for (int i {}; i < nx; ++i)
	for (int j {}; j < ny; ++j)
	for (int k {}; k < nz; ++k)
	{
		EXPECT_NEAR(bkg.get_eX()(0, 1, i, j, k), 0.0, 1e-9);
		EXPECT_NEAR(bkg.get_eY()(0, 1, i, j, k), 0.0, 1e-9);
		EXPECT_NEAR(bkg.get_eZ()(0, 1, i, j, k), E0, 1e-9);

		// u = E x B / B^2 with B = (BR0, Bphi, BZ0) in (R, phi, Z): the
		// R part is -s E0 Bphi / B^2 and the phi part s E0 BR0 / B^2,
		// written with Cartesian unit vectors that already carry s
		const double phi {zones.get_z(0)[k]};
		const double Bphi {Bt0 * R0 / R_center(i)};
		const double B2 {BR0 * BR0 + BZ0 * BZ0 + Bphi * Bphi};
		const auto Rh {Axisym::R_hat(phi, s)};
		const auto Ph {Axisym::phi_hat(phi, s)};
		for (int d {}; d < 3; ++d)
		{
			const double want {(-s * E0 * Bphi * Rh[d]
				+ s * E0 * BR0 * Ph[d]) / B2};
			const auto& u {d == 0 ? bkg.get_uX() : d == 1 ? bkg.get_uY()
				: bkg.get_uZ()};
			EXPECT_NEAR(u(0, 1, i, j, k), want, 1e-9) << "d " << d;
		}
	}
}

INSTANTIATE_TEST_SUITE_P(Handedness, AxisymMeshTest,
	testing::Values(1.0, -1.0));
