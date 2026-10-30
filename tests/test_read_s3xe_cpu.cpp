/*
* Tests for reading a SOLEDGE3X mesh (read_s3xe.cpp): the zones and their
* faces, the geometry, the magnetic field and the wall mask, on a small
* mesh with a circular cross-section written by the test (see
* s3x_test_files.h), where everything is known analytically.
*/

#include <array>
#include <cmath>
#include <filesystem>
#include <gtest/gtest.h>
#include <string>
#include <unistd.h>
#include <vector>

#include "background.h"
#include "read_s3xe.h"
#include "s3x_test_files.h"
#include "zones.h"


namespace
{
	using S3XTest::pi;
	using Vec3 = std::array<double, 3>;

	// Unique path in the temp directory, removed by the fixture
	std::string temp_path(const std::string& name)
	{
		return (std::filesystem::temp_directory_path() / ("flan_test_"
			+ name + "_" + std::to_string(getpid()) + ".h5")).string();
	}

	// Cylindrical unit vectors in the Cartesian frame, with SOLEDGE3X's
	// phi assumed to go clockwise seen from above
	Vec3 R_hat(double phi) {return {std::cos(phi), -std::sin(phi), 0.0};}
	Vec3 phi_hat(double phi) {return {-std::sin(phi), -std::cos(phi), 0.0};}

	class ReadS3XE : public testing::Test
	{
	protected:
		S3XTest::TestMesh m_mesh {};
		std::string m_path {temp_path("s3x_mesh")};

		void TearDown() override {std::filesystem::remove(m_path);}

		// Writes the mesh and loads it into a background with nframes
		// frames
		Background::Background load(int nframes = 2)
		{
			S3XTest::write_mesh(m_path, m_mesh);
			auto mesh {S3X::read_mesh(m_path, m_mesh.nphi)};
			Background::Background bkg {};
			std::vector<double> times(nframes);
			for (int t {}; t < nframes; ++t) times[t] = t;
			bkg.move_into_times(times);
			auto zones {mesh.zones};
			bkg.move_into_zones(zones);
			S3X::calc_geometry(mesh, bkg);
			S3X::calc_bfield(mesh, nframes, bkg);
			S3X::load_chi(mesh, bkg);
			return bkg;
		}

		// Normalized r and theta at the center of a cell, for the analytic
		// expectations
		double r_center(int zone, int i)
		{
			return 0.5 * (S3XTest::corner_r(m_mesh, zone, i)
				+ S3XTest::corner_r(m_mesh, zone, i + 1));
		}

		double theta_center(int zone, int j)
		{
			return 0.5 * (S3XTest::corner_theta(m_mesh, zone, j)
				+ S3XTest::corner_theta(m_mesh, zone, j + 1));
		}
	};
}

// The zones, their sizes and faces come from the file, including how many
// there are
TEST_F(ReadS3XE, Zones)
{
	S3XTest::write_mesh(m_path, m_mesh);
	const auto mesh {S3X::read_mesh(m_path, m_mesh.nphi)};
	const auto& zones {mesh.zones};

	ASSERT_EQ(zones.get_nzones(), 3);
	EXPECT_EQ(zones.get_ghost(), 2);
	EXPECT_EQ(zones.get_nx(0), 12);
	EXPECT_EQ(zones.get_nx(2), 8);
	EXPECT_EQ(zones.get_ny(1), 40);
	EXPECT_EQ(zones.get_nz(0), 6);
	EXPECT_DOUBLE_EQ(zones.get_x_max(2), 8.0);
	EXPECT_DOUBLE_EQ(zones.get_y_max(0), 40.0);
	EXPECT_DOUBLE_EQ(zones.get_z_max(0), 2.0 * pi);
	EXPECT_DOUBLE_EQ(zones.get_z(0)[0], 2.0 * pi / 12.0);

	using Zones::face_absorbing, Zones::face_core;
	EXPECT_EQ(zones.get_face(0, Zones::xmin), face_core);
	EXPECT_EQ(zones.get_face(0, Zones::xmax), 2);
	EXPECT_EQ(zones.get_face(2, Zones::xmin), 0);
	EXPECT_EQ(zones.get_face(1, Zones::xmax), face_absorbing);
	EXPECT_EQ(zones.get_face(0, Zones::ymin), 1);
	EXPECT_EQ(zones.get_face(0, Zones::ymax), 1);
	EXPECT_EQ(zones.get_face(2, Zones::ymin), face_absorbing);
	EXPECT_EQ(zones.get_face(1, Zones::zmin), 1);
	EXPECT_EQ(zones.get_core_zones(), (std::vector<int> {0, 1}));

	EXPECT_DOUBLE_EQ(mesh.R0, 1.7);
	EXPECT_DOUBLE_EQ(mesh.a0, 0.5);
	EXPECT_DOUBLE_EQ(mesh.B0, 2.0);
}

// Cell centers are in meters, with (R, Z, phi) right-handed
TEST_F(ReadS3XE, Positions)
{
	const auto bkg {load()};
	for (const auto [zone, i, j, k] : std::vector<std::array<int, 4>>
		{{0, 0, 0, 0}, {1, 5, 17, 3}, {2, 7, 39, 5}})
	{
		const double R {S3XTest::cell_R(m_mesh, zone, i, j)};
		const double phi {bkg.get_z(zone)[k]};
		double Z {};
		for (int di {}; di < 2; ++di)
		for (int dj {}; dj < 2; ++dj)
		{
			Z += 0.25 * m_mesh.a0 * S3XTest::corner_r(m_mesh, zone, i + di)
				* std::sin(S3XTest::corner_theta(m_mesh, zone, j + dj));
		}
		EXPECT_NEAR(bkg.get_X()(zone, i, j, k), R * std::cos(phi), 1e-12);
		EXPECT_NEAR(bkg.get_Y()(zone, i, j, k), -R * std::sin(phi), 1e-12);
		EXPECT_NEAR(bkg.get_Z()(zone, i, j, k), Z, 1e-12);
	}

	// Cell corners too
	const int nt {m_mesh.ntheta};
	EXPECT_NEAR(bkg.get_grid_X()(0, 0, 0, 0),
		m_mesh.R0 + m_mesh.a0 * m_mesh.r_in, 1e-12);
	EXPECT_NEAR(bkg.get_grid_Z()(2, 8, nt / 2, 0), m_mesh.a0 * m_mesh.r_out,
		1e-12);
}

// The reciprocal basis and Jacobian match the analytic ones of the
// circular cross-section, to within the error of the straight-sided cells:
// grad x = r_hat / (a0 dr), grad y = theta_hat / (a0 r dtheta), grad z =
// phi_hat / R and |J| = (a0 dr) (a0 r dtheta) R
TEST_F(ReadS3XE, BasisAndJacobian)
{
	const auto bkg {load()};
	const double dtheta {pi / m_mesh.ntheta};
	for (const auto [zone, i, j, k] : std::vector<std::array<int, 4>>
		{{0, 3, 7, 0}, {1, 11, 25, 2}, {2, 0, 0, 5}, {2, 4, 31, 1}})
	{
		const double dr {S3XTest::corner_r(m_mesh, zone, 1)
			- S3XTest::corner_r(m_mesh, zone, 0)};
		const double r {r_center(zone, i)};
		const double th {theta_center(zone, j)};
		const double phi {bkg.get_z(zone)[k]};
		const double R {S3XTest::cell_R(m_mesh, zone, i, j)};
		const Vec3 Rh {R_hat(phi)};
		const Vec3 r_hat {std::cos(th) * Rh[0], std::cos(th) * Rh[1],
			std::sin(th)};
		const Vec3 th_hat {-std::sin(th) * Rh[0], -std::sin(th) * Rh[1],
			std::cos(th)};
		const Vec3 ph {phi_hat(phi)};

		const std::array<Vec3, 3> want {{
			{r_hat[0] / (m_mesh.a0 * dr), r_hat[1] / (m_mesh.a0 * dr),
				r_hat[2] / (m_mesh.a0 * dr)},
			{th_hat[0] / (m_mesh.a0 * r * dtheta), th_hat[1] / (m_mesh.a0
				* r * dtheta), th_hat[2] / (m_mesh.a0 * r * dtheta)},
			{ph[0] / R, ph[1] / R, 0.0}}};
		const std::array<std::array<double, 3>, 3> got {{
			{bkg.get_dxdX()(zone, i, j, k), bkg.get_dxdY()(zone, i, j, k),
				bkg.get_dxdZ()(zone, i, j, k)},
			{bkg.get_dydX()(zone, i, j, k), bkg.get_dydY()(zone, i, j, k),
				bkg.get_dydZ()(zone, i, j, k)},
			{bkg.get_dzdX()(zone, i, j, k), bkg.get_dzdY()(zone, i, j, k),
				bkg.get_dzdZ()(zone, i, j, k)}}};

		for (int d {}; d < 3; ++d)
		{
			const double scale {std::hypot(want[d][0], want[d][1],
				want[d][2])};
			for (int c {}; c < 3; ++c)
			{
				// grad phi is exact, the others are within the chord error
				const double tol {(d == 2 ? 1e-12 : 2e-3) * scale};
				EXPECT_NEAR(got[d][c], want[d][c], tol) << "zone " << zone
					<< " cell (" << i << ", " << j << ", " << k << ") grad "
					<< "xyz"[d] << " component " << c;
			}
		}

		const double J {(m_mesh.a0 * dr) * (m_mesh.a0 * r * dtheta) * R};
		EXPECT_NEAR(bkg.get_J()(zone, i, j, k), J, 2e-3 * J);
	}
}

// Geometry ghost cells come from the zone across the face
TEST_F(ReadS3XE, GeometryGhosts)
{
	const auto bkg {load()};

	// Zone 1's psi+ ghost is zone 3's first psi cell
	EXPECT_DOUBLE_EQ(bkg.get_X()(0, 12, 9, 2), bkg.get_X()(2, 0, 9, 2));
	EXPECT_DOUBLE_EQ(bkg.get_dxdZ()(0, 13, 9, 2),
		bkg.get_dxdZ()(2, 1, 9, 2));

	// Around the ring in theta, and periodic in phi
	EXPECT_DOUBLE_EQ(bkg.get_Z()(0, 4, 40, 0), bkg.get_Z()(1, 4, 0, 0));
	EXPECT_DOUBLE_EQ(bkg.get_Y()(1, 4, -1, 0), bkg.get_Y()(0, 4, 39, 0));
	EXPECT_DOUBLE_EQ(bkg.get_Y()(1, 4, 5, -1), bkg.get_Y()(1, 4, 5, 5));
}

// B is B0 R0 / R in the phi direction, the same in every time frame, and
// grad |B| = -B0 R0 / R^2 R_hat. The gradient is a central difference
// inside the zones and across shared faces, and one-sided at the outer
// boundary.
TEST_F(ReadS3XE, MagneticField)
{
	const int nframes {3};
	const auto bkg {load(nframes)};
	const double B0R0 {m_mesh.B0 * m_mesh.R0};

	for (const auto [zone, i, j, k, tol] : std::vector<std::tuple<int, int,
		int, int, double>> {{0, 5, 10, 1, 3e-3}, {0, 11, 3, 4, 3e-3},
		{2, 0, 20, 0, 3e-3}, {1, 6, 39, 2, 3e-3}, {2, 7, 15, 3, 6e-2}})
	{
		const double R {S3XTest::cell_R(m_mesh, zone, i, j)};
		const double phi {bkg.get_z(zone)[k]};
		const Vec3 ph {phi_hat(phi)};
		const Vec3 Rh {R_hat(phi)};
		for (int t {}; t < nframes; ++t)
		{
			EXPECT_NEAR(bkg.get_bX()(zone, t, i, j, k), B0R0 / R * ph[0],
				1e-12);
			EXPECT_NEAR(bkg.get_bY()(zone, t, i, j, k), B0R0 / R * ph[1],
				1e-12);
			EXPECT_NEAR(bkg.get_bZ()(zone, t, i, j, k), 0.0, 1e-12);
			EXPECT_NEAR(bkg.get_bmag()(zone, t, i, j, k), B0R0 / R, 1e-12);

			const double grad {B0R0 / (R * R)};
			EXPECT_NEAR(bkg.get_gradbX()(zone, t, i, j, k), -grad * Rh[0],
				tol * grad) << zone << " " << i << " " << j;
			EXPECT_NEAR(bkg.get_gradbY()(zone, t, i, j, k), -grad * Rh[1],
				tol * grad) << zone << " " << i << " " << j;
			EXPECT_NEAR(bkg.get_gradbZ()(zone, t, i, j, k), 0.0,
				tol * grad) << zone << " " << i << " " << j;
		}
	}
}

// chi is loaded in every toroidal plane
TEST_F(ReadS3XE, Chi)
{
	const auto bkg {load()};
	for (int k {}; k < m_mesh.nphi; ++k)
	{
		EXPECT_DOUBLE_EQ(bkg.get_chi()(2, 7, 3, k), 1.0);
		EXPECT_DOUBLE_EQ(bkg.get_chi()(2, 6, 3, k), 0.0);
		EXPECT_DOUBLE_EQ(bkg.get_chi()(0, 11, 3, k), 0.0);
	}
}

// The number of toroidal planes comes from the plasma file, and every zone
// must have the same
TEST_F(ReadS3XE, Nphi)
{
	const std::string plasma {temp_path("s3x_plasma")};
	S3XTest::write_nphi_only(plasma, m_mesh, m_mesh.nphi);
	EXPECT_EQ(S3X::read_nphi(plasma, 3), m_mesh.nphi);

	S3XTest::write_nphi_only(plasma, m_mesh, m_mesh.nphi + 1);
	EXPECT_EXIT(S3X::read_nphi(plasma, 3),
		testing::ExitedWithCode(EXIT_FAILURE), "zone 3 has Nphi = 7");
	std::filesystem::remove(plasma);
}

// Inconsistent mesh files stop the run with a message saying what's wrong
TEST_F(ReadS3XE, BadMeshes)
{
	auto expect_fail = [&](const S3XTest::TestMesh& bad,
		const std::string& message)
	{
		S3XTest::write_mesh(m_path, bad);
		EXPECT_EXIT(S3X::read_mesh(m_path, bad.nphi),
			testing::ExitedWithCode(EXIT_FAILURE), message);
	};

	S3XTest::TestMesh bad {m_mesh};
	bad.neighbors_cols = 4;
	expect_fail(bad, "/neighbors .* has shape \\(6, 4\\), but should have "
		"shape \\(6, 3\\)");

	bad = m_mesh;
	bad.nzones_written = 4;
	expect_fail(bad, "/neighbors .* should have shape \\(6, 4\\)");

	bad = m_mesh;
	bad.skip_zone3 = true;
	expect_fail(bad, "/zone3/Npsi .* does not exist");

	bad = m_mesh;
	bad.bad_code = -5;
	expect_fail(bad, "zone 3 has neighbor code -5");

	bad = m_mesh;
	bad.mesh_nphi = 4;
	expect_fail(bad, "Only axisymmetric meshes");

	bad = m_mesh;
	bad.corner_extra = 1;
	expect_fail(bad, "/zone1/Rcorners .* has shape \\(14, 41, 1\\)");
}

// The reference values are read from their file and must agree with the
// mesh on R0, a0 and B0
TEST_F(ReadS3XE, Reference)
{
	const std::string ref_path {temp_path("s3x_ref")};
	S3XTest::write_mesh(m_path, m_mesh);
	const auto mesh {S3X::read_mesh(m_path, m_mesh.nphi)};

	S3XTest::write_reference(ref_path, m_mesh);
	const auto ref {S3X::read_reference(ref_path)};
	EXPECT_DOUBLE_EQ(ref.n0, 1e19);
	EXPECT_DOUBLE_EQ(ref.T0, 100.0);
	EXPECT_DOUBLE_EQ(ref.c0, 98227.0);
	EXPECT_DOUBLE_EQ(ref.tau0, 5.23951e-9);
	EXPECT_DOUBLE_EQ(ref.phi0, 100.0);
	S3X::check_reference(ref, mesh);

	S3XTest::write_reference(ref_path, m_mesh, 2.5);
	const auto other {S3X::read_reference(ref_path)};
	EXPECT_EXIT(S3X::check_reference(other, mesh),
		testing::ExitedWithCode(EXIT_FAILURE), "B0 is 2.5.* same case");

	S3XTest::write_reference(ref_path, m_mesh, 0.0, -1.0);
	EXPECT_EXIT(S3X::read_reference(ref_path),
		testing::ExitedWithCode(EXIT_FAILURE), "/T0 .* should be positive");
	std::filesystem::remove(ref_path);
}
