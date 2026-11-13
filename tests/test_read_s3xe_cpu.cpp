/*
* Tests for reading a SOLEDGE3X mesh (read_s3xe.cpp): the zones and their
* faces, the geometry, the magnetic field and the wall mask, on a small
* mesh with a circular cross-section written by the test (see
* s3x_test_files.h), where everything is known analytically.
*/

#include <array>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <gtest/gtest.h>
#include <string>
#include <unistd.h>
#include <vector>

#include "background.h"
#include "constants.h"
#include "options.h"
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
	// phi going clockwise seen from above, since SOLEDGE3X's (R, Z, phi)
	// is right-handed
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

	using Zones::face_absorbing, Zones::face_core, Zones::face_wall;
	EXPECT_EQ(zones.get_face(0, Zones::xmin), face_core);
	EXPECT_EQ(zones.get_face(0, Zones::xmax), 2);
	EXPECT_EQ(zones.get_face(2, Zones::xmin), 0);
	EXPECT_EQ(zones.get_face(1, Zones::xmax), face_absorbing);
	EXPECT_EQ(zones.get_face(0, Zones::ymin), 1);
	EXPECT_EQ(zones.get_face(0, Zones::ymax), 1);
	EXPECT_EQ(zones.get_face(2, Zones::ymin), face_wall);
	EXPECT_EQ(zones.get_face(2, Zones::ymax), face_wall);
	EXPECT_EQ(zones.get_face(2, Zones::xmax), face_absorbing);
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

// The wall outline is read in meters, as written, and a mesh file without
// one has no outlines
TEST_F(ReadS3XE, Walls)
{
	S3XTest::write_mesh(m_path, m_mesh);
	const auto mesh {S3X::read_mesh(m_path, m_mesh.nphi)};
	ASSERT_EQ(mesh.walls.size(), 1u);
	const auto& wall {mesh.walls[0]};
	ASSERT_EQ(wall.R.size(), static_cast<std::size_t>(m_mesh.wall_points + 1));
	ASSERT_EQ(wall.Z.size(), wall.R.size());
	const double r {m_mesh.a0 * S3XTest::wall_r(m_mesh)};
	EXPECT_DOUBLE_EQ(wall.R.front(), m_mesh.R0 + r);
	EXPECT_DOUBLE_EQ(wall.Z.front(), 0.0);
	EXPECT_NEAR(wall.R.back(), wall.R.front(), 1e-12);
	EXPECT_NEAR(wall.Z.back(), wall.Z.front(), 1e-12);

	S3XTest::TestMesh no_walls {m_mesh};
	no_walls.write_walls = false;
	S3XTest::write_mesh(m_path, no_walls);
	EXPECT_TRUE(S3X::read_mesh(m_path, m_mesh.nphi).walls.empty());
}

// The number of toroidal planes comes from the plasma file, and every zone
// must have the same
TEST_F(ReadS3XE, Nphi)
{
	const std::string plasma {temp_path("s3x_plasma")};
	S3XTest::write_nphi_only(plasma, m_mesh, m_mesh.nphi);
	EXPECT_EQ(S3X::read_nphi(plasma), m_mesh.nphi);

	S3XTest::write_nphi_only(plasma, m_mesh, m_mesh.nphi + 1);
	EXPECT_EXIT(S3X::read_nphi(plasma),
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


// -------------------------------------------------------------------------
// Plasma frames and the whole reader
// -------------------------------------------------------------------------

namespace
{
	// A directory with the test mesh, reference file and two plasma frames,
	// read with S3X::read_s3xe
	class ReadS3XEPlasma : public testing::Test
	{
	protected:
		S3XTest::TestMesh m_mesh {};
		std::filesystem::path m_dir {std::filesystem::temp_directory_path()
			/ ("flan_test_s3x_" + std::to_string(getpid()))};
		std::array<S3XTest::TestFrame, 2> m_frames {};

		void SetUp() override
		{
			std::filesystem::create_directories(m_dir);
			m_frames[0].tsimu = 1.0e7;
			m_frames[1].tsimu = 1.0e7 + 1000.0;
		}

		void TearDown() override {std::filesystem::remove_all(m_dir);}

		// Writes the files with the current m_frames
		void write_files()
		{
			const std::string mesh_path {(m_dir / "mesh_S3X.h5").string()};
			S3XTest::write_mesh(mesh_path, m_mesh);
			S3XTest::write_reference((m_dir / "refParam_raptorX.h5")
				.string(), m_mesh);
			const auto mesh {S3X::read_mesh(mesh_path, m_mesh.nphi)};
			for (int t {}; t < 2; ++t)
			{
				char name[32];
				std::snprintf(name, sizeof(name), "plasma_%05d.h5", t);
				S3XTest::write_plasma((m_dir / name).string(), m_mesh,
					mesh.zones, t, m_frames[t]);
			}
		}

		Options::Options options()
		{
			Options::Options opts {};
			opts.set_s3x_dir(m_dir.string());
			opts.set_s3x_frame_start(0);
			opts.set_s3x_frame_end(1);
			return opts;
		}

		Background::Background read()
		{
			write_files();
			return S3X::read_s3xe(options());
		}
	};

	constexpr double n0 {1e19};
	constexpr double T0 {100.0};
	constexpr double c0 {98227.0};
	constexpr double tau0 {5.23951e-9};
	constexpr double phi0 {100.0};
}

// Densities, temperatures and potential are unnormalized with the
// reference values, in every zone and frame, and so are the times. The
// masses come from the plasma file.
TEST_F(ReadS3XEPlasma, Plasma)
{
	const auto bkg {read()};
	ASSERT_EQ(bkg.get_nzones(), 3);
	EXPECT_DOUBLE_EQ(bkg.get_times()[0], 1.0e7 * tau0);
	EXPECT_DOUBLE_EQ(bkg.get_times()[1], (1.0e7 + 1000.0) * tau0);
	EXPECT_DOUBLE_EQ(bkg.get_elec_mass_amu(), 0.00054858);
	EXPECT_DOUBLE_EQ(bkg.get_ion_mass_amu(), 2.014);

	for (const auto [zone, i, j, k] : std::vector<std::array<int, 4>>
		{{0, 0, 0, 0}, {1, 6, 21, 3}, {2, 7, 39, 5}})
	for (int t {}; t < 2; ++t)
	{
		const double R {S3XTest::cell_R(m_mesh, zone, i, j)};
		const double Z {S3XTest::cell_Z(m_mesh, zone, i, j)};
		EXPECT_DOUBLE_EQ(bkg.get_ne()(zone, t, i, j, k),
			S3XTest::test_ne(R, Z, t) * n0);
		EXPECT_DOUBLE_EQ(bkg.get_te()(zone, t, i, j, k),
			S3XTest::test_te(R, Z, t) * T0);
		EXPECT_DOUBLE_EQ(bkg.get_ti()(zone, t, i, j, k),
			S3XTest::test_ti(R, Z, t) * T0);
		EXPECT_DOUBLE_EQ(bkg.get_vp()(zone, t, i, j, k),
			S3XTest::test_phi(R, Z, t) * phi0);
	}

	// Ghost cells come from the zone across the face
	EXPECT_DOUBLE_EQ(bkg.get_ne()(0, 1, 12, 5, 2), bkg.get_ne()(2, 1, 0, 5, 2));
}

// E = -grad vp = -0.5 phi0 Z_hat, to within the error of the curved grid.
// The flow is u_par b + E x B / B^2 + Ti (B x grad B) / B^3, with
// u_par = 0.3 c0.
TEST_F(ReadS3XEPlasma, EfieldAndFlows)
{
	const auto bkg {read()};
	const double E_want {0.5 * phi0};
	for (const auto [zone, i, j, k] : std::vector<std::array<int, 4>>
		{{0, 5, 10, 1}, {0, 11, 3, 4}, {2, 0, 20, 0}, {1, 6, 39, 2},
		{2, 7, 15, 3}, {1, 0, 7, 5}})
	for (int t {}; t < 2; ++t)
	{
		EXPECT_NEAR(bkg.get_eX()(zone, t, i, j, k), 0.0, 5e-3 * E_want)
			<< zone << " " << i << " " << j;
		EXPECT_NEAR(bkg.get_eY()(zone, t, i, j, k), 0.0, 5e-3 * E_want)
			<< zone << " " << i << " " << j;
		EXPECT_NEAR(bkg.get_eZ()(zone, t, i, j, k), -E_want, 5e-3 * E_want)
			<< zone << " " << i << " " << j;

		const Vec3 E {bkg.get_eX()(zone, t, i, j, k),
			bkg.get_eY()(zone, t, i, j, k), bkg.get_eZ()(zone, t, i, j, k)};
		EXPECT_DOUBLE_EQ(bkg.get_emag()(zone, t, i, j, k),
			std::hypot(E[0], E[1], E[2]));

		const Vec3 B {bkg.get_bX()(zone, t, i, j, k),
			bkg.get_bY()(zone, t, i, j, k), bkg.get_bZ()(zone, t, i, j, k)};
		const Vec3 gB {bkg.get_gradbX()(zone, t, i, j, k),
			bkg.get_gradbY()(zone, t, i, j, k),
			bkg.get_gradbZ()(zone, t, i, j, k)};
		const double Ti {bkg.get_ti()(zone, t, i, j, k)};
		const double B2 {B[0] * B[0] + B[1] * B[1] + B[2] * B[2]};
		const double Bm {std::sqrt(B2)};
		const Vec3 ExB {E[1] * B[2] - E[2] * B[1], E[2] * B[0] - E[0] * B[2],
			E[0] * B[1] - E[1] * B[0]};
		const Vec3 BxgB {B[1] * gB[2] - B[2] * gB[1],
			B[2] * gB[0] - B[0] * gB[2], B[0] * gB[1] - B[1] * gB[0]};
		const Vec3 u {bkg.get_uX()(zone, t, i, j, k),
			bkg.get_uY()(zone, t, i, j, k), bkg.get_uZ()(zone, t, i, j, k)};
		for (int c {}; c < 3; ++c)
		{
			const double want {0.3 * c0 * B[c] / Bm + ExB[c] / B2
				+ Ti * BxgB[c] / (B2 * Bm)};
			EXPECT_NEAR(u[c], want, 1e-9 * c0);
		}

		// The drifts are perpendicular to B, so u . b is u_par
		EXPECT_NEAR((u[0] * B[0] + u[1] * B[1] + u[2] * B[2]) / Bm,
			0.3 * c0, 1e-9 * c0);
	}
}

// The ghost cells in the files are the ones fill_ghosts would give, so the
// report shows no differences, along with the reference values
TEST_F(ReadS3XEPlasma, Messages)
{
	write_files();
	testing::internal::CaptureStdout();
	S3X::read_s3xe(options());
	const std::string out {testing::internal::GetCapturedStdout()};
	EXPECT_NE(out.find("n0 = 1e+19 m-3, T0 = 100 eV, c0 = 98227 m/s"),
		std::string::npos);
	EXPECT_NE(out.find("ne: edges 0, corners 0, averaged corners 0, "
		"boundaries 0"), std::string::npos) << out;
	EXPECT_NE(out.find("u_par: edges 0, corners 0, averaged corners 0, "
		"boundaries 0"), std::string::npos) << out;
}

// The ion parallel conductive heat flux comes back from the phi-face
// fluxes: the energy carried by the flow is taken out, the face constant C
// is fitted, and the result is unnormalized with n0 c0 T0 e.
TEST_F(ReadS3XEPlasma, HeatFlux)
{
	// Without face fluxes there's no heat flux, and a message says so
	write_files();
	testing::internal::CaptureStdout();
	auto bkg {S3X::read_s3xe(options())};
	std::string out {testing::internal::GetCapturedStdout()};
	EXPECT_TRUE(bkg.get_q_par_i().get_data().empty());
	EXPECT_NE(out.find("No ion fluxes on the phi faces"), std::string::npos)
		<< out;
	EXPECT_EQ(bkg.interp_q_par_i(2, bkg.get_times()[0], 3.5, 10.5, 1.0),
		0.0);

	// With them, q_par is recovered in every zone and frame, and C is
	// printed next to -2 pi c0 / a0
	for (auto& f : m_frames) f.face_fluxes = true;
	write_files();
	testing::internal::CaptureStdout();
	bkg = S3X::read_s3xe(options());
	out = testing::internal::GetCapturedStdout();
	EXPECT_NE(out.find("Face flux constant C = -1.2"), std::string::npos)
		<< out;
	EXPECT_EQ(out.find("Warning!"), std::string::npos) << out;

	const double q0 {n0 * c0 * T0 * Constants::ev_to_j};
	auto check = [&](const Background::Background& b)
	{
		for (const auto [zone, i, j, k] : std::vector<std::array<int, 4>>
			{{0, 0, 0, 0}, {1, 6, 21, 3}, {2, 3, 39, 5}, {2, 7, 10, 2}})
		for (int t {}; t < 2; ++t)
		{
			const double R {S3XTest::cell_R(m_mesh, zone, i, j)};
			const double Z {S3XTest::cell_Z(m_mesh, zone, i, j)};
			const double want {S3XTest::test_q_par(R, Z, t) * q0};
			EXPECT_NEAR(b.get_q_par_i()(zone, t, i, j, k), want,
				1e-9 * std::abs(want)) << zone << " " << i << " " << j;
		}

		// Ghost cells are filled from the zone across the face
		EXPECT_DOUBLE_EQ(b.get_q_par_i()(0, 1, 12, 5, 2),
			b.get_q_par_i()(2, 1, 0, 5, 2));

		// At a cell center, interpolating gives the cell's value
		const double x {3.5};
		const double y {10.5};
		const double zc {b.get_zones().get_z(2)[1]};
		EXPECT_NEAR(b.interp_q_par_i(2, b.get_times()[1], x, y, zc),
			b.get_q_par_i()(2, 1, 3, 10, 1),
			1e-9 * std::abs(b.get_q_par_i()(2, 1, 3, 10, 1)));
	};
	check(bkg);

	// A C far from -2 pi c0 / a0 gives a warning, but since C is fitted,
	// q_par is still right. So does the opposite sign.
	for (auto& f : m_frames) f.flux_C = 3.0 * 2.0 * S3XTest::pi * c0 / 0.5;
	write_files();
	testing::internal::CaptureStdout();
	bkg = S3X::read_s3xe(options());
	out = testing::internal::GetCapturedStdout();
	EXPECT_NE(out.find("Warning! In"), std::string::npos) << out;
	check(bkg);

	// Every frame needs them, or none
	m_frames[1].face_fluxes = false;
	write_files();
	EXPECT_EXIT(S3X::read_s3xe(options()),
		testing::ExitedWithCode(EXIT_FAILURE), "plasma_00001.h5 doesn't "
		"have the ion fluxes on the phi faces");
}

// Problems with the frames stop the run with a message saying what's
// wrong
TEST_F(ReadS3XEPlasma, BadFrames)
{
	auto expect_fail = [&](const std::string& message)
	{
		write_files();
		EXPECT_EXIT(S3X::read_s3xe(options()),
			testing::ExitedWithCode(EXIT_FAILURE), message);
	};

	// Fewer than two frames
	write_files();
	auto opts {options()};
	opts.set_s3x_frame_end(0);
	EXPECT_EXIT(S3X::read_s3xe(opts), testing::ExitedWithCode(EXIT_FAILURE),
		"at least two frames");

	// A frame that doesn't exist
	opts.set_s3x_frame_end(2);
	EXPECT_EXIT(S3X::read_s3xe(opts), testing::ExitedWithCode(EXIT_FAILURE),
		"plasma_00002.h5 does not exist");

	m_frames[1].tsimu = m_frames[0].tsimu;
	expect_fail("times don't increase");
	SetUp();

	m_frames[1].nzones_written = 4;
	expect_fail("plasma_00001.h5 has 4 zones, but the mesh has 3");
	m_frames[1].nzones_written = 3;

	m_frames[0].nghosts = 3;
	expect_fail("has 3 ghost cells on a side, but Flan expects 2");
	m_frames[0].nghosts = 2;

	m_frames[1].extra_theta = 1;
	expect_fail("/zone2/spec0/n .* has shape \\(16, 45, 10\\), but should "
		"have shape \\(16, 44, 10\\)");
}
