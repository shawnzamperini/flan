/*
* Tests that the MPI communication in Flan works: broadcasting vectors and the
* background from rank 0, and summing the impurity statistics onto rank 0.
* Run with more than one rank, e.g., mpiexec -n 3 ./test_mpi (ctest does
* this, see tests/CMakeLists.txt).
*
* Every rank runs every test. Each test does all of its MPI calls on every
* rank no matter what, and only uses EXPECT (not ASSERT), so a failure on one
* rank can't leave the others waiting forever at a collective call. After
* each test, the MPITest fixture collects which ranks failed, so rank 0's
* summary (the one that's printed) reports a failure on any rank. Don't run
* with --gtest_shuffle, since every rank must run the tests in the same order.
*/

#include <gtest/gtest.h>
#include <cstdio>
#include <string>
#include <vector>

#include "background.h"
#include "impurity_stats.h"
#include "mpi.h"
#include "options.h"
#include "read_test.h"
#include "utilities.h"
#include "vectors.h"


namespace
{
	int mpi_rank()
	{
		int rank {};
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		return rank;
	}

	int mpi_size()
	{
		int nprocs {};
		MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
		return nprocs;
	}

	// The values rank 0 fills its vectors with, different for every element
	template <typename T>
	std::vector<T> root_values(int n)
	{
		std::vector<T> v(n);
		for (int i {}; i < n; ++i) v[i] = static_cast<T>(3 * i + 1) / 2;
		return v;
	}

	// Compare two Vector4Ds or Vector3Ds element by element, adding the name
	// of the array to any failure message. The expected array must not be
	// empty, or the comparison would pass no matter what was broadcast.
	template <typename T>
	void expect_same(const Vectors::Vector4D<T>& got,
		const Vectors::Vector4D<T>& want, const std::string& name)
	{
		EXPECT_FALSE(want.get_data().empty()) << name << " is empty, so the "
			"comparison wouldn't check anything";
		EXPECT_EQ(got.get_dim1(), want.get_dim1()) << name;
		EXPECT_EQ(got.get_dim2(), want.get_dim2()) << name;
		EXPECT_EQ(got.get_dim3(), want.get_dim3()) << name;
		EXPECT_EQ(got.get_dim4(), want.get_dim4()) << name;
		EXPECT_EQ(got.get_data(), want.get_data()) << name;
	}

	template <typename T>
	void expect_same(const Vectors::Vector3D<T>& got,
		const Vectors::Vector3D<T>& want, const std::string& name)
	{
		EXPECT_FALSE(want.get_data().empty()) << name << " is empty, so the "
			"comparison wouldn't check anything";
		EXPECT_EQ(got.get_dim1(), want.get_dim1()) << name;
		EXPECT_EQ(got.get_dim2(), want.get_dim2()) << name;
		EXPECT_EQ(got.get_dim3(), want.get_dim3()) << name;
		EXPECT_EQ(got.get_data(), want.get_data()) << name;
	}

	// Prints only failed checks, tagged with the rank. Used on every rank
	// but 0 so the output isn't repeated once per rank.
	class FailurePrinter : public testing::EmptyTestEventListener
	{
		void OnTestPartResult(const testing::TestPartResult& result) override
		{
			if (result.failed())
			{
				std::printf("[rank %d] %s:%d: Failure\n%s\n", mpi_rank(),
					result.file_name() ? result.file_name() : "unknown",
					result.line_number(), result.summary());
			}
		}
	};
}


// After each test, gathers which ranks failed onto rank 0, and fails the test
// there too if any other rank failed. Rank 0's output is the one printed, so
// without this its summary could say PASSED when another rank failed.
class MPITest : public testing::Test
{
protected:
	void TearDown() override
	{
		const int rank {mpi_rank()};
		const int nprocs {mpi_size()};
		int failed {HasFailure() ? 1 : 0};
		std::vector<int> all_failed(nprocs);
		MPI_Gather(&failed, 1, MPI_INT, all_failed.data(), 1, MPI_INT, 0,
			MPI_COMM_WORLD);

		if (rank == 0)
		{
			std::string bad_ranks {};
			for (int r {1}; r < nprocs; ++r)
				if (all_failed[r]) bad_ranks += " " + std::to_string(r);
			if (!bad_ranks.empty())
			{
				ADD_FAILURE() << "Failed on rank(s)" << bad_ranks
					<< ". See the [rank N] lines for details.";
			}
		}
	}
};


// The other tests only mean something with more than one rank
TEST_F(MPITest, RunningOnMultipleRanks)
{
	EXPECT_GE(mpi_size(), 2) << "Run this test with mpiexec -n 2 or more.";
}


// -------------------------------------------------------------------------
// Broadcasting vectors
// -------------------------------------------------------------------------

// mpi_broadcast_vector resizes the receiving vectors to match rank 0's,
// whether they start out empty, too short or too long
template <typename T>
void check_broadcast_vector()
{
	const int n {37};
	const int rank {mpi_rank()};

	std::vector<T> v {};
	if (rank == 0) v = root_values<T>(n);
	else if (rank % 2 == 1) v.resize(5, T(-1));  // too short
	else v.resize(100, T(-1));                    // too long

	Utilities::mpi_broadcast_vector(v, 0, MPI_COMM_WORLD);

	EXPECT_EQ(v, root_values<T>(n)) << "rank " << rank;
}

TEST_F(MPITest, BroadcastVectorInt) { check_broadcast_vector<int>(); }
TEST_F(MPITest, BroadcastVectorFloat) { check_broadcast_vector<float>(); }
TEST_F(MPITest, BroadcastVectorDouble) { check_broadcast_vector<double>(); }

TEST_F(MPITest, BroadcastVector4D)
{
	const int rank {mpi_rank()};
	const Vectors::Vector4D<double> want {root_values<double>(2 * 3 * 4 * 5),
		2, 3, 4, 5};

	// Other ranks start empty, the same as a Background before the broadcast
	Vectors::Vector4D<double> v {};
	if (rank == 0) v = Vectors::Vector4D<double> {want};
	v.broadcast(MPI_COMM_WORLD);

	expect_same(v, want, "Vector4D");
}

TEST_F(MPITest, BroadcastVector3D)
{
	const int rank {mpi_rank()};
	const Vectors::Vector3D<double> want {root_values<double>(3 * 4 * 5),
		3, 4, 5};

	Vectors::Vector3D<double> v {};
	if (rank == 0) v = Vectors::Vector3D<double> {want};
	v.broadcast(MPI_COMM_WORLD);

	expect_same(v, want, "Vector3D");
}


// -------------------------------------------------------------------------
// Broadcasting the background
// -------------------------------------------------------------------------

// Builds a test background with every broadcast array filled in. It starts
// from read_test, then:
// - Shifts the times and grids by 1 so that none of the simulation bounds
//   are 0. A bound that's really 0 can't be told apart from one that was
//   never broadcast, since the other ranks start at 0.
// - Fills the arrays read_test leaves empty (vp, the B gradient, the
//   physical grid edges, the metric coefficients and the tangent basis
//   vectors), each with its own pattern. Empty arrays would compare equal
//   whether they were broadcast or not.
// It's deterministic, so every rank builds exactly the same background.
Background::Background make_full_background(const std::string& test_opt)
{
	Options::Options opts {};
	opts.set_bkg_source("test");
	opts.set_test_opt(test_opt);
	Background::Background bkg {Test::read_test(opts)};

	auto shift = [](std::vector<BkgFPType> v)
	{
		for (auto& e : v) e += 1.0;
		return v;
	};
	std::vector<BkgFPType> t {shift(bkg.get_times())};
	std::vector<BkgFPType> x {shift(bkg.get_x())};
	std::vector<BkgFPType> y {shift(bkg.get_y())};
	std::vector<BkgFPType> z {shift(bkg.get_z())};
	std::vector<BkgFPType> gx {shift(bkg.get_grid_x())};
	std::vector<BkgFPType> gy {shift(bkg.get_grid_y())};
	std::vector<BkgFPType> gz {shift(bkg.get_grid_z())};
	bkg.move_into_times(t);
	bkg.move_into_x(x);
	bkg.move_into_y(y);
	bkg.move_into_z(z);
	bkg.move_into_grid_x(gx);
	bkg.move_into_grid_y(gy);
	bkg.move_into_grid_z(gz);

	const int nt {bkg.get_dim1()};
	const int nx {bkg.get_dim2()};
	const int ny {bkg.get_dim3()};
	const int nz {bkg.get_dim4()};

	// Value id + 1e-6 * i in element i, different for every array (id) and
	// every element
	auto pattern = [](int id, int n)
	{
		std::vector<BkgFPType> v(n);
		for (int i {}; i < n; ++i) v[i] = id + 1e-6 * i;
		return v;
	};
	auto v4 = [&](int id) {return Vectors::Vector4D<BkgFPType>
		{pattern(id, nt * nx * ny * nz), nt, nx, ny, nz};};
	auto v3 = [&](int id) {return Vectors::Vector3D<BkgFPType>
		{pattern(id, nx * ny * nz), nx, ny, nz};};
	auto v3_grid = [&](int id) {return Vectors::Vector3D<BkgFPType>
		{pattern(id, (nx + 1) * (ny + 1) * (nz + 1)), nx + 1, ny + 1,
		nz + 1};};

	auto vp {v4(1)};           bkg.move_into_vp(vp);
	auto gradbX {v4(2)};       bkg.move_into_gradbX(gradbX);
	auto gradbY {v4(3)};       bkg.move_into_gradbY(gradbY);
	auto gradbZ {v4(4)};       bkg.move_into_gradbZ(gradbZ);
	auto grid_X {v3_grid(5)};  bkg.move_into_grid_X(grid_X);
	auto grid_Y {v3_grid(6)};  bkg.move_into_grid_Y(grid_Y);
	auto grid_Z {v3_grid(7)};  bkg.move_into_grid_Z(grid_Z);
	auto gij_00 {v3(8)};       bkg.move_into_gij_00(gij_00);
	auto gij_01 {v3(9)};       bkg.move_into_gij_01(gij_01);
	auto gij_02 {v3(10)};      bkg.move_into_gij_02(gij_02);
	auto gij_11 {v3(11)};      bkg.move_into_gij_11(gij_11);
	auto gij_12 {v3(12)};      bkg.move_into_gij_12(gij_12);
	auto gij_22 {v3(13)};      bkg.move_into_gij_22(gij_22);
	auto dXdx {v3(14)};        bkg.move_into_dXdx(dXdx);
	auto dYdx {v3(15)};        bkg.move_into_dYdx(dYdx);
	auto dZdx {v3(16)};        bkg.move_into_dZdx(dZdx);
	auto dXdy {v3(17)};        bkg.move_into_dXdy(dXdy);
	auto dYdy {v3(18)};        bkg.move_into_dYdy(dYdy);
	auto dZdy {v3(19)};        bkg.move_into_dZdy(dZdy);
	auto dXdz {v3(20)};        bkg.move_into_dXdz(dXdz);
	auto dYdz {v3(21)};        bkg.move_into_dYdz(dYdz);
	auto dZdz {v3(22)};        bkg.move_into_dZdz(dZdz);

	return bkg;
}

// Rank 0 builds the background and broadcasts it, the same as flan() does
// with the one it reads. Every rank also builds the same background itself
// to compare against.
void check_background_broadcast(const std::string& test_opt)
{
	Background::Background bkg {};
	if (mpi_rank() == 0) bkg = make_full_background(test_opt);
	bkg.broadcast(MPI_COMM_WORLD);

	const Background::Background want {make_full_background(test_opt)};

	// Check the test itself: no bound is 0, so a missing broadcast of any
	// of them is caught
	EXPECT_NE(want.get_t_min(), 0);
	EXPECT_NE(want.get_x_min(), 0);
	EXPECT_NE(want.get_y_min(), 0);
	EXPECT_NE(want.get_z_min(), 0);

	// Dimensions and the simulation bounds. The bounds are used by the
	// boundary conditions and when starting particles, so wrong ones kill
	// every particle on that rank.
	EXPECT_EQ(bkg.get_dim1(), want.get_dim1());
	EXPECT_EQ(bkg.get_dim2(), want.get_dim2());
	EXPECT_EQ(bkg.get_dim3(), want.get_dim3());
	EXPECT_EQ(bkg.get_dim4(), want.get_dim4());
	EXPECT_EQ(bkg.get_t_min(), want.get_t_min());
	EXPECT_EQ(bkg.get_t_max(), want.get_t_max());
	EXPECT_EQ(bkg.get_x_min(), want.get_x_min());
	EXPECT_EQ(bkg.get_x_max(), want.get_x_max());
	EXPECT_EQ(bkg.get_y_min(), want.get_y_min());
	EXPECT_EQ(bkg.get_y_max(), want.get_y_max());
	EXPECT_EQ(bkg.get_z_min(), want.get_z_min());
	EXPECT_EQ(bkg.get_z_max(), want.get_z_max());

	// Times and grid
	EXPECT_EQ(bkg.get_times(), want.get_times());
	EXPECT_EQ(bkg.get_x(), want.get_x());
	EXPECT_EQ(bkg.get_y(), want.get_y());
	EXPECT_EQ(bkg.get_z(), want.get_z());
	EXPECT_EQ(bkg.get_grid_x(), want.get_grid_x());
	EXPECT_EQ(bkg.get_grid_y(), want.get_grid_y());
	EXPECT_EQ(bkg.get_grid_z(), want.get_grid_z());

	// Plasma and fields
	expect_same(bkg.get_ne(), want.get_ne(), "ne");
	expect_same(bkg.get_te(), want.get_te(), "te");
	expect_same(bkg.get_ti(), want.get_ti(), "ti");
	expect_same(bkg.get_vp(), want.get_vp(), "vp");
	expect_same(bkg.get_bX(), want.get_bX(), "bX");
	expect_same(bkg.get_bY(), want.get_bY(), "bY");
	expect_same(bkg.get_bZ(), want.get_bZ(), "bZ");
	expect_same(bkg.get_bmag(), want.get_bmag(), "bmag");
	expect_same(bkg.get_gradbX(), want.get_gradbX(), "gradbX");
	expect_same(bkg.get_gradbY(), want.get_gradbY(), "gradbY");
	expect_same(bkg.get_gradbZ(), want.get_gradbZ(), "gradbZ");
	expect_same(bkg.get_eX(), want.get_eX(), "eX");
	expect_same(bkg.get_eY(), want.get_eY(), "eY");
	expect_same(bkg.get_eZ(), want.get_eZ(), "eZ");
	expect_same(bkg.get_emag(), want.get_emag(), "emag");
	expect_same(bkg.get_uX(), want.get_uX(), "uX");
	expect_same(bkg.get_uY(), want.get_uY(), "uY");
	expect_same(bkg.get_uZ(), want.get_uZ(), "uZ");

	// Geometry
	expect_same(bkg.get_X(), want.get_X(), "X");
	expect_same(bkg.get_Y(), want.get_Y(), "Y");
	expect_same(bkg.get_Z(), want.get_Z(), "Z");
	expect_same(bkg.get_grid_X(), want.get_grid_X(), "grid_X");
	expect_same(bkg.get_grid_Y(), want.get_grid_Y(), "grid_Y");
	expect_same(bkg.get_grid_Z(), want.get_grid_Z(), "grid_Z");
	expect_same(bkg.get_J(), want.get_J(), "J");
	expect_same(bkg.get_gij_00(), want.get_gij_00(), "gij_00");
	expect_same(bkg.get_gij_01(), want.get_gij_01(), "gij_01");
	expect_same(bkg.get_gij_02(), want.get_gij_02(), "gij_02");
	expect_same(bkg.get_gij_11(), want.get_gij_11(), "gij_11");
	expect_same(bkg.get_gij_12(), want.get_gij_12(), "gij_12");
	expect_same(bkg.get_gij_22(), want.get_gij_22(), "gij_22");
	expect_same(bkg.get_dxdX(), want.get_dxdX(), "dxdX");
	expect_same(bkg.get_dxdY(), want.get_dxdY(), "dxdY");
	expect_same(bkg.get_dxdZ(), want.get_dxdZ(), "dxdZ");
	expect_same(bkg.get_dydX(), want.get_dydX(), "dydX");
	expect_same(bkg.get_dydY(), want.get_dydY(), "dydY");
	expect_same(bkg.get_dydZ(), want.get_dydZ(), "dydZ");
	expect_same(bkg.get_dzdX(), want.get_dzdX(), "dzdX");
	expect_same(bkg.get_dzdY(), want.get_dzdY(), "dzdY");
	expect_same(bkg.get_dzdZ(), want.get_dzdZ(), "dzdZ");
	expect_same(bkg.get_dXdx(), want.get_dXdx(), "dXdx");
	expect_same(bkg.get_dYdx(), want.get_dYdx(), "dYdx");
	expect_same(bkg.get_dZdx(), want.get_dZdx(), "dZdx");
	expect_same(bkg.get_dXdy(), want.get_dXdy(), "dXdy");
	expect_same(bkg.get_dYdy(), want.get_dYdy(), "dYdy");
	expect_same(bkg.get_dZdy(), want.get_dZdy(), "dZdy");
	expect_same(bkg.get_dXdz(), want.get_dXdz(), "dXdz");
	expect_same(bkg.get_dYdz(), want.get_dYdz(), "dYdz");
	expect_same(bkg.get_dZdz(), want.get_dZdz(), "dZdz");
}

// The grad-B case has a magnetic field that varies in space, the
// polarization case an electric field that varies in time, and the
// curvature case cylindrical geometry, so between them every array read_test
// sets holds values that differ from cell to cell.
TEST_F(MPITest, BroadcastBackgroundGradB)
	{ check_background_broadcast("gradb"); }
TEST_F(MPITest, BroadcastBackgroundPolarization)
	{ check_background_broadcast("polarization"); }
TEST_F(MPITest, BroadcastBackgroundCurvature)
	{ check_background_broadcast("curvature"); }


// -------------------------------------------------------------------------
// Summing the statistics
// -------------------------------------------------------------------------

// Each rank records statistics scaled by (rank + 1), with a different
// pattern for each array, and reduce_stats sums them onto rank 0. Rank 0
// should get the pattern times 1 + 2 + ... + nprocs. The patterns differ so
// that two arrays swapped in the packing would be caught.
TEST_F(MPITest, ReduceStats)
{
	const int rank {mpi_rank()};
	const int nprocs {mpi_size()};
	const int scale {rank + 1};
	const int total_scale {nprocs * (nprocs + 1) / 2};

	Impurity::Statistics stats {2, 3, 4, 5};
	const int n {static_cast<int>(stats.get_counts().get_data().size())};

	// Each array's value in cell i for a scale of 1
	auto pattern = [](int field, int i) {return 10.0 * field + i % 13 + 1;};

	// Arrays to fill, in no particular order
	std::vector<Vectors::Vector4D<BkgFPType>*> arrays {&stats.get_weights(),
		&stats.get_vX(), &stats.get_vY(), &stats.get_vZ(), &stats.get_vx(),
		&stats.get_vy(), &stats.get_vz(), &stats.get_charge(),
		&stats.get_s()};

	for (int i {}; i < n; ++i)
	{
		stats.get_counts().get_data()[i] = scale * (i % 7 + 1);
		for (std::size_t f {}; f < arrays.size(); ++f)
			arrays[f]->get_data()[i] = scale * pattern(f + 1, i);
	}

	Impurity::Statistics sum {Impurity::reduce_stats(stats)};

	// Only rank 0 gets the sum. The other ranks get an object of the same
	// shape, all zeros.
	const int count_scale {rank == 0 ? total_scale : 0};
	const double array_scale {rank == 0 ? double(total_scale) : 0.0};

	std::vector<Vectors::Vector4D<BkgFPType>*> sums {&sum.get_weights(),
		&sum.get_vX(), &sum.get_vY(), &sum.get_vZ(), &sum.get_vx(),
		&sum.get_vy(), &sum.get_vz(), &sum.get_charge(), &sum.get_s()};
	const std::vector<std::string> names {"weights", "vX", "vY", "vZ", "vx",
		"vy", "vz", "charge", "s"};

	EXPECT_EQ(sum.get_dim1(), 2);
	EXPECT_EQ(sum.get_dim2(), 3);
	EXPECT_EQ(sum.get_dim3(), 4);
	EXPECT_EQ(sum.get_dim4(), 5);
	EXPECT_EQ(static_cast<int>(sum.get_counts().get_data().size()), n);
	for (const auto* a : sums)
		EXPECT_EQ(static_cast<int>(a->get_data().size()), n);

	int bad_counts {};
	std::vector<int> bad (sums.size());
	for (int i {}; i < n; ++i)
	{
		bad_counts += sum.get_counts().get_data()[i]
			!= count_scale * (i % 7 + 1);
		for (std::size_t f {}; f < sums.size(); ++f)
			bad[f] += sums[f]->get_data()[i]
				!= array_scale * pattern(f + 1, i);
	}
	EXPECT_EQ(bad_counts, 0) << "counts on rank " << rank;
	for (std::size_t f {}; f < sums.size(); ++f)
		EXPECT_EQ(bad[f], 0) << names[f] << " on rank " << rank;
}


// gtest_main can't be used, since MPI has to be initialized first
int main(int argc, char** argv)
{
	MPI_Init(&argc, &argv);
	testing::InitGoogleTest(&argc, argv);

	// Only rank 0 prints the full output. The other ranks print just their
	// failures.
	if (mpi_rank() != 0)
	{
		testing::TestEventListeners& listeners
			{testing::UnitTest::GetInstance()->listeners()};
		delete listeners.Release(listeners.default_result_printer());
		listeners.Append(new FailurePrinter);
	}

	int result {RUN_ALL_TESTS()};

	// Every rank exits with failure if any rank failed
	int any_failed {};
	MPI_Allreduce(&result, &any_failed, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

	MPI_Finalize();
	return any_failed;
}
