/*
* Tests for Vectors::ZonedVector, the multi-zone storage used for SOLEDGE3X
* backgrounds (and single-zone Gkeyll backgrounds).
*/

#include <gtest/gtest.h>
#include <vector>

#include "vectors.h"


// Two zones of different sizes with 2 ghost cells, like a SOLEDGE3X plasma
// file
static Vectors::ZonedVector<double> make_two_zone(int dim1 = 3)
{
	return Vectors::ZonedVector<double>(dim1, {4, 2}, {5, 3}, {6, 1}, 2);
}

// Sizes and offsets account for the ghost cells and leading dimension
TEST(ZonedVectorCPU, Layout)
{
	auto v = make_two_zone();

	EXPECT_EQ(v.get_nzones(), 2);
	EXPECT_EQ(v.get_dim1(), 3);
	EXPECT_EQ(v.get_ghost(), 2);
	EXPECT_EQ(v.get_nx(1), 2);
	EXPECT_EQ(v.get_ny(1), 3);
	EXPECT_EQ(v.get_nz(1), 1);

	const std::size_t zone0 {3u * 8 * 9 * 10};
	const std::size_t zone1 {3u * 6 * 7 * 5};
	EXPECT_EQ(v.get_zone_offset(0), 0u);
	EXPECT_EQ(v.get_zone_size(0), zone0);
	EXPECT_EQ(v.get_zone_offset(1), zone0);
	EXPECT_EQ(v.get_zone_size(1), zone1);
	EXPECT_EQ(v.get_data().size(), zone0 + zone1);
}

// The first and last ghost cells of every zone map to the start and end of
// that zone's block, and every element is reached exactly once
TEST(ZonedVectorCPU, GhostIndexingCoversEveryElement)
{
	auto v = make_two_zone();

	EXPECT_EQ(v.calc_index(0, 0, -2, -2, -2), 0u);
	EXPECT_EQ(v.calc_index(0, 2, 5, 6, 7), v.get_zone_size(0) - 1);
	EXPECT_EQ(v.calc_index(1, 0, -2, -2, -2), v.get_zone_offset(1));
	EXPECT_EQ(v.calc_index(1, 2, 3, 4, 2), v.get_data().size() - 1);

	std::vector<int> hits(v.get_data().size(), 0);
	for (int z {}; z < v.get_nzones(); ++z)
	for (int t {}; t < v.get_dim1(); ++t)
	for (int i {-2}; i < v.get_nx(z) + 2; ++i)
	for (int j {-2}; j < v.get_ny(z) + 2; ++j)
	for (int k {-2}; k < v.get_nz(z) + 2; ++k)
	{
		hits[v.calc_index(z, t, i, j, k)]++;
	}
	for (int h : hits) EXPECT_EQ(h, 1);
}

// Writing through one index reads back through the same index and nowhere
// else
TEST(ZonedVectorCPU, ReadWrite)
{
	auto v = make_two_zone();
	v(1, 2, -1, 3, 0) = 7.5;
	EXPECT_DOUBLE_EQ(v(1, 2, -1, 3, 0), 7.5);
	EXPECT_DOUBLE_EQ(v(0, 2, -1, 3, 0), 0.0);
	EXPECT_DOUBLE_EQ(v(1, 1, -1, 3, 0), 0.0);

	// The version without t uses t = 0
	Vectors::ZonedVector<double> g(1, {3}, {3}, {3}, 1);
	g(0, -1, 0, 3) = 2.0;
	EXPECT_DOUBLE_EQ(g(0, 0, -1, 0, 3), 2.0);
}

// A single zone without ghost cells has the same layout as a Vector4D, so it
// can be handed to code (e.g., the GPU copies) expecting one
TEST(ZonedVectorCPU, SingleZoneMatchesVector4D)
{
	Vectors::Vector4D<double> v4(2, 3, 4, 5);
	auto& d = v4.get_data();
	for (std::size_t n {}; n < d.size(); ++n) d[n] = static_cast<double>(n);

	Vectors::ZonedVector<double> zv(2, {3}, {4}, {5});
	for (int t {}; t < 2; ++t)
	for (int i {}; i < 3; ++i)
	for (int j {}; j < 4; ++j)
	for (int k {}; k < 5; ++k)
	{
		EXPECT_EQ(zv.calc_index(0, t, i, j, k),
			static_cast<std::size_t>(v4.calc_index(t, i, j, k)));
	}

	// Moving in a Vector4D keeps every value at the same index
	Vectors::Vector4D<double> ref {v4};
	Vectors::ZonedVector<double> moved {std::move(v4)};
	EXPECT_EQ(moved.get_nzones(), 1);
	EXPECT_EQ(moved.get_dim1(), 2);
	EXPECT_EQ(moved.get_nz(0), 5);
	EXPECT_DOUBLE_EQ(moved(0, 1, 2, 3, 4), ref(1, 2, 3, 4));
	EXPECT_EQ(moved.get_data(), ref.get_data());
}

// Moving in a Vector3D gives a single zone with dim1 = 1
TEST(ZonedVectorCPU, FromVector3D)
{
	Vectors::Vector3D<double> v3(2, 3, 4);
	v3(1, 2, 3) = 9.0;
	Vectors::ZonedVector<double> zv {std::move(v3)};
	EXPECT_EQ(zv.get_dim1(), 1);
	EXPECT_EQ(zv.get_nx(0), 2);
	EXPECT_DOUBLE_EQ(zv(0, 1, 2, 3), 9.0);
}

TEST(ZonedVectorCPU, SameShape)
{
	auto a = make_two_zone();
	auto b = make_two_zone();
	auto c = make_two_zone(2);
	Vectors::ZonedVector<double> d(3, {4, 2}, {5, 3}, {6, 1}, 1);
	EXPECT_TRUE(a.check_same_shape(b));
	EXPECT_FALSE(a.check_same_shape(c));
	EXPECT_FALSE(a.check_same_shape(d));
}
