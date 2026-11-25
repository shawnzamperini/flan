/**
* @file test_zones_cuda.cu
*
* @brief Tests of the zone functions in zones_hd.h on the GPU, against the
* same functions on the host and against Zones and ZonedVector.
*
* The zones are copied with Zones::to_device, so these also check that the
* copy holds the same grids, faces and offsets as the host's view.
*/
#include <cmath>
#include <utility>
#include <vector>

#include <cuda_runtime.h>
#include <gtest/gtest.h>

#include "pcg32.h"
#include "vectors.h"
#include "zones.h"
#include "zones_hd.h"


namespace
{
	constexpr double two_pi {6.283185307179586};

	/// One query: a zone, a dimension and a coordinate in it
	struct Query
	{
		int zone;
		int dim;
		double coord;
	};

	/// What the functions in zones_hd.h give for a Query
	struct Answer
	{
		int cell;			///< cell_index
		int neighbor;		///< neighbor_index from that cell
		int lo;				///< bracket
		int hi;				///< bracket
		double c_cell;		///< center of cell
		double c_neighbor;	///< center of neighbor
		double c_lo;		///< center of lo
		double c_hi;		///< center of hi
	};

	/**
	* Answers a query with the functions in zones_hd.h. The same code runs
	* on the host and in query_kernel.
	*/
	__host__ __device__ Answer answer(const Zones::ZonesView& v,
		const Query& q)
	{
		Answer a {};
		a.cell = Zones::cell_index(v, q.zone, q.dim, q.coord);
		a.neighbor = Zones::neighbor_index(v, q.zone, q.dim, q.coord, a.cell);
		const Zones::IndexPair p {Zones::bracket(v, q.zone, q.dim, q.coord)};
		a.lo = p.lo;
		a.hi = p.hi;
		a.c_cell = Zones::center(v, q.zone, q.dim, a.cell);
		a.c_neighbor = Zones::center(v, q.zone, q.dim, a.neighbor);
		a.c_lo = Zones::center(v, q.zone, q.dim, a.lo);
		a.c_hi = Zones::center(v, q.zone, q.dim, a.hi);
		return a;
	}

	__global__ void query_kernel(const Zones::ZonesView v, const Query* q,
		Answer* out, int n)
	{
		const int i {static_cast<int>(blockIdx.x * blockDim.x
			+ threadIdx.x)};
		if (i < n) out[i] = answer(v, q[i]);
	}

	/// One ZonedVector index: dim1, zone, t, i, j, k
	struct FlatQuery
	{
		int dim1, zone, t, i, j, k;
	};

	__global__ void flat_kernel(const Zones::ZonesView v,
		const FlatQuery* q, std::size_t* out, int n)
	{
		const int m {static_cast<int>(blockIdx.x * blockDim.x
			+ threadIdx.x)};
		if (m < n)
		{
			out[m] = Zones::flat_index(v, q[m].dim1, q[m].zone, q[m].t,
				q[m].i, q[m].j, q[m].k);
		}
	}

	/// Copies a vector to the device, runs a kernel over it, copies back
	template <typename In, typename Out, typename Launch>
	std::vector<Out> run(const std::vector<In>& in, Launch launch)
	{
		const int n {static_cast<int>(in.size())};
		In* in_d {};
		Out* out_d {};
		cudaMalloc(&in_d, n * sizeof(In));
		cudaMalloc(&out_d, n * sizeof(Out));
		cudaMemcpy(in_d, in.data(), n * sizeof(In), cudaMemcpyHostToDevice);
		launch(in_d, out_d, n);
		std::vector<Out> out(n);
		cudaMemcpy(out.data(), out_d, n * sizeof(Out),
			cudaMemcpyDeviceToHost);
		EXPECT_EQ(cudaGetLastError(), cudaSuccess);
		cudaFree(in_d);
		cudaFree(out_d);
		return out;
	}

	/// Grid from its edges, as (centers, edges)
	std::pair<std::vector<double>, std::vector<double>> from_edges(
		std::vector<double> edges)
	{
		std::vector<double> centers(edges.size() - 1);
		for (std::size_t i {}; i < centers.size(); ++i)
			centers[i] = 0.5 * (edges[i] + edges[i + 1]);
		return {centers, edges};
	}

	/// Grid of n uneven cells from lo to hi (widths 1, 2, 1, 2, ...)
	std::pair<std::vector<double>, std::vector<double>> uneven(int n,
		double lo, double hi)
	{
		std::vector<double> w(n);
		double total {};
		for (int i {}; i < n; ++i) total += (w[i] = 1.0 + (i % 2));
		std::vector<double> edges {lo};
		for (int i {}; i < n; ++i)
			edges.push_back(edges.back() + (hi - lo) * w[i] / total);
		edges.back() = hi;
		return from_edges(edges);
	}

	/**
	* Three zones in a row in x with ghost width 2, the middle one a single
	* cell thinner than the ghost width, uneven cells, periodic in z, and
	* zone 2 periodic in y. The other faces are absorbing, core or wall.
	*/
	Zones::Zones make_zones()
	{
		Zones::Zones zones(2);
		const std::array<int, 3> nx {4, 1, 3};
		const std::array<int, 3> ny {3, 3, 3};
		const std::array<int, 3> nz {5, 5, 5};
		for (int z {}; z < 3; ++z)
		{
			auto [x, gx] = uneven(nx[z], 0.0, nx[z]);
			auto [y, gy] = uneven(ny[z], -1.0, 2.0);
			auto [zc, gz] = uneven(nz[z], 0.0, two_pi);
			zones.add_zone(x, y, zc, gx, gy, gz);
			zones.set_face(z, Zones::zmin, z);
			zones.set_face(z, Zones::zmax, z);
		}
		zones.set_face(0, Zones::xmax, 1);
		zones.set_face(1, Zones::xmin, 0);
		zones.set_face(1, Zones::xmax, 2);
		zones.set_face(2, Zones::xmin, 1);
		zones.set_face(0, Zones::xmin, Zones::face_core);
		zones.set_face(2, Zones::xmax, Zones::face_wall);
		zones.set_face(2, Zones::ymin, 2);
		zones.set_face(2, Zones::ymax, 2);
		zones.validate();
		return zones;
	}

	/// Equal, or both NaN (a ghost cell past a boundary face)
	void expect_same(double a, double b, const char* what, int i)
	{
		if (std::isnan(a) || std::isnan(b))
		{
			EXPECT_TRUE(std::isnan(a) && std::isnan(b)) << what << " " << i;
		}
		else
		{
			EXPECT_EQ(a, b) << what << " " << i;
		}
	}
}


// Cells, neighbors, brackets and centers on the device match the host,
// for points inside each zone and up to a cell past its faces
TEST(ZonesCUDA, MatchesHost)
{
	const Zones::Zones zones {make_zones()};
	Zones::ZonesView zones_d {zones.to_device(0)};

	pcg32 rng(2024, 3);
	std::vector<Query> queries {};
	for (int n {}; n < 20000; ++n)
	{
		const int zone {static_cast<int>(rng.next_double() * 3)};
		const int dim {static_cast<int>(rng.next_double() * 3)};
		const double lo {zones.get_min(zone, dim)};
		const double hi {zones.get_max(zone, dim)};
		const double pad {0.2 * (hi - lo)};
		queries.push_back({zone, dim,
			lo - pad + rng.next_double() * (hi - lo + 2 * pad)});
	}

	// The cell centers themselves, where neighbor_index picks a side
	for (int zone {}; zone < 3; ++zone)
	for (int dim {}; dim < 3; ++dim)
	for (const double c : zones.get_centers(zone, dim))
		queries.push_back({zone, dim, c});

	const auto got {run<Query, Answer>(queries,
		[&](const Query* q, Answer* out, int n)
		{query_kernel<<<(n + 127) / 128, 128>>>(zones_d, q, out, n);})};

	for (std::size_t i {}; i < queries.size(); ++i)
	{
		const Answer want {answer(zones.view(), queries[i])};
		const int n {static_cast<int>(i)};
		EXPECT_EQ(got[i].cell, want.cell) << n;
		EXPECT_EQ(got[i].neighbor, want.neighbor) << n;
		EXPECT_EQ(got[i].lo, want.lo) << n;
		EXPECT_EQ(got[i].hi, want.hi) << n;
		expect_same(got[i].c_cell, want.c_cell, "c_cell", n);
		expect_same(got[i].c_neighbor, want.c_neighbor, "c_neighbor", n);
		expect_same(got[i].c_lo, want.c_lo, "c_lo", n);
		expect_same(got[i].c_hi, want.c_hi, "c_hi", n);

		// And the host's view matches Zones' own functions
		const Query& q {queries[i]};
		EXPECT_EQ(want.neighbor,
			zones.neighbor_index(q.zone, q.dim, q.coord, want.cell)) << n;
		EXPECT_EQ(std::make_pair(want.lo, want.hi),
			zones.bracket(q.zone, q.dim, q.coord)) << n;
	}

	Zones::free_zones(zones_d, 0);
	EXPECT_EQ(zones_d.centers[0], nullptr);
}


// Every ghost cell center in the table is the one Zones::center finds by
// stepping across the faces (NaN past a boundary face)
TEST(ZonesCUDA, GhostCentersMatchZones)
{
	const Zones::Zones zones {make_zones()};
	const Zones::ZonesView& v {zones.view()};
	const int g {zones.get_ghost()};
	for (int zone {}; zone < 3; ++zone)
	for (int dim {}; dim < 3; ++dim)
	for (int i {-g}; i < zones.get_n(zone, dim) + g; ++i)
	{
		const double c {Zones::center(v, zone, dim, i)};
		const int face {i < 0 ? 2 * dim : 2 * dim + 1};
		const bool interior {i >= 0 && i < zones.get_n(zone, dim)};
		if (!interior && zones.get_face(zone,
			static_cast<Zones::Face>(face)) < 0)
		{
			EXPECT_TRUE(std::isnan(c)) << zone << " " << dim << " " << i;
		}
		else
		{
			EXPECT_EQ(c, zones.center(zone, dim, i))
				<< zone << " " << dim << " " << i;
		}
	}

	// The thin middle zone: zone 0's second ghost past xmax is in zone 2
	// (zone 0's max is 4, zone 1 is one cell wide, zone 2's first cell
	// is 1 of 4 units of its 3-cell width)
	EXPECT_DOUBLE_EQ(Zones::center(v, 0, 0, 5), 4.0 + 1.0 + 3.0 / 8.0);
}


// flat_index on the device is ZonedVector::calc_index
TEST(ZonesCUDA, FlatIndexMatchesZonedVector)
{
	const Zones::Zones zones {make_zones()};
	Zones::ZonesView zones_d {zones.to_device(0)};
	const int g {zones.get_ghost()};

	pcg32 rng(77, 1);
	auto pick = [&rng](int lo, int hi)
		{return lo + static_cast<int>(rng.next_double() * (hi - lo));};

	std::vector<FlatQuery> queries {};
	std::vector<std::size_t> want {};
	for (const int dim1 : {1, 3})
	{
		const Vectors::ZonedVector<double> v(dim1, zones.get_nx_all(),
			zones.get_ny_all(), zones.get_nz_all(), g);
		for (int n {}; n < 2000; ++n)
		{
			FlatQuery q {};
			q.dim1 = dim1;
			q.zone = pick(0, 3);
			q.t = pick(0, dim1);
			q.i = pick(-g, zones.get_nx(q.zone) + g);
			q.j = pick(-g, zones.get_ny(q.zone) + g);
			q.k = pick(-g, zones.get_nz(q.zone) + g);
			queries.push_back(q);
			want.push_back(v.calc_index(q.zone, q.t, q.i, q.j, q.k));
		}
	}

	const auto got {run<FlatQuery, std::size_t>(queries,
		[&](const FlatQuery* q, std::size_t* out, int n)
		{flat_kernel<<<(n + 127) / 128, 128>>>(zones_d, q, out, n);})};
	for (std::size_t i {}; i < queries.size(); ++i)
	{
		EXPECT_EQ(got[i], want[i]) << i;
		EXPECT_EQ(Zones::flat_index(zones.view(), queries[i].dim1,
			queries[i].zone, queries[i].t, queries[i].i, queries[i].j,
			queries[i].k), want[i]) << i;
	}

	Zones::free_zones(zones_d, 0);
}


// A copy of a Zones has its own view, pointing at its own arrays
TEST(ZonesCUDA, CopyHasItsOwnView)
{
	const Zones::Zones a {make_zones()};
	const Zones::Zones b {a};
	EXPECT_NE(a.view().centers[0], b.view().centers[0]);
	EXPECT_EQ(Zones::center(b.view(), 0, 0, 5), Zones::center(a.view(), 0,
		0, 5));
}
