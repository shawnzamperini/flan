/**
* @file zones_hd.h
*
* @brief The zones' grids as flat arrays (ZonesView), and the index
* functions that use them, shared by the CPU and the GPU.
*
* Zones (zones.h) keeps its grids in nested std::vectors, which the GPU
* can't use. A ZonesView holds the same grids as flat arrays: every zone's
* values one after the other, with per-zone offsets. Zones builds one on the
* host (Zones::view), and Zones::to_device copies its arrays to a GPU and
* returns the same struct with device pointers. The functions here take
* either, so the CPU and the GPU find cells, neighbors and ghost cell
* centers with the same code.
*
* Nothing here includes MPI or the standard containers, so it can be used
* in .cu files.
*/
#pragma once

#include <cmath>
#include <cstddef>

#include "host_device.h"

namespace Zones
{
	/**
	* Faces of a zone. The order matches the rows of the SOLEDGE3X neighbors
	* array (psi-, psi+, theta-, theta+, phi-, phi+), and the opposite face
	* of f is f ^ 1.
	*/
	enum Face : int
	{
		xmin,
		xmax,
		ymin,
		ymax,
		zmin,
		zmax,
		num_faces
	};

	/// Face code for a face where particles are absorbed
	inline constexpr int face_absorbing {-1};

	/// Face code for a face where particles are put back in at the core
	inline constexpr int face_core {-2};

	/**
	* Face code for a material surface, such as a divertor target. Particles
	* are absorbed there the same as at an absorbing face, but are counted
	* separately. -4 rather than -3, since -3 is SOLEDGE3X's code for its
	* outer psi boundary, which is an absorbing face in Flan.
	*/
	inline constexpr int face_wall {-4};

	/**
	* Face code for a z face of a single-zone (Gkeyll) background that is a
	* last closed flux surface boundary (bound type "lcfs"): periodic for
	* particles inside lcfs_x (closed field lines) and absorbing outside it
	* (see Boundary::lcfs_bc). Only z faces can have it.
	*/
	inline constexpr int face_lcfs {-5};

	/// Bits of ZonesView::flags
	inline constexpr int flag_toroidal {1};		///< See Zones::set_toroidal
	inline constexpr int flag_cartesian {2};	///< See Zones::set_cartesian
	inline constexpr int flag_wedge {4};		///< See Zones::is_wedge

	/**
	* The zones' grids and connections as flat arrays, on the host or on a
	* GPU.
	*
	* A plain struct of pointers and sizes, so it can be passed to kernels
	* by value. The pointers are to arrays owned by someone else: Zones on
	* the host (Zones::view), or the GPU copy made by Zones::to_device and
	* freed with free_zones.
	*
	* Cell centers are stored with the ghost cells: zone z's centers in
	* dimension d are centers[d][center_offset[d][z] + i + ghost] for
	* -ghost <= i < n + ghost. A ghost cell's center is the center of the
	* cell it copies, in this zone's frame (Zones::center). Ghost cells past
	* a boundary face have no center and hold NaN.
	*/
	struct ZonesView
	{
		int nzones;		///< Number of zones
		int ghost;		///< Ghost cells on each side of every zone

		/// Interior cells per zone in x, y and z (nzones values each)
		const int* n[3];

		/// Face codes, faces[zone * num_faces + face] (see Zones::get_face)
		const int* faces;

		/// Start of each zone's centers in centers[d] (nzones values)
		const int* center_offset[3];

		/// Cell centers including the ghost cells, n + 2 ghost per zone
		const double* centers[3];

		/// Start of each zone's cell edges in edges[d] (nzones values)
		const int* edge_offset[3];

		/// Cell edges, n + 1 per zone
		const double* edges[3];

		/**
		* Start of each zone in a ZonedVector with dim1 = 1 and this ghost
		* width (nzones + 1 values; the last is the total size). With dim1
		* frames, zone z starts at dim1 * cell_offset[z].
		*/
		const std::size_t* cell_offset;

		/// flag_* bits per zone (nzones values)
		const int* flags;

		/**
		* Start of each zone's zShift table in z_shift, or -1 for a zone
		* without one (nzones values). The table covers the zone's x-y
		* cells including ghosts (see Zones::set_z_shift).
		*/
		const long long* z_shift_offset;

		/// zShift tables of the zones that have one, one after the other
		const double* z_shift;
	};

	/// Indices of the two cells whose centers bracket a coordinate
	struct IndexPair
	{
		int lo;		///< Lower cell, which may be a ghost cell
		int hi;		///< Upper cell, which may be a ghost cell
	};

	/**
	* @param v The zones.
	* @param zone Zone index.
	* @param face Face of the zone.
	* @return The face's code (see Zones::get_face).
	*/
	FLAN_HD inline int face(const ZonesView& v, const int zone,
		const int face)
	{
		return v.faces[zone * num_faces + face];
	}

	/**
	* Coordinate of a cell center in a zone's frame, which may be a ghost
	* cell (see Zones::center).
	*
	* @param v The zones.
	* @param zone Zone index.
	* @param dim Dimension: 0 = x, 1 = y, 2 = z.
	* @param i Cell index, -ghost <= i < n + ghost.
	* @return The center, or NaN for a ghost cell past a boundary face.
	*/
	FLAN_HD inline double center(const ZonesView& v, const int zone,
		const int dim, const int i)
	{
		return v.centers[dim][v.center_offset[dim][zone] + i + v.ghost];
	}

	/**
	* Pointer to a zone's interior cell centers in one dimension (n values).
	*
	* @param v The zones.
	* @param zone Zone index.
	* @param dim Dimension: 0 = x, 1 = y, 2 = z.
	* @return Pointer to the center of cell 0.
	*/
	FLAN_HD inline const double* interior_centers(const ZonesView& v,
		const int zone, const int dim)
	{
		return v.centers[dim] + v.center_offset[dim][zone] + v.ghost;
	}

	/**
	* Pointer to a zone's cell edges in one dimension (n + 1 values).
	*
	* @param v The zones.
	* @param zone Zone index.
	* @param dim Dimension: 0 = x, 1 = y, 2 = z.
	* @return Pointer to the first edge.
	*/
	FLAN_HD inline const double* zone_edges(const ZonesView& v,
		const int zone, const int dim)
	{
		return v.edges[dim] + v.edge_offset[dim][zone];
	}

	/**
	* Index of the cell of a zone containing a coordinate.
	*
	* Same as Indices::get_nearest_cell_index_cpu on the zone's edges: a
	* coordinate on the edge between two cells is in the lower cell, and one
	* outside the zone is in the first or last cell.
	*
	* @param v The zones.
	* @param zone Zone index.
	* @param dim Dimension: 0 = x, 1 = y, 2 = z.
	* @param coord The coordinate.
	* @return Cell index, 0 <= i < n.
	*/
	FLAN_HD inline int cell_index(const ZonesView& v, const int zone,
		const int dim, const double coord)
	{
		const double* e {zone_edges(v, zone, dim)};
		const int n_edges {v.n[dim][zone] + 1};

		// Index of the first edge not less than coord (std::lower_bound,
		// which device code can't call)
		int left {0};
		int right {n_edges};
		while (left < right)
		{
			const int mid {(left + right) >> 1};
			if (e[mid] < coord) left = mid + 1;
			else right = mid;
		}
		if (left == 0) return 0;
		if (left == n_edges) return n_edges - 2;
		return left - 1;
	}

	/**
	* The cell to interpolate towards from a particle's cell, in one
	* dimension (see Zones::neighbor_index).
	*
	* The neighbor on the side of the cell center the particle is on, or
	* the only neighbor at the edge of a zone, as
	* Utilities::get_neighbor_index. With ghost cells, at the edge of a zone
	* whose face is shared with another zone (or periodic), the ghost cell
	* across the face is returned instead.
	*
	* @param v The zones.
	* @param zone Zone the particle is in.
	* @param dim Dimension: 0 = x, 1 = y, 2 = z.
	* @param coord Particle's coordinate in that dimension.
	* @param idx Index of the cell the particle is in.
	* @return Index of the neighboring cell, which may be a ghost cell.
	*/
	FLAN_HD inline int neighbor_index(const ZonesView& v, const int zone,
		const int dim, const double coord, const int idx)
	{
		const int n {v.n[dim][zone]};
		const double dx_from_center {coord - center(v, zone, dim, idx)};
		if (v.ghost > 0)
		{
			if (idx == 0 && dx_from_center <= 0.0
				&& face(v, zone, 2 * dim) >= 0)
			{
				return -1;
			}
			if (idx == n - 1 && dx_from_center > 0.0
				&& face(v, zone, 2 * dim + 1) >= 0)
			{
				return n;
			}
		}

		// Utilities::get_neighbor_index: the side of the center, except at
		// the first or last cell, which have one neighbor
		const int side {2 * (dx_from_center > 0.0) - 1};
		const int is_left_edge {idx == 0};
		const int is_right_edge {idx == n - 1};
		const int offset {side * (1 - is_left_edge - is_right_edge)
			+ is_left_edge - is_right_edge};
		return idx + offset;
	}

	/**
	* The two cells whose centers bracket a coordinate, in one dimension
	* (see Zones::bracket).
	*
	* Utilities::bracket_indices on the zone's centers (the first or last
	* two cells past the outer centers), except that with ghost cells, past
	* the first or last center of a face shared with another zone (or
	* periodic), the pair is the ghost cell and the edge cell.
	*
	* @param v The zones.
	* @param zone Zone the point is in.
	* @param dim Dimension: 0 = x, 1 = y, 2 = z.
	* @param coord Coordinate of the point in that dimension.
	* @return The lower and upper cells.
	*/
	FLAN_HD inline IndexPair bracket(const ZonesView& v, const int zone,
		const int dim, const double coord)
	{
		const int n {v.n[dim][zone]};
		const double* c {interior_centers(v, zone, dim)};
		if (v.ghost > 0)
		{
			if (coord < c[0] && face(v, zone, 2 * dim) >= 0)
				return {-1, 0};
			if (coord > c[n - 1] && face(v, zone, 2 * dim + 1) >= 0)
				return {n - 1, n};
		}

		// Utilities::bracket_indices: i is the first center not less than
		// coord (std::lower_bound)
		int left {0};
		int right {n};
		while (left < right)
		{
			const int mid {(left + right) >> 1};
			if (c[mid] < coord) left = mid + 1;
			else right = mid;
		}
		if (left == 0) return {0, 1};
		const int hi {left < n - 1 ? left : n - 1};
		const int lo {hi - 1 > 0 ? hi - 1 : 0};
		return {lo, hi};
	}

	/**
	* @param v The zones.
	* @param zone Zone to check.
	* @return True if the zone has a zShift table (see Zones::set_z_shift).
	*/
	FLAN_HD inline bool has_z_shift(const ZonesView& v, const int zone)
	{
		return v.z_shift_offset[zone] >= 0;
	}

	/**
	* A zone's zShift at a cell (see Zones::z_shift). The zone must have a
	* table.
	*
	* @param v The zones.
	* @param zone Zone of the cell.
	* @param i x index of the cell, which may be a ghost index.
	* @param j y index of the cell, which may be a ghost index.
	* @return zShift (radians of z).
	*/
	FLAN_HD inline double z_shift(const ZonesView& v, const int zone,
		const int i, const int j)
	{
		const long long py {v.n[1][zone] + 2 * v.ghost};
		return v.z_shift[v.z_shift_offset[zone]
			+ static_cast<long long>(i + v.ghost) * py + j + v.ghost];
	}

	/**
	* Wraps a z coordinate into a zone's z range, for a zone that is
	* periodic in z (see Zones::wrap_z).
	*
	* @param v The zones.
	* @param zone Zone whose range to wrap into.
	* @param z Coordinate to wrap, any number of periods away.
	* @return z plus a whole number of periods, in [z_min, z_min + period).
	*/
	FLAN_HD inline double wrap_z(const ZonesView& v, const int zone,
		const double z)
	{
		const double* e {zone_edges(v, zone, 2)};
		const double lo {e[0]};
		const double period {e[v.n[2][zone]] - lo};
		double w {std::fmod(z - lo, period)};
		if (w < 0.0) w += period;
		return lo + w;
	}

	/**
	* Index into the flat data of a ZonedVector laid out for these zones
	* (same as ZonedVector::calc_index).
	*
	* @param v The zones.
	* @param dim1 Size of the ZonedVector's leading dimension (1 for the
	*        geometry, the number of frames for the plasma).
	* @param zone Zone index.
	* @param t Index in the leading dimension.
	* @param i x index, -ghost <= i < nx + ghost.
	* @param j y index, likewise.
	* @param k z index, likewise.
	* @return Index into the flat data.
	*/
	FLAN_HD inline std::size_t flat_index(const ZonesView& v, const int dim1,
		const int zone, const int t, const int i, const int j, const int k)
	{
		const int g {v.ghost};
		const std::size_t pz {static_cast<std::size_t>(v.n[2][zone] + 2 * g)};
		const std::size_t py {static_cast<std::size_t>(v.n[1][zone] + 2 * g)};
		const std::size_t px {static_cast<std::size_t>(v.n[0][zone] + 2 * g)};
		return static_cast<std::size_t>(dim1) * v.cell_offset[zone]
			+ static_cast<std::size_t>(t) * (px * py * pz)
			+ static_cast<std::size_t>(i + g) * (py * pz)
			+ static_cast<std::size_t>(j + g) * pz
			+ static_cast<std::size_t>(k + g);
	}
}
