#pragma once
/**
* @file interp_hd.h
*
* @brief Field-aligned interpolation on flat zone data, shared by the CPU
*        and the GPU (see host_device.h).
*
* Utilities::interp_aligned (utilities.cpp) on the CPU and the GPU kernels
* (interp_plasma in cuda/interpolate.cuh) both call interp_aligned_flat, so
* a shifted-metric (HERMES) zone is interpolated the same way on both. The
* fields are the flat data of ZonedVectors laid out for the zones (see
* Zones::flat_index), with one entry of the leading dimension per
* background time.
*/

#include <cmath>
#include <cstddef>

#include "host_device.h"
#include "zones_hd.h"


namespace Utilities
{
	/**
	* The two values of a sorted array that bracket x0 (same as
	* Utilities::bracket_indices): the first value not less than x0 and the
	* one before it, or the first or last two past the ends. An array of one
	* value gives {0, 0}.
	*
	* @param x Sorted values.
	* @param n Number of values.
	* @param x0 Value to bracket.
	* @return The lower and upper indices.
	*/
	FLAN_HD inline Zones::IndexPair bracket_sorted(const double* x,
		const int n, const double x0)
	{
		// std::lower_bound, which device code can't call
		int left {0};
		int right {n};
		while (left < right)
		{
			const int mid {(left + right) >> 1};
			if (x[mid] < x0) left = mid + 1;
			else right = mid;
		}
		if (left == 0) return {0, n > 1 ? 1 : 0};
		const int hi {left < n - 1 ? left : n - 1};
		const int lo {hi - 1 > 0 ? hi - 1 : 0};
		return {lo, hi};
	}

	/**
	* Weight of the upper of two bracketing points c0 and c1 for a point c,
	* clamped to [0, 1] so a point past the last centers takes the edge
	* value. Two equal points give 0.
	*/
	FLAN_HD inline double bracket_weight(const double c, const double c0,
		const double c1)
	{
		return c1 == c0 ? 0.0
			: std::fmin(1.0, std::fmax(0.0, (c - c0) / (c1 - c0)));
	}

	/**
	* Field-aligned interpolation of N fields at a point in a zone with a
	* zShift table. See Utilities::interp_aligned for the method: x is
	* interpolated at fixed z, and in each of the two bracketing x columns
	* the field line through the point is followed to the two bracketing y
	* planes, where the field is interpolated linearly in z; the samples are
	* combined bilinearly in (x, y) and linearly in time.
	*
	* With phi_sign nonzero, the fields are the X, Y and Z components of a
	* Cartesian vector (N = 3), and each sample's (X, Y) part is rotated
	* about Z by phi_sign (z0 - zf) to the point's toroidal angle, zf being
	* the sample cell's own angle.
	*
	* @tparam T Data type of the fields.
	* @tparam N Number of fields (1, or 3 for a vector).
	* @param fields Flat data of each field (see Zones::flat_index), with
	*        ntimes entries in the leading dimension.
	* @param zones The zones.
	* @param times Background times, sorted.
	* @param ntimes Number of times.
	* @param zone Zone to sample in. It must have a zShift table.
	* @param t0 Time at which to evaluate the fields (clamped to the times).
	* @param x0 x coordinate of the point.
	* @param y0 y coordinate of the point.
	* @param z0 z coordinate (toroidal angle) of the point.
	* @param phi_sign 0 for scalars; for a vector, +1 if (R, phi, Z) is
	*        right handed and -1 if not.
	* @param[out] out Interpolated value of each field (N values).
	*/
	template <typename T, int N>
	FLAN_HD inline void interp_aligned_flat(const T* const* fields,
		const Zones::ZonesView& zones, const double* times, const int ntimes,
		const int zone, const double t0, const double x0, const double y0,
		const double z0, const double phi_sign, double* out)
	{
		static_assert(N == 1 || N == 3);
		const bool rotate {phi_sign != 0.0};

		const double tc {std::fmin(times[ntimes - 1],
			std::fmax(times[0], t0))};
		const Zones::IndexPair it {bracket_sorted(times, ntimes, tc)};
		const double wt {bracket_weight(tc, times[it.lo], times[it.hi])};
		const Zones::IndexPair ixs {Zones::bracket(zones, zone, 0, x0)};
		const double wx {bracket_weight(x0,
			Zones::center(zones, zone, 0, ixs.lo),
			Zones::center(zones, zone, 0, ixs.hi))};
		const Zones::IndexPair iys {Zones::bracket(zones, zone, 1, y0)};
		const double wy {bracket_weight(y0,
			Zones::center(zones, zone, 1, iys.lo),
			Zones::center(zones, zone, 1, iys.hi))};

		for (int n {}; n < N; ++n) out[n] = 0.0;
		for (int a {}; a < 2; ++a)
		{
			const int ix {a ? ixs.hi : ixs.lo};
			const double w_x {a ? wx : 1.0 - wx};
			if (w_x == 0.0) continue;

			// zShift of the field line through the point, in this column
			const double zs_p {(1.0 - wy) * Zones::z_shift(zones, zone, ix,
				iys.lo) + wy * Zones::z_shift(zones, zone, ix, iys.hi)};
			for (int b {}; b < 2; ++b)
			{
				const int iy {b ? iys.hi : iys.lo};
				const double w_xy {w_x * (b ? wy : 1.0 - wy)};
				if (w_xy == 0.0) continue;

				// Where the field line crosses this plane, and the two
				// cells around it in z (one may be a periodic ghost)
				const double zc {Zones::wrap_z(zones, zone,
					z0 + Zones::z_shift(zones, zone, ix, iy) - zs_p)};
				const Zones::IndexPair izs {Zones::bracket(zones, zone, 2,
					zc)};
				const double wz {bracket_weight(zc,
					Zones::center(zones, zone, 2, izs.lo),
					Zones::center(zones, zone, 2, izs.hi))};

				// Each of the two z cells' values, blended in time, and
				// for a vector rotated from its cell's frame to z0
				double sample[N] {};
				for (int c {}; c < 2; ++c)
				{
					const int iz {c ? izs.hi : izs.lo};
					const double w_z {c ? wz : 1.0 - wz};
					const std::size_t i0 {Zones::flat_index(zones, ntimes,
						zone, it.lo, ix, iy, iz)};
					const std::size_t i1 {Zones::flat_index(zones, ntimes,
						zone, it.hi, ix, iy, iz)};
					double v[N];
					for (int n {}; n < N; ++n)
					{
						v[n] = (1.0 - wt) * static_cast<double>(fields[n][i0])
							+ wt * static_cast<double>(fields[n][i1]);
					}
					if constexpr (N == 3)
					{
						if (rotate)
						{
							const double zf {Zones::center(zones, zone, 2,
								iz)};
							const double ang {phi_sign * (z0 - zf)};
							const double cs {std::cos(ang)};
							const double sn {std::sin(ang)};
							const double vx {v[0]};
							v[0] = cs * vx - sn * v[1];
							v[1] = sn * vx + cs * v[1];
						}
					}
					for (int n {}; n < N; ++n) sample[n] += w_z * v[n];
				}
				for (int n {}; n < N; ++n) out[n] += w_xy * sample[n];
			}
		}
	}

}  // namespace Utilities
