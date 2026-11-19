/**
* @file axisym_mesh.cpp
*
* @brief Geometry, magnetic field, electric field and flows of a multi-zone
* background on an axisymmetric mesh.
*
* Codes like SOLEDGE3X and HERMES describe their grid by the (R, Z) corners
* of the cells of a poloidal plane, split into logically rectangular zones,
* with the field (Br, Bz, Bphi) at the cell centers and the same mesh in
* every toroidal plane. Each zone has index-space coordinates x and y (the
* cell indices across and along the flux surfaces) and z, the toroidal
* angle. Everything here works from an Axisym::Mesh in SI units, so a
* reader only has to fill one in (and the background's zones) to get the
* geometry, B, grad B, E and the flows the impurities need.
*/

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <iostream>
#include <numbers>
#include <string>
#include <utility>
#include <tuple>
#include <vector>

#include "axisym_mesh.h"
#include "background.h"
#include "braginskii.h"
#include "errors.h"
#include "flan_types.h"
#include "vectors.h"
#include "zones.h"

namespace Axisym
{
	/**
	* An empty field shaped like the zones, with their ghost cells.
	*
	* @param zones The zones.
	* @param dim1 Size of the leading dimension (e.g., time frames).
	* @return The field, filled with zeros.
	*/
	Vectors::ZonedVector<BkgFPType> zoned(const Zones::Zones& zones,
		const int dim1)
	{
		return Vectors::ZonedVector<BkgFPType>(dim1, zones.get_nx_all(),
			zones.get_ny_all(), zones.get_nz_all(), zones.get_ghost());
	}

	/**
	* A uniform grid of n cells between lo and hi.
	*
	* @param n Number of cells.
	* @param lo, hi Ends of the grid.
	* @return (centers, edges), with n and n + 1 values.
	*/
	std::pair<std::vector<double>, std::vector<double>> uniform_grid(
		const int n, const double lo, const double hi)
	{
		std::vector<double> centers(n), edges(n + 1);
		for (int i {}; i <= n; ++i) edges[i] = lo + i * (hi - lo) / n;
		for (int i {}; i < n; ++i)
			centers[i] = 0.5 * (edges[i] + edges[i + 1]);
		return {centers, edges};
	}

	/**
	* The cells a derivative at index j of a zone uses in y: j - 1 and
	* j + 1 where they exist (including ghost cells across a face shared
	* with another zone), otherwise j itself on that side.
	*
	* @param zones The zones.
	* @param zone Zone of the cell.
	* @param j y index of the cell.
	* @return {lo, hi}, equal if the zone is one cell wide with no
	*         neighbors.
	*/
	std::pair<int, int> y_stencil(const Zones::Zones& zones, const int zone,
		const int j)
	{
		const bool ghosts {zones.get_ghost() > 0};
		const int ny {zones.get_ny(zone)};
		const bool has_lo {j > 0 || (ghosts
			&& zones.get_face(zone, Zones::ymin) >= 0)};
		const bool has_hi {j < ny - 1 || (ghosts
			&& zones.get_face(zone, Zones::ymax) >= 0)};
		return {has_lo ? j - 1 : j, has_hi ? j + 1 : j};
	}

	/**
	* Derivatives of a field with respect to the computational coordinates
	* x, y and z at an interior cell.
	*
	* A central difference is used where the cells on both sides exist,
	* including ghost cells across faces shared with another zone (or
	* periodic), so ghost cells must be filled. At a face with nothing
	* across it (absorbing, core or wall), the difference is one-sided. Cell
	* spacings come from Zones::Zones::center, so they're right across
	* zones too.
	*
	* @param zones The zones.
	* @param f Field, with filled ghost cells.
	* @param zone Zone of the cell.
	* @param t Index along f's leading dimension.
	* @param i Cell index in x.
	* @param j Cell index in y.
	* @param k Cell index in z.
	* @param dfdy Optional df/dy along the field line from aligned_dfdy.
	*        In zones with a zShift table, df/dy at fixed z is then
	*        dfdy - (dzShift/dy) df/dz, with dzShift/dy over the same two
	*        planes and df/dz the central difference here, instead of a
	*        central difference in y at fixed z. So the gradient is
	*        df/dx grad x + dfdy grad y + df/dz (grad z - dzShift/dy
	*        grad y): the gradient in field-aligned coordinates, written
	*        in the (x, y, z) basis the particle push uses.
	* @return {df/dx, df/dy, df/dz}.
	*/
	std::array<double, 3> derivatives(const Zones::Zones& zones,
		const Vectors::ZonedVector<BkgFPType>& f, const int zone,
		const int t, const int i, const int j, const int k,
		const Vectors::ZonedVector<BkgFPType>* dfdy)
	{
		const std::array<int, 3> idx {i, j, k};
		const bool ghosts {zones.get_ghost() > 0};
		std::array<double, 3> out {};
		for (int d {}; d < 3; ++d)
		{
			const int n {zones.get_n(zone, d)};
			const bool has_lo {idx[d] > 0 || (ghosts
				&& zones.get_face(zone, static_cast<Zones::Face>(2 * d)) >= 0)};
			const bool has_hi {idx[d] < n - 1 || (ghosts
				&& zones.get_face(zone, static_cast<Zones::Face>(2 * d + 1))
				>= 0)};
			const int lo {has_lo ? idx[d] - 1 : idx[d]};
			const int hi {has_hi ? idx[d] + 1 : idx[d]};
			if (lo == hi) continue;

			std::array<int, 3> a {idx};
			std::array<int, 3> b {idx};
			a[d] = lo;
			b[d] = hi;
			out[d] = (f(zone, t, b[0], b[1], b[2]) - f(zone, t, a[0], a[1],
				a[2])) / (zones.center(zone, d, hi) - zones.center(zone, d, lo));
		}

		// Field-aligned df/dy, converted to fixed z
		if (dfdy && zones.has_z_shift(zone))
		{
			const auto [lo, hi] = y_stencil(zones, zone, j);
			if (lo != hi)
			{
				const double dzs_dy {(zones.z_shift(zone, i, hi)
					- zones.z_shift(zone, i, lo)) / (zones.center(zone, 1, hi)
					- zones.center(zone, 1, lo))};
				out[1] = (*dfdy)(zone, t, i, j, k) - dzs_dy * out[2];
			}
		}
		return out;
	}

	/**
	* Cartesian gradient of a field at an interior cell, from its
	* derivatives (see derivatives) and the background's reciprocal basis:
	* grad f = df/dx grad x + df/dy grad y + df/dz grad z.
	*
	* @param bkg Background holding the zones and the reciprocal basis.
	* @param f Field, with filled ghost cells.
	* @param zone Zone of the cell.
	* @param t Index along f's leading dimension.
	* @param i Cell index in x.
	* @param j Cell index in y.
	* @param k Cell index in z.
	* @param dfdy Optional field-aligned df/dy (see derivatives).
	* @return {df/dX, df/dY, df/dZ}.
	*/
	std::array<double, 3> gradient(const Background::Background& bkg,
		const Vectors::ZonedVector<BkgFPType>& f, const int zone,
		const int t, const int i, const int j, const int k,
		const Vectors::ZonedVector<BkgFPType>* dfdy)
	{
		const auto df {derivatives(bkg.get_zones(), f, zone, t, i, j, k,
			dfdy)};
		const std::array<std::array<const Vectors::ZonedVector<BkgFPType>*,
			3>, 3> recp {{
			{&bkg.get_dxdX(), &bkg.get_dxdY(), &bkg.get_dxdZ()},
			{&bkg.get_dydX(), &bkg.get_dydY(), &bkg.get_dydZ()},
			{&bkg.get_dzdX(), &bkg.get_dzdY(), &bkg.get_dzdZ()}}};

		std::array<double, 3> grad {};
		for (int d {}; d < 3; ++d)
		for (int c {}; c < 3; ++c)
		{
			grad[c] += df[d] * (*recp[d][c])(zone, i, j, k);
		}
		return grad;
	}

	/**
	* Parallel gradient b . grad f at a cell of a zone with a zShift table,
	* from df/dy along the field line (aligned_dfdy):
	* b . grad f = (B . grad y / B) df/dy|_zeta, since B . grad zeta = 0.
	*
	* This is HERMES's (BOUT++'s) own parallel gradient. b . grad of the
	* full gradient would be the same if Flan's geometry matched zShift
	* exactly, but it doesn't quite (cell centers and basis vectors from the
	* corners, about 1% per cell on real grids), and then the leftover
	* (B . grad z - dzShift/dy B . grad y) df/dz, from the much larger
	* gradient across the field, swamps the parallel one.
	*
	* @param bkg Background with the field and reciprocal basis.
	* @param dfdy df/dy along the field line, from aligned_dfdy.
	* @param zone Zone of the cell.
	* @param t Index along dfdy's leading dimension (and B's).
	* @param i Cell index in x.
	* @param j Cell index in y.
	* @param k Cell index in z.
	* @return b . grad f.
	*/
	double grad_par(const Background::Background& bkg,
		const Vectors::ZonedVector<BkgFPType>& dfdy, const int zone,
		const int t, const int i, const int j, const int k)
	{
		const Vec3 B {bkg.get_bX()(zone, t, i, j, k), bkg.get_bY()(zone, t,
			i, j, k), bkg.get_bZ()(zone, t, i, j, k)};
		const Vec3 grad_y {bkg.get_dydX()(zone, i, j, k),
			bkg.get_dydY()(zone, i, j, k), bkg.get_dydZ()(zone, i, j, k)};
		return dot(B, grad_y) / std::sqrt(dot(B, B))
			* dfdy(zone, t, i, j, k);
	}

	/**
	* Shifts a periodic row of values in z by a fraction of the row,
	* spectrally: out[k] = f(k + shift), where f is the trigonometric
	* interpolant of the row (period n = row.size()).
	*
	* This is exact for structure resolved by the row, so shifting a
	* field-aligned quantity from one y plane onto the field line through
	* the next doesn't smooth it, the way linear interpolation would. It
	* is how BOUT++'s shifted metric moves data between planes. For an
	* even n, the Nyquist mode is a cosine, which keeps the result real.
	* A plain DFT is used (O(n^2)): rows are a few hundred long at most and
	* this runs once when the background is read.
	*
	* @param row Values at n evenly spaced, periodic points.
	* @param shift Shift in units of the point spacing (any real number).
	* @return The shifted row.
	*/
	std::vector<double> shift_periodic(const std::vector<double>& row,
		const double shift)
	{
		const int n {static_cast<int>(row.size())};
		const double two_pi {2.0 * std::numbers::pi};

		// e^(-2 pi i p / n) for p = 0 to n - 1, so the transforms below
		// need only n sines and cosines instead of n^2
		std::vector<std::complex<double>> w(n);
		for (int p {}; p < n; ++p) w[p] = std::polar(1.0, -two_pi * p / n);

		// Forward transform
		std::vector<std::complex<double>> c(n);
		for (int m {}; m < n; ++m)
		{
			std::complex<double> sum {};
			for (int k {}; k < n; ++k) sum += row[k] * w[(m * k) % n];
			c[m] = sum / static_cast<double>(n);
		}

		// Inverse at the shifted points, with signed frequencies ms so the
		// interpolant is the smoothest one through the points:
		// e^(2 pi i ms (k + shift) / n) = conj(w[ms k mod n]) times a
		// phase for the shift
		for (int m {}; m < n; ++m)
		{
			const int ms {2 * m <= n ? m : m - n};
			c[m] *= std::polar(1.0, two_pi * ms * shift / n);
		}
		std::vector<double> out(n);
		for (int k {}; k < n; ++k)
		{
			double sum {};
			for (int m {}; m < n; ++m)
				sum += (c[m] * std::conj(w[(m * k) % n])).real();
			out[k] = sum;
		}
		return out;
	}

	/**
	* df/dy along the field line, df/dy|_zeta with zeta = z - zShift(x, y)
	* constant on a field line, for every interior cell of the zones with
	* a zShift table (zero elsewhere).
	*
	* On a shifted-metric grid, a field-aligned structure in the y + 1
	* plane is at z + dzShift, often many toroidal cells away, so a central
	* difference in y at fixed z mostly measures unrelated structures.
	* Here the neighboring planes are first shifted onto the field line
	* through each cell (shift_periodic), then differenced centrally (or
	* one-sided at a y face with nothing across it, as in derivatives).
	*
	* This gives both of the y derivatives Flan needs:
	*
	* - df/dy at fixed z, for the gradient, is
	*   df/dy|_zeta - (dzShift/dy) df/dz (see derivatives).
	* - The parallel gradient is (B . grad y / B) df/dy|_zeta (see
	*   grad_par), which is how HERMES itself takes it.
	*
	* The zones' z grids must be evenly spaced and periodic (checked).
	*
	* @param zones The zones.
	* @param f Field, with filled ghost cells.
	* @return df/dy|_zeta, shaped like f.
	*/
	Vectors::ZonedVector<BkgFPType> aligned_dfdy(const Zones::Zones& zones,
		const Vectors::ZonedVector<BkgFPType>& f)
	{
		const int nframes {f.get_dim1()};
		auto out {zoned(zones, nframes)};
		for (int zone {}; zone < zones.get_nzones(); ++zone)
		{
			if (!zones.has_z_shift(zone)) continue;
			const int nx {zones.get_nx(zone)};
			const int ny {zones.get_ny(zone)};
			const int nz {zones.get_nz(zone)};
			const auto& zc {zones.get_z(zone)};
			const double dz {zones.z_period(zone) / nz};
			for (int k {}; k < nz; ++k)
			{
				if (std::abs(zc[k] - zc[0] - k * dz) > 1e-9 * dz * nz)
				{
					Errors::fatal_error("Zone " + std::to_string(zone)
						+ " has a zShift table but its z cells aren't "
						"evenly spaced, so it can't be shifted spectrally.");
				}
			}

			#pragma omp parallel for collapse(2)
			for (int i = 0; i < nx; ++i)
			for (int t = 0; t < nframes; ++t)
			{
				std::vector<double> row(nz);
				auto plane = [&](int j, double shift)
				{
					for (int k {}; k < nz; ++k) row[k] = f(zone, t, i, j, k);
					return shift == 0.0 ? row : shift_periodic(row, shift);
				};
				for (int j {}; j < ny; ++j)
				{
					const auto [lo, hi] = y_stencil(zones, zone, j);
					if (lo == hi) continue;
					const double zs {zones.z_shift(zone, i, j)};

					// The neighboring planes on the field line through
					// each cell of plane j
					const auto f_lo {plane(lo, (zones.z_shift(zone, i, lo)
						- zs) / dz)};
					const auto f_hi {plane(hi, (zones.z_shift(zone, i, hi)
						- zs) / dz)};
					const double dy {zones.center(zone, 1, hi)
						- zones.center(zone, 1, lo)};
					for (int k {}; k < nz; ++k)
					{
						out(zone, t, i, j, k) = static_cast<BkgFPType>(
							(f_hi[k] - f_lo[k]) / dy);
					}
				}
			}
		}
		return out;
	}

	/**
	* Calculates the geometry of every cell from the mesh corners and moves
	* it into the background.
	*
	* The background's zones must already be set to the mesh's zones. In
	* each toroidal plane, a cell's center (R, Z) is the bilinear center of
	* its four corners. The tangent basis vectors along x and y are
	* the differences between the midpoints of opposite cell faces, and
	* along z it's R times the unit vector in phi. From these come the
	* Jacobian J = e_x . (e_y x e_z) and the reciprocal basis, e.g.,
	* grad x = (e_y x e_z) / J.
	*
	* Moved into the background are the cell centers X, Y and Z, the
	* reciprocal basis (dxdX, ..., dzdZ) and |J|, since psi, theta and phi
	* may be a left-handed set but J is used as a cell volume, all with
	* their ghost cells filled, and the cell corners grid_X, grid_Y and
	* grid_Z. The tangent basis and metric aren't kept, to save memory.
	* The background's phi_sign is set to the mesh's.
	*
	* @param mesh The mesh.
	* @param bkg Background to move the geometry into.
	*/
	void calc_geometry(const Mesh& mesh, Background::Background& bkg)
	{
		const Zones::Zones& zones {bkg.get_zones()};
		auto X {zoned(zones, 1)};
		auto Y {zoned(zones, 1)};
		auto Z {zoned(zones, 1)};
		auto J {zoned(zones, 1)};
		std::array<Vectors::ZonedVector<BkgFPType>, 9> recp {};
		for (auto& r : recp) r = zoned(zones, 1);

		// Cell corners have one more value in each dimension and no ghost
		// cells
		auto corners = [&]()
		{
			std::vector<int> nx {}, ny {}, nz {};
			for (int zone {}; zone < zones.get_nzones(); ++zone)
			{
				nx.push_back(zones.get_nx(zone) + 1);
				ny.push_back(zones.get_ny(zone) + 1);
				nz.push_back(zones.get_nz(zone) + 1);
			}
			return Vectors::ZonedVector<BkgFPType>(1, nx, ny, nz, 0);
		};
		auto grid_X {corners()};
		auto grid_Y {corners()};
		auto grid_Z {corners()};
		const double phi_sign {mesh.phi_sign};

		for (int zone {}; zone < zones.get_nzones(); ++zone)
		{
			const int npsi {zones.get_nx(zone)};
			const int ntheta {zones.get_ny(zone)};
			const int nphi {zones.get_nz(zone)};
			const auto& Rc {mesh.R_corners[zone]};
			const auto& Zc {mesh.Z_corners[zone]};

			// Corner (i, j) in meters
			auto R = [&](int i, int j) {return Rc[i * (ntheta + 1) + j];};
			auto Zm = [&](int i, int j) {return Zc[i * (ntheta + 1) + j];};

			for (int i {}; i < npsi; ++i)
			for (int j {}; j < ntheta; ++j)
			{
				const double R_cen {0.25 * (R(i, j) + R(i + 1, j)
					+ R(i, j + 1) + R(i + 1, j + 1))};
				const double Z_cen {0.25 * (Zm(i, j) + Zm(i + 1, j)
					+ Zm(i, j + 1) + Zm(i + 1, j + 1))};
				const double dRdx {0.5 * (R(i + 1, j) + R(i + 1, j + 1)
					- R(i, j) - R(i, j + 1))};
				const double dZdx {0.5 * (Zm(i + 1, j) + Zm(i + 1, j + 1)
					- Zm(i, j) - Zm(i, j + 1))};
				const double dRdy {0.5 * (R(i, j + 1) + R(i + 1, j + 1)
					- R(i, j) - R(i + 1, j))};
				const double dZdy {0.5 * (Zm(i, j + 1) + Zm(i + 1, j + 1)
					- Zm(i, j) - Zm(i + 1, j))};

				for (int k {}; k < nphi; ++k)
				{
					const double phi {zones.get_z(zone)[k]};
					const Vec3 e_x {combine(dRdx, R_hat(phi, phi_sign), dZdx,
						Z_hat)};
					const Vec3 e_y {combine(dRdy, R_hat(phi, phi_sign), dZdy,
						Z_hat)};
					const Vec3 e_z {combine(R_cen, phi_hat(phi, phi_sign), 0.0,
						Z_hat)};
					const double jac {dot(e_x, cross(e_y, e_z))};
					const std::array<Vec3, 3> r {cross(e_y, e_z),
						cross(e_z, e_x), cross(e_x, e_y)};

					X(zone, i, j, k) = R_cen * std::cos(phi);
					Y(zone, i, j, k) = phi_sign * R_cen * std::sin(phi);
					Z(zone, i, j, k) = Z_cen;
					J(zone, i, j, k) = std::abs(jac);
					for (int d {}; d < 3; ++d)
					for (int c {}; c < 3; ++c)
					{
						recp[3 * d + c](zone, i, j, k) = r[d][c] / jac;
					}
				}
			}

			for (int i {}; i <= npsi; ++i)
			for (int j {}; j <= ntheta; ++j)
			for (int k {}; k <= nphi; ++k)
			{
				const double phi {zones.get_grid_z(zone)[k]};
				grid_X(zone, i, j, k) = R(i, j) * std::cos(phi);
				grid_Y(zone, i, j, k) = phi_sign * R(i, j) * std::sin(phi);
				grid_Z(zone, i, j, k) = Zm(i, j);
			}
		}

		// The cell centers and the reciprocal basis vectors are Cartesian
		// vectors, which turn with the toroidal angle in the ghost cells of
		// a wedge
		zones.fill_ghosts(X, Y, phi_sign);
		for (auto* v : {&Z, &J}) zones.fill_ghosts(*v);
		for (int d {}; d < 3; ++d)
		{
			zones.fill_ghosts(recp[3 * d], recp[3 * d + 1], phi_sign);
			zones.fill_ghosts(recp[3 * d + 2]);
		}

		bkg.move_into_X(X);
		bkg.move_into_Y(Y);
		bkg.move_into_Z(Z);
		bkg.move_into_J(J);
		bkg.move_into_grid_X(grid_X);
		bkg.move_into_grid_Y(grid_Y);
		bkg.move_into_grid_Z(grid_Z);
		bkg.move_into_dxdX(recp[0]);
		bkg.move_into_dxdY(recp[1]);
		bkg.move_into_dxdZ(recp[2]);
		bkg.move_into_dydX(recp[3]);
		bkg.move_into_dydY(recp[4]);
		bkg.move_into_dydZ(recp[5]);
		bkg.move_into_dzdX(recp[6]);
		bkg.move_into_dzdY(recp[7]);
		bkg.move_into_dzdZ(recp[8]);
		bkg.set_phi_sign(phi_sign);
	}

	/**
	* Calculates the magnetic field and the gradient of its magnitude, and
	* moves them into the background.
	*
	* The mesh's (Br, Bphi, Bz) are rotated into Cartesian
	* components in each toroidal plane. The field is axisymmetric and
	* doesn't change in time, so every plane and every one of the nframes
	* time frames get the same values. grad |B| is calculated with gradient,
	* so calc_geometry must be called first. Moved into the background are
	* bX, bY, bZ, bmag, gradbX, gradbY and gradbZ, with ghost cells filled.
	*
	* @param mesh The mesh.
	* @param nframes Number of time frames.
	* @param bkg Background to move the field into, with its geometry set.
	*/
	void calc_bfield(const Mesh& mesh, const int nframes,
		Background::Background& bkg)
	{
		if (bkg.get_dxdX().get_data().empty())
		{
			Errors::fatal_error("Axisym::calc_bfield was called before "
				"Axisym::calc_geometry.");
		}

		const Zones::Zones& zones {bkg.get_zones()};
		const double phi_sign {mesh.phi_sign};
		auto bX {zoned(zones, nframes)};
		auto bY {zoned(zones, nframes)};
		auto bZ {zoned(zones, nframes)};
		auto bmag {zoned(zones, nframes)};
		auto gradbX {zoned(zones, nframes)};
		auto gradbY {zoned(zones, nframes)};
		auto gradbZ {zoned(zones, nframes)};

		for (int zone {}; zone < zones.get_nzones(); ++zone)
		{
			const int ntheta {zones.get_ny(zone)};
			for (int i {}; i < zones.get_nx(zone); ++i)
			for (int j {}; j < ntheta; ++j)
			{
				const int c {i * ntheta + j};
				const double BR {mesh.Br[zone][c]};
				const double BZ {mesh.Bz[zone][c]};
				const double Bphi {mesh.Bphi[zone][c]};
				for (int k {}; k < zones.get_nz(zone); ++k)
				{
					const double phi {zones.get_z(zone)[k]};
					const Vec3 B {combine(BR, R_hat(phi, phi_sign), Bphi,
						phi_hat(phi, phi_sign))};
					for (int t {}; t < nframes; ++t)
					{
						bX(zone, t, i, j, k) = B[0];
						bY(zone, t, i, j, k) = B[1];
						bZ(zone, t, i, j, k) = BZ;
						bmag(zone, t, i, j, k) = std::sqrt(B[0] * B[0]
							+ B[1] * B[1] + BZ * BZ);
					}
				}
			}
		}
		zones.fill_ghosts(bX, bY, phi_sign);
		for (auto* v : {&bZ, &bmag}) zones.fill_ghosts(*v);

		// Gradient of |B|, the same in every time frame
		for (int zone {}; zone < zones.get_nzones(); ++zone)
		for (int i {}; i < zones.get_nx(zone); ++i)
		for (int j {}; j < zones.get_ny(zone); ++j)
		for (int k {}; k < zones.get_nz(zone); ++k)
		{
			const auto grad {gradient(bkg, bmag, zone, 0, i, j, k)};
			for (int t {}; t < nframes; ++t)
			{
				gradbX(zone, t, i, j, k) = grad[0];
				gradbY(zone, t, i, j, k) = grad[1];
				gradbZ(zone, t, i, j, k) = grad[2];
			}
		}
		zones.fill_ghosts(gradbX, gradbY, phi_sign);
		zones.fill_ghosts(gradbZ);

		bkg.move_into_bX(bX);
		bkg.move_into_bY(bY);
		bkg.move_into_bZ(bZ);
		bkg.move_into_bmag(bmag);
		bkg.move_into_gradbX(gradbX);
		bkg.move_into_gradbY(gradbY);
		bkg.move_into_gradbZ(gradbZ);
	}

	/**
	* Prints how a code's own ghost cells in a field (e.g., SOLEDGE3X's)
	* differ from what
	* Zones::Zones::fill_ghosts gives, for the first entry of the leading
	* dimension (e.g., the first time frame).
	*
	* The ghost cells are split into those past a single face shared with
	* another zone or periodic (edges), those past two or three such faces
	* (corners), corners averaged from more than one cell (e.g., at the
	* X-point) and those past an absorbing, core or wall face (boundaries). In
	* each group, the largest difference is printed relative to the largest
	* |value| of the field. Codes usually fill boundary ghost cells with
	* their boundary conditions, so those are expected to differ.
	*
	* @param zones The zones.
	* @param v Field holding the code's ghost cells.
	* @param name Name of the field (for the message).
	*/
	void report_ghosts(const Zones::Zones& zones,
		const Vectors::ZonedVector<BkgFPType>& v, const std::string& name)
	{
		Vectors::ZonedVector<BkgFPType> ours {v};
		zones.fill_ghosts(ours);

		double scale {};
		for (const double val : v.get_data())
			scale = std::max(scale, std::abs(val));
		if (scale == 0.0) scale = 1.0;

		// edges, corners, averaged, boundaries
		std::array<double, 4> max_diff {};
		for (const auto& gc : zones.get_ghost_cells())
		{
			const Zones::Cell& c {gc.ghost};
			const std::array<int, 3> idx {c.i, c.j, c.k};
			int nout {};
			bool boundary {};
			for (int d {}; d < 3; ++d)
			{
				const bool below {idx[d] < 0};
				if (!below && idx[d] < zones.get_n(c.zone, d)) continue;
				++nout;
				const auto face {static_cast<Zones::Face>(2 * d + !below)};
				if (zones.get_face(c.zone, face) < 0) boundary = true;
			}
			const int group {boundary ? 3 : gc.sources.size() > 1 ? 2
				: nout > 1 ? 1 : 0};
			const double diff {std::abs(v(c.zone, 0, c.i, c.j, c.k)
				- ours(c.zone, 0, c.i, c.j, c.k))};
			max_diff[group] = std::max(max_diff[group], diff / scale);
		}

		std::cout << "  " << name << ": edges " << max_diff[0]
			<< ", corners " << max_diff[1] << ", averaged corners "
			<< max_diff[2] << ", boundaries " << max_diff[3] << '\n';
	}

	/**
	* Calculates the electric field and the background ion flow, and moves
	* them into the background.
	*
	* The background must already have its geometry and magnetic field,
	* and should have its ion temperature, which is taken as 0 if it
	* hasn't been loaded. E = -grad vp, calculated with gradient. The flow is the parallel flow plus the
	* E x B and grad B drifts, the same as for Gkeyll backgrounds:
	* u = u_par b + E x B / B^2 + Ti (B x grad B) / B^3, where the grad B
	* drift uses v_perp^2 = 2 e Ti / m_i for a Maxwellian, so the mass
	* cancels. Moved into the background are eX, eY, eZ, emag, uX, uY and
	* uZ, with ghost cells filled.
	*
	* @param bkg Background with zones, geometry and magnetic field.
	* @param vp Plasma potential (V), with ghost cells filled.
	* @param upar Parallel ion flow (m/s).
	* @param epar Optional extra electric field along B (V/m), e.g. the
	*        inductive -dApar/dt of an electromagnetic run, added to E's
	*        parallel part. It doesn't change the E x B drift.
	*/
	void calc_efield_and_flows(Background::Background& bkg,
		const Vectors::ZonedVector<BkgFPType>& vp,
		const Vectors::ZonedVector<BkgFPType>& upar,
		const Vectors::ZonedVector<BkgFPType>* epar)
	{
		const Zones::Zones& zones {bkg.get_zones()};
		const int nframes {vp.get_dim1()};
		auto eX {zoned(zones, nframes)};
		auto eY {zoned(zones, nframes)};
		auto eZ {zoned(zones, nframes)};
		auto emag {zoned(zones, nframes)};
		auto uX {zoned(zones, nframes)};
		auto uY {zoned(zones, nframes)};
		auto uZ {zoned(zones, nframes)};

		const auto& bX {bkg.get_bX()};
		const auto& bY {bkg.get_bY()};
		const auto& bZ {bkg.get_bZ()};
		const auto& gX {bkg.get_gradbX()};
		const auto& gY {bkg.get_gradbY()};
		const auto& gZ {bkg.get_gradbZ()};
		const auto& ti {bkg.get_ti()};
		const bool have_ti {!ti.get_data().empty()};

		// Make sure the ghost cell sources are found before the threads
		// start, since that isn't thread safe
		zones.get_ghost_cells();

		// On a shifted-metric grid, y derivatives of the potential are
		// taken along the field line (see aligned_dfdy)
		bool any_shift {};
		for (int zone {}; zone < zones.get_nzones(); ++zone)
			any_shift = any_shift || zones.has_z_shift(zone);
		Vectors::ZonedVector<BkgFPType> dvp_dy {};
		if (any_shift) dvp_dy = aligned_dfdy(zones, vp);
		const auto* dfdy {any_shift ? &dvp_dy : nullptr};

		for (int zone {}; zone < zones.get_nzones(); ++zone)
		{
			#pragma omp parallel for
			for (int i = 0; i < zones.get_nx(zone); ++i)
			for (int t {}; t < nframes; ++t)
			for (int j {}; j < zones.get_ny(zone); ++j)
			for (int k {}; k < zones.get_nz(zone); ++k)
			{
				const auto grad {gradient(bkg, vp, zone, t, i, j, k, dfdy)};
				Vec3 E {-grad[0], -grad[1], -grad[2]};
				const Vec3 B {bX(zone, t, i, j, k), bY(zone, t, i, j, k),
					bZ(zone, t, i, j, k)};

				// On a shifted-metric grid, E's parallel part is
				// -grad_par phi taken purely along the field line (see
				// grad_par), the way HERMES takes it
				if (dfdy && zones.has_z_shift(zone))
				{
					const double Bm {std::sqrt(dot(B, B))};
					const double dE {-grad_par(bkg, *dfdy, zone, t, i, j, k)
						- dot(E, B) / Bm};
					for (int c {}; c < 3; ++c) E[c] += dE * B[c] / Bm;
				}
				if (epar)
				{
					const double Bm {std::sqrt(dot(B, B))};
					const double e {(*epar)(zone, t, i, j, k)};
					for (int c {}; c < 3; ++c) E[c] += e * B[c] / Bm;
				}
				const Vec3 gradB {gX(zone, t, i, j, k), gY(zone, t, i, j, k),
					gZ(zone, t, i, j, k)};
				const double B2 {dot(B, B)};
				const double Bmag {std::sqrt(B2)};
				const double Ti {have_ti ? ti(zone, t, i, j, k) : 0.0};

				const Vec3 ExB {cross(E, B)};
				const Vec3 BxgradB {cross(B, gradB)};
				const double u_par {upar(zone, t, i, j, k)};
				Vec3 u {};
				for (int c {}; c < 3; ++c)
				{
					u[c] = u_par * B[c] / Bmag + ExB[c] / B2
						+ Ti * BxgradB[c] / (B2 * Bmag);
				}

				eX(zone, t, i, j, k) = E[0];
				eY(zone, t, i, j, k) = E[1];
				eZ(zone, t, i, j, k) = E[2];
				emag(zone, t, i, j, k) = std::sqrt(dot(E, E));
				uX(zone, t, i, j, k) = u[0];
				uY(zone, t, i, j, k) = u[1];
				uZ(zone, t, i, j, k) = u[2];
			}
		}
		const double phi_sign {bkg.get_phi_sign()};
		zones.fill_ghosts(eX, eY, phi_sign);
		zones.fill_ghosts(uX, uY, phi_sign);
		for (auto* v : {&eZ, &emag, &uZ}) zones.fill_ghosts(*v);

		bkg.move_into_eX(eX);
		bkg.move_into_eY(eY);
		bkg.move_into_eZ(eZ);
		bkg.move_into_emag(emag);
		bkg.move_into_uX(uX);
		bkg.move_into_uY(uY);
		bkg.move_into_uZ(uZ);
	}

	/**
	* Calculates the parallel conductive heat fluxes of the electrons and
	* the main ions from the background's n and T, and moves them into
	* the background (q_par_e and q_par_i), for backgrounds that don't store
	* them (HERMES).
	*
	* For each species, q_par = -kappa_par grad_par T, with Braginskii's
	* kappa_par = coef n e T tau / m in the form Hermes-3 uses (see
	* Braginskii::kappa_par), and tau = 1 / nu from all of the species'
	* collisions, as in Hermes-3's default ("multispecies") conduction:
	* nu_ee + nu_ei for the electrons and nu_ii + nu_ie for the ions, with
	* nu_ie = nu_ei (me / mi) (ne / ni) from momentum conservation. The
	* optional flux limiter is applied if c.limit_alpha > 0.
	*
	* The plasma is taken to be one ion species with ni = ne / Z
	* (quasineutral; for Z = 1, ni = ne), since the background only holds
	* ne. A run with several ion species (or impurities that matter for
	* the collisions) would need every species' density and would add
	* their collisions to each tau.
	*
	* On zones with a zShift table, grad_par is the parallel gradient along
	* the field line (aligned_dfdy and grad_par), not b . grad T, which
	* picks up part of the much larger gradient across the field. Zones
	* without one (hermes_field_aligned = "off") use b . grad T. q_par is
	* positive along B, the same as the parallel flow and the SOLEDGE3X
	* heat fluxes.
	*
	* At a target (a wall face in y) HERMES conducts nothing through the
	* face: the sheath boundary condition carries the energy there instead.
	* So, as in the run, the target cell's q_par is the average of its two
	* y faces' conducted fluxes: zero at the target and, at the face to its
	* neighbor, -kappa grad_par T with kappa the average of the two cells'
	* and the gradient between them (the one-sided difference). This keeps
	* the thermal force in the target cell from using a gradient that the
	* run sets by its sheath condition, not by conduction.
	*
	* @param bkg Background with zones, geometry, B, ne, te and ti (ghost
	*        cells filled).
	* @param c Conduction coefficients, flux limiter and main ion.
	*/
	void calc_conductive_heat_fluxes(Background::Background& bkg,
		const Conduction& c)
	{
		const Zones::Zones& zones {bkg.get_zones()};

		// Make sure the ghost cell sources are found before the threads
		// start, since that isn't thread safe
		zones.get_ghost_cells();

		const auto& ne {bkg.get_ne()};
		const auto& te {bkg.get_te()};
		const auto& ti {bkg.get_ti()};
		const int nframes {ne.get_dim1()};
		const auto dte_dy {aligned_dfdy(zones, te)};
		const auto dti_dy {aligned_dfdy(zones, ti)};
		auto q_e {zoned(zones, nframes)};
		auto q_i {zoned(zones, nframes)};

		// Each cell's (limited) kappa, kept for the target cells, which
		// use the average of theirs and their neighbor's
		auto kappa_e {zoned(zones, nframes)};
		auto kappa_i {zoned(zones, nframes)};

		const double mi {Braginskii::mp * c.ion_aa};
		const double me {Braginskii::me};

		for (int zone {}; zone < zones.get_nzones(); ++zone)
		{
			#pragma omp parallel for
			for (int i = 0; i < zones.get_nx(zone); ++i)
			for (int t {}; t < nframes; ++t)
			for (int j {}; j < zones.get_ny(zone); ++j)
			for (int k {}; k < zones.get_nz(zone); ++k)
			{
				const double n_e {ne(zone, t, i, j, k)};
				const double n_i {n_e / c.ion_z};	// quasineutral
				const double t_e {te(zone, t, i, j, k)};
				const double t_i {ti(zone, t, i, j, k)};

				const double nu_ei {Braginskii::nu_ei(n_e, n_i, t_e, t_i,
					c.ion_z, c.ion_aa)};
				const double nu_e {Braginskii::nu_ee(n_e, t_e) + nu_ei};
				const double nu_ie {n_i > 0.0
					? nu_ei * (me / mi) * (n_e / n_i) : 0.0};
				const double nu_i {Braginskii::nu_ii(n_i, t_i, c.ion_z,
					c.ion_aa, n_i, t_i, c.ion_z, c.ion_aa) + nu_ie};

				// Along the field line on a shifted-metric zone, otherwise
				// b . grad T
				double g_e {}, g_i {};
				if (zones.has_z_shift(zone))
				{
					g_e = grad_par(bkg, dte_dy, zone, t, i, j, k);
					g_i = grad_par(bkg, dti_dy, zone, t, i, j, k);
				}
				else
				{
					const Vec3 B {bkg.get_bX()(zone, t, i, j, k),
						bkg.get_bY()(zone, t, i, j, k),
						bkg.get_bZ()(zone, t, i, j, k)};
					const double Bm {std::sqrt(dot(B, B))};
					g_e = dot(B, gradient(bkg, te, zone, t, i, j, k)) / Bm;
					g_i = dot(B, gradient(bkg, ti, zone, t, i, j, k)) / Bm;
				}
				const double k_e {Braginskii::limit_kappa(
					Braginskii::kappa_par(c.coef_e, n_e, t_e, nu_e, me),
					g_e, n_e, t_e, me, c.limit_alpha)};
				const double k_i {Braginskii::limit_kappa(
					Braginskii::kappa_par(c.coef_i, n_i, t_i, nu_i, mi),
					g_i, n_i, t_i, mi, c.limit_alpha)};

				q_e(zone, t, i, j, k) = static_cast<BkgFPType>(-k_e * g_e);
				q_i(zone, t, i, j, k) = static_cast<BkgFPType>(-k_i * g_i);
				kappa_e(zone, t, i, j, k) = static_cast<BkgFPType>(k_e);
				kappa_i(zone, t, i, j, k) = static_cast<BkgFPType>(k_i);
			}
		}

		// Target cells: no conduction through the wall face, so q_par is
		// half the flux through the face to the neighbor. That flux uses
		// the one-sided gradient already in q (= -kappa g there) with the
		// average kappa of the two cells. A zone one cell long has no
		// neighbor in y to average with: with walls at both ends it
		// conducts nothing, with one it keeps the one-sided value.
		for (int zone {}; zone < zones.get_nzones(); ++zone)
		{
			const int ny {zones.get_ny(zone)};
			const bool wall_lo {zones.get_face(zone, Zones::ymin)
				== Zones::face_wall};
			const bool wall_hi {zones.get_face(zone, Zones::ymax)
				== Zones::face_wall};
			for (const auto [j, nb, wall] : {std::tuple {0, 1, wall_lo},
				std::tuple {ny - 1, ny - 2, wall_hi}})
			{
				if (!wall || (ny == 1 && !(wall_lo && wall_hi))) continue;
				const bool both {ny == 1};
				#pragma omp parallel for
				for (int i = 0; i < zones.get_nx(zone); ++i)
				for (int t {}; t < nframes; ++t)
				for (int k {}; k < zones.get_nz(zone); ++k)
				{
					if (both)
					{
						q_e(zone, t, i, j, k) = 0;
						q_i(zone, t, i, j, k) = 0;
						continue;
					}
					auto face = [&](Vectors::ZonedVector<BkgFPType>& q,
						const Vectors::ZonedVector<BkgFPType>& kap)
					{
						const double kc {kap(zone, t, i, j, k)};
						const double kn {kap(zone, t, i, nb, k)};
						if (kc == 0.0) return;
						const double g {-q(zone, t, i, j, k) / kc};
						q(zone, t, i, j, k) = static_cast<BkgFPType>(
							-0.5 * 0.5 * (kc + kn) * g);
					};
					face(q_e, kappa_e);
					face(q_i, kappa_i);
				}
			}
		}
		zones.fill_ghosts(q_e);
		zones.fill_ghosts(q_i);
		bkg.move_into_q_par_e(q_e);
		bkg.move_into_q_par_i(q_i);
	}
}
