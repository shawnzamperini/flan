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
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include "axisym_mesh.h"
#include "background.h"
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
	* @return {df/dx, df/dy, df/dz}.
	*/
	std::array<double, 3> derivatives(const Zones::Zones& zones,
		const Vectors::ZonedVector<BkgFPType>& f, const int zone,
		const int t, const int i, const int j, const int k)
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
	* @return {df/dX, df/dY, df/dZ}.
	*/
	std::array<double, 3> gradient(const Background::Background& bkg,
		const Vectors::ZonedVector<BkgFPType>& f, const int zone,
		const int t, const int i, const int j, const int k)
	{
		const auto df {derivatives(bkg.get_zones(), f, zone, t, i, j, k)};
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

		for (auto* v : {&X, &Y, &Z, &J}) zones.fill_ghosts(*v);
		for (auto& r : recp) zones.fill_ghosts(r);

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
		for (auto* v : {&bX, &bY, &bZ, &bmag}) zones.fill_ghosts(*v);

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
		for (auto* v : {&gradbX, &gradbY, &gradbZ}) zones.fill_ghosts(*v);

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
	*/
	void calc_efield_and_flows(Background::Background& bkg,
		const Vectors::ZonedVector<BkgFPType>& vp,
		const Vectors::ZonedVector<BkgFPType>& upar)
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

		for (int zone {}; zone < zones.get_nzones(); ++zone)
		{
			#pragma omp parallel for
			for (int i = 0; i < zones.get_nx(zone); ++i)
			for (int t {}; t < nframes; ++t)
			for (int j {}; j < zones.get_ny(zone); ++j)
			for (int k {}; k < zones.get_nz(zone); ++k)
			{
				const auto grad {gradient(bkg, vp, zone, t, i, j, k)};
				const Vec3 E {-grad[0], -grad[1], -grad[2]};
				const Vec3 B {bX(zone, t, i, j, k), bY(zone, t, i, j, k),
					bZ(zone, t, i, j, k)};
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
		for (auto* v : {&eX, &eY, &eZ, &emag, &uX, &uY, &uZ})
			zones.fill_ghosts(*v);

		bkg.move_into_eX(eX);
		bkg.move_into_eY(eY);
		bkg.move_into_eZ(eZ);
		bkg.move_into_emag(emag);
		bkg.move_into_uX(uX);
		bkg.move_into_uY(uY);
		bkg.move_into_uZ(uZ);
	}
}
