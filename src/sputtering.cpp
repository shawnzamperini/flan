/**
* @file sputtering.cpp
*
* @brief Sputtering of wall surfaces by the background plasma.
*
* There are two parts so far: the sputtering yield, and the wall surfaces
* that get sputtered.
*
* Yield (YieldTable)
* ------------------
* The sputtering yield Y (sputtered atoms per incident ion) depends on the
* incident ion's energy E and angle from the surface normal. It is read
* from a table made with RustBCA by python/rustbca_yields.py (see
* data/sputtering/) and bilinearly interpolated, linear in log10(E) and in
* angle.
*
* Wall surfaces (Surface)
* -----------------------
* SOLEDGE3X doesn't fit its mesh to the wall. Its mesh runs past the wall,
* and a mask, chi, marks the cells inside the wall (chi = 1) and in the
* plasma (chi = 0). So the surface the plasma touches is wherever a plasma
* cell sits next to a wall cell. Seen in the poloidal (R, Z) plane, that's
* a staircase of cell faces, while the real wall, which SOLEDGE3X's mesh
* file also has as an outline (a polygon, see Background::WallPolygon),
* runs smoothly through it:
*
*                      chi = 1 (wall)
*              +----+----+----+----+
*              |    |    | \  |    |          \  = wall outline
*              +----+----+--\-+----+          ## = staircase of faces
*              |    |    |   \|    |               between chi = 0 and
*              +----+----##===\====+               chi = 1 cells
*              |    |    #|    \   |
*              +----#====#+-----\--+
*              |    #     |      \ |
*              +----+-----+-------\+
*                      chi = 0 (plasma)
*
* Neither is good enough alone:
*
* - The staircase faces are where particles can actually be: a particle
*   in a chi = 1 cell is absorbed (s3x_chi_absorb), so a sputtered atom
*   has to start in the plasma cell next to the face. The face also tells
*   us which plasma cell's density, temperature, etc. to use. But the
*   faces point along the grid (each is a psi or theta face), not
*   perpendicular to the wall, and there's more of them than there is
*   wall: a staircase is longer than the slope it follows, up to sqrt(2)
*   times for a wall at 45 degrees to the grid.
*
* - The wall outline has the right direction and length, but it cuts
*   through cells, often through chi = 1 cells, so it doesn't say where in
*   the grid to start a particle or which cell's plasma to use.
*
* So each staircase face becomes a surface element (SurfaceElement) that
* keeps the face's position and plasma cell, but takes its direction from
* the wall outline:
*
* 1. Normal. The element's unit normal n_wall is the normal of the wall
*    outline segment nearest to the middle of the face, pointed into the
*    plasma. This is the normal to use for the ion flux onto the wall
*    (which goes as b . n) and for the direction sputtered atoms leave in.
*
* 2. Area. The face's own area A_face is projected onto the wall:
*    A = A_face |n_face . n_wall|, where n_face is the face's own normal.
*    This makes the areas add up to the wall's, not the longer
*    staircase's. Take a stretch of wall of length L at angle a to one
*    grid direction. The steps under it have faces along that direction
*    adding up to L |cos(a)|, whose normals are at angle a to n_wall, and
*    faces along the other direction adding up to L |sin(a)|, whose
*    normals are at 90 - a to it. Projected, they add up to
*    L cos^2(a) + L sin^2(a) = L.
*
* Where the nearest wall segment is far from a face (more than
* max_wall_dist cell widths), the outline doesn't describe that part of
* the wall. That happens, e.g., where chi marks something the outline
* doesn't include. Then the element falls back to the face's own normal
* and area (SurfaceElement::fallback), which is the best that's known
* there. The same happens for every face if the mesh file has no
* outline.
*
* Everything is done per toroidal cell, so each face appears once in
* every toroidal cell of its zone, with that cell's range of phi.
*
* Sputtered source (Source)
* -------------------------
* For every surface element and background time frame, the rate atoms are
* sputtered off it (atoms/s) is
*
*     S = Gamma x Y_eff x A
*
* using the plasma in the element's plasma cell (no interpolation):
*
* - Gamma is the flux of background ions onto the wall (m^-2 s^-1). Ions
*   move along the field with the parallel flow u_par and across it with
*   the E x B drift v_E, so the velocity into the wall is
*   -(u_par b_hat + v_E) . n_wall (n_wall points into the plasma), and
*   Gamma = n_i max(0, that). The two parts are perpendicular, separate
*   pieces of the motion, so they're kept separately too (flux_par,
*   flux_exb). The background's flow also has the grad-B drift, which
*   isn't used: it's mostly a magnetization flow that doesn't carry
*   particles onto the wall. sput_flux_opt = "parallel" leaves out the
*   E x B part, and sput_bohm_floor = "on" makes the parallel flow onto
*   the wall at least the sound speed.
*
* - Y_eff is the sputtering yield averaged over the energies the ions hit
*   the wall with, at the fixed angle sput_angle_val. Only the ions
*   sputter: electrons are far too light to knock a W atom loose. Te only
*   comes in through the sheath potential, which the electrons set and
*   the ions are accelerated by. An ion falls through
*   the sheath and presheath potential, gaining E_sheath (sheath_energy),
*   on top of its thermal energy. The thermal energy of the ions hitting
*   a wall is that of a Maxwellian weighted by their speed into it, which
*   is a Gamma(2, Ti) distribution, p(eps) = eps / Ti^2 exp(-eps / Ti)
*   (mean 2 Ti). Y_eff = int Y(E_sheath + eps) p(eps) d eps
*   (effective_yield). Averaging matters near threshold, where Y rises
*   steeply: Y at the mean energy can be far from the mean of Y.
*
* - A is the element's area (see Surface).
*
* Starting the sputtered atoms (Source::sample_birth)
* ----------------------------------------------------
* With imp_source_opt = "sputter", every particle is a W atom sputtered off
* the wall:
*
* 1. Which element and when. Over the run, frame t stands for the times
*    nearest to it, a window of length dt_t, so element e sputters
*    S_e(t) dt_t atoms in it. A (frame, element) pair is picked with that
*    weight, and the time uniformly in the frame's window.
* 2. Where. Uniformly along the element's face and its range of phi, just
*    inside the plasma cell (a thousandth of a cell from the face).
* 3. Which ion did it. Only the fast ions sputter near threshold, so the
*    energy of the ion is picked from the ions that sputter, whose
*    energies go as Y(E) p(eps), the same integrand as Y_eff
*    (sample_impact_energy).
* 4. How fast. The sputtered atom's energy is from the Thompson
*    distribution, f(E) ~ E / (E + E_b)^3, with E_b the target's surface
*    binding energy from the yield table. It's cut off at the most energy
*    the ion can give a target atom in one collision,
*    Lambda E_ion, Lambda = 4 m1 m2 / (m1 + m2)^2 (0.043 for D on W)
*    (sample_thompson).
* 5. Which way. A cosine distribution about the element's normal, which
*    is the most common angular distribution of sputtered atoms
*    (cosine_direction).
*
* The atom starts neutral. Its weight, the atoms/s it stands for, is
* source_rate() / imp_num, where source_rate() is the total rate averaged
* over the frames' windows. Atoms are picked in proportion to where they're
* sputtered, so they all have the same weight.
*
* One artifact: the atom starts at a staircase face, but its direction is
* about the wall outline's normal, which can be tilted from the face's. A
* direction that points back through the face goes into the wall cell and
* is absorbed on its first step (with s3x_chi_absorb = "on"). For a cosine
* distribution, that's 0.8% of atoms where the two normals are 10 degrees
* apart, 7% at 30 degrees and 15% at 45 degrees.
*/

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <limits>
#include <map>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>

#include "background.h"
#include "constants.h"
#include "errors.h"
#include "hdf5_utils.h"
#include "options.h"
#include "read_s3xe.h"
#include "sputtering.h"
#include "zones.h"


namespace Sputtering
{
	namespace
	{
		/**
		* Finds the interval of a sorted grid to interpolate in.
		*
		* @param grid Strictly increasing values, at least two.
		* @param value Value to find. Must be within the grid.
		* @return Index i of the lower end of the interval, so the value is
		*         between grid[i] and grid[i + 1], from 0 to size - 2.
		*/
		std::size_t lower_index(const std::vector<double>& grid,
			const double value)
		{
			const auto upper {std::upper_bound(grid.begin(), grid.end(),
				value)};
			const std::size_t i {static_cast<std::size_t>(
				std::max<std::ptrdiff_t>(upper - grid.begin() - 1, 0))};
			return std::min(i, grid.size() - 2);
		}

		/**
		* Checks that a grid has at least two values and is strictly
		* increasing, stopping the run if not.
		*
		* @param grid Values to check.
		* @param name Name of the grid, for the error message.
		* @param path File the grid is from, for the error message.
		*/
		void check_grid(const std::vector<double>& grid,
			const std::string& name, const std::string& path)
		{
			if (grid.size() < 2)
			{
				Errors::fatal_error("The sputtering yield table " + path
					+ " has " + std::to_string(grid.size()) + " " + name
					+ " value(s). At least 2 are needed to interpolate.");
			}
			for (std::size_t i {1}; i < grid.size(); ++i)
			{
				if (!(grid[i] > grid[i - 1]))
				{
					Errors::fatal_error("The " + name + " values in the "
						"sputtering yield table " + path + " are not "
						"strictly increasing.");
				}
			}
		}
	}

	/**
	* Reads a sputtering yield table.
	*
	* The file is an HDF5 file made by python/rustbca_yields.py, with
	* /energy (eV), /angle (degrees from the surface normal), /yield (atoms
	* per ion, shape (energy, angle)), /projectile/Z, /projectile/mass_amu,
	* /target/Z, /target/mass_amu and /target/Es (eV). Its other datasets
	* and attributes aren't needed. The run stops if a dataset is missing,
	* the yield has the wrong shape, there are fewer than two energies or
	* angles, the energies aren't positive and strictly increasing, the
	* angles aren't strictly increasing between 0 and 90, or a yield is
	* negative or not finite.
	*
	* @param path Path to the table.
	*/
	YieldTable::YieldTable(const std::string& path)
		: m_path {path}
	{
		const Hdf5::File file(path);
		const std::string where {"The sputtering yield table " + path};

		const std::vector<double> energy {Hdf5::read_dataset<double>(file,
			"/energy").data};
		m_angle = Hdf5::read_dataset<double>(file, "/angle").data;
		const Hdf5::Dataset<double> yield {Hdf5::read_dataset<double>(file,
			"/yield")};

		check_grid(energy, "energy", path);
		check_grid(m_angle, "angle", path);
		if (energy.front() <= 0.0)
		{
			Errors::fatal_error(where + " has energies <= 0.");
		}
		if (m_angle.front() < 0.0 || m_angle.back() > 90.0)
		{
			Errors::fatal_error(where + " has angles outside 0 to 90 "
				"degrees.");
		}

		const std::vector<hsize_t> shape {energy.size(), m_angle.size()};
		if (yield.shape != shape)
		{
			Errors::fatal_error(where + " has a yield with the wrong shape. "
				"It should be (" + std::to_string(shape[0]) + ", "
				+ std::to_string(shape[1]) + "), the number of energies "
				"and angles.");
		}
		for (const double y : yield.data)
		{
			if (!std::isfinite(y) || y < 0.0)
			{
				Errors::fatal_error(where + " has a negative or non-finite "
					"yield.");
			}
		}
		m_yield = yield.data;

		m_log_energy.reserve(energy.size());
		for (const double e : energy) m_log_energy.push_back(std::log10(e));
		m_energy_min = energy.front();
		m_energy_max = energy.back();

		m_projectile_Z = Hdf5::read_scalar<double>(file, "/projectile/Z");
		m_projectile_mass_amu = Hdf5::read_scalar<double>(file,
			"/projectile/mass_amu");
		m_target_Z = Hdf5::read_scalar<double>(file, "/target/Z");
		m_target_mass_amu = Hdf5::read_scalar<double>(file,
			"/target/mass_amu");
		m_target_Es = Hdf5::read_scalar<double>(file, "/target/Es");
	}

	/**
	* Returns the sputtering yield at an incident energy and angle.
	*
	* The yield is bilinearly interpolated, linear in log10(E) and in
	* angle. Y itself is interpolated, not log(Y), since the yield is 0
	* below threshold. Below the lowest energy in the table the yield is 0,
	* since the table starts below threshold. Above the highest energy, the
	* highest energy is used. The angle is clamped to the table's angles
	* (e.g., 90 uses the largest angle, 89).
	*
	* @param energy Incident energy (eV).
	* @param angle Incident angle from the surface normal (degrees).
	* @return Sputtering yield (atoms per ion).
	*/
	double YieldTable::get_yield(double energy, double angle) const
	{
		if (!(energy >= m_energy_min)) return 0.0;

		const double log_energy {std::min(std::log10(energy),
			m_log_energy.back())};
		angle = std::clamp(angle, m_angle.front(), m_angle.back());

		const std::size_t i {lower_index(m_log_energy, log_energy)};
		const std::size_t j {lower_index(m_angle, angle)};
		const double fe {(log_energy - m_log_energy[i])
			/ (m_log_energy[i + 1] - m_log_energy[i])};
		const double fa {(angle - m_angle[j])
			/ (m_angle[j + 1] - m_angle[j])};

		const std::size_t na {m_angle.size()};
		const double y00 {m_yield[i * na + j]};
		const double y01 {m_yield[i * na + j + 1]};
		const double y10 {m_yield[(i + 1) * na + j]};
		const double y11 {m_yield[(i + 1) * na + j + 1]};

		return (1.0 - fe) * ((1.0 - fa) * y00 + fa * y01)
			+ fe * ((1.0 - fa) * y10 + fa * y11);
	}

	/**
	* @name Accessors
	* The file the table was read from, its energy range (eV), and the
	* projectile and target it's for: atomic numbers, masses (amu) and the
	* target's surface binding energy (eV).
	*/
	///@{
	const std::string& YieldTable::path() const {return m_path;}
	double YieldTable::energy_min() const {return m_energy_min;}
	double YieldTable::energy_max() const {return m_energy_max;}
	double YieldTable::projectile_Z() const {return m_projectile_Z;}
	double YieldTable::projectile_mass_amu() const
		{return m_projectile_mass_amu;}
	double YieldTable::target_Z() const {return m_target_Z;}
	double YieldTable::target_mass_amu() const {return m_target_mass_amu;}
	double YieldTable::target_Es() const {return m_target_Es;}
	///@}

	namespace
	{
		/// A point or vector in the poloidal (R, Z) plane (m)
		struct RZ
		{
			double R {};
			double Z {};
		};

		RZ operator-(const RZ& a, const RZ& b) {return {a.R - b.R, a.Z - b.Z};}
		RZ operator+(const RZ& a, const RZ& b) {return {a.R + b.R, a.Z + b.Z};}
		RZ operator*(const double s, const RZ& a) {return {s * a.R, s * a.Z};}
		double dot(const RZ& a, const RZ& b) {return a.R * b.R + a.Z * b.Z;}
		double norm(const RZ& a) {return std::hypot(a.R, a.Z);}

		/**
		* Returns a unit vector perpendicular to the line from a to b, which
		* must be different points. Which of the two perpendicular
		* directions comes back is arbitrary; callers point it where they
		* need it with orient_like.
		*/
		RZ unit_normal(const RZ& a, const RZ& b)
		{
			const RZ t {b - a};
			const double len {norm(t)};
			return {t.Z / len, -t.R / len};
		}

		/// Returns n, flipped if needed to point the same way as ref
		RZ orient_like(const RZ& n, const RZ& ref)
		{
			return dot(n, ref) < 0.0 ? -1.0 * n : n;
		}

		/// A straight piece of a wall outline, from a to b
		struct Segment
		{
			RZ a {};
			RZ b {};
		};

		/**
		* Splits the wall outlines into their straight segments, leaving
		* out any of zero length (e.g., a closed outline whose last point
		* repeats its first gives no extra segment).
		*/
		std::vector<Segment> wall_segments(
			const std::vector<Background::WallPolygon>& walls)
		{
			std::vector<Segment> segments {};
			for (const auto& wall : walls)
			{
				for (std::size_t p {1}; p < wall.R.size(); ++p)
				{
					const Segment s {{wall.R[p - 1], wall.Z[p - 1]},
						{wall.R[p], wall.Z[p]}};
					if (norm(s.b - s.a) > 0.0) segments.push_back(s);
				}
			}
			return segments;
		}

		/// Distance from point p to the closest point of segment s
		double distance_to_segment(const RZ& p, const Segment& s)
		{
			const RZ ab {s.b - s.a};
			const double t {std::clamp(dot(p - s.a, ab) / dot(ab, ab), 0.0,
				1.0)};
			return norm(p - (s.a + t * ab));
		}

		/**
		* The plasma cell's neighbor across one of its faces, and the two
		* corners of that face, all as (i, j) index offsets from the cell.
		* Corner (di, dj) of cell (i, j) is grid corner (i + di, j + dj).
		* Only the four faces in the poloidal plane are used: walls are
		* outlines in that plane, so they never lie across a z face.
		*/
		struct FaceStencil
		{
			Zones::Face face {};
			int di {};			///< Neighbor is cell (i + di, j + dj)
			int dj {};
			int c0i {};			///< First corner is (i + c0i, j + c0j)
			int c0j {};
			int c1i {};			///< Second corner is (i + c1i, j + c1j)
			int c1j {};
		};

		constexpr std::array<FaceStencil, 4> stencils {{
			{Zones::xmin, -1, 0, 0, 0, 0, 1},
			{Zones::xmax, 1, 0, 1, 0, 1, 1},
			{Zones::ymin, 0, -1, 0, 0, 1, 0},
			{Zones::ymax, 0, 1, 0, 1, 1, 1}}};
	}

	/**
	* Finds the wall surfaces of a SOLEDGE3X background.
	*
	* See the top of this file for the idea. For every cell of every zone
	* that is in the plasma (chi < 0.5), each of its four faces in the
	* poloidal plane is checked. If the cell across that face is in the
	* wall (chi >= 0.5), the face becomes a SurfaceElement:
	*
	* 1. Its two corners and middle in (R, Z) come from the background's
	*    cell corners (grid_X, grid_Y, grid_Z), and its range of phi from
	*    the zone's z grid.
	* 2. Its own area is A_face = (edge length) x (mean R of its corners)
	*    x (phi range). That's exact for a face swept around in phi with a
	*    straight edge in (R, Z).
	* 3. Its own unit normal n_face is perpendicular to the edge, pointed
	*    into the plasma, i.e., toward the plasma cell's center.
	* 4. The wall outline segment nearest to the middle of the face is
	*    found. If it's within max_wall_dist cell widths, its normal,
	*    pointed the same way as n_face, is the element's normal, and the
	*    element's area is A_face |n_face . n_wall|. A cell width is taken
	*    as twice the distance from the face to the cell's center.
	*    Otherwise (or if there's no outline) the element falls back to
	*    n_face and A_face.
	*
	* The cell across a face can be in another zone. Its chi then comes
	* from the ghost cells, which hold copies of the neighboring zone's
	* cells. Faces on a zone boundary with nothing across it (absorbing,
	* core or wall faces) are skipped: there's no cell there to be wall.
	* Each face between a plasma and a wall cell is found exactly once,
	* from its plasma side.
	*
	* @param bkg Background with a wall mask (chi), i.e., from SOLEDGE3X.
	*        Its wall outlines (Background::get_walls) give the normals.
	* @param max_wall_dist How far, in cell widths, the nearest wall
	*        segment can be from a face before the element falls back to
	*        the face's own normal and area. Defaults to 2.
	*/
	Surface::Surface(const Background::Background& bkg,
		const double max_wall_dist)
	{
		const auto& chi {bkg.get_chi()};
		if (chi.get_data().empty())
		{
			Errors::fatal_error("Sputtering surfaces need a background with "
				"a wall mask (chi), which only SOLEDGE3X backgrounds have.");
		}

		const Zones::Zones& zones {bkg.get_zones()};
		const auto& grid_X {bkg.get_grid_X()};
		const auto& grid_Y {bkg.get_grid_Y()};
		const auto& grid_Z {bkg.get_grid_Z()};
		const auto& X {bkg.get_X()};
		const auto& Y {bkg.get_Y()};
		const auto& Z {bkg.get_Z()};
		const std::vector<Segment> segments {wall_segments(bkg.get_walls())};

		for (int zone {}; zone < zones.get_nzones(); ++zone)
		{
			const int nx {zones.get_nx(zone)};
			const int ny {zones.get_ny(zone)};
			const auto& grid_z {zones.get_grid_z(zone)};

			// (R, Z) of grid corner (ci, cj) in toroidal cell k. R is the
			// distance from the axis of the corner's Cartesian position.
			auto corner = [&](int ci, int cj, int k) -> RZ
			{
				return {std::hypot(grid_X(zone, ci, cj, k),
					grid_Y(zone, ci, cj, k)), grid_Z(zone, ci, cj, k)};
			};

			for (int i {}; i < nx; ++i)
			for (int j {}; j < ny; ++j)
			for (int k {}; k < zones.get_nz(zone); ++k)
			{
				// Only plasma cells have faces on the wall (looking from the
				// wall side would find every face a second time)
				if (chi(zone, i, j, k) >= 0.5) continue;

				for (const FaceStencil& s : stencils)
				{
					// Is the cell across this face in the wall? If the
					// neighbor is past the edge of the zone, it's in
					// another zone only if this face connects to one. Its
					// chi is then in this zone's ghost cells.
					const int ni {i + s.di};
					const int nj {j + s.dj};
					const bool past_edge {ni < 0 || ni >= nx || nj < 0
						|| nj >= ny};
					if (past_edge && zones.get_face(zone, s.face) < 0)
						continue;
					if (chi(zone, ni, nj, k) < 0.5) continue;

					SurfaceElement e {};
					e.zone = zone;
					e.i = i;
					e.j = j;
					e.k = k;
					e.face = s.face;
					e.phi_min = grid_z[k];
					e.phi_max = grid_z[k + 1];

					// 1. The face's corners and middle in (R, Z)
					const RZ c0 {corner(i + s.c0i, j + s.c0j, k)};
					const RZ c1 {corner(i + s.c1i, j + s.c1j, k)};
					const RZ mid {0.5 * (c0 + c1)};
					e.R0 = c0.R;
					e.Z0 = c0.Z;
					e.R1 = c1.R;
					e.Z1 = c1.Z;
					e.R_mid = mid.R;
					e.Z_mid = mid.Z;

					// 2. The face's own area
					e.area_face = norm(c1 - c0) * mid.R
						* (e.phi_max - e.phi_min);

					// 3. The face's own normal, pointed toward the plasma
					// cell's center
					const RZ center {std::hypot(X(zone, i, j, k),
						Y(zone, i, j, k)), Z(zone, i, j, k)};
					const RZ to_plasma {center - mid};
					const RZ n_face {orient_like(unit_normal(c0, c1),
						to_plasma)};

					// 4. The nearest wall segment
					const Segment* nearest {nullptr};
					e.wall_dist = std::numeric_limits<double>::infinity();
					for (const Segment& seg : segments)
					{
						const double d {distance_to_segment(mid, seg)};
						if (d < e.wall_dist)
						{
							e.wall_dist = d;
							nearest = &seg;
						}
					}

					// Use it if it's close enough, otherwise fall back to
					// the face itself
					const double cell_width {2.0 * norm(to_plasma)};
					// n_wall is pointed the same way as n_face, i.e., into
					// the plasma. Where a face is perpendicular to the wall
					// (e.g., a short step that straddles the outline), "the
					// same way" is ambiguous and n_wall can come out
					// pointing along the wall or even into it. That hardly
					// matters, since the area, A_face |n_face . n_wall|, is
					// then small, so the element sputters little.
					if (nearest != nullptr
						&& e.wall_dist <= max_wall_dist * cell_width)
					{
						const RZ n_wall {orient_like(unit_normal(nearest->a,
							nearest->b), n_face)};
						e.n_R = n_wall.R;
						e.n_Z = n_wall.Z;
						e.area = e.area_face * std::abs(dot(n_face, n_wall));
					}
					else
					{
						e.n_R = n_face.R;
						e.n_Z = n_face.Z;
						e.area = e.area_face;
						e.fallback = true;
						++m_num_fallback;
					}

					m_total_area += e.area;
					m_elements.push_back(e);
				}
			}
		}
	}

	/**
	* @name Accessors
	* The surface elements, how many there are, how many fell back to their
	* face's own normal and area, and the sum of their areas (m^2).
	*/
	///@{
	const std::vector<SurfaceElement>& Surface::elements() const
		{return m_elements;}
	std::size_t Surface::size() const {return m_elements.size();}
	int Surface::num_fallback() const {return m_num_fallback;}
	double Surface::total_area() const {return m_total_area;}
	///@}

	/**
	* Prints how many surface elements there are and their area, in total
	* and per zone (numbered from 1), and how many fell back to their face's
	* own normal and area because no wall outline segment was close.
	*
	* @param os Stream to print to.
	*/
	void Surface::print_summary(std::ostream& os) const
	{
		struct Totals
		{
			int n {};
			int fallback {};
			double area {};
			double area_face {};
		};
		std::map<int, Totals> per_zone {};
		Totals all {};
		for (const SurfaceElement& e : m_elements)
		{
			for (Totals* t : {&per_zone[e.zone], &all})
			{
				++t->n;
				t->fallback += e.fallback;
				t->area += e.area;
				t->area_face += e.area_face;
			}
		}

		os << "Sputtering surface: " << all.n << " elements, area "
			<< all.area << " m^2 (" << all.area_face << " m^2 of cell faces), "
			<< all.fallback << " without a nearby wall outline\n";
		for (const auto& [zone, t] : per_zone)
		{
			os << "  zone " << std::setw(3) << zone + 1 << ": "
				<< std::setw(7) << t.n << " elements, area " << t.area
				<< " m^2, " << t.fallback << " without a nearby wall "
				"outline\n";
		}
	}

	/**
	* Returns the energy an ion gains falling through the sheath and
	* presheath in front of a floating wall (eV).
	*
	* The sheath potential drop for a floating wall with no secondary
	* electron emission is
	* e dphi_sh = -0.5 Te ln[2 pi (m_e / m_i) (1 + Ti / Te)], which is about
	* 2.84 Te for deuterium at Ti = Te. The presheath adds 0.5 Te. An ion of
	* charge Z gains Z times the total drop. For deuterium at Ti = Te,
	* that's about 3.34 Te, a little more than the 3 Te rule of thumb.
	*
	* @param te Electron temperature (eV), > 0.
	* @param ti Ion temperature (eV), >= 0.
	* @param elec_mass_amu Mass of the electrons (amu).
	* @param ion_mass_amu Mass of the ions (amu).
	* @param charge Charge of the ions (Z).
	* @return Energy gained (eV).
	*/
	double sheath_energy(const double te, const double ti,
		const double elec_mass_amu, const double ion_mass_amu,
		const int charge)
	{
		const double me_mi {elec_mass_amu / ion_mass_amu};
		const double sheath {-0.5 * te * std::log(2.0 * Constants::pi
			* me_mi * (1.0 + ti / te))};
		const double presheath {0.5 * te};
		return charge * (sheath + presheath);
	}

	namespace
	{
		/**
		* Nodes and weights for averaging over the thermal energy
		* distribution of the ions hitting a wall, Gamma(2, Ti).
		*
		* In x = eps / Ti the distribution is x e^-x. The midpoint rule is
		* used on [0, 30] with 300 points, with the weights x e^-x dx
		* normalized to add up to 1. The part beyond x = 30 is about
		* 1e-11. A Gauss-Laguerre rule would need far fewer points for a
		* smooth integrand, but the yield has a kink at threshold (it's 0
		* below), where it converges badly. The midpoint rule doesn't care.
		*
		* @return Pairs of (x, weight), made the first time it's called.
		*/
		const std::vector<std::array<double, 2>>& thermal_quadrature()
		{
			static const std::vector<std::array<double, 2>> nodes {[]
			{
				constexpr int n {300};
				constexpr double x_max {30.0};
				const double dx {x_max / n};
				std::vector<std::array<double, 2>> q(n);
				double total {};
				for (int i {}; i < n; ++i)
				{
					const double x {(i + 0.5) * dx};
					q[i] = {x, x * std::exp(-x) * dx};
					total += q[i][1];
				}
				for (auto& node : q) node[1] /= total;
				return q;
			}()};
			return nodes;
		}
	}

	/**
	* Averages the sputtering yield over the energies ions hit the wall
	* with.
	*
	* The energy is e_sheath plus a thermal energy from Gamma(2, ti), and
	* the average is done with thermal_quadrature.
	*
	* @param table Sputtering yield table.
	* @param e_sheath Energy gained in the sheath and presheath (eV).
	* @param ti Ion temperature (eV), > 0.
	* @param angle Angle of incidence from the surface normal (degrees).
	* @return The averaged yield, and the fraction of the ions whose energy
	*         was above the end of the table (where the yield at the
	*         table's highest energy was used).
	*/
	EffectiveYield effective_yield(const YieldTable& table,
		const double e_sheath, const double ti, const double angle)
	{
		EffectiveYield out {};
		for (const auto& [x, w] : thermal_quadrature())
		{
			const double energy {e_sheath + x * ti};
			out.yield += w * table.get_yield(energy, angle);
			if (energy > table.energy_max()) out.frac_above += w;
		}
		return out;
	}

	/**
	* Calculates the sputtered source on every surface element in every
	* background time frame.
	*
	* See the top of this file for the physics. For element e in frame t,
	* with n_i = ne (a D+ background), b_hat = B / |B| and n_wall the
	* element's normal in Cartesian coordinates at the middle of its phi
	* range, all from the element's plasma cell:
	*
	* 1. u_par = u . b_hat, where u is the background flow. The drifts in u
	*    are perpendicular to B, so this is exactly the parallel flow.
	* 2. The parallel velocity into the wall is -u_par (b_hat . n_wall). With
	*    sput_bohm_floor = "on", it's c_s |b_hat . n_wall| if that's more,
	*    where c_s = sqrt(e (Te + Ti) / m_i).
	* 3. v_E = E x B / B^2, and its velocity into the wall is -v_E . n_wall.
	* 4. flux_par and flux_exb are n_i times those, and are positive into
	*    the wall. The flux is n_i max(0, sum of the two), or
	*    max(0, flux_par) with sput_flux_opt = "parallel".
	* 5. The yield is effective_yield with E_sheath from sheath_energy (with
	*    the background's electron and ion masses) and the angle
	*    sput_angle_val. Te and Ti are held at or above min_te and
	*    min_ti. The mean impact energy saved is E_sheath + 2 Ti.
	* 6. The rate is flux x yield x area.
	*
	* The run stops if the table's sputtered atoms aren't the impurity
	* being followed (imp_atom_num), its projectile isn't singly charged
	* (the background ions are D+), or sput_angle_val isn't between 0 and
	* 90. A projectile whose mass is more than 5% off the background ion
	* mass, and element-frames where more than 1% of the ions hit with more
	* energy than the table covers, are reported by print_summary.
	*
	* @param surface The surface elements.
	* @param bkg The background plasma, with the time frames, plasma and
	*        fields.
	* @param table The sputtering yield table.
	* @param opts Simulation options: imp_atom_num, sput_flux_opt,
	*        sput_bohm_floor, sput_angle_val, min_te and min_ti.
	*/
	Source::Source(const Surface& surface,
		const Background::Background& bkg, const YieldTable& table,
		const Options::Options& opts)
		: m_surface {surface},
		m_nframes {static_cast<int>(bkg.get_times().size())}
	{
		// The table has to be for the impurity being followed, sputtered by
		// the background's singly charged ions
		if (static_cast<int>(std::lround(table.target_Z()))
			!= opts.imp_atom_num())
		{
			Errors::fatal_error("The sputtering yield table " + table.path()
				+ " is for a target with Z = "
				+ std::to_string(table.target_Z()) + ", but the impurity "
				"(imp_atom_num) is " + std::to_string(opts.imp_atom_num())
				+ ".");
		}
		if (std::lround(table.projectile_Z()) != 1)
		{
			Errors::fatal_error("The sputtering yield table " + table.path()
				+ " is for a projectile with Z = "
				+ std::to_string(table.projectile_Z()) + ", but the "
				"background ions are singly charged (D+).");
		}
		const double angle {opts.sput_angle_val()};
		if (!(angle >= 0.0 && angle <= 90.0))
		{
			Errors::fatal_error("sput_angle_val = " + std::to_string(angle)
				+ ", but it must be between 0 and 90 degrees.");
		}
		const double mi_amu {bkg.get_ion_mass_amu()};
		const double me_amu {bkg.get_elec_mass_amu()};
		if (!(mi_amu > 0.0 && me_amu > 0.0))
		{
			Errors::fatal_error("The sputtered source needs the background's "
				"electron and ion masses, but they aren't set.");
		}
		if (std::abs(table.projectile_mass_amu() / mi_amu - 1.0) > 0.05)
		{
			std::ostringstream msg {};
			msg << "The sputtering yield table's projectile mass ("
				<< table.projectile_mass_amu() << " amu) is more than 5% off "
				"the background ion mass (" << mi_amu << " amu).";
			m_warnings.push_back(msg.str());
		}

		const bool use_exb {opts.sput_flux_opt_int() == 0};
		const bool bohm_floor {opts.sput_bohm_floor_int() == 1};
		const double mi_kg {mi_amu * Constants::amu_to_kg};

		const auto& ne {bkg.get_ne()};
		const auto& te {bkg.get_te()};
		const auto& ti {bkg.get_ti()};
		const auto& uX {bkg.get_uX()};
		const auto& uY {bkg.get_uY()};
		const auto& uZ {bkg.get_uZ()};
		const auto& bX {bkg.get_bX()};
		const auto& bY {bkg.get_bY()};
		const auto& bZ {bkg.get_bZ()};
		const auto& eX {bkg.get_eX()};
		const auto& eY {bkg.get_eY()};
		const auto& eZ {bkg.get_eZ()};

		const std::size_t n_elem {m_surface.size()};
		const std::size_t n_all {n_elem * m_nframes};
		for (auto* v : {&m_flux, &m_flux_par, &m_flux_exb, &m_energy,
			&m_yield, &m_rate, &m_e_sheath, &m_ti})
		{
			v->assign(n_all, 0.0);
		}

		// Kept for starting particles
		m_table = table;
		m_angle = angle;
		const double m1 {table.projectile_mass_amu()};
		const double m2 {table.target_mass_amu()};
		m_lambda = 4.0 * m1 * m2 / ((m1 + m2) * (m1 + m2));
		m_imp_mass_kg = opts.imp_mass_amu() * Constants::amu_to_kg;
		m_total_rate.assign(m_nframes, 0.0);

		using Vec3 = std::array<double, 3>;
		auto dot3 = [](const Vec3& a, const Vec3& b)
			{return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];};

		for (int t {}; t < m_nframes; ++t)
		for (std::size_t e {}; e < n_elem; ++e)
		{
			const SurfaceElement& el {m_surface.elements()[e]};
			const int z {el.zone};
			const int i {el.i};
			const int j {el.j};
			const int k {el.k};
			const std::size_t idx {t * n_elem + e};

			// The element's normal in Cartesian coordinates, in the middle
			// of its phi range. Positions are (R cos phi,
			// phi_sign R sin phi, Z), so R_hat is (cos phi,
			// phi_sign sin phi, 0).
			const double phi {0.5 * (el.phi_min + el.phi_max)};
			const Vec3 n_wall {el.n_R * std::cos(phi),
				el.n_R * S3X::phi_sign * std::sin(phi), el.n_Z};

			// Plasma in the element's plasma cell
			const double n_i {std::max(static_cast<double>(
				ne(z, t, i, j, k)), 0.0)};
			const double te_c {std::max(static_cast<double>(
				te(z, t, i, j, k)), opts.min_te())};
			const double ti_c {std::max(static_cast<double>(
				ti(z, t, i, j, k)), opts.min_ti())};
			const Vec3 B {bX(z, t, i, j, k), bY(z, t, i, j, k),
				bZ(z, t, i, j, k)};
			const double B2 {dot3(B, B)};

			// 1-4. Flux onto the wall. Without a field there's no
			// direction for the parallel flow or the drift.
			double flux_par {};
			double flux_exb {};
			if (B2 > 0.0)
			{
				const double Bmag {std::sqrt(B2)};
				const Vec3 b_hat {B[0] / Bmag, B[1] / Bmag, B[2] / Bmag};
				const double b_n {dot3(b_hat, n_wall)};
				const Vec3 u {uX(z, t, i, j, k), uY(z, t, i, j, k),
					uZ(z, t, i, j, k)};

				double v_par_in {-dot3(u, b_hat) * b_n};
				if (bohm_floor)
				{
					const double c_s {std::sqrt(Constants::ev_to_j
						* (te_c + ti_c) / mi_kg)};
					v_par_in = std::max(v_par_in, c_s * std::abs(b_n));
				}

				const Vec3 E {eX(z, t, i, j, k), eY(z, t, i, j, k),
					eZ(z, t, i, j, k)};
				const Vec3 v_E {(E[1] * B[2] - E[2] * B[1]) / B2,
					(E[2] * B[0] - E[0] * B[2]) / B2,
					(E[0] * B[1] - E[1] * B[0]) / B2};

				flux_par = n_i * v_par_in;
				flux_exb = -n_i * dot3(v_E, n_wall);
			}
			const double flux {std::max(0.0,
				use_exb ? flux_par + flux_exb : flux_par)};

			// 5. Effective yield at the fixed angle
			const double e_sheath {sheath_energy(te_c, ti_c, me_amu, mi_amu,
				1)};
			const EffectiveYield y {effective_yield(table, e_sheath, ti_c,
				angle)};
			if (y.frac_above > 0.01) ++m_num_above_table;

			// 6. Rate
			m_flux[idx] = flux;
			m_flux_par[idx] = flux_par;
			m_flux_exb[idx] = flux_exb;
			m_energy[idx] = e_sheath + 2.0 * ti_c;
			m_yield[idx] = y.yield;
			m_rate[idx] = flux * y.yield * el.area;
			m_e_sheath[idx] = e_sheath;
			m_ti[idx] = ti_c;
			m_total_rate[t] += m_rate[idx];
		}

		// Computational x and y edges of each element's plasma cell, for
		// placing atoms on its face
		const Zones::Zones& zones {bkg.get_zones()};
		m_cell_bounds.reserve(n_elem);
		for (const SurfaceElement& el : m_surface.elements())
		{
			const auto& gx {zones.get_grid_x(el.zone)};
			const auto& gy {zones.get_grid_y(el.zone)};
			m_cell_bounds.push_back({gx[el.i], gx[el.i + 1], gy[el.j],
				gy[el.j + 1]});
		}

		// Each frame stands for the times nearer to it than to the other
		// frames (the frame a particle uses is the nearest one), from half
		// way to the frame before to half way to the frame after, and to
		// the ends of the background's time range at the first and last.
		const auto& times {bkg.get_times()};
		m_frame_lo.resize(m_nframes);
		m_frame_hi.resize(m_nframes);
		for (int t {}; t < m_nframes; ++t)
		{
			m_frame_lo[t] = t == 0 ? times[0]
				: 0.5 * (times[t - 1] + times[t]);
			m_frame_hi[t] = t == m_nframes - 1 ? times[t]
				: 0.5 * (times[t] + times[t + 1]);
		}

		// With one frame its window has no length, so give it a length of
		// 1 to weight it
		auto window = [&](int t)
		{
			return m_nframes == 1 ? 1.0 : m_frame_hi[t] - m_frame_lo[t];
		};

		// Atoms sputtered off each element in each frame's window, added
		// up in the order the per-element values are stored, for picking
		// births. The time-averaged rate is the total over the windows
		// divided by their total length.
		m_cdf.resize(n_all);
		double cumulative {};
		double total_time {};
		for (int t {}; t < m_nframes; ++t)
		{
			for (std::size_t e {}; e < n_elem; ++e)
			{
				cumulative += m_rate[t * n_elem + e] * window(t);
				m_cdf[t * n_elem + e] = cumulative;
			}
			total_time += window(t);
		}
		m_source_rate = total_time > 0.0 ? cumulative / total_time : 0.0;
	}

	/**
	* @name Accessors
	* The surface elements, how many there are, the number of time frames,
	* and the per-element values in every frame (index t * size() + e):
	* the ion flux onto the wall and its parallel and E x B parts (m^-2
	* s^-1, positive into the wall), the mean impact energy (eV), the
	* effective yield and the sputtered rate (atoms/s). Also the total rate
	* in each frame, its mean over the frames, and how many element-frames
	* had more than 1% of their ions above the end of the yield table.
	* source_rate is the total rate averaged over the frames' windows
	* (see sample_birth), and frame_lo and frame_hi are the windows.
	*/
	///@{
	const Surface& Source::surface() const {return m_surface;}
	std::size_t Source::size() const {return m_surface.size();}
	int Source::nframes() const {return m_nframes;}
	const std::vector<double>& Source::flux() const {return m_flux;}
	const std::vector<double>& Source::flux_par() const {return m_flux_par;}
	const std::vector<double>& Source::flux_exb() const {return m_flux_exb;}
	const std::vector<double>& Source::energy() const {return m_energy;}
	const std::vector<double>& Source::yield() const {return m_yield;}
	const std::vector<double>& Source::rate() const {return m_rate;}
	const std::vector<double>& Source::total_rate() const
		{return m_total_rate;}
	int Source::num_above_table() const {return m_num_above_table;}
	double Source::source_rate() const {return m_source_rate;}
	const std::vector<double>& Source::frame_lo() const {return m_frame_lo;}
	const std::vector<double>& Source::frame_hi() const {return m_frame_hi;}

	double Source::mean_total_rate() const
	{
		if (m_total_rate.empty()) return 0.0;
		double sum {};
		for (const double r : m_total_rate) sum += r;
		return sum / m_total_rate.size();
	}
	///@}

	/**
	* Prints the sputtered source: the total rate (mean, and the range over
	* the frames), how much of the ion flux onto the wall came from the
	* parallel flow and from the E x B drift, the flux-weighted mean impact
	* energy and yield, and any warnings.
	*
	* @param os Stream to print to.
	*/
	void Source::print_summary(std::ostream& os) const
	{
		double sum_flux_area {};
		double sum_par {};
		double sum_exb {};
		double sum_energy {};
		double sum_yield {};
		const std::size_t n {size()};
		for (int t {}; t < m_nframes; ++t)
		for (std::size_t e {}; e < n; ++e)
		{
			const std::size_t idx {t * n + e};
			const double a {m_surface.elements()[e].area};
			sum_flux_area += m_flux[idx] * a;
			sum_par += m_flux_par[idx] * a;
			sum_exb += m_flux_exb[idx] * a;
			sum_energy += m_energy[idx] * m_flux[idx] * a;
			sum_yield += m_yield[idx] * m_flux[idx] * a;
		}
		const auto [min_it, max_it] {std::minmax_element(m_total_rate.begin(),
			m_total_rate.end())};

		os << "Sputtered source: " << mean_total_rate() << " atoms/s ";
		if (!m_total_rate.empty())
		{
			os << "(" << *min_it << " to " << *max_it << " over "
				<< m_nframes << " frames)\n";
		}
		if (sum_flux_area > 0.0)
		{
			os << "  Ions onto the wall: " << sum_flux_area / m_nframes
				<< " ions/s, net parallel " << sum_par / m_nframes
				<< " and net E x B " << sum_exb / m_nframes << "\n"
				<< "  Mean impact energy " << sum_energy / sum_flux_area
				<< " eV and yield " << sum_yield / sum_flux_area
				<< " (weighted by the ion flux)\n";
		}
		if (m_num_above_table > 0)
		{
			os << "  Warning! " << m_num_above_table << " element-frames "
				"had more than 1% of their ions above the highest energy in "
				"the yield table, where the yield at that energy was used.\n";
		}
		for (const auto& w : m_warnings) os << "  Warning! " << w << "\n";
	}

	/**
	* Picks the energy of an ion that sputtered an atom.
	*
	* The ions hitting the wall have energies E_sheath + eps, with eps from
	* Gamma(2, ti). The ones that sputter an atom are weighted by the
	* yield, so their energies go as Y(E) p(eps). Near threshold that's
	* only the fast tail. The energy is picked from the same points as
	* effective_yield, weighted by w_i Y(E_i), so it's one of those points.
	*
	* @param table Sputtering yield table.
	* @param e_sheath Energy gained in the sheath and presheath (eV).
	* @param ti Ion temperature (eV), > 0.
	* @param angle Angle of incidence from the surface normal (degrees).
	* @param u Uniform random number in [0, 1).
	* @return Energy of the ion (eV), or the mean energy e_sheath + 2 ti if
	*         no energy sputters anything.
	*/
	double sample_impact_energy(const YieldTable& table, const double e_sheath,
		const double ti, const double angle, const double u)
	{
		const auto& nodes {thermal_quadrature()};
		std::vector<double> cdf(nodes.size());
		double total {};
		for (std::size_t n {}; n < nodes.size(); ++n)
		{
			const auto& [x, w] {nodes[n]};
			total += w * table.get_yield(e_sheath + x * ti, angle);
			cdf[n] = total;
		}
		if (!(total > 0.0)) return e_sheath + 2.0 * ti;

		const std::size_t n {static_cast<std::size_t>(std::upper_bound(
			cdf.begin(), cdf.end(), u * total) - cdf.begin())};
		const std::size_t pick {std::min(n, nodes.size() - 1)};
		return e_sheath + nodes[pick][0] * ti;
	}

	/**
	* Picks the energy of a sputtered atom from the Thompson distribution,
	* cut off at e_max.
	*
	* The Thompson distribution is f(E) ~ E / (E + e_b)^3, with e_b the
	* surface binding energy. Its cumulative distribution is
	* F(E) = E^2 / (E + e_b)^2, so with u' = u F(e_max), solving F(E) = u'
	* gives E = e_b sqrt(u') / (1 - sqrt(u')), which is exactly distributed
	* as f(E) between 0 and e_max.
	*
	* @param e_b Surface binding energy (eV), > 0.
	* @param e_max Largest energy (eV). If it isn't > 0, e_b is returned.
	* @param u Uniform random number in [0, 1).
	* @return Energy of the sputtered atom (eV).
	*/
	double sample_thompson(const double e_b, const double e_max,
		const double u)
	{
		if (!(e_max > 0.0)) return e_b;
		const double f_max {e_max * e_max / ((e_max + e_b) * (e_max + e_b))};
		const double s {std::sqrt(u * f_max)};
		return e_b * s / (1.0 - s);
	}

	/**
	* Picks a direction from a cosine distribution about a unit normal.
	*
	* The number of atoms leaving at angle theta from the normal goes as
	* cos(theta) sin(theta), so cos^2(theta) is uniform: cos(theta) =
	* sqrt(u1). The azimuth is 2 pi u2. The direction is put together from
	* the normal and two unit vectors perpendicular to it.
	*
	* @param n Unit normal, which the direction is on the same side of.
	* @param u1 Uniform random number in [0, 1), for the polar angle.
	* @param u2 Uniform random number in [0, 1), for the azimuth.
	* @return Unit direction.
	*/
	std::array<double, 3> cosine_direction(const std::array<double, 3>& n,
		const double u1, const double u2)
	{
		using Vec3 = std::array<double, 3>;
		auto cross = [](const Vec3& a, const Vec3& b) -> Vec3
		{
			return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2],
				a[0] * b[1] - a[1] * b[0]};
		};

		// A perpendicular to n: cross it with the axis it's least along
		Vec3 axis {0.0, 0.0, 0.0};
		const int least {std::abs(n[0]) <= std::abs(n[1])
			? (std::abs(n[0]) <= std::abs(n[2]) ? 0 : 2)
			: (std::abs(n[1]) <= std::abs(n[2]) ? 1 : 2)};
		axis[least] = 1.0;
		Vec3 t1 {cross(n, axis)};
		const double len {std::sqrt(t1[0] * t1[0] + t1[1] * t1[1]
			+ t1[2] * t1[2])};
		for (double& c : t1) c /= len;
		const Vec3 t2 {cross(n, t1)};

		const double cos_th {std::sqrt(u1)};
		const double sin_th {std::sqrt(1.0 - u1)};
		const double az {2.0 * Constants::pi * u2};
		Vec3 d {};
		for (int c {}; c < 3; ++c)
		{
			d[c] = cos_th * n[c] + sin_th * (std::cos(az) * t1[c]
				+ std::sin(az) * t2[c]);
		}
		return d;
	}

	/**
	* Picks where, when and how fast a new sputtered atom starts.
	*
	* See the top of this file for the steps. The atom starts on its
	* element's face, moved a thousandth of a cell into the plasma cell so
	* it isn't on the edge of the wall cell, at a uniformly random place
	* along the face and in the element's range of phi. Its direction is
	* from the element's normal at its own phi.
	*
	* @param rng Random number generator of the calling thread.
	* @return The new atom. Its charge (0) and weight are set by the caller.
	*/
	Birth Source::sample_birth(pcg32& rng) const
	{
		if (m_cdf.empty() || !(m_cdf.back() > 0.0))
		{
			Errors::fatal_error("Can't start sputtered atoms, since nothing "
				"is sputtered (the source is 0 everywhere).");
		}

		// 1. Frame and element, in proportion to the atoms they sputter
		const double r {rng.next_double() * m_cdf.back()};
		const std::size_t idx {std::min(static_cast<std::size_t>(
			std::upper_bound(m_cdf.begin(), m_cdf.end(), r) - m_cdf.begin()),
			m_cdf.size() - 1)};
		const std::size_t n_elem {size()};
		Birth b {};
		b.frame = static_cast<int>(idx / n_elem);
		b.element = idx % n_elem;
		const SurfaceElement& el {m_surface.elements()[b.element]};
		b.zone = el.zone;

		// Time, uniform in the frame's window
		b.t = m_frame_lo[b.frame] + rng.next_double()
			* (m_frame_hi[b.frame] - m_frame_lo[b.frame]);

		// 2. Place: uniform along the face, just inside the plasma cell
		constexpr double inset {1e-3};
		const auto& [x_lo, x_hi, y_lo, y_hi] {m_cell_bounds[b.element]};
		const double ux {rng.next_double()};
		switch (el.face)
		{
			case Zones::xmin:
				b.x = x_lo + inset * (x_hi - x_lo);
				b.y = y_lo + ux * (y_hi - y_lo);
				break;
			case Zones::xmax:
				b.x = x_hi - inset * (x_hi - x_lo);
				b.y = y_lo + ux * (y_hi - y_lo);
				break;
			case Zones::ymin:
				b.x = x_lo + ux * (x_hi - x_lo);
				b.y = y_lo + inset * (y_hi - y_lo);
				break;
			default:	// ymax
				b.x = x_lo + ux * (x_hi - x_lo);
				b.y = y_hi - inset * (y_hi - y_lo);
				break;
		}
		b.z = el.phi_min + rng.next_double() * (el.phi_max - el.phi_min);

		// 3. Energy of the ion that sputtered it
		b.e_impact = sample_impact_energy(m_table, m_e_sheath[idx],
			m_ti[idx], m_angle, rng.next_double());

		// 4. Energy of the atom, Thompson cut off at Lambda E_ion
		b.e_atom = sample_thompson(m_table.target_Es(),
			m_lambda * b.e_impact, rng.next_double());
		const double speed {std::sqrt(2.0 * b.e_atom * Constants::ev_to_j
			/ m_imp_mass_kg)};

		// 5. Direction: cosine about the normal at the atom's phi
		const std::array<double, 3> n {el.n_R * std::cos(b.z),
			el.n_R * S3X::phi_sign * std::sin(b.z), el.n_Z};
		const double u1 {rng.next_double()};
		const double u2 {rng.next_double()};
		const auto d {cosine_direction(n, u1, u2)};
		b.vX = speed * d[0];
		b.vY = speed * d[1];
		b.vZ = speed * d[2];
		return b;
	}
}
