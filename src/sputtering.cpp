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
*/

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <limits>
#include <map>
#include <ostream>
#include <string>
#include <vector>

#include "background.h"
#include "errors.h"
#include "hdf5_utils.h"
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
}
