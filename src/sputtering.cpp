/**
* @file sputtering.cpp
*
* @brief Sputtering of wall surfaces by the background plasma.
*
* The parts: the sputtering yield, the wall surfaces that get sputtered,
* the sputtered source and its atoms, and where they land again.
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
* Grids fitted to the targets (HERMES) are simpler. A hypnotoad grid ends
* at the divertor targets, so the targets are the zone faces marked as wall
* faces (Zones::face_wall), there's no wall mask, and the plasma cells
* along a target are the ones whose plasma hits it. Every cell face on a
* wall face is an element with its own normal and area, which are the
* target's, so nothing falls back and no outline is needed. The ion flux
* onto it is then the flow along the field through the target face,
* n u_par b . n, the same as for a staircase face with the true normal.
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
*
* Redeposition (Deposition)
* -------------------------
* An atom that moves into a wall cell (chi = 1) is absorbed by the boundary
* check, which then puts it on the surface element between that wall cell
* and the plasma cell it came from (each element knows its wall cell). An
* atom's weight is the atoms/s it stands for, so adding up the weights
* landing on an element gives the rate atoms are redeposited there, and
* subtracting it from the element's source gives its net erosion. How
* long after being sputtered each atom landed, and its charge then, are
* tallied too. On a grid fitted to the targets, an atom absorbed at a wall
* face is put on the element at the cell of the face it crossed. On a
* SOLEDGE3X background the wall faces (the ends of the zones) are inside
* the material, with no elements on them, so atoms absorbed there are only
* counted in total.
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
#include <type_traits>
#include <vector>

#include "mpi.h"

#include "background.h"
#include "constants.h"
#include "errors.h"
#include "hdf5_utils.h"
#include "indices.h"
#include "options.h"
#include "sheath.h"
#include "sputtering.h"
#include "zones.h"

#ifdef USE_CUDA
#include <cuda_runtime.h>
#endif


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
		return Sputtering::get_yield(view(), energy, angle);
	}

	/**
	* The table as flat arrays, for the routines shared with the GPU
	* (sputtering_hd.h). The pointers are into this table's own vectors.
	*
	* @return The view.
	*/
	YieldView YieldTable::view() const
	{
		YieldView v {};
		v.n_energy = static_cast<int>(m_log_energy.size());
		v.n_angle = static_cast<int>(m_angle.size());
		v.log_energy = m_log_energy.data();
		v.angle = m_angle.data();
		v.yield = m_yield.data();
		v.energy_min = m_energy_min;
		v.target_Es = m_target_Es;
		return v;
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
	* Finds the wall surfaces of a background.
	*
	* A background with a wall mask (SOLEDGE3X) has its mesh run past the
	* wall, so the surface is the staircase of faces between plasma and
	* wall cells (find_staircase). Its wall faces (the ends of the zones)
	* are inside the material and aren't used. A background without one
	* whose grid is fitted to the targets (HERMES) has the targets at the
	* wall faces of its zones (find_wall_faces).
	*
	* A slab background with a wall (slab_wall) is Cartesian, and its wall
	* face is found the same way (find_wall_faces).
	*
	* The run stops if neither gives any surface, e.g., for a Gkeyll
	* background, or a slab without a wall.
	*
	* @param bkg The background, with its geometry.
	* @param max_wall_dist How far, in cell widths, the nearest wall
	*        outline segment can be from a staircase face before the
	*        element falls back to the face's own normal and area. Only
	*        used with a wall mask. Defaults to 2.
	*/
	Surface::Surface(const Background::Background& bkg,
		const double max_wall_dist)
	{
		if (!bkg.get_chi().get_data().empty())
			find_staircase(bkg, max_wall_dist);
		else
			find_wall_faces(bkg);

		if (m_elements.empty())
		{
			Errors::fatal_error("Sputtering surfaces need a background with "
				"a wall mask (chi), as SOLEDGE3X backgrounds have, with "
				"wall faces on a toroidal grid fitted to the targets, as "
				"HERMES backgrounds have, or a slab background with a wall "
				"(slab_wall). This background has none of these.");
		}
		build_flat();
	}

	/**
	* Copies the elements into flat arrays, one per field, for
	* ElementsView (the routines shared with the GPU take plain arrays).
	*/
	void Surface::build_flat()
	{
		m_flat = {};
		for (const SurfaceElement& el : m_elements)
		{
			m_flat.zone.push_back(el.zone);
			m_flat.i.push_back(el.i);
			m_flat.j.push_back(el.j);
			m_flat.k.push_back(el.k);
			m_flat.face.push_back(static_cast<int>(el.face));
			m_flat.n_R.push_back(el.n_R);
			m_flat.n_Z.push_back(el.n_Z);
			m_flat.phi_min.push_back(el.phi_min);
			m_flat.phi_max.push_back(el.phi_max);
			m_flat.cartesian.push_back(el.cartesian ? 1 : 0);
		}
	}

	/**
	* The elements as flat arrays (see build_flat), pointing into this
	* surface's own vectors.
	*
	* @return The view.
	*/
	ElementsView Surface::view() const
	{
		ElementsView v {};
		v.n = m_elements.size();
		v.zone = m_flat.zone.data();
		v.i = m_flat.i.data();
		v.j = m_flat.j.data();
		v.k = m_flat.k.data();
		v.face = m_flat.face.data();
		v.n_R = m_flat.n_R.data();
		v.n_Z = m_flat.n_Z.data();
		v.phi_min = m_flat.phi_min.data();
		v.phi_max = m_flat.phi_max.data();
		v.cartesian = m_flat.cartesian.data();
		return v;
	}

	/**
	* Finds the wall surfaces of a grid fitted to the targets (HERMES).
	*
	* Hypnotoad grids end at the divertor targets, so the targets are the
	* zone faces with the code Zones::face_wall, and every plasma cell on
	* one has a face that is part of the target. Each becomes a
	* SurfaceElement with that face's own corners, area and normal (pointed
	* into the plasma), exactly as find_staircase works them out before
	* looking for a wall outline: the face is the wall here, so there's
	* nothing to correct and no element falls back. There's no wall cell
	* across the face (wall_zone = -1), so redeposited atoms are put on the
	* element by where they cross the face (Deposition::add_at_wall_face).
	*
	* Only zones whose z is the toroidal angle (Zones::is_toroidal) or that
	* are Cartesian (Zones::is_cartesian, the straight slab) are used,
	* since the area and the normal's direction need to know what z is.
	* Wall faces in z are skipped, as they aren't in the poloidal plane (or
	* the X-Y plane).
	*
	* On a Cartesian zone the face's corners are its cell edges in (X, Y)
	* directly, held in the R and Z fields of the element, and its area is
	* the edge length times the cell's length in z (cartesian = true).
	*
	* @param bkg Background whose zones have wall faces.
	*/
	void Surface::find_wall_faces(const Background::Background& bkg)
	{
		m_wall_faces = true;
		const Zones::Zones& zones {bkg.get_zones()};
		const auto& grid_X {bkg.get_grid_X()};
		const auto& grid_Y {bkg.get_grid_Y()};
		const auto& grid_Z {bkg.get_grid_Z()};
		const auto& X {bkg.get_X()};
		const auto& Y {bkg.get_Y()};
		const auto& Z {bkg.get_Z()};

		for (int zone {}; zone < zones.get_nzones(); ++zone)
		{
			const bool cartesian {zones.is_cartesian(zone)};
			if (!zones.is_toroidal(zone) && !cartesian) continue;
			const int nx {zones.get_nx(zone)};
			const int ny {zones.get_ny(zone)};
			const auto& grid_x {zones.get_grid_x(zone)};
			const auto& grid_y {zones.get_grid_y(zone)};
			const auto& grid_z {zones.get_grid_z(zone)};

			// A cell corner in the poloidal plane, (R, Z), or in the X-Y
			// plane on a Cartesian zone, where x and y are X and Y
			auto corner = [&](int ci, int cj, int k) -> RZ
			{
				if (cartesian) return {grid_x[ci], grid_y[cj]};
				return {std::hypot(grid_X(zone, ci, cj, k),
					grid_Y(zone, ci, cj, k)), grid_Z(zone, ci, cj, k)};
			};
			// A cell's center in the same plane
			auto center_of = [&](int i, int j, int k) -> RZ
			{
				if (cartesian)
				{
					return {0.5 * (grid_x[i] + grid_x[i + 1]),
						0.5 * (grid_y[j] + grid_y[j + 1])};
				}
				return {std::hypot(X(zone, i, j, k), Y(zone, i, j, k)),
					Z(zone, i, j, k)};
			};

			for (const FaceStencil& s : stencils)
			{
				if (zones.get_face(zone, s.face) != Zones::face_wall)
					continue;

				// The cells along the face: one row or column of the zone
				const bool x_face {s.face == Zones::xmin
					|| s.face == Zones::xmax};
				const int i_lo {s.face == Zones::xmax ? nx - 1 : 0};
				const int i_hi {s.face == Zones::xmin ? 1 : nx};
				const int j_lo {s.face == Zones::ymax ? ny - 1 : 0};
				const int j_hi {s.face == Zones::ymin ? 1 : ny};
				for (int i {x_face ? i_lo : 0}; i < (x_face ? i_hi : nx); ++i)
				for (int j {x_face ? 0 : j_lo}; j < (x_face ? ny : j_hi); ++j)
				for (int k {}; k < zones.get_nz(zone); ++k)
				{
					SurfaceElement e {};
					e.zone = zone;
					e.i = i;
					e.j = j;
					e.k = k;
					e.face = s.face;
					e.wall_zone = -1;
					e.wall_i = -1;
					e.wall_j = -1;
					e.phi_min = grid_z[k];
					e.phi_max = grid_z[k + 1];

					const RZ c0 {corner(i + s.c0i, j + s.c0j, k)};
					const RZ c1 {corner(i + s.c1i, j + s.c1j, k)};
					const RZ mid {0.5 * (c0 + c1)};
					e.R0 = c0.R;
					e.Z0 = c0.Z;
					e.R1 = c1.R;
					e.Z1 = c1.Z;
					e.R_mid = mid.R;
					e.Z_mid = mid.Z;
					// Swept around in phi (R dphi), or straight along Z
					e.cartesian = cartesian;
					e.area_face = norm(c1 - c0) * (e.phi_max - e.phi_min)
						* (cartesian ? 1.0 : mid.R);
					e.area = e.area_face;

					const RZ center {center_of(i, j, k)};
					const RZ n_face {orient_like(unit_normal(c0, c1),
						center - mid)};
					e.n_R = n_face.R;
					e.n_Z = n_face.Z;

					m_total_area += e.area;
					m_elements.push_back(e);
				}
			}
		}
	}

	/**
	* Finds the staircase of wall surfaces of a background with a wall
	* mask (SOLEDGE3X).
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
	*        the face's own normal and area.
	*/
	void Surface::find_staircase(const Background::Background& bkg,
		const double max_wall_dist)
	{
		const auto& chi {bkg.get_chi()};

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

					// The wall cell, in its own zone if it's across a
					// zone face (resolve_cell finds the cell a ghost index
					// is a copy of)
					e.wall_zone = zone;
					e.wall_i = ni;
					e.wall_j = nj;
					if (past_edge)
					{
						const Zones::Cell wall {zones.resolve_cell(
							Zones::Cell {zone, ni, nj, k}).front()};
						e.wall_zone = wall.zone;
						e.wall_i = wall.i;
						e.wall_j = wall.j;
					}
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
	* face's own normal and area, the sum of their areas (m^2), and whether
	* they are the zones' wall faces (targets) rather than a staircase.
	*/
	///@{
	const std::vector<SurfaceElement>& Surface::elements() const
		{return m_elements;}
	std::size_t Surface::size() const {return m_elements.size();}
	int Surface::num_fallback() const {return m_num_fallback;}
	double Surface::total_area() const {return m_total_area;}
	bool Surface::from_wall_faces() const {return m_wall_faces;}
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

		// Target faces are the wall itself, so there's no outline to fall
		// back from
		if (m_wall_faces)
		{
			os << "Sputtering surface (target faces): " << all.n
				<< " elements, area " << all.area << " m^2\n";
			for (const auto& [zone, t] : per_zone)
			{
				os << "  zone " << std::setw(3) << zone + 1 << ": "
					<< std::setw(7) << t.n << " elements, area " << t.area
					<< " m^2\n";
			}
			return;
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
	* Whether a background has wall surfaces a Surface can be made from:
	* a wall mask (chi), or a wall face on the x or y side of a toroidal or
	* Cartesian zone. The Surface constructor stops the run when there are
	* none, so this lets options that are on by default (sheath_model)
	* skip backgrounds without walls instead.
	*
	* @param bkg The background.
	* @return True if Surface(bkg) would find elements to make.
	*/
	bool has_surfaces(const Background::Background& bkg)
	{
		if (!bkg.get_chi().get_data().empty()) return true;
		const Zones::Zones& zones {bkg.get_zones()};
		for (int zone {}; zone < zones.get_nzones(); ++zone)
		{
			if (!zones.is_toroidal(zone) && !zones.is_cartesian(zone))
				continue;
			for (const auto f : {Zones::xmin, Zones::xmax, Zones::ymin,
				Zones::ymax})
			{
				if (zones.get_face(zone, f) == Zones::face_wall) return true;
			}
		}
		return false;
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
		* @return The nodes x and weights, made the first time it's called.
		*/
		struct QuadratureNodes
		{
			std::vector<double> x {};	///< Energy over Ti
			std::vector<double> w {};	///< Weight
		};

		const QuadratureNodes& thermal_quadrature()
		{
			static const QuadratureNodes nodes {[]
			{
				constexpr int n {300};
				constexpr double x_max {30.0};
				const double dx {x_max / n};
				QuadratureNodes q {};
				double total {};
				for (int i {}; i < n; ++i)
				{
					const double x {(i + 0.5) * dx};
					q.x.push_back(x);
					q.w.push_back(x * std::exp(-x) * dx);
					total += q.w.back();
				}
				for (double& w : q.w) w /= total;
				return q;
			}()};
			return nodes;
		}

		/// The thermal quadrature as a Quadrature view
		Quadrature quadrature_view()
		{
			const QuadratureNodes& q {thermal_quadrature()};
			return {static_cast<int>(q.x.size()), q.x.data(), q.w.data()};
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
		const QuadratureNodes& q {thermal_quadrature()};
		for (std::size_t n {}; n < q.x.size(); ++n)
		{
			const double x {q.x[n]};
			const double w {q.w[n]};
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
		const Options::Options& opts, const Sheath::Model impact_model)
		: m_surface {surface},
		m_nframes {static_cast<int>(bkg.get_times().size())},
		m_phi_sign {bkg.get_phi_sign()}
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
			&m_yield, &m_rate, &m_e_sheath, &m_ti, &m_phi_wall, &m_te,
			&m_impact_angle})
		{
			v->assign(n_all, 0.0);
		}
		m_point.assign(n_all, -1);

		// D+ impacts from the sheath model: each element and frame gets a
		// point of the table, sampled once they're all known
		if (impact_model != Sheath::Model::off)
		{
			m_impacts = std::make_shared<Sheath::ImpactTable>(impact_model,
				me_amu, mi_amu, opts.sheath_fd(), opts.sheath_drop());
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
			// phi_sign R sin phi, Z), with the background's phi_sign, so
			// R_hat is (cos phi, phi_sign sin phi, 0).
			const Vec3 n_wall {el.normal_xyz(
				0.5 * (el.phi_min + el.phi_max), m_phi_sign)};

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
			double cos_theta {};	// |b . n|, for the sheath
			if (B2 > 0.0)
			{
				const double Bmag {std::sqrt(B2)};
				const Vec3 b_hat {B[0] / Bmag, B[1] / Bmag, B[2] / Bmag};
				const double b_n {dot3(b_hat, n_wall)};
				cos_theta = std::abs(b_n);
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

			// 5. Effective yield at the fixed angle, or (with the sheath
			// model) a point of the impact table, whose yield is found
			// once the table is sampled, below
			const double e_sheath {sheath_energy(te_c, ti_c, me_amu, mi_amu,
				1)};
			if (m_impacts && B2 > 0.0)
			{
				m_point[idx] = m_impacts->add(Sheath::make_params(te_c, ti_c,
					n_i, std::sqrt(B2), cos_theta, me_amu, mi_amu,
					opts.sheath_fd(), opts.sheath_drop()));
			}
			m_flux[idx] = flux;
			m_flux_par[idx] = flux_par;
			m_flux_exb[idx] = flux_exb;
			m_e_sheath[idx] = e_sheath;
			m_ti[idx] = ti_c;
			m_te[idx] = te_c;
			if (m_point[idx] < 0)
			{
				const EffectiveYield y {effective_yield(table, e_sheath,
					ti_c, angle)};
				if (y.frac_above > 0.01) ++m_num_above_table;
				m_energy[idx] = e_sheath + 2.0 * ti_c;
				m_yield[idx] = y.yield;
				m_impact_angle[idx] = angle;
			}

			// 6. Rate
			m_rate[idx] = flux * m_yield[idx] * el.area;

			// The floating drop alone (no presheath), as the sheath model
			// has it, or sheath_drop Te. For self-sputtering without the
			// sheath field (SelfSource).
			m_phi_wall[idx] = opts.sheath_drop() > 0.0
				? opts.sheath_drop() * te_c : e_sheath - 0.5 * te_c;
			m_total_rate[t] += m_rate[idx];
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

		// With the sheath model, the yield averaged over the impact
		// samples of each element's point, weighted by their share of the
		// flux, with energies over Te times the element's Te. The energies
		// and angles of the ions hitting the wall, and of those weighted by
		// their yield (the atoms they sputter), are histogrammed at the
		// same time, averaged over the frames' windows like the rates
		// (ions/s and atoms/s in each bin). Each element and frame is
		// independent, so they're split between the threads, each with
		// histograms of its own that are added up at the end.
		if (m_impacts)
		{
			m_impacts->build();

			// The samples flattened, for SourceView
			m_sample_offset.assign(1, 0);
			for (std::size_t p {}; p < m_impacts->size(); ++p)
			{
				for (const auto& smp : m_impacts->samples(
					static_cast<int>(p)))
				{
					m_sample_energy.push_back(smp.energy);
					m_sample_angle.push_back(smp.angle);
					m_sample_weight.push_back(smp.weight);
				}
				m_sample_offset.push_back(static_cast<int>(
					m_sample_energy.size()));
			}

			// Energy bins of 20 a decade from 0.01 eV to 100 keV
			constexpr int per_decade {20};
			constexpr double log_lo {-2.0};
			const int n_ebin {per_decade * 7};
			for (int b {}; b <= n_ebin; ++b)
			{
				m_ion_energy_edges.push_back(std::pow(10.0, log_lo
					+ static_cast<double>(b) / per_decade));
			}
			for (auto* v : {&m_ion_energy, &m_sput_ion_energy})
				v->assign(n_ebin, 0.0);
			for (auto* v : {&m_ion_angle, &m_sput_ion_angle})
				v->assign(90, 0.0);
			double total_window {};
			for (int t {}; t < m_nframes; ++t) total_window += window(t);

			int num_above {};
			#pragma omp parallel reduction(+:num_above)
			{
				std::vector<double> ion_e(n_ebin, 0.0), sput_e(n_ebin, 0.0);
				std::vector<double> ion_a(90, 0.0), sput_a(90, 0.0);

				#pragma omp for schedule(dynamic, 256)
				for (std::size_t idx = 0; idx < n_all; ++idx)
				{
					const int point {m_point[idx]};
					if (point < 0) continue;
					const std::size_t e {idx % n_elem};
					const int t {static_cast<int>(idx / n_elem)};
					const double ions {m_flux[idx]
						* m_surface.elements()[e].area * window(t)
						/ total_window};
					double y {}, energy {}, ang {}, above {};
					for (const auto& smp : m_impacts->samples(point))
					{
						const double en {smp.energy * m_te[idx]};
						const double ys {table.get_yield(en, smp.angle)};
						y += smp.weight * ys;
						energy += smp.weight * en;
						ang += smp.weight * smp.angle;
						if (en > table.energy_max()) above += smp.weight;
						if (!(ions > 0.0)) continue;

						const double w {ions * smp.weight};
						const int bin {static_cast<int>(std::floor(
							(std::log10(en) - log_lo) * per_decade))};
						if (bin >= 0 && bin < n_ebin)
						{
							ion_e[bin] += w;
							sput_e[bin] += w * ys;
						}
						const int abin {std::clamp(static_cast<int>(
							smp.angle), 0, 89)};
						ion_a[abin] += w;
						sput_a[abin] += w * ys;
					}
					if (above > 0.01) ++num_above;
					m_yield[idx] = y;
					m_energy[idx] = energy;
					m_impact_angle[idx] = ang;
					m_rate[idx] = m_flux[idx] * y
						* m_surface.elements()[e].area;
				}

				#pragma omp critical
				{
					for (int b {}; b < n_ebin; ++b)
					{
						m_ion_energy[b] += ion_e[b];
						m_sput_ion_energy[b] += sput_e[b];
					}
					for (int b {}; b < 90; ++b)
					{
						m_ion_angle[b] += ion_a[b];
						m_sput_ion_angle[b] += sput_a[b];
					}
				}
			}
			m_num_above_table += num_above;

			std::fill(m_total_rate.begin(), m_total_rate.end(), 0.0);
			for (std::size_t idx {}; idx < n_all; ++idx)
				m_total_rate[idx / n_elem] += m_rate[idx];
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

		// Each element's rate averaged the same way, to compare with how
		// much is redeposited on it (see Deposition)
		m_mean_rate.assign(n_elem, 0.0);
		for (int t {}; t < m_nframes; ++t)
		for (std::size_t e {}; e < n_elem; ++e)
		{
			m_mean_rate[e] += m_rate[t * n_elem + e] * window(t)
				/ total_time;
		}
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
	* (see sample_birth), mean_rate is each element's rate averaged the
	* same way, and frame_lo and frame_hi are the windows.
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
	const std::vector<double>& Source::mean_rate() const {return m_mean_rate;}
	const std::vector<double>& Source::frame_lo() const {return m_frame_lo;}
	const std::vector<double>& Source::frame_hi() const {return m_frame_hi;}

	/**
	* The D+ yield table the source was made with.
	*/
	const YieldTable& Source::table() const {return m_table;}

	/**
	* Computational x and y edges of each element's plasma cell,
	* {x_lo, x_hi, y_lo, y_hi}, one per element.
	*/
	const std::vector<std::array<double, 4>>& Source::cell_bounds() const
		{return m_cell_bounds;}

	/**
	* The floating sheath drop in front of each element (V), without the
	* presheath's 0.5 Te: what the sheath model's field holds, or
	* sheath_drop Te if that's set. Per frame and element, like rate().
	*/
	const std::vector<double>& Source::phi_wall() const {return m_phi_wall;}

	/**
	* The background's toroidal handedness, for the elements' normals
	* (SurfaceElement::normal_xyz).
	*/
	double Source::phi_sign() const {return m_phi_sign;}

	/**
	* Mean angle from the normal the ions hit each element at (degrees),
	* weighted by their share of the flux: from the sheath model's
	* samples with sput_impact_opt = "sheath", or sput_angle_val. Per
	* frame and element, like rate().
	*/
	const std::vector<double>& Source::impact_angle() const
	{
		return m_impact_angle;
	}

	/// The D+ impact table, or nullptr without the sheath model
	const Sheath::ImpactTable* Source::impacts() const
	{
		return m_impacts.get();
	}

	/**
	* @name Ion impact histograms
	* With the sheath model only (empty otherwise): the edges of the
	* energy bins (eV, 20 a decade from 0.01 eV to 100 keV), and the rate
	* ions hit the whole wall in each energy bin and each 1 degree bin of
	* angle from the normal (ions/s), and the same weighted by their
	* yield, the rate the ions in each bin sputter atoms (atoms/s). All
	* averaged over the frames' windows, like source_rate.
	*/
	///@{
	const std::vector<double>& Source::ion_energy_edges() const
		{return m_ion_energy_edges;}
	const std::vector<double>& Source::ion_energy() const
		{return m_ion_energy;}
	const std::vector<double>& Source::ion_angle() const
		{return m_ion_angle;}
	const std::vector<double>& Source::sputtering_ion_energy() const
		{return m_sput_ion_energy;}
	const std::vector<double>& Source::sputtering_ion_angle() const
		{return m_sput_ion_angle;}
	///@}

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
		double sum_angle {};
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
			sum_angle += m_impact_angle[idx] * m_flux[idx] * a;
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
				<< " eV, angle " << sum_angle / sum_flux_area
				<< " degrees from the normal, and yield " << sum_yield
				/ sum_flux_area << " (weighted by the ion flux)\n";
		}
		if (m_impacts)
		{
			os << "  Ion impacts from the sheath model ("
				<< (m_impacts->model() == Sheath::Model::brooks ? "brooks"
				: "borodkina") << "): " << m_impacts->size() << " points, "
				<< m_impacts->num_ions() << " test ions each, "
				<< m_impacts->lost() << " lost, energy conserved to "
				<< m_impacts->max_energy_error() << " of the drop\n";
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
		return Sputtering::sample_impact_energy(table.view(),
			quadrature_view(), e_sheath, ti, angle, u);
	}

	// sample_thompson is in sputtering_hd.h, shared with the GPU

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
		std::array<double, 3> d {};
		Sputtering::cosine_direction(n.data(), u1, u2, d.data());
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
		return Sputtering::sample_birth(view(), rng);
	}

	/**
	* What sample_birth needs, as flat arrays pointing into this source's
	* own vectors (and its surface's and yield table's), for the routine
	* shared with the GPU.
	*
	* @return The view.
	*/
	SourceView Source::view() const
	{
		SourceView v {};
		v.el = m_surface.view();
		v.table = m_table.view();
		v.quad = quadrature_view();
		v.n_idx = m_cdf.size();
		v.cdf = m_cdf.data();
		v.frame_lo = m_frame_lo.data();
		v.frame_hi = m_frame_hi.data();
		static_assert(sizeof(std::array<double, 4>) == 4 * sizeof(double));
		v.cell_bounds = m_cell_bounds.empty() ? nullptr
			: m_cell_bounds.front().data();
		v.e_sheath = m_e_sheath.data();
		v.ti = m_ti.data();
		v.te = m_te.data();
		v.energy = m_energy.data();
		v.impact_angle = m_impact_angle.data();
		v.point = m_point.data();
		v.sample_offset = m_sample_offset.data();
		v.sample_energy = m_sample_energy.data();
		v.sample_angle = m_sample_angle.data();
		v.sample_weight = m_sample_weight.data();
		v.angle = m_angle;
		v.lambda = m_lambda;
		v.imp_mass_kg = m_imp_mass_kg;
		v.phi_sign = m_phi_sign;
		return v;
	}

	namespace
	{
		/// Adds to a CPU tally from inside OpenMP loops, for the routines
		/// shared with the GPU (which add with atomicAdd)
		struct OmpAdd
		{
			void operator()(double* p, const double value) const
			{
				#pragma omp atomic
				*p += value;
			}
		};
	}

	/**
	* Sets up empty tallies for a surface.
	*
	* Each wall cell that borders the plasma gets a list of the surface
	* elements on it (usually one, or two at a corner of the staircase), so
	* an atom absorbed in that wall cell can be put on one of them. Each
	* element on a wall face of its zone (a target) is found instead by its
	* zone, face, cell along the face and toroidal cell, so an atom
	* absorbed at that face can be put on it (add_at_wall_face). Both are
	* kept as sorted keys (cell_key, wall_face_key) searched by bisection,
	* which works the same on the GPU. The time since birth is binned in 10
	* bins per decade from 1e-10 s to 1e-1 s (shorter and longer times go
	* in the first and last bins), and the charge from 0 to the atomic
	* number.
	*
	* @param surface The surface elements. It has to outlive this object.
	* @param atomic_number Atomic number of the impurity.
	* @param phi_sign Handedness of (R, phi, Z), for the impact angles.
	*/
	Deposition::Deposition(const Surface& surface, const int atomic_number,
		const double phi_sign)
		: m_surface {&surface},
		m_redeposited(surface.size(), 0.0),
		m_charge(atomic_number + 1, 0.0),
		m_prompt(surface.size(), 0.0),
		m_prompt_charge(atomic_number + 1, 0.0),
		m_phi_sign {phi_sign},
		m_angle(90, 0.0),
		m_element_energy(surface.size(), 0.0)
	{
		// (key, element) pairs, sorted by key and then element, so each
		// wall cell's elements stay in element order (the first is the
		// default pick, as before)
		std::vector<std::pair<std::uint64_t, int>> cells {};
		std::vector<std::pair<std::uint64_t, int>> faces {};
		for (std::size_t e {}; e < surface.size(); ++e)
		{
			const SurfaceElement& el {surface.elements()[e]};
			const int ie {static_cast<int>(e)};
			if (el.on_wall_face())
			{
				faces.emplace_back(wall_face_key(el.zone,
					static_cast<int>(el.face), el.i, el.j, el.k), ie);
				continue;
			}
			cells.emplace_back(cell_key(el.wall_zone, el.wall_i, el.wall_j,
				el.k), ie);
		}
		std::sort(cells.begin(), cells.end());
		std::sort(faces.begin(), faces.end());

		for (const auto& [key, e] : cells)
		{
			if (m_wall_cell_keys.empty() || m_wall_cell_keys.back() != key)
			{
				m_wall_cell_keys.push_back(key);
				m_wall_cell_start.push_back(static_cast<int>(
					m_wall_cell_elems.size()));
			}
			m_wall_cell_elems.push_back(e);
		}
		m_wall_cell_start.push_back(static_cast<int>(
			m_wall_cell_elems.size()));

		// One element per wall face cell; a later one with the same key
		// replaced the earlier in the map this used to be, so keep the last
		for (const auto& [key, e] : faces)
		{
			if (!m_wall_face_keys.empty() && m_wall_face_keys.back() == key)
			{
				m_wall_face_elems.back() = e;
				continue;
			}
			m_wall_face_keys.push_back(key);
			m_wall_face_elems.push_back(e);
		}

		constexpr int per_decade {10};
		constexpr int decades {9};
		for (int b {}; b <= per_decade * decades; ++b)
		{
			m_lifetime_edges.push_back(std::pow(10.0,
				-10.0 + static_cast<double>(b) / per_decade));
		}
		m_lifetime.assign(m_lifetime_edges.size() - 1, 0.0);

		// Impact energies from 0.01 eV to 100 keV, 20 bins a decade
		for (int b {}; b <= 20 * 7; ++b)
			m_energy_edges.push_back(std::pow(10.0, -2.0 + b / 20.0));
		m_energy.assign(m_energy_edges.size() - 1, 0.0);
	}

	/**
	* The lookups and tallies as flat arrays, pointing into this object's
	* own vectors, for the routines shared with the GPU
	* (sputtering_hd.h). Adding through the view adds to this object.
	*
	* @return The view.
	*/
	DepositionView Deposition::view()
	{
		DepositionView d {};
		d.el = m_surface->view();
		d.phi_sign = m_phi_sign;
		d.n_wall_cells = m_wall_cell_keys.size();
		d.wall_cell_keys = m_wall_cell_keys.data();
		d.wall_cell_start = m_wall_cell_start.data();
		d.wall_cell_elems = m_wall_cell_elems.data();
		d.n_wall_faces = m_wall_face_keys.size();
		d.wall_face_keys = m_wall_face_keys.data();
		d.wall_face_elems = m_wall_face_elems.data();
		d.n_lifetime = static_cast<int>(m_lifetime.size());
		d.lifetime_edges = m_lifetime_edges.data();
		d.n_energy = static_cast<int>(m_energy.size());
		d.energy_edges = m_energy_edges.data();
		d.n_angle = static_cast<int>(m_angle.size());
		d.n_charge = static_cast<int>(m_charge.size());
		d.redeposited = m_redeposited.data();
		d.prompt = m_prompt.data();
		d.element_energy = m_element_energy.data();
		d.lifetime = m_lifetime.data();
		d.charge = m_charge.data();
		d.prompt_charge = m_prompt_charge.data();
		d.energy = m_energy.data();
		d.angle = m_angle.data();
		d.scalars = m_scalars.data();
		return d;
	}

	/**
	* Adds tallies kept elsewhere (a GPU's, copied back to the host) to
	* these. They must be for the same surface and impurity, so the arrays
	* have the same sizes.
	*
	* @param other The tallies to add (only the tally arrays are read).
	*/
	void Deposition::add_tallies(const DepositionView& other)
	{
		auto add = [](std::vector<double>& to, const double* from)
		{
			for (std::size_t n {}; n < to.size(); ++n) to[n] += from[n];
		};
		add(m_redeposited, other.redeposited);
		add(m_prompt, other.prompt);
		add(m_element_energy, other.element_energy);
		add(m_lifetime, other.lifetime);
		add(m_charge, other.charge);
		add(m_prompt_charge, other.prompt_charge);
		add(m_energy, other.energy);
		add(m_angle, other.angle);
		for (int n {}; n < num_dep_scalars; ++n)
			m_scalars[n] += other.scalars[n];
	}

	/**
	* Adds an atom absorbed after moving into a wall cell (chi = 1).
	*
	* The atom is put on a surface element bordering that wall cell: the one
	* whose plasma cell it just came from if there is one, otherwise the
	* first. If no element borders the wall cell (the atom went more than a
	* cell into the wall in one step), it's counted as unassigned. Its time
	* since birth and charge are tallied, and its impact energy and angle if
	* given. The work is done by the routine shared with the GPU
	* (Sputtering::add_into_wall_cell in sputtering_hd.h). Called from
	* inside OpenMP loops, so the adds are atomic.
	*
	* @param zone Zone of the wall cell.
	* @param i x index of the wall cell.
	* @param j y index of the wall cell.
	* @param k z index of the wall cell.
	* @param prev_zone Zone of the cell the atom was in before.
	* @param prev_i x index of the cell it was in before.
	* @param prev_j y index of the cell it was in before.
	* @param pt Time since the atom was born (s).
	* @param q Its charge.
	* @param weight Its weight (atoms/s).
	* @param prompt_by The time since birth up to which its return is
	*        prompt (Slots::Slots::prompt_by). Below 0 (the default) it
	*        isn't.
	* @param impact Its velocity and mass, for the impact tallies, or
	*        nullptr.
	* @return Index of the element it was put on, or -1 if none.
	*/
	long Deposition::add_into_wall_cell(const int zone, const int i,
		const int j, const int k, const int prev_zone, const int prev_i,
		const int prev_j, const double pt, const int q, const double weight,
		const double prompt_by, const Impact* impact)
	{
		return Sputtering::add_into_wall_cell(view(), zone, i, j, k,
			prev_zone, prev_i, prev_j, pt, q, weight, prompt_by,
			impact != nullptr ? impact->v.data() : nullptr,
			impact != nullptr ? impact->mass : 0.0, OmpAdd {});
	}

	/**
	* Adds an atom absorbed at a wall face (e.g., a divertor target at the
	* end of a zone).
	*
	* If the surface has elements on wall faces (a grid fitted to the
	* targets), the atom is put on the one at the cell of the face it
	* crossed: its position along the face and in z picks the cell. z is
	* wrapped first if the zone is periodic in z, since the z faces are
	* checked after the y faces, so z can still be just past the zone's
	* end. Otherwise (SOLEDGE3X, whose wall faces are inside the material,
	* or a face with no element) it's only counted in wall_face. The work
	* is done by the routine shared with the GPU
	* (Sputtering::add_at_wall_face in sputtering_hd.h). Called from inside
	* OpenMP loops, so the adds are atomic.
	*
	* @param zones The background's zones.
	* @param zone Zone whose face the atom crossed.
	* @param face The face it crossed.
	* @param pos Its computational position (x, y, z) in that zone, past
	*        the face.
	* @param pt Time since the atom was born (s).
	* @param q Its charge.
	* @param weight Its weight (atoms/s).
	* @param prompt_by The time since birth up to which its return is
	*        prompt (Slots::Slots::prompt_by). Below 0 (the default) it
	*        isn't.
	* @param impact Its velocity and mass, for the impact tallies, or
	*        nullptr.
	* @return Index of the element it was put on, or -1 if none.
	*/
	long Deposition::add_at_wall_face(const Zones::Zones& zones,
		const int zone, const Zones::Face face,
		const std::array<double, 3>& pos, const double pt, const int q,
		const double weight, const double prompt_by, const Impact* impact)
	{
		return Sputtering::add_at_wall_face(view(), zones.view(), zone,
			static_cast<int>(face), pos.data(), pt, q, weight, prompt_by,
			impact != nullptr ? impact->v.data() : nullptr,
			impact != nullptr ? impact->mass : 0.0, OmpAdd {});
	}

	/**
	* Turns on self-sputtering: the boundary check then passes each atom
	* put on an element to it (SelfSource::rebirth).
	*
	* @param self The self-sputtering source, which has to outlive this
	*        object, or nullptr to turn it off.
	*/
	void Deposition::set_self(SelfSource* self) {m_self = self;}

	/**
	* The self-sputtering source, or nullptr if self-sputtering is off.
	*/
	SelfSource* Deposition::self() const {return m_self;}

	/**
	* Adds up the tallies of every MPI rank onto one of them. Must be
	* called by every rank. The other ranks' tallies are left as they were.
	*
	* @param root Rank to add them up on.
	*/
	void Deposition::reduce(const int root)
	{
		int rank {};
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		auto sum = [&](std::vector<double>& v)
		{
			std::vector<double> total(v.size(), 0.0);
			MPI_Reduce(v.data(), total.data(), static_cast<int>(v.size()),
				MPI_DOUBLE, MPI_SUM, root, MPI_COMM_WORLD);
			if (rank == root) v = total;
		};
		std::vector<double> scalars(m_scalars.begin(), m_scalars.end());
		sum(scalars);
		std::copy(scalars.begin(), scalars.end(), m_scalars.begin());
		sum(m_energy);
		sum(m_angle);
		sum(m_element_energy);
		sum(m_redeposited);
		sum(m_lifetime);
		sum(m_charge);
		sum(m_prompt);
		sum(m_prompt_charge);
	}

	/**
	* @name Accessors
	* The rate atoms are redeposited on each surface element, into wall
	* cells with no element (unassigned) and at wall faces, and their total
	* (atoms/s). Also the edges of the time since birth bins (s), the rate
	* in each bin, and the rate at each charge state, for every redeposited
	* atom. The prompt returns onto each element and at each charge state,
	* and their total (every prompt return, including those that land with
	* no element under them). The mean time since birth of the redeposited
	* atoms (s). The impact energy bin edges (eV), and the rate in each
	* energy bin and each 1 degree bin of the angle from the element's
	* normal, their means (eV, degrees), and the sum of weight * energy on
	* each element (eV atoms/s), for the particles added with an Impact.
	*/
	///@{
	const std::vector<double>& Deposition::energy_edges() const
		{return m_energy_edges;}
	const std::vector<double>& Deposition::impact_energy() const
		{return m_energy;}
	const std::vector<double>& Deposition::impact_angle() const
		{return m_angle;}
	const std::vector<double>& Deposition::element_energy() const
		{return m_element_energy;}
	double Deposition::mean_lifetime() const
	{
		const double tot {total()};
		return tot > 0.0 ? m_scalars[dep_lifetime_sum] / tot : 0.0;
	}
	double Deposition::mean_energy() const
	{
		return m_scalars[dep_energy_weight] > 0.0 ? m_scalars[dep_energy_sum] / m_scalars[dep_energy_weight] : 0.0;
	}
	double Deposition::mean_angle() const
	{
		return m_scalars[dep_angle_weight] > 0.0 ? m_scalars[dep_angle_sum] / m_scalars[dep_angle_weight] : 0.0;
	}
	const std::vector<double>& Deposition::redeposited() const
		{return m_redeposited;}
	double Deposition::unassigned() const {return m_scalars[dep_unassigned];}
	double Deposition::wall_face() const {return m_scalars[dep_wall_face];}
	double Deposition::total() const
	{
		double sum {m_scalars[dep_unassigned] + m_scalars[dep_wall_face]};
		for (const double r : m_redeposited) sum += r;
		return sum;
	}
	const std::vector<double>& Deposition::lifetime_edges() const
		{return m_lifetime_edges;}
	const std::vector<double>& Deposition::lifetime() const
		{return m_lifetime;}
	const std::vector<double>& Deposition::charge() const {return m_charge;}
	const std::vector<double>& Deposition::prompt() const {return m_prompt;}
	const std::vector<double>& Deposition::prompt_charge() const
		{return m_prompt_charge;}
	double Deposition::prompt_total() const
	{
		double sum {};
		for (const double r : m_prompt_charge) sum += r;
		return sum;
	}
	///@}

	/**
	* Prints how much of the sputtered source comes back to the wall, how
	* quickly, and in what charge state.
	*
	* @param os Stream to print to.
	* @param source_rate The sputtered source (atoms/s), to compare with.
	*/
	void Deposition::print_summary(std::ostream& os,
		const double source_rate) const
	{
		const double tot {total()};
		os << "Redeposition: " << tot << " atoms/s back on the wall";
		if (source_rate > 0.0)
			os << ", " << 100.0 * tot / source_rate << "% of the source";
		os << "\n";
		if (!(tot > 0.0)) return;

		// Fraction back within 1 us, 10 us, 100 us and 1 ms of being born.
		// The bin edges are powers of 10^0.1, so these are bin edges.
		os << "  Back within";
		const std::vector<double> times {1e-6, 1e-5, 1e-4, 1e-3};
		for (std::size_t n {}; n < times.size(); ++n)
		{
			double within {};
			for (std::size_t b {}; b < m_lifetime.size(); ++b)
			{
				if (m_lifetime_edges[b + 1] <= times[n] * (1.0 + 1e-9))
					within += m_lifetime[b];
			}
			os << " " << times[n] * 1e6 << " us: " << 100.0 * within / tot
				<< "%" << (n + 1 < times.size() ? "," : "\n");
		}
		double mean_q {};
		for (std::size_t q {}; q < m_charge.size(); ++q)
			mean_q += q * m_charge[q] / tot;
		os << "  Mean charge when they land " << mean_q << ", "
			<< 100.0 * m_charge[0] / tot << "% as neutrals\n";
		os << "  Mean time since birth " << mean_lifetime() * 1e6
			<< " us";
		if (m_scalars[dep_energy_weight] > 0.0)
			os << ", impact energy " << mean_energy() << " eV";
		if (m_scalars[dep_angle_weight] > 0.0)
		{
			os << ", angle from the normal " << mean_angle()
				<< " degrees";
		}
		os << "\n";
		const double prompt {prompt_total()};
		os << "  Prompt (within one singly charged gyroperiod of first "
			"ionizing): " << prompt << " atoms/s";
		if (source_rate > 0.0)
			os << ", " << 100.0 * prompt / source_rate << "% of the source";
		os << "\n";
		if (m_scalars[dep_wall_face] > 0.0 || m_scalars[dep_unassigned] > 0.0)
		{
			os << "  " << 100.0 * m_scalars[dep_wall_face] / tot << "% at wall faces and "
				<< 100.0 * m_scalars[dep_unassigned] / tot << "% into wall cells with no "
				"surface element\n";
		}
	}
	/**
	* Works out the energy and angle an impurity hits a wall element with,
	* for the self-sputtering yield.
	*
	* With the sheath field in the push (phi_wall = 0 here), the particle
	* has already been accelerated through the sheath, so both come from
	* its velocity. Without it, an ion of charge q still falls through the
	* sheath before it hits, gaining q phi_wall, the same assumption as for
	* the D+ ions of the source. That energy is added to its motion along
	* the normal, v_n' = sqrt(v_n^2 + 2 q e phi_wall / m), since the
	* sheath's field is along the normal, and the angle is from there. A
	* particle moving away from the wall when it was absorbed (possible
	* when the staircase's face and the element's normal differ) is taken
	* to hit it with its speed along the normal.
	*
	* @param v Cartesian velocity (m/s).
	* @param mass Mass (kg).
	* @param n The element's unit normal, into the plasma.
	* @param q Charge.
	* @param phi_wall Sheath drop to add (V), 0 with the sheath field on.
	* @return Impact energy (eV) and angle from the normal (degrees).
	*/
	SelfImpact self_impact(const std::array<double, 3>& v, const double mass,
		const std::array<double, 3>& n, const int q, const double phi_wall)
	{
		return self_impact(v.data(), mass, n.data(), q, phi_wall);
	}

	// roulette is in sputtering_hd.h, shared with the GPU

	/**
	* Sets up self-sputtering off the source's surface elements.
	*
	* The table has to be for the impurity on the wall material: its
	* projectile the impurity (atomic number imp_atom_num, mass within 5%
	* of imp_mass_amu) and its target the D+ table's. The largest fraction
	* of its energy the impurity can give a wall atom in one collision is
	* Lambda = 4 m1 m2 / (m1 + m2)^2, 1 for W on W. Roulette is played
	* below sput_self_min_weight times the source particles' weight,
	* source_rate / imp_num.
	*
	* @param source The sputtered source: its elements, their plasma cells,
	*        and the sheath drop in front of them. It has to outlive this
	*        object.
	* @param table Yield table of the impurity on the wall material.
	* @param opts Options (imp_atom_num, imp_mass_amu, imp_num,
	*        sput_self_min_weight, sput_self_max_gen).
	* @param sheath_on Whether the push has the sheath field. If not, an
	*        ion's impact energy gets the sheath drop added (self_impact).
	*/
	SelfSource::SelfSource(const Source& source, const YieldTable& table,
		const Options::Options& opts, const bool sheath_on)
		: m_table {table},
		m_source {&source},
		m_sheath_on {sheath_on},
		m_imp_mass_kg {opts.imp_mass_amu() * Constants::amu_to_kg},
		m_max_gen {opts.sput_self_max_gen()},
		m_self(source.size(), 0.0),
		m_gen_weight(opts.sput_self_max_gen() + 1, 0.0),
		m_gen_landed(opts.sput_self_max_gen() + 1, 0.0)
	{
		if (std::lround(table.projectile_Z()) != opts.imp_atom_num())
		{
			Errors::fatal_error("The self-sputtering table "
				+ table.path() + " is for a projectile of atomic number "
				+ std::to_string(table.projectile_Z()) + ", but the impurity "
				"(imp_atom_num) is " + std::to_string(opts.imp_atom_num())
				+ ".");
		}
		if (std::abs(table.projectile_mass_amu() / opts.imp_mass_amu() - 1.0)
			> 0.05)
		{
			Errors::fatal_error("The projectile mass in the self-sputtering "
				"table " + table.path() + " (" + std::to_string(
				table.projectile_mass_amu()) + " amu) is more than 5% off "
				"the impurity's (imp_mass_amu, " + std::to_string(
				opts.imp_mass_amu()) + " amu).");
		}
		if (std::lround(table.target_Z()) != std::lround(
			source.table().target_Z()))
		{
			Errors::fatal_error("The self-sputtering table " + table.path()
				+ " is for a target of atomic number " + std::to_string(
				table.target_Z()) + ", but the sputtering table "
				+ source.table().path() + " is for " + std::to_string(
				source.table().target_Z()) + ". They need the same wall.");
		}

		const double m1 {table.projectile_mass_amu()};
		const double m2 {table.target_mass_amu()};
		m_lambda = 4.0 * m1 * m2 / ((m1 + m2) * (m1 + m2));
		m_min_weight = opts.sput_self_min_weight() * source.source_rate()
			/ opts.imp_num();
		m_phi_sign = source.phi_sign();
	}

	/**
	* What a particle that landed on element e does next: sputters more
	* atoms and carries on as one of them, or ends.
	*
	* 1. Its impact energy and angle (self_impact), with the sheath drop
	*    added if the push has no sheath field, from the element's plasma
	*    cell at the particle's frame.
	* 2. The yield Y at those from the table. The self-sputtered rate off
	*    the element goes up by w Y, its expected value.
	* 3. Its new weight w Y, through roulette (roulette).
	* 4. If it carries on, it's a neutral atom at the place it hit, moved
	*    a thousandth of a cell back into the element's plasma cell, with
	*    a Thompson energy cut off at Lambda E (sample_thompson) and a
	*    cosine direction about the element's normal (cosine_direction),
	*    as for the D+ source. If that would be generation sput_self_max_gen
	*    + 1, it ends instead and avalanche is set: the chain should have
	*    ended long before (see SelfSource::print_summary).
	*
	* Called from inside OpenMP loops, so the adds are atomic.
	*
	* @param e Index of the element it landed on.
	* @param frame Its background frame (Slots::Slots::tidx).
	* @param pos Its computational position (x, y, z) when absorbed, z
	*        already wrapped into its zone if periodic.
	* @param impact Its velocity and mass.
	* @param q Its charge.
	* @param weight Its weight (atoms/s).
	* @param generation How many times it was self-sputtered before.
	* @param rng Random number generator of the calling thread.
	* @return What it does next.
	*/
	Rebirth SelfSource::rebirth(const std::size_t e, const int frame,
		const std::array<double, 3>& pos, const Impact& impact, const int q,
		const double weight, const int generation, pcg32& rng)
	{
		return Sputtering::rebirth(view(), e, frame, pos.data(),
			impact.v.data(), impact.mass, q, weight, generation, rng,
			OmpAdd {});
	}

	/**
	* Self-sputtering as flat arrays, pointing into this object's and the
	* source's own data, for the routine shared with the GPU (rebirth in
	* sputtering_hd.h).
	*
	* @return The view.
	*/
	SelfView SelfSource::view()
	{
		const SourceView src {m_source->view()};
		SelfView v {};
		v.el = src.el;
		v.cell_bounds = src.cell_bounds;
		v.phi_wall = m_source->phi_wall().data();
		v.nframes = m_source->nframes();
		v.table = m_table.view();
		v.energy_max = m_table.energy_max();
		v.lambda = m_lambda;
		v.imp_mass_kg = m_imp_mass_kg;
		v.phi_sign = m_phi_sign;
		v.min_weight = m_min_weight;
		v.max_gen = m_max_gen;
		v.sheath_on = m_sheath_on;
		v.self = m_self.data();
		v.gen_weight = m_gen_weight.data();
		v.gen_landed = m_gen_landed.data();
		v.scalars = m_scalars.data();
		return v;
	}

	/**
	* Adds up the tallies of every MPI rank onto one of them. Must be
	* called by every rank. The other ranks' tallies are left as they were.
	*
	* @param root Rank to add them up on.
	*/
	void SelfSource::reduce(const int root)
	{
		int rank {};
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		auto sum = [&](std::vector<double>& v)
		{
			std::vector<double> total(v.size(), 0.0);
			MPI_Reduce(v.data(), total.data(), static_cast<int>(v.size()),
				MPI_DOUBLE, MPI_SUM, root, MPI_COMM_WORLD);
			if (rank == root) v = total;
		};
		sum(m_self);
		sum(m_gen_weight);
		sum(m_gen_landed);
		std::vector<double> scalars(m_scalars.begin(), m_scalars.end());
		sum(scalars);
		if (rank != root) return;
		std::copy(scalars.begin(), scalars.end(), m_scalars.begin());
	}

	/**
	* @name Accessors
	* The table, the roulette weight (atoms/s) and the generation limit;
	* the self-sputtered rate off each element and in total (sum of w Y
	* over the landings, atoms/s); the weight carried on into each
	* generation (index g, 1 to sput_self_max_gen; 0 is left for the
	* source) and the weight landing on elements from each generation
	* (0 to sput_self_max_gen); the rate landing on elements and the weight-averaged yield,
	* impact energy (eV) and angle (degrees) of those landings; the rate
	* ended by roulette (its expected weight, w Y); and how many landings
	* were below the threshold or past the table's highest energy.
	*/
	///@{
	const YieldTable& SelfSource::table() const {return m_table;}
	double SelfSource::min_weight() const {return m_min_weight;}
	int SelfSource::max_gen() const {return m_max_gen;}
	const std::vector<double>& SelfSource::self_sputtered() const
		{return m_self;}
	const std::vector<double>& SelfSource::generation_weight() const
		{return m_gen_weight;}
	const std::vector<double>& SelfSource::generation_landed() const
		{return m_gen_landed;}
	double SelfSource::self_total() const
		{return m_scalars[Sputtering::self_total];}
	double SelfSource::landed() const {return m_scalars[self_landed];}
	double SelfSource::mean_yield() const
	{
		return landed() > 0.0 ? m_scalars[Sputtering::self_total] / landed()
			: 0.0;
	}
	double SelfSource::mean_energy() const
	{
		return landed() > 0.0 ? m_scalars[self_energy_sum] / landed() : 0.0;
	}
	double SelfSource::mean_angle() const
	{
		return landed() > 0.0 ? m_scalars[self_angle_sum] / landed() : 0.0;
	}
	double SelfSource::roulette_lost() const
		{return m_scalars[self_roulette_lost];}
	long long SelfSource::num_below_threshold() const
		{return std::llround(m_scalars[self_num_below]);}
	long long SelfSource::num_above_table() const
		{return std::llround(m_scalars[self_num_above]);}
	///@}

	/**
	* The self-sputtering multiplication measured along the chains: the
	* weight carried into generations 2 and up over that into 1 and up.
	* If every generation comes back and sputters the same way, the weight
	* in generation g is S M^g and this is M. It doesn't need the source
	* rate, so it also works on one rank mid-run (for the avalanche
	* message).
	*
	* @return The multiplication, 0 if nothing was self-sputtered.
	*/
	double SelfSource::chain_multiplication() const
	{
		return Sputtering::chain_multiplication(m_gen_weight);
	}

	/**
	* The multiplication along the chains from the weight carried into
	* each generation (see SelfSource::chain_multiplication).
	*
	* @param gen_weight Weight carried into each generation, 0 to
	*        sput_self_max_gen.
	* @return The multiplication, 0 if nothing was self-sputtered.
	*/
	double chain_multiplication(const std::vector<double>& gen_weight)
	{
		double from1 {};
		for (std::size_t g {1}; g < gen_weight.size(); ++g)
			from1 += gen_weight[g];
		if (!(from1 > 0.0)) return 0.0;
		return (from1 - gen_weight[1]) / from1;
	}

	/**
	* Prints the self-sputtered rate, the multiplication and the mean
	* yield, energy and angle of the impurities landing on the wall.
	*
	* The multiplication M (the self-sputtered atoms per atom sputtered,
	* from any source, that are followed) is T / (S + T), with S the D+
	* source and T the self-sputtered rate, since T = S M / (1 - M) for a
	* chain that sputters the same way each generation. It is printed with
	* the chain estimate (chain_multiplication).
	*
	* @param os Stream to print to.
	* @param source_rate The D+ sputtered source (atoms/s).
	*/
	void SelfSource::print_summary(std::ostream& os,
		const double source_rate) const
	{
		const double total {self_total()};
		os << "Self-sputtering: " << total << " atoms/s";
		if (source_rate > 0.0)
		{
			os << ", " << total / source_rate << " times the D+ "
				"source\n  Multiplication " << total / (source_rate
				+ total) << " (along the chains "
				<< chain_multiplication() << ")";
		}
		os << "\n  Landing on elements " << landed() << " atoms/s, mean "
			"yield " << mean_yield() << ", impact energy " << mean_energy()
			<< " eV, angle from the normal " << mean_angle() << " degrees"
			<< (m_sheath_on ? "" : " (with the sheath drop added)") << "\n";
		double carried {};
		int deepest {};
		for (std::size_t g {1}; g < m_gen_weight.size(); ++g)
		{
			carried += m_gen_weight[g];
			if (m_gen_weight[g] > 0.0) deepest = static_cast<int>(g);
		}
		os << "  Carried on " << carried << " atoms/s, down to generation "
			<< deepest << "; " << roulette_lost() << " atoms/s ended by "
			"roulette\n";
		if (num_above_table() > 0)
		{
			os << "  Warning! " << num_above_table() << " landings were past the "
				"self-sputtering table's highest energy (" << m_table
				.energy_max() << " eV), where its yield there was used.\n";
		}
	}

	// -----------------------------------------------------------------
	// Copies for the GPU
	// -----------------------------------------------------------------

#ifdef USE_CUDA
	namespace
	{
		/// Copies n values to the current GPU, nullptr for none
		template <typename T>
		T* copy_to_gpu(const T* src, const std::size_t n, const char* name)
		{
			if (n == 0 || src == nullptr) return nullptr;
			T* dst {};
			Errors::check_cuda(cudaMalloc(&dst, n * sizeof(T)), name);
			Errors::check_cuda(cudaMemcpy(dst, src, n * sizeof(T),
				cudaMemcpyHostToDevice), name);
			return dst;
		}

		/// n zeros on the current GPU, nullptr for none
		double* zeros_on_gpu(const std::size_t n, const char* name)
		{
			if (n == 0) return nullptr;
			double* dst {};
			Errors::check_cuda(cudaMalloc(&dst, n * sizeof(double)), name);
			Errors::check_cuda(cudaMemset(dst, 0, n * sizeof(double)), name);
			return dst;
		}

		/// Copies n values from a GPU (zeros for none)
		std::vector<double> copy_from_gpu(const double* src,
			const std::size_t n, const char* name)
		{
			std::vector<double> out(n, 0.0);
			if (n == 0 || src == nullptr) return out;
			Errors::check_cuda(cudaMemcpy(out.data(), src, n * sizeof(double),
				cudaMemcpyDeviceToHost), name);
			return out;
		}

		/// Frees a device array a view points to (views hold const
		/// pointers, but these were allocated by copy_to_gpu)
		template <typename T>
		void free_gpu(T*& p)
		{
			cudaFree(const_cast<std::remove_const_t<T>*>(p));
			p = nullptr;
		}

		/// Copies the elements to the current GPU
		ElementsView elements_to_gpu(const ElementsView& h)
		{
			ElementsView d {};
			d.n = h.n;
			d.zone = copy_to_gpu(h.zone, h.n, "elements");
			d.i = copy_to_gpu(h.i, h.n, "elements");
			d.j = copy_to_gpu(h.j, h.n, "elements");
			d.k = copy_to_gpu(h.k, h.n, "elements");
			d.face = copy_to_gpu(h.face, h.n, "elements");
			d.n_R = copy_to_gpu(h.n_R, h.n, "elements");
			d.n_Z = copy_to_gpu(h.n_Z, h.n, "elements");
			d.phi_min = copy_to_gpu(h.phi_min, h.n, "elements");
			d.phi_max = copy_to_gpu(h.phi_max, h.n, "elements");
			d.cartesian = copy_to_gpu(h.cartesian, h.n, "elements");
			return d;
		}

		/// Frees elements copied by elements_to_gpu
		void free_elements(ElementsView& d)
		{
			free_gpu(d.zone);
			free_gpu(d.i);
			free_gpu(d.j);
			free_gpu(d.k);
			free_gpu(d.face);
			free_gpu(d.n_R);
			free_gpu(d.n_Z);
			free_gpu(d.phi_min);
			free_gpu(d.phi_max);
			free_gpu(d.cartesian);
			d.n = 0;
		}
	}
#endif

	/**
	* Copies the source to a GPU, for starting sputtered atoms there with
	* the same routine as the CPU (sample_birth in sputtering_hd.h).
	*
	* @param device_id GPU to copy to.
	* @param imp_num Number of particles, so each stands for
	*        source_rate / imp_num atoms/s, as on the CPU.
	* @return The copy, to free with free_source. Without CUDA, an empty
	*         one.
	*/
	SourceDevice Source::to_device(const int device_id, const int imp_num)
		const
	{
		SourceDevice out {};
#ifdef USE_CUDA
		if (m_cdf.empty() || !(m_cdf.back() > 0.0))
		{
			Errors::fatal_error("Can't start sputtered atoms, since nothing "
				"is sputtered (the source is 0 everywhere).");
		}
		Errors::check_cuda(cudaSetDevice(device_id), "Source::to_device");
		const SourceView h {view()};
		SourceView& d {out.view};
		d = h;
		d.el = elements_to_gpu(h.el);
		d.table.log_energy = copy_to_gpu(h.table.log_energy,
			h.table.n_energy, "yield table");
		d.table.angle = copy_to_gpu(h.table.angle, h.table.n_angle,
			"yield table");
		d.table.yield = copy_to_gpu(h.table.yield, static_cast<std::size_t>(
			h.table.n_energy) * h.table.n_angle, "yield table");
		d.quad.x = copy_to_gpu(h.quad.x, h.quad.n, "quadrature");
		d.quad.w = copy_to_gpu(h.quad.w, h.quad.n, "quadrature");
		d.cdf = copy_to_gpu(h.cdf, h.n_idx, "source cdf");
		d.frame_lo = copy_to_gpu(h.frame_lo, m_frame_lo.size(), "frames");
		d.frame_hi = copy_to_gpu(h.frame_hi, m_frame_hi.size(), "frames");
		d.cell_bounds = copy_to_gpu(h.cell_bounds, 4 * h.el.n,
			"cell bounds");
		d.e_sheath = copy_to_gpu(h.e_sheath, h.n_idx, "source");
		d.ti = copy_to_gpu(h.ti, h.n_idx, "source");
		d.te = copy_to_gpu(h.te, h.n_idx, "source");
		d.energy = copy_to_gpu(h.energy, h.n_idx, "source");
		d.impact_angle = copy_to_gpu(h.impact_angle, h.n_idx, "source");
		d.point = copy_to_gpu(h.point, h.n_idx, "source");
		d.sample_offset = copy_to_gpu(h.sample_offset,
			m_sample_offset.size(), "impact samples");
		d.sample_energy = copy_to_gpu(h.sample_energy,
			m_sample_energy.size(), "impact samples");
		d.sample_angle = copy_to_gpu(h.sample_angle,
			m_sample_angle.size(), "impact samples");
		d.sample_weight = copy_to_gpu(h.sample_weight,
			m_sample_weight.size(), "impact samples");
		out.device_id = device_id;
		out.weight = m_source_rate / imp_num;
#else
		(void)device_id;
		(void)imp_num;
#endif
		return out;
	}

	/**
	* Frees a source copied to a GPU (Source::to_device).
	*
	* @param source_d The copy. Left empty.
	*/
	void free_source(SourceDevice& source_d)
	{
#ifdef USE_CUDA
		if (source_d.device_id < 0) return;
		cudaSetDevice(source_d.device_id);
		SourceView& d {source_d.view};
		free_elements(d.el);
		free_gpu(d.table.log_energy);
		free_gpu(d.table.angle);
		free_gpu(d.table.yield);
		free_gpu(d.quad.x);
		free_gpu(d.quad.w);
		free_gpu(d.cdf);
		free_gpu(d.frame_lo);
		free_gpu(d.frame_hi);
		free_gpu(d.cell_bounds);
		free_gpu(d.e_sheath);
		free_gpu(d.ti);
		free_gpu(d.te);
		free_gpu(d.energy);
		free_gpu(d.impact_angle);
		free_gpu(d.point);
		free_gpu(d.sample_offset);
		free_gpu(d.sample_energy);
		free_gpu(d.sample_angle);
		free_gpu(d.sample_weight);
#endif
		source_d = {};
	}

	/**
	* Makes empty tallies on a GPU with this object's lookups and bins, for
	* the boundary kernel to add landings to with the same routines as the
	* CPU (add_at_wall_face and add_into_wall_cell in sputtering_hd.h).
	*
	* @param device_id GPU to put them on.
	* @return The tallies, to add back with add_from_device and free with
	*         free_deposition. Without CUDA, an empty one.
	*/
	DepositionDevice Deposition::to_device(const int device_id)
	{
		DepositionDevice out {};
#ifdef USE_CUDA
		Errors::check_cuda(cudaSetDevice(device_id), "Deposition::to_device");
		const DepositionView h {view()};
		DepositionView& d {out.view};
		d = h;
		d.el = elements_to_gpu(h.el);
		d.wall_cell_keys = copy_to_gpu(h.wall_cell_keys, h.n_wall_cells,
			"wall cells");
		d.wall_cell_start = copy_to_gpu(h.wall_cell_start,
			m_wall_cell_start.size(), "wall cells");
		d.wall_cell_elems = copy_to_gpu(h.wall_cell_elems,
			m_wall_cell_elems.size(), "wall cells");
		d.wall_face_keys = copy_to_gpu(h.wall_face_keys, h.n_wall_faces,
			"wall faces");
		d.wall_face_elems = copy_to_gpu(h.wall_face_elems, h.n_wall_faces,
			"wall faces");
		d.lifetime_edges = copy_to_gpu(h.lifetime_edges,
			static_cast<std::size_t>(h.n_lifetime) + 1, "lifetime bins");
		d.energy_edges = copy_to_gpu(h.energy_edges,
			static_cast<std::size_t>(h.n_energy) + 1, "energy bins");
		d.redeposited = zeros_on_gpu(m_redeposited.size(), "deposition");
		d.prompt = zeros_on_gpu(m_prompt.size(), "deposition");
		d.element_energy = zeros_on_gpu(m_element_energy.size(),
			"deposition");
		d.lifetime = zeros_on_gpu(m_lifetime.size(), "deposition");
		d.charge = zeros_on_gpu(m_charge.size(), "deposition");
		d.prompt_charge = zeros_on_gpu(m_prompt_charge.size(), "deposition");
		d.energy = zeros_on_gpu(m_energy.size(), "deposition");
		d.angle = zeros_on_gpu(m_angle.size(), "deposition");
		d.scalars = zeros_on_gpu(num_dep_scalars, "deposition");
		out.device_id = device_id;
#else
		(void)device_id;
#endif
		return out;
	}

	/**
	* Adds a GPU's tallies (Deposition::to_device, added to by the
	* boundary kernel) to these. Blocks until the GPU's kernels are done.
	*
	* @param dep_d The GPU's tallies.
	*/
	void Deposition::add_from_device(const DepositionDevice& dep_d)
	{
#ifdef USE_CUDA
		if (dep_d.device_id < 0) return;
		Errors::check_cuda(cudaSetDevice(dep_d.device_id),
			"Deposition::add_from_device");
		const DepositionView& d {dep_d.view};
		auto from = [](const double* p, const std::size_t n)
			{return copy_from_gpu(p, n, "Deposition::add_from_device");};
		std::vector<double> redeposited {from(d.redeposited,
			m_redeposited.size())};
		std::vector<double> prompt {from(d.prompt, m_prompt.size())};
		std::vector<double> element_energy {from(d.element_energy,
			m_element_energy.size())};
		std::vector<double> lifetime {from(d.lifetime, m_lifetime.size())};
		std::vector<double> charge {from(d.charge, m_charge.size())};
		std::vector<double> prompt_charge {from(d.prompt_charge,
			m_prompt_charge.size())};
		std::vector<double> energy {from(d.energy, m_energy.size())};
		std::vector<double> angle {from(d.angle, m_angle.size())};
		std::vector<double> scalars {from(d.scalars, num_dep_scalars)};
		DepositionView h {};
		h.redeposited = redeposited.data();
		h.prompt = prompt.data();
		h.element_energy = element_energy.data();
		h.lifetime = lifetime.data();
		h.charge = charge.data();
		h.prompt_charge = prompt_charge.data();
		h.energy = energy.data();
		h.angle = angle.data();
		h.scalars = scalars.data();
		add_tallies(h);
#else
		(void)dep_d;
#endif
	}

	/**
	* Frees tallies on a GPU (Deposition::to_device).
	*
	* @param dep_d The tallies. Left empty.
	*/
	void free_deposition(DepositionDevice& dep_d)
	{
#ifdef USE_CUDA
		if (dep_d.device_id < 0) return;
		cudaSetDevice(dep_d.device_id);
		DepositionView& d {dep_d.view};
		free_elements(d.el);
		free_gpu(d.wall_cell_keys);
		free_gpu(d.wall_cell_start);
		free_gpu(d.wall_cell_elems);
		free_gpu(d.wall_face_keys);
		free_gpu(d.wall_face_elems);
		free_gpu(d.lifetime_edges);
		free_gpu(d.energy_edges);
		free_gpu(d.redeposited);
		free_gpu(d.prompt);
		free_gpu(d.element_energy);
		free_gpu(d.lifetime);
		free_gpu(d.charge);
		free_gpu(d.prompt_charge);
		free_gpu(d.energy);
		free_gpu(d.angle);
		free_gpu(d.scalars);
#endif
		dep_d = {};
	}

	/**
	* Copies self-sputtering to a GPU, with empty tallies, for the boundary
	* kernel to self-sputter with the same routine as the CPU (rebirth in
	* sputtering_hd.h).
	*
	* @param device_id GPU to put it on.
	* @return The copy, to add back with add_from_device and free with
	*         free_self. Without CUDA, an empty one.
	*/
	SelfDevice SelfSource::to_device(const int device_id)
	{
		SelfDevice out {};
#ifdef USE_CUDA
		Errors::check_cuda(cudaSetDevice(device_id),
			"SelfSource::to_device");
		// Only the lookups and settings are copied; the tallies start at 0
		const SelfView h {view()};
		const std::size_t n_elem {h.el.n};
		SelfView& d {out.view};
		d = h;
		d.el = elements_to_gpu(h.el);
		d.cell_bounds = copy_to_gpu(h.cell_bounds, 4 * n_elem,
			"self-sputtering");
		d.phi_wall = copy_to_gpu(h.phi_wall, m_source->phi_wall().size(),
			"self-sputtering");
		d.table.log_energy = copy_to_gpu(h.table.log_energy,
			h.table.n_energy, "self-sputtering table");
		d.table.angle = copy_to_gpu(h.table.angle, h.table.n_angle,
			"self-sputtering table");
		d.table.yield = copy_to_gpu(h.table.yield, static_cast<std::size_t>(
			h.table.n_energy) * h.table.n_angle, "self-sputtering table");
		d.self = zeros_on_gpu(m_self.size(), "self-sputtering");
		d.gen_weight = zeros_on_gpu(m_gen_weight.size(), "self-sputtering");
		d.gen_landed = zeros_on_gpu(m_gen_landed.size(), "self-sputtering");
		d.scalars = zeros_on_gpu(num_self_scalars, "self-sputtering");
		Errors::check_cuda(cudaMalloc(&out.avalanche, sizeof(int)),
			"SelfSource::to_device");
		Errors::check_cuda(cudaMemset(out.avalanche, 0, sizeof(int)),
			"SelfSource::to_device");
		out.device_id = device_id;
#else
		(void)device_id;
#endif
		return out;
	}

	/**
	* Adds a GPU's tallies (SelfSource::to_device, added to by the boundary
	* kernel) to these. Blocks until the GPU's kernels are done.
	*
	* @param self_d The GPU's copy.
	*/
	void SelfSource::add_from_device(const SelfDevice& self_d)
	{
#ifdef USE_CUDA
		if (self_d.device_id < 0) return;
		Errors::check_cuda(cudaSetDevice(self_d.device_id),
			"SelfSource::add_from_device");
		const SelfView& d {self_d.view};
		auto add = [](std::vector<double>& to, const double* from)
		{
			const std::vector<double> v {copy_from_gpu(from, to.size(),
				"SelfSource::add_from_device")};
			for (std::size_t n {}; n < to.size(); ++n) to[n] += v[n];
		};
		add(m_self, d.self);
		add(m_gen_weight, d.gen_weight);
		add(m_gen_landed, d.gen_landed);
		const std::vector<double> scalars {copy_from_gpu(d.scalars,
			num_self_scalars, "SelfSource::add_from_device")};
		for (std::size_t n {}; n < num_self_scalars; ++n)
			m_scalars[n] += scalars[n];
#else
		(void)self_d;
#endif
	}

	/**
	* Frees self-sputtering on a GPU (SelfSource::to_device).
	*
	* @param self_d The copy. Left empty.
	*/
	void free_self(SelfDevice& self_d)
	{
#ifdef USE_CUDA
		if (self_d.device_id < 0) return;
		cudaSetDevice(self_d.device_id);
		SelfView& d {self_d.view};
		free_elements(d.el);
		free_gpu(d.cell_bounds);
		free_gpu(d.phi_wall);
		free_gpu(d.table.log_energy);
		free_gpu(d.table.angle);
		free_gpu(d.table.yield);
		free_gpu(d.self);
		free_gpu(d.gen_weight);
		free_gpu(d.gen_landed);
		free_gpu(d.scalars);
		cudaFree(self_d.avalanche);
#endif
		self_d = {};
	}

	/**
	* Whether a chain on a GPU reached sput_self_max_gen (the boundary
	* kernel sets the flag). Read after each boundary step, as the CPU
	* checks after its own, so the run stops at the same point. Blocks
	* until the GPU's kernels are done.
	*
	* @param self_d Self-sputtering on the GPU.
	* @return True if there was an avalanche.
	*/
	/**
	* The multiplication along the chains, mean yield and mean impact
	* energy of a GPU's self-sputtering so far, for the avalanche message.
	* Blocks until the GPU's kernels are done.
	*
	* @param self_d Self-sputtering on the GPU.
	* @return Those, as SelfSource's accessors would give them.
	*/
	SelfStats self_stats_from_device(const SelfDevice& self_d)
	{
		SelfStats out {};
#ifdef USE_CUDA
		if (self_d.device_id < 0) return out;
		Errors::check_cuda(cudaSetDevice(self_d.device_id),
			"Sputtering::self_stats_from_device");
		const SelfView& d {self_d.view};
		const std::vector<double> gen {copy_from_gpu(d.gen_weight,
			static_cast<std::size_t>(d.max_gen) + 1,
			"Sputtering::self_stats_from_device")};
		const std::vector<double> sc {copy_from_gpu(d.scalars,
			num_self_scalars, "Sputtering::self_stats_from_device")};
		out.chain_multiplication = chain_multiplication(gen);
		if (sc[self_landed] > 0.0)
		{
			out.mean_yield = sc[self_total] / sc[self_landed];
			out.mean_energy = sc[self_energy_sum] / sc[self_landed];
		}
#else
		(void)self_d;
#endif
		return out;
	}

	bool avalanche_on_device(const SelfDevice& self_d)
	{
#ifdef USE_CUDA
		if (self_d.device_id < 0 || self_d.avalanche == nullptr) return false;
		int flag {};
		Errors::check_cuda(cudaMemcpy(&flag, self_d.avalanche, sizeof(int),
			cudaMemcpyDeviceToHost), "Sputtering::avalanche_on_device");
		return flag != 0;
#else
		(void)self_d;
		return false;
#endif
	}
}
