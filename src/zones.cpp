/**
* @file zones.cpp
*
* @brief Zones describe how the simulation volume is split into logically
*        rectangular blocks and how those blocks connect.
*
* Gkeyll backgrounds are a single zone. SOLEDGE3X splits its grid into zones
* to handle X-point geometries (e.g., core, SOL, private flux region and
* divertor legs), and a particle moving between them crosses a zone face.
* Each zone has its own computational grid, and each of its six faces is
* either shared with another zone (or itself, for a periodic face) or is a
* boundary (absorbing or core).
*/

#include <algorithm>
#include <array>
#include <cmath>
#include <string>
#include <utility>
#include <vector>

#include "errors.h"
#include "mpi.h"
#include "utilities.h"
#include "vectors.h"
#include "zones.h"

namespace Zones
{
	/**
	* Creates an empty Zones with no zones and no ghost cells.
	*/
	Zones::Zones()
	{}

	/**
	* Creates an empty Zones whose zones will have ghost cells.
	*
	* @param ghost Number of ghost cells on each side of every zone. This
	*        is the ghost width of the background arrays.
	*/
	Zones::Zones(int ghost)
		: m_ghost {ghost}
	{
		if (ghost < 0)
		{
			Errors::fatal_error("Zones was given a negative number of ghost "
				"cells (" + std::to_string(ghost) + ").");
		}
	}

	/**
	* Adds a zone with the given grid.
	*
	* The cell centers and edges are moved in. Every face of the new zone
	* is absorbing until set with set_face.
	*
	* @param x x cell centers.
	* @param y y cell centers.
	* @param z z cell centers.
	* @param grid_x x cell edges, one more than x.
	* @param grid_y y cell edges, one more than y.
	* @param grid_z z cell edges, one more than z.
	* @return Index of the new zone.
	*/
	int Zones::add_zone(std::vector<double> x, std::vector<double> y,
		std::vector<double> z, std::vector<double> grid_x,
		std::vector<double> grid_y, std::vector<double> grid_z)
	{
		m_x.push_back(std::move(x));
		m_y.push_back(std::move(y));
		m_z.push_back(std::move(z));
		m_grid_x.push_back(std::move(grid_x));
		m_grid_y.push_back(std::move(grid_y));
		m_grid_z.push_back(std::move(grid_z));
		for (auto& faces : m_faces) faces.push_back(face_absorbing);
		m_ghost_cells.clear();

		return get_nzones() - 1;
	}

	/**
	* Sets the code of one face of a zone.
	*
	* @param zone Zone index.
	* @param face Face to set.
	* @param code Zone across the face (>= 0, the zone itself for
	*        periodic), face_absorbing or face_core.
	*/
	void Zones::set_face(const int zone, const Face face, const int code)
	{
		m_faces[face][zone] = code;
		m_ghost_cells.clear();
	}

	/**
	* Returns the code of one face of a zone.
	*
	* @param zone Zone index.
	* @param face Face to get.
	* @return Zone across the face (>= 0, the zone itself for periodic),
	*         face_absorbing or face_core.
	*/
	int Zones::get_face(const int zone, const Face face) const
	{
		return m_faces[face][zone];
	}

	/**
	* @name Size accessors
	*
	* Return the number of zones, the ghost cell width, or the number of
	* cells of a zone (or of every zone) in x, y or z.
	*/
	///@{
	int Zones::get_nzones() const {return std::ssize(m_x);}

	int Zones::get_ghost() const {return m_ghost;}

	int Zones::get_nx(const int zone) const {return std::ssize(m_x[zone]);}

	int Zones::get_ny(const int zone) const {return std::ssize(m_y[zone]);}

	int Zones::get_nz(const int zone) const {return std::ssize(m_z[zone]);}

	std::vector<int> Zones::get_nx_all() const
	{
		std::vector<int> n {};
		for (const auto& v : m_x) n.push_back(std::ssize(v));
		return n;
	}

	std::vector<int> Zones::get_ny_all() const
	{
		std::vector<int> n {};
		for (const auto& v : m_y) n.push_back(std::ssize(v));
		return n;
	}

	std::vector<int> Zones::get_nz_all() const
	{
		std::vector<int> n {};
		for (const auto& v : m_z) n.push_back(std::ssize(v));
		return n;
	}
	///@}

	/**
	* @name Grid accessors
	*
	* Return a zone's cell centers or cell edges in x, y or z.
	*/
	///@{
	const std::vector<double>& Zones::get_x(const int zone) const
		{return m_x[zone];}

	const std::vector<double>& Zones::get_y(const int zone) const
		{return m_y[zone];}

	const std::vector<double>& Zones::get_z(const int zone) const
		{return m_z[zone];}

	const std::vector<double>& Zones::get_grid_x(const int zone) const
		{return m_grid_x[zone];}

	const std::vector<double>& Zones::get_grid_y(const int zone) const
		{return m_grid_y[zone];}

	const std::vector<double>& Zones::get_grid_z(const int zone) const
		{return m_grid_z[zone];}
	///@}

	/**
	* @name Bound accessors
	*
	* Return the first or last cell edge of a zone in x, y or z.
	*/
	///@{
	double Zones::get_x_min(const int zone) const
		{return m_grid_x[zone].front();}

	double Zones::get_x_max(const int zone) const
		{return m_grid_x[zone].back();}

	double Zones::get_y_min(const int zone) const
		{return m_grid_y[zone].front();}

	double Zones::get_y_max(const int zone) const
		{return m_grid_y[zone].back();}

	double Zones::get_z_min(const int zone) const
		{return m_grid_z[zone].front();}

	double Zones::get_z_max(const int zone) const
		{return m_grid_z[zone].back();}
	///@}

	/**
	* @name Accessors by dimension
	*
	* Return the number of cells, or the first or last cell edge, of a zone
	* in dimension dim, where 0 = x, 1 = y and 2 = z. These let code loop
	* over the dimensions (e.g., over the faces, where a face's dimension is
	* face / 2) instead of repeating itself for x, y and z.
	*/
	///@{
	int Zones::get_n(const int zone, const int dim) const
	{
		const std::array<const std::vector<std::vector<double>>*, 3> centers
			{&m_x, &m_y, &m_z};
		return std::ssize((*centers[dim])[zone]);
	}

	double Zones::get_min(const int zone, const int dim) const
	{
		const std::array<const std::vector<std::vector<double>>*, 3> edges
			{&m_grid_x, &m_grid_y, &m_grid_z};
		return (*edges[dim])[zone].front();
	}

	double Zones::get_max(const int zone, const int dim) const
	{
		const std::array<const std::vector<std::vector<double>>*, 3> edges
			{&m_grid_x, &m_grid_y, &m_grid_z};
		return (*edges[dim])[zone].back();
	}
	///@}

	/**
	* Returns the zones that have at least one core face, in increasing
	* order.
	*
	* These are the zones a particle can be put back into when it hits a
	* core boundary.
	*
	* @return Indices of the core zones.
	*/
	std::vector<int> Zones::get_core_zones() const
	{
		std::vector<int> core {};
		for (int zone {}; zone < get_nzones(); ++zone)
		{
			for (int face {}; face < num_faces; ++face)
			{
				if (m_faces[face][zone] == face_core)
				{
					core.push_back(zone);
					break;
				}
			}
		}
		return core;
	}

	/**
	* Returns a zone's cell centers in one dimension.
	*
	* @param zone Zone index.
	* @param dim Dimension: 0 = x, 1 = y, 2 = z.
	* @return The cell centers, without ghost cells.
	*/
	const std::vector<double>& Zones::get_centers(const int zone,
		const int dim) const
	{
		const std::array<const std::vector<std::vector<double>>*, 3> centers
			{&m_x, &m_y, &m_z};
		return (*centers[dim])[zone];
	}

	/**
	* Returns the coordinate of a cell center in a zone's frame, including
	* ghost cells.
	*
	* An interior cell's center is just its center. A ghost cell's center is
	* the center of the cell it's a copy of in the zone across the face,
	* shifted into this zone's frame, e.g., past the max face it's
	* max(this zone) + (center - min(other zone)). This works for any cell
	* widths, and steps through neighbors thinner than the ghost width the
	* same way resolve_cell does. A ghost cell past an absorbing or core
	* face has no center, which stops the run with Errors::fatal_error.
	*
	* @param zone Zone index.
	* @param dim Dimension: 0 = x, 1 = y, 2 = z.
	* @param i Cell index, which may be a ghost cell index.
	* @return Coordinate of the cell center.
	*/
	double Zones::center(const int zone, const int dim, int i) const
	{
		// Interior cell
		if (i >= 0 && i < get_n(zone, dim)) return get_centers(zone, dim)[i];

		// Step across faces, adding up how far each zone's frame is
		// shifted from the next
		double shift {};
		int z {zone};
		while (i < 0 || i >= get_n(z, dim))
		{
			const bool past_max {i >= get_n(z, dim)};
			const int code {m_faces[2 * dim + past_max][z]};
			if (code < 0)
			{
				Errors::fatal_error("Zone " + std::to_string(zone) + " has "
					"no cell center at ghost index " + std::to_string(i)
					+ ", since it's past an absorbing or core face.");
			}

			if (past_max)
			{
				shift += get_max(z, dim) - get_min(code, dim);
				i -= get_n(z, dim);
			}
			else
			{
				shift += get_min(z, dim) - get_max(code, dim);
				i += get_n(code, dim);
			}
			z = code;
		}
		return get_centers(z, dim)[i] + shift;
	}

	/**
	* Returns the cell to interpolate towards from a particle's cell, in one
	* dimension.
	*
	* This is Utilities::get_neighbor_index (the neighbor on the side of
	* the cell center the particle is on, or the inward neighbor at the
	* edge of the zone), except at the edge of a zone with ghost cells. If
	* the particle is on the outer side of an edge cell's center and that
	* face is shared with another zone (or periodic), the ghost cell (-1 or
	* the number of cells) is returned, so interpolation continues into the
	* zone across the face. Without ghost cells this is exactly
	* Utilities::get_neighbor_index.
	*
	* @param zone Zone the particle is in.
	* @param dim Dimension: 0 = x, 1 = y, 2 = z.
	* @param coord Particle's coordinate in that dimension.
	* @param idx Index of the cell the particle is in.
	* @return Index of the neighboring cell, which may be a ghost cell.
	*/
	int Zones::neighbor_index(const int zone, const int dim,
		const double coord, const int idx) const
	{
		const auto& centers {get_centers(zone, dim)};
		if (m_ghost > 0)
		{
			const int n {static_cast<int>(centers.size())};
			const double dx_from_center {coord - centers[idx]};
			if (idx == 0 && dx_from_center <= 0.0
				&& m_faces[2 * dim][zone] >= 0)
			{
				return -1;
			}
			if (idx == n - 1 && dx_from_center > 0.0
				&& m_faces[2 * dim + 1][zone] >= 0)
			{
				return n;
			}
		}
		return Utilities::get_neighbor_index(coord, centers, idx);
	}

	/**
	* Returns the two cells whose centers bracket a coordinate, in one
	* dimension, for interpolating at it.
	*
	* This is Utilities::bracket_indices on the zone's cell centers, which
	* holds at the first or last two cells past the outer centers, except
	* at the edge of a zone with ghost cells. Past the first or last center,
	* if that face is shared with another zone (or periodic), the pair is
	* the ghost cell and the edge cell, so interpolation continues into the
	* zone across the face.
	*
	* @param zone Zone the point is in.
	* @param dim Dimension: 0 = x, 1 = y, 2 = z.
	* @param coord Coordinate of the point in that dimension.
	* @return Indices of the lower and upper cells, either of which may be
	*         a ghost cell.
	*/
	std::pair<int, int> Zones::bracket(const int zone, const int dim,
		const double coord) const
	{
		const auto& centers {get_centers(zone, dim)};
		if (m_ghost > 0)
		{
			const int n {static_cast<int>(centers.size())};
			if (coord < centers.front() && m_faces[2 * dim][zone] >= 0)
				return {-1, 0};
			if (coord > centers.back() && m_faces[2 * dim + 1][zone] >= 0)
				return {n - 1, n};
		}
		return Utilities::bracket_indices(centers, coord);
	}

	/**
	* Returns the zones whose given face is a core face, in increasing
	* order.
	*
	* A particle that hits a core face is put back on one of the core faces
	* of the same kind (e.g., any xmin core face), which these are.
	*
	* @param face Face to check.
	* @return Indices of the zones.
	*/
	std::vector<int> Zones::get_core_zones(const Face face) const
	{
		std::vector<int> core {};
		for (int zone {}; zone < get_nzones(); ++zone)
		{
			if (m_faces[face][zone] == face_core) core.push_back(zone);
		}
		return core;
	}

	/**
	* Checks that the zones are consistent, stopping the run with
	* Errors::fatal_error if not.
	*
	* For every zone, the cell centers and edges must be the right sizes and
	* the edges increasing. Every face code must be a valid zone or face_*
	* code. A face shared with another zone must be shared both ways (the
	* opposite face of the zone across it points back), and the two zones
	* must have the same cell edges along the face. A periodic face (the
	* zone itself) may have a non-periodic opposite face. That's what lets a particle cross the face
	* by only shifting the coordinate normal to it.
	*/
	void Zones::validate() const
	{
		const std::array<const std::vector<std::vector<double>>*, 3> centers
			{&m_x, &m_y, &m_z};
		const std::array<const std::vector<std::vector<double>>*, 3> edges
			{&m_grid_x, &m_grid_y, &m_grid_z};
		const std::array<char, 3> names {'x', 'y', 'z'};
		const std::array<std::string, num_faces> face_names
			{"xmin", "xmax", "ymin", "ymax", "zmin", "zmax"};

		if (get_nzones() == 0) Errors::fatal_error("There are no zones.");

		for (int zone {}; zone < get_nzones(); ++zone)
		{
			const std::string zname {"Zone " + std::to_string(zone)};

			// Grid sizes and ordering
			for (int d {}; d < 3; ++d)
			{
				const auto& c = (*centers[d])[zone];
				const auto& e = (*edges[d])[zone];
				if (c.empty() || e.size() != c.size() + 1)
				{
					Errors::fatal_error(zname + " has "
						+ std::to_string(c.size()) + " " + names[d]
						+ " cell centers and " + std::to_string(e.size())
						+ " edges. There must be at least one cell and one "
						"more edge than centers.");
				}
				for (std::size_t n {1}; n < e.size(); ++n)
				{
					if (e[n] <= e[n - 1])
					{
						Errors::fatal_error(zname + " " + names[d]
							+ " cell edges are not increasing.");
					}
				}
			}

			// Face codes and connectivity
			for (int face {}; face < num_faces; ++face)
			{
				const int code {m_faces[face][zone]};
				const std::string fname {zname + " " + face_names[face]};
				if (code == face_absorbing || code == face_core) continue;
				if (code < 0 || code >= get_nzones())
				{
					Errors::fatal_error(fname + " has invalid code "
						+ std::to_string(code) + ".");
				}

				// The neighbor must point back across its opposite face. A
				// face connecting the zone to itself (periodic) doesn't
				// need its opposite face to be periodic too, since Gkeyll
				// backgrounds have always allowed a periodic bound on one
				// side only.
				const int opposite {face ^ 1};
				if (code != zone && m_faces[opposite][code] != zone)
				{
					Errors::fatal_error(fname + " connects to zone "
						+ std::to_string(code) + ", but that zone's "
						+ face_names[opposite] + " connects to "
						+ std::to_string(m_faces[opposite][code]) + ".");
				}

				// Same cells along the face, i.e., in the two dimensions
				// other than the one normal to the face.
				const int normal {face / 2};
				for (int d {}; d < 3; ++d)
				{
					if (d == normal) continue;
					if ((*centers[d])[zone].size()
						!= (*centers[d])[code].size())
					{
						Errors::fatal_error(fname + " connects to zone "
							+ std::to_string(code) + ", but they have "
							"different numbers of cells in " + names[d]
							+ ".");
					}

					const auto& e = (*edges[d])[zone];
					const auto& e_other = (*edges[d])[code];
					for (std::size_t n {}; n < e.size(); ++n)
					{
						const double tol {1e-9 * std::max(1.0,
							std::abs(e[n]))};
						if (std::abs(e[n] - e_other[n]) > tol)
						{
							Errors::fatal_error(fname + " connects to zone "
								+ std::to_string(code) + ", but they have "
								"different " + names[d] + " cell edges "
								"along it.");
						}
					}
				}
			}
		}
	}

	/**
	* Checks whether a cell's index in one dimension is inside its zone,
	* i.e., not a ghost cell index.
	*
	* @param cell Cell to check.
	* @param dim Dimension to check: 0 = x, 1 = y, 2 = z.
	* @return True if the index is between 0 and the number of cells.
	*/
	bool Zones::in_range(const Cell& cell, const int dim) const
	{
		const std::array<int, 3> idx {cell.i, cell.j, cell.k};
		return idx[dim] >= 0 && idx[dim] < get_n(cell.zone, dim);
	}

	/**
	* Moves a cell across zone faces until its indices are inside a zone,
	* one dimension at a time in the given order.
	*
	* For each dimension in order, while the cell's index is past one of
	* the zone's faces in that dimension, the cell moves into the zone
	* across that face and the index is shifted into that zone's frame
	* (e.g., i < 0 becomes nx(neighbor) + i). Several steps are needed
	* when a neighbor is thinner than the distance past the face. The other
	* two indices don't change, since zones that share a face have the same
	* cells along it.
	*
	* @param[in,out] cell Cell to move. Only valid if true is returned.
	* @param order Dimensions to move in, in order (0 = x, 1 = y, 2 = z).
	* @param ndims Number of entries of order to use.
	* @return False if an absorbing or core face was hit, which has nothing
	*         across it.
	*/
	bool Zones::step_across(Cell& cell, const std::array<int, 3>& order,
		const int ndims) const
	{
		// Far more steps than any real zone layout needs
		constexpr int max_steps {100};

		for (int n {}; n < ndims; ++n)
		{
			const int dim {order[n]};
			std::array<int*, 3> idx {&cell.i, &cell.j, &cell.k};
			int steps {};
			while (!in_range(cell, dim))
			{
				const bool past_max {*idx[dim] >= get_n(cell.zone, dim)};
				const int face {2 * dim + past_max};
				const int code {m_faces[face][cell.zone]};
				if (code < 0) return false;

				// Shift into the frame of the zone across the face
				if (past_max) *idx[dim] -= get_n(cell.zone, dim);
				else *idx[dim] += get_n(code, dim);
				cell.zone = code;

				if (++steps > max_steps)
				{
					Errors::fatal_error("A ghost cell of zone "
						+ std::to_string(cell.zone) + " didn't reach the "
						"inside of a zone after " + std::to_string(max_steps)
						+ " steps across zone faces.");
				}
			}
		}
		return true;
	}

	/**
	* Finds the interior cells that a cell's value comes from.
	*
	* An interior cell is its own source. A ghost cell past one face is a
	* copy of the cell it lands on after stepping across the faces (see
	* step_across). A ghost cell past faces in two or three dimensions (a
	* corner ghost) is stepped across in every order of those dimensions,
	* and its sources are the different cells reached:
	*
	* - Normally every order reaches the same cell, which is the source.
	* - If some orders hit an absorbing or core face, the cells the others
	*   reach are used.
	* - If orders reach different cells, as at an X-point where the cells
	*   meeting there aren't a regular 2x2 block, the value is the average
	*   of them.
	*
	* If no order reaches an interior cell (e.g., past an absorbing face),
	* the source is the nearest interior cell of the cell's own zone, so
	* the ghost copies the edge of the zone.
	*
	* @param cell Cell to resolve.
	* @return The interior cells whose average is the cell's value.
	*/
	std::vector<Cell> Zones::resolve_cell(const Cell& cell) const
	{
		// Dimensions the cell is outside its zone in
		std::array<int, 3> out {};
		int nout {};
		for (int dim {}; dim < 3; ++dim)
		{
			if (!in_range(cell, dim)) out[nout++] = dim;
		}
		if (nout == 0) return {cell};

		// Step across the faces in every order of those dimensions,
		// keeping each different cell reached
		std::vector<Cell> sources {};
		do
		{
			Cell c {cell};
			if (step_across(c, out, nout)
				&& std::ranges::find(sources, c) == sources.end())
			{
				sources.push_back(c);
			}
		} while (std::next_permutation(out.begin(), out.begin() + nout));

		// Nothing across the faces, so use the zone's own edge
		if (sources.empty())
		{
			Cell c {cell};
			c.i = std::clamp(c.i, 0, get_nx(c.zone) - 1);
			c.j = std::clamp(c.j, 0, get_ny(c.zone) - 1);
			c.k = std::clamp(c.k, 0, get_nz(c.zone) - 1);
			sources.push_back(c);
		}
		return sources;
	}

	/**
	* Finds the sources of every ghost cell of every zone with
	* resolve_cell and stores them in m_ghost_cells.
	*/
	void Zones::build_ghost_cells() const
	{
		m_ghost_cells.clear();
		const int g {m_ghost};
		for (int zone {}; zone < get_nzones(); ++zone)
		{
			for (int i {-g}; i < get_nx(zone) + g; ++i)
			for (int j {-g}; j < get_ny(zone) + g; ++j)
			for (int k {-g}; k < get_nz(zone) + g; ++k)
			{
				const Cell cell {zone, i, j, k};
				if (in_range(cell, 0) && in_range(cell, 1)
					&& in_range(cell, 2))
				{
					continue;
				}
				m_ghost_cells.push_back({cell, resolve_cell(cell)});
			}
		}
	}

	/**
	* Fills the ghost cells of a ZonedVector from the zones around them.
	*
	* Each ghost cell gets the value (or average of values) of its sources
	* from resolve_cell, for every entry of the leading dimension (e.g.,
	* every time frame). Copying is correct for scalars, Cartesian vector
	* components and the basis vectors, since crossing a face only shifts a
	* computational coordinate by a constant. The sources are found the
	* first time this is called and reused after that.
	*
	* v must have the same zones and ghost width as this, otherwise the
	* run is stopped with Errors::fatal_error. Without ghost cells (e.g., a
	* Gkeyll background) this does nothing.
	*
	* @tparam T Type of the values.
	* @param v Vector whose ghost cells are filled.
	*/
	template <typename T>
	void Zones::fill_ghosts(Vectors::ZonedVector<T>& v) const
	{
		bool same {v.get_nzones() == get_nzones()
			&& v.get_ghost() == m_ghost};
		for (int zone {}; same && zone < get_nzones(); ++zone)
		{
			same = v.get_nx(zone) == get_nx(zone)
				&& v.get_ny(zone) == get_ny(zone)
				&& v.get_nz(zone) == get_nz(zone);
		}
		if (!same)
		{
			Errors::fatal_error("fill_ghosts was given a ZonedVector whose "
				"zones or ghost width don't match the zones.");
		}
		if (m_ghost == 0) return;

		if (m_ghost_cells.empty()) build_ghost_cells();

		for (int t {}; t < v.get_dim1(); ++t)
		{
			for (const auto& gc : m_ghost_cells)
			{
				const Cell& g {gc.ghost};
				if (gc.sources.size() == 1)
				{
					const Cell& s {gc.sources[0]};
					v(g.zone, t, g.i, g.j, g.k) = v(s.zone, t, s.i, s.j, s.k);
				}
				else
				{
					double sum {};
					for (const auto& s : gc.sources)
						sum += v(s.zone, t, s.i, s.j, s.k);
					v(g.zone, t, g.i, g.j, g.k)
						= static_cast<T>(sum / gc.sources.size());
				}
			}
		}
	}

	/**
	* Broadcasts the zones from rank 0 to every other rank.
	*
	* Must be called by every rank in comm. The other ranks' zones are
	* overwritten with rank 0's.
	*
	* @param comm MPI communicator to broadcast over.
	*/
	void Zones::broadcast(MPI_Comm comm)
	{
		MPI_Bcast(&m_ghost, 1, MPI_INT, 0, comm);

		int nzones {get_nzones()};
		MPI_Bcast(&nzones, 1, MPI_INT, 0, comm);

		for (auto* grids : {&m_x, &m_y, &m_z, &m_grid_x, &m_grid_y,
			&m_grid_z})
		{
			grids->resize(nzones);
			for (auto& v : *grids) Utilities::mpi_broadcast_vector(v, 0, comm);
		}
		for (auto& faces : m_faces)
			Utilities::mpi_broadcast_vector(faces, 0, comm);
		m_ghost_cells.clear();
	}

	/**
	* Converts a Gkeyll-style *bound_type option into a face code.
	*
	* The bound types are 0 = absorbing, 1 = periodic and 2 = core. A
	* periodic face connects the zone to itself. Any other value stops the
	* run with Errors::fatal_error.
	*
	* @param bound_type Value of a *bound_type_int option.
	* @param zone Zone the face belongs to.
	* @return Face code for set_face.
	*/
	int face_from_bound_type(const int bound_type, const int zone)
	{
		switch (bound_type)
		{
			case 0: return face_absorbing;
			case 1: return zone;
			case 2: return face_core;
		}
		Errors::fatal_error("Bound type " + std::to_string(bound_type)
			+ " has no matching face code.");
	}
}

// Instantiate the templates since the definitions are in this file
template void Zones::Zones::fill_ghosts<int>(
	Vectors::ZonedVector<int>& v) const;
template void Zones::Zones::fill_ghosts<float>(
	Vectors::ZonedVector<float>& v) const;
template void Zones::Zones::fill_ghosts<double>(
	Vectors::ZonedVector<double>& v) const;
