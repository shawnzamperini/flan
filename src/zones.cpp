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
* boundary (absorbing, core or wall).
*/

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <numbers>
#include <string>
#include <utility>
#include <vector>

#include "errors.h"
#include "mpi.h"
#include "utilities.h"
#include "vectors.h"
#include "zones.h"

#ifdef USE_CUDA
#include <cuda_runtime.h>
#endif

namespace Zones
{
	/**
	* Creates an empty Zones with no zones and no ghost cells.
	*/
	Zones::Zones()
	{
		rebuild_flat();
	}

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
		rebuild_flat();
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
		rebuild_flat();

		return get_nzones() - 1;
	}

	/**
	* Sets the code of one face of a zone.
	*
	* @param zone Zone index.
	* @param face Face to set.
	* @param code Zone across the face (>= 0, the zone itself for
	*        periodic), face_absorbing, face_core, face_wall or face_lcfs.
	*/
	void Zones::set_face(const int zone, const Face face, const int code)
	{
		m_faces[face][zone] = code;
		m_ghost_cells.clear();
		rebuild_flat();
	}

	/**
	* Returns the code of one face of a zone.
	*
	* @param zone Zone index.
	* @param face Face to get.
	* @return Zone across the face (>= 0, the zone itself for periodic),
	*         face_absorbing, face_core, face_wall or face_lcfs.
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
	* same way resolve_cell does. A ghost cell past a boundary (absorbing,
	* core or wall) face has no center, which stops the run with
	* Errors::fatal_error.
	*
	* Within the ghost width the center comes from the table in the flat
	* arrays (Zones::center in zones_hd.h, which the GPU uses too). Further
	* out it's found by stepping across faces (step_to_center), which is
	* also how the table is filled.
	*
	* @param zone Zone index.
	* @param dim Dimension: 0 = x, 1 = y, 2 = z.
	* @param i Cell index, which may be a ghost cell index.
	* @return Coordinate of the cell center.
	*/
	double Zones::center(const int zone, const int dim, int i) const
	{
		const bool in_table {i >= -m_ghost && i < get_n(zone, dim) + m_ghost};
		const double c {in_table ? ::Zones::center(m_flat.view(), zone, dim, i)
			: step_to_center(zone, dim, i)};
		if (std::isnan(c))
		{
			Errors::fatal_error("Zone " + std::to_string(zone) + " has "
				"no cell center at ghost index " + std::to_string(i)
				+ ", since it's past a boundary face.");
		}
		return c;
	}

	/**
	* Finds a cell center in a zone's frame by stepping across the zone's
	* faces (see center).
	*
	* @param zone Zone index.
	* @param dim Dimension: 0 = x, 1 = y, 2 = z.
	* @param i Cell index, which may be a ghost cell index.
	* @return Coordinate of the cell center, or NaN if it's past a boundary
	*         face (or a face to a zone that hasn't been added yet).
	*/
	double Zones::step_to_center(const int zone, const int dim, int i) const
	{
		// Step across faces, adding up how far each zone's frame is
		// shifted from the next
		double shift {};
		int z {zone};
		while (i < 0 || i >= get_n(z, dim))
		{
			const bool past_max {i >= get_n(z, dim)};
			const int code {m_faces[2 * dim + past_max][z]};
			if (code < 0 || code >= get_nzones())
				return std::numeric_limits<double>::quiet_NaN();

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
	* Utilities::get_neighbor_index. The work is done by
	* Zones::neighbor_index in zones_hd.h, which the GPU uses too.
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
		return ::Zones::neighbor_index(m_flat.view(), zone, dim, coord, idx);
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
	* zone across the face. The work is done by Zones::bracket in
	* zones_hd.h, which the GPU uses too.
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
		const IndexPair p {::Zones::bracket(m_flat.view(), zone, dim, coord)};
		return {p.lo, p.hi};
	}

	/**
	* The zones as flat arrays on the host, for the functions in zones_hd.h.
	* The view is valid until the zones next change.
	*
	* @return The view.
	*/
	const ZonesView& Zones::view() const {return m_flat.view();}

	/**
	* Rebuilds the flat arrays (m_flat) from the nested vectors. Called
	* whenever the zones change, so view() is always current.
	*
	* Ghost cell centers come from step_to_center, so they're the same
	* values center() gave before there was a table. A face to a zone that
	* hasn't been added yet gives NaN until the zones are complete.
	*/
	void Zones::rebuild_flat()
	{
		FlatZones f {};
		const int nzones {get_nzones()};
		f.nzones = nzones;
		f.ghost = m_ghost;

		const std::array<const std::vector<std::vector<double>>*, 3> edges
			{&m_grid_x, &m_grid_y, &m_grid_z};
		for (int dim {}; dim < 3; ++dim)
		{
			for (int zone {}; zone < nzones; ++zone)
			{
				const int n {get_n(zone, dim)};
				f.n[dim].push_back(n);
				f.center_offset[dim].push_back(std::ssize(f.centers[dim]));
				for (int i {-m_ghost}; i < n + m_ghost; ++i)
					f.centers[dim].push_back(step_to_center(zone, dim, i));
				f.edge_offset[dim].push_back(std::ssize(f.edges[dim]));
				const auto& e {(*edges[dim])[zone]};
				f.edges[dim].insert(f.edges[dim].end(), e.begin(), e.end());
			}
		}

		f.cell_offset.push_back(0);
		for (int zone {}; zone < nzones; ++zone)
		{
			for (int face {}; face < num_faces; ++face)
				f.faces.push_back(m_faces[face][zone]);

			std::size_t cells {1};
			for (int dim {}; dim < 3; ++dim)
			{
				cells *= static_cast<std::size_t>(get_n(zone, dim)
					+ 2 * m_ghost);
			}
			f.cell_offset.push_back(f.cell_offset.back() + cells);

			f.flags.push_back((is_toroidal(zone) ? flag_toroidal : 0)
				| (is_cartesian(zone) ? flag_cartesian : 0)
				| (is_wedge(zone) ? flag_wedge : 0));

			if (has_z_shift(zone))
			{
				f.z_shift_offset.push_back(std::ssize(f.z_shift));
				f.z_shift.insert(f.z_shift.end(), m_z_shift[zone].begin(),
					m_z_shift[zone].end());
			}
			else
			{
				f.z_shift_offset.push_back(-1);
			}
		}

		// The zones with a core face of each kind, for putting particles
		// back in at the core (Boundary::core_face)
		f.core_offset.push_back(0);
		for (int face {}; face < num_faces; ++face)
		{
			for (const int zone : get_core_zones(static_cast<Face>(face)))
				f.core_zones.push_back(zone);
			f.core_offset.push_back(std::ssize(f.core_zones));
		}

		m_flat = std::move(f);
		m_flat.refresh_view();
	}

	/**
	* Copies the flat arrays and points this copy's view at them.
	*
	* @param other FlatZones to copy.
	*/
	FlatZones::FlatZones(const FlatZones& other)
		: nzones {other.nzones}
		, ghost {other.ghost}
		, n {other.n}
		, faces {other.faces}
		, center_offset {other.center_offset}
		, centers {other.centers}
		, edge_offset {other.edge_offset}
		, edges {other.edges}
		, cell_offset {other.cell_offset}
		, flags {other.flags}
		, z_shift_offset {other.z_shift_offset}
		, z_shift {other.z_shift}
		, core_offset {other.core_offset}
		, core_zones {other.core_zones}
	{
		refresh_view();
	}

	/**
	* Copies the flat arrays and points this view at them.
	*
	* @param other FlatZones to copy.
	* @return This.
	*/
	FlatZones& FlatZones::operator=(const FlatZones& other)
	{
		if (this != &other)
		{
			FlatZones copy {other};
			*this = std::move(copy);
		}
		return *this;
	}

	/**
	* Points the view at the arrays. Needed after the arrays are built or
	* copied, since a std::vector's buffer moves when it's copied or
	* reallocated.
	*/
	void FlatZones::refresh_view()
	{
		m_view.nzones = nzones;
		m_view.ghost = ghost;
		for (int dim {}; dim < 3; ++dim)
		{
			m_view.n[dim] = n[dim].data();
			m_view.center_offset[dim] = center_offset[dim].data();
			m_view.centers[dim] = centers[dim].data();
			m_view.edge_offset[dim] = edge_offset[dim].data();
			m_view.edges[dim] = edges[dim].data();
		}
		m_view.faces = faces.data();
		m_view.cell_offset = cell_offset.data();
		m_view.flags = flags.data();
		m_view.z_shift_offset = z_shift_offset.data();
		m_view.z_shift = z_shift.data();
		m_view.core_offset = core_offset.data();
		m_view.core_zones = core_zones.data();
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
				if (code == face_absorbing || code == face_core
					|| code == face_wall) continue;
				if (code == face_lcfs)
				{
					if (face == zmin || face == zmax) continue;
					Errors::fatal_error(fname + " is an lcfs face, which "
						"only z faces can be.");
				}
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
	* @return False if a boundary (absorbing, core or wall) face was hit,
	*         which has nothing across it.
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
	* - If some orders hit a boundary (absorbing, core or wall) face, the
	*   cells the others reach are used.
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
	* Returns every ghost cell of every zone and its sources (see
	* resolve_cell), the same ones fill_ghosts uses.
	*
	* @return The ghost cells.
	*/
	const std::vector<GhostCell>& Zones::get_ghost_cells() const
	{
		if (m_ghost_cells.empty()) build_ghost_cells();
		return m_ghost_cells;
	}

	/**
	* Fills the ghost cells of a ZonedVector from the zones around them.
	*
	* Each ghost cell gets the value (or average of values) of its sources
	* from resolve_cell, for every entry of the leading dimension (e.g.,
	* every time frame). Copying is correct for scalars, and for Cartesian
	* vectors as long as no zone is a toroidal wedge (see is_wedge), since
	* crossing a face only shifts a computational coordinate by a
	* constant. On a wedge, the (X, Y) parts of Cartesian vectors need the
	* overload that rotates them. The sources are found the first time this
	* is called and reused after that.
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
	* Fills the ghost cells of the X and Y components of a Cartesian vector
	* (e.g., B, E or the reciprocal basis vectors), rotating them on a
	* toroidal wedge.
	*
	* On a zone whose z is the toroidal angle (see set_toroidal), a ghost
	* cell's source can be at a different toroidal angle than the ghost
	* itself: past a periodic z face of a wedge 2 pi / 5 wide, the ghost
	* is the source moved a fifth of the way around the torus. The source's
	* Cartesian vector is then rotated about Z by phi_sign (z_ghost -
	* z_source) before it's copied, so the ghost holds the vector at its
	* own position. Sources are rotated one at a time before averaging.
	* The angles are the cell centers along z, continued across faces for
	* the ghost (see center).
	*
	* Without a wedge (any_wedge), the rotation would be by 0 or a whole
	* turn, so the components are copied as fill_ghosts(v) does, bit for
	* bit. The Z component is a scalar under the rotation and uses
	* fill_ghosts(v).
	*
	* @tparam T Type of the values.
	* @param vX X component of the vector, ghost cells filled in place.
	* @param vY Y component of the vector, ghost cells filled in place.
	* @param phi_sign +1 if (R, z, Z) is right handed, -1 if not (see
	*        Background::Background::get_phi_sign).
	*/
	template <typename T>
	void Zones::fill_ghosts(Vectors::ZonedVector<T>& vX,
		Vectors::ZonedVector<T>& vY, const double phi_sign) const
	{
		if (!any_wedge())
		{
			fill_ghosts(vX);
			fill_ghosts(vY);
			return;
		}

		// Checks the shapes
		fill_ghosts(vX);
		fill_ghosts(vY);

		for (const auto& gc : m_ghost_cells)
		{
			const Cell& g {gc.ghost};
			if (!is_wedge(g.zone)) continue;

			// Rotation of each source, the same in every time frame
			std::vector<double> cs {}, sn {};
			const double zg {center(g.zone, 2, g.k)};
			for (const auto& s : gc.sources)
			{
				const double ang {phi_sign * (zg - center(s.zone, 2, s.k))};
				cs.push_back(std::cos(ang));
				sn.push_back(std::sin(ang));
			}

			for (int t {}; t < vX.get_dim1(); ++t)
			{
				double sum_x {}, sum_y {};
				for (std::size_t n {}; n < gc.sources.size(); ++n)
				{
					const Cell& s {gc.sources[n]};
					const double x {vX(s.zone, t, s.i, s.j, s.k)};
					const double y {vY(s.zone, t, s.i, s.j, s.k)};
					sum_x += cs[n] * x - sn[n] * y;
					sum_y += sn[n] * x + cs[n] * y;
				}
				const double ns {static_cast<double>(gc.sources.size())};
				vX(g.zone, t, g.i, g.j, g.k) = static_cast<T>(sum_x / ns);
				vY(g.zone, t, g.i, g.j, g.k) = static_cast<T>(sum_y / ns);
			}
		}
	}

	/**
	* Sets a zone's field-line shift table, for field-aligned interpolation
	* and derivatives on a shifted-metric grid (e.g., BOUT++/HERMES).
	*
	* On such a grid z is the toroidal angle, and the field line through
	* cell (i, j) at angle z passes through cell (i, j + 1) at angle
	* z + zShift(i, j + 1) - zShift(i, j). Turbulence is stretched along
	* the field, so a quantity at (i, j + 1) is far better matched by its
	* value along the field line than at the same z, which may be many
	* toroidal cells (or several turns of a toroidal wedge) away.
	*
	* The table covers the zone's cells including the ghost cells, as
	* (nx + 2 ghost) x (ny + 2 ghost) values in C order with ghost index
	* -ghost first. A ghost cell's value must continue the field line from
	* the zone's own cells (e.g., include the core's twist across the
	* branch cut), since the ghost cells' plasma values are copied at the
	* same z from the zone they come from. The zone's z faces must be
	* periodic to itself, since shifted points are wrapped in z.
	*
	* @param zone Zone to set the table for.
	* @param table zShift (radians of z) on the zone's cells and ghosts.
	*/
	void Zones::set_z_shift(const int zone, std::vector<double> table)
	{
		const std::size_t want {static_cast<std::size_t>(
			(get_nx(zone) + 2 * m_ghost) * (get_ny(zone) + 2 * m_ghost))};
		if (table.size() != want)
		{
			Errors::fatal_error("The zShift table of zone "
				+ std::to_string(zone) + " has " + std::to_string(
				table.size()) + " values, but the zone has "
				+ std::to_string(want) + " cells including ghost cells.");
		}
		if (m_faces[zmin][zone] != zone || m_faces[zmax][zone] != zone)
		{
			Errors::fatal_error("Zone " + std::to_string(zone) + " was "
				"given a zShift table, but it isn't periodic in z.");
		}
		if (std::ssize(m_z_shift) < get_nzones())
			m_z_shift.resize(get_nzones());
		m_z_shift[zone] = std::move(table);
		rebuild_flat();
	}

	/**
	* @param zone Zone to check.
	* @return True if the zone has a zShift table (see set_z_shift).
	*/
	bool Zones::has_z_shift(const int zone) const
	{
		return zone < std::ssize(m_z_shift) && !m_z_shift[zone].empty();
	}

	/**
	* A zone's zShift at a cell (see set_z_shift). The zone must have a
	* table.
	*
	* @param zone Zone of the cell.
	* @param i x index of the cell, which may be a ghost index.
	* @param j y index of the cell, which may be a ghost index.
	* @return zShift (radians of z).
	*/
	double Zones::z_shift(const int zone, const int i, const int j) const
	{
		// The shared version (zones_hd.h), which the GPU uses
		return ::Zones::z_shift(m_flat.view(), zone, i, j);
	}

	/**
	* @param zone Zone to get the period of.
	* @return The zone's extent in z, which is its period in z when its z
	*         faces are periodic.
	*/
	double Zones::z_period(const int zone) const
	{
		return m_grid_z[zone].back() - m_grid_z[zone].front();
	}

	/**
	* Wraps a z coordinate into a zone's z range, for a zone that is
	* periodic in z.
	*
	* @param zone Zone whose range to wrap into.
	* @param z Coordinate to wrap, any number of periods away.
	* @return z plus a whole number of periods, in [z_min, z_min + period).
	*/
	double Zones::wrap_z(const int zone, const double z) const
	{
		// The shared version (zones_hd.h), which the GPU uses
		return ::Zones::wrap_z(m_flat.view(), zone, z);
	}

	/**
	* Marks a zone's z as the toroidal angle in radians, as on SOLEDGE3X
	* and HERMES grids (not Gkeyll or the slab, where z is a length).
	*
	* Cartesian vectors then turn with z: a vector at angle z moved to
	* angle z' is rotated about Z by phi_sign (z' - z). This matters when
	* the zone is a wedge (is_wedge), whose periodic z faces join points
	* that are a wedge's width apart around the torus, so ghost cells
	* (fill_ghosts(vX, vY, phi_sign)) and particles that wrap in z
	* (Boundary::check_bounds_cpu) rotate their vectors.
	*
	* @param zone Zone whose z is the toroidal angle.
	*/
	void Zones::set_toroidal(const int zone)
	{
		if (zone < 0 || zone >= get_nzones())
		{
			Errors::fatal_error("Zones::set_toroidal: there is no zone "
				+ std::to_string(zone) + ".");
		}
		if (std::ssize(m_toroidal) < get_nzones())
			m_toroidal.resize(get_nzones());
		m_toroidal[zone] = 1;
		rebuild_flat();
	}

	/**
	* @param zone Zone to check.
	* @return True if the zone's z is the toroidal angle (set_toroidal).
	*/
	bool Zones::is_toroidal(const int zone) const
	{
		return zone < std::ssize(m_toroidal) && m_toroidal[zone] == 1;
	}

	/**
	* Marks a zone as Cartesian: its x, y and z are the Cartesian X, Y and
	* Z in meters, as on the straight slab. Wall surfaces on its faces
	* (Sputtering::Surface) can then be found without a toroidal angle.
	* Gkeyll zones aren't Cartesian, even though their z is a length, since
	* x and y there are flux coordinates.
	*
	* @param zone Zone whose coordinates are Cartesian.
	*/
	void Zones::set_cartesian(const int zone)
	{
		if (zone < 0 || zone >= get_nzones())
		{
			Errors::fatal_error("Zones::set_cartesian: there is no zone "
				+ std::to_string(zone) + ".");
		}
		if (std::ssize(m_cartesian) < get_nzones())
			m_cartesian.resize(get_nzones());
		m_cartesian[zone] = 1;
		rebuild_flat();
	}

	/**
	* @param zone Zone to check.
	* @return True if the zone's coordinates are Cartesian (set_cartesian).
	*/
	bool Zones::is_cartesian(const int zone) const
	{
		return zone < std::ssize(m_cartesian) && m_cartesian[zone] == 1;
	}

	/**
	* Whether a zone is a toroidal wedge: its z is the toroidal angle, it's
	* periodic in z, and its period isn't a whole turn (to within 1e-9 of
	* 2 pi). Only then do vectors need rotating across its z faces, since
	* on a whole torus the rotation is by 2 pi.
	*
	* @param zone Zone to check.
	* @return True if the zone is a wedge of the torus.
	*/
	bool Zones::is_wedge(const int zone) const
	{
		if (!is_toroidal(zone)) return false;
		if (get_face(zone, zmin) != zone || get_face(zone, zmax) != zone)
			return false;
		const double two_pi {2.0 * std::numbers::pi};
		return std::abs(z_period(zone) - two_pi) > 1e-9 * two_pi;
	}

	/**
	* @return True if any zone is a toroidal wedge (is_wedge).
	*/
	bool Zones::any_wedge() const
	{
		for (int zone {}; zone < get_nzones(); ++zone)
			if (is_wedge(zone)) return true;
		return false;
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

		// zShift tables, sized on rank 0 (empty if no zone has one)
		int ntables {static_cast<int>(m_z_shift.size())};
		MPI_Bcast(&ntables, 1, MPI_INT, 0, comm);
		m_z_shift.resize(ntables);
		for (auto& v : m_z_shift) Utilities::mpi_broadcast_vector(v, 0, comm);
		Utilities::mpi_broadcast_vector(m_toroidal, 0, comm);
		Utilities::mpi_broadcast_vector(m_cartesian, 0, comm);
		m_ghost_cells.clear();
		rebuild_flat();
	}

#ifdef USE_CUDA
	namespace
	{
		/**
		* Copies a vector to a new array on the current GPU.
		*
		* @param v Values to copy.
		* @param what What the values are, for the error message.
		* @return Device pointer, or nullptr if v is empty.
		*/
		template <typename T>
		const T* copy_to_device(const std::vector<T>& v, const char* what)
		{
			if (v.empty()) return nullptr;
			T* d {};
			Errors::check_cuda(cudaMalloc(&d, v.size() * sizeof(T)), what);
			Errors::check_cuda(cudaMemcpy(d, v.data(), v.size() * sizeof(T),
				cudaMemcpyHostToDevice), what);
			return d;
		}
	}
#endif

	/**
	* Copies the flat arrays (view()) to a GPU.
	*
	* @param device_id ID of the GPU to copy to.
	* @return The same view with pointers to the GPU's copies, to be freed
	*         with free_zones. Stops the run if GPU support wasn't compiled
	*         in.
	*/
	ZonesView Zones::to_device(int device_id) const
	{
		ZonesView d {m_flat.view()};
#ifdef USE_CUDA
		Errors::check_cuda(cudaSetDevice(device_id), "Zones::to_device");
		for (int dim {}; dim < 3; ++dim)
		{
			d.n[dim] = copy_to_device(m_flat.n[dim], "Zones::to_device: n");
			d.center_offset[dim] = copy_to_device(m_flat.center_offset[dim],
				"Zones::to_device: center offsets");
			d.centers[dim] = copy_to_device(m_flat.centers[dim],
				"Zones::to_device: centers");
			d.edge_offset[dim] = copy_to_device(m_flat.edge_offset[dim],
				"Zones::to_device: edge offsets");
			d.edges[dim] = copy_to_device(m_flat.edges[dim],
				"Zones::to_device: edges");
		}
		d.faces = copy_to_device(m_flat.faces, "Zones::to_device: faces");
		d.cell_offset = copy_to_device(m_flat.cell_offset,
			"Zones::to_device: cell offsets");
		d.flags = copy_to_device(m_flat.flags, "Zones::to_device: flags");
		d.z_shift_offset = copy_to_device(m_flat.z_shift_offset,
			"Zones::to_device: zShift offsets");
		d.z_shift = copy_to_device(m_flat.z_shift,
			"Zones::to_device: zShift");
		d.core_offset = copy_to_device(m_flat.core_offset,
			"Zones::to_device: core offsets");
		d.core_zones = copy_to_device(m_flat.core_zones,
			"Zones::to_device: core zones");
#else
		(void)device_id;
		Errors::fatal_error("Zones::to_device was called, but GPU support "
			"wasn't compiled in.");
#endif
		return d;
	}

	/**
	* Frees the GPU arrays of a view from Zones::to_device and sets its
	* pointers to nullptr. Does nothing without GPU support.
	*
	* @param v View from Zones::to_device.
	* @param device_id ID of the GPU it's on.
	*/
	void free_zones(ZonesView& v, int device_id)
	{
#ifdef USE_CUDA
		cudaSetDevice(device_id);
		auto release = [](auto*& p)
		{
			cudaFree(const_cast<void*>(static_cast<const void*>(p)));
			p = nullptr;
		};
		for (int dim {}; dim < 3; ++dim)
		{
			release(v.n[dim]);
			release(v.center_offset[dim]);
			release(v.centers[dim]);
			release(v.edge_offset[dim]);
			release(v.edges[dim]);
		}
		release(v.faces);
		release(v.cell_offset);
		release(v.flags);
		release(v.z_shift_offset);
		release(v.z_shift);
		release(v.core_offset);
		release(v.core_zones);
#else
		(void)v;
		(void)device_id;
#endif
	}

	/**
	* Converts a Gkeyll-style *bound_type option into a face code.
	*
	* The bound types are 0 = absorbing, 1 = periodic, 2 = core and 3 =
	* lcfs (z faces only). A periodic face connects the zone to itself. Any other value stops the
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
			case 3: return face_lcfs;
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
template void Zones::Zones::fill_ghosts<float>(
	Vectors::ZonedVector<float>& vX, Vectors::ZonedVector<float>& vY,
	const double phi_sign) const;
template void Zones::Zones::fill_ghosts<double>(
	Vectors::ZonedVector<double>& vX, Vectors::ZonedVector<double>& vY,
	const double phi_sign) const;
