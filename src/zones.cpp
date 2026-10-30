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

#include <string>
#include <utility>
#include <vector>

#include "errors.h"
#include "mpi.h"
#include "utilities.h"
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
	* Checks that the zones are consistent, stopping the run with
	* Errors::fatal_error if not.
	*
	* For every zone, the cell centers and edges must be the right sizes and
	* the edges increasing. Every face code must be a valid zone or face_*
	* code. A shared face must be shared both ways (the opposite face of
	* the zone across it points back), and the two zones must have the same
	* number of cells along the face.
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

				// The neighbor must point back across its opposite face
				const int opposite {face ^ 1};
				if (m_faces[opposite][code] != zone)
				{
					Errors::fatal_error(fname + " connects to zone "
						+ std::to_string(code) + ", but that zone's "
						+ face_names[opposite] + " connects to "
						+ std::to_string(m_faces[opposite][code]) + ".");
				}

				// Same number of cells along the face, i.e., in the two
				// dimensions other than the one normal to the face.
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
