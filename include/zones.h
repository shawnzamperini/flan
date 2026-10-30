/**
* @file zones.h
*
* @brief Header file for zones.cpp
*/

#pragma once

#include <array>
#include <vector>

#ifndef __CUDACC__
#include "mpi.h"
#endif

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
	* Describes how the simulation volume is split into logically
	* rectangular zones and how the zones connect.
	*
	* Each zone has its own computational grid (cell centers and edges in x,
	* y and z). Each face of a zone has a code: a value >= 0 is the zone
	* across that face, where the zone itself means the face is periodic,
	* and a negative value is one of the face_* boundary codes. Zones that
	* share a face share it entirely and have the same cell indices along
	* it, so a particle crossing a face only needs its coordinate shifted.
	*
	* A Gkeyll background is a single zone whose faces come from the
	* *bound_type options. A SOLEDGE3X background has one zone per
	* SOLEDGE3X zone. See zones.cpp for the member function docs.
	*/
	class Zones
	{
	private:
		int m_ghost {};		///< Ghost cells on each side of every zone
		std::vector<std::vector<double>> m_x {};		///< x cell centers
		std::vector<std::vector<double>> m_y {};		///< y cell centers
		std::vector<std::vector<double>> m_z {};		///< z cell centers
		std::vector<std::vector<double>> m_grid_x {};	///< x cell edges
		std::vector<std::vector<double>> m_grid_y {};	///< y cell edges
		std::vector<std::vector<double>> m_grid_z {};	///< z cell edges
		std::array<std::vector<int>, num_faces> m_faces {};	///< Face codes

	public:

		Zones();

		explicit Zones(int ghost);

		int add_zone(std::vector<double> x, std::vector<double> y,
			std::vector<double> z, std::vector<double> grid_x,
			std::vector<double> grid_y, std::vector<double> grid_z);

		void set_face(const int zone, const Face face, const int code);

		int get_face(const int zone, const Face face) const;

		int get_nzones() const;

		int get_ghost() const;

		int get_nx(const int zone) const;

		int get_ny(const int zone) const;

		int get_nz(const int zone) const;

		std::vector<int> get_nx_all() const;

		std::vector<int> get_ny_all() const;

		std::vector<int> get_nz_all() const;

		const std::vector<double>& get_x(const int zone) const;

		const std::vector<double>& get_y(const int zone) const;

		const std::vector<double>& get_z(const int zone) const;

		const std::vector<double>& get_grid_x(const int zone) const;

		const std::vector<double>& get_grid_y(const int zone) const;

		const std::vector<double>& get_grid_z(const int zone) const;

		double get_x_min(const int zone) const;

		double get_x_max(const int zone) const;

		double get_y_min(const int zone) const;

		double get_y_max(const int zone) const;

		double get_z_min(const int zone) const;

		double get_z_max(const int zone) const;

		std::vector<int> get_core_zones() const;

		void validate() const;

#ifndef __CUDACC__
		void broadcast(MPI_Comm comm);
#endif
	};

	int face_from_bound_type(const int bound_type, const int zone);
}
