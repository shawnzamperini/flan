/**
* @file background.cpp
*
* @brief Contains class definition for a loaded in Gkeyll background plasma. 
*/
#include <string>
#include <vector>

#include "background.h"
#include "background_device.h"
#include "errors.h"
#include "flan_types.h"
#include "mpi.h"
#include "utilities.h"
#include "vectors.h"
#include "zones.h"


/**
* @namespace Background
* @brief Namespace to contain anything related to the background plasma
*/
namespace Background
{
	/**
	* @name Time and computational grid accessors
	*
	* Return the background times, the zones and each zone's computational
	* (curvilinear) coordinates. The x, y and z vectors hold a zone's cell
	* centers, and the grid_x, grid_y and grid_z vectors hold its cell
	* edges. The background times are stored at "time centers", so there is
	* no corresponding time grid.
	*/
	///@{
	const std::vector<BkgFPType>& Background::get_times() const
		{return m_times;}
	const Zones::Zones& Background::get_zones() const
		{return m_zones;}
	int Background::get_nzones() const
		{return m_zones.get_nzones();}
	const std::vector<BkgFPType>& Background::get_x(const int zone) const
		{return m_zones.get_x(zone);}
	const std::vector<BkgFPType>& Background::get_y(const int zone) const
		{return m_zones.get_y(zone);}
	const std::vector<BkgFPType>& Background::get_z(const int zone) const
		{return m_zones.get_z(zone);}
	const std::vector<BkgFPType>& Background::get_grid_x(const int zone) const
		{return m_zones.get_grid_x(zone);}
	const std::vector<BkgFPType>& Background::get_grid_y(const int zone) const
		{return m_zones.get_grid_y(zone);}
	const std::vector<BkgFPType>& Background::get_grid_z(const int zone) const
		{return m_zones.get_grid_z(zone);}
	///@}

	/**
	* @name Plasma and field accessors
	*
	* Return the plasma and field data, each stored as a function of
	* (zone, t, x, y, z). These are the electron density (ne), electron and ion
	* temperatures (te, ti), plasma potential (vp), magnetic field
	* components and magnitude (bX, bY, bZ, bmag), magnetic field gradient
	* (gradbX, gradbY, gradbZ), electric field components and magnitude
	* (eX, eY, eZ, emag), electric field gradient (gradeX, gradeY, gradeZ)
	* and background flow (uX, uY, uZ). Vector components are in physical
	* (Cartesian) X, Y and Z. bR is only kept for debugging.
	*
	* Most have a const version and a non-const version that allows the
	* data to be modified in place.
	*/
	///@{
	const Vectors::ZonedVector<BkgFPType>& Background::get_ne() const
		{return m_ne;}
	Vectors::ZonedVector<BkgFPType>& Background::get_ne() 
		{return m_ne;}
	const Vectors::ZonedVector<BkgFPType>& Background::get_te() const
		{return m_te;}
	const Vectors::ZonedVector<BkgFPType>& Background::get_ti() const
		{return m_ti;}
	const Vectors::ZonedVector<BkgFPType>& Background::get_vp() const
		{return m_vp;}
	const Vectors::ZonedVector<BkgFPType>& Background::get_bX() const
		{return m_bX;}
	const Vectors::ZonedVector<BkgFPType>& Background::get_bY() const
		{return m_bY;}
	const Vectors::ZonedVector<BkgFPType>& Background::get_bZ() const
		{return m_bZ;}
	const Vectors::ZonedVector<BkgFPType>& Background::get_bR() const
		{return m_bR;}
	const Vectors::ZonedVector<BkgFPType>& Background::get_bmag() const
		{return m_bmag;}
	Vectors::ZonedVector<BkgFPType>& Background::get_te() 
		{return m_te;}
	Vectors::ZonedVector<BkgFPType>& Background::get_ti() 
		{return m_ti;}
	Vectors::ZonedVector<BkgFPType>& Background::get_vp() 
		{return m_vp;}
	Vectors::ZonedVector<BkgFPType>& Background::get_bX() 
		{return m_bX;}
	Vectors::ZonedVector<BkgFPType>& Background::get_bY() 
		{return m_bY;}
	Vectors::ZonedVector<BkgFPType>& Background::get_bZ() 
		{return m_bZ;}
	Vectors::ZonedVector<BkgFPType>& Background::get_bR() 
		{return m_bR;}
	Vectors::ZonedVector<BkgFPType>& Background::get_bmag() 
		{return m_bmag;}
	const Vectors::ZonedVector<BkgFPType>& Background::get_gradbX() const
		{return m_gradbX;}
	const Vectors::ZonedVector<BkgFPType>& Background::get_gradbY() const
		{return m_gradbY;}
	const Vectors::ZonedVector<BkgFPType>& Background::get_gradbZ() const
		{return m_gradbZ;}
	const Vectors::ZonedVector<BkgFPType>& Background::get_eX() const
		{return m_eX;}
	const Vectors::ZonedVector<BkgFPType>& Background::get_eY() const
		{return m_eY;}
	const Vectors::ZonedVector<BkgFPType>& Background::get_eZ() const
		{return m_eZ;}
	const Vectors::ZonedVector<BkgFPType>& Background::get_emag() const
		{return m_emag;}
	Vectors::ZonedVector<BkgFPType>& Background::get_eX() 
		{return m_eX;}
	Vectors::ZonedVector<BkgFPType>& Background::get_eY() 
		{return m_eY;}
	Vectors::ZonedVector<BkgFPType>& Background::get_eZ() 
		{return m_eZ;}
	Vectors::ZonedVector<BkgFPType>& Background::get_emag() 
		{return m_emag;}
	const Vectors::ZonedVector<BkgFPType>& Background::get_gradeX() const
		{return m_gradeX;}
	const Vectors::ZonedVector<BkgFPType>& Background::get_gradeY() const
		{return m_gradeY;}
	const Vectors::ZonedVector<BkgFPType>& Background::get_gradeZ() const
		{return m_gradeZ;}
	const Vectors::ZonedVector<BkgFPType>& Background::get_uX() const
		{return m_uX;}
	const Vectors::ZonedVector<BkgFPType>& Background::get_uY() const
		{return m_uY;}
	const Vectors::ZonedVector<BkgFPType>& Background::get_uZ() const
		{return m_uZ;}
	Vectors::ZonedVector<BkgFPType>& Background::get_uX() 
		{return m_uX;}
	Vectors::ZonedVector<BkgFPType>& Background::get_uY() 
		{return m_uY;}
	Vectors::ZonedVector<BkgFPType>& Background::get_uZ() 
		{return m_uZ;}
	///@}

	/**
	* @name Physical coordinate and Jacobian accessors
	*
	* Return the physical (Cartesian) coordinates of the grid, each stored
	* as a function of (zone, x, y, z). X, Y and Z hold the coordinates of the
	* cell centers, and grid_X, grid_Y and grid_Z hold the coordinates of
	* the cell edges. J is the Jacobian of the computational coordinate
	* system.
	*/
	///@{
	const Vectors::ZonedVector<BkgFPType>& Background::get_X() const
		{return m_X;}
	const Vectors::ZonedVector<BkgFPType>& Background::get_Y() const
		{return m_Y;}
	const Vectors::ZonedVector<BkgFPType>& Background::get_Z() const
		{return m_Z;}
	const Vectors::ZonedVector<BkgFPType>& Background::get_grid_X() const
		{return m_grid_X;}
	const Vectors::ZonedVector<BkgFPType>& Background::get_grid_Y() const
		{return m_grid_Y;}
	const Vectors::ZonedVector<BkgFPType>& Background::get_grid_Z() const
		{return m_grid_Z;}
	const Vectors::ZonedVector<BkgFPType>& Background::get_J() const
		{return m_J;}
	///@}

	/**
	* @name Metric coefficient accessors
	*
	* Return the six independent components of the metric tensor g_ij of
	* the computational coordinate system, each stored as a function of
	* (zone, x, y, z). The indices 0, 1 and 2 correspond to x, y and z.
	*/
	///@{
	const Vectors::ZonedVector<BkgFPType>& Background::get_gij_00() const
		{return m_gij_00;}
	const Vectors::ZonedVector<BkgFPType>& Background::get_gij_01() const
		{return m_gij_01;}
	const Vectors::ZonedVector<BkgFPType>& Background::get_gij_02() const
		{return m_gij_02;}
	const Vectors::ZonedVector<BkgFPType>& Background::get_gij_11() const
		{return m_gij_11;}
	const Vectors::ZonedVector<BkgFPType>& Background::get_gij_12() const
		{return m_gij_12;}
	const Vectors::ZonedVector<BkgFPType>& Background::get_gij_22() const
		{return m_gij_22;}
	///@}

	/**
	* @name Basis vector accessors
	*
	* Return the components of the basis vectors that convert between the
	* computational (x, y, z) and physical (X, Y, Z) coordinates, each
	* stored as a function of (zone, x, y, z).
	*
	* The reciprocal (contravariant) basis vectors are the derivatives of
	* the computational coordinates with respect to the physical ones
	* (dxdX, dxdY, ..., dzdZ). The tangent (covariant) basis vectors are
	* the derivatives of the physical coordinates with respect to the
	* computational ones (dXdx, dYdx, ..., dZdz). Only the reciprocal basis
	* vectors have non-const versions.
	*/
	///@{
	const Vectors::ZonedVector<BkgFPType>& Background::get_dxdX() const
		{return m_dxdX;}
	const Vectors::ZonedVector<BkgFPType>& Background::get_dxdY() const
		{return m_dxdY;}
	const Vectors::ZonedVector<BkgFPType>& Background::get_dxdZ() const
		{return m_dxdZ;}
	const Vectors::ZonedVector<BkgFPType>& Background::get_dydX() const
		{return m_dydX;}
	const Vectors::ZonedVector<BkgFPType>& Background::get_dydY() const
		{return m_dydY;}
	const Vectors::ZonedVector<BkgFPType>& Background::get_dydZ() const
		{return m_dydZ;}
	const Vectors::ZonedVector<BkgFPType>& Background::get_dzdX() const
		{return m_dzdX;}
	const Vectors::ZonedVector<BkgFPType>& Background::get_dzdY() const
		{return m_dzdY;}
	const Vectors::ZonedVector<BkgFPType>& Background::get_dzdZ() const
		{return m_dzdZ;}

	const Vectors::ZonedVector<BkgFPType>& Background::get_dXdx() const
		{return m_dXdx;}
	const Vectors::ZonedVector<BkgFPType>& Background::get_dYdx() const
		{return m_dYdx;}
	const Vectors::ZonedVector<BkgFPType>& Background::get_dZdx() const
		{return m_dZdx;}
	const Vectors::ZonedVector<BkgFPType>& Background::get_dXdy() const
		{return m_dXdy;}
	const Vectors::ZonedVector<BkgFPType>& Background::get_dYdy() const
		{return m_dYdy;}
	const Vectors::ZonedVector<BkgFPType>& Background::get_dZdy() const
		{return m_dZdy;}
	const Vectors::ZonedVector<BkgFPType>& Background::get_dXdz() const
		{return m_dXdz;}
	const Vectors::ZonedVector<BkgFPType>& Background::get_dYdz() const
		{return m_dYdz;}
	const Vectors::ZonedVector<BkgFPType>& Background::get_dZdz() const
		{return m_dZdz;}

	Vectors::ZonedVector<BkgFPType>& Background::get_dxdX() {return m_dxdX;}
	Vectors::ZonedVector<BkgFPType>& Background::get_dxdY() {return m_dxdY;}
	Vectors::ZonedVector<BkgFPType>& Background::get_dxdZ() {return m_dxdZ;}
	Vectors::ZonedVector<BkgFPType>& Background::get_dydX() {return m_dydX;}
	Vectors::ZonedVector<BkgFPType>& Background::get_dydY() {return m_dydY;}
	Vectors::ZonedVector<BkgFPType>& Background::get_dydZ() {return m_dydZ;}
	Vectors::ZonedVector<BkgFPType>& Background::get_dzdX() {return m_dzdX;}
	Vectors::ZonedVector<BkgFPType>& Background::get_dzdY() {return m_dzdY;}
	Vectors::ZonedVector<BkgFPType>& Background::get_dzdZ() {return m_dzdZ;}
	///@}

	/**
	* @name Covariant magnetic field accessors
	*
	* Return the covariant components of the magnetic field (b_x, b_y,
	* b_z) in the computational coordinate system, each stored as a
	* function of (zone, x, y, z). These are only kept for debugging and post
	* analysis.
	*/
	///@{
	const Vectors::ZonedVector<BkgFPType>& Background::get_b_x() const
		{return m_b_x;}
	const Vectors::ZonedVector<BkgFPType>& Background::get_b_y() const
		{return m_b_y;}
	const Vectors::ZonedVector<BkgFPType>& Background::get_b_z() const
		{return m_b_z;}
	///@}

	/**
	* Returns the ion parallel conductive heat flux q_par (W/m^2), stored
	* as a function of (zone, t, x, y, z). It is the heat conducted along
	* B in the frame moving with the ion flow (not the energy carried by
	* the flow), positive along B like the parallel flow. It sets how
	* skewed the ion distribution is, eps = q_par / (n T v_th), which is
	* what gives the thermal force in the collision model. Only SOLEDGE3X
	* backgrounds whose files have the face fluxes have it, otherwise it
	* is empty.
	*/
	const Vectors::ZonedVector<BkgFPType>& Background::get_q_par_i() const
		{return m_q_par_i;}

	/**
	* Returns the electron parallel conductive heat flux q_par (W/m^2),
	* like get_q_par_i but for the electrons. It skews the electron
	* distribution, which gives the electron thermal force.
	*/
	const Vectors::ZonedVector<BkgFPType>& Background::get_q_par_e() const
		{return m_q_par_e;}

	/**
	* Returns the wall mask, stored as a function of (zone, x, y, z). It is
	* 1 in cells inside the wall and 0 in the plasma. Only SOLEDGE3X
	* backgrounds have one, otherwise it is empty.
	*/
	const Vectors::ZonedVector<BkgFPType>& Background::get_chi() const
		{return m_chi;}

	/**
	* Returns the wall outlines in the poloidal plane. Only SOLEDGE3X
	* backgrounds whose mesh file has a /walls group have them, otherwise
	* it is empty.
	*/
	const std::vector<WallPolygon>& Background::get_walls() const
		{return m_walls;}

	/**
	* Sets the wall outlines in the poloidal plane.
	*
	* @param walls Wall outlines, each a polyline of (R, Z) points (m).
	*/
	void Background::set_walls(const std::vector<WallPolygon>& walls)
		{m_walls = walls;}

	/**
	* @name Species masses
	*
	* Get or set the masses of the background electrons and main ions
	* (amu), which the impurities collide with. The reader sets them, e.g.,
	* from the gkyl_*_mass_amu options for Gkeyll or from the plasma files
	* for SOLEDGE3X.
	*/
	///@{
	double Background::get_elec_mass_amu() const {return m_elec_mass_amu;}

	double Background::get_ion_mass_amu() const {return m_ion_mass_amu;}

	void Background::set_elec_mass_amu(const double elec_mass_amu)
		{m_elec_mass_amu = elec_mass_amu;}

	void Background::set_ion_mass_amu(const double ion_mass_amu)
		{m_ion_mass_amu = ion_mass_amu;}
	///@}

	/**
	* Stops the run with Errors::fatal_error if the background does not
	* have exactly one zone.
	*
	* Used by the parts of the code that only handle a single zone so far.
	*
	* @param caller Name of the function or feature that needs a single
	*        zone (for the message).
	*/
	void Background::check_single_zone(const std::string_view caller) const
	{
		if (m_zones.get_nzones() != 1)
		{
			Errors::fatal_error(std::string(caller) + " only works with a "
				"single zone, but the background has "
				+ std::to_string(m_zones.get_nzones()) + ".");
		}
	}

	/**
	* @name Dimension accessors
	*
	* Return the dimensions of the 4D plasma data. dim1 is the number of
	* times, and dim2, dim3 and dim4 are the number of cells in x, y and z.
	* A dim1 of 0 means it has not been set yet. dim2, dim3 and dim4 are
	* only defined for a single-zone background, and stop the run with
	* Errors::fatal_error otherwise. Use get_zones for the per-zone sizes.
	*/
	///@{
	int Background::get_dim1() const {return m_dim1;}

	int Background::get_dim2() const
	{
		check_single_zone("Background::get_dim2");
		return m_zones.get_nx(0);
	}

	int Background::get_dim3() const
	{
		check_single_zone("Background::get_dim3");
		return m_zones.get_ny(0);
	}

	int Background::get_dim4() const
	{
		check_single_zone("Background::get_dim4");
		return m_zones.get_nz(0);
	}
	///@}

	/**
	* Checks that a new dimension matches the dimension already stored.
	*
	* A stored dimension of 0 means it has not been set yet, so any new
	* dimension is accepted. Otherwise the two must be equal.
	*
	* @param m_dim Dimension currently stored in the Background.
	* @param in_dim Dimension of the data being loaded.
	* @param dim_num Which dimension (1-4) is being checked (for the
	*                message).
	* @return Empty if the dimensions are consistent. Otherwise a line
	*         describing the mismatch, for set_dims to put in its error
	*         message.
	*/
	std::string Background::check_dim(const int m_dim, const int in_dim,
		const int dim_num)
	{
		// Assume a zero dimension value just means it's not set.
		if (m_dim != 0 && m_dim != in_dim)
		{
			return "\n  dim" + std::to_string(dim_num) + ": already set to "
				+ std::to_string(m_dim) + ", but this data has " 
				+ std::to_string(in_dim);
		}
		return {};
	}

	/**
	* Checks the dimensions of a single-zone Vector4D against the times and
	* the zone, and stores the number of times if they are consistent.
	*
	* The zone must already be set with move_into_zones, and there must be
	* exactly one, without ghost cells. The first dimension of v is checked against the number
	* of times already loaded with check_dim, and the others against the
	* zone's number of cells. If any don't match, the run is stopped with
	* Errors::fatal_error, listing every dimension that doesn't match. A
	* mismatch means the data is not being read in correctly or is from
	* different Gkeyll runs.
	*
	* @tparam T Data type stored in the Vector4D.
	* @param v Vector4D whose dimensions are checked.
	* @param data Name of the data being loaded (for the error message).
	*/
	template <typename T>
	void Background::set_dims(const Vectors::Vector4D<T>& v,
		const std::string_view data)
	{
		check_single_zone("Loading " + std::string(data)
			+ " as a Vector4D");
		if (m_zones.get_ghost() != 0)
		{
			Errors::fatal_error("In Background, " + std::string(data)
				+ " was loaded as a Vector4D, which has no ghost cells, but "
				"the zones have ghost cells.");
		}

		// Check dimensions first, collecting a line for each one that
		// doesn't match.
		const std::string mismatches {check_dim(m_dim1, v.get_dim1(), 1)
			+ check_dim(m_zones.get_nx(0), v.get_dim2(), 2)
			+ check_dim(m_zones.get_ny(0), v.get_dim3(), 3)
			+ check_dim(m_zones.get_nz(0), v.get_dim4(), 4)};

		if (!mismatches.empty())
		{
			Errors::fatal_error("In Background, " + std::string(data) 
				+ " has different dimensions from what has already been "
				"set. This indicates the data is not being read in "
				"correctly or is from different Gkeyll runs." + mismatches);
		}

		// Then overwrite the number of times. This either sets it from 0
		// to its correct value, or overwrites it with the same value.
		m_dim1 = v.get_dim1();
	}

	/**
	* Checks the shape of a ZonedVector against the zones and the times,
	* and stores the number of times if they are consistent.
	*
	* The zones must already be set with move_into_zones. The data must
	* have the same number of zones and the same number of cells in each
	* as the zones, and the zones' ghost width. Data at the cell corners
	* (grid_X, grid_Y and grid_Z) has one more value in each dimension
	* instead, and no ghost cells. Time-dependent data must have the same
	* number of times as data already loaded (checked with check_dim),
	* and other data a leading dimension of 1. If anything doesn't match,
	* the run is stopped with Errors::fatal_error, listing every mismatch.
	*
	* @param v Data whose shape is checked.
	* @param data Name of the data being loaded (for the error message).
	* @param time_dependent True for data stored as (zone, t, x, y, z),
	*        false for data stored as (zone, x, y, z).
	* @param corners True for data at the cell corners.
	*/
	void Background::set_dims(const Vectors::ZonedVector<BkgFPType>& v,
		const std::string_view data, const bool time_dependent,
		const bool corners)
	{
		if (m_zones.get_nzones() == 0)
		{
			Errors::fatal_error("In Background, " + std::string(data)
				+ " was loaded before the zones. Load the zones first with "
				"move_into_zones.");
		}

		std::string mismatches {};
		const int extra {corners ? 1 : 0};
		const int ghost {corners ? 0 : m_zones.get_ghost()};
		if (v.get_nzones() != m_zones.get_nzones())
		{
			mismatches += "\n  " + std::to_string(v.get_nzones())
				+ " zones, but there are "
				+ std::to_string(m_zones.get_nzones());
		}
		else
		{
			for (int z {}; z < v.get_nzones(); ++z)
			{
				if (v.get_nx(z) != m_zones.get_nx(z) + extra
					|| v.get_ny(z) != m_zones.get_ny(z) + extra
					|| v.get_nz(z) != m_zones.get_nz(z) + extra)
				{
					mismatches += "\n  zone " + std::to_string(z)
						+ " has the wrong number of cells";
				}
			}
		}
		if (v.get_ghost() != ghost)
		{
			mismatches += "\n  ghost width is "
				+ std::to_string(v.get_ghost()) + ", but should be "
				+ std::to_string(ghost);
		}
		if (time_dependent) mismatches += check_dim(m_dim1, v.get_dim1(), 1);
		else if (v.get_dim1() != 1)
		{
			mismatches += "\n  dim1 is " + std::to_string(v.get_dim1())
				+ ", but should be 1";
		}

		if (!mismatches.empty())
		{
			Errors::fatal_error("In Background, " + std::string(data)
				+ " doesn't match the zones or the data already loaded."
				+ mismatches);
		}

		if (time_dependent) m_dim1 = v.get_dim1();
	}

	/**
	* @name Range accessors
	*
	* Return the minimum and maximum values of the background times, or of
	* a zone's x, y and z grid edges. The time range is computed when the
	* times are set with move_into_times.
	*/
	///@{
	BkgFPType Background::get_t_min() const {return m_t_min;}
	BkgFPType Background::get_t_max() const {return m_t_max;}
	BkgFPType Background::get_x_min(const int zone) const
		{return m_zones.get_x_min(zone);}
	BkgFPType Background::get_x_max(const int zone) const
		{return m_zones.get_x_max(zone);}
	BkgFPType Background::get_y_min(const int zone) const
		{return m_zones.get_y_min(zone);}
	BkgFPType Background::get_y_max(const int zone) const
		{return m_zones.get_y_max(zone);}
	BkgFPType Background::get_z_min(const int zone) const
		{return m_zones.get_z_min(zone);}
	BkgFPType Background::get_z_max(const int zone) const
		{return m_zones.get_z_max(zone);}
	///@}

	/**
	* @name Time and zone setters
	*
	* Move the background times and the zones (each with its computational
	* grid) into the Background with move semantics, so no copy is made.
	* move_into_times also stores the minimum and maximum times, which are
	* returned by get_t_min and get_t_max. move_into_zones checks the zones
	* with Zones::Zones::validate first.
	*/
	///@{
	void Background::move_into_times(std::vector<BkgFPType>& times)
	{
		m_times = std::move(times);

		// Store min/max values
		auto [min_it, max_it] = std::minmax_element(m_times.begin(),
			m_times.end());
		m_t_min = *min_it;
		m_t_max = *max_it;
	}

	void Background::move_into_zones(Zones::Zones& zones)
	{
		zones.validate();
		m_zones = std::move(zones);
	}
	///@}

	/**
	* @name Plasma and field setters
	*
	* Move plasma and field data into the Background with move semantics.
	* The data comes in as a Vector4D for a single zone, which is moved
	* into a ZonedVector without copying. Before each move, set_dims checks
	* that the data has the same dimensions as the zone and any data
	* already loaded, and stops the run if not.
	* See the plasma and field accessors for what each quantity is.
	*/
	///@{
	void Background::move_into_ne(Vectors::Vector4D<BkgFPType>& ne)
	{
		set_dims(ne, "ne");
		m_ne = Vectors::ZonedVector<BkgFPType>(std::move(ne));
	}

	void Background::move_into_te(Vectors::Vector4D<BkgFPType>& te)
	{
		set_dims(te, "te");
		m_te = Vectors::ZonedVector<BkgFPType>(std::move(te));
	}

	void Background::move_into_ti(Vectors::Vector4D<BkgFPType>& ti)
	{
		set_dims(ti, "ti");
		m_ti = Vectors::ZonedVector<BkgFPType>(std::move(ti));
	}

	void Background::move_into_q_par_i(Vectors::Vector4D<BkgFPType>& q_par_i)
	{
		set_dims(q_par_i, "q_par_i");
		m_q_par_i = Vectors::ZonedVector<BkgFPType>(std::move(q_par_i));
	}

	void Background::move_into_q_par_e(Vectors::Vector4D<BkgFPType>& q_par_e)
	{
		set_dims(q_par_e, "q_par_e");
		m_q_par_e = Vectors::ZonedVector<BkgFPType>(std::move(q_par_e));
	}

	void Background::move_into_vp(Vectors::Vector4D<BkgFPType>& vp)
	{
		set_dims(vp, "vp");
		m_vp = Vectors::ZonedVector<BkgFPType>(std::move(vp));
	}

	void Background::move_into_bX(Vectors::Vector4D<BkgFPType>& bX)
	{
		set_dims(bX, "bX");
		m_bX = Vectors::ZonedVector<BkgFPType>(std::move(bX));
	}

	void Background::move_into_bY(Vectors::Vector4D<BkgFPType>& bY)
	{
		set_dims(bY, "bY");
		m_bY = Vectors::ZonedVector<BkgFPType>(std::move(bY));
	}

	void Background::move_into_bZ(Vectors::Vector4D<BkgFPType>& bZ)
	{
		set_dims(bZ, "bZ");
		m_bZ = Vectors::ZonedVector<BkgFPType>(std::move(bZ));
	}

	void Background::move_into_bR(Vectors::Vector4D<BkgFPType>& bR)
	{
		set_dims(bR, "bR");
		m_bR = Vectors::ZonedVector<BkgFPType>(std::move(bR));
	}

	void Background::move_into_bmag(Vectors::Vector4D<BkgFPType>& bmag)
	{
		set_dims(bmag, "bmag");
		m_bmag = Vectors::ZonedVector<BkgFPType>(std::move(bmag));
	}

	void Background::move_into_gradbX(Vectors::Vector4D<BkgFPType>& gradbX)
	{
		set_dims(gradbX, "gradbX");
		m_gradbX = Vectors::ZonedVector<BkgFPType>(std::move(gradbX));
	}

	void Background::move_into_gradbY(Vectors::Vector4D<BkgFPType>& gradbY)
	{
		set_dims(gradbY, "gradbY");
		m_gradbY = Vectors::ZonedVector<BkgFPType>(std::move(gradbY));
	}

	void Background::move_into_gradbZ(Vectors::Vector4D<BkgFPType>& gradbZ)
	{
		set_dims(gradbZ, "gradbZ");
		m_gradbZ = Vectors::ZonedVector<BkgFPType>(std::move(gradbZ));
	}

	void Background::move_into_eX(Vectors::Vector4D<BkgFPType>& eX)
	{
		set_dims(eX, "eX");
		m_eX = Vectors::ZonedVector<BkgFPType>(std::move(eX));
	}

	void Background::move_into_eY(Vectors::Vector4D<BkgFPType>& eY)
	{
		set_dims(eY, "eY");
		m_eY = Vectors::ZonedVector<BkgFPType>(std::move(eY));
	}

	void Background::move_into_eZ(Vectors::Vector4D<BkgFPType>& eZ)
	{
		set_dims(eZ, "eZ");
		m_eZ = Vectors::ZonedVector<BkgFPType>(std::move(eZ));
	}

	void Background::move_into_emag(Vectors::Vector4D<BkgFPType>& emag)
	{
		set_dims(emag, "emag");
		m_emag = Vectors::ZonedVector<BkgFPType>(std::move(emag));
	}

	void Background::move_into_gradeX(Vectors::Vector4D<BkgFPType>& gradeX)
	{
		set_dims(gradeX, "gradeX");
		m_gradeX = Vectors::ZonedVector<BkgFPType>(std::move(gradeX));
	}

	void Background::move_into_gradeY(Vectors::Vector4D<BkgFPType>& gradeY)
	{
		set_dims(gradeY, "gradeY");
		m_gradeY = Vectors::ZonedVector<BkgFPType>(std::move(gradeY));
	}

	void Background::move_into_gradeZ(Vectors::Vector4D<BkgFPType>& gradeZ)
	{
		set_dims(gradeZ, "gradeZ");
		m_gradeZ = Vectors::ZonedVector<BkgFPType>(std::move(gradeZ));
	}

	void Background::move_into_uX(Vectors::Vector4D<BkgFPType>& uX)
	{
		set_dims(uX, "uX");
		m_uX = Vectors::ZonedVector<BkgFPType>(std::move(uX));
	}

	void Background::move_into_uY(Vectors::Vector4D<BkgFPType>& uY)
	{
		set_dims(uY, "uY");
		m_uY = Vectors::ZonedVector<BkgFPType>(std::move(uY));
	}

	void Background::move_into_uZ(Vectors::Vector4D<BkgFPType>& uZ)
	{
		set_dims(uZ, "uZ");
		m_uZ = Vectors::ZonedVector<BkgFPType>(std::move(uZ));
	}
	///@}

	/**
	* @name Geometry setters
	*
	* Move geometry data (physical coordinates, Jacobian, metric
	* coefficients, basis vectors and covariant magnetic field components)
	* for a single zone into the Background with move semantics, without
	* copying. Unlike the plasma and field setters, these do not check
	* dimensions, because set_dims only handles Vector4D data. See the geometry accessors for what each
	* quantity is.
	*/
	///@{
	void Background::move_into_X(Vectors::Vector3D<BkgFPType>& X)
	{
		// Need a Vector3D version of set_dims
		//set_dims(X, "X");
		m_X = Vectors::ZonedVector<BkgFPType>(std::move(X));
	}

	void Background::move_into_Y(Vectors::Vector3D<BkgFPType>& Y)
	{
		// Need a Vector3D version of set_dims
		//set_dims(Y, "Y");
		m_Y = Vectors::ZonedVector<BkgFPType>(std::move(Y));
	}

	void Background::move_into_Z(Vectors::Vector3D<BkgFPType>& Z)
	{
		// Need a Vector3D version of set_dims
		//set_dims(Z, "Z");
		m_Z = Vectors::ZonedVector<BkgFPType>(std::move(Z));
	}

	void Background::move_into_grid_X(Vectors::Vector3D<BkgFPType>& grid_X)
	{
		// Need a Vector3D version of set_dims
		//set_dims(grid_X, "grid_X");
		m_grid_X = Vectors::ZonedVector<BkgFPType>(std::move(grid_X));
	}

	void Background::move_into_grid_Y(Vectors::Vector3D<BkgFPType>& grid_Y)
	{
		// Need a Vector3D version of set_dims
		//set_dims(grid_Y, "grid_Y");
		m_grid_Y = Vectors::ZonedVector<BkgFPType>(std::move(grid_Y));
	}

	void Background::move_into_grid_Z(Vectors::Vector3D<BkgFPType>& grid_Z)
	{
		// Need a Vector3D version of set_dims
		//set_dims(grid_Z, "grid_Z");
		m_grid_Z = Vectors::ZonedVector<BkgFPType>(std::move(grid_Z));
	}

	void Background::move_into_J(Vectors::Vector3D<BkgFPType>& J)
	{
		// Need a Vector3D version of set_dims
		//set_dims(J, "J");
		m_J = Vectors::ZonedVector<BkgFPType>(std::move(J));
	}

	void Background::move_into_gij_00(Vectors::Vector3D<BkgFPType>& gij_00)
	{
		// Need a Vector3D version of set_dims
		//set_dims(gij_00, "gij_00");
		m_gij_00 = Vectors::ZonedVector<BkgFPType>(std::move(gij_00));
	}

	void Background::move_into_gij_01(Vectors::Vector3D<BkgFPType>& gij_01)
	{
		// Need a Vector3D version of set_dims
		//set_dims(gij_01, "gij_01");
		m_gij_01 = Vectors::ZonedVector<BkgFPType>(std::move(gij_01));
	}

	void Background::move_into_gij_02(Vectors::Vector3D<BkgFPType>& gij_02)
	{
		// Need a Vector3D version of set_dims
		//set_dims(gij_02, "gij_02");
		m_gij_02 = Vectors::ZonedVector<BkgFPType>(std::move(gij_02));
	}

	void Background::move_into_gij_11(Vectors::Vector3D<BkgFPType>& gij_11)
	{
		// Need a Vector3D version of set_dims
		//set_dims(gij_11, "gij_11");
		m_gij_11 = Vectors::ZonedVector<BkgFPType>(std::move(gij_11));
	}

	void Background::move_into_gij_12(Vectors::Vector3D<BkgFPType>& gij_12)
	{
		// Need a Vector3D version of set_dims
		//set_dims(gij_12, "gij_12");
		m_gij_12 = Vectors::ZonedVector<BkgFPType>(std::move(gij_12));
	}

	void Background::move_into_gij_22(Vectors::Vector3D<BkgFPType>& gij_22)
	{
		// Need a Vector3D version of set_dims
		//set_dims(gij_22, "gij_22");
		m_gij_22 = Vectors::ZonedVector<BkgFPType>(std::move(gij_22));
	}

	void Background::move_into_dxdX(Vectors::Vector3D<BkgFPType>& dxdX)
	{
		m_dxdX = Vectors::ZonedVector<BkgFPType>(std::move(dxdX));
	}

	void Background::move_into_dxdY(Vectors::Vector3D<BkgFPType>& dxdY)
	{
		m_dxdY = Vectors::ZonedVector<BkgFPType>(std::move(dxdY));
	}

	void Background::move_into_dxdZ(Vectors::Vector3D<BkgFPType>& dxdZ)
	{
		m_dxdZ = Vectors::ZonedVector<BkgFPType>(std::move(dxdZ));
	}

	void Background::move_into_dydX(Vectors::Vector3D<BkgFPType>& dydX)
	{
		m_dydX = Vectors::ZonedVector<BkgFPType>(std::move(dydX));
	}

	void Background::move_into_dydY(Vectors::Vector3D<BkgFPType>& dydY)
	{
		m_dydY = Vectors::ZonedVector<BkgFPType>(std::move(dydY));
	}

	void Background::move_into_dydZ(Vectors::Vector3D<BkgFPType>& dydZ)
	{
		m_dydZ = Vectors::ZonedVector<BkgFPType>(std::move(dydZ));
	}

	void Background::move_into_dzdX(Vectors::Vector3D<BkgFPType>& dzdX)
	{
		m_dzdX = Vectors::ZonedVector<BkgFPType>(std::move(dzdX));
	}

	void Background::move_into_dzdY(Vectors::Vector3D<BkgFPType>& dzdY)
	{
		m_dzdY = Vectors::ZonedVector<BkgFPType>(std::move(dzdY));
	}

	void Background::move_into_dzdZ(Vectors::Vector3D<BkgFPType>& dzdZ)
	{
		m_dzdZ = Vectors::ZonedVector<BkgFPType>(std::move(dzdZ));
	}

	void Background::move_into_dXdx(Vectors::Vector3D<BkgFPType>& dXdx)
	{
		m_dXdx = Vectors::ZonedVector<BkgFPType>(std::move(dXdx));
	}

	void Background::move_into_dYdx(Vectors::Vector3D<BkgFPType>& dYdx)
	{
		m_dYdx = Vectors::ZonedVector<BkgFPType>(std::move(dYdx));
	}

	void Background::move_into_dZdx(Vectors::Vector3D<BkgFPType>& dZdx)
	{
		m_dZdx = Vectors::ZonedVector<BkgFPType>(std::move(dZdx));
	}

	void Background::move_into_dXdy(Vectors::Vector3D<BkgFPType>& dXdy)
	{
		m_dXdy = Vectors::ZonedVector<BkgFPType>(std::move(dXdy));
	}

	void Background::move_into_dYdy(Vectors::Vector3D<BkgFPType>& dYdy)
	{
		m_dYdy = Vectors::ZonedVector<BkgFPType>(std::move(dYdy));
	}

	void Background::move_into_dZdy(Vectors::Vector3D<BkgFPType>& dZdy)
	{
		m_dZdy = Vectors::ZonedVector<BkgFPType>(std::move(dZdy));
	}

	void Background::move_into_dXdz(Vectors::Vector3D<BkgFPType>& dXdz)
	{
		m_dXdz = Vectors::ZonedVector<BkgFPType>(std::move(dXdz));
	}

	void Background::move_into_dYdz(Vectors::Vector3D<BkgFPType>& dYdz)
	{
		m_dYdz = Vectors::ZonedVector<BkgFPType>(std::move(dYdz));
	}

	void Background::move_into_dZdz(Vectors::Vector3D<BkgFPType>& dZdz)
	{
		m_dZdz = Vectors::ZonedVector<BkgFPType>(std::move(dZdz));
	}

	void Background::move_into_b_x(Vectors::Vector3D<BkgFPType>& b_x)
	{
		m_b_x = Vectors::ZonedVector<BkgFPType>(std::move(b_x));
	}

	void Background::move_into_b_y(Vectors::Vector3D<BkgFPType>& b_y)
	{
		m_b_y = Vectors::ZonedVector<BkgFPType>(std::move(b_y));
	}

	void Background::move_into_b_z(Vectors::Vector3D<BkgFPType>& b_z)
	{
		m_b_z = Vectors::ZonedVector<BkgFPType>(std::move(b_z));
	}
	///@}

	/**
	* @name Setters for data with any number of zones
	*
	* Move plasma, field and geometry data that is already a ZonedVector,
	* with any number of zones and ghost cells, into the Background with
	* move semantics. Before each move, set_dims checks that the data
	* matches the zones and any data already loaded, and stops the run if
	* not. See the accessors for what each quantity is.
	*/
	///@{
	void Background::move_into_ne(Vectors::ZonedVector<BkgFPType>& ne)
	{
		set_dims(ne, "ne", true);
		m_ne = std::move(ne);
	}

	void Background::move_into_te(Vectors::ZonedVector<BkgFPType>& te)
	{
		set_dims(te, "te", true);
		m_te = std::move(te);
	}

	void Background::move_into_ti(Vectors::ZonedVector<BkgFPType>& ti)
	{
		set_dims(ti, "ti", true);
		m_ti = std::move(ti);
	}

	void Background::move_into_q_par_i(
		Vectors::ZonedVector<BkgFPType>& q_par_i)
	{
		set_dims(q_par_i, "q_par_i", true);
		m_q_par_i = std::move(q_par_i);
	}

	void Background::move_into_q_par_e(
		Vectors::ZonedVector<BkgFPType>& q_par_e)
	{
		set_dims(q_par_e, "q_par_e", true);
		m_q_par_e = std::move(q_par_e);
	}

	void Background::move_into_vp(Vectors::ZonedVector<BkgFPType>& vp)
	{
		set_dims(vp, "vp", true);
		m_vp = std::move(vp);
	}

	void Background::move_into_bX(Vectors::ZonedVector<BkgFPType>& bX)
	{
		set_dims(bX, "bX", true);
		m_bX = std::move(bX);
	}

	void Background::move_into_bY(Vectors::ZonedVector<BkgFPType>& bY)
	{
		set_dims(bY, "bY", true);
		m_bY = std::move(bY);
	}

	void Background::move_into_bZ(Vectors::ZonedVector<BkgFPType>& bZ)
	{
		set_dims(bZ, "bZ", true);
		m_bZ = std::move(bZ);
	}

	void Background::move_into_bR(Vectors::ZonedVector<BkgFPType>& bR)
	{
		set_dims(bR, "bR", true);
		m_bR = std::move(bR);
	}

	void Background::move_into_bmag(Vectors::ZonedVector<BkgFPType>& bmag)
	{
		set_dims(bmag, "bmag", true);
		m_bmag = std::move(bmag);
	}

	void Background::move_into_gradbX(Vectors::ZonedVector<BkgFPType>& gradbX)
	{
		set_dims(gradbX, "gradbX", true);
		m_gradbX = std::move(gradbX);
	}

	void Background::move_into_gradbY(Vectors::ZonedVector<BkgFPType>& gradbY)
	{
		set_dims(gradbY, "gradbY", true);
		m_gradbY = std::move(gradbY);
	}

	void Background::move_into_gradbZ(Vectors::ZonedVector<BkgFPType>& gradbZ)
	{
		set_dims(gradbZ, "gradbZ", true);
		m_gradbZ = std::move(gradbZ);
	}

	void Background::move_into_eX(Vectors::ZonedVector<BkgFPType>& eX)
	{
		set_dims(eX, "eX", true);
		m_eX = std::move(eX);
	}

	void Background::move_into_eY(Vectors::ZonedVector<BkgFPType>& eY)
	{
		set_dims(eY, "eY", true);
		m_eY = std::move(eY);
	}

	void Background::move_into_eZ(Vectors::ZonedVector<BkgFPType>& eZ)
	{
		set_dims(eZ, "eZ", true);
		m_eZ = std::move(eZ);
	}

	void Background::move_into_emag(Vectors::ZonedVector<BkgFPType>& emag)
	{
		set_dims(emag, "emag", true);
		m_emag = std::move(emag);
	}

	void Background::move_into_gradeX(Vectors::ZonedVector<BkgFPType>& gradeX)
	{
		set_dims(gradeX, "gradeX", true);
		m_gradeX = std::move(gradeX);
	}

	void Background::move_into_gradeY(Vectors::ZonedVector<BkgFPType>& gradeY)
	{
		set_dims(gradeY, "gradeY", true);
		m_gradeY = std::move(gradeY);
	}

	void Background::move_into_gradeZ(Vectors::ZonedVector<BkgFPType>& gradeZ)
	{
		set_dims(gradeZ, "gradeZ", true);
		m_gradeZ = std::move(gradeZ);
	}

	void Background::move_into_uX(Vectors::ZonedVector<BkgFPType>& uX)
	{
		set_dims(uX, "uX", true);
		m_uX = std::move(uX);
	}

	void Background::move_into_uY(Vectors::ZonedVector<BkgFPType>& uY)
	{
		set_dims(uY, "uY", true);
		m_uY = std::move(uY);
	}

	void Background::move_into_uZ(Vectors::ZonedVector<BkgFPType>& uZ)
	{
		set_dims(uZ, "uZ", true);
		m_uZ = std::move(uZ);
	}

	void Background::move_into_X(Vectors::ZonedVector<BkgFPType>& X)
	{
		set_dims(X, "X", false);
		m_X = std::move(X);
	}

	void Background::move_into_Y(Vectors::ZonedVector<BkgFPType>& Y)
	{
		set_dims(Y, "Y", false);
		m_Y = std::move(Y);
	}

	void Background::move_into_Z(Vectors::ZonedVector<BkgFPType>& Z)
	{
		set_dims(Z, "Z", false);
		m_Z = std::move(Z);
	}

	void Background::move_into_J(Vectors::ZonedVector<BkgFPType>& J)
	{
		set_dims(J, "J", false);
		m_J = std::move(J);
	}

	void Background::move_into_gij_00(Vectors::ZonedVector<BkgFPType>& gij_00)
	{
		set_dims(gij_00, "gij_00", false);
		m_gij_00 = std::move(gij_00);
	}

	void Background::move_into_gij_01(Vectors::ZonedVector<BkgFPType>& gij_01)
	{
		set_dims(gij_01, "gij_01", false);
		m_gij_01 = std::move(gij_01);
	}

	void Background::move_into_gij_02(Vectors::ZonedVector<BkgFPType>& gij_02)
	{
		set_dims(gij_02, "gij_02", false);
		m_gij_02 = std::move(gij_02);
	}

	void Background::move_into_gij_11(Vectors::ZonedVector<BkgFPType>& gij_11)
	{
		set_dims(gij_11, "gij_11", false);
		m_gij_11 = std::move(gij_11);
	}

	void Background::move_into_gij_12(Vectors::ZonedVector<BkgFPType>& gij_12)
	{
		set_dims(gij_12, "gij_12", false);
		m_gij_12 = std::move(gij_12);
	}

	void Background::move_into_gij_22(Vectors::ZonedVector<BkgFPType>& gij_22)
	{
		set_dims(gij_22, "gij_22", false);
		m_gij_22 = std::move(gij_22);
	}

	void Background::move_into_dxdX(Vectors::ZonedVector<BkgFPType>& dxdX)
	{
		set_dims(dxdX, "dxdX", false);
		m_dxdX = std::move(dxdX);
	}

	void Background::move_into_dxdY(Vectors::ZonedVector<BkgFPType>& dxdY)
	{
		set_dims(dxdY, "dxdY", false);
		m_dxdY = std::move(dxdY);
	}

	void Background::move_into_dxdZ(Vectors::ZonedVector<BkgFPType>& dxdZ)
	{
		set_dims(dxdZ, "dxdZ", false);
		m_dxdZ = std::move(dxdZ);
	}

	void Background::move_into_dydX(Vectors::ZonedVector<BkgFPType>& dydX)
	{
		set_dims(dydX, "dydX", false);
		m_dydX = std::move(dydX);
	}

	void Background::move_into_dydY(Vectors::ZonedVector<BkgFPType>& dydY)
	{
		set_dims(dydY, "dydY", false);
		m_dydY = std::move(dydY);
	}

	void Background::move_into_dydZ(Vectors::ZonedVector<BkgFPType>& dydZ)
	{
		set_dims(dydZ, "dydZ", false);
		m_dydZ = std::move(dydZ);
	}

	void Background::move_into_dzdX(Vectors::ZonedVector<BkgFPType>& dzdX)
	{
		set_dims(dzdX, "dzdX", false);
		m_dzdX = std::move(dzdX);
	}

	void Background::move_into_dzdY(Vectors::ZonedVector<BkgFPType>& dzdY)
	{
		set_dims(dzdY, "dzdY", false);
		m_dzdY = std::move(dzdY);
	}

	void Background::move_into_dzdZ(Vectors::ZonedVector<BkgFPType>& dzdZ)
	{
		set_dims(dzdZ, "dzdZ", false);
		m_dzdZ = std::move(dzdZ);
	}

	void Background::move_into_dXdx(Vectors::ZonedVector<BkgFPType>& dXdx)
	{
		set_dims(dXdx, "dXdx", false);
		m_dXdx = std::move(dXdx);
	}

	void Background::move_into_dYdx(Vectors::ZonedVector<BkgFPType>& dYdx)
	{
		set_dims(dYdx, "dYdx", false);
		m_dYdx = std::move(dYdx);
	}

	void Background::move_into_dZdx(Vectors::ZonedVector<BkgFPType>& dZdx)
	{
		set_dims(dZdx, "dZdx", false);
		m_dZdx = std::move(dZdx);
	}

	void Background::move_into_dXdy(Vectors::ZonedVector<BkgFPType>& dXdy)
	{
		set_dims(dXdy, "dXdy", false);
		m_dXdy = std::move(dXdy);
	}

	void Background::move_into_dYdy(Vectors::ZonedVector<BkgFPType>& dYdy)
	{
		set_dims(dYdy, "dYdy", false);
		m_dYdy = std::move(dYdy);
	}

	void Background::move_into_dZdy(Vectors::ZonedVector<BkgFPType>& dZdy)
	{
		set_dims(dZdy, "dZdy", false);
		m_dZdy = std::move(dZdy);
	}

	void Background::move_into_dXdz(Vectors::ZonedVector<BkgFPType>& dXdz)
	{
		set_dims(dXdz, "dXdz", false);
		m_dXdz = std::move(dXdz);
	}

	void Background::move_into_dYdz(Vectors::ZonedVector<BkgFPType>& dYdz)
	{
		set_dims(dYdz, "dYdz", false);
		m_dYdz = std::move(dYdz);
	}

	void Background::move_into_dZdz(Vectors::ZonedVector<BkgFPType>& dZdz)
	{
		set_dims(dZdz, "dZdz", false);
		m_dZdz = std::move(dZdz);
	}

	void Background::move_into_b_x(Vectors::ZonedVector<BkgFPType>& b_x)
	{
		set_dims(b_x, "b_x", false);
		m_b_x = std::move(b_x);
	}

	void Background::move_into_b_y(Vectors::ZonedVector<BkgFPType>& b_y)
	{
		set_dims(b_y, "b_y", false);
		m_b_y = std::move(b_y);
	}

	void Background::move_into_b_z(Vectors::ZonedVector<BkgFPType>& b_z)
	{
		set_dims(b_z, "b_z", false);
		m_b_z = std::move(b_z);
	}

	void Background::move_into_chi(Vectors::ZonedVector<BkgFPType>& chi)
	{
		set_dims(chi, "chi", false);
		m_chi = std::move(chi);
	}

	void Background::move_into_grid_X(Vectors::ZonedVector<BkgFPType>& grid_X)
	{
		set_dims(grid_X, "grid_X", false, true);
		m_grid_X = std::move(grid_X);
	}

	void Background::move_into_grid_Y(Vectors::ZonedVector<BkgFPType>& grid_Y)
	{
		set_dims(grid_Y, "grid_Y", false, true);
		m_grid_Y = std::move(grid_Y);
	}

	void Background::move_into_grid_Z(Vectors::ZonedVector<BkgFPType>& grid_Z)
	{
		set_dims(grid_Z, "grid_Z", false, true);
		m_grid_Z = std::move(grid_Z);
	}
	///@}

	/**
	* Interpolates 4D background data at a given time and position in a
	* zone.
	*
	* This is a convenience wrapper around Utilities::interp_vec4d that
	* passes in the background times and the zones, so callers only need
	* to supply the data and the point to evaluate it at.
	*
	* @param vec Data to interpolate, stored as a function of
	*            (zone, t, x, y, z).
	* @param zone Zone the point is in.
	* @param t0 Time at which to evaluate the data.
	* @param x0 x-coordinate at which to evaluate the data.
	* @param y0 y-coordinate at which to evaluate the data.
	* @param z0 z-coordinate at which to evaluate the data.
	* @return Interpolated value at (t0, x0, y0, z0).
	*/
	BkgFPType Background::interp_4d(
		const Vectors::ZonedVector<BkgFPType>& vec, const int zone,
		const double t0, const double x0, const double y0,
		const double z0) const
	{
		// This is really just a convenient wrapper to handle passing in
		// times and the x, y, z vectors.
		return Utilities::interp_vec4d(vec, zone, m_zones, m_times, t0, x0,
			y0, z0);
	}

	/**
	* @name Interpolation at a point
	*
	* Interpolate one plasma quantity at the time and position (t0, x0,
	* y0, z0) in a zone using interp_4d. The quantities are the electron
	* temperature (te), ion temperature (ti), electron density (ne) and the
	* background flow components (uX, uY, uZ).
	*/
	///@{
	BkgFPType Background::interp_te(const int zone, const double t0,
		const double x0, const double y0, const double z0) const
	{
		return interp_4d(m_te, zone, t0, x0, y0, z0);
	}

	BkgFPType Background::interp_ti(const int zone, const double t0,
		const double x0, const double y0, const double z0) const
	{
		return interp_4d(m_ti, zone, t0, x0, y0, z0);
	}

	BkgFPType Background::interp_ne(const int zone, const double t0,
		const double x0, const double y0, const double z0) const
	{
		return interp_4d(m_ne, zone, t0, x0, y0, z0);
	}

	/**
	* Interpolates the ion parallel conductive heat flux q_par (W/m^2,
	* positive along B) at a point, like interp_ti. A background without a
	* heat flux (see get_q_par_i) gives 0, i.e., no skew of the ion
	* distribution.
	*/
	BkgFPType Background::interp_q_par_i(const int zone, const double t0,
		const double x0, const double y0, const double z0) const
	{
		if (m_q_par_i.get_data().empty()) return 0.0;
		return interp_4d(m_q_par_i, zone, t0, x0, y0, z0);
	}

	/**
	* Interpolates the electron parallel conductive heat flux, like
	* interp_q_par_i. A background without one gives 0.
	*/
	BkgFPType Background::interp_q_par_e(const int zone, const double t0,
		const double x0, const double y0, const double z0) const
	{
		if (m_q_par_e.get_data().empty()) return 0.0;
		return interp_4d(m_q_par_e, zone, t0, x0, y0, z0);
	}

	BkgFPType Background::interp_uX(const int zone, const double t0,
		const double x0, const double y0, const double z0) const
	{
		return interp_4d(m_uX, zone, t0, x0, y0, z0);
	}

	BkgFPType Background::interp_uY(const int zone, const double t0,
		const double x0, const double y0, const double z0) const
	{
		return interp_4d(m_uY, zone, t0, x0, y0, z0);
	}

	BkgFPType Background::interp_uZ(const int zone, const double t0,
		const double x0, const double y0, const double z0) const
	{
		return interp_4d(m_uZ, zone, t0, x0, y0, z0);
	}
	///@}

	/**
	* Broadcasts the background from rank 0 to every other MPI rank.
	*
	* Every rank must call this. After it returns, every rank holds a copy
	* of the number of times, the time bounds, the species masses, the
	* times, the zones (with
	* their computational grids), the geometry data and the plasma and
	* field data loaded on rank 0. Not broadcast are bR, the covariant
	* magnetic field components (b_x, b_y, b_z) and the electric field
	* gradient (gradeX, gradeY, gradeZ), which are only saved to the output
	* file on rank 0.
	*
	* @param comm MPI communicator to broadcast over.
	*/
	void Background::broadcast(MPI_Comm comm)
	{
		// Broadcast scalars
		MPI_Bcast(&m_dim1, 1, MPI_INT, 0, comm);

		// Broadcast the time bounds. Rank 0 sets these in move_into_times
		// while loading the background, so the other ranks only get them
		// from here.
		const MPI_Datatype bkg_type {mpi_type<BkgFPType>::type};
		MPI_Bcast(&m_t_min, 1, bkg_type, 0, comm);
		MPI_Bcast(&m_t_max, 1, bkg_type, 0, comm);

		// Species masses
		MPI_Bcast(&m_elec_mass_amu, 1, MPI_DOUBLE, 0, comm);
		MPI_Bcast(&m_ion_mass_amu, 1, MPI_DOUBLE, 0, comm);

		// Broadcast the times and the zones, which hold the computational
		// grids
		Utilities::mpi_broadcast_vector(m_times, 0, comm);
		m_zones.broadcast(comm);

		// Broadcast geometry
		m_X.broadcast(comm);
		m_Y.broadcast(comm);
		m_Z.broadcast(comm);
		m_grid_X.broadcast(comm);
		m_grid_Y.broadcast(comm);
		m_grid_Z.broadcast(comm);
		m_J.broadcast(comm);
		m_gij_00.broadcast(comm);
		m_gij_01.broadcast(comm);
		m_gij_02.broadcast(comm);
		m_gij_11.broadcast(comm);
		m_gij_12.broadcast(comm);
		m_gij_22.broadcast(comm);
		m_dxdX.broadcast(comm);
		m_dxdY.broadcast(comm);
		m_dxdZ.broadcast(comm);
		m_dydX.broadcast(comm);
		m_dydY.broadcast(comm);
		m_dydZ.broadcast(comm);
		m_dzdX.broadcast(comm);
		m_dzdY.broadcast(comm);
		m_dzdZ.broadcast(comm);
		m_dXdx.broadcast(comm);
		m_dYdx.broadcast(comm);
		m_dZdx.broadcast(comm);
		m_dXdy.broadcast(comm);
		m_dYdy.broadcast(comm);
		m_dZdy.broadcast(comm);
		m_dXdz.broadcast(comm);
		m_dYdz.broadcast(comm);
		m_dZdz.broadcast(comm);
		m_chi.broadcast(comm);

		// Wall outlines: how many, then each one's points
		int nwalls {static_cast<int>(m_walls.size())};
		MPI_Bcast(&nwalls, 1, MPI_INT, 0, comm);
		m_walls.resize(nwalls);
		for (auto& wall : m_walls)
		{
			Utilities::mpi_broadcast_vector(wall.R, 0, comm);
			Utilities::mpi_broadcast_vector(wall.Z, 0, comm);
		}

		// Broadcast plasma and fields
		m_ne.broadcast(comm);
		m_te.broadcast(comm);
		m_ti.broadcast(comm);
		m_q_par_i.broadcast(comm);
		m_q_par_e.broadcast(comm);
		m_vp.broadcast(comm);
		m_bX.broadcast(comm);
		m_bY.broadcast(comm);
		m_bZ.broadcast(comm);
		m_bmag.broadcast(comm);
		m_gradbX.broadcast(comm);
		m_gradbY.broadcast(comm);
		m_gradbZ.broadcast(comm);
		m_eX.broadcast(comm);
		m_eY.broadcast(comm);
		m_eZ.broadcast(comm);
		m_emag.broadcast(comm);
		m_uX.broadcast(comm);
		m_uY.broadcast(comm);
		m_uZ.broadcast(comm);
	}

#ifndef USE_CUDA

	// If these definitions gets called, it means you're in the CPU-only version
	// of the code trying to call GPU-specific code. The real definitions are
	// in src/cuda/background.cu since they call CUDA code.

	/**
	* Placeholder for copying the background to a GPU in CPU-only builds.
	*
	* The real definition is in src/cuda/background.cu. This version stops
	* the run with Errors::fatal_error, since GPU support was not compiled
	* in.
	*
	* @param device_id ID of the GPU to copy to (unused).
	*/
	BackgroundDevice Background::to_device(int device_id)
	{
		Errors::fatal_error("Background::to_device was called but GPU "
			"support was not compiled in.");
	}

	/**
	* Placeholder for freeing a device-side background in CPU-only builds.
	*
	* The real definition is in src/cuda/background.cu. This version stops
	* the run with Errors::fatal_error, since GPU support was not compiled
	* in.
	*
	* @param bkg_d Device-side background to free (unused).
	* @param device_id ID of the GPU the background is on (unused).
	*/
	void free_bkg(BackgroundDevice& bkg_d, int device_id)
	{
		Errors::fatal_error("Background::free_bkg was called but GPU "
			"support was not compiled in.");
	}

#endif

}  // namespace Background

// This is a peculiarity of splitting class declarations and definitions
// between a header and source file. We need to instantiate the templates with
// the needed definitions so the linker can see them. Hurts flexibility, but
// not an issue since it's pretty straightforward.
template void Background::Background::set_dims<int>(
	const Vectors::Vector4D<int>& v, const std::string_view data);
template void Background::Background::set_dims<BkgFPType>(
	const Vectors::Vector4D<BkgFPType>& v, const std::string_view data);
