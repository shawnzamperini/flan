/**
* @file background.h
* @brief Header file for background.cpp
*/
#pragma once

#include <string>
#include <string_view>
#include <vector>

#include "background_device.h"
#include "flan_types.h"
#include "vectors.h"
#include "zones.h"

#ifndef __CUDACC__
#include "mpi.h"
#endif

namespace Background
{
	/**
	* @brief Contains time-dependent background plasma for simulation
	*
	* The simulation volume is split into one or more zones (see
	* Zones::Zones), each with its own computational grid. Every plasma,
	* field and geometry array is stored per zone in a ZonedVector, indexed
	* as (zone, t, x, y, z) or (zone, x, y, z).
	*/
	class Background
	{
	private:
		std::vector<BkgFPType> m_times {};
		Zones::Zones m_zones {};
		Vectors::ZonedVector<BkgFPType> m_ne {};
		Vectors::ZonedVector<BkgFPType> m_te {};
		Vectors::ZonedVector<BkgFPType> m_ti {};
		Vectors::ZonedVector<BkgFPType> m_vp {};
		Vectors::ZonedVector<BkgFPType> m_bX {};
		Vectors::ZonedVector<BkgFPType> m_bY {};
		Vectors::ZonedVector<BkgFPType> m_bZ {};
		Vectors::ZonedVector<BkgFPType> m_bmag {};
		Vectors::ZonedVector<BkgFPType> m_gradbX {};
		Vectors::ZonedVector<BkgFPType> m_gradbY {};
		Vectors::ZonedVector<BkgFPType> m_gradbZ {};
		Vectors::ZonedVector<BkgFPType> m_eX {};
		Vectors::ZonedVector<BkgFPType> m_eY {};
		Vectors::ZonedVector<BkgFPType> m_eZ {};
		Vectors::ZonedVector<BkgFPType> m_emag {};
		Vectors::ZonedVector<BkgFPType> m_uX {};
		Vectors::ZonedVector<BkgFPType> m_uY {};
		Vectors::ZonedVector<BkgFPType> m_uZ {};
		Vectors::ZonedVector<BkgFPType> m_X {};
		Vectors::ZonedVector<BkgFPType> m_Y {};
		Vectors::ZonedVector<BkgFPType> m_Z {};
		Vectors::ZonedVector<BkgFPType> m_grid_X {};
		Vectors::ZonedVector<BkgFPType> m_grid_Y {};
		Vectors::ZonedVector<BkgFPType> m_grid_Z {};
		Vectors::ZonedVector<BkgFPType> m_J {};
		Vectors::ZonedVector<BkgFPType> m_gij_00 {};
		Vectors::ZonedVector<BkgFPType> m_gij_01 {};
		Vectors::ZonedVector<BkgFPType> m_gij_02 {};
		Vectors::ZonedVector<BkgFPType> m_gij_11 {};
		Vectors::ZonedVector<BkgFPType> m_gij_12 {};
		Vectors::ZonedVector<BkgFPType> m_gij_22 {};
		Vectors::ZonedVector<BkgFPType> m_dxdX {};
		Vectors::ZonedVector<BkgFPType> m_dxdY {};
		Vectors::ZonedVector<BkgFPType> m_dxdZ {};
		Vectors::ZonedVector<BkgFPType> m_dydX {};
		Vectors::ZonedVector<BkgFPType> m_dydY {};
		Vectors::ZonedVector<BkgFPType> m_dydZ {};
		Vectors::ZonedVector<BkgFPType> m_dzdX {};
		Vectors::ZonedVector<BkgFPType> m_dzdY {};
		Vectors::ZonedVector<BkgFPType> m_dzdZ {};
		Vectors::ZonedVector<BkgFPType> m_dXdx {};
		Vectors::ZonedVector<BkgFPType> m_dYdx {};
		Vectors::ZonedVector<BkgFPType> m_dZdx {};
		Vectors::ZonedVector<BkgFPType> m_dXdy {};
		Vectors::ZonedVector<BkgFPType> m_dYdy {};
		Vectors::ZonedVector<BkgFPType> m_dZdy {};
		Vectors::ZonedVector<BkgFPType> m_dXdz {};
		Vectors::ZonedVector<BkgFPType> m_dYdz {};
		Vectors::ZonedVector<BkgFPType> m_dZdz {};

		// The following are just for debugging and post analysis
		Vectors::ZonedVector<BkgFPType> m_bR {};
		Vectors::ZonedVector<BkgFPType> m_b_x {};
		Vectors::ZonedVector<BkgFPType> m_b_y {};
		Vectors::ZonedVector<BkgFPType> m_b_z {};
		Vectors::ZonedVector<BkgFPType> m_gradeX {};
		Vectors::ZonedVector<BkgFPType> m_gradeY {};
		Vectors::ZonedVector<BkgFPType> m_gradeZ {};

		int m_dim1 {};

		BkgFPType m_t_min {};
		BkgFPType m_t_max {};

		// Masses of the background electrons and main ions (amu)
		double m_elec_mass_amu {};
		double m_ion_mass_amu {};

		void check_single_zone(const std::string_view caller) const;

	public:

		// Accessors
		const std::vector<BkgFPType>& get_times() const;
		const Zones::Zones& get_zones() const;
		int get_nzones() const;
		const std::vector<BkgFPType>& get_x(const int zone) const;
		const std::vector<BkgFPType>& get_y(const int zone) const;
		const std::vector<BkgFPType>& get_z(const int zone) const;
		const std::vector<BkgFPType>& get_grid_x(const int zone) const;
		const std::vector<BkgFPType>& get_grid_y(const int zone) const;
		const std::vector<BkgFPType>& get_grid_z(const int zone) const;
		const Vectors::ZonedVector<BkgFPType>& get_ne() const;
		Vectors::ZonedVector<BkgFPType>& get_ne();
		const Vectors::ZonedVector<BkgFPType>& get_te() const;
		Vectors::ZonedVector<BkgFPType>& get_te();
		const Vectors::ZonedVector<BkgFPType>& get_ti() const;
		Vectors::ZonedVector<BkgFPType>& get_ti();
		const Vectors::ZonedVector<BkgFPType>& get_vp() const;
		Vectors::ZonedVector<BkgFPType>& get_vp();
		const Vectors::ZonedVector<BkgFPType>& get_bX() const;
		Vectors::ZonedVector<BkgFPType>& get_bX();
		const Vectors::ZonedVector<BkgFPType>& get_bY() const;
		Vectors::ZonedVector<BkgFPType>& get_bY();
		const Vectors::ZonedVector<BkgFPType>& get_bZ() const;
		Vectors::ZonedVector<BkgFPType>& get_bZ();
		const Vectors::ZonedVector<BkgFPType>& get_bR() const;
		Vectors::ZonedVector<BkgFPType>& get_bR();
		const Vectors::ZonedVector<BkgFPType>& get_bmag() const;
		Vectors::ZonedVector<BkgFPType>& get_bmag();
		const Vectors::ZonedVector<BkgFPType>& get_gradbX() const;
		const Vectors::ZonedVector<BkgFPType>& get_gradbY() const;
		const Vectors::ZonedVector<BkgFPType>& get_gradbZ() const;
		const Vectors::ZonedVector<BkgFPType>& get_eX() const;
		const Vectors::ZonedVector<BkgFPType>& get_eY() const;
		const Vectors::ZonedVector<BkgFPType>& get_eZ() const;
		const Vectors::ZonedVector<BkgFPType>& get_emag() const;
		Vectors::ZonedVector<BkgFPType>& get_eX();
		Vectors::ZonedVector<BkgFPType>& get_eY();
		Vectors::ZonedVector<BkgFPType>& get_eZ();
		Vectors::ZonedVector<BkgFPType>& get_emag();
		const Vectors::ZonedVector<BkgFPType>& get_gradeX() const;
		const Vectors::ZonedVector<BkgFPType>& get_gradeY() const;
		const Vectors::ZonedVector<BkgFPType>& get_gradeZ() const;
		const Vectors::ZonedVector<BkgFPType>& get_uX() const;
		const Vectors::ZonedVector<BkgFPType>& get_uY() const;
		const Vectors::ZonedVector<BkgFPType>& get_uZ() const;
		Vectors::ZonedVector<BkgFPType>& get_uX();
		Vectors::ZonedVector<BkgFPType>& get_uY();
		Vectors::ZonedVector<BkgFPType>& get_uZ();
		const Vectors::ZonedVector<BkgFPType>& get_X() const;
		const Vectors::ZonedVector<BkgFPType>& get_Y() const;
		const Vectors::ZonedVector<BkgFPType>& get_Z() const;
		const Vectors::ZonedVector<BkgFPType>& get_grid_X() const;
		const Vectors::ZonedVector<BkgFPType>& get_grid_Y() const;
		const Vectors::ZonedVector<BkgFPType>& get_grid_Z() const;
		const Vectors::ZonedVector<BkgFPType>& get_J() const;
		const Vectors::ZonedVector<BkgFPType>& get_gij_00() const;
		const Vectors::ZonedVector<BkgFPType>& get_gij_01() const;
		const Vectors::ZonedVector<BkgFPType>& get_gij_02() const;
		const Vectors::ZonedVector<BkgFPType>& get_gij_11() const;
		const Vectors::ZonedVector<BkgFPType>& get_gij_12() const;
		const Vectors::ZonedVector<BkgFPType>& get_gij_22() const;
		const Vectors::ZonedVector<BkgFPType>& get_dxdX() const;
		const Vectors::ZonedVector<BkgFPType>& get_dxdY() const;
		const Vectors::ZonedVector<BkgFPType>& get_dxdZ() const;
		const Vectors::ZonedVector<BkgFPType>& get_dydX() const;
		const Vectors::ZonedVector<BkgFPType>& get_dydY() const;
		const Vectors::ZonedVector<BkgFPType>& get_dydZ() const;
		const Vectors::ZonedVector<BkgFPType>& get_dzdX() const;
		const Vectors::ZonedVector<BkgFPType>& get_dzdY() const;
		const Vectors::ZonedVector<BkgFPType>& get_dzdZ() const;
		Vectors::ZonedVector<BkgFPType>& get_dxdX();
		Vectors::ZonedVector<BkgFPType>& get_dxdY();
		Vectors::ZonedVector<BkgFPType>& get_dxdZ();
		Vectors::ZonedVector<BkgFPType>& get_dydX();
		Vectors::ZonedVector<BkgFPType>& get_dydY();
		Vectors::ZonedVector<BkgFPType>& get_dydZ();
		Vectors::ZonedVector<BkgFPType>& get_dzdX();
		Vectors::ZonedVector<BkgFPType>& get_dzdY();
		Vectors::ZonedVector<BkgFPType>& get_dzdZ();
		const Vectors::ZonedVector<BkgFPType>& get_dXdx() const;
		const Vectors::ZonedVector<BkgFPType>& get_dYdx() const;
		const Vectors::ZonedVector<BkgFPType>& get_dZdx() const;
		const Vectors::ZonedVector<BkgFPType>& get_dXdy() const;
		const Vectors::ZonedVector<BkgFPType>& get_dYdy() const;
		const Vectors::ZonedVector<BkgFPType>& get_dZdy() const;
		const Vectors::ZonedVector<BkgFPType>& get_dXdz() const;
		const Vectors::ZonedVector<BkgFPType>& get_dYdz() const;
		const Vectors::ZonedVector<BkgFPType>& get_dZdz() const;
		const Vectors::ZonedVector<BkgFPType>& get_b_x() const;
		const Vectors::ZonedVector<BkgFPType>& get_b_y() const;
		const Vectors::ZonedVector<BkgFPType>& get_b_z() const;
		double get_elec_mass_amu() const;
		double get_ion_mass_amu() const;
		void set_elec_mass_amu(const double elec_mass_amu);
		void set_ion_mass_amu(const double ion_mass_amu);
		int get_dim1() const;
		int get_dim2() const;
		int get_dim3() const;
		int get_dim4() const;

		// Psuedo-accessors to get the bound of the simulation spacetime.
		BkgFPType get_t_min() const;
		BkgFPType get_t_max() const;
		BkgFPType get_x_min(const int zone) const;
		BkgFPType get_x_max(const int zone) const;
		BkgFPType get_y_min(const int zone) const;
		BkgFPType get_y_max(const int zone) const;
		BkgFPType get_z_min(const int zone) const;
		BkgFPType get_z_max(const int zone) const;

		// Setters using move semantics. No idea if this is the proper
		// way to do it but it works.
		void move_into_times(std::vector<BkgFPType>& times);
		void move_into_zones(Zones::Zones& zones);
		void move_into_ne(Vectors::Vector4D<BkgFPType>& ne);
		void move_into_te(Vectors::Vector4D<BkgFPType>& te);
		void move_into_ti(Vectors::Vector4D<BkgFPType>& ti);
		void move_into_vp(Vectors::Vector4D<BkgFPType>& vp);
		void move_into_bX(Vectors::Vector4D<BkgFPType>& bX);
		void move_into_bY(Vectors::Vector4D<BkgFPType>& bY);
		void move_into_bZ(Vectors::Vector4D<BkgFPType>& bZ);
		void move_into_bR(Vectors::Vector4D<BkgFPType>& bR);
		void move_into_bmag(Vectors::Vector4D<BkgFPType>& bmag);
		void move_into_gradbX(Vectors::Vector4D<BkgFPType>& gradbX);
		void move_into_gradbY(Vectors::Vector4D<BkgFPType>& gradbY);
		void move_into_gradbZ(Vectors::Vector4D<BkgFPType>& gradbZ);
		void move_into_eX(Vectors::Vector4D<BkgFPType>& eX); 
		void move_into_eY(Vectors::Vector4D<BkgFPType>& eY);
		void move_into_eZ(Vectors::Vector4D<BkgFPType>& eZ);
		void move_into_emag(Vectors::Vector4D<BkgFPType>& emag);
		void move_into_gradeX(Vectors::Vector4D<BkgFPType>& gradeX); 
		void move_into_gradeY(Vectors::Vector4D<BkgFPType>& gradeY);
		void move_into_gradeZ(Vectors::Vector4D<BkgFPType>& gradeZ);
		void move_into_uX(Vectors::Vector4D<BkgFPType>& uX); 
		void move_into_uY(Vectors::Vector4D<BkgFPType>& uY);
		void move_into_uZ(Vectors::Vector4D<BkgFPType>& uZ);
		void move_into_X(Vectors::Vector3D<BkgFPType>& X);
		void move_into_Y(Vectors::Vector3D<BkgFPType>& Y);
		void move_into_Z(Vectors::Vector3D<BkgFPType>& Z);
		void move_into_grid_X(Vectors::Vector3D<BkgFPType>& grid_X);
		void move_into_grid_Y(Vectors::Vector3D<BkgFPType>& grid_Y);
		void move_into_grid_Z(Vectors::Vector3D<BkgFPType>& grid_Z);
		void move_into_J(Vectors::Vector3D<BkgFPType>& J);
		void move_into_gij_00(Vectors::Vector3D<BkgFPType>& gij_00);
		void move_into_gij_01(Vectors::Vector3D<BkgFPType>& gij_01);
		void move_into_gij_02(Vectors::Vector3D<BkgFPType>& gij_02);
		void move_into_gij_11(Vectors::Vector3D<BkgFPType>& gij_11);
		void move_into_gij_12(Vectors::Vector3D<BkgFPType>& gij_12);
		void move_into_gij_22(Vectors::Vector3D<BkgFPType>& gij_22);
		void move_into_dxdX(Vectors::Vector3D<BkgFPType>& dxdX);
		void move_into_dxdY(Vectors::Vector3D<BkgFPType>& dxdY);
		void move_into_dxdZ(Vectors::Vector3D<BkgFPType>& dxdZ);
		void move_into_dydX(Vectors::Vector3D<BkgFPType>& dydX);
		void move_into_dydY(Vectors::Vector3D<BkgFPType>& dydY);
		void move_into_dydZ(Vectors::Vector3D<BkgFPType>& dydZ);
		void move_into_dzdX(Vectors::Vector3D<BkgFPType>& dzdX);
		void move_into_dzdY(Vectors::Vector3D<BkgFPType>& dzdY);
		void move_into_dzdZ(Vectors::Vector3D<BkgFPType>& dzdZ);
		void move_into_dXdx(Vectors::Vector3D<BkgFPType>& dXdx);
		void move_into_dYdx(Vectors::Vector3D<BkgFPType>& dYdx);
		void move_into_dZdx(Vectors::Vector3D<BkgFPType>& dZdx);
		void move_into_dXdy(Vectors::Vector3D<BkgFPType>& dXdy);
		void move_into_dYdy(Vectors::Vector3D<BkgFPType>& dYdy);
		void move_into_dZdy(Vectors::Vector3D<BkgFPType>& dZdy);
		void move_into_dXdz(Vectors::Vector3D<BkgFPType>& dXdz);
		void move_into_dYdz(Vectors::Vector3D<BkgFPType>& dYdz);
		void move_into_dZdz(Vectors::Vector3D<BkgFPType>& dZdz);
		void move_into_b_x(Vectors::Vector3D<BkgFPType>& b_x);
		void move_into_b_y(Vectors::Vector3D<BkgFPType>& b_y);
		void move_into_b_z(Vectors::Vector3D<BkgFPType>& b_z);

		// Setters for data with any number of zones
		void move_into_ne(Vectors::ZonedVector<BkgFPType>& ne);
		void move_into_te(Vectors::ZonedVector<BkgFPType>& te);
		void move_into_ti(Vectors::ZonedVector<BkgFPType>& ti);
		void move_into_vp(Vectors::ZonedVector<BkgFPType>& vp);
		void move_into_bX(Vectors::ZonedVector<BkgFPType>& bX);
		void move_into_bY(Vectors::ZonedVector<BkgFPType>& bY);
		void move_into_bZ(Vectors::ZonedVector<BkgFPType>& bZ);
		void move_into_bR(Vectors::ZonedVector<BkgFPType>& bR);
		void move_into_bmag(Vectors::ZonedVector<BkgFPType>& bmag);
		void move_into_gradbX(Vectors::ZonedVector<BkgFPType>& gradbX);
		void move_into_gradbY(Vectors::ZonedVector<BkgFPType>& gradbY);
		void move_into_gradbZ(Vectors::ZonedVector<BkgFPType>& gradbZ);
		void move_into_eX(Vectors::ZonedVector<BkgFPType>& eX);
		void move_into_eY(Vectors::ZonedVector<BkgFPType>& eY);
		void move_into_eZ(Vectors::ZonedVector<BkgFPType>& eZ);
		void move_into_emag(Vectors::ZonedVector<BkgFPType>& emag);
		void move_into_gradeX(Vectors::ZonedVector<BkgFPType>& gradeX);
		void move_into_gradeY(Vectors::ZonedVector<BkgFPType>& gradeY);
		void move_into_gradeZ(Vectors::ZonedVector<BkgFPType>& gradeZ);
		void move_into_uX(Vectors::ZonedVector<BkgFPType>& uX);
		void move_into_uY(Vectors::ZonedVector<BkgFPType>& uY);
		void move_into_uZ(Vectors::ZonedVector<BkgFPType>& uZ);
		void move_into_X(Vectors::ZonedVector<BkgFPType>& X);
		void move_into_Y(Vectors::ZonedVector<BkgFPType>& Y);
		void move_into_Z(Vectors::ZonedVector<BkgFPType>& Z);
		void move_into_J(Vectors::ZonedVector<BkgFPType>& J);
		void move_into_gij_00(Vectors::ZonedVector<BkgFPType>& gij_00);
		void move_into_gij_01(Vectors::ZonedVector<BkgFPType>& gij_01);
		void move_into_gij_02(Vectors::ZonedVector<BkgFPType>& gij_02);
		void move_into_gij_11(Vectors::ZonedVector<BkgFPType>& gij_11);
		void move_into_gij_12(Vectors::ZonedVector<BkgFPType>& gij_12);
		void move_into_gij_22(Vectors::ZonedVector<BkgFPType>& gij_22);
		void move_into_dxdX(Vectors::ZonedVector<BkgFPType>& dxdX);
		void move_into_dxdY(Vectors::ZonedVector<BkgFPType>& dxdY);
		void move_into_dxdZ(Vectors::ZonedVector<BkgFPType>& dxdZ);
		void move_into_dydX(Vectors::ZonedVector<BkgFPType>& dydX);
		void move_into_dydY(Vectors::ZonedVector<BkgFPType>& dydY);
		void move_into_dydZ(Vectors::ZonedVector<BkgFPType>& dydZ);
		void move_into_dzdX(Vectors::ZonedVector<BkgFPType>& dzdX);
		void move_into_dzdY(Vectors::ZonedVector<BkgFPType>& dzdY);
		void move_into_dzdZ(Vectors::ZonedVector<BkgFPType>& dzdZ);
		void move_into_dXdx(Vectors::ZonedVector<BkgFPType>& dXdx);
		void move_into_dYdx(Vectors::ZonedVector<BkgFPType>& dYdx);
		void move_into_dZdx(Vectors::ZonedVector<BkgFPType>& dZdx);
		void move_into_dXdy(Vectors::ZonedVector<BkgFPType>& dXdy);
		void move_into_dYdy(Vectors::ZonedVector<BkgFPType>& dYdy);
		void move_into_dZdy(Vectors::ZonedVector<BkgFPType>& dZdy);
		void move_into_dXdz(Vectors::ZonedVector<BkgFPType>& dXdz);
		void move_into_dYdz(Vectors::ZonedVector<BkgFPType>& dYdz);
		void move_into_dZdz(Vectors::ZonedVector<BkgFPType>& dZdz);
		void move_into_b_x(Vectors::ZonedVector<BkgFPType>& b_x);
		void move_into_b_y(Vectors::ZonedVector<BkgFPType>& b_y);
		void move_into_b_z(Vectors::ZonedVector<BkgFPType>& b_z);
		void move_into_grid_X(Vectors::ZonedVector<BkgFPType>& grid_X);
		void move_into_grid_Y(Vectors::ZonedVector<BkgFPType>& grid_Y);
		void move_into_grid_Z(Vectors::ZonedVector<BkgFPType>& grid_Z);

		BkgFPType interp_4d(const Vectors::ZonedVector<BkgFPType>& vec,
			const int zone, const double t0, const double x0, const double y0,
			const double z0) const;

		BkgFPType interp_te(const int zone, const double t0, const double x0,
			const double y0, const double z0) const;

		BkgFPType interp_ti(const int zone, const double t0, const double x0,
			const double y0, const double z0) const;

		BkgFPType interp_ne(const int zone, const double t0, const double x0,
			const double y0, const double z0) const;

		BkgFPType interp_uX(const int zone, const double t0, const double x0,
			const double y0, const double z0) const;

		BkgFPType interp_uY(const int zone, const double t0, const double x0,
			const double y0, const double z0) const;

		BkgFPType interp_uZ(const int zone, const double t0, const double x0,
			const double y0, const double z0) const;

		std::string check_dim(const int m_dim, const int in_dim, 
			const int dim_num);

		template <typename T>
		void set_dims(const Vectors::Vector4D<T>& v,
			const std::string_view data);

		void set_dims(const Vectors::ZonedVector<BkgFPType>& v,
			const std::string_view data, const bool time_dependent,
			const bool corners = false);

#ifndef __CUDACC__
		void broadcast(MPI_Comm comm);
#endif

		BackgroundDevice to_device(int device_id = 0);

		BkgFPType interp_te();
	};

	void free_bkg(BackgroundDevice& bkg_d, int device_id = 0);
}
