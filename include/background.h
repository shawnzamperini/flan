/**
* @file background.h
* @brief Header file for background.cpp
*/
#ifndef BACKGROUND_H
#define BACKGROUND_H

#include <vector>

#include "flan_types.h"
#include "vectors.h"

/**
* @namespace Background
* @brief Namespace to contain anything related to the background plasma
*/
namespace Background
{
	/**
	* @brief Contains time-dependent background plasma for simulation
	*/
	class Background
	{
	private:
		std::vector<BkgFPType> m_times {};
		std::vector<BkgFPType> m_x {};
		std::vector<BkgFPType> m_y {};
		std::vector<BkgFPType> m_z {};
		std::vector<BkgFPType> m_grid_x {};
		std::vector<BkgFPType> m_grid_y {};
		std::vector<BkgFPType> m_grid_z {};
		Vectors::Vector4D<BkgFPType> m_ne {};
		Vectors::Vector4D<BkgFPType> m_te {};
		Vectors::Vector4D<BkgFPType> m_ti {};
		Vectors::Vector4D<BkgFPType> m_vp {};
		Vectors::Vector4D<BkgFPType> m_bX {};
		Vectors::Vector4D<BkgFPType> m_bY {};
		Vectors::Vector4D<BkgFPType> m_bZ {};
		Vectors::Vector4D<BkgFPType> m_gradbX {};
		Vectors::Vector4D<BkgFPType> m_gradbY {};
		Vectors::Vector4D<BkgFPType> m_gradbZ {};
		Vectors::Vector4D<BkgFPType> m_eX {};
		Vectors::Vector4D<BkgFPType> m_eY {};
		Vectors::Vector4D<BkgFPType> m_eZ {};
		Vectors::Vector4D<BkgFPType> m_uX {};
		Vectors::Vector4D<BkgFPType> m_uY {};
		Vectors::Vector4D<BkgFPType> m_uZ {};
		Vectors::Vector3D<BkgFPType> m_X {};
		Vectors::Vector3D<BkgFPType> m_Y {};
		Vectors::Vector3D<BkgFPType> m_Z {};
		Vectors::Vector3D<BkgFPType> m_grid_X {};
		Vectors::Vector3D<BkgFPType> m_grid_Y {};
		Vectors::Vector3D<BkgFPType> m_grid_Z {};
		Vectors::Vector3D<BkgFPType> m_J {};
		Vectors::Vector3D<BkgFPType> m_gij_00 {};
		Vectors::Vector3D<BkgFPType> m_gij_01 {};
		Vectors::Vector3D<BkgFPType> m_gij_02 {};
		Vectors::Vector3D<BkgFPType> m_gij_11 {};
		Vectors::Vector3D<BkgFPType> m_gij_12 {};
		Vectors::Vector3D<BkgFPType> m_gij_22 {};

		// The following are just for debugging
		Vectors::Vector4D<BkgFPType> m_bR {};
		Vectors::Vector3D<BkgFPType> m_b_x {};
		Vectors::Vector3D<BkgFPType> m_b_y {};
		Vectors::Vector3D<BkgFPType> m_b_z {};

		int m_dim1 {};
		int m_dim2 {};
		int m_dim3 {};
		int m_dim4 {};

	public:

		// Accessors
		const std::vector<BkgFPType>& get_times() const;
		const std::vector<BkgFPType>& get_x() const;
		const std::vector<BkgFPType>& get_y() const;
		const std::vector<BkgFPType>& get_z() const;
		const std::vector<BkgFPType>& get_grid_x() const;
		const std::vector<BkgFPType>& get_grid_y() const;
		const std::vector<BkgFPType>& get_grid_z() const;
		const Vectors::Vector4D<BkgFPType>& get_ne() const;
		const Vectors::Vector4D<BkgFPType>& get_te() const;
		const Vectors::Vector4D<BkgFPType>& get_ti() const;
		const Vectors::Vector4D<BkgFPType>& get_vp() const;
		const Vectors::Vector4D<BkgFPType>& get_bX() const;
		const Vectors::Vector4D<BkgFPType>& get_bY() const;
		const Vectors::Vector4D<BkgFPType>& get_bZ() const;
		const Vectors::Vector4D<BkgFPType>& get_bR() const;
		const Vectors::Vector4D<BkgFPType>& get_gradbX() const;
		const Vectors::Vector4D<BkgFPType>& get_gradbY() const;
		const Vectors::Vector4D<BkgFPType>& get_gradbZ() const;
		const Vectors::Vector4D<BkgFPType>& get_eX() const;
		const Vectors::Vector4D<BkgFPType>& get_eY() const;
		const Vectors::Vector4D<BkgFPType>& get_eZ() const;
		const Vectors::Vector4D<BkgFPType>& get_uX() const;
		const Vectors::Vector4D<BkgFPType>& get_uY() const;
		const Vectors::Vector4D<BkgFPType>& get_uZ() const;
		const Vectors::Vector3D<BkgFPType>& get_X() const;
		const Vectors::Vector3D<BkgFPType>& get_Y() const;
		const Vectors::Vector3D<BkgFPType>& get_Z() const;
		const Vectors::Vector3D<BkgFPType>& get_grid_X() const;
		const Vectors::Vector3D<BkgFPType>& get_grid_Y() const;
		const Vectors::Vector3D<BkgFPType>& get_grid_Z() const;
		const Vectors::Vector3D<BkgFPType>& get_J() const;
		const Vectors::Vector3D<BkgFPType>& get_gij_00() const;
		const Vectors::Vector3D<BkgFPType>& get_gij_01() const;
		const Vectors::Vector3D<BkgFPType>& get_gij_02() const;
		const Vectors::Vector3D<BkgFPType>& get_gij_11() const;
		const Vectors::Vector3D<BkgFPType>& get_gij_12() const;
		const Vectors::Vector3D<BkgFPType>& get_gij_22() const;
		const Vectors::Vector3D<BkgFPType>& get_b_x() const;
		const Vectors::Vector3D<BkgFPType>& get_b_y() const;
		const Vectors::Vector3D<BkgFPType>& get_b_z() const;
		int get_dim1() const;
		int get_dim2() const;
		int get_dim3() const;
		int get_dim4() const;

		/**
		* @brief Safety check that any data that is read into Background is 
		* consistent with the already-set dimension.
		*
		* @param m_dim Current dimension value
		* @param in_dim New dimension value
		* @param data Name of data (for print out)
		* @param dim_num Dimension unmber (for print out)
		* @return Returns true if dimensions match, false (and a error message) if
		* not
		*/
		bool check_dim(const int m_dim, const int in_dim, 
			const std::string_view data, const int dim_num);

		/**
		* @brief Function to ensure dimensions are consistent across the different
		* vectors.
		*
		* Right now this just prints the appropriate error message and then lets
		* the program crash. I'd like to work in some kind of intentional 
		* terminate command, but not sure what the safest way is for that so I
		* let it crash.
		*
		* @param v Vector4D containing the dimensions we're checking against
		* @param data Name of data being loaded (for print out)
		*/
		template <typename T>
		void set_dims(Vectors::Vector4D<T>& v, const std::string_view data);

		// Helper functions to get the start/end times of the simulation.
		BkgFPType get_t_min() const;
		BkgFPType get_t_max() const;
		BkgFPType get_x_min() const;
		BkgFPType get_x_max() const;
		BkgFPType get_y_min() const;
		BkgFPType get_y_max() const;
		BkgFPType get_z_min() const;
		BkgFPType get_z_max() const;

		// Setters using move semantics. No idea if this is the proper
		// way to do it but it works.
		void move_into_times(std::vector<BkgFPType>& times);
		void move_into_x(std::vector<BkgFPType>& x);
		void move_into_y(std::vector<BkgFPType>& y);
		void move_into_z(std::vector<BkgFPType>& z);
		void move_into_grid_x(std::vector<BkgFPType>& grid_x);
		void move_into_grid_y(std::vector<BkgFPType>& grid_y);
		void move_into_grid_z(std::vector<BkgFPType>& grid_z);
		void move_into_ne(Vectors::Vector4D<BkgFPType>& ne);
		void move_into_te(Vectors::Vector4D<BkgFPType>& te);
		void move_into_ti(Vectors::Vector4D<BkgFPType>& ti);
		void move_into_vp(Vectors::Vector4D<BkgFPType>& vp);
		void move_into_bX(Vectors::Vector4D<BkgFPType>& bX);
		void move_into_bY(Vectors::Vector4D<BkgFPType>& bY);
		void move_into_bZ(Vectors::Vector4D<BkgFPType>& bZ);
		void move_into_bR(Vectors::Vector4D<BkgFPType>& bR);
		void move_into_gradbX(Vectors::Vector4D<BkgFPType>& gradbX);
		void move_into_gradbY(Vectors::Vector4D<BkgFPType>& gradbY);
		void move_into_gradbZ(Vectors::Vector4D<BkgFPType>& gradbZ);
		void move_into_eX(Vectors::Vector4D<BkgFPType>& eX); 
		void move_into_eY(Vectors::Vector4D<BkgFPType>& eY);
		void move_into_eZ(Vectors::Vector4D<BkgFPType>& eZ);
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
		void move_into_b_x(Vectors::Vector3D<BkgFPType>& b_x);
		void move_into_b_y(Vectors::Vector3D<BkgFPType>& b_y);
		void move_into_b_z(Vectors::Vector3D<BkgFPType>& b_z);
	};
}

#endif
