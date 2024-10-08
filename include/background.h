/**
* @file background.h
* @brief Header file for background.cpp
*/
#ifndef BACKGROUND_H
#define BACKGROUND_H

#include <vector>
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
		std::vector<double> m_times {};
		std::vector<double> m_x {};
		std::vector<double> m_y {};
		std::vector<double> m_z {};
		std::vector<double> m_grid_x {};
		std::vector<double> m_grid_y {};
		std::vector<double> m_grid_z {};
		Vectors::Vector4D<double> m_ne {};
		Vectors::Vector4D<double> m_te {};
		Vectors::Vector4D<double> m_ti {};
		Vectors::Vector4D<double> m_vp {};
		Vectors::Vector4D<double> m_b {};
		Vectors::Vector4D<double> m_ex {};
		Vectors::Vector4D<double> m_ey {};
		Vectors::Vector4D<double> m_ez {};
		int m_dim1 {};
		int m_dim2 {};
		int m_dim3 {};
		int m_dim4 {};

	public:

		// Accessors
		const std::vector<double>& get_times() const;
		const std::vector<double>& get_x() const;
		const std::vector<double>& get_y() const;
		const std::vector<double>& get_z() const;
		const std::vector<double>& get_grid_x() const;
		const std::vector<double>& get_grid_y() const;
		const std::vector<double>& get_grid_z() const;
		const Vectors::Vector4D<double>& get_ne() const;
		const Vectors::Vector4D<double>& get_te() const;
		const Vectors::Vector4D<double>& get_ti() const;
		const Vectors::Vector4D<double>& get_vp() const;
		const Vectors::Vector4D<double>& get_b() const;
		const Vectors::Vector4D<double>& get_ex() const;
		const Vectors::Vector4D<double>& get_ey() const;
		const Vectors::Vector4D<double>& get_ez() const;
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
		double get_t_min() const;
		double get_t_max() const;
		double get_x_min() const;
		double get_x_max() const;
		double get_y_min() const;
		double get_y_max() const;
		double get_z_min() const;
		double get_z_max() const;

		// Setters using move semantics. No idea if this is the proper
		// way to do it but it works.
		void move_into_times(std::vector<double>& times);
		void move_into_x(std::vector<double>& x);
		void move_into_y(std::vector<double>& y);
		void move_into_z(std::vector<double>& z);
		void move_into_grid_x(std::vector<double>& grid_x);
		void move_into_grid_y(std::vector<double>& grid_y);
		void move_into_grid_z(std::vector<double>& grid_z);
		void move_into_ne(Vectors::Vector4D<double>& ne);
		void move_into_te(Vectors::Vector4D<double>& te);
		void move_into_ti(Vectors::Vector4D<double>& ti);
		void move_into_vp(Vectors::Vector4D<double>& vp);
		void move_into_b(Vectors::Vector4D<double>& b);
		void move_into_ex(Vectors::Vector4D<double>& ex); 
		void move_into_ey(Vectors::Vector4D<double>& ey);
		void move_into_ez(Vectors::Vector4D<double>& ez);
	};
}

#endif
