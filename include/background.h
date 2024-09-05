#ifndef BACKGROUND_H
#define BACKGROUND_H

#include <vector>
#include "vectors.h"

namespace Background
{
	// Encapsulating class containing the Gkeyll background.
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
		int m_dim1 {};
		int m_dim2 {};
		int m_dim3 {};
		int m_dim4 {};

	public:

		// Accessors
		std::vector<double> get_times() const {return m_times;}
		std::vector<double> get_x() const {return m_x;}
		std::vector<double> get_y() const {return m_y;}
		std::vector<double> get_z() const {return m_z;}
		std::vector<double> get_grid_x() const {return m_grid_x;}
		std::vector<double> get_grid_y() const {return m_grid_y;}
		std::vector<double> get_grid_z() const {return m_grid_z;}
		const Vectors::Vector4D<double>& get_ne() const {return m_ne;}
		const Vectors::Vector4D<double>& get_te() const {return m_te;}
		const Vectors::Vector4D<double>& get_ti() const {return m_ti;}
		const Vectors::Vector4D<double>& get_vp() const {return m_vp;}
		const Vectors::Vector4D<double>& get_b() const {return m_b;}
		int get_dim1() const {return m_dim1;}
		int get_dim2() const {return m_dim2;}
		int get_dim3() const {return m_dim3;}
		int get_dim4() const {return m_dim4;}

		bool check_dim(const int m_dim, const int in_dim, 
			const std::string_view data, const int dim_num)
		{
			// Assume a zero dimension value just means it's not set.
			if (m_dim != 0)
			{
				if (m_dim != in_dim)
				{
					std::cerr << "Error! In Background, " << data << " has" 
						<< " different dimensions from what already has been" 
						<< " set. This indicates the data is not being read" 
					  	<< "	in correctly or is from different Gkeyll runs."
					  	<< '\n';
					std::cerr << "  Previous dim: m_dim" << dim_num << " = " 
						<< m_dim << '\n';
					std::cerr << "  New dim:     in_dim" << dim_num << " = " 
						<< in_dim << '\n';

					return false;
				}
			}
			return true;
		}

		// Function to ensure dimensions are consistent across the different
		// vectors.
		template <typename T>
		void set_dims(Vectors::Vector4D<T>& v, const std::string_view data)
		{
			// Check dimensions first, printing errors if something ain't
			// right.
			bool dim_good1 {check_dim(m_dim1, v.get_dim1(), data, 1)};
			bool dim_good2 {check_dim(m_dim2, v.get_dim2(), data, 2)};
			bool dim_good3 {check_dim(m_dim3, v.get_dim3(), data, 3)};
			bool dim_good4 {check_dim(m_dim4, v.get_dim4(), data, 4)};

			// Then overwrite the dimensions. This will either:
			//  - Set the dimensions from 0 to their correct value (good)
			//  - Overwrite the dimension with the same value (good)
			//  - Overwrite the dimension with a different value (bad, error
			//    (message from check_dim should be yelling at you before
			//    this).
			if (dim_good1 && dim_good2 && dim_good3 && dim_good4)
			{
				m_dim1 = v.get_dim1();
				m_dim2 = v.get_dim2();
				m_dim3 = v.get_dim3();
				m_dim4 = v.get_dim4();
			}
			else
			{
				std::cerr << "Fix dimension error before continuing\n";
				// Best way to terminate program?
			}

		}

		// Helper functions to get the start/end times of the simulation.
		double get_t_min() const {return m_times.front();}
		double get_t_max() const {return m_times.back();}
		double get_x_min() const {return m_x.front();}
		double get_x_max() const {return m_x.back();}
		double get_y_min() const {return m_y.front();}
		double get_y_max() const {return m_y.back();}
		double get_z_min() const {return m_z.front();}
		double get_z_max() const {return m_z.back();}

		// Setters using move semantics. No idea if this is the proper
		// way to do it but it works.
		void move_into_times(std::vector<double>& times)
		{
			m_times = std::move(times);	
		}
		void move_into_x(std::vector<double>& x)
		{
			m_x = std::move(x);	
		}
		void move_into_y(std::vector<double>& y)
		{
			m_y = std::move(y);	
		}
		void move_into_z(std::vector<double>& z)
		{
			m_z = std::move(z);	
		}
		void move_into_grid_x(std::vector<double>& grid_x)
		{
			m_grid_x = std::move(grid_x);	
		}
		void move_into_grid_y(std::vector<double>& grid_y)
		{
			m_grid_y = std::move(grid_y);	
		}
		void move_into_grid_z(std::vector<double>& grid_z)
		{
			m_grid_z = std::move(grid_z);	
		}
		void move_into_ne(Vectors::Vector4D<double>& ne) 
		{
			set_dims(ne, "ne");
			m_ne.move_into_data(ne);
		}
		void move_into_te(Vectors::Vector4D<double>& te) 
		{
			set_dims(te, "te");
			m_te.move_into_data(te);
		}
		void move_into_ti(Vectors::Vector4D<double>& ti) 
		{
			set_dims(ti, "ti");
			m_ti.move_into_data(ti);
		}
		void move_into_vp(Vectors::Vector4D<double>& vp) 
		{
			set_dims(vp, "vp");
			m_vp.move_into_data(vp);
		}
		void move_into_b(Vectors::Vector4D<double>& b) 
		{
			// This is skipped for now since elctrostatic so the time dimension
			// is always just one long.
			//set_dims(b, "b");
			m_b.move_into_data(b);
		}
	};
}

#endif
