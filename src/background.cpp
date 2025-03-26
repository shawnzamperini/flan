/**
* @file background.cpp
* 
* Contains class definition for a loaded in Gkeyll background plasma. This will
* likely be expanded to included other codes one day. At the bottom of the file
* are template instantiations needed for it to compile. 
*/
#include <vector>
#include "vectors.h"
#include "background.h"

namespace Background
{

	// Accessors... I'll get around doxygen commenting these one day...
	const std::vector<double>& Background::get_times() const {return m_times;}
	const std::vector<double>& Background::get_x() const {return m_x;}
	const std::vector<double>& Background::get_y() const {return m_y;}
	const std::vector<double>& Background::get_z() const {return m_z;}
	const std::vector<double>& Background::get_grid_x() const {return m_grid_x;}
	const std::vector<double>& Background::get_grid_y() const {return m_grid_y;}
	const std::vector<double>& Background::get_grid_z() const {return m_grid_z;}
	const Vectors::Vector4D<double>& Background::get_ne() const {return m_ne;}
	const Vectors::Vector4D<double>& Background::get_te() const {return m_te;}
	const Vectors::Vector4D<double>& Background::get_ti() const {return m_ti;}
	const Vectors::Vector4D<double>& Background::get_vp() const {return m_vp;}
	const Vectors::Vector4D<double>& Background::get_bX() const {return m_bX;}
	const Vectors::Vector4D<double>& Background::get_bY() const {return m_bY;}
	const Vectors::Vector4D<double>& Background::get_bZ() const {return m_bZ;}
	const Vectors::Vector4D<double>& Background::get_eX() const {return m_eX;}
	const Vectors::Vector4D<double>& Background::get_eY() const {return m_eY;}
	const Vectors::Vector4D<double>& Background::get_eZ() const {return m_eZ;}
	const Vectors::Vector3D<double>& Background::get_X() const {return m_X;}
	const Vectors::Vector3D<double>& Background::get_Y() const {return m_Y;}
	const Vectors::Vector3D<double>& Background::get_Z() const {return m_Z;}
	const Vectors::Vector3D<double>& Background::get_grid_X() const {return m_grid_X;}
	const Vectors::Vector3D<double>& Background::get_grid_Y() const {return m_grid_Y;}
	const Vectors::Vector3D<double>& Background::get_grid_Z() const {return m_grid_Z;}
	const Vectors::Vector3D<double>& Background::get_J() const {return m_J;}
	int Background::get_dim1() const {return m_dim1;}
	int Background::get_dim2() const {return m_dim2;}
	int Background::get_dim3() const {return m_dim3;}
	int Background::get_dim4() const {return m_dim4;}

	bool Background::check_dim(const int m_dim, const int in_dim, 
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
					<< " in correctly or is from different Gkeyll runs."
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

	template <typename T>
	void Background::set_dims(Vectors::Vector4D<T>& v, 
		const std::string_view data)
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
	double Background::get_t_min() const {return m_times.front();}
	double Background::get_t_max() const {return m_times.back();}
	double Background::get_x_min() const {return m_x.front();}
	double Background::get_x_max() const {return m_x.back();}
	double Background::get_y_min() const {return m_y.front();}
	double Background::get_y_max() const {return m_y.back();}
	double Background::get_z_min() const {return m_z.front();}
	double Background::get_z_max() const {return m_z.back();}

	// Setters using move semantics. No idea if this is the proper
	// way to do it but it works.
	void Background::move_into_times(std::vector<double>& times)
	{
		m_times = std::move(times);	
	}
	void Background::move_into_x(std::vector<double>& x)
	{
		m_x = std::move(x);	
	}
	void Background::move_into_y(std::vector<double>& y)
	{
		m_y = std::move(y);	
	}
	void Background::move_into_z(std::vector<double>& z)
	{
		m_z = std::move(z);	
	}
	void Background::move_into_grid_x(std::vector<double>& grid_x)
	{
		m_grid_x = std::move(grid_x);	
	}
	void Background::move_into_grid_y(std::vector<double>& grid_y)
	{
		m_grid_y = std::move(grid_y);	
	}
	void Background::move_into_grid_z(std::vector<double>& grid_z)
	{
		m_grid_z = std::move(grid_z);	
	}
	void Background::move_into_ne(Vectors::Vector4D<double>& ne) 
	{
		set_dims(ne, "ne");
		m_ne.move_into_data(ne);
	}
	void Background::move_into_te(Vectors::Vector4D<double>& te) 
	{
		set_dims(te, "te");
		m_te.move_into_data(te);
	}
	void Background::move_into_ti(Vectors::Vector4D<double>& ti) 
	{
		set_dims(ti, "ti");
		m_ti.move_into_data(ti);
	}
	void Background::move_into_vp(Vectors::Vector4D<double>& vp) 
	{
		set_dims(vp, "vp");
		m_vp.move_into_data(vp);
	}
	void Background::move_into_bX(Vectors::Vector4D<double>& bX) 
	{
		set_dims(bX, "bX");
		m_bX.move_into_data(bX);
	}
	void Background::move_into_bY(Vectors::Vector4D<double>& bY) 
	{
		set_dims(bY, "bY");
		m_bY.move_into_data(bY);
	}
	void Background::move_into_bZ(Vectors::Vector4D<double>& bZ) 
	{
		set_dims(bZ, "bZ");
		m_bZ.move_into_data(bZ);
	}
	void Background::move_into_eX(Vectors::Vector4D<double>& eX) 
	{
		set_dims(eX, "eX");
		m_eX.move_into_data(eX);
	}
	void Background::move_into_eY(Vectors::Vector4D<double>& eY) 
	{
		set_dims(eY, "eY");
		m_eY.move_into_data(eY);
	}
	void Background::move_into_eZ(Vectors::Vector4D<double>& eZ) 
	{
		set_dims(eZ, "eZ");
		m_eZ.move_into_data(eZ);
	}
	void Background::move_into_X(Vectors::Vector3D<double>& X) 
	{
		// Need a Vector3D version of set_dims
		//set_dims(X, "X");
		m_X.move_into_data(X);
	}
	void Background::move_into_Y(Vectors::Vector3D<double>& Y) 
	{
		// Need a Vector3D version of set_dims
		//set_dims(Y, "Y");
		m_Y.move_into_data(Y);
	}
	void Background::move_into_Z(Vectors::Vector3D<double>& Z) 
	{
		// Need a Vector3D version of set_dims
		//set_dims(Z, "Z");
		m_Z.move_into_data(Z);
	}
	void Background::move_into_grid_X(Vectors::Vector3D<double>& grid_X) 
	{
		// Need a Vector3D version of set_dims
		//set_dims(grid_X, "grid_X");
		m_grid_X.move_into_data(grid_X);
	}
	void Background::move_into_grid_Y(Vectors::Vector3D<double>& grid_Y) 
	{
		// Need a Vector3D version of set_dims
		//set_dims(grid_Y, "grid_Y");
		m_grid_Y.move_into_data(grid_Y);
	}
	void Background::move_into_grid_Z(Vectors::Vector3D<double>& grid_Z) 
	{
		// Need a Vector3D version of set_dims
		//set_dims(grid_Z, "grid_Z");
		m_grid_Z.move_into_data(grid_Z);
	}
	void Background::move_into_J(Vectors::Vector3D<double>& J) 
	{
		// Need a Vector3D version of set_dims
		//set_dims(J, "J");
		m_J.move_into_data(J);
	}
}

// This is a pecularity of splitting a class declarations and definitions 
// a header and source file. We need to instatiate the templates with the 
// needed definitions so the linker can see them. Hurts flexibility, but not
// an issue since it's pretty straightforward.
template void Background::Background::set_dims<int>(
	Vectors::Vector4D<int>& v, const std::string_view data);
template void Background::Background::set_dims<double>(
	Vectors::Vector4D<double>& v, const std::string_view data);
