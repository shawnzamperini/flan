/**
* @file background.cpp
* 
* Contains class definition for a loaded in Gkeyll background plasma. This will
* likely be expanded to included other codes one day. At the bottom of the file
* are template instantiations needed for it to compile. 
*/
#include <vector>

#include "background.h"
#include "flan_types.h"
#include "vectors.h"

namespace Background
{

	// Accessors... I'll get around doxygen commenting these one day...
	const std::vector<BkgFPType>& Background::get_times() const {return m_times;}
	const std::vector<BkgFPType>& Background::get_x() const {return m_x;}
	const std::vector<BkgFPType>& Background::get_y() const {return m_y;}
	const std::vector<BkgFPType>& Background::get_z() const {return m_z;}
	const std::vector<BkgFPType>& Background::get_grid_x() const {return m_grid_x;}
	const std::vector<BkgFPType>& Background::get_grid_y() const {return m_grid_y;}
	const std::vector<BkgFPType>& Background::get_grid_z() const {return m_grid_z;}
	const Vectors::Vector4D<BkgFPType>& Background::get_ne() const {return m_ne;}
	const Vectors::Vector4D<BkgFPType>& Background::get_te() const {return m_te;}
	const Vectors::Vector4D<BkgFPType>& Background::get_ti() const {return m_ti;}
	const Vectors::Vector4D<BkgFPType>& Background::get_vp() const {return m_vp;}
	const Vectors::Vector4D<BkgFPType>& Background::get_bX() const {return m_bX;}
	const Vectors::Vector4D<BkgFPType>& Background::get_bY() const {return m_bY;}
	const Vectors::Vector4D<BkgFPType>& Background::get_bZ() const {return m_bZ;}
	const Vectors::Vector4D<BkgFPType>& Background::get_bR() const {return m_bR;}
	const Vectors::Vector4D<BkgFPType>& Background::get_gradbX() const {return m_gradbX;}
	const Vectors::Vector4D<BkgFPType>& Background::get_gradbY() const {return m_gradbY;}
	const Vectors::Vector4D<BkgFPType>& Background::get_gradbZ() const {return m_gradbZ;}
	const Vectors::Vector4D<BkgFPType>& Background::get_eX() const {return m_eX;}
	const Vectors::Vector4D<BkgFPType>& Background::get_eY() const {return m_eY;}
	const Vectors::Vector4D<BkgFPType>& Background::get_eZ() const {return m_eZ;}
	const Vectors::Vector4D<BkgFPType>& Background::get_gradeX() const 
		{return m_gradeX;}
	const Vectors::Vector4D<BkgFPType>& Background::get_gradeY() const 
		{return m_gradeY;}
	const Vectors::Vector4D<BkgFPType>& Background::get_gradeZ() const 
		{return m_gradeZ;}
	const Vectors::Vector4D<BkgFPType>& Background::get_uX() const {return m_uX;}
	const Vectors::Vector4D<BkgFPType>& Background::get_uY() const {return m_uY;}
	const Vectors::Vector4D<BkgFPType>& Background::get_uZ() const {return m_uZ;}
	const Vectors::Vector3D<BkgFPType>& Background::get_X() const {return m_X;}
	const Vectors::Vector3D<BkgFPType>& Background::get_Y() const {return m_Y;}
	const Vectors::Vector3D<BkgFPType>& Background::get_Z() const {return m_Z;}
	const Vectors::Vector3D<BkgFPType>& Background::get_grid_X() const {return m_grid_X;}
	const Vectors::Vector3D<BkgFPType>& Background::get_grid_Y() const {return m_grid_Y;}
	const Vectors::Vector3D<BkgFPType>& Background::get_grid_Z() const {return m_grid_Z;}
	const Vectors::Vector3D<BkgFPType>& Background::get_J() const {return m_J;}
	const Vectors::Vector3D<BkgFPType>& Background::get_gij_00() const {
        return m_gij_00;
    }
	const Vectors::Vector3D<BkgFPType>& Background::get_gij_01() const {
        return m_gij_01;
    }
	const Vectors::Vector3D<BkgFPType>& Background::get_gij_02() const {
        return m_gij_02;
    }
	const Vectors::Vector3D<BkgFPType>& Background::get_gij_11() const {
        return m_gij_11;
    }
	const Vectors::Vector3D<BkgFPType>& Background::get_gij_12() const {
        return m_gij_12;
    }
	const Vectors::Vector3D<BkgFPType>& Background::get_gij_22() const {
        return m_gij_22;
    }
	const Vectors::Vector3D<BkgFPType>& Background::get_b_x() const {
        return m_b_x;
    }
	const Vectors::Vector3D<BkgFPType>& Background::get_b_y() const {
        return m_b_y;
    }
	const Vectors::Vector3D<BkgFPType>& Background::get_b_z() const {
        return m_b_z;
    }

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
	BkgFPType Background::get_t_min() const {return m_times.front();}
	BkgFPType Background::get_t_max() const {return m_times.back();}
	BkgFPType Background::get_x_min() const {return m_x.front();}
	BkgFPType Background::get_x_max() const {return m_x.back();}
	BkgFPType Background::get_y_min() const {return m_y.front();}
	BkgFPType Background::get_y_max() const {return m_y.back();}
	BkgFPType Background::get_z_min() const {return m_z.front();}
	BkgFPType Background::get_z_max() const {return m_z.back();}

	// Setters using move semantics. No idea if this is the proper
	// way to do it but it works.
	void Background::move_into_times(std::vector<BkgFPType>& times)
	{
		m_times = std::move(times);	
	}
	void Background::move_into_x(std::vector<BkgFPType>& x)
	{
		m_x = std::move(x);	
	}
	void Background::move_into_y(std::vector<BkgFPType>& y)
	{
		m_y = std::move(y);	
	}
	void Background::move_into_z(std::vector<BkgFPType>& z)
	{
		m_z = std::move(z);	
	}
	void Background::move_into_grid_x(std::vector<BkgFPType>& grid_x)
	{
		m_grid_x = std::move(grid_x);	
	}
	void Background::move_into_grid_y(std::vector<BkgFPType>& grid_y)
	{
		m_grid_y = std::move(grid_y);	
	}
	void Background::move_into_grid_z(std::vector<BkgFPType>& grid_z)
	{
		m_grid_z = std::move(grid_z);	
	}
	void Background::move_into_ne(Vectors::Vector4D<BkgFPType>& ne) 
	{
		set_dims(ne, "ne");
		m_ne.move_into_data(ne);
	}
	void Background::move_into_te(Vectors::Vector4D<BkgFPType>& te) 
	{
		set_dims(te, "te");
		m_te.move_into_data(te);
	}
	void Background::move_into_ti(Vectors::Vector4D<BkgFPType>& ti) 
	{
		set_dims(ti, "ti");
		m_ti.move_into_data(ti);
	}
	void Background::move_into_vp(Vectors::Vector4D<BkgFPType>& vp) 
	{
		set_dims(vp, "vp");
		m_vp.move_into_data(vp);
	}
	void Background::move_into_bX(Vectors::Vector4D<BkgFPType>& bX) 
	{
		set_dims(bX, "bX");
		m_bX.move_into_data(bX);
	}
	void Background::move_into_bY(Vectors::Vector4D<BkgFPType>& bY) 
	{
		set_dims(bY, "bY");
		m_bY.move_into_data(bY);
	}
	void Background::move_into_bZ(Vectors::Vector4D<BkgFPType>& bZ) 
	{
		set_dims(bZ, "bZ");
		m_bZ.move_into_data(bZ);
	}
	void Background::move_into_bR(Vectors::Vector4D<BkgFPType>& bR) 
	{
		set_dims(bR, "bR");
		m_bR.move_into_data(bR);
	}
	void Background::move_into_gradbX(Vectors::Vector4D<BkgFPType>& gradbX) 
	{
		set_dims(gradbX, "gradbX");
		m_gradbX.move_into_data(gradbX);
	}
	void Background::move_into_gradbY(Vectors::Vector4D<BkgFPType>& gradbY) 
	{
		set_dims(gradbY, "gradbY");
		m_gradbY.move_into_data(gradbY);
	}
	void Background::move_into_gradbZ(Vectors::Vector4D<BkgFPType>& gradbZ) 
	{
		set_dims(gradbZ, "gradbZ");
		m_gradbZ.move_into_data(gradbZ);
	}
	void Background::move_into_eX(Vectors::Vector4D<BkgFPType>& eX) 
	{
		set_dims(eX, "eX");
		m_eX.move_into_data(eX);
	}
	void Background::move_into_eY(Vectors::Vector4D<BkgFPType>& eY) 
	{
		set_dims(eY, "eY");
		m_eY.move_into_data(eY);
	}
	void Background::move_into_eZ(Vectors::Vector4D<BkgFPType>& eZ) 
	{
		set_dims(eZ, "eZ");
		m_eZ.move_into_data(eZ);
	}
	void Background::move_into_gradeX(Vectors::Vector4D<BkgFPType>& gradeX) 
	{
		set_dims(gradeX, "gradeX");
		m_gradeX.move_into_data(gradeX);
	}
	void Background::move_into_gradeY(Vectors::Vector4D<BkgFPType>& gradeY) 
	{
		set_dims(gradeY, "gradeY");
		m_gradeY.move_into_data(gradeY);
	}
	void Background::move_into_gradeZ(Vectors::Vector4D<BkgFPType>& gradeZ) 
	{
		set_dims(gradeZ, "gradeZ");
		m_gradeZ.move_into_data(gradeZ);
	}
	void Background::move_into_uX(Vectors::Vector4D<BkgFPType>& uX) 
	{
		set_dims(uX, "uX");
		m_uX.move_into_data(uX);
	}
	void Background::move_into_uY(Vectors::Vector4D<BkgFPType>& uY) 
	{
		set_dims(uY, "uY");
		m_uY.move_into_data(uY);
	}
	void Background::move_into_uZ(Vectors::Vector4D<BkgFPType>& uZ) 
	{
		set_dims(uZ, "uZ");
		m_uZ.move_into_data(uZ);
	}
	void Background::move_into_X(Vectors::Vector3D<BkgFPType>& X) 
	{
		// Need a Vector3D version of set_dims
		//set_dims(X, "X");
		m_X.move_into_data(X);
	}
	void Background::move_into_Y(Vectors::Vector3D<BkgFPType>& Y) 
	{
		// Need a Vector3D version of set_dims
		//set_dims(Y, "Y");
		m_Y.move_into_data(Y);
	}
	void Background::move_into_Z(Vectors::Vector3D<BkgFPType>& Z) 
	{
		// Need a Vector3D version of set_dims
		//set_dims(Z, "Z");
		m_Z.move_into_data(Z);
	}
	void Background::move_into_grid_X(Vectors::Vector3D<BkgFPType>& grid_X) 
	{
		// Need a Vector3D version of set_dims
		//set_dims(grid_X, "grid_X");
		m_grid_X.move_into_data(grid_X);
	}
	void Background::move_into_grid_Y(Vectors::Vector3D<BkgFPType>& grid_Y) 
	{
		// Need a Vector3D version of set_dims
		//set_dims(grid_Y, "grid_Y");
		m_grid_Y.move_into_data(grid_Y);
	}
	void Background::move_into_grid_Z(Vectors::Vector3D<BkgFPType>& grid_Z) 
	{
		// Need a Vector3D version of set_dims
		//set_dims(grid_Z, "grid_Z");
		m_grid_Z.move_into_data(grid_Z);
	}
	void Background::move_into_J(Vectors::Vector3D<BkgFPType>& J) 
	{
		// Need a Vector3D version of set_dims
		//set_dims(J, "J");
		m_J.move_into_data(J);
	}
	void Background::move_into_gij_00(Vectors::Vector3D<BkgFPType>& gij_00) 
	{
		// Need a Vector3D version of set_dims
		//set_dims(gij_00, "gij_00");
		m_gij_00.move_into_data(gij_00);
	}
	void Background::move_into_gij_01(Vectors::Vector3D<BkgFPType>& gij_01) 
	{
		// Need a Vector3D version of set_dims
		//set_dims(gij_01, "gij_01");
		m_gij_01.move_into_data(gij_01);
	}
	void Background::move_into_gij_02(Vectors::Vector3D<BkgFPType>& gij_02) 
	{
		// Need a Vector3D version of set_dims
		//set_dims(gij_02, "gij_02");
		m_gij_02.move_into_data(gij_02);
	}
	void Background::move_into_gij_11(Vectors::Vector3D<BkgFPType>& gij_11) 
	{
		// Need a Vector3D version of set_dims
		//set_dims(gij_11, "gij_11");
		m_gij_11.move_into_data(gij_11);
	}
	void Background::move_into_gij_12(Vectors::Vector3D<BkgFPType>& gij_12) 
	{
		// Need a Vector3D version of set_dims
		//set_dims(gij_12, "gij_12");
		m_gij_12.move_into_data(gij_12);
	}
	void Background::move_into_gij_22(Vectors::Vector3D<BkgFPType>& gij_22) 
	{
		// Need a Vector3D version of set_dims
		//set_dims(gij_22, "gij_22");
		m_gij_22.move_into_data(gij_22);
	}
	void Background::move_into_b_x(Vectors::Vector3D<BkgFPType>& b_x) 
	{
		m_b_x.move_into_data(b_x);
	}
	void Background::move_into_b_y(Vectors::Vector3D<BkgFPType>& b_y) 
	{
		m_b_y.move_into_data(b_y);
	}
	void Background::move_into_b_z(Vectors::Vector3D<BkgFPType>& b_z) 
	{
		m_b_z.move_into_data(b_z);
	}
}

// This is a pecularity of splitting a class declarations and definitions 
// a header and source file. We need to instatiate the templates with the 
// needed definitions so the linker can see them. Hurts flexibility, but not
// an issue since it's pretty straightforward.
template void Background::Background::set_dims<int>(
	Vectors::Vector4D<int>& v, const std::string_view data);
template void Background::Background::set_dims<BkgFPType>(
	Vectors::Vector4D<BkgFPType>& v, const std::string_view data);
