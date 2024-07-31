#ifndef VECTORS_H
#define VECTORS_H

#include <vector>

// This file contains some simple implementations of multidimensional vectors.
// Instead of storing them as actual mutlidimensional vectors, they are stored
// as flatten vectors and then () is overloaded to act as indexing. 
namespace Vectors
{
	class Vector3D
	{
	private:
		std::vector<double> m_data;
		int m_dim1 {};
		int m_dim2 {};
		int m_dim3 {};
		int m_size {};

	public:
		// Default constructor
		Vector3D()
		{}

		// Constructor to allow creating an empty vector
		Vector3D(int dim1, int dim2, int dim3)
			: m_dim1 {dim1}
			, m_dim2 {dim2}
			, m_dim3 {dim3}
		{
			m_size = dim1 * dim2 * dim3;
			m_data.resize(m_size);
		}

		// Create vector and fill with the data passed on
		Vector3D(std::vector<double> data_in, int dim1, int dim2, int dim3)
			: Vector3D(dim1, dim2, dim3)
		{
			// Copy over into the class
			// Note: Expensive! Could we use move semantics instead?
			m_data = data_in;
		}
	
		// Accessors
		int get_dim1() {return m_dim1;}
		int get_dim2() {return m_dim2;}
		int get_dim3() {return m_dim3;}
		int get_size() {return m_size;}
		std::vector<double> get_data() {return m_data;}

		// Overload parentheses to act as indexing.
		double& operator()(int i, int j, int k)
		{
			int index = i * (m_dim2 * m_dim3) + j * m_dim3 + k;
			return m_data[index];
		}
	};

	class Vector4D
	{
	private:
		std::vector<double> m_data;
		int m_dim1 {};
		int m_dim2 {};
		int m_dim3 {};
		int m_dim4 {};
		int m_size {};

	public:
		// Default constructor
		Vector4D()
		{}

		// Constructor to allow creating an empty vector
		Vector4D(int dim1, int dim2, int dim3, int dim4)
			: m_dim1 {dim1}
			, m_dim2 {dim2}
			, m_dim3 {dim3}
			, m_dim4 {dim4}
		{
			m_size = dim1 * dim2 * dim3 * dim4;
			m_data.resize(m_size);
		}

		// Create vector and fill with the data passed on
		Vector4D(std::vector<double> data_in, int dim1, int dim2, int dim3, int dim4)
			: Vector4D(dim1, dim2, dim3, dim4)
		{
			// Copy over into the class
			// Note: Expensive! Could we use move semantics instead?
			m_data = data_in;
		}
	
		// Accessors
		int get_dim1() {return m_dim1;}
		int get_dim2() {return m_dim2;}
		int get_dim3() {return m_dim3;}
		int get_dim4() {return m_dim4;}
		int get_size() {return m_size;}
		std::vector<double> get_data() {return m_data;}	

		// Overload parentheses to act as indexing.
		double& operator()(int i, int j, int k, int l)
		{
			int index = i * (m_dim2 * m_dim3 * m_dim4) + j * (m_dim3 * m_dim4)
				+ k * m_dim4 + l;
			return m_data[index];
		}

		// Slice a single dimension (which returns a vector with one 
		// lower dimensionality).
		Vector3D slice_dim4(int dim_index)
		{
			Vector3D vec {m_dim1, m_dim2, m_dim3};
			for (int i {}; i < m_dim1; ++i)
			{
				for (int j {}; j < m_dim2; ++j)
				{
					for (int k {}; k < m_dim3; ++k)
					{
						int index = i * (m_dim2 * m_dim3 * m_dim4) + j 
							* (m_dim3 * m_dim4) + k * m_dim4 + dim_index;
		   				vec(i, j, k) = m_data[index];
					}
				}
			}
			return vec;
		}
	};
}

#endif
