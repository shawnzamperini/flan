#ifndef VECTORS_H
#define VECTORS_H

#include <vector>
#include <iostream>

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

		// Constructor to create vector and fill with the data passed on
		Vector3D(std::vector<double> data_in, int dim1, int dim2, int dim3)
			: Vector3D(dim1, dim2, dim3)
		{
			// Copy over into the class
			// Note: Expensive! Could we use move semantics instead?
			m_data = data_in;
		}

		// Artificial constructor to allow us to compile read_gkyl.load_values.
		// That function selectively returns a brace-enclosed initializer list
		// for either a Vector3D or Vector4D based on the template parameter.
		// The compiler nonetheless needs to see a constructor for Vector3D
		// that accepts 5 inputs, since that is what Vector4D uses. This
		// constructor should never get used.
		Vector3D(std::vector<double> data_in, int dim1, int dim2, int dim3, 
				int dim4)
		{
			std::cerr << "Error! Vector3D was called with 4 dimensions. This"
				<< " is a programming error. Please fix\n";
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

		// Move the passed in vector to data. Hmmm... this is really
		// some sort of abstraction of operator=...
		void move_into_data(Vectors::Vector3D& vec)
		{
			// Copy over the dimensions
			m_dim1 = vec.get_dim1();
			m_dim2 = vec.get_dim2();
			m_dim3 = vec.get_dim3();

			// Move the data into this
			m_data = std::move(vec.m_data);
		}
	};

	class Vector4D
	{
	private:
		std::vector<double> m_data {};
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

		// Constructor to create vector and fill with the data passed on
		Vector4D(std::vector<double> data_in, int dim1, int dim2, int dim3, 
			int dim4)
			: Vector4D(dim1, dim2, dim3, dim4)
		{
			// Copy over into the class
			// Note: Expensive! Could we use move semantics instead?
			m_data = data_in;
		}

		// Artificial constructor to allow us to compile read_gkyl.load_values.
		// That function selectively returns a brace-enclosed initializer list
		// for either a Vector3D or Vector4D based on the template parameter.
		// The compiler nonetheless needs to see a constructor for Vector4D
		// that accepts these inputs, since that is what Vector3D uses. This
		// constructor should never get used.
		Vector4D(std::vector<double> data_in, int dim1, int dim2, int dim3)
		{
			std::cerr << "Error! Vector4D was called with 3 dimensions. This"
				<< " is a programming error. Please fix\n";
		}

		// We are making the copy copy constructor yell at us since we 
		// probably never want to make a copy of a Vector4D. Can comment
		// this out if we ever find this is not the case.
		Vector4D(const Vector4D& v)
			: m_data {v.m_data}
			, m_dim1 {v.m_dim1}
			, m_dim2 {v.m_dim2}
			, m_dim3 {v.m_dim3}
			, m_dim4 {v.m_dim4}
			, m_size {v.m_size}
		{
			std::cout << "Warning! Why are you making an expensive copy"
				<< " of a Vector4D?\n";
		}

		// We could alternatively prevent copies with delete.
		//Vector4D(const Vector4D& v) = delete;
	
		// Accessors
		int get_dim1() const {return m_dim1;}
		int get_dim2() const {return m_dim2;}
		int get_dim3() const {return m_dim3;}
		int get_dim4() const {return m_dim4;}

		// Currently not working because we aren't recalculating m_size when
		// using the move semantics functions below. Use get_data().size()
		// instead.
		//int get_size() const {return m_size;}
		std::vector<double> get_data() const {return m_data;}

		// Convert from 4D index to the 1D one used in the underlying m_data.
		int calc_index(int i, int j, int k, int l) const
		{
			return i * (m_dim2 * m_dim3 * m_dim4) + j 
				* (m_dim3 * m_dim4) + k * m_dim4 + l;
		}

		// Overload parentheses to act as indexing.
		double& operator()(int i, int j, int k, int l)
		{
			return m_data[calc_index(i, j, k, l)];
		}

		// Assignment operator moves data from one Vector4D into this one.
		Vector4D& operator=(Vector4D&& other) noexcept
		{
			// Self-assignment check
			if (this != &other)
			{	
				m_data = std::move(other.m_data);
			}
			return *this;
		}
			
		// Move the passed in vector to data. Hmmm... this is really
		// some sort of abstraction of operator=...
		void move_into_data(Vectors::Vector4D& vec)
		{
			// Copy over the dimensions
			m_dim1 = vec.get_dim1();
			m_dim2 = vec.get_dim2();
			m_dim3 = vec.get_dim3();
			m_dim4 = vec.get_dim4();

			// Move the data into this
			m_data = std::move(vec.m_data);
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

		// Represent that data as an actual 4D vector. Unsure if this is
		// needed but doesn't hurt to include.
		std::vector<std::vector<std::vector<std::vector<double>>>> as_4d() const
		{

			// Size check
			if (std::ssize(m_data) != m_dim1 * m_dim2 * m_dim3 * m_dim4)
			{
				std::cerr << "Error! The number of elements does not match" 
				<< " the specified dimensions.\n";	
			}

			// Create an actual 4D vector, fill in one element at a time. Must
			// initialize vector with the dimensions. Incredibly messy.
			std::vector<std::vector<std::vector<std::vector<double>>>> vec(
        		m_dim1, std::vector<std::vector<std::vector<double>>>(
            	m_dim2, std::vector<std::vector<double>>(
                m_dim3, std::vector<double>(m_dim4, 0.0))));
			std::cout << "as_4d dim1 = " << vec.size() << "\n";
			std::cout << "as_4d dim2 = " << vec[0].size() << "\n";
			std::cout << "as_4d dim3 = " << vec[0][0].size() << "\n";
			std::cout << "as_4d dim4 = " << vec[0][0][0].size() << "\n";
			for (int i {}; i < m_dim1; ++i)
			{
				for (int j {}; j < m_dim2; ++j)
				{
					for (int k {}; k < m_dim3; ++k)
					{
						for (int l {}; l < m_dim4; ++l)
						{
							vec[i][j][k][l] = m_data[calc_index(i,j,k,l)];
						}
					}
				}
			}

			// Return
			return vec;	
		}
	};
}

#endif