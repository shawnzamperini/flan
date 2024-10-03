/**
* @file vectors.cpp
*
* @brief Implementations of multidimensional vectors
*
* This file contains some simple implementations of multidimensional vectors.
* Instead of storing them as actual multidimensional vectors, they are stored
* as flatten vectors and then () is overloaded to act as indexing. 
*
* The bottom of this file includes instatiations of Vector3D and Vector4D for
* the different types needed (e.g., int, double, ...). This is a pecularity
* needed when splitting the declarations and definitions between a header and
* source file, needed so the linker can find the classes for each template
* type.
*/

#include <vector>

#include "vectors.h"

namespace Vectors
{

	// ************************
	// * Vector3D definitions *
	// ************************

	template <typename T>
	Vector3D<T>::Vector3D()
	{}

	template <typename T>
	Vector3D<T>::Vector3D(int dim1, int dim2, int dim3)
		: m_dim1 {dim1}
		, m_dim2 {dim2}
		, m_dim3 {dim3}
	{
		m_size = dim1 * dim2 * dim3;
		m_data.resize(m_size);
	}

	template <typename T>
	Vector3D<T>::Vector3D(std::vector<T> data_in, int dim1, int dim2, 
		int dim3)
		: Vector3D<T>::Vector3D(dim1, dim2, dim3)
	{
		// Copy over into the class
		// Note: Expensive! Could we use move semantics instead?
		m_data = data_in;
	}

	template <typename T>
	Vector3D<T>::Vector3D(std::vector<T> data_in, int dim1, int dim2, int dim3, 
			int dim4)
	{
		std::cerr << "Error! Vector3D was called with 4 dimensions. This"
			<< " is a programming error. Please fix\n";
	}

	template <typename T>
	int Vector3D<T>::get_dim1() {return m_dim1;}

	template <typename T>
	int Vector3D<T>::get_dim2() {return m_dim2;}

	template <typename T>
	int Vector3D<T>::get_dim3() {return m_dim3;}

	template <typename T>
	int Vector3D<T>::get_size() {return m_size;}

	template <typename T>
	std::vector<T> Vector3D<T>::get_data() {return m_data;}
	
	template <typename T>
	int Vector3D<T>::calc_index(const int i, const int j, const int k) const
	{
		return i * (m_dim2 * m_dim3) + j * m_dim3 + k;
	}
		
	template <typename T>
	T& Vector3D<T>::operator()(int i, int j, int k)
	{
		return m_data[calc_index(i, j, k)];
	}

	template <typename T>
	const T& Vector3D<T>::operator()(const int i, const int j, 
		const int k) const
	{
		return m_data[calc_index(i, j, k)];
	}

	template <typename T>
	void Vector3D<T>::move_into_data(Vectors::Vector3D<T>& vec)
	{
		// Copy over the dimensions
		m_dim1 = vec.get_dim1();
		m_dim2 = vec.get_dim2();
		m_dim3 = vec.get_dim3();

		// Move the data into this
		m_data = std::move(vec.m_data);
	}
	
	// ************************
	// * Vector4D definitions *
	// ************************

	template <typename T>
	Vector4D<T>::Vector4D()
	{}

	template <typename T>
	Vector4D<T>::Vector4D(int dim1, int dim2, int dim3, int dim4)
		: m_dim1 {dim1}
		, m_dim2 {dim2}
		, m_dim3 {dim3}
		, m_dim4 {dim4}
	{
		m_size = dim1 * dim2 * dim3 * dim4;
		m_data.resize(m_size);
	}

	template <typename T>
	Vector4D<T>::Vector4D(std::vector<T> data_in, int dim1, int dim2, int dim3, 
		int dim4)
		: Vector4D<T>::Vector4D(dim1, dim2, dim3, dim4)
	{
		// Copy over into the class
		// Note: Expensive! Could we use move semantics instead?
		m_data = data_in;
	}

	template <typename T>
	Vector4D<T>::Vector4D(std::vector<T> data_in, int dim1, int dim2, int dim3)
	{
		std::cerr << "Error! Vector4D was called with 3 dimensions. This"
			<< " is a programming error. Please fix\n";
	}

	template <typename T>
	Vector4D<T>::Vector4D(const Vector4D<T>& v)
		: m_data {v.m_data}
		, m_dim1 {v.m_dim1}
		, m_dim2 {v.m_dim2}
		, m_dim3 {v.m_dim3}
		, m_dim4 {v.m_dim4}
		, m_size {v.m_size}
	{
		// This is a useful warning to check. Generally don't want to do
		// this, but it is needed when creating private copies of vectors
		// for each OpenMP thread.
		//std::cout << "Warning! Why are you making an expensive copy"
		//	<< " of a Vector4D?\n";
	}
	
	template <typename T>
	int Vector4D<T>::get_dim1() const {return m_dim1;}

	template <typename T>
	int Vector4D<T>::get_dim2() const {return m_dim2;}

	template <typename T>
	int Vector4D<T>::get_dim3() const {return m_dim3;}

	template <typename T>
	int Vector4D<T>::get_dim4() const {return m_dim4;}

	template <typename T>
	const std::vector<T>& Vector4D<T>::get_data() const {return m_data;}

	template <typename T>
	int Vector4D<T>::calc_index(const int i, const int j, const int k, 
		const int l) const
	{
		return i * (m_dim2 * m_dim3 * m_dim4) + j 
			* (m_dim3 * m_dim4) + k * m_dim4 + l;
	}

	template <typename T>
	T& Vector4D<T>::operator()(const int i, const int j, const int k, 
		const int l)
	{
		return m_data[calc_index(i, j, k, l)];
	}

	template <typename T>
	const T& Vector4D<T>::operator()(const int i, const int j, 
		const int k, const int l) const
	{
		return m_data[calc_index(i, j, k, l)];
	}

	template <typename T>
	Vector4D<T>& Vector4D<T>::operator=(Vector4D<T>&& other) 
		noexcept
	{
		// Self-assignment check
		if (this != &other)
		{	
			m_data = std::move(other.m_data);
		}
		return *this;
	}

	template <typename T>
	bool Vector4D<T>::check_same_shape(const Vector4D<T>& other) const
	{
		if (m_dim1 != other.m_dim1 || m_dim2 != other.m_dim2 || 
			m_dim3 != other.m_dim3 || m_dim4 != other.m_dim4)
		{
			std::cerr << "Error! Vector4D's are not of the same shape.\n";
			std::cerr << "  this   other\n";
			std::cerr << std::setw(6) << m_dim1 << std::setw(8) 
				<< other.m_dim1 << '\n'; 
			std::cerr << std::setw(6) << m_dim2 << std::setw(8) 
				<< other.m_dim2 << '\n'; 
			std::cerr << std::setw(6) << m_dim3 << std::setw(8) 
				<< other.m_dim3 << '\n'; 
			std::cerr << std::setw(6) << m_dim4 << std::setw(8) 
				<< other.m_dim4 << '\n'; 
			
			return false;
		}
		return true;
	}

	template <typename T>
	Vector4D<T> Vector4D<T>::operator+(const Vector4D& other) const
	{
		// Safety check to make sure the two Vector4D's are the same 
		// shape.
		bool same_shape {check_same_shape(other)};

		if (same_shape)
		{
			// Use the standard library transform function to apply the
			// addition function and store into a new vector (ret_data). 
			std::vector<T> ret_data (m_data.size());
			std::transform(m_data.begin(), m_data.end(), 
				other.m_data.begin(), ret_data.begin(), std::plus<T>());
			return {ret_data, m_dim1, m_dim2, m_dim3, m_dim4};	
		}
		else
		{
			std::cerr << "Error! The two Vector4D's are not the " 
				<< "same size."
				<< " See previous error message for details.\n";
			return {};
		}
	}

	template <typename T>
	void Vector4D<T>::move_into_data(Vector4D<T>& vec)
	{
		// Copy over the dimensions
		m_dim1 = vec.get_dim1();
		m_dim2 = vec.get_dim2();
		m_dim3 = vec.get_dim3();
		m_dim4 = vec.get_dim4();

		// Move the data into this
		m_data = std::move(vec.m_data);
	}

	template <typename T>
	Vector3D<T> Vector4D<T>::slice_dim4(int dim_index)
	{
		Vector3D<T> vec {m_dim1, m_dim2, m_dim3};
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

	template <typename T>
	std::vector<std::vector<std::vector<std::vector<T>>>> 
		Vector4D<T>::as_4d() const
	{

		// Size check
		if (std::ssize(m_data) != m_dim1 * m_dim2 * m_dim3 * m_dim4)
		{
			std::cerr << "Error! The number of elements does not match" 
			<< " the specified dimensions.\n";	
		}

		// Create an actual 4D vector, fill in one element at a time. Must
		// initialize vector with the dimensions. Incredibly messy.
		std::vector<std::vector<std::vector<std::vector<T>>>> vec(
			m_dim1, std::vector<std::vector<std::vector<T>>>(
			m_dim2, std::vector<std::vector<T>>(
			m_dim3, std::vector<T>(m_dim4, 0.0))));
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
}

// This is a pecularity of splitting a class declarations and definitions 
// a header and source file. We need to instatiate the templates with the 
// needed definitions so the linker can see them. Hurts flexibility, but not
// an issue since it's pretty straightforward.
template class Vectors::Vector3D<int>;
template class Vectors::Vector3D<double>;
template class Vectors::Vector4D<int>;
template class Vectors::Vector4D<double>;
