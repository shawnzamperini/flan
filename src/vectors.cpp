/**
* @file vectors.cpp
*
* @brief Implementations of multidimensional vectors
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
	/**
	* @class Vector3D
	* @brief Implementation of a 3D vector
	*
	* This class provides an intuitive interface to 3D data. A value at a 
	* particular index is accessed by indexing via the parentheses operator,
	* e.g., vec(i,j,k). The underlying data is stored as a contiguous 1D
	* vector. 
	*/
	
	/**
	* @brief Default empty constructor
	*/
	template <typename T>
	Vector3D<T>::Vector3D()
	{}

	/**
	* @brief Constructor to create an empty vector of indicated size
	*/
	template <typename T>
	Vector3D<T>::Vector3D(int dim1, int dim2, int dim3)
		: m_dim1 {dim1}
		, m_dim2 {dim2}
		, m_dim3 {dim3}
	{
		m_size = dim1 * dim2 * dim3;
		m_data.resize(m_size);
	}

	/**
	* @brief Constructor to create vector and fill with data passed in
	*/
	template <typename T>
	Vector3D<T>::Vector3D(std::vector<T> data_in, int dim1, int dim2, 
		int dim3)
		: Vector3D<T>::Vector3D(dim1, dim2, dim3)
	{
		// Copy over into the class
		// Note: Expensive! Could we use move semantics instead?
		m_data = data_in;
	}

	/**
	* @brief Artificial constructor that should not be directly called.
	*
	* Artificial constructor to allow us to compile read_gkyl::load_values.
	* That function selectively returns a brace-enclosed initializer list
	* for either a Vector3D or Vector4D based on the template parameter.
	* The compiler nonetheless needs to see a constructor for Vector3D
	* that accepts 5 inputs, since that is what Vector4D uses. This
	* constructor should never get used.
	*/
	template <typename T>
	Vector3D<T>::Vector3D(std::vector<T> data_in, int dim1, int dim2, int dim3, 
			int dim4)
	{
		std::cerr << "Error! Vector3D was called with 4 dimensions. This"
			<< " is a programming error. Please fix\n";
	}

	/**
	* @brief Get first dimension size
	* @return Returns dimension as an int
	*/
	template <typename T>
	int Vector3D<T>::get_dim1() {return m_dim1;}

	/**
	* @brief Get second dimension size
	* @return Returns dimension as an int
	*/
	template <typename T>
	int Vector3D<T>::get_dim2() {return m_dim2;}

	/**
	* @brief Get third dimension size
	* @return Returns dimension as an int
	*/
	template <typename T>
	int Vector3D<T>::get_dim3() {return m_dim3;}

	/**
	* @brief Get size of underlying data
	* @return Returns size as an int
	*/
	template <typename T>
	int Vector3D<T>::get_size() {return m_size;}

	/**
	* @brief Get the 1D vector containing all the data
	* @return Returns a 1D vector of underlying data type
	*/
	template <typename T>
	std::vector<T> Vector3D<T>::get_data() {return m_data;}
	
	/**
	* @brief Overload parentheses to act as indexing.
	* @return Returns reference to value represented by this index of 
	* indicated data type
	*/
	template <typename T>
	T& Vector3D<T>::operator()(int i, int j, int k)
	{
		int index = i * (m_dim2 * m_dim3) + j * m_dim3 + k;
		return m_data[index];
	}

	/**
	* @brief Move the passed in vector to data. 
	* 
	* This uses move semantics. Hmmm... this is really some sort of 
	* abstraction of operator=...
	*/
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
	

	/**
	* @class Vector4D
	* @brief Implementation of a 4D vector
	*
	* This class provides an intuitive interface to 4D data. A value at a 
	* particular index is accessed by indexing via the parentheses operator,
	* e.g., vec(i,j,k,l). The underlying data is stored as a contiguous 1D
	* vector. 
	*/

	/**
	* @brief Default constructor
	*/
	template <typename T>
	Vector4D<T>::Vector4D()
	{}

	/**
	* @brief Constructor to allow creating an empty vector
	*/
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

	/**
	* @brief Constructor to create vector and fill with the data passed on
	*/
	template <typename T>
	Vector4D<T>::Vector4D(std::vector<T> data_in, int dim1, int dim2, int dim3, 
		int dim4)
		: Vector4D<T>::Vector4D(dim1, dim2, dim3, dim4)
	{
		// Copy over into the class
		// Note: Expensive! Could we use move semantics instead?
		m_data = data_in;
	}

	/**
	* @brief Artificial constructor that should not be directly called.
	* 
	* Artificial constructor to allow us to compile read_gkyl.load_values.
	* That function selectively returns a brace-enclosed initializer list
	* for either a Vector3D or Vector4D based on the template parameter.
	* The compiler nonetheless needs to see a constructor for Vector4D
	* that accepts these inputs, since that is what Vector3D uses. This
	* constructor should never get used.
	*/
	template <typename T>
	Vector4D<T>::Vector4D(std::vector<T> data_in, int dim1, int dim2, int dim3)
	{
		std::cerr << "Error! Vector4D was called with 3 dimensions. This"
			<< " is a programming error. Please fix\n";
	}

	/**
	* @brief Copy constructor
	*
	* We can make the copy constructor yell at us since we probably never 
	* want to make a copy of a Vector4D, with the exception of when a Vector4D
	* is involved in an OpenMP reduction (in that case you need to make a copy
	* for each thread. Useful to uncomment out this print statement to check 
	* for inefficiencies. 
	*/
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
	
	/**
	* @brief Get first dimension size
	* @return Returns dimension as an int
	*/
	template <typename T>
	int Vector4D<T>::get_dim1() const {return m_dim1;}

	/**
	* @brief Get second dimension size
	* @return Returns dimension as an int
	*/
	template <typename T>
	int Vector4D<T>::get_dim2() const {return m_dim2;}

	/**
	* @brief Get third dimension size
	* @return Returns dimension as an int
	*/
	template <typename T>
	int Vector4D<T>::get_dim3() const {return m_dim3;}

	/**
	* @brief Get fourth dimension size
	* @return Returns dimension as an int
	*/
	template <typename T>
	int Vector4D<T>::get_dim4() const {return m_dim4;}

	/**
	* @brief Get a reference to the underlying 1D data
	* @return Returns a reference to the underlying 1D vector
	*/
	template <typename T>
	const std::vector<T>& Vector4D<T>::get_data() const {return m_data;}
}

// This is a pecularity of splitting a class declarations and definitions 
// a header and source file. We need to instatiate the templates with the 
// needed definitions so the linker can see them. Hurts flexibility, but not
// an issue since it's pretty straightforward.
template class Vectors::Vector3D<int>;
template class Vectors::Vector3D<double>;
template class Vectors::Vector4D<int>;
template class Vectors::Vector4D<double>;
