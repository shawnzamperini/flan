/**
* @file vectors.cpp
*
* @brief Implementations of multidimensional vectors
*
* This file contains some simple implementations of multidimensional vectors.
* Instead of storing them as actual multidimensional vectors, they are stored
* as flattened vectors and then () is overloaded to act as indexing. The last
* index varies fastest (row-major order), the same as the flattened arrays
* used on the GPU (see calc_3d_index_cuda and calc_4d_index_cuda).
*
* The bottom of this file includes instantiations of Vector3D and Vector4D
* for the different types needed (e.g., int, double, ...). This is a
* peculiarity needed when splitting the declarations and definitions between
* a header and source file, needed so the linker can find the classes for
* each template type.
*/

#include <algorithm>
#include <iomanip>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

#include "errors.h"
#include "mpi.h"
#include "utilities.h"
#include "vectors.h"

namespace Vectors
{

	// ************************
	// * Vector3D definitions *
	// ************************

	/**
	* Creates an empty Vector3D with all dimensions 0.
	*/
	template <typename T>
	Vector3D<T>::Vector3D()
	{}

	/**
	* Creates a Vector3D of the given shape, with every element value
	* initialized (zero for numbers).
	*
	* @param dim1 Size of the first dimension.
	* @param dim2 Size of the second dimension.
	* @param dim3 Size of the third dimension.
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
	* Creates a Vector3D of the given shape, filled with data_in.
	*
	* The run is stopped with Errors::fatal_error if data_in doesn't have
	* dim1 * dim2 * dim3 values, since indexing would then go out of
	* bounds.
	*
	* @param data_in Values in flattened (row-major) order.
	* @param dim1 Size of the first dimension.
	* @param dim2 Size of the second dimension.
	* @param dim3 Size of the third dimension.
	*/
	template <typename T>
	Vector3D<T>::Vector3D(std::vector<T> data_in, int dim1, int dim2, 
		int dim3)
		: Vector3D<T>::Vector3D(dim1, dim2, dim3)
	{
		if (std::ssize(data_in) != m_size)
		{
			Errors::fatal_error("Vector3D was given " 
				+ std::to_string(data_in.size()) + " values for a shape of " 
				+ std::to_string(dim1) + " x " + std::to_string(dim2) + " x " 
				+ std::to_string(dim3) + ".");
		}

		// Copy over into the class
		// Note: Expensive! Could we use move semantics instead?
		m_data = data_in;
	}

	/**
	* Artificial constructor that should never be called.
	*
	* It lets read_gkyl's load_values compile. That function returns a
	* brace-enclosed initializer list for either a Vector3D or a Vector4D
	* depending on its template parameter, so the compiler needs to see a
	* Vector3D constructor that takes the 5 arguments a Vector4D uses. If
	* it's ever called, the run is stopped with Errors::fatal_error.
	*/
	template <typename T>
	Vector3D<T>::Vector3D(std::vector<T> data_in, int dim1, int dim2, int dim3, 
			int dim4)
	{
		Errors::fatal_error("Vector3D was called with 4 dimensions. This "
			"is a programming error. Please fix.");
	}

	/**
	* @name Dimension accessors
	*
	* Return the size of the first, second or third dimension.
	*/
	///@{
	template <typename T>
	int Vector3D<T>::get_dim1() const {return m_dim1;}

	template <typename T>
	int Vector3D<T>::get_dim2() const {return m_dim2;}

	template <typename T>
	int Vector3D<T>::get_dim3() const {return m_dim3;}
	///@}

	/**
	* Returns the number of elements set when the vector was constructed.
	*
	* Not updated by resize or move_into_data, so it can be out of date
	* after either. get_data().size() is always correct.
	*
	* @return Number of elements at construction.
	*/
	template <typename T>
	int Vector3D<T>::get_size() const {return m_size;}

	/**
	* @name Data accessors
	*
	* Return the underlying flattened data. The non-const version lets the
	* data be filled in or changed in place, e.g., element by element with
	* calc_index.
	*/
	///@{
	template <typename T>
	const std::vector<T>& Vector3D<T>::get_data() const {return m_data;}

	template <typename T>
	std::vector<T>& Vector3D<T>::get_data() {return m_data;}
	///@}
	
	/**
	* Converts a 3D index to the index into the flattened data.
	*
	* @param i Index along the first dimension.
	* @param j Index along the second dimension.
	* @param k Index along the third dimension.
	* @return Index into get_data().
	*/
	template <typename T>
	int Vector3D<T>::calc_index(const int i, const int j, const int k) const
	{
		return i * (m_dim2 * m_dim3) + j * m_dim3 + k;
	}

	/**
	* Converts an index into the flattened data back to a 3D index. The
	* inverse of calc_index.
	*
	* @param idx Index into get_data().
	* @return Tuple of the i, j and k indices.
	*/
	template <typename T>
	std::tuple<int, int, int> Vector3D<T>::get_ijk(const int idx)
	{
		// Calculate each index
		int i {idx / (m_dim2 * m_dim3)};  // int division
		int remainder {idx % (m_dim2 * m_dim3)};
		int j {remainder / m_dim3};
		int k {remainder % m_dim3};

		// Return as tuple
		return {i, j, k};
	}

	/**
	* @name Single index from a flattened index
	*
	* Return the i, j or k index for an index into the flattened data, i.e.,
	* one part of get_ijk.
	*
	* @param idx Index into get_data().
	*/
	///@{
	template <typename T>
	int Vector3D<T>::get_i(const int idx) const
	{
		return static_cast<int>(idx / (m_dim2 * m_dim3));
	}

	template <typename T>
	int Vector3D<T>::get_j(const int idx) const
	{
		return static_cast<int>(idx / m_dim3) % m_dim2;
	}

	template <typename T>
	int Vector3D<T>::get_k(const int idx) const
	{
		return idx % m_dim3;
	}
	///@}
		
	/**
	* @name Element access
	*
	* Return a reference to the element at (i, j, k). There is no bounds
	* checking.
	*
	* @param i Index along the first dimension.
	* @param j Index along the second dimension.
	* @param k Index along the third dimension.
	*/
	///@{
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
	///@}

	/**
	* Moves the data and dimensions of vec into this vector, without
	* copying the data.
	*
	* vec's data is left empty, but its dimensions aren't changed, so it
	* shouldn't be used afterwards.
	*
	* @param vec Vector to move from.
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
	* Resizes the vector to a new shape, updating the dimensions along with
	* the data.
	*
	* Safer than resizing get_data() directly, since that wouldn't update
	* the dimensions. The data is not rearranged, so existing values are
	* only still in the right place if just dim1 changed.
	*
	* @param dim1 New size of the first dimension.
	* @param dim2 New size of the second dimension.
	* @param dim3 New size of the third dimension.
	*/
	template <typename T>
	void Vector3D<T>::resize(const int dim1, const int dim2, const int dim3)
	{
		// When we resize, we need to make sure to update the dimension sizes
		m_data.resize(dim1 * dim2 * dim3);
		m_dim1 = dim1;
		m_dim2 = dim2;
		m_dim3 = dim3;
	}

	/**
	* Broadcasts the vector from rank 0 to every other rank.
	*
	* Must be called by every rank in comm. The other ranks' dimensions and
	* data are overwritten with rank 0's.
	*
	* @param comm MPI communicator to broadcast over.
	*/
	template <typename T>
	void Vector3D<T>::broadcast(MPI_Comm comm)
	{
		// Broadcast dimensions
		MPI_Bcast(&m_dim1, 1, MPI_INT, 0, comm);
		MPI_Bcast(&m_dim2, 1, MPI_INT, 0, comm);
		MPI_Bcast(&m_dim3, 1, MPI_INT, 0, comm);

		// Broadcast underlying 1D data from root (assumed to be rank 0) to 
		// other ranks
		Utilities::mpi_broadcast_vector(m_data, 0, comm);
	}
	
	// ************************
	// * Vector4D definitions *
	// ************************

	/**
	* Creates an empty Vector4D with all dimensions 0.
	*/
	template <typename T>
	Vector4D<T>::Vector4D()
	{}

	/**
	* Creates a Vector4D of the given shape, with every element value
	* initialized (zero for numbers).
	*
	* @param dim1 Size of the first dimension.
	* @param dim2 Size of the second dimension.
	* @param dim3 Size of the third dimension.
	* @param dim4 Size of the fourth dimension.
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
	* Creates a Vector4D of the given shape, filled with data_in.
	*
	* The run is stopped with Errors::fatal_error if data_in doesn't have
	* dim1 * dim2 * dim3 * dim4 values, since indexing would then go out of
	* bounds.
	*
	* @param data_in Values in flattened (row-major) order.
	* @param dim1 Size of the first dimension.
	* @param dim2 Size of the second dimension.
	* @param dim3 Size of the third dimension.
	* @param dim4 Size of the fourth dimension.
	*/
	template <typename T>
	Vector4D<T>::Vector4D(std::vector<T> data_in, int dim1, int dim2, int dim3, 
		int dim4)
		: Vector4D<T>::Vector4D(dim1, dim2, dim3, dim4)
	{
		if (std::ssize(data_in) != m_size)
		{
			Errors::fatal_error("Vector4D was given " 
				+ std::to_string(data_in.size()) + " values for a shape of " 
				+ std::to_string(dim1) + " x " + std::to_string(dim2) + " x " 
				+ std::to_string(dim3) + " x " + std::to_string(dim4) + ".");
		}

		// Copy over into the class
		// Note: Expensive! Could we use move semantics instead?
		m_data = data_in;
	}

	/**
	* Artificial constructor that should never be called.
	*
	* It lets read_gkyl's load_values compile. That function returns a
	* brace-enclosed initializer list for either a Vector3D or a Vector4D
	* depending on its template parameter, so the compiler needs to see a
	* Vector4D constructor that takes the 4 arguments a Vector3D uses. If
	* it's ever called, the run is stopped with Errors::fatal_error.
	*/
	template <typename T>
	Vector4D<T>::Vector4D(std::vector<T> data_in, int dim1, int dim2, int dim3)
	{
		Errors::fatal_error("Vector4D was called with 3 dimensions. This "
			"is a programming error. Please fix.");
	}

	/**
	* Copy constructor. Copies the data and dimensions of v.
	*
	* Copies of a Vector4D are expensive and usually not wanted, except for
	* the private copies made for each thread in an OpenMP reduction.
	* Uncommenting the print statement inside is a quick way to find
	* unwanted copies.
	*
	* @param v Vector to copy.
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
	* @name Dimension accessors
	*
	* Return the size of the first, second, third or fourth dimension.
	*/
	///@{
	template <typename T>
	int Vector4D<T>::get_dim1() const {return m_dim1;}

	template <typename T>
	int Vector4D<T>::get_dim2() const {return m_dim2;}

	template <typename T>
	int Vector4D<T>::get_dim3() const {return m_dim3;}

	template <typename T>
	int Vector4D<T>::get_dim4() const {return m_dim4;}
	///@}

	/**
	* @name Data accessors
	*
	* Return the underlying flattened data. The non-const version lets the
	* data be filled in or changed in place, e.g., element by element with
	* calc_index, or through a raw pointer from get_data().data().
	*/
	///@{
	template <typename T>
	const std::vector<T>& Vector4D<T>::get_data() const {return m_data;}

	template <typename T>
	std::vector<T>& Vector4D<T>::get_data() {return m_data;}
	///@}

	/**
	* Converts a 4D index to the index into the flattened data.
	*
	* @param i Index along the first dimension.
	* @param j Index along the second dimension.
	* @param k Index along the third dimension.
	* @param l Index along the fourth dimension.
	* @return Index into get_data().
	*/
	template <typename T>
	int Vector4D<T>::calc_index(const int i, const int j, const int k, 
		const int l) const
	{
		return i * (m_dim2 * m_dim3 * m_dim4) + j 
			* (m_dim3 * m_dim4) + k * m_dim4 + l;
	}

	/**
	* @name Element access
	*
	* Return a reference to the element at (i, j, k, l). There is no bounds
	* checking.
	*
	* @param i Index along the first dimension.
	* @param j Index along the second dimension.
	* @param k Index along the third dimension.
	* @param l Index along the fourth dimension.
	*/
	///@{
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
	///@}

	/**
	* Move assignment. Moves the data and dimensions of other into this
	* vector, without copying the data.
	*
	* other's data is left empty, but its dimensions aren't changed, so it
	* shouldn't be used afterwards.
	*
	* @param other Vector to move from.
	* @return This vector.
	*/
	template <typename T>
	Vector4D<T>& Vector4D<T>::operator=(Vector4D<T>&& other) 
		noexcept
	{
		// Self-assignment check
		if (this != &other)
		{	
			m_dim1 = other.m_dim1;
			m_dim2 = other.m_dim2;
			m_dim3 = other.m_dim3;
			m_dim4 = other.m_dim4;
			m_data = std::move(other.m_data);
		}
		return *this;
	}

	/**
	* Checks if other has the same shape as this vector.
	*
	* @param other Vector to compare to.
	* @return True if all four dimensions match, false otherwise.
	*/
	template <typename T>
	bool Vector4D<T>::check_same_shape(const Vector4D<T>& other) const
	{
		return m_dim1 == other.m_dim1 && m_dim2 == other.m_dim2 
			&& m_dim3 == other.m_dim3 && m_dim4 == other.m_dim4;
	}

	/**
	* Stops the run with Errors::fatal_error because a and b are not the
	* same shape, showing both shapes side by side.
	*
	* Used by operator+ and operator*.
	*
	* @param a First vector.
	* @param b Second vector.
	* @param op Operator that was called (e.g., "+"), for the message.
	*/
	template <typename T>
	[[noreturn]] void shape_mismatch_error(const Vector4D<T>& a, 
		const Vector4D<T>& b, const std::string& op)
	{
		std::ostringstream msg {};
		msg << "Vector4D's are not the same shape for operator" << op 
			<< ".\n  this   other\n"
			<< std::setw(6) << a.get_dim1() << std::setw(8) << b.get_dim1() 
			<< '\n'
			<< std::setw(6) << a.get_dim2() << std::setw(8) << b.get_dim2() 
			<< '\n'
			<< std::setw(6) << a.get_dim3() << std::setw(8) << b.get_dim3() 
			<< '\n'
			<< std::setw(6) << a.get_dim4() << std::setw(8) << b.get_dim4();
		Errors::fatal_error(msg.str());
	}

	/**
	* Adds two vectors element by element.
	*
	* The run is stopped with Errors::fatal_error if they are not the same
	* shape.
	*
	* @param other Vector to add to this one.
	* @return New vector holding the sum.
	*/
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
		else shape_mismatch_error(*this, other, "+");
	}

	/**
	* Multiplies two vectors element by element.
	*
	* The run is stopped with Errors::fatal_error if they are not the same
	* shape.
	*
	* @param other Vector to multiply this one by.
	* @return New vector holding the product.
	*/
	template <typename T>
	Vector4D<T> Vector4D<T>::operator*(const Vector4D& other) const
	{
		// Safety check to make sure the two Vector4D's are the same 
		// shape.
		bool same_shape {check_same_shape(other)};

		if (same_shape)
		{
			// Use the standard library transform function to apply the
			// multiply function and store into a new vector (ret_data). 
			std::vector<T> ret_data (m_data.size());
			std::transform(m_data.begin(), m_data.end(), 
				other.m_data.begin(), ret_data.begin(), std::multiplies<T>());
			return {ret_data, m_dim1, m_dim2, m_dim3, m_dim4};	
		}
		else shape_mismatch_error(*this, other, "*");
	}

	/**
	* @name Move from another vector
	*
	* Move the data and dimensions of vec into this vector, without copying
	* the data. One overload takes an rvalue (e.g., a temporary) and the
	* other a normal reference; both move. vec's data is left empty, but
	* its dimensions aren't changed, so it shouldn't be used afterwards.
	*
	* @param vec Vector to move from.
	*/
	///@{
	template <typename T>
	void Vector4D<T>::move_into_data(Vector4D<T>&& vec)
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
	///@}

	/**
	* Returns the 3D slice at one index of the first dimension, e.g., one
	* time frame of a (t, x, y, z) array.
	*
	* @param dim_index Index along the first dimension to take.
	* @return Copy of the slice, with dimensions (dim2, dim3, dim4).
	*/
	template <typename T>
	Vector3D<T> Vector4D<T>::slice_dim1(int dim_index)
	{
		// Sticking to traditional indexing variables to help make it clear
		// i = dim1, j = dim2, k = dim3, l = dim4
		Vector3D<T> vec {m_dim2, m_dim3, m_dim4};
		for (int j {}; j < m_dim2; ++j)
		{
			for (int k {}; k < m_dim3; ++k)
			{
				for (int l {}; l < m_dim4; ++l)
				{
					// from calc_index just i = dim_index
					int index = dim_index * (m_dim2 * m_dim3 * m_dim4) + j 
						* (m_dim3 * m_dim4) + k * m_dim4 + l;
					vec(j, k, l) = m_data[index];
				}
			}
		}
		return vec;
	}

	/**
	* Returns the 3D slice at one index of the fourth dimension.
	*
	* @param dim_index Index along the fourth dimension to take.
	* @return Copy of the slice, with dimensions (dim1, dim2, dim3).
	*/
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

	/**
	* Copies the data into a nested std::vector, indexed as
	* vec[i][j][k][l]. Currently unused.
	*
	* Prints the dimensions of the result to stdout. The run is stopped
	* with Errors::fatal_error if the number of elements does not match the
	* dimensions.
	*
	* @return Nested 4D copy of the data.
	*/
	template <typename T>
	std::vector<std::vector<std::vector<std::vector<T>>>> 
		Vector4D<T>::as_4d() const
	{

		// Size check
		if (std::ssize(m_data) != m_dim1 * m_dim2 * m_dim3 * m_dim4)
		{
			Errors::fatal_error("The number of elements in a Vector4D ("
				+ std::to_string(m_data.size()) + ") does not match its "
				"dimensions.");
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

	/**
	* Resizes the vector to a new shape, updating the dimensions along with
	* the data.
	*
	* Safer than resizing get_data() directly, since that wouldn't update
	* the dimensions. The data is not rearranged, so existing values are
	* only still in the right place if just dim1 changed.
	*
	* @param dim1 New size of the first dimension.
	* @param dim2 New size of the second dimension.
	* @param dim3 New size of the third dimension.
	* @param dim4 New size of the fourth dimension.
	*/
	template <typename T>
	void Vector4D<T>::resize(const int dim1, const int dim2, const int dim3,
        const int dim4)
	{
		// When we resize, we need to make sure to update the dimension sizes
		m_data.resize(dim1 * dim2 * dim3 * dim4);
		m_dim1 = dim1;
		m_dim2 = dim2;
		m_dim3 = dim3;
		m_dim4 = dim4;
	}

	/**
	* Broadcasts the vector from rank 0 to every other rank.
	*
	* Must be called by every rank in comm. The other ranks' dimensions and
	* data are overwritten with rank 0's.
	*
	* @param comm MPI communicator to broadcast over.
	*/
	template <typename T>
	void Vector4D<T>::broadcast(MPI_Comm comm)
	{
		// Broadcast dimensions
		MPI_Bcast(&m_dim1, 1, MPI_INT, 0, comm);
		MPI_Bcast(&m_dim2, 1, MPI_INT, 0, comm);
		MPI_Bcast(&m_dim3, 1, MPI_INT, 0, comm);
		MPI_Bcast(&m_dim4, 1, MPI_INT, 0, comm);

		// Broadcast underlying 1D data from root (assumed to be rank 0) to
		// other ranks
		Utilities::mpi_broadcast_vector(m_data, 0, comm);
	}

	// ***************************
	// * ZonedVector definitions *
	// ***************************

	/**
	* Creates an empty ZonedVector with no zones.
	*/
	template <typename T>
	ZonedVector<T>::ZonedVector()
	{}

	/**
	* Creates a ZonedVector with the given zone sizes, with every element
	* value initialized (zero for numbers).
	*
	* nx, ny and nz must have one entry per zone. The run is stopped with
	* Errors::fatal_error if they don't or if any size is negative.
	*
	* @param dim1 Size of the leading dimension, e.g., the number of
	*        background times, or 1 for time-independent values.
	* @param nx Number of interior cells in x for each zone.
	* @param ny Number of interior cells in y for each zone.
	* @param nz Number of interior cells in z for each zone.
	* @param ghost Number of ghost cells on each side of every zone.
	*/
	template <typename T>
	ZonedVector<T>::ZonedVector(int dim1, const std::vector<int>& nx,
		const std::vector<int>& ny, const std::vector<int>& nz, int ghost)
		: m_dim1 {dim1}
		, m_ghost {ghost}
		, m_nx {nx}
		, m_ny {ny}
		, m_nz {nz}
	{
		if (nx.size() != ny.size() || nx.size() != nz.size())
		{
			Errors::fatal_error("ZonedVector was given zone sizes for "
				+ std::to_string(nx.size()) + ", " + std::to_string(ny.size())
				+ " and " + std::to_string(nz.size()) + " zones in x, y and "
				"z. These must be equal.");
		}
		if (dim1 < 0 || ghost < 0
			|| std::ranges::any_of(nx, [](int n) {return n < 0;})
			|| std::ranges::any_of(ny, [](int n) {return n < 0;})
			|| std::ranges::any_of(nz, [](int n) {return n < 0;}))
		{
			Errors::fatal_error("ZonedVector was given a negative size.");
		}

		calc_layout();
		m_data.resize(m_offset.back());
	}

	/**
	* Creates a single-zone ZonedVector without ghost cells by moving in
	* the data of a Vector4D, without copying it.
	*
	* The dimensions of vec are taken as (dim1, nx, ny, nz). The layout is
	* identical, so the data is used as-is. vec's data is left empty.
	*
	* @param vec Vector to move from.
	*/
	template <typename T>
	ZonedVector<T>::ZonedVector(Vector4D<T>&& vec)
		: m_dim1 {vec.get_dim1()}
		, m_ghost {0}
		, m_nx {vec.get_dim2()}
		, m_ny {vec.get_dim3()}
		, m_nz {vec.get_dim4()}
	{
		calc_layout();
		m_data = std::move(vec.get_data());
	}

	/**
	* Creates a single-zone ZonedVector without ghost cells and with
	* dim1 = 1 by moving in the data of a Vector3D, without copying it.
	*
	* The dimensions of vec are taken as (nx, ny, nz). vec's data is left
	* empty.
	*
	* @param vec Vector to move from.
	*/
	template <typename T>
	ZonedVector<T>::ZonedVector(Vector3D<T>&& vec)
		: m_dim1 {1}
		, m_ghost {0}
		, m_nx {vec.get_dim1()}
		, m_ny {vec.get_dim2()}
		, m_nz {vec.get_dim3()}
	{
		calc_layout();
		m_data = std::move(vec.get_data());
	}

	/**
	* Calculates the offset of each zone into the flattened data and the
	* strides used for indexing within each zone.
	*
	* m_offset gets one more entry than there are zones, with the last
	* being the total number of elements. Called whenever the zone sizes
	* change.
	*/
	template <typename T>
	void ZonedVector<T>::calc_layout()
	{
		const std::size_t nzones {m_nx.size()};
		m_offset.assign(nzones + 1, 0);
		m_stride_t.resize(nzones);
		m_stride_i.resize(nzones);
		m_stride_j.resize(nzones);

		for (std::size_t z {}; z < nzones; ++z)
		{
			// Size of each dimension, including ghost cells
			const std::size_t px {static_cast<std::size_t>(m_nx[z]
				+ 2 * m_ghost)};
			const std::size_t py {static_cast<std::size_t>(m_ny[z]
				+ 2 * m_ghost)};
			const std::size_t pz {static_cast<std::size_t>(m_nz[z]
				+ 2 * m_ghost)};

			m_stride_j[z] = pz;
			m_stride_i[z] = py * pz;
			m_stride_t[z] = px * py * pz;
			m_offset[z + 1] = m_offset[z]
				+ static_cast<std::size_t>(m_dim1) * m_stride_t[z];
		}
	}

	/**
	* @name Dimension accessors
	*
	* Return the number of zones, the size of the leading dimension, the
	* ghost cell width, or the number of interior cells of a zone in x, y
	* or z.
	*/
	///@{
	template <typename T>
	int ZonedVector<T>::get_nzones() const {return std::ssize(m_nx);}

	template <typename T>
	int ZonedVector<T>::get_dim1() const {return m_dim1;}

	template <typename T>
	int ZonedVector<T>::get_ghost() const {return m_ghost;}

	template <typename T>
	int ZonedVector<T>::get_nx(const int zone) const {return m_nx[zone];}

	template <typename T>
	int ZonedVector<T>::get_ny(const int zone) const {return m_ny[zone];}

	template <typename T>
	int ZonedVector<T>::get_nz(const int zone) const {return m_nz[zone];}
	///@}

	/**
	* Returns the index into get_data() where a zone's values start.
	*
	* @param zone Zone index.
	* @return Offset of the zone into the flattened data.
	*/
	template <typename T>
	std::size_t ZonedVector<T>::get_zone_offset(const int zone) const
	{
		return m_offset[zone];
	}

	/**
	* Returns the number of values stored for a zone, including ghost
	* cells and every entry of the leading dimension.
	*
	* @param zone Zone index.
	* @return Number of values in the zone.
	*/
	template <typename T>
	std::size_t ZonedVector<T>::get_zone_size(const int zone) const
	{
		return m_offset[zone + 1] - m_offset[zone];
	}

	/**
	* @name Data accessors
	*
	* Return the underlying flattened data. The non-const version lets the
	* data be filled in or changed in place.
	*/
	///@{
	template <typename T>
	const std::vector<T>& ZonedVector<T>::get_data() const {return m_data;}

	template <typename T>
	std::vector<T>& ZonedVector<T>::get_data() {return m_data;}
	///@}

	/**
	* Converts a zone and cell index to the index into the flattened data.
	*
	* i, j and k are interior cell indices and may extend into the ghost
	* cells, e.g., i = -1 is the first ghost cell below the zone in x. In
	* Debug builds, out of range indices stop the run with
	* Errors::fatal_error. Otherwise there is no bounds checking.
	*
	* @param zone Zone index.
	* @param t Index along the leading dimension.
	* @param i Cell index in x.
	* @param j Cell index in y.
	* @param k Cell index in z.
	* @return Index into get_data().
	*/
	template <typename T>
	std::size_t ZonedVector<T>::calc_index(const int zone, const int t,
		const int i, const int j, const int k) const
	{
#ifdef DEBUG
		if (zone < 0 || zone >= get_nzones() || t < 0 || t >= m_dim1
			|| i < -m_ghost || i >= m_nx[zone] + m_ghost
			|| j < -m_ghost || j >= m_ny[zone] + m_ghost
			|| k < -m_ghost || k >= m_nz[zone] + m_ghost)
		{
			Errors::fatal_error("ZonedVector index (" + std::to_string(zone)
				+ ", " + std::to_string(t) + ", " + std::to_string(i) + ", "
				+ std::to_string(j) + ", " + std::to_string(k)
				+ ") is out of range.");
		}
#endif
		return m_offset[zone] + t * m_stride_t[zone]
			+ (i + m_ghost) * m_stride_i[zone]
			+ (j + m_ghost) * m_stride_j[zone] + (k + m_ghost);
	}

	/**
	* @name Element access
	*
	* Return a reference to the element at (zone, t, i, j, k). The versions
	* without t are for vectors with dim1 = 1 and use t = 0. See calc_index
	* for the allowed index ranges.
	*
	* @param zone Zone index.
	* @param t Index along the leading dimension.
	* @param i Cell index in x.
	* @param j Cell index in y.
	* @param k Cell index in z.
	*/
	///@{
	template <typename T>
	T& ZonedVector<T>::operator()(const int zone, const int t, const int i,
		const int j, const int k)
	{
		return m_data[calc_index(zone, t, i, j, k)];
	}

	template <typename T>
	const T& ZonedVector<T>::operator()(const int zone, const int t,
		const int i, const int j, const int k) const
	{
		return m_data[calc_index(zone, t, i, j, k)];
	}

	template <typename T>
	T& ZonedVector<T>::operator()(const int zone, const int i, const int j,
		const int k)
	{
		return m_data[calc_index(zone, 0, i, j, k)];
	}

	template <typename T>
	const T& ZonedVector<T>::operator()(const int zone, const int i,
		const int j, const int k) const
	{
		return m_data[calc_index(zone, 0, i, j, k)];
	}
	///@}

	/**
	* Checks whether another ZonedVector has the same number of zones, zone
	* sizes, ghost width and leading dimension.
	*
	* @param other Vector to compare against.
	* @return True if the shapes are the same.
	*/
	template <typename T>
	bool ZonedVector<T>::check_same_shape(const ZonedVector<T>& other) const
	{
		return m_dim1 == other.m_dim1 && m_ghost == other.m_ghost
			&& m_nx == other.m_nx && m_ny == other.m_ny
			&& m_nz == other.m_nz;
	}

	/**
	* Broadcasts the vector from rank 0 to every other rank.
	*
	* Must be called by every rank in comm. The other ranks' zone sizes and
	* data are overwritten with rank 0's, and their layout is recalculated.
	*
	* @param comm MPI communicator to broadcast over.
	*/
	template <typename T>
	void ZonedVector<T>::broadcast(MPI_Comm comm)
	{
		// Broadcast dimensions
		MPI_Bcast(&m_dim1, 1, MPI_INT, 0, comm);
		MPI_Bcast(&m_ghost, 1, MPI_INT, 0, comm);
		Utilities::mpi_broadcast_vector(m_nx, 0, comm);
		Utilities::mpi_broadcast_vector(m_ny, 0, comm);
		Utilities::mpi_broadcast_vector(m_nz, 0, comm);

		// The offsets and strides follow from the dimensions
		calc_layout();

		// Broadcast underlying 1D data from root (assumed to be rank 0) to
		// other ranks
		Utilities::mpi_broadcast_vector(m_data, 0, comm);
	}
}  // namespace Vectors

// This is a peculiarity of splitting a class's declarations and definitions 
// between a header and source file. We need to instantiate the templates with
// the needed definitions so the linker can see them. Hurts flexibility, but not
// an issue since it's pretty straightforward.
template class Vectors::Vector3D<int>;
template class Vectors::Vector3D<float>;
template class Vectors::Vector3D<double>;
template class Vectors::Vector4D<int>;
template class Vectors::Vector4D<float>;
template class Vectors::Vector4D<double>;
template class Vectors::ZonedVector<int>;
template class Vectors::ZonedVector<float>;
template class Vectors::ZonedVector<double>;
