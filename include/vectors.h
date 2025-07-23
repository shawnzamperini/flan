/**
* @file vectors.h
*
* @brief Header file for vectors.cpp
*/
#ifndef VECTORS_H
#define VECTORS_H

#include <vector>
#include <iostream>
#include <iomanip>

namespace Vectors
{
	/**
	* @brief Implementation of a 3D vector
	*
	* This class provides an intuitive interface to 3D data. A value at a 
	* particular index is accessed by indexing via the parentheses operator,
	* e.g., vec(i,j,k). The underlying data is stored as a contiguous 1D
	* vector. 
	*/
	template <typename T>
	class Vector3D
	{
	private:
		std::vector<T> m_data;
		int m_dim1 {};
		int m_dim2 {};
		int m_dim3 {};
		int m_size {};

	public:

		/**
		* @brief Default empty constructor
		*/
		Vector3D();

		/**
		* @brief Constructor to create an empty vector of indicated size
		*/
		Vector3D(int dim1, int dim2, int dim3);

		/**
		* @brief Constructor to create vector and fill with data passed in
		*/
		Vector3D(std::vector<T> data_in, int dim1, int dim2, int dim3);

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
		Vector3D(std::vector<T> data_in, int dim1, int dim2, int dim3, 
				int dim4);
	
		/**
		* @brief Get first dimension size
		* @return Returns dimension as an int
		*/
		int get_dim1() const;

		/**
		* @brief Get second dimension size
		* @return Returns dimension as an int
		*/
		int get_dim2() const;

		/**
		* @brief Get third dimension size
		* @return Returns dimension as an int
		*/
		int get_dim3() const;

		/**
		* @brief Get size of underlying data
		* @return Returns size as an int
		*/
		int get_size() const;

		/**
		* @brief Get the 1D vector containing all the data (const)
		* @return Returns a 1D vector of underlying data type
		*/
		const std::vector<T>& get_data() const;

		/**
		* @brief Get the 1D vector containing all the data (non-const)
		* @return Returns a 1D vector of underlying data type
		*/
		std::vector<T>& get_data();

		/**
		* @brief Convert from 3D index to the 1D one used in the underlying 
		* m_data.
		* @return Returns the index used to index the underlying 1D data 
		* vector.
		*/
		int calc_index(const int i, const int j, const int k) const;
	
		/**
		* @brief Get the i,j,k indices from a given index in the underlying
		* 1D data (i.e., inverse of calc_index).
		*/
		std::tuple<int, int, int> get_ijk(const int idx);

		/**
		* @brief Get the index along the first dimension from an index
		* calculated by calc_index.
		*/
		int get_i(const int idx) const;

		/**
		* @brief Get the index along the second dimension from an index
		* calculated by calc_index.
		*/
		int get_j(const int idx) const;

		/**
		* @brief Get the index along the third dimension from an index
		* calculated by calc_index.
		*/
		int get_k(const int idx) const;

		/**
		* @brief Overload parentheses to act as indexing.
		* @return Returns reference to value represented by this index of 
		* indicated data type
		*/
		T& operator()(int i, int j, int k);

		/**
		* @brief Overload parentheses to act as indexing (for const objects).
		* @return Returns reference to value represented by this index of 
		* indicated data type
		*/
		const T& operator()(const int i, const int j, 
			const int k) const;

		/**
		* @brief Move the passed in vector to data. 
		* 
		* This uses move semantics. Hmmm... this is really some sort of 
		* abstraction of operator=...
		*/
		void move_into_data(Vectors::Vector3D<T>& vec);

		/**
		* @brief Resize Vector3D to the passed in dimensions
		*
		* This is needed to add a layer of safety. One could resize m_data
		* themselves, but then they would also need to set dim1, dim2 and dim3.
		* This function takes care of that to prevent that accident.
		*/
		void resize(const int dim1, const int dim2, const int dim3);
	};

	/**
	* @brief Implementation of a 4D vector
	*
	* This class provides an intuitive interface to 4D data. A value at a 
	* particular index is accessed by indexing via the parentheses operator,
	* e.g., vec(i,j,k,l). The underlying data is stored as a contiguous 1D
	* vector. 
	*/
	template <typename T>
	class Vector4D
	{
	private:
		std::vector<T> m_data {};
		int m_dim1 {};
		int m_dim2 {};
		int m_dim3 {};
		int m_dim4 {};
		int m_size {};

	public:

		/**
		* @brief Default constructor
		*/
		Vector4D();

		/**
		* @brief Constructor to allow creating an empty vector
		*/
		Vector4D(int dim1, int dim2, int dim3, int dim4);

		/**
		* @brief Constructor to create vector and fill with the data passed on
		*/
		Vector4D(std::vector<T> data_in, int dim1, int dim2, int dim3, 
			int dim4);

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
		Vector4D(std::vector<T> data_in, int dim1, int dim2, int dim3);

		/**
		* @brief Copy constructor
		*
		* We can make the copy constructor yell at us since we probably never 
		* want to make a copy of a Vector4D, with the exception of when a Vector4D
		* is involved in an OpenMP reduction (in that case you need to make a copy
		* for each thread. Useful to uncomment out this print statement to check 
		* for inefficiencies. 
		*/
		Vector4D(const Vector4D& v);

		// We could alternatively prevent copies with delete.
		//Vector4D(const Vector4D& v) = delete;
	
		/**
		* @brief Get first dimension size
		* @return Returns dimension as an int
		*/
		int get_dim1() const;

		/**
		* @brief Get second dimension size
		* @return Returns dimension as an int
		*/
		int get_dim2() const;

		/**
		* @brief Get third dimension size
		* @return Returns dimension as an int
		*/
		int get_dim3() const;

		/**
		* @brief Get fourth dimension size
		* @return Returns dimension as an int
		*/
		int get_dim4() const;

		/**
		* @brief Get a reference to the underlying 1D data
		* @return Returns a reference to the underlying 1D vector
		*/
		const std::vector<T>& get_data() const;

		/**
		* @brief Convert from 4D index to the 1D one used in the underlying 
		* m_data.
		* @return Returns the index used to index the underlying 1D data 
		* vector.
		*/
		int calc_index(const int i, const int j, const int k, 
			const int l) const;

		/**
		* @brief Overload parentheses to act as indexing.
		* @return Returns value represented by this 4D index
		*/
		T& operator()(const int i, const int j, const int k, const int l);

		/**
		* @brief I guess we need a separate identical overload of the 
		* parentheses overload (?) to handle const objects?
		* @return Returns value represented by this 4D index
		*/
		const T& operator()(const int i, const int j, const int k, const int l) 
			const;

		/**
		* @brief Assignment operator moves data from one Vector4D into this 
		* one.
		* @return Returns vector with other's data moved into it
		*/
		Vector4D& operator=(Vector4D&& other) noexcept;

		/**
		* @brief Helper function that checks if other is of the same shape as 
		* this
		* Vector4D.
		* @return Return true if they're the same shape, false (and an error 
		* message) if not
		*/
		bool check_same_shape(const Vector4D& other) const;

		/**
		* @brief Addition operator just adds two vectors together.
		* @return Returns a new Vector4D that is the two summed together
		*/
		Vector4D operator+(const Vector4D& other) const;
			
		/**
		* @brief Move the passed in vector to data. 
		*
		* Hmmm... this is really some sort of abstraction of operator=...
		*/
		//void move_into_data(Vectors::Vector4D<T>& vec);
		void move_into_data(Vector4D<T>&& vec);
		void move_into_data(Vector4D<T>& vec);

		/**
		* @brief Slice a single dimension (which returns a vector with one 
		* lower dimensionality).
		* @return Returns a Vector3D
		*/
		Vector3D<T> slice_dim1(int dim_index);

		/**
		* @brief Slice a single dimension (which returns a vector with one 
		* lower dimensionality).
		* @return Returns a Vector3D
		*/
		Vector3D<T> slice_dim4(int dim_index);

		/**
		* @brief Represent that data as an actual 4D vector. Unsure if this is
		* needed but doesn't hurt to include.
		*/
		std::vector<std::vector<std::vector<std::vector<T>>>> as_4d() const;

        void resize(const int dim1, const int dim2, const int dim3,
            const int dim4);
	};
}

#endif
