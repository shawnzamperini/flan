/**
* @file vectors.h
*
* @brief Header file for vectors.cpp
*/
#ifndef VECTORS_H
#define VECTORS_H

#include <iomanip>
#include <iostream>
#include <vector>

#ifndef __CUDACC__
#include "mpi.h"
#endif

namespace Vectors
{
	/**
	* 3D array of values, stored as a flattened 1D vector.
	*
	* A value is accessed with the parentheses operator, e.g., vec(i, j, k).
	* The data is stored in row-major order (k varies fastest), so it can be
	* copied straight to and from GPU arrays. See vectors.cpp for the
	* member function docs.
	*/
	template <typename T>
	class Vector3D
	{
	private:
		std::vector<T> m_data;	///< Values in flattened row-major order
		int m_dim1 {};			///< Size of the first dimension
		int m_dim2 {};			///< Size of the second dimension
		int m_dim3 {};			///< Size of the third dimension
		int m_size {};			///< Number of elements at construction. Not
								///< updated by resize or move_into_data.

	public:

		Vector3D();

		Vector3D(int dim1, int dim2, int dim3);

		Vector3D(std::vector<T> data_in, int dim1, int dim2, int dim3);

		Vector3D(std::vector<T> data_in, int dim1, int dim2, int dim3,
				int dim4);

		int get_dim1() const;

		int get_dim2() const;

		int get_dim3() const;

		int get_size() const;

		const std::vector<T>& get_data() const;

		std::vector<T>& get_data();

		int calc_index(const int i, const int j, const int k) const;

		std::tuple<int, int, int> get_ijk(const int idx);

		int get_i(const int idx) const;

		int get_j(const int idx) const;

		int get_k(const int idx) const;

		T& operator()(int i, int j, int k);

		const T& operator()(const int i, const int j,
			const int k) const;

		void move_into_data(Vectors::Vector3D<T>& vec);

		void resize(const int dim1, const int dim2, const int dim3);

#ifndef __CUDACC__
		void broadcast(MPI_Comm comm);
#endif

	};

	/**
	* 4D array of values, stored as a flattened 1D vector.
	*
	* A value is accessed with the parentheses operator, e.g.,
	* vec(i, j, k, l). The data is stored in row-major order (l varies
	* fastest), so it can be copied straight to and from GPU arrays. The
	* background and statistics arrays are Vector4Ds with dimensions
	* (t, x, y, z). See vectors.cpp for the member function docs.
	*/
	template <typename T>
	class Vector4D
	{
	private:
		std::vector<T> m_data {};	///< Values in flattened row-major order
		int m_dim1 {};				///< Size of the first dimension
		int m_dim2 {};				///< Size of the second dimension
		int m_dim3 {};				///< Size of the third dimension
		int m_dim4 {};				///< Size of the fourth dimension
		int m_size {};				///< Number of elements at construction.
									///< Not updated by resize,
									///< move_into_data or move assignment.

	public:

		Vector4D();

		Vector4D(int dim1, int dim2, int dim3, int dim4);

		Vector4D(std::vector<T> data_in, int dim1, int dim2, int dim3,
			int dim4);

		Vector4D(std::vector<T> data_in, int dim1, int dim2, int dim3);

		Vector4D(const Vector4D& v);

		// We could alternatively prevent copies with delete.
		//Vector4D(const Vector4D& v) = delete;

		int get_dim1() const;

		int get_dim2() const;

		int get_dim3() const;

		int get_dim4() const;

		const std::vector<T>& get_data() const;

		std::vector<T>& get_data();

		int calc_index(const int i, const int j, const int k,
			const int l) const;

		T& operator()(const int i, const int j, const int k, const int l);

		const T& operator()(const int i, const int j, const int k, const int l)
			const;

		Vector4D& operator=(Vector4D&& other) noexcept;

		bool check_same_shape(const Vector4D& other) const;

		Vector4D operator+(const Vector4D& other) const;

		Vector4D operator*(const Vector4D& other) const;

		//void move_into_data(Vectors::Vector4D<T>& vec);
		void move_into_data(Vector4D<T>&& vec);
		void move_into_data(Vector4D<T>& vec);

		Vector3D<T> slice_dim1(int dim_index);

		Vector3D<T> slice_dim4(int dim_index);

		std::vector<std::vector<std::vector<std::vector<T>>>> as_4d() const;

        void resize(const int dim1, const int dim2, const int dim3,
            const int dim4);

#ifndef __CUDACC__
		void broadcast(MPI_Comm comm);
#endif
	};
}

#endif
