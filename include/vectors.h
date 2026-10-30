/**
* @file vectors.h
*
* @brief Header file for vectors.cpp
*/
#ifndef VECTORS_H
#define VECTORS_H

#include <cstddef>
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

	/**
	* Values on a grid split into zones, stored as a single flattened 1D
	* vector.
	*
	* Each zone is a logically rectangular block of nx x ny x nz cells
	* surrounded by a layer of ghost cells ghost wide on every side, and has
	* a leading dimension of size dim1 (the number of background times, or 1
	* for time-independent quantities like the geometry). A value is
	* accessed with vec(zone, t, i, j, k), or vec(zone, i, j, k) when dim1
	* is 1, where i, j and k are interior cell indices that may extend into
	* the ghost cells, i.e., -ghost <= i < nx + ghost.
	*
	* The zones are stored one after the other, each in row-major order
	* (k varies fastest). With a single zone and no ghost cells the layout
	* is identical to a Vector4D(dim1, nx, ny, nz). See vectors.cpp for the
	* member function docs.
	*/
	template <typename T>
	class ZonedVector
	{
	private:
		std::vector<T> m_data {};			///< Values, zone after zone
		int m_dim1 {};						///< Size of the leading dimension
		int m_ghost {};						///< Ghost cells on each side
		std::vector<int> m_nx {};			///< Interior cells in x per zone
		std::vector<int> m_ny {};			///< Interior cells in y per zone
		std::vector<int> m_nz {};			///< Interior cells in z per zone
		std::vector<std::size_t> m_offset {};	///< Start of each zone
		std::vector<std::size_t> m_stride_t {};	///< Stride of dim1 per zone
		std::vector<std::size_t> m_stride_i {};	///< Stride of i per zone
		std::vector<std::size_t> m_stride_j {};	///< Stride of j per zone

		void calc_layout();

	public:

		ZonedVector();

		ZonedVector(int dim1, const std::vector<int>& nx,
			const std::vector<int>& ny, const std::vector<int>& nz,
			int ghost = 0);

		explicit ZonedVector(Vector4D<T>&& vec);

		explicit ZonedVector(Vector3D<T>&& vec);

		int get_nzones() const;

		int get_dim1() const;

		int get_ghost() const;

		int get_nx(const int zone) const;

		int get_ny(const int zone) const;

		int get_nz(const int zone) const;

		std::size_t get_zone_offset(const int zone) const;

		std::size_t get_zone_size(const int zone) const;

		const std::vector<T>& get_data() const;

		std::vector<T>& get_data();

		std::size_t calc_index(const int zone, const int t, const int i,
			const int j, const int k) const;

		T& operator()(const int zone, const int t, const int i, const int j,
			const int k);

		const T& operator()(const int zone, const int t, const int i,
			const int j, const int k) const;

		T& operator()(const int zone, const int i, const int j, const int k);

		const T& operator()(const int zone, const int i, const int j,
			const int k) const;

		bool check_same_shape(const ZonedVector& other) const;

#ifndef __CUDACC__
		void broadcast(MPI_Comm comm);
#endif
	};
}

#endif
