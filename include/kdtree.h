/**
* @file kdtree.h
*
* This file contains a class KDTree that is a wrapped for interfacing with
* the nanoflann (name is a coincidence) library. It uses a KDTree nearest
* neighbor search to find the nearest point to a given point in a set of
* of coordinates. 
*/
#ifndef KDTREE_H
#define KDTREE_H

#include <memory>
#include <vector>

#include "nanoflann.hpp"
#include "vectors.h"


namespace KDTree
{
	/**
	* @brief 3D PointCloud structure
	*/
	template <typename T>
	struct PointCloud {

		// Pointers to existing X,Y,Z data in vectors
		const std::vector<T>* X {};
		const std::vector<T>* Y {};
		const std::vector<T>* Z {};

		// Return number of data points. They should all be the same, so
		// just use X here.
		inline size_t kdtree_get_point_count() const {return (*X).size();}

		// Return the dim'th dimension of the idx'th point
		inline T kdtree_get_pt(const size_t idx, const size_t dim) const 
		{
			if (dim == 0) return (*X)[idx];
			else if (dim == 1) return (*Y)[idx];
			else return (*Z)[idx];
		}

		// Optional: Bounding box (not used here)
		template <class BBOX>
		bool kdtree_get_bbox(BBOX & /*bb*/) const { return false; }

		// Constructor to assign pointers
		PointCloud(std::vector<T>* Xptr, std::vector<T>* Yptr, 
			std::vector<T>* Zptr)
			: X {Xptr}
			, Y {Yptr}
			, Z {Zptr}
		{}
	};

	/**
	* @brief Alias for KD-Tree with three dimensions (X,Y,Z)
	*/
	using KDTree_t = nanoflann::KDTreeSingleIndexAdaptor<
		nanoflann::L2_Simple_Adaptor<double, PointCloud<double>>,
		PointCloud<double>,
		3 /* dimension */
	>;

	/**
	* @brief Build KDTree from the given X,Y,Z coordinates
	*/
	std::unique_ptr<KDTree_t> build_tree(const Vectors::Vector3D<double>& X,
		const Vectors::Vector3D<double>& Y, 
		const Vectors::Vector3D<double>& Z);

	/**
	* @brief Perform nearest neighbor search to X,Y,Z
	*/
	std::size_t nearest_neighbor(std::unique_ptr<KDTree_t>& tree, 
		const double X_pt, const double Y_pt, const double Z_pt);
}

#endif
