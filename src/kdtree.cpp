/**
* @file kdtree.cpp
*/
#include <memory>

#include "kdtree.h"
#include "vectors.h"

namespace KDTree
{
	// The KDTree stores a reference to the PointCloud, and so it must remain
	// in scope the entire time. So we define it here with empty vectors and
	// it'll be modified later.
	std::vector<double> empty_vec {0};
	PointCloud<double> cloud {&empty_vec, &empty_vec, &empty_vec};

	// Build KDTree and return it for use in nearest_neighbor later
	std::unique_ptr<KDTree_t> build_tree(const Vectors::Vector3D<double>& X,
		const Vectors::Vector3D<double>& Y, 
		const Vectors::Vector3D<double>& Z)
	{
		// Create a point cloud, passing in the addresses to the underlying
		// flattened 1D data.
		cloud.X = &(X.get_data());
		cloud.Y = &(Y.get_data()); 
		cloud.Z = &(Z.get_data());

		// Build KD-Tree. {dims, cloud, max leaf}. Not sure what max leaf does
		// yet. I am using index here because that's what the examples called
		// it. Maybe I'm mixing terminology up here in a bad way, idk.
		constexpr int max_leaf_size {10};
		auto index = std::make_unique<KDTree_t>(3, cloud, 
			nanoflann::KDTreeSingleIndexAdaptorParams(max_leaf_size));
		index->buildIndex();
		return index;
	}

	// Find nearest neighbor to a point in tree.
	std::size_t nearest_neighbor(std::unique_ptr<KDTree_t>& tree, 
		const double X_pt, const double Y_pt, const double Z_pt)
	{

		// From pointcloud_kdd_radius.cpp
		double query_pt[3] = {X_pt, Y_pt, Z_pt};
		size_t                num_results = 1;
        std::vector<uint32_t> ret_index(num_results);
        std::vector<double>    out_dist_sqr(num_results);

        num_results = tree->knnSearch(
            &query_pt[0], num_results, &ret_index[0], &out_dist_sqr[0]);

        // In case of less points in the tree than requested:
        //ret_index.resize(num_results);
        //out_dist_sqr.resize(num_results);

		// Return
		//return nearest_index;
		return ret_index[0];
	}
}
