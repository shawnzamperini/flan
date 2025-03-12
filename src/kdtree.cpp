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
		auto index = std::make_unique<KDTree_t>(3, cloud, 
			nanoflann::KDTreeSingleIndexAdaptorParams(10));
		index->buildIndex();
		return index;
	}

	// Find nearest neighbor to a point in tree.
	std::size_t nearest_neighbor(std::unique_ptr<KDTree_t>& tree, 
		const double X_pt, const double Y_pt, const double Z_pt)
	{
		// Variables for nearest neighbor search
		std::size_t nearest_index;
		double distance_sq;
		double query_pt[3] = {X_pt, Y_pt, Z_pt};

		// Perform nearest neighbor search
		std::cout << "result_set...\n";
		nanoflann::KNNResultSet<double> result_set {1}; // 1 nearest neighbor
		std::cout << "init...\n";
		result_set.init(&nearest_index, &distance_sq);
		std::cout << "findNeighbors... " << X_pt << " " << Y_pt << " " 
			<< Z_pt << "\n";
		tree->findNeighbors(result_set, &query_pt[0]);

		// Return
		return nearest_index;
	}
}
