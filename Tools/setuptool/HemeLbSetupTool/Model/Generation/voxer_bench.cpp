#include <memory>
#include <chrono>

#include "TestResources/Meshes.hpp"
#include "TriangleSorter.h"
#include "SurfaceVoxeliser.h"
#include "Neighbours.h"

#include "range.hpp"

int main( int argc, char **argv)
{
	TriTree::Int levels = 5;
	TriTree::Int tri_level = 3;

	auto sphere = SimpleMeshFactory::MkSphere();
	auto tree = TrianglesToTreeSerial(levels, tri_level, sphere->points,
			sphere->triangles);
	SurfaceVoxeliser voxer(1 << tri_level, sphere->points, sphere->triangles,
			       sphere->normals, sphere->labels, sphere->iolets);

	auto t0 = std::chrono::high_resolution_clock::now();
	for (auto i: range(10))
		auto edge_tree = voxer(tree, tri_level);

	auto t1 = std::chrono::high_resolution_clock::now();
	auto dt = std::chrono::duration_cast<std::chrono::milliseconds>(t1-t0);
	//std::chrono::milliseconds dt(t1-t0);
	std::cout << dt.count() << " ms" << std::endl;
}
