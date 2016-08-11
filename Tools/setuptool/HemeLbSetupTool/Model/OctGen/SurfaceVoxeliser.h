// -*- mode: c++; -*-
#ifndef HEMELBSETUPTOOL_SURFACEVOXELISER_H
#define HEMELBSETUPTOOL_SURFACEVOXELISER_H

#include <array>

#include "Vector.h"
#include "Oct.h"
#include "TriTree.h"
#include "Cgal.h"
#include "CgalSearch.h"

struct Cut {
	Cut() : dist(std::numeric_limits<double>::infinity()), id(-1) {
	}
	double dist;
	int id;
};
struct Vox {
	std::array<Cut, 26> closest_cut;
};

typedef std::shared_ptr<Vox> VoxPtr;
typedef Octree<VoxPtr> VoxTree;

class SurfaceVoxeliser {
public:
  SurfaceVoxeliser(const int node_size, const std::vector<Vector>& points, const std::vector<Index>& triangles,
		  const std::vector<Vector>&normals, const std::vector<int>& labels);

  // We maybe need a destructor because the order of destruction of mesh and searcher seems to matter?
  //~SurfaceVoxeliser();
  
  VoxTree::NodePtr ComputeIntersectionsForRegion(const TriTree::Node& node);
private:
  const static int NDIR = 13;
  const std::vector<Vector>& Points;
  const std::vector<Index>& Triangles;
  const std::vector<Vector>& Normals;
  const std::vector<int>& Labels;

  CgalMeshPtr mesh;
  std::unique_ptr<CgalSearchTree> searcher;

  typedef std::pair<Index, Index> Seg;

  std::vector<int> direction_indices;
  std::vector<Index> directions;
  std::vector<std::vector<Seg>> relative_segments;
};
#endif
