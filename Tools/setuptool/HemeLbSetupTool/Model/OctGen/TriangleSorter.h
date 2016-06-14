// -*- mode: c++; -*-
#ifndef HEMELBSETUPTOOL_TRIANGLESORTER_H
#define HEMELBSETUPTOOL_TRIANGLESORTER_H

#include "Oct.h"
#include <vector>
#include "Vector.h"
#include <boost/container/flat_set.hpp>

typedef boost::container::flat_set<int> IdList;
typedef Octree<IdList> TriTree;

// Create an Octree with n_levels potential number of levels.
// Nodes will be created down to the level tri_level but only if voxel-voxel
// links originating within that node could potentially intersect the triangles
// specified by (points, triangles), where:
// points is a numpy array of doubles with shape == (nPoints, 3) specifying the vertices.
// triangles is a numpy array of ints with shape == (nTris, 3) specifiying the IDs of the verticies comprising the triangle.
// The create nodes at tri_level will have an attribute 'triIds' added which is
// a list of the IDs of the triangles that might intersect its voxel-voxel links.

TriTree TrianglesToTreeSerial(const int n_levels, const int tri_level, const std::vector<Vector>& points, const std::vector<Index>& triangles);

TriTree TrianglesToTreeParallel(const int n_levels, const int tri_level, const std::vector<Vector>& points, const std::vector<Index>& triangles, int nprocs);

TriTree TrianglesToTree_Worker(const int n_levels, const int tri_level,
			       const std::vector<Vector>& points,
			       std::vector<Index>::const_iterator triStart, std::vector<Index>::const_iterator triEnd,
			       const int tri_index_start);


// This guy accumulates multiple trees together, merging their id lists
class TreeSummer {
public:
  TreeSummer(TriTree::Int levels, TriTree::Int tri_level);
  void Add(TriTree& source);
  TriTree& GetTree();
  
private:
  const TriTree::Int levels;
  const TriTree::Int tri_level;
  // This will be the output
  TriTree tree;
  
};

#endif
