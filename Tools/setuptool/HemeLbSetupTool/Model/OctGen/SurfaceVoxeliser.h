// -*- mode: c++; -*-
#ifndef HEMELBSETUPTOOL_SURFACEVOXELISER_H
#define HEMELBSETUPTOOL_SURFACEVOXELISER_H

#include "Vector.h"
#include "Oct.h"
#include "TriTree.h"

typedef std::pair<Vector, Vector> BBox;

// This class creates an octree with leaf nodes that make a minimal,
// 26-separable voxelisation of the surface represented by the
// triangle mesh, following Huang et al
// (http://dx.doi.org/10.1109/SVV.1998.729593)

// The leaf nodes have added to them a set listing the triangles that
// they could intersect (triIds)
class SurfaceVoxeliser {
public:
  typedef size_t IndexT;
  SurfaceVoxeliser(const std::vector<Vector>& points, const std::vector<Index>& triangles, const std::vector<Vector>& normals);
  BBox AABB_Tri(IndexT iTri);
  
  // Flag all voxels lying within Rc = sqrt(3)/2 of the point
  void FilterPoint(IndexT iPt, const std::vector<Index>& voxels, std::vector<bool>& mask);
  // Flag all voxels lying within Rc = sqrt(3)/2 of the line segment
  void FilterEdge(IndexT iPt, IndexT jPt, const std::vector<Index>& voxels, std::vector<bool>& mask);
  // Mark as inside all points within the triangular prism defined by
  // the following 5 planes (see Huang Fig. 12)
  void FilterTriangle(IndexT iTri, const std::vector<Index>& voxels, std::vector<bool>& mask);

  // For the triangles attached at tri_level on the inTree create the
  // voxel nodes to the outTree
  void DoSubTree(TriTree::Node& inTree, TriTree::Node& outTree);
  
  TriTree operator()(TriTree& inTree, const int tri_level);
  
private:
  const std::vector<Vector>& points;
  const std::vector<Index>& triangles;
  const std::vector<Vector>& normals;
};

#endif
