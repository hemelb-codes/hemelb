// -*- mode: c++; -*-
#ifndef HEMELBSETUPTOOL_SURFACEVOXELISER_H
#define HEMELBSETUPTOOL_SURFACEVOXELISER_H

#include <array>
#include <memory>

#include "Vector.h"
#include "Oct.h"
#include "TriTree.h"
#include "Cgal.h"
#include "CgalSearch.h"
#include "FluidSiteTree.h"
#include "Iolet.h"

struct Cut {
  Cut() : dist(std::numeric_limits<double>::infinity()), id(-1) {
  }
  double dist;
  int id;
};
// This represents a site (whether fluid or solid) that has one or more
// links cutting the surface.
struct EdgeSite {
  std::array<Cut, 26> closest_cut;
};

typedef std::shared_ptr<EdgeSite> EdgeSitePtr;
typedef Octree<EdgeSitePtr> EdgeSiteTree;

class SurfaceVoxeliser {
public:
  SurfaceVoxeliser(const int node_size, const std::vector<Vector>& points, const std::vector<Index>& triangles,
		   const std::vector<Vector>&normals, const std::vector<int>& labels, const std::vector<Iolet>& iolets);

  // We maybe need a destructor because the order of destruction of mesh and searcher seems to matter?
  //~SurfaceVoxeliser();
  
  EdgeSiteTree::NodePtr ComputeIntersectionsForRegion(TriTree::ConstNodePtr node) const;
  FluidTree::NodePtr ClassifyRegion(EdgeSiteTree::ConstNodePtr node) const;

  FluidTree operator()(const TriTree& inTree, const TriTree::Int tri_level);
private:
  const static int NDIR = 13;
  const std::vector<Vector>& Points;
  const std::vector<Index>& Triangles;
  const std::vector<Vector>& Normals;
  const std::vector<int>& Labels;
  const std::vector<Iolet>& Iolets;
  
  CgalMeshPtr mesh;
  std::unique_ptr<CgalSearchTree> searcher;

  typedef std::pair<Index, Index> Seg;

  std::vector<int> direction_indices;
  std::vector<Index> directions;
  std::vector<std::vector<Seg>> relative_segments;
};
#endif
