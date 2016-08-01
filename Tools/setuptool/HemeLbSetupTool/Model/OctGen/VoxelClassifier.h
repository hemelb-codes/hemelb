// -*- mode: c++; -*-
#ifndef HEMELBSETUPTOOL_VOXELCLASSIFIER_H
#define HEMELBSETUPTOOL_VOXELCLASSIFIER_H

#include "Vector.h"
#include "Oct.h"
#include "TriTree.h"

struct LinkBase {
	typedef std::array<LinkBase*, 26> LinkArray;
	LinkArray links;
};
class Vox : public LinkBase {

};
class Cut : public LinkBase {

};


typedef Octree<Vox> VoxTree;

class VoxelClassifier {
public:
  VoxelClassifier(const std::vector<Vector>& points, const std::vector<Index>& triangles,
		  const std::vector<Vector>&normals, const std::vector<int>& labels);
  
  double IntersectLinkWithTriangle(const Index& coord, const Index& disp, const int iTri) const;
  VoxTree::NodePtr DoRegion(const TriTree::NodePtr node);
private:
  const std::vector<Vector>& Points;
  const std::vector<Index>& Triangles;
  const std::vector<Vector>& Normals;
  const std::vector<int>& Labels;
};
#endif
