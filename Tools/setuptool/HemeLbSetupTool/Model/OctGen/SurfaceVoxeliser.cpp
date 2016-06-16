#include "SurfaceVoxeliser.h"
#include <algorithm>
#include "range.hpp"
#include "enumerate.hpp"

SurfaceVoxeliser::SurfaceVoxeliser(const std::vector<Vector>& p, const std::vector<Index>& t, const std::vector<Vector>& n) :
  points(p), triangles(t), normals(n) {
}

// Flag all voxels lying within Rc = sqrt(3)/2 of the point
void SurfaceVoxeliser::FilterPoint(IndexT iPt, const std::vector<Index>& voxels, std::vector<bool>& mask) {
  auto& p = points[iPt];
  
  //auto mask_i = mask.begin();
  for (auto i: range(voxels.size())) {
    auto dr = Vector(voxels[i]) - p;
    auto dr2 = dr.GetMagnitudeSquared();
    mask[i] = mask[i] | (dr2 <= 0.75);
  }
}

// Flag all voxels lying within Rc = sqrt(3)/2 of the line segment
void SurfaceVoxeliser::FilterEdge(IndexT iPt, IndexT jPt, const std::vector<Index>& voxels, std::vector<bool>& mask) {
  /*
   * x .              . b
   *
   *
   *       .  p       
   *
   * a .
   *
   * Our point of interest is x
   * line is a -> b
   * it's equation is a + lambda n
   * where 0 <= lambda <= 1 
   * and n = b - a
   * p is the closest point on the line to x
   * define r = x - a
   *
   * Since xp is perpendicular to ab:
   *     p.n = x.n
   * hence
   *     lambda = (x - a).n / n.n
   * and
   *     p = a + n*(r.n) / (n.n)
   *
   * if rho = |x - p|, Pythagoras tells us:
   *     r.r = rho^2 + (r.n)^2 / (n.n)
   */
  auto& a = points[iPt];
  auto& b = points[jPt];
  
  auto n = b - a;
  auto n2 = n.GetMagnitudeSquared();
  
  for (auto i: range(voxels.size())) {
    auto r = Vector(voxels[i]) - a;
    auto r_n = Vector::Dot(r, n);
    auto lambda = r_n / n2;

    auto r2 = r.GetMagnitudeSquared();
    auto rho2 = r2 - r_n*r_n / n2;
    
    // Inside points have 0 <= lamdba <= 1 and rho^2 < Rc^2
    mask[i] = mask[i] || ((0 <= lambda) && (lambda <= 1) && (rho2 <= 0.75));
  }
}

// Mark as inside all points within the triangular prism defined by
// the following 5 planes (see Huang Fig. 12)
void SurfaceVoxeliser::FilterTriangle(IndexT iTri, const std::vector<Index>& voxels, std::vector<bool>& mask) {
  
  auto& ids = triangles[iTri];
  auto& norm = normals[iTri];
  const std::vector<Vector> tri_points = {points[ids.x],
					  points[ids.y],
					  points[ids.z]};
  // t_26 is defined in Huang figure 10
  double t26 = 0;
  for (auto n: Vector(norm))
    t26 += std::abs(n) / 2;

  // These will hold the data on our 5 planes.
  // NB: we want the *inward* normals
  std::vector<Vector> plane_normals(5);
  std::vector<double> plane_offsets(5);
  
  // These 6 points define the prism
  std::vector<Vector> upper_points(3);
  std::vector<Vector> lower_points(3);
  for (auto i: range(3)) {
    upper_points[i] = tri_points[i] + norm * t26;
    lower_points[i] = tri_points[i] - norm * t26;
  }
  
  // Top and bottom are easy
  plane_normals[0] = norm;
  plane_offsets[0] = Vector::Dot(plane_normals[0], lower_points[0]);
  
  plane_normals[1] = -norm;
  plane_offsets[1] = Vector::Dot(plane_normals[1], upper_points[0]);
  
  // Think about the edges now
  // Will need cross products to get their normals. Our normal can be
  // used as one of the inputs.
  auto v10 = tri_points[1] - tri_points[0];
  auto v21 = tri_points[2] - tri_points[1];
  auto v02 = tri_points[0] - tri_points[2];
  
  plane_normals[2] = Vector::Cross(norm, v10);
  plane_offsets[2] = Vector::Dot(plane_normals[2], upper_points[0]);
  
  plane_normals[3] = Vector::Cross(norm, v21);
  plane_offsets[3] = Vector::Dot(plane_normals[3], upper_points[1]);
  
  plane_normals[4] = Vector::Cross(norm, v02);
  plane_offsets[4] = Vector::Dot(plane_normals[4], upper_points[2]);
  
  // greater than because the normals are inwards
  for (auto pair: enumerate(voxels)) {
    auto i = pair.first;
    Vector vox(pair.second);
    bool tmp = true;
    for (auto j: range(5))
      tmp &= Vector::Dot(vox, plane_normals[j]) >= plane_offsets[j];
    
    mask[i] = mask[i] | tmp;
  }
}

// Return the axis-aligned bounding box for a triangle.
BBox SurfaceVoxeliser::AABB_Tri(IndexT iTri) {
  auto ids = triangles[iTri];
  // bounding box of the triangle
  Vector min(std::numeric_limits<double>::infinity());
  Vector max(0);
  for (auto iPt: ids) {
    min.UpdatePointwiseMin(points[iPt]);
    max.UpdatePointwiseMax(points[iPt]);
  }

  return BBox(min, max);
}

void SurfaceVoxeliser::DoSubTree(TriTree::Node& inTree, TriTree::Node& outTree) {
  const int box_size = 1 << inTree.Level();
  
  for (auto iTri: inTree.Data()) {
    auto ids = triangles[iTri];
    auto bbox = AABB_Tri(iTri);
    // voxels that could in principle intersect
    Index vlo(bbox.first);
    // +2 for use as upper bound in range statement
    Index vhi(bbox.second);
    vhi += 2;
    // now clip this against the node's BB
    vlo.x = std::max(vlo.x, int(inTree.X()));
    vlo.y = std::max(vlo.y, int(inTree.Y()));
    vlo.z = std::max(vlo.z, int(inTree.Z()));

    vhi.x = std::min(vhi.x, inTree.X() + box_size);
    vhi.y = std::min(vhi.y, inTree.Y() + box_size);
    vhi.z = std::min(vhi.z, inTree.Z() + box_size);

    auto vshape = vhi - vlo;
    auto vsize = vshape.x * vshape.y * vshape.z;
    
    // Get voxels
    std::vector<Index> voxels(vsize);
    auto cursor = voxels.begin();
    for (auto i: range(vlo.x, vhi.x))
      for (auto j: range(vlo.y, vhi.y))
	for (auto k: range(vlo.z, vhi.z)) {
	  *cursor = Index(i,j,k);
	  ++cursor;
	}
    // Now apply the three tests:
    // 1) within point sphere?
    // 2) within edge cylinder?
    // 3) within the planes?
    std::vector<bool> mask(vsize);
    
    for (auto i: range(3)) {
      auto iPt = ids[i];
      auto jPt = ids[(i + 1) % 3];
      FilterPoint(iPt, voxels, mask);
      FilterEdge(iPt, jPt, voxels, mask);
    }
    
    FilterTriangle(iTri, voxels, mask);

    // Now add iTri to the included voxels
    for (auto i: range(vsize))
      if (mask[i]) {
	auto& vox = voxels[i];
	auto vNode = outTree.GetCreate(vox.x, vox.y, vox.z, 0);
	vNode->Data().insert(iTri);
      }
	
  }
}


TriTree SurfaceVoxeliser::operator()(TriTree& inTree, const int tri_level) {
  TriTree outTree(inTree.Level());
  inTree.IterDepthFirst(tri_level, tri_level,
			[&](TriTree::Node& inNode){
			  TriTree::NodePtr outNode = outTree.GetCreate(inNode.X(),
								       inNode.Y(),
								       inNode.Z(),
								       inNode.Level());
			  DoSubTree(inNode, *outNode);
			});
  return outTree;
}
