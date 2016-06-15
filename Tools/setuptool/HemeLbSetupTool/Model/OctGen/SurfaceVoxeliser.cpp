#include "SurfaceVoxeliser.h"
#include "range.hpp"
#include "enumerate.hpp"

SurfaceVoxeliser::SurfaceVoxeliser(const std::vector<Vector>& p, const std::vector<Index>& t, const std::vector<Vector>& n) :
  points(p), triangles(t), normals(n) {
}

// Flag all voxels lying within Rc = sqrt(3)/2 of the point
std::vector<bool> SurfaceVoxeliser::FilterPoint(IndexT iPt, const std::vector<Index>& voxels) {
  auto& p = points[iPt];
  
  std::vector<bool> mask(voxels.size());
  auto mask_i = mask.begin();
  for (auto vox: voxels) {
    auto dr = Vector(vox) - p;
    auto dr2 = dr.GetMagnitudeSquared();
    *mask_i = dr2 < 0.75;
    ++mask_i;
  }
  return mask;
}

// Flag all voxels lying within Rc = sqrt(3)/2 of the line segment
std::vector<bool> SurfaceVoxeliser::FilterEdge(IndexT iPt, IndexT jPt, const std::vector<Index>& voxels) {
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
  
  std::vector<bool> mask(voxels.size());
  auto mask_i = mask.begin();
  for (auto vox: voxels) {
    auto r = Vector(vox) - a;
    auto r_n = Vector::Dot(r, n);
    auto lambda = r_n / n2;

    auto r2 = r.GetMagnitudeSquared();
    auto rho2 = r2 - r_n*r_n / n2;
    
    // Inside points have 0 <= lamdba <= 1 and rho^2 < Rc^2
    *mask_i = (0 <= lambda) && (lambda <= 1) && (rho2 <= 0.75);
    ++mask_i;
  }
  return mask;    
}
// Mark as inside all points within the triangular prism defined by
// the following 5 planes (see Huang Fig. 12)
std::vector<bool> SurfaceVoxeliser::FilterTriangle(IndexT iTri, const std::vector<Index>& voxels) {
  
  std::vector<bool> mask(voxels.size());
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
    
    mask[i] = tmp;
  }
  return mask;
}
