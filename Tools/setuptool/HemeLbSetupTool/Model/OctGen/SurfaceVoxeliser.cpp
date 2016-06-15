#include "SurfaceVoxeliser.h"

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
