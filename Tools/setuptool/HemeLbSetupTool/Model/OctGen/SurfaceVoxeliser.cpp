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
    auto dr = vox - p;
    auto dr2 = dr.GetMagnitudeSquared();
    *mask_i = dr2 < 0.75;
    ++mask_i;
  }
  return mask;
}

