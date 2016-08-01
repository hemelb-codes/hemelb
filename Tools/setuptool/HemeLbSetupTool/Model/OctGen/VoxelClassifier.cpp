#include "VoxelClassifier.h"

VoxelClassifier::VoxelClassifier(const std::vector<Vector>& p,
				 const std::vector<Index>& t,
				 const std::vector<Vector>&n,
				 const std::vector<int>& l) :
  Points(p), Triangles(t), Normals(n), Labels(l) {
}

double VoxelClassifier::IntersectLinkWithTriangle(const Index& coord, const Index& disp, const int iTri) const {
  auto tri_pt_ids = Triangles[iTri];
  // # Work relative to the voxel coord
  const std::vector<const Vector> tri_points = {
    Points[tri_pt_ids[0]] - coord,
    Points[tri_pt_ids[1]] - coord,
    Points[tri_pt_ids[2]] - coord
  };
  auto& norm = Normals[iTri];
  // assume that the line intersects the triangle's plane at a fraction t 
  // along the displacement vector
  double t = Vector::Dot(tri_points[0], norm) / Vector::Dot(disp, norm);
  // To intersect, t must be in (0,1)
  if (t < 0 or t > 1)
    return std::numeric_limits<double>::infinity();
  
  // so the intersection is on our line segment, is it in the triangle?
  // first find the point of intersections
  auto r = disp * t;
  // get coords relative to one point of the tri
  const auto v10 = tri_points[1] - tri_points[0];
  const auto v20 = tri_points[2] - tri_points[0];
  const auto vr0 =  r - tri_points[0];

  const auto v10_v10 = v10.GetMagnitudeSquared();
  const auto v20_v20 = v20.GetMagnitudeSquared();
  const auto v10_v20 = Vector::Dot(v10, v20);
        
  const auto denom = v10_v10 * v20_v20 - v10_v20 * v10_v20;
  
  const auto vr0_v10 = Vector::Dot(vr0, v10);
  const auto vr0_v20 = Vector::Dot(vr0, v20);
  // compute the barycentric coords of the intersection
  const auto u = (v20_v20 * vr0_v10 - v10_v20 * vr0_v20) / denom;
  const auto v = (v10_v10 * vr0_v20 - v10_v20 * vr0_v10) / denom;
  
  if (u>=0 && v>=0 && u+v<=1)
    return t;
  
  return std::numeric_limits<double>::infinity();
}
