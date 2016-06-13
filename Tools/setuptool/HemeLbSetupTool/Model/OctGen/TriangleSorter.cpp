#include "TriangleSorter.h"

TriTree TrianglesToTree(const int n_levels, const int tri_level, const std::vector<Vector>& points, const std::vector<Index>& triangles) {
  return TrianglesToTree_Worker(n_levels, tri_level, points, triangles.begin(), triangles.end(), 0);
}
  
TriTree TrianglesToTree_Worker(const int n_levels, const int tri_level,
			       const std::vector<Vector>& points,
			       std::vector<Index>::const_iterator triStart, std::vector<Index>::const_iterator triEnd,
			       const int tri_index_start) {
  typedef TriTree::Int Int;
  Int cube_size = 1 << n_levels;
  Int tri_box_size = 1 << tri_level;
  
  std::vector<Int> non_halo_edges, lowers, uppers;
  for (Int i = 0; i< cube_size; i += tri_box_size) {
    non_halo_edges.push_back(i);
    // i==0 use lower=0 cos we know that no tri can be outside our bounding box.
    // The upper and lower bounds of the regions of influence of the tri_level nodes.
    // +/- 1 because we have to consider the halo.
    lowers.push_back(i ? i - 1 : 0);
    uppers.push_back(i + tri_box_size);
  }
  
  TriTree tree(n_levels);
  //for (auto i_tri = 0; i_tri < ntris; ++i_tri) {
  for (auto tri_it = triStart; tri_it != triEnd; ++tri_it) {
    auto triPtIds = *tri_it;
    
    // bounding box of the triangle
    Vector min(std::numeric_limits<double>::infinity());
    Vector max(0);
    for (auto iPt: triPtIds) {
      min.UpdatePointwiseMin(points[iPt]);
      max.UpdatePointwiseMax(points[iPt]);
    }
    
    Index range_min;
    Index range_max;
    for (auto iDim=0; iDim<3; ++iDim) {
      // iterator pointing to the first element > min, then go back one
      auto lower_ptr = std::upper_bound(lowers.cbegin(), lowers.cend(), min[iDim]) - 1;
      // iterator pointing to the first element > max, then go forward one
      auto upper_ptr = std::upper_bound(uppers.cbegin(), uppers.cend(), max[iDim]) + 1;
      // convert to indices
      range_min[iDim] = lower_ptr - lowers.cbegin();
      range_max[iDim] = upper_ptr - uppers.cbegin();
    }
    
    // This bit is the only non-thread-safe part.
    // We will later merge the trees if this was run in parallel
    for (auto i = range_min[0]; i < range_max[0]; ++i)
      for (auto j = range_min[1]; j < range_max[1]; ++j)
	for (auto k = range_min[2]; k < range_max[2]; ++k) {
	  auto node = tree.GetCreate(non_halo_edges[i], non_halo_edges[j], non_halo_edges[k], tri_level);
	  node->Data().push_back(tri_it - triStart + tri_index_start);
	}
  }
  return tree;
}
