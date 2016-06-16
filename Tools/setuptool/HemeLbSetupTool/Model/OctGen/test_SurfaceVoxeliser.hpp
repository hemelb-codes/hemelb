#ifndef HEMELBSETUPTOOL_TEST_SURFACEVOXELISER_HPP
#define HEMELBSETUPTOOL_TEST_SURFACEVOXELISER_HPP

#include <cppunit/extensions/HelperMacros.h>
#include "TestResources/Meshes.hpp"
#include "SurfaceVoxeliser.h"
#include "TriangleSorter.h"
#include <array>
#include <memory>
#include <stack>

#include "range.hpp"
#include "enumerate.hpp"

class SurfaceVoxeliserTests : public CppUnit::TestFixture {
  CPPUNIT_TEST_SUITE(SurfaceVoxeliserTests);
  
  CPPUNIT_TEST(TrivialPoints);
  CPPUNIT_TEST(TrivialEdges);
  CPPUNIT_TEST(TrivialPlane);
  CPPUNIT_TEST(Trivial);
  CPPUNIT_TEST(TrivialBig);
  CPPUNIT_TEST(SphereTri90);
  CPPUNIT_TEST(Connected1);
  CPPUNIT_TEST(Connected2);
  CPPUNIT_TEST(Sphere);
  CPPUNIT_TEST(DuctPoint);
  CPPUNIT_TEST(DuctEdge);
  CPPUNIT_TEST(Duct);
  CPPUNIT_TEST_SUITE_END();
  
public:
  void TrivialPoints() {
    // 8 cube
    auto levels = 3;
    auto tri_level = 2;
    auto n = 1 << levels;
    auto triv = SimpleMeshFactory::MkTrivial();
    auto tree = TrianglesToTreeSerial(levels, tri_level, triv->points, triv->triangles);
    SurfaceVoxeliser voxer( triv->points, triv->triangles, triv->normals);
    
    // Get full grid of points
    std::vector<Index> coords(n * n * n);
    auto cursor = coords.begin();
    for (auto i: range(n))
      for (auto j: range(n))
	for (auto k: range(n)) {
	  *cursor = Index(i,j,k);
	  ++cursor;
	}
    // sanity check 
    CPPUNIT_ASSERT(coords[3*n*n + 4*n + 5] == Index(3,4,5));

    for (auto iPt: range(triv->points.size())) {
      std::vector<bool> mask(coords.size());
      voxer.FilterPoint(iPt, coords, mask);
      // innies
      for (auto i: range(coords.size())) {
	auto dr2 = (Vector(coords[i]) - triv->points[iPt]).GetMagnitudeSquared();
	if (mask[i]) {
	  // innie
	  CPPUNIT_ASSERT(dr2 <= 3.0/4.0);
	} else {
	  // outie
	  CPPUNIT_ASSERT(dr2 > 3.0/4.0);
	}
      }
      
    }
  }
  
  template <class T>
  T hypot2(T x, T y) {
    return x*x + y*y;
  }
  
  void TrivialEdges() {
    // 8 cube
    auto levels = 3;
    auto tri_level = 2;
    auto n = 1 << levels;
    auto triv = SimpleMeshFactory::MkTrivial();
    auto tree = TrianglesToTreeSerial(levels, tri_level, triv->points, triv->triangles);
    SurfaceVoxeliser voxer( triv->points, triv->triangles, triv->normals);
    
    // Get full grid of points
    std::vector<Index> coords(n * n * n);
    auto cursor = coords.begin();
    for (auto i: range(n))
      for (auto j: range(n))
	for (auto k: range(n)) {
	  *cursor = Index(i,j,k);
	  ++cursor;
	}
    // sanity check 
    CPPUNIT_ASSERT(coords[3*n*n + 4*n + 5] == Index(3,4,5));
    const std::array<std::array<int, 2>,3> edges = {{{0,1}, {1,2}, {2,0}}};
    for (auto pair: enumerate(edges)) {
      auto n = pair.first;
      auto i = pair.second[0];
      auto j = pair.second[1];

      std::vector<bool> mask(coords.size());
      voxer.FilterEdge(i, j, coords, mask);
      if (n == 0) {
	// cylinder along z axis
	for (auto i: range(coords.size())) {
	  auto v = coords[i];
	  auto r2 = hypot2(v.x - 1.2, v.y - 1.2);
	  if (mask[i]) {
	    // innie
	    CPPUNIT_ASSERT(v.z >= 1.2);
	    CPPUNIT_ASSERT(v.z <= 2.2);
	    CPPUNIT_ASSERT(r2 <= 3.0/4.0);
	  } else {
	    // outie
	    CPPUNIT_ASSERT(v.z < 1.2 || v.z > 2.2 || r2 > 0.75);
	  }
	}
      }

      if (n == 2) {
	// y cylinder
	for (auto i: range(coords.size())) {
	  auto v = coords[i];
	  auto r2 = hypot2(v.x - 1.2, v.z - 1.2);
	  if (mask[i]) {
	    // innies
            CPPUNIT_ASSERT(v.y >= 1.2);
            CPPUNIT_ASSERT(v.y <= 2.2);
	    CPPUNIT_ASSERT(r2 <= 3.0/4.0);
	  } else {
	    // outies
	    CPPUNIT_ASSERT(v.y < 1.2 || v.y > 2.2 || r2 > 3.0/4.0);
	  }
	}
      }
    }
  }

  void TrivialPlane() {
    // 8 cube
    auto levels = 3;
    auto tri_level = 2;
    auto n = 1 << levels;
    auto triv = SimpleMeshFactory::MkTrivial();
    auto tree = TrianglesToTreeSerial(levels, tri_level, triv->points, triv->triangles);
    SurfaceVoxeliser voxer( triv->points, triv->triangles, triv->normals);
    
    // Get full grid of points
    std::vector<Index> coords(n * n * n);
    auto cursor = coords.begin();
    for (auto i: range(n))
      for (auto j: range(n))
	for (auto k: range(n)) {
	  *cursor = Index(i,j,k);
	  ++cursor;
	}
    // sanity check 
    CPPUNIT_ASSERT(coords[3*n*n + 4*n + 5] == Index(3,4,5));
    
    // This tri shouldn't get any points
    std::vector<bool> mask(coords.size());
    voxer.FilterTriangle(0, coords, mask);
    for (auto flag: mask)
      CPPUNIT_ASSERT(!flag);

    // This one should get (1,2,2)
    std::vector<bool> pt122_mask(coords.size());
    voxer.FilterTriangle(1, coords, pt122_mask);
    for (auto i: range(coords.size())) {
      if (coords[i].x == 1 && coords[i].y == 2 && coords[i].z == 2)
	CPPUNIT_ASSERT(pt122_mask[i]);
      else
	CPPUNIT_ASSERT(!pt122_mask[i]);
    }
  }

  void Trivial() {
    // 16 cube
    auto levels = 4;
    auto tri_level = 2;
    auto triv = SimpleMeshFactory::MkTrivial();
    auto tree = TrianglesToTreeSerial(levels, tri_level, triv->points, triv->triangles);
    SurfaceVoxeliser voxer(triv->points, triv->triangles, triv->normals);
    auto vox_tree = voxer(tree, tri_level);

    std::vector<Index> expected_nodes = {{1,1,1},
  					 {1,1,2},
  					 {1,2,1},
  					 {1,2,2},
					 
  					 {2,1,1},
  					 {2,1,2},
  					 {2,2,1},
  					 {2,2,2},
					 
  					 {1,1,3},
  					 {1,2,3},
  					 {1,3,1},
  					 {1,3,2}};

    vox_tree.IterDepthFirst([&](TriTree::Node& node) {
	Index ijk(node.X(), node.Y(), node.Z());
        if (node.Level() == 0) {
    	  // leaf
	  auto ijk_it = std::find(expected_nodes.begin(), expected_nodes.end(), ijk);
    	  CPPUNIT_ASSERT(ijk_it != expected_nodes.end());
    	  expected_nodes.erase(ijk_it);
    	} else if (node.Level() == 1) {
	  CPPUNIT_ASSERT(ijk.x == 0 || ijk.x == 2);
	  CPPUNIT_ASSERT(ijk.y == 0 || ijk.y == 2);
	  CPPUNIT_ASSERT(ijk.z == 0 || ijk.z == 2);
	} else {
	  // All higher levels only in the zero octant
	  CPPUNIT_ASSERT(ijk.x == 0);
	  CPPUNIT_ASSERT(ijk.y == 0);
	  CPPUNIT_ASSERT(ijk.z == 0);
    	} 
      });
    CPPUNIT_ASSERT(expected_nodes.empty());
  }

  void TrivialBig() {
    // We want to check that this copes with triangles much bigger
    // than a tri_level bounding box.
    
    // 16 cube
    auto levels = 4;
    auto tri_level = 2;
    auto triv = SimpleMeshFactory::MkTrivial();
    
    // tri_level has nodes 4 voxels to a side
    // trivial's BBox is (1.2, 1.2, 1.2) - (1.2, 2.2, 2.2)
    for (auto& p: triv->points) {
      p -= 1.2;
      p *= 4.0;
      p += 2.2;
    }
    // Now (2.2, 2.2, 2.2) - (2.2, 6.2, 6.2)
    
    auto tree = TrianglesToTreeSerial(levels, tri_level, triv->points, triv->triangles);
    SurfaceVoxeliser voxer(triv->points, triv->triangles, triv->normals);
    auto vox_tree = voxer(tree, tri_level);

    vox_tree.IterDepthFirst(0, 0,
			    [](TriTree::Node& vox) {
			      CPPUNIT_ASSERT(vox.X() == 2 || vox.X() == 3);
			    });
  }
  void SphereTri90() {
    // Initial testing showed that this fails for some triangles, including this
    TriTree::Int levels = 5;
    TriTree::Int tri_level = 3;
    
    auto sphere = SimpleMeshFactory::MkSphere();
    auto tree = TrianglesToTreeSerial(levels, tri_level, sphere->points, sphere->triangles);
    SurfaceVoxeliser voxer(sphere->points, sphere->triangles, sphere->normals);
    
    auto iTri = 90;
    auto tri_pt_ids = sphere->triangles[iTri];
    // tri_pts = points[tri_pt_ids]
    auto norm = sphere->normals[iTri];
    
    auto bb = voxer.AABB_Tri(iTri);
    // Round down
    auto lo = Index(bb.first) - 1;
    // 1 to round up, 1 to be safe, 1 to get range upper bound
    auto hi = Index(bb.second) + 3;
    
    // This region of interest bounds the triangle safely
    // auto shape = hi - lo;
    // auto size = shape.x * shape.y * shape.z;
    
    // std::vector<Index> voxels(size);
    // auto flat_i = 0;
    std::vector<Index> voxels;
    for (auto i: range(lo.x, hi.x))
      for (auto j: range(lo.y, hi.y))
	for (auto k: range(lo.z, hi.z))
	  voxels.push_back({i, j, k});
	  
    std::vector<bool> union_mask(voxels.size());
    
    std::vector<std::vector<Index>> expected_voxels = {
      {{22,8,18}, {22,9,18}, {23,9,18}},
      {{22,8,13}, {22,9,13}, {23,9,13}},
      {{25,15,13}, {25,16,13}}
    };
    for (auto i: range(3)) {
      auto iPt = tri_pt_ids[i];
      std::vector<bool> inside_mask(voxels.size());
      voxer.FilterPoint(iPt, voxels, inside_mask);
      auto& ev = expected_voxels[i];
      for(auto iVox: range(voxels.size())) {
	auto ptr = std::find(ev.begin(), ev.end(), voxels[iVox]);
	if(inside_mask[iVox]) {
	  // Flagged as included, so must be in expected_voxels
	  CPPUNIT_ASSERT(ptr != ev.end());
	  // Remove
	  ev.erase(ptr);
	} else {
	  // Not included, so it must NOT be in expected_voxels
	  CPPUNIT_ASSERT(ptr == ev.end());
	}
      }
      // Must be no leftover voxels
      CPPUNIT_ASSERT(ev.empty());
      // OR this point's flags onto the union
      for(auto iFlag: range(inside_mask.size()))
	union_mask[i] = union_mask[i] | inside_mask[i];
      
    }
    // # Can't face working this out by hand!
    // for i in xrange(3):
    //     inside_mask[:] = False
    //     iPt = tri_pt_ids[i]
    //     jPt = tri_pt_ids[(i+1)%3]
    //     voxer.FilterEdge(iPt, jPt, voxels, inside_mask)
    //     union |= inside_mask
    
    // # or this...
    // inside_mask[:] = False
    // voxer.FilterTriangle(iTri, voxels, inside_mask)
    // union |= inside_mask
  }

  template<class T>
  struct Ar3 {
    // quick n dirty 3d array
      
    const Index shape;
    const Index strides;
    const int size;

    // shared_ptr => rule of zero
    std::shared_ptr<T> data;

    Ar3(const Index& shp) : shape(shp),
			    strides(shp.y*shp.z, shp.z, 1),
			    size(shp.x*shp.y*shp.z) {
      data = std::shared_ptr<T>(new T[size], std::default_delete<T[]>());
    }
    Ar3(int x, int y, int z) : Ar3(Index(x,y,z)) {
      
    }
    
    void fill(const T& val) {
      std::fill(data.get(), data.get() + size, val);
    }
    
    T& operator()(int i, int j, int k) {
      return operator[](Index(i,j,k));
    }
    T& operator[](const Index& ijk) {
      auto i = Index::Dot(ijk, strides);
      return data.get()[i];
    }
  };

  // Helper function that does a 6-connected flood fill to mark the
  // continuous region with the same value as idx.
  Ar3<char> ConnectedRegion(Ar3<char>& array, Index idx) {
    const auto& shape = array.shape;
    Ar3<char> ans(shape);
    Ar3<char> checked(shape);
    ans.fill(0);
    checked.fill(0);
          
    auto target = array[idx];
    
    const std::array<Index, 6> deltas = {{
	{ 1, 0, 0},
	{ 0, 1, 0},
	{ 0, 0, 1},
	{-1, 0, 0},
	{ 0,-1, 0},
	{ 0, 0,-1}      
      }};
    std::stack<Index> stack;
    stack.push(idx);
    
    while (!stack.empty()) {
      auto xyz = stack.top();
      stack.pop();
      ans[xyz] = 1;
      checked[xyz] = 1;
      for (auto& delta: deltas) {
	auto uvw = xyz + delta;
	if (uvw.x < 0 || uvw.x >= shape.x)
	  continue;
	if (uvw.y < 0 || uvw.y >= shape.y)
	  continue;
	if (uvw.z < 0 || uvw.z >= shape.z)
	  continue;
	if (!checked[uvw] && array[uvw] == target)
	  stack.push(uvw);
      }

    }

    return ans;
  }

  void Connected1() {
    Ar3<char> im(3,5,5);
    const char raw[] = {
                   0,0,0,0,0,
                   0,1,1,1,0,
                   0,1,1,1,0,
                   0,1,1,1,0,
                   0,0,0,0,0,
                   
                   0,0,0,0,0,
		   0,1,1,1,0,
		   0,1,1,1,0,
		   0,1,1,1,0,
		   0,0,0,0,0,
		   
                   0,0,0,0,0,
		   0,1,1,1,0,
		   0,1,1,1,0,
		   0,1,1,1,0,
		   0,0,0,0,0
    };
    std::memcpy(im.data.get(), raw, im.size);
    auto ff = ConnectedRegion(im, Index(0,2,2));
    for (auto i: range(3))
      for (auto j: range(5))
	for (auto k: range(5)) {
	  CPPUNIT_ASSERT(ff(i,j,k) == im(i,j,k));
	}
  }

  void Connected2() {
    Ar3<char> im(1,5,5);
    const char raw[] = {
      1,1,0,0,0,
      1,1,0,0,0,
      0,0,0,0,0,
      0,0,0,1,1,
      0,0,0,1,1
    };
    std::memcpy(im.data.get(), raw, im.size);
    auto ff = ConnectedRegion(im, Index(0,0,0));

    Ar3<char> expected(1,5,5);
    const char expected_raw[] = {
      1,1,0,0,0,
      1,1,0,0,0,
      0,0,0,0,0,
      0,0,0,0,0,
      0,0,0,0,0
    };
    std::memcpy(expected.data.get(), expected_raw, expected.size);
    for (auto i: range(1))
      for (auto j: range(5))
	for (auto k: range(5)) {
	  CPPUNIT_ASSERT(ff(i,j,k) == expected(i,j,k));
	}
  }
  
  Ar3<char> TreeToMaskArray(TriTree& tree) {
    auto cube_size = 1 << tree.Level();
    Ar3<char> cube(cube_size, cube_size, cube_size);
    cube.fill(0);
    
    tree.IterDepthFirst(0,0, [&](TriTree::Node& node) {
	// Index ijk(node.X(), node.Y(), node.Z());
	// std::cout << ijk << "," << std::endl;
        cube(node.X(), node.Y(), node.Z()) = 1;
      });
    
    return cube;
  }
  
  void Sphere() {
    TriTree::Int levels = 5;
    auto box = 1 << levels;
    TriTree::Int tri_level = 3;
    
    auto sphere = SimpleMeshFactory::MkSphere();
    auto tree = TrianglesToTreeSerial(levels, tri_level, sphere->points, sphere->triangles);
    SurfaceVoxeliser voxer(sphere->points, sphere->triangles, sphere->normals);

    auto vox_tree = voxer(tree, tri_level);
    
    auto edge_mask = TreeToMaskArray(vox_tree);
    auto interior = ConnectedRegion(edge_mask, Index(15,15,15));
    
    for (auto i: range(box))
      for (auto j: range(box))
	for (auto k: range(box)) {
	  Index ijk(i,j,k);
	  if (interior[ijk]) {
	    auto r2 = (Vector(ijk) - 15.5).GetMagnitudeSquared();
	    CPPUNIT_ASSERT(r2 < 100.0);
	  }
	}

    // Make sure we've seen all the triangles
    // First gather the list of them
    IdList seen_tris;
    vox_tree.IterDepthFirst(0,0, [&seen_tris](const TriTree::Node& node) {
	seen_tris.insert(boost::container::ordered_unique_range_t(),
			 node.Data().begin(), node.Data().end());
      });

    // Now assert that the list == 0...nTris-1
    auto nTri = sphere->triangles.size();
    int i = 0;
    for (auto triId: seen_tris)
      CPPUNIT_ASSERT(triId == i++);
  }
  
  void DuctPoint() {
    // 16 cube
    TriTree::Int levels = 4;
    TriTree::Int tri_level = 2;
    auto duct = SimpleMeshFactory::MkDuct();
    auto tree = TrianglesToTreeSerial(levels, tri_level, duct->points, duct->triangles);
    SurfaceVoxeliser voxer(duct->points, duct->triangles, duct->normals);
    
    // first 4 points all on lower face, x = 0.5/3.5, y = 0.5/3.5, z = 2.5
    // Box is (0, 0, 2) -- (4,4,3) inclusive
    std::vector<Index> voxels;
    for (auto i: range(0,5))
      for (auto j: range(0,5))
	for (auto k: range(2,4))
	  voxels.push_back({i,j,k});
    
    std::vector<bool> mask(voxels.size());
    voxer.FilterPoint(0, voxels, mask);
    // All 8 voxels around the point should be included
    CPPUNIT_ASSERT(mask[10*0 + 2*0 + 0]);
    CPPUNIT_ASSERT(mask[10*0 + 2*0 + 1]);
    CPPUNIT_ASSERT(mask[10*0 + 2*1 + 0]);
    CPPUNIT_ASSERT(mask[10*0 + 2*1 + 1]);
    CPPUNIT_ASSERT(mask[10*1 + 2*0 + 0]);
    CPPUNIT_ASSERT(mask[10*1 + 2*0 + 1]);
    CPPUNIT_ASSERT(mask[10*1 + 2*1 + 0]);
    CPPUNIT_ASSERT(mask[10*1 + 2*1 + 1]);

    std::fill(mask.begin(), mask.end(), false);
    voxer.FilterPoint(2, voxels, mask);
    // All 8 voxels around the point should be included
    CPPUNIT_ASSERT(mask[10*3 + 2*3 + 0]);
    CPPUNIT_ASSERT(mask[10*3 + 2*3 + 1]);
    CPPUNIT_ASSERT(mask[10*3 + 2*4 + 0]);
    CPPUNIT_ASSERT(mask[10*3 + 2*4 + 1]);
    CPPUNIT_ASSERT(mask[10*4 + 2*3 + 0]);
    CPPUNIT_ASSERT(mask[10*4 + 2*3 + 1]);
    CPPUNIT_ASSERT(mask[10*4 + 2*4 + 0]);
    CPPUNIT_ASSERT(mask[10*4 + 2*4 + 1]);
  }

  void DuctEdge() {
    // 16 cube
    TriTree::Int levels = 4;
    TriTree::Int tri_level = 2;
    auto duct = SimpleMeshFactory::MkDuct();
    auto tree = TrianglesToTreeSerial(levels, tri_level, duct->points, duct->triangles);
    SurfaceVoxeliser voxer(duct->points, duct->triangles, duct->normals);

    // Seen some problems with the edge (3.5, 3.5, 2.5) - (3.5, 3.5, 14.5)
    // or edge (2, 6)
    std::vector<Index> voxels;
    for (auto i: range(3,5))
      for (auto j: range(3,5))
	for (auto k: range(2,16))
	  voxels.push_back({i,j,k});
    
    std::vector<bool> mask(voxels.size());
    voxer.FilterEdge(2, 6, voxels, mask);

    for (auto i: range(3,5))
      for (auto j: range(3,5))
	for (auto k: range(2,16)) {
	  auto flag = mask[(i-3)*28 + (j-3)*14 + (k-2)];
	  if (k == 2 || k == 15)
	    CPPUNIT_ASSERT(flag == false);
	  else
	    CPPUNIT_ASSERT(flag == true);
	}
  }
  
  void Duct() {
    // 16 cube
    TriTree::Int levels = 4;
    TriTree::Int tri_level = 2;
    auto duct = SimpleMeshFactory::MkDuct();
    auto tree = TrianglesToTreeSerial(levels, tri_level, duct->points, duct->triangles);
    SurfaceVoxeliser voxer(duct->points, duct->triangles, duct->normals);
    auto vox_tree = voxer(tree, tri_level);
    
    auto edge_mask = TreeToMaskArray(vox_tree);
    CPPUNIT_ASSERT(edge_mask(0,1,2));
    auto interior = ConnectedRegion(edge_mask, Index(2,2,10));

    auto box = 1 << levels;
    for (auto i: range(box))
      for (auto j: range(box))
	for (auto k: range(box)) {
	  Index ijk(i,j,k);
	  if (interior[ijk]) {
	    CPPUNIT_ASSERT(i == 2);
	    CPPUNIT_ASSERT(j == 2);
	  }
	}
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION(SurfaceVoxeliserTests);

#endif
