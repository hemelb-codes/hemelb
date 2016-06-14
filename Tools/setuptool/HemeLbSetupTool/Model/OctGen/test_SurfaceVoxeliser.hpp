#ifndef HEMELBSETUPTOOL_TEST_SURFACEVOXELISER_HPP
#define HEMELBSETUPTOOL_TEST_SURFACEVOXELISER_HPP

#include <cppunit/extensions/HelperMacros.h>
#include "TestResources/Meshes.hpp"
#include "SurfaceVoxeliser.h"
#include "range.hpp"

class SurfaceVoxeliserTests : public CppUnit::TestFixture {
  CPPUNIT_TEST_SUITE(SurfaceVoxeliserTests);
  CPPUNIT_TEST(TrivialPoints);
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
      auto mask = voxer.FilterPoint(iPt, coords);
      // innies
      for (auto i: range(coords.size())) {
	auto dr2 = (coords[i] - triv->points[iPt]).GetMagnitudeSquared();
	if (dr2 > 0.75) {
	  // outie
	  CPPUNIT_ASSERT(!mask[i]);
	} else {
	  // innie
	  CPPUNIT_ASSERT(mask[i]);
	}
      }
      
    }
  }
};
CPPUNIT_TEST_SUITE_REGISTRATION(SurfaceVoxeliserTests);

#endif
