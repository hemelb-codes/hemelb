#ifndef HEMELBSETUPTOOL_TEST_TRIANGLESORTER_HPP
#define HEMELBSETUPTOOL_TEST_TRIANGLESORTER_HPP

#include <cppunit/extensions/HelperMacros.h>
#include <set>
#include "TestResources/Meshes.hpp"
#include "TriangleSorter.h"


class TriangleSorterTests : public CppUnit::TestFixture {
  CPPUNIT_TEST_SUITE(TriangleSorterTests);
  CPPUNIT_TEST(Trivial);
  CPPUNIT_TEST(TrivialDeep);
  CPPUNIT_TEST(Sphere);
  CPPUNIT_TEST_SUITE_END();
public:
  
  void Trivial() {
    // 16 cube
    TriTree::Int levels = 4;
    // put triangles onto the 8 cube level
    TriTree::Int tri_level = 3;
    auto triv = SimpleMeshFactory::MkTrivial();
    // Sort!
    auto tree = TrianglesToTree(levels, tri_level, triv->points, triv->triangles);

    CPPUNIT_ASSERT(tree.Level() == levels);
    TriTree::Int i = tri_level;
    tree.IterDepthFirst([&i](TriTree::Node& node) mutable {
	CPPUNIT_ASSERT(node.X() == 0);
	CPPUNIT_ASSERT(node.Y() == 0);
	CPPUNIT_ASSERT(node.Z() == 0);
	CPPUNIT_ASSERT(node.Level() == i);
	++i;
      });
    
    auto node = tree.Get(0,0,0, tri_level);
    auto ids = node->Data();
    std::sort(ids.begin(), ids.end());
    CPPUNIT_ASSERT(ids.size() == 2);
    CPPUNIT_ASSERT(ids[0] == 0);
    CPPUNIT_ASSERT(ids[1] == 1);

  }

  void TrivialDeep() {
    // 1024 cube
    TriTree::Int levels = 10;
    // put triangles onto the 8 cube level
    TriTree::Int tri_level = 3;
    auto triv = SimpleMeshFactory::MkTrivial();
    // Sort!
    auto tree = TrianglesToTree(levels, tri_level, triv->points, triv->triangles);
    
    CPPUNIT_ASSERT(tree.Level() == levels);
    TriTree::Int i = tri_level;
    tree.IterDepthFirst([&i](TriTree::Node& node) mutable {
	CPPUNIT_ASSERT(node.X() == 0);
	CPPUNIT_ASSERT(node.Y() == 0);
	CPPUNIT_ASSERT(node.Z() == 0);
	CPPUNIT_ASSERT(node.Level() == i);
	++i;
      });
    
    auto node = tree.Get(0,0,0, tri_level);
    auto ids = node->Data();
    std::sort(ids.begin(), ids.end());
    CPPUNIT_ASSERT(ids.size() == 2);
    CPPUNIT_ASSERT(ids[0] == 0);
    CPPUNIT_ASSERT(ids[1] == 1);
  }

  template <class T>
  bool overlap1d(T aMin, T aMax, T bMin, T bMax) {
    return (aMin < bMax) && (aMax > bMin);
  }

  bool overlap3d(const Vector& aMin, const Vector& aMax,
		 const Vector& bMin, const Vector& bMax) {
    bool ans = true;
    for (auto d =0; d<3; ++d)
      ans &= overlap1d(aMin[d], aMax[d], bMin[d], bMax[d]);
    
    return ans;
  }
  
  void Sphere() {
    TriTree::Int levels = 5;
    TriTree::Int tri_level = 3;
    
    auto sphere = SimpleMeshFactory::MkSphere();
    auto tree = TrianglesToTree(levels, tri_level, sphere->points, sphere->triangles);

    std::set<int> seen_tris;

    tree.IterDepthFirst(tri_level, tri_level, [&](TriTree::Node& node) {
	const auto& triIds = node.Data();
	const Vector offset(node.X(), node.Y(), node.Z());
	auto node_min = offset - 1;
	auto node_max = offset + (1 << tri_level) + 1;
	
	for (auto iTri: triIds) {
	  Vector tri_min(std::numeric_limits<double>::infinity());
	  Vector tri_max(0);
	  for (auto iPt: sphere->triangles[iTri]) {
	    tri_min.UpdatePointwiseMin(sphere->points[iPt]);
	    tri_max.UpdatePointwiseMax(sphere->points[iPt]);
	  }

	  CPPUNIT_ASSERT(overlap3d(tri_min, tri_max, node_min, node_max));
	  seen_tris.insert(triIds.begin(), triIds.end());
	}
      });

    // Assert that we have seen all the triangles in the input
    CPPUNIT_ASSERT(seen_tris.size() == sphere->triangles.size());
    int i = 0;
    for (auto seen_id : seen_tris) {
      CPPUNIT_ASSERT(seen_id == i);
      ++i;
    }
  }
    
};
CPPUNIT_TEST_SUITE_REGISTRATION(TriangleSorterTests);

#endif
