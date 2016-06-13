#ifndef HEMELBSETUPTOOL_TEST_TRIANGLESORTER_HPP
#define HEMELBSETUPTOOL_TEST_TRIANGLESORTER_HPP

#include <cppunit/extensions/HelperMacros.h>
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

  void Sphere() {
    // levels = 5
    // tri_level = 3
    // points, triangles, normals = GetSphereNumpy()
    // tree = TriangleSorter.TrianglesToTree(levels, tri_level, points, triangles)
    
    // seen_tris = set()
    // for node in tree.IterDepthFirst(tri_level, tri_level):
    //     n = len(node.triIds)
    //     node_points = np.zeros((n, 3, 3))
    //     node_point_ids = triangles[node.triIds]
    //     for iPt in xrange(3):
    //         node_points[:,iPt,:] = points[node_point_ids[:,iPt]]
        
    //     tri_min = node_points.min(axis=1)
    //     tri_max = node_points.max(axis=1)
        
    //     nd_min = (node.offset - 1)[np.newaxis,:]
    //     nd_max = (node.offset + 2**tri_level + 1)[np.newaxis,:]
        
    //     assert np.all(overlap3d(tri_min, tri_max, nd_min, nd_max))
            
    //     seen_tris.update(node.triIds)
    // allids = sorted(seen_tris)
    // assert np.all(np.equal(allids, np.arange(len(triangles))))
  }
    
};
CPPUNIT_TEST_SUITE_REGISTRATION(TriangleSorterTests);

#endif
