#ifndef HEMELBSETUPTOOL_TEST_SURFACEVOXELISER_HPP
#define HEMELBSETUPTOOL_TEST_SURFACEVOXELISER_HPP

#include <cppunit/extensions/HelperMacros.h>
#include "TestResources/Meshes.hpp"
#include "TriangleSorter.h"
#include "SurfaceVoxeliser.h"
#include "Neighbours.h"

#include "range.hpp"
#include "enumerate.hpp"

class SurfaceVoxeliserTests : public CppUnit::TestFixture {
  CPPUNIT_TEST_SUITE(SurfaceVoxeliserTests);
  
  CPPUNIT_TEST(NeighbourInverses);
  CPPUNIT_TEST(Intersections);
  CPPUNIT_TEST(Trivial);
  
  CPPUNIT_TEST_SUITE_END();
  
public:
  void NeighbourInverses() {
    auto& vectors = Neighbours::GetDisplacements();
    auto& opps = Neighbours::GetOpposites();
    CPPUNIT_ASSERT(vectors.size() == 26);
    CPPUNIT_ASSERT(opps.size() == 26);
    for (auto i: range(26)) {
      CPPUNIT_ASSERT(vectors[i] + vectors[opps[i]] == Index::Zero());
    }
  }
  
  void Intersections() {
    auto triv = SimpleMeshFactory::MkTrivial();
    SurfaceVoxeliser voxer(4, triv->points, triv->triangles, triv->normals, triv->labels);

//    double t;
//    // Well away from the square so def no intersection
//    t = voxer.IntersectLinkWithTriangle({2,3,4}, {1,0,0}, 0);
//    CPPUNIT_ASSERT(t > 1);
//
//    t = voxer.IntersectLinkWithTriangle({1,2,2}, {1,0,0}, 0);
//    CPPUNIT_ASSERT(t > 1);
//    t = voxer.IntersectLinkWithTriangle({1,2,2}, {1,0,0}, 1);
//    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.2, t, 1e-6);
//    t = vc.IntersectLinkWithTriangle({2,2,2}, {-1,0,0}, 1);
//    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.8, t, 1e-6);
  }

  void Trivial() {
	    // 16 cube
	    auto levels = 4;
	    auto tri_level = 2;
	    auto triv = SimpleMeshFactory::MkTrivial();
	    auto tree = TrianglesToTreeSerial(levels, tri_level, triv->points, triv->triangles);

	    SurfaceVoxeliser voxer(4, triv->points, triv->triangles, triv->normals, triv->labels);
	    VoxTree edge_tree(levels);
	    std::list<VoxTree::NodePtr> edges;
	    auto append = std::back_inserter(edges);
	    tree.IterDepthFirst(tri_level, tri_level,
	    		[&voxer, &append](TriTree::Node& node) mutable {
	    	append = voxer.ComputeIntersectionsForRegion(node);
	    });
//
//    vc = SurfaceVoxeliser(points, triangles, normals, labels)
//    new_tree = vc(tree, tri_level)
//
//    for voxel in new_tree.IterDepthFirst(0,0):
//        # Only nodes with x = 2 should exist
//        assert voxel.offset[0] == 2
//        # The cut distances should be ~0.8
//        cut_mask = voxel.cut_tri >= 0
//        real_cd = voxel.cut_dist[cut_mask]
//        assert np.all((real_cd - 0.8)**2 < 1e-6)
  }

};

CPPUNIT_TEST_SUITE_REGISTRATION(SurfaceVoxeliserTests);

#endif
