// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.
#ifndef HEMELBSETUPTOOL_TEST_FLOODFILL_HPP
#define HEMELBSETUPTOOL_TEST_FLOODFILL_HPP

#include <cppunit/extensions/HelperMacros.h>
#include "TestResources/Meshes.hpp"
#include "FloodFill.h"

class FloodFillTests : public CppUnit::TestFixture {
  CPPUNIT_TEST_SUITE(FloodFillTests);
  
  CPPUNIT_TEST(GetStartSingle);
  CPPUNIT_TEST(GetStartMulti);
  CPPUNIT_TEST(FillSphere);
  CPPUNIT_TEST_SUITE_END();

  using Idx = FloodFill::Idx;
public:

  void GetStartSingle() {
    // Make a FluidTree
    FluidTree tree(4);
    auto node = tree.GetCreate(1,2,3,0);
    FloodFill ff(tree);
    auto start_node = ff.GetStart();

    CPPUNIT_ASSERT(Idx({1,2,3,0}) == start_node);
  }

  void GetStartMulti() {
    // Make a FluidTree
    FluidTree tree(4);
    // Cube of sites
    for (auto i: {3,4,5})
      for (auto j: {3,4,5})
	for (auto k: {3,4,5})
	  auto node = tree.GetCreate(i,j,k,0);

    FloodFill ff(tree);
    auto start_node = ff.GetStart();

    for (auto i: range(3))
      CPPUNIT_ASSERT(start_node[i] >= 3 && start_node[i] <= 5);
    CPPUNIT_ASSERT(start_node[3] == 0);
  }

  void FillSphere() {
    TriTree::Int levels = 5;
    TriTree::Int tri_level = 3;
    
    Vector sphere_centre(15.5);
    double sphere_radius = 10.0;
    auto r2 = sphere_radius*sphere_radius;
    
    auto sphere = SimpleMeshFactory::MkSphere();
    auto tree = TrianglesToTreeSerial(levels, tri_level, sphere->points,
				      sphere->triangles);
    
    SurfaceVoxeliser voxer(1 << tri_level, sphere->points,
			   sphere->triangles, sphere->normals, sphere->labels,
			   sphere->iolets);
    auto fluid_tree = voxer(tree, tri_level);

    // Fill the thing
    FloodFill ff(fluid_tree);
    auto mask = ff();
    
    // Assert things...
    // (24, 17, 20) is known to be outside so shouldn't exist
    CPPUNIT_ASSERT(!mask.Get(24, 17, 20, 0));
    
    mask.IterDepthFirst([&](MaskTree::ConstNodePtr node) {
	if (node->Level() == 0) {
	  // Get the node's coordinate
	  Vector posn(node->X(), node->Y(), node->Z());
	  auto dr = posn - sphere_centre;
	  CPPUNIT_ASSERT(dr.GetMagnitudeSquared() < r2);
	}
      });
      
    
  }
};
CPPUNIT_TEST_SUITE_REGISTRATION(FloodFillTests);
#endif
