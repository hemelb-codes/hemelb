#ifndef HEMELBSETUPTOOL_TEST_FLOODFILL_HPP
#define HEMELBSETUPTOOL_TEST_FLOODFILL_HPP

#include <cppunit/extensions/HelperMacros.h>
#include "TestResources/Meshes.hpp"
#include "FloodFill.h"
#include "FloodFillImpl.h"

class FloodFillTests : public CppUnit::TestFixture {
	  CPPUNIT_TEST_SUITE(FloodFillTests);

	  CPPUNIT_TEST(GetStartSingle);
	  CPPUNIT_TEST(GetStartMulti);

	  CPPUNIT_TEST(BuildMaskSingle);
	  CPPUNIT_TEST_SUITE_END();
public:

	  void GetStartSingle() {
		  // Make a FluidTree
		  FluidTree tree(4);
		  auto node = tree.GetCreate(1,2,3,0);
		  auto start_node = GetStart(tree);

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

		  auto start_node = GetStart(tree);

		  for (auto i: range(3))
			  CPPUNIT_ASSERT(start_node[i] >= 3 && start_node[i] <= 5);
		  CPPUNIT_ASSERT(start_node[3] == 0);
	  }

	  void BuildMaskSingle() {
		  // Make a FluidTree
		  FluidTree tree(4);
		  auto node = tree.GetCreate(1,2,3,0);

		  MaskBuilder b(tree);
		  auto output = b();

		  // One node at each level
		  CPPUNIT_ASSERT(output.indices[3].size() == 8);
		  CPPUNIT_ASSERT(output.indices[2].size() == 8);
		  CPPUNIT_ASSERT(output.indices[1].size() == 8);
		  // bottom level always empty
		  CPPUNIT_ASSERT(output.indices[0].size() == 0);
		  // because we use the mask instead
		  CPPUNIT_ASSERT(output.mask.size() == 8);

		  // top level only has octant 0,0,0
		  CPPUNIT_ASSERT_EQUAL(0UL, output.indices[3][0]);
		  for (auto i: range(1, 8))
			  CPPUNIT_ASSERT_EQUAL(~0UL, output.indices[3][i]);

		  // same for second level
		  CPPUNIT_ASSERT_EQUAL(0UL, output.indices[2][0]);
		  for (auto i: range(1, 8))
			  CPPUNIT_ASSERT_EQUAL(~0UL, output.indices[2][i]);

		  // Now 0, 1, 1 == 0 + 2 + 1
		  for (auto i: range(8))
			  CPPUNIT_ASSERT_EQUAL(i==3 ? 0 : ~0UL, output.indices[1][i]);

		  // Bottom level mask is true at 1,0,1 == 5
		  for (auto i: range(8))
			  CPPUNIT_ASSERT((i==5) == output.mask[i]);

		  CPPUNIT_ASSERT(output.Get(1,2,3));
		  CPPUNIT_ASSERT(output.Get(1,2,4) == false);

	  }
};
CPPUNIT_TEST_SUITE_REGISTRATION(FloodFillTests);
#endif
