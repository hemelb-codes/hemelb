#ifndef HEMELBSETUPTOOL_TEST_BUILDCGALPOLYGON_HPP
#define HEMELBSETUPTOOL_TEST_BUILDCGALPOLYGON_HPP

#include <cppunit/extensions/HelperMacros.h>
#include "TestResources/Meshes.hpp"
#include "BuildCGALPolygon.h"

//#include <iostream>

class BuildCGALPolygonTests : public CppUnit::TestFixture {
	  CPPUNIT_TEST_SUITE(BuildCGALPolygonTests);
	  CPPUNIT_TEST(Sphere);
	  CPPUNIT_TEST_SUITE_END();
public:
	  void Sphere() {
		  auto sphere = SimpleMeshFactory::MkSphere();
		  auto surface = MkCGALMesh(sphere->points, sphere->triangles, sphere->labels);
		  CPPUNIT_ASSERT(surface->is_closed());
		  CPPUNIT_ASSERT(surface->is_pure_triangle());
		  CPPUNIT_ASSERT(surface->is_valid());
	  }
};
CPPUNIT_TEST_SUITE_REGISTRATION(BuildCGALPolygonTests);
#endif
