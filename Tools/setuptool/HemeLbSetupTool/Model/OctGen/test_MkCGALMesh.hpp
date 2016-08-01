#ifndef HEMELBSETUPTOOL_TEST_MKCGALMESH_HPP
#define HEMELBSETUPTOOL_TEST_MKCGALMESH_HPP

#include <cppunit/extensions/HelperMacros.h>
#include "TestResources/Meshes.hpp"
#include "MkCGALMesh.h"

//#include <iostream>

class MkCGALMeshTests : public CppUnit::TestFixture {
	  CPPUNIT_TEST_SUITE(MkCGALMeshTests);
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
CPPUNIT_TEST_SUITE_REGISTRATION(MkCGALMeshTests);
#endif
