// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.
#ifndef HEMELBSETUPTOOL_TEST_MKCGALMESHTESTS_HPP
#define HEMELBSETUPTOOL_TEST_MKCGALMESHTESTS_HPP

#include <cppunit/extensions/HelperMacros.h>
#include "TestResources/Meshes.hpp"
#include "MkCgalMesh.h"

//#include <iostream>

class MkCgalMeshTests : public CppUnit::TestFixture {
	  CPPUNIT_TEST_SUITE(MkCgalMeshTests);
	  CPPUNIT_TEST(Sphere);
	  CPPUNIT_TEST_SUITE_END();
public:
	  void Sphere() {
		  auto sphere = SimpleMeshFactory::MkSphere();
		  auto surface = MkCgalMesh(sphere->points, sphere->triangles);
		  CPPUNIT_ASSERT(surface->is_closed());
		  CPPUNIT_ASSERT(surface->is_pure_triangle());
		  CPPUNIT_ASSERT(surface->is_valid());
	  }
};
CPPUNIT_TEST_SUITE_REGISTRATION(MkCgalMeshTests);
#endif
