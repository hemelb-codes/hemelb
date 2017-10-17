// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>

#include "test_Oct.hpp"
#include "test_MkCgalMesh.hpp"
#include "test_TriangleSorter.hpp"
#include "test_SegmentFactory.hpp"
#include "test_SurfaceVoxeliser.hpp"
#include "test_FloodFill.hpp"
#include "test_SectionTree.hpp"

#include "test.hpp"

bool run_all_tests() {
  CppUnit::TextUi::TestRunner runner;
  CppUnit::TestFactoryRegistry &registry = CppUnit::TestFactoryRegistry::getRegistry();
  runner.addTest( registry.makeTest() );
  return runner.run("", false);
}
