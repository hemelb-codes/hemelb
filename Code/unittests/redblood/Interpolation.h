//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELB_UNITTESTS_REDBLOOD_INTERPOLATION_H
#define HEMELB_UNITTESTS_REDBLOOD_INTERPOLATION_H

#include <cppunit/TestFixture.h>
#include "redblood/interpolation.h"
#include "unittests/redblood/Fixtures.h"

namespace hemelb { namespace unittests {

class InterpolationTests : public CppUnit::TestFixture, public Comparisons {
    CPPUNIT_TEST_SUITE(InterpolationTests);
    CPPUNIT_TEST(testRegionIterator);
    CPPUNIT_TEST_SUITE_END();

public:

    void testRegionIterator() {
      using hemelb::redblood::RegionIterator;

      LatticeVector vectors[] = {
        LatticeVector(4, 3, 2),
        LatticeVector(4, 3, 3),
        LatticeVector(4, 3, 4),
        LatticeVector(4, 4, 2),
        LatticeVector(4, 5, 2),
        LatticeVector(5, 3, 2),
        LatticeVector(6, 3, 2),
        LatticeVector(6, 5, 4),
      };
      size_t incs[] = {
        0, 1, 1, 1, 3, 3, 9, 8,
        666 // break
      };

      // Checks iteration goes through correct sequence
      RegionIterator iterator(LatticeVector(5, 4, 3), 1);
      for(size_t i(0); incs[i] < 666; ++i)  {
        for(size_t j(0); j < incs[i]; ++j, ++iterator);
        CPPUNIT_ASSERT(is_zero(*iterator - vectors[i]));
        CPPUNIT_ASSERT(iterator.isValid());
      }
      // Checks iterator becomes invalid
      ++iterator;
      CPPUNIT_ASSERT(not iterator.isValid());
    }
};

CPPUNIT_TEST_SUITE_REGISTRATION(InterpolationTests);
}}

#endif // ONCE

