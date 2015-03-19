//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELB_UNITTESTS_REDBLOOD_FLOWEXTENSION_TESTS_H
#define HEMELB_UNITTESTS_REDBLOOD_FLOWEXTENSION_TESTS_H

#include <cppunit/TestFixture.h>
#include "redblood/FlowExtension.h"
#include "unittests/redblood/Fixtures.h"
#include "util/Vector3D.h"
#include "units.h"

namespace hemelb
{
  namespace unittests
  {
    namespace redblood
    {
      class FlowExtensionTests : public FlowExtensionFixture {
        CPPUNIT_TEST_SUITE (FlowExtensionTests);
        CPPUNIT_TEST (testAxis);
        CPPUNIT_TEST (testCircumference);
        CPPUNIT_TEST (testLinearWeight);
        CPPUNIT_TEST_SUITE_END();

        typedef util::Vector3D<LatticeDistance> Point;

        public:

          void testAxis() {
            // Check points along the axis from just outside the start of the cylinder (x = 0) to just inside
            CPPUNIT_ASSERT_MESSAGE("Point (-0.001, 0, 0) is outside the cylinder", !contains(flowExt, Point(-0.001, 0, 0)));
            CPPUNIT_ASSERT_MESSAGE("Point ( 0.000, 0, 0) is inside the cylinder",   contains(flowExt, Point( 0.000, 0, 0)));
            CPPUNIT_ASSERT_MESSAGE("Point ( 0.001, 0, 0) is inside the cylinder",   contains(flowExt, Point( 0.001, 0, 0)));

            // Check points along the axis from just inside the end of the cylinder (x = 10) to just outside
            CPPUNIT_ASSERT_MESSAGE("Point ( 9.999, 0, 0) is inside the cylinder",   contains(flowExt, Point( 9.999, 0, 0)));
            CPPUNIT_ASSERT_MESSAGE("Point (10.000, 0, 0) is inside the cylinder",   contains(flowExt, Point(10.000, 0, 0)));
            CPPUNIT_ASSERT_MESSAGE("Point (10.001, 0, 0) is outside the cylinder", !contains(flowExt, Point(10.001, 0, 0)));
          }

          void testCircumference() {
            // Points around the circumference should be inside the cylinder
            CPPUNIT_ASSERT_MESSAGE("Point (0.0,    1.000,   0.000 ) is inside the cylinder", contains(flowExt, Point(0.0, 1.0, 0.0)));
            CPPUNIT_ASSERT_MESSAGE("Point (0.0,  cos(45),  sin(45)) is inside the cylinder", contains(flowExt, Point(0.0, cos(45), sin(45))));
            CPPUNIT_ASSERT_MESSAGE("Point (0.0,    0.000,   1.000 ) is inside the cylinder", contains(flowExt, Point(0.0, 0.0, 1.0)));
            CPPUNIT_ASSERT_MESSAGE("Point (0.0,  cos(45), -sin(45)) is inside the cylinder", contains(flowExt, Point(0.0, cos(45), -sin(45))));
            CPPUNIT_ASSERT_MESSAGE("Point (0.0,   -1.000,   0.000 ) is inside the cylinder", contains(flowExt, Point(0.0, -1.0, 0.0)));
            CPPUNIT_ASSERT_MESSAGE("Point (0.0, -cos(45), -sin(45)) is inside the cylinder", contains(flowExt, Point(0.0, -cos(45), -sin(45))));
            CPPUNIT_ASSERT_MESSAGE("Point (0.0,    0.000,  -1.000 ) is inside the cylinder", contains(flowExt, Point(0.0, 0.0, -1.0)));
            CPPUNIT_ASSERT_MESSAGE("Point (0.0, -cos(45),  sin(45)) is inside the cylinder", contains(flowExt, Point(0.0, -cos(45), sin(45))));

            // Points forming a square around the axis should be outside the cylinder
            CPPUNIT_ASSERT_MESSAGE("Point(0.0,  1.0,  1.0) is outside the cylinder", !contains(flowExt, Point(0.0,  1.0, 1.0)));
            CPPUNIT_ASSERT_MESSAGE("Point(0.0,  1.0, -1.0) is outside the cylinder", !contains(flowExt, Point(0.0,  1.0, -1.0)));
            CPPUNIT_ASSERT_MESSAGE("Point(0.0, -1.0, -1.0) is outside the cylinder", !contains(flowExt, Point(0.0, -1.0, -1.0)));
            CPPUNIT_ASSERT_MESSAGE("Point(0.0, -1.0,  1.0) is outside the cylinder", !contains(flowExt, Point(0.0, -1.0, 1.0)));
          }

          void testLinearWeight()
          {
            FlowExtension const flow = {{-1.0, 0, 0}, {0.5, 0.5, 0.5}, 2.0, 0.5, 1.5};
            for(auto y: {0.5, 0.7})
            {
              CPPUNIT_ASSERT_DOUBLES_EQUAL(0e0, linearWeight(flow, {2.4, y, 0.5}), 1e-8);
              CPPUNIT_ASSERT_DOUBLES_EQUAL(1e0, linearWeight(flow, {0.5, y, 0.5}), 1e-8);
              for(auto epsilon: {0.3, 0.5, 0.7, 1.0})
              {
                auto const pos = flow.origin + flow.normal * flow.fadeLength * epsilon
                  + LatticePosition{0, y - 0.5, 0};
                CPPUNIT_ASSERT_DOUBLES_EQUAL(1e0 - epsilon, linearWeight(flow, pos), 1e-8);
              }
            }
            CPPUNIT_ASSERT_DOUBLES_EQUAL(0e0, linearWeight(flow, {0.5, 0.5, 1.0 + 1e-8}), 1e-8);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(1e0, linearWeight(flow, {0.5, 0.5, 1.0 - 1e-8}), 1e-8);
          }
      };

      CPPUNIT_TEST_SUITE_REGISTRATION (FlowExtensionTests);

    } // namespace: redblood
  } // namespace: unittests
} // namespace: hemelb

#endif // HEMELB_UNITTESTS_REDBLOOD_FLOWEXTENSION_TESTS_H
