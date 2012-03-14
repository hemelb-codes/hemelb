#ifndef HEMELB_UNITTESTS_LBTESTS_BOUNDARIES_BOUNDARYTESTS_H
#define HEMELB_UNITTESTS_LBTESTS_BOUNDARIES_BOUNDARYTESTS_H
#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>

#include "unittests/helpers/FourCubeBasedTestFixture.h"
#include "lb/boundaries/BoundaryValues.h"

namespace hemelb
{
  namespace unittests
  {
    namespace lbtests
    {
      namespace boundaries
      {
        using namespace hemelb::lb::boundaries;
        /**
         * Class to test the collision operators. These tests are for the functions involved in
         * calculating the post-collision values, specifically CalculatePreCollision and Collide.
         * For each collision operator, we test that these functions return the expected values of
         * hydrodynamic quantities and f-distribution values.
         *
         * Note that we are only testing collision operators here, so we
         * can assume that the kernel objects work perfectly.
         */
        class BoundaryTests : public helpers::FourCubeBasedTestFixture
        {
          CPPUNIT_TEST_SUITE(BoundaryTests);
              CPPUNIT_TEST(TestConstruct);
            CPPUNIT_TEST_SUITE_END();
          public:
            void setUp()
            {
              FourCubeBasedTestFixture::setUp();
            }
            void tearDown()
            {
              FourCubeBasedTestFixture::tearDown();
            }
          private:
            void TestConstruct()
            {
              BoundaryValues inlets(hemelb::geometry::INLET_TYPE, latDat, simConfig->Inlets, simState, unitConverter);
            }
        };
        //BoundaryTests
        CPPUNIT_TEST_SUITE_REGISTRATION(BoundaryTests);
      } // boundaries
    } //lbtests
  } //unittests
} //hemelb
#endif // HEMELB_UNITTESTS_LBTESTS_BOUNDARIES_BOUNDARYTESTS_H
