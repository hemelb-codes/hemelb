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
              CPPUNIT_TEST(TestUpdate);
            CPPUNIT_TEST_SUITE_END();
          public:
            void setUp()
            {
              FourCubeBasedTestFixture::setUp();

              inlets = new BoundaryValues(hemelb::geometry::INLET_TYPE,
                                          latDat,
                                          simConfig->Inlets,
                                          simState,
                                          unitConverter);
            }
            void tearDown()
            {
              delete inlets;
              FourCubeBasedTestFixture::tearDown();
            }
          private:
            void TestConstruct()
            {
              double targetStartDensity = pressureToDensity(80.0-1.0);
              CPPUNIT_ASSERT_EQUAL(targetStartDensity, inlets->GetBoundaryDensity(0));
            }
            void TestUpdate()
            {
              inlets->RequestComms();
              double targetStartDensity = pressureToDensity(80.0 - 1.0);
              CPPUNIT_ASSERT_EQUAL(targetStartDensity, inlets->GetBoundaryDensity(0));
              for (unsigned int step = 0; step < simState->GetTotalTimeSteps() / 2; step++)
              {
                simState->Increment();
              }
              inlets->RequestComms();
              double targetMidDensity = pressureToDensity(80.0 + 1.0);
              CPPUNIT_ASSERT_EQUAL(targetMidDensity, inlets->GetBoundaryDensity(0));
              for (unsigned int step = 0; step < simState->GetTotalTimeSteps() / 2; step++)
              {
                simState->Increment();
              }
              inlets->RequestComms();
              double targetEndDensity = pressureToDensity(80.0 - 1.0);
              CPPUNIT_ASSERT_EQUAL(targetEndDensity, inlets->GetBoundaryDensity(0));
            }
            double pressureToDensity(double pressure)
            {
              double inverse_velocity = simState->GetTimeStepLength() / latDat->GetVoxelSize();
              return 1
                  + ( pressure - REFERENCE_PRESSURE_mmHg) * mmHg_TO_PASCAL * inverse_velocity * inverse_velocity
                      / (Cs2 * BLOOD_DENSITY_Kg_per_m3);
            }
            BoundaryValues *inlets;
        };
        //BoundaryTests
        CPPUNIT_TEST_SUITE_REGISTRATION(BoundaryTests);
      } // boundaries
    } //lbtests
  } //unittests
} //hemelb
#endif // HEMELB_UNITTESTS_LBTESTS_BOUNDARIES_BOUNDARYTESTS_H
