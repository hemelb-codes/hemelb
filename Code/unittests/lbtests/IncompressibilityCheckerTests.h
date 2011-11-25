#ifndef HEMELB_UNITTESTS_LBTESTS_IMCOMPRESSIBILITYCHECKERTESTS_H
#define HEMELB_UNITTESTS_LBTESTS_IMCOMPRESSIBILITYCHECKERTESTS_H

#include <cppunit/TestFixture.h>

#include "lb/IncompressibilityChecker.h"
#include "unittests/FourCubeLatticeData.h"
#include "net/BroadcastMock.h"

namespace hemelb
{
  namespace unittests
  {
    namespace lbtests
    {
      class IncompressibilityCheckerTests : public CppUnit::TestFixture
      {
          CPPUNIT_TEST_SUITE( IncompressibilityCheckerTests);
          CPPUNIT_TEST( TestIncompressibilityChecker);
          CPPUNIT_TEST_SUITE_END();

        public:
          void setUp()
          {
            eps = 1e-5;
          }

          void tearDown()
          {
          }

          void AdvanceActorOneTimeStep(net::IteratedAction& actor)
          {
            actor.RequestComms();
            actor.PreSend();
            actor.PreReceive();
            actor.PostReceive();
            actor.EndIteration();
          }

          void TestIncompressibilityChecker()
          {
            // The following objects will be ignored since we are using a broadcast mock, no need
            // to initialise them properly
            net::Net net;
            lb::SimulationState simulationState(1000, 1);

            FourCubeLatticeData latticeData;
            LbTestsHelper::InitialiseAnisotropicTestData(&latticeData);
            latticeData.SwapOldAndNew();
            lb::IncompressibilityChecker<net::BroadcastMock> incompChecker(&latticeData,
                                                                           &net,
                                                                           &simulationState);

            // These are the smallest and largest density values in FourCubeLatticeData by default
            AdvanceActorOneTimeStep(incompChecker);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(10, incompChecker.GetGlobalSmallestDensity(), eps);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(21.45, incompChecker.GetGlobalLargestDensity(), eps);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(11.45, incompChecker.GetMaxDensityDifference(), eps);

            // The broadcast mock injects some smaller and larger densities coming from one of the children
            AdvanceActorOneTimeStep(incompChecker);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, incompChecker.GetGlobalSmallestDensity(), eps);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(100.0, incompChecker.GetGlobalLargestDensity(), eps);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(99.0, incompChecker.GetMaxDensityDifference(), eps);

            // The previous values are not reported by any children anymore. Testing that the checker remembers them
            AdvanceActorOneTimeStep(incompChecker);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, incompChecker.GetGlobalSmallestDensity(), eps);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(100.0, incompChecker.GetGlobalLargestDensity(), eps);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(99.0, incompChecker.GetMaxDensityDifference(), eps);
          }

        private:
          distribn_t eps;
      };

      CPPUNIT_TEST_SUITE_REGISTRATION( IncompressibilityCheckerTests);
    }
  }
}

#endif /* HEMELB_UNITTESTS_LBTESTS_IMCOMPRESSIBILITYCHECKERTESTS_H */
