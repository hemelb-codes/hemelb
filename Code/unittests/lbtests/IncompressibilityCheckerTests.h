#ifndef HEMELB_UNITTESTS_LBTESTS_INCOMPRESSIBILITYCHECKERTESTS_H
#define HEMELB_UNITTESTS_LBTESTS_INCOMPRESSIBILITYCHECKERTESTS_H

#include <cppunit/TestFixture.h>

#include "lb/IncompressibilityChecker.hpp"
#include "unittests/FourCubeLatticeData.h"
#include "unittests/lbtests/BroadcastMock.h"
#include "unittests/reporting/Mocks.h"

#include "unittests/helpers/FourCubeBasedTestFixture.h"

namespace hemelb
{
  namespace unittests
  {
    namespace lbtests
    {
      class IncompressibilityCheckerTests : public helpers::FourCubeBasedTestFixture
      {
          CPPUNIT_TEST_SUITE (IncompressibilityCheckerTests);
          CPPUNIT_TEST (TestIncompressibilityChecker);CPPUNIT_TEST_SUITE_END();

        public:
          void setUp()
          {
            FourCubeBasedTestFixture::setUp();
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
            LbTestsHelper::InitialiseAnisotropicTestData<lb::lattices::D3Q15>(latDat);
            lb::MacroscopicPropertyCache* cache = new lb::MacroscopicPropertyCache(*simState, *latDat);

            cache->densityCache.SetRefreshFlag();
            lbtests::LbTestsHelper::UpdatePropertyCache<lb::lattices::D3Q15>(*latDat, *cache, *simState);

            hemelb::reporting::Timers timings;
            net::Net net;
            hemelb::lb::IncompressibilityChecker<net::BroadcastMock> incompChecker(latDat,
                                                                                   &net,
                                                                                   simState,
                                                                                   *cache,
                                                                                   timings,
                                                                                   10.0); // Will accept a max/min of (21.45, 12) but not (100,1)

            // These are the smallest and largest density values in FourCubeLatticeData by default
            //! @23 The lattice class below must be consistent with the one used in FourCubeLatticeData. Consider templating FourCubeLatticeData over lattice class, so both can be controlled from the test.
            distribn_t numDirections = (distribn_t) lb::lattices::D3Q15::NUMVECTORS;
            distribn_t numSites = (distribn_t) latDat->GetLocalFluidSiteCount();
            distribn_t smallestDefaultDensity = numDirections * (numDirections + 1) / 20; // = sum_{j=1}^{numDirections} j/10 = 12 with current configuration of FourCubeLatticeData
            distribn_t largestDefaultDensity = (numDirections * (numDirections + 1) / 20)
                + ( (numSites - 1) * numDirections / 100); // = sum_{j=1}^{numDirections} j/10 + (numSites-1)/100 = 21.45 with current configuration of FourCubeLatticeData

            // Not available until the first broadcast has finished
            CPPUNIT_ASSERT(!incompChecker.AreDensitiesAvailable());

            cache->densityCache.SetRefreshFlag();
            LbTestsHelper::UpdatePropertyCache<lb::lattices::D3Q15>(*latDat, *cache, *simState);
            AdvanceActorOneTimeStep(incompChecker);

            CPPUNIT_ASSERT_DOUBLES_EQUAL(smallestDefaultDensity, incompChecker.GetGlobalSmallestDensity(), eps);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(largestDefaultDensity, incompChecker.GetGlobalLargestDensity(), eps);
            CPPUNIT_ASSERT_DOUBLES_EQUAL( (largestDefaultDensity - smallestDefaultDensity)
                                             / hemelb::lb::REFERENCE_DENSITY,
                                         incompChecker.GetMaxRelativeDensityDifference(),
                                         eps);
            CPPUNIT_ASSERT(incompChecker.IsDensityDiffWithinRange());

            // The broadcast mock injects some smaller and larger densities (1,100) coming from one of the children
            CPPUNIT_ASSERT(incompChecker.AreDensitiesAvailable());

            cache->densityCache.SetRefreshFlag();
            LbTestsHelper::UpdatePropertyCache<lb::lattices::D3Q15>(*latDat, *cache, *simState);
            AdvanceActorOneTimeStep(incompChecker);

            CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, incompChecker.GetGlobalSmallestDensity(), eps);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(100.0, incompChecker.GetGlobalLargestDensity(), eps);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(99.0, incompChecker.GetMaxRelativeDensityDifference(), eps);
            CPPUNIT_ASSERT(!incompChecker.IsDensityDiffWithinRange());

            // The previous values are not reported by any children anymore. Testing that the checker remembers them
            CPPUNIT_ASSERT(incompChecker.AreDensitiesAvailable());

            LbTestsHelper::UpdatePropertyCache<lb::lattices::D3Q15>(*latDat, *cache, *simState);
            AdvanceActorOneTimeStep(incompChecker);

            CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, incompChecker.GetGlobalSmallestDensity(), eps);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(100.0, incompChecker.GetGlobalLargestDensity(), eps);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(99.0, incompChecker.GetMaxRelativeDensityDifference(), eps);
            CPPUNIT_ASSERT(!incompChecker.IsDensityDiffWithinRange());
          }

        private:
          distribn_t eps;
      };

      CPPUNIT_TEST_SUITE_REGISTRATION (IncompressibilityCheckerTests);
    }
  }
}

#endif /* HEMELB_UNITTESTS_LBTESTS_IMCOMPRESSIBILITYCHECKERTESTS_H */
