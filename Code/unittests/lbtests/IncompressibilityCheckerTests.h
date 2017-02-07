
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_UNITTESTS_LBTESTS_INCOMPRESSIBILITYCHECKERTESTS_H
#define HEMELB_UNITTESTS_LBTESTS_INCOMPRESSIBILITYCHECKERTESTS_H

#include <cppunit/TestFixture.h>

#include "lb/IncompressibilityChecker.h"
#include "net/phased/steps.h"
#include "timestep/TimeStepManager.h"

#include "unittests/FourCubeLatticeData.h"
#include "unittests/reporting/Mocks.h"

#include "unittests/helpers/FourCubeBasedTestFixture.h"
#include "debug/Debugger.h"

namespace hemelb
{
  namespace unittests
  {
    namespace lbtests
    {
      class IncompressibilityCheckerTests : public helpers::FourCubeBasedTestFixture
      {
          CPPUNIT_TEST_SUITE (IncompressibilityCheckerTests);
          CPPUNIT_TEST (TestIncompressibilityCheckerRootNode);
          CPPUNIT_TEST (TestIncompressibilityCheckerLeafNode);CPPUNIT_TEST_SUITE_END();

        public:
          void setUp()
          {
            FourCubeBasedTestFixture::setUp();
            LbTestsHelper::InitialiseAnisotropicTestData<lb::lattices::D3Q15>(latDat);
            cache = new lb::MacroscopicPropertyCache(*simState, *latDat);

            cache->densityCache.SetRefreshFlag();
            cache->velocityCache.SetRefreshFlag();
            lbtests::LbTestsHelper::UpdatePropertyCache<lb::lattices::D3Q15>(*latDat,
                                                                             *cache,
                                                                             *simState);

            // These are the smallest and largest density values in FourCubeLatticeData by default
            //! @23 The lattice class below must be consistent with the one used in FourCubeLatticeData. Consider templating FourCubeLatticeData over lattice class, so both can be controlled from the test.
            distribn_t numDirections = (distribn_t) lb::lattices::D3Q15::NUMVECTORS;
            distribn_t numSites = (distribn_t) latDat->GetLocalFluidSiteCount();
            smallestDefaultDensity = numDirections * (numDirections + 1) / 20; // = sum_{j=1}^{numDirections} j/10 = 12 with current configuration of FourCubeLatticeData
            largestDefaultDensity = (numDirections * (numDirections + 1) / 20)
                + ( (numSites - 1) * numDirections / 100); // = sum_{j=1}^{numDirections} j/10 + (numSites-1)/100 = 21.45 with current configuration of FourCubeLatticeData
            largestDefaultVelocityMagnitude = 0.0433012701892219; // with current configuration of FourCubeLatticeData

            eps = 1e-9;

            timings = new hemelb::reporting::Timers(Comms());
          }

          void tearDown()
          {
            delete cache;
            FourCubeBasedTestFixture::tearDown();
          }

  	  void AdvanceActorOneTimeStep(timestep::Actor& actor)
          {
	    timestep::TimeStepManager tsm(1);
	    tsm.AddToPhase(0, &actor);
	    tsm.DoStep();
	  }

          void TestIncompressibilityCheckerRootNode()
          {
            hemelb::lb::IncompressibilityChecker incompChecker(latDat,
                                                               Comms(),
                                                               simState,
                                                               *cache,
                                                               *timings,
                                                               10.0); // Will accept a max/min of (21.45, 12) but not (100,1)

            // Not available until the first broadcast has finished
            AdvanceActorOneTimeStep(incompChecker);

            CPPUNIT_ASSERT_DOUBLES_EQUAL(smallestDefaultDensity,
                                         incompChecker.GetGlobalSmallestDensity(),
                                         eps);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(largestDefaultDensity,
                                         incompChecker.GetGlobalLargestDensity(),
                                         eps);
            CPPUNIT_ASSERT_DOUBLES_EQUAL( (largestDefaultDensity - smallestDefaultDensity)
                                             / hemelb::lb::REFERENCE_DENSITY,
                                         incompChecker.GetMaxRelativeDensityDifference(),
                                         eps);
            CPPUNIT_ASSERT(incompChecker.IsDensityDiffWithinRange());
            CPPUNIT_ASSERT_DOUBLES_EQUAL(largestDefaultVelocityMagnitude,
                                         incompChecker.GetGlobalLargestVelocityMagnitude(),
                                         eps);

            // Insert values at 1 & 100
            debug::Debugger::Get()->BreakHere();
            cache->densityCache.Put(0, 1.0);
            cache->densityCache.Put(1, 100.0);
            cache->velocityCache.Put(2, LatticeVelocity(10.0, 0., 0.0));

            AdvanceActorOneTimeStep(incompChecker);

            CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, incompChecker.GetGlobalSmallestDensity(), eps);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(100.0, incompChecker.GetGlobalLargestDensity(), eps);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(99.0,
                                         incompChecker.GetMaxRelativeDensityDifference(),
                                         eps);
            CPPUNIT_ASSERT(!incompChecker.IsDensityDiffWithinRange());
            CPPUNIT_ASSERT_DOUBLES_EQUAL(10.0,
                                         incompChecker.GetGlobalLargestVelocityMagnitude(),
                                         eps);

            // Reset the previous values. Testing that the checker remembers the extrema
            cache->densityCache.SetRefreshFlag();
            lbtests::LbTestsHelper::UpdatePropertyCache<lb::lattices::D3Q15>(*latDat,
                                                                             *cache,
                                                                             *simState);
            AdvanceActorOneTimeStep(incompChecker);

            CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, incompChecker.GetGlobalSmallestDensity(), eps);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(100.0, incompChecker.GetGlobalLargestDensity(), eps);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(99.0,
                                         incompChecker.GetMaxRelativeDensityDifference(),
                                         eps);
            CPPUNIT_ASSERT(!incompChecker.IsDensityDiffWithinRange());
            CPPUNIT_ASSERT_DOUBLES_EQUAL(10.0,
                                         incompChecker.GetGlobalLargestVelocityMagnitude(),
                                         eps);
          }

          void TestIncompressibilityCheckerLeafNode()
          {
            hemelb::lb::IncompressibilityChecker incompChecker(latDat,
                                                               Comms(),
                                                               simState,
                                                               *cache,
                                                               *timings,
                                                               10.0); // Will accept a max/min of (21.45, 12) but not (100,1)

            AdvanceActorOneTimeStep(incompChecker);

            // This time the current leaf node reported global minimum and maximum values
            CPPUNIT_ASSERT_DOUBLES_EQUAL(smallestDefaultDensity,
                                         incompChecker.GetGlobalSmallestDensity(),
                                         eps);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(largestDefaultDensity,
                                         incompChecker.GetGlobalLargestDensity(),
                                         eps);
            CPPUNIT_ASSERT_DOUBLES_EQUAL( (largestDefaultDensity - smallestDefaultDensity)
                                             / hemelb::lb::REFERENCE_DENSITY,
                                         incompChecker.GetMaxRelativeDensityDifference(),
                                         eps);
            CPPUNIT_ASSERT(incompChecker.IsDensityDiffWithinRange());
            CPPUNIT_ASSERT_DOUBLES_EQUAL(largestDefaultVelocityMagnitude,
                                         incompChecker.GetGlobalLargestVelocityMagnitude(),
                                         eps);
          }

        private:
          lb::MacroscopicPropertyCache* cache;
          distribn_t smallestDefaultDensity;
          distribn_t largestDefaultDensity;
          distribn_t largestDefaultVelocityMagnitude;
          hemelb::reporting::Timers* timings;
          distribn_t eps;
      };

      CPPUNIT_TEST_SUITE_REGISTRATION (IncompressibilityCheckerTests);
    }
  }
}

#endif /* HEMELB_UNITTESTS_LBTESTS_IMCOMPRESSIBILITYCHECKERTESTS_H */
