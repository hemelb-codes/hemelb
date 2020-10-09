
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <catch2/catch.hpp>

#include "lb/IncompressibilityChecker.h"
#include "timestep/TimeStepManager.h"

#include "tests/lb/LbTestsHelper.h"
#include "tests/reporting/Mocks.h"
#include "tests/helpers/FourCubeLatticeData.h"
#include "tests/helpers/FourCubeBasedTestFixture.h"

namespace hemelb
{
  namespace tests
  {
    TEST_CASE_METHOD(helpers::FourCubeBasedTestFixture, "IncompressibilityCheckerTests") {
      using LATTICE = lb::lattices::D3Q15;

      LbTestsHelper::InitialiseAnisotropicTestData<LATTICE>(latDat);

      auto cache = std::make_unique<lb::MacroscopicPropertyCache>(*simState, *latDat);
      cache->densityCache.SetRefreshFlag();
      cache->velocityCache.SetRefreshFlag();
      LbTestsHelper::UpdatePropertyCache<LATTICE>(*latDat, *cache, *simState);

      // These are the smallest and largest density values in
      // FourCubeLatticeData by default. The lattice class below must
      // be consistent with the one used in
      // FourCubeLatticeData.
      // 
      // TODO: Consider templating FourCubeLatticeData over lattice
      // class, so both can be controlled from the test.
      distribn_t numDirections = (distribn_t) LATTICE::NUMVECTORS;
      distribn_t numSites = (distribn_t) latDat->GetLocalFluidSiteCount();

      // = sum_{j=1}^{numDirections} j/10 = 12 with current configuration of FourCubeLatticeData
      distribn_t smallestDefaultDensity = numDirections * (numDirections + 1) / 20;
      // = sum_{j=1}^{numDirections} j/10 + (numSites-1)/100 = 21.45 with current configuration of FourCubeLatticeData
      distribn_t largestDefaultDensity = (numDirections * (numDirections + 1) / 20) + ( (numSites - 1) * numDirections / 100);
      distribn_t largestDefaultVelocityMagnitude = 0.0433012701892219; // with current configuration of FourCubeLatticeData
      auto apprx = [](double x) {
	return Approx(x).margin(1e-9);
      };

      auto timings = std::make_unique<hemelb::reporting::Timers>(Comms());

      auto AdvanceActorOneTimeStep = [&](timestep::Actor& actor) {
	timestep::TimeStepManager tsm(1);
	tsm.AddToPhase(0, &actor);
	tsm.DoStep();
      };
      
      SECTION("IncompressibilityCheckerRootNode") {
	hemelb::lb::IncompressibilityChecker incompChecker(latDat,
							   Comms(),
							   simState.get(),
							   *cache,
							   *timings,
							   10.0); // Will accept a max/min of (21.45, 12) but not (100,1)

	// Not available until the first broadcast has finished
	AdvanceActorOneTimeStep(incompChecker);

	REQUIRE(apprx(smallestDefaultDensity) == incompChecker.GetGlobalSmallestDensity());
	REQUIRE(apprx(largestDefaultDensity) == incompChecker.GetGlobalLargestDensity());
	REQUIRE(apprx( (largestDefaultDensity - smallestDefaultDensity)
		       / hemelb::lb::REFERENCE_DENSITY) ==
		incompChecker.GetMaxRelativeDensityDifference());
	REQUIRE(incompChecker.IsDensityDiffWithinRange());
	REQUIRE(apprx(largestDefaultVelocityMagnitude) == incompChecker.GetGlobalLargestVelocityMagnitude());

	// Insert values at 1 & 100
	cache->densityCache.Put(0, 1.0);
	cache->densityCache.Put(1, 100.0);
	cache->velocityCache.Put(2, LatticeVelocity(10.0, 0., 0.0));
	AdvanceActorOneTimeStep(incompChecker);

	REQUIRE(apprx(1.0) == incompChecker.GetGlobalSmallestDensity());
	REQUIRE(apprx(100.0) == incompChecker.GetGlobalLargestDensity());
	REQUIRE(apprx(99.0) == incompChecker.GetMaxRelativeDensityDifference());
	REQUIRE(!incompChecker.IsDensityDiffWithinRange());
	REQUIRE(apprx(10.0) == incompChecker.GetGlobalLargestVelocityMagnitude());

	// The previous values are not reported by any children anymore. Testing that the checker remembers them
	AdvanceActorOneTimeStep(incompChecker);

	REQUIRE(apprx(1.0) == incompChecker.GetGlobalSmallestDensity());
	REQUIRE(apprx(100.0) == incompChecker.GetGlobalLargestDensity());
	REQUIRE(apprx(99.0) == incompChecker.GetMaxRelativeDensityDifference());
	REQUIRE(!incompChecker.IsDensityDiffWithinRange());
	REQUIRE(apprx(10.0) == incompChecker.GetGlobalLargestVelocityMagnitude());
      }

      SECTION("IncompressibilityCheckerLeafNode") {
	lb::IncompressibilityChecker incompChecker(latDat,
						   Comms(),
						   simState.get(),
						   *cache,
						   *timings,
						   10.0); // Will accept a max/min of (21.45, 12) but not (100,1)

	AdvanceActorOneTimeStep(incompChecker);

	// This time the current leaf node reported global minimum and maximum values
	REQUIRE(apprx(smallestDefaultDensity) == incompChecker.GetGlobalSmallestDensity());
	REQUIRE(apprx(largestDefaultDensity) == incompChecker.GetGlobalLargestDensity());
	REQUIRE(apprx( (largestDefaultDensity - smallestDefaultDensity)
		       / hemelb::lb::REFERENCE_DENSITY) ==
				      incompChecker.GetMaxRelativeDensityDifference());
	REQUIRE(incompChecker.IsDensityDiffWithinRange());
	REQUIRE(apprx(largestDefaultVelocityMagnitude) ==
				     incompChecker.GetGlobalLargestVelocityMagnitude());

      }
    }
  }
}
