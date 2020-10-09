
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <catch2/catch.hpp>

#include "lb/kernels/Kernels.h"
#include "util/UnitConverter.h"

#include "tests/helpers/FourCubeBasedTestFixture.h"
#include "tests/lb/LbTestsHelper.h"

namespace hemelb
{
  namespace tests
  {
    // Class to test the collision operators. These tests are for the
    // functions involved in calculating the post-collision values,
    // specifically CalculatePreCollision and Collide.  For each
    // collision operator, we test that these functions return the
    // expected values of hydrodynamic quantities and f-distribution
    // values.
    // 
    // Note that we are only testing collision operators here, so we
    // can assume that the kernel objects work perfectly.
    TEST_CASE_METHOD(public helpers::FourCubeBasedTestFixture, "Collisions") {
      const distribn_t allowedError = 1e-10;
      auto apprx = [&](double x) {
	return Approx(x).margin(allowedError);
      };

      // Initialise the fOld and the hydro vars.
      distribn_t fOld[lb::lattices::D3Q15::NUMVECTORS];
      LbTestsHelper::InitialiseAnisotropicTestData<lb::lattices::D3Q15>(0, fOld);
      lb::kernels::HydroVars<lb::kernels::LBGK<lb::lattices::D3Q15> > hydroVars(fOld);

      SECTION("TestNormal") {
	lb::collisions::Normal<lb::kernels::LBGK<lb::lattices::D3Q15> > normal(initParams);

	// Test the pre-collision step, which should calculate the
	// correct post-collisional density, velocity and equilibrium
	// distribution.
	normal.CalculatePreCollision(hydroVars, latDat->GetSite(0));

	// Calculate the expected density, velocity and f_eq.
	distribn_t expectedRho;
	util::Vector3D<distribn_t> expectedMomentum;

	LbTestsHelper::CalculateRhoMomentum<lb::lattices::D3Q15>(fOld, expectedRho, expectedMomentum);

	distribn_t expectedFeq[lb::lattices::D3Q15::NUMVECTORS];
	LbTestsHelper::CalculateLBGKEqmF<lb::lattices::D3Q15>(expectedRho,
							      expectedMomentum,
							      expectedFeq);

	// Compare.
	LbTestsHelper::CompareHydros(expectedRho,
				     expectedMomentum,
				     expectedFeq,
				     hydroVars,
				     allowedError);

	// Next, compare the collision function itself. The result
	// should be the equilibrium distribution.  Make a copy for
	// the second collision to compare against.
	lb::kernels::HydroVars<lb::kernels::LBGK<lb::lattices::D3Q15> > hydroVarsCopy(hydroVars);

	lb::kernels::LBGK<lb::lattices::D3Q15> lbgk(initParams);
	lbgk.Collide(lbmParams, hydroVars);
	normal.Collide(lbmParams, hydroVarsCopy);

	for (unsigned int ii = 0; ii < lb::lattices::D3Q15::NUMVECTORS; ++ii) {
	  REQUIRE(hydroVars.GetFPostCollision()[ii] == apprx(hydroVarsCopy.GetFPostCollision()[ii]));
	}
      }
    }	
  }
}

