// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <catch2/catch.hpp>

#include "lb/kernels/Kernels.h"
#include "util/UnitConverter.h"

#include "tests/helpers/FourCubeBasedTestFixture.h"
#include "tests/lb/LbTestsHelper.h"

namespace hemelb::tests
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
    TEST_CASE_METHOD(public helpers::FourCubeBasedTestFixture<>, "Collisions") {
      const distribn_t allowedError = 1e-10;
      auto apprx = [&](double x) {
	return Approx(x).margin(allowedError);
      };

      // Initialise the fOld and the hydro vars.
      distribn_t fOld[lb::D3Q15::NUMVECTORS];
      LbTestsHelper::InitialiseAnisotropicTestData<lb::D3Q15>(0, fOld);
      lb::HydroVars<lb::LBGK<lb::D3Q15> > hydroVars(fOld);

      SECTION ("TestNonZeroVelocityEquilibriumFixedDensity") {
	lb::BoundaryValues inletBoundary = BuildIolets(geometry::INLET_TYPE);
	initParams.boundaryObject = &inletBoundary;

	lb::NonZeroVelocityEquilibriumFixedDensity<lb::LBGK<lb::D3Q15> >
	  nonZeroVFixedDensityILet(initParams);

	// Test the pre-collision step, which should calculate the correct
	// post-collisional density, velocity and equilibrium distribution.
	geometry::Site<geometry::FieldData> dummySite(0, *latDat);
	dom->SetIoletId(0, 0);

	nonZeroVFixedDensityILet.CalculatePreCollision(hydroVars, dummySite);

	// Calculate the expected density, velocity and f_eq.
	distribn_t expectedRho = inletBoundary.GetBoundaryDensity(0);
	util::Vector3D<distribn_t> expectedMomentum;

	distribn_t originalRho;
	LbTestsHelper::CalculateRhoMomentum<lb::D3Q15>(fOld, originalRho, expectedMomentum);

	// Now need to scale the momentum, expectedV, to account for the difference between
	// original and enforced densities.
	expectedMomentum *= (expectedRho / originalRho);

	distribn_t expectedFeq[lb::D3Q15::NUMVECTORS];
	LbTestsHelper::CalculateLBGKEqmF<lb::D3Q15>(expectedRho,
							      expectedMomentum,
							      expectedFeq);

	// Compare.
	LbTestsHelper::CompareHydros(expectedRho,
				     expectedMomentum,
				     expectedFeq,
				     hydroVars,
				     allowedError);

	// Next, compare the collision function itself. The result should be the equilibrium
	// distribution.
	nonZeroVFixedDensityILet.Collide(&lbmParams, hydroVars);

	for (unsigned int ii = 0; ii < lb::D3Q15::NUMVECTORS; ++ii) {
	  REQUIRE(apprx(expectedFeq[ii]) == hydroVars.GetFPostCollision()[ii]);
	}
      }

      SECTION("TestZeroVelocityEquilibriumFixedDensity") {
	lb::BoundaryValues outletBoundary = BuildIolets(geometry::OUTLET_TYPE);
	initParams.boundaryObject = &outletBoundary;

	lb::ZeroVelocityEquilibriumFixedDensity<lb::LBGK<lb::D3Q15> >
	  zeroVFixedDensityOLet(initParams);

	// Test the pre-collision step, which should calculate the
	// correct post-collisional density, velocity and equilibrium
	// distribution.
	dom->SetIoletId(0, 0);
	zeroVFixedDensityOLet.CalculatePreCollision(hydroVars, latDat->GetSite(0));

	// Calculate the expected density, velocity and f_eq.
	distribn_t expectedRho = outletBoundary.GetBoundaryDensity(0);
	auto expectedMomentum = util::Vector3D<distribn_t>::Zero();

	distribn_t expectedFeq[lb::D3Q15::NUMVECTORS];
	LbTestsHelper::CalculateLBGKEqmF<lb::D3Q15>(expectedRho,
							      expectedMomentum,
							      expectedFeq);

	// Compare.
	LbTestsHelper::CompareHydros(expectedRho,
				     expectedMomentum,
				     expectedFeq,
				     hydroVars,
				     allowedError);

	// Next, compare the collision function itself. The result
	// should be the equilibrium distribution.
	zeroVFixedDensityOLet.Collide(&lbmParams, hydroVars);

	for (unsigned int ii = 0; ii < lb::D3Q15::NUMVECTORS; ++ii)
	{
	  REQUIRE(apprx(expectedFeq[ii]) == hydroVars.GetFPostCollision()[ii]);
	}
      }

      SECTION("TestZeroVelocityEquilibrium") {
	lb::ZeroVelocityEquilibrium<lb::LBGK<lb::D3Q15> > zeroVEqm(initParams);

	// Test the pre-collision step, which should calculate the
	// correct post-collisional density, velocity and equilibrium
	// distribution.
	zeroVEqm.CalculatePreCollision(hydroVars, latDat->GetSite(0));

	// Calculate the expected density, velocity and f_eq.
	distribn_t expectedRho = 0.0;
	auto expectedMomentum = util::Vector3D<distribn_t>::Zero();

	for (unsigned int ii = 0; ii < lb::D3Q15::NUMVECTORS; ++ii) {
	  expectedRho += fOld[ii];
	}

	distribn_t expectedFeq[lb::D3Q15::NUMVECTORS];
	LbTestsHelper::CalculateLBGKEqmF<lb::D3Q15>(expectedRho,
							      expectedMomentum,
							      expectedFeq);

	// Compare.
	LbTestsHelper::CompareHydros(expectedRho,
				     expectedMomentum,
				     expectedFeq,
				     hydroVars,
				     allowedError);

	// Next, compare the collision function itself. The result
	// should be the equilibrium distribution.
	zeroVEqm.Collide(&lbmParams, hydroVars);

	for (unsigned int ii = 0; ii < lb::D3Q15::NUMVECTORS; ++ii)
	{
	  REQUIRE(apprx(expectedFeq[ii]) == hydroVars.GetFPostCollision()[ii]);
	}
      }

      SECTION("TestNormal") {
	lb::Normal<lb::LBGK<lb::D3Q15> > normal(initParams);

	// Test the pre-collision step, which should calculate the
	// correct post-collisional density, velocity and equilibrium
	// distribution.
	normal.CalculatePreCollision(hydroVars, latDat->GetSite(0));

	// Calculate the expected density, velocity and f_eq.
	distribn_t expectedRho;
	util::Vector3D<distribn_t> expectedMomentum;

	LbTestsHelper::CalculateRhoMomentum<lb::D3Q15>(fOld, expectedRho, expectedMomentum);

	distribn_t expectedFeq[lb::D3Q15::NUMVECTORS];
	LbTestsHelper::CalculateLBGKEqmF<lb::D3Q15>(expectedRho,
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
	lb::HydroVars<lb::LBGK<lb::D3Q15> > hydroVarsCopy(hydroVars);

	lb::LBGK<lb::D3Q15> lbgk(initParams);
	lbgk.Collide(&lbmParams, hydroVars);
	normal.Collide(&lbmParams, hydroVarsCopy);

	for (unsigned int ii = 0; ii < lb::D3Q15::NUMVECTORS; ++ii) {
	  REQUIRE(hydroVars.GetFPostCollision()[ii] == apprx(hydroVarsCopy.GetFPostCollision()[ii]));
	}
      }
    }
}

