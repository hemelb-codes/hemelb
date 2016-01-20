
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_UNITTESTS_LBTESTS_COLLISIONTESTS_H
#define HEMELB_UNITTESTS_LBTESTS_COLLISIONTESTS_H

#include <cppunit/TestFixture.h>

#include "unittests/helpers/FourCubeBasedTestFixture.h"
#include "util/UnitConverter.h"

namespace hemelb
{
  namespace unittests
  {
    namespace lbtests
    {
      /**
       * Class to test the collision operators. These tests are for the functions involved in
       * calculating the post-collision values, specifically CalculatePreCollision and Collide.
       * For each collision operator, we test that these functions return the expected values of
       * hydrodynamic quantities and f-distribution values.
       *
       * Note that we are only testing collision operators here, so we
       * can assume that the kernel objects work perfectly.
       */
      class CollisionTests : public helpers::FourCubeBasedTestFixture
      {
          CPPUNIT_TEST_SUITE( CollisionTests);
          CPPUNIT_TEST( TestNonZeroVelocityEquilibriumFixedDensity);
          CPPUNIT_TEST( TestZeroVelocityEquilibriumFixedDensity);
          CPPUNIT_TEST( TestZeroVelocityEquilibrium);
          CPPUNIT_TEST( TestNormal);CPPUNIT_TEST_SUITE_END();
        public:
          void TestNonZeroVelocityEquilibriumFixedDensity()
          {
            lb::iolets::BoundaryValues inletBoundary(geometry::INLET_TYPE,
                                                         latDat,
                                                         simConfig->GetInlets(),
                                                         simState,
                                                         Comms(),
                                                         *unitConverter);

            initParams.boundaryObject = &inletBoundary;

            lb::collisions::NonZeroVelocityEquilibriumFixedDensity<lb::kernels::LBGK<lb::lattices::D3Q15> >
                nonZeroVFixedDensityILet(initParams);

            distribn_t allowedError = 1e-10;

            // Initialise the fOld and the hydro vars.
            distribn_t fOld[lb::lattices::D3Q15::NUMVECTORS];

            LbTestsHelper::InitialiseAnisotropicTestData<lb::lattices::D3Q15>(0, fOld);

            lb::kernels::HydroVars<lb::kernels::LBGK<lb::lattices::D3Q15> > hydroVars(fOld);

            // Test the pre-collision step, which should calculate the correct
            // post-collisional density, velocity and equilibrium distribution.
            geometry::Site<geometry::LatticeData> dummySite(0, *latDat);
            latDat->SetIoletId(0, 0);

            nonZeroVFixedDensityILet.CalculatePreCollision(hydroVars, dummySite);

            // Calculate the expected density, velocity and f_eq.
            distribn_t expectedRho = inletBoundary.GetBoundaryDensity(0);
            distribn_t expectedMomentum[3];

            distribn_t originalRho;
            LbTestsHelper::CalculateRhoMomentum<lb::lattices::D3Q15>(fOld, originalRho, expectedMomentum);

            // Now need to scale the momentum, expectedV, to account for the difference between
            // original and enforced densities.
            for (unsigned axis = 0; axis < 3; ++axis)
              expectedMomentum[axis] *= (expectedRho / originalRho);

            distribn_t expectedFeq[lb::lattices::D3Q15::NUMVECTORS];
            LbTestsHelper::CalculateLBGKEqmF<lb::lattices::D3Q15>(expectedRho,
                                                                  expectedMomentum[0],
                                                                  expectedMomentum[1],
                                                                  expectedMomentum[2],
                                                                  expectedFeq);

            // Compare.
            LbTestsHelper::CompareHydros(expectedRho,
                                         expectedMomentum[0],
                                         expectedMomentum[1],
                                         expectedMomentum[2],
                                         expectedFeq,
                                         "Non-0 velocity eqm fixed density, calculate pre collision",
                                         hydroVars,
                                         allowedError);

            // Next, compare the collision function itself. The result should be the equilibrium
            // distribution.
            nonZeroVFixedDensityILet.Collide(lbmParams, hydroVars);

            for (unsigned int ii = 0; ii < lb::lattices::D3Q15::NUMVECTORS; ++ii)
            {
              CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Non-0 velocity eqm fixed density, collide",
                                                   expectedFeq[ii],
                                                   hydroVars.GetFPostCollision()[ii],
                                                   allowedError);
            }
          }

          void TestZeroVelocityEquilibriumFixedDensity()
          {
            lb::iolets::BoundaryValues outletBoundary(geometry::OUTLET_TYPE,
                                                          latDat,
                                                          simConfig->GetOutlets(),
                                                          simState,
                                                          Comms(),
                                                          *unitConverter);
            initParams.boundaryObject = &outletBoundary;

            lb::collisions::ZeroVelocityEquilibriumFixedDensity<lb::kernels::LBGK<lb::lattices::D3Q15> >
                zeroVFixedDensityOLet(initParams);

            distribn_t allowedError = 1e-10;

            // Initialise the fOld and the hydro vars.
            distribn_t fOld[lb::lattices::D3Q15::NUMVECTORS];

            LbTestsHelper::InitialiseAnisotropicTestData<lb::lattices::D3Q15>(0, fOld);

            lb::kernels::HydroVars<lb::kernels::LBGK<lb::lattices::D3Q15> > hydroVars(fOld);

            // Test the pre-collision step, which should calculate the correct
            // post-collisional density, velocity and equilibrium distribution.
            latDat->SetIoletId(0, 0);
            zeroVFixedDensityOLet.CalculatePreCollision(hydroVars, latDat->GetSite(0));

            // Calculate the expected density, velocity and f_eq.
            distribn_t expectedRho = outletBoundary.GetBoundaryDensity(0);
            distribn_t expectedMomentum[3] = { 0., 0., 0. };

            distribn_t expectedFeq[lb::lattices::D3Q15::NUMVECTORS];
            LbTestsHelper::CalculateLBGKEqmF<lb::lattices::D3Q15>(expectedRho,
                                                                  expectedMomentum[0],
                                                                  expectedMomentum[1],
                                                                  expectedMomentum[2],
                                                                  expectedFeq);

            // Compare.
            LbTestsHelper::CompareHydros(expectedRho,
                                         expectedMomentum[0],
                                         expectedMomentum[1],
                                         expectedMomentum[2],
                                         expectedFeq,
                                         "0 velocity eqm fixed density, calculate pre collision",
                                         hydroVars,
                                         allowedError);

            // Next, compare the collision function itself. The result should be the equilibrium
            // distribution.
            zeroVFixedDensityOLet.Collide(lbmParams, hydroVars);

            for (unsigned int ii = 0; ii < lb::lattices::D3Q15::NUMVECTORS; ++ii)
            {
              CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("0 velocity eqm fixed density, collide",
                                                   expectedFeq[ii],
                                                   hydroVars.GetFPostCollision()[ii],
                                                   allowedError);
            }
          }

          void TestZeroVelocityEquilibrium()
          {
            lb::collisions::ZeroVelocityEquilibrium<lb::kernels::LBGK<lb::lattices::D3Q15> > zeroVEqm(initParams);

            distribn_t allowedError = 1e-10;

            // Initialise the fOld and the hydro vars.
            distribn_t fOld[lb::lattices::D3Q15::NUMVECTORS];

            LbTestsHelper::InitialiseAnisotropicTestData<lb::lattices::D3Q15>(0, fOld);

            lb::kernels::HydroVars<lb::kernels::LBGK<lb::lattices::D3Q15> > hydroVars(fOld);

            // Test the pre-collision step, which should calculate the correct
            // post-collisional density, velocity and equilibrium distribution.
            zeroVEqm.CalculatePreCollision(hydroVars, latDat->GetSite(0));

            // Calculate the expected density, velocity and f_eq.
            distribn_t expectedRho = 0.0;
            distribn_t expectedMomentum[3] = { 0., 0., 0. };

            for (unsigned int ii = 0; ii < lb::lattices::D3Q15::NUMVECTORS; ++ii)
            {
              expectedRho += fOld[ii];
            }

            distribn_t expectedFeq[lb::lattices::D3Q15::NUMVECTORS];
            LbTestsHelper::CalculateLBGKEqmF<lb::lattices::D3Q15>(expectedRho,
                                                                  expectedMomentum[0],
                                                                  expectedMomentum[1],
                                                                  expectedMomentum[2],
                                                                  expectedFeq);

            // Compare.
            LbTestsHelper::CompareHydros(expectedRho,
                                         expectedMomentum[0],
                                         expectedMomentum[1],
                                         expectedMomentum[2],
                                         expectedFeq,
                                         "0 velocity eqm, calculate pre collision",
                                         hydroVars,
                                         allowedError);

            // Next, compare the collision function itself. The result should be the equilibrium
            // distribution.
            zeroVEqm.Collide(lbmParams, hydroVars);

            for (unsigned int ii = 0; ii < lb::lattices::D3Q15::NUMVECTORS; ++ii)
            {
              CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("0 velocity eqm, collide",
                                                   expectedFeq[ii],
                                                   hydroVars.GetFPostCollision()[ii],
                                                   allowedError);
            }
          }

          void TestNormal()
          {
            lb::collisions::Normal<lb::kernels::LBGK<lb::lattices::D3Q15> > normal(initParams);

            distribn_t allowedError = 1e-10;

            // Initialise the fOld and the hydro vars.
            distribn_t fOld[lb::lattices::D3Q15::NUMVECTORS];

            LbTestsHelper::InitialiseAnisotropicTestData<lb::lattices::D3Q15>(0, fOld);

            lb::kernels::HydroVars<lb::kernels::LBGK<lb::lattices::D3Q15> > hydroVars(fOld);

            // Test the pre-collision step, which should calculate the correct
            // post-collisional density, velocity and equilibrium distribution.
            normal.CalculatePreCollision(hydroVars, latDat->GetSite(0));

            // Calculate the expected density, velocity and f_eq.
            distribn_t expectedRho;
            distribn_t expectedMomentum[3];

            LbTestsHelper::CalculateRhoMomentum<lb::lattices::D3Q15>(fOld, expectedRho, expectedMomentum);

            distribn_t expectedFeq[lb::lattices::D3Q15::NUMVECTORS];
            LbTestsHelper::CalculateLBGKEqmF<lb::lattices::D3Q15>(expectedRho,
                                                                  expectedMomentum[0],
                                                                  expectedMomentum[1],
                                                                  expectedMomentum[2],
                                                                  expectedFeq);

            // Compare.
            LbTestsHelper::CompareHydros(expectedRho,
                                         expectedMomentum[0],
                                         expectedMomentum[1],
                                         expectedMomentum[2],
                                         expectedFeq,
                                         "Normal, calculate pre collision",
                                         hydroVars,
                                         allowedError);

            // Next, compare the collision function itself. The result should be the equilibrium
            // distribution.
            // Make a copy for the second collision to compare against.
            lb::kernels::HydroVars<lb::kernels::LBGK<lb::lattices::D3Q15> > hydroVarsCopy(hydroVars);

            lb::kernels::LBGK<lb::lattices::D3Q15> lbgk(initParams);
            lbgk.Collide(lbmParams, hydroVars);
            normal.Collide(lbmParams, hydroVarsCopy);

            for (unsigned int ii = 0; ii < lb::lattices::D3Q15::NUMVECTORS; ++ii)
            {
              CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Normal, collide",
                                                   hydroVars.GetFPostCollision()[ii],
                                                   hydroVarsCopy.GetFPostCollision()[ii],
                                                   allowedError);
            }
          }
      };
      CPPUNIT_TEST_SUITE_REGISTRATION( CollisionTests);
    }
  }
}

#endif /* HEMELB_UNITTESTS_LBTESTS_COLLISIONTESTS_H_ */
