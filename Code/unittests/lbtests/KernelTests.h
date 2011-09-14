#ifndef HEMELB_UNITTESTS_LBTESTS_KERNELTESTS_H
#define HEMELB_UNITTESTS_LBTESTS_KERNELTESTS_H

#include <cppunit/TestFixture.h>
#include <cstring>
#include <sstream>

#include "lb/kernels/Entropic.h"
#include "lb/kernels/LBGK.h"
#include "unittests/lbtests/LbTestsHelper.h"

namespace hemelb
{
  namespace unittests
  {
    namespace lbtests
    {

      /**
       * Class containing tests for the functionality of the lattice-Boltzmann kernels.
       */
      class KernelTests : public CppUnit::TestFixture
      {
        public:
          void setUp()
          {
            // Initialise both kernels.
            lb::kernels::InitParams initParams;

            initParams.siteCount = 1;

            entropic = new lb::kernels::Entropic(initParams);
            lbgk = new lb::kernels::LBGK(initParams);

            // Initialise the LBM parameters.
            int timeStepsPerCycle = 1000;
            double voxelSize = 0.01;

            lbmParams = new lb::LbmParameters(PULSATILE_PERIOD_s / (distribn_t) timeStepsPerCycle,
                                              voxelSize);
          }

          void tearDown()
          {
            delete entropic;
            delete lbgk;
            delete lbmParams;
          }

          void TestEntropic()
          {
            // Initialise the original f distribution to something asymmetric.
            distribn_t f_original[D3Q15::NUMVECTORS];

            for (unsigned int ii = 0; ii < D3Q15::NUMVECTORS; ++ii)
            {
              f_original[ii] = ((float) (1 + ii)) / 10.0;
            }

            lb::kernels::HydroVars<lb::kernels::Entropic> hydroVars0(f_original);
            lb::kernels::HydroVars<lb::kernels::Entropic> hydroVars1(f_original);

            // Calculate density, velocity, equilibrium f.
            entropic->CalculateDensityVelocityFeq(hydroVars0, 0);

            // Manually set density and velocity and calculate eqm f.
            hydroVars1.density = 1.0;
            hydroVars1.v_x = 0.4;
            hydroVars1.v_y = 0.5;
            hydroVars1.v_z = 0.6;

            entropic->CalculateFeq(hydroVars1, 1);

            // Calculate expected values.
            distribn_t expectedDensity0 = 12.0; // (sum 1 to 15) / 10
            distribn_t expectedDensity1 = 1.0; // Unchanged

            distribn_t expectedVelocity0[3];
            LbTestsHelper::CalculateVelocity<D3Q15>(hydroVars0.f, expectedVelocity0);
            distribn_t expectedVelocity1[3] = { 0.4, 0.5, 0.6 };

            distribn_t expectedFEq0[D3Q15::NUMVECTORS];
            LbTestsHelper::CalculateEntropicEqmF<D3Q15>(expectedDensity0,
                                                            expectedVelocity0[0],
                                                            expectedVelocity0[1],
                                                            expectedVelocity0[2],
                                                            expectedFEq0);
            distribn_t expectedFEq1[D3Q15::NUMVECTORS];
            LbTestsHelper::CalculateEntropicEqmF<D3Q15>(expectedDensity1,
                                                            expectedVelocity1[0],
                                                            expectedVelocity1[1],
                                                            expectedVelocity1[2],
                                                            expectedFEq1);

            // Now compare the expected and actual values.
            distribn_t allowedError = 1e-10;

            LbTestsHelper::CompareHydros(expectedDensity0,
                                             expectedVelocity0[0],
                                             expectedVelocity0[1],
                                             expectedVelocity0[2],
                                             expectedFEq0,
                                             "Entropic, case 0",
                                             hydroVars0,
                                             allowedError);
            LbTestsHelper::CompareHydros(expectedDensity1,
                                             expectedVelocity1[0],
                                             expectedVelocity1[1],
                                             expectedVelocity1[2],
                                             expectedFEq1,
                                             "Entropic, case 1",
                                             hydroVars1,
                                             allowedError);

            // Do the collision and test the result.
            distribn_t postCollision0[D3Q15::NUMVECTORS];
            distribn_t postCollision1[D3Q15::NUMVECTORS];

            // Set the values in f_neq.
            for (unsigned int ii = 0; ii < D3Q15::NUMVECTORS; ++ii)
            {
              hydroVars0.f_neq[ii] = f_original[ii] - hydroVars0.f_eq[ii];
              hydroVars1.f_neq[ii] = f_original[ii] - hydroVars1.f_eq[ii];
            }

            for (unsigned int ii = 0; ii < D3Q15::NUMVECTORS; ++ii)
            {
              postCollision0[ii] = entropic->DoCollide(lbmParams, hydroVars0, ii);
              postCollision1[ii] = entropic->DoCollide(lbmParams, hydroVars1, ii);
            }

            // Get the expected post-collision densities.
            distribn_t expectedPostCollision0[D3Q15::NUMVECTORS];
            distribn_t expectedPostCollision1[D3Q15::NUMVECTORS];

            LbTestsHelper::CalculateEntropicCollision<D3Q15>(f_original,
                                                                 hydroVars0.f_eq,
                                                                 lbmParams->Tau(),
                                                                 lbmParams->Beta(),
                                                                 expectedPostCollision0);

            LbTestsHelper::CalculateEntropicCollision<D3Q15>(f_original,
                                                                 hydroVars1.f_eq,
                                                                 lbmParams->Tau(),
                                                                 lbmParams->Beta(),
                                                                 expectedPostCollision1);

            // Compare.
            for (unsigned int ii = 0; ii < D3Q15::NUMVECTORS; ++ii)
            {
              std::stringstream message("Post-collision ");
              message << ii;

              CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE(message.str(),
                                                   postCollision0[ii],
                                                   expectedPostCollision0[ii],
                                                   allowedError);

              CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE(message.str(),
                                                   postCollision1[ii],
                                                   expectedPostCollision1[ii],
                                                   allowedError);
            }
          }

          void TestLBGK()
          {
            // Initialise the original f distribution to something asymmetric.
            distribn_t f_original[D3Q15::NUMVECTORS];

            for (unsigned int ii = 0; ii < D3Q15::NUMVECTORS; ++ii)
            {
              f_original[ii] = ((float) (1 + ii)) / 10.0;
            }

            lb::kernels::HydroVars<lb::kernels::LBGK> hydroVars0(f_original);
            lb::kernels::HydroVars<lb::kernels::LBGK> hydroVars1(f_original);

            // Calculate density, velocity, equilibrium f.
            lbgk->CalculateDensityVelocityFeq(hydroVars0, 0);

            // Manually set density and velocity and calculate eqm f.
            hydroVars1.density = 1.0;
            hydroVars1.v_x = 0.4;
            hydroVars1.v_y = 0.5;
            hydroVars1.v_z = 0.6;

            lbgk->CalculateFeq(hydroVars1, 1);

            // Calculate expected values.
            distribn_t expectedDensity0 = 12.0; // (sum 1 to 15) / 10
            distribn_t expectedDensity1 = 1.0; // Unchanged

            distribn_t expectedVelocity0[3];
            LbTestsHelper::CalculateVelocity<D3Q15>(hydroVars0.f, expectedVelocity0);
            distribn_t expectedVelocity1[3] = { 0.4, 0.5, 0.6 };

            distribn_t expectedFEq0[D3Q15::NUMVECTORS];
            LbTestsHelper::CalculateLBGKEqmF<D3Q15>(expectedDensity0,
                                                        expectedVelocity0[0],
                                                        expectedVelocity0[1],
                                                        expectedVelocity0[2],
                                                        expectedFEq0);
            distribn_t expectedFEq1[D3Q15::NUMVECTORS];
            LbTestsHelper::CalculateLBGKEqmF<D3Q15>(expectedDensity1,
                                                        expectedVelocity1[0],
                                                        expectedVelocity1[1],
                                                        expectedVelocity1[2],
                                                        expectedFEq1);

            // Now compare the expected and actual values.
            distribn_t allowedError = 1e-10;

            LbTestsHelper::CompareHydros(expectedDensity0,
                                             expectedVelocity0[0],
                                             expectedVelocity0[1],
                                             expectedVelocity0[2],
                                             expectedFEq0,
                                             "LBGK, case 0",
                                             hydroVars0,
                                             allowedError);
            LbTestsHelper::CompareHydros(expectedDensity1,
                                             expectedVelocity1[0],
                                             expectedVelocity1[1],
                                             expectedVelocity1[2],
                                             expectedFEq1,
                                             "LBGK, case 1",
                                             hydroVars1,
                                             allowedError);

            // Do the collision and test the result.
            distribn_t postCollision0[D3Q15::NUMVECTORS];
            distribn_t postCollision1[D3Q15::NUMVECTORS];

            // Set the values in f_neq.
            for (unsigned int ii = 0; ii < D3Q15::NUMVECTORS; ++ii)
            {
              hydroVars0.f_neq[ii] = f_original[ii] - hydroVars0.f_eq[ii];
              hydroVars1.f_neq[ii] = f_original[ii] - hydroVars1.f_eq[ii];
            }

            for (unsigned int ii = 0; ii < D3Q15::NUMVECTORS; ++ii)
            {
              postCollision0[ii] = lbgk->DoCollide(lbmParams, hydroVars0, ii);
              postCollision1[ii] = lbgk->DoCollide(lbmParams, hydroVars1, ii);
            }

            // Get the expected post-collision densities.
            distribn_t expectedPostCollision0[D3Q15::NUMVECTORS];
            distribn_t expectedPostCollision1[D3Q15::NUMVECTORS];

            LbTestsHelper::CalculateLBGKCollision<D3Q15>(f_original,
                                                             hydroVars0.f_eq,
                                                             lbmParams->Omega(),
                                                             expectedPostCollision0);

            LbTestsHelper::CalculateLBGKCollision<D3Q15>(f_original,
                                                             hydroVars1.f_eq,
                                                             lbmParams->Omega(),
                                                             expectedPostCollision1);

            // Compare.
            for (unsigned int ii = 0; ii < D3Q15::NUMVECTORS; ++ii)
            {
              std::stringstream message("Post-collision ");
              message << ii;

              CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE(message.str(),
                                                   postCollision0[ii],
                                                   expectedPostCollision0[ii],
                                                   allowedError);

              CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE(message.str(),
                                                   postCollision1[ii],
                                                   expectedPostCollision1[ii],
                                                   allowedError);
            }
          }

        private:
          lb::LbmParameters* lbmParams;
          lb::kernels::Entropic* entropic;
          lb::kernels::LBGK* lbgk;
      };

    }
  }
}

#endif /* HEMELB_UNITTESTS_LBTESTS_KERNELTESTS_H */
