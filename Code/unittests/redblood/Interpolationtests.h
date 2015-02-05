//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELB_UNITTESTS_REDBLOOD_INTERPOLATION_H
#define HEMELB_UNITTESTS_REDBLOOD_INTERPOLATION_H

#include <cppunit/TestFixture.h>
#include "redblood/interpolation.h"
#include "redblood/VelocityInterpolation.h"
#include "lb/lattices/D3Q15.h"
#include "lb/kernels/LBGK.h"
#include "lb/kernels/GuoForcingLBGK.h"
#include "unittests/redblood/Fixtures.h"
#include "unittests/helpers/Comparisons.h"
#include "unittests/helpers/LatticeDataAccess.h"

namespace hemelb
{
  namespace unittests
  {
    namespace redblood
    {
      class InterpolationTests : public CppUnit::TestFixture
      {
          CPPUNIT_TEST_SUITE (InterpolationTests);
          CPPUNIT_TEST (testIndexIterator);
          CPPUNIT_TEST (testOffLattice);
          CPPUNIT_TEST (testOffLatticeZeroOutsideStencil);
          CPPUNIT_TEST (testInterpolateLinearFunction);
          CPPUNIT_TEST (testInterpolateQuadraticFunction);CPPUNIT_TEST_SUITE_END();

          struct PlanarFunction
          {
              LatticePosition operator()(LatticePosition const &pos) const
              {
                return LatticePosition(LatticePosition(1, 1, 1).Dot(pos),
                                       LatticePosition(-1, 2, 1).Dot(pos),
                                       LatticePosition(0, 0, 1).Dot(pos));
              }
              LatticePosition operator()(Dimensionless x, Dimensionless y, Dimensionless z) const
              {
                return operator()(LatticePosition(x, y, z));
              }
          };

          struct QuadraticFunction
          {
              LatticePosition operator()(LatticePosition const &pos) const
              {
                Dimensionless const offset(0);
                return LatticePosition( (LatticePosition(1, 1, 1).Dot(pos) - offset)
                                           * (LatticePosition(1, 1, 1).Dot(pos) - offset),
                                       (LatticePosition(0, 1, 0).Dot(pos) - offset)
                                           * (LatticePosition(0, 1, 0).Dot(pos) - offset),
                                       (LatticePosition(0, 0, 1).Dot(pos) - offset)
                                           * (LatticePosition(0, 0, 1).Dot(pos) - offset));
              }
              LatticePosition operator()(Dimensionless x, Dimensionless y, Dimensionless z) const
              {
                return operator()(LatticePosition(x, y, z));
              }
          };

        public:
          void testIndexIterator()
          {
            LatticeVector vectors[] = { LatticeVector(4, 3, 2),
                                        LatticeVector(4, 3, 3),
                                        LatticeVector(4, 3, 4),
                                        LatticeVector(4, 4, 2),
                                        LatticeVector(4, 5, 2),
                                        LatticeVector(5, 3, 2),
                                        LatticeVector(6, 3, 2),
                                        LatticeVector(6, 5, 4), };
            size_t incs[] = { 0, 1, 1, 1, 3, 3, 9, 8, 666 // break
                };

            // Checks iteration goes through correct sequence
            IndexIterator iterator(LatticeVector(5, 4, 3), 1);

            for (size_t i(0); incs[i] < 666; ++i)
            {
              for (size_t j(0); j < incs[i]; ++j, ++iterator)
                ;

              CPPUNIT_ASSERT(helpers::is_zero(*iterator - vectors[i]));
              CPPUNIT_ASSERT(iterator.IsValid());
            }

            // Checks iterator becomes invalid
            ++iterator;
            CPPUNIT_ASSERT(not iterator.IsValid());
          }

          void testOffLattice()
          {
            LatticePosition const pos(56.51, 52.9, 15.2);
            stencil::FourPoint stencil;
            InterpolationIterator iterator(pos, stencil);

            LatticeVector vectors[] = { LatticeVector(55, 51, 14),
                                        LatticeVector(55, 51, 15),
                                        LatticeVector(55, 51, 16),
                                        LatticeVector(55, 51, 17),
                                        LatticeVector(55, 52, 14),
                                        LatticeVector(56, 51, 14),
                                        LatticeVector(58, 54, 17), };
            size_t incs[] = { 0, 1, 1, 1, 1, 12, 47, 666 // break
                };

            // Checks iteration goes through correct sequence
            for (size_t i(0); incs[i] < 666; ++i)
            {
              for (size_t j(0); j < incs[i]; ++j, ++iterator)
                ;

              CPPUNIT_ASSERT(helpers::is_zero(*iterator - vectors[i]));
              CPPUNIT_ASSERT(helpers::is_zero(stencil(pos - vectors[i]) - iterator.weight()));
              CPPUNIT_ASSERT(iterator.IsValid());
            }

            ++iterator;
            CPPUNIT_ASSERT(not iterator.IsValid());
          }

          void testOffLatticeZeroOutsideStencil()
          {
            LatticePosition const pos(56.51, 52.9, 15.2);
            stencil::FourPoint stencil;
            InterpolationIterator iterator(pos, stencil);
            // Checks that outside iteration box, weights are zero
            LatticeVector zero_vecs[] = { LatticeVector(57, 53, 13),
                                          LatticeVector(57, 53, 18),
                                          LatticeVector(57, 50, 17),
                                          LatticeVector(57, 55, 17),
                                          LatticeVector(54, 53, 17),
                                          LatticeVector(59, 53, 17) };

            for (size_t i(0); i < 6; ++i)
            {
              LatticeVector const dx(1, 0, 0), dy(0, 1, 0), dz(0, 0, 1);
              CPPUNIT_ASSERT(helpers::is_zero(stencil(pos - zero_vecs[i])));
              // checks we are one step outside the iteration box only.
              // this is really a test on the zero_vecs data, eg a test of the test.
              size_t const one_non_zero = size_t(not helpers::is_zero(stencil(pos + dx
                  - zero_vecs[i]))) + size_t(not helpers::is_zero(stencil(pos - dx - zero_vecs[i])))
                  + size_t(not helpers::is_zero(stencil(pos + dy - zero_vecs[i])))
                  + size_t(not helpers::is_zero(stencil(pos - dy - zero_vecs[i])))
                  + size_t(not helpers::is_zero(stencil(pos + dz - zero_vecs[i])))
                  + size_t(not helpers::is_zero(stencil(pos - dz - zero_vecs[i])));
              CPPUNIT_ASSERT(one_non_zero == 1);
            }
          }

          template<class FUNCTION>
          void check(Dimensionless x, Dimensionless y, Dimensionless z, Dimensionless tolerance =
                         1e-8)
          {
            using hemelb::redblood::stencil::FOUR_POINT;

            FUNCTION func;
            LatticePosition expected(func(x, y, z));
            LatticePosition actual(interpolate(func, x, y, z, FOUR_POINT));
            CPPUNIT_ASSERT(helpers::is_zero(actual - expected, tolerance));
          }

          // Test interpolation when the point is on the grid
          void testInterpolateLinearFunction()
          {
            check<PlanarFunction>(0, 0, 0);
            check<PlanarFunction>(0.1, 0.5, 0.6);
            check<PlanarFunction>(-5.1, 0.5, 8.7);
            check<PlanarFunction>(-5, 0, -1);
          }

          void testInterpolateQuadraticFunction()
          {
            QuadraticFunction quad;
            // Error depends on variation on scale larger than stencil
            Dimensionless const tolerance( (quad(0, 0, 0) - quad(12, 12, 12)).GetMagnitude()
                * 1e-2);
            check<QuadraticFunction>(0, 0, 0, tolerance);
            check<QuadraticFunction>(0.1, 0.5, 0.6, tolerance);
            check<QuadraticFunction>(-5.1, 0.5, 8.7, tolerance);
            check<QuadraticFunction>(-5, 0, -1, tolerance);
          }
      };

      class VelocityInterpolationTests : public helpers::FourCubeBasedTestFixture
      {
          typedef lb::lattices::D3Q15 D3Q15;
          typedef lb::kernels::LBGK<D3Q15> LBGK;
          typedef lb::kernels::GuoForcingLBGK<D3Q15> GuoForcingLBGK;
          CPPUNIT_TEST_SUITE (VelocityInterpolationTests);
          // CPPUNIT_TEST(testLatticeDataFunctor);
          CPPUNIT_TEST (testVelocityDataFromLatticeWithForces);
          CPPUNIT_TEST (testVelocityDataFromLatticeWithoutForces);CPPUNIT_TEST_SUITE_END();

        public:
          void setUp()
          {
            FourCubeBasedTestFixture::setUp();

            LatticeVector const min(latDat->GetGlobalSiteMins());
            LatticeVector const max(latDat->GetGlobalSiteMaxes());

            for (size_t i(min[0]); i <= max[0]; ++i)
              for (size_t j(min[1]); j <= max[1]; ++j)
                for (size_t k(min[2]); k <= max[2]; ++k)
                {
                  distribn_t fold[D3Q15::NUMVECTORS];

                  for (size_t u(0); u < D3Q15::NUMVECTORS; ++u)
                  {
                    fold[u] = 0e0;
                  }

                  fold[2] = i;
                  fold[4] = j;
                  fold[6] = k;
                  site_t const local(latDat->GetContiguousSiteId(LatticeVector(i, j, k)));
                  latDat->SetFOld<D3Q15>(local, fold);
                }
          }
          void tearDown()
          {
            FourCubeBasedTestFixture::tearDown();
          }

          template<class KERNEL>
          void velocityFromLatticeDataTester(bool doforce)
          {
            using hemelb::redblood::details::VelocityFromLatticeData;
            using hemelb::redblood::details::HasForce;

            // Check that type traits works as expected
            CPPUNIT_ASSERT(HasForce<KERNEL>::value == doforce);

            KERNEL kernel(initParams);
            VelocityFromLatticeData<KERNEL> velocityFunctor(*latDat);
            size_t const N(latDat->GetMidDomainCollisionCount(0));

            CPPUNIT_ASSERT(N > 0);

            for (size_t index(0); index < N; ++index)
            {
              geometry::Site<geometry::LatticeData const> const site(index, *latDat);

              // Value to test
              LatticeVelocity const actual = velocityFunctor(index);

              // Check against directly computed values:
              // We know how the lattice was setup
              LatticeVector const position(site.GetGlobalSiteCoords());
              LatticeVelocity const momentum = LatticeVelocity(D3Q15::CX[2],
                                                               D3Q15::CY[2],
                                                               D3Q15::CZ[2])
                  * Dimensionless(position[0])
                  + LatticeVelocity(D3Q15::CX[4], D3Q15::CY[4], D3Q15::CZ[4])
                      * Dimensionless(position[1])
                  + LatticeVelocity(D3Q15::CX[6], D3Q15::CY[6], D3Q15::CZ[6])
                      * Dimensionless(position[2]) + (doforce ?
                    site.GetForce() * 0.5 :
                    LatticeVelocity(0, 0, 0));
              LatticeDensity const density(position[0] + position[1] + position[2]);

              CPPUNIT_ASSERT(helpers::is_zero(actual - momentum / density));

              // Check against Kernel + HydroVars implementation
              typename KERNEL::KHydroVars hydroVars(site);
              kernel.DoCalculateDensityMomentumFeq(hydroVars, site.GetIndex());

              CPPUNIT_ASSERT(helpers::is_zero(actual - hydroVars.velocity));
            }
          }

          void testVelocityDataFromLatticeWithoutForces()
          {
            helpers::LatticeDataAccess(latDat).ZeroOutForces();
            velocityFromLatticeDataTester<LBGK>(false);
          }
          void testVelocityDataFromLatticeWithForces()
          {
            size_t const N(latDat->GetMidDomainCollisionCount(0));

            for (size_t i(0); i < N; ++i)
              latDat->GetSite(i).SetForce(LatticeForceVector(i, 2 * i, double(i * i) * 0.0001));

            velocityFromLatticeDataTester<GuoForcingLBGK>(true);
          }
      };

      CPPUNIT_TEST_SUITE_REGISTRATION (InterpolationTests);
    }
  }
}

#endif  // ONCE
