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
#include "redblood/Interpolation.h"
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
          CPPUNIT_TEST (testOffLattice<stencil::FourPoint> );
          CPPUNIT_TEST (testOffLattice<stencil::CosineApprox> );
          CPPUNIT_TEST (testOffLattice<stencil::ThreePoint> );
          CPPUNIT_TEST (testOffLattice<stencil::TwoPoint> );
          CPPUNIT_TEST (testOffLatticeZeroOutsideStencil<stencil::FourPoint> );
          CPPUNIT_TEST (testOffLatticeZeroOutsideStencil<stencil::CosineApprox> );
          CPPUNIT_TEST (testOffLatticeZeroOutsideStencil<stencil::ThreePoint> );
          CPPUNIT_TEST (testOffLatticeZeroOutsideStencil<stencil::TwoPoint> );
          CPPUNIT_TEST (testInterpolateLinearFunction<stencil::FourPoint> );
          CPPUNIT_TEST (testInterpolateLinearFunction<stencil::CosineApprox> );
          CPPUNIT_TEST (testInterpolateLinearFunction<stencil::ThreePoint> );
          CPPUNIT_TEST (testInterpolateLinearFunction<stencil::TwoPoint> );
          CPPUNIT_TEST (testInterpolateQuadraticFunction<stencil::FourPoint> );
          CPPUNIT_TEST (testInterpolateQuadraticFunction<stencil::CosineApprox> );
          CPPUNIT_TEST (testInterpolateQuadraticFunction<stencil::ThreePoint> );
          CPPUNIT_TEST (testInterpolateQuadraticFunction<stencil::TwoPoint> );CPPUNIT_TEST_SUITE_END();

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

          std::vector<std::pair<LatticeVector, size_t>> offLatticeData(stencil::FourPoint const &)
          {
            return
            {
              { LatticeVector(55, 51, 14), 0},
              { LatticeVector(55, 51, 15), 1},
              { LatticeVector(55, 51, 16), 1},
              { LatticeVector(55, 51, 17), 1},
              { LatticeVector(55, 52, 14), 1},
              { LatticeVector(56, 51, 14), 12},
              { LatticeVector(58, 54, 17), 47}
            };
          }
          std::vector<std::pair<LatticeVector, size_t>> offLatticeData(
              stencil::CosineApprox const &)
          {
            return offLatticeData(stencil::FourPoint());
          }
          std::vector<std::pair<LatticeVector, size_t>> offLatticeData(stencil::ThreePoint const &)
          {
            return
            {
              { LatticeVector(56, 52, 14), 0},
              { LatticeVector(56, 52, 15), 1},
              { LatticeVector(56, 52, 16), 1},
              { LatticeVector(56, 53, 14), 1},
              { LatticeVector(57, 52, 14), 6},
              { LatticeVector(58, 54, 16), 17}
            };
          }
          std::vector<std::pair<LatticeVector, size_t>> offLatticeData(stencil::TwoPoint const &)
          {
            return
            {
              { LatticeVector(56, 52, 15), 0},
              { LatticeVector(56, 52, 16), 1},
              { LatticeVector(56, 53, 15), 1},
              { LatticeVector(57, 52, 15), 2},
              { LatticeVector(57, 53, 16), 3}
            };
          }
          template<class STENCIL> void testOffLattice()
          {
            LatticePosition const pos(56.51, 52.9, 15.2);
            InterpolationIterator<STENCIL> iterator(pos);

            // Expected point and number of times to increment iterator
            auto const expected = offLatticeData(STENCIL());

            // Checks iteration goes through correct sequence
            for (auto const &item : expected)
            {
              for (size_t j(0); j < item.second; ++j, ++iterator)
                ;

              CPPUNIT_ASSERT(iterator.IsValid());
              CPPUNIT_ASSERT_EQUAL(item.first.x, iterator->x);
              CPPUNIT_ASSERT_EQUAL(item.first.y, iterator->y);
              CPPUNIT_ASSERT_EQUAL(item.first.z, iterator->z);
              CPPUNIT_ASSERT_DOUBLES_EQUAL(iterator.weight(),
                                           STENCIL::stencil(pos - item.first),
                                           1e-8);
            }

            ++iterator;
            CPPUNIT_ASSERT(not iterator.IsValid());
          }

          std::vector<LatticeVector> outsidePoints(stencil::FourPoint const&) const
          {
            return
            {
              LatticeVector(57, 53, 13), LatticeVector(57, 53, 18), LatticeVector(57, 50, 17),
              LatticeVector(57, 55, 17), LatticeVector(54, 53, 17), LatticeVector(59, 53, 17)
            };
          }
          std::vector<LatticeVector> outsidePoints(stencil::CosineApprox const&) const
          {
            return outsidePoints(stencil::FourPoint());
          }
          std::vector<LatticeVector> outsidePoints(stencil::ThreePoint const&) const
          {
            return
            {
              LatticeVector(56, 52, 13), LatticeVector(56, 52, 17), LatticeVector(56, 51, 16),
              LatticeVector(56, 55, 16), LatticeVector(55, 52, 16), LatticeVector(59, 52, 16)
            };
          }
          std::vector<LatticeVector> outsidePoints(stencil::TwoPoint const&) const
          {
            return
            {
              LatticeVector(56, 52, 14), LatticeVector(56, 52, 17), LatticeVector(56, 51, 15),
              LatticeVector(56, 54, 16), LatticeVector(55, 52, 16), LatticeVector(58, 52, 16)
            };
          }
          template<class STENCIL> void testOffLatticeZeroOutsideStencil()
          {
            LatticePosition const pos(56.51, 52.9, 15.2);
            InterpolationIterator<STENCIL> iterator(pos);

            // Checks that outside iteration box, weights are zero
            for (auto const& vec : outsidePoints(STENCIL()))
            {
              LatticeVector const dx(1, 0, 0), dy(0, 1, 0), dz(0, 0, 1);
              CPPUNIT_ASSERT_DOUBLES_EQUAL(0e0, STENCIL::stencil(pos - vec), 1e-8);
              // checks we are one step outside the iteration box only.
              // this is really a test on the zero_vecs data, eg a test of the test.
              size_t const one_non_zero = 0
                  + size_t(not helpers::is_zero(STENCIL::stencil(pos + dx - vec)))
                  + size_t(not helpers::is_zero(STENCIL::stencil(pos - dx - vec)))
                  + size_t(not helpers::is_zero(STENCIL::stencil(pos + dy - vec)))
                  + size_t(not helpers::is_zero(STENCIL::stencil(pos - dy - vec)))
                  + size_t(not helpers::is_zero(STENCIL::stencil(pos + dz - vec)))
                  + size_t(not helpers::is_zero(STENCIL::stencil(pos - dz - vec)));
              CPPUNIT_ASSERT_EQUAL(size_t(1), one_non_zero);
            }
          }

          template<class FUNCTION, class STENCIL>
          void check(Dimensionless x, Dimensionless y, Dimensionless z, Dimensionless tolerance =
                         1e-8)
          {
            FUNCTION func;
            LatticePosition expected(func(x, y, z));
            LatticePosition actual(interpolate<FUNCTION, STENCIL>(func, x, y, z));
            CPPUNIT_ASSERT_DOUBLES_EQUAL(expected.x, actual.x, tolerance);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(expected.y, actual.y, tolerance);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(expected.z, actual.z, tolerance);
          }

          // Test interpolation when the point is on the grid
          template<class STENCIL> void testInterpolateLinearFunction()
          {
            auto const tolerance = std::is_same<STENCIL, stencil::CosineApprox>::value ?
              5e-2 :
              1e-8;
            check<PlanarFunction, STENCIL>(0, 0, 0, tolerance);
            check<PlanarFunction, STENCIL>(0.1, 0.5, 0.6, tolerance);
            check<PlanarFunction, STENCIL>(-5.1, 0.5, 8.7, tolerance);
            check<PlanarFunction, STENCIL>(-5, 0, -1, tolerance);
          }

          template<class STENCIL> void testInterpolateQuadraticFunction()
          {
            QuadraticFunction quad;
            // Error depends on variation on scale larger than stencil
            Dimensionless const tolerance( (quad(0, 0, 0) - quad(12, 12, 12)).GetMagnitude()
                * 1e-2);
            check<QuadraticFunction, STENCIL>(0, 0, 0, tolerance);
            check<QuadraticFunction, STENCIL>(0.1, 0.5, 0.6, tolerance);
            check<QuadraticFunction, STENCIL>(-5.1, 0.5, 8.7, tolerance);
            check<QuadraticFunction, STENCIL>(-5, 0, -1, tolerance);
          }
      };

      class VelocityInterpolationTests : public helpers::FourCubeBasedTestFixture
      {
          typedef lb::lattices::D3Q15 D3Q15;
          typedef lb::kernels::LBGK<D3Q15> LBGK;
          typedef lb::kernels::GuoForcingLBGK<D3Q15> GuoForcingLBGK;
          CPPUNIT_TEST_SUITE (VelocityInterpolationTests);
          CPPUNIT_TEST (testVelocityDataFromLatticeWithForces);
          CPPUNIT_TEST (testVelocityDataFromLatticeWithoutForces);CPPUNIT_TEST_SUITE_END();

        public:
          void setUp()
          {
            FourCubeBasedTestFixture::setUp();

            LatticeVector const min(latDat->GetGlobalSiteMins());
            LatticeVector const max(latDat->GetGlobalSiteMaxes());

            for (std::size_t i(min[0]); i <= std::size_t(max[0]); ++i)
              for (std::size_t j(min[1]); j <= std::size_t(max[1]); ++j)
                for (std::size_t k(min[2]); k <= std::size_t(max[2]); ++k)
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
