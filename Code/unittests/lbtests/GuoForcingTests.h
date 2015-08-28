//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELB_UNITTESTS_GUOFORCING_H
#define HEMELB_UNITTESTS_GUOFORCING_H

#include <random>
#include <numeric>
#include <cppunit/TestFixture.h>
#include "lb/lattices/D3Q15.h"
#include "lb/lattices/Lattice.h"
#include "lb/streamers/Streamers.h"
#include "lb/kernels/GuoForcingLBGK.h"
#include "unittests/helpers/FourCubeBasedTestFixture.h"
#include "unittests/helpers/Comparisons.h"
#include "unittests/helpers/LatticeDataAccess.h"
#include <iostream>

namespace hemelb
{
  namespace unittests
  {

    class GuoForcingTests : public helpers::FourCubeBasedTestFixture
    {

        CPPUNIT_TEST_SUITE (GuoForcingTests);
        CPPUNIT_TEST (testHydroVarsGetForce);
        CPPUNIT_TEST (testCalculatDensityAndMomentum);
        CPPUNIT_TEST (testCalculateDensityMomentumFEq);
        CPPUNIT_TEST (testForceDistributionRegression);
        CPPUNIT_TEST (testForceDistributionBadTau);
        CPPUNIT_TEST (testForceDistributionColinear);
        CPPUNIT_TEST (testForceDistributionZeroVelocity);
        CPPUNIT_TEST (testZerothOrderMomentOfForceIsNull);
        CPPUNIT_TEST (testDoCollide);
        CPPUNIT_TEST (testNoForceMeansNoForceDistribution);
        CPPUNIT_TEST (testGuoEquation);
        CPPUNIT_TEST (testSimpleCollideAndStream);
        CPPUNIT_TEST (testSimpleBounceBack);CPPUNIT_TEST_SUITE_END();

        typedef lb::lattices::D3Q19 LatticeType;
        typedef lb::kernels::GuoForcingLBGK<LatticeType> Kernel;
        typedef Kernel::KHydroVars HydroVars;
        typedef lb::lattices::Lattice<LatticeType> Lattice;
      public:
        void setUp()
        {
          FourCubeBasedTestFixture::setUp();
          // Set forces equal to site index
          size_t const nFluidSites(latDat->GetLocalFluidSiteCount());
          for (size_t site(0); site < nFluidSites; ++site)
            latDat->GetSite(site).SetForce(GetMeAForce(site));
          propertyCache = new lb::MacroscopicPropertyCache(*simState, *latDat);
        }

        void tearDown()
        {
          delete propertyCache;
          FourCubeBasedTestFixture::tearDown();
        }

        // Tests the correct force is assigned to HydroVars, and that the right
        // specialization is used as well.
        void testHydroVarsGetForce()
        {
          size_t const nFluidSites(latDat->GetLocalFluidSiteCount());
          for (size_t i(1); i < nFluidSites; i <<= 1)
            CPPUNIT_ASSERT(helpers::is_zero(HydroVars(latDat->GetSite(i)).force - GetMeAForce(i)));
        }

        // Checks equation 19 of Guo paper (doi: 10.1103/PhysRevE.65.046308)
        void testCalculatDensityAndMomentum()
        {
          // Create fake density and compute resulting momentum and velocity
          // velocity is not yet corrected for forces
          distribn_t f[LatticeType::NUMVECTORS];
          PhysicalVelocity expected_momentum(0, 0, 0);
          PhysicalDensity expected_density(0);
          for (size_t i(0); i < LatticeType::NUMVECTORS; ++i)
          {
            f[i] = distribn_t(std::rand()) / distribn_t(RAND_MAX);
            expected_density += f[i];
            expected_momentum += PhysicalVelocity(
                LatticeType::CX[i], LatticeType::CY[i], LatticeType::CZ[i]) * f[i];
          }
          // ensures that functions take const value
          distribn_t const * const const_f = f;

          PhysicalVelocity const force(1, 2, 3);
          PhysicalVelocity momentum(0, 0, 0);
          PhysicalDensity density(0);

          // Check when forces are zero
          Lattice::CalculateDensityAndMomentum(const_f,
                                               0,
                                               0,
                                               0,
                                               density,
                                               momentum[0],
                                               momentum[1],
                                               momentum[2]);
          CPPUNIT_ASSERT_DOUBLES_EQUAL(density, expected_density, 1e-8);
          CPPUNIT_ASSERT(helpers::is_zero(momentum - expected_momentum));

          // Check when forces are not zero
          Lattice::CalculateDensityAndMomentum(const_f,
                                               force[0],
                                               force[1],
                                               force[2],
                                               density,
                                               momentum[0],
                                               momentum[1],
                                               momentum[2]);
          CPPUNIT_ASSERT_DOUBLES_EQUAL(density, expected_density, 1e-8);
          CPPUNIT_ASSERT(helpers::is_zero(momentum - expected_momentum - force * 0.5));
        }

        // Checks Guo Correction is included in FEq and velocity
        void testCalculateDensityMomentumFEq()
        {
          LatticeForceVector const force(1, 2, 3);
          distribn_t f[LatticeType::NUMVECTORS];
          for (size_t i(0); i < LatticeType::NUMVECTORS; ++i)
            f[i] = distribn_t(std::rand()) / distribn_t(RAND_MAX);
          // ensures that functions take const value
          distribn_t const * const const_f = f;

          PhysicalVelocity noforce_momentum(0, 0, 0), force_momentum(0, 0, 0);
          PhysicalVelocity noforce_velocity(0, 0, 0), force_velocity(0, 0, 0);
          PhysicalDensity noforce_density(0), force_density(0);
          distribn_t noforce_feq[LatticeType::NUMVECTORS], force_feq[LatticeType::NUMVECTORS],
              feq[LatticeType::NUMVECTORS];

          // Compute without force
          Lattice::CalculateDensityMomentumFEq(const_f,
                                               noforce_density,
                                               noforce_momentum[0],
                                               noforce_momentum[1],
                                               noforce_momentum[2],
                                               noforce_velocity[0],
                                               noforce_velocity[1],
                                               noforce_velocity[2],
                                               noforce_feq);
          // Compute with forces
          Lattice::CalculateDensityMomentumFEq(const_f,
                                               force[0],
                                               force[1],
                                               force[2],
                                               force_density,
                                               force_momentum[0],
                                               force_momentum[1],
                                               force_momentum[2],
                                               force_velocity[0],
                                               force_velocity[1],
                                               force_velocity[2],
                                               force_feq);
          CPPUNIT_ASSERT_DOUBLES_EQUAL(noforce_density, force_density, 1e-8);
          CPPUNIT_ASSERT(helpers::is_zero(noforce_momentum + force * 0.5 - force_momentum));
          CPPUNIT_ASSERT(helpers::is_zero(force_velocity - force_momentum / force_density));

          // Compute feq with forces directly
          Lattice::CalculateFeq(force_density,
                                force_momentum[0],
                                force_momentum[1],
                                force_momentum[2],
                                feq);
          for (size_t i(0); i < LatticeType::NUMVECTORS; ++i)
            CPPUNIT_ASSERT_DOUBLES_EQUAL(force_feq[i], feq[i], 1e-8);
        }

        void ForceDistributionTestInstance(LatticeVelocity const &_velocity, const distribn_t _tau,
                                           const LatticeForceVector &_force)
        {
          distribn_t expected_Fi[LatticeType::NUMVECTORS];
          distribn_t actual_Fi[LatticeType::NUMVECTORS];

          ForceDistribution(_velocity, _tau, _force, expected_Fi);
          Lattice::CalculateForceDistribution(_tau,
                                              _velocity[0],
                                              _velocity[1],
                                              _velocity[2],
                                              _force[0],
                                              _force[1],
                                              _force[2],
                                              actual_Fi);

          for (size_t i(0); i < LatticeType::NUMVECTORS; ++i)
            CPPUNIT_ASSERT_DOUBLES_EQUAL(actual_Fi[i], expected_Fi[i], 1e-8);
        }

        // Equation 20 of Guo paper (doi: 10.1103/PhysRevE.65.046308)
        void testForceDistributionRegression()
        {
          ForceDistributionTestInstance(LatticeVelocity(1, 0, 0),
                                        0.25,
                                        LatticeForceVector(1.5, 0, 0));
          ForceDistributionTestInstance(LatticeVelocity(0, 1, 0),
                                        0.25,
                                        LatticeForceVector(1.5, 0, 0));
          ForceDistributionTestInstance(LatticeVelocity(0, 0, 1),
                                        0.25,
                                        LatticeForceVector(3.0, 0, 0));
        }
        // tau = 0.5 => zero force: look at prefactor
        void testForceDistributionBadTau()
        {
          distribn_t Fi[LatticeType::NUMVECTORS];
          Lattice::CalculateForceDistribution(0.5, 1.0, 10.0, 100.0, 1, 1, 1, Fi);
          for (size_t i(0); i < LatticeType::NUMVECTORS; ++i)
            CPPUNIT_ASSERT_DOUBLES_EQUAL(Fi[i], 0, 1e-8);
        }

        // velocity and force colinear with CX, CY, CZ:
        // tests second term in eq 20
        void testForceDistributionColinear()
        {
          distribn_t Fi[LatticeType::NUMVECTORS];
          for (size_t i(0); i < LatticeType::NUMVECTORS; ++i)
          {
            LatticeVelocity const ei(LatticeType::CX[i], LatticeType::CY[i], LatticeType::CZ[i]);
            Lattice::CalculateForceDistribution(0.25, ei[0], ei[1], ei[2], ei[0], ei[1], ei[2], Fi);
            distribn_t ei_norm(ei.Dot(ei));
            distribn_t const expected = (1. - 0.5 / 0.25) * LatticeType::EQMWEIGHTS[i]
                * (ei_norm * ei_norm * 9.0);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(Fi[i], expected, 1e-8);
          }
        }
        // force colinear with CX, CY, CZ, and zero force:
        // tests first term in eq 20
        void testForceDistributionZeroVelocity()
        {
          distribn_t Fi[LatticeType::NUMVECTORS];
          for (size_t i(0); i < LatticeType::NUMVECTORS; ++i)
          {
            LatticeVelocity const ei(LatticeType::CX[i], LatticeType::CY[i], LatticeType::CZ[i]);
            Lattice::CalculateForceDistribution(0.25, 0, 0, 0, ei[0], ei[1], ei[2], Fi);
            distribn_t ei_norm(ei.Dot(ei));
            distribn_t const expected = (1. - 0.5 / 0.25) * LatticeType::EQMWEIGHTS[i] * (ei_norm * 3.0);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(Fi[i], expected, 1e-8);
          }
        }

        void testNoForceMeansNoForceDistribution()
        {
          auto const N = 10;
          std::random_device rd;
          std::default_random_engine e1(rd());
          std::uniform_real_distribution<distribn_t> rdm(-10, 10);

          distribn_t Fi[LatticeType::NUMVECTORS];
          for(size_t i(0); i < N; ++i)
          {
            Lattice::CalculateForceDistribution(1e0, rdm(e1), rdm(e1), rdm(e1), 0e0, 0e0, 0e0, Fi);
            for(size_t j(0); j < LatticeType::NUMVECTORS; ++j)
            {
              CPPUNIT_ASSERT_DOUBLES_EQUAL(0e0, Fi[j], 1e-12);
            }
          }
        }

        void testZerothOrderMomentOfForceIsNull()
        {
          auto const N = 10;
          std::random_device rd;
          std::default_random_engine e1(rd());
          std::uniform_real_distribution<distribn_t> rdm(-10, 10);

          distribn_t Fi[LatticeType::NUMVECTORS];
          // No velocity
          for(size_t i(0); i < N; ++i)
          {
            Lattice::CalculateForceDistribution(
                1e0, 0e0, 0e0, 0e0, rdm(e1), rdm(e1), rdm(e1), Fi);
            auto const moment = std::accumulate(Fi, Fi+LatticeType::NUMVECTORS, 0e0);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(0e0, moment, 1e-8);
          }
          // Force and Velocity are perpendicular
          for(size_t i(0); i < N; ++i)
          {
            auto const a = rdm(e1), b = rdm(e1);
            Lattice::CalculateForceDistribution(1e0, rdm(e1), a, b, 0, -b, a, Fi);
            auto const moment0 = std::accumulate(Fi, Fi+LatticeType::NUMVECTORS, 0e0);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(0e0, moment0, 1e-8);

            Lattice::CalculateForceDistribution(1e0, 0, a, b, rdm(e1), -b, a, Fi);
            auto const moment1 = std::accumulate(Fi, Fi+LatticeType::NUMVECTORS, 0e0);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(0e0, moment1, 1e-8);
          }

          // Force and velocity tied via pops
          // populations
          for(size_t i(0); i < 1; ++i)
          {
            LatticeForceVector const F(rdm(e1), rdm(e1), rdm(e1));
            LatticeVelocity v = F * 0.5;
            LatticeDensity rho(0);
            for(size_t j(0); j < LatticeType::NUMVECTORS; ++j)
            {
              auto const fi = rdm(e1);
              rho += fi;
              v += LatticeVelocity(LatticeType::CX[j], LatticeType::CY[j], LatticeType::CZ[j]) * fi;
            }
            v = v / rho;
            Lattice::CalculateForceDistribution(1e0, v.x, v.y, v.z, F.x, F.y, F.z, Fi);
            auto const moment = std::accumulate(Fi, Fi+LatticeType::NUMVECTORS, 0e0);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(0e0, moment, 1e-8);
          }
        }

        void testGuoEquation()
        {
          std::random_device rd;
          std::default_random_engine e1(rd());
          std::uniform_real_distribution<distribn_t> rdm(-10, 10);

          typedef lb::lattices::D3Q19 Lattice;
          LatticeVelocity const v(rdm(e1), rdm(e1), rdm(e1));
          LatticeForceVector const F(rdm(e1), rdm(e1), rdm(e1));
          LatticeForce sumFi(0);
          for(size_t j(0); j < Lattice::NUMVECTORS; ++j)
          {
            LatticeVelocity const ei(Lattice::CXD[j], Lattice::CYD[j], Lattice::CZD[j]);
            sumFi += Lattice::EQMWEIGHTS[j] * ((ei - v) * 3e0 + ei * ei.Dot(v) * 9e0).Dot(F);
          }
          CPPUNIT_ASSERT_DOUBLES_EQUAL(0e0, sumFi, 1e-8);
        }

        // Checks collision process includes forcing
        // Assumes that equilibrium and non-equilibrium populations, velocity,
        // and forces are computed elsewhere.
        // This just makes sure the force distribution included correctly.
        void testDoCollide()
        {
          distribn_t f[LatticeType::NUMVECTORS], Fi[LatticeType::NUMVECTORS];
          LatticeForceVector const force(1, 4, 8);
          HydroVars hydroVars(f, force);
          for (size_t i(0); i < LatticeType::NUMVECTORS; ++i)
          {
            f[i] = distribn_t(std::rand()) / distribn_t(RAND_MAX);
            hydroVars.SetFNeq(i, distribn_t(std::rand()) / distribn_t(RAND_MAX));
          }

          hydroVars.velocity = LatticeVelocity(1, 2, 3);
          ForceDistribution(hydroVars.velocity, lbmParams->GetTau(), hydroVars.force, Fi);

          Kernel kernel(initParams);
          kernel.DoCollide(lbmParams, hydroVars);
          for (size_t i(0); i < LatticeType::NUMVECTORS; ++i)
          {
            CPPUNIT_ASSERT_DOUBLES_EQUAL(hydroVars.GetFPostCollision()[i],
                                         f[i] + hydroVars.GetFNeq()[i] * lbmParams->GetOmega()
                                             + Fi[i],
                                         1e-8);
          }
        }

        void testSimpleCollideAndStream()
        {
          // Lattice with a single non-zero distribution in the middle
          // Same for forces
          LatticeVector const position(2, 2, 2);
          geometry::Site<geometry::LatticeData> const site(latDat->GetSite(position));
          helpers::allZeroButOne<LatticeType>(latDat, position);

          // Get collided distributions at this site
          // Will compare zero and non-zero forces, to make sure they are different
          // Assumes collides works, since tested in TestDoCollide
          distribn_t withForce[LatticeType::NUMVECTORS], withoutForce[LatticeType::NUMVECTORS];
          FPostCollision(site.GetFOld<LatticeType>(), site.GetForce(), withForce);
          FPostCollision(site.GetFOld<LatticeType>(), LatticeForceVector(0, 0, 0), withoutForce);

          // Stream that site
          using lb::streamers::SimpleCollideAndStream;
          using lb::collisions::Normal;
          SimpleCollideAndStream<Normal<Kernel> > streamer(initParams);
          streamer.StreamAndCollide<false>(site.GetIndex(), 1, lbmParams, latDat, *propertyCache);

          // Now check streaming worked correctly
          for (size_t i(0); i < LatticeType::NUMVECTORS; ++i)
          {
            auto const index = site.GetStreamedIndex<LatticeType>(i);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(withForce[i], *helpers::GetFNew(latDat, index), 1e-8);
            // And that forces from streaming site were used
            CPPUNIT_ASSERT(not helpers::is_zero(withForce[i] - withoutForce[i]));
          }
        }

        void testSimpleBounceBack()
        {
          // Lattice with a single non-zero distribution in the middle
          // Same for forces
          LatticeVector const position(1, 1, 1);
          geometry::Site<geometry::LatticeData> const site(latDat->GetSite(position));
          helpers::allZeroButOne<LatticeType>(latDat, position);

          // Get collided distributions at this site
          // Will compare zero and non-zero forces, to make sure they are different
          // Assumes collides works, since tested in TestDoCollide
          distribn_t withForce[LatticeType::NUMVECTORS], withoutForce[LatticeType::NUMVECTORS];
          FPostCollision(site.GetFOld<LatticeType>(), site.GetForce(), withForce);
          FPostCollision(site.GetFOld<LatticeType>(), LatticeForceVector(0, 0, 0), withoutForce);

          // Stream that site
          using lb::streamers::SimpleBounceBack;
          using lb::collisions::Normal;
          SimpleBounceBack<Normal<Kernel> >::Type streamer(initParams);
          streamer.StreamAndCollide<false>(site.GetIndex(), 1, lbmParams, latDat, *propertyCache);

          distribn_t const * const actual = helpers::GetFNew<LatticeType>(latDat, position);
          bool paranoia(false);
          for (size_t i(0); i < LatticeType::NUMVECTORS; ++i)
          {
            if (not site.HasWall(i))
              continue;
            paranoia = true;
            CPPUNIT_ASSERT_DOUBLES_EQUAL(withForce[i], actual[LatticeType::INVERSEDIRECTIONS[i]], 1e-8);
            // And that forces from streaming site were used
            CPPUNIT_ASSERT(not helpers::is_zero(withForce[i] - withoutForce[i]));
          }
          CPPUNIT_ASSERT(paranoia);
        }

      protected:
        LatticeForceVector GetMeAForce(size_t site)
        {
          return LatticeForceVector(LatticeForce(site),
                                    LatticeForce(site) + 0.1,
                                    LatticeForce(site) + 0.3);
        }

        // Sets up _force distribution from actual 3d vector
        void ForceDistribution(LatticeVelocity const &_velocity, const distribn_t _tau,
                               const LatticeForceVector &_force, distribn_t Fi[])
        {
          const distribn_t inv_cs2(3e0);
          const distribn_t inv_cs4(9e0);
          const distribn_t prefactor(1. - 0.5 / _tau);
          LatticeVelocity result(0, 0, 0);
          for (size_t i(0); i < LatticeType::NUMVECTORS; ++i)
          {
            LatticeVelocity const ei(LatticeType::CX[i], LatticeType::CY[i], LatticeType::CZ[i]);
            LatticeForceVector const forcing( (ei - _velocity) * inv_cs2
                + ei * (ei.Dot(_velocity) * inv_cs4));
            Fi[i] = forcing.Dot(_force) * LatticeType::EQMWEIGHTS[i] * prefactor;
          }
        }

        void FPostCollision(distribn_t const *_f, LatticeForceVector const &_force,
                            distribn_t _fout[])
        {
          HydroVars hydroVars(_f, _force);
          Kernel kernel(initParams);
          // Index is never used ... so give a fake value
          kernel.DoCalculateDensityMomentumFeq(hydroVars, 999999);
          kernel.DoCollide(lbmParams, hydroVars);

          std::copy(hydroVars.GetFPostCollision().f,
                    hydroVars.GetFPostCollision().f + LatticeType::NUMVECTORS,
                    _fout);
        }

      private:
        lb::MacroscopicPropertyCache* propertyCache;
    };

  }
} // namespace

CPPUNIT_TEST_SUITE_REGISTRATION (hemelb::unittests::GuoForcingTests);

#endif
