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

#include <cppunit/TestFixture.h>
#include "lb/lattices/D3Q15.h"
#include "lb/lattices/Lattice.h"
#include "lb/kernels/GuoForcingLBGK.h"
#include "unittests/helpers/FourCubeBasedTestFixture.h"
#include "unittests/helpers/Comparisons.h"
#include <iostream>

namespace hemelb { namespace unittests {

  class GuoForcingTests
        : public helpers::FourCubeBasedTestFixture {
    CPPUNIT_TEST_SUITE(GuoForcingTests);
      CPPUNIT_TEST(testHydroVarsGetForce);
      CPPUNIT_TEST(testCalculatDensityAndMomentum);
      CPPUNIT_TEST(testCalculateDensityMomentumFEq);
      CPPUNIT_TEST(testForceDistributionRegression);
      CPPUNIT_TEST(testForceDistributionBadTau);
      CPPUNIT_TEST(testForceDistributionColinear);
      CPPUNIT_TEST(testForceDistributionZeroVelocity);
      CPPUNIT_TEST(testDoCollide);
    CPPUNIT_TEST_SUITE_END();
    typedef lb::lattices::D3Q15 D3Q15;
    typedef lb::kernels::GuoForcingLBGK<lb::lattices::D3Q15> Kernel;
    typedef Kernel::KHydroVars HydroVars;
    typedef lb::lattices::Lattice<D3Q15> Lattice;
    public:
      void setUp() {
        FourCubeBasedTestFixture::setUp();
        // Set forces equal to site index
        size_t const nFluidSites(latDat->GetLocalFluidSiteCount());
        for(size_t site(0); site < nFluidSites; ++site)
          latDat->GetSite(site).SetForce( GetMeAForce(site) );
      }

      // Tests the correct force is assigned to HydroVars, and that the right
      // specialization is used as well.
      void testHydroVarsGetForce() {
        size_t const nFluidSites(latDat->GetLocalFluidSiteCount());
        for(size_t i(1); i < nFluidSites; i <<= 1)
          CPPUNIT_ASSERT(helpers::is_zero(
            HydroVars(latDat->GetSite(i)).force - GetMeAForce(i)
          ));
      }

      // Checks equation 19 of Guo paper (doi: 10.1103/PhysRevE.65.046308)
      void testCalculatDensityAndMomentum() {
        // Create fake density and compute resulting momentum and velocity
        // velocity is not yet corrected for forces
        distribn_t f[D3Q15::NUMVECTORS];
        PhysicalVelocity expected_momentum(0, 0, 0);
        PhysicalDensity expected_density(0);
        for(size_t i(0); i < lb::lattices::D3Q15::NUMVECTORS; ++i) {
          f[i] = distribn_t(std::rand()) / distribn_t(RAND_MAX);
          expected_density += f[i];
          expected_momentum +=
            PhysicalVelocity(D3Q15::CX[i], D3Q15::CY[i], D3Q15::CZ[i]) * f[i];
        }
        // ensures that functions take const value
        distribn_t const * const const_f = f;

        PhysicalVelocity const force(1, 2, 3);
        PhysicalVelocity momentum(0, 0, 0);
        PhysicalDensity density(0);

        // Check when forces are zero
        Lattice :: CalculateDensityAndMomentum(const_f, 0, 0, 0, density,
            momentum[0], momentum[1], momentum[2]);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(density, expected_density, 1e-8);
        CPPUNIT_ASSERT(helpers::is_zero(momentum - expected_momentum));

        // Check when forces are not zero
        Lattice :: CalculateDensityAndMomentum(const_f, force[0], force[1],
            force[2], density, momentum[0], momentum[1], momentum[2]);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(density, expected_density, 1e-8);
        CPPUNIT_ASSERT(helpers::is_zero(
              momentum - expected_momentum - force * 0.5));
      }

      // Checks Guo Correction is included in FEq and velocity
      void testCalculateDensityMomentumFEq() {
        LatticeForceVector const force(1, 2, 3);
        distribn_t f[D3Q15::NUMVECTORS];
        for(size_t i(0); i < lb::lattices::D3Q15::NUMVECTORS; ++i)
          f[i] = distribn_t(std::rand()) / distribn_t(RAND_MAX);
        // ensures that functions take const value
        distribn_t const * const const_f = f;

        PhysicalVelocity noforce_momentum(0, 0, 0), force_momentum(0, 0, 0);
        PhysicalVelocity noforce_velocity(0, 0, 0), force_velocity(0, 0, 0);
        PhysicalDensity noforce_density(0), force_density(0);
        distribn_t noforce_feq[D3Q15::NUMVECTORS],
          force_feq[D3Q15::NUMVECTORS], feq[D3Q15::NUMVECTORS];

        // Compute without force
        Lattice :: CalculateDensityMomentumFEq(
            const_f,
            noforce_density,
            noforce_momentum[0], noforce_momentum[1], noforce_momentum[2],
            noforce_velocity[0], noforce_velocity[1], noforce_velocity[2],
            noforce_feq
        );
        // Compute with forces
        Lattice :: CalculateDensityMomentumFEq(
            const_f, force[0], force[1], force[2],
            force_density,
            force_momentum[0], force_momentum[1], force_momentum[2],
            force_velocity[0], force_velocity[1], force_velocity[2],
            force_feq
        );
        CPPUNIT_ASSERT_DOUBLES_EQUAL(noforce_density, force_density, 1e-8);
        CPPUNIT_ASSERT(helpers::is_zero(
          noforce_momentum + force * 0.5 - force_momentum
        ));
        CPPUNIT_ASSERT(helpers::is_zero(
          force_velocity - force_momentum / force_density
        ));

        // Compute feq with forces directly
        Lattice :: CalculateFeq(
            force_density,
            force_momentum[0], force_momentum[1], force_momentum[2],
            feq
        );
        for(size_t i(0); i < D3Q15::NUMVECTORS; ++i)
          CPPUNIT_ASSERT_DOUBLES_EQUAL(force_feq[i], feq[i], 1e-8);
      }

      void ForceDistribution(
          LatticeVelocity const &_velocity, const distribn_t _tau,
          const LatticeForceVector &_force,
          distribn_t Fi[]
      ) {
        const distribn_t inv_cs2(1e0 / 3e0);
        const distribn_t inv_cs4(1e0 / 9e0);
        const distribn_t prefactor(1. - 0.5 / _tau);
        LatticeVelocity result(0, 0, 0);
        for(size_t i(0); i < D3Q15::NUMVECTORS; ++i) {
          LatticeVelocity const ei(D3Q15::CX[i], D3Q15::CY[i], D3Q15::CZ[i]);
          LatticeForceVector const forcing(
            (ei - _velocity) * inv_cs2 + ei * (ei.Dot(_velocity) * inv_cs4)
          );
          Fi[i] = forcing.Dot(_force) * D3Q15::EQMWEIGHTS[i] * prefactor;
        }
      }
      void ForceDistributionTest(
            LatticeVelocity const &_velocity, const distribn_t _tau,
            const LatticeForceVector &_force
      ) {
        distribn_t expected_Fi[D3Q15::NUMVECTORS];
        distribn_t actual_Fi[D3Q15::NUMVECTORS];

        ForceDistribution(_velocity, _tau, _force, expected_Fi);
        Lattice :: CalculateForceDistribution(
            _tau, _velocity[0], _velocity[1], _velocity[2],
            _force[0], _force[1], _force[2], actual_Fi
        );

        for(size_t i(0); i < D3Q15::NUMVECTORS; ++i)
          CPPUNIT_ASSERT_DOUBLES_EQUAL(actual_Fi[i], expected_Fi[i], 1e-8);
      }

      // Equation 20 of Guo paper (doi: 10.1103/PhysRevE.65.046308)
      void testForceDistributionRegression() {
        ForceDistributionTest(
            LatticeVelocity(1, 0, 0), 0.25, LatticeForceVector(1.5, 0, 0));
        ForceDistributionTest(
            LatticeVelocity(0, 1, 0), 0.25, LatticeForceVector(1.5, 0, 0));
        ForceDistributionTest(
            LatticeVelocity(0, 0, 1), 0.25, LatticeForceVector(3.0, 0, 0));
      }
      // tau = 0.5 => zero force: look at prefactor
      void testForceDistributionBadTau() {
        distribn_t Fi[D3Q15::NUMVECTORS];
        Lattice::CalculateForceDistribution(
            0.5, 1.0, 10.0, 100.0, 1, 1, 1, Fi);
        for(size_t i(0); i < D3Q15::NUMVECTORS; ++i)
          CPPUNIT_ASSERT_DOUBLES_EQUAL(Fi[i], 0, 1e-8);
      }

      // velocity and force colinear with CX, CY, CZ:
      // tests second term in eq 20
      void testForceDistributionColinear() {
        distribn_t Fi[D3Q15::NUMVECTORS];
        for(size_t i(0); i < D3Q15::NUMVECTORS; ++i) {
          LatticeVelocity const ei(D3Q15::CX[i], D3Q15::CY[i], D3Q15::CZ[i]);
          Lattice::CalculateForceDistribution(
              0.25, ei[0], ei[1], ei[2], ei[0], ei[1], ei[2], Fi);
          distribn_t ei_norm(ei.Dot(ei));
          distribn_t const expected
            = (1. - 0.5/0.25) * D3Q15::EQMWEIGHTS[i] * (ei_norm * ei_norm/9.0);
          CPPUNIT_ASSERT_DOUBLES_EQUAL(Fi[i], expected, 1e-8);
        }
      }
      // force colinear with CX, CY, CZ, and zero force:
      // tests first term in eq 20
      void testForceDistributionZeroVelocity() {
        distribn_t Fi[D3Q15::NUMVECTORS];
        for(size_t i(0); i < D3Q15::NUMVECTORS; ++i) {
          LatticeVelocity const ei(D3Q15::CX[i], D3Q15::CY[i], D3Q15::CZ[i]);
          Lattice::CalculateForceDistribution(
              0.25, 0, 0, 0, ei[0], ei[1], ei[2], Fi);
          distribn_t ei_norm(ei.Dot(ei));
          distribn_t const expected
            = (1. - 0.5/0.25) * D3Q15::EQMWEIGHTS[i] * (ei_norm/3.0);
          CPPUNIT_ASSERT_DOUBLES_EQUAL(Fi[i], expected, 1e-8);
        }
      }

      // Checks collision process includes forcing
      // Assumes that equilibrium and non-equilibrium populations, velocity,
      // and forces are computed elsewhere.
      // This just makes sure the force distribution included correctly.
      void testDoCollide() {
        distribn_t f[D3Q15::NUMVECTORS], Fi[D3Q15::NUMVECTORS];
        LatticeForceVector const force(1, 4, 8);
        HydroVars hydroVars(f, force);
        for(size_t i(0); i < lb::lattices::D3Q15::NUMVECTORS; ++i) {
          f[i] = distribn_t(std::rand()) / distribn_t(RAND_MAX);
          hydroVars.SetFNeq(i, distribn_t(std::rand()) / distribn_t(RAND_MAX));
        }

        hydroVars.velocity = LatticeVelocity(1, 2, 3);
        ForceDistribution(hydroVars.velocity, lbmParams->GetTau(),
            hydroVars.force, Fi);

        Kernel kernel(initParams);
        kernel.DoCollide(lbmParams, hydroVars);
        for(size_t i(0); i < D3Q15::NUMVECTORS; ++i) {
          CPPUNIT_ASSERT_DOUBLES_EQUAL(
              hydroVars.GetFPostCollision()[i],
              f[i] + hydroVars.GetFNeq()[i] * lbmParams->GetOmega() + Fi[i],
              1e-8
          );
        }
      }

    protected:
      LatticeForceVector GetMeAForce(size_t site) {
        return LatticeForceVector(
            LatticeForce(site),
            LatticeForce(site) + 0.1,
            LatticeForce(site) + 0.3
        );
      }
  };

}}
  CPPUNIT_TEST_SUITE_REGISTRATION(hemelb::unittests::GuoForcingTests);

#endif
