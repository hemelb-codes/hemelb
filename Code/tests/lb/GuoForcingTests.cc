// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <random>
#include <numeric>
#include <iostream>

#include <catch2/catch.hpp>

#include "lb/lattices/D3Q15.h"
#include "lb/lattices/Lattice.h"
#include "lb/streamers/Streamers.h"
#include "lb/kernels/GuoForcingLBGK.h"
#include "tests/helpers/FourCubeBasedTestFixture.h"
#include "unittests/helpers/Comparisons.h"
#include "unittests/helpers/LatticeDataAccess.h"

namespace hemelb
{
  namespace tests
  {
    // Temporary alias for test migration
    namespace helpers
    {
      using namespace ::hemelb::unittests::helpers;
    }

    TEST_CASE_METHOD(helpers::FourCubeBasedTestFixture<>, "GuoForcingTests") {

        typedef lb::lattices::D3Q19 LatticeType;
        typedef lb::kernels::GuoForcingLBGK<LatticeType> Kernel;
        typedef Kernel::KHydroVars HydroVars;
        typedef lb::lattices::Lattice<LatticeType> Lattice;
	
	// Set forces equal to site index
	size_t const nFluidSites(latDat->GetLocalFluidSiteCount());

        auto GetMeAForce = [&](size_t site) -> LatticeForceVector {
          return LatticeForceVector(LatticeForce(site),
                                    LatticeForce(site) + 0.1,
                                    LatticeForce(site) + 0.3);
        };

        // Sets up _force distribution from actual 3d vector
        auto ForceDistribution = [&](LatticeVelocity const &_velocity, const distribn_t _tau,
				     const LatticeForceVector &_force, distribn_t Fi[]) {
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
        };

        auto FPostCollision = [&](distribn_t const *_f, LatticeForceVector const &_force,
                            distribn_t _fout[]) {
          HydroVars hydroVars(_f, _force);
          Kernel kernel(initParams);
          // Index is never used ... so give a fake value
          kernel.DoCalculateDensityMomentumFeq(hydroVars, 999999);
          kernel.DoCollide(lbmParams, hydroVars);

          std::copy(hydroVars.GetFPostCollision().f,
                    hydroVars.GetFPostCollision().f + LatticeType::NUMVECTORS,
                    _fout);
        };
	auto ForceDistributionTestInstance = [&](LatticeVelocity const &_velocity, const distribn_t _tau,
						 const LatticeForceVector &_force) {
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
            REQUIRE(actual_Fi[i] == Approx(expected_Fi[i]).margin(1e-8));
        };

	for (size_t site(0); site < nFluidSites; ++site)
	  latDat->GetSite(site).SetForce(GetMeAForce(site));
	auto propertyCache = std::make_unique<lb::MacroscopicPropertyCache>(*simState, *latDat);

        SECTION("testHydroVarsGetForce",
		"Tests the correct force is assigned to HydroVars, and that the right specialization is used as well.") {
          size_t const nFluidSites(latDat->GetLocalFluidSiteCount());
          for (size_t i(1); i < nFluidSites; i <<= 1)
            REQUIRE(helpers::is_zero(HydroVars(latDat->GetSite(i)).force - GetMeAForce(i)));
        }

        SECTION("testCalculatDensityAndMomentum",
		"Checks equation 19 of Guo paper (doi: 10.1103/PhysRevE.65.046308)") {
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
          REQUIRE(density == Approx(expected_density).margin(1e-8));
	  REQUIRE(helpers::is_zero(momentum - expected_momentum));

          // Check when forces are not zero
          Lattice::CalculateDensityAndMomentum(const_f,
                                               force[0],
                                               force[1],
                                               force[2],
                                               density,
                                               momentum[0],
                                               momentum[1],
                                               momentum[2]);
          REQUIRE(density == Approx(expected_density).margin(1e-8));
          REQUIRE(helpers::is_zero(momentum - expected_momentum - force * 0.5));
        }

        SECTION("testCalculateDensityMomentumFEq",
		"Checks Guo Correction is included in FEq and velocity") {
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
	  REQUIRE(noforce_density == Approx(force_density).margin(1e-8));
          REQUIRE(helpers::is_zero(noforce_momentum + force * 0.5 - force_momentum));
          REQUIRE(helpers::is_zero(force_velocity - force_momentum / force_density));

          // Compute feq with forces directly
          Lattice::CalculateFeq(force_density,
                                force_momentum[0],
                                force_momentum[1],
                                force_momentum[2],
                                feq);
          for (size_t i(0); i < LatticeType::NUMVECTORS; ++i)
            REQUIRE(force_feq[i] == Approx(feq[i]).margin(1e-8));
        }

        

        SECTION("testForceDistributionRegression",
		"Equation 20 of Guo paper (doi: 10.1103/PhysRevE.65.046308)") {
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

        SECTION("testForceDistributionBadTau",
		"tau = 0.5 => zero force: look at prefactor") {
          distribn_t Fi[LatticeType::NUMVECTORS];
          Lattice::CalculateForceDistribution(0.5, 1.0, 10.0, 100.0, 1, 1, 1, Fi);
          for (size_t i(0); i < LatticeType::NUMVECTORS; ++i)
            REQUIRE(Fi[i] == Approx(0.).margin(1e-8));
        }

        SECTION("testForceDistributionColinear",
		"velocity and force colinear with CX, CY, CZ:"
		" tests second term in eq 20") {
          distribn_t Fi[LatticeType::NUMVECTORS];
          for (size_t i(0); i < LatticeType::NUMVECTORS; ++i)
          {
            LatticeVelocity const ei(LatticeType::CX[i], LatticeType::CY[i], LatticeType::CZ[i]);
            Lattice::CalculateForceDistribution(0.25, ei[0], ei[1], ei[2], ei[0], ei[1], ei[2], Fi);
            distribn_t ei_norm(ei.Dot(ei));
            distribn_t const expected = (1. - 0.5 / 0.25) * LatticeType::EQMWEIGHTS[i]
                * (ei_norm * ei_norm * 9.0);
            REQUIRE(Fi[i] == Approx(expected).margin(1e-8));
          }
        }
        SECTION("testForceDistributionZeroVelocity",
		"force colinear with CX, CY, CZ, and zero force:"
		" tests first term in eq 20") {
          distribn_t Fi[LatticeType::NUMVECTORS];
          for (size_t i(0); i < LatticeType::NUMVECTORS; ++i)
          {
            LatticeVelocity const ei(LatticeType::CX[i], LatticeType::CY[i], LatticeType::CZ[i]);
            Lattice::CalculateForceDistribution(0.25, 0, 0, 0, ei[0], ei[1], ei[2], Fi);
            distribn_t ei_norm(ei.Dot(ei));
            distribn_t const expected = (1. - 0.5 / 0.25) * LatticeType::EQMWEIGHTS[i] * (ei_norm * 3.0);
	    REQUIRE(Fi[i] == Approx(expected).margin(1e-8));
          }
        }

        SECTION("testNoForceMeansNoForceDistribution") {
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
	      REQUIRE(Approx(0e0).margin(1e-12) == Fi[j]);
            }
          }
        }

	SECTION("testZerothOrderMomentOfForceIsNull")
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
	    REQUIRE(Approx(0e0).margin(1e-8) == moment);
          }
          // Force and Velocity are perpendicular
          for(size_t i(0); i < N; ++i)
          {
            auto const a = rdm(e1), b = rdm(e1);
            Lattice::CalculateForceDistribution(1e0, rdm(e1), a, b, 0, -b, a, Fi);
            auto const moment0 = std::accumulate(Fi, Fi+LatticeType::NUMVECTORS, 0e0);
	    REQUIRE(Approx(0e0).margin(1e-8) == moment0);

            Lattice::CalculateForceDistribution(1e0, 0, a, b, rdm(e1), -b, a, Fi);
            auto const moment1 = std::accumulate(Fi, Fi+LatticeType::NUMVECTORS, 0e0);
	    REQUIRE(Approx(0e0).margin(1e-8) == moment1);
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
	    REQUIRE(Approx(0e0).margin(1e-8) == moment);
          }
        }

        SECTION("testGuoEquation")
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
	  REQUIRE(Approx(0e0).margin(1e-8) == sumFi);
        }

        SECTION("testDoCollide",
		"Checks collision process includes forcing."
		"Assumes that equilibrium and non-equilibrium populations, velocity, and forces are computed elsewhere."
		"This just makes sure the force distribution included correctly.") {
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
	    REQUIRE(hydroVars.GetFPostCollision()[i]
		    == Approx(f[i]
			      + hydroVars.GetFNeq()[i] * lbmParams->GetOmega()
			      + Fi[i]).margin(1e-8));
          }
        }

        SECTION("testSimpleCollideAndStream")
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
            REQUIRE(withForce[i] == Approx(*helpers::GetFNew(latDat, index)).margin(1e-8));
            // And that forces from streaming site were used
	    REQUIRE(not helpers::is_zero(withForce[i] - withoutForce[i]));
          }
        }

        SECTION("testSimpleBounceBack")
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
            REQUIRE(withForce[i] == Approx(actual[LatticeType::INVERSEDIRECTIONS[i]]).margin(1e-8));
            // And that forces from streaming site were used
	    REQUIRE(not helpers::is_zero(withForce[i] - withoutForce[i]));
          }
	  REQUIRE(paranoia);
        }

    }

  }
}

