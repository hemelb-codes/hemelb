// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_TESTS_LB_LBTESTSHELPER_H
#define HEMELB_TESTS_LB_LBTESTSHELPER_H

#include <cmath>
#include <string>

#include <catch2/catch.hpp>

#include "constants.h"
#include "lb/HFunction.h"
#include "lb/kernels/BaseKernel.h"
#include "lb/MacroscopicPropertyCache.h"

#include "tests/helpers/ApproxVector.h"
#include "tests/helpers/FourCubeLatticeData.h"

namespace hemelb
{
  namespace tests
  {
    namespace LbTestsHelper
    {
      using Vec = util::Vector3D<distribn_t>;

      template<typename Lattice>
      void CalculateRhoMomentum(const distribn_t f[Lattice::NUMVECTORS], distribn_t& rho,  Vec& v)
      {
	rho = 0.0;

	v.x = v.y = v.z = 0.0;
	for (unsigned int ii = 0; ii < Lattice::NUMVECTORS; ++ii)
	  {
	    rho += f[ii];
	    v.x += f[ii] * Lattice::CXD[ii];
	    v.y += f[ii] * Lattice::CYD[ii];
	    v.z += f[ii] * Lattice::CZD[ii];
	  }
      }

      template<typename Lattice>
      Vec CalculateMomentum(const distribn_t f[Lattice::NUMVECTORS])
      {
	auto v = Vec::Zero();
	for (unsigned int ii = 0; ii < Lattice::NUMVECTORS; ++ii) {
	  v += f[ii] * Vec{Lattice::CXD[ii], Lattice::CYD[ii], Lattice::CZD[ii]};
	}
	return v;
      }

      template<typename Lattice>
      void CalculateAnsumaliEntropicEqmF(distribn_t density,
					 const Vec& momentum,
					 distribn_t f[Lattice::NUMVECTORS])
      {
	// Calculate velocity.
	auto velocity = momentum / density;

	distribn_t B[3];
	distribn_t term1[3];
	distribn_t term2[3];
	for (int ii = 0; ii < 3; ++ii)
	  {
	    B[ii] = sqrt(1.0 + 3.0 * velocity[ii] * velocity[ii]);
	    term1[ii] = 2.0 - B[ii];
	    term2[ii] = (2.0 * velocity[ii] + B[ii]) / (1.0 - velocity[ii]);
	  }

	for (unsigned int ii = 0; ii < Lattice::NUMVECTORS; ++ii)
	  {
	    f[ii] = density * Lattice::EQMWEIGHTS[ii];

	    f[ii] *= term1[0] * pow(term2[0], (double) Lattice::CX[ii]);
	    f[ii] *= term1[1] * pow(term2[1], (double) Lattice::CY[ii]);
	    f[ii] *= term1[2] * pow(term2[2], (double) Lattice::CZ[ii]);
	  }
      }

      template<typename Lattice>
      void CalculateLBGKEqmF(distribn_t density,
			     const util::Vector3D<distribn_t>& momentum,
			     distribn_t f[Lattice::NUMVECTORS])
      {
	for (unsigned int ii = 0; ii < Lattice::NUMVECTORS; ++ii)
	  {
	    // Calculate the dot-product of the velocity with the direction vector.
	    distribn_t vSum = momentum.Dot(util::Vector3D<double>{Lattice::CXD[ii], Lattice::CYD[ii], Lattice::CZD[ii]});

	    // Calculate the squared magnitude of the velocity.
	    distribn_t momentum2Sum = momentum.GetMagnitudeSquared();

	    // F eqm = density proportional component...
	    f[ii] = density;

	    // ... - v^2 component...
	    f[ii] -= ( (3.0 / 2.0) * momentum2Sum / density);

	    // ... + v^1 component
	    f[ii] += 3.0 * vSum + (9.0 / 2.0) * vSum * vSum / density;

	    // Multiply by eqm weight.
	    f[ii] *= Lattice::EQMWEIGHTS[ii];
	  }
      }

      template<typename Lattice>
      void CalculateEntropicCollision(const distribn_t f[Lattice::NUMVECTORS],
					     const distribn_t f_eq[Lattice::NUMVECTORS],
					     distribn_t tau,
					     distribn_t beta,
					     distribn_t f_collided[Lattice::NUMVECTORS])
      {
	lb::HFunction<Lattice> HFunc(f, f_eq);

	distribn_t alpha = hemelb::util::NumericalMethods::NewtonRaphson(&HFunc, 2.0, 1.0E-100);

	for (unsigned int ii = 0; ii < Lattice::NUMVECTORS; ++ii)
	  {
	    f_collided[ii] = f[ii] + alpha * beta * (f[ii] - f_eq[ii]);
	  }
      }

      template<typename Lattice>
      void CalculateLBGKCollision(const distribn_t f[Lattice::NUMVECTORS],
					 const distribn_t f_eq[Lattice::NUMVECTORS],
					 distribn_t omega,
					 distribn_t f_collided[Lattice::NUMVECTORS])
      {
	for (unsigned int ii = 0; ii < Lattice::NUMVECTORS; ++ii)
	  {
	    f_collided[ii] = f[ii] + omega * (f[ii] - f_eq[ii]);
	  }
      }

      //using DISTS = distribn_t[lb::lattices::D3Q15::NUMVECTORS];
      template<typename Kernel>
      void CompareHydros(const distribn_t expectedDensity,
			 const util::Vector3D<distribn_t>& expectedMomentum,
			 const distribn_t (&expectedFEq)[Kernel::LatticeType::NUMVECTORS],
			 lb::kernels::HydroVars<Kernel> &hydroVars,
			 distribn_t allowedError)
      {
	// Compare density
	REQUIRE(Approx(expectedDensity).margin(allowedError) == hydroVars.density);
	// Compare momentum
	REQUIRE(ApproxVector<distribn_t>(expectedMomentum).Margin(allowedError) == hydroVars.momentum);

	// Compare equilibrium f
	for (unsigned int ii = 0; ii < lb::lattices::D3Q15::NUMVECTORS; ++ii) {
	  REQUIRE(Approx(expectedFEq[ii]).margin(allowedError) == hydroVars.GetFEq().f[ii]);
	}
      }

      template<typename LatticeType>
      void InitialiseAnisotropicTestData(site_t site, distribn_t* distribution)
      {
	for (unsigned int direction = 0; direction < LatticeType::NUMVECTORS; ++direction)
	  {
	    distribution[direction] = ((distribn_t) (direction + 1)) / 10.0 + ((distribn_t) (site)) / 100.0;
	  }
      }

      template<typename LatticeType>
      void InitialiseAnisotropicTestData(FourCubeLatticeData* latticeData)
      {
	for (site_t site = 0; site < latticeData->GetLocalFluidSiteCount(); ++site)
	  {
	    distribn_t fOld[LatticeType::NUMVECTORS];
	    InitialiseAnisotropicTestData<LatticeType> (site, (distribn_t*)fOld);
	    latticeData->SetFOld<LatticeType>(site, fOld);
	  }
      }


          /**
           * Updates a property cache for the macroscopic variables selected. This should have
           * identical behaviour to in UpdateSiteMinsAndMaxes in BaseStreamer where the cache
           * is normally populated.
           * @param latDat
           * @param cache
           * @param simState
           */
          template<typename Lattice>
          void UpdatePropertyCache(geometry::LatticeData& latDat,
                                          lb::MacroscopicPropertyCache& cache,
                                          lb::SimulationState& simState)
          {
            for (site_t site = 0; site < latDat.GetLocalFluidSiteCount(); ++site)
            {
              distribn_t density, feq[Lattice::NUMVECTORS];
              util::Vector3D<distribn_t> momentum;
              util::Vector3D<distribn_t> velocity;

              Lattice::CalculateDensityMomentumFEq(latDat.GetSite(site).GetFOld<Lattice> (),
                                                   density,
                                                   momentum[0],
                                                   momentum[1],
                                                   momentum[2],
                                                   velocity[0],
                                                   velocity[1],
                                                   velocity[2],
                                                   feq);

              if (cache.densityCache.RequiresRefresh())
              {
                cache.densityCache.Put(site, density);
              }
              if (cache.velocityCache.RequiresRefresh())
              {
                cache.velocityCache.Put(site, velocity);
              }

              // TODO stress cache filling not yet implemented.
            }
          }

          /**
           * For every site in latticeData and every link crossing a boundary, set the distance to the
           * boundary to be wallDistance lattice length units
           *
           * @param latticeData Lattice object
           * @param wallDistance Distance to the wall in lattice units and in [0,1)
           */
          template<typename Lattice>
          void SetWallAndIoletDistances(FourCubeLatticeData& latticeData, distribn_t wallDistance)
          {
            assert(wallDistance >= 0);
            assert(wallDistance < 1);

            for (site_t siteIndex = 0; siteIndex < latticeData.GetTotalFluidSites(); siteIndex++)
            {
              for (Direction direction = 1; direction < Lattice::NUMVECTORS; direction++)
              {
                // -1 means that the a given link does not cross any boundary
                if (latticeData.GetSite(siteIndex).template GetWallDistance<Lattice> (direction) != (distribn_t) -1)
                {
                  latticeData.SetBoundaryDistance(siteIndex, direction, wallDistance);
                }
              }
            }
          }

    }
  }
}

#endif // ONCE
