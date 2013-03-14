// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_LB_STREAMERS_LADDIOLET_H
#define HEMELB_LB_STREAMERS_LADDIOLET_H

#include "lb/kernels/BaseKernel.h"
#include "lb/streamers/BaseStreamer.h"

namespace hemelb
{
  namespace lb
  {
    namespace streamers
    {
      template<typename CollisionImpl>
      class LaddIolet : public BaseStreamer<LaddIolet<CollisionImpl> >
      {
        public:
          typedef CollisionImpl CollisionType;

        private:
          CollisionType collider;
          typedef typename CollisionType::CKernel::LatticeType LatticeType;

          static inline site_t GetBBIndex(site_t siteIndex, int direction)
          {
            return (siteIndex * LatticeType::NUMVECTORS) + LatticeType::INVERSEDIRECTIONS[direction];
          }

        public:
          LaddIolet(kernels::InitParams& initParams) :
            collider(initParams), boundaries(initParams.boundaryObject)
          {

          }

          template<bool tDoRayTracing>
          inline void DoStreamAndCollide(const site_t firstIndex,
                                         const site_t siteCount,
                                         const LbmParameters* lbmParams,
                                         geometry::LatticeData* latDat,
                                         lb::MacroscopicPropertyCache& propertyCache)
          {
            for (site_t siteIndex = firstIndex; siteIndex < (firstIndex + siteCount); siteIndex++)
            {
              geometry::Site<geometry::LatticeData> site = latDat->GetSite(siteIndex);

              const distribn_t* fOld = site.GetFOld<LatticeType> ();

              kernels::HydroVars<typename CollisionType::CKernel> hydroVars(fOld);

              ///< @todo #126 This value of tau will be updated by some kernels within the collider code (e.g. LBGKNN). It would be nicer if tau is handled in a single place.
              hydroVars.tau = lbmParams->GetTau();

              collider.CalculatePreCollision(hydroVars, site);

              collider.Collide(lbmParams, hydroVars);

              int boundaryId = site.GetBoundaryId();
              boundaries::iolets::InOutLetParabolicVelocity* iolet =
                  dynamic_cast<boundaries::iolets::InOutLetParabolicVelocity*> (boundaries->GetLocalIolet(boundaryId));
              LatticePosition sitePos(site.GetGlobalSiteCoords());

              for (Direction ii = 0; ii < LatticeType::NUMVECTORS; ii++)
              {
                // Translating from Ladd, J. Fluid Mech. "Numerical simulations
                // of particulate suspensions via a discretized Boltzmann
                // equation. Part 1. Theoretical foundation", 1994
                // Eq (3.2) -- simple bounce-back -- becomes:
                //   f_i'(r, t+1) = f_i(r, t*)
                // Eq (3.3) --- modified BB -- becomes:
                //   f_i'(r, t+1) = f_i(r, t*) - 2 a1_i \rho u . c_i
                // where u is the velocity of the boundary half way along the
                // link and a1_i = w_1 / cs2

                // The actual bounce-back lines, including streaming and collision. Basically swap
                // the non-equilibrium components of f in each of the opposing pairs of directions.
                // We always bounce back into the same component of f, irrespective of whether it's
                // wall or iolet.
                // NOTE: site.HasBoundary only returns true if that link intersects a wall, not an
                // iolet, so we have to call GetBBIndex below if it is an iolet.
                site_t streamingDestination = site.HasBoundary(ii)
                  ? GetBBIndex(siteIndex, ii)
                  : site.GetStreamedIndex<LatticeType> (ii);

                // We now must compute the Ladd correction, which is zero in
                // the case of BB at a wall.
                distribn_t correction = 0.;
                if (site.HasIolet(ii))
                {
                  streamingDestination = GetBBIndex(siteIndex, ii);
                  LatticePosition halfWay(sitePos);
                  halfWay.x += 0.5 * LatticeType::CX[ii];
                  halfWay.y += 0.5 * LatticeType::CY[ii];
                  halfWay.z += 0.5 * LatticeType::CZ[ii];

                  LatticeVelocity wallVel(iolet->GetVelocityAtPosition(halfWay));

                  correction = 2. * LatticeType::EQMWEIGHTS[ii] * hydroVars.density * (wallVel.x * LatticeType::CX[ii]
                      + wallVel.y * LatticeType::CY[ii] + wallVel.z * LatticeType::CZ[ii]) / Cs2;
                }

                // Remember, oFNeq currently hold the equilibrium distribution. We
                // simultaneously use this and correct it, here.
                * (latDat->GetFNew(streamingDestination)) = hydroVars.GetFPostCollision()[ii] - correction;
              }

              //TODO: Necessary to specify sub-class?
              BaseStreamer<LaddIolet>::template UpdateMinsAndMaxes<tDoRayTracing>(site,
                                                                                  hydroVars,
                                                                                  lbmParams,
                                                                                  propertyCache);
            }
          }

          template<bool tDoRayTracing>
          inline void DoPostStep(const site_t iFirstIndex,
                                 const site_t iSiteCount,
                                 const LbmParameters* iLbmParams,
                                 geometry::LatticeData* bLatDat,
                                 lb::MacroscopicPropertyCache& propertyCache)
          {

          }

        private:
          boundaries::BoundaryValues* boundaries;
      };
    }
  }
}

#endif /* HEMELB_LB_STREAMERS_LADDIOLET_H */
