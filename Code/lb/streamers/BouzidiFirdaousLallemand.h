// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_LB_STREAMERS_BOUZIDIFIRDAOUSLALLEMAND_H
#define HEMELB_LB_STREAMERS_BOUZIDIFIRDAOUSLALLEMAND_H

#include "lb/streamers/BaseStreamer.h"
#include "lb/streamers/SimpleCollideAndStreamDelegate.h"
#include "lb/streamers/BouzidiFirdaousLallemandDelegate.h"

namespace hemelb
{
  namespace lb
  {
    namespace streamers
    {
      /**
       * Implement the boundary condition described in Bouzidi, Firdaous and
       * Lallemand "Momentum transfer of a Boltzmann-lattice fluid with
       * boundaries", Phys. Fluids 13, 3452-3459 (2001).
       *
       * This is based on the idea of doing interpolated bounce-back.
       *
       * Note that since the method requires data from neighbouring sites (in
       * some circumstances), it has a DoPostStep method.
       */
      template<typename CollisionImpl>
      class BouzidiFirdaousLallemand : public BaseStreamer<BouzidiFirdaousLallemand<CollisionImpl> >
      {
        public:
          typedef CollisionImpl CollisionType;

        private:
          CollisionType collider;
          SimpleCollideAndStreamDelegate<CollisionType> bulkLinkDelegate;
          BouzidiFirdaousLallemandDelegate<CollisionType> wallLinkDelegate;
          typedef typename CollisionType::CKernel::LatticeType LatticeType;

        public:
          BouzidiFirdaousLallemand(kernels::InitParams& initParams) :
            collider(initParams), bulkLinkDelegate(collider, initParams), wallLinkDelegate(collider, initParams)
          {

          }

          template<bool tDoRayTracing>
          inline void DoStreamAndCollide(const site_t firstIndex,
                                         const site_t siteCount,
                                         const LbmParameters* const lbmParams,
                                         geometry::LatticeData* const latticeData,
                                         lb::MacroscopicPropertyCache& propertyCache)
          {
            for (site_t siteIndex = firstIndex; siteIndex < (firstIndex + siteCount); ++siteIndex)
            {
              geometry::Site<geometry::LatticeData> site = latticeData->GetSite(siteIndex);

              const distribn_t* distribution = site.GetFOld<LatticeType> ();

              kernels::HydroVars<typename CollisionType::CKernel> hydroVars(distribution);

              ///< @todo #126 This value of tau will be updated by some kernels within the collider code (e.g. LBGKNN). It would be nicer if tau is handled in a single place.
              hydroVars.tau = lbmParams->GetTau();

              // In the first step, we stream and collide as we would for the SimpleCollideAndStream
              // streamer.
              collider.CalculatePreCollision(hydroVars, site);

              collider.Collide(lbmParams, hydroVars);

              for (unsigned int direction = 0; direction < LatticeType::NUMVECTORS; direction++)
              {
                if (site.HasBoundary(direction))
                {
                  wallLinkDelegate.StreamLink(latticeData, site, hydroVars, direction);
                }
                else
                {
                  // This is a standard link to another fluid site.
                  bulkLinkDelegate.StreamLink(latticeData, site, hydroVars, direction);
                }
              }

              BaseStreamer<BouzidiFirdaousLallemand>::template UpdateMinsAndMaxes<tDoRayTracing>(site,
                                                                                                 hydroVars,
                                                                                                 lbmParams,
                                                                                                 propertyCache);
            }
          }

          template<bool tDoRayTracing>
          inline void DoPostStep(const site_t firstIndex,
                                 const site_t siteCount,
                                 const LbmParameters* lbmParameters,
                                 geometry::LatticeData* latticeData,
                                 lb::MacroscopicPropertyCache& propertyCache)
          {
            for (site_t siteIndex = firstIndex; siteIndex < (firstIndex + siteCount); siteIndex++)
            {
              geometry::Site<geometry::LatticeData> site = latticeData->GetSite(siteIndex);
              for (unsigned int direction = 0; direction < LatticeType::NUMVECTORS; direction++)
              {
                if (site.HasBoundary(direction))
                {
                  wallLinkDelegate.PostStepLink(latticeData, site, direction);
                }
              }
            }
          }
      };
    }
  }
}

#endif /* HEMELB_LB_STREAMERS_BOUZIDIFIRDAOUSLALLEMAND_H */
