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

namespace hemelb
{
  namespace lb
  {
    namespace streamers
    {

      template<typename CollisionImpl>
      class BouzidiFirdaousLallemand : public BaseStreamer<BouzidiFirdaousLallemand<CollisionImpl> >
      {
        public:
          typedef CollisionImpl CollisionType;

        private:
          CollisionType collider;
          typedef typename CollisionType::CKernel::LatticeType LatticeType;

        public:
          BouzidiFirdaousLallemand(kernels::InitParams& initParams) :
            collider(initParams)
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
                  site_t invDirection = LatticeType::INVERSEDIRECTIONS[direction];
                  site_t bbDestination = (siteIndex * LatticeType::NUMVECTORS) + invDirection;
                  distribn_t q = site.GetWallDistance<LatticeType> (direction);

                  if (site.HasBoundary(invDirection) || q < 0.5)
                  {
                    // If there IS NO fluid site in the opposite direction, fall back to SBB.
                    // If there IS such a site, we have to wait for the site in the opposite
                    // direction to finish in order to complete this update. So just bounce-back
                    // the post collision f that we would have otherwise thrown away (to avoid
                    // having to collide twice).
                    * (latticeData->GetFNew(bbDestination)) = hydroVars.GetFPostCollision()[direction];
                  }
                  else
                  {
                    // We have a fluid site and have all the data needed to complete this direction!
                    // Implement Eq (5b) from Bouzidi et al.
                    * (latticeData->GetFNew(bbDestination)) = (hydroVars.GetFPostCollision()[direction] + (2.0 * q - 1)
                        * hydroVars.GetFPostCollision()[invDirection]) / (2.0 * q);
                  }
                }
                else
                {
                  // This is a standard link to another fluid site.
                  * (latticeData->GetFNew(site.GetStreamedIndex<LatticeType> (direction)))
                      = hydroVars.GetFPostCollision()[direction];
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
              distribn_t* fNew = latticeData->GetFNew(siteIndex * LatticeType::NUMVECTORS);

              for (unsigned int direction = 0; direction < LatticeType::NUMVECTORS; direction++)
              {
                if (site.HasBoundary(direction))
                {
                  site_t invDirection = LatticeType::INVERSEDIRECTIONS[direction];
                  distribn_t q = site.GetWallDistance<LatticeType> (direction);
                  // If there is no fluid site in the opposite direction, fall back to simple
                  // bounce back, which has been done above.

                  // If q >= 0.5, then we handled that fully above also.
                  if (!site.HasBoundary(invDirection) && q < 0.5)
                  {
                    // So, we have a fluid site and all the data needed to complete this direction!
                    // Implement Eq (5a) from Bouzidi et al.

                    // Note that:
                    // - fNew[direction] is the newly-arrived fPostColl[direction] from the neighbouring site
                    // - fNew[invDirection] is the above-bounced-back fPostColl[direction] for this site.
                    fNew[invDirection] = 2.0 * q * fNew[invDirection] + (1.0 - 2.0 * q) * fNew[direction];
                  }
                }
              }
            }
          }
      };
    }
  }
}

#endif /* HEMELB_LB_STREAMERS_BOUZIDIFIRDAOUSLALLEMAND_H */
