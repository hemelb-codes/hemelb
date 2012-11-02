// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_LB_STREAMERS_FINTERPOLATION_H
#define HEMELB_LB_STREAMERS_FINTERPOLATION_H

#include "lb/streamers/BaseStreamer.h"

namespace hemelb
{
  namespace lb
  {
    namespace streamers
    {

      template<typename CollisionImpl>
      class FInterpolation : public BaseStreamer<FInterpolation<CollisionImpl> >
      {
        public:
          typedef CollisionImpl CollisionType;

        private:
          CollisionType collider;
          typedef typename CollisionType::CKernel::LatticeType LatticeType;

        public:
          FInterpolation(kernels::InitParams& initParams) :
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
            for (site_t index = firstIndex; index < (firstIndex + siteCount); index++)
            {
              geometry::Site<geometry::LatticeData> site = latticeData->GetSite(index);

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
                * (latticeData->GetFNew(site.GetStreamedIndex<LatticeType> (direction)))
                    = hydroVars.GetFPostCollision()[direction];
              }

              BaseStreamer<FInterpolation>::template UpdateMinsAndMaxes<tDoRayTracing>(site,
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

              kernels::HydroVars<typename CollisionType::CKernel> hydroVars(site.GetFOld<LatticeType> ());

              ///< @todo #126 This value of tau will be updated by some kernels within the collider code (e.g. LBGKNN). It would be nicer if tau is handled in a single place.
              hydroVars.tau = lbmParameters->GetTau();

              // In the first step, we stream and collide as we would for the SimpleCollideAndStream
              // streamer.
              collider.CalculatePreCollision(hydroVars, site);

              collider.Collide(lbmParameters, hydroVars);

              // Iterate over the direction indices.
              for (unsigned int direction = 1; direction < LatticeType::NUMVECTORS; direction++)
              {
                // If there's a boundary in that direction and none in the other direction, do the
                // f-interpolation.
                if (site.HasBoundary(direction))
                {
                  int inverseDirection = LatticeType::INVERSEDIRECTIONS[direction];

                  if (!site.HasBoundary(inverseDirection))
                  {
                    // Calculate 2 x the distance to the boundary.
                    distribn_t twoQ = 2.0 * site.GetWallDistance<LatticeType> (direction);

                    distribn_t thisDirectionNew = *latticeData->GetFNew(siteIndex
                        * CollisionType::CKernel::LatticeType::NUMVECTORS + direction);
                    distribn_t thisDirectionOld = hydroVars.GetFPostCollision()[direction];
                    distribn_t oppDirectionOld = hydroVars.GetFPostCollision()[inverseDirection];

                    // Interpolate between the values of the f direction to work out a new streamed value.
                    distribn_t streamed = (twoQ < 1.0)
                      ? (thisDirectionNew + twoQ * (thisDirectionOld - thisDirectionNew))
                      : (oppDirectionOld + (1. / twoQ) * (thisDirectionOld - oppDirectionOld));

                    // This streamed value is assigned to the f-distribution in the direction facing
                    // away from the boundary.
                    * (latticeData->GetFNew(siteIndex * LatticeType::NUMVECTORS + inverseDirection)) = streamed;
                  }
                  // If there are boundaries in both directions perform simple bounce-back using the
                  // post-collision values in f_old.

                  else
                  {
                    * (latticeData->GetFNew(siteIndex * LatticeType::NUMVECTORS + inverseDirection))
                        = hydroVars.GetFPostCollision()[direction];
                  }
                }
              }
            }
          }

      };
    }
  }
}

#endif /* HEMELB_LB_STREAMERS_FINTERPOLATION_H */
