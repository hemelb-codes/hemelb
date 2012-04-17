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
              const geometry::Site site = latticeData->GetSite(index);

              distribn_t* distribution = site.GetFOld<LatticeType> ();

              kernels::HydroVars<typename CollisionType::CKernel> hydroVars(distribution);

              ///< @todo #126 This value of tau will be updated by some kernels within the collider code (e.g. LBGKNN). It would be nicer if tau is handled in a single place.
              hydroVars.tau = lbmParams->GetTau();

              // In the first step, we stream and collide as we would for the SimpleCollideAndStream
              // streamer.
              collider.CalculatePreCollision(hydroVars, site);

              collider.Collide(lbmParams, hydroVars);

              for (unsigned int direction = 0; direction < LatticeType::NUMVECTORS; direction++)
              {
                // Note that the post-step of this boundary condition relies on the post-collsion
                // value being written over f_old.
                distribution[direction] = hydroVars.GetFPostCollision()[direction];

                * (latticeData->GetFNew(site.GetStreamedIndex<LatticeType> (direction))) = distribution[direction];
              }

              BaseStreamer<FInterpolation>::template UpdateMinsAndMaxes<tDoRayTracing>(hydroVars.v_x,
                                                                                       hydroVars.v_y,
                                                                                       hydroVars.v_z,
                                                                                       site,
                                                                                       hydroVars.GetFNeq().f,
                                                                                       hydroVars.density,
                                                                                       hydroVars.tau,
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
              geometry::Site site = latticeData->GetSite(siteIndex);

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
                    distribn_t thisDirectionOld = site.GetFOld<LatticeType> ()[direction];
                    distribn_t oppDirectionOld = site.GetFOld<LatticeType> ()[inverseDirection];

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
                    * (latticeData->GetFNew(siteIndex * LatticeType::NUMVECTORS + inverseDirection)) = site.GetFOld<
                        LatticeType> ()[direction];
                  }
                }
              }
            }
          }

          inline void DoReset(kernels::InitParams* init)
          {
            collider.Reset(init);
          }

      };
    }
  }
}

#endif /* HEMELB_LB_STREAMERS_FINTERPOLATION_H */
