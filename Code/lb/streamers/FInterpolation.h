#ifndef HEMELB_LB_STREAMERS_FINTERPOLATION_H
#define HEMELB_LB_STREAMERS_FINTERPOLATION_H

#include "lb/streamers/BaseStreamer.h"

namespace hemelb
{
  namespace lb
  {
    namespace streamers
    {

      template<typename CollisionType>
      class FInterpolation : public BaseStreamer<FInterpolation<CollisionType> >
      {
        private:
          CollisionType collider;

        public:
          FInterpolation(kernels::InitParams& initParams) :
              collider(initParams)
          {

          }

          template<bool tDoRayTracing>
          inline void DoStreamAndCollide(const site_t firstIndex,
                                         const site_t siteCount,
                                         const LbmParameters* const iLbmParams,
                                         geometry::LatticeData* const latticeData,
                                         hemelb::vis::Control *visControl)
          {
            for (site_t index = firstIndex; index < (firstIndex + siteCount); index++)
            {
              distribn_t* distribution = latticeData->GetFOld(index * D3Q15::NUMVECTORS);

              kernels::HydroVars<typename CollisionType::CKernel> hydroVars(distribution);

              // In the first step, we stream and collide as we would for the SimpleCollideAndStream
              // streamer.
              collider.CalculatePreCollision(hydroVars, index - firstIndex);

              collider.Collide(iLbmParams, hydroVars);

              for (unsigned int direction = 0; direction < D3Q15::NUMVECTORS; direction++)
              {
                // Note that the post-step of this boundary condition relies on the post-collsion
                // value being written over f_old.
                distribution[direction] = hydroVars.GetFPostCollision()[direction];

                * (latticeData->GetFNew(latticeData->GetStreamedIndex(index, direction))) =
                    distribution[direction];
              }

              BaseStreamer<FInterpolation>::template UpdateMinsAndMaxes<tDoRayTracing>(hydroVars.v_x,
                                                                                       hydroVars.v_y,
                                                                                       hydroVars.v_z,
                                                                                       index,
                                                                                       hydroVars.GetFNeq().f,
                                                                                       hydroVars.density,
                                                                                       latticeData,
                                                                                       iLbmParams,
                                                                                       visControl);
            }
          }

          template<bool tDoRayTracing>
          inline void DoPostStep(const site_t firstIndex,
                                 const site_t siteCount,
                                 const LbmParameters* lbmParameters,
                                 geometry::LatticeData* latticeData,
                                 hemelb::vis::Control *visControl)
          {
            for (site_t siteIndex = firstIndex; siteIndex < (firstIndex + siteCount); siteIndex++)
            {
              // Iterate over the direction indices.
              for (unsigned int direction = 1; direction < D3Q15::NUMVECTORS; direction++)
              {
                // If there's a boundary in that direction and none in the other direction, do the
                // f-interpolation.
                if (latticeData->HasBoundary(siteIndex, direction))
                {
                  int inverseDirection = D3Q15::INVERSEDIRECTIONS[direction];

                  if (!latticeData->HasBoundary(siteIndex, inverseDirection))
                  {
                    // Calculate 2 x the distance to the boundary.
                    distribn_t twoQ = 2.0 * latticeData->GetCutDistance(siteIndex, direction);

                    distribn_t thisDirectionNew = *latticeData->GetFNew(siteIndex
                        * D3Q15::NUMVECTORS + direction);
                    distribn_t thisDirectionOld = *latticeData->GetFOld(siteIndex
                        * D3Q15::NUMVECTORS + direction);
                    distribn_t oppDirectionOld = *latticeData->GetFOld(siteIndex * D3Q15::NUMVECTORS
                        + inverseDirection);

                    // Interpolate between the values of the f direction to work out a new streamed value.
                    distribn_t streamed = (twoQ < 1.0) ?
                      (thisDirectionNew + twoQ * (thisDirectionOld - thisDirectionNew)) :
                      (oppDirectionOld + (1. / twoQ) * (thisDirectionOld - oppDirectionOld));

                    // This streamed value is assigned to the f-distribution in the direction facing
                    // away from the boundary.
                    * (latticeData->GetFNew(siteIndex * D3Q15::NUMVECTORS + inverseDirection)) =
                        streamed;
                  }
                  // If there are boundaries in both directions perform simple bounce-back using the
                  // post-collision values in f_old.
                  else
                  {
                    * (latticeData->GetFNew(siteIndex * D3Q15::NUMVECTORS + inverseDirection)) =
                        *latticeData->GetFOld(siteIndex * D3Q15::NUMVECTORS + direction);
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
