#ifndef HEMELB_LB_STREAMERS_SIMPLECOLLIDEANDSTREAM_H
#define HEMELB_LB_STREAMERS_SIMPLECOLLIDEANDSTREAM_H

#include "lb/streamers/BaseStreamer.h"
#include "lb/kernels/BaseKernel.h"
#include "lb/HFunction.h"

namespace hemelb
{
  namespace lb
  {
    namespace streamers
    {

      template<typename CollisionType>
      class SimpleCollideAndStream : public BaseStreamer<SimpleCollideAndStream<CollisionType> >
      {
        private:
          CollisionType collider;

        public:
          SimpleCollideAndStream(kernels::InitParams& initParams) :
            collider(initParams)
          {

          }

          template<bool tDoRayTracing>
          void DoStreamAndCollide(const site_t iFirstIndex,
                                  const site_t iSiteCount,
                                  const LbmParameters* iLbmParams,
                                  geometry::LatticeData* bLatDat,
                                  hemelb::vis::Control *iControl)
          {
            for (site_t iIndex = iFirstIndex; iIndex < (iFirstIndex + iSiteCount); iIndex++)
            {
              distribn_t* lFOld = bLatDat->GetFOld(iIndex * D3Q15::NUMVECTORS);

              kernels::HydroVars<typename CollisionType::CKernel> hydroVars(lFOld);

              collider.CalculatePreCollision(hydroVars, iIndex - iFirstIndex);

              for (unsigned int ii = 0; ii < D3Q15::NUMVECTORS; ii++)
              {
                * (bLatDat->GetFNew(bLatDat->GetStreamedIndex(iIndex, ii))) = lFOld[ii]
                    = collider.Collide(iLbmParams, ii, hydroVars);
              }

              BaseStreamer<SimpleCollideAndStream>::template UpdateMinsAndMaxes<tDoRayTracing>(hydroVars.v_x,
                                                                                               hydroVars.v_y,
                                                                                               hydroVars.v_z,
                                                                                               iIndex,
                                                                                               hydroVars.GetFNeq().f,
                                                                                               hydroVars.density,
                                                                                               bLatDat,
                                                                                               iLbmParams,
                                                                                               iControl);
            }
          }

          template<bool tDoRayTracing>
          void DoPostStep(const site_t iFirstIndex,
                          const site_t iSiteCount,
                          const LbmParameters* iLbmParams,
                          geometry::LatticeData* bLatDat,
                          hemelb::vis::Control *iControl)
          {

          }

          void DoReset(kernels::InitParams* init)
          {
            collider.Reset(init);
          }

      };
    }
  }
}

#endif /* HEMELB_LB_STREAMERS_SIMPLECOLLIDEANDSTREAM */
