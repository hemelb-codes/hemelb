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

      template<typename CollisionImpl>
      class SimpleCollideAndStream : public BaseStreamer<SimpleCollideAndStream<CollisionImpl> >
      {
        public:
          typedef CollisionImpl CollisionType;

        private:
          CollisionType collider;

        public:
          SimpleCollideAndStream(kernels::InitParams& initParams) :
              collider(initParams)
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
              const geometry::Site site = latDat->GetSite(siteIndex);

              distribn_t* lFOld = site.GetFOld();

              kernels::HydroVars<typename CollisionType::CKernel> hydroVars(lFOld);

              collider.CalculatePreCollision(hydroVars, site);

              collider.Collide(lbmParams, hydroVars);

              for (unsigned int ii = 0; ii < CollisionType::CKernel::LatticeType::NUMVECTORS; ii++)
              {
                * (latDat->GetFNew(site.GetStreamedIndex(ii))) = lFOld[ii] = hydroVars.GetFPostCollision()[ii];

              }

              BaseStreamer<SimpleCollideAndStream>::template UpdateMinsAndMaxes<tDoRayTracing>(hydroVars.v_x,
                                                                                               hydroVars.v_y,
                                                                                               hydroVars.v_z,
                                                                                               site,
                                                                                               hydroVars.GetFNeq().f,
                                                                                               hydroVars.density,
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

          inline void DoReset(kernels::InitParams* init)
          {
            collider.Reset(init);
          }

      };
    }
  }
}

#endif /* HEMELB_LB_STREAMERS_SIMPLECOLLIDEANDSTREAM */
