#ifndef HEMELB_LB_STREAMERS_NASHZEROTHORDERPRESSUREBB_H
#define HEMELB_LB_STREAMERS_NASHZEROTHORDERPRESSUREBB_H

#include "lb/streamers/BaseStreamer.h"
#include "lb/streamers/NashZerothOrderPressureDelegate.h"
#include "util/utilityFunctions.h"

namespace hemelb
{
  namespace lb
  {
    namespace streamers
    {
      template<typename CollisionImpl>
      class NashZerothOrderPressureBB : public BaseStreamer<NashZerothOrderPressureBB<CollisionImpl> >
      {
        public:
          typedef CollisionImpl CollisionType;
          typedef typename CollisionType::CKernel::LatticeType LatticeType;

        private:
          CollisionType collider;
          SimpleCollideAndStreamDelegate<CollisionType> bulkLinkDelegate;
          SimpleBounceBackDelegate<CollisionType> wallLinkDelegate;
          NashZerothOrderPressureDelegate<CollisionType> ioletLinkDelegate;

        public:
          NashZerothOrderPressureBB(kernels::InitParams& initParams) :
            collider(initParams), bulkLinkDelegate(collider, initParams), wallLinkDelegate(collider, initParams), ioletLinkDelegate(collider, initParams)
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

              const distribn_t* f = site.GetFOld<LatticeType> ();

              kernels::HydroVars<typename CollisionType::CKernel> hydroVars(f);

              // First calculate the density and macro-velocity
              collider.CalculatePreCollision(hydroVars, site);

              collider.Collide(lbmParams, hydroVars);

              // Start by doing the normal stream and collide operation.
              for (Direction direction = 0; direction < LatticeType::NUMVECTORS; ++direction)
              {
                if (site.HasIolet(direction))
                {
                  ioletLinkDelegate.StreamLink(latDat, site, hydroVars, direction);
                }
                else if (site.HasBoundary(direction))
                {
                  wallLinkDelegate.StreamLink(latDat, site, hydroVars, direction);
                }
                else
                {
                  bulkLinkDelegate.StreamLink(latDat, site, hydroVars, direction);
                }
              }

              ///< @todo #126 It would be nicer if tau is handled in a single place.
              hydroVars.tau = lbmParams->GetTau();

              BaseStreamer<NashZerothOrderPressureBB>::template UpdateMinsAndMaxes<tDoRayTracing>(site,
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

          void DoReset(kernels::InitParams* init)
          {
            collider.Reset(init);
          }
      };
    }
  }
}

#endif /* HEMELB_LB_STREAMERS_NASHZEROTHORDERPRESSUREBB_H */
