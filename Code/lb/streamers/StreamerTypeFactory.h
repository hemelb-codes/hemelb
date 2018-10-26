
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_STREAMERS_STREAMERTYPEFACTORY_H
#define HEMELB_LB_STREAMERS_STREAMERTYPEFACTORY_H

#include "lb/kernels/BaseKernel.h"
#include "lb/streamers/BaseStreamer.h"
#include "lb/streamers/SimpleCollideAndStreamDelegate.h"

namespace hemelb
{
  namespace lb
  {
    namespace streamers
    {
      void DoStreamAndCollideGPU(
        site_t firstIndex,
        site_t siteCount,
        distribn_t lbmParams_tau,
        distribn_t lbmParams_omega,
        const site_t* neighbourIndices,
        const unsigned* wallIntersections,
        const distribn_t* fOld,
        distribn_t* fNew
      );

      /**
       * Template to produce Streamers that can cope with fluid-fluid and
       * fluid-wall links. Requires two classes as arguments: 1) the Collision
       * class and 2) a StreamerDelegate class that will handle the wall links.
       *
       * It is intended that a simpler metafunction partially specialise this
       * template on WallLinkImpl.
       */
      template<typename CollisionImpl, typename WallLinkImpl>
      class WallStreamerTypeFactory : public BaseStreamer<WallStreamerTypeFactory<CollisionImpl, WallLinkImpl> >
      {
        public:
          typedef CollisionImpl CollisionType;

        private:
          CollisionType collider;
          SimpleCollideAndStreamDelegate<CollisionType> bulkLinkDelegate;
          WallLinkImpl wallLinkDelegate;

          typedef typename CollisionType::CKernel::LatticeType LatticeType;

        public:
          WallStreamerTypeFactory(kernels::InitParams& initParams) :
            collider(initParams), bulkLinkDelegate(collider, initParams), wallLinkDelegate(collider, initParams)
          {

          }

          template<bool tDoRayTracing>
          inline void DoStreamAndCollide(const site_t firstIndex,
                                         const site_t siteCount,
                                         const LbmParameters* lbmParams,
                                         geometry::LatticeData* latDat,
                                         lb::MacroscopicPropertyCache& propertyCache)
          {
            unsigned localFluidSites = latDat->GetLocalFluidSiteCount();
            site_t sharedFs = latDat->GetNumSharedFs();

            if ( siteCount == 0 )
            {
              return;
            }

            if (lbmParams->UseGPU() && !propertyCache.RequiresRefresh())
            {
            // copy fOld from host to device
            CUDA_SAFE_CALL(cudaMemcpyAsync(
              latDat->GetFOldGPU(),
              latDat->GetSite(0).GetFOld<LatticeType>(),
              (localFluidSites * LatticeType::NUMVECTORS + sharedFs + 1) * sizeof(distribn_t),
              cudaMemcpyHostToDevice
            ));
            CUDA_SAFE_CALL(cudaMemcpyAsync(
              latDat->GetFNewGPU(),
              latDat->GetFNew(0),
              (localFluidSites * LatticeType::NUMVECTORS + sharedFs + 1) * sizeof(distribn_t),
              cudaMemcpyHostToDevice
            ));

            // launch WallStreamer_DoStreamAndCollide kernel
            DoStreamAndCollideGPU(
              firstIndex,
              siteCount,
              lbmParams->GetTau(),
              lbmParams->GetOmega(),
              latDat->GetNeighbourIndicesGPU(),
              latDat->GetWallIntersectionsGPU(),
              latDat->GetFOldGPU(),
              latDat->GetFNewGPU()
            );
            CUDA_SAFE_CALL(cudaGetLastError());

            // copy fNew from device to host
            CUDA_SAFE_CALL(cudaMemcpy(
              latDat->GetFNew(0),
              latDat->GetFNewGPU(),
              (localFluidSites * LatticeType::NUMVECTORS + sharedFs + 1) * sizeof(distribn_t),
              cudaMemcpyDeviceToHost
            ));
            }

            else
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

              for (Direction ii = 0; ii < LatticeType::NUMVECTORS; ii++)
              {
                if (site.HasWall(ii))
                {
                  wallLinkDelegate.StreamLink(lbmParams, latDat, site, hydroVars, ii);
                }
                else
                {
                  bulkLinkDelegate.StreamLink(lbmParams, latDat, site, hydroVars, ii);
                }
              }

              //TODO: Necessary to specify sub-class?
              BaseStreamer<WallStreamerTypeFactory>::template UpdateMinsAndMaxes<tDoRayTracing>(site,
                                                                                                hydroVars,
                                                                                                lbmParams,
                                                                                                propertyCache);
            }
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
                if (site.HasWall(direction))
                {
                  wallLinkDelegate.PostStepLink(latticeData, site, direction);
                }
              }
            }
          }

      };

      /**
       * Template to produce Streamers that can cope with fluid-fluid and
       * fluid-iolet links. Requires two classes as arguments: 1) the Collision
       * class and 2) a StreamerDelegate class that will handle the iolet links.
       *
       * It is intended that a simpler metafunction partially specialise this
       * template on IoletLinkImpl.
       */
      template<typename CollisionImpl, typename IoletLinkImpl>
      class IoletStreamerTypeFactory : public BaseStreamer<IoletStreamerTypeFactory<CollisionImpl, IoletLinkImpl> >
      {
        public:
          typedef CollisionImpl CollisionType;

        private:
          CollisionType collider;
          SimpleCollideAndStreamDelegate<CollisionType> bulkLinkDelegate;
          IoletLinkImpl ioletLinkDelegate;

          typedef typename CollisionType::CKernel::LatticeType LatticeType;

        public:
          IoletStreamerTypeFactory(kernels::InitParams& initParams) :
            collider(initParams), bulkLinkDelegate(collider, initParams), ioletLinkDelegate(collider, initParams)
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

              for (Direction ii = 0; ii < LatticeType::NUMVECTORS; ii++)
              {
                if (site.HasIolet(ii))
                {
                  ioletLinkDelegate.StreamLink(lbmParams, latDat, site, hydroVars, ii);
                }
                else
                {
                  bulkLinkDelegate.StreamLink(lbmParams, latDat, site, hydroVars, ii);
                }
              }

              //TODO: Necessary to specify sub-class?
              BaseStreamer<IoletStreamerTypeFactory>::template UpdateMinsAndMaxes<tDoRayTracing>(site,
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
                if (site.HasIolet(direction))
                {
                  ioletLinkDelegate.PostStepLink(latticeData, site, direction);
                }
              }
            }
          }
      };

      /**
       * Template to produce Streamers that can cope with fluid-fluid,
       * fluid-wall and fluid-iolet links. Requires three classes as arguments:
       * 1) the Collision class,
       * 2) a StreamerDelegate class that will handle the wall links, and
       * 3) a StreamerDelegate class that will handle the iolet links.
       *
       * It is intended that a simpler metafunction partially specialise this
       * template on WallLinkImpl and IoletLinkImpl.
       */
      template<typename CollisionImpl, typename WallLinkImpl, typename IoletLinkImpl>
      class WallIoletStreamerTypeFactory : public BaseStreamer<WallIoletStreamerTypeFactory<CollisionImpl,
          WallLinkImpl, IoletLinkImpl> >
      {
        public:
          typedef CollisionImpl CollisionType;
          typedef typename CollisionType::CKernel::LatticeType LatticeType;

        private:
          CollisionType collider;
          SimpleCollideAndStreamDelegate<CollisionType> bulkLinkDelegate;
          WallLinkImpl wallLinkDelegate;
          IoletLinkImpl ioletLinkDelegate;

        public:
          WallIoletStreamerTypeFactory(kernels::InitParams& initParams) :
            collider(initParams), bulkLinkDelegate(collider, initParams), wallLinkDelegate(collider, initParams),
                ioletLinkDelegate(collider, initParams)
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

              for (Direction ii = 0; ii < LatticeType::NUMVECTORS; ii++)
              {
                if (site.HasIolet(ii))
                {
                  ioletLinkDelegate.StreamLink(lbmParams, latDat, site, hydroVars, ii);
                }
                else if (site.HasWall(ii))
                {
                  wallLinkDelegate.StreamLink(lbmParams, latDat, site, hydroVars, ii);
                }
                else
                {
                  bulkLinkDelegate.StreamLink(lbmParams, latDat, site, hydroVars, ii);
                }
              }

              //TODO: Necessary to specify sub-class?
              BaseStreamer<WallIoletStreamerTypeFactory>::template UpdateMinsAndMaxes<tDoRayTracing>(site,
                                                                                                     hydroVars,
                                                                                                     lbmParams,
                                                                                                     propertyCache);
            }
          }

          template<bool tDoRayTracing>
          inline void DoPostStep(const site_t firstIndex,
                                 const site_t siteCount,
                                 const LbmParameters* lbmParams,
                                 geometry::LatticeData* latticeData,
                                 lb::MacroscopicPropertyCache& propertyCache)
          {
            for (site_t siteIndex = firstIndex; siteIndex < (firstIndex + siteCount); siteIndex++)
            {
              geometry::Site<geometry::LatticeData> site = latticeData->GetSite(siteIndex);
              for (unsigned int direction = 0; direction < LatticeType::NUMVECTORS; direction++)
              {
                if (site.HasWall(direction))
                {
                  wallLinkDelegate.PostStepLink(latticeData, site, direction);
                }
                else if (site.HasIolet(direction))
                {
                  ioletLinkDelegate.PostStepLink(latticeData, site, direction);
                }
              }
            }

          }
      };
    }
  }
}
#endif // HEMELB_LB_STREAMERS_STREAMERTYPEFACTORY_H
