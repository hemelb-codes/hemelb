
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_STREAMERS_STREAMERTYPEFACTORY_H
#define HEMELB_LB_STREAMERS_STREAMERTYPEFACTORY_H

#include "lb/kernels/BaseKernel.h"
#include "lb/streamers/BaseStreamer.h"
#include "lb/streamers/SimpleCollideAndStreamDelegate.h"
#include <cuda_runtime.h>

#define CUDA_SAFE_CALL(x)                           \
{                                                   \
    cudaError_t error = x;                          \
    if ( error != cudaSuccess ) {                   \
      const char *name = cudaGetErrorName(error);   \
      const char *str = cudaGetErrorString(error);  \
      std::cerr << "\n"                             \
                << "CUDA Error at " #x "\n"         \
                << name << ": " << error << "\n";   \
      exit(1);                                      \
    }                                               \
}

namespace hemelb
{
  namespace lb
  {
    namespace streamers
    {
      __global__ void WallStreamerTypeFactory_DoStreamAndCollideKernel(
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
            static bool init = true;
            static site_t* neighbourIndices_dev;
            static unsigned* wallIntersections_dev;
            static distribn_t* fOld_dev;
            static distribn_t* fNew_dev;

            unsigned numSites = latDat->GetLocalFluidSiteCount();

            if ( init )
            {
              CUDA_SAFE_CALL(cudaMalloc(&neighbourIndices_dev, numSites * LatticeType::NUMVECTORS * sizeof(site_t)));
              CUDA_SAFE_CALL(cudaMalloc(&wallIntersections_dev, numSites * sizeof(unsigned)));
              CUDA_SAFE_CALL(cudaMalloc(&fOld_dev, (numSites * LatticeType::NUMVECTORS + 1) * sizeof(distribn_t)));
              CUDA_SAFE_CALL(cudaMalloc(&fNew_dev, (numSites * LatticeType::NUMVECTORS + 1) * sizeof(distribn_t)));

              std::vector<unsigned> wallIntersections(numSites);

              for ( site_t siteIndex = 0; siteIndex < numSites; ++siteIndex )
              {
                wallIntersections[siteIndex] = latDat->GetSite(siteIndex).GetSiteData().GetWallIntersectionData();
              }

              CUDA_SAFE_CALL(cudaMemcpy(neighbourIndices_dev, latDat->GetNeighbourIndices().data(), numSites * LatticeType::NUMVECTORS * sizeof(site_t), cudaMemcpyHostToDevice));
              CUDA_SAFE_CALL(cudaMemcpy(wallIntersections_dev, wallIntersections.data(), numSites * sizeof(unsigned), cudaMemcpyHostToDevice));

              init = false;
            }

            if (lbmParams->UseGPU() && !propertyCache.RequiresRefresh())
            {
            // copy fOld from host to device
            CUDA_SAFE_CALL(cudaMemcpy(fOld_dev, latDat->GetSite(0).GetFOld<LatticeType>(), numSites * LatticeType::NUMVECTORS * sizeof(distribn_t), cudaMemcpyHostToDevice));

            // launch WallStreamer_DoStreamAndCollide kernel
            const int BLOCK_SIZE = 256;
            const int GRID_SIZE = (siteCount + BLOCK_SIZE - 1) / BLOCK_SIZE;

            WallStreamerTypeFactory_DoStreamAndCollideKernel<<<GRID_SIZE, BLOCK_SIZE>>>(
              firstIndex,
              siteCount,
              lbmParams->GetTau(),
              lbmParams->GetOmega(),
              neighbourIndices_dev,
              wallIntersections_dev,
              fOld_dev,
              fNew_dev
            );
            CUDA_SAFE_CALL(cudaGetLastError());

            CUDA_SAFE_CALL(cudaDeviceSynchronize());

            // copy fNew from device to host
            CUDA_SAFE_CALL(cudaMemcpy(latDat->GetFNew(0), fNew_dev, numSites * LatticeType::NUMVECTORS * sizeof(distribn_t), cudaMemcpyDeviceToHost));
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
