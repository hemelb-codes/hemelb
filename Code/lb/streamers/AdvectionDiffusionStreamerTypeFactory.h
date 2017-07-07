
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_STREAMERS_ADVECTIONDIFFUSIONSTREAMERTYPEFACTORY_H
#define HEMELB_LB_STREAMERS_ADVECTIONDIFFUSIONSTREAMERTYPEFACTORY_H

#include "lb/kernels/BaseKernel.h"
#include "lb/streamers/BaseStreamer.h"
#include "lb/streamers/SimpleCollideAndStreamDelegate.h"

namespace hemelb
{
  namespace lb
  {
    namespace streamers
    {
      /**
       * Template to produce Streamers that can cope with fluid-fluid and
       * fluid-wall links. Requires two classes as arguments: 1) the Collision
       * class and 2) a StreamerDelegate class that will handle the wall links.
       *
       * It is intended that a simpler metafunction partially specialise this
       * template on WallLinkImpl.
       */
      template<typename CollisionImpl, typename WallLinkImpl>
      class VesselWallStreamerTypeFactory : public AdvectionDiffusionBaseStreamer<VesselWallStreamerTypeFactory<CollisionImpl, WallLinkImpl> >
      {
        public:
          typedef CollisionImpl CollisionType;

        private:
          CollisionType collider;
          SimpleCollideAndStreamDelegate<CollisionType> bulkLinkDelegate;
          WallLinkImpl wallLinkDelegate;

          typedef typename CollisionType::CKernel::LatticeType LatticeType;

        public:
          VesselWallStreamerTypeFactory(kernels::InitParams& initParams) :
            collider(initParams), bulkLinkDelegate(collider, initParams), wallLinkDelegate(collider, initParams)
          {

          }

          template<bool tDoRayTracing>
          inline void DoStreamAndCollide(const site_t firstIndex,
                                         const site_t siteCount,
                                         const LbmParameters* lbmParams,
                                         geometry::LatticeData* latDat,
                                         lb::MacroscopicPropertyCache& propertyCache,
                                         lb::MacroscopicPropertyCache& coupledPropertyCache)
          {
            for (site_t siteIndex = firstIndex; siteIndex < (firstIndex + siteCount); siteIndex++)
            {
              geometry::Site<geometry::LatticeData> site = latDat->GetSite(siteIndex);

              const distribn_t* fOld = site.GetFOld<LatticeType> ();

              kernels::HydroVars<typename CollisionType::CKernel> hydroVars(fOld);

              ///< @todo #126 This value of tau will be updated by some kernels within the collider code (e.g. LBGKNN). It would be nicer if tau is handled in a single place.
              hydroVars.tau = lbmParams->GetTau();

              collider.CalculatePreCollision(hydroVars, coupledPropertyCache, site);

              collider.Collide(lbmParams, hydroVars);

              for (Direction ii = 0; ii < LatticeType::NUMVECTORS; ii++)
              {
                if (site.HasVesselWall(ii))
                {
                  wallLinkDelegate.StreamLink(lbmParams, latDat, site, hydroVars, ii);
                }
                else
                {
                  bulkLinkDelegate.StreamLink(lbmParams, latDat, site, hydroVars, ii);
                }
              }

              //TODO: Necessary to specify sub-class?
              AdvectionDiffusionBaseStreamer<VesselWallStreamerTypeFactory>::template UpdateMinsAndMaxes<tDoRayTracing>(site,
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
                if (site.HasVesselWall(direction))
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
      class AdvectionDiffusionIoletStreamerTypeFactory : public AdvectionDiffusionBaseStreamer<AdvectionDiffusionIoletStreamerTypeFactory<CollisionImpl, IoletLinkImpl> >
      {
        public:
          typedef CollisionImpl CollisionType;

        private:
          CollisionType collider;
          SimpleCollideAndStreamDelegate<CollisionType> bulkLinkDelegate;
          IoletLinkImpl ioletLinkDelegate;

          typedef typename CollisionType::CKernel::LatticeType LatticeType;

        public:
          AdvectionDiffusionIoletStreamerTypeFactory(kernels::InitParams& initParams) :
            collider(initParams), bulkLinkDelegate(collider, initParams), ioletLinkDelegate(collider, initParams)
          {

          }

          template<bool tDoRayTracing>
          inline void DoStreamAndCollide(const site_t firstIndex,
                                         const site_t siteCount,
                                         const LbmParameters* lbmParams,
                                         geometry::LatticeData* latDat,
                                         lb::MacroscopicPropertyCache& propertyCache,
                                         lb::MacroscopicPropertyCache& coupledPropertyCache)
          {
            for (site_t siteIndex = firstIndex; siteIndex < (firstIndex + siteCount); siteIndex++)
            {
              geometry::Site<geometry::LatticeData> site = latDat->GetSite(siteIndex);

              const distribn_t* fOld = site.GetFOld<LatticeType> ();

              kernels::HydroVars<typename CollisionType::CKernel> hydroVars(fOld);

              ///< @todo #126 This value of tau will be updated by some kernels within the collider code (e.g. LBGKNN). It would be nicer if tau is handled in a single place.
              hydroVars.tau = lbmParams->GetTau();

              collider.CalculatePreCollision(hydroVars, coupledPropertyCache, site);

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
              AdvectionDiffusionBaseStreamer<AdvectionDiffusionIoletStreamerTypeFactory>::template UpdateMinsAndMaxes<tDoRayTracing>(site,
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
      class VesselWallIoletStreamerTypeFactory : public AdvectionDiffusionBaseStreamer<VesselWallIoletStreamerTypeFactory<CollisionImpl,
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
          VesselWallIoletStreamerTypeFactory(kernels::InitParams& initParams) :
            collider(initParams), bulkLinkDelegate(collider, initParams), wallLinkDelegate(collider, initParams),
                ioletLinkDelegate(collider, initParams)
          {

          }

          template<bool tDoRayTracing>
          inline void DoStreamAndCollide(const site_t firstIndex,
                                         const site_t siteCount,
                                         const LbmParameters* lbmParams,
                                         geometry::LatticeData* latDat,
                                         lb::MacroscopicPropertyCache& propertyCache,
                                         lb::MacroscopicPropertyCache& coupledPropertyCache)
          {
            for (site_t siteIndex = firstIndex; siteIndex < (firstIndex + siteCount); siteIndex++)
            {
              geometry::Site<geometry::LatticeData> site = latDat->GetSite(siteIndex);

              const distribn_t* fOld = site.GetFOld<LatticeType> ();

              kernels::HydroVars<typename CollisionType::CKernel> hydroVars(fOld);

              ///< @todo #126 This value of tau will be updated by some kernels within the collider code (e.g. LBGKNN). It would be nicer if tau is handled in a single place.
              hydroVars.tau = lbmParams->GetTau();

              collider.CalculatePreCollision(hydroVars, coupledPropertyCache, site);

              collider.Collide(lbmParams, hydroVars);

              for (Direction ii = 0; ii < LatticeType::NUMVECTORS; ii++)
              {
                if (site.HasIolet(ii))
                {
                  ioletLinkDelegate.StreamLink(lbmParams, latDat, site, hydroVars, ii);
                }
                else if (site.HasVesselWall(ii))
                {
                  wallLinkDelegate.StreamLink(lbmParams, latDat, site, hydroVars, ii);
                }
                else
                {
                  bulkLinkDelegate.StreamLink(lbmParams, latDat, site, hydroVars, ii);
                }
              }

              //TODO: Necessary to specify sub-class?
              AdvectionDiffusionBaseStreamer<VesselWallIoletStreamerTypeFactory>::template UpdateMinsAndMaxes<tDoRayTracing>(site,
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
                if (site.HasVesselWall(direction))
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

      template<typename CollisionImpl, typename WallLinkImpl>
      class StentWallStreamerTypeFactory : public AdvectionDiffusionBaseStreamer<StentWallStreamerTypeFactory<CollisionImpl, WallLinkImpl> >
      {
        public:
          typedef CollisionImpl CollisionType;

        private:
          CollisionType collider;
          SimpleCollideAndStreamDelegate<CollisionType> bulkLinkDelegate;
          WallLinkImpl wallLinkDelegate;

          typedef typename CollisionType::CKernel::LatticeType LatticeType;

        public:
          StentWallStreamerTypeFactory(kernels::InitParams& initParams) :
            collider(initParams), bulkLinkDelegate(collider, initParams), wallLinkDelegate(collider, initParams)
          {

          }

          template<bool tDoRayTracing>
          inline void DoStreamAndCollide(const site_t firstIndex,
                                         const site_t siteCount,
                                         const LbmParameters* lbmParams,
                                         geometry::LatticeData* latDat,
                                         lb::MacroscopicPropertyCache& propertyCache,
                                         lb::MacroscopicPropertyCache& coupledPropertyCache)
          {
            for (site_t siteIndex = firstIndex; siteIndex < (firstIndex + siteCount); siteIndex++)
            {
              geometry::Site<geometry::LatticeData> site = latDat->GetSite(siteIndex);

              const distribn_t* fOld = site.GetFOld<LatticeType> ();

              kernels::HydroVars<typename CollisionType::CKernel> hydroVars(fOld);

              ///< @todo #126 This value of tau will be updated by some kernels within the collider code (e.g. LBGKNN). It would be nicer if tau is handled in a single place.
              hydroVars.tau = lbmParams->GetTau();

              collider.CalculatePreCollision(hydroVars, coupledPropertyCache, site);

              collider.Collide(lbmParams, hydroVars);

              for (Direction ii = 0; ii < LatticeType::NUMVECTORS; ii++)
              {
                if (site.HasStentWall(ii))
                {
                  wallLinkDelegate.StreamLink(lbmParams, latDat, site, hydroVars, ii);
                }
                else
                {
                  bulkLinkDelegate.StreamLink(lbmParams, latDat, site, hydroVars, ii);
                }
              }

              //TODO: Necessary to specify sub-class?
              AdvectionDiffusionBaseStreamer<StentWallStreamerTypeFactory>::template UpdateMinsAndMaxes<tDoRayTracing>(site,
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
                if (site.HasStentWall(direction))
                {
                  wallLinkDelegate.PostStepLink(latticeData, site, direction);
                }
              }
            }
          }

      };

      template<typename CollisionImpl, typename WallLinkImpl, typename IoletLinkImpl>
      class StentWallIoletStreamerTypeFactory : public AdvectionDiffusionBaseStreamer<StentWallIoletStreamerTypeFactory<CollisionImpl,
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
          StentWallIoletStreamerTypeFactory(kernels::InitParams& initParams) :
            collider(initParams), bulkLinkDelegate(collider, initParams), wallLinkDelegate(collider, initParams),
                ioletLinkDelegate(collider, initParams)
          {

          }

          template<bool tDoRayTracing>
          inline void DoStreamAndCollide(const site_t firstIndex,
                                         const site_t siteCount,
                                         const LbmParameters* lbmParams,
                                         geometry::LatticeData* latDat,
                                         lb::MacroscopicPropertyCache& propertyCache,
                                         lb::MacroscopicPropertyCache& coupledPropertyCache)
          {
            for (site_t siteIndex = firstIndex; siteIndex < (firstIndex + siteCount); siteIndex++)
            {
              geometry::Site<geometry::LatticeData> site = latDat->GetSite(siteIndex);

              const distribn_t* fOld = site.GetFOld<LatticeType> ();

              kernels::HydroVars<typename CollisionType::CKernel> hydroVars(fOld);

              ///< @todo #126 This value of tau will be updated by some kernels within the collider code (e.g. LBGKNN). It would be nicer if tau is handled in a single place.
              hydroVars.tau = lbmParams->GetTau();

              collider.CalculatePreCollision(hydroVars, coupledPropertyCache, site);

              collider.Collide(lbmParams, hydroVars);

              for (Direction ii = 0; ii < LatticeType::NUMVECTORS; ii++)
              {
                if (site.HasIolet(ii))
                {
                  ioletLinkDelegate.StreamLink(lbmParams, latDat, site, hydroVars, ii);
                }
                else if (site.HasStentWall(ii))
                {
                  wallLinkDelegate.StreamLink(lbmParams, latDat, site, hydroVars, ii);
                }
                else
                {
                  bulkLinkDelegate.StreamLink(lbmParams, latDat, site, hydroVars, ii);
                }
              }

              //TODO: Necessary to specify sub-class?
              AdvectionDiffusionBaseStreamer<StentWallIoletStreamerTypeFactory>::template UpdateMinsAndMaxes<tDoRayTracing>(site,
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
                if (site.HasStentWall(direction))
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

      template<typename CollisionImpl, typename WallLinkImpl, typename VesselWallLinkImpl>
      class StentWallVesselWallStreamerTypeFactory : public AdvectionDiffusionBaseStreamer<StentWallVesselWallStreamerTypeFactory<CollisionImpl,
          WallLinkImpl, VesselWallLinkImpl> >
      {
        public:
          typedef CollisionImpl CollisionType;
          typedef typename CollisionType::CKernel::LatticeType LatticeType;

        private:
          CollisionType collider;
          SimpleCollideAndStreamDelegate<CollisionType> bulkLinkDelegate;
          WallLinkImpl wallLinkDelegate;
          VesselWallLinkImpl vesselWallLinkDelegate;

        public:
          StentWallVesselWallStreamerTypeFactory(kernels::InitParams& initParams) :
            collider(initParams), bulkLinkDelegate(collider, initParams), wallLinkDelegate(collider, initParams),
                vesselWallLinkDelegate(collider, initParams)
          {

          }

          template<bool tDoRayTracing>
          inline void DoStreamAndCollide(const site_t firstIndex,
                                         const site_t siteCount,
                                         const LbmParameters* lbmParams,
                                         geometry::LatticeData* latDat,
                                         lb::MacroscopicPropertyCache& propertyCache,
                                         lb::MacroscopicPropertyCache& coupledPropertyCache)
          {
            for (site_t siteIndex = firstIndex; siteIndex < (firstIndex + siteCount); siteIndex++)
            {
              geometry::Site<geometry::LatticeData> site = latDat->GetSite(siteIndex);

              const distribn_t* fOld = site.GetFOld<LatticeType> ();

              kernels::HydroVars<typename CollisionType::CKernel> hydroVars(fOld);

              ///< @todo #126 This value of tau will be updated by some kernels within the collider code (e.g. LBGKNN). It would be nicer if tau is handled in a single place.
              hydroVars.tau = lbmParams->GetTau();

              collider.CalculatePreCollision(hydroVars, coupledPropertyCache, site);

              collider.Collide(lbmParams, hydroVars);

              for (Direction ii = 0; ii < LatticeType::NUMVECTORS; ii++)
              {
                if (site.HasVesselWall(ii))
                {
                  vesselWallLinkDelegate.StreamLink(lbmParams, latDat, site, hydroVars, ii);
                }
                else if (site.HasStentWall(ii))
                {
                  wallLinkDelegate.StreamLink(lbmParams, latDat, site, hydroVars, ii);
                }
                else
                {
                  bulkLinkDelegate.StreamLink(lbmParams, latDat, site, hydroVars, ii);
                }
              }

              //TODO: Necessary to specify sub-class?
              AdvectionDiffusionBaseStreamer<StentWallVesselWallStreamerTypeFactory>::template UpdateMinsAndMaxes<tDoRayTracing>(site,
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
                if (site.HasStentWall(direction))
                {
                  wallLinkDelegate.PostStepLink(latticeData, site, direction);
                }
                else if (site.HasVesselWall(direction))
                {
                  vesselWallLinkDelegate.PostStepLink(latticeData, site, direction);
                }
              }
            }

          }
      };

      template<typename CollisionImpl, typename WallLinkImpl, typename VesselWallLinkImpl, typename IoletLinkImpl>
      class StentWallVesselWallIoletStreamerTypeFactory : public AdvectionDiffusionBaseStreamer<StentWallVesselWallIoletStreamerTypeFactory<CollisionImpl,
          WallLinkImpl, VesselWallLinkImpl, IoletLinkImpl> >
      {
        public:
          typedef CollisionImpl CollisionType;
          typedef typename CollisionType::CKernel::LatticeType LatticeType;

        private:
          CollisionType collider;
          SimpleCollideAndStreamDelegate<CollisionType> bulkLinkDelegate;
          WallLinkImpl wallLinkDelegate;
          VesselWallLinkImpl vesselWallLinkDelegate;
          IoletLinkImpl ioletLinkDelegate;

        public:
          StentWallVesselWallIoletStreamerTypeFactory(kernels::InitParams& initParams) :
            collider(initParams), bulkLinkDelegate(collider, initParams), wallLinkDelegate(collider, initParams),
                vesselWallLinkDelegate(collider, initParams), ioletLinkDelegate(collider, initParams)
          {

          }

          template<bool tDoRayTracing>
          inline void DoStreamAndCollide(const site_t firstIndex,
                                         const site_t siteCount,
                                         const LbmParameters* lbmParams,
                                         geometry::LatticeData* latDat,
                                         lb::MacroscopicPropertyCache& propertyCache,
                                         lb::MacroscopicPropertyCache& coupledPropertyCache)
          {
            for (site_t siteIndex = firstIndex; siteIndex < (firstIndex + siteCount); siteIndex++)
            {
              geometry::Site<geometry::LatticeData> site = latDat->GetSite(siteIndex);

              const distribn_t* fOld = site.GetFOld<LatticeType> ();

              kernels::HydroVars<typename CollisionType::CKernel> hydroVars(fOld);

              ///< @todo #126 This value of tau will be updated by some kernels within the collider code (e.g. LBGKNN). It would be nicer if tau is handled in a single place.
              hydroVars.tau = lbmParams->GetTau();

              collider.CalculatePreCollision(hydroVars, coupledPropertyCache, site);

              collider.Collide(lbmParams, hydroVars);

              for (Direction ii = 0; ii < LatticeType::NUMVECTORS; ii++)
              {
                if (site.HasVesselWall(ii))
                {
                  vesselWallLinkDelegate.StreamLink(lbmParams, latDat, site, hydroVars, ii);
                }
                else if (site.HasStentWall(ii))
                {
                  wallLinkDelegate.StreamLink(lbmParams, latDat, site, hydroVars, ii);
                }
                else if (site.HasIolet(ii))
                {
                  ioletLinkDelegate.StreamLink(lbmParams, latDat, site, hydroVars, ii);
                }
                else
                {
                  bulkLinkDelegate.StreamLink(lbmParams, latDat, site, hydroVars, ii);
                }
              }

              //TODO: Necessary to specify sub-class?
              AdvectionDiffusionBaseStreamer<StentWallVesselWallIoletStreamerTypeFactory>::template UpdateMinsAndMaxes<tDoRayTracing>(site,
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
                if (site.HasStentWall(direction))
                {
                  wallLinkDelegate.PostStepLink(latticeData, site, direction);
                }
                else if (site.HasVesselWall(direction))
                {
                  vesselWallLinkDelegate.PostStepLink(latticeData, site, direction);
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
#endif // HEMELB_LB_STREAMERS_ADVECTIONDIFFUSIONSTREAMERTYPEFACTORY_H
