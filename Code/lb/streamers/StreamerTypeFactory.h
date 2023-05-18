// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_STREAMERS_STREAMERTYPEFACTORY_H
#define HEMELB_LB_STREAMERS_STREAMERTYPEFACTORY_H

#include "lb/streamers/BaseStreamer.h"
#include "lb/streamers/SimpleCollideAndStreamDelegate.h"

namespace hemelb::lb::streamers
{
      /**
       * Template to produce Streamers that can cope with fluid-fluid and
       * fluid-wall links. Requires two classes as arguments: 1) the Collision
       * class and 2) a StreamerDelegate class that will handle the wall links.
       *
       * It is intended that a simpler metafunction partially specialise this
       * template on WallLinkImpl.
       */
      template<typename CollisionImpl, LinkStreamer WallLinkImpl>
      class WallStreamerTypeFactory : public BaseStreamer
      {
        public:
          using CollisionType = CollisionImpl;

        private:
          CollisionType collider;
          SimpleCollideAndStreamDelegate<CollisionType> bulkLinkDelegate;
          WallLinkImpl wallLinkDelegate;

          using LatticeType = typename CollisionType::CKernel::LatticeType;

        public:
          WallStreamerTypeFactory(InitParams& initParams) :
              collider(initParams), bulkLinkDelegate(collider, initParams),
                  wallLinkDelegate(collider, initParams)
          {

          }

          inline void StreamAndCollide(const site_t firstIndex, const site_t siteCount,
                                         const LbmParameters* lbmParams,
                                         geometry::FieldData& latDat,
                                         lb::MacroscopicPropertyCache& propertyCache)
          {
            for (site_t siteIndex = firstIndex; siteIndex < (firstIndex + siteCount); siteIndex++)
            {
              geometry::Site<geometry::FieldData> site = latDat.GetSite(siteIndex);
              HydroVars<typename CollisionType::CKernel> hydroVars(site);

              ///< @todo #126 This value of tau will be updated by some kernels
              //within the collider code (e.g. LBGKNN). It would be nicer if
              //tau is handled in a single place.
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
              UpdateMinsAndMaxes(site,
                                                                                                hydroVars,
                                                                                                lbmParams,
                                                                                                propertyCache);
            }
          }

          inline void PostStep(const site_t firstIndex, const site_t siteCount,
                                 const LbmParameters* lbmParameters,
                                 geometry::FieldData& latticeData,
                                 lb::MacroscopicPropertyCache& propertyCache)
          {
            for (site_t siteIndex = firstIndex; siteIndex < (firstIndex + siteCount); siteIndex++)
            {
              auto&& site = latticeData.GetSite(siteIndex);
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

      template <typename Collision, template<typename> typename WallLink>
      using WallStreamer = WallStreamerTypeFactory<Collision, WallLink<Collision>>;
      /**
       * Template to produce Streamers that can cope with fluid-fluid and
       * fluid-iolet links. Requires two classes as arguments: 1) the Collision
       * class and 2) a StreamerDelegate class that will handle the iolet links.
       *
       * It is intended that a simpler metafunction partially specialise this
       * template on IoletLinkImpl.
       */
      template<typename CollisionImpl, LinkStreamer IoletLinkImpl>
      class IoletStreamerTypeFactory : public BaseStreamer
      {
        public:
          using CollisionType = CollisionImpl;

        private:
          CollisionType collider;
          SimpleCollideAndStreamDelegate<CollisionType> bulkLinkDelegate;
          IoletLinkImpl ioletLinkDelegate;

          using LatticeType = typename CollisionType::CKernel::LatticeType;

        public:
          IoletStreamerTypeFactory(InitParams& initParams) :
              collider(initParams), bulkLinkDelegate(collider, initParams),
                  ioletLinkDelegate(collider, initParams)
          {

          }

          void StreamAndCollide(const site_t firstIndex, const site_t siteCount,
                                         const LbmParameters* lbmParams,
                                         geometry::FieldData& latDat,
                                         lb::MacroscopicPropertyCache& propertyCache)
          {
            for (site_t siteIndex = firstIndex; siteIndex < (firstIndex + siteCount); siteIndex++)
            {
              geometry::Site<geometry::FieldData> site = latDat.GetSite(siteIndex);
              HydroVars<typename CollisionType::CKernel> hydroVars(site);

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

              UpdateMinsAndMaxes(site,
                                                                                                 hydroVars,
                                                                                                 lbmParams,
                                                                                                 propertyCache);
            }
          }

          void PostStep(const site_t firstIndex, const site_t siteCount,
                                 const LbmParameters* lbmParameters,
                                 geometry::FieldData& latticeData,
                                 lb::MacroscopicPropertyCache& propertyCache)
          {
            for (site_t siteIndex = firstIndex; siteIndex < (firstIndex + siteCount); siteIndex++)
            {
              geometry::Site<geometry::FieldData> site = latticeData.GetSite(siteIndex);
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
      template<typename CollisionImpl, LinkStreamer WallLinkImpl, LinkStreamer IoletLinkImpl>
      class WallIoletStreamerTypeFactory : public BaseStreamer
      {
        public:
          using CollisionType = CollisionImpl;
          using LatticeType = typename CollisionType::CKernel::LatticeType;

        private:
          CollisionType collider;
          SimpleCollideAndStreamDelegate<CollisionType> bulkLinkDelegate;
          WallLinkImpl wallLinkDelegate;
          IoletLinkImpl ioletLinkDelegate;

        public:
          WallIoletStreamerTypeFactory(InitParams& initParams) :
              collider(initParams), bulkLinkDelegate(collider, initParams),
                  wallLinkDelegate(collider, initParams), ioletLinkDelegate(collider, initParams)
          {

          }

          void StreamAndCollide(const site_t firstIndex, const site_t siteCount,
                                         const LbmParameters* lbmParams,
                                         geometry::FieldData& latDat,
                                         lb::MacroscopicPropertyCache& propertyCache)
          {
            for (site_t siteIndex = firstIndex; siteIndex < (firstIndex + siteCount); siteIndex++)
            {
              geometry::Site<geometry::FieldData> site = latDat.GetSite(siteIndex);
              HydroVars<typename CollisionType::CKernel> hydroVars(site);

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
              UpdateMinsAndMaxes(site,
                                                                                                     hydroVars,
                                                                                                     lbmParams,
                                                                                                     propertyCache);
            }
          }

          void PostStep(const site_t firstIndex, const site_t siteCount,
                                 const LbmParameters* lbmParams, geometry::FieldData& latticeData,
                                 lb::MacroscopicPropertyCache& propertyCache)
          {
            for (site_t siteIndex = firstIndex; siteIndex < (firstIndex + siteCount); siteIndex++)
            {
              geometry::Site<geometry::FieldData> site = latticeData.GetSite(siteIndex);
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
#endif // HEMELB_LB_STREAMERS_STREAMERTYPEFACTORY_H
