// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_LB_STREAMERS_LADDIOLET_H
#define HEMELB_LB_STREAMERS_LADDIOLET_H

#include "lb/kernels/BaseKernel.h"
#include "lb/streamers/BaseStreamer.h"
#include "lb/streamers/LaddIoletDelegate.h"
#include "lb/streamers/SimpleCollideAndStreamDelegate.h"

namespace hemelb
{
  namespace lb
  {
    namespace streamers
    {

      template<typename CollisionImpl>
      class LaddIolet : public BaseStreamer<LaddIolet<CollisionImpl> >
      {
        public:
          typedef CollisionImpl CollisionType;

        private:
          CollisionType collider;
          class SimpleCollideAndStreamDelegate<CollisionType> bulkLinkDelegate;
          class LaddIoletDelegate<CollisionType> ioletLinkDelegate;

          typedef typename CollisionType::CKernel::LatticeType LatticeType;

        public:
          LaddIolet(kernels::InitParams& initParams) :
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
                  ioletLinkDelegate.StreamLink(latDat, site, hydroVars, ii);
                }
                else
                {
                  bulkLinkDelegate.StreamLink(latDat, site, hydroVars, ii);
                }
              }

              //TODO: Necessary to specify sub-class?
              BaseStreamer<LaddIolet>::template UpdateMinsAndMaxes<tDoRayTracing>(site,
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
      };

      template<typename CollisionImpl, typename WallLinkImpl>
      class LaddIoletWall : public BaseStreamer<LaddIoletWall<CollisionImpl, WallLinkImpl> >
      {
        public:
          typedef CollisionImpl CollisionType;

        private:
          CollisionType collider;
          SimpleCollideAndStreamDelegate<CollisionType> bulkLinkDelegate;
          WallLinkImpl wallLinkDelegate;
          LaddIoletDelegate<CollisionType> ioletLinkDelegate;

          typedef typename CollisionType::CKernel::LatticeType LatticeType;

        public:
          LaddIoletWall(kernels::InitParams& initParams) :
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
                  ioletLinkDelegate.StreamLink(latDat, site, hydroVars, ii);
                }
                else if (site.HasBoundary(ii))
                {
                  wallLinkDelegate.StreamLink(latDat, site, hydroVars, ii);
                }
                else
                {
                  bulkLinkDelegate.StreamLink(latDat, site, hydroVars, ii);
                }
              }

              //TODO: Necessary to specify sub-class?
              BaseStreamer<LaddIoletWall>::template UpdateMinsAndMaxes<tDoRayTracing>(site,
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
      };

      template<class T>
      struct LaddIoletSBB
      {
          typedef LaddIoletWall<T, SimpleBounceBackDelegate<T> > Type;
      };

      template<class T>
      struct LaddIoletBFL
      {
          typedef LaddIoletWall<T, BouzidiFirdaousLallemandDelegate<T> > Type;
      };

    }
  }
}

#endif /* HEMELB_LB_STREAMERS_LADDIOLET_H */
