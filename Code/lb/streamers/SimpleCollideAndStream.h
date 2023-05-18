// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_STREAMERS_SIMPLECOLLIDEANDSTREAM_H
#define HEMELB_LB_STREAMERS_SIMPLECOLLIDEANDSTREAM_H

#include "lb/streamers/BaseStreamer.h"
#include "lb/streamers/SimpleCollideAndStreamDelegate.h"
#include "lb/HFunction.h"

namespace hemelb::lb::streamers
{

    template<typename CollisionImpl>
    class SimpleCollideAndStream : public BaseStreamer
    {
    public:
        using CollisionType = CollisionImpl;

    private:
        CollisionType collider;
        SimpleCollideAndStreamDelegate<CollisionType> bulkLinkDelegate;
        static_assert(LinkStreamer<SimpleCollideAndStreamDelegate<CollisionType>>);
        using LatticeType = typename CollisionType::CKernel::LatticeType;

    public:
        SimpleCollideAndStream(InitParams& initParams) :
                collider(initParams), bulkLinkDelegate(collider, initParams)
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

                for (unsigned int ii = 0; ii < LatticeType::NUMVECTORS; ii++)
                {
                    bulkLinkDelegate.StreamLink(lbmParams, latDat, site, hydroVars, ii);
                }

                UpdateMinsAndMaxes(site, hydroVars, lbmParams, propertyCache);
            }
        }

        void PostStep(const site_t iFirstIndex, const site_t iSiteCount,
                      const LbmParameters* iLbmParams, geometry::FieldData& bLatDat,
                      lb::MacroscopicPropertyCache& propertyCache)
        {
        }

    };
}
#endif
