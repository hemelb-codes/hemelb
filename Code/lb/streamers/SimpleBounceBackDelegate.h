// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_STREAMERS_SIMPLEBOUNCEBACKDELEGATE_H
#define HEMELB_LB_STREAMERS_SIMPLEBOUNCEBACKDELEGATE_H

#include "lb/streamers/LinkStreamer.h"
#include "lb/streamers/SimpleCollideAndStream.h"

namespace hemelb::lb::streamers
{

    template<typename CollisionImpl>
    class SimpleBounceBackDelegate
    {
    public:
        using CollisionType = CollisionImpl;
        using LatticeType = typename CollisionType::CKernel::LatticeType;

        static site_t GetBBIndex(site_t siteIndex, int direction)
        {
            return (siteIndex * LatticeType::NUMVECTORS) + LatticeType::INVERSEDIRECTIONS[direction];
        }

        SimpleBounceBackDelegate(CollisionType& delegatorCollider,
                                 InitParams& initParams)
        {
        }

        void StreamLink(const LbmParameters* lbmParams,
                        geometry::FieldData& latticeData,
                        const geometry::Site<geometry::Domain>& site,
                        HydroVars<typename CollisionType::CKernel>& hydroVars,
                        const Direction& direction)
        {
            // Propagate the outgoing post-collisional f into the opposite direction.
            * (latticeData.GetFNew(GetBBIndex(site.GetIndex(), direction))) =
                    hydroVars.GetFPostCollision()[direction];
        }
        void PostStepLink(geometry::FieldData& latticeData,
                          const geometry::Site<geometry::FieldData>& site,
                          const Direction& direction)
        {
            // Nothing to do :)
        }
    };


}

#endif /* HEMELB_LB_STREAMERS_SIMPLEBOUNCEBACKDELEGATE_H */
