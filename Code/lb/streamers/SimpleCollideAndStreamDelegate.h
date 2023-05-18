// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_STREAMERS_SIMPLECOLLIDEANDSTREAMDELEGATE_H
#define HEMELB_LB_STREAMERS_SIMPLECOLLIDEANDSTREAMDELEGATE_H

#include "geometry/FieldData.h"
#include "lb/streamers/LinkStreamer.h"

namespace hemelb::lb::streamers
{

    /// Do bulk streaming along a link
    /// Instantiations satisfy LinkStreamer
    template<typename CollisionImpl>
    class SimpleCollideAndStreamDelegate
    {
    public:
        using CollisionType = CollisionImpl;
        using LatticeType = typename CollisionType::CKernel::LatticeType;

        SimpleCollideAndStreamDelegate(CollisionType& delegatorCollider,
                                       InitParams& initParams)
        {
        }

        void StreamLink(const LbmParameters* lbmParams,
                        geometry::FieldData& latticeData,
                        const geometry::Site<geometry::FieldData>& site,
                        HydroVars<typename CollisionType::CKernel>& hydroVars,
                        const Direction& direction)
        {
            * (latticeData.GetFNew(site.GetStreamedIndex<LatticeType>(direction))) =
                    hydroVars.GetFPostCollision()[direction];
        }

        void PostStepLink(geometry::FieldData& latticeData,
                          const geometry::Site<geometry::FieldData>& site,
                          const Direction& direction) {

        }
    };

}

#endif /* HEMELB_LB_STREAMERS_SIMPLECOLLIDEANDSTREAMDELEGATE_H */
