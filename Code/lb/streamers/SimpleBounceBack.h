// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_STREAMERS_SIMPLEBOUNCEBACK_H
#define HEMELB_LB_STREAMERS_SIMPLEBOUNCEBACK_H

#include "lb/concepts.h"
#include "lb/streamers/BulkStreamer.h"

namespace hemelb::lb
{

    template<collision_type C>
    class BounceBackLink
    {
    public:
        using CollisionType = C;
        using VarsType = typename CollisionType::VarsType;
        using LatticeType = typename CollisionType::LatticeType;

        static site_t GetBBIndex(site_t siteIndex, int direction)
        {
            return (siteIndex * LatticeType::NUMVECTORS) + LatticeType::INVERSEDIRECTIONS[direction];
        }

        BounceBackLink(CollisionType& delegatorCollider,
                       InitParams& initParams)
        {
        }

        void StreamLink(const LbmParameters* lbmParams,
                        geometry::FieldData& latticeData,
                        const geometry::Site<geometry::Domain>& site,
                        VarsType& hydroVars,
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

#endif
