
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_STREAMERS_SIMPLECOLLIDEANDSTREAMDELEGATE_H
#define HEMELB_LB_STREAMERS_SIMPLECOLLIDEANDSTREAMDELEGATE_H

namespace hemelb
{
  namespace lb
  {
    namespace streamers
    {

      template<typename CollisionImpl>
      class SimpleCollideAndStreamDelegate : public BaseStreamerDelegate<CollisionImpl>
      {
        public:
          typedef CollisionImpl CollisionType;
          typedef typename CollisionType::CKernel::LatticeType LatticeType;

          SimpleCollideAndStreamDelegate(CollisionType& delegatorCollider, kernels::InitParams& initParams)
          {
          }

          inline void StreamLink(const LbmParameters* lbmParams,
                                 geometry::LatticeData* const latticeData,
                                 const geometry::Site<geometry::LatticeData>& site,
                                 kernels::HydroVars<typename CollisionType::CKernel>& hydroVars,
                                 const Direction& direction)
          {
            * (latticeData->GetFNew(site.GetStreamedIndex<LatticeType> (direction)))
                = hydroVars.GetFPostCollision()[direction];
          }

      };

    }
  }
}

#endif /* HEMELB_LB_STREAMERS_SIMPLECOLLIDEANDSTREAMDELEGATE_H */
