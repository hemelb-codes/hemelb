// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_LB_STREAMERS_SIMPLECOLLIDEANDSTREAMDELEGATE_H
#define HEMELB_LB_STREAMERS_SIMPLECOLLIDEANDSTREAMDELEGATE_H

namespace hemelb
{
  namespace lb
  {
    namespace streamers
    {

      template<typename CollisionImpl>
      class SimpleCollideAndStreamDelegate
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
