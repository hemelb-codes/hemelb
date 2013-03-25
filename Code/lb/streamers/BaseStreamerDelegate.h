// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_LB_STREAMERS_BASESTREAMERDELEGATE_H
#define HEMELB_LB_STREAMERS_BASESTREAMERDELEGATE_H

#include "geometry/LatticeData.h"
#include "lb/kernels/BaseKernel.h"

namespace hemelb
{
  namespace lb
  {
    namespace streamers
    {
      /**
       * Base class for Streamer delegates.
       *
       * Sets out the interface that streamers should implement and typedefs
       * that they should make.
       *
       * Unfortunately C++ template rules cause a difficulty: since the base
       * class is a dependent name (i.e. only know on instantiation) and
       * unqualified names are assumed to be non-dependent we cannot use
       * CollisionType and LatticeType in subclasses without using the (rather
       * lond) name of the appropriate instantiation of this template.
       *
       * Subclasses should therefore redeclare CollisionType and LatticeType
       * for their convenience.
       */
      template<class CollisionImpl>
      class BaseStreamerDelegate
      {
        public:
          typedef CollisionImpl CollisionType;
          typedef typename CollisionType::CKernel::LatticeType LatticeType;

          /**
           * Constructor
           * @param delegatorCollider
           * @param initParams
           * @return
           */
          //BaseStreamerDelegate(CollisionType& delegatorCollider, kernels::InitParams& initParams);
        protected:
          /**
           * Protected default ctor to make life a bit easier for subclasses.
           */
          BaseStreamerDelegate()
          {
          }
        public:
          /**
           * Perform the streaming operation from site along direction
           *
           * hydroVars must be post-collision
           *
           * @param latticeData
           * @param site
           * @param hydroVars
           * @param direction
           */
          inline void StreamLink(const LbmParameters* lbmParams,
                                 geometry::LatticeData* const latticeData,
                                 const geometry::Site<geometry::LatticeData>& site,
                                 kernels::HydroVars<typename CollisionType::CKernel>& hydroVars,
                                 const Direction& direction)
          {
          }

          /**
           * Perform any post-step operations for the link from site along direction
           *
           * @param latticeData
           * @param site
           * @param direction
           */
          inline void PostStepLink(geometry::LatticeData* const latticeData,
                                   const geometry::Site<geometry::LatticeData>& site,
                                   const Direction& direction)
          {
          }
      };

    }
  }
}

#endif /* HEMELB_LB_STREAMERS_BASESTREAMERDELEGATE_H */
