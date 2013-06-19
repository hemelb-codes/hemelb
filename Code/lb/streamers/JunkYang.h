// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_LB_STREAMERS_JUNKYANG_H
#define HEMELB_LB_STREAMERS_JUNKYANG_H

#include "lb/streamers/JunkYangFactory.h"

namespace hemelb
{
  namespace lb
  {
    namespace streamers
    {
      /**
       * Null implementation of an iolet link delegate.
       */
      template<typename CollisionImpl>
      struct NoIoletLink : BaseStreamerDelegate<CollisionImpl>
      {
          NoIoletLink(CollisionImpl& collider, kernels::InitParams& initParams)
          {
          }
      };

      /**
       * Metafunction returning a streamer that does JY at walls only.
       */
      template<typename CollisionImpl>
      struct JunkYang
      {
          typedef JunkYangFactory<CollisionImpl, NoIoletLink<CollisionImpl> > Type;
      };

    }
  }
}

#endif /* HEMELB_LB_STREAMERS_JUNKYANG_H */
