// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_LB_STREAMERS_SIMPLEBOUNCEBACK_H
#define HEMELB_LB_STREAMERS_SIMPLEBOUNCEBACK_H

#include "lb/streamers/SimpleBounceBackDelegate.h"
#include "lb/streamers/StreamerTypeFactory.h"

namespace hemelb
{
  namespace lb
  {
    namespace streamers
    {

      template<typename CollisionImpl>
      struct SimpleBounceBack
      {
          typedef WallStreamerTypeFactory<CollisionImpl, SimpleBounceBackDelegate<CollisionImpl> > Type;
      };

    }
  }
}

#endif /* HEMELB_LB_STREAMERS_SIMPLEBOUNCEBACK_H */
