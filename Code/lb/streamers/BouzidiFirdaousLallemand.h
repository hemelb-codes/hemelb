// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_LB_STREAMERS_BOUZIDIFIRDAOUSLALLEMAND_H
#define HEMELB_LB_STREAMERS_BOUZIDIFIRDAOUSLALLEMAND_H

#include "lb/streamers/BouzidiFirdaousLallemandDelegate.h"
#include "lb/streamers/StreamerTypeFactory.h"

namespace hemelb
{
  namespace lb
  {
    namespace streamers
    {

      template<typename CollisionImpl>
      struct BouzidiFirdaousLallemand
      {
          typedef WallStreamerTypeFactory<CollisionImpl, BouzidiFirdaousLallemandDelegate<CollisionImpl> > Type;
      };

    }
  }
}

#endif /* HEMELB_LB_STREAMERS_BOUZIDIFIRDAOUSLALLEMAND_H */
