// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_LB_STREAMERS_LADDIOLET_H
#define HEMELB_LB_STREAMERS_LADDIOLET_H

#include "lb/streamers/StreamerTypeFactory.h"
#include "lb/streamers/LaddIoletDelegate.h"
#include "lb/streamers/SimpleBounceBackDelegate.h"

namespace hemelb
{
  namespace lb
  {
    namespace streamers
    {

      template<class T>
      struct LaddIolet
      {
          typedef IoletStreamerTypeFactory<T, LaddIoletDelegate<T> > Type;
      };

      template<class T>
      struct LaddIoletSBB
      {
          typedef WallIoletStreamerTypeFactory<T, SimpleBounceBackDelegate<T> , LaddIoletDelegate<T> > Type;
      };

      template<class T>
      struct LaddIoletBFL
      {
          typedef WallIoletStreamerTypeFactory<T, BouzidiFirdaousLallemandDelegate<T> , LaddIoletDelegate<T> > Type;
      };

    }
  }
}

#endif /* HEMELB_LB_STREAMERS_LADDIOLET_H */
