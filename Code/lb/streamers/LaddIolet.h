
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_STREAMERS_LADDIOLET_H
#define HEMELB_LB_STREAMERS_LADDIOLET_H

#include "lb/streamers/StreamerTypeFactory.h"
#include "lb/streamers/LaddIoletDelegate.h"
#include "lb/streamers/SimpleBounceBackDelegate.h"
#include "lb/streamers/BouzidiFirdaousLallemandDelegate.h"
#include "lb/streamers/GuoZhengShiDelegate.h"
#include "lb/streamers/JunkYangFactory.h"

namespace hemelb
{
  namespace lb
  {
    namespace streamers
    {

      template<class CollisionType>
      struct LaddIolet
      {
          typedef IoletStreamerTypeFactory<CollisionType, LaddIoletDelegate<CollisionType> > Type;
      };

      template<class CollisionType>
      struct LaddIoletSBB
      {
          typedef WallIoletStreamerTypeFactory<CollisionType,
              SimpleBounceBackDelegate<CollisionType>, LaddIoletDelegate<CollisionType> > Type;
      };

      template<class CollisionType>
      struct LaddIoletBFL
      {
          typedef WallIoletStreamerTypeFactory<CollisionType,
              BouzidiFirdaousLallemandDelegate<CollisionType>, LaddIoletDelegate<CollisionType> > Type;
      };

      template<class CollisionType>
      struct LaddIoletGZS
      {
          typedef WallIoletStreamerTypeFactory<CollisionType, GuoZhengShiDelegate<CollisionType>,
              LaddIoletDelegate<CollisionType> > Type;
      };

      template<class CollisionType>
      struct LaddIoletJY
      {
          typedef JunkYangFactory<CollisionType, LaddIoletDelegate<CollisionType> > Type;
      };
    }
  }
}

#endif /* HEMELB_LB_STREAMERS_LADDIOLET_H */
