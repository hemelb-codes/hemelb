
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.
#ifndef HEMELB_LB_STREAMERS_NASHZEROTHORDERPRESSUREIOLET_H
#define HEMELB_LB_STREAMERS_NASHZEROTHORDERPRESSUREIOLET_H

#include "lb/streamers/StreamerTypeFactory.h"
#include "lb/streamers/NashZerothOrderPressureDelegate.h"
#include "lb/streamers/SimpleBounceBackDelegate.h"
#include "lb/streamers/BouzidiFirdaousLallemandDelegate.h"
#include "lb/streamers/JunkYangFactory.h"

namespace hemelb
{
  namespace lb
  {
    namespace streamers
    {

      template<class CollisionType>
      struct NashZerothOrderPressureIoletSBB
      {
          typedef StreamerTypeFactory<CollisionType, SimpleBounceBackDelegate<CollisionType> ,
              NashZerothOrderPressureDelegate<CollisionType> > Type;
      };

      template<class CollisionType>
      struct NashZerothOrderPressureIoletBFL
      {
          typedef StreamerTypeFactory<CollisionType, BouzidiFirdaousLallemandDelegate<CollisionType> ,
              NashZerothOrderPressureDelegate<CollisionType> > Type;
      };

      template<class CollisionType>
      struct NashZerothOrderPressureIoletGZS
      {
          typedef StreamerTypeFactory<CollisionType, GuoZhengShiDelegate<CollisionType> ,
              NashZerothOrderPressureDelegate<CollisionType> > Type;
      };

      template<class CollisionType>
      struct NashZerothOrderPressureIoletJY
      {
          typedef JunkYangFactory<CollisionType, NashZerothOrderPressureDelegate<CollisionType> > Type;
      };

    }
  }
}

#endif /* HEMELB_LB_STREAMERS_NASHZEROTHORDERPRESSUREIOLET_H */
