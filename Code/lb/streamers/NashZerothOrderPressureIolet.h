#ifndef HEMELB_LB_STREAMERS_NASHZEROTHORDERPRESSUREIOLET_H
#define HEMELB_LB_STREAMERS_NASHZEROTHORDERPRESSUREIOLET_H

#include "lb/streamers/StreamerTypeFactory.h"
#include "lb/streamers/NashZerothOrderPressureDelegate.h"
#include "lb/streamers/SimpleBounceBackDelegate.h"
#include "lb/streamers/BouzidiFirdaousLallemandDelegate.h"

namespace hemelb
{
  namespace lb
  {
    namespace streamers
    {

      template<class T>
      struct NashZerothOrderPressureIolet
      {
          typedef IoletStreamerTypeFactory<T, NashZerothOrderPressureDelegate<T> > Type;
      };

      template<class T>
      struct NashZerothOrderPressureIoletSBB
      {
          typedef WallIoletStreamerTypeFactory<T, SimpleBounceBackDelegate<T> , NashZerothOrderPressureDelegate<T> >
              Type;
      };

      template<class T>
      struct NashZerothOrderPressureIoletBFL
      {
          typedef WallIoletStreamerTypeFactory<T, BouzidiFirdaousLallemandDelegate<T> ,
              NashZerothOrderPressureDelegate<T> > Type;
      };
    }
  }
}

#endif /* HEMELB_LB_STREAMERS_NASHZEROTHORDERPRESSUREIOLET_H */
