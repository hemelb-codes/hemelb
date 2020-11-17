
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

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
