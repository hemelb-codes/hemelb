
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

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
