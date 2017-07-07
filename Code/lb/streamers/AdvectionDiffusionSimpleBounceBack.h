
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_STREAMERS_ADVECTIONDIFFUSIONSIMPLEBOUNCEBACK_H
#define HEMELB_LB_STREAMERS_ADVECTIONDIFFUSIONSIMPLEBOUNCEBACK_H

#include "lb/streamers/SimpleBounceBackDelegate.h"
#include "lb/streamers/AdvectionDiffusionStreamerTypeFactory.h"

namespace hemelb
{
  namespace lb
  {
    namespace streamers
    {

      template<typename CollisionImpl>
      struct AdvectionDiffusionSimpleBounceBack
      {
          typedef VesselWallStreamerTypeFactory<CollisionImpl, SimpleBounceBackDelegate<CollisionImpl> > Type;
      };

      template<typename CollisionImpl>
      struct AdvectionDiffusionSimpleBounceBackIolet
      {
          typedef AdvectionDiffusionIoletStreamerTypeFactory<CollisionImpl, SimpleBounceBackDelegate<CollisionImpl> > Type;
      };

    }
  }
}

#endif /* HEMELB_LB_STREAMERS_ADVECTIONDIFFUSIONSIMPLEBOUNCEBACK_H */
