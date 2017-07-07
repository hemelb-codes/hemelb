
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_STREAMERS_ADVECTIONDIFFUSIONGUOZHENGSHI_H
#define HEMELB_LB_STREAMERS_ADVECTIONDIFFUSIONGUOZHENGSHI_H

#include "lb/streamers/AdvectionDiffusionStreamerTypeFactory.h"
#include "lb/streamers/GuoZhengShiDelegate.h"

namespace hemelb
{
  namespace lb
  {
    namespace streamers
    {

      template<typename CollisionImpl>
      struct AdvectionDiffusionGuoZhengShi
      {
          typedef VesselWallStreamerTypeFactory<CollisionImpl, GuoZhengShiDelegate<CollisionImpl> > Type;
      };

      template<typename CollisionImpl>
      struct AdvectionDiffusionGuoZhengShiIolet
      {
          typedef AdvectionDiffusionIoletStreamerTypeFactory<CollisionImpl, GuoZhengShiDelegate<CollisionImpl> > Type;
      };

    }
  }
}

#endif /* HEMELB_LB_STREAMERS_ADVECTIONDIFFUSIONGUOZHENGSHI_H */
