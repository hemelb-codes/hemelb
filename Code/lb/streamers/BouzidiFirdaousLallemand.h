
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

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
