// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_GEOMETRY_NEIGHBOURING_NEIGHBOURINGSITE_H
#define HEMELB_GEOMETRY_NEIGHBOURING_NEIGHBOURINGSITE_H

#include "geometry/Site.h"


namespace hemelb{
  namespace geometry{
    namespace neighbouring{
      class NeighbouringLatticeData;
      typedef BaseSite<NeighbouringLatticeData> NeighbouringSite;
      typedef BaseSite<const NeighbouringLatticeData> ConstNeighbouringSite;
    }
  }
}

#endif
