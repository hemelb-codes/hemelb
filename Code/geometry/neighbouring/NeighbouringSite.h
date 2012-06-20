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
