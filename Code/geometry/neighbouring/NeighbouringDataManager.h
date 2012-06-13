#ifndef HEMELB_GEOMETRY_NEIGHBOURING_NEIGHBOURINGDATAMANAGER_H
#define HEMELB_GEOMETRY_NEIGHBOURING_NEIGHBOURINGDATAMANAGER_H

#include "geometry/neighbouring/NeighbouringLatticeData.h"
#include "net/net.h"

namespace hemelb
{
  namespace geometry
  {
    namespace neighbouring
    {

      class NeighbouringDataManager
      {
        public:
          NeighbouringDataManager(const LatticeData & localLatticeData,
                                  NeighbouringLatticeData & neighbouringLatticeData,net::BaseNet & net);
        private:
          const LatticeData & localLatticeData;
          const NeighbouringLatticeData & neighbouringLatticeData;
          const net::BaseNet & net;
      };

    }
  }
}

#endif
