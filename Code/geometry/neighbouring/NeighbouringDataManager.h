#ifndef HEMELB_GEOMETRY_NEIGHBOURING_NEIGHBOURINGDATAMANAGER_H
#define HEMELB_GEOMETRY_NEIGHBOURING_NEIGHBOURINGDATAMANAGER_H

#include "geometry/neighbouring/NeighbouringLatticeData.h"

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
                                  NeighbouringLatticeData & neighbouringLatticeData);
        private:
          const LatticeData & localLatticeData;
          const NeighbouringLatticeData & neighbouringLatticeData;
      };

    }
  }
}

#endif
