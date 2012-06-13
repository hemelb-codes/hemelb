#ifndef HEMELB_GEOMETRY_NEIGHBOURING_NEIGHBOURINGDATAMANAGER_H
#define HEMELB_GEOMETRY_NEIGHBOURING_NEIGHBOURINGDATAMANAGER_H

#include "geometry/neighbouring/NeighbouringLatticeData.h"
#include "net/net.h"
#include <vector>
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
                                  NeighbouringLatticeData & neighbouringLatticeData,net::InterfaceDelegationNet & net);
          void RegisterNeededSite(site_t globalId);
        private:
          const LatticeData & localLatticeData;
          const NeighbouringLatticeData & neighbouringLatticeData;
          const net::InterfaceDelegationNet & net;

          std::vector<site_t> neededSites;
      };

    }
  }
}

#endif
