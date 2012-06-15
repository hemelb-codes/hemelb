#ifndef HEMELB_GEOMETRY_NEIGHBOURING_NEIGHBOURINGDATAMANAGER_H
#define HEMELB_GEOMETRY_NEIGHBOURING_NEIGHBOURINGDATAMANAGER_H

#include "geometry/neighbouring/NeighbouringLatticeData.h"
#include "net/net.h"
#include <vector>
#include <map>
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
                                  NeighbouringLatticeData & neighbouringLatticeData,
                                  net::InterfaceDelegationNet & net);
          void RegisterNeededSite(site_t globalId);
          void ShareNeeds();
          std::vector<proc_t> & GetProcsNeedingSite(site_t site){return procsNeedingEachSite[site];}
          std::vector<site_t> & GetNeededSites(){return neededSites;}

        private:
          const LatticeData & localLatticeData;
          NeighbouringLatticeData & neighbouringLatticeData;
          net::InterfaceDelegationNet & net;

          std::vector<site_t> neededSites;
          std::map<site_t, std::vector<proc_t> > procsNeedingEachSite;

      };

    }
  }
}

#endif
