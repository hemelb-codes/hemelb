#include "geometry/neighbouring/NeighbouringDataManager.h"

namespace hemelb
{
  namespace geometry
  {
    namespace neighbouring
    {

      NeighbouringDataManager::NeighbouringDataManager(const LatticeData & localLatticeData,
                                                       NeighbouringLatticeData & neighbouringLatticeData,net::InterfaceDelegationNet & net) :
          localLatticeData(localLatticeData), neighbouringLatticeData(neighbouringLatticeData),net(net)
      {
      }
      void NeighbouringDataManager::RegisterNeededSite(site_t globalId){
        neededSites.push_back(globalId);
      }
    }
  }
}
