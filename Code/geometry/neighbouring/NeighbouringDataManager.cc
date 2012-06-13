#include "geometry/neighbouring/NeighbouringDataManager.h"

namespace hemelb
{
  namespace geometry
  {
    namespace neighbouring
    {

      NeighbouringDataManager::NeighbouringDataManager(const LatticeData & localLatticeData,
                                                       NeighbouringLatticeData & neighbouringLatticeData,net::BaseNet & net) :
          localLatticeData(localLatticeData), neighbouringLatticeData(neighbouringLatticeData),net(net)
      {
      }

    }
  }
}
