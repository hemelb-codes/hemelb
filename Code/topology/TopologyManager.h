#ifndef HEMELB_TOPOLOGY_TOPOLOGYMANAGER_H
#define HEMELB_TOPOLOGY_TOPOLOGYMANAGER_H

#include "lb/LocalLatticeData.h"
#include "lb/GlobalLatticeData.h"
#include "topology/NetworkTopology.h"

namespace hemelb
{
  namespace topology
  {
    class TopologyManager
    {
      public:
        void FindTopology(NetworkTopology &bNetworkTopology,
                          bool & oWasSuccessful);

        void
        AssignFluidSitesToProcessors(int & proc_count,
                                     int & fluid_sites_per_unit,
                                     int & unvisited_fluid_sites,
                                     const int iCurrentProcId,
                                     const int unitLevel,
                                     lb::LocalLatticeData * iLocalLatDat,
                                     lb::GlobalLatticeData &iGlobLatDat,
                                     NetworkTopology * bNetTopology);

      private:
        // Site coordinates.
        struct SiteLocation
        {
            short int i, j, k;
        };
    };
  }
}

#endif /* HEMELB_TOPOLOGY_TOPOLOGYMANAGER_H */
