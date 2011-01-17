#ifndef HEMELB_TOPOLOGY_TOPOLOGYMANAGER_H
#define HEMELB_TOPOLOGY_TOPOLOGYMANAGER_H

#include "lb/GlobalLatticeData.h"
#include "topology/NetworkTopology.h"

namespace hemelb
{
  namespace topology
  {
    class TopologyManager
    {
      public:
        TopologyManager(NetworkTopology* bNetworkTopology, bool *oSuccess);

        void DecomposeDomain(int iTotalFluidSites,
                             NetworkTopology* bNetTop,
                             const lb::GlobalLatticeData & bGlobLatDat);

      private:
        // Site coordinates.
        struct SiteLocation
        {
            short int i, j, k;
        };

        void
        AssignFluidSitesToProcessors(int & proc_count,
                                     int & fluid_sites_per_unit,
                                     int & unvisited_fluid_sites,
                                     const int iCurrentProcId,
                                     const bool iIsMachineLevel,
                                     const lb::GlobalLatticeData &iGlobLatDat,
                                     NetworkTopology * bNetTopology);
    };
  }
}

#endif /* HEMELB_TOPOLOGY_CONCRETETOPOLOGYMANAGER_H */
