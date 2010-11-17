#ifndef HEMELB_TOPOLOGY_TOPOLOGYMANAGER_H
#define HEMELB_TOPOLOGY_TOPOLOGYMANAGER_H

#include "lb/LocalLatticeData.h"
#include "lb/GlobalLatticeData.h"
#include "topology/NetworkTopology.h"

#include <cmath>

namespace hemelb
{
  namespace topology
  {
    class BaseTopologyManager
    {
      public:
        void
        AssignFluidSitesToProcessors(int & proc_count,
                                     int & fluid_sites_per_unit,
                                     int & unvisited_fluid_sites,
                                     const int iCurrentProcId,
                                     const bool iIsMachineLevel,
                                     lb::LocalLatticeData * iLocalLatDat,
                                     const lb::GlobalLatticeData &iGlobLatDat,
                                     NetworkTopology * bNetTopology);

      protected:

        // This is here to prevent instantiation.
        BaseTopologyManager();

        // Site coordinates.
        struct SiteLocation
        {
            short int i, j, k;
        };
    };
  }
}

#endif /* HEMELB_TOPOLOGY_TOPOLOGYMANAGER_H */
