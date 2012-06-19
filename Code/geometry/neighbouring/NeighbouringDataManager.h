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
          std::vector<site_t> &GetNeedsForProc(proc_t proc){return needsEachProcHasFromMe[proc];}
          std::vector<site_t> & GetNeededSites(){return neededSites;}
          void TransferNonFieldDependentInformation();
          void TransferFieldDependentInformation();
          virtual proc_t ProcForSite(site_t site); // virtual to make this class testable

        private:
          const LatticeData & localLatticeData;
          NeighbouringLatticeData & neighbouringLatticeData;
          net::InterfaceDelegationNet & net;

          std::vector<site_t> neededSites;
          std::vector<std::vector<site_t> > needsEachProcHasFromMe;

      };

    }
  }
}

#endif
