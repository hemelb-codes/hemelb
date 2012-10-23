// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_GEOMETRY_NEIGHBOURING_NEIGHBOURINGDATAMANAGER_H
#define HEMELB_GEOMETRY_NEIGHBOURING_NEIGHBOURINGDATAMANAGER_H

#include "geometry/neighbouring/NeighbouringLatticeData.h"
#include "geometry/neighbouring/RequiredSiteInformation.h"
#include "net/net.h"
#include "net/IteratedAction.h"
#include <vector>
#include <map>
namespace hemelb
{
  namespace geometry
  {
    namespace neighbouring
    {

      class NeighbouringDataManager : public net::IteratedAction
      {
        public:
          NeighbouringDataManager(const LatticeData & localLatticeData,
                                  NeighbouringLatticeData & neighbouringLatticeData,
                                  net::InterfaceDelegationNet & net);
          // Initially, the required site information will not be used -- we just transfer everything.
          // This considerably simplifies matters.
          // Nevertheless, we provide the interface here in its final form
          void RegisterNeededSite(site_t globalId,
                                  RequiredSiteInformation requirements = RequiredSiteInformation(true));
          void ShareNeeds();
          std::vector<site_t> &GetNeedsForProc(proc_t proc)
          {
            return needsEachProcHasFromMe[proc];
          }
          std::vector<site_t> & GetNeededSites()
          {
            return neededSites;
          }
          void TransferNonFieldDependentInformation();
          void TransferFieldDependentInformation();
          virtual proc_t ProcForSite(site_t site); // virtual to make this class testable
        protected:
          void RequestComms();
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
