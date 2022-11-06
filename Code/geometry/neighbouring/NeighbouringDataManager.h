// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_GEOMETRY_NEIGHBOURING_NEIGHBOURINGDATAMANAGER_H
#define HEMELB_GEOMETRY_NEIGHBOURING_NEIGHBOURINGDATAMANAGER_H

#include "geometry/Domain.h"
#include "geometry/neighbouring/NeighbouringDomain.h"
#include "geometry/neighbouring/RequiredSiteInformation.h"
#include "net/net.h"
#include "net/IteratedAction.h"
#include <vector>
#include <map>

namespace hemelb
{
  namespace geometry
  {
    class Domain;
    class FieldData;

    namespace neighbouring
    {

      class NeighbouringDataManager : public net::IteratedAction
      {
        public:
          NeighbouringDataManager(const FieldData& localLatticeData,
                                  NeighbouringFieldData& neighbouringFieldData,
                                  net::InterfaceDelegationNet & net);
          // Initially, the required site information will not be used -- we just transfer everything.
          // This considerably simplifies matters.
          // Nevertheless, we provide the interface here in its final form
          void RegisterNeededSite(site_t globalId, RequiredSiteInformation requirements =
                                      RequiredSiteInformation(true));
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
          // NB this is virtual so that the class can be tested.
          virtual proc_t ProcForSite(site_t site);
        protected:
          void RequestComms();
        private:
          const FieldData& localFieldData;
          NeighbouringFieldData& neighbouringFieldData;
          net::InterfaceDelegationNet & net;

          std::vector<site_t> neededSites;
          std::vector<std::vector<site_t> > needsEachProcHasFromMe;

          bool needsHaveBeenShared;

      };

    }
  }
}

#endif
