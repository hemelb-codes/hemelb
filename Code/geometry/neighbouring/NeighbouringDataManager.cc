// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <algorithm>

#include "geometry/neighbouring/NeighbouringDataManager.h"
#include "geometry/Domain.h"
#include "geometry/FieldData.h"

#include "log/Logger.h"
namespace hemelb
{
  namespace geometry
  {
    namespace neighbouring
    {

      NeighbouringDataManager::NeighbouringDataManager(
              const FieldData& localLatticeData, NeighbouringFieldData& neighbouringLatticeData,
              net::InterfaceDelegationNet & net) :
              localFieldData(localLatticeData), neighbouringFieldData(neighbouringLatticeData),
              net(net), needsEachProcHasFromMe(net.Size()), needsHaveBeenShared(false)
      {
      }
      void NeighbouringDataManager::RegisterNeededSite(site_t globalId,
                                                       RequiredSiteInformation requirements)
      {
        //ignore the requirements, we require everying.
        if (std::find(neededSites.begin(), neededSites.end(), globalId) == neededSites.end())
        {
          neededSites.push_back(globalId);
        }
        else
        {
          // Merge requirements.
        }
      }

      proc_t NeighbouringDataManager::ProcForSite(site_t site)
      {
        return localFieldData.GetDomain().ProcProvidingSiteByGlobalNoncontiguousId(site);
      }

      void NeighbouringDataManager::TransferNonFieldDependentInformation()
      {
        auto&& neigh_dom = neighbouringFieldData.GetDomain();
        auto&& local_dom = localFieldData.GetDomain();
        // Ordering is important here, to ensure the requests are registered in the same order
        // on the sending and receiving procs.
        // But, the needsEachProcHasFromMe is always ordered,
        // by the same order, as the neededSites, so this should be OK.
        for (std::vector<site_t>::iterator localNeed = neededSites.begin();
            localNeed != neededSites.end(); localNeed++)
        {
          proc_t source = ProcForSite(*localNeed);
          auto site = neigh_dom.GetSite(*localNeed);

          net.RequestReceiveR(site.GetSiteData().GetWallIntersectionData(), source);
          net.RequestReceiveR(site.GetSiteData().GetIoletIntersectionData(), source);
          net.RequestReceiveR(site.GetSiteData().GetIoletId(), source);
          net.RequestReceiveR(site.GetSiteData().GetSiteType(), source);
          net.RequestReceive(site.GetWallDistances(),
                             local_dom.GetLatticeInfo().GetNumVectors() - 1,
                             source);
          net.RequestReceiveR(site.GetWallNormal(), source);
        }
        for (proc_t other = 0; other < net.Size(); other++)
        {
          for (std::vector<site_t>::iterator needOnProcFromMe =
              needsEachProcHasFromMe[other].begin();
              needOnProcFromMe != needsEachProcHasFromMe[other].end(); needOnProcFromMe++)
          {
            site_t localContiguousId =
                local_dom.GetLocalContiguousIdFromGlobalNoncontiguousId(*needOnProcFromMe);

//            Site<Domain> site =
//                const_cast<Domain&>(localFieldData).GetSite(localContiguousId);
auto const site = localFieldData.GetSite(localContiguousId);
            // have to cast away the const, because no respect for const-ness for sends in MPI
            net.RequestSendR(site.GetSiteData().GetWallIntersectionData(), other);
            net.RequestSendR(site.GetSiteData().GetIoletIntersectionData(), other);
            net.RequestSendR(site.GetSiteData().GetIoletId(), other);
            net.RequestSendR(site.GetSiteData().GetSiteType(), other);
            net.RequestSend(site.GetWallDistances(),
                            local_dom.GetLatticeInfo().GetNumVectors() - 1,
                            other);
            net.RequestSendR(site.GetWallNormal(), other);
          }
        }
        net.Dispatch();
      }

      void NeighbouringDataManager::TransferFieldDependentInformation()
      {
        RequestComms();
        net.Dispatch();
      }

      void NeighbouringDataManager::RequestComms()
      {
        /*if (needsHaveBeenShared == false)
         {
         hemelb::log::Logger::Log<hemelb::log::Debug, hemelb::log::OnePerCore>("NDM needs are shared now.");
         ShareNeeds();
         }*/ ///TODO: Re-enable!
        // Ordering is important here, to ensure the requests are registered in the same order
        // on the sending and receiving procs.
        // But, the needsEachProcHasFromMe is always ordered,
        // by the same order, as the neededSites, so this should be OK.
        auto&& local_dom = localFieldData.GetDomain();
        auto const NV = local_dom.GetLatticeInfo().GetNumVectors();
        for (std::vector<site_t>::iterator localNeed = neededSites.begin();
            localNeed != neededSites.end(); localNeed++)
        {
          proc_t source = ProcForSite(*localNeed);
          auto site = neighbouringFieldData.GetSite(*localNeed);
          net.RequestReceive(site.GetFOld(NV),
                             NV,
                             source);

        }
        for (proc_t other = 0; other < net.Size(); other++)
        {
          for (std::vector<site_t>::iterator needOnProcFromMe =
              needsEachProcHasFromMe[other].begin();
              needOnProcFromMe != needsEachProcHasFromMe[other].end(); needOnProcFromMe++)
          {
            site_t localContiguousId =
                local_dom.GetLocalContiguousIdFromGlobalNoncontiguousId(*needOnProcFromMe);
//            Site<Domain> site =
//                const_cast<Domain&>(localFieldData).GetSite(localContiguousId);
            auto const site = localFieldData.GetSite(localContiguousId);
            // have to cast away the const, because no respect for const-ness for sends in MPI
            net.RequestSend(site.GetFOld(NV),
                            NV,
                            other);

          }
        }
      }

      void NeighbouringDataManager::ShareNeeds()
      {
        hemelb::log::Logger::Log<hemelb::log::Debug, hemelb::log::OnePerCore>("NDM ShareNeeds().");
        //if (needsHaveBeenShared == true)
        //  return; //TODO: Fix!
        
        // build a table of which procs needs can be achieved from which proc
        std::vector<std::vector<site_t> > needsIHaveFromEachProc(net.Size());
        std::vector<int> countOfNeedsIHaveFromEachProc(net.Size(), 0);
        for (std::vector<site_t>::iterator localNeed = neededSites.begin();
            localNeed != neededSites.end(); localNeed++)
        {
          needsIHaveFromEachProc[ProcForSite(*localNeed)].push_back(*localNeed);
          countOfNeedsIHaveFromEachProc[ProcForSite(*localNeed)]++;

        }

        // every proc must send to all procs, how many it needs from that proc
        net.RequestAllToAllSend(countOfNeedsIHaveFromEachProc);

        // every proc must receive from all procs, how many it needs to give that proc
        std::vector<int> countOfNeedsOnEachProcFromMe(net.Size(), 0);
        net.RequestAllToAllReceive(countOfNeedsOnEachProcFromMe);
        net.Dispatch();

        const int netSize = net.Size(); // avoid calling e.g. MPI_Comm_size many times; it is not inlined
        for (proc_t other = 0; other < netSize; other++)
        {

          // now, for every proc, which I need something from,send the ids of those
          net.RequestSendV(std::span<const site_t>(needsIHaveFromEachProc[other]), other);
          // and, for every proc, which needs something from me, receive those ids
          needsEachProcHasFromMe[other].resize(countOfNeedsOnEachProcFromMe[other]);
          net.RequestReceiveV(std::span<site_t>(needsEachProcHasFromMe[other]), other);
          // In principle, this bit could have been implemented as a separate GatherV onto every proc
          // However, in practice, we expect the needs to be basically local
          // so using point-to-point will be more efficient.
        }

        net.Dispatch();
        needsHaveBeenShared = true;
      }
    }
  }
}
