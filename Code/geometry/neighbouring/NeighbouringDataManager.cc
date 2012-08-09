#include <algorithm>

#include "geometry/neighbouring/NeighbouringDataManager.h"
#include "geometry/LatticeData.h"

#include "log/Logger.h"
namespace hemelb
{
  namespace geometry
  {
    namespace neighbouring
    {

      NeighbouringDataManager::NeighbouringDataManager(const LatticeData & localLatticeData,
                                                       NeighbouringLatticeData & neighbouringLatticeData,
                                                       net::InterfaceDelegationNet & net) :
          localLatticeData(localLatticeData), neighbouringLatticeData(neighbouringLatticeData), net(net), needsEachProcHasFromMe(net.GetCommunicator().GetSize())
      {
      }
      void NeighbouringDataManager::RegisterNeededSite(site_t globalId, RequiredSiteInformation requirements)
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
        return localLatticeData.ProcProvidingSiteByGlobalNoncontiguousId(site);
      }

      void NeighbouringDataManager::TransferNonFieldDependentInformation()
      {
        // Ordering is important here, to ensure the requests are registered in the same order
        // on the sending and receiving procs.
        // But, the needsEachProcHasFromMe is always ordered,
        // by the same order, as the neededSites, so this should be OK.
        for (std::vector<site_t>::iterator localNeed = neededSites.begin(); localNeed != neededSites.end(); localNeed++)
        {
          proc_t source = ProcForSite(*localNeed);
          NeighbouringSite site = neighbouringLatticeData.GetSite(*localNeed);

          net.RequestReceiveR(site.GetSiteData().GetIntersectionData(), source);
          net.RequestReceiveR(site.GetSiteData().GetOtherRawData(), source);
          net.RequestReceive(site.GetWallDistances(), localLatticeData.GetLatticeInfo().GetNumVectors() - 1, source);
          net.RequestReceiveR(site.GetWallNormal(), source);
        }
        for (proc_t other = 0; other < net.GetCommunicator().GetSize(); other++)
        {
          for (std::vector<site_t>::iterator needOnProcFromMe = needsEachProcHasFromMe[other].begin();
              needOnProcFromMe != needsEachProcHasFromMe[other].end(); needOnProcFromMe++)
          {
            site_t localContiguousId =
                localLatticeData.GetLocalContiguousIdFromGlobalNoncontiguousId(*needOnProcFromMe);


            Site site = const_cast<LatticeData&>(localLatticeData).GetSite(localContiguousId);
            // have to cast away the const, because no respect for const-ness for sends in MPI
            net.RequestSendR(site.GetSiteData().GetIntersectionData(), other);
            net.RequestSendR(site.GetSiteData().GetOtherRawData(), other);
            net.RequestSend(site.GetWallDistances(), localLatticeData.GetLatticeInfo().GetNumVectors() - 1, other);
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
        // Ordering is important here, to ensure the requests are registered in the same order
        // on the sending and receiving procs.
        // But, the needsEachProcHasFromMe is always ordered,
        // by the same order, as the neededSites, so this should be OK.
        for (std::vector<site_t>::iterator localNeed = neededSites.begin(); localNeed != neededSites.end(); localNeed++)
        {
          proc_t source = ProcForSite(*localNeed);
          NeighbouringSite site = neighbouringLatticeData.GetSite(*localNeed);
          net.RequestReceive(site.GetFOld(localLatticeData.GetLatticeInfo().GetNumVectors()),
                             localLatticeData.GetLatticeInfo().GetNumVectors(),
                             source);

        }
        for (proc_t other = 0; other < net.GetCommunicator().GetSize(); other++)
        {
          for (std::vector<site_t>::iterator needOnProcFromMe = needsEachProcHasFromMe[other].begin();
              needOnProcFromMe != needsEachProcHasFromMe[other].end(); needOnProcFromMe++)
          {
            site_t localContiguousId =
                localLatticeData.GetLocalContiguousIdFromGlobalNoncontiguousId(*needOnProcFromMe);
            Site site = const_cast<LatticeData&>(localLatticeData).GetSite(localContiguousId);
            // have to cast away the const, because no respect for const-ness for sends in MPI
            net.RequestSend(const_cast<distribn_t*>(site.GetFOld(localLatticeData.GetLatticeInfo().GetNumVectors())),
                            localLatticeData.GetLatticeInfo().GetNumVectors(),
                            other);

          }
        }
      }

      void NeighbouringDataManager::ShareNeeds()
      {
        // build a table of which procs needs can be achieved from which proc
        std::vector<std::vector<site_t> > needsIHaveFromEachProc(net.GetCommunicator().GetSize());
        std::vector<int> countOfNeedsIHaveFromEachProc(net.GetCommunicator().GetSize(), 0);
        for (std::vector<site_t>::iterator localNeed = neededSites.begin(); localNeed != neededSites.end(); localNeed++)
        {
          needsIHaveFromEachProc[ProcForSite(*localNeed)].push_back(*localNeed);
          countOfNeedsIHaveFromEachProc[ProcForSite(*localNeed)]++;

        }
        // every proc must send to all procs, how many it needs from that proc
        net.RequestAllToAllSend(countOfNeedsIHaveFromEachProc);

        // every proc must receive from all procs, how many it needs to give that proc
        std::vector<int> countOfNeedsOnEachProcFromMe(net.GetCommunicator().GetSize(), 0);
        net.RequestAllToAllReceive(countOfNeedsOnEachProcFromMe);
        net.Dispatch();

        for (proc_t other = 0; other < net.GetCommunicator().GetSize(); other++)
        {

          // now, for every proc, which I need something from,send the ids of those
          net.RequestSendV(needsIHaveFromEachProc[other], other);
          // and, for every proc, which needs something from me, receive those ids
          needsEachProcHasFromMe[other].resize(countOfNeedsOnEachProcFromMe[other]);
          net.RequestReceiveV(needsEachProcHasFromMe[other], other);
          // In principle, this bit could have been implemented as a separate GatherV onto every proc
          // However, in practice, we expect the needs to be basically local
          // so using point-to-point will be more efficient.
        }

        net.Dispatch();

      }
    }
  }
}
