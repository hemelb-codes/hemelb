
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <algorithm>

#include "geometry/neighbouring/NeighbouringDataManager.h"
#include "geometry/LatticeData.h"
#include "log/Logger.h"
#include "comm/MapAllToAll.h"

namespace hemelb
{
  namespace geometry
  {
    namespace neighbouring
    {

      NeighbouringDataManager::NeighbouringDataManager(
          const LatticeData & localLatticeData, NeighbouringLatticeData & neighbouringLatticeData,
          net::InterfaceDelegationNet & net) :
          localLatticeData(localLatticeData), neighbouringLatticeData(neighbouringLatticeData),
              net(net), needsHaveBeenShared(false)
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
        return localLatticeData.ProcProvidingSiteByGlobalNoncontiguousId(site);
      }

      void NeighbouringDataManager::TransferNonFieldDependentInformation()
      {
        // Ordering is important here, to ensure the requests are registered in the same order
        // on the sending and receiving procs.
        // But, the needsEachProcHasFromMe is always ordered,
        // by the same order, as the neededSites, so this should be OK.
        for (IdVec::iterator localNeed = neededSites.begin();
            localNeed != neededSites.end(); localNeed++)
        {
          proc_t source = ProcForSite(*localNeed);
          NeighbouringSite site = neighbouringLatticeData.GetSite(*localNeed);

          net.RequestReceiveR(site.GetSiteData().GetWallIntersectionData(), source);
          net.RequestReceiveR(site.GetSiteData().GetIoletIntersectionData(), source);
          net.RequestReceiveR(site.GetSiteData().GetIoletId(), source);
          net.RequestReceiveR(site.GetSiteData().GetSiteType(), source);
          net.RequestReceive(site.GetWallDistances(),
                             localLatticeData.GetLatticeInfo().GetNumVectors() - 1,
                             source);
          net.RequestReceiveR(site.GetWallNormal(), source);
        }

        for (IdsMap::const_iterator iter = needsEachProcHasFromMe.begin();
            iter != needsEachProcHasFromMe.end();
            ++iter)
        {
          proc_t other = iter->first;
          const IdVec& neededIds = iter->second;
          for (IdVec::const_iterator needOnProcFromMe = neededIds.begin();
              needOnProcFromMe != neededIds.end();
              ++needOnProcFromMe)
          {
            site_t localContiguousId =
                            localLatticeData.GetLocalContiguousIdFromGlobalNoncontiguousId(*needOnProcFromMe);

            const Site<const LatticeData> site = localLatticeData.GetSite(localContiguousId);
            const SiteData& sd = site.GetSiteData();

            net.RequestSendR(sd.GetWallIntersectionData(), other);
            net.RequestSendR(sd.GetIoletIntersectionData(), other);
            net.RequestSendR(sd.GetIoletId(), other);
            net.RequestSendR(sd.GetSiteType(), other);
            net.RequestSend(site.GetWallDistances(),
                            localLatticeData.GetLatticeInfo().GetNumVectors() - 1,
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
        for (IdVec::iterator localNeed = neededSites.begin();
            localNeed != neededSites.end(); localNeed++)
        {
          proc_t source = ProcForSite(*localNeed);
          NeighbouringSite site = neighbouringLatticeData.GetSite(*localNeed);
          net.RequestReceive(site.GetFOld(localLatticeData.GetLatticeInfo().GetNumVectors()),
                             localLatticeData.GetLatticeInfo().GetNumVectors(),
                             source);

        }

        const unsigned Q = localLatticeData.GetLatticeInfo().GetNumVectors();
        for (IdsMap::const_iterator iter = needsEachProcHasFromMe.begin();
            iter != needsEachProcHasFromMe.end();
            ++iter)
        {
          proc_t other = iter->first;
          const IdVec& neededIds = iter->second;
          for (IdVec::const_iterator needOnProcFromMe = neededIds.begin();
              needOnProcFromMe != neededIds.end(); ++needOnProcFromMe)
          {
            site_t localContiguousId =
                localLatticeData.GetLocalContiguousIdFromGlobalNoncontiguousId(*needOnProcFromMe);
            const Site<const LatticeData> site = localLatticeData.GetSite(localContiguousId);
            // have to cast away the const, because no respect for const-ness for sends in MPI
            net.RequestSend(site.GetFOld(Q),
                            Q,
                            other);
          }
        }
      }

      void NeighbouringDataManager::ShareNeeds()
      {
        hemelb::log::Logger::Log<hemelb::log::Debug, hemelb::log::OnePerCore>("NDM ShareNeeds().");
        //if (needsHaveBeenShared == true)
        //  return; //TODO: Fix!

        // build a table of which sites are needed by this rank, by other rank.
        IdsMap needsIHaveFromEachProc;
        // This map will count the number per-rank
        CountMap countOfNeedsIHaveFromEachProc;
        for (IdVec::iterator localNeed = neededSites.begin();
            localNeed != neededSites.end(); localNeed++)
        {
          site_t needId = *localNeed;
          int needRank = ProcForSite(needId);
          needsIHaveFromEachProc[needRank].push_back(needId);
          countOfNeedsIHaveFromEachProc[needRank]++;
        }

        // This will store the numbers of sites other ranks need from this rank
        CountMap countOfNeedsOnEachProcFromMe;
        // This is collective
        comm::Communicator::ConstPtr comms = net.GetCommunicator();
        comm::MapAllToAll(comms, countOfNeedsIHaveFromEachProc, countOfNeedsOnEachProcFromMe, 1234);

	//comm::Request::ReqVec requestQueue;
	auto requestQueue = comms->MakeRequestList();
	
        // Now, for every rank, which I need something from, send the ids of those
        for (CountMap::const_iterator countIt = countOfNeedsIHaveFromEachProc.begin();
            countIt != countOfNeedsIHaveFromEachProc.end();
            ++countIt)
        {
          int other = countIt->first;
          requestQueue->push_back(comms->Isend(needsIHaveFromEachProc[other], other));
        }

        // And for every rank, which needs something from me, receive those ids
        for (CountMap::const_iterator countIt = countOfNeedsOnEachProcFromMe.begin();
            countIt != countOfNeedsOnEachProcFromMe.end();
            ++countIt)
        {
          int other = countIt->first;
          int size = countIt->second;
          IdVec& otherNeeds = needsEachProcHasFromMe[other];
          otherNeeds.resize(size);
          requestQueue->push_back(comms->Irecv(otherNeeds, other));
        }

	requestQueue->WaitAll();
        needsHaveBeenShared = true;
      }
    }
  }
}
