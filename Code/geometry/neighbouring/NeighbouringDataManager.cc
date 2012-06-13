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
          localLatticeData(localLatticeData), neighbouringLatticeData(neighbouringLatticeData), net(net)
      {
      }
      void NeighbouringDataManager::RegisterNeededSite(site_t globalId)
      {
        neededSites.push_back(globalId);
      }

      void NeighbouringDataManager::ShareNeeds()
      {
        // build a table of which procs needs can be achieved from which proc
        std::vector<std::vector<site_t> > needsIHaveFromEachProc(net.GetCommunicator().GetSize());
        std::vector<int> countOfNeedsIHaveFromEachProc(net.GetCommunicator().GetSize(), 0);
        for (std::vector<site_t>::iterator localNeed = neededSites.begin(); localNeed != neededSites.end(); localNeed++)
        {
          log::Logger::Log<log::Debug, log::OnePerCore>("Considering source for need %i",*localNeed);
          log::Logger::Log<log::Debug, log::OnePerCore>("Source should be: %d",localLatticeData.ProcProvidingSiteByGlobalNoncontiguousId(*localNeed));
          needsIHaveFromEachProc[localLatticeData.ProcProvidingSiteByGlobalNoncontiguousId(*localNeed)].push_back(*localNeed);
          countOfNeedsIHaveFromEachProc[localLatticeData.ProcProvidingSiteByGlobalNoncontiguousId(*localNeed)]++;
        }
        // every proc must send to all procs, how many it needs from that proc
        net.RequestAllToAllSend(countOfNeedsIHaveFromEachProc);

        // every proc must receive from all procs, how many it needs to give that proc
        std::vector<int> countOfNeedsOnEachProcFromMe(net.GetCommunicator().GetSize(), 0);
        std::vector<std::vector<site_t> > needsEachProcHasFromMe(net.GetCommunicator().GetSize());
        net.RequestAllToAllReceive(countOfNeedsOnEachProcFromMe);

        net.Dispatch();
        for (proc_t other = 0; other < net.GetCommunicator().GetSize(); other++)
        {
          // now, for every proc, which I need something from,send the ids of those
          net.RequestSendV(needsIHaveFromEachProc[other], other);
          // and, for every proc, which needs something from me, receive those ids
          net.RequestReceiveV(needsEachProcHasFromMe[other], other);
          // In principle, this bit could have been implemented as a separate GatherV onto every proc
          // However, in practice, we expect the needs to be basically local
          // so using point-to-point will be more efficient.
        }
        net.Dispatch();
        // Now, every proc knows the needs, of which other procs it must send sites to
        // so iterate that and invert the matrix into the map we use for storing this
        for (proc_t other = 0; other < net.GetCommunicator().GetSize(); other++)
        {
          for (std::vector<site_t>::iterator neededSite = needsEachProcHasFromMe[other].begin();
              neededSite != needsEachProcHasFromMe[other].end(); neededSite++)
          {
            procsNeedingEachSite.insert(std::pair<site_t, proc_t>(*neededSite, other));
          }
        }
      }
    }
  }
}
