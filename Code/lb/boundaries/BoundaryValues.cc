#include "lb/boundaries/BoundaryValues.h"
#include "lb/boundaries/BoundaryComms.h"
#include "util/utilityFunctions.h"
#include "util/fileutils.h"
#include <algorithm>
#include <fstream>

namespace hemelb
{
  namespace lb
  {
    namespace boundaries
    {

      BoundaryValues::BoundaryValues(geometry::SiteType ioletType,
                                      geometry::LatticeData* latticeData,
                                      const std::vector<iolets::InOutLet*> &incoming_iolets,
                                      SimulationState* simulationState,
                                      util::UnitConverter* units):
          net::IteratedAction(),totalIoletCount(incoming_iolets.size()), localIoletCount(0), state(simulationState),  unitConverter(units)
      {

        std::vector<int> *procsList = new std::vector<int>[totalIoletCount];

        // Determine which iolets need comms and create them
        for (int ioletIndex = 0; ioletIndex < totalIoletCount; ioletIndex++)
        {
          // First create a copy of all iolets
          iolets::InOutLet* iolet = incoming_iolets[ioletIndex]->Clone();

          iolet->Initialise(unitConverter);

          iolets.push_back(iolet);

          bool isIOletOnThisProc = IsIOletOnThisProc(ioletType, latticeData, ioletIndex);
          procsList[ioletIndex] = GatherProcList(isIOletOnThisProc);

          // With information on whether a proc has an IOlet and the list of procs for each IOlte
          // on the BC task we can create the comms
          if (isIOletOnThisProc || IsCurrentProcTheBCProc())
          {
            localIoletCount++;

            localIoletIDs.push_back(ioletIndex);
            if (iolet->IsCommsRequired())
            {
              iolet->SetComms(new BoundaryComms(state, procsList[ioletIndex], isIOletOnThisProc));
            }
          }
        }

        // Send out initial values
        Reset();

        // Clear up
        delete[] procsList;

      }

      BoundaryValues::~BoundaryValues()
      {

        for (int i = 0; i <totalIoletCount; i++)
        {
          delete iolets[i];
        }
      }

      bool BoundaryValues::IsIOletOnThisProc(geometry::SiteType ioletType, geometry::LatticeData* latticeData, int boundaryId)
      {
        for (site_t i = 0; i < latticeData->GetLocalFluidSiteCount(); i++)
        {
          const geometry::Site site = latticeData->GetSite(i);

          if (site.GetSiteType() == ioletType && site.GetBoundaryId() == boundaryId)
          {
            return true;
          }
        }

        return false;
      }

      std::vector<int> BoundaryValues::GatherProcList(bool hasBoundary)
      {
        std::vector<int> processorsNeedingIoletList(0);

        // This is where the info about whether a proc contains the given inlet/outlet is sent
        // If it does contain the given inlet/outlet it sends a true value, else it sends a false.
        int isIOletOnThisProc = hasBoundary; // true if inlet i is on this proc

        // These should be bool, but MPI only supports MPI_INT
        // For each inlet/outlet there is an array of length equal to total number of procs.
        // Each stores true/false value. True if proc of rank equal to the index contains
        // the given inlet/outlet.
        proc_t processorCount = topology::NetworkTopology::Instance()->GetProcessorCount();
        int *processorsNeedingIoletFlags = new int[processorCount];

        MPI_Gather(&isIOletOnThisProc,
                   1,
                   hemelb::MpiDataType(isIOletOnThisProc),
                   processorsNeedingIoletFlags,
                   1,
                   hemelb::MpiDataType(processorsNeedingIoletFlags[0]),
                   GetBCProcRank(),
                   MPI_COMM_WORLD);

        if (IsCurrentProcTheBCProc())
        {
          // Now we have an array for each IOlet with true (1) at indices corresponding to
          // processes that are members of that group. We have to convert this into arrays
          // of ints which store a list of processor ranks.
          for (proc_t process = 0; process < processorCount; process++)
          {
            if (processorsNeedingIoletFlags[process])
            {
              processorsNeedingIoletList.push_back(process);
            }
          }
        }

        // Clear up
        delete[] processorsNeedingIoletFlags;

        return processorsNeedingIoletList; // return by copy
      }

      bool BoundaryValues::IsCurrentProcTheBCProc()
      {
        return topology::NetworkTopology::Instance()->GetLocalRank() == GetBCProcRank();
      }

      proc_t BoundaryValues::GetBCProcRank()
      {
        return 0;
      }

      void BoundaryValues::RequestComms()
      {
        for (int i = 0; i < localIoletCount; i++)
        {
          HandleComms(GetLocalIolet(i));
        }
      }

      void BoundaryValues::HandleComms(iolets::InOutLet* iolet)
      {

        if (iolet->IsCommsRequired())
        {
          iolet->DoComms(IsCurrentProcTheBCProc(),state->GetTimeStep());
        }

      }

      void BoundaryValues::EndIteration()
      {
        for (int i = 0; i < localIoletCount; i++)
        {
          if (GetLocalIolet(i)->IsCommsRequired())
          {
            GetLocalIolet(i)->GetComms()->FinishSend();
          }
        }
      }

      void BoundaryValues::FinishReceive()
      {
        for (int i = 0; i < localIoletCount; i++)
        {
          if (GetLocalIolet(i)->IsCommsRequired())
          {
            GetLocalIolet(i)->GetComms()->Wait();
          }
        }
      }

      void BoundaryValues::Reset()
      {
        for (int i = 0; i < localIoletCount; i++)
        {
          GetLocalIolet(i)->Reset(*state);
          if (GetLocalIolet(i)->IsCommsRequired())
          {
            GetLocalIolet(i)->GetComms()->WaitAllComms();

          }
        }
      }

      // This assumes the program has already waited for comms to finish before
      distribn_t BoundaryValues::GetBoundaryDensity(const int index)
      {
        return iolets[index]->GetDensity(state->Get0IndexedTimeStep());
      }

      distribn_t BoundaryValues::GetDensityMin(int iBoundaryId)
      {
        return iolets[iBoundaryId]->GetDensityMin();
      }

      distribn_t BoundaryValues::GetDensityMax(int iBoundaryId)
      {
        return iolets[iBoundaryId]->GetDensityMax();
      }

    }
  }
}
