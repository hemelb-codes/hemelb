
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "lb/stents/BoundaryValues.h"
#include "lb/iolets/BoundaryComms.h"
#include "util/utilityFunctions.h"
#include "util/fileutils.h"
#include <algorithm>
#include <fstream>

namespace hemelb
{
  namespace lb
  {
    namespace stents
    {
      BoundaryValues::BoundaryValues(geometry::LatticeData* latticeData,
                                     const std::vector<stents::Stent*> &incoming_stents,
                                     SimulationState* simulationState,
                                     const net::MpiCommunicator& comms,
                                     const util::UnitConverter& units) :
        net::IteratedAction(), totalStentCount(incoming_stents.size()), localStentCount(0),
            state(simulationState), unitConverter(units), bcComms(comms)
      {
        std::vector<int> *procsList = new std::vector<int>[totalStentCount];

        // Determine which iolets need comms and create them
        for (int stentIndex = 0; stentIndex < totalStentCount; stentIndex++)
        {
          // First create a copy of all iolets
          stents::Stent* stent = (incoming_stents[stentIndex])->Clone();

          stent->Initialise(&unitConverter);

          stents.push_back(stent);

          bool isStentOnThisProc = IsStentOnThisProc(latticeData, stentIndex);
          hemelb::log::Logger::Log<hemelb::log::Debug, hemelb::log::OnePerCore>("BOUNDARYVALUES.CC - isstentonthisproc? : %d", isStentOnThisProc);
          procsList[stentIndex] = GatherProcList(isStentOnThisProc);

          // With information on whether a proc has an IOlet and the list of procs for each IOlte
          // on the BC task we can create the comms
          if (isStentOnThisProc || bcComms.IsCurrentProcTheBCProc())
          {
            localStentCount++;
            localStentIDs.push_back(stentIndex);
//            hemelb::log::Logger::Log<hemelb::log::Warning, hemelb::log::OnePerCore>("BOUNDARYVALUES.H - stentIndex: %d", stentIndex);

//            if (stent->IsCommsRequired()) //DEREK: POTENTIAL MULTISCALE ISSUE (this if-statement)
//            {
              stent->SetComms(new iolets::BoundaryComms(state, procsList[stentIndex], bcComms, isStentOnThisProc));
//            }
          }
        }

        // Send out initial values
        Reset();

        // Clear up
        delete[] procsList;

        hemelb::log::Logger::Log<hemelb::log::Debug, hemelb::log::OnePerCore>("BOUNDARYVALUES.H - stentCount: %d, first stent ID %d", localStentCount, localStentIDs[0]);

      }

      BoundaryValues::~BoundaryValues()
      {

        for (int i = 0; i < totalStentCount; i++)
        {
          delete stents[i];
        }
      }

      bool BoundaryValues::IsStentOnThisProc(geometry::LatticeData* latticeData,
                                             int boundaryId)
      {
        for (site_t i = 0; i < latticeData->GetLocalFluidSiteCount(); i++)
        {
          const geometry::Site<geometry::LatticeData> site = latticeData->GetSite(i);

          if (site.GetStentId() == boundaryId)
          {
            return true;
          }
        }

        return true;
      }

      std::vector<int> BoundaryValues::GatherProcList(bool hasBoundary)
      {
        std::vector<int> processorsNeedingStentList(0);

        // This is where the info about whether a proc contains the given inlet/outlet is sent
        // If it does contain the given inlet/outlet it sends a true value, else it sends a false.
        int isStentOnThisProc = hasBoundary; // true if inlet i is on this proc

        // These should be bool, but MPI only supports MPI_INT
        // For each inlet/outlet there is an array of length equal to total number of procs.
        // Each stores true/false value. True if proc of rank equal to the index contains
        // the given inlet/outlet.

        std::vector<int> processorsNeedingStentFlags = bcComms.Gather(isStentOnThisProc, bcComms.GetBCProcRank());

        if (bcComms.IsCurrentProcTheBCProc())
        {
          // Now we have an array for each Stent with true (1) at indices corresponding to
          // processes that are members of that group. We have to convert this into arrays
          // of ints which store a list of processor ranks.
          for (proc_t process = 0; process < processorsNeedingStentFlags.size(); ++process)
          {
            if (processorsNeedingStentFlags[process])
            {
              processorsNeedingStentList.push_back(process);
            }
          }
        }

        return processorsNeedingStentList; // return by copy
      }

      void BoundaryValues::RequestComms()
      {
        for (int i = 0; i < localStentCount; i++)
        {
          HandleComms(GetLocalStent(i));
        }
      }

      void BoundaryValues::HandleComms(stents::Stent* stent)
      {

        if (stent->IsCommsRequired())
        {
          stent->DoComms(bcComms, state->GetTimeStep());
        }

      }

      void BoundaryValues::EndIteration()
      {
        for (int i = 0; i < localStentCount; i++)
        {
          if (GetLocalStent(i)->IsCommsRequired())
          {
            GetLocalStent(i)->GetComms()->FinishSend();
          }
        }
      }

      void BoundaryValues::FinishReceive()
      {
        for (int i = 0; i < localStentCount; i++)
        {
          if (GetLocalStent(i)->IsCommsRequired())
          {
            GetLocalStent(i)->GetComms()->Wait();
          }
        }
      }

      void BoundaryValues::Reset()
      {
        for (int i = 0; i < localStentCount; i++)
        {
          GetLocalStent(i)->Reset(*state);
          if (GetLocalStent(i)->IsCommsRequired())
          {
            GetLocalStent(i)->GetComms()->WaitAllComms();

          }
        }
      }

      // This assumes the program has already waited for comms to finish before
      LatticeDensity BoundaryValues::GetBoundaryDensity(const int index)
      {
        return stents[index]->GetDensity(state->Get0IndexedTimeStep());
      }

      LatticeDensity BoundaryValues::GetDensityMin(int iBoundaryId)
      {
        return stents[iBoundaryId]->GetDensityMin();
      }

      LatticeDensity BoundaryValues::GetDensityMax(int iBoundaryId)
      {
        return stents[iBoundaryId]->GetDensityMax();
      }

    }
  }
}
