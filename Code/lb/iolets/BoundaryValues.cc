
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "lb/iolets/BoundaryValues.h"
#include "util/utilityFunctions.h"
#include "util/fileutils.h"
#include <algorithm>
#include <fstream>

namespace hemelb
{
  namespace lb
  {
    namespace iolets
    {
      BoundaryValues::BoundaryValues(geometry::SiteType ioletType,
                                     geometry::LatticeData* latticeData,
                                     const std::vector<iolets::InOutLet*> &incoming_iolets,
                                     SimulationState* simulationState,
                                     comm::Async::Ptr commQ,
                                     const util::UnitConverter& units) :
        ioletType(ioletType), totalIoletCount(incoming_iolets.size()), localIoletCount(0),
	procsForIolet(incoming_iolets.size()),
	state(simulationState), unitConverter(units), asyncCommsQ(commQ)
      {
        // Determine which iolets need comms and create them
        for (int ioletIndex = 0; ioletIndex < totalIoletCount; ioletIndex++)
        {
          // First create a copy of all iolets
          iolets::InOutLet* iolet = (incoming_iolets[ioletIndex])->Clone();

          iolet->Initialise(&unitConverter);

          iolets.push_back(iolet);

          bool isIOletOnThisProc = IsIOletOnThisProc(ioletType, latticeData, ioletIndex);
          hemelb::log::Logger::Log<hemelb::log::Debug, hemelb::log::OnePerCore>("BOUNDARYVALUES.CC - isioletonthisproc? : %d", isIOletOnThisProc);
          procsForIolet[ioletIndex] = GatherProcList(isIOletOnThisProc);

          // With information on whether a proc has an IOlet and the list of procs for each IOlte
          // on the BC task we can create the comms
          if (isIOletOnThisProc || IsCurrentProcTheBCProc())
          {
            localIoletCount++;
            localIoletIDs.push_back(ioletIndex);
          }
        }

        // Send out initial values
        Reset();


        hemelb::log::Logger::Log<hemelb::log::Debug, hemelb::log::OnePerCore>("BOUNDARYVALUES.H - ioletCount: %d, first iolet ID %d", localIoletCount, localIoletIDs[0]);

      }

      BoundaryValues::~BoundaryValues()
      {

        for (int i = 0; i < totalIoletCount; i++)
        {
          delete iolets[i];
        }
      }

      bool BoundaryValues::IsIOletOnThisProc(geometry::SiteType ioletType,
                                             geometry::LatticeData* latticeData,
                                             int boundaryId)
      {
        for (site_t i = 0; i < latticeData->GetLocalFluidSiteCount(); i++)
        {
          const geometry::Site<geometry::LatticeData> site = latticeData->GetSite(i);

          if (site.GetSiteType() == ioletType && site.GetIoletId() == boundaryId)
          {
            return true;
          }
        }

        return true;
      }

      std::vector<int> BoundaryValues::GatherProcList(bool hasBoundary)
      {
        std::vector<int> processorsNeedingIoletList(0);

        // This is where the info about whether a proc contains the
        // given inlet/outlet is sent If it does contain the given
        // inlet/outlet it sends a true value, else it sends a false.
        // Ideally bools, but MPI doesn't support them :(
        char isIOletOnThisProc = hasBoundary; // true if inlet i is on this proc

        // For each inlet/outlet there is an array of length equal to
        // total number of procs.  Each stores true/false value. True
        // if proc of rank equal to the index contains the given
        // inlet/outlet.
	
        auto processorsNeedingIoletFlags = asyncCommsQ->GetComm()->Gather(isIOletOnThisProc, GetBCProcRank());

        if (IsCurrentProcTheBCProc())
        {
          // Now we have an array for each IOlet with true (1) at indices corresponding to
          // processes that are members of that group. We have to convert this into arrays
          // of ints which store a list of processor ranks.
          for (proc_t process = 0; process < processorsNeedingIoletFlags.size(); ++process)
          {
            if (processorsNeedingIoletFlags[process])
            {
              processorsNeedingIoletList.push_back(process);
            }
          }
        }

        return processorsNeedingIoletList; // return by move
      }


      void BoundaryValues::BeginAll() { /* not needed */ }
      void BoundaryValues::PreSend() { /* not needed */ }
      void BoundaryValues::PreWait() { /* not needed */ }
      void BoundaryValues::Wait() { /* not needed */ }
      void BoundaryValues::EndAll() { /* not needed */ }

      void BoundaryValues::Begin()
      {
	for (int i = 0; i < localIoletCount; i++)
        {
	  auto iolet = GetLocalIolet(i);
	  if (iolet->IsCommsRequired())
	  {
	    iolet->Begin(this);
	  }
        }
      }
      void BoundaryValues::Receive()
      {
	for (int i = 0; i < localIoletCount; i++)
        {
	  auto iolet = GetLocalIolet(i);
	  if (iolet->IsCommsRequired())
	  {
	    iolet->Receive(this, asyncCommsQ);
	  }
        }
      }
      void BoundaryValues::Send()
      {
	for (int i = 0; i < localIoletCount; i++)
        {
          auto iolet = GetLocalIolet(i);
	  if (iolet->IsCommsRequired())
	  {
	    iolet->Send(this, asyncCommsQ);
	  }
        }
      }
      // Wait done by AsyncActor
      void BoundaryValues::End()
      {
	for (int i = 0; i < localIoletCount; i++)
        {
          auto iolet = GetLocalIolet(i);
	  if (iolet->IsCommsRequired())
	  {
	    iolet->CommsComplete(this);
	  }
        }
      }
      
      void BoundaryValues::ForceCommunication()
      {
	Begin();
	{
	  auto tmpQ = comm::Async::New(asyncCommsQ->GetComm());
	  for (int i = 0; i < localIoletCount; i++)
	  {
	    auto iolet = GetLocalIolet(i);
	    if (iolet->IsCommsRequired())
	    {
	      iolet->Receive(this, tmpQ);
	      iolet->Send(this, tmpQ);
	    }
	  }
	}
	End();
      }

      void BoundaryValues::Reset()
      {
        for (int i = 0; i < localIoletCount; i++)
        {
          GetLocalIolet(i)->Reset(*state);
        }
	// TODO: Are we sure we need to force a comms cycle?
	ForceCommunication();
      }

      // This assumes the program has already waited for comms to finish before
      LatticeDensity BoundaryValues::GetBoundaryDensity(const int index)
      {
        return iolets[index]->GetDensity(state->Get0IndexedTimeStep());
      }

    }
  }
}
