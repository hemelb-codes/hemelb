#include "lb/boundaries/BoundaryValues.h"
#include "util/utilityFunctions.h"
#include "util/fileutils.h"
#include <algorithm>
#include <fstream>
#include <math.h>

namespace hemelb
{
  namespace lb
  {
    namespace boundaries
    {

      BoundaryValues::BoundaryValues(geometry::LatticeData::SiteType IOtype,
                                     geometry::LatticeData* iLatDat,
                                     std::vector<iolets::InOutLet*> &iiolets,
                                     SimulationState* iSimState) :
        net::IteratedAction(), mState(iSimState)
      {
        nTotIOlets = (int) iiolets.size();

        FindBCProcRank();

        std::vector<int> *procsList = new std::vector<int>[nTotIOlets];

        nIOlets = 0;
        ioletIDs.resize(0);
        iolets.resize(0);
        mComms.resize(0);

        for (int i = 0; i < nTotIOlets; i++)
        {
          iolets.push_back(iiolets[i]->Clone());

          bool IOletOnThisProc = IsIOletOnThisProc(IOtype, iLatDat, i);

          procsList[i] = GatherProcList(IOletOnThisProc);

          if (IOletOnThisProc || IsCurrentProcTheBCProc())
          {
            nIOlets++;
            ioletIDs.push_back(i);
            if (iolets[i]->DoComms())
            {
              mComms.push_back(new BoundaryComms(mState, procsList[i], IOletOnThisProc, BCproc));
            }
            else
            {
              mComms.push_back(NULL);
            }
          }
        }

        density_cycle = new std::vector<distribn_t>[nIOlets];

        Reset();

        // Clear up
        delete[] procsList;

      }

      BoundaryValues::~BoundaryValues()
      {

        delete[] density_cycle;

        for (int i = 0; i < nIOlets; i++)
        {
          delete iolets[i];
        }

        for (int i = 0; i < nIOlets; i++)
        {
          if (mComms[i] != NULL)
          {
            delete mComms[i];
          }
        }
      }

      void BoundaryValues::FindBCProcRank()
      {
        proc_t BCrank = 0;

        if (IsCurrentProcTheBCProc())
          BCrank = topology::NetworkTopology::Instance()->GetLocalRank();

        // Since only one proc will update BCrank, the sum of all BCrank is the BCproc
        MPI_Allreduce(&BCrank, &BCproc, 1, hemelb::MpiDataType(BCrank), MPI_SUM, MPI_COMM_WORLD);
      }

      bool BoundaryValues::IsIOletOnThisProc(geometry::LatticeData::SiteType IOtype,
                                             geometry::LatticeData* iLatDat,
                                             int iBoundaryId)
      {
        for (site_t i = 0; i < iLatDat->GetLocalFluidSiteCount(); i++)
        {
          if (iLatDat->GetSiteType(i) == IOtype && iLatDat->GetBoundaryId(i) == iBoundaryId)
          {
            return true;
          }
        }

        return false;
      }

      std::vector<int> BoundaryValues::GatherProcList(bool hasBoundary)
      {
        std::vector<int> procsList(0);

        // This is where the info about whether a proc contains the given inlet/outlet is sent
        // If it does contain the given inlet/outlet it sends a true value, else it sends a false.
        int IOletOnThisProc = hasBoundary; // true if inlet i is on this proc

        // These should be bool, but MPI only supports MPI_INT
        // For each inlet/outlet there is an array of length equal to total number of procs.
        // Each stores true/false value. True if proc of rank equal to the index contains
        // the given inlet/outlet.
        int nTotProcs = topology::NetworkTopology::Instance()->GetProcessorCount();
        int *boolList = new int[nTotProcs];

        MPI_Gather(&IOletOnThisProc,
                   1,
                   hemelb::MpiDataType(IOletOnThisProc),
                   boolList,
                   1,
                   hemelb::MpiDataType(boolList[0]),
                   BCproc,
                   MPI_COMM_WORLD);

        if (IsCurrentProcTheBCProc())
        {
          // Now we have an array for each IOlet with true (1) at indices corresponding to
          // processes that are members of that group. We have to convert this into arrays
          // of ints which store a list of processor ranks.
          for (int j = 0; j < nTotProcs; j++)
          {
            if (boolList[j])
            {
              procsList.push_back(j);
            }
          }
        }

        // Clear up
        delete[] boolList;

        return procsList;
      }

      bool BoundaryValues::IsCurrentProcTheBCProc()
      {
        return topology::NetworkTopology::Instance()->IsCurrentProcTheIOProc();
      }

      void BoundaryValues::RequestComms()
      {
        if (IsCurrentProcTheBCProc())
        {
          for (int i = 0; i < nIOlets; i++)
          {
            unsigned long time_step = (mState->GetTimeStep() - 1) % density_cycle[i].size();

            if (iolets[ioletIDs[i]]->DoComms())
            {
              iolets[ioletIDs[i]]->UpdateCycle(density_cycle[i], mState);
              mComms[i]->Send(&density_cycle[i][time_step]);
            }
            else
            {
              iolets[ioletIDs[i]]->UpdateCycle(density_cycle[i], mState);
              iolets[ioletIDs[i]]->density = density_cycle[i][time_step];
            }
          }
        }
        else
        {
          for (int i = 0; i < nIOlets; i++)
          {
            if (iolets[ioletIDs[i]]->DoComms())
            {
              mComms[i]->Receive(&iolets[ioletIDs[i]]->density);
            }
            else
            {
              unsigned long time_step = (mState->GetTimeStep() - 1) % density_cycle[i].size();
              iolets[ioletIDs[i]]->UpdateCycle(density_cycle[i], mState);
              iolets[ioletIDs[i]]->density = density_cycle[i][time_step];
            }
          }
        }
      }

      void BoundaryValues::EndIteration()
      {
        for (int i = 0; i < nIOlets; i++)
        {
          if (iolets[ioletIDs[i]]->DoComms())
          {
            mComms[i]->FinishSend();
          }
        }
      }

      void BoundaryValues::FinishReceive()
      {
        for (int i = 0; i < nIOlets; i++)
        {
          if (iolets[ioletIDs[i]]->DoComms())
          {
            mComms[i]->Wait();
          }
        }
      }

      void BoundaryValues::Reset()
      {
        for (int i = 0; i < nTotIOlets; i++)
        {
          iolets[i]->ResetValues();
        }

        for (int i = 0; i < nIOlets; i++)
        {
          if (iolets[ioletIDs[i]]->DoComms())
          {
            if (IsCurrentProcTheBCProc())
            {
              iolets[ioletIDs[i]]->InitialiseCycle(density_cycle[i], mState);
              mComms[i]->Send(&density_cycle[i][0]);
            }
            else
            {
              mComms[i]->Receive(&iolets[ioletIDs[i]]->density);
            }
          }
          else
          {
            iolets[ioletIDs[i]]->InitialiseCycle(density_cycle[i], mState);
            iolets[ioletIDs[i]]->density = density_cycle[i][0];
          }
        }

        for (int i = 0; i < nIOlets; i++)
        {
          if (iolets[ioletIDs[i]]->DoComms())
          {
            mComms[i]->WaitAllComms();
          }
        }
      }

      // This assumes the program has already waited for comms to finish before
      distribn_t BoundaryValues::GetBoundaryDensity(const int index)
      {
        return iolets[index]->density;
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
