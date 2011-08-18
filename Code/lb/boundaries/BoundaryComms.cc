#include "lb/boundaries/BoundaryComms.h"
#include "topology/NetworkTopology.h"
#include "util/utilityFunctions.h"
#include <math.h>

namespace hemelb
{
  namespace lb
  {
    namespace boundaries
    {

      BoundaryComms::BoundaryComms(SimulationState* iSimState,
                                   int iTotIOlets,
                                   geometry::LatticeData::SiteType IOtype,
                                   geometry::LatticeData* iLatDat) :
        nTotIOlets(iTotIOlets), mState(iSimState)
      {
        FindBCProcRank();

        FindIOlets(IOtype, iLatDat);

        if (IsCurrentProcTheBCProc())
        {
          nProcs = new int[nTotIOlets];
          procsList = new int*[nTotIOlets];

          for (int i = 0; i < nTotIOlets; i++)
            nProcs[i] = 0;
        }

        // Now BC process must find out which process belongs to what group
        GatherProcList();

        if (IsCurrentProcTheBCProc())
        {
          int nRequests = 0;

          requestOffset = new int[nTotIOlets];

          for (int i = 0; i < nTotIOlets; i++)
          {
            requestOffset[i] = nRequests;
            nRequests += nProcs[i];
          }

          request = new MPI_Request[nRequests];
          status = new MPI_Status[nRequests];
        }
        else
        {
          request = new MPI_Request[nTotIOlets];
          status = new MPI_Status[nTotIOlets];
        }
      }

      void BoundaryComms::FindBCProcRank()
      {
        proc_t BCrank = 0;

        if (IsCurrentProcTheBCProc())
          BCrank = topology::NetworkTopology::Instance()->GetLocalRank();

        // Since only one proc will update BCrank, the sum of all BCrank is the BCproc
        MPI_Allreduce(&BCrank, &BCproc, 1, hemelb::MpiDataType(BCrank), MPI_SUM, MPI_COMM_WORLD);
      }

      void BoundaryComms::FindIOlets(geometry::LatticeData::SiteType IOtype,
                                     geometry::LatticeData* iLatDat)
      {
        nIOlets = 0;
        IOlets = std::vector<int>(0);

        for (site_t i = 0; i < iLatDat->GetLocalFluidSiteCount(); i++)
        {
          if (iLatDat->GetSiteType(i) == IOtype
              && !util::VectorFunctions::member(IOlets, iLatDat->GetBoundaryId(i)))
          {
            nIOlets++;
            IOlets.push_back(iLatDat->GetBoundaryId(i));
          }
        }
      }

      void BoundaryComms::GatherProcList()
      {
        if (IsCurrentProcTheBCProc())
        {
          // These should be bool, but MPI only supports MPI_INT
          // For each inlet/outlet there is an array of length equal to total number of procs.
          // Each stores true/false value. True if proc of rank equal to the index contains
          // the given inlet/outlet.
          int **boolList = new int*[nTotIOlets];

          int nTotProcs = topology::NetworkTopology::Instance()->GetProcessorCount();

          for (int i = 0; i < nTotIOlets; i++)
          {
            boolList[i] = new int[nTotProcs];
          }

          MPI_Status tempStat;

          // Non-BC procs will be sending at this point
          for (int i = 0; i < nTotIOlets; i++)
          {
            for (int proc = 0; proc < nTotProcs; proc++)
            {
              if (proc != BCproc)
              {
                MPI_Recv(&boolList[i][proc], 1, MPI_INT, proc, 100, MPI_COMM_WORLD, &tempStat);
              }
              else
              {
                // All of them should be false, but in case we ever want BCproc to also have some sites on board
                boolList[i][proc] = util::VectorFunctions::member(IOlets, i);
              }
            }
          }

          // Now we have an array for each IOlet with true (1) at indices corresponding to
          // processes that are members of that group. We have to convert this into arrays
          // of ints which store a list of processor ranks.

          for (int i = 0; i < nTotIOlets; i++)
          {
            for (int j = 0; j < nTotProcs; j++)
            {
              if (boolList[i][j])
                nProcs[i]++;
            }

            procsList[i] = new int[nProcs[i]];

            int memberIndex = 0;

            for (int j = 0; j < nTotProcs; j++)
            {
              if (boolList[i][j])
              {
                procsList[i][memberIndex] = j;
                memberIndex++;
              }
            }
          }

          // Clear up
          for (int i = 0; i < nTotIOlets; i++)
          {
            delete[] boolList[i];
          }

          delete[] boolList;

        }
        else
        {
          // This is where the info about whether a proc contains a given inlet/outlet is sent
          // If it does contain the given inlet/outlet it sends a true value, else it sends a false.
          for (int i = 0; i < nTotIOlets; i++)
          {
            int IOletOnThisProc = util::VectorFunctions::member(IOlets, i); // true if inlet i is on this proc

            MPI_Ssend(&IOletOnThisProc, 1, MPI_INT, BCproc, 100, MPI_COMM_WORLD);
          }
        }
      }

      BoundaryComms::~BoundaryComms()
      {

        if (IsCurrentProcTheBCProc())
        {
          for (int i = 0; i < nTotIOlets; i++)
            delete[] procsList[i];

          delete[] requestOffset;
          delete[] procsList;
          delete[] nProcs;
        }

        // Communicators and groups
        delete[] request;
        delete[] status;
      }

      inline bool BoundaryComms::IsCurrentProcTheBCProc()
      {
        return topology::NetworkTopology::Instance()->IsCurrentProcTheIOProc();
      }

      proc_t BoundaryComms::GetBCProcRank()
      {
        return BCproc;
      }

      // This assumes that the BCproc never receives
      void BoundaryComms::Wait(const int index)
      {
        MPI_Wait(&request[index], &status[index]);
      }

      void BoundaryComms::WaitAllComms()
      {
        // Now wait for all to complete
        if (IsCurrentProcTheBCProc())
        {
          for (int i = 0; i < nTotIOlets; i++)
          {
            MPI_Waitall(nProcs[i], &request[requestOffset[i]], &status[requestOffset[i]]);
          }
        }
        else
        {
          for (int i = 0; i < nIOlets; i++)
          {
            MPI_Wait(&request[IOlets[i]], &status[IOlets[i]]);
          }
        }
      }

      void BoundaryComms::SendAndReceive(distribn_t* density)
      {
        if (IsCurrentProcTheBCProc())
        {
          int message = 0;

          for (int i = 0; i < nTotIOlets; i++)
          {
            for (int proc = 0; proc < nProcs[i]; proc++)
            {
              MPI_Isend(&density[i],
                        1,
                        hemelb::MpiDataType(density[0]),
                        procsList[i][proc],
                        100,
                        MPI_COMM_WORLD,
                        &request[message++]);
            }
          }
        }
        else
        {
          for (int i = 0; i < nIOlets; i++)
          {
            MPI_Irecv(&density[IOlets[i]],
                      1,
                      hemelb::MpiDataType(density[0]),
                      BCproc,
                      100,
                      MPI_COMM_WORLD,
                      &request[IOlets[i]]);
          }
        }
      }

      void BoundaryComms::FinishSend()
      {
        // Don't move on to next step with BC proc until all messages have been sent
        // Precautionary measure to make sure proc doesn't overwrite, before message is sent
        if (IsCurrentProcTheBCProc())
        {
          for (int i = 0; i < nTotIOlets; i++)
            MPI_Waitall(nProcs[i], &request[requestOffset[i]], &status[requestOffset[i]]);
        }
      }

    }
  }
}
