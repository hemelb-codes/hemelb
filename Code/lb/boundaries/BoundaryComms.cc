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
                                   int iProcs,
                                   int* iProcsList,
                                   bool iHasBoundary,
                                   proc_t iBCproc) :
        BCproc(iBCproc), hasBoundary(iHasBoundary), nProcs(iProcs), procsList(iProcsList),
            mState(iSimState)
      {
        if (IsCurrentProcTheBCProc())
        {
          sendRequest = new MPI_Request[nProcs];
          sendStatus = new MPI_Status[nProcs];
        }
      }

      BoundaryComms::~BoundaryComms()
      {

        if (IsCurrentProcTheBCProc())
        {
          delete[] sendRequest;
          delete[] sendStatus;
        }
      }

      inline bool BoundaryComms::IsCurrentProcTheBCProc()
      {
        return topology::NetworkTopology::Instance()->GetLocalRank() == BCproc;
      }

      void BoundaryComms::Wait()
      {
        MPI_Wait(&request, &status);
      }

      void BoundaryComms::WaitAllComms()
      {
        // Now wait for all to complete
        if (IsCurrentProcTheBCProc())
        {
          MPI_Waitall(nProcs, sendRequest, sendStatus);

          if (hasBoundary)
            MPI_Wait(&request, &status);
        }
        else
        {
          MPI_Wait(&request, &status);
        }

      }

      void BoundaryComms::SendAndReceive(distribn_t* density)
      {
        if (IsCurrentProcTheBCProc())
        {
          for (int proc = 0; proc < nProcs; proc++)
          {
            MPI_Isend(density,
                      1,
                      hemelb::MpiDataType(*density),
                      procsList[proc],
                      100,
                      MPI_COMM_WORLD,
                      &sendRequest[proc]);

            if (hasBoundary)
              MPI_Irecv(density,
                        1,
                        hemelb::MpiDataType(*density),
                        BCproc,
                        100,
                        MPI_COMM_WORLD,
                        &request);
          }
        }
        else
        {
          MPI_Irecv(density,
                    1,
                    hemelb::MpiDataType(*density),
                    BCproc,
                    100,
                    MPI_COMM_WORLD,
                    &request);
        }
      }

      void BoundaryComms::FinishSend()
      {
        // Don't move on to next step with BC proc until all messages have been sent
        // Precautionary measure to make sure proc doesn't overwrite, before message is sent
        if (IsCurrentProcTheBCProc())
        {
          MPI_Waitall(nProcs, sendRequest, sendStatus);
        }
      }

    }
  }
}
