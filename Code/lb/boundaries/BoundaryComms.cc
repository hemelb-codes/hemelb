#include "lb/boundaries/BoundaryComms.h"
#include "lb/boundaries/BoundaryValues.h"
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
                                   std::vector<int> &iProcsList,
                                   bool iHasBoundary,
                                   proc_t iBCproc) :
          BCproc(iBCproc), hasBoundary(iHasBoundary), nProcs((int) iProcsList.size()), procsList(iProcsList), mState(iSimState)
      {
        if (BoundaryValues::IsCurrentProcTheBCProc())
        {
          sendRequest = new MPI_Request[nProcs];
          sendStatus = new MPI_Status[nProcs];
        }
      }

      BoundaryComms::~BoundaryComms()
      {

        if (BoundaryValues::IsCurrentProcTheBCProc())
        {
          delete[] sendRequest;
          delete[] sendStatus;
        }
      }

      void BoundaryComms::Wait()
      {
        if (hasBoundary || !BoundaryValues::IsCurrentProcTheBCProc())
        {
          MPI_Wait(&receiveRequest, &receiveStatus);
        }
      }

      void BoundaryComms::WaitAllComms()
      {
        // Now wait for all to complete
        if (BoundaryValues::IsCurrentProcTheBCProc())
        {
          MPI_Waitall(nProcs, sendRequest, sendStatus);

          if (hasBoundary)
            MPI_Wait(&receiveRequest, &receiveStatus);
        }
        else
        {
          MPI_Wait(&receiveRequest, &receiveStatus);
        }

      }

      void BoundaryComms::SendAndReceive(distribn_t* density)
      {
        if (BoundaryValues::IsCurrentProcTheBCProc())
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
            {
              MPI_Irecv(density,
                        1,
                        hemelb::MpiDataType(*density),
                        BCproc,
                        100,
                        MPI_COMM_WORLD,
                        &receiveRequest);
            }
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
                    &receiveRequest);
        }
      }

      void BoundaryComms::FinishSend()
      {
        // Don't move on to next step with BC proc until all messages have been sent
        // Precautionary measure to make sure proc doesn't overwrite, before message is sent
        if (BoundaryValues::IsCurrentProcTheBCProc())
        {
          MPI_Waitall(nProcs, sendRequest, sendStatus);
        }
      }

    }
  }
}
