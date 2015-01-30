// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#include "lb/iolets/BoundaryComms.h"
#include "lb/iolets/BoundaryValues.h"
#include "net/IOCommunicator.h"
#include "util/utilityFunctions.h"

namespace hemelb
{
  namespace lb
  {
    namespace iolets
    {

      BoundaryComms::BoundaryComms(SimulationState* iSimState, std::vector<int> &iProcsList, const BoundaryCommunicator& boundaryComm, bool iHasBoundary) :
          hasBoundary(iHasBoundary), nProcs((int) iProcsList.size()), procsList(iProcsList), bcComm(boundaryComm)
      {
        /* iProcsList contains the procs containing said Boundary/iolet, but NOT proc 0 (the BoundaryControlling/BC proc)! */

        // Only BC proc sends
        if (bcComm.IsCurrentProcTheBCProc())
        {
          sendRequest = new MPI_Request[nProcs];
          sendStatus = new MPI_Status[nProcs];
        }
        else
        {
          sendRequest = nullptr;
          sendStatus = nullptr;
        }
      }

      BoundaryComms::~BoundaryComms()
      {

        if (bcComm.IsCurrentProcTheBCProc())
        {
          delete[] sendRequest;
          delete[] sendStatus;
        }
      }

      void BoundaryComms::Wait()
      {
        if (hasBoundary)
        {
          HEMELB_MPI_CALL(
              MPI_Wait, (&receiveRequest, &receiveStatus)
          );
        }
      }

      void BoundaryComms::WaitAllComms()
      {
        // Now wait for all to complete
        if (bcComm.IsCurrentProcTheBCProc())
        {
          HEMELB_MPI_CALL(
              MPI_Waitall, (nProcs, sendRequest, sendStatus)
          );

          if (hasBoundary)
            HEMELB_MPI_CALL(
                MPI_Wait, (&receiveRequest, &receiveStatus)
            );
        }
        else
        {
          HEMELB_MPI_CALL(
              MPI_Wait, (&receiveRequest, &receiveStatus)
          );
        }

      }

      // It is up to the caller to make sure only BCproc calls send
      void BoundaryComms::Send(distribn_t* density)
      {
        for (int proc = 0; proc < nProcs; proc++)
        {
          HEMELB_MPI_CALL(
              MPI_Isend, (
                  density,
                  1,
                  net::MpiDataType(*density),
                  procsList[proc],
                  100,
                  bcComm,
                  &sendRequest[proc]
              ));
        }
      }

      void BoundaryComms::Receive(distribn_t* density)
      {
        if (hasBoundary)
        {
          HEMELB_MPI_CALL(
              MPI_Irecv, (
                  density,
                  1,
                  net::MpiDataType(*density),
                  bcComm.GetBCProcRank(),
                  100,
                  bcComm,
                  &receiveRequest
              ));
        }
      }

      void BoundaryComms::FinishSend()
      {
        // Don't move on to next step with BC proc until all messages have been sent
        // Precautionary measure to make sure proc doesn't overwrite, before message is sent
        if (bcComm.IsCurrentProcTheBCProc())
        {
          HEMELB_MPI_CALL(
              MPI_Waitall, (nProcs, sendRequest, sendStatus)
          );
        }
      }

    }
  }
}
