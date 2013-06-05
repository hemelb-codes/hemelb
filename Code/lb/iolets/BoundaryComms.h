// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_LB_BOUNDARIES_BOUNDARYCOMMS_H
#define HEMELB_LB_BOUNDARIES_BOUNDARYCOMMS_H

#include "geometry/LatticeData.h"
#include "lb/SimulationState.h"
#include "util/UnitConverter.h"

namespace hemelb
{
  namespace lb
  {
    namespace iolets
    {

      class BoundaryComms
      {
        public:
          BoundaryComms(SimulationState* iSimState, std::vector<int> &iProcsList, bool iHasBoundary);
          ~BoundaryComms();

          void Wait();

          // It is up to the caller to make sure only BCproc calls send
          void Send(distribn_t* density);
          void Receive(distribn_t* density);

          const std::vector<int>& GetListOfProcs() const
          {
            return procsList;
          }

          void ReceiveDoubles(double* double_array, int size);
          void WaitAllComms();
          void FinishSend();

        private:
          // This is necessary to support BC proc having fluid sites
          bool hasBoundary;

          // These are only assigned on the BC proc as it is the only one that needs to know
          // which proc has which IOlet
          int nProcs;
          std::vector<int> procsList;

          MPI_Request *sendRequest;
          MPI_Status *sendStatus;

          MPI_Request receiveRequest;
          MPI_Status receiveStatus;

          SimulationState* mState;
      };

    }
  }
}

#endif /* HEMELB_LB_BOUNDARIES_BOUNDARYCOMMS_H */
