#ifndef HEMELB_LB_BOUNDARIES_BOUNDARYCOMMS_H
#define HEMELB_LB_BOUNDARIES_BOUNDARYCOMMS_H

#include "geometry/LatticeData.h"
#include "lb/SimulationState.h"
#include "util/UnitConverter.h"

namespace hemelb
{
  namespace lb
  {
    namespace boundaries
    {

      class BoundaryComms
      {
        public:
          BoundaryComms(SimulationState* iSimState,
                        std::vector<int> &iProcsList,
                        bool iHasBoundary,
                        proc_t iBCproc);
          ~BoundaryComms();

          void Wait();

          // It is up to the caller to make sure only BCproc calls send
          void Send(distribn_t* density);
          void Receive(distribn_t* density);
          void WaitAllComms();
          void FinishSend();

        private:
          proc_t BCproc; // Process responsible for sending out BC info

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
