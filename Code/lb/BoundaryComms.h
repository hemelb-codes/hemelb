#ifndef HEMELB_LB_BOUNDARYCOMMS_H
#define HEMELB_LB_BOUNDARYCOMMS_H

#include "net/IteratedAction.h"
#include "geometry/LatticeData.h"
#include "lb/SimulationState.h"
#include "util/UnitConverter.h"

namespace hemelb
{
  namespace lb
  {

    /*
     * This class deals with in/outlet boundary conditions. This class is not concerned with updating
     * the densities. It assumes that it has already done when Broadcast is called. Every in/outlet will
     * have a corresponding communicator which includes the BCproc and all processes that contain the
     * given in/outlet. Once per cycle BroadcastBoundaryDensities needs to be called so that BCproc
     * updates the relevant processes with the new values.
     */
    class BoundaryComms : public net::IteratedAction
    {
      public:
        BoundaryComms(SimulationState* iSimState, int iTotIOlets);
        ~BoundaryComms();

        void Initialise(geometry::LatticeData::SiteType IOtype,
                        geometry::LatticeData* iLatDat,
                        std::vector<distribn_t>* iDensityCycleVector);

        void RequestComms();
        void EndIteration();
        void Reset();

        distribn_t GetBoundaryDensity(const int index);
        void WaitAllComms();

      private:
        proc_t BCproc; // Process responsible for sending out BC info
        std::vector<distribn_t>* density_cycle_vector;
        distribn_t *density_cycle;
        distribn_t *density;

        // Total number of inlets/outlets in simulation
        int nTotIOlets;

        // Number of inlets/outlets on this process
        int nIOlets;

        // List of indices of inlets/outlets on this process
        std::vector<int> IOlets;

        // These are only assigned on the BC proc as it is the only one that needs to know
        // which proc has which IOlet
        int* requestOffset;
        int* nProcs;
        int** procsList;

        MPI_Request* request;
        MPI_Status* status;

        SimulationState* mState;
    };

  }
}

#endif /* HEMELB_LB_BOUNDARYCOMMS_H */
