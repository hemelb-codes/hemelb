#ifndef HEMELB_LB_BOUNDARYCOMMS_H
#define HEMELB_LB_BOUNDARYCOMMS_H

#include "geometry/LatticeData.h"

namespace hemelb
{
  namespace lb
  {

    /*
     * This class deals with in/outlet boundary conditions. UpdateBoundaryConditions will calculate/obtain
     * all inlet/outlet boundary values. It should only be called by BCproc. Every in/outlet will
     * have a corresponding communicator which includes the BCproc and all processes that contain the
     * given in/outlet. Once per cycle BroadcastBoundaryDensities needs to be called so that BCproc
     * updates the relevant processes with the new values.
     */
    class BoundaryComms
    {
      public:
        BoundaryComms(const geometry::LatticeData* iLatDat, const SimConfig* iSimConfig);
        ~BoundaryComms();

        void UpdateBoundaryDensities(unsigned long time_step,
                                     unsigned long timeStepsPerCycle,
                                     distribn_t* inlet_density,
                                     distribn_t* inlet_density_avg,
                                     distribn_t* inlet_density_amp,
                                     distribn_t* inlet_density_phs,
                                     distribn_t* outlet_density,
                                     distribn_t* outlet_density_avg,
                                     distribn_t* outlet_density_amp,
                                     distribn_t* outlet_density_phs);

        void BroadcastBoundaryDensities(distribn_t* inlet_density, distribn_t* outlet_density);

      private:
        proc_t BCproc; // Process responsible for sending out BC info

        // Total number of inlets/outlets in simulation
        int nTotInlets;
        int nTotOutlets;

        // Number of inlets/outlets on this process
        int nInlets;
        int nOutlets;

        // List of indices of inlets/outlets on this process
        std::vector<int> inlets;
        std::vector<int> outlets;

        // Communicators and groups
        MPI_Comm* outlet_comms;
        MPI_Comm* inlet_comms;
        MPI_Group* outlet_groups;
        MPI_Group* inlet_groups;

        // Just for testing ATM
        void printStuff();
    };

  }
}

#endif /* HEMELB_LB_BOUNDARYCOMMS_H */
