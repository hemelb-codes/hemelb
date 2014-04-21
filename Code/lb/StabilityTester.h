
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_STABILITYTESTER_H
#define HEMELB_LB_STABILITYTESTER_H

#include "net/phased/Concern.h"
#include "net/phased/steps.h"
#include "geometry/LatticeData.h"

namespace hemelb
{
  namespace lb
  {
    /**
     * Class to repeatedly assess the stability of the simulation, using non-blocking collective
     */
    template<class LatticeType>
    class StabilityTester : public net::phased::Concern
    {
      public:
        StabilityTester(const geometry::LatticeData * iLatDat, net::Net* net,
                        SimulationState* simState, reporting::Timers& timings,
                        bool checkForConvergence, double relativeTolerance);


        /**
         * Override the reset method in the base class, to reset the stability variables.
         */
        void Reset();

        bool CallAction(int action);

      protected:
        /**
         * Compute the local stability/convergence state.
         */
        void PreSend(void);

        /**
         * Initiate the collective.
         */
        void Send(void);

        /**
         * Wait on the collectives to finish.
         */
        void Wait(void);

        /**
         * Computes the relative difference between the densities at the beginning and end of a
         * timestep, i.e. |(rho_new - rho_old) / (rho_old - rho_0)|.
         *
         * @param fNew Distribution function after stream and collide, i.e. solution of the current timestep.
         * @param fOld Distribution function at the end of the previous timestep.
         * @return relative difference between the densities computed from fNew and fOld.
         */
        double ComputeRelativeDifference(const distribn_t* fNew,
                                         const distribn_t* fOld) const;

        /**
         * Apply the stability value sent by the root node to the simulation logic.
         */
        void Effect();

      private:

        const geometry::LatticeData * mLatDat;

        /**
         * Local and global stability.
         */
        int localStability;
        int globalStability;

        /**
         * Pointer to the simulation state used in the rest of the simulation.
         */
        lb::SimulationState* mSimState;

        /** Timing object. */
        reporting::Timers& timings;

        /** Whether to check for steady flow simulation convergence */
        bool checkForConvergence;

        /** Relative error tolerance in convergence check */
        double relativeTolerance;
        /**
         * Private communicator for non-blocking collectives.
         */
        net::MpiCommunicator collectiveComm;
        /**
         * Request object for the collective
         */
        MPI_Request collectiveReq;
    };
  }
}

#include "lb/StabilityTester.hpp"
#endif /* HEMELB_LB_STABILITYTESTER_H */
