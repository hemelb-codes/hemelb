
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_STABILITYTESTER_H
#define HEMELB_LB_STABILITYTESTER_H

#include "timestep/CollectiveActor.h"
#include "geometry/LatticeData.h"

namespace hemelb
{
  namespace lb
  {
    /**
     * Class to repeatedly assess the stability of the simulation, using non-blocking collective
     */
    template<class LatticeType>
    class StabilityTester : public timestep::CollectiveActor
    {
    public:
        StabilityTester(const geometry::LatticeData * iLatDat, comm::Communicator::ConstPtr comms,
                        SimulationState* simState, reporting::Timers& timings,
                        bool checkForConvergence, double relativeTolerance);

        ~StabilityTester() = default;

        /**
         * Override the reset method in the base class, to reset the stability variables.
         */
        void Reset();

    private:
      // Actor interface
      // (This is a template so inline is implied)

      // Nothing to do
      virtual void BeginAll() override {}
      // Nothing to do
      virtual void Begin() override {}
      // Receive is defined in CollectiveActor
      // virtual void Receive();
      // Compute the local stability/convergence state
      virtual void PreSend() override;
      // Initiate the collective
      virtual void Send() override;
      // Nothing to do
      virtual void PreWait() override {}
      // Wait is defined in CollectiveActor
      // virtual void Wait();
      // Apply the stability value sent by the root node to the simulation logic.
      virtual void End() override;
      // Nothing to do
      virtual void EndAll() override {}

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

        const geometry::LatticeData* mLatDat;

        /**
         * Local and global stability.
         */
        int localStability;
        int globalStability;

        /**
         * Pointer to the simulation state used in the rest of the simulation.
         */
        lb::SimulationState* mSimState;

        /** Whether to check for steady flow simulation convergence */
        bool checkForConvergence;

        /** Relative error tolerance in convergence check */
        double relativeTolerance;

        reporting::Timer& workTimer;
    };
  }
}

#include "lb/StabilityTester.hpp"
#endif /* HEMELB_LB_STABILITYTESTER_H */
