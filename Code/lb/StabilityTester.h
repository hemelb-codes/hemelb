// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_LB_STABILITYTESTER_H
#define HEMELB_LB_STABILITYTESTER_H

#include "net/PhasedBroadcastRegular.h"
#include "geometry/LatticeData.h"

namespace hemelb
{
  namespace lb
  {
    /**
     * Class to repeatedly assess the stability of the simulation, using the PhasedBroadcast
     * interface.
     *
     * PhasedBroadcastRegular is used because we know in advance which iterations we will
     * want to communicate on. The default parameters suffice: no initial action is necessary
     * because we can assess the stability just before communicating (it doesn't have to happen
     * at the same time on all nodes), only one communication is needed between depths, which
     * can't overlap. We go down the tree to pass the overall stability to all nodes, and we go up
     * the tree to compose the local stability for all nodes to discover whether the simulation as
     * a whole is stable.
     */
    template<class LatticeType>
    class StabilityTester : public net::PhasedBroadcastRegular<>
    {
      public:
        StabilityTester(const geometry::LatticeData * iLatDat, net::Net* net,
                        SimulationState* simState, reporting::Timers& timings,
                        bool checkForConvergence, double relativeTolerance) :
            net::PhasedBroadcastRegular<>(net, simState, SPREADFACTOR), mLatDat(iLatDat),
                mSimState(simState), timings(timings), checkForConvergence(checkForConvergence),
                relativeTolerance(relativeTolerance)
        {

          Reset();
        }

        /**
         * Override the reset method in the base class, to reset the stability variables.
         */
        void Reset()
        {
          mUpwardsStability = UndefinedStability;
          mDownwardsStability = UndefinedStability;

          mSimState->SetStability(UndefinedStability);

          for (unsigned int ii = 0; ii < SPREADFACTOR; ii++)
          {
            mChildrensStability[ii] = UndefinedStability;
          }
        }

      protected:
        /**
         * Override the methods from the base class to propagate data from the root, and
         * to send data about this node and its childrens' stabilities up towards the root.
         */
        void ProgressFromChildren(unsigned long splayNumber)
        {
          ReceiveFromChildren<int>(mChildrensStability, 1);
        }

        void ProgressFromParent(unsigned long splayNumber)
        {
          ReceiveFromParent<int>(&mDownwardsStability, 1);
        }

        void ProgressToChildren(unsigned long splayNumber)
        {
          SendToChildren<int>(&mDownwardsStability, 1);
        }

        void ProgressToParent(unsigned long splayNumber)
        {
          SendToParent<int>(&mUpwardsStability, 1);
        }

        /**
         * The algorithm that checks distribution function convergence must be run in this
         * method rather than in ProgressToParent to make sure that the current timestep has
         * finished streaming.
         *
         * @param splayNumber
         */
        void PostSendToParent(unsigned long splayNumber)
        {
          timings[hemelb::reporting::Timers::monitoring].Start();

          // No need to bother testing out local lattice points if we're going to be
          // sending up a 'Unstable' value anyway.
          if (mUpwardsStability != Unstable)
          {
            bool unconvergedSitePresent = false;

            for (site_t i = 0; i < mLatDat->GetLocalFluidSiteCount(); i++)
            {
              for (unsigned int l = 0; l < LatticeType::NUMVECTORS; l++)
              {
                distribn_t value = *mLatDat->GetFNew(i * LatticeType::NUMVECTORS + l);

                // Note that by testing for value > 0.0, we also catch stray NaNs.
                if (! (value > 0.0))
                {
                  mUpwardsStability = Unstable;
                  break;
                }
              }
              ///@todo: If we refactor the previous loop out, we can get away with a single break statement
              if (mUpwardsStability == Unstable)
              {
                break;
              }

              if (checkForConvergence)
              {
                distribn_t relativeDifference =
                    ComputeRelativeDifference(mLatDat->GetFNew(i * LatticeType::NUMVECTORS),
                                              mLatDat->GetSite(i).GetFOld<LatticeType>());

                if (relativeDifference > relativeTolerance)
                {
                  // The simulation is stable but hasn't converged in the whole domain yet.
                  unconvergedSitePresent = true;
                }
              }
            }

            switch (mUpwardsStability)
            {
              case UndefinedStability:
              case Stable:
              case StableAndConverged:
                mUpwardsStability = (checkForConvergence && !unconvergedSitePresent) ?
                  StableAndConverged :
                  Stable;
                break;
              case Unstable:
                break;
            }
          }

          timings[hemelb::reporting::Timers::monitoring].Stop();
        }

        /**
         * Computes the relative difference between the densities at the beginning and end of a
         * timestep, i.e. |(rho_new - rho_old) / (rho_old - rho_0)|.
         *
         * @param fNew Distribution function after stream and collide, i.e. solution of the current timestep.
         * @param fOld Distribution function at the end of the previous timestep.
         * @return relative difference between the densities computed from fNew and fOld.
         */
        inline double ComputeRelativeDifference(const distribn_t* fNew,
                                                const distribn_t* fOld) const
        {
          distribn_t new_density = 0.;
          distribn_t old_density = 0.;
          for (unsigned int l = 0; l < LatticeType::NUMVECTORS; l++)
          {
            new_density += fNew[l];
            old_density += fOld[l];
          }

          // This is equivalent to REFERENCE_PRESSURE_mmHg in lattice units
          distribn_t ref_density = 1.0;

          if (old_density == ref_density)
          {
            // We want to avoid returning inf if the site is at pressure = REFERENCE_PRESSURE_mmHg
            return 0.0;
          }

          return fabs( (new_density - old_density) / (old_density - ref_density));
        }

        /**
         * Take the combined stability information (an int, with a value of hemelb::lb::Unstable
         * if any child node is unstable) and start passing it back down the tree.
         */
        void TopNodeAction()
        {
          mDownwardsStability = mUpwardsStability;
        }

        /**
         * Override the method from the base class to use the data from child nodes.
         */
        void PostReceiveFromChildren(unsigned long splayNumber)
        {
          timings[hemelb::reporting::Timers::monitoring].Start();

          // No need to test children's stability if this node is already unstable.
          if (mUpwardsStability != Unstable)
          {
            for (int ii = 0; ii < (int) SPREADFACTOR; ii++)
            {
              if (mChildrensStability[ii] == Unstable)
              {
                mUpwardsStability = Unstable;
                break;
              }
            }

            // If the simulation wasn't found to be unstable and we need to check for convergence, do it now.
            if ( (mUpwardsStability != Unstable) && checkForConvergence)
            {
              bool anyStableNotConverged = false;
              bool anyConverged = false;

              // mChildrensStability will contain UndefinedStability for non-existent children
              for (int ii = 0; ii < (int) SPREADFACTOR; ii++)
              {
                if (mChildrensStability[ii] == StableAndConverged)
                {
                  anyConverged = true;
                }

                if (mChildrensStability[ii] == Stable)
                {
                  anyStableNotConverged = true;
                }
              }

              // With the current configuration the root node of the tree won't own any fluid sites. Its
              // state only depends on children nodes not on local state.
              if (anyConverged
                  && (mUpwardsStability == StableAndConverged || GetParent() == NOPARENT))
              {
                mUpwardsStability = StableAndConverged;
              }

              if (anyStableNotConverged)
              {
                mUpwardsStability = Stable;
              }
            }
          }

          timings[hemelb::reporting::Timers::monitoring].Stop();
        }

        /**
         * Apply the stability value sent by the root node to the simulation logic.
         */
        void Effect()
        {
          mSimState->SetStability((Stability) mDownwardsStability);
        }

      private:
        /**
         * Slightly arbitrary spread factor for the tree.
         */
        static const unsigned int SPREADFACTOR = 10;

        const geometry::LatticeData * mLatDat;

        /**
         * Stability value of this node and its children to propagate upwards.
         */
        int mUpwardsStability;
        /**
         * Stability value as understood by the root node, to pass downwards.
         */
        int mDownwardsStability;
        /**
         * Array for storing the passed-up stability values from child nodes.
         */
        int mChildrensStability[SPREADFACTOR];
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
    };
  }
}

#endif /* HEMELB_LB_STABILITYTESTER_H */
