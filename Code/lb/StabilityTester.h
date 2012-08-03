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
        StabilityTester(const geometry::LatticeData * iLatDat,
                        net::Net* net,
                        SimulationState* simState,
                        reporting::Timers& timings) :
            net::PhasedBroadcastRegular<>(net, simState, SPREADFACTOR), mLatDat(iLatDat), mSimState(simState), timings(timings)
        {
          Reset();
        }

        /**
         * Override the reset method in the base class, to reset the stability variables.
         */
        void Reset()
        {
          // Re-initialise all values to be Stable.
          mUpwardsStability = Stable;
          mDownwardsStability = Stable;

          mSimState->SetStability(Stable);

          for (unsigned int ii = 0; ii < SPREADFACTOR; ii++)
          {
            mChildrensStability[ii] = Stable;
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
          timings[hemelb::reporting::Timers::monitoring].Start();

          // No need to bother testing out local lattice points if we're going to be
          // sending up a 'Unstable' value anyway.
          if (mUpwardsStability != Unstable)
          {
            for (site_t i = 0; i < mLatDat->GetLocalFluidSiteCount(); i++)
            {
              for (unsigned int l = 0; l < LatticeType::NUMVECTORS; l++)
              {
                distribn_t value = *mLatDat->GetFNew(i * LatticeType::NUMVECTORS + l);

                // Note that by testing for value > 0.0, we also catch stray NaNs.
                if (! (value > 0.0))
                {
                  mUpwardsStability = Unstable;
                }
              }
            }

            timings[hemelb::reporting::Timers::monitoring].Stop();
          }

          SendToParent<int>(&mUpwardsStability, 1);
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
          if (mUpwardsStability == Stable)
          {
            for (int ii = 0; ii < (int) SPREADFACTOR; ii++)
            {
              if (mChildrensStability[ii] == Unstable)
              {
                mUpwardsStability = Unstable;
                break;
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
    };
  }
}

#endif /* HEMELB_LB_STABILITYTESTER_H */
