#ifndef HEMELB_LB_STABILITYCHECKER_H
#define HEMELB_LB_STABILITYCHECKER_H

#include "net/PhasedBroadcast.h"
#include "geometry/LatticeData.h"

namespace hemelb
{
  namespace lb
  {
    class StabilityTester : public net::PhasedBroadcast
    {
      public:
        StabilityTester(const geometry::LatticeData * iLatDat,
                        net::Net* net,
                        topology::NetworkTopology *iNetTop,
                        SimulationState* simState);

        /**
         * Override the reset method in the base class, to reset the stability variables.
         */
        void Reset();

      protected:
        /**
         * Override the methods from the base class to propagate data from the root, and
         * to send data about this node and its childrens' stabilities up towards the root.
         */
        void ProgressFromChildren();
        void ProgressFromParent();
        void ProgressToChildren();
        void ProgressToParent();

        /**
         * Take the combined stability information (an int, with a value of hemelb::lb::Unstable
         * if any child node is unstable) and start passing it back down the tree.
         */
        void TopNodeAction();

        /**
         * Override the method from the base class to use the data from child nodes.
         */
        void PostReceiveFromChildren();

        /**
         * Apply the stability value sent by the root node to the simulation logic.
         */
        void Effect();

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
         * Pointer to the stability variable used in the rest of the simulation.
         */
        int * mPublicSimulationStability;
    };
  }
}

#endif /* HEMELB_LB_STABILITYCHECKER_H */
