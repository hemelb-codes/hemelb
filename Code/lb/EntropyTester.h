#ifndef HEMELB_LB_ENTROPYTESTER_H
#define HEMELB_LB_ENTROPYTESTER_H

#include "net/PhasedBroadcastRegular.h"
#include "geometry/LatticeData.h"

namespace hemelb
{
  namespace lb
  {

    class EntropyTester : public net::PhasedBroadcastRegular<false, 1, 1, false, true>
    {
      public:
        EntropyTester(int* collisionTypes,
                      unsigned int typesTested,
                      const geometry::LatticeData * iLatDat,
                      net::Net* net,
                      SimulationState* simState);

        ~EntropyTester();

        void PreReceive();

        /**
         * Override the reset method in the base class, to reset the stability variables.
         */
        void Reset();

      protected:
        /**
         * Override the methods from the base class to propagate data from the root, and
         * to send data about this node and its childrens' stabilities up towards the root.
         */
        void ProgressFromChildren(unsigned long splayNumber);
        void ProgressToParent(unsigned long splayNumber);

        /**
         * Take the combined stability information (an int, with a value of hemelb::lb::Unstable
         * if any child node is unstable) and start passing it back down the tree.
         */
        void TopNodeAction();

        /**
         * Override the method from the base class to use the data from child nodes.
         */
        void PostReceiveFromChildren(unsigned long splayNumber);

      private:
        enum HTHEOREM
        {
          OBEYED,
          DISOBEYED
        };

        /**
         * Slightly arbitrary spread factor for the tree.
         */
        static const unsigned int SPREADFACTOR = 10;

        const geometry::LatticeData * mLatDat;

        /**
         * Stability value of this node and its children to propagate upwards.
         */
        int mUpwardsValue;
        /**
         * Array for storing the passed-up stability values from child nodes.
         */
        int mChildrensValues[SPREADFACTOR];

        bool mCollisionTypesTested[COLLISION_TYPES];
        double* mHPreCollision;
    };

  }
}

#endif /* HEMELB_LB_ENTROPYTESTER_H */
