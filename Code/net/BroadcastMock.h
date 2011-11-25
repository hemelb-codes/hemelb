#ifndef HEMELB_NET_BROADCASTMOCK_H
#define HEMELB_NET_BROADCASTMOCK_H

#include "net/PhasedBroadcastRegular.h"
#include "lb/SimulationState.h"
#include "net/net.h"

namespace hemelb
{
  namespace net
  {
    class BroadcastMock : public net::PhasedBroadcastRegular<>
    {
        /*
         * In this mock, we pretend that the current process is the root node of a phased
         * broadcast and that a pair of values is going up the tree. The mock simulates
         * three rounds of communications:
         *  1) All children report (10,15).
         *  2) One child reports (1,100) and the rest (10,15).
         *  3) All back to (10,15).
         */

      public:
        BroadcastMock(net::Net * net,
                      const lb::SimulationState * simState,
                      unsigned int spreadFactor);

        virtual ~BroadcastMock();

        /*
         *  Overwritten IteraredAction methods that implement the mock.
         */
        void RequestComms();
        void PostReceive();
        void EndIteration();

      protected:
        /**
         * Receives data from each child. This is a set length per child, and each child's
         * data is inserted contiguously into the provided array.
         *
         * @param dataStart Pointer to the start of the array.
         * @param countPerChild Number of elements to receive per child.
         */
        template<class T>
        void ReceiveFromChildren(T* dataStart, int countPerChild);

      private:
        /** Number of children nodes in the tree */
        unsigned spreadFactor;

        /** Number of times that the phased broadcast has been completed */
        unsigned callCounter;
    };
  }
}

#endif /* HEMELB_NET_BROADCASTMOCK_H */

