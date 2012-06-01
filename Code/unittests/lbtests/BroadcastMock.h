#ifndef HEMELB_UNITTESTS_LBTESTS_BROADCASTMOCK_H
#define HEMELB_UNITTESTS_LBTESTS_BROADCASTMOCK_H

#include "net/PhasedBroadcastRegular.h"
#include "lb/SimulationState.h"
#include "net/net.h"

namespace hemelb
{
  namespace net
  {
    /**
     * In this mock, we pretend that the current process is the root node of a phased
     * broadcast and that a pair of values is going up the tree. The mock simulates
     * three rounds of communications:
     *  1) All children report (14,15).
     *  2) One child reports (1,100) and the rest (14,15).
     *  3) All back to (14,15).
     */
    class BroadcastMock : public net::PhasedBroadcastRegular<>
    {

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

    /// TODO: Split file into .h and .cc ?
    BroadcastMock::BroadcastMock(net::Net * net,
                                 const lb::SimulationState * simState,
                                 unsigned int spreadFactor) :
      net::PhasedBroadcastRegular<>(net, simState, spreadFactor), spreadFactor(spreadFactor),
          callCounter(0)
    {
    }

    BroadcastMock::~BroadcastMock()
    {
    }

    void BroadcastMock::RequestComms()
    {
      // Action taken when a node has to receive from its children in the tree. Calls ReceiveFromChildren mock below.
      ProgressFromChildren(0);
      // Action taken when a node has to send to its parent in the tree.
      ProgressToParent(0);
    }

    void BroadcastMock::PostReceive()
    {
      // Action taken after data has been received from its children
      PostReceiveFromChildren(0);

      // Action taken after data has been sent to its parent
      PostSendToParent(0);
    }

    void BroadcastMock::EndIteration()
    {
      // Action taken when upwards-travelling data reaches the top node
      TopNodeAction();

      // Action taken by all nodes when downwards-travelling data has been sent to every node.
      Effect();
    }

    template<class T>
    void BroadcastMock::ReceiveFromChildren(T* dataStart, int countPerChild)
    {
      assert(countPerChild == 2);

      switch (callCounter++)
      {
        case 0:
          for (unsigned childIndex = 0; childIndex < spreadFactor; childIndex++)
          {
            dataStart[childIndex * countPerChild] = 14.0;
            dataStart[childIndex * countPerChild + 1] = 15.0;
          }
          break;

        case 1:
          for (unsigned childIndex = 0; childIndex < spreadFactor - 1; childIndex++)
          {
            dataStart[childIndex * countPerChild] = 14.0;
            dataStart[childIndex * countPerChild + 1] = 15.0;
          }
          dataStart[ (spreadFactor - 1) * countPerChild] = 1.0;
          dataStart[ (spreadFactor - 1) * countPerChild + 1] = 100.0;
          break;

        case 2:
          for (unsigned childIndex = 0; childIndex < spreadFactor; childIndex++)
          {
            dataStart[childIndex * countPerChild] = 14.0;
            dataStart[childIndex * countPerChild + 1] = 15.0;
          }
          break;

        default:
          /// TODO: control should never reach this branch. Exception?
          assert(false);
      }
    }

  }
}

#endif /* HEMELB_UNITTESTS_LBTESTS_BROADCASTMOCK_H */
