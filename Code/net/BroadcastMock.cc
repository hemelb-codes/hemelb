#include "BroadcastMock.h"

namespace hemelb
{

  namespace net
  {

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
            dataStart[childIndex * countPerChild] = 10.0;
            dataStart[childIndex * countPerChild + 1] = 15.0;
          }
          break;

        case 1:
          for (unsigned childIndex = 0; childIndex < spreadFactor - 1; childIndex++)
          {
            dataStart[childIndex * countPerChild] = 10.0;
            dataStart[childIndex * countPerChild + 1] = 15.0;
          }
          dataStart[ (spreadFactor - 1) * countPerChild] = 1.0;
          dataStart[ (spreadFactor - 1) * countPerChild + 1] = 100.0;
          break;

        case 2:
          for (unsigned childIndex = 0; childIndex < spreadFactor; childIndex++)
          {
            dataStart[childIndex * countPerChild] = 10.0;
            dataStart[childIndex * countPerChild + 1] = 15.0;
          }
          break;

        default:
          /// TODO: control should never reach this branch. Exception?
          assert(false);
      }
    }

    // Explicit instantiation
    template void BroadcastMock::ReceiveFromChildren<distribn_t>(distribn_t*, int);
  }
}
