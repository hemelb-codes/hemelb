// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "tests/lb/BroadcastMocks.h"

// This file defines two PhasedBroadcastRegular mocks. One that
// behaves like the root node of the broadcast tree and another for a
// leaf node.

namespace hemelb
{
  namespace net
  {

    BroadcastMockRootNode::BroadcastMockRootNode(net::Net * net,
                                                 const lb::SimulationState * simState,
                                                 unsigned int spreadFactor) :
        net::PhasedBroadcastRegular<>(net, simState, spreadFactor), spreadFactor(spreadFactor), callCounter(0)
    {
    }

    BroadcastMockRootNode::~BroadcastMockRootNode()
    {
    }

    void BroadcastMockRootNode::RequestComms()
    {
      // Action taken when a node has to receive from its children in the tree. Calls ReceiveFromChildren mock below.
      ProgressFromChildren(0);
    }

    void BroadcastMockRootNode::PostReceive()
    {
      // Action taken after data has been received from its children
      PostReceiveFromChildren(0);

      // Action taken after data has been sent to its parent
      PostSendToParent(0);
    }

    void BroadcastMockRootNode::EndIteration()
    {
      // Action taken when upwards-travelling data reaches the top node
      TopNodeAction();

      // Action taken by all nodes when downwards-travelling data has been sent to every node.
      Effect();
    }

    
    BroadcastMockLeafNode::BroadcastMockLeafNode(net::Net * net,
                                                 const lb::SimulationState * simState,
                                                 unsigned int spreadFactor) :
        net::PhasedBroadcastRegular<>(net, simState, spreadFactor), /*spreadFactor(spreadFactor),*/ iterationCounter(0)
    {
    }

    BroadcastMockLeafNode::~BroadcastMockLeafNode()
    {
    }

    void BroadcastMockLeafNode::RequestComms()
    {
      if (DownwardPass(iterationCounter))
      {
        // Action taken when a node has to receive from its parent in the tree. Calls ReceiveFromParent mock below.
        ProgressFromParent(0);
      }
      else
      {
        // Action taken when a node has to send to its parent in the tree. Calls SendToParent mock below
        ProgressToParent(0);
      }

      iterationCounter++;
    }

    void BroadcastMockLeafNode::PostReceive()
    {
      // Action taken after data has been sent to its parent
      PostSendToParent(0);
    }

    void BroadcastMockLeafNode::EndIteration()
    {
      // Action taken by all nodes when downwards-travelling data has been sent to every node.
      Effect();
    }

    bool BroadcastMockLeafNode::DownwardPass(unsigned iterationCounter)
    {
      return (iterationCounter % 2 == 0);
    }

  }
}

