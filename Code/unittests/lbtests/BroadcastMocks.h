
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_UNITTESTS_LBTESTS_BROADCASTMOCKS_H
#define HEMELB_UNITTESTS_LBTESTS_BROADCASTMOCKS_H

#include "net/PhasedBroadcastRegular.h"
#include "lb/SimulationState.h"
#include "net/net.h"

/**
 * This file defines two PhasedBroadcastRegular mocks. One that behaves like the root
 * node of the broadcast tree and another for a leaf node.
 *
 * @todo Split file into .h and .cc ?
 */

namespace hemelb
{
  namespace net
  {
    /**
     * In this mock, we pretend that the current process is the root node of a phased
     * broadcast and that a pair of values is going up the tree. The mock simulates
     * three cycles of communications (i.e. three messages travelling up the tree and
     * reaching the root node):
     *  1) All children report (14,15).
     *  2) One child reports (1,100) and the rest (14,15).
     *  3) All back to (14,15).
     */
    class BroadcastMockRootNode : public net::PhasedBroadcastRegular<>
    {

      public:
        BroadcastMockRootNode(net::Net * net, const lb::SimulationState * simState, unsigned int spreadFactor);

        virtual ~BroadcastMockRootNode();

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

    template<class T>
    void BroadcastMockRootNode::ReceiveFromChildren(T* dataStart, int countPerChild)
    {
      assert(countPerChild == 3);

      switch (callCounter++)
      {
        case 0:
          for (unsigned childIndex = 0; childIndex < spreadFactor; childIndex++)
          {
            dataStart[childIndex * countPerChild] = 14.0;
            dataStart[childIndex * countPerChild + 1] = 15.0;
            dataStart[childIndex * countPerChild + 2] = 0.01;
          }
          break;

        case 1:
          for (unsigned childIndex = 0; childIndex < spreadFactor - 1; childIndex++)
          {
            dataStart[childIndex * countPerChild] = 14.0;
            dataStart[childIndex * countPerChild + 1] = 15.0;
            dataStart[childIndex * countPerChild + 2] = 1.0;
          }
          dataStart[ (spreadFactor - 1) * countPerChild] = 1.0;
          dataStart[ (spreadFactor - 1) * countPerChild + 1] = 100.0;
          dataStart[ (spreadFactor - 1) * countPerChild + 2] = 10.0;
          break;

        case 2:
          for (unsigned childIndex = 0; childIndex < spreadFactor; childIndex++)
          {
            dataStart[childIndex * countPerChild] = 14.0;
            dataStart[childIndex * countPerChild + 1] = 15.0;
            dataStart[childIndex * countPerChild + 2] = 0.01;
          }
          break;

        default:
          // Sanity check. Control should never reach this branch
          assert(false);
      }
    }

    /**
     * In this mock, we pretend that the current process is a leaf node of a phased
     * broadcast and that a pair of values is going down and up the tree. The mock
     * simulates three complete executions of the phased broadcast algorithm. For all
     * three executions, the current leaf node reports (12, 21.45) to its parent (upwards
     * tree pass). For the downwards pass, the current leaf receives:
     *  1) unitialised data.
     *  2) the min and max values reported by the current leaf node itself (12, 21.45).
     *     No other hypothetical leaf node has reported any smaller min or larger max densities.
     *  3) (1,100). Another hypothetical leaf node reported these values during the previous
     *     upwards pass.
     */
    class BroadcastMockLeafNode : public net::PhasedBroadcastRegular<>
    {

      public:
        BroadcastMockLeafNode(net::Net * net, const lb::SimulationState * simState, unsigned int spreadFactor);

        virtual ~BroadcastMockLeafNode();

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
        void ReceiveFromParent(T* dataStart, int countPerChild);

        /**
         * Helper function for sending data to parent nodes.
         */
        template<class T>
        void SendToParent(T* data, int count);

      private:
        /** Number of children nodes in the tree */
        // unsigned spreadFactor;

        /** Number of iterations of the phased broadcast algorithm */
        unsigned iterationCounter;

        /** Minimum and maximum values sent up by the leaf node */
        distribn_t minSentUp, maxSentUp, maxVelSentUp;

        /**
         * Helper method to find out whether we are in a downward pass
         *
         * @param iterationCounter iteration number in the phased broadcast algorithm
         * @return whether we are in a downward pass
         */
        bool DownwardPass(unsigned iterationCounter);
    };

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

    template<class T>
    void BroadcastMockLeafNode::ReceiveFromParent(T* dataStart, int countPerChild)
    {
      assert(countPerChild == 3);

      switch (iterationCounter)
      {
        case 0:
          // The first pass down contains rubbish
          dataStart[0] = DBL_MAX;
          dataStart[1] = -DBL_MAX;
          dataStart[1] = 0;
          break;

        case 2:
          dataStart[0] = minSentUp;
          dataStart[1] = maxSentUp;
          dataStart[2] = maxVelSentUp;
          break;

        case 4:
          dataStart[0] = 1.0;
          dataStart[1] = 100.0;
          dataStart[2] = 10.0;
          break;

        default:
          // Sanity check. Control should never reach this branch
          assert(false);
      }

    }

    template<class T>
    void BroadcastMockLeafNode::SendToParent(T* data, int count)
    {
      assert(count == 3);
      minSentUp = data[0];
      maxSentUp = data[1];
      maxVelSentUp = data[2];
    }

    bool BroadcastMockLeafNode::DownwardPass(unsigned iterationCounter)
    {
      return (iterationCounter % 2 == 0);
    }

  }
}

#endif /* HEMELB_UNITTESTS_LBTESTS_BROADCASTMOCKS_H */
