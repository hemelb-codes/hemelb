// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_NET_PHASEDBROADCAST_H
#define HEMELB_NET_PHASEDBROADCAST_H

#include <vector>

#include "debug/Debugger.h"
#include "net/IteratedAction.h"
#include "net/net.h"
#include "lb/SimulationState.h"
#include "net/IOCommunicator.h"

namespace hemelb
{
  namespace net
  {
    /*
     * PROGRAMME
     * ---------
     *
     * For 0-indexed cycle counts and 0-indexed depths, tree depth = N.
     *
     * Iterations 0 to (N-1) are nodes at depth (it) passing down to nodes at depth (it + 1).
     * Action happens on all nodes at the end of iteration (N-1).
     * Iterations N to (2N-1) are nodes at depth (2N - it) passing up to nodes at depth (2N - it - 1).
     */

    /**
     * PhasedBroadcast - a class to control the general process of communication between
     * all nodes over multiple iterations. By using a common interface, we can ensure that
     * communication happens asynchronously at a single point in each iteration - only one
     * communication from any node to any other node will be performed, giving efficient
     * performance.
     *
     * Communication uses a tree structure with a single top-level node, and a constant number of
     * children per node down the tree. All nodes at the same depth communicate on the same
     * iteration with nodes either above or below them (depending on the position through the
     * communication programme defined above). The class supports actions that have to be performed
     * before any communication begins, followed by potentially multiple, overlapping communication
     * stages between related nodes in each consecutive pair of depths in the tree. Communication
     * can be up the tree towards the topmost node, down the tree (beginning at the topmost node)
     * or both (down the tree then up). An action to be performed on all nodes, when communication
     * from top to bottom of the tree has been completed, is also supported.
     *
     * The class is called via the IteratedAction interface. Classes that use this interface should
     * derive from either PhasedBroadcastRegular (for communication that is to happen at regular
     * intervals) or PhasedBroadcastIrregular (for communication that is to happen at irregular
     * intervals). The derived class should override one or more virtual functions in the chosen
     * base class, and should perform communication using the templated
     * [SendTo|ReceiveFrom][Children|Parent] functions defined in this class.
     *
     * This class is made general using template parameters:
     *
     * initialAction = if true, an extra iteration occurs at the start of each broadcast cycle
     * splay = the number of consecutive iterations communication between a pair of nodes needs to
     *   go on for. Useful if the passed data is an array of variable length; one node can spend an
     *   iteration telling the other how many elements will be passed then the next iteration
     *   sending them.
     * ovrlp = the number of iterations where a node both receives and sends. overlap <= splay.
     * down = true if parents communicate to their child nodes in the pattern.
     * up = true if children communicate to their parent nodes in the pattern.
     */
    template<bool initialAction = false, unsigned splay = 1, unsigned ovrlp = 0, bool down = true,
        bool up = true>
    class PhasedBroadcast : public IteratedAction
    {
      public:
        PhasedBroadcast(Net * iNet,
                        const lb::SimulationState * iSimState,
                        unsigned int spreadFactor) :
                          mSimState(iSimState), mMyDepth(0), mTreeDepth(0), mNet(iNet)
        {
          // Calculate the correct values for the depth variables.
          proc_t noSeenToThisDepth = 1;
          proc_t noAtCurrentDepth = 1;

          const MpiCommunicator& netTop = mNet->GetCommunicator();

          while (noSeenToThisDepth < netTop.Size())
          {
            // Go down a level. I.e. increase the depth of the tree, to a new level which has M times
            // as many nodes on it.
            ++mTreeDepth;
            noAtCurrentDepth *= spreadFactor;
            noSeenToThisDepth += noAtCurrentDepth;

            // If this node is at the current depth, it must have a rank lower than or equal to the highest
            // rank at the current depth but greater than the highest rank at the previous depth.
            if (noSeenToThisDepth > netTop.Rank() && ( (noSeenToThisDepth
                - noAtCurrentDepth) <= netTop.Rank()))
            {
              mMyDepth = mTreeDepth;
            }
          }

          // In a M-tree, with a root of 0, each node N's parent has rank floor((N-1) / M)
          if (netTop.Rank() == 0)
          {
            mParent = NOPARENT;
          }
          else
          {
            mParent = (netTop.Rank() - 1) / spreadFactor;
          }

          // The children of a node N in a M-tree with root 0 are those in the range (M*N)+1,...,(M*N) + M
          for (unsigned int child = (spreadFactor * netTop.Rank()) + 1; child
              <= spreadFactor * (1 + netTop.Rank()); ++child)
          {
            if (child < (unsigned int) netTop.Size())
            {
              mChildren.push_back(child);
            }
          }
        }

      protected:

        const lb::SimulationState * mSimState;

        /**
         * Helper function for sending data to child nodes.
         */
        template<class T>
        void SendToChildren(T* data, int count)
        {
          for (std::vector<int>::iterator it = mChildren.begin(); it != mChildren.end(); ++it)
          {
            mNet->RequestSend<T> (data, count, *it);
          }
        }

        /**
         * Helper function for receiving data from parent nodes.
         */
        template<class T>
        void ReceiveFromParent(T* data, int count)
        {
          if (mParent != NOPARENT)
          {
            mNet ->RequestReceive<T> (data, count, mParent);
          }
        }

        /**
         * Receives data from each child. This is a set length per child, and each child's
         * data is inserted contiguously into the provided array.
         *
         * The user must handle the case where the number of children is smaller than the
         * spread factor and the array pointed by dataStart is not completely filled in.
         *
         * @param dataStart Pointer to the start of the array.
         * @param countPerChild Number of elements to receive per child.
         */
        template<class T>
        void ReceiveFromChildren(T* dataStart, int countPerChild)
        {
          T* data = dataStart;
          for (std::vector<int>::iterator it = mChildren.begin(); it != mChildren.end(); ++it)
          {
            mNet->RequestReceive<T> (data, countPerChild, *it);
            data += countPerChild;
          }
        }

        /**
         * Receives data from each child. This is a variable length per child, and each child's
         * data is inserted into the appropriate array.
         *
         * The user must handle the case where the number of children is smaller than the
         * spread factor and the array pointed by dataStart is not completely filled in.
         *
         * @param dataStart Pointer to the start of the array.
         * @param countPerChild Number of elements to receive per child.
         */
        template<class T>
        void ReceiveFromChildren(T** dataStart, unsigned int* countPerChild)
        {
          unsigned int childNum = 0;
          for (std::vector<int>::iterator it = mChildren.begin(); it != mChildren.end(); ++it)
          {
            mNet->RequestReceive<T> (dataStart[childNum], countPerChild[childNum], *it);
            ++childNum;
          }
        }

        /**
         * Helper function for sending data to parent nodes.
         */
        template<class T>
        void SendToParent(T* data, int count)
        {
          if (mParent != NOPARENT)
          {
            mNet->RequestSend<T> (data, count, mParent);
          }
        }

        /**
         * Returns the total number of iterations spent doing a complete traversal -- all the way
         * down the tree and back up again.
         *
         * @return
         */
        unsigned long GetRoundTripLength() const
        {
          unsigned long delayTime = initialAction
            ? 1
            : 0;

          unsigned long multiplier = (down
            ? 1
            : 0) + (up
            ? 1
            : 0);

          return delayTime + multiplier * GetTraverseTime();
        }

        /**
         * Get the index of the iteration when we first start descending the tree.
         *
         * @return
         */
        unsigned long GetFirstDescending() const
        {
          return (initialAction
            ? 1
            : 0);
        }

        /**
         * Get the index of the iteration when we first start ascending the tree.
         *
         * @return
         */
        unsigned long GetFirstAscending() const
        {
          return GetFirstDescending() + (down
            ? GetTraverseTime()
            : 0);
        }

        /**
         * Get the 0-indexed depth of this rank.
         *
         * @return
         */
        unsigned long GetMyDepth() const
        {
          return mMyDepth;
        }

        /**
         * Get the 0-indexed depth of the whole tree.
         *
         * @return
         */
        unsigned long GetTreeDepth() const
        {
          return mTreeDepth;
        }

        /**
         * Get the number of iterations required for a half-cycle (messages going either all the
         * way up the tree or all the way down)
         *
         * @return
         */
        unsigned long GetTraverseTime() const
        {
          return mTreeDepth * (splay - ovrlp) + ovrlp;
        }

        /**
         * Gets the overlap value for a send to a parent node given the number of iterations
         * through the 'upwards' phase. Returns true if this node should send to is parent.
         */
        bool GetSendParentOverlap(unsigned long subCycleNumber, unsigned long* sendOverlap)
        {
          // The first cycle we're sending to parents.
          unsigned long firstSendCycle = (mTreeDepth - mMyDepth) * (splay - ovrlp);

          // If we're either not far enough or too far through the cycle, don't send.
          if (subCycleNumber < firstSendCycle || subCycleNumber >= (firstSendCycle + splay))
          {
            return false;
          }
          // Otherwise calculate the overlap value and return true.
          else
          {
            *sendOverlap = subCycleNumber - firstSendCycle;
            return true;
          }
        }

        /**
         * Gets the overlap value for a send to child nodes given the number of iterations
         * through the 'downwards' phase. Returns true if this node should send to is children.
         */
        bool GetSendChildrenOverlap(unsigned long subCycleNumber, unsigned long* sendOverlap)
        {
          // The first cycle we're sending to children.
          unsigned long firstSendCycle = mMyDepth * (splay - ovrlp);

          // If we're either not far enough or too far through the cycle, don't send.
          if (subCycleNumber < firstSendCycle || subCycleNumber >= (firstSendCycle + splay))
          {
            return false;
          }
          // Otherwise calculate the overlap value and return true.
          else
          {
            *sendOverlap = subCycleNumber - firstSendCycle;
            return true;
          }
        }

        /**
         * Gets the overlap value for a receive from a parent node given the number of iterations
         * through the 'downwards' phase. Returns true if this node should receive from its parent.
         */
        bool GetReceiveParentOverlap(unsigned long subCycleNumber, unsigned long* receiveOverlap)
        {
          // The first cycle we're receiving from parents.
          unsigned long firstReceiveCycle = (mMyDepth - 1) * (splay - ovrlp);

          // If we're either not far enough or too far through the cycle, don't receive.
          if (subCycleNumber < firstReceiveCycle || subCycleNumber >= (firstReceiveCycle + splay))
          {
            return false;
          }
          // Otherwise calculate the overlap value and return true.
          else
          {
            *receiveOverlap = subCycleNumber - firstReceiveCycle;
            return true;
          }
        }

        /**
         * Gets the overlap value for a receive from child node given the number of iterations
         * through the 'upwards' phase. Returns true if this node should receive from its children.
         */
        bool GetReceiveChildrenOverlap(unsigned long subCycleNumber, unsigned long* receiveOverlap)
        {
          // The first cycle we're receiving from parents.
          unsigned long firstReceiveCycle = (mTreeDepth - (mMyDepth + 1)) * (splay - ovrlp);

          // If we're either not far enough or too far through the cycle, don't receive.
          if (subCycleNumber < firstReceiveCycle || subCycleNumber >= (firstReceiveCycle + splay))
          {
            return false;
          }
          // Otherwise calculate the overlap value and return true.
          else
          {
            *receiveOverlap = subCycleNumber - firstReceiveCycle;
            return true;
          }
        }

        /**
         * Returns the number of the parent node.
         * @return
         */
        int GetParent() const
        {
          return mParent;
        }

        const std::vector<int>& GetChildren() const
        {
          return mChildren;
        }

      protected:
        static const int NOPARENT = -1;

      private:
        /**
         * Note that depths are 0-indexed. I.e. a tree with a single node has
         * a depth of zero, and the node itself is considered to be at depth 0.
         */
        unsigned int mMyDepth;
        unsigned int mTreeDepth;

        /**
         * This node's parent rank.
         */
        int mParent;
        /**
         * This node's child ranks.
         */
        std::vector<int> mChildren;
      protected:
        Net * mNet;
    };
  }
}

#endif /* HEMELB_NET_PHASEDBROADCAST_H */
