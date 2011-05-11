#ifndef HEMELB_NET_PHASEDBROADCAST_H
#define HEMELB_NET_PHASEDBROADCAST_H

#include <vector>
#include "net/IteratedAction.h"
#include "net/net.h"
#include "lb/SimulationState.h"
#include "topology/NetworkTopology.h"

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
     * This class is made general using template parameters.
     *
     * initialAction = if true, an extra iteration occurs at the start of each broadcast cycle
     * splay = the number of consecutive iterations communication between a pair of nodes needs to
     *   go on for. Useful if the passed data is an array of variable length; one node can spend an
     *   iteration telling the other how many elements will be passed then the next iteration
     *   sending them.
     * overlap = the number of iterations where a node both receives and sends. overlap <= splay.
     * goDown = true if parents communicate to their child nodes in the pattern.
     * goUp = true if children communicate to their parent nodes in the pattern.
     */
    template<bool initialAction = false, int splay = 1, int overlap = 0, bool goDown = true,
        bool goUp = true>
    class PhasedBroadcast : public IteratedAction
    {
      public:
        PhasedBroadcast(Net * iNet,
                        const lb::SimulationState * iSimState,
                        unsigned int spreadFactor)
        {
          mNet = iNet;
          mSimState = iSimState;

          // Initialise member variables.
          mTreeDepth = 0;
          mMyDepth = 0;
          mSpreadFactor = spreadFactor;

          // Calculate the correct values for the depth variables.
          proc_t noSeenToThisDepth = 1;
          proc_t noAtCurrentDepth = 1;

          hemelb::topology::NetworkTopology* netTop = hemelb::topology::NetworkTopology::Instance();

          while (noSeenToThisDepth < netTop->GetProcessorCount())
          {
            // Go down a level. I.e. increase the depth of the tree, to a new level which has M times
            // as many nodes on it.
            ++mTreeDepth;
            noAtCurrentDepth *= spreadFactor;
            noSeenToThisDepth += noAtCurrentDepth;

            // If this node is at the current depth, it must have a rank lower than or equal to the highest
            // rank at the current depth but greater than the highest rank at the previous depth.
            if (noSeenToThisDepth > netTop->GetLocalRank() && ( (noSeenToThisDepth
                - noAtCurrentDepth) <= netTop->GetLocalRank()))
            {
              mMyDepth = mTreeDepth;
            }
          }

          // In a M-tree, with a root of 0, each node N's parent has rank floor((N-1) / M)
          if (netTop->GetLocalRank() == 0)
          {
            mParent = NOPARENT;
          }
          else
          {
            mParent = (netTop->GetLocalRank() - 1) / spreadFactor;
          }

          // The children of a node N in a M-tree with root 0 are those in the range (M*N)+1,...,(M*N) + M
          for (unsigned int child = (spreadFactor * netTop->GetLocalRank()) + 1; child
              <= spreadFactor * (1 + netTop->GetLocalRank()); ++child)
          {
            if (child < (unsigned int) netTop->GetProcessorCount())
            {
              mChildren.push_back(child);
            }
          }
        }

        /**
         * Overridable function to reset the state of the broadcaster.
         */
        virtual void Reset()
        {

        }

        /**
         * Function that requests all the communications from the Net object.
         */
        void RequestComms()
        {
          unsigned int iCycleNumber = Get0IndexedIterationNumber();

          // Passing down the tree.
          if (iCycleNumber < mTreeDepth)
          {
            if (mMyDepth == iCycleNumber)
            {
              ProgressToChildren();
            }
            else if (mMyDepth == (iCycleNumber + 1))
            {
              ProgressFromParent();
            }
          }
          // Passing up the tree.
          else
          {
            if (mMyDepth == (2 * mTreeDepth - iCycleNumber))
            {
              ProgressToParent();
            }
            else if (mMyDepth == (2 * mTreeDepth - iCycleNumber - 1))
            {
              ProgressFromChildren();
            }
          }
        }

        /**
         * Function to be called after the Receives have completed, where the
         * data is used.
         */
        void PostReceive()
        {
          unsigned int iCycleNumber = Get0IndexedIterationNumber();
          // Passing down the tree.
          if (mMyDepth == (iCycleNumber + 1))
          {
            PostReceiveFromParent();
          }
          // Passing up the tree to this node.
          else if (mMyDepth == (2 * mTreeDepth - iCycleNumber - 1))
          {
            PostReceiveFromChildren();

            // If this node is the root of the tree, it must act.
            if (topology::NetworkTopology::Instance()->GetLocalRank() == 0)
            {
              TopNodeAction();
            }
          }

          // If we're halfway through the programme, all top-down changes have occurred and
          // can be applied on all nodes at once safely.
          if (iCycleNumber == (mTreeDepth - 1))
          {
            Effect();
          }
        }

      protected:

        /**
         * Overridable function for the initial action performed by a node at the beginning of the
         * cycle. Only has an effect if the template paramter initialAction is true.
         */
        virtual void InitialAction()
        {

        }

        /**
         * Overridable function for when a node has to receive from its children in the tree.
         *
         * Use ReceiveFromChildren to do this.
         */
        virtual void ProgressFromChildren()
        {

        }

        /**
         * Overridable function for when a node has to receive from its parent in the tree.
         *
         * Use ReceiveFromParent to do this.
         */
        virtual void ProgressFromParent()
        {

        }

        /**
         * Overridable function for when a node has to send to its children in the tree.
         *
         * Use SendToChildren to do this.
         */
        virtual void ProgressToChildren()
        {

        }

        /**
         * Overridable function for when a node has to send to its parent in the tree.
         *
         * Use SendToParent to do this.
         */
        virtual void ProgressToParent()
        {

        }

        /**
         * Overridable function, called by a node after data has been received from its children.
         */
        virtual void PostReceiveFromChildren()
        {

        }

        /**
         * Overridable function, called by a node after data has been received from its parent.
         */
        virtual void PostReceiveFromParent()
        {

        }

        /**
         * Action taken when upwards-travelling data reaches the top node.
         */
        virtual void TopNodeAction()
        {

        }

        /**
         * Action taken by all nodes when downwards-travelling data has been sent to every node.
         */
        virtual void Effect()
        {

        }

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
        unsigned int GetRoundTripLength() const
        {
          return 2 * mTreeDepth;
        }

      private:
        static const int NOPARENT = -1;

        /**
         * Returns the number of the iteration, as an integer between inclusive-0 and
         * exclusive-2 * (the tree depth)
         */
        unsigned long Get0IndexedIterationNumber() const
        {
          if (mTreeDepth > 0)
          {
            return (mSimState->TimeStep - 1) % GetRoundTripLength();
          }
          else
          {
            return 0;
          }
        }

        /**
         * Note that depths are 0-indexed. I.e. a tree with a single node has
         * a depth of zero, and the node itself is considered to be at depth 0.
         */
        unsigned int mMyDepth;
        unsigned int mTreeDepth;

        /**
         * The maximum number of children and node has.
         */
        unsigned int mSpreadFactor;

        /**
         * This node's parent rank.
         */
        int mParent;
        /**
         * This node's child ranks.
         */
        std::vector<int> mChildren;

        const lb::SimulationState * mSimState;
        Net * mNet;
    };
  }
}

#endif /* HEMELB_NET_PHASEDBROADCAST_H */
