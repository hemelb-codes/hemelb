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
    class PhasedBroadcast : public IteratedAction
    {
      public:
        PhasedBroadcast(Net * iNet,
                        const topology::NetworkTopology *iNetTop,
                        const lb::SimulationState * iSimState,
                        unsigned int spreadFactor);

        /**
         * Overridable function to reset the state of the broadcaster.
         */
        virtual void Reset();

        /**
         * Function that requests all the communications from the Net object.
         */
        void RequestComms();

        /**
         * Function to be called after the Receives have completed, where the
         * data is used.
         */
        void PostReceive();

      protected:
        /**
         * Overridable function for when a node has to receive from its children in the tree.
         *
         * Use ReceiveFromChildren to do this.
         */
        virtual void ProgressFromChildren();

        /**
         * Overridable function for when a node has to receive from its parent in the tree.
         *
         * Use ReceiveFromParent to do this.
         */
        virtual void ProgressFromParent();

        /**
         * Overridable function for when a node has to send to its children in the tree.
         *
         * Use SendToChildren to do this.
         */
        virtual void ProgressToChildren();

        /**
         * Overridable function for when a node has to send to its parent in the tree.
         *
         * Use SendToParent to do this.
         */
        virtual void ProgressToParent();

        /**
         * Overridable function, called by a node after data has been received from its children.
         */
        virtual void PostReceiveFromChildren();

        /**
         * Overridable function, called by a node after data has been received from its parent.
         */
        virtual void PostReceiveFromParent();

        /**
         * Action taken when upwards-travelling data reaches the top node.
         */
        virtual void TopNodeAction();

        /**
         * Action taken by all nodes when downwards-travelling data has been sent to every node.
         */
        virtual void Effect();

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
        unsigned int GetRoundTripLength() const;

      private:
        static const int NOPARENT = -1;

        /**
         * Returns the number of the iteration, as an integer between inclusive-0 and
         * exclusive-2 * (the tree depth)
         */
        unsigned long Get0IndexedIterationNumber() const;

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
        const topology::NetworkTopology * mNetTop;
    };
  }
}

#endif /* HEMELB_NET_PHASEDBROADCAST_H */
