#ifndef HEMELB_NET_PHASEDBROADCAST_H
#define HEMELB_NET_PHASEDBROADCAST_H

#include <vector>

#include "debug/Debugger.h"
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
    template<bool initialAction = false, unsigned int splay = 1, unsigned int overlap = 0,
        bool goDown = true, bool goUp = true>
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
          const unsigned int iCycleNumber = Get0IndexedIterationNumber();
          const unsigned int firstAscent = GetFirstAscending();
          const unsigned int firstDescent = GetFirstDescending();
          const unsigned int traversalLength = GetTraverseTime();

          // Nothing to do for initial action case.

          // Next, deal with the case of a cycle with an initial pass-down the tree.
          if (goDown)
          {
            if (iCycleNumber >= firstDescent && iCycleNumber < firstAscent)
            {
              unsigned int sendOverlap;
              unsigned int receiveOverlap;

              if (GetSendChildrenOverlap(iCycleNumber - firstDescent, &sendOverlap))
              {
                ProgressToChildren(sendOverlap);
              }

              if (GetReceiveParentOverlap(iCycleNumber - firstDescent, &receiveOverlap))
              {
                ProgressFromParent(receiveOverlap);
              }
            }
          }

          if (goUp)
          {
            if (iCycleNumber >= firstAscent)
            {
              unsigned int sendOverlap;
              unsigned int receiveOverlap;

              if (GetSendParentOverlap(iCycleNumber - firstAscent, &sendOverlap))
              {
                ProgressToParent(sendOverlap);
              }

              if (GetReceiveChildrenOverlap(iCycleNumber - firstAscent, &receiveOverlap))
              {
                ProgressFromChildren(receiveOverlap);
              }
            }
          }
        }

        void PreReceive()
        {
          // The only thing to do while waiting is the initial action.
          if (initialAction)
          {
            if (Get0IndexedIterationNumber() == 0)
            {
              InitialAction();
            }
          }
        }

        /**
         * Function to be called after the Receives have completed, where the
         * data is used.
         */
        void PostReceive()
        {
          const unsigned int iCycleNumber = Get0IndexedIterationNumber();
          const unsigned int firstAscent = GetFirstAscending();
          const unsigned int traversalLength = GetTraverseTime();
          const unsigned int cycleLength = GetRoundTripLength();

          // Deal with the case of a cycle with an initial pass-down the tree.
          if (goDown)
          {
            const unsigned int firstDescent = GetFirstDescending();

            if (iCycleNumber >= firstDescent && iCycleNumber < firstAscent)
            {
              unsigned int receiveOverlap;

              if (GetReceiveParentOverlap(iCycleNumber - firstDescent, &receiveOverlap))
              {
                PostReceiveFromParent(receiveOverlap);
              }

              // If we're halfway through the programme, all top-down changes have occurred and
              // can be applied on all nodes at once safely.
              if ( (iCycleNumber - firstDescent) == (traversalLength - 1))
              {
                Effect();
              }
            }
          }

          if (goUp)
          {
            if (iCycleNumber >= firstAscent)
            {
              unsigned int receiveOverlap;

              if (GetReceiveChildrenOverlap(iCycleNumber - firstAscent, &receiveOverlap))
              {
                PostReceiveFromChildren(receiveOverlap);
              }
            }
          }

          // If this node is the root of the tree and we've just finished the upwards half, it
          // must act.
          if (iCycleNumber == (GetRoundTripLength() - 1)
              && topology::NetworkTopology::Instance()->GetLocalRank() == 0)
          {
            TopNodeAction();
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
         * Use ReceiveFromChildren to do this. The parameter splayNumber is 0 indexed and less
         * than splay.
         */
        virtual void ProgressFromChildren(unsigned int splayNumber)
        {

        }

        /**
         * Overridable function for when a node has to receive from its parent in the tree.
         *
         * Use ReceiveFromParent to do this. The parameter splayNumber is 0 indexed and less
         * than splay.
         */
        virtual void ProgressFromParent(unsigned int splayNumber)
        {

        }

        /**
         * Overridable function for when a node has to send to its children in the tree.
         *
         * Use SendToChildren to do this. The parameter splayNumber is 0 indexed and less
         * than splay.
         */
        virtual void ProgressToChildren(unsigned int splayNumber)
        {

        }

        /**
         * Overridable function for when a node has to send to its parent in the tree.
         *
         * Use SendToParent to do this. The parameter splayNumber is 0 indexed and less
         * than splay.
         */
        virtual void ProgressToParent(unsigned int splayNumber)
        {

        }

        /**
         * Overridable function, called by a node after data has been received from its children.
         * The parameter splayNumber is 0 indexed and less than splay.
         */
        virtual void PostReceiveFromChildren(unsigned int splayNumber)
        {

        }

        /**
         * Overridable function, called by a node after data has been received from its parent. The
         * parameter splayNumber is 0 indexed and less than splay.
         */
        virtual void PostReceiveFromParent(unsigned int splayNumber)
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
          unsigned int delayTime = initialAction
            ? 1
            : 0;

          unsigned int multiplier = (goDown
            ? 1
            : 0) + (goUp
            ? 1
            : 0);

          return delayTime + multiplier * GetTraverseTime();
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
            unsigned long stepsPassed = (mSimState->CycleId - 1) * mSimState->TimeStepsPerCycle
                + mSimState->TimeStep - 1;

            return stepsPassed % GetRoundTripLength();
          }
          else
          {
            return 0;
          }
        }

        /**
         * Get the number of iterations required for a half-cycle (messages going either all the
         * way up the tree or all the way down)
         *
         * @return
         */
        unsigned int GetTraverseTime() const
        {
          return mTreeDepth * (splay - overlap) + overlap;
        }

        /*
         * Gets the overlap value for a send to a parent node given the number of iterations
         * through the 'upwards' phase. Returns true if this node should send to is parent.
         */
        bool GetSendParentOverlap(unsigned int subCycleNumber, unsigned int* sendOverlap)
        {
          // The first cycle we're sending to parents.
          unsigned int firstSendCycle = (mTreeDepth - mMyDepth) * (splay - overlap);

          // If we're either not far enough or too far through the cycle, don't send.
          if (subCycleNumber < firstSendCycle || subCycleNumber > (firstSendCycle + overlap))
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

        /*
         * Gets the overlap value for a send to child nodes given the number of iterations
         * through the 'downwards' phase. Returns true if this node should send to is children.
         */
        bool GetSendChildrenOverlap(unsigned int subCycleNumber, unsigned int* sendOverlap)
        {
          // The first cycle we're sending to children.
          unsigned int firstSendCycle = mMyDepth * (splay - overlap);

          // If we're either not far enough or too far through the cycle, don't send.
          if (subCycleNumber < firstSendCycle || subCycleNumber > (firstSendCycle + overlap))
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

        /*
         * Gets the overlap value for a receive from a parent node given the number of iterations
         * through the 'downwards' phase. Returns true if this node should receive from its parent.
         */
        bool GetReceiveParentOverlap(unsigned int subCycleNumber, unsigned int* receiveOverlap)
        {
          // The first cycle we're receiving from parents.
          unsigned int firstReceiveCycle = (mMyDepth - 1) * (splay - overlap);

          // If we're either not far enough or too far through the cycle, don't receive.
          if (subCycleNumber < firstReceiveCycle || subCycleNumber > (firstReceiveCycle + overlap))
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

        /*
         * Gets the overlap value for a receive from child node given the number of iterations
         * through the 'upwards' phase. Returns true if this node should receive from its children.
         */
        bool GetReceiveChildrenOverlap(unsigned int subCycleNumber, unsigned int* receiveOverlap)
        {
          // The first cycle we're receiving from parents.
          unsigned int firstReceiveCycle = (mTreeDepth - (mMyDepth + 1)) * (splay - overlap);

          // If we're either not far enough or too far through the cycle, don't receive.
          if (subCycleNumber < firstReceiveCycle || subCycleNumber > (firstReceiveCycle + overlap))
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

        unsigned int GetFirstDescending() const
        {
          return (initialAction
            ? 1
            : 0);
        }

        unsigned int GetFirstAscending() const
        {
          return GetFirstDescending() + (goDown
            ? GetTraverseTime()
            : 0);
        }

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

        const lb::SimulationState * mSimState;
        Net * mNet;
    };
  }
}

#endif /* HEMELB_NET_PHASEDBROADCAST_H */
