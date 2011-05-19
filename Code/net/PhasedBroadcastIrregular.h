#ifndef HEMELB_NET_PHASEDBROADCASTIRREGULAR_H
#define HEMELB_NET_PHASEDBROADCASTIRREGULAR_H

#include <queue>

#include "net/PhasedBroadcast.h"

namespace hemelb
{
  namespace net
  {
    /**
     * PhasedBroadcastIrregular - a class for performing phased broadcasts that are not restricted
     * to starting at regular intervals.
     */
    template<bool initAction = false, unsigned splay = 1, unsigned ovrlp = 0, bool down = true,
        bool up = true>
    class PhasedBroadcastIrregular : public PhasedBroadcast<initAction, splay, ovrlp, down, up>
    {
        // Typedef for the base class's type
        typedef PhasedBroadcast<initAction, splay, ovrlp, down, up> base;

      public:
        /**
         * Constructor that initialises the base class.
         *
         * @param iNet
         * @param iSimState
         * @param spreadFactor
         * @return
         */
        PhasedBroadcastIrregular(Net * iNet,
                                 const lb::SimulationState * iSimState,
                                 unsigned int spreadFactor) :
          base(iNet, iSimState, spreadFactor)
        {

        }

        /**
         * Request the broadcasting to start on this iteration. Returns the number of the iteration
         * on which it will complete.
         *
         * @return
         */
        unsigned int Start()
        {
          unsigned int currentTimeStep = base::mSimState->GetTimeStepsPassed();
          if (!startIterations.empty() && startIterations.back() == currentTimeStep)
          {
            startIterations.push_back(currentTimeStep);
          }
          return currentTimeStep + base::GetRoundTripLength() - 1;
        }

        virtual void Reset()
        {
          startIterations.clear();

          base::Reset();
        }

        /**
         * Function that requests all the communications from the Net object.
         */
        void RequestComms()
        {
          const unsigned int iCycleNumber = Get0IndexedIterationNumber();
          const unsigned int firstAscent = base::GetFirstAscending();
          const unsigned int firstDescent = base::GetFirstDescending();
          const unsigned int traversalLength = base::GetTraverseTime();

          // Nothing to do for initial action case.

          // Next, deal with the case of a cycle with an initial pass down the tree.
          if (down)
          {
            if (iCycleNumber >= firstDescent && iCycleNumber < firstAscent)
            {
              unsigned int sendOverlap;
              unsigned int receiveOverlap;

              if (base::GetSendChildrenOverlap(iCycleNumber - firstDescent, &sendOverlap))
              {
                ProgressToChildren(sendOverlap);
              }

              if (base::GetReceiveParentOverlap(iCycleNumber - firstDescent, &receiveOverlap))
              {
                ProgressFromParent(receiveOverlap);
              }
            }
          }

          // Now deal with the case of a pass up the tree.
          if (up)
          {
            if (iCycleNumber >= firstAscent)
            {
              unsigned int sendOverlap;
              unsigned int receiveOverlap;

              if (base::GetSendParentOverlap(iCycleNumber - firstAscent, &sendOverlap))
              {
                ProgressToParent(sendOverlap);
              }

              if (base::GetReceiveChildrenOverlap(iCycleNumber - firstAscent, &receiveOverlap))
              {
                ProgressFromChildren(receiveOverlap);
              }
            }
          }
        }

        /**
         * Function that acts while waiting for data to be received. This performs the initial
         * action on the first iteration only.
         */
        void PreReceive()
        {
          // The only thing to do while waiting is the initial action.
          if (initAction)
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
          const unsigned int firstAscent =
              PhasedBroadcast<initAction, splay, ovrlp, down, up>::GetFirstAscending();
          const unsigned int traversalLength =
              PhasedBroadcast<initAction, splay, ovrlp, down, up>::GetTraverseTime();
          const unsigned int cycleLength =
              PhasedBroadcast<initAction, splay, ovrlp, down, up>::GetRoundTripLength();

          // Deal with the case of a cycle with an initial pass down the tree.
          if (down)
          {
            const unsigned int firstDescent =
                PhasedBroadcast<initAction, splay, ovrlp, down, up>::GetFirstDescending();

            if (iCycleNumber >= firstDescent && iCycleNumber < firstAscent)
            {
              unsigned int receiveOverlap;

              if (base::GetReceiveParentOverlap(iCycleNumber - firstDescent, &receiveOverlap))
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

          // Deal with the case of a pass up the tree.
          if (up)
          {
            if (iCycleNumber >= firstAscent)
            {
              unsigned int receiveOverlap;

              if (base::GetReceiveChildrenOverlap(iCycleNumber - firstAscent, &receiveOverlap))
              {
                PostReceiveFromChildren(receiveOverlap);
              }
            }
          }

          // If this node is the root of the tree and we've just finished the upwards half, it
          // must act.
          if (iCycleNumber == (base::GetRoundTripLength() - 1)
              && topology::NetworkTopology::Instance()->GetLocalRank() == 0)
          {
            TopNodeAction();
          }
        }

        /**
         * Returns the number of the iteration, as an integer between inclusive-0 and
         * exclusive-2 * (the tree depth)
         */
        unsigned long Get0IndexedIterationNumber() const
        {
          if (base::GetTreeDepth() > 0)
          {
            unsigned long stepsPassed = (base::mSimState->CycleId - 1)
                * base::mSimState->TimeStepsPerCycle + base::mSimState->TimeStep - 1;

            return stepsPassed % base::GetRoundTripLength();
          }
          else
          {
            return 0;
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

      private:
        queue<unsigned int> startIterations;
    };
  }
}

#endif /* HEMELB_NET_PHASEDBROADCASTIRREGULAR_H */
