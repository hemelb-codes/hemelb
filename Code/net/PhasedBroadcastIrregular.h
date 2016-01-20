
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_NET_PHASEDBROADCASTIRREGULAR_H
#define HEMELB_NET_PHASEDBROADCASTIRREGULAR_H

#include <list>

#include "net/PhasedBroadcast.h"

namespace hemelb
{
  namespace net
  {
    /**
     * PhasedBroadcastIrregular - a class for performing phased broadcasts that are not restricted
     * to starting at regular intervals. A longer description is given in PhasedBroadcast.h.
     */
    template<bool initAction = false, unsigned splay = 1, unsigned ovrlp = 0, bool down = true,
        bool up = true>
    class PhasedBroadcastIrregular : public PhasedBroadcast<initAction, splay, ovrlp, down, up>
    {
        // Typedef for the base class's type
        typedef PhasedBroadcast<initAction, splay, ovrlp, down, up> base;
        typedef std::list<unsigned long> storeType;

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
          performInstantBroadcast = false;
        }

        /**
         * Request the broadcasting to start on this iteration. Returns the number of the iteration
         * on which it will complete.
         *
         * @return
         */
        unsigned long Start()
        {
          unsigned long currentTimeStep = base::mSimState->GetTimeStep();
          unsigned long finishTime = currentTimeStep + base::GetRoundTripLength() - 1;

          // If it's going to take longer than is available for the simulation, assume the
          // implementing class will do its thing instantaneously this iteration.
          if (finishTime > base::mSimState->GetTotalTimeSteps())
          {
            performInstantBroadcast = true;
            return currentTimeStep;
          }
          else
          {
            if (startIterations.empty() || startIterations.back() != currentTimeStep)
            {
              startIterations.push_back(currentTimeStep);
            }

            return finishTime;
          }
        }

        /**
         * Function that requests all the communications from the Net object.
         */
        void RequestComms()
        {
          const unsigned long currentIt = base::mSimState->GetTimeStep();
          const unsigned long firstAscent = base::GetFirstAscending();
          const unsigned long firstDescent = base::GetFirstDescending();

          for (storeType::const_iterator it = startIterations.begin(); it != startIterations.end(); it++)
          {
            unsigned long progress = currentIt - *it;

            // Next, deal with the case of a cycle with an initial pass down the tree.
            if (down)
            {
              if (progress >= firstDescent && progress < firstAscent)
              {
                unsigned long sendOverlap;
                unsigned long receiveOverlap;

                if (base::GetSendChildrenOverlap(progress - firstDescent, &sendOverlap)
                    && !base::GetChildren().empty())
                {
                  ProgressToChildren(*it, sendOverlap);
                }

                if (base::GetReceiveParentOverlap(progress - firstDescent, &receiveOverlap)
                    && base::GetParent() >= 0)
                {
                  ProgressFromParent(*it, receiveOverlap);
                }
              }
            }

            // Now deal with the case of a pass up the tree.
            if (up)
            {
              if (progress >= firstAscent)
              {
                unsigned long sendOverlap;
                unsigned long receiveOverlap;

                if (base::GetSendParentOverlap(progress - firstAscent, &sendOverlap)
                    && base::GetParent() >= 0)
                {
                  ProgressToParent(*it, sendOverlap);
                }

                if (base::GetReceiveChildrenOverlap(progress - firstAscent, &receiveOverlap)
                    && !base::GetChildren().empty())
                {
                  ProgressFromChildren(*it, receiveOverlap);
                }
              }
            }
          }

          unsigned long searchValue = currentIt - (base::GetRoundTripLength() + 1);

          // Use this time to clear out the array. We do this once every iteration so it
          // suffices to get rid of one value.
          for (storeType::iterator it = startIterations.begin(); it != startIterations.end(); it++)
          {
            if (*it == searchValue)
            {
              ClearOut(searchValue);
              startIterations.erase(it);
              break;
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
            if (IsInitialAction())
            {
              InitialAction(base::mSimState->GetTimeStep());
            }
          }
        }

        /**
         * Function to be called after the Receives have completed, where the
         * data is used.
         */
        void PostReceive()
        {
          const unsigned long firstAscent = base::GetFirstAscending();
          const unsigned long traversalLength = base::GetTraverseTime();
          const unsigned long cycleLength = base::GetRoundTripLength();
          const unsigned long currentIt = base::mSimState->GetTimeStep();

          for (storeType::const_iterator it = startIterations.begin(); it != startIterations.end(); it++)
          {
            const unsigned long progress = currentIt - *it;

            // Deal with the case of a cycle with an initial pass down the tree.
            if (down)
            {
              const unsigned long firstDescent = base::GetFirstDescending();

              if (progress >= firstDescent && progress < firstAscent)
              {
                unsigned long receiveOverlap;

                if (base::GetReceiveParentOverlap(progress - firstDescent, &receiveOverlap))
                {
                  PostReceiveFromParent(*it, receiveOverlap);
                }

                // If we're halfway through the programme, all top-down changes have occurred and
                // can be applied on all nodes at once safely.
                if (progress == (traversalLength - 1))
                {
                  Effect(currentIt);
                }
              }
            }

            // Deal with the case of a pass up the tree.
            if (up)
            {
              if (progress >= firstAscent && progress < cycleLength)
              {
                unsigned long receiveOverlap;

                if (base::GetReceiveChildrenOverlap(progress - firstAscent, &receiveOverlap)
                    && !base::GetChildren().empty())
                {
                  PostReceiveFromChildren(*it, receiveOverlap);
                }

                unsigned long sendOverlap;

                if (base::GetSendParentOverlap(progress - firstAscent, &sendOverlap)
                    && base::GetParent() >= 0)
                {
                  PostSendToParent(*it, sendOverlap);
                }
              }
            }

            // If this node is the root of the tree and we've just finished the upwards half, it
            // must act.
            if (progress == (base::GetRoundTripLength() - 1)
                && this->mNet->Rank() == 0)
            {
              TopNodeAction(*it);
            }
          }

          if (performInstantBroadcast)
          {
            InstantBroadcast(currentIt);
            performInstantBroadcast = false;
          }
        }

      protected:
        /**
         * Returns true if we are performing an instant broadcast on the current iteration.
         *
         * @return
         */
        bool IsInstantBroadcast() const
        {
          return performInstantBroadcast;
        }

        /**
         * Returns true if we are performing the initial action for a phased broadcast on
         * the current iteration.
         *
         * @return
         */
        bool IsInitialAction() const
        {
          for (storeType::const_iterator it = startIterations.begin(); it != startIterations.end(); it++)
          {
            if (*it == base::mSimState->GetTimeStep())
            {
              return true;
            }
          }

          return false;
        }

        /**
         * Overridable function for the initial action performed by a node at the beginning of the
         * cycle. Only has an effect if the template paramter initialAction is true.
         */
        virtual void InitialAction(unsigned long startIteration)
        {

        }

        /**
         * Overridable function for when a node has to receive from its children in the tree.
         *
         * Use ReceiveFromChildren to do this. The parameter splayNumber is 0 indexed and less
         * than splay.
         */
        virtual void ProgressFromChildren(unsigned long startIteration, unsigned long splayNumber)
        {

        }

        /**
         * Overridable function for when a node has to receive from its parent in the tree.
         *
         * Use ReceiveFromParent to do this. The parameter splayNumber is 0 indexed and less
         * than splay.
         */
        virtual void ProgressFromParent(unsigned long startIteration, unsigned long splayNumber)
        {

        }

        /**
         * Overridable function for when a node has to send to its children in the tree.
         *
         * Use SendToChildren to do this. The parameter splayNumber is 0 indexed and less
         * than splay.
         */
        virtual void ProgressToChildren(unsigned long startIteration, unsigned long splayNumber)
        {

        }

        /**
         * Overridable function for when a node has to send to its parent in the tree.
         *
         * Use SendToParent to do this. The parameter splayNumber is 0 indexed and less
         * than splay.
         */
        virtual void ProgressToParent(unsigned long startIteration, unsigned long splayNumber)
        {

        }

        /**
         * Overridable function, called by a node after data has been received from its children.
         * The parameter splayNumber is 0 indexed and less than splay.
         */
        virtual void PostReceiveFromChildren(unsigned long startIteration,
                                             unsigned long splayNumber)
        {

        }

        /**
         * Overridable function, called by a node after data has been received from its parent. The
         * parameter splayNumber is 0 indexed and less than splay.
         */
        virtual void PostReceiveFromParent(unsigned long startIteration, unsigned long splayNumber)
        {

        }

        /**
         * Overridable function, called by a node after data has been sent to its parent.
         * @param startIteration The iteration on which this phased send began.
         * @param splayNumber The number of steps we have passed through this phase of the sending.
         */
        virtual void PostSendToParent(unsigned long startIteration, unsigned long splayNumber)
        {

        }

        /**
         * Action taken when upwards-travelling data reaches the top node.
         */
        virtual void TopNodeAction(unsigned long startIteration)
        {

        }

        /**
         * Action taken by all nodes when downwards-travelling data has been sent to every node.
         */
        virtual void Effect(unsigned long startIteration)
        {

        }

        /**
         * For the case where we don't have enough iterations for a phased broadcast to complete.
         */
        virtual void InstantBroadcast(unsigned long startIteration)
        {

        }

        /**
         * For clearing out the results of an iteration after completion.
         */
        virtual void ClearOut(unsigned long startIteration)
        {

        }

      private:
        storeType startIterations;
        bool performInstantBroadcast;
    };
  }
}

#endif /* HEMELB_NET_PHASEDBROADCASTIRREGULAR_H */
