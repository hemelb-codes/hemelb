#ifndef HEMELB_LB_BOUNDARIES_INOUTLETCYCLE_H
#define HEMELB_LB_BOUNDARIES_INOUTLETCYCLE_H

#include "lb/boundaries/iolets/InOutLet.h"

namespace hemelb
{
  namespace lb
  {
    namespace boundaries
    {
      namespace iolets
      {
        /*
         * This class exists to allow templating the update period as well as the bool which determines
         * whether the density values are to be updated locally or centrally through the BCproc. These
         * depend on the update rule of the InOutLet so it was chosen to be templated and it also allows
         * for some IOlets to be updated at regular intervals from the BCproc whilst allowing others
         * with simpler updating rules to update locally and save on communications.
         *
         * tUpdatePeriod - densityCycle is recalculated every tUpdatePeriod time steps. If
         *                 it is 0 then the densityCycle is only initialised at the begining
         *                 and updated only when simulation is reset
         *
         * tComms - if true BCproc updates values and sends them out. If false relevant procs update
         *          locally.
         */
        template<unsigned long tUpdatePeriod, bool tComms>
        class InOutLetCycle : public InOutLet
        {
          public:
            virtual void InitialiseCycle(std::vector<distribn_t> &densityCycle,
                                         const SimulationState *iState);
            virtual void UpdateCycle(std::vector<distribn_t> &densityCycle,
                                     const SimulationState *iState);
            virtual bool DoComms();

          protected:
            InOutLetCycle();
            virtual ~InOutLetCycle();
        };

        template<unsigned long tUpdatePeriod, bool tComms>
        InOutLetCycle<tUpdatePeriod, tComms>::InOutLetCycle() :
            InOutLet()
        {

        }

        template<unsigned long tUpdatePeriod, bool tComms>
        InOutLetCycle<tUpdatePeriod, tComms>::~InOutLetCycle()
        {

        }

        template<unsigned long tUpdatePeriod, bool tComms>
        void InOutLetCycle<tUpdatePeriod, tComms>::InitialiseCycle(std::vector<distribn_t> &densityCycle,
                                                                   const SimulationState *iState)
        {
          // Currently this is unnecessary, but it may be only done here in the future
          // Only called when initialising so doesn't impact performance anyway
          ResetValues();

          if (tUpdatePeriod == 0)
          {
            densityCycle.resize(iState->GetTimeStepsPerCycle());
          }
          else
          {
            densityCycle.resize(tUpdatePeriod);
          }

          CalculateCycle(densityCycle, iState);
        }

        template<unsigned long tUpdatePeriod, bool tComms>
        void InOutLetCycle<tUpdatePeriod, tComms>::UpdateCycle(std::vector<distribn_t> &densityCycle,
                                                               const SimulationState *iState)
        {
          if (tUpdatePeriod == 0)
          {
            return;
          }
          else if (iState->Get0IndexedTimeStep() % tUpdatePeriod == 0)
          {
            CalculateCycle(densityCycle, iState);
          }
        }

        template<unsigned long tUpdatePeriod, bool tComms>
        bool InOutLetCycle<tUpdatePeriod, tComms>::DoComms()
        {
          return tComms;
        }

      }
    }
  }
}

#endif /* HEMELB_LB_BOUNDARIES_INOUTLETCYCLE_H */
