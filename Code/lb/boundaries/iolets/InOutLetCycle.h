#ifndef HEMELB_LB_BOUNDARIES_IOLETS_INOUTLETCYCLE_H
#define HEMELB_LB_BOUNDARIES_IOLETS_INOUTLETCYCLE_H

#include "lb/boundaries/iolets/InOutLet.h"

namespace hemelb
{
  namespace lb
  {
    namespace boundaries
    {
      namespace iolets
      {

        /***
         * InOutLet with varying update period and communications strategy.
         * This class exists to allow templating the update period as well as the bool which determines
         * whether the density values are to be updated locally or centrally through the BCproc. These
         * depend on the update rule of the InOutLet so it was chosen to be templated and it also allows
         * for some IOlets to be updated at regular intervals from the BCproc whilst allowing others
         * with simpler updating rules to update locally and save on communications.
         *
         * @tparam updatePeriod densityCycle is recalculated every updatePeriod time steps. If
         *                 it is 0 then the densityCycle is only initialised at the begining
         *                 and updated only when simulation is reset
         * @tparam comms if true BCproc updates values and sends them out. If false relevant procs update
         *          locally.
         */
        template<unsigned long updatePeriod, bool comms>
        class InOutLetCycle : public InOutLet
        {
          public:

            virtual void UpdateCycle(std::vector<distribn_t> &densityCycle, const SimulationState *state);
            virtual bool DoComms();
            unsigned int GetUpdatePeriod()
            {
              return updatePeriod;
            }
          protected:
            InOutLetCycle();
            virtual ~InOutLetCycle();
        };

        template<unsigned long updatePeriod, bool comms>
        InOutLetCycle<updatePeriod, comms>::InOutLetCycle() :
            InOutLet()
        {

        }

        template<unsigned long updatePeriod, bool comms>
        InOutLetCycle<updatePeriod, comms>::~InOutLetCycle()
        {

        }

        template<unsigned long updatePeriod, bool comms>
        void InOutLetCycle<updatePeriod, comms>::UpdateCycle(std::vector<distribn_t> &densityCycle,
                                                             const SimulationState *state)
        {
          log::Logger::Log<log::Debug, log::OnePerCore>("Update cycle for iolet %d %d",updatePeriod,state->Get0IndexedTimeStep());
          if ( (updatePeriod != 0 && state->Get0IndexedTimeStep() % updatePeriod == 0)
              || (updatePeriod == 0 && state->Get0IndexedTimeStep() == 0))
          {
            log::Logger::Log<log::Debug, log::OnePerCore>("Calculating iolet cycle");
            CalculateCycle(densityCycle, state);
          }
        }

        template<unsigned long updatePeriod, bool comms>
        bool InOutLetCycle<updatePeriod, comms>::DoComms()
        {
          return comms;
        }
      }
    }
  }
}

#endif /* HEMELB_LB_BOUNDARIES_INOUTLETCYCLE_H */
