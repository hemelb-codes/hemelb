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
         * Generic base InOutLet with varying update period and communications strategy.
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

            virtual bool GetIsCommsRequired()
            {
              return comms;
            }
            LatticeDensity GetDensity(unsigned long time_step)
            {
              if (!GetIsCommsRequired())
              {
                UpdateCycle(time_step);
              }
              return densityBuffer[time_step % updatePeriod];
            }

            virtual void Reset(SimulationState &state)
            {
              densityBuffer.resize(state.GetTotalTimeSteps());
            }
            virtual void DoComms(bool is_io_proc, unsigned long time_step);
          protected:
            InOutLetCycle() :
                InOutLet(), densityBuffer()
            {
            }
            virtual ~InOutLetCycle()
            {
              delete comms;
            }
            ;
            /***
             * Fill in a vector of values of the density for this IOLet.
             * The Iolet will not subsequently calculate it's density values, but will read them from this buffer.
             * This array is NOT typically one value for each point in the pulsatile cycle.
             * The length of the densityCycle argument is an arbitrary choice of how often to calculate the densities for this and some subsequent steps.
             * This length is usually one, indicating that the calculation is done every time step,
             * or zero, indicating that the calculation is done once for the whole simulation.
             * @param densityCycle An array of densities for this iolet
             * @param iState Simulation state for the iolet
             */
            virtual void CalculateCycle() = 0; // fill the density buffer

          private:
            void UpdateCycle(unsigned long time_step)
            {
              if (shouldUpdateBuffer(time_step))
              {
                CalculateCycle();
              }
            }
            bool shouldUpdateBuffer(unsigned int time_step)
            {
              if (updatePeriod == 0)
              {
                return time_step == 0;
              }
              return time_step % updatePeriod == 0;
            }
            std::vector<LatticeDensity> densityBuffer;
        }
        ;
      }
    }
  }
}
#endif /* HEMELB_LB_BOUNDARIES_INOUTLETCYCLE_H */
