// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

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
         * @tparam updatePeriod densityBuffer is recalculated every updatePeriod time steps. If
         *                 it is 0 then the densityBuffer is only initialised at the begining
         *                 and updated only when simulation is reset
         * @tparam comms if true BCproc updates values and sends them out. If false relevant procs update
         *          locally.
         */
        template<unsigned long updatePeriod, bool comms>
        class InOutLetCycle : public InOutLet
        {
          public:
            //@override
            virtual bool IsCommsRequired() const
            {
              return comms;
            }
            LatticeDensity GetDensity(LatticeTime time_step) const
            {
              if (!IsCommsRequired())
              {
                // logically constant method with a non-const implementation.
                UpdateCycle(time_step);
              }
              return densityBuffer[time_step % updatePeriod];
            }

            virtual void Reset(SimulationState &state)
            {
              densityBuffer.resize(state.GetTotalTimeSteps());
            }
            void DoComms(bool is_io_proc, LatticeTime time_step);
          protected:
            InOutLetCycle() :
                InOutLet(), densityBuffer()
            {
            }
            virtual ~InOutLetCycle()
            {
              delete comms;
            }
            /***
             * Fill in a vector of values of the density for this IOLet.
             * The Iolet will not subsequently calculate it's density values, but will read them from this buffer.
             * This array is NOT typically one value for each point in the pulsatile cycle.
             * The length of the densityCycle argument is an arbitrary choice of how often to calculate the densities for this and some subsequent steps.
             * This length is usually one, indicating that the calculation is done every time step,
             * or zero, indicating that the calculation is done once for the whole simulation.
             */
            virtual void CalculateCycle() const = 0; // fill the density buffer

          private:
            void UpdateCycle(LatticeTime time_step) const
            {
              if (ShouldUpdateBuffer(time_step))
              {
                CalculateCycle();
              }
            }
            bool ShouldUpdateBuffer(LatticeTime time_step) const
            {
              if (updatePeriod == 0)
              {
                return time_step == 0;
              }
              return time_step % updatePeriod == 0;
            }
            // a cache, logically constant methods can modify this.
            // (Perfect use-case for mutable.)
            mutable std::vector<LatticeDensity> densityBuffer;
        }
        ;
      }
    }
  }
}
#endif /* HEMELB_LB_BOUNDARIES_INOUTLETCYCLE_H */
