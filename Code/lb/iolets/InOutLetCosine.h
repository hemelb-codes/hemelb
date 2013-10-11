// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_LB_IOLETS_INOUTLETCOSINE_H
#define HEMELB_LB_IOLETS_INOUTLETCOSINE_H

#include "lb/iolets/InOutLet.h"

namespace hemelb
{
  namespace lb
  {
    namespace iolets
    {

      /*
       * Template values are chosen to be tUpdatePeriod = 1 and tComms = false
       * The cosine pressure trace can be easily calculated locally by any proc and
       * there is no need to store values for any time step beyond the current one
       */
      class InOutLetCosine : public InOutLet
      {
        public:
          InOutLetCosine();
          virtual ~InOutLetCosine();
          virtual InOutLet* Clone() const;
          virtual void Reset(SimulationState &state)
          {
            //pass;
          }
          LatticeDensity GetDensityMean() const;
          LatticeDensity GetDensityAmp() const;
          LatticeDensity GetDensity(unsigned long time_step) const;
          PhysicalPressure GetPressureMin() const
          {
            return pressureMeanPhysical - pressureAmpPhysical;
          }
          PhysicalPressure GetPressureMax() const
          {
            return pressureMeanPhysical + pressureAmpPhysical;
          }
          // TODO I do not like returning references, this method should be const and we should have a setter.
          // but, the way the IO code in SimConfig is currently set up prevents this for now.
          PhysicalPressure & GetPressureMean()
          {
            return pressureMeanPhysical;
          }
          void SetPressureMean(PhysicalPressure pressure)
          {
            pressureMeanPhysical = pressure;
          }
          PhysicalPressure & GetPressureAmp()
          {
            return pressureAmpPhysical;
          }
          void SetPressureAmp(PhysicalPressure pressure)
          {
            pressureAmpPhysical = pressure;
          }
          double & GetPhase()
          {
            return phase;
          }
          void SetPhase(double aPhase)
          {
            phase = aPhase;
          }
          PhysicalTime & GetPeriod()
          {
            return period;
          }
          void SetPeriod(PhysicalTime aPeriod)
          {
            period = aPeriod;
          }
          void SetWarmup(unsigned int warmup)
          {
            warmUpLength = warmup;
          }
        private:

          PhysicalPressure pressureMeanPhysical;
          PhysicalPressure pressureAmpPhysical;

          double phase;
          PhysicalTime period;
          unsigned int warmUpLength;
      };

    }
  }
}

#endif /* HEMELB_LB_IOLETS_INOUTLETCOSINE_H */
