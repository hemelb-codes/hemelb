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

          LatticeDensity GetDensity(unsigned long time_step) const;

          const LatticeDensity& GetDensityMean() const
          {
            return densityMean;
          }
          void SetDensityMean(const LatticeDensity& rho)
          {
            densityMean = rho;
          }

          const LatticeDensity& GetDensityAmp() const
          {
            return densityAmp;
          }
          void SetDensityAmp(const LatticeDensity& rho)
          {
            densityAmp = rho;
          }

          LatticeDensity GetDensityMin() const
          {
            return (densityMean - densityAmp);
          }
          LatticeDensity GetDensityMax() const
          {
            return (densityMean + densityAmp);
          }

          LatticePressure GetPressureMean() const
          {
            return densityMean * Cs2;
          }
          void SetPressureMean(const LatticePressure& pressure)
          {
            densityMean = pressure / Cs2;
          }

          LatticePressure GetPressureAmp() const
          {
            return densityAmp * Cs2;
          }
          void SetPressureAmp(const LatticePressure& pressure)
          {
            densityAmp = pressure / Cs2;
          }

          const Angle& GetPhase() const
          {
            return phase;
          }
          void SetPhase(const Angle& aPhase)
          {
            phase = aPhase;
          }

          const LatticeTime& GetPeriod() const
          {
            return period;
          }
          void SetPeriod(const LatticeTime& aPeriod)
          {
            period = aPeriod;
          }

          void SetWarmup(unsigned int warmup)
          {
            warmUpLength = warmup;
          }
        private:

          LatticeDensity densityMean;
          LatticeDensity densityAmp;

          Angle phase;
          LatticeTime period;
          unsigned int warmUpLength;
      };

    }
  }
}

#endif /* HEMELB_LB_IOLETS_INOUTLETCOSINE_H */
