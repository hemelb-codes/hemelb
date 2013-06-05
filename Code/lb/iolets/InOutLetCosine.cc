// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#include "lb/iolets/InOutLetCosine.h"
#include "configuration/SimConfig.h"
#include "topology/NetworkTopology.h"

namespace hemelb
{
  namespace lb
  {
    namespace iolets
    {

      InOutLetCosine::InOutLetCosine() :
        InOutLet(), pressureMeanPhysical(0.0), pressureAmpPhysical(0.0), phase(0.0), period(1.0),
            warmUpLength(0)
      {

      }

      void InOutLetCosine::DoIO(TiXmlElement *iParent, bool iIsLoading,
                                configuration::SimConfig* iSimConfig)
      {
        iSimConfig->DoIOForCosineInOutlet(iParent, iIsLoading, this);
      }

      InOutLet* InOutLetCosine::Clone() const
      {
        InOutLetCosine* copy = new InOutLetCosine(*this);

        return copy;
      }

      InOutLetCosine::~InOutLetCosine()
      {

      }

      LatticeDensity InOutLetCosine::GetDensity(unsigned long time_step) const
      {
        double w = 2.0 * PI / period;

        // Calculate the target
        LatticeDensity target = GetDensityMean() + GetDensityAmp() * cos(w
            * units->ConvertTimeStepToPhysicalUnits(time_step) + phase);

        // If we're past the warm-up phase, just use the target.
        if (time_step >= warmUpLength)
        {
          return target;
        }

        // Otherwise, interpolate between the minimum in the simulation and the target.
        double interpolationFactor = ((double) time_step) / ((double) warmUpLength);

        return interpolationFactor * target + (1. - interpolationFactor) * minimumSimulationDensity;
      }

      LatticeDensity InOutLetCosine::GetDensityMean() const
      {
        return units->ConvertPressureToLatticeUnits(pressureMeanPhysical) / Cs2;
      }

      LatticeDensity InOutLetCosine::GetDensityAmp() const
      {
        return units->ConvertPressureDifferenceToLatticeUnits(pressureAmpPhysical) / Cs2;
      }

    }
  }
}
