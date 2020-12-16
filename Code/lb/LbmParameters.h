// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_LBMPARAMETERS_H
#define HEMELB_LB_LBMPARAMETERS_H

#include <cmath>
#include "constants.h"
#include <vector>
#include <cassert>

namespace hemelb
{
  namespace lb
  {
    enum StressTypes
    {
      VonMises = 0,
      ShearStress = 1,
      IgnoreStress = 2
    };

    class LbmParameters
    {
      public:
        inline LbmParameters(PhysicalTime timeStepSeconds,
		    PhysicalDistance voxelSizeMetres,
		    PhysicalDensity fluidDensityKgm3 = DEFAULT_FLUID_DENSITY_Kg_per_m3) :
	  timeStep(timeStepSeconds),
	  voxelSize(voxelSizeMetres),
	  fluidDensity(fluidDensityKgm3),
	  tau(0.5 + (timeStep * BLOOD_VISCOSITY_Pa_s / fluidDensity) / (Cs2 * voxelSize * voxelSize)),
	  omega(-1.0 / tau),
	  stressParameter((1.0 - 1.0 / (2.0 * tau)) / std::sqrt(2.0)),
	  beta(-1.0 / (2.0 * tau))
        {
        }

        inline const PhysicalTime& GetTimeStep() const
        {
          return timeStep;
        }

        inline const PhysicalDistance& GetVoxelSize() const
        {
          return voxelSize;
        }

        inline const PhysicalDensity& GetFluidDensity() const
        {
          return fluidDensity;
        }

        inline const distribn_t& GetOmega() const
        {
          return omega;
        }

        inline const distribn_t& GetTau() const
        {
          return tau;
        }

        inline const distribn_t& GetStressParameter() const
        {
          return stressParameter;
        }

        inline const distribn_t& GetBeta() const
        {
          return beta;
        }

        StressTypes StressType;

      private:
        PhysicalTime timeStep;
        PhysicalDistance voxelSize;
        PhysicalDensity fluidDensity;
        distribn_t tau;
        distribn_t omega;
        distribn_t stressParameter;
        distribn_t beta; ///< Viscous dissipation in ELBM
    };
  }
}

#endif //HEMELB_LB_LBMPARAMETERS_H
